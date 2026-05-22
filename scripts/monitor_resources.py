#!/usr/bin/env python3
"""Sample CPU / RAM / GPU / disk usage during an ssign run.

Standalone, no ssign imports. Run alongside an ssign invocation, point
``--log-file`` at the tee'd stdout, and the monitor will tag each sample
with the current pipeline step.

GPU sampling is single-device (index 0); fine for Imperial CX3 (one A40
per node) but will under-report on a multi-GPU host.

Install (CX3 venv):
    pip install psutil nvidia-ml-py
(`nvidia-ml-py` is the maintained replacement for the inactive `pynvml`
PyPI package; both expose the same `pynvml` module.)

The `step` column is populated by tailing the log for the line emitted
at `core/runner.py` (search for "Starting step"). If that format string
changes, this regex needs to follow.

Example:
    python scripts/monitor_resources.py \\
        --out resources.csv --log-file run.log --interval 5 &
    MON_PID=$!
    ssign run input.gbff --outdir results 2>&1 | tee run.log
    kill $MON_PID
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import signal
import sys
import time
from pathlib import Path

import psutil

_GPU = None
_GPU_TOTAL_GB = 0.0
try:
    import pynvml

    pynvml.nvmlInit()
    _GPU = pynvml.nvmlDeviceGetHandleByIndex(0)
    _GPU_TOTAL_GB = pynvml.nvmlDeviceGetMemoryInfo(_GPU).total / 2**30
except (ImportError, Exception) as _gpu_err:  # NVMLError lives on pynvml
    if not isinstance(_gpu_err, ImportError):
        print(f"[monitor] GPU sampling disabled: {_gpu_err}", file=sys.stderr)

STEP_RE = re.compile(r"Starting step (\d+)\s*/\s*(\d+)\s*:\s*(.+?)(?:\s*\(|$)")


def gpu_sample() -> tuple[str, str]:
    if _GPU is None:
        return ("", "")
    util = pynvml.nvmlDeviceGetUtilizationRates(_GPU).gpu
    used_gb = pynvml.nvmlDeviceGetMemoryInfo(_GPU).used / 2**30
    return (f"{util:.0f}", f"{used_gb:.2f}")


def tail_step(log_path: Path | None, state: dict) -> str:
    if log_path is None or not log_path.exists():
        return state.get("step", "")
    size = log_path.stat().st_size
    if size < state.get("pos", 0):  # log rotated/truncated
        state["pos"] = 0
    with log_path.open("rb") as fh:
        fh.seek(state["pos"])
        chunk = fh.read()
        state["pos"] = fh.tell()
    for line in chunk.decode(errors="replace").splitlines():
        m = STEP_RE.search(line)
        if m:
            state["step"] = f"{m.group(1)}/{m.group(2)}:{m.group(3).strip()}"
    return state.get("step", "")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="Output CSV path")
    ap.add_argument("--log-file", help="ssign stdout log to tail for step boundaries")
    ap.add_argument("--interval", type=float, default=5.0, help="Seconds between samples (default 5)")
    args = ap.parse_args()

    log_path = Path(args.log_file) if args.log_file else None
    state: dict = {"pos": 0, "step": ""}

    psutil.cpu_percent(percpu=True)  # prime the per-call delta
    disk_prev = psutil.disk_io_counters()
    t_prev = time.monotonic()

    stop = {"flag": False}
    signal.signal(signal.SIGTERM, lambda *_: stop.update(flag=True))
    signal.signal(signal.SIGINT, lambda *_: stop.update(flag=True))

    print(
        f"[monitor] pid={os.getpid()} writing to {args.out} every {args.interval}s (gpu={'yes' if _GPU else 'no'})",
        file=sys.stderr,
    )

    ram_total_gb = psutil.virtual_memory().total / 2**30

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "timestamp",
                "cpu_overall_pct",
                "cpu_per_core",
                "ram_used_gb",
                "ram_total_gb",
                "gpu_util_pct",
                "gpu_mem_gb",
                "gpu_mem_total_gb",
                "disk_read_mb_s",
                "disk_write_mb_s",
                "step",
            ]
        )
        fh.flush()
        while not stop["flag"]:
            time.sleep(args.interval)
            now = time.monotonic()
            dt = max(now - t_prev, 1e-6)
            per_core = psutil.cpu_percent(percpu=True)
            overall = sum(per_core) / len(per_core) if per_core else 0.0
            ram_used_gb = psutil.virtual_memory().used / 2**30
            disk_now = psutil.disk_io_counters()
            r_mb = (disk_now.read_bytes - disk_prev.read_bytes) / 2**20 / dt
            wr_mb = (disk_now.write_bytes - disk_prev.write_bytes) / 2**20 / dt
            disk_prev, t_prev = disk_now, now
            gu, gm = gpu_sample()
            step = tail_step(log_path, state)
            w.writerow(
                [
                    time.strftime("%Y-%m-%dT%H:%M:%S"),
                    f"{overall:.1f}",
                    ";".join(f"{c:.0f}" for c in per_core),
                    f"{ram_used_gb:.2f}",
                    f"{ram_total_gb:.2f}",
                    gu,
                    gm,
                    f"{_GPU_TOTAL_GB:.2f}" if _GPU is not None else "",
                    f"{r_mb:.1f}",
                    f"{wr_mb:.1f}",
                    step,
                ]
            )
            fh.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
