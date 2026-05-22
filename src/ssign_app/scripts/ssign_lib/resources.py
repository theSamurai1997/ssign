"""Cgroup/scheduler-aware CPU and RAM detection.

`os.cpu_count()` and `psutil.virtual_memory().total` read host-wide values
and ignore the job's actual allocation under PBS/SLURM cgroups. On a shared
HPC node this leads to massive thread oversubscription (24 CPUs requested
on a 4-CPU allocation thrashes at ~2% efficiency) and OOM-by-default for
memory-heavy steps. These helpers query the real allocation.

Lookup priority for both is: explicit env override → scheduler env vars →
scheduler CLI (qstat/scontrol) → cgroup v2 → cgroup v1 → host total.
Take the smallest reasonable value, since under-using cores is far cheaper
than thrashing them.
"""

from __future__ import annotations

import logging
import os
import re
import subprocess

logger = logging.getLogger(__name__)

# cgroup v1 "no limit" sentinel: math.inf, encoded as 2**63 - 1 rounded
# down to the kernel page size. Anything within ~1 page of this is the
# "unlimited" marker, not a real limit.
_CGROUP_UNLIMITED_BYTES = 2**63 - 4096


def effective_cpu_count() -> int:
    """Number of CPUs the current process is allowed to schedule on.

    Honors cgroup CPU pinning (PBS cpuset, SLURM `--cpus-per-task`). Falls
    back to host total when `sched_getaffinity` isn't available (non-Linux).
    """
    if hasattr(os, "sched_getaffinity"):
        try:
            return max(1, len(os.sched_getaffinity(0)))
        except OSError:
            pass
    return os.cpu_count() or 4


def _parse_size_to_gb(value: str) -> float | None:
    """Parse a size string like '32gb', '128GiB', '4096mb' → GB (decimal)."""
    m = re.match(r"\s*([\d.]+)\s*([kmgt]i?b?)?\s*$", value.lower())
    if not m:
        return None
    n = float(m.group(1))
    unit = (m.group(2) or "").rstrip("b")
    mult = {
        "": 1 / 2**30,
        "k": 1 / 2**20,
        "ki": 1 / 2**20,
        "m": 1 / 1024,
        "mi": 1 / 1024,
        "g": 1.0,
        "gi": 1.0,
        "t": 1024.0,
        "ti": 1024.0,
    }.get(unit)
    return n * mult if mult else None


def _pbs_alloc_gb() -> float | None:
    """Parse Resource_List.mem from `qstat -f $PBS_JOBID`."""
    job = os.environ.get("PBS_JOBID")
    if not job:
        return None
    try:
        out = subprocess.run(["qstat", "-f", job], capture_output=True, text=True, timeout=5)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return None
    if out.returncode != 0:  # job vanished, auth error, etc.
        return None
    m = re.search(r"Resource_List\.mem\s*=\s*(\S+)", out.stdout)
    return _parse_size_to_gb(m.group(1)) if m else None


def _cgroup_mem_gb() -> float | None:
    """Read the memory cgroup limit, skipping the 'unlimited' sentinels."""
    # cgroup v2 reports the literal string "max" when no limit is set.
    for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            with open(path) as fh:
                raw = fh.read().strip()
        except OSError:
            continue
        if raw == "max" or not raw.isdigit():
            continue
        n = int(raw)
        if n >= _CGROUP_UNLIMITED_BYTES:
            continue
        return n / 2**30
    return None


def effective_ram_gb() -> float:
    """Smallest of every job-allocation signal we can find, falling back to host.

    Sources tried (returns the minimum of all that produce a positive value):
    `SSIGN_MAX_RAM_GB` env → `SLURM_MEM_PER_NODE` env → PBS qstat parse →
    cgroup v2/v1 limit → psutil host total. The min handles HPCs that
    enforce RAM externally (PBS) while the kernel cgroup still reports
    "unlimited" — the qstat parse catches what /sys misses.
    """
    candidates: list[float] = []

    override = os.environ.get("SSIGN_MAX_RAM_GB")
    if override:
        try:
            candidates.append(float(override))
        except ValueError:
            logger.warning(f"SSIGN_MAX_RAM_GB={override!r} not a number; ignored")

    slurm_mb = os.environ.get("SLURM_MEM_PER_NODE")
    if slurm_mb and slurm_mb.isdigit():
        candidates.append(int(slurm_mb) / 1024)

    pbs = _pbs_alloc_gb()
    if pbs is not None:
        candidates.append(pbs)

    cg = _cgroup_mem_gb()
    if cg is not None:
        candidates.append(cg)

    try:
        import psutil

        candidates.append(psutil.virtual_memory().total / 2**30)
    except Exception:
        pass

    return min(c for c in candidates if c > 0) if candidates else 0.0
