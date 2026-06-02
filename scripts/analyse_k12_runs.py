#!/usr/bin/env python3
"""Post-run analysis for the K-12 dual-run validation.

Compares two ssign runs of the same genome (e.g. RTX6000 interactive vs L40S
PBS) to:

1. Confirm both finished and produced the expected outputs.
2. Diff predictions to check determinism across GPU models.
3. Extract per-tool wallclocks from the runner stdout log and emit
   JSON-Lines rows compatible with calibration/runs.jsonl.
4. Pull PLM-Effector configuration (batch_size, dtype, per-type timings)
   from streamed subprocess output, when present in the log.

Outputs Markdown summary on stdout, JSONL rows to a file when
``--calibration-out`` is given.

Run after both ssign runs complete:

    scripts/analyse_k12_runs.py \\
        --run-a results-a-dir --machine-a CX3-RTX6000 --log-a runA.log \\
        --run-b results-b-dir --machine-b CX3-L40S    --log-b runB.log \\
        --calibration-out /tmp/k12_runs.jsonl
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

# Runner stdout emits "— " (U+2014) separators. Match U+2014 / U+2013 / "--"
# so the parser survives a future runner change to en-dash or ASCII fallback.
_DASH = r"(?:—|–|--)"

# Lines look like:
#   [ 10%] Extracting proteins — Done: Bakta (re-annotated): ... | 4m 8s elapsed
#   [ 30%] Running in parallel: ... — 3 tools running simultaneously | 4m 40s elapsed
_STEP_DONE_RE = re.compile(
    rf"\[\s*(?P<pct>\d+)%\]\s+(?P<label>.+?)\s+{_DASH}\s+Done:\s+(?P<msg>.+?)\s+\|\s+(?P<elapsed>[\dhms\s]+)\s+elapsed",
)
_PARALLEL_START_RE = re.compile(
    rf"\[\s*(?P<pct>\d+)%\]\s+Running in parallel:\s+(?P<label>.+?)\s+{_DASH}\s+\d+\s+tools running simultaneously\s+\|\s+(?P<elapsed>[\dhms\s]+)\s+elapsed",
)

# Streamed PLM-E subprocess lines (only present with the runner's stream_stderr=True path):
#   [run_plm_effector.py] PLM-Effector: extracting features for ... (batch_size=64, chunk_size=256, dtype=bf16, PLMs=...)
#   [run_plm_effector.py] PLM-Effector: T1SE — wrote 4314 predictions (3 passing threshold) to ...
_PLME_CONFIG_RE = re.compile(
    r"\[run_plm_effector\.py\]\s+PLM-Effector:\s+extracting features.+?batch_size=(?P<bs>\d+),\s+chunk_size=(?P<cs>\d+),\s+dtype=(?P<dt>\w+)"
)
_PLME_TYPE_DONE_RE = re.compile(
    rf"\[run_plm_effector\.py\]\s+PLM-Effector:\s+(?P<type>T[1-9][a-zA-Z]?SE)\s+{_DASH}\s+wrote\s+(?P<n>\d+)\s+predictions\s+\((?P<pos>\d+)\s+passing threshold\)"
)

# One row per (unit, multiplier) so we don't rebuild this on every call.
_ELAPSED_PART_RE = re.compile(r"(\d+)\s*([hms])")
_ELAPSED_UNIT_S = {"h": 3600, "m": 60, "s": 1}

# Friendly step label → canonical tool name in calibration/runs.jsonl.
# Steps that don't map to a single tool (or are bookkeeping) are dropped.
_LABEL_TO_TOOL = {
    "Extracting proteins": "bakta",
    "Running MacSyFinder": "macsyfinder",
    "Predicting localization (DeepLocPro)": "deeplocpro",
    "Predicting secretion type (DeepSecE)": "deepsece",
    "Predicting signal peptides (SignalP)": "signalp",
    "Predicting effectors (PLM-Effector, 5 types)": "plm_effector",
    "Running InterProScan": "interproscan",
    "Running BLASTp": "blastp",
    "Running HH-suite": "hh_suite",
    "Running EggNOG-mapper": "eggnog",
    "Running pLM-BLAST": "plm_blast",
    "Computing ProtParam": "protparam",
}


@dataclass
class StepRecord:
    label: str
    pct: int
    elapsed_s: int
    message: str
    parallel_start_s: int | None = None


@dataclass
class RunSummary:
    machine: str
    run_dir: Path
    log_path: Path | None
    steps: list[StepRecord] = field(default_factory=list)
    plme_batch_size: int | None = None
    plme_chunk_size: int | None = None
    plme_dtype: str | None = None
    plme_type_counts: dict[str, tuple[int, int]] = field(default_factory=dict)
    final_elapsed_s: int | None = None


def parse_elapsed(token: str) -> int:
    """Convert "4m 8s" / "1h 23m 45s" / "19s" to total seconds."""
    return sum(int(n) * _ELAPSED_UNIT_S[u] for n, u in _ELAPSED_PART_RE.findall(token))


def parse_log(path: Path) -> tuple[list[StepRecord], dict]:
    """Walk a runner stdout log and pull step-done lines + PLM-E config.

    Returns (step_records, plme_meta). `plme_meta` keys: batch_size, chunk_size,
    dtype, type_counts (dict of {effector_type: (total, positive)}).
    """
    steps: list[StepRecord] = []
    plme_meta: dict = {"batch_size": None, "chunk_size": None, "dtype": None, "type_counts": {}}
    parallel_start_s: int | None = None
    with path.open() as f:
        for raw in f:
            line = raw.rstrip("\n")
            m = _PARALLEL_START_RE.search(line)
            if m:
                parallel_start_s = parse_elapsed(m.group("elapsed"))
                continue
            m = _STEP_DONE_RE.search(line)
            if m:
                label = m.group("label").strip()
                elapsed_s = parse_elapsed(m.group("elapsed"))
                steps.append(
                    StepRecord(
                        label=label,
                        pct=int(m.group("pct")),
                        elapsed_s=elapsed_s,
                        message=m.group("msg").strip(),
                        parallel_start_s=parallel_start_s if _is_parallel_member(label) else None,
                    )
                )
                continue
            m = _PLME_CONFIG_RE.search(line)
            if m:
                plme_meta["batch_size"] = int(m.group("bs"))
                plme_meta["chunk_size"] = int(m.group("cs"))
                plme_meta["dtype"] = m.group("dt")
                continue
            m = _PLME_TYPE_DONE_RE.search(line)
            if m:
                plme_meta["type_counts"][m.group("type")] = (int(m.group("n")), int(m.group("pos")))
    return steps, plme_meta


def _is_parallel_member(label: str) -> bool:
    """True for the three predictor steps that run as a parallel group.

    The parallel group is fixed in `_step_predict_parallel`; if it ever
    grows, extend this set rather than parsing the announcement line.
    """
    return label in {
        "Predicting localization (DeepLocPro)",
        "Predicting secretion type (DeepSecE)",
        "Predicting signal peptides (SignalP)",
    }


def compute_durations(steps: list[StepRecord]) -> dict[str, int]:
    """Per-step wallclock from cumulative-elapsed values.

    Sequential steps: end_n - end_{n-1} where prev_end is the previous
    non-parallel-member step (so parallel members don't all attribute the
    wallclock of the parallel block to the one that finished first).
    Parallel members: end - parallel_group_start_elapsed.
    """
    durations: dict[str, int] = {}
    prev_sequential_end: int = 0
    for step in steps:
        if step.parallel_start_s is not None:
            dur = max(0, step.elapsed_s - step.parallel_start_s)
        else:
            dur = max(0, step.elapsed_s - prev_sequential_end)
            prev_sequential_end = step.elapsed_s
        durations[step.label] = dur
    return durations


def load_run(run_dir: Path, log_path: Path | None, machine: str) -> RunSummary:
    summary = RunSummary(machine=machine, run_dir=run_dir, log_path=log_path)
    if log_path and log_path.exists():
        summary.steps, plme = parse_log(log_path)
        summary.plme_batch_size = plme["batch_size"]
        summary.plme_chunk_size = plme["chunk_size"]
        summary.plme_dtype = plme["dtype"]
        summary.plme_type_counts = plme["type_counts"]
        if summary.steps:
            summary.final_elapsed_s = summary.steps[-1].elapsed_s
    return summary


def compare_results(run_a_dir: Path, run_b_dir: Path, sample_id: str = "ecoli_k12") -> dict:
    """Diff the two runs' final results.csv on locus_tag + key prediction columns.

    `keep_default_na=False` keeps strings like "no_hits" / "None" / "NA" as the
    literal strings the consensus script wrote, instead of letting pandas
    coerce them to NaN and changing equality semantics downstream.
    """
    import pandas as pd

    out: dict = {}
    path_a = run_a_dir / f"{sample_id}_results.csv"
    path_b = run_b_dir / f"{sample_id}_results.csv"
    out["path_a"] = str(path_a)
    out["path_b"] = str(path_b)
    out["exists_a"] = path_a.exists()
    out["exists_b"] = path_b.exists()
    if not (path_a.exists() and path_b.exists()):
        return out

    df_a = pd.read_csv(path_a, keep_default_na=False, na_values=[""])
    df_b = pd.read_csv(path_b, keep_default_na=False, na_values=[""])
    out["shape_a"] = df_a.shape
    out["shape_b"] = df_b.shape
    out["columns_match"] = list(df_a.columns) == list(df_b.columns)

    if "locus_tag" in df_a.columns and "locus_tag" in df_b.columns:
        # Duplicate locus_tags would Cartesian-explode the merges below and
        # silently inflate the {col}_total counts. Surface as a warning, not
        # an assert, so the comparison still produces partial info.
        if not df_a["locus_tag"].is_unique:
            out["warning_a_duplicate_locus_tags"] = True
        if not df_b["locus_tag"].is_unique:
            out["warning_b_duplicate_locus_tags"] = True

        set_a = set(df_a["locus_tag"])
        set_b = set(df_b["locus_tag"])
        common = set_a & set_b
        out["n_locus_common"] = len(common)
        out["n_locus_only_a"] = len(set_a - common)
        out["n_locus_only_b"] = len(set_b - common)

        for col in ("broad_annotation", "confidence_tier", "ss_type"):
            if col in df_a.columns and col in df_b.columns:
                joined = df_a[["locus_tag", col]].merge(
                    df_b[["locus_tag", col]],
                    on="locus_tag",
                    suffixes=("_a", "_b"),
                )
                if joined.empty:
                    out[f"{col}_match"] = None
                    continue
                same = (joined[f"{col}_a"] == joined[f"{col}_b"]) | (
                    joined[f"{col}_a"].isna() & joined[f"{col}_b"].isna()
                )
                out[f"{col}_match"] = int(same.sum())
                out[f"{col}_total"] = len(joined)
    return out


def calibration_rows(
    summary: RunSummary, tier: str = "extended", genome: str = "ecoli_k12", n_proteins: int | None = None
) -> list[dict]:
    """Map parsed step durations to runs.jsonl rows for the canonical tools."""
    durations = compute_durations(summary.steps)
    today = datetime.date.today().isoformat()
    run_id = f"{today}-{summary.machine.lower()}-k12-{tier}"
    base_input = {"genome": genome, "tier": tier, "use_input_annotations": False}
    if n_proteins is not None:
        base_input["n_proteins"] = n_proteins

    rows: list[dict] = []
    for label, dur in durations.items():
        tool = _LABEL_TO_TOOL.get(label)
        if tool is None:
            continue
        row = {
            "run_id": run_id,
            "timestamp": today,
            "machine": summary.machine,
            "input": dict(base_input),
            "tool": tool,
            "wallclock_s": dur,
            "success": True,
            "notes": "",
        }
        if tool == "plm_effector" and summary.plme_batch_size:
            row["notes"] = (
                f"batch_size={summary.plme_batch_size}, dtype={summary.plme_dtype}, "
                f"types={','.join(summary.plme_type_counts.keys()) or 'unknown'}"
            )
        rows.append(row)
    if summary.final_elapsed_s is not None:
        rows.append(
            {
                "run_id": run_id,
                "timestamp": today,
                "machine": summary.machine,
                "input": dict(base_input),
                "tool": "_pipeline_total",
                "wallclock_s": summary.final_elapsed_s,
                "success": True,
                "notes": "",
            }
        )
    return rows


def render_markdown(summary_a: RunSummary, summary_b: RunSummary, compare: dict) -> str:
    """Render side-by-side comparison as Markdown."""
    durs_a = compute_durations(summary_a.steps)
    durs_b = compute_durations(summary_b.steps)
    all_labels = list(durs_a.keys()) + [k for k in durs_b.keys() if k not in durs_a]

    lines = []
    lines.append(f"# K-12 dual-run analysis ({datetime.date.today().isoformat()})")
    lines.append("")
    lines.append("## Run summaries")
    lines.append("")
    lines.append(f"| | {summary_a.machine} | {summary_b.machine} |")
    lines.append("| --- | --- | --- |")
    lines.append(f"| run dir | `{summary_a.run_dir}` | `{summary_b.run_dir}` |")
    lines.append(f"| log    | `{summary_a.log_path or 'n/a'}` | `{summary_b.log_path or 'n/a'}` |")
    lines.append(
        f"| total wallclock | {_fmt_secs(summary_a.final_elapsed_s)} | {_fmt_secs(summary_b.final_elapsed_s)} |"
    )
    lines.append(f"| PLM-E batch_size | {summary_a.plme_batch_size or '—'} | {summary_b.plme_batch_size or '—'} |")
    lines.append(f"| PLM-E dtype | {summary_a.plme_dtype or '—'} | {summary_b.plme_dtype or '—'} |")
    lines.append("")
    lines.append("## Per-tool wallclock (seconds)")
    lines.append("")
    lines.append(f"| step | {summary_a.machine} | {summary_b.machine} | speedup |")
    lines.append("| --- | ---: | ---: | ---: |")
    for label in all_labels:
        a = durs_a.get(label)
        b = durs_b.get(label)
        ratio = f"{a / b:.2f}×" if (a and b) else "—"
        lines.append(f"| {label} | {a if a is not None else '—'} | {b if b is not None else '—'} | {ratio} |")
    lines.append("")
    lines.append("## Results agreement")
    lines.append("")
    if compare.get("exists_a") and compare.get("exists_b"):
        lines.append(f"- Shapes: {compare['shape_a']} vs {compare['shape_b']}")
        lines.append(f"- Columns identical: {compare['columns_match']}")
        if "n_locus_common" in compare:
            lines.append(
                f"- locus_tag overlap: {compare['n_locus_common']} common, "
                f"{compare['n_locus_only_a']} only-A, {compare['n_locus_only_b']} only-B"
            )
            for col in ("broad_annotation", "confidence_tier", "ss_type"):
                if f"{col}_match" in compare:
                    lines.append(f"- {col} matches: {compare[f'{col}_match']}/{compare[f'{col}_total']}")
    else:
        lines.append(f"- Missing result CSV: A={compare.get('exists_a')}, B={compare.get('exists_b')}")
    lines.append("")
    if summary_a.plme_type_counts or summary_b.plme_type_counts:
        lines.append("## PLM-Effector per-type predictions")
        lines.append("")
        lines.append(f"| type | {summary_a.machine} (n / passing) | {summary_b.machine} (n / passing) |")
        lines.append("| --- | --- | --- |")
        types = sorted(set(summary_a.plme_type_counts) | set(summary_b.plme_type_counts))
        for t in types:
            a = summary_a.plme_type_counts.get(t)
            b = summary_b.plme_type_counts.get(t)
            a_cell = f"{a[0]}/{a[1]}" if a else "—"
            b_cell = f"{b[0]}/{b[1]}" if b else "—"
            lines.append(f"| {t} | {a_cell} | {b_cell} |")
        lines.append("")
    return "\n".join(lines)


def _fmt_secs(s: int | None) -> str:
    if s is None:
        return "—"
    h, rem = divmod(s, 3600)
    m, ss = divmod(rem, 60)
    if h:
        return f"{h}h {m}m {ss}s"
    if m:
        return f"{m}m {ss}s"
    return f"{ss}s"


def main(argv: list | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-a", required=True, type=Path, help="results dir of run A")
    parser.add_argument("--log-a", required=True, type=Path, help="runner stdout log of run A")
    parser.add_argument("--machine-a", required=True, help="machine label for run A (e.g. CX3-RTX6000)")
    parser.add_argument("--run-b", required=True, type=Path, help="results dir of run B")
    parser.add_argument("--log-b", required=True, type=Path, help="runner stdout log of run B")
    parser.add_argument("--machine-b", required=True, help="machine label for run B (e.g. CX3-L40S)")
    parser.add_argument("--tier", default="extended", help="tier label (default: extended)")
    parser.add_argument("--genome", default="ecoli_k12")
    parser.add_argument(
        "--n-proteins", type=int, default=None, help="optional override of n_proteins for calibration rows"
    )
    parser.add_argument(
        "--calibration-out", type=Path, default=None, help="path to write JSONL rows; omit to skip emission"
    )
    args = parser.parse_args(argv)

    summary_a = load_run(args.run_a, args.log_a, args.machine_a)
    summary_b = load_run(args.run_b, args.log_b, args.machine_b)
    compare = compare_results(args.run_a, args.run_b, sample_id=args.genome)
    print(render_markdown(summary_a, summary_b, compare))

    if args.calibration_out is not None:
        rows = calibration_rows(
            summary_a, tier=args.tier, genome=args.genome, n_proteins=args.n_proteins
        ) + calibration_rows(summary_b, tier=args.tier, genome=args.genome, n_proteins=args.n_proteins)
        with args.calibration_out.open("w") as f:
            for row in rows:
                f.write(json.dumps(row) + "\n")
        print(f"\nWrote {len(rows)} calibration rows to {args.calibration_out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
