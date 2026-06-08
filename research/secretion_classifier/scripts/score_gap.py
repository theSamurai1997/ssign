#!/usr/bin/env python3
"""Score gap-validation runs against literature ground truth.

For each of the 6 validation genomes, two ssign runs were submitted:
one with the default proximity filter and one with the four
whole-genome flags enabled. This script joins each run's predicted
substrate list against `data/ground_truth/<genome>.tsv` (literature-
curated known effectors with `locus_tag` + `ss_type` columns) and
computes per-(genome, SS-type, condition) recall.

SS-type normalization: ground-truth and ssign use different label
flavors ("T3SS-1", "T6SS_H1" in ground truth; "T3SS_typeF", "T6SSi"
in ssign output). Both sides reduce to a canonical "T1SS" / "T2SS"
... "T6SS" prefix for the recall comparison. Instance-level recall
(SPI-1 vs SPI-2, H1 vs H2 vs H3) is a Tier-2 question and is NOT
scored here — this script measures only "did ssign find the
substrate" at the SS-type granularity.

Usage:
    score_gap.py --results-dir /path/to/runs/ \\
                 --ground-truth-dir research/secretion_classifier/data/ground_truth/ \\
                 --output gap_results.tsv

The results-dir is expected to contain subdirectories named like
`gap_<sample>_<proxy|whole>_<jobid>` (output of the
submit_gap_validation.sh PBS script). The script finds the
`<sample_id>_substrates_filtered.tsv` inside each run and parses it.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# Genomes we expect in the experiment. Keys must match the substrate-TSV
# `sample_id` (which ssign derives from the input filename stem).
GENOMES = [
    "legionella_pneumophila",
    "coxiella_rsa493",
    "salmonella_lt2",
    "yersinia_pestis_co92",
    "pseudomonas_pao1",
    "vibrio_cholerae_n16961",
]

CONDITIONS = ["proxy", "whole"]

# Reduce any SS-type label to a canonical "T<N>SS" prefix.
# Tolerates T5 subtype letters (T5aSS, T5bSS, T5cSS) and T6 lineage
# suffixes (T6SSi, T6SSii) used in MacSyFinder model names, plus the
# ground-truth annotations like T6SS_H1 / T3SS-1.
_SS_PREFIX = re.compile(r"(T[1-6])[A-Za-z0-9]*SS", re.IGNORECASE)


def normalize_ss_type(label: str) -> str | None:
    """Reduce 'T4SS_typeT' / 'T3SS-1' / 'T5aSS_PF03797' / 'T6SSi' to 'T<N>SS'.

    Returns None if no T1-T6 prefix is found.
    """
    if not label:
        return None
    m = _SS_PREFIX.match(label.strip())
    return f"{m.group(1).upper()}SS" if m else None


def load_ground_truth(gt_path: Path) -> dict[str, set[str]]:
    """Read ground_truth/<genome>.tsv. Return {canonical_ss_type: {locus_tag, ...}}."""
    by_type: dict[str, set[str]] = {}
    with gt_path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            ss = normalize_ss_type(row["ss_type"])
            tag = (row.get("locus_tag") or "").strip()
            if ss and tag and tag != "-":
                by_type.setdefault(ss, set()).add(tag)
    return by_type


def load_predicted_substrates(substrates_path: Path) -> dict[str, set[str]]:
    """Read ssign's <sample>_substrates_filtered.tsv. Return {canonical_ss_type: {locus_tag, ...}}.

    A single prediction can be associated with multiple SS types via the
    comma-joined `nearby_ss_types` column; we explode the row so that
    locus_tag X with nearby_ss_types='T3SS_typeF,T6SSi' counts toward
    recall for both T3SS and T6SS.
    """
    by_type: dict[str, set[str]] = {}
    with substrates_path.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tag = (row.get("locus_tag") or "").strip()
            if not tag:
                continue
            nearby = (row.get("nearby_ss_types") or "").strip()
            for chunk in nearby.split(","):
                ss = normalize_ss_type(chunk)
                if ss:
                    by_type.setdefault(ss, set()).add(tag)
    return by_type


def _jobid_from_dir(p: Path) -> int:
    """Extract the trailing numeric jobid from a RUN_DIR name.

    PBS template builds RUN_DIR=$HOME/runs/<sample>_<gpu>_<datetime>_<jobid>
    so the jobid is the last underscore-separated token.
    """
    try:
        return int(p.name.rsplit("_", 1)[-1])
    except (ValueError, IndexError):
        return -1


def find_run_pair(results_dir: Path, sample: str) -> tuple[Path | None, Path | None]:
    """Find the (proxy, whole) RUN_DIR pair for a sample.

    The PBS template builds dirs as <sample>_<gpu>_<datetime>_<jobid>
    with no condition suffix. The submit_gap_validation.sh script calls
    qsub twice per sample, proxy first then whole, so the proxy run gets
    a lower jobid than the whole run. We pair by sorting the per-sample
    runs by jobid and taking the lowest two.

    If there are fewer than 2 candidate directories, returns (None, None).
    Limitation: if the user submitted a sample more than once (reruns),
    this heuristic picks the two earliest jobids, which may not be what
    the user wants. In that case, run the script after manually cleaning
    up stale RUN_DIRs.
    """
    candidates = [p for p in results_dir.glob(f"{sample}_*") if p.is_dir() and _jobid_from_dir(p) > 0]
    if len(candidates) < 2:
        return None, None
    by_jobid = sorted(candidates, key=_jobid_from_dir)
    return by_jobid[0], by_jobid[1]


def find_substrates_tsv(run_dir: Path, sample: str) -> Path | None:
    """Locate <sample>_substrates_filtered.tsv anywhere under run_dir."""
    matches = list(run_dir.rglob(f"{sample}_substrates_filtered.tsv"))
    return matches[0] if matches else None


def score_genome(
    sample: str,
    results_dir: Path,
    ground_truth_dir: Path,
) -> list[dict]:
    """Compute per-SS-type recall for one genome under both conditions."""
    gt_path = ground_truth_dir / f"{sample}.tsv"
    if not gt_path.is_file():
        sys.stderr.write(f"WARN: no ground truth for {sample} at {gt_path}\n")
        return []

    truth = load_ground_truth(gt_path)
    rows: list[dict] = []

    proxy_dir, whole_dir = find_run_pair(results_dir, sample)
    predicted: dict[str, dict[str, set[str]]] = {}
    for cond, run_dir in (("proxy", proxy_dir), ("whole", whole_dir)):
        if run_dir is None:
            sys.stderr.write(f"WARN: no run dir for {sample}/{cond}\n")
            predicted[cond] = {}
            continue
        sub_tsv = find_substrates_tsv(run_dir, sample)
        if sub_tsv is None:
            sys.stderr.write(f"WARN: no substrates_filtered.tsv in {run_dir}\n")
            predicted[cond] = {}
            continue
        predicted[cond] = load_predicted_substrates(sub_tsv)

    all_ss_types = sorted(set(truth) | set(predicted["proxy"]) | set(predicted["whole"]))
    for ss in all_ss_types:
        truth_set = truth.get(ss, set())
        proxy_set = predicted["proxy"].get(ss, set())
        whole_set = predicted["whole"].get(ss, set())
        n_truth = len(truth_set)
        n_proxy = len(proxy_set)
        n_whole = len(whole_set)
        hit_proxy = len(truth_set & proxy_set)
        hit_whole = len(truth_set & whole_set)
        recall_proxy = hit_proxy / n_truth if n_truth else float("nan")
        recall_whole = hit_whole / n_truth if n_truth else float("nan")
        gap = recall_whole - recall_proxy if n_truth else float("nan")
        rows.append(
            {
                "genome": sample,
                "ss_type": ss,
                "n_truth": n_truth,
                "n_pred_proxy": n_proxy,
                "n_pred_whole": n_whole,
                "hits_proxy": hit_proxy,
                "hits_whole": hit_whole,
                "recall_proxy": f"{recall_proxy:.3f}" if n_truth else "NA",
                "recall_whole": f"{recall_whole:.3f}" if n_truth else "NA",
                "gap": f"{gap:.3f}" if n_truth else "NA",
            }
        )
    return rows


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--results-dir",
        required=True,
        type=Path,
        help="Directory containing gap_<sample>_<cond>_<jobid>/ subdirs from CX3",
    )
    p.add_argument(
        "--ground-truth-dir",
        type=Path,
        default=Path("research/secretion_classifier/data/ground_truth"),
        help="Directory containing per-genome ground_truth/<sample>.tsv files",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("gap_results.tsv"),
        help="Output TSV path",
    )
    args = p.parse_args()

    if not args.results_dir.is_dir():
        sys.stderr.write(f"FATAL: results-dir {args.results_dir} not a directory\n")
        return 2
    if not args.ground_truth_dir.is_dir():
        sys.stderr.write(f"FATAL: ground-truth-dir {args.ground_truth_dir} not a directory\n")
        return 2

    all_rows: list[dict] = []
    for sample in GENOMES:
        all_rows.extend(score_genome(sample, args.results_dir, args.ground_truth_dir))

    if not all_rows:
        sys.stderr.write("FATAL: no rows scored. Check directory layout.\n")
        return 1

    fieldnames = [
        "genome",
        "ss_type",
        "n_truth",
        "n_pred_proxy",
        "n_pred_whole",
        "hits_proxy",
        "hits_whole",
        "recall_proxy",
        "recall_whole",
        "gap",
    ]
    with args.output.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)

    # Stdout summary
    print(f"Wrote {len(all_rows)} rows to {args.output}\n")
    print(f"{'genome':<28}{'ss':<8}{'truth':>6}{'proxy':>7}{'whole':>7}{'r_p':>7}{'r_w':>7}{'gap':>7}")
    for r in all_rows:
        print(
            f"{r['genome']:<28}{r['ss_type']:<8}"
            f"{r['n_truth']:>6}{r['n_pred_proxy']:>7}{r['n_pred_whole']:>7}"
            f"{r['recall_proxy']:>7}{r['recall_whole']:>7}{r['gap']:>7}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
