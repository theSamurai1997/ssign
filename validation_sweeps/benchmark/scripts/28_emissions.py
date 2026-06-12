#!/usr/bin/env python3
"""Phase 2 task 8.1: assemble ssign's full substrate-emission set (the precision denominator).

Recall (scripts 24/25) asks "of known effectors, how many did ssign emit?". Precision asks the
opposite: "of everything ssign emits as a substrate, how much is real?". This script builds that
emission set: one row per protein in every panel run's `# Secreted Proteins` chunk (ssign's substrate
output), tagged with how it was called (`substrate_source`: proximity vs T5SS-self), the SS type it
was assigned to, its annotation, the per-protein tool signals, its sequence, and whether it is one of
the 51 gold effectors (the only known-true emissions). Tiers 8.2/8.3 then bound the unlabelled rest.

Inputs : data/phase2/runs/<tag>/<unit>/results/<unit>_results.csv   (# Secreted Proteins chunk)
         data/phase2/actual_per_effector.<tag>.tsv                  (the 51 gold emissions -> is_gold)
Outputs: data/phase2/emissions.<tag>.tsv                            (one row per emission)
Run:     .venv/bin/python scripts/28_emissions.py --run-tag panel_genbank_default
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import bench_runout as bo  # noqa: E402
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
RUNS = BENCH / "data" / "phase2" / "runs"

ANNOT_COLS = ("broad_annotation", "detailed_annotation", "gbff_annotation", "product")
# substrate_source + nearby_ss_types already come from SIGNAL_COLS — don't list them twice.
OUT_COLS = (
    "unit_id",
    "locus_tag",
    "aa_length",
    *ANNOT_COLS,
    *bo.SIGNAL_COLS,
    "is_gold",
    "gold_gene",
    "sequence",
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-tag", required=True)
    args = ap.parse_args()
    run_root = RUNS / args.run_tag
    if not run_root.is_dir():
        sys.exit(f"no run dir: {run_root.relative_to(BENCH)}")

    # The 51 gold emissions: emitted_secreted rows in the actual-call table, keyed by (unit, ssign_locus).
    gold = {}
    for r in read_tsv(BENCH / "data" / "phase2" / f"actual_per_effector.{args.run_tag}.tsv"):
        if r["ssign_call"] == bo.CALL_EMITTED:
            gold[(r["unit_id"], r["ssign_locus"])] = r["gene"]

    rows = []
    for unit_dir in sorted(run_root.iterdir()):
        rc = unit_dir / "results" / f"{unit_dir.name}_results.csv"
        if not rc.exists():
            continue
        for r in bo._read_secreted_chunk(rc):
            locus = (r.get("locus_tag") or "").strip()
            if not locus:
                continue
            key = (unit_dir.name, locus)
            rows.append(
                {
                    "unit_id": unit_dir.name,
                    "locus_tag": locus,
                    "aa_length": (r.get("aa_length") or "").strip(),
                    **{c: (r.get(c) or "").strip() for c in ANNOT_COLS},
                    **{c: (r.get(c) or "").strip() for c in bo.SIGNAL_COLS},
                    "is_gold": "yes" if key in gold else "no",
                    "gold_gene": gold.get(key, ""),
                    "sequence": (r.get("sequence") or "").strip().upper(),
                }
            )

    out = BENCH / "data" / "phase2" / f"emissions.{args.run_tag}.tsv"
    write_tsv(out, list(OUT_COLS), rows)
    src = Counter(r["substrate_source"] for r in rows)
    ngold = sum(r["is_gold"] == "yes" for r in rows)
    noseq = sum(not r["sequence"] for r in rows)
    print(f"wrote {out.relative_to(BENCH)}  ({len(rows)} emissions)")
    print(f"  by substrate_source: {dict(src)}")
    print(f"  is_gold: {ngold}   |  proximity subset: {src.get('proximity', 0)}   |  no-sequence: {noseq}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
