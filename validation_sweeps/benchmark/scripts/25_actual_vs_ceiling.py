#!/usr/bin/env python3
"""
25_actual_vs_ceiling.py  (Phase 2 task 6.5: actual recall vs the Phase 1 ceiling)

Per SS type, on the testable denominator (same as Phase 1):
  ceiling@N  = fraction reachable@N (could the +/-N rule reach the effector's own machinery)
  actual     = fraction ssign actually emitted as secreted
and the concordance 2x2 at N=3:
  emitted & reachable@3      true positive (ssign got a reachable effector)
  emitted & NOT reachable@3  ssign emitted it via machinery MacSyFinder detected that the
                             literature answer key didn't anchor -> actual can exceed ceiling
                             here; these rows are listed for inspection, NOT treated as an error.
  not-emitted & reachable@3  ssign MISS within reach -> the interesting false negatives
  not-emitted & NOT reachable@3  out of reach anyway

Input : data/phase2/actual_per_effector.<run_tag>.tsv  (from script 24)
Output: data/phase2/actual_vs_ceiling.<run_tag>.tsv
Run:    .venv/bin/python scripts/25_actual_vs_ceiling.py --run-tag pilot_genbank_default
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import bench_index as bi  # noqa: E402
import bench_runout as bo  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
load = bi.load_tsv
NS = (3, 5, 7)
SS_ORDER = ("T1SS", "T2SS", "T3SS", "T4SS", "T6SS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-tag", required=True)
    args = ap.parse_args()
    per_eff = BENCH / "data" / "phase2" / f"actual_per_effector.{args.run_tag}.tsv"
    if not per_eff.exists():
        sys.exit(f"missing {per_eff.relative_to(BENCH)} -- run 24_actual_call.py --run-tag {args.run_tag} first")
    rows = load(per_eff)

    def reach(r, n):
        return r.get(f"reachable_n{n}") == "true"

    def emitted(r):
        return r["ssign_call"] == bo.CALL_EMITTED

    out_rows = []
    for ss in (*SS_ORDER, "ALL"):
        sel = [r for r in rows if r["testable"] == "yes" and (ss == "ALL" or r["ss_type"] == ss)]
        n = len(sel)
        if not n:
            continue
        rec = {"ss_type": ss, "n_testable": n}
        for k in NS:
            rec[f"ceiling_n{k}"] = sum(reach(r, k) for r in sel)
        rec["actual_emitted"] = sum(emitted(r) for r in sel)
        # concordance 2x2 at N=3
        rec["emit_reach3"] = sum(emitted(r) and reach(r, 3) for r in sel)
        rec["emit_unreach3"] = sum(emitted(r) and not reach(r, 3) for r in sel)
        rec["miss_reach3"] = sum(not emitted(r) and reach(r, 3) for r in sel)
        rec["miss_unreach3"] = sum(not emitted(r) and not reach(r, 3) for r in sel)
        out_rows.append(rec)

    cols = [
        "ss_type",
        "n_testable",
        "ceiling_n3",
        "ceiling_n5",
        "ceiling_n7",
        "actual_emitted",
        "emit_reach3",
        "emit_unreach3",
        "miss_reach3",
        "miss_unreach3",
    ]
    out_path = BENCH / "data" / "phase2" / f"actual_vs_ceiling.{args.run_tag}.tsv"
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"wrote {out_path.relative_to(BENCH)}  (run_tag={args.run_tag})")
    print(f"\n  {'SS':5} {'test':>5} {'ceil@3':>7} {'ceil@5':>7} {'ceil@7':>7} {'actual':>7}   (emitted of testable)")
    for r in out_rows:
        n = r["n_testable"]
        pct = lambda x: f"{100 * x / n:3.0f}%"
        print(
            f"  {r['ss_type']:5} {n:5d} {pct(r['ceiling_n3']):>7} {pct(r['ceiling_n5']):>7} "
            f"{pct(r['ceiling_n7']):>7} {pct(r['actual_emitted']):>7}"
        )
    # emitted-but-unreachable: ssign found machinery the answer key didn't anchor (not an error)
    discord = [
        r
        for r in rows
        if r["testable"] == "yes" and r["ssign_call"] == bo.CALL_EMITTED and r.get("reachable_n7") != "true"
    ]
    print(f"\n  emitted but ceiling-unreachable@7 (ssign machinery != literature answer key): {len(discord)}")
    for r in discord[:12]:
        print(
            f"    {r['gene']}/{r['uniprot']} ({r['ss_type']}) nearby={r.get('nearby_ss_types', '')} reason={r.get('ceiling_reason', '')}"
        )
    if len(discord) > 12:
        print(f"    ... +{len(discord) - 12} more")


if __name__ == "__main__":
    raise SystemExit(main())
