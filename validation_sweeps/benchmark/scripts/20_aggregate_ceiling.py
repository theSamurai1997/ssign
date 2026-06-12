#!/usr/bin/env python3
"""
20_aggregate_ceiling.py  (Phase 1 task 5.4: aggregate ceiling / impossible)

Roll up ceiling_per_effector.tsv to the headline numbers: for each SS type (and each
genome), what fraction of TESTABLE verified effectors the +/-N proximity rule could reach,
at N = 3, 5, 7.

Denominator = testable effectors (own-instance machinery anchored AND effector locatable).
Untestable effectors are counted and reported alongside but NEVER enter the fraction, so the
ceiling is an honest "of the effectors we can actually evaluate, this many are reachable".
  ceiling@N    = reachable@N / testable
  impossible@N = 1 - ceiling@N   (testable but machinery is >N genes away or on another replicon)

Input : data/phase1/ceiling_per_effector.tsv
Output: data/phase1/ceiling_summary.tsv      (per SS type + an ALL row)
        data/phase1/ceiling_by_genome.tsv     (per genome x SS type)
Run:    .venv/bin/python scripts/20_aggregate_ceiling.py
"""

from __future__ import annotations

import csv
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
P1 = BENCH / "data" / "phase1"
IN = P1 / "ceiling_per_effector.tsv"
SUMMARY = P1 / "ceiling_summary.tsv"
BY_GENOME = P1 / "ceiling_by_genome.tsv"
NS = (3, 5, 7)


def load(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def summarize(rows):
    """rows for one group -> dict of counts + ceiling/impossible fractions per N."""
    n_gold = len(rows)
    test = [r for r in rows if r["testable"] == "yes"]
    nt = len(test)
    out = {"n_gold": n_gold, "n_testable": nt, "n_untestable": n_gold - nt}
    for n in NS:
        reach = sum(r[f"reachable_n{n}"] == "true" for r in test)
        out[f"reachable_n{n}"] = reach
        out[f"ceiling_n{n}"] = round(reach / nt, 4) if nt else ""
        out[f"impossible_n{n}"] = round(1 - reach / nt, 4) if nt else ""
    return out


def write(path, key_field, groups):
    cols = [key_field, "n_gold", "n_testable", "n_untestable"]
    for n in NS:
        cols += [f"reachable_n{n}", f"ceiling_n{n}", f"impossible_n{n}"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for key, rows in groups:
            w.writerow({key_field: key, **summarize(rows)})


def main():
    rows = load(IN)
    order = ["T1SS", "T2SS", "T3SS", "T4SS", "T6SS"]

    by_type = [(ss, [r for r in rows if r["ss_type"] == ss]) for ss in order]
    by_type = [(ss, rs) for ss, rs in by_type if rs] + [("ALL", rows)]
    write(SUMMARY, "ss_type", by_type)

    # per genome x type, testable genomes first, ordered by type then accession
    genomes = sorted(
        {(r["ss_type"], r["refseq_genome"]) for r in rows if r["refseq_genome"].strip() not in ("", "-")},
        key=lambda k: (order.index(k[0]) if k[0] in order else 9, k[1]),
    )
    by_geno = [(f"{ss}:{g}", [r for r in rows if r["ss_type"] == ss and r["refseq_genome"] == g]) for ss, g in genomes]
    write(BY_GENOME, "type_genome", by_geno)

    print(f"wrote {SUMMARY.relative_to(BENCH)} and {BY_GENOME.relative_to(BENCH)}\n")
    print(f"  {'type':5s} {'gold':>4s} {'test':>4s} {'untest':>6s}   ceiling@3 / @5 / @7")
    for ss, rs in by_type:
        s = summarize(rs)
        c = "  ".join(f"{s[f'ceiling_n{n}']:.0%}" if s[f"ceiling_n{n}"] != "" else "  -" for n in NS)
        print(f"  {ss:5s} {s['n_gold']:4d} {s['n_testable']:4d} {s['n_untestable']:6d}   {c}")


if __name__ == "__main__":
    raise SystemExit(main())
