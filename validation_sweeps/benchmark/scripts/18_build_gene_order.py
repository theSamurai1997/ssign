#!/usr/bin/env python3
"""
18_build_gene_order.py  (Phase 1 task 5.1: per-genome gene-order index)

Parse every cached genome once and emit a flat, alias-aware gene-order index that Phase 1
(ceiling) and Phase 2 (recall bridge) both read. One row per CDS:
  record_acc, ordinal (0-based position on that replicon, sorted by start), locus_tag, gene,
  aliases (gene/old_locus_tag/gene_synonym/locus_tag, ';'-joined), start, end, strand.

The heavy lifting (feature merge, ordinal assignment, drift-tolerant accession + locus_tag
matching) lives in bench_index.py so downstream scripts get identical behaviour.

Input : data/refseq_cache/*.gb
Output: data/phase1/gene_order_index.tsv
Run:    .venv/bin/python scripts/18_build_gene_order.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import bench_index as bi  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
GOLD = BENCH / "data" / "phase1" / "effector_gold_set_phase1.tsv"
OUT = bi.INDEX_TSV


def panel_accessions():
    """Every genome accession a gold effector is placed in (for a coverage report)."""
    accs = set()
    with open(GOLD) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["refseq_genome"].strip()
            if g and g != "-":
                accs.add(g)
    return accs


def main():
    idx = bi.build_from_genbank()  # all cached records; index is reused by Phase 2 too
    OUT.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=bi.INDEX_FIELDS, delimiter="\t")
        w.writeheader()
        for rec, cds in idx.records.items():
            for d in cds:
                w.writerow(
                    {
                        "record_acc": rec,
                        "ordinal": d["ordinal"],
                        "locus_tag": d["locus_tag"],
                        "gene": d["gene"],
                        "aliases": ";".join(sorted(a.lower() for a in d["aliases"])),
                        "start": d["start"],
                        "end": d["end"],
                        "strand": d["strand"],
                    }
                )
                n += 1

    # coverage: does every gold-effector genome resolve in the index?
    panel = panel_accessions()
    unresolved = sorted(a for a in panel if idx.resolve_record(a) is None)
    print(f"wrote {OUT.relative_to(BENCH)}")
    print(f"  replicons indexed : {len(idx.records)}")
    print(f"  CDS rows          : {n}")
    print(f"  gold genomes      : {len(panel)} distinct; resolved {len(panel) - len(unresolved)}/{len(panel)}")
    if unresolved:
        print("  UNRESOLVED gold genomes (need fetch):")
        for a in unresolved:
            print(f"    {a}")


if __name__ == "__main__":
    raise SystemExit(main())
