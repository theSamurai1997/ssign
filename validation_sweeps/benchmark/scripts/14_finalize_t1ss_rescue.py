#!/usr/bin/env python3
"""
14_finalize_t1ss_rescue.py  (Phase 0a augmentation: merge the two rescue passes)

Combines the IPG pass (script 12, identical-protein placements) and the BLAST pass
(script 13, >=95% representative-strain placements) into one rescue table for the 19
T1SS effectors that the corpus left without a genome. Adds a `species_match` flag
(exact = same genus+species as the characterized strain; genus_only = same genus,
representative species; mismatch) so the slightly weaker placements are explicit, not
hidden -- per Teo's representative_strain decision.

This does NOT yet mutate the gold set. It produces the reviewable rescue table; folding
the placements (and fetching their genomes + locating the cognate T1SS transporter for
the ceiling measurement) happens in Phase 1, where the gene-order machinery already lives.

Inputs:
  data/t1ss_rescue/ipg_placements.tsv
  data/t1ss_rescue/blast_placements.tsv
Output:
  data/t1ss_rescue/t1ss_rescued.tsv     one row per effector, final placement + flags

Run:
  .venv/bin/python scripts/14_finalize_t1ss_rescue.py
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
RESCUE = BENCH / "data" / "t1ss_rescue"
IPG_TSV = RESCUE / "ipg_placements.tsv"
BLAST_TSV = RESCUE / "blast_placements.tsv"
OUT_TSV = RESCUE / "t1ss_rescued.tsv"

FIELDS = [
    "gene",
    "uniprot",
    "corpus_organism",
    "final_status",
    "tier",
    "placement_assembly",
    "placement_nuc",
    "start",
    "stop",
    "strand",
    "placement_protein",
    "placement_organism",
    "species_match",
    "pident",
    "qcov",
    "n_genome_options",
    "primary_ref",
]


def species_tokens(org: str):
    """(genus, species) lowercased from an organism string, ignoring parens/strain."""
    head = org.split("(")[0].strip().split()
    genus = head[0].lower() if head else ""
    species = head[1].lower() if len(head) > 1 else ""
    return genus, species


def match_level(corpus_org, placement_org):
    if not placement_org:
        return ""
    cg, cs = species_tokens(corpus_org)
    pg, ps = species_tokens(placement_org)
    if cg and cg == pg and cs and cs == ps:
        return "exact"
    if cg and cg == pg:
        return "genus_only"
    return "mismatch"


def main():
    ipg = {r["uniprot"]: r for r in csv.DictReader(open(IPG_TSV), delimiter="\t")}
    blast = {r["uniprot"]: r for r in csv.DictReader(open(BLAST_TSV), delimiter="\t")} if BLAST_TSV.exists() else {}

    rows = []
    for acc, ir in ipg.items():
        out = {f: "" for f in FIELDS}
        out.update(gene=ir["gene"], uniprot=acc, corpus_organism=ir["organism"], primary_ref=ir["primary_ref"])
        if ir["status"] == "placed_ipg":
            src = ir
            out["final_status"] = "placed_ipg"
        elif acc in blast and blast[acc]["status"] == "placed_blast":
            src = blast[acc]
            out["final_status"] = "placed_blast"
            out["pident"], out["qcov"] = src.get("pident", ""), src.get("qcov", "")
        else:
            b = blast.get(acc, {})
            out["final_status"] = b.get("status", "unplaceable")
            out["pident"], out["qcov"] = b.get("pident", ""), b.get("qcov", "")
            rows.append(out)
            continue
        out.update(
            tier=src.get("tier", ""),
            placement_assembly=src.get("placement_assembly", ""),
            placement_nuc=src.get("placement_nuc", ""),
            start=src.get("start", ""),
            stop=src.get("stop", ""),
            strand=src.get("strand", ""),
            placement_protein=src.get("placement_protein", ""),
            placement_organism=src.get("placement_organism", ""),
            species_match=match_level(ir["organism"], src.get("placement_organism", "")),
            n_genome_options=src.get("n_genome_options", ""),
        )
        rows.append(out)

    rows.sort(key=lambda r: (r["final_status"] != "placed_ipg", r["final_status"] != "placed_blast", r["gene"]))
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    placed = [r for r in rows if r["final_status"].startswith("placed")]
    print(f"wrote {OUT_TSV.relative_to(BENCH)}")
    print("  19 unplaced T1SS effectors")
    print(
        f"  RESCUED (placed): {len(placed)}  "
        f"[ipg {sum(1 for r in placed if r['final_status'] == 'placed_ipg')} / "
        f"blast {sum(1 for r in placed if r['final_status'] == 'placed_blast')}]"
    )
    print(f"  still unplaced  : {len(rows) - len(placed)}")
    print(f"  species_match   : {dict(Counter(r['species_match'] for r in placed))}")
    print()
    for r in rows:
        flag = r["species_match"] or "-"
        if r["final_status"] == "placed_ipg":
            idcov = "100% (IPG)"
        elif r["final_status"] == "placed_blast":
            idcov = f"id={r['pident']}/cov={r['qcov']}"
        else:
            idcov = r["final_status"]
        print(
            f"  {r['gene']:10s} {r['uniprot']:8s} {r['final_status']:13s} {flag:10s} "
            f"{r['placement_nuc']:18s} {idcov:14s} {r['placement_organism'][:34]}"
        )


if __name__ == "__main__":
    raise SystemExit(main())
