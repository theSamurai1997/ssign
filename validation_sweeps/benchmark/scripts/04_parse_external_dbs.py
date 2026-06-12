#!/usr/bin/env python3
"""Phase 0a, step 2b.2: parse the external effector DBs to one common schema.

Sources (downloaded 2026-06-10, see data/external_dbs/MANIFEST.md):
  SecReT4  data/external_dbs/secret4/verified_effectors.fas   (540 verified T4SS effectors)
           header: >gi|<gi>|ref|<protein_acc>| <description> [<organism>]
  SecReT6  data/external_dbs/secret6/effector_exp_protein.fasta (331 experimental T6SS effectors)
           header: >EFF##### <genome_acc>:<coords> [<organism>]

We keep only the experimentally-supported subset of each DB (SecReT4 'verified',
SecReT6 '_exp_'); the prediction tiers are not used. EffectiveDB is excluded
entirely (it is a prediction engine with no per-effector experimental citations).

Common schema (one row per external effector):
  source_db, source_id, ss_type, organism, refseq_genome, protein_acc, coords, description

protein_acc / coords are whatever the DB gives us; the next step (05) resolves
these to a UniProt/locus key for dedup against the corpus gold set.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data" / "external_dbs"
OUT = EXT / "parsed_external.tsv"

HEADER = [
    "source_db",
    "source_id",
    "ss_type",
    "organism",
    "refseq_genome",
    "protein_acc",
    "coords",
    "description",
]

ORG_RE = re.compile(r"\[([^\]]+)\]\s*$")


def parse_secret4(path: Path) -> list[dict]:
    """SecReT4 verified effectors: gi|ref header carries a RefSeq protein acc."""
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(">"):
            continue
        h = line[1:].strip()
        org_m = ORG_RE.search(h)
        organism = org_m.group(1).strip() if org_m else ""
        body = ORG_RE.sub("", h).strip()
        protein_acc, gi, desc = "", "", body
        # >gi|52840385|ref|YP_094184.1| hypothetical protein lpg0130
        m = re.match(r"gi\|(\d+)\|ref\|([^|]+)\|\s*(.*)$", body)
        if m:
            gi, protein_acc, desc = m.group(1), m.group(2).strip(), m.group(3).strip()
        rows.append(
            {
                "source_db": "SecReT4",
                "source_id": protein_acc or f"gi:{gi}",
                "ss_type": "T4SS",
                "organism": organism,
                "refseq_genome": "",
                "protein_acc": protein_acc,
                "coords": "",
                "description": desc,
            }
        )
    return rows


def parse_secret6(path: Path) -> list[dict]:
    """SecReT6 experimental effectors: EFF id + genome_acc:coords header."""
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(">"):
            continue
        h = line[1:].strip()
        org_m = ORG_RE.search(h)
        organism = org_m.group(1).strip() if org_m else ""
        body = ORG_RE.sub("", h).strip()
        # >EFF01497 NC_016776:c2339019-2336050
        m = re.match(r"(\S+)\s+(\S+):(\S+)$", body)
        eff_id, genome, coords = "", "", ""
        if m:
            eff_id, genome, coords = m.group(1), m.group(2), m.group(3)
        else:
            eff_id = body.split()[0] if body else ""
        rows.append(
            {
                "source_db": "SecReT6",
                "source_id": eff_id,
                "ss_type": "T6SS",
                "organism": organism,
                "refseq_genome": genome,
                "protein_acc": "",
                "coords": coords,
                "description": "",
            }
        )
    return rows


def main() -> int:
    s4 = parse_secret4(EXT / "secret4" / "verified_effectors.fas")
    s6 = parse_secret6(EXT / "secret6" / "effector_exp_protein.fasta")
    allrows = s4 + s6

    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t")
        w.writeheader()
        w.writerows(allrows)

    from collections import Counter

    print(f"SecReT4: {len(s4)} | SecReT6: {len(s6)} | total {len(allrows)} -> {OUT}")
    print("SecReT4 with protein_acc:", sum(1 for r in s4 if r["protein_acc"]))
    print("SecReT6 with genome+coords:", sum(1 for r in s6 if r["refseq_genome"] and r["coords"]))
    print("\nSecReT4 top organisms:")
    for o, c in Counter(r["organism"][:45] for r in s4).most_common(10):
        print(f"  {c:4} {o}")
    print("\nSecReT6 top organisms:")
    for o, c in Counter(r["organism"][:45] for r in s6).most_common(10):
        print(f"  {c:4} {o}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
