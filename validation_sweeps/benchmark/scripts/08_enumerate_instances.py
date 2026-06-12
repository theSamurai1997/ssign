#!/usr/bin/env python3
"""Phase 0b, step 3.1: enumerate the distinct secretion-system instances to curate.

An "instance" is one secretion system in one genome: the key is
(ss_type, refseq_genome, sys_instance_id). sys_instance_id is the corpus's
within-genome label (a gene-range like "Xcp (xcpP-Z)" for T2SS, a family name
like "LEE"/"Hrp" for T3SS, or a representative locus for T1SS). Each instance
becomes one literature-curation job in step 3.4 and gets its own apparatus answer.

We enumerate from the corpus gold set only: those rows carry sys_instance_id and a
genome. The 59 folded external net-new rows have neither (external DBs don't group
effectors into instances), so they are NOT enumerated here; assigning them to an
instance is a Phase 1 concern (task 5.2, by genome+ss_type+proximity).

Rows with no refseq_genome cannot be placed (no genome to fetch/resolve against) and
are logged as a coverage gap, not curated. This hits T1SS hardest (it is spread thin
and mostly predicted, so few validated rows retain a genome).

Input : data/gold_build/gold_set_corpus.tsv
Output: data/machinery/instances.tsv          one row per instance (~92)
        data/machinery/unplaceable_effectors.tsv  rows dropped for no genome, with reason
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "data" / "gold_build" / "gold_set_corpus.tsv"
OUTDIR = ROOT / "data" / "machinery"


def slug(s: str) -> str:
    """Compact, filesystem-safe token from a sys_instance label ('Xcp (xcpP-Z)' -> 'Xcp')."""
    head = re.split(r"[\s(]", s.strip(), maxsplit=1)[0]
    return re.sub(r"[^A-Za-z0-9]", "", head) or "unlabeled"


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(GOLD.open(), delimiter="\t"))

    def has_genome(r: dict) -> bool:
        # "-" is the corpus placeholder for "no accession"; treat it as missing.
        return r["refseq_genome"].strip() not in ("", "-")

    placed = [r for r in rows if has_genome(r)]
    unplaceable = [r for r in rows if not has_genome(r)]

    # group effectors by instance key
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in placed:
        key = (r["ss_type"], r["refseq_genome"].strip(), r["sys_instance_id"].strip())
        groups[key].append(r)

    # stable ordering: ss_type, then genome, then label
    keys = sorted(groups, key=lambda k: (k[0], k[1], k[2]))
    counters: dict[str, int] = defaultdict(int)
    instances = []
    for ss_type, genome, label in keys:
        counters[ss_type] += 1
        members = groups[(ss_type, genome, label)]
        instances.append(
            {
                "instance_id": f"{ss_type}_{counters[ss_type]:02d}",
                "ss_type": ss_type,
                "refseq_genome": genome,
                "organism": members[0]["organism"],
                "sys_instance_label": label,
                "sys_slug": slug(label) if label else "unlabeled",
                "n_effectors": len(members),
                "effector_loci": ",".join(sorted(m["locus_tag"].strip() for m in members if m["locus_tag"].strip())),
            }
        )

    with (OUTDIR / "instances.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(instances[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(instances)

    with (OUTDIR / "unplaceable_effectors.tsv").open("w", newline="") as f:
        cols = ["gene", "uniprot", "locus_tag", "organism", "ss_type", "sys_instance_id", "primary_ref"]
        w = csv.DictWriter(f, fieldnames=cols + ["reason"], delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in unplaceable:
            w.writerow({**r, "reason": "no refseq_genome -> cannot fetch/resolve machinery; excluded from Phase 0b/1"})

    by_type: dict[str, int] = defaultdict(int)
    for it in instances:
        by_type[it["ss_type"]] += 1
    print(f"instances: {len(instances)}  by type: {dict(sorted(by_type.items()))}")
    print(f"genomes:   {len({it['refseq_genome'] for it in instances})}")
    print(f"unplaceable effectors (no genome): {len(unplaceable)}", end="  ")
    ut: dict[str, int] = defaultdict(int)
    for r in unplaceable:
        ut[r["ss_type"]] += 1
    print(f"by type: {dict(sorted(ut.items()))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
