#!/usr/bin/env python3
"""
17_fold_t1ss_rescue.py  (Phase 1 prep: fold the rescued T1SS placements into the gold set)

19 of 28 validated T1SS effectors had no genome in the corpus, so they were untestable
for the proximity ceiling. Scripts 12-16 placed 16 of them into representative RefSeq
genomes and verified, by reading the gene order, that the T1SS transporter operon sits
1-3 genes away (the adjacency scan = an ssign-independent, literature-family ceiling).
This script folds those 16 placements back into the effector gold set so Phase 1 can treat
them on the same footing as the curated-machinery effectors.

What it does, per gold-set row (keyed by UniProt accession):
  - rescued + placed (16): fill refseq_genome with the placement replicon, give it a
    synthetic per-effector instance id (T1SS_Rnn -- T1SS is one-effector-one-operon, so
    each placed effector is its own system instance), and carry the placement coordinates,
    tier (ipg_identical / representative_strain), species_match, and
    ceiling_source='adjacency_scan' (its ceiling comes from t1ss_ceiling.tsv, not the
    literature answer key).
  - rescued but unplaceable (3): testable='no' with the reason; kept in the gold set,
    excluded from the Phase 1 denominator, documented.
  - everyone else (the 9 curated-machinery T1SS effectors + all T2-T6SS): untouched
    except ceiling_source='answer_key'; testability is decided in script 19 by whether the
    effector's own instance has an anchored machinery locus.

Inputs : data/gold_build/effector_gold_set.tsv  (582 rows)
         data/t1ss_rescue/t1ss_rescued.tsv       (19 rows: 16 placed + 3 unplaceable)
         data/t1ss_rescue/t1ss_ceiling.tsv       (16 placed: adjacency + nearest dist)
Outputs: data/phase1/effector_gold_set_phase1.tsv  (582 rows + Phase 1 columns)
         data/phase1/fold_provenance.tsv           (every rescued row's fate)

Run:
  .venv/bin/python scripts/17_fold_t1ss_rescue.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
GOLD = BENCH / "data" / "gold_build" / "effector_gold_set.tsv"
RESCUED = BENCH / "data" / "t1ss_rescue" / "t1ss_rescued.tsv"
CEILING = BENCH / "data" / "t1ss_rescue" / "t1ss_ceiling.tsv"
OUT = BENCH / "data" / "phase1" / "effector_gold_set_phase1.tsv"
PROV = BENCH / "data" / "phase1" / "fold_provenance.tsv"

# Phase 1 columns appended to every gold-set row
NEW_COLS = [
    "ceiling_source",  # 'answer_key' | 'adjacency_scan'
    "placement_tier",  # ipg_identical | representative_strain | '' (curated)
    "species_match",  # exact | genus_only | '' (curated)
    "placement_start",  # effector CDS coords on the placement replicon (rescued only)
    "placement_stop",
    "placement_strand",
    "placement_effector_locus",
    "testable",  # 'yes' | 'no' | '' (decided in script 19 for answer_key rows)
    "testable_reason",
]

# the 3 T1SS effectors the rescue could not place, by UniProt accession (see tasks.md 4.3)
UNPLACEABLE_REASON = {
    "Q07162": "best genome match only 61% identity (genuinely divergent)",
    "P82115": "no RefSeq genome for the exact strain (Photorhabdus sp. Az29)",
    "P55123": "no RefSeq genome for the exact strain (Pasteurella haemolytica-like 5943B)",
}


def load(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    gold = load(GOLD)
    rescued = {r["uniprot"]: r for r in load(RESCUED)}
    ceiling = {r["uniprot"]: r for r in load(CEILING)}

    placed = [u for u, r in rescued.items() if r["final_status"].startswith("placed")]
    unplaceable = [u for u, r in rescued.items() if not r["final_status"].startswith("placed")]
    if len(placed) != len(ceiling):
        sys.exit(f"placed={len(placed)} but ceiling has {len(ceiling)} rows; expected equal")

    # stable synthetic instance ids, ordered by gene then accession for reproducibility
    rnn = {u: f"T1SS_R{i:02d}" for i, u in enumerate(sorted(placed, key=lambda u: (rescued[u]["gene"], u)), 1)}

    fieldnames = list(gold[0].keys()) + [c for c in NEW_COLS if c not in gold[0]]
    prov = []
    folded = unplaced = 0
    for r in gold:
        for c in NEW_COLS:
            r.setdefault(c, "")
        u = r["uniprot"]
        if u in rnn:  # rescued + placed
            res, cel = rescued[u], ceiling[u]
            r["refseq_genome"] = res["placement_nuc"]
            r["sys_instance_id"] = rnn[u]
            r["ceiling_source"] = "adjacency_scan"
            r["placement_tier"] = res["tier"]
            r["species_match"] = res["species_match"]
            r["placement_start"] = res["start"]
            r["placement_stop"] = res["stop"]
            r["placement_strand"] = res["strand"]
            r["placement_effector_locus"] = cel["effector_locus"]
            r["testable"] = "yes"
            r["testable_reason"] = ""
            folded += 1
            prov.append(
                {
                    "uniprot": u,
                    "gene": r["gene"],
                    "fate": "folded_placed",
                    "instance_id": rnn[u],
                    "placement_nuc": res["placement_nuc"],
                    "tier": res["tier"],
                    "species_match": res["species_match"],
                    "adjacency": cel["adjacency"],
                    "nearest_dist": cel["nearest_component_dist"],
                }
            )
        elif u in unplaceable:
            reason = UNPLACEABLE_REASON.get(u, "no genome placement")
            r["ceiling_source"] = "adjacency_scan"
            r["testable"] = "no"
            r["testable_reason"] = reason
            unplaced += 1
            prov.append(
                {
                    "uniprot": u,
                    "gene": r["gene"],
                    "fate": "unplaceable",
                    "instance_id": "",
                    "placement_nuc": "",
                    "tier": "",
                    "species_match": "",
                    "adjacency": "",
                    "nearest_dist": "",
                    "reason": reason,
                }
            )
        else:  # curated-machinery effector (9 T1SS + all T2-T6SS)
            r["ceiling_source"] = "answer_key"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(gold)

    prov_fields = [
        "uniprot",
        "gene",
        "fate",
        "instance_id",
        "placement_nuc",
        "tier",
        "species_match",
        "adjacency",
        "nearest_dist",
        "reason",
    ]
    with open(PROV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=prov_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(prov)

    print(f"wrote {OUT.relative_to(BENCH)}  ({len(gold)} rows, +{len(NEW_COLS)} cols)")
    print(f"  folded placed     : {folded}  (instance ids {rnn[placed[0]]}..{rnn[placed[-1]]} after sort)")
    print(f"  flagged unplaceable: {unplaced}")
    print(f"wrote {PROV.relative_to(BENCH)}")


if __name__ == "__main__":
    raise SystemExit(main())
