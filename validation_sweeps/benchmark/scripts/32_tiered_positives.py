#!/usr/bin/env python3
"""Dataset group 1, task 1.3: union validated + audited-predicted into a tiered table.

Combines the frozen Phase-1 gold set (582 validated, instance-tagged) with the
audited-predicted survivors (predicted_audited.tsv) into one positive table carrying
an `evidence_tier` column (validated | predicted). The tier is what the model-training
loss later uses to down-weight predicted relative to validated; this change only emits
the column, it does not pick a weight ratio.

Schema is the union of both inputs: the shared corpus core, plus tier-specific columns
(validated carries Phase-1 placement/testable fields; predicted carries citation_status,
audit_tier). Missing cells are blank, never invented.

Inputs:
  data/phase1/effector_gold_set_phase1.tsv   (validated backbone)
  data/dataset/predicted_audited.tsv         (predicted survivors, from step 31)
Output:
  data/dataset/positives_tiered.tsv
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from bench_io import by_type, read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "dataset"
VALIDATED = ROOT / "data" / "phase1" / "effector_gold_set_phase1.tsv"
PREDICTED = OUT / "predicted_audited.tsv"

# Canonical key order; any input column not listed is appended in first-seen order.
LEAD = [
    "gene",
    "uniprot",
    "locus_tag",
    "organism",
    "refseq_genome",
    "ss_type",
    "sys_instance_id",
    "evidence_level",
    "evidence_tier",
    "primary_ref",
    "verification_status",
    "verification_notes",
    "citation_status",
    "audit_tier",
]


def main() -> int:
    validated = read_tsv(VALIDATED)
    predicted = read_tsv(PREDICTED)
    for r in validated:
        r["evidence_tier"] = "validated"
    for r in predicted:
        r["evidence_tier"] = "predicted"

    rows = validated + predicted
    header = list(LEAD)
    for r in rows:  # append any remaining columns in first-seen order
        for k in r:
            if k not in header:
                header.append(k)

    write_tsv(OUT / "positives_tiered.tsv", header, rows)

    print(f"validated: {len(validated):>4}   {by_type(validated)}")
    print(f"predicted: {len(predicted):>4}   {by_type(predicted)}")
    print(f"total:     {len(rows):>4}   {by_type(rows)}")
    cites = Counter(r.get("citation_status", "") for r in predicted)
    print(f"predicted citation_status: {dict(cites)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
