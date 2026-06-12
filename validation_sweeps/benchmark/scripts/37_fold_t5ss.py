"""Dataset group 3, task 3.5 (+ group 2 reconciliation): build the final positive table.

Unions the three label-side branches into positives_all.tsv:
  - validated gold set  (positives_tiered.tsv, evidence_tier=validated, instance from Phase 1)
  - predicted corpus    (positives_tiered.tsv predicted rows, with the group-2 instance
                         assignment layered in from positives_instanced.tsv)
  - sourced T5SS        (t5ss_effectors.tsv, resolved + verified in step 36)

Two things this fold reconciles:

  1. Instance assignment for predicted rows. Script 32 wrote the predicted rows BEFORE the
     group-2 instance assignment (33/35). This fold layers positives_instanced.tsv's
     instance_source / type_level / provenance back onto those rows by (gene, locus_tag,
     refseq_genome, primary_ref) key, so the deliverable carries the group-2 work.

  2. One canonical instance column. The validated branch keys the benchmark instance id in
     sys_instance_id (e.g. T1SS_R06); the predicted branch kept it in instance_id (the
     predicted sys_instance_id held an unreliable name). This fold unifies on sys_instance_id
     for every tier, so group 4's feature join keys one column.

instance_source is harmonized across tiers: gold (validated, Phase 1) | auto | literature
(predicted) | self (T5a/c/d/e autotransporter) | none (unresolved / type-level). type_level
is "yes" when no specific instance is assigned. self_secreted/subtype are T5SS-only.

Inputs:  data/dataset/positives_tiered.tsv, positives_instanced.tsv, t5ss_effectors.tsv
Output:  data/dataset/positives_all.tsv
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from bench_io import by_type, read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "dataset"

# predicted-only provenance carried over from the group-2 instance assignment
INSTANCE_PROV = ["instance_source", "type_level", "instance_candidates", "instance_quote", "instance_source_doi"]


def _key(r: dict) -> tuple:
    # uniprot is in the key because (gene, locus, genome, ref) is NOT unique: a few rows
    # share a blank locus+genome (e.g. lktA cited twice for different UniProt accessions).
    return (
        r.get("gene", ""),
        r.get("locus_tag", ""),
        r.get("refseq_genome", ""),
        r.get("primary_ref", ""),
        r.get("uniprot", ""),
    )


def main() -> int:
    tiered = read_tsv(OUT / "positives_tiered.tsv")
    instanced = read_tsv(OUT / "positives_instanced.tsv")
    t5_path = OUT / "t5ss_effectors.tsv"
    t5 = read_tsv(t5_path) if t5_path.exists() else []

    inst_by = {_key(r): r for r in instanced}
    assert len(inst_by) == len(instanced), "predicted instance key collision -- _key is not unique"

    # canonical benchmark instance for validated effectors: the gold-set sys_instance_id is a
    # messy literature label (names, loci, blanks); ceiling_per_effector.instance_id is the
    # benchmark's authoritative assignment (same id space as the machinery answer key + the
    # predicted tier), so group 4's feature join keys one consistent column.
    ceiling = read_tsv(ROOT / "data" / "phase1" / "ceiling_per_effector.tsv")
    ceil_by = {(r["gene"], r["uniprot"], r["refseq_genome"]): r for r in ceiling}

    for r in tiered:
        if r["evidence_tier"] == "validated":
            r["instance_source"] = "gold"
            # take the instance straight from the ceiling (every validated row has one): canonical
            # id where the benchmark placed it, blank where it couldn't -> strictly canonical-or-blank,
            # never a stray gold label that no machinery key would match.
            c = ceil_by.get((r["gene"], r["uniprot"], r["refseq_genome"]))
            r["sys_instance_id"] = c["instance_id"].strip() if c else ""
            r["type_level"] = "no" if r["sys_instance_id"] else "yes"
        else:  # predicted: layer in the group-2 instance assignment, unify id into sys_instance_id
            m = inst_by.get(_key(r))
            if m:
                r["sys_instance_id"] = m.get("instance_id", "").strip()
                for c in INSTANCE_PROV:
                    r[c] = m.get(c, "")
            else:  # no instance row (should not happen; keep as type-level rather than guess)
                r["instance_source"], r["type_level"] = "none", "yes"

    for r in t5:  # map T5SS rows onto the tiered schema
        r["evidence_level"] = "validated"
        r["evidence_tier"] = "validated" if r.get("verified") == "yes" else "predicted"
        r["sys_instance_id"] = ""
        r["citation_status"] = "RESOLVED" if r.get("doi_resolves") == "yes" else "UNRESOLVED"
        # self-secreted autotransporter is its own system; T5b TpsA is a type-level substrate
        self_sec = r.get("self_secreted") == "true"
        r["instance_source"] = "self" if self_sec else "none"
        r["type_level"] = "no" if self_sec else "yes"

    rows = tiered + t5
    header = list(tiered[0].keys()) if tiered else []
    for r in rows:  # append instance + T5SS-only columns in first-seen order
        for k in r:
            if k not in header:
                header.append(k)

    write_tsv(OUT / "positives_all.tsv", header, rows)

    print(f"validated:  {sum(1 for r in tiered if r['evidence_tier'] == 'validated')}")
    print(
        f"predicted:  {sum(1 for r in tiered if r['evidence_tier'] == 'predicted')}   "
        f"instance_source={dict(Counter(r['instance_source'] for r in tiered if r['evidence_tier'] == 'predicted'))}"
    )
    print(
        f"T5SS folded: {len(t5)}   {by_type(t5, 'subtype')}   "
        f"self_secreted=true: {sum(1 for r in t5 if r.get('self_secreted') == 'true')}"
    )
    print(f"total positives_all: {len(rows)}   {by_type(rows)}")
    print(f"  instance_source: {dict(sorted(Counter(r.get('instance_source', '') for r in rows).items()))}")
    print(f"  type_level=yes:  {sum(1 for r in rows if r.get('type_level') == 'yes')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
