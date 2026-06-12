#!/usr/bin/env python3
"""Dataset group 6, task 6.3: apply the pass-2 citation audit to positives_all.tsv.

The pass-2 re-audit (task 6.2, `pass2_results.tsv`) produced four kinds of decision for the 30
ssign-found effectors it re-sourced. This script applies them to the training-label table,
deterministically and reversibly:

  RESOLVED     -> overwrite the row's sourcing DOI with the verified corrected DOI.
  CONFIRMED    -> no change (the existing DOI was already correct); logged for the record.
  MISASSIGNED  -> the protein is real but the SS-type label is wrong (BopA/BopE are Bsa T3SS, not
                  T6SS). Re-label ss_type, drop the now-invalid T6SS instance assignment (-> type-level
                  T3SS positive), and correct the DOI. The label, not the row, was the defect.
  NOT_FOUND    -> no primary paper supports the row as a secreted effector of this system; REMOVE it
                  from the training positives (quarantined to positives_removed_citation.tsv, never
                  silently deleted).

Plus one curated de-duplication: idx 24 (TplE_alias_Tle4) == idx 15 (Tle4) == P. aeruginosa PA1510;
drop the alias row.

Targets are located in positives_all.tsv by the exact (sys_instance_id, gene, locus_tag) identity
that script 41 recorded when it matched each found effector (pos_* columns), so only the audited
instance is touched, never another instance of the same gene. Verified 30/30 unique matches.

Inputs : data/dataset/citation_consistency_found.tsv  (pos_* identity + old sourcing_doi)
         data/dataset/pass2_results.tsv                (per-idx decision + final_doi)
         data/dataset/positives_all.tsv                (the table to correct)
Outputs: data/dataset/positives_all.tsv                (corrected, in place)
         data/dataset/positives_all.pre_citation_audit.tsv   (backup of the original)
         data/dataset/positives_removed_citation.tsv   (the removed rows + reason)
         data/dataset/citation_corrections_log.tsv     (every field change)
Run:     .venv/bin/python scripts/43_apply_citation_corrections.py
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from bench_io import norm_doi, read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "data" / "dataset"
POS = DATASET / "positives_all.tsv"

# Curated de-duplication: the alias row to drop -> the row it duplicates (both = PA1510).
DUPLICATE_IDX = {"24": "duplicate_of PA1510 / Tle4 (idx 15)"}


def main() -> int:
    cc = {(r["gene"], r["ss_type"], r["organism"]): r for r in read_tsv(DATASET / "citation_consistency_found.tsv")}
    decisions = read_tsv(DATASET / "pass2_results.tsv")
    pos = read_tsv(POS)
    pos_header = list(pos[0].keys())

    # Locate each audited effector's exact positives_all row index. Key on (sys_instance_id, gene) —
    # the same instance key script 41 resolved every found effector by — NOT on locus_tag: some rows
    # have a legitimately blank locus_tag (e.g. prtB, placed by coordinates) and gating on it would
    # drop the match. locus_tag / real UniProt only break a tie if (instance, gene) is non-unique.
    def find_row(gene, uni, locus, inst):
        uni = uni if uni and uni != "-" else ""
        cands = [
            i
            for i, p in enumerate(pos)
            if (p.get("gene") or "").strip() == gene and (p.get("sys_instance_id") or "").strip() == inst
        ]
        if len(cands) <= 1:
            return cands[0] if cands else None
        for i in cands:  # tie-break only when (instance, gene) is ambiguous
            p = pos[i]
            if (locus and (p.get("locus_tag") or "").strip() == locus) or (
                uni and (p.get("uniprot") or "").strip() == uni
            ):
                return i
        return None  # ambiguous and untie-breakable -> skip (logged), never edit a guessed row

    log = []
    remove_idx = {}  # positives_all row index -> reason
    for d in decisions:
        c = cc.get((d["gene"], d["ss_type"], d["organism"]))
        if not c:
            log.append(
                {"idx": d["idx"], "gene": d["gene"], "action": "SKIP", "field": "", "old": "", "new": "no cc identity"}
            )
            continue
        inst, puni, plocus = c["pos_sys_instance_id"].strip(), c["pos_uniprot"].strip(), c["pos_locus_tag"].strip()
        old_doi = norm_doi(c["sourcing_doi"])
        ri = find_row(d["gene"], puni, plocus, inst)
        if ri is None:
            log.append(
                {"idx": d["idx"], "gene": d["gene"], "action": "SKIP", "field": "", "old": "", "new": "row not found"}
            )
            continue
        row = pos[ri]
        status = d["status"]

        if d["idx"] in DUPLICATE_IDX:
            remove_idx[ri] = DUPLICATE_IDX[d["idx"]]
            log.append(
                {
                    "idx": d["idx"],
                    "gene": d["gene"],
                    "action": "REMOVE",
                    "field": "",
                    "old": "",
                    "new": DUPLICATE_IDX[d["idx"]],
                }
            )
            continue
        if status == "NOT_FOUND":
            reason = "unsupported_no_primary_paper"
            remove_idx[ri] = reason
            log.append({"idx": d["idx"], "gene": d["gene"], "action": "REMOVE", "field": "", "old": "", "new": reason})
            continue

        def set_field(field, new):
            old = row.get(field, "")
            if (old or "") != (new or ""):
                row[field] = new
                log.append(
                    {"idx": d["idx"], "gene": d["gene"], "action": status, "field": field, "old": old, "new": new}
                )

        new_doi = norm_doi(d["final_doi"])
        if status in ("RESOLVED", "MISASSIGNED") and new_doi:
            # Update whichever column actually held the wrong DOI (primary_ref by default).
            doi_field = (
                "instance_source_doi" if norm_doi(row.get("instance_source_doi", "")) == old_doi else "primary_ref"
            )
            set_field(doi_field, new_doi)
        if status == "MISASSIGNED" and d["corrected_ss_type"]:
            set_field("ss_type", d["corrected_ss_type"])
            set_field("sys_instance_id", "")  # the old T6SS instance no longer applies
            set_field("type_level", "yes")
            set_field("instance_source", "audit_reassign")
        # CONFIRMED -> record, no field change.
        if status == "CONFIRMED":
            log.append(
                {
                    "idx": d["idx"],
                    "gene": d["gene"],
                    "action": "CONFIRMED",
                    "field": "",
                    "old": "",
                    "new": norm_doi(row.get("primary_ref", "")),
                }
            )

    removed = [{**pos[i], "removal_reason": r} for i, r in sorted(remove_idx.items())]
    kept = [p for i, p in enumerate(pos) if i not in remove_idx]

    # Snapshot the pristine table once; never overwrite it on a re-run (the script is idempotent:
    # removals are gone and DOI edits are no-ops the second time, so re-running must not lose the original).
    backup = DATASET / "positives_all.pre_citation_audit.tsv"
    if not backup.exists():
        backup.write_text(POS.read_text())
    write_tsv(POS, pos_header, kept)
    write_tsv(DATASET / "positives_removed_citation.tsv", pos_header + ["removal_reason"], removed)
    write_tsv(DATASET / "citation_corrections_log.tsv", ["idx", "gene", "action", "field", "old", "new"], log)

    print(
        f"positives_all.tsv: {len(pos)} -> {len(kept)} rows  ({len(removed)} removed, backup = positives_all.pre_citation_audit.tsv)"
    )
    print("decision tally   :", dict(Counter(d["status"] for d in decisions)))
    print("removed          :", [(r["gene"], r["removal_reason"]) for r in removed])
    print("field edits      :", dict(Counter(l["field"] for l in log if l["field"])))
    print(f"full change log  : data/dataset/citation_corrections_log.tsv ({len(log)} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
