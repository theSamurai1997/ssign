#!/usr/bin/env python3
"""Apply the recall-figure answer-key audit corrections (task #74-79).

The 4-pass audit + fresh blind-agent fix-finding produced
`data/phase2/verification/reconciled_fixes.tsv` (101 rows). Every replacement DOI was
Crossref/PubMed-confirmed this run; every UniProt accession verified live. This script applies
those corrections to the source answer-key tables, deterministically and reversibly:

  effector_gold_set.tsv   -> T1/T2/T3/T4/T6 rows (primary_ref, uniprot, locus_tag, gene)
  positives_all.tsv       -> T5SS rows (primary_ref [+ instance_source_doi if it held the old DOI])

Decisions baked in (from the user, 2026-06-16):
  DROP (quarantined, not deleted): R207 BopC, R208 BopE (misclassified T6SS dups; correct T3SS rows
        already exist), R214 YPK_3548 (duplicate of R215/YezP), R179+R186 TcpB/Btp1 (no VirB
        translocation assay exists; T4SS-dependence unestablished).
  R159: gene MavF -> SdbC, locus lpg2391 kept (lpg2391 IS sdbC, a Dot/Icm substrate; "MavF" was a
        mislabel), uniprot -> Q5ZSX5, primary_ref -> Huang 2011 (10.1111/j.1462-5822.2010.01531.x).
  keep_current: R060 lip, R066 zmpA (no better primary exists).

Each row_id is located in the source table by the exact (ss_type, gene, +locus_tag/uniprot/organism)
identity recorded in the verification batch files (dry-run confirmed 101/101 unique matches).

Inputs : data/phase2/verification/reconciled_fixes.tsv   (the corrections)
         data/phase2/verification/batch_*.tsv            (row_id -> source identity)
         data/gold_build/effector_gold_set.tsv
         data/dataset/positives_all.tsv
Outputs: the two source tables (corrected in place)
         *.pre_audit_fix.tsv                              (one-time pristine backups)
         data/gold_build/effector_gold_set.removed_audit.tsv  (quarantined drops + reason)
         data/phase2/verification/audit_fix_changelog.tsv (every field change)
Run:     .venv/bin/python scripts/55_apply_audit_corrections.py
"""

from __future__ import annotations

import csv
import glob
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "data" / "gold_build" / "effector_gold_set.tsv"
POS = ROOT / "data" / "dataset" / "positives_all.tsv"
VER = ROOT / "data" / "phase2" / "verification"

DROP = {"R207", "R208", "R214", "R179", "R186"}
KEEP = {"R060", "R066"}
DROP_REASON = {
    "R207": "misclassified T6SS dup; correct T3SS/Bsa BopC already in gold",
    "R208": "misclassified T6SS dup; correct T3SS/Bsa BopE already in gold",
    "R214": "duplicate of R215 (YezP/YPK_3549); YPK_3548 has no secretion evidence",
    "R179": "TcpB/Btp1: no VirB-translocation assay exists; T4SS-dependence unestablished",
    "R186": "TcpB/Btp1 (B. abortus): no VirB-translocation assay exists; T4SS-dependence unestablished",
}
# R159 special-case overrides (gene rename + accession + DOI)
R159 = {"gene": "SdbC", "uniprot": "Q5ZSX5", "primary_ref": "10.1111/j.1462-5822.2010.01531.x"}


def rd(p):
    with open(p) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return rows, list(rows[0].keys())


def norm(s):
    return (s or "").strip()


def locate(rows, ssype, gene, locus, uni, organism):
    """Return the unique source-row index for this identity, or None if not unique."""
    cands = [i for i, r in enumerate(rows) if norm(r.get("ss_type")) == ssype and norm(r.get("gene")) == gene]
    if len(cands) == 1:
        return cands[0]
    if locus:
        c = [i for i in cands if norm(rows[i].get("locus_tag")) == locus]
        if len(c) == 1:
            return c[0]
    if uni and uni != "-":
        c = [i for i in cands if norm(rows[i].get("uniprot")) == uni]
        if len(c) == 1:
            return c[0]
    if organism:
        tok = organism.split()[0].rstrip(".").lower()
        c = [i for i in cands if tok in norm(rows[i].get("organism")).lower()]
        if len(c) == 1:
            return c[0]
    return None


def main() -> int:
    gold, gold_hdr = rd(GOLD)
    pos, pos_hdr = rd(POS)
    ident = {}
    for bf in glob.glob(str(VER / "batch_*.tsv")):
        for r in csv.DictReader(open(bf), delimiter="\t"):
            ident[r["row_id"]] = r
    fixes, _ = rd(VER / "reconciled_fixes.tsv")

    log = []
    drop_gold = {}  # gold index -> reason

    def set_field(table, ri, field, new, rid, gene):
        old = table[ri].get(field, "")
        if (old or "") != (new or ""):
            table[ri][field] = new
            log.append({"row_id": rid, "gene": gene, "field": field, "old": old, "new": new})

    for fx in fixes:
        rid = fx["row_id"]
        idr = ident[rid]
        ssype, gene = norm(idr.get("ss_type")), norm(idr.get("gene"))
        locus, uni, org = norm(idr.get("locus_tag")), norm(idr.get("uniprot")), norm(idr.get("organism"))
        table, _ = (pos, "positives") if ssype == "T5SS" else (gold, "gold")
        ri = locate(table, ssype, gene, locus, uni, org)
        if ri is None:
            log.append({"row_id": rid, "gene": gene, "field": "", "old": "", "new": "ROW NOT FOUND - SKIPPED"})
            continue

        if rid in DROP:
            drop_gold[ri] = DROP_REASON[rid]  # all 5 drops are gold rows (T4/T6)
            log.append({"row_id": rid, "gene": gene, "field": "", "old": "", "new": f"DROP: {DROP_REASON[rid]}"})
            continue
        if rid in KEEP:
            log.append({"row_id": rid, "gene": gene, "field": "", "old": "", "new": "keep_current (no better primary)"})
            continue

        if rid == "R159":
            set_field(table, ri, "gene", R159["gene"], rid, gene)
            set_field(table, ri, "uniprot", R159["uniprot"], rid, gene)
            set_field(table, ri, "primary_ref", R159["primary_ref"], rid, gene)
            continue

        new_doi = norm(fx.get("new_primary_ref"))
        new_uni = norm(fx.get("new_uniprot"))
        blank_uni = norm(fx.get("set_uniprot_blank")) == "yes"
        new_loc = norm(fx.get("new_locus_tag"))

        if new_doi:
            # mirror script 43: if instance_source_doi held the old primary, update it too
            old_doi = norm(table[ri].get("primary_ref"))
            set_field(table, ri, "primary_ref", new_doi, rid, gene)
            if "instance_source_doi" in table[ri] and norm(table[ri].get("instance_source_doi")) == old_doi and old_doi:
                set_field(table, ri, "instance_source_doi", new_doi, rid, gene)
        if blank_uni:
            set_field(table, ri, "uniprot", "", rid, gene)
        elif new_uni:
            set_field(table, ri, "uniprot", new_uni, rid, gene)
        if new_loc:
            set_field(table, ri, "locus_tag", new_loc, rid, gene)

    # quarantine drops
    removed = [{**gold[i], "removal_reason": r} for i, r in sorted(drop_gold.items())]
    gold_kept = [r for i, r in enumerate(gold) if i not in drop_gold]

    # one-time pristine backups (idempotent: never overwrite)
    for src, bak in [(GOLD, GOLD.with_suffix(".pre_audit_fix.tsv")), (POS, POS.with_suffix(".pre_audit_fix.tsv"))]:
        if not bak.exists():
            bak.write_text(src.read_text())

    def write_tsv(p, hdr, rows):
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    write_tsv(GOLD, gold_hdr, gold_kept)
    write_tsv(POS, pos_hdr, pos)
    write_tsv(GOLD.parent / "effector_gold_set.removed_audit.tsv", gold_hdr + ["removal_reason"], removed)
    write_tsv(VER / "audit_fix_changelog.tsv", ["row_id", "gene", "field", "old", "new"], log)

    print(f"gold: {len(gold)} -> {len(gold_kept)} rows ({len(removed)} quarantined)")
    print(f"positives: {len(pos)} rows (edited in place)")
    print("field edits:", dict(Counter(l["field"] for l in log if l["field"])))
    print("drops:", [(r["gene"], r["ss_type"], r["removal_reason"][:40]) for r in removed])
    print(f"changelog: data/phase2/verification/audit_fix_changelog.tsv ({len(log)} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
