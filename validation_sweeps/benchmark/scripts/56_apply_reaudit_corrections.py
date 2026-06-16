#!/usr/bin/env python3
"""Apply the re-audit (acceptance-test) corrections to the answer key -- round 2.

After script 55 applied the first audit pass, the corrected recall-figure proteins (240) were
re-audited under a STRICT same-species-secretion-evidence bar (manual main-session pass; the blind
agents tripped a content filter on the pathogen batches). Every DOI was resolved live this run
(crossref/europepmc) and every UniProt accession/locus verified against the live index. The verdicts
+ verified replacement refs live in `data/phase2/verification/reaudit/` (reconciled_*.tsv) and are
distilled into the machine-readable `reaudit/fixes_verified.tsv`, which drives this script.

This is round 2 -- it edits the ALREADY-round-1-corrected gold set, reversibly and independently of
round 1:

  action=drop          -> quarantine the gold row (9 rows: evidence too weak / structural / incoherent)
  action=fix           -> set primary_ref / uniprot / locus_tag on the gold row (18 field edits)
  action=dedup_figure  -> NOT handled here (it is a figure-build instance-collapse done in script 54;
                          listed in fixes_verified.tsv only for the record). T3013/T3025.

Rows are located in the gold set by genome identity (gene + locus_tag, then uniprot, then organism).
ss_type is deliberately NOT used to pre-filter: the verification batch files mislabel one row's
ss_type (T3005 BopA is tagged T3SS in the batch but is a T6SS row in the gold set), which would
mis-target a same-gene row. Locus/accession identity is authoritative. A located row whose organism
disagrees with the batch is reported and skipped rather than edited.

Inputs : data/phase2/verification/reaudit/fixes_verified.tsv   (corrections, source of truth)
         data/phase2/verification/reaudit/batch_*.tsv          (row_id -> source identity)
         data/gold_build/effector_gold_set.tsv                 (round-1-corrected)
Outputs: data/gold_build/effector_gold_set.tsv                 (round-2-corrected, in place)
         data/gold_build/effector_gold_set.pre_reaudit2.tsv    (one-time round-2 pristine backup)
         data/gold_build/effector_gold_set.removed_reaudit2.tsv (round-2 drops only; script 54 reads
                                                                  this alongside round-1's removed_audit)
         data/phase2/verification/audit_fix_changelog_reaudit2.tsv (every round-2 field change)
Run:     .venv/bin/python scripts/56_apply_reaudit_corrections.py
"""

from __future__ import annotations

import csv
import glob
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "data" / "gold_build" / "effector_gold_set.tsv"
VER = ROOT / "data" / "phase2" / "verification"
REAUDIT = VER / "reaudit"
FIXES = REAUDIT / "fixes_verified.tsv"
REMOVED2 = GOLD.parent / "effector_gold_set.removed_reaudit2.tsv"
BACKUP = GOLD.parent / "effector_gold_set.pre_reaudit2.tsv"
CHANGELOG = VER / "audit_fix_changelog_reaudit2.tsv"


def rd(p):
    with open(p) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return rows, list(rows[0].keys())


def norm(s):
    return (s or "").strip()


def org_tok(organism):
    return organism.split()[0].rstrip(".").lower() if organism else ""


def locate(rows, gene, locus, uni, organism):
    """Genome-identity match: gene + (locus_tag | uniprot | organism). ss_type is NOT used."""
    by_gene = [i for i, r in enumerate(rows) if norm(r.get("gene")) == gene]
    for key, val in (("locus_tag", locus), ("uniprot", uni)):
        if val and val != "-":
            c = [i for i in by_gene if norm(rows[i].get(key)) == val]
            if len(c) == 1:
                return c[0]
            if len(c) > 1:
                return None  # ambiguous on a strong key -> refuse
    if organism:
        c = [i for i in by_gene if org_tok(organism) in norm(rows[i].get("organism")).lower()]
        if len(c) == 1:
            return c[0]
    return by_gene[0] if len(by_gene) == 1 else None


def main() -> int:
    gold, gold_hdr = rd(GOLD)
    ident = {}
    for bf in glob.glob(str(REAUDIT / "batch_*.tsv")):
        for r in rd(bf)[0]:
            ident[r["row_id"]] = r
    fixes, _ = rd(FIXES)

    # group corrections by row_id so a multi-field row is located once (by its ORIGINAL identity)
    by_row = defaultdict(list)
    for fx in fixes:
        by_row[fx["row_id"]].append(fx)

    log = []
    drop_gold = {}  # gold index -> reason

    for rid, group in by_row.items():
        action = norm(group[0]["action"])
        if action == "dedup_figure":
            continue  # handled in script 54
        idr = ident[rid]
        gene, locus = norm(idr.get("gene")), norm(idr.get("locus_tag"))
        uni, org = norm(idr.get("uniprot")), norm(idr.get("organism"))
        ri = locate(gold, gene, locus, uni, org)
        if ri is None:
            log.append(
                {"row_id": rid, "gene": gene, "field": "", "old": "", "new": "ROW NOT FOUND/AMBIGUOUS - SKIPPED"}
            )
            continue
        if org and org_tok(org) not in norm(gold[ri].get("organism")).lower():
            log.append(
                {
                    "row_id": rid,
                    "gene": gene,
                    "field": "",
                    "old": norm(gold[ri].get("organism")),
                    "new": f"ORGANISM MISMATCH (batch={org}) - SKIPPED",
                }
            )
            continue

        if action == "drop":
            drop_gold[ri] = norm(group[0]["note"])
            log.append({"row_id": rid, "gene": gene, "field": "", "old": "", "new": "DROP: " + norm(group[0]["note"])})
        else:  # fix (possibly multiple fields)
            for fx in group:
                field, new = norm(fx["field"]), norm(fx["new_value"])
                old = gold[ri].get(field, "")
                if (old or "") != (new or ""):
                    gold[ri][field] = new
                    log.append({"row_id": rid, "gene": gene, "field": field, "old": old, "new": new})

    if not BACKUP.exists():
        BACKUP.write_text(GOLD.read_text())

    removed2 = [{**gold[i], "removal_reason": r} for i, r in sorted(drop_gold.items())]
    gold_kept = [r for i, r in enumerate(gold) if i not in drop_gold]

    def write_tsv(p, hdr, rows):
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=hdr, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    write_tsv(REMOVED2, gold_hdr + ["removal_reason"], removed2)
    write_tsv(GOLD, gold_hdr, gold_kept)
    write_tsv(CHANGELOG, ["row_id", "gene", "field", "old", "new"], log)

    skipped = [l for l in log if "SKIPPED" in l["new"]]
    print(f"gold: {len(gold)} -> {len(gold_kept)} rows ({len(removed2)} quarantined this round)")
    print("field edits:", dict(Counter(l["field"] for l in log if l["field"])))
    print("drops:", [(r["gene"], r["ss_type"], r["locus_tag"], r["organism"][:22]) for r in removed2])
    if skipped:
        print("!! SKIPPED:", [(l["row_id"], l["new"]) for l in skipped])
    print(f"changelog: {CHANGELOG.relative_to(ROOT)} ({len(log)} entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
