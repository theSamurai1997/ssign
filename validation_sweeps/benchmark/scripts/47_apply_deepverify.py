#!/usr/bin/env python3
"""Apply the full-table citation audit to positives_all.tsv (strict policy: drop every refuted row).

Two removal sources:
  - pass-1 deterministic drops (FLAG_WRONG_TOPIC / FLAG_GENE_ABSENT / DOI_UNRESOLVED, and FETCH_ERROR
    if any): 252 rows whose cited paper is wrong-topic, doesn't name the gene, or has a dead DOI.
    FETCH_ERROR means CrossRef was unreachable after retries (transient); 0 occurred this run because
    the caches were warm, but if a future run has any, re-run script 44 to settle them before applying.
  - deep-verify REFUTED (wrong_organism / wrong_protein / wrong_system / no_effector_evidence): 215 rows
    where reading the actual paper gave positive counter-evidence the citation does not support the claim.
Survivors that read clean keep their row and gain two columns:
  citation_trust  - verified_paper | verified_external | unverifiable | fallback_consistent
  citation_quote  - the verbatim sentence from the cited paper that supports the claim (verified rows)

Row identity uses index alignment with citation_consistency_full.tsv (script 44 emitted it 1:1 in
positives order); the deep-verify trust tier is keyed by (gene, sys_instance_id, organism).

Backup  : data/dataset/positives_all.pre_deepverify.tsv  (written once)
Log     : data/dataset/deepverify_removed.tsv             (every removed row + reason)
Run     : python3 scripts/47_apply_deepverify.py
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "data" / "dataset"
sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

POS = DATASET / "positives_all.tsv"
BACKUP = DATASET / "positives_all.pre_deepverify.tsv"
DROP_PASS1 = {"FLAG_WRONG_TOPIC", "FLAG_GENE_ABSENT", "DOI_UNRESOLVED", "FETCH_ERROR"}


def key(r: dict) -> tuple:
    return ((r.get("gene") or "").strip(), (r.get("sys_instance_id") or "").strip(), (r.get("organism") or "").strip())


def main() -> int:
    pos = read_tsv(POS)
    c44 = read_tsv(DATASET / "citation_consistency_full.tsv")
    dv = {key(r): r for r in read_tsv(DATASET / "deepverify_results_full.tsv")}
    if len(pos) != len(c44):
        print(f"ABORT: positives_all ({len(pos)}) and citation_consistency_full ({len(c44)}) differ in length")
        return 1

    kept, removed = [], []
    for prow, crow in zip(pos, c44):
        if prow.get("gene", "").strip() != crow.get("gene", "").strip():
            print(f"ABORT: row misalignment at gene {prow.get('gene')} vs {crow.get('gene')}")
            return 1
        v1 = crow["verdict"]
        if v1 in DROP_PASS1:
            removed.append({**_log(prow), "stage": "pass1", "reason": v1})
            continue
        d = dv.get(key(prow))
        tier = d["trust_tier"] if d else "fallback_unknown"
        if tier.startswith("refuted_"):
            removed.append({**_log(prow), "stage": "deepverify", "reason": tier})
            continue
        prow["citation_trust"] = tier
        prow["citation_quote"] = (d or {}).get("quote", "")
        kept.append(prow)

    if not kept:
        print("ABORT: every row was removed; refusing to overwrite positives_all")
        return 1
    if not BACKUP.exists():
        BACKUP.write_text(POS.read_text())
    cols = list(kept[0].keys())  # original cols + the two appended
    write_tsv(POS, cols, kept)
    write_tsv(DATASET / "deepverify_removed.tsv", list(removed[0].keys()), removed)

    print(f"positives_all: {len(pos)} -> {len(kept)} kept, {len(removed)} removed")
    print(f"  removed by stage: {dict(Counter(r['stage'] for r in removed))}")
    print(f"  pass1 reasons:    {dict(Counter(r['reason'] for r in removed if r['stage'] == 'pass1'))}")
    print(f"  deepverify reasons: {dict(Counter(r['reason'] for r in removed if r['stage'] == 'deepverify'))}")
    print("\nkept by citation_trust:")
    for t, c in Counter(r["citation_trust"] for r in kept).most_common():
        print(f"  {c:4d}  {t}")
    print("\nkept by ss_type:")
    for ss, c in sorted(Counter(r["ss_type"] for r in kept).items()):
        print(f"  {ss}: {c}")
    return 0


def _log(prow: dict) -> dict:
    return {
        k: prow.get(k, "")
        for k in (
            "gene",
            "uniprot",
            "locus_tag",
            "sys_instance_id",
            "ss_type",
            "organism",
            "evidence_tier",
            "primary_ref",
        )
    }


if __name__ == "__main__":
    raise SystemExit(main())
