#!/usr/bin/env python3
"""Dataset group 1, task 1.2: verify the predicted rows to the gold-set bar.

Holds the predicted tier to the same standard scripts 02/03 applied to the validated
gold set, reusing their network code verbatim (single source of truth):

  - VERIFIED-clean rows (171) pass straight through; the 2026-06-08 audit already
    confirmed their DOI resolves and the paper names the protein.
  - PARTIAL repair-queue rows (154): rebuild the UniProt accession from the trusted
    locus_tag (script 02 `rebuild_uniprot`) and lift an explicitly-named corrected DOI
    from the audit notes (script 02 `CORRECT_DOI`). Never guess a DOI.
  - Hard drop only what the validated pipeline hard-dropped: unplaceable rows (no usable
    locus_tag AND no usable UniProt accession). Phase-1/feature placement needs one.
  - Independent DOI.org Handle recheck over every kept row (script 03 `doi_resolves`),
    recording citation_status. As in the validated path, a broken DOI is REPORTED, not
    dropped: the protein/locus/instance label is load-bearing, the citation is metadata.
    Predicted rows are the down-weighted tier; unresolved citations are flagged for an
    optional re-citation pass, not silently lost.

Inputs (from data/dataset/, produced by 30_extract_predicted.py):
  30_predicted_verified_clean.tsv, 30_predicted_repair_queue.tsv
Outputs (to data/dataset/):
  predicted_audited.tsv             kept rows + uniprot_status, citation_status, audit_tier
  predicted_audit_provenance.tsv    every input row + fate + fate_reason
  .doi_cache.json                   DOI.org resolution cache (re-runs are offline)
"""

from __future__ import annotations

import importlib.util
import json
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "dataset"
DOI_CACHE = OUT / ".doi_cache.json"


def _load(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).parent / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rg = _load("repair_gold_set", "02_repair_gold_set.py")
fz = _load("finalize_gold_set", "03_finalize_gold_set.py")
from bench_io import by_type, norm_doi, read_tsv, write_tsv  # noqa: E402  (scripts/ is on sys.path)

KEEP_COLS = [
    "gene",
    "uniprot",
    "locus_tag",
    "organism",
    "refseq_genome",
    "ss_type",
    "sys_instance_id",
    "evidence_level",
    "primary_ref",
    "verification_status",
    "verification_notes",
]
AUDIT_COLS = ["uniprot_status", "citation_status", "audit_tier"]


def main() -> int:
    verified = read_tsv(OUT / "30_predicted_verified_clean.tsv")
    partial = read_tsv(OUT / "30_predicted_repair_queue.tsv")
    uni_cache = rg.load_cache()

    kept: list[dict] = []
    provenance: list[dict] = []

    def consider(r: dict, audit_tier: str) -> None:
        locus = (r.get("locus_tag") or "").strip()
        old_uni = (r.get("uniprot") or "").strip()
        notes = r.get("verification_notes") or ""
        r = dict(r)
        r["audit_tier"] = audit_tier

        if audit_tier == "PARTIAL":
            uni_flagged = "uniprot" in notes.lower()
            rebuilt = rg.rebuild_uniprot(locus, uni_cache) if (uni_flagged or not rg.has_locus(old_uni)) else None
            if uni_flagged:
                r["uniprot_status"] = "CORRECTED" if rebuilt else "UNRESOLVED"
            elif not rg.has_locus(old_uni):
                r["uniprot_status"] = "FILLED" if rebuilt else "MISSING"
            else:
                r["uniprot_status"] = "KEPT"
            if rebuilt:
                r["uniprot"] = rebuilt
            m = rg.CORRECT_DOI.search(notes)
            if m:
                r["primary_ref"] = m.group(1).rstrip(".;,)")
        else:
            r["uniprot_status"] = "VERIFIED_AUDIT"

        unplaceable = not rg.has_locus(r.get("uniprot", "")) and not rg.has_locus(locus)
        if unplaceable:
            provenance.append(
                {
                    **r,
                    "fate": "DROP_UNPLACEABLE",
                    "fate_reason": "no usable locus_tag and no resolvable UniProt accession",
                }
            )
            return
        kept.append(r)
        provenance.append({**r, "fate": "KEEP", "fate_reason": f"audit_tier={audit_tier}"})

    for r in verified:
        consider(r, "VERIFIED")
    for r in partial:
        consider(r, "PARTIAL")

    rg.save_cache(uni_cache)

    # Fold in agent re-citations (step 34, Crossref-verified) before the recheck, so a
    # corrected DOI replaces the broken one in primary_ref. Idempotent: re-runnable.
    recite_path = OUT / "doi_recite.jsonl"
    recite: dict[str, str] = {}
    if recite_path.exists():
        for line in recite_path.read_text().splitlines():
            if line.strip():
                rec = json.loads(line)
                recite[norm_doi(rec["broken_doi"])] = rec["corrected_doi"].strip()
    if recite:
        n = 0
        for r in kept + provenance:
            cur = norm_doi(r.get("primary_ref", ""))
            if cur in recite:
                r["primary_ref"] = recite[cur]
                n += 1
        print(f"applied {len(recite)} re-citations to {n} rows")

    # Independent DOI.org recheck over kept rows (report, never drop).
    doi_cache: dict = json.loads(DOI_CACHE.read_text()) if DOI_CACHE.exists() else {}
    distinct = sorted({(r.get("primary_ref") or "").strip() for r in kept if (r.get("primary_ref") or "").strip()})
    print(f"DOI.org recheck of {len(distinct)} distinct kept DOIs ...")
    for d in distinct:
        fz.doi_resolves(d, doi_cache)
        time.sleep(0.05)
    DOI_CACHE.write_text(json.dumps(doi_cache, indent=0))

    def citation_status(raw: str) -> str:
        # doi_resolves caches by the cleaned DOI, so normalize the same way or a
        # "doi:"-prefixed ref would miss the cache and false-fail as UNRESOLVED.
        doi = norm_doi(raw)
        if not doi:
            return "NONE"
        return "RESOLVED" if doi_cache.get(doi) else "UNRESOLVED"

    for r in kept:
        r["citation_status"] = citation_status(r.get("primary_ref", ""))
    for p in provenance:  # mirror citation_status onto kept provenance rows
        if p["fate"] == "KEEP":
            p["citation_status"] = citation_status(p.get("primary_ref", ""))

    write_tsv(OUT / "predicted_audited.tsv", KEEP_COLS + AUDIT_COLS, kept)
    write_tsv(OUT / "predicted_audit_provenance.tsv", KEEP_COLS + AUDIT_COLS + ["fate", "fate_reason"], provenance)

    fates = Counter(p["fate"] for p in provenance)
    cites = Counter(r["citation_status"] for r in kept)
    print(f"\ninputs: {len(verified)} verified-clean + {len(partial)} partial = {len(verified) + len(partial)}")
    print(f"kept:   {len(kept)}   {by_type(kept)}")
    print("fates:", dict(fates))
    print("citation_status (kept):", dict(cites))
    unresolved = [r for r in kept if r["citation_status"] == "UNRESOLVED"]
    print(f"  UNRESOLVED DOIs (flagged, kept): {len(unresolved)}   {by_type(unresolved)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
