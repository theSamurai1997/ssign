#!/usr/bin/env python3
"""Phase 0a, steps 2.5(apply)-2.7: merge repairs and emit the corpus-derived gold set.

Inputs (all under data/gold_build/):
  03_verified_clean.tsv     349 validated+VERIFIED effectors, clean (pass straight through)
  04_repaired.tsv           183 PARTIAL rows with UniProt rebuilt (uniprot_status, doi_status)
  doi_repair_pilot.jsonl    5  DOI-repair records (gene, uniprot, status, corrected_doi, quote)
  doi_repair_workflow.jsonl 49 DOI-repair records
  doi_repair_manual.jsonl   5  DOI-repair records (filter-blocked rows done by hand)

What it does:
  - Applies the DOI repairs to the 183 PARTIAL rows (corrected DOI -> primary_ref).
  - Applies Teo's adjudications on contradicted/identity-error rows:
      DROP_CONTRADICTED : celY, pemB, paeX, CagF, Sca4 (real/secreted but NOT a substrate
                          of the labelled system, route unresolved) and Lem4 (wrong protein).
      RECLASSIFY_T1SS   : TRP47, TRP32 (experimentally shown to be T1SS, not T4SS substrates).
      KEEP (T4SS)       : TcpB/Btp1 x2 (host-translocated; Brucella's only route is VirB).
  - Drops rows that cannot be placed in Phase 1 (no usable locus_tag AND no usable UniProt).
  - 2.6 independent check: re-confirms every corrected DOI resolves via Crossref (does not
    trust the agent's word), and records DOIs that fail.
  - 2.7 emits gold_set_corpus.tsv (the gold set) and gold_set_provenance.tsv (every row's fate).

Decision rows are keyed by UniProt accession (unique per protein), never by gene name.
"""

from __future__ import annotations

import csv
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "gold_build"
UA = "ssign-benchmark/0.1 (teoreid@gmail.com)"

# --- Teo's adjudications (2026-06-10), keyed by UniProt accession ----------
DROP_CONTRADICTED = {
    "P27032": "celY: Out/T2SS-independent (Zhou & Ingram 2000, 'not affected by out genes')",
    "Q47474": "pemB: outer-membrane lipoprotein, 'not released into the extracellular medium'",
    "E0SDD1": "paeX: periplasmic, 'cannot be considered true Out-dependent secretion'",
    "O25276": "CagF: a T4SS chaperone for CagA, not itself translocated (machinery, not cargo)",
    "Q52658": "Sca4: secretion proven but route unresolved (paper: 'T1SS or T4SS'); T4SS unsupported",
    "Q5ZX05": "Lem4: stored accession is PilN (type IV pilus machinery), not the Lem4 effector",
}
RECLASSIFY_T1SS = {
    "Q2GHU2": "TRP47: negative in T4SS assay, TolC-dependent T1SS substrate (Wakeel 2011)",
    "Q2GHT8": "TRP32: negative in T4SS assay, TolC-dependent T1SS substrate (Wakeel 2011)",
}

PASS_HEADER = [
    "gene",
    "uniprot",
    "locus_tag",
    "organism",
    "refseq_genome",
    "ss_type",
    "sys_instance_id",
    "evidence_level",
    "primary_ref",
    "family",
    "length",
    "proteome",
    "verification_status",
    "verification_notes",
]


def has_val(s: str) -> bool:
    return bool((s or "").strip()) and (s or "").strip() != "-"


def read_tsv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_doi_records() -> dict:
    """Map UniProt accession -> DOI-repair record, across all three passes."""
    by_uni: dict[str, dict] = {}
    for name in ("doi_repair_pilot.jsonl", "doi_repair_workflow.jsonl", "doi_repair_manual.jsonl"):
        p = OUT / name
        if not p.exists():
            continue
        for line in p.read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if r.get("uniprot"):
                by_uni[r["uniprot"].strip()] = r
    return by_uni


def doi_resolves(doi: str, cache: dict) -> bool:
    """True if the DOI is registered, via the authoritative DOI.org Handle API.

    This tests resolution itself (any registration agency), not Crossref coverage, so
    it does not false-fail on older or non-Crossref DOIs the way /works does. The Handle
    API returns responseCode 1 for a known handle, 100 for not-found.
    """
    doi = doi.strip().removeprefix("doi:").strip()
    if not doi:
        return False
    if doi in cache:
        return cache[doi]
    url = "https://doi.org/api/handles/" + urllib.parse.quote(doi, safe="")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    ok = False
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.load(resp)
            ok = data.get("responseCode") == 1
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                ok = False
                break
            time.sleep(2 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            time.sleep(2 * (attempt + 1))
    cache[doi] = ok
    return ok


def main() -> int:
    verified_clean = read_tsv(OUT / "03_verified_clean.tsv")
    repaired = read_tsv(OUT / "04_repaired.tsv")
    doi_by_uni = load_doi_records()

    gold: list[dict] = []
    provenance: list[dict] = []

    def emit(row: dict, fate: str, reason: str, *, keep: bool) -> None:
        rec = {k: row.get(k, "") for k in PASS_HEADER}
        rec["fate"] = fate
        rec["fate_reason"] = reason
        provenance.append(rec)
        if keep:
            gold.append({k: row.get(k, "") for k in PASS_HEADER})

    # 1) verified-clean rows pass straight through
    for r in verified_clean:
        emit(r, "VERIFIED_CLEAN", "validated+VERIFIED, no repair needed", keep=True)

    # 2) repaired (PARTIAL) rows: apply UniProt fix, DOI repair, and adjudications
    missing_doi: list[str] = []
    for r in repaired:
        uni = (r.get("uniprot") or "").strip()  # carried-forward accession from 02_*
        rec = doi_by_uni.get(uni)

        # placeability: Phase 1 needs coordinates from locus_tag in refseq_genome
        placeable = has_val(r.get("locus_tag", "")) or has_val(uni)

        if uni in DROP_CONTRADICTED:
            emit(r, "DROP_CONTRADICTED", DROP_CONTRADICTED[uni], keep=False)
            continue
        if not placeable:
            emit(r, "DROP_UNPLACEABLE", "no usable locus_tag or UniProt -> cannot position", keep=False)
            continue
        if r.get("uniprot_status") == "UNRESOLVED":
            # audit said the accession was wrong and we could not rebuild it; if the row is
            # still placeable by locus we keep it but flag the soft-unverified accession
            note = "audit-flagged UniProt could not be rebuilt; kept on locus_tag"
        else:
            note = ""

        # DOI: prefer a repaired record; else the stored DOI if the audit left it OK
        if rec:
            if rec.get("status") == "RESOLVED" and rec.get("corrected_doi"):
                r["primary_ref"] = rec["corrected_doi"]
            elif rec.get("status") in ("CONTRADICTED", "NOT_FOUND"):
                # contradicted rows not in DROP_CONTRADICTED were reclassified instead;
                # NOT_FOUND should not occur (none returned) -> guard by dropping
                if uni not in RECLASSIFY_T1SS:
                    emit(r, "DROP_NO_CITATION", f"DOI status {rec.get('status')}", keep=False)
                    continue
        elif r.get("doi_status") == "BROKEN_RESIDUAL":
            missing_doi.append(f"{r.get('gene')}/{uni}")
            emit(r, "DROP_NO_CITATION", "broken DOI with no repair record found", keep=False)
            continue

        if uni in RECLASSIFY_T1SS:
            r["ss_type"] = "T1SS"
            if rec and rec.get("corrected_doi"):
                r["primary_ref"] = rec["corrected_doi"]
            emit(r, "RECLASSIFIED_T1SS", RECLASSIFY_T1SS[uni], keep=True)
            continue

        fate_reason = "PARTIAL repaired (UniProt + DOI verified)"
        if note:
            fate_reason += f"; {note}"
        emit(r, "REPAIRED_VERIFIED", fate_reason, keep=True)

    # 2.6) independent DOI recheck via Crossref over the final gold set
    cache: dict[str, bool] = {}
    dois = sorted({(g.get("primary_ref") or "").strip() for g in gold if has_val(g.get("primary_ref", ""))})
    doi_fail: list[str] = []
    print(f"independent DOI.org recheck of {len(dois)} distinct DOIs ...", file=sys.stderr)
    for d in dois:
        if not doi_resolves(d, cache):
            doi_fail.append(d)
        time.sleep(0.05)

    # 2.7) emit
    def write(path: Path, header: list[str], rows: list[dict]) -> None:
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in header})

    write(OUT / "gold_set_corpus.tsv", PASS_HEADER, gold)
    write(OUT / "gold_set_provenance.tsv", PASS_HEADER + ["fate", "fate_reason"], provenance)

    # report
    from collections import Counter

    fates = Counter(p["fate"] for p in provenance)
    by_type = Counter(g["ss_type"] for g in gold)
    print(
        f"\ninputs: {len(verified_clean)} verified-clean + {len(repaired)} repaired = {len(verified_clean) + len(repaired)}"
    )
    print("fates:")
    for k, v in sorted(fates.items()):
        print(f"  {k:20} {v}")
    print(f"\nGOLD SET: {len(gold)} effectors")
    print(f"  by ss_type: {dict(sorted(by_type.items()))}")
    print(f"\nDOI recheck: {len(dois)} distinct, {len(doi_fail)} failed to resolve")
    for d in doi_fail:
        print(f"  FAILED: {d}")
    if missing_doi:
        print(f"\nWARN: {len(missing_doi)} broken-DOI rows had no repair record: {missing_doi}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
