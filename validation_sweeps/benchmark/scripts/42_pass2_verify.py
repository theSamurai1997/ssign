#!/usr/bin/env python3
"""Dataset group 6, task 6.2: merge the pass-2 agent re-audit + deterministically verify its DOIs.

Pass-2 sent the 30 non-CONSISTENT found effectors to literature agents (5 batches) + a main-session
manual pass (batch C, AUP-blocked on B. pseudomallei). Each returned a status + a corrected DOI +
a verbatim secretion quote. This script does NOT trust those returns blindly: for every DOI an agent
hands back, it independently confirms the DOI is real (registered, on CrossRef with a title) and that
the effector's gene or genus actually appears in the CrossRef metadata — the same token check pass-1
used. A fabricated or wrong-paper DOI therefore cannot pass silently.

Per merged row -> verify_status:
  VERIFIED          - DOI is real AND gene/genus/SS-topic appears in the CrossRef record.
  RESOLVES_NOABSTRACT - DOI is real but CrossRef has no abstract (can't token-confirm); the agent's
                        verbatim quote is the remaining evidence, accept with a flag.
  UNVERIFIED_DOI    - the returned DOI does not resolve / is not on CrossRef -> do NOT trust.
  NO_DOI            - status NOT_FOUND: the agent found no qualifying paper (row to drop/hold).

Inputs : data/dataset/pass2_raw/batch_*.json   (agent + manual returns)
         data/dataset/pass2_input.json          (organism/genus/ss per idx)
Outputs: data/dataset/pass2_results.tsv         (merged + verified, one row per idx)
Run:     .venv/bin/python scripts/42_pass2_verify.py
"""

from __future__ import annotations

import importlib.util
import json
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "data" / "dataset"
RAW = DATASET / "pass2_raw"

_spec = importlib.util.spec_from_file_location("cc", Path(__file__).parent / "41_citation_consistency.py")
_cc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_cc)
from bench_io import norm_doi, write_tsv  # noqa: E402

OUT_COLS = [
    "idx",
    "gene",
    "ss_type",
    "organism",
    "status",
    "final_doi",
    "verify_status",
    "doi_real",
    "gene_or_genus_in_paper",
    "corrected_ss_type",
    "confidence",
    "verbatim_quote",
    "quote_source_doi",
    "notes",
]


def main() -> int:
    meta = {r["idx"]: r for r in json.loads((DATASET / "pass2_input.json").read_text())}
    results = []
    for bf in sorted(RAW.glob("batch_*.json")):
        results.extend(json.loads(bf.read_text()))
    results.sort(key=lambda r: r["idx"])

    doi_cache = json.loads(_cc.DOI_CACHE.read_text()) if _cc.DOI_CACHE.exists() else {}
    cr_cache = json.loads(_cc.CROSSREF_CACHE.read_text()) if _cc.CROSSREF_CACHE.exists() else {}

    rows = []
    for r in results:
        m = meta.get(r["idx"], {})
        gene = r.get("gene") or m.get("gene", "")  # meta holds the authoritative gene
        genus = (m.get("organism") or "").split()[0] if m.get("organism") else ""
        doi = norm_doi(r.get("final_doi", ""))
        status = r.get("status", "")

        if not doi:
            verify, doi_real, hit = "NO_DOI", False, False
        else:
            needs_net = doi not in doi_cache or doi not in cr_cache
            resolves = _cc._fz.doi_resolves(doi, doi_cache)
            content = _cc.crossref_content(doi, cr_cache)
            if needs_net:
                time.sleep(0.05)
            in_cr = bool(content and content["in_crossref"])
            has_abs = bool(content and content["abstract"])
            doi_real = resolves or (in_cr and bool(content["title"]))
            hit = _cc.effector_hit(gene, genus, r.get("ss_type", ""), _cc.crossref_haystack(content))
            if not doi_real:
                verify = "UNVERIFIED_DOI"
            elif hit:
                verify = "VERIFIED"
            elif not has_abs:
                verify = "RESOLVES_NOABSTRACT"
            else:
                verify = "UNVERIFIED_DOI"  # real DOI, abstract present, but effector absent -> suspect

        rows.append(
            {
                "idx": r["idx"],
                "gene": gene,
                "ss_type": r.get("ss_type", ""),
                "organism": m.get("organism", ""),
                "status": status,
                "final_doi": doi,
                "verify_status": verify,
                "doi_real": "yes" if doi_real else "no",
                "gene_or_genus_in_paper": "yes" if hit else "no",
                "corrected_ss_type": r.get("corrected_ss_type") or "",
                "confidence": r.get("confidence", ""),
                "verbatim_quote": r.get("verbatim_quote", ""),
                "quote_source_doi": r.get("quote_source_doi", ""),
                "notes": r.get("notes", ""),
            }
        )

    _cc.DOI_CACHE.write_text(json.dumps(doi_cache, indent=0))
    _cc.CROSSREF_CACHE.write_text(json.dumps(cr_cache, indent=0))
    write_tsv(DATASET / "pass2_results.tsv", OUT_COLS, rows)

    print(f"wrote data/dataset/pass2_results.tsv  ({len(rows)} rows)")
    print("status   :", dict(Counter(r["status"] for r in rows)))
    print("verify   :", dict(Counter(r["verify_status"] for r in rows)))
    mis = [r for r in rows if r["status"] == "MISASSIGNED"]
    nf = [r for r in rows if r["status"] == "NOT_FOUND"]
    bad = [r for r in rows if r["verify_status"] == "UNVERIFIED_DOI"]
    print(
        f"\nMISASSIGNED ss_type ({len(mis)}):", [(r["gene"], r["ss_type"] + "->" + r["corrected_ss_type"]) for r in mis]
    )
    print(f"NOT_FOUND / drop-or-hold ({len(nf)}):", [r["gene"] for r in nf])
    print(f"UNVERIFIED returned DOI ({len(bad)}):", [(r["gene"], r["final_doi"]) for r in bad])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
