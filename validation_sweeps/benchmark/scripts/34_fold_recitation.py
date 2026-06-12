#!/usr/bin/env python3
"""Dataset group 1, task 1.2 (finish): validate agent re-citations, emit doi_recite.jsonl.

The 31_verify_predicted.py audit left 61 predicted rows with broken DOIs (15 distinct,
mostly RTX-toxin / Nle-effector families). A re-citation agent proposed corrected DOIs
(data/dataset/recite_results.json). This script independently verifies each proposed DOI
actually resolves (DOI.org Handle API, via script 03 `doi_resolves` -- never trust the
agent's word, same as the validated pipeline's 2.6 check) and emits the verified
broken->corrected map as doi_recite.jsonl.

Script 31 loads doi_recite.jsonl and applies the corrections, so re-running 31 -> 32 folds
the recitations in idempotently (mirrors how script 03 consumes doi_repair_*.jsonl).

Input:  data/dataset/recite_results.json   (agent output)
Output: data/dataset/doi_recite.jsonl      (verified corrections only)
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "dataset"

_spec = importlib.util.spec_from_file_location("fz", Path(__file__).parent / "03_finalize_gold_set.py")
fz = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fz)


def norm_doi(d: str) -> str:
    return (d or "").strip().removeprefix("doi:").strip()


def main() -> int:
    results = json.loads((OUT / "recite_results.json").read_text())
    doi_cache: dict = json.loads((OUT / ".doi_cache.json").read_text()) if (OUT / ".doi_cache.json").exists() else {}

    records, kept, rejected = [], 0, 0
    for r in results:
        broken = norm_doi(r.get("broken_doi", ""))
        corrected = norm_doi(r.get("corrected_doi", ""))
        if not broken or not corrected:
            rejected += 1
            print(f"  SKIP (no correction): {r.get('broken_doi')!r} genes={r.get('genes')}")
            continue
        if not fz.doi_resolves(corrected, doi_cache):
            rejected += 1
            print(f"  REJECT (corrected DOI does not resolve): {corrected} for {r.get('genes')}")
            continue
        kept += 1
        records.append(
            {
                "broken_doi": broken,
                "corrected_doi": corrected,
                "ss_type": r.get("ss_type", ""),
                "genes": r.get("genes", []),
                "quote": r.get("quote", ""),
                "source": r.get("source", ""),
                "note": r.get("note", ""),
            }
        )

    (OUT / ".doi_cache.json").write_text(json.dumps(doi_cache, indent=0))
    with (OUT / "doi_recite.jsonl").open("w") as f:
        for rec in records:
            f.write(json.dumps(rec) + "\n")

    print(f"\nverified corrections: {kept}   rejected/empty: {rejected}")
    print(f"wrote {OUT / 'doi_recite.jsonl'}")
    print("next: re-run 31_verify_predicted.py then 32_tiered_positives.py to fold them in")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
