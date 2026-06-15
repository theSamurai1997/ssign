#!/usr/bin/env python3
"""Full-table deterministic citation-consistency sweep (all 925 positives, not just ssign-found).

Script 41 ran the deterministic CrossRef check over the 51 ssign-FOUND gold effectors only. This
script runs the identical check over EVERY row of positives_all.tsv, so the whole training-label
table gets the same provenance verdict, not just the recall-headline subset.

Simpler than 41: positives_all carries the sourcing DOI (primary_ref, else instance_source_doi) and
organism on the row itself, so there is no found->ceiling->positives join bridge. We reuse 41's
CrossRef fetch + JATS/entity-stripping + word-boundary token matching verbatim (single source of
truth) and apply the same verdict ladder per row.

The fetch cost is per DISTINCT DOI (330 of them, most already cached from the found run), not per row.

Verdict per row (same ladder as 41):
  CONSISTENT       - row's gene or genus appears in the resolved paper's title/abstract/journal.
  FLAG_WRONG_TOPIC - abstract present, no gene/genus AND no SS-type phrase -> wrong-field DOI.
  FLAG_GENE_ABSENT - abstract present, SS-type topic present, gene/genus not named -> imprecise.
  DOI_UNRESOLVED   - DOI.org 404 / no DOI on the row.
  FETCH_ERROR      - DOI resolves but CrossRef unreachable after retries -> re-run to settle.
  INDETERMINATE    - resolves but no CrossRef abstract -> title-only can't refute; queue pass-2.

Inputs : data/dataset/positives_all.tsv
Outputs: data/dataset/citation_consistency_full.tsv   (one row per positive)
         data/dataset/.doi_cache.json / .crossref_cache.json   (shared with 41; re-runs offline)
Run:     .venv/bin/python scripts/44_citation_consistency_full.py
"""

from __future__ import annotations

import importlib.util
import json
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "data" / "dataset"

# Reuse 41's CrossRef machinery verbatim (module name can't start with a digit -> load by spec).
_spec = importlib.util.spec_from_file_location("cc41", Path(__file__).parent / "41_citation_consistency.py")
_cc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_cc)
from bench_io import by_type, norm_doi, read_tsv, write_tsv  # noqa: E402  (scripts/ on sys.path)


def verdict_for(gene: str, genus: str, ss: str, doi: str, doi_cache: dict, cr_cache: dict) -> tuple[str, dict]:
    """Run the deterministic ladder for one row. Returns (verdict, evidence_fields)."""
    needs_net = bool(doi) and (doi not in doi_cache or doi not in cr_cache)
    resolves = _cc._fz.doi_resolves(doi, doi_cache) if doi else False
    content = _cc.crossref_content(doi, cr_cache) if doi else None
    if needs_net:
        time.sleep(0.05)

    in_crossref = bool(content and content["in_crossref"])
    has_abstract = bool(content and content["abstract"])
    haystack = _cc.crossref_haystack(content)
    gene_hit = bool(gene) and _cc._token_hit(gene, haystack)
    genus_hit = bool(genus) and _cc._token_hit(genus, haystack)
    ss_hit = any(ph in haystack for ph in _cc.SS_PHRASES.get(ss, []))

    if not resolves:
        verdict = "DOI_UNRESOLVED"
    elif doi and content is None:
        verdict = "FETCH_ERROR"
    elif gene_hit or genus_hit:
        verdict = "CONSISTENT"
    elif not in_crossref or not has_abstract:
        verdict = "INDETERMINATE"
    elif ss_hit:
        verdict = "FLAG_GENE_ABSENT"
    else:
        verdict = "FLAG_WRONG_TOPIC"

    return verdict, {
        "doi_resolves": "yes" if resolves else "no",
        "in_crossref": "yes" if in_crossref else "no",
        "has_abstract": "yes" if has_abstract else "no",
        "gene_in_paper": "yes" if gene_hit else "no",
        "genus_in_paper": "yes" if genus_hit else "no",
        "ss_phrase_in_paper": "yes" if ss_hit else "no",
        "crossref_title": content["title"] if content else "",
    }


def main() -> int:
    pos = read_tsv(DATASET / "positives_all.tsv")
    doi_cache: dict = json.loads(_cc.DOI_CACHE.read_text()) if _cc.DOI_CACHE.exists() else {}
    cr_cache: dict = json.loads(_cc.CROSSREF_CACHE.read_text()) if _cc.CROSSREF_CACHE.exists() else {}

    rows: list[dict] = []
    for p in pos:
        gene = (p.get("gene") or "").strip()
        organism = (p.get("organism") or "").strip()
        genus = organism.split()[0] if organism else ""
        ss = (p.get("ss_type") or "").strip()
        doi = norm_doi(p.get("primary_ref", "")) or norm_doi(p.get("instance_source_doi", ""))
        verdict, ev = verdict_for(gene, genus, ss, doi, doi_cache, cr_cache)
        rows.append(
            {
                "gene": gene,
                "uniprot": (p.get("uniprot") or "").strip(),
                "locus_tag": (p.get("locus_tag") or "").strip(),
                "sys_instance_id": (p.get("sys_instance_id") or "").strip(),
                "ss_type": ss,
                "organism": organism,
                "evidence_tier": (p.get("evidence_tier") or "").strip(),
                "sourcing_doi": doi,
                **ev,
                "verdict": verdict,
            }
        )

    _cc.DOI_CACHE.write_text(json.dumps(doi_cache, indent=0))
    _cc.CROSSREF_CACHE.write_text(json.dumps(cr_cache, indent=0))
    cols = list(rows[0].keys())
    write_tsv(DATASET / "citation_consistency_full.tsv", cols, rows)

    verdicts = Counter(r["verdict"] for r in rows)
    print(f"wrote data/dataset/citation_consistency_full.tsv  ({len(rows)} positives)")
    print("verdicts:", dict(verdicts))
    print("by ss_type:", by_type(rows))
    print("\nby evidence_tier:")
    for tier in sorted({r["evidence_tier"] for r in rows}):
        sub = [r for r in rows if r["evidence_tier"] == tier]
        vc = Counter(r["verdict"] for r in sub)
        print(f"  {tier or '(blank)':12s} n={len(sub):3d}  {dict(vc)}")

    flag_types = ("FLAG_WRONG_TOPIC", "FLAG_GENE_ABSENT", "DOI_UNRESOLVED", "FETCH_ERROR", "INDETERMINATE")
    flagged = [r for r in rows if r["verdict"] in flag_types]
    print(f"\nqueued for pass-2 (flagged + indeterminate): {len(flagged)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
