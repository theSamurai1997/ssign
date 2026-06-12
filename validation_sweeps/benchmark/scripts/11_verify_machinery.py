#!/usr/bin/env python3
"""
11_verify_machinery.py  (Phase 0b, task 3.6)

Independent verification pass over the resolved machinery genes from 3.5.
Three checks, then a per-gene `verified` flag and per-instance `anchor_status`:

  1. DOI resolves   -- HEAD https://doi.org/<doi> WITHOUT following redirects:
                       a 3xx (redirect to the publisher) means the DOI is registered;
                       404 means it is dead. Cached to .doi_cache.json.
  2. PMID valid     -- batch NCBI esummary over every distinct PMID; a PMID that
                       returns a PubMed record is a real, indexed citation.
  3. quote names    -- the gene's native symbol or unified family token appears in
                       the verbatim quote curated in 3.4 (else it is an operon-level
                       "collective" quote, recorded as such, not a failure).

verified (per gene) := match_tier == "alias"  (locus_tag is the genome's own name)
                       AND coordinates present
                       AND (doi_resolves OR pmid_valid)
Product-tier and ambiguous matches are never auto-verified -> verdict "review".
Unresolved genes -> "unverified".

anchor_status (per instance):
  ANCHORED      >=1 verified machinery gene (operon is locatable for Phase 1)
  REVIEW        genes resolved but none alias-verified (only product/ambiguous)
  NEEDS_ANCHOR  0 machinery genes resolved

Inputs:
  data/machinery/machinery_resolved.tsv   (from 3.5)
  data/machinery_raw/*.json               (quotes, pmids from 3.4)
Outputs:
  data/machinery/machinery_answer_key.tsv
  data/machinery/.doi_cache.json          (resolution cache)

Run:
  .venv/bin/python scripts/11_verify_machinery.py
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
RESOLVED = BENCH / "data" / "machinery" / "machinery_resolved.tsv"
RAW_DIR = BENCH / "data" / "machinery_raw"
OUT_TSV = BENCH / "data" / "machinery" / "machinery_answer_key.tsv"
DOI_CACHE = BENCH / "data" / "machinery" / ".doi_cache.json"

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
UA = "ssign-benchmark/0.1 (teoreid@gmail.com)"
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")


def doi_resolves(doi: str) -> bool:
    """True iff doi.org has the DOI registered (3xx redirect, not 404)."""
    url = "https://doi.org/" + urllib.parse.quote(doi, safe="/:")
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})

    class _NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, *a, **k):
            return None

    opener = urllib.request.build_opener(_NoRedirect)
    try:
        opener.open(req, timeout=20)
        return True  # 2xx (rare for doi.org) also counts as registered
    except urllib.error.HTTPError as e:
        return e.code in (301, 302, 303, 307, 308)
    except Exception:  # noqa: BLE001 - network hiccup -> treat as unknown=False
        return False


def resolve_dois(dois):
    cache = json.loads(DOI_CACHE.read_text()) if DOI_CACHE.exists() else {}
    todo = [d for d in dois if d not in cache]
    for i, d in enumerate(sorted(todo), 1):
        cache[d] = doi_resolves(d)
        if i % 20 == 0:
            print(f"  doi.org {i}/{len(todo)}...", file=sys.stderr)
            DOI_CACHE.write_text(json.dumps(cache, indent=0))
        time.sleep(0.15)
    DOI_CACHE.write_text(json.dumps(cache, indent=0))
    return cache


def valid_pmids(pmids):
    """Batch NCBI esummary; return the set of PMIDs that return a real record."""
    valid = set()
    pmids = sorted(p for p in pmids if p.isdigit())
    for i in range(0, len(pmids), 180):
        chunk = pmids[i : i + 180]
        params = urllib.parse.urlencode(
            {
                "db": "pubmed",
                "id": ",".join(chunk),
                "retmode": "json",
                "tool": "ssign-benchmark",
                "email": "teoreid@gmail.com",
            }
        )
        # Fail loud on a persistent error: a swallowed failure would mark every
        # PMID in the chunk invalid (up to 180 genes), silently corrupting the
        # verdict. Retry once, then raise.
        for attempt in range(2):
            try:
                req = urllib.request.Request(f"{ESUMMARY}?{params}", headers={"User-Agent": UA})
                data = json.loads(urllib.request.urlopen(req, timeout=60).read())
                res = data.get("result", {})
                for uid in res.get("uids", []):
                    rec = res.get(uid, {})
                    if rec.get("title") and "error" not in rec:
                        valid.add(uid)
                break
            except Exception as e:  # noqa: BLE001 - network, retry then raise
                if attempt == 1:
                    raise RuntimeError(f"esummary failed for PMID chunk {i}-{i + len(chunk)}: {e}")
                print(f"  esummary retry ({e})", file=sys.stderr)
                time.sleep(3)
        time.sleep(0.4)
    return valid


def quote_match(gene, family, quote):
    """gene_named | collective | absent."""
    if not quote:
        return "absent"
    tokens = [gene]
    m = re.match(r"^([A-Za-z][A-Za-z0-9]{1,7})\b", family or "")
    if m:
        tokens.append(m.group(1))
    for t in tokens:
        if re.search(rf"(?<![A-Za-z0-9]){re.escape(t)}(?![A-Za-z0-9])", quote, re.IGNORECASE):
            return "gene_named"
    return "collective"


def main():
    # quotes / families from 3.4, keyed (instance_id, gene)
    quotes = {}
    for jf in sorted(RAW_DIR.glob("*.json")):
        with open(jf) as fh:
            d = json.load(fh)
        for g in d["machinery_genes"]:
            quotes[(d["instance_id"], g["gene"].strip())] = (g.get("quote", ""), g.get("family", ""))

    with open(RESOLVED) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    doi_cache = resolve_dois({r["doi"] for r in rows if r["doi"]})
    pmid_ok = valid_pmids({r["pmid"] for r in rows if r["pmid"]})

    per_inst = {}  # instance_id -> {verified, resolved, status, ss_type}
    for r in rows:
        doi = r["doi"]
        r["doi_wellformed"] = str(bool(DOI_RE.match(doi))).lower()
        r["doi_resolves"] = str(bool(doi_cache.get(doi))).lower()
        r["pmid_valid"] = str(r["pmid"] in pmid_ok).lower()
        q, fam = quotes.get((r["instance_id"], r["gene"]), (r.get("quote", ""), r["family"]))
        r["quote_match"] = quote_match(r["gene"], fam, q)

        cited = r["doi_resolves"] == "true" or r["pmid_valid"] == "true"
        has_coords = bool(r["start"])
        if r["match_tier"] == "alias" and has_coords and cited:
            r["verdict"] = "verified"
        elif r["match_tier"] in ("alias", "product") and has_coords:
            r["verdict"] = "review"
        else:
            r["verdict"] = "unverified"

        s = per_inst.setdefault(
            r["instance_id"], {"verified": 0, "resolved": 0, "status": r["status"], "ss_type": r["ss_type"]}
        )
        if r["verdict"] == "verified":
            s["verified"] += 1
        if r["match_tier"] != "none":
            s["resolved"] += 1

    for r in rows:
        s = per_inst[r["instance_id"]]
        r["anchor_status"] = "ANCHORED" if s["verified"] else "REVIEW" if s["resolved"] else "NEEDS_ANCHOR"

    out_fields = list(rows[0].keys())
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # ---- report ----
    n = len(rows)
    nv = sum(1 for r in rows if r["verdict"] == "verified")
    nr = sum(1 for r in rows if r["verdict"] == "review")
    nu = sum(1 for r in rows if r["verdict"] == "unverified")
    print(f"\nwrote {OUT_TSV.relative_to(BENCH)}")
    print(f"machinery genes        : {n}")
    print(f"  verified             : {nv} ({nv / n:.0%})")
    print(f"  review (product/no-cite): {nr} ({nr / n:.0%})")
    print(f"  unverified (no locus) : {nu} ({nu / n:.0%})")
    print(f"  DOI resolves          : {sum(1 for r in rows if r['doi_resolves'] == 'true')}/{n}")
    print(f"  PMID valid            : {sum(1 for r in rows if r['pmid_valid'] == 'true')}/{n}")
    print(
        f"  quote names gene      : {sum(1 for r in rows if r['quote_match'] == 'gene_named')}"
        f" | collective {sum(1 for r in rows if r['quote_match'] == 'collective')}"
    )

    anchored = [i for i, s in per_inst.items() if s["verified"]]
    review = [i for i, s in per_inst.items() if not s["verified"] and s["resolved"]]
    need = [i for i, s in per_inst.items() if not s["resolved"]]
    print(f"\ninstances ({len(per_inst)}): ANCHORED {len(anchored)} | REVIEW {len(review)} | NEEDS_ANCHOR {len(need)}")
    for label, ids in (("REVIEW", review), ("NEEDS_ANCHOR", need)):
        for i in sorted(ids):
            s = per_inst[i]
            print(f"    {label:12s} {i:10s} {s['ss_type']:5s} {s['status']}")

    # genes whose stored DOI does NOT resolve (curation-integrity flags for hand review)
    bad_doi = sorted(
        {(r["doi"], r["pmid"]) for r in rows if r["doi_resolves"] == "false" and r["pmid_valid"] == "false"}
    )
    if bad_doi:
        print(f"\nDOIs that neither resolve nor have a valid PMID ({len(bad_doi)}):")
        for doi, pmid in bad_doi:
            print(f"    doi={doi}  pmid={pmid}")


if __name__ == "__main__":
    sys.exit(main())
