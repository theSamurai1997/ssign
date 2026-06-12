#!/usr/bin/env python3
"""Phase 0a, steps 2.4-2.5 (scriptable part): repair the PARTIAL gold-set rows.

The repair queue (03_repair_queue.tsv) holds 183 real validated effectors whose
*metadata* the 2026-06-08 audit flagged: a wrong UniProt accession (the stored ID
maps to a different gene), and/or a broken/wrong primary-reference DOI. The
locus_tag is trusted, so we rebuild from it.

What this script does (deterministic + one network source, UniProt REST):
  2.4  Rebuild the UniProt accession from the trusted locus_tag. For each row we
       query UniProt for the protein whose ordered-locus-name (OLN) equals the
       locus_tag, and take that accession. OLN exact-match is the disambiguator,
       so a multi-organism hit list collapses to the one real protein.
         uniprot_status: CONFIRMED (rebuilt == stored) | CORRECTED (differs)
                         | UNRESOLVED (no OLN match) | NO_LOCUS (locus_tag empty)
  2.5  Lift a corrected DOI when the audit named one verbatim ("Correct DOI: ...").
       Ambiguous "bad DOI, candidates are X or Y" notes are NOT auto-picked; they
       and any UNRESOLVED UniProt go to the residual file for the agent pass (2.5)
       and the independent verification pass (2.6).

UniProt is used here strictly to look up an ID/coordinate from a known locus, never
to discover effectors. Responses are cached so re-runs are offline and free.

Outputs (to data/gold_build/):
  04_repaired.tsv          every repair-queue row + repair columns
  04_repair_residual.tsv   rows still needing the agent/verification pass
  .uniprot_cache.json      query cache (locus -> hit list)
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

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "gold_build"
IN = OUT / "03_repair_queue.tsv"
CACHE = OUT / ".uniprot_cache.json"

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
UA = "ssign-benchmark/0.1 (teoreid@gmail.com)"

# audit phrasing that means "the stored DOI is wrong/unreachable"
BAD_DOI = re.compile(
    r"\b(bad doi|wrong doi|404|not findable|does not resolve|doesn'?t resolve"
    r"|invalid doi|wrong reference|placeholder)\b",
    re.IGNORECASE,
)
# audit phrasing that hands us the replacement verbatim
CORRECT_DOI = re.compile(r"correct doi[:\s]*((?:10\.\d{4,9}/)\S+)", re.IGNORECASE)


def norm_locus(s: str) -> str:
    return (s or "").upper().replace("_", "").strip()


def load_cache() -> dict:
    if CACHE.exists():
        return json.loads(CACHE.read_text())
    return {}


def save_cache(cache: dict) -> None:
    CACHE.write_text(json.dumps(cache, indent=0))


def uniprot_query(term: str, cache: dict) -> list[dict]:
    """Return [{acc, oln:[...], organism}] for `gene:<term>`. Cached, with retry."""
    if term in cache:
        return cache[term]
    params = {
        "query": f"gene:{term}",
        "fields": "accession,gene_oln,organism_name,reviewed",
        "format": "json",
        "size": "10",
    }
    url = UNIPROT + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                data = json.load(r)
            hits = []
            for res in data.get("results", []):
                oln = [x.get("value", "") for g in res.get("genes", []) for x in g.get("orderedLocusNames", [])]
                hits.append(
                    {
                        "acc": res["primaryAccession"],
                        "oln": oln,
                        "organism": res.get("organism", {}).get("scientificName", ""),
                        "reviewed": res.get("entryType", "").startswith("UniProtKB reviewed"),
                    }
                )
            cache[term] = hits
            return hits
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            last = e
            time.sleep(2 * (attempt + 1))
    print(f"  WARN uniprot query failed for {term}: {last}", file=sys.stderr)
    cache[term] = []
    return []


def has_locus(locus: str) -> bool:
    """True if the locus_tag field holds a real value (not blank or a '-' placeholder)."""
    return bool(locus.strip()) and locus.strip() != "-"


def rebuild_uniprot(locus: str, cache: dict) -> str | None:
    """Accession whose OLN == locus (underscore/case-insensitive), or None.

    Among OLN-exact matches, prefer a reviewed Swiss-Prot entry over TrEMBL so a
    correction lands on the curated record when one exists.
    """
    if not has_locus(locus):
        return None
    target = norm_locus(locus)
    matches: list[dict] = []
    for term in (locus.strip(), locus.replace("_", "").strip()):
        for hit in uniprot_query(term, cache):
            if any(norm_locus(o) == target for o in hit["oln"]):
                matches.append(hit)
        if matches:
            break
        time.sleep(0.2)  # be polite between distinct queries
    if not matches:
        return None
    matches.sort(key=lambda h: not h["reviewed"])  # reviewed first
    return matches[0]["acc"]


def main() -> int:
    rows = list(csv.DictReader(IN.open(), delimiter="\t"))
    header = list(rows[0].keys()) if rows else []
    cache = load_cache()

    extra = ["uniprot_old", "uniprot_rebuilt", "uniprot_status", "doi_old", "doi_new", "doi_status"]
    repaired: list[dict] = []
    residual: list[dict] = []

    for i, r in enumerate(rows, 1):
        locus = (r.get("locus_tag") or "").strip()
        old_uni = (r.get("uniprot") or "").strip()
        notes = r.get("verification_notes") or ""
        uni_flagged = bool(re.search(r"uniprot", notes, re.IGNORECASE))

        # 2.4 UniProt: the audit is the authority on what is wrong. Only *replace*
        # an accession the audit flagged; for the rest, cross-check and keep.
        rebuilt = rebuild_uniprot(locus, cache)
        r["uniprot_old"] = old_uni
        r["uniprot_rebuilt"] = rebuilt or ""
        stored_empty = not has_locus(old_uni)  # blank or '-' placeholder accession
        if uni_flagged:
            # audit said the stored accession is wrong -> must repair from locus
            r["uniprot_status"] = "CORRECTED" if rebuilt else "UNRESOLVED"
        elif stored_empty:
            # no accession was stored -> fill from locus (nothing to conflict with)
            r["uniprot_status"] = "FILLED" if rebuilt else "MISSING"
        else:
            # audit did not question the accession -> keep it, OLN match is a check
            if rebuilt and rebuilt == old_uni:
                r["uniprot_status"] = "CONFIRMED"
            elif rebuilt and rebuilt != old_uni:
                r["uniprot_status"] = "CONFLICT"  # stored kept; rebuilt differs, verify
            elif not has_locus(locus):
                r["uniprot_status"] = "KEPT_NO_LOCUS"  # gene-symbol/blank locus, kept
            else:
                r["uniprot_status"] = "KEPT_UNVERIFIED"  # locus not OLN-indexed, kept

        # carry the rebuilt accession forward when we repaired or filled; else keep stored
        r["uniprot"] = rebuilt if (rebuilt and (uni_flagged or stored_empty)) else old_uni

        # 2.5 DOI: lift only an explicitly-named correction; never guess
        r["doi_old"] = (r.get("primary_ref") or "").strip()
        m = CORRECT_DOI.search(notes)
        if m:
            r["doi_new"] = m.group(1).rstrip(".;,)")
            r["doi_status"] = "CORRECTED"
        elif BAD_DOI.search(notes):
            r["doi_new"] = ""
            r["doi_status"] = "BROKEN_RESIDUAL"
        else:
            r["doi_new"] = r["doi_old"]
            r["doi_status"] = "OK"

        repaired.append(r)
        # residual = real blocking work for the agent pass (2.5) / verify pass (2.6):
        # an audit-flagged accession we could not rebuild; a broken DOI with no named
        # correction; or an unplaceable row (neither a usable accession nor locus, so
        # Phase 1 cannot position it). KEPT_*/CONFLICT are non-blocking (verify eyeballs).
        unplaceable = not has_locus(r["uniprot"]) and not has_locus(locus)
        if r["uniprot_status"] == "UNRESOLVED" or r["doi_status"] == "BROKEN_RESIDUAL" or unplaceable:
            residual.append(r)

        if i % 25 == 0:
            save_cache(cache)
            print(f"  ...{i}/{len(rows)}")

    save_cache(cache)

    full_header = header + extra
    for path, data in (("04_repaired.tsv", repaired), ("04_repair_residual.tsv", residual)):
        with (OUT / path).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=full_header, delimiter="\t")
            w.writeheader()
            for r in data:
                w.writerow({k: r.get(k, "") for k in full_header})

    def tally(key):
        c = {}
        for r in repaired:
            c[r[key]] = c.get(r[key], 0) + 1
        return dict(sorted(c.items()))

    print(f"\nrepair queue rows:     {len(repaired)}")
    print(f"  uniprot_status:      {tally('uniprot_status')}")
    print(f"  doi_status:          {tally('doi_status')}")
    print(f"residual (-> agent 2.5 + verify 2.6): {len(residual)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
