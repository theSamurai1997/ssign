#!/usr/bin/env python3
"""Phase 0a, step 2b.2/2b.3 helper: pull full metadata for each SecReT6 effector.

The bulk FASTA only carries genome+coords. Each effector's record page
(`effector.php?id=EFF#####`) additionally carries the UniProt ID, locus_tag, gene
name, effector type, PMID and DOI, so one fetch per effector gives us both the
dedup key (UniProt/locus) and the citation in a single pass.

Input : data/external_dbs/parsed_external.tsv (SecReT6 rows -> EFF ids)
Output: data/external_dbs/secret6/secret6_records.jsonl  (one JSON object per EFF id)
        .secret6_html_cache/  (raw pages, so re-runs are offline)
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data" / "external_dbs"
CACHE = EXT / "secret6" / ".secret6_html_cache"
OUT = EXT / "secret6" / "secret6_records.jsonl"
BASE = "https://bioinfo-mml.sjtu.edu.cn/SecReT6/effector.php?id="
UA = "ssign-benchmark/0.1 (teoreid@gmail.com)"


def strip(html: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html)).strip()


# Field labels on the record page. When a cell is empty the stripped text runs the
# label straight into the next label, so any value equal to one of these is a bleed,
# not real data, and is discarded.
FIELD_LABELS = {
    "Coordinates",
    "RefSeq",
    "GenBank",
    "Uniprot",
    "KEGG",
    "PDB",
    "Pfam",
    "NCBI",
    "Strain",
    "Replicon",
    "Reference",
    "References",
    "Description",
    "Effector",
    "Cognate",
    "Related",
    "Sequence",
    "Summary",
    "Taxonomy",
}


def after(text: str, label: str, pat: str = r"(\S+)") -> str:
    """First token after `label`, or '' if absent or if the token is itself a label
    (which means the labelled cell was empty and the next label bled in)."""
    m = re.search(re.escape(label) + r"\s+" + pat, text)
    if not m:
        return ""
    v = m.group(1).strip()
    return "" if v in FIELD_LABELS else v


def parse(html: str) -> dict:
    t = strip(html)
    # "Locus tag (Gene)  Bf638R_1979 (Bfe2)  Coordinates ..."
    lt = re.search(r"Locus tag \(Gene\)\s+(\S+)(?:\s+\(([^)]+)\))?", t)
    locus = lt.group(1) if lt and lt.group(1) not in FIELD_LABELS else ""
    gene = lt.group(2) if lt and lt.group(2) else ""
    uni = after(t, "Uniprot ID")
    # genome accession is labelled "RefSeq" on most records but "GenBank" on some
    refseq = after(t, "RefSeq") or after(t, "GenBank")
    eff_type = after(t, "Effector type")
    pm = re.search(r"PMID:\s*(\d+)", t)
    pmid = pm.group(1) if pm else ""
    doi_m = re.search(r"\b(10\.\d{4,9}/[^\s\"<>]+)", t)
    doi = doi_m.group(1).rstrip(".;,)") if doi_m else ""
    rel = re.search(r"Related T6SS \(Type\)\s+(\S+)", t)
    related = rel.group(1) if rel else ""
    return {
        "locus_tag": locus,
        "gene": gene,
        "uniprot": uni,
        "refseq_genome": refseq,
        "effector_type": eff_type,
        "pmid": pmid,
        "doi": doi,
        "related_t6ss": related,
    }


def fetch(eff_id: str) -> str:
    cp = CACHE / f"{eff_id}.html"
    if cp.exists():
        return cp.read_text(encoding="utf-8", errors="replace")
    req = urllib.request.Request(BASE + eff_id, headers={"User-Agent": UA})
    for attempt in range(3):
        try:
            html = urllib.request.urlopen(req, timeout=45).read().decode("utf-8", "replace")
            cp.write_text(html, encoding="utf-8")
            time.sleep(0.4)
            return html
        except Exception as e:  # noqa: BLE001 - network, retry
            print(f"  retry {eff_id}: {e}", file=sys.stderr)
            time.sleep(3 * (attempt + 1))
    return ""


def main() -> int:
    CACHE.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((EXT / "parsed_external.tsv").open(), delimiter="\t"))
    eff_ids = [r["source_id"] for r in rows if r["source_db"] == "SecReT6" and r["source_id"].startswith("EFF")]
    print(f"scraping {len(eff_ids)} SecReT6 effector records ...", file=sys.stderr)

    out = []
    for i, eid in enumerate(eff_ids, 1):
        html = fetch(eid)
        rec = {"source_id": eid, **parse(html)} if html else {"source_id": eid, "error": "fetch_failed"}
        out.append(rec)
        if i % 50 == 0:
            print(f"  ...{i}/{len(eff_ids)}", file=sys.stderr)

    with OUT.open("w") as f:
        for r in out:
            f.write(json.dumps(r) + "\n")

    have_uni = sum(1 for r in out if r.get("uniprot"))
    have_loc = sum(1 for r in out if r.get("locus_tag"))
    have_doi = sum(1 for r in out if r.get("doi") or r.get("pmid"))
    print(f"\n{len(out)} records -> {OUT}")
    print(f"  with uniprot: {have_uni} | with locus: {have_loc} | with doi/pmid: {have_doi}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
