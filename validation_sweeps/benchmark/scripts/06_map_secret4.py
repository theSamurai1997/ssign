#!/usr/bin/env python3
"""Phase 0a, step 2b.3 helper: resolve every SecReT4 verified effector to a locus_tag.

SecReT4's FASTA header carries a RefSeq protein accession (YP_/NP_/WP_) but no
locus_tag, which is the key we dedup on against the corpus gold set. NCBI efetch
(GenPept) returns the CDS /locus_tag and source /organism for each accession, so we
batch-fetch and attach those. RefSeq is used here strictly as a coordinate/ID
lookup, never to discover effectors.

We also tag each entry with its T4SS subtype and a detection gate:
  gate = GATED  if the organism's T4SS is Dot/Icm (type-IVB, MPF_I): Legionella,
                Coxiella, Rickettsiella. These effectors are genome-dispersed and
                ssign's detection of type-IVB is unconfirmed -> held pending Phase 2.
  gate = ACTIVE otherwise (type-IVA VirB/D4, which ssign detects).

Input : data/external_dbs/parsed_external.tsv (SecReT4 rows with protein_acc)
Output: data/external_dbs/secret4/secret4_mapped.json
        .ncbi_cache.json (acc -> {locus, organism}; re-runs are offline)
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data" / "external_dbs"
CACHE = EXT / "secret4" / ".ncbi_cache.json"
OUT = EXT / "secret4" / "secret4_mapped.json"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA = "ssign-benchmark/0.1 (teoreid@gmail.com)"

GATED_GENERA = {"Legionella", "Coxiella", "Rickettsiella"}


def genus(org: str) -> str:
    return org.split()[0] if org.split() else ""


def load_cache() -> dict:
    return json.loads(CACHE.read_text()) if CACHE.exists() else {}


def efetch_batch(accs: list[str], cache: dict) -> None:
    """Fill cache[acc] = {'locus':..., 'organism':...} for any uncached accs."""
    todo = [a for a in accs if a not in cache]
    for i in range(0, len(todo), 120):
        chunk = todo[i : i + 120]
        data = urllib.parse.urlencode(
            {
                "db": "protein",
                "id": ",".join(chunk),
                "rettype": "gp",
                "retmode": "text",
                "tool": "ssign-benchmark",  # NCBI E-utilities policy: identify tool + email
                "email": "teoreid@gmail.com",
            }
        ).encode()
        txt = ""
        for attempt in range(4):
            try:
                req = urllib.request.Request(EFETCH, data=data, headers={"User-Agent": UA})
                txt = urllib.request.urlopen(req, timeout=120).read().decode("utf-8", "replace")
                break
            except Exception as e:  # noqa: BLE001 - network, retry
                print(f"  retry batch {i}: {e}", file=sys.stderr)
                time.sleep(3 * (attempt + 1))
        for rec in txt.split("\nLOCUS")[1:]:
            vm = re.search(r"VERSION\s+(\S+)", rec)
            lm = re.search(r'/locus_tag="([^"]+)"', rec)
            om = re.search(r'/organism="([^"]+)"', rec)
            if vm:
                cache[vm.group(1)] = {"locus": lm.group(1) if lm else "", "organism": om.group(1) if om else ""}
        CACHE.write_text(json.dumps(cache))
        time.sleep(0.4)


def main() -> int:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        r
        for r in csv.DictReader((EXT / "parsed_external.tsv").open(), delimiter="\t")
        if r["source_db"] == "SecReT4" and r["protein_acc"]
    ]
    accs = sorted({r["protein_acc"] for r in rows})
    cache = load_cache()
    print(f"mapping {len(accs)} SecReT4 accessions ({len(rows)} rows) ...", file=sys.stderr)
    efetch_batch(accs, cache)

    out = []
    for r in rows:
        info = cache.get(r["protein_acc"], {})
        g = genus(r["organism"])
        out.append(
            {
                "source_db": "SecReT4",
                "source_id": r["protein_acc"],
                "ss_type": "T4SS",
                "organism": r["organism"],
                "ncbi_organism": info.get("organism", ""),
                "locus_tag": info.get("locus", ""),
                "uniprot": "",  # SecReT4 FASTA has no UniProt; dedup on locus
                "refseq_genome": "",
                "doi": "",  # attached later from SecReT4 record pages if needed
                "description": r["description"],
                "t4ss_subtype": "type-IVB (Dot/Icm)" if g in GATED_GENERA else "type-IVA (VirB/D4)",
                "gate": "GATED" if g in GATED_GENERA else "ACTIVE",
            }
        )

    OUT.write_text(json.dumps(out, indent=1))
    mapped = sum(1 for r in out if r["locus_tag"])
    from collections import Counter

    print(f"\n{len(out)} rows -> {OUT}; with locus: {mapped}")
    print("gate:", dict(Counter(r["gate"] for r in out)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
