#!/usr/bin/env python3
"""Phase 0b, step 3.2: fetch the RefSeq/INSDC GenBank flatfile for every instance genome.

instances.tsv lists one replicon accession per instance (chromosome or plasmid where the
system sits), 62 distinct. We pull each as a GenBank flatfile with full feature table, so
later steps can resolve paper-named machinery genes to a locus_tag + coordinates (step 3.5)
and build the gene-order index (Phase 1, task 5.1). RefSeq is used strictly for coordinate
and ID lookup here, never to discover effectors or machinery.

rettype=gbwithparts forces CON (contig) records to expand to full features+sequence, so
finished replicons that NCBI stores as assembly pointers still come back complete.

One file per accession under data/refseq_cache/<acc>.gb; existing files are skipped, so a
re-run only fetches what is missing (network failures are resumable).

Input : data/machinery/instances.tsv  (refseq_genome column)
Output: data/refseq_cache/<acc>.gb
"""

from __future__ import annotations

import csv
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INSTANCES = ROOT / "data" / "machinery" / "instances.tsv"
CACHE = ROOT / "data" / "refseq_cache"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA = "ssign-benchmark/0.1 (teoreid@gmail.com)"


def fetch(acc: str) -> bool:
    out = CACHE / f"{acc}.gb"
    if out.exists() and out.stat().st_size > 0:
        return True
    params = urllib.parse.urlencode(
        {
            "db": "nuccore",
            "id": acc,
            "rettype": "gbwithparts",
            "retmode": "text",
            "tool": "ssign-benchmark",
            "email": "teoreid@gmail.com",
        }
    )
    for attempt in range(4):
        try:
            req = urllib.request.Request(f"{EFETCH}?{params}", headers={"User-Agent": UA})
            data = urllib.request.urlopen(req, timeout=300).read()
            if not data.startswith(b"LOCUS"):
                raise ValueError(f"unexpected response head: {data[:60]!r}")
            out.write_bytes(data)
            time.sleep(0.4)  # NCBI: <=3 req/s without an API key
            return True
        except Exception as e:  # noqa: BLE001 - network, retry
            print(f"  retry {acc} ({attempt + 1}/4): {e}", file=sys.stderr)
            time.sleep(3 * (attempt + 1))
    return False


def main() -> int:
    CACHE.mkdir(parents=True, exist_ok=True)
    accs = sorted(
        {
            r["refseq_genome"].strip()
            for r in csv.DictReader(INSTANCES.open(), delimiter="\t")
            if r["refseq_genome"].strip()
        }
    )
    print(f"fetching {len(accs)} replicon GenBank records -> {CACHE}", file=sys.stderr)
    ok = 0
    for i, acc in enumerate(accs, 1):
        if fetch(acc):
            ok += 1
        else:
            print(f"  FAILED: {acc}", file=sys.stderr)
        if i % 10 == 0:
            print(f"  ...{i}/{len(accs)} ({ok} ok)", file=sys.stderr)
    print(f"done: {ok}/{len(accs)} cached")
    return 0 if ok == len(accs) else 1


if __name__ == "__main__":
    raise SystemExit(main())
