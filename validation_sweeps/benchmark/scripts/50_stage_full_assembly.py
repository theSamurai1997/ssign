#!/usr/bin/env python3
"""Re-stage the T1SS placement genomes as their FULL assembly (all replicons), not one replicon.

Root cause of the T1SS false negatives: placed orphan effectors (hlyA, apxIA, ...) were staged as
the single replicon they sit on (a plasmid or a short WGS contig). The T1SS OMF (TolC) is a shared,
chromosomally-encoded channel and is a `loner` in the TXSScan model, so it lives on a DIFFERENT
replicon. Staging one replicon leaves TolC out, MacSyFinder finds only 2 of 3 mandatory genes, and no
T1SS is called. This script unions every replicon of each genome's GCF assembly into one multi-record
GenBank input, so TolC is present and the loner rule can find it (proven for hlyA: chromosome
NZ_CP024997.2 carries TolC at CWB37_RS19950).

Method per accession:
  1. read its cached GenBank, take the assembly (GCF) from the DBLINK record;
  2. eutils: assembly -> all RefSeq nuccore replicons;
  3. fetch any missing replicon (gbwithparts) into refseq_cache;
  4. concatenate all replicon GenBank flatfiles into inputs_gb_fullasm/<acc>.gbff (chromosome first).

Matches the panel's annotation path (input GenBank annotations preserved, no Bakta), so the rerun is
apples-to-apples with the existing recall numbers.

Inputs : data/refseq_cache/<acc>.gb  (the staged single replicon, for the DBLINK)
Outputs: inputs_gb_fullasm/<acc>.gbff (full assembly, multi-record), data/refseq_cache/<new>.gb
Run    : .venv/bin/python scripts/50_stage_full_assembly.py
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "refseq_cache"
OUT = ROOT / "inputs_gb_fullasm"
EUT = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA = "ssign-benchmark/0.1 (teoreid@gmail.com)"

# placement T1SS genomes whose staged unit is a single partial replicon (missing the OMF replicon)
TARGETS = {
    "NZ_CP031766.1": "hlyA",
    "NZ_CBDBTK010000022.1": "apxIA",
    "NZ_JABJZG010000001.1": "ltxA",
    "NZ_SMAM01000003.1": "lktA",
    "NZ_JBCGCZ010000007.1": "Serralysin",
    "NC_010939.1": "apxIIA",
}


def _get(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode()


def assembly_of(acc: str) -> str:
    """Pull the GCF assembly accession from the cached replicon's DBLINK line."""
    for line in (CACHE / f"{acc}.gb").read_text().splitlines():
        if "Assembly:" in line:
            return line.split("Assembly:")[1].strip()
        if line.startswith("FEATURES"):
            break
    return ""


def replicons_of(asm: str) -> list[str]:
    """assembly accession -> its RefSeq nuccore replicon accessions (chromosome + plasmids)."""
    d = json.loads(_get(EUT + "esearch.fcgi?db=assembly&term=" + asm + "&retmode=json"))
    auid = d["esearchresult"]["idlist"][0]
    time.sleep(0.4)
    d = json.loads(
        _get(EUT + f"elink.fcgi?dbfrom=assembly&db=nuccore&id={auid}&retmode=json&linkname=assembly_nuccore_refseq")
    )
    uids = [link for ls in d["linksets"] for db in ls.get("linksetdbs", []) for link in db["links"]]
    time.sleep(0.4)
    d = json.loads(_get(EUT + "esummary.fcgi?db=nuccore&id=" + ",".join(uids) + "&retmode=json"))
    res = d["result"]
    # chromosome first (largest), then plasmids -> stable ordered input
    accs = sorted(((res[u]["accessionversion"], int(res[u].get("slen", 0))) for u in res["uids"]), key=lambda x: -x[1])
    return [a for a, _ in accs]


def fetch(acc: str) -> Path:
    out = CACHE / f"{acc}.gb"
    if out.exists() and out.stat().st_size > 0:
        return out
    params = urllib.parse.urlencode({"db": "nuccore", "id": acc, "rettype": "gbwithparts", "retmode": "text"})
    out.write_text(_get(EFETCH + "?" + params))
    time.sleep(0.4)
    return out


def main() -> int:
    OUT.mkdir(exist_ok=True)
    for acc, gene in TARGETS.items():
        asm = assembly_of(acc)
        if not asm:
            print(f"  {acc} ({gene}): NO assembly in DBLINK -- skipped")
            continue
        reps = replicons_of(asm)
        parts = [fetch(r).read_text() for r in reps]
        (OUT / f"{acc}.gbff").write_text("".join(parts))
        n_cds = sum(p.count("/translation=") for p in parts)
        print(f"  {acc} ({gene}): {asm} -> {len(reps)} replicons ({', '.join(reps)}), {n_cds} CDS")
    print(f"\nfull-assembly inputs in {OUT}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
