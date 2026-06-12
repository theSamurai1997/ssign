#!/usr/bin/env python3
"""Fetch RefSeq GenBank files for the 6 benchmarking genomes from NCBI.

Caches to validation_sweeps/analysis/refseq_cache/. Idempotent.

These are the same input files ssign was run on. We re-fetch them locally
so we can build a coordinate-based bridge from Bakta locus_tag (in ssign
output) back to the original RefSeq locus_tag.
"""

from __future__ import annotations

import time
from pathlib import Path

from Bio import Entrez

ANALYSIS = Path(__file__).resolve().parents[1] / "analysis"
CACHE = ANALYSIS / "refseq_cache"
CACHE.mkdir(exist_ok=True)

Entrez.email = "teoreid@gmail.com"  # NCBI etiquette: identify yourself

# Each sample -> list of RefSeq nucleotide accessions (some genomes have
# multiple replicons / plasmids).
GENOMES = {
    "legionella_pneumophila": ["NC_002942.5"],
    "coxiella_rsa493": ["NC_002971.4"],
    "salmonella_lt2": ["NC_003197.2", "NC_003277.2"],  # chromosome + pSLT
    "yersinia_pestis_co92": ["NC_003143.1", "NC_003131.1", "NC_003134.1", "NC_003132.1"],
    "pseudomonas_pao1": ["NC_002516.2"],
    "vibrio_cholerae_n16961": ["NC_002505.1", "NC_002506.1"],
}


def fetch(sample: str, accs: list[str]) -> Path:
    out = CACHE / f"{sample}.gb"
    if out.exists() and out.stat().st_size > 0:
        return out
    print(f"  fetching {sample} ({', '.join(accs)})")
    with out.open("w") as f:
        for acc in accs:
            try:
                handle = Entrez.efetch(db="nuccore", id=acc, rettype="gbwithparts", retmode="text")
                f.write(handle.read())
                handle.close()
            except Exception as e:
                print(f"    SKIPPED {acc}: {e}")
            time.sleep(1)
    return out


def main() -> int:
    for sample, accs in GENOMES.items():
        path = fetch(sample, accs)
        print(f"  {sample}: {path.stat().st_size // 1024} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
