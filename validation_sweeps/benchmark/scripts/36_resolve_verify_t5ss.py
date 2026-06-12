#!/usr/bin/env python3
"""Dataset group 3, tasks 3.3 + 3.4: resolve + verify agent-sourced T5SS effectors.

The 3.2 sourcing agents proposed T5SS examples (gene, organism, refseq_genome, locus_tag,
subtype, self_secreted, DOI, verbatim quote) in data/dataset/t5ss_raw/*.json. This script
holds each to the same bar as the benchmark machinery answer key (script 11):

  3.3 resolve  -- ensure the cited genome is cached (fetch via 09_fetch_refseq if not),
                  then locate the locus_tag in that genome's NCBI annotation via
                  bench_index, pulling contig/start/end/strand. Falls back to a UNIQUE
                  /gene symbol match when the locus_tag scheme is absent. No coordinates
                  found -> unplaceable (never fabricated).
  3.4 verify   -- DOI resolves (DOI.org Handle API, script 03 doi_resolves) + locus exists
                  in genome + gene named in the verbatim quote. verified := placed AND
                  doi_resolves (gene-in-quote recorded but a collective/operon quote is not
                  a hard fail, mirrors script 11).

Inputs:
  data/dataset/t5ss_raw/*.json        (agent output)
Outputs:
  data/dataset/t5ss_effectors.tsv     (placed rows + verified flag, gold-set-compatible cols)
  data/dataset/t5ss_unplaceable.tsv   (no resolvable locus_tag/coords, with reason)
"""

from __future__ import annotations

import importlib.util
import json
import time
from collections import Counter
from pathlib import Path

from bench_io import norm_doi, write_tsv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "dataset"
RAW = OUT / "t5ss_raw"
CACHE = ROOT / "data" / "refseq_cache"
DOI_CACHE = OUT / ".doi_cache.json"


def _load(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).parent / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bi = _load("bench_index", "bench_index.py")
nine = _load("fetch_refseq9", "09_fetch_refseq.py")
fz = _load("finalize_gold_set", "03_finalize_gold_set.py")

from Bio import Entrez  # noqa: E402

Entrez.email = "teoreid@gmail.com"  # NCBI etiquette; without it esearch/elink get throttled


def replicons_for(genome_acc: str, cache: dict) -> list[str]:
    """Nucleotide RefSeq accessions for a cited genome. A GCF_/GCA_ assembly accession is
    expanded to its replicons via Entrez (bench_index/09 key on nucleotide accessions, not
    assembly); a nucleotide accession passes through unchanged.

    Uses assembly->elink->nuccore (canonical) and retries: a bare nuccore "[Assembly]"
    search is silently throttled by NCBI without an API key, which left some assemblies
    unexpanded on the first pass."""

    if not (genome_acc.startswith("GCF_") or genome_acc.startswith("GCA_")):
        return [genome_acc]
    if genome_acc in cache:
        return cache[genome_acc]

    accs: list[str] = []
    for attempt in range(3):
        try:
            time.sleep(0.4)  # stay under NCBI's 3 req/s without an API key
            # field-tag the accession so a bare term can't match a second (e.g. suppressed) assembly
            h = Entrez.esearch(db="assembly", term=f"{genome_acc}[Assembly Accession]")
            uids = Entrez.read(h)["IdList"]
            h.close()
            if not uids:
                break
            if len(uids) > 1:  # ambiguous -> don't silently expand the wrong assembly's replicons
                print(f"  WARN {genome_acc}: {len(uids)} assembly UIDs, expected 1; skipping expansion")
                break
            time.sleep(0.4)
            h = Entrez.elink(dbfrom="assembly", db="nuccore", id=uids[0], linkname="assembly_nuccore_refseq")
            links = Entrez.read(h)
            h.close()
            nuc_uids = [lk["Id"] for ls in links[0].get("LinkSetDb", []) for lk in ls["Link"]]
            if not nuc_uids:
                break
            time.sleep(0.4)
            h = Entrez.esummary(db="nuccore", id=",".join(nuc_uids))
            accs = [d["AccessionVersion"] for d in Entrez.read(h)]
            h.close()
            break
        except Exception as e:  # noqa: BLE001  network/parse — retry, then leave empty
            print(f"  WARN assembly->nuccore attempt {attempt + 1} failed for {genome_acc}: {e}")
            time.sleep(1.5 * (attempt + 1))
    cache[genome_acc] = accs
    return accs


OUT_COLS = [
    "gene",
    "uniprot",
    "locus_tag",
    "organism",
    "refseq_genome",
    "ss_type",
    "subtype",
    "self_secreted",
    "contig",
    "start",
    "stop",
    "strand",
    "primary_ref",
    "quote",
    "uniprot_note",
    "locus_match",
    "doi_resolves",
    "gene_in_quote",
    "verified",
]


def main() -> int:
    raw: list[dict] = []
    for f in sorted(RAW.glob("*.json")):
        try:
            raw.extend(json.loads(f.read_text()))
        except (json.JSONDecodeError, FileNotFoundError):
            print(f"  WARN could not parse {f.name}, skipping")
    if not raw:
        raise SystemExit(f"no entries in {RAW}/*.json -- run the 3.2 sourcing agents first")

    # 3.3a resolve each cited genome to nucleotide replicons (GCF assembly -> NC_/NZ_), fetch uncached
    ass_cache: dict = {}
    cited = sorted({(r.get("refseq_genome") or "").strip() for r in raw if (r.get("refseq_genome") or "").strip()})
    nuc_for: dict[str, list[str]] = {g: replicons_for(g, ass_cache) for g in cited}
    all_nuc = sorted({a for accs in nuc_for.values() for a in accs})
    missing = [a for a in all_nuc if not (CACHE / f"{a}.gb").exists()]
    if missing:
        print(f"fetching {len(missing)} uncached T5SS replicons ...")
        for a in missing:
            print(f"  {'ok ' if nine.fetch(a) else 'FAIL'} {a}")

    idx = bi.build_from_genbank(all_nuc)
    doi_cache: dict = json.loads(DOI_CACHE.read_text()) if DOI_CACHE.exists() else {}

    placed, unplaceable = [], []
    for r in raw:
        genome = (r.get("refseq_genome") or "").strip()
        locus = (r.get("locus_tag") or "").strip()
        gene = (r.get("gene") or "").strip()

        hit, match = None, "none"
        for nuc in nuc_for.get(genome, [genome]):  # try each replicon of the cited genome
            if locus and (by_locus := idx.find(nuc, locus)):
                hit, match = by_locus, "locus_tag"
                break
            if by_gene := idx.find_by_gene(nuc, gene):
                hit, match = by_gene, "gene_symbol"
                break

        doi = norm_doi(r.get("primary_ref", ""))
        doi_ok = bool(doi) and fz.doi_resolves(doi, doi_cache)
        quote = r.get("quote") or ""
        gene_in_quote = bool(gene) and bi.normalize(gene) in bi.normalize(quote)

        if hit is None:
            unplaceable.append({**r, "reason": f"locus_tag {locus!r}/gene {gene!r} not found in {genome}"})
            continue

        rec_acc, _ordinal, cds = hit
        placed.append(
            {
                "gene": gene,
                "uniprot": r.get("uniprot", ""),
                "locus_tag": cds["locus_tag"] or locus,
                "organism": r.get("organism", ""),
                "refseq_genome": genome,
                "ss_type": "T5SS",
                "subtype": r.get("subtype", ""),
                "self_secreted": str(r.get("self_secreted", "")).lower(),
                "contig": rec_acc,
                "start": cds["start"],
                "stop": cds["end"],
                "strand": cds["strand"],
                "primary_ref": r.get("primary_ref", ""),
                "quote": quote,
                "uniprot_note": r.get("note", ""),
                "locus_match": match,
                "doi_resolves": "yes" if doi_ok else "no",
                "gene_in_quote": "yes" if gene_in_quote else "no",
                "verified": "yes" if (doi_ok and match == "locus_tag") else "review",
            }
        )

    DOI_CACHE.write_text(json.dumps(doi_cache, indent=0))
    write_tsv(OUT / "t5ss_effectors.tsv", OUT_COLS, placed)
    write_tsv(
        OUT / "t5ss_unplaceable.tsv",
        ["gene", "organism", "refseq_genome", "locus_tag", "subtype", "primary_ref", "reason"],
        unplaceable,
    )

    print(f"\nsourced: {len(raw)}   placed: {len(placed)}   unplaceable: {len(unplaceable)}")
    print(f"  by subtype (placed): {dict(sorted(Counter(p['subtype'] for p in placed).items()))}")
    print(
        f"  verified: {sum(1 for p in placed if p['verified'] == 'yes')}   "
        f"review: {sum(1 for p in placed if p['verified'] == 'review')}"
    )
    print(f"  self_secreted=true: {sum(1 for p in placed if p['self_secreted'] == 'true')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
