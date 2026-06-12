#!/usr/bin/env python3
"""
bench_index.py  (shared: alias-aware gene-order index over the cached genomes)

The Phase 1 ceiling analysis and the Phase 2 recall bridge both need the same primitive:
"given a genome accession and a locus_tag, where does that gene sit in gene order, and how
many genes away is another locus on the same replicon?" This module is that primitive.

Two drift problems it absorbs (both documented in tasks.md Checkpoint A 4.0):
  - locus_tag underscore drift: the corpus drops the underscore RefSeq uses
    (corpus `ECs4550` vs RefSeq `ECs_4550`, `CBU0041` vs `CBU_0041`). normalize() folds
    case and removes underscores/spaces so the two forms collide. Within one replicon
    locus_tags stay unique under that fold, so no spurious merging.
  - genome-accession version + RefSeq-prefix drift: corpus `NC_002516` vs cache
    `NC_002516.2`; corpus INSDC `HG326223` vs cache RefSeq `NZ_HG326223.1`. Each cached
    record is registered under its full id, its version-stripped id, and (for NZ_/NC_
    records) the prefix+version-stripped INSDC form, so any of those resolves it.

Gene order is per replicon (one GenBank record): all CDS sorted by start coordinate get a
0-based ordinal. Distance between two loci on the same replicon = |ordinal difference|,
which is exactly ssign's "+/-N genes" notion. Loci on different replicons have no
gene-order distance (return None -> the caller treats that as out-of-reach).

Build once from GenBank (slow, ~92 files) into a flat TSV via 18_build_gene_order.py;
downstream scripts load_from_tsv() (fast).
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
CACHE_DIR = BENCH / "data" / "refseq_cache"
INDEX_TSV = BENCH / "data" / "phase1" / "gene_order_index.tsv"

NAME_QUALS = ("gene", "gene_synonym", "old_locus_tag", "locus_tag")
INDEX_FIELDS = ["record_acc", "ordinal", "locus_tag", "gene", "aliases", "start", "end", "strand"]


def load_tsv(path) -> list[dict]:
    """Read a tab-separated file into a list of dict rows (shared by the Phase 1/2 scripts)."""
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def normalize(tag: str) -> str:
    """Fold a locus_tag/gene name for drift-tolerant matching: lowercase, drop _ and spaces."""
    return re.sub(r"[_\s]", "", (tag or "")).lower()


def strip_version(acc: str) -> str:
    return acc.split(".")[0].strip()


def strip_prefix(acc: str) -> str:
    """Drop a leading RefSeq prefix (NZ_/NC_/NT_/NW_) so the bare INSDC core remains."""
    return re.sub(r"^(NZ|NC|NT|NW)_", "", acc)


def accession_base(acc: str) -> str:
    """Version+prefix-stripped, lowercased accession base (NC_002516.2 -> 002516) for
    drift-tolerant grouping of replicon spellings."""
    return strip_prefix(strip_version(acc or "")).lower()


def _acc_keys(acc: str) -> set[str]:
    """All forms under which a record accession should be findable."""
    v = strip_version(acc)
    return {acc, v, strip_prefix(v)}


class GeneOrderIndex:
    """Per-replicon ordered CDS lists + alias->ordinal maps + a drift-tolerant accession resolver."""

    def __init__(self):
        self.records: dict[str, list[dict]] = {}  # record_acc -> ordered CDS dicts
        self.alias2ord: dict[str, dict[str, int]] = {}  # record_acc -> {normalized alias: ordinal}
        self.gene2ords: dict[str, dict[str, list[int]]] = {}  # record_acc -> {normalized /gene: [ordinals]}
        self._resolve: dict[str, str] = {}  # any acc form -> canonical record_acc

    def _register(self, record_acc: str):
        for k in _acc_keys(record_acc):
            # full/exact ids win over shared prefix-stripped collisions
            self._resolve.setdefault(k, record_acc)
        self._resolve[record_acc] = record_acc

    def resolve_record(self, genome_acc: str) -> str | None:
        g = genome_acc.strip()
        for k in (g, strip_version(g), strip_prefix(strip_version(g))):
            if k in self._resolve:
                return self._resolve[k]
        return None

    def add_record(self, record_acc: str, cds: list[dict]):
        cds = sorted(cds, key=lambda d: (d["start"], d["end"]))
        a2o: dict[str, int] = {}
        g2o: dict[str, list[int]] = defaultdict(list)
        for i, d in enumerate(cds):
            d["ordinal"] = i
            for a in d["aliases"]:
                a2o.setdefault(normalize(a), i)
            if d["gene"]:
                g2o[normalize(d["gene"])].append(i)
        self.records[record_acc] = cds
        self.alias2ord[record_acc] = a2o
        self.gene2ords[record_acc] = g2o
        self._register(record_acc)

    def find(self, genome_acc: str, locus_tag: str):
        """Return (record_acc, ordinal, cds_dict) for a locus, or None if not found."""
        rec = self.resolve_record(genome_acc)
        if rec is None:
            return None
        o = self.alias2ord[rec].get(normalize(locus_tag))
        if o is None:
            return None
        return rec, o, self.records[rec][o]

    def find_by_gene(self, genome_acc: str, gene: str):
        """Fallback when a locus_tag scheme is absent from the assembly: locate by /gene
        symbol, but ONLY when that symbol is unique on the replicon (no paralog ambiguity).
        Returns (record_acc, ordinal, cds_dict) or None.
        """
        rec = self.resolve_record(genome_acc)
        if rec is None or not gene:
            return None
        ords = self.gene2ords[rec].get(normalize(gene))
        if not ords or len(ords) != 1:
            return None
        o = ords[0]
        return rec, o, self.records[rec][o]

    def gene_distance(self, genome_acc: str, locus_a: str, locus_b: str):
        """|ordinal difference| if both loci sit on the same resolved replicon, else None."""
        fa, fb = self.find(genome_acc, locus_a), self.find(genome_acc, locus_b)
        if not fa or not fb or fa[0] != fb[0]:
            return None
        return abs(fa[1] - fb[1])


def build_from_genbank(accessions=None) -> GeneOrderIndex:
    """Parse cached GenBank records into an index. accessions=None -> every cached .gb."""
    from Bio import SeqIO

    idx = GeneOrderIndex()
    files = sorted(CACHE_DIR.glob("*.gb"))
    want = {strip_version(a) for a in accessions} if accessions else None
    for gb in files:
        for rec in SeqIO.parse(str(gb), "genbank"):
            if want is not None and strip_version(rec.id) not in want and strip_version(rec.name) not in want:
                continue
            # merge gene+CDS features sharing a locus_tag so aliases from both land together
            by_locus: dict[str, dict] = {}
            for feat in rec.features:
                if feat.type not in ("gene", "CDS"):
                    continue
                q = feat.qualifiers
                lt = q.get("locus_tag", [None])[0]
                key = lt or f"{rec.id}:{int(feat.location.start)}"
                e = by_locus.setdefault(
                    key,
                    {
                        "locus_tag": lt or "",
                        "gene": q.get("gene", [""])[0],
                        "start": int(feat.location.start) + 1,
                        "end": int(feat.location.end),
                        "strand": feat.location.strand or 0,
                        "aliases": set(),
                    },
                )
                if feat.type == "CDS" and not e["gene"]:
                    e["gene"] = q.get("gene", [""])[0]
                for ql in NAME_QUALS:
                    for val in q.get(ql, []):
                        e["aliases"].add(val)
            # keep only loci with a real CDS-bearing locus_tag (drop bare-coordinate keys)
            cds = [d for d in by_locus.values() if d["locus_tag"]]
            if cds:
                idx.add_record(rec.id, cds)
    return idx


def load_from_tsv(path: Path = INDEX_TSV) -> GeneOrderIndex:
    """Reconstruct the index from the flat TSV emitted by 18_build_gene_order.py."""
    idx = GeneOrderIndex()
    by_rec: dict[str, list[dict]] = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_rec[row["record_acc"]].append(
                {
                    "locus_tag": row["locus_tag"],
                    "gene": row["gene"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "strand": int(row["strand"]),
                    "aliases": set(a for a in row["aliases"].split(";") if a),
                }
            )
    for rec, cds in by_rec.items():
        idx.add_record(rec, cds)  # re-sorts by start; ordinals already consistent
    return idx
