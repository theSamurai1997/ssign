#!/usr/bin/env python3
"""
12_rescue_t1ss_ipg.py  (Phase 0a augmentation: rescue unplaced T1SS effectors)

The 19 T1SS gold-set effectors with no `refseq_genome` carry a UniProt accession
but were curated from 1980s-90s EMBL gene-clone deposits, never from a genome
assembly -- so there is no gene-order context to measure proximity against.

Teo's decision (Checkpoint A follow-up, 2026-06-10): place them into a
representative same-species genome, sequence-confirmed (>=95% identity), flagged
as a `representative_strain` tier. This script does the cheap, exact first pass:

  NCBI Identical Protein Groups (IPG). For each effector, IPG lists every genome
  + coordinate where a *100%-identical* protein occurs. If that identical protein
  already sits in a RefSeq complete/chromosome-level genome, that is the strongest
  possible placement (`ipg_identical` tier, no BLAST needed). RTX toxins / HlyA
  vary across strains, so IPG returns only their old clones -> those fall through
  to the >=95% BLAST pass (next script), listed here as `needs_blast`.

RefSeq is used ONLY as a coordinate lookup for an already-verified effector, never
to discover effectors (same locked rule as the machinery resolver). The effector's
original characterization DOI is preserved; the genome assembly is recorded as the
coordinate source.

Inputs:
  data/gold_build/effector_gold_set.tsv   (the 19 are ss_type=T1SS, empty refseq_genome)
Outputs:
  data/t1ss_rescue/ipg_placements.tsv      one row per effector (placed or not)
  data/t1ss_rescue/sequences.fasta         UniProt protein seq per effector (for the BLAST pass)
  data/t1ss_rescue/.uniprot_cache/<acc>.txt  cached UniProt entries

Run:
  .venv/bin/python scripts/12_rescue_t1ss_ipg.py
"""

from __future__ import annotations

import csv
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
GOLD = BENCH / "data" / "gold_build" / "effector_gold_set.tsv"
OUT_DIR = BENCH / "data" / "t1ss_rescue"
UNIPROT_CACHE = OUT_DIR / ".uniprot_cache"
OUT_TSV = OUT_DIR / "ipg_placements.tsv"
OUT_FASTA = OUT_DIR / "sequences.fasta"

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
UA = "ssign-benchmark/0.1 (teoreid@gmail.com)"

FIELDS = [
    "gene",
    "uniprot",
    "organism",
    "taxid",
    "query_protein_id",
    "status",
    "tier",
    "placement_assembly",
    "placement_nuc",
    "start",
    "stop",
    "strand",
    "placement_protein",
    "placement_strain",
    "placement_organism",
    "assembly_level",
    "n_genome_options",
    "other_assemblies",
    "primary_ref",
]


def _get(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    return urllib.request.urlopen(req, timeout=timeout).read().decode("utf-8", "replace")


def uniprot_entry(acc: str) -> str:
    cache = UNIPROT_CACHE / f"{acc}.txt"
    if cache.exists():
        return cache.read_text()
    txt = _get(f"https://rest.uniprot.org/uniprotkb/{acc}.txt")
    cache.write_text(txt)
    time.sleep(0.2)
    return txt


def parse_uniprot(txt: str):
    """Return (taxid, first_embl_protein_id, sequence)."""
    taxid, prot_id, seq, in_seq = "", "", [], False
    for line in txt.splitlines():
        if line.startswith("OX   NCBI_TaxID="):
            taxid = line.split("=", 1)[1].split(maxsplit=1)[0].rstrip(";").strip()
        elif line.startswith("DR   EMBL;") and not prot_id:
            # DR   EMBL; <nuc>; <protein_id>; <status>; <moltype>.
            parts = [p.strip() for p in line[len("DR   EMBL;") :].split(";")]
            if len(parts) >= 2 and parts[1] not in ("-", ""):
                prot_id = parts[1]
        elif line.startswith("SQ   "):
            in_seq = True
        elif in_seq and line.startswith("     "):
            seq.append(line.strip().replace(" ", ""))
        elif line.startswith("//"):
            in_seq = False
    return taxid, prot_id, "".join(seq)


def ipg_rows(query_id: str):
    """Fetch the IPG report for a protein; return list of dict rows."""
    params = urllib.parse.urlencode(
        {
            "db": "ipg",
            "id": query_id,
            "rettype": "ipg",
            "retmode": "text",
            "tool": "ssign-benchmark",
            "email": "teoreid@gmail.com",
        }
    )
    txt = _get(f"{EFETCH}?{params}")
    time.sleep(0.4)
    lines = [ln for ln in txt.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith("Id\t"):
        return []
    hdr = lines[0].split("\t")
    out = []
    for ln in lines[1:]:
        vals = ln.split("\t")
        if len(vals) < len(hdr):
            vals += [""] * (len(hdr) - len(vals))
        out.append(dict(zip(hdr, vals)))
    return out


def assembly_level(nuc: str) -> tuple[int, str]:
    """Rank a nucleotide accession by assembly completeness (lower = better)."""
    if nuc.startswith("NC_"):
        return (0, "refseq_complete")
    if nuc.startswith(("NZ_CP", "NZ_LR", "NZ_LT", "NZ_LN", "NZ_OU")):
        return (1, "refseq_complete_wgs")
    if nuc.startswith("NZ_CM"):
        return (2, "refseq_chromosome_from_wgs")
    if nuc.startswith("NZ_"):
        return (3, "refseq_wgs_contig")
    if nuc.startswith(("CP", "LR", "LT")):
        return (4, "insdc_complete")
    return (5, "other")


def pick_placement(rows, uniprot_strain_hint: str):
    """Choose the best RefSeq genome row carrying the identical protein."""
    refseq = [
        r
        for r in rows
        if r.get("Source") == "RefSeq" and r.get("Assembly", "").strip() and r.get("Nucleotide Accession", "").strip()
    ]
    if not refseq:
        return None, []
    hint = (uniprot_strain_hint or "").lower()

    def key(r):
        lvl, _ = assembly_level(r["Nucleotide Accession"])
        strain_match = 0 if hint and hint in (r.get("Strain", "") + " " + r.get("Organism", "")).lower() else 1
        return (strain_match, lvl, r["Nucleotide Accession"])

    refseq.sort(key=key)
    best = refseq[0]
    other = sorted({r["Assembly"] for r in refseq if r["Assembly"] != best["Assembly"]})
    return best, other


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    UNIPROT_CACHE.mkdir(parents=True, exist_ok=True)

    strip = lambda a: (a or "").split(".")[0].strip()
    gold = list(csv.DictReader(open(GOLD), delimiter="\t"))
    unplaced = [
        r
        for r in gold
        if r["ss_type"] == "T1SS" and not (strip(r["refseq_genome"]) and strip(r["refseq_genome"]) != "-")
    ]

    rows, fasta = [], []
    for r in unplaced:
        acc = r["uniprot"].strip()
        txt = uniprot_entry(acc)
        taxid, prot_id, seq = parse_uniprot(txt)
        query = prot_id or acc
        if seq:
            fasta.append(f">{acc}|{r['gene']}|taxid{taxid}\n{seq}")

        ipg = ipg_rows(query)
        best, other = pick_placement(ipg, r["organism"])

        row = {f: "" for f in FIELDS}
        row.update(
            gene=r["gene"],
            uniprot=acc,
            organism=r["organism"],
            taxid=taxid,
            query_protein_id=query,
            primary_ref=r["primary_ref"],
        )
        if best:
            lvl_n, lvl_s = assembly_level(best["Nucleotide Accession"])
            row.update(
                status="placed_ipg",
                tier="ipg_identical",
                placement_assembly=best["Assembly"],
                placement_nuc=best["Nucleotide Accession"],
                start=best.get("Start", ""),
                stop=best.get("Stop", ""),
                strand=best.get("Strand", ""),
                placement_protein=best.get("Protein", ""),
                placement_strain=best.get("Strain", "") or best.get("Organism", ""),
                placement_organism=best.get("Organism", ""),
                assembly_level=lvl_s,
                n_genome_options=len(
                    {rr["Assembly"] for rr in ipg if rr.get("Source") == "RefSeq" and rr.get("Assembly", "").strip()}
                ),
                other_assemblies=";".join(other),
            )
        else:
            row.update(status="needs_blast", tier="", n_genome_options=0)
        rows.append(row)
        print(
            f"  {r['gene']:10s} {acc:8s} -> {row['status']:12s} "
            f"{row['placement_assembly']:18s} {row['placement_strain'][:30]}",
            file=sys.stderr,
        )

    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    OUT_FASTA.write_text("\n".join(fasta) + "\n")

    placed = [r for r in rows if r["status"] == "placed_ipg"]
    need = [r for r in rows if r["status"] == "needs_blast"]
    print(f"\nwrote {OUT_TSV.relative_to(BENCH)}")
    print("  19 unplaced T1SS effectors")
    print(f"  placed by IPG (identical protein in a RefSeq genome): {len(placed)}")
    print(f"  need >=95% BLAST pass                                : {len(need)}")
    if placed:
        print("\n  IPG placements:")
        for r in placed:
            print(
                f"    {r['gene']:10s} {r['uniprot']:8s} -> {r['placement_nuc']:18s} "
                f"({r['assembly_level']}) {r['placement_strain'][:34]}"
            )
    if need:
        print("\n  fall through to BLAST:")
        for r in need:
            print(f"    {r['gene']:10s} {r['uniprot']:8s} {r['organism'][:50]}")


if __name__ == "__main__":
    sys.exit(main())
