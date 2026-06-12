#!/usr/bin/env python3
"""
13_rescue_t1ss_blast.py  (Phase 0a augmentation: rescue unplaced T1SS effectors, pass 2)

The 10 effectors that pass 1 (IPG, script 12) could not place at 100% identity --
their characterized sequence is a strain variant not present verbatim in any genome.
Per Teo's >=95% representative-strain decision, place them by sequence similarity:

  1. Remote NCBI blastp of the UniProt protein against refseq_protein, restricted to
     the effector's own species (ENTREZ_QUERY txid<taxid>[ORGN]).
  2. Keep the best hit with pident >= 95 AND query coverage >= 90 (a near-full-length,
     near-identical match = the same gene in a representative genome).
  3. Resolve that matched RefSeq protein to a genome + coordinates via IPG (the same
     pick_placement logic script 12 already uses), preferring complete assemblies.

RefSeq stays a coordinate lookup for an already-verified effector. The original
characterization DOI is preserved; %identity / %coverage / matched accession are all
recorded so the placement is independently checkable. No >=95% full-length hit in the
species -> status=unplaceable (documented, not invented).

Inputs:
  data/t1ss_rescue/ipg_placements.tsv   (rows with status=needs_blast)
  data/t1ss_rescue/sequences.fasta      (UniProt protein sequences)
Outputs:
  data/t1ss_rescue/blast_placements.tsv
  data/t1ss_rescue/.blast_rids.json     (RID cache -> resumable)

Run:
  .venv/bin/python scripts/13_rescue_t1ss_blast.py
"""

from __future__ import annotations

import csv
import importlib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
ipg_mod = importlib.import_module("12_rescue_t1ss_ipg")  # reuse ipg_rows / pick_placement / assembly_level

RESCUE = BENCH / "data" / "t1ss_rescue"
IPG_TSV = RESCUE / "ipg_placements.tsv"
FASTA = RESCUE / "sequences.fasta"
OUT_TSV = RESCUE / "blast_placements.tsv"
RID_CACHE = RESCUE / ".blast_rids.json"

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
UA = "ssign-benchmark/0.1 (teoreid@gmail.com)"

# Floor lowered 95 -> 90 (Teo, 2026-06-11) to admit prtA/Q07295, a same-species 93%
# ortholog. The true per-row identity is recorded, so the bulk (98-100%) stays distinguishable.
PIDENT_MIN = 90.0
QCOV_MIN = 90.0

FIELDS = [
    "gene",
    "uniprot",
    "organism",
    "taxid",
    "status",
    "tier",
    "matched_protein",
    "pident",
    "qcov",
    "evalue",
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


def load_fasta(path):
    seqs, acc, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if acc:
                seqs[acc] = "".join(buf)
            acc = line[1:].split("|")[0]
            buf = []
        else:
            buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def _open_retry(req_or_url, timeout, tries=4):
    """urlopen with backoff; the BLAST queue intermittently drops connections."""
    last = None
    for i in range(tries):
        try:
            return urllib.request.urlopen(req_or_url, timeout=timeout).read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001 - transient network, retry
            last = e
            time.sleep(5 * (i + 1))
    raise last


def blast_submit(seq, taxid):
    data = urllib.parse.urlencode(
        {
            "CMD": "Put",
            "PROGRAM": "blastp",
            "DATABASE": "refseq_protein",
            "QUERY": seq,
            "ENTREZ_QUERY": f"txid{taxid}[ORGN]",
            "HITLIST_SIZE": "50",
            "EXPECT": "1e-20",
        }
    ).encode()
    req = urllib.request.Request(BLAST_URL, data=data, headers={"User-Agent": UA})
    txt = _open_retry(req, 60)
    rid = None
    for line in txt.splitlines():
        if line.strip().startswith("RID = "):
            rid = line.split("=", 1)[1].strip()
            break
    return rid


def blast_ready(rid):
    q = urllib.parse.urlencode({"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"})
    try:
        txt = _open_retry(f"{BLAST_URL}?{q}", 60)
    except Exception as e:  # noqa: BLE001 - persistent failure: keep polling, don't crash the run
        print(f"  poll error (treating as WAITING): {e}", file=sys.stderr)
        return "WAITING"
    for line in txt.splitlines():
        if "Status=" in line:
            return line.split("Status=", 1)[1].strip()
    return "UNKNOWN"


def blast_report(rid):
    # FORMAT_TYPE=Tabular is broken via this API (returns an empty wrapper even when
    # ThereAreHits=yes); the classic Text report is what actually carries the hits.
    q = urllib.parse.urlencode(
        {"CMD": "Get", "RID": rid, "FORMAT_TYPE": "Text", "ALIGNMENTS": "50", "DESCRIPTIONS": "50"}
    )
    return _open_retry(f"{BLAST_URL}?{q}", 120)


def alignment_stats(lines, acc, qlen):
    """For a subject's alignment block: (best-HSP true %identity, union query coverage).

    Reads the real `Identities = a/b` per HSP (unrounded, not the rounded summary
    integer) and unions the `Query <start> ... <end>` ranges across all HSPs, so a
    multi-HSP RTX toxin isn't undercounted. Returns (None, 0.0) if no block found.
    """
    i = next((j for j, l in enumerate(lines) if l.startswith(">") and acc in l), None)
    if i is None:
        return None, 0.0
    best_pid, qlo, qhi = None, None, None
    j = i + 1
    while j < len(lines) and not lines[j].startswith(">"):
        l = lines[j]
        mi = re.search(r"Identities\s*=\s*(\d+)/(\d+)", l)
        if mi:
            pid = 100.0 * int(mi.group(1)) / int(mi.group(2))
            best_pid = pid if best_pid is None else max(best_pid, pid)
        mq = re.match(r"Query\s+(\d+)\s+\S+\s+(\d+)", l)
        if mq:
            s, e = sorted((int(mq.group(1)), int(mq.group(2))))
            qlo = s if qlo is None else min(qlo, s)
            qhi = e if qhi is None else max(qhi, e)
        j += 1
    if best_pid is None:
        return None, 0.0
    qcov = 100.0 * (qhi - qlo + 1) / qlen if (qlo and qhi and qlen) else 0.0
    return best_pid, qcov


def parse_best_hit(report, qlen):
    """Best hit passing the real >=PIDENT_MIN identity and >=QCOV_MIN coverage gates.

    The summary table ('Sequences producing significant alignments') is used only to
    RANK candidates by its (rounded) max-ident column; the actual gate is applied on
    each candidate's alignment block via alignment_stats, so a 94.6%-rounded-to-95%
    summary value cannot sneak through.
    """
    lines = report.splitlines()
    try:
        start = next(i for i, l in enumerate(lines) if l.startswith("Sequences producing significant"))
    except StopIteration:
        return None

    cands = []  # (rounded ident%, accession)
    for l in lines[start + 1 :]:
        if not l.strip():
            if cands:
                break
            continue
        m = re.match(r"^(\S+)\s+.*?\s(\d+)%\s*$", l)
        if m:
            cands.append((int(m.group(2)), m.group(1)))
    # widen the pre-filter by 1% so a true 95.x rounded down to 95 (or 94 shown) is not lost
    cands = sorted((c for c in cands if c[0] >= PIDENT_MIN - 1), reverse=True)

    for _, acc in cands:
        pid, qcov = alignment_stats(lines, acc, qlen)
        if pid is not None and pid >= PIDENT_MIN and qcov >= QCOV_MIN:
            return (pid, qcov, acc, "")
    return None


def main():
    rows_in = [r for r in csv.DictReader(open(IPG_TSV), delimiter="\t") if r["status"] == "needs_blast"]
    seqs = load_fasta(FASTA)
    rids = json.loads(RID_CACHE.read_text()) if RID_CACHE.exists() else {}

    # submit (or reuse) one RID per effector
    for r in rows_in:
        acc = r["uniprot"]
        if rids.get(acc):
            continue
        seq = seqs.get(acc, "")
        if not seq or not r["taxid"]:
            rids[acc] = ""
            continue
        rid = blast_submit(seq, r["taxid"])
        rids[acc] = rid or ""
        print(f"  submitted {acc} ({r['gene']}) taxid{r['taxid']} -> RID {rid}", file=sys.stderr)
        RID_CACHE.write_text(json.dumps(rids, indent=0))
        time.sleep(3)  # be polite to the BLAST queue

    # poll all to completion
    pending = {r["uniprot"] for r in rows_in if rids.get(r["uniprot"])}
    waited = 0
    while pending and waited < 1800:
        for acc in list(pending):
            st = blast_ready(rids[acc])
            if st == "READY":
                pending.discard(acc)
            elif st in ("FAILED", "UNKNOWN"):
                print(f"  RID for {acc} status={st}", file=sys.stderr)
                pending.discard(acc)
            time.sleep(2)
        if pending:
            print(f"  waiting on {len(pending)} BLAST jobs... ({waited}s)", file=sys.stderr)
            time.sleep(20)
            waited += 20

    # collect results
    rows = []
    for r in rows_in:
        acc = r["uniprot"]
        row = {f: "" for f in FIELDS}
        row.update(
            gene=r["gene"],
            uniprot=acc,
            organism=r["organism"],
            taxid=r["taxid"],
            primary_ref=r["primary_ref"],
            status="unplaceable",
            tier="",
        )
        rid = rids.get(acc)
        seq = seqs.get(acc, "")
        if not rid or not seq:
            row["status"] = "blast_error"  # never submitted (no RID / no sequence)
        if rid and seq:
            # network fetch can fail transiently -> a distinct status, NOT "no hit";
            # parsing is pure and left to crash loudly if a real bug appears.
            try:
                report = blast_report(rid)
            except Exception as e:  # noqa: BLE001 - network, distinguish from biology
                print(f"  fetch error {acc}: {e}", file=sys.stderr)
                report = None
            best = parse_best_hit(report, len(seq)) if report is not None else None
            if report is None:
                row["status"] = "blast_error"
            if best:
                pident, qcov, subj_acc, evalue = best
                placement, other = ipg_mod.pick_placement(ipg_mod.ipg_rows(subj_acc), r["organism"])
                row.update(matched_protein=subj_acc, pident=f"{pident:.1f}", qcov=f"{qcov:.0f}", evalue=evalue)
                if placement:
                    lvl_n, lvl_s = ipg_mod.assembly_level(placement["Nucleotide Accession"])
                    row.update(
                        status="placed_blast",
                        tier="representative_strain",
                        placement_assembly=placement["Assembly"],
                        placement_nuc=placement["Nucleotide Accession"],
                        start=placement.get("Start", ""),
                        stop=placement.get("Stop", ""),
                        strand=placement.get("Strand", ""),
                        placement_protein=placement.get("Protein", ""),
                        placement_strain=placement.get("Strain", "") or placement.get("Organism", ""),
                        placement_organism=placement.get("Organism", ""),
                        assembly_level=lvl_s,
                        other_assemblies=";".join(other),
                    )
                else:
                    row["status"] = "hit_no_genome"  # >=95% hit exists but only in a clone, not a genome
            time.sleep(0.4)
        rows.append(row)
        print(
            f"  {r['gene']:10s} {acc:8s} -> {row['status']:14s} "
            f"id={row['pident']} cov={row['qcov']} {row['placement_nuc']}",
            file=sys.stderr,
        )

    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    placed = [r for r in rows if r["status"] == "placed_blast"]
    print(f"\nwrote {OUT_TSV.relative_to(BENCH)}")
    print(f"  needs_blast effectors : {len(rows)}")
    print(f"  placed (>={PIDENT_MIN:.0f}% id)     : {len(placed)}")
    for r in rows:
        print(
            f"    {r['gene']:10s} {r['uniprot']:8s} {r['status']:14s} "
            f"id={r['pident'] or '-':>5s} cov={r['qcov'] or '-':>3s} "
            f"{r['placement_nuc']:18s} {r['placement_strain'][:28]}"
        )


if __name__ == "__main__":
    sys.exit(main())
