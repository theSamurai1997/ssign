#!/usr/bin/env python3
"""
16_t1ss_replace_fragmented.py  (T1SS rescue: re-place effectors stuck on fragmented contigs)

Script 15 flagged some rescued effectors as adjacency NOT_FOUND. Some are genuine biology
(apxIIA: trans-secreted, no operon transporter -- left as-is). Others are assembly
artifacts: the effector landed on a short WGS contig that truncates the operon (hlyA on a
3 kb single-CDS contig), so the adjacent HlyB/HlyD are simply on a neighbouring contig.

This re-places ONLY the artifact cases into a COMPLETE genome where the operon is intact,
then re-verifies adjacency. Method per effector:
  1. From its cached blastp hits (script 13 RIDs), take every subject >= PIDENT_MIN.
  2. IPG each subject -> genome rows; keep those on a RefSeq *complete* replicon
     (NC_ reference, or a chromosome-level NZ_CP/NZ_LR... record).
  3. Fetch each candidate, locate the effector by coordinates, run the script-15 adjacency
     check. Take the first genome where ABC+MFP are CONFIRMED adjacent.
If no complete genome confirms adjacency, the effector stays NOT_FOUND (a genuine
non-operonic case, like apxIIA), not silently forced.

Only operates on rows the caller lists in TARGETS (artifact cases); genuine exceptions are
excluded by name so this never overrides real biology.

Inputs : data/t1ss_rescue/{t1ss_rescued.tsv, .blast_rids.json, sequences.fasta}
Outputs: rewrites the affected rows in t1ss_rescued.tsv, logs every candidate considered
         to data/t1ss_rescue/replacement_log.tsv. Re-run script 15 afterwards to refresh
         t1ss_ceiling.tsv from the updated placements (15 stays the single source).

Run:
  .venv/bin/python scripts/16_t1ss_replace_fragmented.py
  .venv/bin/python scripts/15_t1ss_adjacency.py   # refresh the ceiling table
"""

from __future__ import annotations

import csv
import importlib
import sys
from pathlib import Path

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
blast = importlib.import_module("13_rescue_t1ss_blast")
ipg = importlib.import_module("12_rescue_t1ss_ipg")
adj = importlib.import_module("15_t1ss_adjacency")
finalize = importlib.import_module("14_finalize_t1ss_rescue")  # match_level()

RESCUE = BENCH / "data" / "t1ss_rescue"
RESCUED_TSV = RESCUE / "t1ss_rescued.tsv"
LOG_TSV = RESCUE / "replacement_log.tsv"
CACHE = BENCH / "data" / "refseq_cache"

# effectors to re-place (assembly-artifact NOT_FOUND); genuine exceptions are NOT listed
TARGETS = {"P08715": "hlyA", "P23694": "Serralysin"}


def is_complete_replicon(nuc: str) -> bool:
    # NC_ = curated RefSeq reference (complete); NZ_CP/LR/LT/LN/OU = complete genome WGS.
    return nuc.startswith("NC_") or nuc.startswith(("NZ_CP", "NZ_LR", "NZ_LT", "NZ_LN", "NZ_OU"))


def all_blast_subjects(rid, qlen, pid_floor):
    """Every summary-table subject accession with rounded ident >= floor (best first)."""
    report = blast.blast_report(rid)
    lines = report.splitlines()
    try:
        s = next(i for i, l in enumerate(lines) if l.startswith("Sequences producing significant"))
    except StopIteration:
        return []
    out = []
    for l in lines[s + 1 :]:
        if not l.strip():
            if out:
                break
            continue
        m = blast.re.match(r"^(\S+)\s+.*?\s(\d+)%\s*$", l)
        if m and int(m.group(2)) >= pid_floor:
            out.append(m.group(1))
    return out


# adjacency_for(nuc, start, stop) is reused from script 15 (single source for the scan rule):
adjacency_for = adj.adjacency_for


def main():
    rescued = list(csv.DictReader(open(RESCUED_TSV), delimiter="\t"))
    rids = blast.json.loads((RESCUE / ".blast_rids.json").read_text())
    seqs = blast.load_fasta(RESCUE / "sequences.fasta")

    log = []
    new_placement = {}  # uniprot -> (nuc, start, stop, assembly, organism, ceiling_fields)
    for acc, gene in TARGETS.items():
        rid, seq = rids.get(acc), seqs.get(acc, "")
        if not rid or not seq:
            print(f"  {gene}: no RID/seq, skip", file=sys.stderr)
            continue
        subjects = all_blast_subjects(rid, len(seq), int(blast.PIDENT_MIN))
        print(f"\n{gene} ({acc}): {len(subjects)} blast subjects >= {blast.PIDENT_MIN:.0f}%", file=sys.stderr)

        # gather complete-replicon candidate placements from each subject's IPG
        seen, candidates = set(), []
        for subj in subjects[:25]:
            if len(candidates) >= 15:
                break
            for row in ipg.ipg_rows(subj):
                if row.get("Source") != "RefSeq":
                    continue
                nuc = row.get("Nucleotide Accession", "")
                if not nuc or not is_complete_replicon(nuc) or nuc in seen:
                    continue
                seen.add(nuc)
                candidates.append(
                    {
                        "subject": subj,
                        "assembly": row.get("Assembly", ""),
                        "nuc": nuc,
                        "start": row.get("Start", ""),
                        "stop": row.get("Stop", ""),
                        "organism": row.get("Organism", ""),
                        "strain": row.get("Strain", ""),
                    }
                )
        # NC_ (reference) first
        candidates.sort(key=lambda c: (not c["nuc"].startswith("NC_"), c["nuc"]))
        candidates = candidates[:12]
        print(f"  {len(candidates)} complete-replicon candidates (capped)", file=sys.stderr)

        chosen = None
        for c in candidates:
            status, fields = adjacency_for(c["nuc"], c["start"], c["stop"])
            log.append(
                {
                    **c,
                    "gene": gene,
                    "uniprot": acc,
                    "adjacency": status,
                    "nearest": fields.get("nearest_component_dist", ""),
                    "note": fields.get("note", ""),
                }
            )
            print(
                f"    {c['nuc']:18s} {c['organism'][:34]:34s} {status:10s} "
                f"near={fields.get('nearest_component_dist', '')} [{fields.get('note', '')}]",
                file=sys.stderr,
            )
            if status == "CONFIRMED":
                chosen = (c, fields)
                break
        if chosen:
            new_placement[acc] = chosen
            print(f"  -> RE-PLACED {gene} into {chosen[0]['nuc']} ({chosen[0]['organism'][:30]})", file=sys.stderr)
        else:
            print(f"  -> NO complete genome confirms adjacency for {gene}; left as genuine exception", file=sys.stderr)

    # apply updates
    for r in rescued:
        acc = r["uniprot"]
        if acc in new_placement:
            c, _ = new_placement[acc]
            r["placement_assembly"] = c["assembly"]
            r["placement_nuc"] = c["nuc"]
            r["start"], r["stop"] = c["start"], c["stop"]
            r["placement_organism"] = c["organism"]
            r["species_match"] = finalize.match_level(r["corpus_organism"], c["organism"])
            r["pident"], r["qcov"] = "", ""  # now a complete-genome relocation; identity re-derivable from subject
    if new_placement:
        with open(RESCUED_TSV, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=rescued[0].keys(), delimiter="\t")
            w.writeheader()
            w.writerows(rescued)
        # the ceiling table (t1ss_ceiling.tsv) is regenerated by re-running script 15,
        # which reads the updated placements here -- script 15 stays the single source.

    with open(LOG_TSV, "w", newline="") as fh:
        if log:
            w = csv.DictWriter(fh, fieldnames=list(log[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(log)

    print(f"\nre-placed {len(new_placement)}/{len(TARGETS)} artifact effectors")
    for acc in TARGETS:
        if acc in new_placement:
            c, f = new_placement[acc]
            print(
                f"  {TARGETS[acc]:10s} -> {c['nuc']:16s} {c['organism'][:34]:34s} "
                f"adjacency CONFIRMED, nearest={f['nearest_component_dist']}"
            )
        else:
            print(f"  {TARGETS[acc]:10s} -> not relocated (genuine exception)")


if __name__ == "__main__":
    raise SystemExit(main())
