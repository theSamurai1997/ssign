#!/usr/bin/env python3
"""
24_actual_call.py  (Phase 2 tasks 6.3 + 6.4: per-effector ssign call for one run)

For each gold effector, decide what ssign DID with it in a given run:
  emitted_secreted : ssign emitted the effector's protein in the "# Secreted Proteins" list.
  not_emitted      : the protein is in ssign's input (found in results_raw) but NOT emitted.
  not_in_input     : ssign's run has no protein matching the effector (ORF not called, or the
                     effector's replicon was absent) -- should be rare given the staged inputs.
  no_run           : no ssign output for the effector's unit (genome not run / failed).

Matching ssign's protein to the gold effector (the "bridge", task 6.3, folded in here):
  1. locus_tag, drift-normalised (exact in GenBank --use-input-annotations mode, where ssign's
     locus_tag IS the RefSeq tag).
  2. else exact protein-sequence identity (FASTA/Bakta mode: Bakta assigns its own locus_tags
     but calls the same ORF -> identical translation). Effector protein sequence comes from the
     cached RefSeq GenBank by effector_locus.
  3. else unmatched -> not_in_input.

For a miss (not_emitted) the matched protein's per-tool signals are carried through so we can
see WHY ssign didn't emit it (not near a detected SS, low DLP prob, etc.) -- and those columns
double as the labelled feature set for the secretion-classifier model.

Inputs : data/phase1/ceiling_per_effector.tsv   (effector_locus + testable + reachable_n*)
         data/phase2/effector_unit_map.tsv        (uniprot -> unit_id)
         data/phase2/runs/<run_tag>/<unit>/results/<unit>_results{,_raw}.csv
         data/refseq_cache/<acc>.gb               (effector protein sequence, for the seq bridge)
Output : data/phase2/actual_per_effector.<run_tag>.tsv
Run:     .venv/bin/python scripts/24_actual_call.py --run-tag pilot_genbank_default
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
import bench_index as bi  # noqa: E402
import bench_runout as bo  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
CEILING = BENCH / "data" / "phase1" / "ceiling_per_effector.tsv"
MANIFEST = BENCH / "data" / "phase2" / "panel_manifest.tsv"
RUNS = BENCH / "data" / "phase2" / "runs"
CACHE = BENCH / "data" / "refseq_cache"

norm = bi.normalize  # lowercase + strip _/spaces; folds locus_tag underscore drift
load = bi.load_tsv


accbase = bi.accession_base


OUT_COLS = [
    "uniprot",
    "gene",
    "ss_type",
    "unit_id",
    "effector_locus",
    "ssign_locus",
    "match_method",
    "match_identity",
    "ssign_call",
    "testable",
    "reachable_n3",
    "reachable_n5",
    "reachable_n7",
    "ceiling_reason",
    *bo.SIGNAL_COLS,
]


def cache_file(acc):
    """Drift-tolerant cache lookup: version/prefix-stripped base -> the cached .gb stem."""
    base = accbase(acc)
    hits = sorted(
        (p for p in CACHE.glob("*.gb") if accbase(p.stem) == base),
        key=lambda p: len(p.stem),
        reverse=True,
    )
    return hits[0] if hits else None


class ProteinSeqs:
    """Lazy genome-accession -> {normalised locus_tag/old_locus_tag: protein sequence}."""

    def __init__(self):
        self._g = {}

    def get(self, genome_acc, locus_tag):
        key = accbase(genome_acc)
        if key not in self._g:
            self._g[key] = self._build(genome_acc)
        return self._g[key].get(norm(locus_tag))

    def _build(self, genome_acc):
        f = cache_file(genome_acc)
        out = {}
        if not f:
            return out
        for rec in SeqIO.parse(f, "genbank"):
            for feat in rec.features:
                if feat.type != "CDS":
                    continue
                tr = feat.qualifiers.get("translation", [""])[0].strip().upper()
                if not tr:
                    continue
                for q in ("locus_tag", "old_locus_tag"):
                    for tag in feat.qualifiers.get(q, []):
                        out[norm(tag)] = tr
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-tag", required=True, help="subdir under data/phase2/runs/ to score")
    args = ap.parse_args()
    run_root = RUNS / args.run_tag
    if not run_root.is_dir():
        sys.exit(f"no run dir: {run_root.relative_to(BENCH)} (rsync results back first)")

    # One row per gold effector (582). Do NOT key by uniprot: ~188 rows share an (often empty)
    # uniprot, and the run unit is fixed by the genome anyway, not the effector.
    ceiling = load(CEILING)
    # genome accession base -> run unit, from the staged panel manifest.
    unit_of_base = {}
    for m in load(MANIFEST):
        for acc in m["replicons"].split(","):
            unit_of_base[accbase(acc)] = m["unit_id"]
    seqs = ProteinSeqs()
    gene_order = bi.load_from_tsv()  # effector RefSeq coordinates, for the FASTA coordinate bridge

    runs = {}  # unit_id -> RunOutput | None (None = no/failed output)

    def run_for(unit):
        if unit not in runs:
            d = run_root / unit / "results"
            rc, raw = d / f"{unit}_results.csv", d / f"{unit}_results_raw.csv"
            runs[unit] = bo.RunOutput.load(rc, raw) if (rc.exists() and raw.exists()) else None
        return runs[unit]

    rows = []
    for c in ceiling:
        genome = c.get("refseq_genome", "").strip()
        unit = unit_of_base.get(accbase(genome), "") if genome and genome != "-" else ""
        out = {k: "" for k in OUT_COLS}
        out.update(
            uniprot=c.get("uniprot", ""),
            gene=c.get("gene", ""),
            ss_type=c["ss_type"],
            unit_id=unit,
            effector_locus=c.get("effector_locus", ""),
            testable=c.get("testable", ""),
            reachable_n3=c.get("reachable_n3", ""),
            reachable_n5=c.get("reachable_n5", ""),
            reachable_n7=c.get("reachable_n7", ""),
            ceiling_reason=c.get("reason", ""),
        )
        run = run_for(unit) if unit else None
        if run is None:
            out["ssign_call"] = bo.CALL_NO_RUN
            rows.append(out)
            continue

        eff_locus = c.get("effector_locus", "").strip()
        sl = mm = None
        ident = ""
        # 1. locus_tag match (normalised; exact in GenBank --use-input-annotations mode)
        if eff_locus:
            nl = norm(eff_locus)
            sl = next((lt for lt in run.by_locus if norm(lt) == nl), None)
            if sl:
                mm, ident = "locus_tag", "1.000"
        # 2. coordinate bridge (FASTA/Bakta mode: same ORF at the same contig+strand+3'-stop,
        #    even though Bakta renamed the locus_tag). The effector's RefSeq coords come from the
        #    gene-order index; results_raw carries the run's coords.
        if sl is None and eff_locus:
            ghit = gene_order.find(c.get("refseq_genome", ""), eff_locus)
            if ghit:
                rec_acc, _ordinal, cds = ghit
                scoord = run.find_by_coord(rec_acc, cds["start"], cds["end"], cds["strand"])
                if scoord:
                    sl, mm, ident = scoord, "coordinate", "coord"
        # 3. exact protein-sequence match, then 4. >=90% reciprocal identity+coverage
        #    (only fires if results_raw carries a `sequence` column; currently it does not).
        if sl is None and eff_locus:
            seq = seqs.get(c.get("refseq_genome", ""), eff_locus)
            if seq:
                sl = run.by_seq.get(seq)
                if sl:
                    mm, ident = "sequence", "1.000"
                else:
                    hit = run.fuzzy_find(seq, min_frac=0.90)
                    if hit:
                        sl, mm, ident = hit[0], "sequence_90", f"{hit[1]:.3f}"
        if sl is None:
            out.update(ssign_call=bo.CALL_NOT_IN_INPUT, match_method="unmatched")
            rows.append(out)
            continue

        out["ssign_locus"] = sl
        out["match_method"] = mm
        out["match_identity"] = ident
        out.update(run.by_locus[sl])  # sequence + signals (sequence overwritten below to drop bulk)
        out["sequence"] = ""  # don't carry the AA string into the per-effector table
        out["ssign_call"] = bo.CALL_EMITTED if sl in run.secreted else bo.CALL_NOT_EMITTED
        rows.append(out)

    out_path = BENCH / "data" / "phase2" / f"actual_per_effector.{args.run_tag}.tsv"
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    from collections import Counter

    calls = Counter(r["ssign_call"] for r in rows)
    matched = [r for r in rows if r["match_method"] in ("locus_tag", "coordinate", "sequence", "sequence_90")]
    mm = Counter(r["match_method"] for r in matched)
    print(f"wrote {out_path.relative_to(BENCH)}  ({len(rows)} effectors)")
    print(f"  ssign calls : {dict(calls)}")
    print(f"  bridge      : {dict(mm)}  (unmatched={sum(1 for r in rows if r['match_method'] == 'unmatched')})")
    # recall among testable
    test = [r for r in rows if r["testable"] == "yes"]
    emit = sum(r["ssign_call"] == bo.CALL_EMITTED for r in test)
    if test:
        print(f"  testable    : {len(test)}  emitted_secreted={emit} ({100 * emit / len(test):.0f}%)")


if __name__ == "__main__":
    raise SystemExit(main())
