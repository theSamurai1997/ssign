#!/usr/bin/env python3
"""
bench_runout.py  (shared lib: read one ssign run's output for Phase 2 scoring)

ssign writes two CSVs per genome:
  <sid>_results.csv      multi-chunk; the first chunk "# Secreted Proteins" is the
                         post-filter list of proteins ssign EMITTED as secreted.
  <sid>_results_raw.csv  one row per protein in the genome (secreted or not), with the
                         per-tool signals and `contig/start/end/strand` coordinates.

Phase 2 matches a gold effector to ssign's protein by locus_tag (exact, after
drift-normalisation), which holds in GenBank `--use-input-annotations` mode where ssign keeps
the RefSeq tags. In FASTA mode Bakta re-annotates with its own locus_tags, so the fallback is
a **coordinate bridge**: Bakta calls the same ORF at the same genomic position as RefSeq, so
matching on (contig, strand, 3'-stop) recovers it. A protein-sequence identity bridge also
exists for when results_raw carries a `sequence` column (it currently does not, so coordinates
are the active FASTA path).

Signals kept (explain a miss AND double as features for the secretion-classifier model):
  nearby_ss_types, substrate_source, predicted_localization, dlp_extracellular_prob,
  dse_ss_type, signalp_prediction, signalp_probability, plm_effector_secreted,
  plm_effector_type, n_tools_agreeing, confidence_tier.
"""

from __future__ import annotations

import csv
import difflib
import io
from pathlib import Path

# ssign-call vocabulary (shared by scripts 24 + 25 so a typo can't silently break the join).
CALL_EMITTED = "emitted_secreted"
CALL_NOT_EMITTED = "not_emitted"
CALL_NOT_IN_INPUT = "not_in_input"
CALL_NO_RUN = "no_run"

SIGNAL_COLS = (
    "nearby_ss_types",
    "substrate_source",
    "predicted_localization",
    "dlp_extracellular_prob",
    "dse_ss_type",
    "signalp_prediction",
    "signalp_probability",
    "plm_effector_secreted",
    "plm_effector_type",
    "n_tools_agreeing",
    "confidence_tier",
)


def _read_secreted_chunk(path: Path):
    """Rows of the '# Secreted Proteins' chunk of a results.csv (up to the first blank line)."""
    with open(path, newline="") as fh:
        lines = fh.read().split("\n")
    # find the secreted-proteins header line (the comment), data starts on the next line
    start = next((i for i, ln in enumerate(lines) if ln.strip().lower().startswith("# secreted")), None)
    if start is None:
        return []
    end = next((i for i in range(start + 1, len(lines)) if lines[i].strip() == ""), len(lines))
    block = lines[start + 1 : end]
    if not block:
        return []
    # rejoin and parse as real CSV so a quoted field with an embedded newline survives
    return list(csv.DictReader(io.StringIO("\n".join(block))))


def _signals(row):
    return {c: (row.get(c, "") or "") for c in SIGNAL_COLS}


def _strand_norm(s) -> int:
    """'+'/'1' -> 1, '-'/'-1' -> -1, else 0."""
    s = str(s).strip()
    return 1 if s in ("+", "1") else -1 if s in ("-", "-1") else 0


def _contig_base(acc: str) -> str:
    """Strip a trailing .version so NC_002516.2 and NC_002516 compare equal."""
    acc = (acc or "").strip()
    head, _, tail = acc.rpartition(".")
    return head if (head and tail.isdigit()) else acc


def _three_prime(start: int, end: int, strand: int) -> int:
    """The stop-codon coordinate (most annotator-stable across Bakta vs RefSeq)."""
    return end if strand == 1 else start


class RunOutput:
    """One genome's ssign result: which proteins were emitted secreted, and every protein's
    sequence + signals for matching and miss-diagnosis."""

    def __init__(self, secreted, by_locus, by_seq, by_coord=None):
        self.secreted = secreted  # set of locus_tags emitted as secreted
        self.by_locus = by_locus  # locus_tag -> {sequence, **signals}
        self.by_seq = by_seq  # sequence -> locus_tag (first wins; sequence bridge, if a seq col exists)
        self.by_coord = by_coord or {}  # (contig_base, strand, 3'-stop) -> locus_tag (FASTA coord bridge)

    def find_by_coord(self, contig, start, end, strand, tol=3):
        """ssign locus whose 3'-stop sits within `tol` bp of the effector's, same contig+strand.
        Bakta may shift the start codon but rarely the stop, so the 3' end is the stable anchor."""
        cb, st = _contig_base(contig), _strand_norm(strand)
        if st == 0:
            return None
        tp = _three_prime(int(start), int(end), st)
        for d in range(-tol, tol + 1):
            lt = self.by_coord.get((cb, st, tp + d))
            if lt:
                return lt
        return None

    def fuzzy_find(self, seq, min_frac=0.90):
        """Best ssign protein whose IDENTICAL-residue count is >= min_frac of BOTH its length
        and the effector's. Requiring the reciprocal fraction enforces >=90% identity AND
        >=90% coverage at once, so a Bakta ORF with a shifted start/stop (a few residues
        trimmed or added) still matches its RefSeq effector. Returns (locus_tag, frac) or None.
        Length-prefiltered so this only aligns plausible candidates (fast: only runs for the
        few effectors that miss the exact locus + exact sequence tiers)."""
        if not seq:
            return None
        L = len(seq)
        best = None
        for lt, rec in self.by_locus.items():
            s = rec["sequence"]
            if not s or not (0.85 <= len(s) / L <= 1.18):
                continue
            matched = sum(b.size for b in difflib.SequenceMatcher(None, seq, s, autojunk=False).get_matching_blocks())
            frac = min(matched / L, matched / len(s))
            if frac >= min_frac and (best is None or frac > best[1]):
                best = (lt, frac)
        return best

    @classmethod
    def load(cls, results_csv: Path, results_raw_csv: Path):
        secreted_rows = _read_secreted_chunk(results_csv)
        secreted = {r["locus_tag"] for r in secreted_rows if r.get("locus_tag")}

        by_locus, by_seq, by_coord = {}, {}, {}
        with open(results_raw_csv, newline="") as fh:
            for r in csv.DictReader(fh):
                lt = r.get("locus_tag", "")
                if not lt:
                    continue
                seq = (r.get("sequence", "") or "").strip().upper()
                rec = {"sequence": seq, **_signals(r)}
                by_locus[lt] = rec
                if seq:
                    by_seq.setdefault(seq, lt)
                st = _strand_norm(r.get("strand", ""))
                start, end = r.get("start", ""), r.get("end", "")
                if st and str(start).strip().lstrip("-").isdigit() and str(end).strip().lstrip("-").isdigit():
                    key = (_contig_base(r.get("contig", "")), st, _three_prime(int(start), int(end), st))
                    by_coord.setdefault(key, lt)
        # secreted-chunk signals are richer for emitted proteins; fold them in
        for r in secreted_rows:
            lt = r.get("locus_tag", "")
            if lt and lt in by_locus:
                by_locus[lt].update({c: (r.get(c, "") or "") for c in SIGNAL_COLS if r.get(c)})
        return cls(secreted, by_locus, by_seq, by_coord)
