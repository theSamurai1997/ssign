"""Localization-correctness gate for MacSyFinder-validated secretion systems.

For each validated SS, compare every machinery component's DeepLocPro
``predicted_localization`` against the literature-derived expected set
in ``ssign_lib/data/component_localizations.tsv``. Components whose
DLP confidence is below ``dlp_confidence_threshold`` are excluded from
the calculation. A system passes when at least
``required_fraction_correct`` of its confidently-localized components
land in an acceptable compartment. Systems below the cutoff are
treated as "not real" and their substrates are dropped downstream.

See task #58 in the task list and the head of
``component_localizations.tsv`` for the data schema.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

# Rule key: (ss_type, txsscan_gene_name) -> set of acceptable DLP localization labels.
# Empty set is the sentinel for "no defensible expected localization" (component
# is skipped, contributing to neither numerator nor denominator).
LocalizationRules = dict[tuple[str, str], frozenset[str]]


@dataclass(frozen=True)
class SystemVerdict:
    """Per-system scoring outcome.

    n_correct + n_incorrect == n_scored (the gate's denominator).
    Components excluded for low DLP confidence or missing rules are
    NOT counted in any of these.
    """

    sys_id: str
    ss_type: str
    n_scored: int
    n_correct: int
    fraction_correct: float
    passed: bool


def default_localization_table_path() -> Path:
    """Path to the shipped literature-derived localization table."""
    return Path(__file__).parent / "data" / "component_localizations.tsv"


def load_component_localizations(path: str | Path | None = None) -> LocalizationRules:
    """Parse the TSV into ``{(ss_type, gene_name): frozenset(acceptable_loc)}``.

    Comment lines (starting with ``#``) and blank lines are skipped.
    Multi-compartment cells use ``|`` as the separator. Rows with an
    empty ``acceptable_localizations`` value land in the dict with an
    empty frozenset — callers should treat that as "skip this component".
    """
    if path is None:
        path = default_localization_table_path()
    path = Path(path)

    rules: LocalizationRules = {}
    with open(path) as f:
        non_comment = (line for line in f if line.strip() and not line.startswith("#"))
        reader = csv.DictReader(non_comment, delimiter="\t")
        for row in reader:
            ss_type = (row.get("ss_type") or "").strip()
            gene = (row.get("txsscan_gene_name") or "").strip()
            if not ss_type or not gene:
                continue
            cell = (row.get("acceptable_localizations") or "").strip()
            accept = frozenset(p.strip() for p in cell.split("|") if p.strip())
            rules[(ss_type, gene)] = accept
    return rules


def evaluate_system(
    sys_id: str,
    ss_type: str,
    components: list[tuple[str, str]],
    predictions: dict[str, dict[str, str]],
    rules: LocalizationRules,
    dlp_confidence_threshold: float,
    required_fraction_correct: float,
) -> SystemVerdict:
    """Score one SS and return whether it passes the gate.

    ``components`` is a list of ``(locus_tag, gene_name)`` for the
    components of this system. ``predictions`` maps locus_tag to a
    cross-validate output row (must contain ``predicted_localization``
    and the five ``*_prob`` columns).

    Pass rule: ``n_scored == 0`` is **fail-open** (we have no evidence
    to drop the system, so keep it). Otherwise pass iff
    ``n_correct / n_scored >= required_fraction_correct``.
    """
    n_correct = 0
    n_scored = 0
    for locus_tag, gene_name in components:
        acceptable = rules.get((ss_type, gene_name))
        if not acceptable:
            continue  # no rule (or sentinel empty) -> skip component
        row = predictions.get(locus_tag)
        if not row:
            continue  # protein not in DLP output (e.g. failed extraction)
        max_prob = _max_dlp_prob(row)
        if max_prob < dlp_confidence_threshold:
            continue  # low-confidence DLP call: skip from both numerator and denominator
        predicted = (row.get("predicted_localization") or "").strip()
        n_scored += 1
        if predicted in acceptable:
            n_correct += 1

    if n_scored == 0:
        return SystemVerdict(sys_id, ss_type, 0, 0, 0.0, True)  # fail-open
    fraction = n_correct / n_scored
    return SystemVerdict(sys_id, ss_type, n_scored, n_correct, fraction, fraction >= required_fraction_correct)


def aggregate_failed_ss_types(verdicts: list[SystemVerdict]) -> set[str]:
    """Return the set of ss_types where every system of that type failed.

    Substrate-level filtering downstream operates at ss_type granularity
    (``nearby_ss_types`` carries types, not sys_ids), so a partial failure
    is permissive: as long as one system of a type passes, substrates
    tagged with that type stay. An ss_type only "fails" when none of its
    systems pass — at that point the whole class is suspect and the
    substrates near them get dropped.
    """
    by_type: dict[str, list[bool]] = {}
    for v in verdicts:
        by_type.setdefault(v.ss_type, []).append(v.passed)
    return {ss_type for ss_type, results in by_type.items() if not any(results)}


# Column names in the cross_validate predictions TSV. The extracellular
# one is `dlp_`-prefixed (historical, see cross_validate_predictions.py:190);
# the other four are the bare names from DLP. Hard-coding the prefix in one
# place rather than five.
_PROB_COLS = (
    "dlp_extracellular_prob",
    "periplasmic_prob",
    "outer_membrane_prob",
    "cytoplasmic_prob",
    "cytoplasmic_membrane_prob",
)


def _max_dlp_prob(row: dict[str, str]) -> float:
    """Highest probability across the five DLP classes in this row.

    Tolerant of missing columns (treated as 0) and non-numeric strings
    (the empty string in particular, which is what cross_validate
    writes for proteins outside the DLP output).
    """
    best = 0.0
    for col in _PROB_COLS:
        try:
            v = float(row.get(col, 0) or 0)
        except (TypeError, ValueError):
            continue
        if v > best:
            best = v
    return best
