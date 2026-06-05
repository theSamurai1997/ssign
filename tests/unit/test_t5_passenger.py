"""Tests for the passenger-FASTA extractor (#75a)."""

from __future__ import annotations

import pytest
from _helpers import read_tsv_rows, write_tsv

from ssign_app.scripts.ssign_lib.constants import LINKER_LENGTH, MIN_PASSENGER_FOR_ANNOTATION
from ssign_app.scripts.ssign_lib.fasta_io import read_fasta, write_fasta
from ssign_app.scripts.ssign_lib.t5_passenger import (
    FULL,
    PASSENGER,
    _passenger_slice,
    _routing_decision,
    build_passenger_substituted_fasta,
    load_t5_classifications,
)

T5_CLASSIFICATION_FIELDS = [
    "locus_tag",
    "sample_id",
    "ss_type",
    "domain_group",
    "t5_quality_flag",
    "passenger_length",
    "sp_end",
    "barrel_start",
]


def _make_classification(
    locus,
    ss_type="T5aSS",
    domain_group="Classical AT",
    t5_quality_flag="",
    passenger_length=300,
    sp_end=22,
    barrel_start=352,
):
    """Build one classifications-TSV row. Default = clean Classical AT.

    Default geometry: sp_end=22, barrel_start=352, LINKER_LENGTH=30 →
    passenger_length = 352 - 30 - 22 = 300, matching t5ss_handler.classify_t5a's
    formula (barrel_start - LINKER_LENGTH - sp_end).
    """
    return {
        "locus_tag": locus,
        "sample_id": "test",
        "ss_type": ss_type,
        "domain_group": domain_group,
        "t5_quality_flag": t5_quality_flag,
        "passenger_length": passenger_length,
        "sp_end": sp_end,
        "barrel_start": barrel_start,
    }


# ---------------------------------------------------------------------------
# _passenger_slice
# ---------------------------------------------------------------------------


def test_passenger_slice_basic():
    """sp_end=22, barrel_start=352 → passenger_length 352-30-22=300, slice [22:322].

    Matches t5ss_handler.classify_t5a's formula. In this convention,
    LINKER_LENGTH effectively counts 1 position inside the barrel envelope,
    so barrel_start should be set to the position where the slice should
    stop + LINKER_LENGTH (not the literal first barrel residue).
    """
    seq = "M" * 22 + "P" * 300 + "L" * 30 + "B" * 100
    passenger = _passenger_slice(seq, sp_end=22, barrel_start=352)
    assert len(passenger) == 300
    assert passenger == "P" * 300


def test_passenger_slice_boundary_immediately_after_sp():
    """If barrel_start == sp_end + LINKER_LENGTH + 1, passenger is 0 aa."""
    seq = "M" * 22 + "L" * 30 + "B" * 100
    # sp_end=22, barrel_start=53 → end=53-30=23, start=22 → 1 aa passenger
    out = _passenger_slice(seq, sp_end=22, barrel_start=53)
    assert len(out) == 1


def test_passenger_slice_returns_empty_on_invalid_geometry():
    """Negative or impossible coordinates → empty string (caller falls back to full)."""
    seq = "M" * 200
    assert _passenger_slice(seq, sp_end=-1, barrel_start=100) == ""
    assert _passenger_slice(seq, sp_end=100, barrel_start=50) == ""  # barrel before SP
    assert _passenger_slice(seq, sp_end=10, barrel_start=10000) == ""  # past seq end


# ---------------------------------------------------------------------------
# _routing_decision
# ---------------------------------------------------------------------------


def test_routing_full_when_no_classification():
    """Non-T5aSS substrates have no classifications-TSV entry → full protein."""
    assert _routing_decision("X", None, min_passenger_length=25) == FULL


def test_routing_full_when_not_t5ass():
    """T5bSS/T5cSS rows live in classifications.tsv but get full protein."""
    cls = {
        "ss_type": "T5bSS",
        "t5_quality_flag": "",
        "passenger_length": 0,
        "sp_end": 22,
        "barrel_start": 353,
    }
    assert _routing_decision("X", cls, min_passenger_length=25) == FULL


def test_routing_full_when_quality_flag_set():
    """Any non-empty t5_quality_flag triggers fallback to full protein."""
    for flag in ("barrel_only", "no_signalp", "no_sec_signal", "omp_porin_no_at", "unclassified"):
        cls = {
            "ss_type": "T5aSS",
            "t5_quality_flag": flag,
            "passenger_length": 300,
            "sp_end": 22,
            "barrel_start": 353,
        }
        assert _routing_decision("X", cls, min_passenger_length=25) == FULL, flag


def test_routing_full_when_passenger_too_short():
    """passenger_length below threshold → full protein (cpc-90 floor)."""
    cls = {
        "ss_type": "T5aSS",
        "t5_quality_flag": "",
        "passenger_length": 24,
        "sp_end": 22,
        "barrel_start": 77,
    }
    assert _routing_decision("X", cls, min_passenger_length=25) == FULL


def test_routing_passenger_when_clean():
    cls = {
        "ss_type": "T5aSS",
        "t5_quality_flag": "",
        "passenger_length": 300,
        "sp_end": 22,
        "barrel_start": 353,
    }
    assert _routing_decision("X", cls, min_passenger_length=25) == PASSENGER


def test_routing_full_when_geometry_missing():
    """sp_end or barrel_start absent (None after _parse_int on empty TSV cell) → full."""
    cls_no_sp = {
        "ss_type": "T5aSS",
        "t5_quality_flag": "",
        "passenger_length": 300,
        "sp_end": None,
        "barrel_start": 353,
    }
    assert _routing_decision("X", cls_no_sp, min_passenger_length=25) == FULL

    cls_no_barrel = dict(cls_no_sp, sp_end=22, barrel_start=None)
    assert _routing_decision("X", cls_no_barrel, min_passenger_length=25) == FULL


# ---------------------------------------------------------------------------
# load_t5_classifications
# ---------------------------------------------------------------------------


def test_load_classifications_reads_tsv(tmp_path):
    path = tmp_path / "classifications.tsv"
    write_tsv(
        path,
        T5_CLASSIFICATION_FIELDS,
        [
            _make_classification("KO_001"),
            _make_classification("KO_002", t5_quality_flag="barrel_only", passenger_length=0),
        ],
    )
    out = load_t5_classifications(path)
    assert set(out) == {"KO_001", "KO_002"}
    assert out["KO_001"]["sp_end"] == 22
    assert out["KO_001"]["barrel_start"] == 352
    assert out["KO_001"]["t5_quality_flag"] == ""
    assert out["KO_002"]["t5_quality_flag"] == "barrel_only"
    assert out["KO_002"]["passenger_length"] == 0


def test_load_classifications_missing_file_returns_empty(tmp_path):
    """Pipeline runs without any T5aSS substrates produce no classifications TSV."""
    assert load_t5_classifications(tmp_path / "nonexistent.tsv") == {}


def test_load_classifications_handles_empty_cells(tmp_path):
    """sp_end='' (SignalP didn't run) parses to None, not a crash."""
    path = tmp_path / "classifications.tsv"
    write_tsv(
        path,
        T5_CLASSIFICATION_FIELDS,
        [_make_classification("KO_001", sp_end="", barrel_start="", t5_quality_flag="no_signalp")],
    )
    out = load_t5_classifications(path)
    assert out["KO_001"]["sp_end"] is None
    assert out["KO_001"]["barrel_start"] is None


# ---------------------------------------------------------------------------
# build_passenger_substituted_fasta — end-to-end
# ---------------------------------------------------------------------------


@pytest.fixture
def passenger_test_dir(tmp_path):
    """Three-protein fixture: clean T5aSS, flagged T5aSS, non-T5SS."""
    proteins = {
        "AT_CLEAN": "M" * 22 + "P" * 300 + "L" * 30 + "B" * 100,  # 452 aa
        "AT_FLAGGED": "M" * 22 + "P" * 300 + "L" * 30 + "B" * 100,
        "T2SS_SUB": "M" * 250,
        "AT_SHORT": "M" * 22 + "P" * 24 + "L" * 30 + "B" * 100,  # passenger=24, below cutoff
    }
    fasta = tmp_path / "proteins.fasta"
    write_fasta(proteins, fasta)

    classifications = tmp_path / "t5_classifications.tsv"
    write_tsv(
        classifications,
        T5_CLASSIFICATION_FIELDS,
        [
            _make_classification("AT_CLEAN"),
            _make_classification("AT_FLAGGED", t5_quality_flag="barrel_only"),
            _make_classification("AT_SHORT", passenger_length=24, sp_end=22, barrel_start=77),
        ],
    )
    return tmp_path, fasta, classifications


def test_build_routes_clean_t5ass_to_passenger(passenger_test_dir):
    tmp, fasta, cls = passenger_test_dir
    out_fasta = tmp / "out.fasta"
    out_source = tmp / "source.tsv"
    routing = build_passenger_substituted_fasta(fasta, cls, out_fasta, out_source)

    sequences = read_fasta(out_fasta)
    assert routing["AT_CLEAN"] == PASSENGER
    assert len(sequences["AT_CLEAN"]) == 300
    assert sequences["AT_CLEAN"] == "P" * 300


def test_build_falls_back_to_full_on_flag(passenger_test_dir):
    tmp, fasta, cls = passenger_test_dir
    out_fasta = tmp / "out.fasta"
    out_source = tmp / "source.tsv"
    routing = build_passenger_substituted_fasta(fasta, cls, out_fasta, out_source)

    sequences = read_fasta(out_fasta)
    assert routing["AT_FLAGGED"] == FULL
    assert len(sequences["AT_FLAGGED"]) == 452


def test_build_falls_back_when_passenger_below_cutoff(passenger_test_dir):
    tmp, fasta, cls = passenger_test_dir
    out_fasta = tmp / "out.fasta"
    out_source = tmp / "source.tsv"
    routing = build_passenger_substituted_fasta(fasta, cls, out_fasta, out_source)

    sequences = read_fasta(out_fasta)
    assert routing["AT_SHORT"] == FULL
    assert len(sequences["AT_SHORT"]) == 176  # full protein length


def test_build_non_t5ass_passes_through_unchanged(passenger_test_dir):
    tmp, fasta, cls = passenger_test_dir
    out_fasta = tmp / "out.fasta"
    out_source = tmp / "source.tsv"
    routing = build_passenger_substituted_fasta(fasta, cls, out_fasta, out_source)

    sequences = read_fasta(out_fasta)
    assert "T2SS_SUB" not in routing  # only T5aSS rows recorded
    assert sequences["T2SS_SUB"] == "M" * 250


def test_build_writes_source_sidecar(passenger_test_dir):
    tmp, fasta, cls = passenger_test_dir
    out_fasta = tmp / "out.fasta"
    out_source = tmp / "source.tsv"
    build_passenger_substituted_fasta(fasta, cls, out_fasta, out_source)

    rows = read_tsv_rows(out_source)
    sources = {r["locus_tag"]: r["t5_annotation_source"] for r in rows}
    assert sources["AT_CLEAN"] == PASSENGER
    assert sources["AT_FLAGGED"] == FULL
    assert sources["AT_SHORT"] == FULL
    assert "T2SS_SUB" not in sources  # non-T5aSS excluded from sidecar


def test_build_respects_min_passenger_argument(passenger_test_dir):
    """Raising the threshold flips clean Classical AT to fallback."""
    tmp, fasta, cls = passenger_test_dir
    out_fasta = tmp / "out.fasta"
    out_source = tmp / "source.tsv"
    routing = build_passenger_substituted_fasta(fasta, cls, out_fasta, out_source, min_passenger_length=500)
    assert routing["AT_CLEAN"] == FULL  # passenger=300 < 500
    assert routing["AT_SHORT"] == FULL


def test_default_min_passenger_is_25():
    """Pin the documented default. If this changes, README + design_decisions need updating."""
    assert MIN_PASSENGER_FOR_ANNOTATION == 25


def test_passenger_slice_matches_t5ss_handler_formula():
    """
    The geometry formula in t5ss_handler.classify_t5a:
        passenger_length = (barrel_start - LINKER_LENGTH) - (sp_end + 1) + 1
    must match our slice length exactly. If LINKER_LENGTH ever changes, this
    test pins both modules to the same definition.
    """
    sp_end = 22
    barrel_start = 353
    expected_length = (barrel_start - LINKER_LENGTH) - (sp_end + 1) + 1
    seq = "M" * 500
    actual = _passenger_slice(seq, sp_end=sp_end, barrel_start=barrel_start)
    assert len(actual) == expected_length


def test_build_no_t5ass_substrates_writes_full_fasta_unchanged(tmp_path):
    """Pipeline with no T5aSS substrates: every entry passes through, no sidecar rows."""
    proteins = {"GENE1": "M" * 100, "GENE2": "M" * 200}
    fasta = tmp_path / "proteins.fasta"
    write_fasta(proteins, fasta)

    empty_cls = tmp_path / "empty_classifications.tsv"
    write_tsv(empty_cls, T5_CLASSIFICATION_FIELDS, [])

    out_fasta = tmp_path / "out.fasta"
    out_source = tmp_path / "source.tsv"
    routing = build_passenger_substituted_fasta(fasta, empty_cls, out_fasta, out_source)

    assert routing == {}
    out = read_fasta(out_fasta)
    assert out == proteins
    # Sidecar exists but has only the header row.
    assert read_tsv_rows(out_source) == []


def test_build_handles_missing_classifications_tsv(tmp_path):
    """If t5ss_handler step was skipped entirely, helper must not crash."""
    proteins = {"GENE1": "M" * 100}
    fasta = tmp_path / "proteins.fasta"
    write_fasta(proteins, fasta)

    out_fasta = tmp_path / "out.fasta"
    out_source = tmp_path / "source.tsv"
    routing = build_passenger_substituted_fasta(fasta, tmp_path / "does_not_exist.tsv", out_fasta, out_source)
    assert routing == {}
    assert read_fasta(out_fasta) == proteins
