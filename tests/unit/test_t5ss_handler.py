"""Tests for t5ss_handler.py.

T5SS components are self-secreting — the component IS its own substrate.
Pin: only ss_type values starting with "T5" are picked up; non-T5
components are ignored; predictions are looked up by locus_tag and merged
in; the substrate output tags everything with tool="T5SS-self".

Geometric Pfam filter (applies to T5aSS only):
- PF03797 + passenger >= 100aa → Classical AT       → KEEP
- PF03797 + passenger 1-99aa   → Minimal passenger  → KEEP
- PF03797 + passenger == 0     → Barrel-only        → DROP
- PF13505 only, no PF03797     → OMP/Porin          → DROP
- Neither HMM hit              → Unclassified-AT    → KEEP (lenient)

Wrapper-level tests use placeholder protein sequences that don't hit either
HMM, so every T5aSS component falls through to "Unclassified-AT" and is kept.
The keep/drop branches are exercised via direct calls to ``classify_t5a``.
"""

import os

import pytest
from _helpers import (
    PREDICTIONS_FIELDS,
    SS_COMPONENT_FIELDS,
    make_prediction_row,
    make_ss_component_row,
    read_tsv_rows,
    run_script_main,
    write_tsv,
)
from t5ss_handler import _parse_sp_end, classify_t5a
from t5ss_handler import main as t5ss_main

from ssign_app.scripts.ssign_lib.fasta_io import write_fasta


def _write_proteins_fasta(tmp_dir, loci):
    path = os.path.join(tmp_dir, "proteins.faa")
    # 20-aa placeholder — too short to hit PF03797 (LENG 255) or PF13505 (LENG 174)
    write_fasta({locus: "MKTLLLTLLCAFSVAQAVDL" for locus in loci}, path)
    return path


def _run_t5ss(monkeypatch, tmp_dir, components, predictions, pfam_hits=None):
    """Wrapper-test driver. Skips the real pyhmmer scan (the geometric classifier
    is covered directly by test_classify_*). Pass ``pfam_hits`` to inject hits
    for a wrapper-level check that drops barrel-only/OMP entries from the TSV."""
    import t5ss_handler

    monkeypatch.setattr(t5ss_handler, "scan_bundled_pfams", lambda _fasta: pfam_hits or {})

    ss_path = write_tsv(os.path.join(tmp_dir, "ss_components.tsv"), SS_COMPONENT_FIELDS, components)
    pred_path = write_tsv(os.path.join(tmp_dir, "predictions.tsv"), PREDICTIONS_FIELDS, predictions)
    proteins_path = _write_proteins_fasta(tmp_dir, [c["locus_tag"] for c in components])
    out_substrates = os.path.join(tmp_dir, "t5ss_substrates.tsv")
    out_domains = os.path.join(tmp_dir, "t5ss_domains.tsv")
    run_script_main(
        monkeypatch,
        t5ss_main,
        [
            "t5ss_handler",
            "--ss-components",
            ss_path,
            "--predictions",
            pred_path,
            "--proteins",
            proteins_path,
            "--sample",
            "test_sample",
            "--out-substrates",
            out_substrates,
            "--out-domains",
            out_domains,
        ],
    )
    return read_tsv_rows(out_substrates), read_tsv_rows(out_domains)


# --- Direct tests of the geometric classifier --------------------------------


def test_classify_classical_at():
    """PF03797 at residue 1300, signal peptide ends at 22 → passenger ~1247 aa."""
    group, length = classify_t5a({"PF03797": (1300, 1500)}, sp_end=22)
    assert group == "Classical AT"
    assert length >= 1247


def test_classify_minimal_passenger():
    """PF03797 at residue 100, SP ends at 22 → passenger 23-69 = 47 aa."""
    group, length = classify_t5a({"PF03797": (100, 350)}, sp_end=22)
    assert group == "Minimal passenger"
    assert 1 <= length < 100


def test_classify_barrel_only():
    """PF03797 at residue 30, SP ends at 22 → passenger ~0 → barrel-only."""
    group, length = classify_t5a({"PF03797": (30, 280)}, sp_end=22)
    assert group == "Barrel-only"
    assert length == 0


def test_classify_omp_porin():
    """PF13505 hit but no PF03797 → OMP/Porin, dropped."""
    group, length = classify_t5a({"PF13505": (10, 200)}, sp_end=22)
    assert group == "OMP/Porin (no AT barrel)"
    assert length == 0


def test_classify_unclassified_when_no_hmm_hit():
    """No HMM hits → lenient default keeps the protein for downstream review."""
    group, length = classify_t5a({}, sp_end=22)
    assert group == "Unclassified-AT"
    assert length == 0


def test_classify_missing_signalp_uses_position_1():
    """A real AT with SignalP miss → passenger computed from position 1.

    The lenient default protects against false-negatives when SignalP fails:
    barrel at 1300 → passenger ~= 1300 - 30 - 1 - 1 = 1268 → Classical AT.
    """
    group, length = classify_t5a({"PF03797": (1300, 1500)}, sp_end=None)
    assert group == "Classical AT"
    assert length >= 1265


def test_passenger_length_at_threshold():
    """Boundary check: passenger of exactly MIN_PASSENGER_LENGTH (100) is Classical."""
    # passenger_length = (barrel_start - 30) - (sp_end + 1) + 1
    # = (152 - 30) - (22 + 1) + 1 = 100 → exactly at threshold
    group, length = classify_t5a({"PF03797": (152, 400)}, sp_end=22)
    assert group == "Classical AT"
    assert length == 100
    # One residue shorter → Minimal
    group2, length2 = classify_t5a({"PF03797": (151, 400)}, sp_end=22)
    assert group2 == "Minimal passenger"
    assert length2 == 99


# --- _parse_sp_end --------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("22", 22),
        ("22-23", 22),  # SignalP CS can be reported as range; take left
        ("22.0", 22),
        ("", None),
        ("not-a-number", None),
    ],
)
def test_parse_sp_end(raw, expected):
    assert _parse_sp_end(raw) == expected


# --- Wrapper-level tests (use placeholder FASTA, lenient "Unclassified-AT") ---


@pytest.mark.parametrize("ss_type", ["T5aSS", "T5bSS", "T5cSS"])
def test_every_t5_subtype_recognised(monkeypatch, tmp_dir, ss_type):
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", ss_type)],
        [make_prediction_row("AT_1", dlp_ext=0.9)],
    )
    assert {r["locus_tag"] for r in substrates} == {"AT_1"}
    assert substrates[0]["nearby_ss_types"] == ss_type


def test_non_t5_component_ignored(monkeypatch, tmp_dir):
    substrates, domains = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [
            make_ss_component_row("T1_1", "T1SS"),
            make_ss_component_row("T2_1", "T2SS"),
            make_ss_component_row("T5_1", "T5aSS"),
            make_ss_component_row("T6_1", "T6SSi"),
        ],
        [make_prediction_row(loc, dlp_ext=0.9) for loc in ["T1_1", "T2_1", "T5_1", "T6_1"]],
    )
    assert {r["locus_tag"] for r in substrates} == {"T5_1"}
    assert {r["locus_tag"] for r in domains} == {"T5_1"}


def test_substrate_tool_is_t5ss_self(monkeypatch, tmp_dir):
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [make_prediction_row("AT_1", dlp_ext=0.9)],
    )
    assert substrates[0]["tool"] == "T5SS-self"


def test_prediction_fields_merged_into_substrate(monkeypatch, tmp_dir):
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [
            make_prediction_row(
                "AT_1",
                dlp_ext=0.95,
                dse_type="T1SS",
                dse_prob=0.7,
                signalp_pred="OTHER",
                signalp_prob=0.05,
                product="autotransporter passenger",
            )
        ],
    )
    s = substrates[0]
    assert s["dlp_extracellular_prob"] == "0.95"
    assert s["dse_ss_type"] == "T1SS"
    assert s["dse_max_prob"] == "0.7"
    assert s["product"] == "autotransporter passenger"


def test_t5ss_kept_even_when_dlp_low(monkeypatch, tmp_dir):
    # T5aSS median OM probability is ~0.47; the wrapper does NOT apply DLP
    # thresholds — components pass through regardless of DLP score.
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [make_prediction_row("AT_1", dlp_ext=0.05)],
    )
    assert {r["locus_tag"] for r in substrates} == {"AT_1"}
    assert substrates[0]["dlp_extracellular_prob"] == "0.05"


def test_missing_prediction_yields_empty_fields(monkeypatch, tmp_dir):
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [],
    )
    assert {r["locus_tag"] for r in substrates} == {"AT_1"}
    s = substrates[0]
    assert s["dlp_extracellular_prob"] in ("0", "0.0")
    assert s["dse_ss_type"] == ""
    assert s["product"] == ""


def test_domain_output_records_classification(monkeypatch, tmp_dir):
    """No pfam_hits injected → T5aSS = Unclassified-AT (kept).
    T5bSS / T5cSS bypass the geometric filter → tagged as {ss_type}-component."""
    _, domains = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [
            make_ss_component_row("T5a_1", "T5aSS"),
            make_ss_component_row("T5b_1", "T5bSS"),
        ],
        [],
    )
    by_locus = {r["locus_tag"]: r for r in domains}
    assert by_locus["T5a_1"]["domain_group"] == "Unclassified-AT"
    assert by_locus["T5b_1"]["domain_group"] == "T5bSS-component"
    assert by_locus["T5a_1"]["ss_type"] == "T5aSS"


def test_wrapper_drops_barrel_only_entries(monkeypatch, tmp_dir):
    """Inject a PF03797 hit at residue 30 (passenger=0) → Barrel-only → drop.
    Inject one at residue 1300 (passenger>>100) → Classical AT → keep."""
    substrates, domains = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [
            make_ss_component_row("KEEP_1", "T5aSS"),
            make_ss_component_row("DROP_1", "T5aSS"),
        ],
        [make_prediction_row(loc, dlp_ext=0.9) for loc in ["KEEP_1", "DROP_1"]],
        pfam_hits={
            "KEEP_1": {"PF03797": (1300, 1555)},
            "DROP_1": {"PF03797": (30, 285)},
        },
    )
    assert {r["locus_tag"] for r in substrates} == {"KEEP_1"}
    by_locus = {r["locus_tag"]: r for r in domains}
    assert by_locus["KEEP_1"]["domain_group"] == "Classical AT"
    assert by_locus["DROP_1"]["domain_group"] == "Barrel-only"


def test_no_t5_components_writes_empty_outputs(monkeypatch, tmp_dir):
    substrates, domains = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("T2_1", "T2SS")],
        [make_prediction_row("T2_1", dlp_ext=0.9)],
    )
    assert substrates == []
    assert domains == []


def test_invalid_dlp_prob_falls_back_to_zero(monkeypatch, tmp_dir):
    bad_pred = make_prediction_row("AT_1")
    bad_pred["dlp_extracellular_prob"] = "not-a-number"
    substrates, _ = _run_t5ss(
        monkeypatch,
        tmp_dir,
        [make_ss_component_row("AT_1", "T5aSS")],
        [bad_pred],
    )
    assert substrates[0]["dlp_extracellular_prob"] in ("0", "0.0")
