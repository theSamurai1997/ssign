"""Tests for the K-12 run-analysis script.

Focuses on the parsing logic (elapsed-time regex, parallel-group handling,
PLM-E config extraction). The pandas-based results comparison is exercised
via a small synthetic-CSV test.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_SCRIPT_DIR = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(_SCRIPT_DIR))

import analyse_k12_runs as analyse  # noqa: E402


class TestParseElapsed:
    @pytest.mark.parametrize(
        "token,expected",
        [
            ("19s", 19),
            ("4m 8s", 4 * 60 + 8),
            ("1h 23m 45s", 3600 + 23 * 60 + 45),
            ("2m", 120),
            ("1h", 3600),
        ],
    )
    def test_known_forms(self, token, expected):
        assert analyse.parse_elapsed(token) == expected


class TestParseLog:
    def _write(self, tmp_path, body):
        p = tmp_path / "run.log"
        p.write_text(body)
        return p

    def test_sequential_steps_recorded(self, tmp_path):
        log = self._write(
            tmp_path,
            "ssign — running on /tmp/foo.gbff\n"
            "  [  5%] Detecting input format — Done: Format: genbank | 19s elapsed\n"
            "  [ 10%] Extracting proteins — Done: Bakta (re-annotated): extracted 4314 proteins | 4m 8s elapsed\n"
            "  [ 15%] Running MacSyFinder — Done: MacSyFinder v2 complete | 4m 37s elapsed\n",
        )
        steps, plme = analyse.parse_log(log)
        assert [s.label for s in steps] == [
            "Detecting input format",
            "Extracting proteins",
            "Running MacSyFinder",
        ]
        assert [s.elapsed_s for s in steps] == [19, 248, 277]
        assert all(s.parallel_start_s is None for s in steps)
        assert plme["batch_size"] is None  # no PLM-E lines present

    def test_parallel_group_marks_members(self, tmp_path):
        log = self._write(
            tmp_path,
            "  [ 25%] Extracting SS neighborhood — Done: 128 neighborhood proteins | 4m 40s elapsed\n"
            "  [ 30%] Running in parallel: Predicting localization (DeepLocPro), Predicting secretion type (DeepSecE), Predicting signal peptides (SignalP) — 3 tools running simultaneously | 4m 40s elapsed\n"
            "  [ 35%] Predicting secretion type (DeepSecE) — Done: DeepSecE complete | 5m 22s elapsed\n"
            "  [ 35%] Predicting localization (DeepLocPro) — Done: DeepLocPro complete | 6m 4s elapsed\n"
            "  [ 40%] Predicting signal peptides (SignalP) — Done: SignalP complete | 7m 11s elapsed\n",
        )
        steps, _ = analyse.parse_log(log)
        labels = [s.label for s in steps]
        assert "Predicting localization (DeepLocPro)" in labels
        dse = next(s for s in steps if s.label.endswith("(DeepSecE)"))
        assert dse.parallel_start_s == 4 * 60 + 40
        neighbor = next(s for s in steps if s.label == "Extracting SS neighborhood")
        assert neighbor.parallel_start_s is None

    def test_plme_config_and_per_type_lines(self, tmp_path):
        log = self._write(
            tmp_path,
            "[run_plm_effector.py] PLM-Effector: extracting features for T1SE,T2SE on cuda (batch_size=64, chunk_size=256, dtype=bf16, PLMs=ESM-1b,ProtT5)\n"
            "[run_plm_effector.py] PLM-Effector: T1SE — wrote 4314 predictions (12 passing threshold) to /tmp/T1SE.tsv\n"
            "[run_plm_effector.py] PLM-Effector: T2SE — wrote 4314 predictions (5 passing threshold) to /tmp/T2SE.tsv\n",
        )
        _, plme = analyse.parse_log(log)
        assert plme["batch_size"] == 64
        assert plme["chunk_size"] == 256
        assert plme["dtype"] == "bf16"
        assert plme["type_counts"]["T1SE"] == (4314, 12)
        assert plme["type_counts"]["T2SE"] == (4314, 5)


class TestComputeDurations:
    def test_sequential_only(self):
        steps = [
            analyse.StepRecord(label="A", pct=5, elapsed_s=10, message=""),
            analyse.StepRecord(label="B", pct=10, elapsed_s=30, message=""),
            analyse.StepRecord(label="C", pct=15, elapsed_s=60, message=""),
        ]
        assert analyse.compute_durations(steps) == {"A": 10, "B": 20, "C": 30}

    def test_parallel_group(self):
        steps = [
            analyse.StepRecord(label="Extracting SS neighborhood", pct=25, elapsed_s=280, message=""),
            analyse.StepRecord(
                label="Predicting secretion type (DeepSecE)", pct=35, elapsed_s=322, message="", parallel_start_s=280
            ),
            analyse.StepRecord(
                label="Predicting localization (DeepLocPro)", pct=35, elapsed_s=364, message="", parallel_start_s=280
            ),
            analyse.StepRecord(
                label="Predicting signal peptides (SignalP)", pct=40, elapsed_s=431, message="", parallel_start_s=280
            ),
            # Step right after the parallel block. prev_sequential_end was last
            # set by "Extracting SS neighborhood" (280), so this step's duration
            # is 451 - 280 = 171 — i.e. the parallel block's wallclock 151s plus
            # whatever sequential work followed, exactly what cumulative-elapsed
            # gives us.
            analyse.StepRecord(label="Predicting effectors (PLM-Effector, 5 types)", pct=45, elapsed_s=451, message=""),
        ]
        durs = analyse.compute_durations(steps)
        assert durs["Predicting secretion type (DeepSecE)"] == 42
        assert durs["Predicting localization (DeepLocPro)"] == 84
        assert durs["Predicting signal peptides (SignalP)"] == 151
        # The PLM-E step is sequential, but prev_sequential_end is still 280
        # (the pre-parallel anchor) because parallel members don't update it.
        assert durs["Predicting effectors (PLM-Effector, 5 types)"] == 171


class TestCompareResults:
    def test_identical_results(self, tmp_path):
        pd = pytest.importorskip("pandas")
        run_a = tmp_path / "a"
        run_b = tmp_path / "b"
        run_a.mkdir()
        run_b.mkdir()
        df = pd.DataFrame(
            {
                "locus_tag": ["L1", "L2", "L3"],
                "broad_annotation": ["Adhesin", "Toxin", "Adhesin"],
                "confidence_tier": ["High", "Medium", "Low"],
                "ss_type": ["T1SS", "T2SS", "T1SS"],
            }
        )
        df.to_csv(run_a / "ecoli_k12_results.csv", index=False)
        df.to_csv(run_b / "ecoli_k12_results.csv", index=False)
        out = analyse.compare_results(run_a, run_b)
        assert out["exists_a"] and out["exists_b"]
        assert out["shape_a"] == (3, 4)
        assert out["columns_match"] is True
        assert out["broad_annotation_match"] == 3
        assert out["confidence_tier_match"] == 3
        assert out["ss_type_match"] == 3

    def test_one_disagreement_flagged(self, tmp_path):
        pd = pytest.importorskip("pandas")
        run_a = tmp_path / "a"
        run_b = tmp_path / "b"
        run_a.mkdir()
        run_b.mkdir()
        df_a = pd.DataFrame({"locus_tag": ["L1", "L2"], "confidence_tier": ["High", "Low"]})
        df_b = pd.DataFrame({"locus_tag": ["L1", "L2"], "confidence_tier": ["High", "Medium"]})
        df_a.to_csv(run_a / "ecoli_k12_results.csv", index=False)
        df_b.to_csv(run_b / "ecoli_k12_results.csv", index=False)
        out = analyse.compare_results(run_a, run_b)
        assert out["confidence_tier_match"] == 1
        assert out["confidence_tier_total"] == 2

    def test_missing_file(self, tmp_path):
        out = analyse.compare_results(tmp_path / "a", tmp_path / "b")
        assert out["exists_a"] is False
        assert out["exists_b"] is False


class TestCalibrationRows:
    def test_emits_one_row_per_known_tool(self):
        # Steps in execution order: compute_durations diffs cumulative elapsed
        # against the previous sequential step, so the input ordering matters.
        steps = [
            analyse.StepRecord(label="Detecting input format", pct=5, elapsed_s=19, message=""),
            analyse.StepRecord(label="Extracting proteins", pct=10, elapsed_s=248, message=""),
            analyse.StepRecord(label="Running MacSyFinder", pct=15, elapsed_s=277, message=""),
        ]
        summary = analyse.RunSummary(
            machine="CX3-L40S",
            run_dir=Path("/tmp/runB"),
            log_path=None,
            steps=steps,
            final_elapsed_s=277,
        )
        rows = analyse.calibration_rows(summary, tier="extended", genome="ecoli_k12", n_proteins=4314)
        tools = [r["tool"] for r in rows]
        assert "bakta" in tools  # mapped from "Extracting proteins"
        assert "macsyfinder" in tools
        assert "_pipeline_total" in tools
        bakta = next(r for r in rows if r["tool"] == "bakta")
        assert bakta["wallclock_s"] == 248 - 19  # 229s
        assert bakta["machine"] == "CX3-L40S"
        assert bakta["input"]["n_proteins"] == 4314
        assert bakta["success"] is True

    def test_plme_row_includes_batch_and_dtype(self):
        steps = [
            analyse.StepRecord(label="Predicting effectors (PLM-Effector, 5 types)", pct=45, elapsed_s=900, message=""),
        ]
        summary = analyse.RunSummary(
            machine="CX3-L40S",
            run_dir=Path("/tmp/runB"),
            log_path=None,
            steps=steps,
            plme_batch_size=64,
            plme_dtype="bf16",
            plme_type_counts={"T1SE": (4314, 5), "T2SE": (4314, 3)},
        )
        rows = analyse.calibration_rows(summary)
        plme_rows = [r for r in rows if r["tool"] == "plm_effector"]
        assert len(plme_rows) == 1
        notes = plme_rows[0]["notes"]
        assert "batch_size=64" in notes
        assert "dtype=bf16" in notes
        assert "T1SE" in notes
