"""Tests for src/ssign_app/core/runner.py.

Surface tested here:

- PipelineConfig dataclass — defaults, mutable-default safety, factories
- StepResult dataclass — basic shape
- run_script() — subprocess wrapper: missing script, exit codes, timeout
- PipelineRunner.__init__ — callback default, semaphore wiring, state
- check_dependencies() — warning detection via shutil.which monkeypatch
- _elapsed_str() — pre-start vs post-start formatting
- _save_progress() + load_progress() — JSON round-trip
- _try_resume() — completed-step set computation, missing/stale work_dir handling

Deep orchestration (the run() main loop, parallel groups, per-step _step_*
methods) is out of unit-test scope; tests/integration/test_pipeline_fixture.py
exercises that on the real T1SS fixture.
"""

import json
import os
import shutil
import sys
import threading
import time
from pathlib import Path

import pytest

from ssign_app.core import runner

# runner.py is the orchestrator under src/ssign_app/core/. Import via the
# installed package — editable install puts src/ on sys.path.
from ssign_app.core.runner import (
    PipelineConfig,
    PipelineRunner,
    StepResult,
    run_script,
)

# ---------------------------------------------------------------------------
# PipelineConfig
# ---------------------------------------------------------------------------


class TestPipelineConfig:
    def test_defaults_are_sane(self):
        c = PipelineConfig()
        assert c.wholeness_threshold == 0.8
        assert c.conf_threshold == 0.8
        assert c.proximity_window == 3
        assert c.required_fraction_correct == 0.8
        # signalp_mode / deeplocpro_mode are auto-detected at __post_init__:
        # in a clean dev env neither binary is on PATH, so both fall back to
        # "remote" with a logger warning.
        assert c.deeplocpro_mode in ("local", "remote")
        assert c.signalp_mode in ("local", "remote")

    def test_excluded_systems_default(self):
        c = PipelineConfig()
        assert c.excluded_systems == ["Flagellum", "Tad", "T3SS"]

    def test_excluded_systems_isolated_per_instance(self):
        # Mutable defaults via field(default_factory=...) — each instance
        # must get its own list. Pin this so a refactor to plain `= [...]`
        # fails the test.
        a = PipelineConfig()
        b = PipelineConfig()
        a.excluded_systems.append("MUTATED")
        assert "MUTATED" not in b.excluded_systems

    def test_plm_effector_types_default(self):
        c = PipelineConfig()
        assert c.plm_effector_types == ["T1SE", "T2SE", "T3SE", "T4SE", "T6SE"]

    def test_plm_effector_types_isolated_per_instance(self):
        a = PipelineConfig()
        b = PipelineConfig()
        a.plm_effector_types.append("EXTRA")
        assert "EXTRA" not in b.plm_effector_types

    def test_whole_genome_flags_default_false(self):
        # All four prediction tools must default to neighborhood-only.
        # plme_whole_genome was missing until 2026-06-02 — PLM-E silently
        # ran on the full proteome (~34× the intended work on K-12).
        c = PipelineConfig()
        assert c.dlp_whole_genome is False
        assert c.dse_whole_genome is False
        assert c.sp_whole_genome is False
        assert c.plme_whole_genome is False

    def test_skip_flags_align_with_install_tier(self):
        # At tier='base' (small install, no annotation DBs), heavy DB-bound
        # tools are off and the prediction core is on.
        c = PipelineConfig(tier="base")
        assert c.skip_blastp is True
        assert c.skip_eggnog is True
        assert c.skip_hhsuite is True
        assert c.skip_interproscan is True
        assert c.skip_plmblast is True
        # base-tier predictors all run
        assert c.skip_deeplocpro is False
        assert c.skip_signalp is False
        assert c.skip_deepsece is False
        assert c.skip_plm_effector is False
        assert c.skip_protparam is False

    def test_skip_flags_extended_tier_enables_annotation_tools(self):
        # At tier='extended', EggNOG / IPS / pLM-BLAST come on because
        # the extended DB bundle ships them. HH-suite stays off — its
        # hhblits step needs UniRef30 which is full-tier only.
        # BLAST NR is also off (full-tier only).
        c = PipelineConfig(tier="extended")
        assert c.skip_eggnog is False
        assert c.skip_interproscan is False
        assert c.skip_plmblast is False
        assert c.skip_hhsuite is True
        assert c.skip_blastp is True

    def test_skip_flags_full_tier_enables_blast_and_hhsuite(self):
        c = PipelineConfig(tier="full")
        assert c.skip_blastp is False
        assert c.skip_hhsuite is False

    def test_explicit_skip_overrides_tier_default(self):
        # CLI --skip-eggnog at extended-tier should still skip.
        c = PipelineConfig(tier="extended", skip_eggnog=True)
        assert c.skip_eggnog is True
        # EggNOG is the override; pLM-BLAST default-on sticks.
        assert c.skip_plmblast is False

    def test_unknown_tier_raises(self):
        import pytest

        with pytest.raises(ValueError, match="Unknown tier"):
            PipelineConfig(tier="ludicrous")

    def test_cpu_per_genome_factory_yields_positive_int(self):
        c = PipelineConfig()
        assert isinstance(c.cpu_per_genome, int)
        assert c.cpu_per_genome >= 1


class TestDTUModeAutoDetect:
    """signalp_mode / deeplocpro_mode resolve from None → local / remote."""

    def test_explicit_local_passes_through(self):
        c = PipelineConfig(signalp_mode="local")
        assert c.signalp_mode == "local"

    def test_explicit_remote_passes_through(self):
        c = PipelineConfig(deeplocpro_mode="remote")
        assert c.deeplocpro_mode == "remote"

    def test_auto_falls_back_to_remote_when_no_binary(self, monkeypatch, tmp_path):
        # Empty PATH guarantees shutil.which returns None; HOME redirected so
        # the conda-env scan can't accidentally pick up the dev machine's envs.
        monkeypatch.setenv("PATH", "")
        monkeypatch.setenv("HOME", str(tmp_path))
        c = PipelineConfig(signalp_path="", deeplocpro_path="")
        assert c.signalp_mode == "remote"
        assert c.deeplocpro_mode == "remote"

    def test_auto_picks_local_when_conda_env_holds_binary(self, monkeypatch, tmp_path):
        # Stage ~/.conda/envs/signalp6/bin/signalp6 under a fake HOME.
        bin_dir = tmp_path / ".conda" / "envs" / "signalp6" / "bin"
        bin_dir.mkdir(parents=True)
        stub = bin_dir / "signalp6"
        stub.write_text("#!/bin/sh\n")
        stub.chmod(0o755)
        monkeypatch.setenv("PATH", "")
        monkeypatch.setenv("HOME", str(tmp_path))
        c = PipelineConfig(signalp_path="", deeplocpro_path="")
        assert c.signalp_mode == "local"
        assert c.signalp_path == str(bin_dir)
        # The other tool isn't in any conda env → still remote.
        assert c.deeplocpro_mode == "remote"

    def test_configured_path_wins_over_conda_env(self, monkeypatch, tmp_path):
        # Put a stub in a conda env AND in a user-configured dir; the
        # configured dir should be the one that ends up on signalp_path.
        conda_bin = tmp_path / ".conda" / "envs" / "signalp6" / "bin"
        conda_bin.mkdir(parents=True)
        (conda_bin / "signalp6").write_text("#!/bin/sh\n")
        (conda_bin / "signalp6").chmod(0o755)
        configured = tmp_path / "custom"
        configured.mkdir()
        (configured / "signalp6").write_text("#!/bin/sh\n")
        (configured / "signalp6").chmod(0o755)
        monkeypatch.setenv("PATH", "")
        monkeypatch.setenv("HOME", str(tmp_path))
        c = PipelineConfig(signalp_path=str(configured))
        assert c.signalp_mode == "local"
        assert c.signalp_path == str(configured)

    def test_auto_picks_local_when_path_holds_executable(self, monkeypatch, tmp_path):
        # Drop a stub executable into tmp_path and pass that dir as the path.
        stub = tmp_path / "signalp6"
        stub.write_text("#!/bin/sh\n")
        stub.chmod(0o755)
        monkeypatch.setenv("PATH", "")  # ensure PATH lookup also fails
        c = PipelineConfig(signalp_path=str(tmp_path))
        assert c.signalp_mode == "local"

    def test_auto_picks_local_when_binary_on_path(self, monkeypatch, tmp_path):
        stub = tmp_path / "deeplocpro"
        stub.write_text("#!/bin/sh\n")
        stub.chmod(0o755)
        monkeypatch.setenv("PATH", str(tmp_path))
        c = PipelineConfig()
        assert c.deeplocpro_mode == "local"

    def test_signalp_path_env_var_fallback_drives_local(self, monkeypatch, tmp_path):
        stub = tmp_path / "signalp6"
        stub.write_text("#!/bin/sh\n")
        stub.chmod(0o755)
        monkeypatch.setenv("PATH", "")
        monkeypatch.setenv("SSIGN_SIGNALP_PATH", str(tmp_path))
        c = PipelineConfig()
        assert c.signalp_mode == "local"
        assert c.signalp_path == str(tmp_path)


class TestPipelineConfigMarkerFallback:
    """When ~/.ssign/db_root is set, PipelineConfig auto-resolves DB paths.

    Mirrors what doctor does — same DatabasePath.resolve_path source of
    truth — so they never disagree about where Bakta / EggNOG / pLM-BLAST
    live. Per-DB SSIGN_* env vars still override the marker.
    """

    def _stage_dbs(self, tmp_path, monkeypatch):
        """Stage a fake fetch_databases.sh layout and point the marker at it."""
        root = tmp_path / "dbs"
        (root / "bakta" / "db-light").mkdir(parents=True)
        (root / "hhsuite" / "pfam" / "PfamA_v38_2").mkdir(parents=True)
        (root / "plm_blast" / "ECOD30").mkdir(parents=True)
        (root / "plm_effector_weights").mkdir(parents=True)
        (root / "bakta" / "db-light" / "version.json").touch()
        (root / "hhsuite" / "pfam" / "PfamA_v38_2" / "PfamA_v38_2_a3m.ffdata").touch()
        (root / "plm_blast" / "ECOD30" / "0.emb").touch()
        home = tmp_path / "fake-home"
        (home / ".ssign").mkdir(parents=True)
        (home / ".ssign" / "db_root").write_text(str(root))
        monkeypatch.setenv("HOME", str(home))
        return root

    def test_marker_resolves_bakta_to_inner_db_dir(self, tmp_path, monkeypatch):
        root = self._stage_dbs(tmp_path, monkeypatch)
        c = PipelineConfig(input_path="x.gbff")
        assert c.bakta_db == str(root / "bakta" / "db-light")

    def test_marker_resolves_pfam_to_versioned_subdir(self, tmp_path, monkeypatch):
        root = self._stage_dbs(tmp_path, monkeypatch)
        c = PipelineConfig(input_path="x.gbff")
        assert c.hhsuite_pfam_db == str(root / "hhsuite" / "pfam" / "PfamA_v38_2")

    def test_marker_resolves_ecod(self, tmp_path, monkeypatch):
        root = self._stage_dbs(tmp_path, monkeypatch)
        c = PipelineConfig(input_path="x.gbff")
        assert c.plmblast_db == str(root / "plm_blast" / "ECOD30")

    def test_marker_resolves_plm_effector_weights(self, tmp_path, monkeypatch):
        root = self._stage_dbs(tmp_path, monkeypatch)
        c = PipelineConfig(input_path="x.gbff")
        assert c.plm_effector_weights_dir == str(root / "plm_effector_weights")

    def test_env_var_wins_over_marker(self, tmp_path, monkeypatch):
        # Env var pointing at a valid alternate layout wins. Lay down a
        # second Bakta DB in a separate dir and point BAKTA_DB at it.
        root = self._stage_dbs(tmp_path, monkeypatch)
        alt = tmp_path / "alt-bakta"
        (alt / "db-light").mkdir(parents=True)
        (alt / "db-light" / "version.json").touch()
        monkeypatch.setenv("BAKTA_DB", str(alt))
        c = PipelineConfig(input_path="x.gbff")
        assert c.bakta_db == str(alt / "db-light")
        # Other DBs still resolve via marker
        assert c.plmblast_db == str(root / "plm_blast" / "ECOD30")

    def test_env_var_pointing_at_parent_dir_auto_descends(self, tmp_path, monkeypatch):
        # The key bug-fix: BAKTA_DB=<root>/bakta (the parent of db-light)
        # used to crash Bakta with "version.json not readable". Now the
        # runner descends into db-light/ via the sentinel `db*/version.json`.
        root = self._stage_dbs(tmp_path, monkeypatch)
        monkeypatch.setenv("BAKTA_DB", str(root / "bakta"))
        c = PipelineConfig(input_path="x.gbff")
        assert c.bakta_db == str(root / "bakta" / "db-light")

    def test_bogus_env_var_falls_through_to_marker(self, tmp_path, monkeypatch):
        # If the env var points at a non-existent dir, runner ignores it and
        # falls back to the marker layout — better than handing the broken
        # path to the wrapper and dying with a tool-specific error.
        root = self._stage_dbs(tmp_path, monkeypatch)
        monkeypatch.setenv("BAKTA_DB", "/does/not/exist")
        c = PipelineConfig(input_path="x.gbff")
        assert c.bakta_db == str(root / "bakta" / "db-light")


class TestPipelineConfigHHsuiteEnvFallback:
    """SSIGN_HHSUITE_* env vars fill in empty fields; explicit paths win.

    Documented in --hhsuite-*-db CLI help and docs/how-to/install.md.
    """

    @pytest.mark.parametrize(
        "env_name, attr",
        [
            ("SSIGN_HHSUITE_PFAM", "hhsuite_pfam_db"),
            ("SSIGN_HHSUITE_PDB70", "hhsuite_pdb70_db"),
            ("SSIGN_HHSUITE_UNICLUST", "hhsuite_uniclust_db"),
        ],
    )
    def test_empty_field_picks_up_env_var(self, monkeypatch, env_name, attr):
        monkeypatch.setenv(env_name, "/tmp/fake_db")
        c = PipelineConfig()
        assert getattr(c, attr) == "/tmp/fake_db"

    @pytest.mark.parametrize(
        "env_name, attr",
        [
            ("SSIGN_HHSUITE_PFAM", "hhsuite_pfam_db"),
            ("SSIGN_HHSUITE_PDB70", "hhsuite_pdb70_db"),
            ("SSIGN_HHSUITE_UNICLUST", "hhsuite_uniclust_db"),
        ],
    )
    def test_explicit_path_overrides_env_var(self, monkeypatch, env_name, attr):
        monkeypatch.setenv(env_name, "/tmp/from_env")
        c = PipelineConfig(**{attr: "/explicit/path"})
        assert getattr(c, attr) == "/explicit/path"

    def test_no_env_var_set_keeps_empty_default(self, monkeypatch):
        for env in ("SSIGN_HHSUITE_PFAM", "SSIGN_HHSUITE_PDB70", "SSIGN_HHSUITE_UNICLUST"):
            monkeypatch.delenv(env, raising=False)
        c = PipelineConfig()
        assert c.hhsuite_pfam_db == ""
        assert c.hhsuite_pdb70_db == ""
        assert c.hhsuite_uniclust_db == ""

    def test_env_var_fallback_logs_audit_trail(self, monkeypatch, caplog):
        # Stale env vars in a dev shell should be visible — log INFO when the
        # fallback fires so the user can spot they're picking up a leftover.
        monkeypatch.setenv("SSIGN_HHSUITE_PFAM", "/tmp/fake_pfam")
        with caplog.at_level("INFO", logger="ssign_app.core.runner"):
            PipelineConfig()
        assert any("hhsuite_pfam_db" in rec.message and "/tmp/fake_pfam" in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# StepResult
# ---------------------------------------------------------------------------


class TestStepResult:
    def test_default_output_files_is_empty_dict(self):
        r = StepResult(name="step", success=True, message="ok")
        assert r.output_files == {}

    def test_output_files_isolated_per_instance(self):
        a = StepResult(name="a", success=True, message="ok")
        b = StepResult(name="b", success=True, message="ok")
        a.output_files["k"] = "v"
        assert b.output_files == {}


# ---------------------------------------------------------------------------
# run_script — subprocess wrapper
# ---------------------------------------------------------------------------


class TestRunScript:
    def test_missing_script_returns_minus_one(self):
        rc, stdout, stderr = run_script("does_not_exist.py", [])
        assert rc == -1
        assert "not found" in stderr.lower()

    def test_real_exit_code_returned(self, tmp_dir, monkeypatch):
        # Drop a tiny Python script into BIN_DIR so run_script picks it up

        script_path = os.path.join(tmp_dir, "exit_code.py")
        with open(script_path, "w") as f:
            f.write("import sys; sys.exit(int(sys.argv[1]))\n")
        monkeypatch.setattr(runner, "BIN_DIR", Path(tmp_dir))
        rc, _, _ = run_script("exit_code.py", ["7"])
        assert rc == 7

    def test_stdout_captured(self, tmp_dir, monkeypatch):

        script_path = os.path.join(tmp_dir, "echo.py")
        with open(script_path, "w") as f:
            f.write("print('hello from script')\n")
        monkeypatch.setattr(runner, "BIN_DIR", Path(tmp_dir))
        rc, stdout, _ = run_script("echo.py", [])
        assert rc == 0
        assert "hello from script" in stdout

    def test_timeout_yields_minus_one(self, tmp_dir, monkeypatch):

        script_path = os.path.join(tmp_dir, "slow.py")
        with open(script_path, "w") as f:
            f.write("import time; time.sleep(10)\n")
        monkeypatch.setattr(runner, "BIN_DIR", Path(tmp_dir))
        rc, _, stderr = run_script("slow.py", [], timeout=1)
        assert rc == -1
        assert "timeout" in stderr.lower()


# ---------------------------------------------------------------------------
# PipelineRunner __init__ + small accessors
# ---------------------------------------------------------------------------


class TestRunnerInit:
    def test_default_progress_callback_is_noop(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="t")
        r = PipelineRunner(c)
        # Default callback must be safely invokable without args
        r.progress("step", 50, "msg")  # no exception

    def test_results_files_start_empty(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="t")
        r = PipelineRunner(c)
        assert r.results == []
        assert r.files == {}

    def test_api_semaphores_default_to_empty_dict(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="t")
        r = PipelineRunner(c)
        assert r.api_sem == {}

    def test_api_semaphores_passed_through(self, tmp_dir):

        c = PipelineConfig(outdir=tmp_dir, sample_id="t")
        sem = threading.Semaphore(2)
        r = PipelineRunner(c, api_semaphores={"dtu": sem})
        assert r.api_sem == {"dtu": sem}


class TestProgressMonotonic:
    """`progress(...)` must never go backwards. Parallel-group `as_completed`
    is non-deterministic, so a step that finishes after a higher-ordinal
    sibling shouldn't reset the displayed bar. Holds the highest pct seen."""

    def test_lower_pct_is_clamped_to_max_seen(self, tmp_dir):
        seen = []
        r = PipelineRunner(PipelineConfig(outdir=tmp_dir), progress_callback=lambda s, p, m: seen.append(p))
        r.progress("a", 30, "")
        r.progress("b", 10, "")  # lower; should be clamped
        r.progress("c", 50, "")
        assert seen == [30, 30, 50]

    def test_equal_pct_passes_through(self, tmp_dir):
        seen = []
        r = PipelineRunner(PipelineConfig(outdir=tmp_dir), progress_callback=lambda s, p, m: seen.append(p))
        r.progress("a", 40, "")
        r.progress("b", 40, "")
        assert seen == [40, 40]


class TestRunScriptExceptionLogging:
    """`run_script()` catches generic Exception and returns (-1, "", str(e)).
    Audit found the traceback was being thrown away — must be logged so the
    underlying cause is recoverable from the log output."""

    def test_unexpected_exception_is_logged_with_traceback(self, monkeypatch, caplog):
        def boom(*args, **kwargs):
            raise PermissionError("denied")

        monkeypatch.setattr(runner.subprocess, "run", boom)
        # Use any script name that resolves under BIN_DIR
        any_script = next(iter(runner.BIN_DIR.glob("*.py")), None)
        if any_script is None:
            pytest.skip("BIN_DIR has no scripts")
        with caplog.at_level("ERROR", logger="ssign_app.core.runner"):
            rc, _, stderr = runner.run_script(any_script.name, [])
        assert rc == -1
        assert "denied" in stderr
        assert any("raised unexpectedly" in rec.message for rec in caplog.records)


class TestSemaphoreAcquireTimeout:
    """`DTU_SEMAPHORE_TIMEOUT_S` is the upper bound on rate-limit-hold wait."""

    def test_timeout_constant_is_a_positive_int(self):
        assert isinstance(runner.DTU_SEMAPHORE_TIMEOUT_S, int)
        assert runner.DTU_SEMAPHORE_TIMEOUT_S > 0


class TestElapsedStr:
    def test_returns_empty_before_start(self, tmp_dir):
        r = PipelineRunner(PipelineConfig(outdir=tmp_dir))
        assert r._elapsed_str() == ""

    def test_returns_seconds_only_under_a_minute(self, tmp_dir):
        r = PipelineRunner(PipelineConfig(outdir=tmp_dir))
        r.start_time = time.monotonic() - 5
        assert r._elapsed_str() == "5s"

    def test_returns_minutes_and_seconds(self, tmp_dir):
        r = PipelineRunner(PipelineConfig(outdir=tmp_dir))
        r.start_time = time.monotonic() - 125  # 2m 5s
        assert r._elapsed_str() == "2m 5s"


# ---------------------------------------------------------------------------
# check_dependencies
# ---------------------------------------------------------------------------


class TestCheckDependencies:
    @pytest.mark.parametrize("tool", ["hmmsearch", "macsyfinder"])
    def test_warns_when_core_tool_missing(self, tmp_dir, monkeypatch, tool):
        c = PipelineConfig(outdir=tmp_dir, skip_deepsece=True)
        r = PipelineRunner(c)
        monkeypatch.setattr(shutil, "which", lambda name: None)
        warnings = r.check_dependencies()
        assert any(tool in w for w in warnings)

    def test_no_warnings_when_all_present(self, tmp_dir, monkeypatch):
        c = PipelineConfig(outdir=tmp_dir, skip_deepsece=True)
        r = PipelineRunner(c)
        monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")
        warnings = r.check_dependencies()
        assert warnings == []

    def test_skip_deepsece_suppresses_torch_check(self, tmp_dir, monkeypatch):
        # Even with torch absent, skip_deepsece=True must produce no torch warning
        c = PipelineConfig(outdir=tmp_dir, skip_deepsece=True)
        r = PipelineRunner(c)
        monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")
        # Force the torch import to fail by clobbering it in sys.modules
        monkeypatch.setitem(sys.modules, "torch", None)
        warnings = r.check_dependencies()
        # No torch-related warning should appear when DSE is skipped
        assert not any("torch" in w for w in warnings)


# ---------------------------------------------------------------------------
# _save_progress + load_progress + _try_resume
# ---------------------------------------------------------------------------


class TestProgressRoundTrip:
    def test_save_then_load_recovers_steps(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r = PipelineRunner(c)
        r.work_dir = tmp_dir
        r.results = [
            StepResult(name="extract_proteins", success=True, message="ok"),
            StepResult(name="macsyfinder", success=True, message="ok"),
            StepResult(name="deeplocpro", success=False, message="failed"),
        ]
        r.files = {"proteins": os.path.join(tmp_dir, "proteins.faa")}
        r._save_progress()

        loaded_results, loaded_files, loaded_work, _ = PipelineRunner.load_progress(
            tmp_dir,
            sample_id="my_sample",
        )
        assert [s.name for s in loaded_results] == [
            "extract_proteins",
            "macsyfinder",
            "deeplocpro",
        ]
        assert {s.name: s.success for s in loaded_results} == {
            "extract_proteins": True,
            "macsyfinder": True,
            "deeplocpro": False,
        }
        assert loaded_files["proteins"] == r.files["proteins"]
        assert loaded_work == tmp_dir

    def test_load_returns_none_when_no_progress_file(self, tmp_dir):
        results, files, work, cfg = PipelineRunner.load_progress(
            tmp_dir,
            sample_id="nope",
        )
        assert results is None
        assert files == {}
        assert work == ""
        assert cfg == {}

    def test_save_creates_hidden_ssign_dir(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r = PipelineRunner(c)
        r.work_dir = tmp_dir
        r._save_progress()
        progress_file = os.path.join(tmp_dir, ".ssign", "my_sample_progress.json")
        assert os.path.exists(progress_file)
        with open(progress_file) as f:
            data = json.load(f)
        # Includes sample_id, work_dir, steps, files, config
        assert data["sample_id"] == "my_sample"
        assert "config" in data


class TestTryResume:
    def test_no_prior_progress_returns_empty_set(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="fresh")
        r = PipelineRunner(c)
        assert r._try_resume() == set()

    def test_resumes_only_successful_steps(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        r1.results = [
            StepResult(name="extract_proteins", success=True, message="ok"),
            StepResult(name="macsyfinder", success=True, message="ok"),
            StepResult(name="deeplocpro", success=False, message="failed"),
        ]
        r1._save_progress()

        r2 = PipelineRunner(c)
        completed = r2._try_resume()
        assert "extract_proteins" in completed
        assert "macsyfinder" in completed
        assert "deeplocpro" not in completed

    def test_resume_restores_work_dir_and_files(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        # Use a real file path so the resume validator keeps it
        proteins = os.path.join(tmp_dir, "proteins.faa")
        with open(proteins, "w") as f:
            f.write(">x\nMKT\n")
        r1.results = [StepResult(name="extract_proteins", success=True, message="ok")]
        r1.files = {"proteins": proteins}
        r1._save_progress()

        r2 = PipelineRunner(c)
        r2._try_resume()
        assert r2.work_dir == tmp_dir
        assert r2.files["proteins"] == proteins

    def test_resume_skipped_if_work_dir_gone(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        # Reference a work_dir that doesn't exist
        r1.work_dir = os.path.join(tmp_dir, "deleted_work")
        os.makedirs(r1.work_dir)
        r1.results = [StepResult(name="extract_proteins", success=True, message="ok")]
        r1._save_progress()
        # Wipe the work_dir before resume
        os.rmdir(r1.work_dir)

        r2 = PipelineRunner(c)
        completed = r2._try_resume()
        assert completed == set()


class TestReportFiguresUpstreamGuards:
    """`_step_report` and `_step_figures` must skip cleanly when their upstream
    `integrated` CSV is missing — otherwise they invoke generate_report.py /
    generate_figures.py with an empty path and produce a cryptic subprocess
    error."""

    @pytest.mark.parametrize("step_name", ["report", "figures"])
    def test_skip_when_integrated_missing(self, tmp_dir, step_name):
        c = PipelineConfig(outdir=tmp_dir, sample_id="x")
        r = PipelineRunner(c)
        r.work_dir = tmp_dir
        r.files = {}
        result = r._step_report() if step_name == "report" else r._step_figures()
        assert result.success is False
        assert result.name == step_name
        assert "skipping" in result.message.lower()

    @pytest.mark.parametrize("step_name", ["report", "figures"])
    def test_skip_when_integrated_path_does_not_exist(self, tmp_dir, step_name):
        c = PipelineConfig(outdir=tmp_dir, sample_id="x")
        r = PipelineRunner(c)
        r.work_dir = tmp_dir
        r.files = {"integrated": os.path.join(tmp_dir, "never_written.csv")}
        result = r._step_report() if step_name == "report" else r._step_figures()
        assert result.success is False
        assert "skipping" in result.message.lower()

    def test_resume_aborted_if_persisted_output_missing(self, tmp_dir):
        # Manifest says step succeeded with file X; X was deleted between
        # runs. Resume must refuse — better fresh than corrupt.
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        proteins = os.path.join(tmp_dir, "proteins.faa")
        with open(proteins, "w") as f:
            f.write(">x\nMKT\n")
        r1.results = [StepResult(name="extract_proteins", success=True, message="ok")]
        r1.files = {"proteins": proteins}
        r1._save_progress()

        os.remove(proteins)

        r2 = PipelineRunner(c)
        assert r2._try_resume() == set()

    def test_resume_aborted_if_persisted_output_empty(self, tmp_dir):
        # Zero-byte TSV looks like exists==True but breaks downstream
        # parsers. Must not be silently accepted.
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        empty_path = os.path.join(tmp_dir, "ss_components.tsv")
        open(empty_path, "w").close()
        r1.results = [StepResult(name="macsyfinder", success=True, message="ok")]
        r1.files = {"ss_components": empty_path}
        r1._save_progress()

        r2 = PipelineRunner(c)
        assert r2._try_resume() == set()

    def test_resume_ignores_intentionally_blank_paths(self, tmp_dir):
        # Some steps record an empty-string path to mean "no output produced
        # this run" (e.g., skipped optional tool). Validator must skip those.
        c = PipelineConfig(outdir=tmp_dir, sample_id="my_sample")
        r1 = PipelineRunner(c)
        r1.work_dir = tmp_dir
        proteins = os.path.join(tmp_dir, "proteins.faa")
        with open(proteins, "w") as f:
            f.write(">x\nMKT\n")
        r1.results = [StepResult(name="extract_proteins", success=True, message="ok")]
        r1.files = {"proteins": proteins, "blastp_results": ""}
        r1._save_progress()

        r2 = PipelineRunner(c)
        completed = r2._try_resume()
        assert "extract_proteins" in completed
        assert r2.files["blastp_results"] == ""


# ---------------------------------------------------------------------------
# Constants integration
# ---------------------------------------------------------------------------


def test_hhsuite_min_prob_default_matches_constants():
    """PipelineConfig.hhsuite_min_prob must mirror the canonical value in
    ssign_lib.constants — single source of truth."""
    from ssign_app.scripts.ssign_lib.constants import HHSUITE_MIN_PROB

    assert PipelineConfig().hhsuite_min_prob == HHSUITE_MIN_PROB


class TestBuildRawCsv:
    """Raw CSV must carry every per-tool intermediate column, not just
    duplicate the substrate-filtered master CSV (CX3 K-12 run b2060a9
    surfaced this -- raw and non-raw were ~identical 38-column files)."""

    def _make_runner(self, tmp_path, files):
        cfg = PipelineConfig(sample_id="test")
        r = PipelineRunner(cfg)
        r.files = files
        return r

    def test_raw_csv_left_joins_every_per_tool_intermediate(self, tmp_path):
        import pandas as pd

        # Base: 4 proteins from gene_info
        gi = tmp_path / "gene_info.tsv"
        gi.write_text("locus_tag\tgbff_annotation\nP1\tprotein 1\nP2\tprotein 2\nP3\tprotein 3\nP4\tprotein 4\n")

        # Per-tool intermediates with non-overlapping per-tool columns
        pred = tmp_path / "predictions.tsv"
        pred.write_text("locus_tag\tdlp_extracellular_prob\tplm_effector_secreted\nP1\t0.95\tTrue\nP2\t0.10\tFalse\n")
        blastp = tmp_path / "blastp.csv"
        blastp.write_text("locus_tag,blastp_hit_description\nP1,hemolysin family\n")
        ips = tmp_path / "ips.tsv"
        ips.write_text("locus_tag\tinterpro_descriptions\nP3\tPF03797 autotransporter\n")

        r = self._make_runner(
            tmp_path,
            {
                "gene_info": str(gi),
                "predictions": str(pred),
                "blastp": str(blastp),
                "interproscan": str(ips),
            },
        )

        out = tmp_path / "raw.csv"
        r._build_raw_csv(out)

        df = pd.read_csv(out)
        # Every protein from gene_info is present
        assert sorted(df["locus_tag"].tolist()) == ["P1", "P2", "P3", "P4"]
        # Per-tool columns are preserved
        for col in (
            "gbff_annotation",
            "dlp_extracellular_prob",
            "plm_effector_secreted",
            "blastp_hit_description",
            "interpro_descriptions",
        ):
            assert col in df.columns, f"raw CSV missing {col}"
        # And carry per-protein data through
        row_p1 = df[df["locus_tag"] == "P1"].iloc[0]
        assert row_p1["blastp_hit_description"] == "hemolysin family"
        assert row_p1["dlp_extracellular_prob"] == 0.95
        # Tool didn't cover P4 -- column is NaN, not dropped
        row_p4 = df[df["locus_tag"] == "P4"].iloc[0]
        assert pd.isna(row_p4["blastp_hit_description"])

    def test_raw_csv_disambiguates_overlapping_columns(self, tmp_path):
        # If two intermediates both define a column with the same name
        # (e.g. both have "product"), we must keep both rather than have
        # the second silently overwrite the first via pandas merge.
        import pandas as pd

        gi = tmp_path / "gene_info.tsv"
        gi.write_text("locus_tag\tproduct\nP1\tgi_product_value\n")
        blastp = tmp_path / "blastp.csv"
        blastp.write_text("locus_tag,product\nP1,blastp_product_value\n")

        r = self._make_runner(tmp_path, {"gene_info": str(gi), "blastp": str(blastp)})
        out = tmp_path / "raw.csv"
        r._build_raw_csv(out)
        df = pd.read_csv(out)
        # Base keeps its name; collision is prefix-labelled.
        assert "product" in df.columns
        assert "blastp__product" in df.columns
        assert df.iloc[0]["product"] == "gi_product_value"
        assert df.iloc[0]["blastp__product"] == "blastp_product_value"


class TestStepSampleNullProteins:
    """The null-sampling step runs `sample_null_proteins.py` and concatenates
    its output FASTA with the neighborhood FASTA so DLP/DSE pick up both in
    a single tool invocation. SignalP and PLM-Effector continue to read
    `neighborhood_proteins`."""

    def _setup(self, tmp_path):
        # 15-protein two-contig fixture, mirrors conftest's two_contig_genes.
        from ssign_app.scripts.ssign_lib.fasta_io import write_fasta

        proteins = tmp_path / "proteins.faa"
        seqs = {f"GENE_{i:04d}": "M" + ("A" * 30) for i in range(10)}
        seqs.update({f"GENEB_{i:04d}": "M" + ("L" * 30) for i in range(5)})
        write_fasta(seqs, proteins)

        # 5 SS neighborhood proteins (a subset of the proteome -- the upstream
        # extract_neighborhood step would write exactly these).
        neighborhood = tmp_path / "neighborhood.faa"
        write_fasta({f"GENE_{i:04d}": seqs[f"GENE_{i:04d}"] for i in range(2, 7)}, neighborhood)

        gene_order = tmp_path / "gene_order.tsv"
        with open(gene_order, "w") as f:
            f.write("contig\tgene_index\tlocus_tag\n")
            for i in range(10):
                f.write(f"contig_A\t{i}\tGENE_{i:04d}\n")
            for i in range(5):
                f.write(f"contig_B\t{i}\tGENEB_{i:04d}\n")

        ss_components = tmp_path / "ss_components.tsv"
        with open(ss_components, "w") as f:
            f.write("locus_tag\tss_type\n")
            f.write("GENE_0004\tT2SS\n")
            f.write("GENE_0005\tT2SS\n")

        return proteins, neighborhood, gene_order, ss_components

    def test_writes_three_files_and_concat_contains_both(self, tmp_path):
        proteins, neighborhood, gene_order, ss_components = self._setup(tmp_path)

        config = PipelineConfig(
            outdir=str(tmp_path),
            sample_id="t",
            enrichment_stats=True,
            n_null_proteins=3,
            null_seed=42,
            proximity_window=1,  # tight window to leave a large null pool
        )
        r = PipelineRunner(config)
        r.work_dir = str(tmp_path)
        r.files = {
            "proteins": str(proteins),
            "gene_order": str(gene_order),
            "ss_components": str(ss_components),
            "neighborhood_proteins": str(neighborhood),
        }

        result = r._step_sample_null_proteins()
        assert result.success, result.message
        assert "null_proteins_fasta" in r.files
        assert "null_proteins_ids" in r.files
        assert "dlp_dse_input" in r.files

        # Concat must contain every neighborhood ID + every null ID, no duplicates
        from ssign_app.scripts.ssign_lib.fasta_io import read_fasta as _rf

        concat_ids = set(_rf(r.files["dlp_dse_input"]).keys())
        neigh_ids = set(_rf(str(neighborhood)).keys())
        null_ids = {line.strip() for line in open(r.files["null_proteins_ids"]) if line.strip()}
        assert neigh_ids.issubset(concat_ids)
        assert null_ids.issubset(concat_ids)
        assert null_ids.isdisjoint(neigh_ids)
        assert len(concat_ids) == len(neigh_ids) + len(null_ids)

    def test_missing_upstream_files_returns_failure(self, tmp_path):
        config = PipelineConfig(outdir=str(tmp_path), sample_id="t", enrichment_stats=True)
        r = PipelineRunner(config)
        r.work_dir = str(tmp_path)
        r.files = {}  # nothing upstream
        result = r._step_sample_null_proteins()
        assert not result.success
        assert "missing" in result.message.lower()


class TestWriteStepTimings:
    """`_write_step_timings` produces outdir/step_timings.csv after every
    run. Companion to resources.csv (sampler thread output) so the user
    can see per-tool wallclocks without grepping the log."""

    def test_writes_one_row_per_step_with_duration(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="x")
        r = PipelineRunner(c)
        r.results = [
            StepResult("macsyfinder", True, "ok", duration_s=82.0),
            StepResult("deeplocpro", True, "DLP complete", duration_s=111.4),
            StepResult("eggnog", False, "emapper.py missing", duration_s=0.5),
        ]
        r._write_step_timings()

        path = os.path.join(tmp_dir, "step_timings.csv")
        assert os.path.exists(path)
        with open(path) as fh:
            import csv as _csv

            rows = list(_csv.DictReader(fh))
        assert len(rows) == 3
        assert {r["step"] for r in rows} == {"macsyfinder", "deeplocpro", "eggnog"}
        eggnog = next(r for r in rows if r["step"] == "eggnog")
        assert eggnog["success"] == "0"
        assert eggnog["message"].startswith("emapper.py")
        macsy = next(r for r in rows if r["step"] == "macsyfinder")
        assert macsy["duration_s"] == "82.00"

    def test_skips_when_results_empty(self, tmp_dir):
        c = PipelineConfig(outdir=tmp_dir, sample_id="x")
        r = PipelineRunner(c)
        r.results = []
        r._write_step_timings()
        assert not os.path.exists(os.path.join(tmp_dir, "step_timings.csv"))

    def test_incremental_write_survives_midpipeline_kill(self, tmp_dir):
        # If the pipeline is killed (SIGKILL, OOM, PBS walltime) before the
        # last step finishes, per-step writes leave complete timing data
        # for every step that did complete instead of an empty file.
        c = PipelineConfig(outdir=tmp_dir, sample_id="x")
        r = PipelineRunner(c)
        path = os.path.join(tmp_dir, "step_timings.csv")
        completed = [
            StepResult("detect_format", True, "ok", duration_s=0.5),
            StepResult("extract_proteins", True, "ok", duration_s=12.0),
            StepResult("macsyfinder", True, "ok", duration_s=82.0),
        ]
        for step in completed:
            r._record_result(step)
        assert os.path.exists(path)
        with open(path) as fh:
            import csv as _csv

            rows = list(_csv.DictReader(fh))
        assert [row["step"] for row in rows] == ["detect_format", "extract_proteins", "macsyfinder"]

    def test_atomic_replace_leaves_no_tmp_files(self, tmp_dir):
        # PBS walltime kills the process with SIGKILL — any in-flight
        # write should not leave a partially-written step_timings.csv.
        # write-then-os.replace guarantees the file is either pre-write
        # state or fully-new state, never half-written. After a write,
        # the .tmp.<pid> sidecar must be gone.
        c = PipelineConfig(outdir=tmp_dir, sample_id="x")
        r = PipelineRunner(c)
        r.results = [StepResult("step1", True, "ok", duration_s=1.0)]
        r._write_step_timings()
        leftover = [p for p in os.listdir(tmp_dir) if p.startswith("step_timings.csv.tmp")]
        assert leftover == []


class TestResolveStepInputFasta:
    """Centralises FASTA selection for the four prediction steps so they
    cannot drift. DLP/DSE see the null-pool concat when enrichment-stats
    is on; SignalP/PLM-E never do. Whole-genome mode forces the full
    proteome regardless."""

    def _runner(self, tmp_dir, files):
        c = PipelineConfig(outdir=tmp_dir, sample_id="x")
        r = PipelineRunner(c)
        r.files = files
        return r

    def test_whole_genome_forces_full_proteome(self, tmp_dir):
        r = self._runner(tmp_dir, {"proteins": "/p", "neighborhood_proteins": "/n", "dlp_dse_input": "/c"})
        assert r._resolve_step_input_fasta(whole_genome=True) == "/p"
        assert r._resolve_step_input_fasta(whole_genome=True, include_null_concat=True) == "/p"

    def test_default_returns_neighborhood(self, tmp_dir):
        r = self._runner(tmp_dir, {"proteins": "/p", "neighborhood_proteins": "/n"})
        assert r._resolve_step_input_fasta(whole_genome=False) == "/n"

    def test_null_concat_preferred_for_dlp_dse(self, tmp_dir):
        r = self._runner(tmp_dir, {"proteins": "/p", "neighborhood_proteins": "/n", "dlp_dse_input": "/c"})
        assert r._resolve_step_input_fasta(whole_genome=False, include_null_concat=True) == "/c"

    def test_null_concat_ignored_for_signalp_plme(self, tmp_dir):
        # PLM-E and SignalP pass include_null_concat=False because their
        # models are too expensive to run on the null pool.
        r = self._runner(tmp_dir, {"proteins": "/p", "neighborhood_proteins": "/n", "dlp_dse_input": "/c"})
        assert r._resolve_step_input_fasta(whole_genome=False, include_null_concat=False) == "/n"

    def test_falls_back_to_full_proteome_when_no_neighborhood(self, tmp_dir):
        # MacSyFinder found no systems → no neighborhood FASTA staged.
        r = self._runner(tmp_dir, {"proteins": "/p"})
        assert r._resolve_step_input_fasta(whole_genome=False) == "/p"
        assert r._resolve_step_input_fasta(whole_genome=False, include_null_concat=True) == "/p"

    def test_raises_when_nothing_staged(self, tmp_dir):
        # No silent --input "" — surface "input-processing never ran"
        # before launching a subprocess that would crash on it anyway.
        r = self._runner(tmp_dir, {})
        with pytest.raises(RuntimeError, match="input-processing"):
            r._resolve_step_input_fasta(whole_genome=False)
        with pytest.raises(RuntimeError, match="input-processing"):
            r._resolve_step_input_fasta(whole_genome=True)


class TestStepPlmEffectorInput:
    """`_step_plm_effector` must read the SS neighborhood by default and
    only fall back to the full proteome when plme_whole_genome=True. The
    bug fixed in 2026-06-02 (commit X) had it always read the full
    proteome — 34× more work on K-12, 42m wallclock on an L40S GPU."""

    def _stub_run_script(self, captured):
        def _fake(script_name, args, **kwargs):
            captured.append((script_name, list(args)))
            return (0, "", "")

        return _fake

    def _make_runner(self, tmp_path, plme_whole_genome, files):
        weights = tmp_path / "weights"
        weights.mkdir()
        config = PipelineConfig(
            outdir=str(tmp_path),
            sample_id="t",
            plm_effector_weights_dir=str(weights),
            plme_whole_genome=plme_whole_genome,
        )
        r = PipelineRunner(config)
        r.work_dir = str(tmp_path)
        r.files = files
        return r

    def _input_arg(self, captured):
        # First captured call is run_plm_effector.py; "--input" is the
        # second positional in our args list.
        _, args = captured[0]
        i = args.index("--input")
        return args[i + 1]

    def test_default_reads_neighborhood(self, tmp_path, monkeypatch):
        neigh = tmp_path / "neighborhood.faa"
        neigh.write_text(">a\nMKT\n")
        proteins = tmp_path / "proteins.faa"
        proteins.write_text(">a\nMKT\n>b\nMKL\n")
        r = self._make_runner(
            tmp_path,
            plme_whole_genome=False,
            files={"proteins": str(proteins), "neighborhood_proteins": str(neigh)},
        )
        captured = []
        monkeypatch.setattr(runner, "run_script", self._stub_run_script(captured))
        result = r._step_plm_effector()
        assert result.success, result.message
        assert self._input_arg(captured) == str(neigh)

    def test_whole_genome_flag_reads_full_proteome(self, tmp_path, monkeypatch):
        neigh = tmp_path / "neighborhood.faa"
        neigh.write_text(">a\nMKT\n")
        proteins = tmp_path / "proteins.faa"
        proteins.write_text(">a\nMKT\n>b\nMKL\n")
        r = self._make_runner(
            tmp_path,
            plme_whole_genome=True,
            files={"proteins": str(proteins), "neighborhood_proteins": str(neigh)},
        )
        captured = []
        monkeypatch.setattr(runner, "run_script", self._stub_run_script(captured))
        result = r._step_plm_effector()
        assert result.success, result.message
        assert self._input_arg(captured) == str(proteins)

    def test_ignores_dlp_dse_input_concat(self, tmp_path, monkeypatch):
        # The enrichment-stats dual-fasta (neighborhood + null sample) is
        # for DLP/DSE only. PLM-E must NOT pick it up — running ensembles
        # over the null pool is too expensive for the marginal info gained.
        neigh = tmp_path / "neighborhood.faa"
        neigh.write_text(">a\nMKT\n")
        concat = tmp_path / "dlp_dse_input.faa"
        concat.write_text(">a\nMKT\n>null1\nMKQ\n")
        r = self._make_runner(
            tmp_path,
            plme_whole_genome=False,
            files={
                "proteins": str(tmp_path / "proteins.faa"),
                "neighborhood_proteins": str(neigh),
                "dlp_dse_input": str(concat),
            },
        )
        captured = []
        monkeypatch.setattr(runner, "run_script", self._stub_run_script(captured))
        result = r._step_plm_effector()
        assert result.success, result.message
        assert self._input_arg(captured) == str(neigh)


class TestCheckRequiredExecutables:
    """Pre-flight gate that hard-fails on missing emapper.py / plmblast.py
    for enabled steps. Added 2026-06-02 after a CX3 PBS job wasted 56
    minutes before crashing on the missing emapper.py binary."""

    def test_passes_when_steps_skipped(self, tmp_dir, monkeypatch):
        # No optional steps enabled → no required executables.
        monkeypatch.setattr(shutil, "which", lambda _: None)
        monkeypatch.delenv("SSIGN_PLMBLAST_SCRIPT", raising=False)
        c = PipelineConfig(outdir=tmp_dir, sample_id="x", skip_eggnog=True, skip_plmblast=True)
        r = PipelineRunner(c)
        assert r.check_required_executables() == []

    def test_flags_missing_emapper(self, tmp_dir, monkeypatch):
        monkeypatch.setattr(shutil, "which", lambda _: None)
        monkeypatch.delenv("SSIGN_PLMBLAST_SCRIPT", raising=False)
        c = PipelineConfig(outdir=tmp_dir, sample_id="x", skip_eggnog=False, skip_plmblast=True)
        r = PipelineRunner(c)
        errors = r.check_required_executables()
        assert len(errors) == 1
        assert "emapper.py" in errors[0]
        assert "pip install --no-deps eggnog-mapper" in errors[0]

    def test_flags_missing_plmblast(self, tmp_dir, monkeypatch):
        monkeypatch.setattr(shutil, "which", lambda _: None)
        monkeypatch.delenv("SSIGN_PLMBLAST_SCRIPT", raising=False)
        c = PipelineConfig(outdir=tmp_dir, sample_id="x", skip_eggnog=True, skip_plmblast=False)
        r = PipelineRunner(c)
        errors = r.check_required_executables()
        assert len(errors) == 1
        assert "plmblast.py" in errors[0]
        assert "SSIGN_PLMBLAST_SCRIPT" in errors[0]

    def test_plmblast_env_var_satisfies_check(self, tmp_dir, monkeypatch, tmp_path):
        # Existing file at SSIGN_PLMBLAST_SCRIPT is enough — even without PATH.
        monkeypatch.setattr(shutil, "which", lambda _: None)
        fake_script = tmp_path / "plmblast.py"
        fake_script.write_text("# stub")
        monkeypatch.setenv("SSIGN_PLMBLAST_SCRIPT", str(fake_script))
        c = PipelineConfig(outdir=tmp_dir, sample_id="x", skip_eggnog=True, skip_plmblast=False)
        r = PipelineRunner(c)
        assert r.check_required_executables() == []


class TestPoolEnrichmentStats:
    """Cross-genome pooling: sum M and k across genomes per (broad_type, tool),
    weighted-average p_bg by n_null, re-run binomial test, BH FDR on pooled."""

    def _write_per_genome_tsv(self, path, rows):
        from ssign_app.scripts.enrichment_testing import OUT_FIELDS

        with open(path, "w", newline="") as f:
            import csv as _csv

            writer = _csv.DictWriter(f, fieldnames=OUT_FIELDS, delimiter="\t")
            writer.writeheader()
            for r in rows:
                writer.writerow(r)

    def _row(self, **kwargs):
        # Defaults match the OUT_FIELDS schema so the per-genome TSV is valid.
        base = {
            "sample_id": "g",
            "scope_kind": "system",
            "scope_id": "sys_1",
            "ss_type": "T2SS",
            "tool": "DLP",
            "M": 5,
            "k": 2,
            "p_bg": 0.1,
            "fold_enrich": 4.0,
            "pvalue": 0.01,
            "qvalue": 0.01,
            "significant": True,
            "n_null": 100,
        }
        base.update(kwargs)
        return base

    def test_sums_M_and_k_across_genomes(self, tmp_path):
        from ssign_app.core.runner import pool_enrichment_stats

        a = tmp_path / "a_enrichment_stats.tsv"
        b = tmp_path / "b_enrichment_stats.tsv"
        # Each genome: 1 T2SS system, single tool DLP row.
        self._write_per_genome_tsv(a, [self._row(scope_kind="system", ss_type="T2SS", M=5, k=2, p_bg=0.1, n_null=100)])
        self._write_per_genome_tsv(b, [self._row(scope_kind="system", ss_type="T2SS", M=10, k=4, p_bg=0.2, n_null=100)])

        out = tmp_path / "pooled.tsv"
        n = pool_enrichment_stats([str(a), str(b)], str(out))
        assert n == 1  # one (broad_type, tool) row

        import csv as _csv

        rows = list(_csv.DictReader(open(out), delimiter="\t"))
        assert len(rows) == 1
        r = rows[0]
        assert r["scope_kind"] == "broad_type_pool"
        assert r["scope_id"] == "T2SS"
        assert int(r["M"]) == 15
        assert int(r["k"]) == 6
        # Weighted p_bg: (100*0.1 + 100*0.2) / 200 = 0.15
        assert abs(float(r["p_bg"]) - 0.15) < 1e-9

    def test_prefers_broad_type_row_over_per_system(self, tmp_path):
        from ssign_app.core.runner import pool_enrichment_stats

        # Genome has 2 T2SS systems → per-genome script emitted both per-system
        # rows AND the broad_type aggregate. Pool should use the aggregate.
        tsv = tmp_path / "g_enrichment_stats.tsv"
        rows_in = [
            self._row(scope_kind="system", scope_id="s1", ss_type="T2SS", M=4, k=1, p_bg=0.1, n_null=100),
            self._row(scope_kind="system", scope_id="s2", ss_type="T2SS", M=4, k=1, p_bg=0.1, n_null=100),
            self._row(scope_kind="broad_type", scope_id="T2SS", ss_type="T2SS", M=7, k=2, p_bg=0.1, n_null=100),
        ]
        self._write_per_genome_tsv(tsv, rows_in)

        out = tmp_path / "pooled.tsv"
        pool_enrichment_stats([str(tsv)], str(out))

        import csv as _csv

        pooled = list(_csv.DictReader(open(out), delimiter="\t"))
        # Only one pooled row (T2SS, DLP); M comes from the broad_type aggregate
        # (which dedupes overlapping neighborhoods), not the sum of per-system.
        assert len(pooled) == 1
        assert int(pooled[0]["M"]) == 7
        assert int(pooled[0]["k"]) == 2

    def test_collapses_subtypes_to_broad_type(self, tmp_path):
        from ssign_app.core.runner import pool_enrichment_stats

        # T5aSS + T5bSS rows pool under T5SS.
        a = tmp_path / "a_enrichment_stats.tsv"
        b = tmp_path / "b_enrichment_stats.tsv"
        self._write_per_genome_tsv(a, [self._row(scope_kind="system", ss_type="T5aSS", M=3, k=1, p_bg=0.1, n_null=100)])
        self._write_per_genome_tsv(b, [self._row(scope_kind="system", ss_type="T5bSS", M=4, k=2, p_bg=0.1, n_null=100)])

        out = tmp_path / "pooled.tsv"
        pool_enrichment_stats([str(a), str(b)], str(out))

        import csv as _csv

        rows = list(_csv.DictReader(open(out), delimiter="\t"))
        assert len(rows) == 1
        assert rows[0]["scope_id"] == "T5SS"
        assert int(rows[0]["M"]) == 7
        assert int(rows[0]["k"]) == 3

    def test_missing_input_file_skipped(self, tmp_path):
        from ssign_app.core.runner import pool_enrichment_stats

        present = tmp_path / "present.tsv"
        self._write_per_genome_tsv(present, [self._row()])
        missing = tmp_path / "does_not_exist.tsv"
        out = tmp_path / "pooled.tsv"
        n = pool_enrichment_stats([str(present), str(missing)], str(out))
        assert n == 1  # missing file is silently skipped, present one still pools

    def test_empty_inputs_writes_header_only(self, tmp_path):
        from ssign_app.core.runner import pool_enrichment_stats

        out = tmp_path / "pooled.tsv"
        n = pool_enrichment_stats([], str(out))
        assert n == 0
        # File exists with just the header row
        assert os.path.exists(out)
        with open(out) as f:
            assert len(f.readlines()) == 1
