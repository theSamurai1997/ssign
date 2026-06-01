"""Tests for ssign_lib.subprocess_diag.dump_failure_log."""

import logging
import subprocess
import types

from ssign_lib.subprocess_diag import dump_failure_log


def _result(returncode=1, stdout="", stderr=""):
    return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr=stderr)


class TestDumpFailureLog:
    def test_writes_log_with_both_streams(self, tmp_path):
        r = _result(returncode=1, stdout="ouch", stderr="real cause here")
        err = dump_failure_log("Bakta", r, ["bakta", "--db", "x"], str(tmp_path))
        assert isinstance(err, RuntimeError)
        log = tmp_path / "bakta_failure.log"
        assert log.is_file()
        body = log.read_text()
        assert "Bakta exit code: 1" in body
        assert "bakta --db x" in body
        assert "ouch" in body
        assert "real cause here" in body

    def test_empty_streams_emit_placeholder(self, tmp_path):
        r = _result(returncode=2, stdout="", stderr="")
        dump_failure_log("X", r, ["x"], str(tmp_path))
        body = (tmp_path / "x_failure.log").read_text()
        # Both stream sections must still appear, with a sentinel so the user
        # knows the tool wrote nothing rather than the log being truncated.
        assert body.count("(empty)") == 2

    def test_tool_name_with_spaces_dashes_safe_in_filename(self, tmp_path):
        r = _result(returncode=1, stdout="", stderr="")
        dump_failure_log("pLM-BLAST embedding", r, ["x"], str(tmp_path))
        # spaces / dashes → underscores; case-folded; everything else preserved.
        assert (tmp_path / "plm_blast_embedding_failure.log").is_file()

    def test_returned_runtime_error_embeds_log_path(self, tmp_path):
        r = _result(returncode=42, stdout="", stderr="")
        err = dump_failure_log("InterProScan", r, ["ips"], str(tmp_path))
        msg = str(err)
        assert "exit code 42" in msg
        assert str(tmp_path / "interproscan_failure.log") in msg

    def test_write_failure_does_not_mask_tool_failure(self, tmp_path, caplog):
        # Pointing output_dir at a non-existent subdirectory makes the open()
        # raise OSError. The caller still gets a RuntimeError back so the
        # tool failure surfaces, with a warning in the log about the disk issue.
        bad_dir = str(tmp_path / "does_not_exist")
        r = _result(returncode=1, stdout="stack trace", stderr="")
        with caplog.at_level(logging.WARNING):
            err = dump_failure_log("Bakta", r, ["bakta"], bad_dir)
        assert isinstance(err, RuntimeError)
        assert "exit code 1" in str(err)
        assert any("Could not write" in r.message for r in caplog.records)

    def test_accepts_real_completed_process(self, tmp_path):
        # The helper should accept an actual subprocess.CompletedProcess
        # (not just SimpleNamespace stubs) — this catches attribute mismatches
        # that pure-stub tests would miss.
        real = subprocess.run(["false"], capture_output=True, text=True)
        err = dump_failure_log("false", real, ["false"], str(tmp_path))
        assert isinstance(err, RuntimeError)
        assert "false exit code 1" in str(err)
        assert (tmp_path / "false_failure.log").is_file()
