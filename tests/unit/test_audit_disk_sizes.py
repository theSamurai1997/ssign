"""Unit tests for scripts/audit_disk_sizes.py.

The script's main job is "walk the manifest, du each path, present
three views". We don't test du itself (system tool); we pin the
view-formatting logic + the manifest-field completeness invariant
(every manifest entry has a non-empty `tool` field, which the audit
relies on for its per-tool rollup).
"""

import os
import sys

import pytest

_REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(_REPO_ROOT, "scripts"))


@pytest.fixture
def audit():
    """Import the script as a module so we can call its helpers."""
    import importlib

    if "audit_disk_sizes" in sys.modules:
        del sys.modules["audit_disk_sizes"]
    return importlib.import_module("audit_disk_sizes")


class TestHumaniseReexport:
    """audit_disk_sizes pulls humanise_bytes from ssign_lib.resources;
    this test just confirms the shared helper behaves the way the audit
    output assumes."""

    def test_imported_helper_handles_common_sizes(self):
        from ssign_lib.resources import humanise_bytes

        assert humanise_bytes(0) == "0 B"
        assert humanise_bytes(1024).startswith("1.0 KB")
        assert humanise_bytes(1024 * 1024).startswith("1.0 MB")
        assert humanise_bytes(1024**3).startswith("1.0 GB")
        assert humanise_bytes(1024**4).endswith("TB")


class TestManifestToolFieldCoverage:
    """Every database + weights entry in the live manifest must carry a
    non-empty `tool` field. Catches the case where someone adds a new
    entry but forgets to categorise it for audit/reporting consumers."""

    def test_all_database_paths_have_tool(self):
        from ssign_lib.dependency_manifest import DATABASE_PATHS

        missing = [e.name for e in DATABASE_PATHS if not e.tool]
        assert not missing, f"set .tool on these DatabasePath entries: {missing}"

    def test_all_model_weights_have_tool(self):
        from ssign_lib.dependency_manifest import MODEL_WEIGHTS

        missing = [e.name for e in MODEL_WEIGHTS if not e.tool]
        assert not missing, f"set .tool on these ModelWeights entries: {missing}"


class TestPerTierRollup:
    """Tier rollup is cumulative — extended includes base, full includes
    both. Assertions parse the printed rows by tier label rather than
    counting tokens so they don't break on padding changes."""

    @staticmethod
    def _tier_row(out: str, tier: str) -> str:
        for line in out.splitlines():
            stripped = line.strip()
            if stripped.startswith(tier):
                return stripped
        raise AssertionError(f"no row for tier={tier!r} in:\n{out}")

    def test_base_only_counted_in_all_three(self, audit, capsys):
        m = audit.Measurement("X", "DeepSecE", "base", "/p", 100, "")
        audit._print_per_tier([m], sys.stdout)
        out = capsys.readouterr().out
        for tier in ("base", "extended", "full"):
            assert "100 B" in self._tier_row(out, tier)

    def test_extended_only_counted_in_extended_and_full(self, audit, capsys):
        m = audit.Measurement("X", "Bakta", "extended", "/p", 100, "")
        audit._print_per_tier([m], sys.stdout)
        out = capsys.readouterr().out
        assert "0 B" in self._tier_row(out, "base")
        assert "100 B" in self._tier_row(out, "extended")
        assert "100 B" in self._tier_row(out, "full")

    def test_missing_entries_ignored(self, audit, capsys):
        m = audit.Measurement("X", "Bakta", "extended", None, 0, "not fetched")
        audit._print_per_tier([m], sys.stdout)
        out = capsys.readouterr().out
        for tier in ("base", "extended", "full"):
            assert "0 B" in self._tier_row(out, tier)


class TestPerToolRollup:
    def test_sums_across_multiple_entries_same_tool(self, audit, capsys):
        ms = [
            audit.Measurement("HH-suite Pfam", "HH-suite", "full", "/p1", 1024, ""),
            audit.Measurement("HH-suite PDB70", "HH-suite", "full", "/p2", 2048, ""),
        ]
        audit._print_per_tool(ms, sys.stdout)
        out = capsys.readouterr().out
        # Combined 3 KB total
        assert "HH-suite" in out
        assert "3.0 KB" in out


class TestMarkdownEmit:
    def test_emits_header_and_one_row_per_tool_tier(self, audit, capsys):
        ms = [
            audit.Measurement("Bakta DB", "Bakta", "extended", "/p", 1024, ""),
            audit.Measurement("BLAST NR", "BLAST+", "full", "/p", 2048, ""),
        ]
        audit._emit_markdown(ms, sys.stdout)
        out = capsys.readouterr().out
        assert "| Tool | Tier | Size |" in out
        assert "| Bakta | extended |" in out
        assert "| BLAST+ | full |" in out


class TestDuPortableFlag:
    """`du -sk` (POSIX kilobyte summary) is what the script uses, NOT
    `du -sb` (GNU-only). Pinning this so a future "let's go back to -b
    for byte precision" regression doesn't silently break the macOS run."""

    def test_du_command_uses_sk_not_sb(self, audit, monkeypatch, tmp_path):
        # Stub subprocess.run, capture the argv.
        captured = {}

        class _Result:
            returncode = 0
            stdout = "1024\t/path\n"
            stderr = ""

        def fake_run(cmd, **kwargs):
            captured["cmd"] = list(cmd)
            return _Result()

        monkeypatch.setattr("subprocess.run", fake_run)
        size, note = audit._du_bytes(str(tmp_path))
        assert captured["cmd"][:2] == ["du", "-sk"]
        # 1024 KB -> 1 MB
        assert size == 1024 * 1024
        assert note == ""
