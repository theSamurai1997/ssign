"""Unit tests for scripts/audit_disk_sizes.py.

The script's main job is "walk the manifest, du each path, present
three views". We don't test du itself (system tool); we pin the
view-formatting logic + the manifest-categorisation completeness
invariant (every manifest entry maps to a tool name).
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


class TestHumanise:
    def test_bytes(self, audit):
        assert audit._humanise(0) == "0 B"
        assert audit._humanise(512) == "512 B"

    def test_kilo_mega_giga(self, audit):
        assert audit._humanise(1024).startswith("1.0 KB")
        assert audit._humanise(1024 * 1024).startswith("1.0 MB")
        assert audit._humanise(1024**3).startswith("1.0 GB")

    def test_terabyte_ceiling(self, audit):
        assert audit._humanise(1024**4).endswith("TB")


class TestToolMapCoverage:
    """Every entry in the live manifest must map to a tool name in
    _TOOL_FOR_ENTRY. Catches the case where someone adds a new database
    or weights entry but forgets to update the audit categorisation."""

    def test_all_database_paths_categorised(self, audit):
        from ssign_lib.dependency_manifest import DATABASE_PATHS

        uncategorised = [e.name for e in DATABASE_PATHS if e.name not in audit._TOOL_FOR_ENTRY]
        assert not uncategorised, f"add to _TOOL_FOR_ENTRY: {uncategorised}"

    def test_all_model_weights_categorised(self, audit):
        from ssign_lib.dependency_manifest import MODEL_WEIGHTS

        uncategorised = [e.name for e in MODEL_WEIGHTS if e.name not in audit._TOOL_FOR_ENTRY]
        assert not uncategorised, f"add to _TOOL_FOR_ENTRY: {uncategorised}"


class TestPerTierRollup:
    """Tier rollup is cumulative — extended includes base, full
    includes both. Mirrors how a user actually plans an install."""

    def test_base_only_counted_in_all_three(self, audit, capsys):
        m = audit.Measurement("X", "DeepSecE", "base", "/p", 100, "")
        audit._print_per_tier([m], sys.stdout)
        out = capsys.readouterr().out
        # All three tiers should show 100 bytes.
        assert out.count("100 B") == 3

    def test_extended_only_counted_in_extended_and_full(self, audit, capsys):
        m = audit.Measurement("X", "Bakta", "extended", "/p", 100, "")
        audit._print_per_tier([m], sys.stdout)
        out = capsys.readouterr().out
        # 100 B appears in extended + full rows; base row is "0 B".
        assert "0 B" in out
        assert out.count("100 B") == 2

    def test_missing_entries_ignored(self, audit, capsys):
        m = audit.Measurement("X", "Bakta", "extended", None, 0, "not fetched")
        audit._print_per_tier([m], sys.stdout)
        out = capsys.readouterr().out
        # All tiers are 0 B; entry isn't fetched so nothing gets summed.
        assert out.count("0 B") == 3


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
