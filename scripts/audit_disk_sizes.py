#!/usr/bin/env python3
"""Measure on-disk size of every ssign database + model bundle, grouped
by tool and tier.

Reads the dependency manifest (single source of truth for "what ssign
expects on disk") and walks each path with `du`. Skips entries that
aren't actually fetched on this machine. Produces three views:

  1. Per-tool table (one row per database/weights bundle, with size)
  2. Per-tier rollup (cumulative install footprint at base, extended, full)
  3. Markdown snippet you can paste into docs/how-to/install.md

Usage:
  python scripts/audit_disk_sizes.py
  python scripts/audit_disk_sizes.py --data-root ~/.ssign
  python scripts/audit_disk_sizes.py --markdown > /tmp/sizes.md

The script does not require ssign to be on PATH; it imports the
manifest from src/ssign_app/scripts/ssign_lib/ directly so you can run
it from a fresh clone on CX3 before activating any venv.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from dataclasses import dataclass

# Add src/ to path so we can import the manifest without installing ssign.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "src", "ssign_app", "scripts"))

from ssign_lib.dependency_manifest import (  # noqa: E402
    DATABASE_PATHS,
    MODEL_WEIGHTS,
    DatabasePath,
    ModelWeights,
    Tier,
)

# Map manifest entry names to the tool the user thinks of them as
# belonging to. Single source of truth for the per-tool rollup; if a
# new manifest entry appears that isn't in this map, the audit will
# warn rather than silently miscategorise.
_TOOL_FOR_ENTRY: dict[str, str] = {
    "Bakta DB": "Bakta",
    "EggNOG DB": "EggNOG-mapper",
    "InterProScan DB": "InterProScan",
    "HH-suite Pfam": "HH-suite",
    "HH-suite PDB70": "HH-suite",
    "HH-suite UniRef30": "HH-suite",
    "pLM-BLAST ECOD70": "pLM-BLAST",
    "BLAST NR": "BLAST+",
    "DeepSecE checkpoint": "DeepSecE",
    "PLM-Effector ensemble weights": "PLM-Effector",
}


@dataclass
class Measurement:
    name: str  # manifest entry name
    tool: str  # canonical tool name (Bakta, EggNOG-mapper, ...)
    tier: Tier
    resolved_path: str | None
    bytes_used: int  # 0 when not fetched / unreachable
    note: str = ""  # populated when the path resolved but du failed, etc.


def _du_bytes(path: str) -> tuple[int, str]:
    """Return (bytes, note). bytes==0 with non-empty note on failure."""
    try:
        out = subprocess.run(
            ["du", "-sb", path],
            capture_output=True,
            text=True,
            timeout=300,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        return (0, f"du failed: {e.stderr.strip()[:120]}")
    except subprocess.TimeoutExpired:
        return (0, "du timed out (>5 min)")
    except FileNotFoundError:
        return (0, "du not on PATH")
    first = out.stdout.split()
    if not first or not first[0].isdigit():
        return (0, f"unexpected du output: {out.stdout[:120]}")
    return (int(first[0]), "")


def _humanise(n: int) -> str:
    """Pretty-print bytes as B / KB / MB / GB. Uses 1024-base."""
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024
    return f"{n:.1f} TB"


def _measure_database(entry: DatabasePath, db_root: str) -> Measurement:
    tool = _TOOL_FOR_ENTRY.get(entry.name, "OTHER")
    path = entry.resolve_path(db_root)
    if not path or not os.path.isdir(path):
        return Measurement(entry.name, tool, entry.tier, None, 0, "not fetched")
    size, note = _du_bytes(path)
    return Measurement(entry.name, tool, entry.tier, path, size, note)


def _measure_weights(entry: ModelWeights, data_root: str, db_root: str) -> Measurement:
    tool = _TOOL_FOR_ENTRY.get(entry.name, "OTHER")
    root = db_root if entry.under_db_root else data_root
    path = os.path.join(root, entry.default_subpath)
    if not os.path.exists(path):
        return Measurement(entry.name, tool, entry.tier, None, 0, "not fetched")
    if os.path.isdir(path):
        size, note = _du_bytes(path)
    else:
        # Individual file (e.g., DeepSecE checkpoint)
        try:
            size = os.path.getsize(path)
            note = ""
        except OSError as e:
            size, note = 0, str(e)
    return Measurement(entry.name, tool, entry.tier, path, size, note)


def _resolve_db_root(data_root: str) -> str:
    """Read the db_root marker fetch_databases.sh writes, else fall
    back to <data_root>/databases. Mirrors doctor.resolve_db_root."""
    marker = os.path.join(data_root, "db_root")
    if os.path.isfile(marker):
        try:
            with open(marker) as fh:
                value = fh.read().strip()
                if value:
                    return value
        except OSError:
            pass
    return os.path.join(data_root, "databases")


def _print_per_entry(measurements: list[Measurement], stream) -> None:
    print("\nPer-entry  (path resolution + actual on-disk size)", file=stream)
    print("─" * 60, file=stream)
    print(f"  {'Entry':<32} {'Tier':<10} {'Size':>10}  Path", file=stream)
    for m in measurements:
        size_str = _humanise(m.bytes_used) if m.resolved_path else "(missing)"
        path_str = m.resolved_path or m.note
        print(f"  {m.name:<32} {m.tier:<10} {size_str:>10}  {path_str}", file=stream)


def _print_per_tool(measurements: list[Measurement], stream) -> None:
    by_tool: dict[str, int] = {}
    for m in measurements:
        if m.bytes_used:
            by_tool[m.tool] = by_tool.get(m.tool, 0) + m.bytes_used
    print("\nPer-tool  (sum across all that tool's databases and weights)", file=stream)
    print("─" * 60, file=stream)
    print(f"  {'Tool':<24} {'Size':>10}", file=stream)
    for tool in sorted(by_tool, key=lambda t: by_tool[t], reverse=True):
        print(f"  {tool:<24} {_humanise(by_tool[tool]):>10}", file=stream)


def _print_per_tier(measurements: list[Measurement], stream) -> None:
    """Cumulative footprint: tier=extended includes all base+extended;
    tier=full includes all base+extended+full. Matches how a user
    actually plans an install."""
    by_tier: dict[Tier, int] = {"base": 0, "extended": 0, "full": 0}
    for m in measurements:
        if not m.bytes_used:
            continue
        if m.tier == "base":
            by_tier["base"] += m.bytes_used
            by_tier["extended"] += m.bytes_used
            by_tier["full"] += m.bytes_used
        elif m.tier == "extended":
            by_tier["extended"] += m.bytes_used
            by_tier["full"] += m.bytes_used
        else:
            by_tier["full"] += m.bytes_used
    print("\nCumulative install size per tier", file=stream)
    print("─" * 60, file=stream)
    for t in ("base", "extended", "full"):
        print(f"  {t:<12} {_humanise(by_tier[t]):>10}", file=stream)


def _emit_markdown(measurements: list[Measurement], stream) -> None:
    """Markdown table you can paste into docs/how-to/install.md."""
    print("| Tool | Tier | Size |", file=stream)
    print("|---|---|---|", file=stream)
    by_tool_tier: dict[tuple[str, str], int] = {}
    for m in measurements:
        if m.bytes_used:
            by_tool_tier[(m.tool, m.tier)] = by_tool_tier.get((m.tool, m.tier), 0) + m.bytes_used
    for (tool, tier), size in sorted(by_tool_tier.items()):
        print(f"| {tool} | {tier} | {_humanise(size)} |", file=stream)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--data-root",
        default=os.path.expanduser("~/.ssign"),
        help="ssign data root (default: ~/.ssign).",
    )
    p.add_argument(
        "--markdown",
        action="store_true",
        help="Emit a markdown table suitable for pasting into install.md.",
    )
    args = p.parse_args(argv)

    db_root = _resolve_db_root(args.data_root)
    print(f"Data root: {args.data_root}")
    print(f"DB root:   {db_root}")

    measurements: list[Measurement] = []
    for entry in DATABASE_PATHS:
        measurements.append(_measure_database(entry, db_root))
    for entry in MODEL_WEIGHTS:
        measurements.append(_measure_weights(entry, args.data_root, db_root))

    # Warn about manifest entries the tool-map doesn't know about.
    unknown = [m.name for m in measurements if m.tool == "OTHER"]
    if unknown:
        print(f"\nWARNING: manifest entries not categorised in _TOOL_FOR_ENTRY: {unknown}", file=sys.stderr)

    _print_per_entry(measurements, sys.stdout)
    _print_per_tool(measurements, sys.stdout)
    _print_per_tier(measurements, sys.stdout)

    if args.markdown:
        print()
        _emit_markdown(measurements, sys.stdout)

    return 0


if __name__ == "__main__":
    sys.exit(main())
