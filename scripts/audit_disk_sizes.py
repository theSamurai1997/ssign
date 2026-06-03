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
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import get_args

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
from ssign_lib.resources import humanise_bytes  # noqa: E402

# 30 min per path. UniRef30 + nr can take ~10-15 min on gpfs metadata
# walks; 5 min was tight. Per-path, so total wallclock is bounded by the
# slowest path, not the sum (we run in parallel below).
_DU_TIMEOUT_S = 1800


@dataclass
class Measurement:
    name: str  # manifest entry name
    tool: str  # canonical tool name from the manifest's .tool field
    tier: Tier
    resolved_path: str | None
    bytes_used: int  # 0 when not fetched / unreachable
    note: str = ""  # populated when the path resolved but du failed, etc.


def _du_bytes(path: str) -> tuple[int, str]:
    """Return (bytes, note). bytes==0 with non-empty note on failure.

    Uses `du -sk` (kilobyte summary) instead of GNU's `-sb` — `-k` is
    POSIX, so the script runs on macOS as well as Linux. Multiplies the
    KB count by 1024 for byte accuracy.
    """
    try:
        out = subprocess.run(
            ["du", "-sk", path],
            capture_output=True,
            text=True,
            timeout=_DU_TIMEOUT_S,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        return (0, f"du failed: {e.stderr.strip()[:120]}")
    except subprocess.TimeoutExpired:
        return (0, f"du timed out (>{_DU_TIMEOUT_S // 60} min)")
    except FileNotFoundError:
        return (0, "du not on PATH")
    first = out.stdout.split()
    if not first or not first[0].isdigit():
        return (0, f"unexpected du output: {out.stdout[:120]}")
    return (int(first[0]) * 1024, "")


def _measure_database(entry: DatabasePath, db_root: str) -> Measurement:
    path = entry.resolve_path(db_root)
    if not path or not os.path.isdir(path):
        return Measurement(entry.name, entry.tool, entry.tier, None, 0, "not fetched")
    size, note = _du_bytes(path)
    return Measurement(entry.name, entry.tool, entry.tier, path, size, note)


def _measure_weights(entry: ModelWeights, data_root: str, db_root: str) -> Measurement:
    root = db_root if entry.under_db_root else data_root
    path = os.path.join(root, entry.default_subpath)
    if not os.path.exists(path):
        return Measurement(entry.name, entry.tool, entry.tier, None, 0, "not fetched")
    if os.path.isdir(path):
        size, note = _du_bytes(path)
    else:
        # Individual file (e.g., DeepSecE checkpoint)
        try:
            size = os.path.getsize(path)
            note = ""
        except OSError as e:
            size, note = 0, str(e)
    return Measurement(entry.name, entry.tool, entry.tier, path, size, note)


def _resolve_db_root(data_root: str) -> str:
    """Single source of truth for the db_root lookup: defer to
    `doctor.resolve_db_root` so we don't drift from its marker handling.
    """
    from ssign_app.scripts.doctor import resolve_db_root

    return resolve_db_root(data_root)


def _print_per_entry(measurements: list[Measurement], stream) -> None:
    print("\nPer-entry  (path resolution + actual on-disk size)", file=stream)
    print("─" * 60, file=stream)
    print(f"  {'Entry':<32} {'Tier':<10} {'Size':>10}  Path", file=stream)
    for m in measurements:
        size_str = humanise_bytes(m.bytes_used) if m.resolved_path else "(missing)"
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
        print(f"  {tool:<24} {humanise_bytes(by_tool[tool]):>10}", file=stream)


# Tier ordering for the cumulative rollup: a base entry counts toward
# every higher tier. Derived from the Tier Literal so adding a new tier
# only requires touching the manifest.
_TIER_ORDER: tuple[Tier, ...] = get_args(Tier)


def _print_per_tier(measurements: list[Measurement], stream) -> None:
    """Cumulative footprint: tier=extended includes all base+extended;
    tier=full includes all base+extended+full. Matches how a user
    actually plans an install."""
    by_tier: dict[Tier, int] = {t: 0 for t in _TIER_ORDER}
    for m in measurements:
        if not m.bytes_used:
            continue
        # Add to this entry's tier and every tier above it.
        idx = _TIER_ORDER.index(m.tier)
        for higher in _TIER_ORDER[idx:]:
            by_tier[higher] += m.bytes_used
    print("\nCumulative install size per tier", file=stream)
    print("─" * 60, file=stream)
    for t in _TIER_ORDER:
        print(f"  {t:<12} {humanise_bytes(by_tier[t]):>10}", file=stream)


def _emit_markdown(measurements: list[Measurement], stream) -> None:
    """Markdown table you can paste into docs/how-to/install.md."""
    print("| Tool | Tier | Size |", file=stream)
    print("|---|---|---|", file=stream)
    by_tool_tier: dict[tuple[str, str], int] = {}
    for m in measurements:
        if m.bytes_used:
            by_tool_tier[(m.tool, m.tier)] = by_tool_tier.get((m.tool, m.tier), 0) + m.bytes_used
    for (tool, tier), size in sorted(by_tool_tier.items()):
        print(f"| {tool} | {tier} | {humanise_bytes(size)} |", file=stream)


def _gather_measurements(data_root: str, db_root: str) -> list[Measurement]:
    """Fan out du calls across paths.

    Each measurement is dominated by gpfs metadata walk, which is
    I/O-bound and benefits from concurrency. 10 paths × ~10 min serial
    on gpfs was the original bottleneck; ~3 min wallclock at width=10.
    """
    tasks: list = []
    for entry in DATABASE_PATHS:
        tasks.append((_measure_database, entry, (db_root,)))
    for entry in MODEL_WEIGHTS:
        tasks.append((_measure_weights, entry, (data_root, db_root)))

    width = max(1, len(tasks))
    measurements: list[Measurement] = []
    with ThreadPoolExecutor(max_workers=width) as pool:
        futures = [pool.submit(fn, entry, *args) for fn, entry, args in tasks]
        measurements = [f.result() for f in futures]
    return measurements


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

    measurements = _gather_measurements(args.data_root, db_root)

    # The TestToolMapCoverage equivalent (every manifest entry has a
    # non-empty .tool) lives in tests/. Here we just surface the
    # uncategorised ones at runtime as a safety net.
    unknown = [m.name for m in measurements if not m.tool]
    if unknown:
        print(f"\nWARNING: manifest entries with empty .tool: {unknown}", file=sys.stderr)

    _print_per_entry(measurements, sys.stdout)
    _print_per_tool(measurements, sys.stdout)
    _print_per_tier(measurements, sys.stdout)

    if args.markdown:
        print()
        _emit_markdown(measurements, sys.stdout)

    return 0


if __name__ == "__main__":
    sys.exit(main())
