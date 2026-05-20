"""``ssign doctor`` — verify a fresh ssign install on the user's machine.

Three classes of dependency, three sections of report:

1. Python packages + lazy symbols (the same list the integration test asserts).
   These should always pass after ``pip install ssign[extended]``; if any fail
   here, ``pyproject.toml`` is incomplete and CI should catch it before release.
2. External binaries that pip cannot install (HH-suite, BLAST+, InterProScan,
   EggNOG-mapper, optional DTU SignalP/DeepLocPro).
3. Databases + model weights on disk (gigabytes, fetched separately).

Each failure prints the exact fix command. Exit non-zero on any required
failure so ``ssign doctor && ssign run …`` works in scripts.
"""

from __future__ import annotations

import argparse
import importlib
import os
import shutil
import sys
from dataclasses import dataclass

from ssign_app.scripts.ssign_lib.dependency_manifest import (
    DatabasePath,
    ExternalBinary,
    ModelWeights,
    PythonDep,
    Tier,
    binaries_for_tier,
    databases_for_tier,
    deps_for_tier,
    weights_for_tier,
)

DEFAULT_DATA_ROOT = os.path.expanduser("~/.ssign")
DEFAULT_TIER: Tier = "extended"


@dataclass
class CheckResult:
    name: str
    ok: bool
    detail: str
    fix: str = ""


# ---------------------------------------------------------------------------
# Section 1 — Python imports
# ---------------------------------------------------------------------------


def _resolve_symbol(dotted: str) -> None:
    """Import ``a.b.c.D`` by importing ``a.b.c`` then resolving attribute ``D``.

    Raises ``ImportError`` / ``AttributeError`` on failure (the caller catches).
    """
    parts = dotted.split(".")
    # walk back until we hit something importable; the rest are attributes
    for split in range(len(parts), 0, -1):
        mod_name = ".".join(parts[:split])
        try:
            mod = importlib.import_module(mod_name)
        except ImportError:
            continue
        obj = mod
        for attr in parts[split:]:
            obj = getattr(obj, attr)
        return
    raise ImportError(f"could not import any prefix of {dotted!r}")


def check_python_dep(dep: PythonDep) -> CheckResult:
    try:
        importlib.import_module(dep.module)
    except Exception as e:  # noqa: BLE001 — surface anything from native ext loaders
        return CheckResult(
            name=dep.module,
            ok=False,
            detail=f"{type(e).__name__}: {e}",
            fix=f"pip install {dep.pip_name}",
        )
    for sym in dep.symbols:
        try:
            _resolve_symbol(sym)
        except Exception as e:  # noqa: BLE001
            return CheckResult(
                name=f"{dep.module} (symbol {sym})",
                ok=False,
                detail=f"{type(e).__name__}: {e}",
                fix=f"pip install --upgrade {dep.pip_name}  # symbol moved/missing in installed version",
            )
    return CheckResult(name=dep.module, ok=True, detail="")


# ---------------------------------------------------------------------------
# Section 2 — External binaries
# ---------------------------------------------------------------------------


def check_external_binary(b: ExternalBinary) -> CheckResult:
    found = shutil.which(b.binary)
    if found:
        return CheckResult(name=b.name, ok=True, detail=found)
    return CheckResult(
        name=b.name + (" [optional]" if b.optional else ""),
        ok=b.optional,  # optional missing = treated as OK for exit code
        detail=f"`{b.binary}` not on PATH",
        fix=b.install_hint,
    )


# ---------------------------------------------------------------------------
# Section 3 — Databases + model weights
# ---------------------------------------------------------------------------


def _resolve_db_path(d: DatabasePath, data_root: str) -> str:
    env_value = os.environ.get(d.env_var, "").strip()
    if env_value:
        return env_value
    return os.path.join(data_root, "databases", d.default_subpath)


def check_database(d: DatabasePath, data_root: str) -> CheckResult:
    resolved = _resolve_db_path(d, data_root)
    sentinel = os.path.join(resolved, d.sentinel_file)
    if os.path.isfile(sentinel):
        return CheckResult(name=d.name, ok=True, detail=resolved)
    if os.path.isdir(resolved):
        return CheckResult(
            name=d.name,
            ok=False,
            detail=f"{resolved} exists but sentinel {d.sentinel_file!r} is missing",
            fix=d.install_hint,
        )
    return CheckResult(
        name=d.name,
        ok=False,
        detail=f"not found at {resolved} (env: ${d.env_var})",
        fix=d.install_hint,
    )


def check_weights(w: ModelWeights, data_root: str) -> CheckResult:
    resolved = os.path.join(data_root, w.default_subpath)
    if os.path.exists(resolved):
        return CheckResult(name=w.name, ok=True, detail=resolved)
    return CheckResult(
        name=w.name,
        ok=False,
        detail=f"not found at {resolved}",
        fix=w.install_hint,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def _render(section: str, results: list[CheckResult], stream) -> tuple[int, int]:
    ok_count = sum(1 for r in results if r.ok)
    total = len(results)
    print(f"\n{section}  ({ok_count}/{total} OK)", file=stream)
    print("─" * len(section), file=stream)
    for r in results:
        marker = "OK  " if r.ok else "FAIL"
        print(f"  {marker}  {r.name}", file=stream)
        if r.detail:
            print(f"           {r.detail}", file=stream)
        if not r.ok and r.fix:
            print(f"           fix: {r.fix}", file=stream)
    return ok_count, total


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def run(
    tier: Tier = DEFAULT_TIER,
    imports_only: bool = False,
    data_root: str = DEFAULT_DATA_ROOT,
    stream=sys.stdout,
) -> int:
    """Run all checks. Returns 0 if every required check passed, 1 otherwise.

    ``imports_only=True`` skips binaries / DBs / weights — used by CI where
    only the Python environment can be verified.
    """
    print(f"ssign doctor — checking tier '{tier}'", file=stream)
    print(f"  data root: {data_root}  (override sub-paths via SSIGN_* env vars)", file=stream)

    failures = 0

    py_results = [check_python_dep(d) for d in deps_for_tier(tier)]
    ok, total = _render("Python packages", py_results, stream)
    failures += total - ok

    if imports_only:
        print(f"\nResult: {failures} failures (imports-only).", file=stream)
        return 1 if failures else 0

    bin_results = [check_external_binary(b) for b in binaries_for_tier(tier)]
    ok, total = _render("External binaries", bin_results, stream)
    # optional-missing was already counted as ok=True in check_external_binary
    failures += sum(1 for r in bin_results if not r.ok)

    db_results = [check_database(d, data_root) for d in databases_for_tier(tier)]
    ok, total = _render("Databases", db_results, stream)
    failures += total - ok

    w_results = [check_weights(w, data_root) for w in weights_for_tier(tier)]
    ok, total = _render("Model weights", w_results, stream)
    failures += total - ok

    print(file=stream)
    if failures:
        print(f"{failures} check(s) failed.  Address the `fix:` lines above and re-run.", file=stream)
        return 1
    print("All checks passed.", file=stream)
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="ssign doctor",
        description="Verify a fresh ssign install: Python packages, external binaries, databases, model weights.",
    )
    p.add_argument(
        "--tier",
        choices=("base", "extended", "full"),
        default=DEFAULT_TIER,
        help=f"Install tier to verify against (default: {DEFAULT_TIER}).",
    )
    p.add_argument(
        "--imports-only",
        action="store_true",
        help="Only check Python imports (skip binaries, DBs, weights). Used by CI.",
    )
    p.add_argument(
        "--data-root",
        default=DEFAULT_DATA_ROOT,
        help=f"Root directory for databases + models (default: {DEFAULT_DATA_ROOT}).",
    )
    args = p.parse_args(argv)
    return run(tier=args.tier, imports_only=args.imports_only, data_root=args.data_root)


if __name__ == "__main__":
    sys.exit(main())
