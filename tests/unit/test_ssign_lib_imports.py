"""Guard against bare ``from ssign_lib.X import ...`` imports inside ssign_lib/.

Why this matters: every module under
``src/ssign_app/scripts/ssign_lib/`` is imported in two different
contexts -- as ``ssign_app.scripts.ssign_lib.X`` (in-process by the
runner) and as ``ssign_lib.X`` (by scripts that prepend ``scripts/``
to sys.path before running standalone). Bare absolute imports work
in the second context but fail in the first when an earlier in-process
step has not already leaked ``scripts/`` onto sys.path. This caused a
multi-genome ``build_passenger_fasta`` crash on 2026-06-05
(CX3 job 2934238). Relative imports (``from .X import ...``) work in
both contexts.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

SSIGN_LIB_DIR = Path(__file__).parent.parent.parent / "src" / "ssign_app" / "scripts" / "ssign_lib"


def _list_ssign_lib_python_files() -> list[Path]:
    return sorted(p for p in SSIGN_LIB_DIR.glob("*.py") if p.name != "__init__.py")


def _bare_ssign_lib_imports(path: Path) -> list[str]:
    """Real import statements importing ssign_lib as a top-level package.

    AST-based so docstring snippets like ``from ssign_lib.foo import ...``
    inside Usage examples don't trip the check.
    """
    tree = ast.parse(path.read_text(), filename=str(path))
    offenders: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.level == 0:
            mod = node.module or ""
            if mod == "ssign_lib" or mod.startswith("ssign_lib."):
                offenders.append(f"line {node.lineno}: from {mod} import ...")
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "ssign_lib" or alias.name.startswith("ssign_lib."):
                    offenders.append(f"line {node.lineno}: import {alias.name}")
    return offenders


@pytest.mark.parametrize("path", _list_ssign_lib_python_files(), ids=lambda p: p.name)
def test_no_bare_ssign_lib_imports(path: Path) -> None:
    offenders = _bare_ssign_lib_imports(path)
    assert not offenders, (
        f"{path.name} uses bare 'ssign_lib...' imports:\n  "
        + "\n  ".join(offenders)
        + "\nUse 'from .X import ...' (relative) so the module loads both via "
        "the ssign_app namespace and as a standalone script. See "
        "src/ssign_app/scripts/ssign_lib/t5_passenger.py commit a934e48 "
        "for the original bug and fix."
    )


def test_guard_covers_at_least_one_file() -> None:
    """Sanity: parameterised test must run on at least one file.

    If ssign_lib/ ever empties (move/rename), this fails loudly rather
    than the parametrised test silently passing on zero inputs.
    """
    assert _list_ssign_lib_python_files(), f"No .py files found under {SSIGN_LIB_DIR}"
