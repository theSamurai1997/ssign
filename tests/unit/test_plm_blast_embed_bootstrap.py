"""Unit tests for _plm_blast_embed_bootstrap.py.

Covers the torch.distributed / torch.multiprocessing stub-installation
logic. Does NOT exercise the runpy hand-off to upstream embeddings.py
(that requires a real pLM-BLAST clone + ProtT5 weights — covered by
the integration test).
"""

from __future__ import annotations

import builtins
import sys
import types

import _plm_blast_embed_bootstrap as bootstrap
import pytest


class TestInstallTorchSubmoduleStubs:
    def test_real_torch_distributed_present_returns_empty(self):
        """When the host torch wheel ships the real submodules, no stubbing."""
        try:
            import torch.distributed  # noqa: F401
            import torch.multiprocessing  # noqa: F401
        except (ImportError, ModuleNotFoundError):
            pytest.skip("Host torch lacks distributed/multiprocessing; can't test happy path")
        assert bootstrap.install_torch_submodule_stubs() == []

    def test_missing_torch_returns_empty(self, monkeypatch):
        """When torch itself is absent, the bootstrap is a no-op — the
        downstream embeddings.py import will fail loudly, which is correct."""
        monkeypatch.delitem(sys.modules, "torch", raising=False)
        monkeypatch.delitem(sys.modules, "torch.distributed", raising=False)
        monkeypatch.delitem(sys.modules, "torch.multiprocessing", raising=False)
        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "torch" or name.startswith("torch."):
                raise ImportError(f"{name} not installed (test stub)")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        assert bootstrap.install_torch_submodule_stubs() == []

    def test_stubs_missing_distributed(self, monkeypatch):
        """When torch.distributed is missing, install a stub and report it."""
        try:
            import torch  # noqa: F401
        except (ImportError, ModuleNotFoundError):
            pytest.skip("torch not installed")

        monkeypatch.delitem(sys.modules, "torch.distributed", raising=False)
        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "torch.distributed":
                raise ModuleNotFoundError("No module named 'torch.distributed' (test stub)")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        stubbed = bootstrap.install_torch_submodule_stubs()
        assert "torch.distributed" in stubbed
        assert isinstance(sys.modules["torch.distributed"], types.ModuleType)

    def test_stub_has_no_init_process_group(self, monkeypatch):
        """The stub must NOT silently allow multi-GPU calls — accessing
        init_process_group on the stub should raise AttributeError so the
        failure point is the real call site, not an obscure runtime hang."""
        try:
            import torch  # noqa: F401
        except (ImportError, ModuleNotFoundError):
            pytest.skip("torch not installed")

        monkeypatch.delitem(sys.modules, "torch.distributed", raising=False)
        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "torch.distributed":
                raise ModuleNotFoundError("stripped wheel (test stub)")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", fake_import)
        bootstrap.install_torch_submodule_stubs()
        stub = sys.modules["torch.distributed"]
        with pytest.raises(AttributeError):
            _ = stub.init_process_group  # type: ignore[attr-defined]


class TestMain:
    def test_no_args_returns_2(self, capsys):
        assert bootstrap.main([]) == 2
        captured = capsys.readouterr()
        assert "Usage:" in captured.err

    def test_missing_script_returns_2(self, tmp_path, capsys):
        missing = str(tmp_path / "nope_embeddings.py")
        assert bootstrap.main([missing]) == 2
        captured = capsys.readouterr()
        assert "embeddings.py not found" in captured.err
