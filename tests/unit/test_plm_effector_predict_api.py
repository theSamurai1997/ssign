"""Tests for plm_effector.predict_api.

`write_predictions_tsv` is already covered by the existing
`tests/unit/test_run_plm_effector.py::TestWritePredictionsTsv` class.
This file fills in the rest:

- `_resolve_device`: torch.device coercion + CUDA-availability gate.
- `predict`: argument validation paths (effector_types / missing input /
  missing weights). The orchestration body (extract_all_features →
  run_ensemble → write) requires GPU + 15 GB of weights and is
  exercised by `tests/integration/test_run_plm_effector_integration.py`.
"""

import importlib.util
import os
import sys
from types import SimpleNamespace

import pytest

# predict_api imports torch lazily — this module-level import is fine
# without torch installed. Tests that exercise _resolve_device or the
# full orchestration mark themselves needs_torch below.
from plm_effector.predict_api import (
    _VALID_EFFECTOR_TYPES,
    _resolve_device,
    predict,
)

needs_torch = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None,
    reason="torch not installed (covered by tests/integration/)",
)


@needs_torch
class TestResolveDevice:
    def test_cpu_returns_cpu_device(self):
        import torch

        result = _resolve_device("cpu")
        assert result == torch.device("cpu")

    def test_passthrough_torch_device(self):
        import torch

        device = torch.device("cpu")
        assert _resolve_device(device) is device

    def test_unrecognised_string_raises_value_error(self):
        with pytest.raises(ValueError, match="Unrecognised device"):
            _resolve_device("tpu")

    def test_cuda_when_unavailable_raises_runtime_error(self, monkeypatch):
        import torch

        monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
        with pytest.raises(RuntimeError, match="no CUDA device"):
            _resolve_device("cuda")

    def test_default_none_treated_as_cuda(self, monkeypatch):
        import torch

        # Same gate as "cuda" — default value picks GPU when available
        monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
        with pytest.raises(RuntimeError, match="no CUDA device"):
            _resolve_device(None)

    def test_cuda_when_available_returns_cuda_device(self, monkeypatch):
        import torch

        monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
        result = _resolve_device("cuda")
        assert result.type == "cuda"


class TestPredictValidation:
    def test_invalid_effector_type_raises_value_error(self, tmp_dir):
        input_fasta = os.path.join(tmp_dir, "in.faa")
        with open(input_fasta, "w") as f:
            f.write(">P1\nMKT\n")
        weights = os.path.join(tmp_dir, "weights")
        os.makedirs(weights)
        with pytest.raises(ValueError, match="unrecognised effector_types"):
            predict(
                proteins_fasta=input_fasta,
                weights_dir=weights,
                effector_types=["T9SE"],
                out_dir=os.path.join(tmp_dir, "out"),
            )

    def test_empty_effector_types_raises(self, tmp_dir):
        input_fasta = os.path.join(tmp_dir, "in.faa")
        with open(input_fasta, "w") as f:
            f.write(">P1\nMKT\n")
        weights = os.path.join(tmp_dir, "weights")
        os.makedirs(weights)
        with pytest.raises(ValueError, match="effector_types is empty"):
            predict(
                proteins_fasta=input_fasta,
                weights_dir=weights,
                effector_types=[],
                out_dir=os.path.join(tmp_dir, "out"),
            )

    def test_missing_input_raises_file_not_found(self, tmp_dir):
        weights = os.path.join(tmp_dir, "weights")
        os.makedirs(weights)
        with pytest.raises(FileNotFoundError, match="Input FASTA not found"):
            predict(
                proteins_fasta=os.path.join(tmp_dir, "does_not_exist.faa"),
                weights_dir=weights,
                effector_types=["T1SE"],
                out_dir=os.path.join(tmp_dir, "out"),
            )

    def test_missing_weights_dir_raises_file_not_found(self, tmp_dir):
        input_fasta = os.path.join(tmp_dir, "in.faa")
        with open(input_fasta, "w") as f:
            f.write(">P1\nMKT\n")
        with pytest.raises(FileNotFoundError, match="weights directory not found"):
            predict(
                proteins_fasta=input_fasta,
                weights_dir=os.path.join(tmp_dir, "no_weights"),
                effector_types=["T1SE"],
                out_dir=os.path.join(tmp_dir, "out"),
            )

    @needs_torch
    @pytest.mark.parametrize("effector_type", _VALID_EFFECTOR_TYPES)
    def test_every_canonical_effector_type_accepted(self, tmp_dir, effector_type, monkeypatch):
        """All five canonical effector types must pass validation when
        requested singly. We mock the downstream extract → ensemble
        pipeline and verify validation lets the call through to
        write_predictions_tsv with the correct effector_type tag."""
        import numpy as np
        import plm_effector.predict_api as predict_api

        input_fasta = os.path.join(tmp_dir, "in.faa")
        with open(input_fasta, "w") as f:
            f.write(">P1\nMKT\n")
        weights = os.path.join(tmp_dir, "weights")
        os.makedirs(weights)

        # Stub feature_extraction and ensemble so predict() doesn't need a real PLM.
        plms = ["esm1", "esm2_t33", "ProtT5", "ProtBert"] if effector_type == "T4SE" else ["esm1", "esm2_t33", "ProtT5"]

        def fake_iter_chunks(chunk_paths, **_kw):
            n_chunks = len(next(iter(chunk_paths.values())))
            for _ in range(n_chunks):
                yield {pt: {"features": "irrelevant"} for pt in chunk_paths}

        def fake_pretrained_types_for(effector_types):
            if isinstance(effector_types, str):
                effector_types = [effector_types]
            out = ["esm1", "esm2_t33", "ProtT5"]
            if "T4SE" in effector_types:
                out.append("ProtBert")
            return out

        fake_extract = SimpleNamespace(
            extract_all_features=lambda **_kw: {pt: ["/fake/chunk0.npz"] for pt in plms},
            iter_chunk_features=fake_iter_chunks,
            pretrained_types_for=fake_pretrained_types_for,
        )
        n_base = 8 if effector_type == "T4SE" else 6
        fake_ensemble = SimpleNamespace(
            run_ensemble=lambda **_kw: (
                np.array([">P1"]),
                np.zeros((1, n_base)),
                np.array([0.5]),
                np.array([True]),
            ),
        )
        monkeypatch.setitem(sys.modules, "plm_effector.feature_extraction", fake_extract)
        monkeypatch.setitem(sys.modules, "plm_effector.ensemble", fake_ensemble)
        import torch

        monkeypatch.setattr(predict_api, "_resolve_device", lambda _d: torch.device("cpu"))
        monkeypatch.setattr(predict_api.os, "remove", lambda _p: None)

        out_dir = os.path.join(tmp_dir, "out")
        summary = predict(
            proteins_fasta=input_fasta,
            weights_dir=weights,
            effector_types=[effector_type],
            out_dir=out_dir,
            device="cpu",
        )
        assert effector_type in summary
        n_positive, out_path = summary[effector_type]
        assert n_positive == 1
        with open(out_path) as f:
            lines = f.readlines()
        assert lines[1].rstrip("\n").split("\t")[-1] == effector_type


@needs_torch
class TestPredictMultiType:
    """Verify multi-type predict extracts features once and runs ensemble
    per type, producing one TSV per requested effector type."""

    def test_multi_type_extracts_features_once(self, tmp_dir, monkeypatch):
        import numpy as np
        import plm_effector.predict_api as predict_api

        input_fasta = os.path.join(tmp_dir, "in.faa")
        with open(input_fasta, "w") as f:
            f.write(">P1\nMKT\n")
        weights = os.path.join(tmp_dir, "weights")
        os.makedirs(weights)

        extract_calls = []
        ensemble_calls = []

        def fake_extract_all_features(**kw):
            extract_calls.append(kw["pretrained_types"])
            return {pt: ["/fake/chunk0.npz"] for pt in kw["pretrained_types"]}

        def fake_iter_chunks(chunk_paths, **_kw):
            n_chunks = len(next(iter(chunk_paths.values())))
            for _ in range(n_chunks):
                yield {pt: {"features": "irrelevant"} for pt in chunk_paths}

        def fake_pretrained_types_for(effector_types):
            out = ["esm1", "esm2_t33", "ProtT5"]
            if "T4SE" in effector_types:
                out.append("ProtBert")
            return out

        def fake_run_ensemble(**kw):
            ensemble_calls.append(kw["effector_type"])
            n_base = 8 if kw["effector_type"] == "T4SE" else 6
            return (
                np.array([">P1"]),
                np.zeros((1, n_base)),
                np.array([0.5]),
                np.array([True]),
            )

        fake_extract = SimpleNamespace(
            extract_all_features=fake_extract_all_features,
            iter_chunk_features=fake_iter_chunks,
            pretrained_types_for=fake_pretrained_types_for,
        )
        fake_ensemble = SimpleNamespace(run_ensemble=fake_run_ensemble)
        monkeypatch.setitem(sys.modules, "plm_effector.feature_extraction", fake_extract)
        monkeypatch.setitem(sys.modules, "plm_effector.ensemble", fake_ensemble)
        import torch

        monkeypatch.setattr(predict_api, "_resolve_device", lambda _d: torch.device("cpu"))
        monkeypatch.setattr(predict_api.os, "remove", lambda _p: None)

        out_dir = os.path.join(tmp_dir, "out")
        summary = predict(
            proteins_fasta=input_fasta,
            weights_dir=weights,
            effector_types=["T1SE", "T2SE", "T3SE", "T4SE", "T6SE"],
            out_dir=out_dir,
            device="cpu",
        )

        # The key savings claim: features extracted ONCE across all types.
        assert len(extract_calls) == 1
        assert extract_calls[0] == ["esm1", "esm2_t33", "ProtT5", "ProtBert"]

        # One ensemble run per type (one chunk × 5 types).
        assert ensemble_calls == ["T1SE", "T2SE", "T3SE", "T4SE", "T6SE"]

        # One TSV per type.
        assert set(summary) == {"T1SE", "T2SE", "T3SE", "T4SE", "T6SE"}
        for eff_type, (n_positive, out_path) in summary.items():
            assert n_positive == 1
            assert os.path.exists(out_path)
            assert os.path.basename(out_path) == f"{eff_type}.tsv"
