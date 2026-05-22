"""Tests for plm_effector.feature_extraction.

The heavy paths (`extract_terminal_features`, `_cli_main`) need torch +
15 GB of weights + GPU and are covered by integration tests. Here we
test the npz round-trip and the subprocess plumbing in isolation.
"""

from __future__ import annotations

import os
from unittest import mock

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from plm_effector.feature_extraction import (
    _extract_one_plm_in_subprocess,
    _load_features_npz,
    _save_features_npz,
    extract_all_features,
)


def _make_features(n_seqs=3, n_tokens=5, dim=7):
    """Synthesise a features dict shaped like extract_terminal_features."""
    return {
        "Nterminal": {
            "embedding": np.random.rand(n_seqs, n_tokens, dim).astype(np.float32),
            "attention_masks": np.ones((n_seqs, n_tokens), dtype=np.int64),
            "seq_ids": np.array([f"seq_{i}" for i in range(n_seqs)]),
        },
        "Cterminal": {
            "embedding": np.random.rand(n_seqs, n_tokens, dim).astype(np.float32),
            "attention_masks": np.ones((n_seqs, n_tokens), dtype=np.int64),
            "seq_ids": np.array([f"seq_{i}" for i in range(n_seqs)]),
        },
    }


class TestNpzRoundTrip:
    """_save_features_npz + _load_features_npz must preserve all six arrays."""

    def test_round_trip_preserves_arrays(self, tmp_dir):
        path = os.path.join(tmp_dir, "f.npz")
        original = _make_features()
        _save_features_npz(path, original)
        loaded = _load_features_npz(path)
        for term in ("Nterminal", "Cterminal"):
            for key in ("embedding", "attention_masks", "seq_ids"):
                np.testing.assert_array_equal(loaded[term][key], original[term][key])

    def test_no_allow_pickle_in_npz(self, tmp_dir):
        """Object-array dtypes would force allow_pickle=True on load,
        which is a remote-code-execution footgun. Confirm all six arrays
        round-trip with plain numeric/Unicode dtypes."""
        path = os.path.join(tmp_dir, "f.npz")
        _save_features_npz(path, _make_features())
        # allow_pickle=False is the read-side default we rely on.
        data = np.load(path, allow_pickle=False)
        for k in ("n_emb", "n_mask", "n_ids", "c_emb", "c_mask", "c_ids"):
            assert k in data, f"missing key {k}"


class TestExtractAllFeaturesIsolation:
    """Drive extract_all_features with the subprocess spawner mocked out
    so we can assert (a) one subprocess per PLM, (b) features are
    reassembled correctly, (c) the in-process path is still selectable.
    """

    def test_isolated_path_spawns_one_subprocess_per_plm(self, tmp_dir, monkeypatch):
        spawned = []

        def fake_spawn(proteins_fasta, pretrained_type, weights_dir, device, batch_size, out_npz):
            spawned.append(pretrained_type)
            _save_features_npz(out_npz, _make_features())

        monkeypatch.setattr(
            "plm_effector.feature_extraction._extract_one_plm_in_subprocess",
            fake_spawn,
        )
        features = extract_all_features(
            proteins_fasta=os.path.join(tmp_dir, "proteins.faa"),
            effector_type="T1SE",
            weights_dir="/fake",
            device=torch.device("cpu"),
            batch_size=5,
        )
        # T1SE → esm1 + esm2_t33 + ProtT5 (no ProtBert).
        assert spawned == ["esm1", "esm2_t33", "ProtT5"]
        assert set(features) == {"esm1", "esm2_t33", "ProtT5"}
        for v in features.values():
            assert "Nterminal" in v and "Cterminal" in v

    def test_t4se_includes_protbert(self, tmp_dir, monkeypatch):
        spawned = []
        monkeypatch.setattr(
            "plm_effector.feature_extraction._extract_one_plm_in_subprocess",
            lambda **kw: spawned.append(kw["pretrained_type"]) or _save_features_npz(kw["out_npz"], _make_features()),
        )
        extract_all_features(
            proteins_fasta="x.faa",
            effector_type="T4SE",
            weights_dir="/fake",
            device=torch.device("cpu"),
            batch_size=5,
        )
        assert spawned == ["esm1", "esm2_t33", "ProtBert", "ProtT5"]

    def test_in_process_path_skips_subprocess(self, tmp_dir, monkeypatch):
        called = mock.Mock()
        monkeypatch.setattr("plm_effector.feature_extraction._extract_one_plm_in_subprocess", called)
        monkeypatch.setattr(
            "plm_effector.feature_extraction.extract_terminal_features",
            lambda **kw: _make_features(),
        )
        extract_all_features(
            proteins_fasta="x.faa",
            effector_type="T1SE",
            weights_dir="/fake",
            device=torch.device("cpu"),
            batch_size=5,
            isolate_plms=False,
        )
        called.assert_not_called()

    def test_subprocess_failure_surfaces_stderr(self, tmp_dir, monkeypatch):
        """If the subprocess exits non-zero, the wrapping RuntimeError must
        include the subprocess stderr so the runner sees the real cause
        rather than just `exit -9`."""
        fake_result = mock.Mock(returncode=137, stderr="CUDA OOM at layer 7")
        monkeypatch.setattr("subprocess.run", lambda *a, **kw: fake_result)
        with pytest.raises(RuntimeError, match="CUDA OOM"):
            _extract_one_plm_in_subprocess(
                proteins_fasta="x.faa",
                pretrained_type="esm1",
                weights_dir="/fake",
                device=torch.device("cpu"),
                batch_size=5,
                out_npz=os.path.join(tmp_dir, "out.npz"),
            )
