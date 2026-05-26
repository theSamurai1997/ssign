"""Tests for plm_effector.feature_extraction.

The heavy paths (`iter_terminal_feature_chunks`, `_cli_main`) need torch +
15 GB of weights + GPU and are covered by integration tests. Here we
test the npz round-trip, the chunk plumbing, and the subprocess wiring
in isolation.
"""

from __future__ import annotations

import os
from unittest import mock

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from plm_effector.feature_extraction import (  # noqa: F401  (pretrained_types_for used by class body)
    _chunk_path,
    _discover_chunk_paths,
    _extract_one_plm_in_subprocess,
    _load_features_npz,
    _save_features_npz,
    extract_all_features,
    extract_terminal_features,
    pretrained_types_for,
)


def _make_features(n_seqs=3, n_tokens=5, dim=7, seed=0):
    """Synthesise a features dict shaped like one chunk's worth of output."""
    rng = np.random.default_rng(seed)
    return {
        "Nterminal": {
            "embedding": rng.random((n_seqs, n_tokens, dim), dtype=np.float32),
            "attention_masks": np.ones((n_seqs, n_tokens), dtype=np.int64),
            "seq_ids": np.array([f"seq_{seed}_{i}" for i in range(n_seqs)]),
        },
        "Cterminal": {
            "embedding": rng.random((n_seqs, n_tokens, dim), dtype=np.float32),
            "attention_masks": np.ones((n_seqs, n_tokens), dtype=np.int64),
            "seq_ids": np.array([f"seq_{seed}_{i}" for i in range(n_seqs)]),
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


class TestChunkPaths:
    def test_chunk_path_zero_padding(self, tmp_dir):
        # Sort-friendly zero-padded chunk index — 4 digits handles up to 10k
        # chunks before lexical order breaks.
        assert _chunk_path(tmp_dir, "esm1", 0).endswith("esm1_chunk0000.npz")
        assert _chunk_path(tmp_dir, "esm1", 17).endswith("esm1_chunk0017.npz")

    def test_discover_returns_sorted(self, tmp_dir):
        for idx in (2, 0, 1):
            path = _chunk_path(tmp_dir, "esm1", idx)
            _save_features_npz(path, _make_features(seed=idx))
        paths = _discover_chunk_paths(tmp_dir, "esm1")
        assert [os.path.basename(p) for p in paths] == [
            "esm1_chunk0000.npz",
            "esm1_chunk0001.npz",
            "esm1_chunk0002.npz",
        ]

    def test_discover_filters_by_prefix(self, tmp_dir):
        # Two PLMs writing into the same directory should not see each other's chunks.
        _save_features_npz(_chunk_path(tmp_dir, "esm1", 0), _make_features())
        _save_features_npz(_chunk_path(tmp_dir, "esm2_t33", 0), _make_features())
        assert len(_discover_chunk_paths(tmp_dir, "esm1")) == 1
        assert len(_discover_chunk_paths(tmp_dir, "esm2_t33")) == 1


class TestPretrainedTypesFor:
    """pretrained_types_for unions PLM requirements across effector types."""

    def test_single_non_t4se_excludes_protbert(self):
        assert pretrained_types_for("T1SE") == ["esm1", "esm2_t33", "ProtT5"]

    def test_t4se_adds_protbert(self):
        assert pretrained_types_for("T4SE") == ["esm1", "esm2_t33", "ProtT5", "ProtBert"]

    def test_union_adds_protbert_when_any_type_needs_it(self):
        plms = pretrained_types_for(["T1SE", "T2SE", "T3SE", "T4SE", "T6SE"])
        assert plms == ["esm1", "esm2_t33", "ProtT5", "ProtBert"]

    def test_union_omits_protbert_without_t4se(self):
        plms = pretrained_types_for(["T1SE", "T2SE", "T3SE", "T6SE"])
        assert plms == ["esm1", "esm2_t33", "ProtT5"]


class TestExtractAllFeaturesIsolation:
    """Drive extract_all_features with the subprocess spawner mocked out
    so we can assert (a) one subprocess per PLM, (b) chunk paths are
    returned correctly, (c) the in-process path is still selectable.
    """

    def test_isolated_path_spawns_one_subprocess_per_plm(self, tmp_dir, monkeypatch):
        spawned = []

        def fake_spawn(*, proteins_fasta, pretrained_type, weights_dir, device, batch_size, chunk_size, out_dir):
            spawned.append(pretrained_type)
            # Simulate the subprocess writing two chunks
            for chunk_idx in range(2):
                _save_features_npz(
                    _chunk_path(out_dir, pretrained_type, chunk_idx),
                    _make_features(seed=chunk_idx),
                )
            return _discover_chunk_paths(out_dir, pretrained_type)

        monkeypatch.setattr(
            "plm_effector.feature_extraction._extract_one_plm_in_subprocess",
            fake_spawn,
        )
        chunk_paths = extract_all_features(
            proteins_fasta=os.path.join(tmp_dir, "proteins.faa"),
            pretrained_types=["esm1", "esm2_t33", "ProtT5"],
            weights_dir="/fake",
            device=torch.device("cpu"),
            feature_cache_dir=tmp_dir,
            batch_size=5,
            chunk_size=8,
        )
        assert spawned == ["esm1", "esm2_t33", "ProtT5"]
        assert set(chunk_paths) == {"esm1", "esm2_t33", "ProtT5"}
        for pt, paths in chunk_paths.items():
            assert len(paths) == 2, f"{pt}: expected 2 chunks"
            for p in paths:
                assert os.path.exists(p)

    def test_extracts_protbert_when_requested(self, tmp_dir, monkeypatch):
        spawned = []

        def fake_spawn(*, proteins_fasta, pretrained_type, weights_dir, device, batch_size, chunk_size, out_dir):
            spawned.append(pretrained_type)
            _save_features_npz(_chunk_path(out_dir, pretrained_type, 0), _make_features())
            return _discover_chunk_paths(out_dir, pretrained_type)

        monkeypatch.setattr(
            "plm_effector.feature_extraction._extract_one_plm_in_subprocess",
            fake_spawn,
        )
        extract_all_features(
            proteins_fasta="x.faa",
            pretrained_types=["esm1", "esm2_t33", "ProtT5", "ProtBert"],
            weights_dir="/fake",
            device=torch.device("cpu"),
            feature_cache_dir=tmp_dir,
            batch_size=5,
            chunk_size=8,
        )
        assert spawned == ["esm1", "esm2_t33", "ProtT5", "ProtBert"]

    def test_in_process_path_skips_subprocess(self, tmp_dir, monkeypatch):
        called = mock.Mock()
        monkeypatch.setattr("plm_effector.feature_extraction._extract_one_plm_in_subprocess", called)

        def fake_iter(*, proteins_fasta, pretrained_type, weights_dir, device, batch_size, chunk_size):
            # Two synthetic chunks per PLM
            yield _make_features(seed=0)
            yield _make_features(seed=1)

        monkeypatch.setattr(
            "plm_effector.feature_extraction.iter_terminal_feature_chunks",
            fake_iter,
        )
        chunk_paths = extract_all_features(
            proteins_fasta="x.faa",
            pretrained_types=["esm1", "esm2_t33", "ProtT5"],
            weights_dir="/fake",
            device=torch.device("cpu"),
            feature_cache_dir=tmp_dir,
            batch_size=5,
            chunk_size=8,
            isolate_plms=False,
        )
        called.assert_not_called()
        for pt, paths in chunk_paths.items():
            assert len(paths) == 2

    def _fake_popen(self, returncode: int, stderr_lines: list):
        """Build a Popen-like mock that yields stderr lines then exits with returncode."""
        proc = mock.Mock()
        proc.stderr = iter(stderr_lines)
        proc.wait = mock.Mock(return_value=returncode)
        return proc

    def test_subprocess_failure_surfaces_stderr(self, tmp_dir, monkeypatch):
        """If the subprocess exits non-zero, the wrapping RuntimeError must
        include the subprocess stderr tail so the runner sees the real cause
        rather than just `exit -9`."""
        fake_proc = self._fake_popen(137, ["loading ESM-1b...\n", "CUDA OOM at layer 7\n"])
        monkeypatch.setattr("subprocess.Popen", lambda *a, **kw: fake_proc)
        with pytest.raises(RuntimeError, match="CUDA OOM"):
            _extract_one_plm_in_subprocess(
                proteins_fasta="x.faa",
                pretrained_type="esm1",
                weights_dir="/fake",
                device=torch.device("cpu"),
                batch_size=5,
                chunk_size=8,
                out_dir=tmp_dir,
            )

    def test_subprocess_zero_exit_but_no_chunks_is_error(self, tmp_dir, monkeypatch):
        """Subprocess returning 0 but writing nothing means something silent
        went wrong (e.g. empty FASTA filtered upstream). Must error rather
        than return empty chunk list and produce a downstream KeyError."""
        fake_proc = self._fake_popen(0, [])
        monkeypatch.setattr("subprocess.Popen", lambda *a, **kw: fake_proc)
        with pytest.raises(RuntimeError, match="wrote no chunk files"):
            _extract_one_plm_in_subprocess(
                proteins_fasta="x.faa",
                pretrained_type="esm1",
                weights_dir="/fake",
                device=torch.device("cpu"),
                batch_size=5,
                chunk_size=8,
                out_dir=tmp_dir,
            )


class TestChunkingEquivalence:
    """Per-protein outputs must be bit-identical between chunked and
    non-chunked feature extraction. Mocks the model/tokenizer at the
    `_load_model_and_tokenizer` boundary with a deterministic stand-in
    so we can test the chunk-slicing logic without a real PLM.
    """

    def _install_fake_model(self, monkeypatch, dim=4):
        """Replace the model with one that returns deterministic per-token
        features derived from input_ids, so we can check that chunking
        preserves per-protein outputs.
        """
        import plm_effector.feature_extraction as fe

        class FakeOutputs:
            def __init__(self, last_hidden_state):
                self.last_hidden_state = last_hidden_state

        class FakeModel:
            def __init__(self, dim):
                self.dim = dim

            def __call__(self, input_ids=None, attention_mask=None, **kw):
                # Per-token features = input_ids broadcast over a deterministic
                # ramp, so different proteins land on different feature values.
                ids = input_ids.float().unsqueeze(-1)
                ramp = torch.arange(self.dim, dtype=torch.float32).view(1, 1, -1)
                return FakeOutputs(ids * 0.01 + ramp)

            def to(self, device):
                return self

            def eval(self):
                return self

        class FakeTokenizer:
            padding_side = "right"
            pad_token = "X"

            def __call__(self, sequences, return_tensors="pt", padding="max_length", truncation=True, max_length=8):
                # Token id = position-based hash of residue char so tokens
                # depend on the actual sequence content, not just length.
                ids = []
                masks = []
                for s in sequences:
                    s_clean = s.replace(" ", "")[:max_length]
                    tok = [1 + (ord(c) % 20) for c in s_clean]
                    pad_len = max_length - len(tok)
                    if self.padding_side == "left":
                        tok = [0] * pad_len + tok
                        mask = [0] * pad_len + [1] * (max_length - pad_len)
                    else:
                        tok = tok + [0] * pad_len
                        mask = [1] * (max_length - pad_len) + [0] * pad_len
                    ids.append(tok)
                    masks.append(mask)
                return {
                    "input_ids": torch.tensor(ids, dtype=torch.long),
                    "attention_mask": torch.tensor(masks, dtype=torch.long),
                }

        monkeypatch.setattr(
            fe,
            "_load_model_and_tokenizer",
            lambda *a, **k: (FakeModel(dim), FakeTokenizer()),
        )
        # Use short sequences so max_length fits comfortably
        monkeypatch.setattr(fe, "_MODEL_MAX_LENGTH", {**fe._MODEL_MAX_LENGTH, "esm1": 8})
        monkeypatch.setattr(fe, "_MODEL_REAL_SEQUENCE_LEN", {**fe._MODEL_REAL_SEQUENCE_LEN, "esm1": 6})

    def test_chunked_concat_matches_all_at_once(self, tmp_dir, monkeypatch):
        self._install_fake_model(monkeypatch)
        fasta = os.path.join(tmp_dir, "p.faa")
        with open(fasta, "w") as f:
            for i in range(10):
                # Vary residue composition so per-protein outputs differ
                seq = "ACDE"[i % 4] * (3 + (i % 4))
                f.write(f">prot_{i}\n{seq}\n")

        device = torch.device("cpu")
        # All-at-once: chunk_size larger than n_proteins → one chunk
        all_at_once = extract_terminal_features(
            proteins_fasta=fasta,
            pretrained_type="esm1",
            weights_dir="/fake",
            device=device,
            batch_size=2,
            chunk_size=100,
        )
        # Chunked: 10 proteins in chunks of 3 → 4 chunks
        chunked = extract_terminal_features(
            proteins_fasta=fasta,
            pretrained_type="esm1",
            weights_dir="/fake",
            device=device,
            batch_size=2,
            chunk_size=3,
        )

        for terminal in ("Nterminal", "Cterminal"):
            np.testing.assert_array_equal(chunked[terminal]["seq_ids"], all_at_once[terminal]["seq_ids"])
            np.testing.assert_array_equal(
                chunked[terminal]["attention_masks"],
                all_at_once[terminal]["attention_masks"],
            )
            np.testing.assert_allclose(
                chunked[terminal]["embedding"],
                all_at_once[terminal]["embedding"],
                rtol=0,
                atol=0,
            )
