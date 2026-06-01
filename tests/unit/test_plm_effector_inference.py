"""Regression tests for the ``models.*`` pickle-path shim in inference.py.

Upstream PLM-Effector trained its checkpoints with ``models.py`` as a
top-level import, so the saved .pth files pickle their classes as
``models.SimpleMLP`` etc. After vendoring, the same file lives at
``ssign_app.scripts.plm_effector.models``, which made torch.load raise
``ModuleNotFoundError: No module named 'models'`` until inference.py
started aliasing the vendored module into ``sys.modules['models']``.
These tests pin that alias in place.
"""

import sys

import torch

from ssign_app.scripts.plm_effector import inference
from ssign_app.scripts.plm_effector.models import SimpleMLP


def test_models_alias_registered_at_import():
    assert sys.modules.get("models") is sys.modules["ssign_app.scripts.plm_effector.models"]


def test_loadmodel_4test_handles_upstream_models_qualified_pickle(tmp_path, monkeypatch):
    # Forcing __module__ to "models" reproduces upstream's pickle layout.
    # Without the sys.modules alias, the load below would raise
    # ModuleNotFoundError("No module named 'models'").
    inference.clear_model_cache()
    monkeypatch.setattr(SimpleMLP, "__module__", "models")
    model = SimpleMLP(input_dim=4, hidden_dim=8, hidden_layer=1, dropout=0.0).eval()
    pth = tmp_path / "upstream.pth"
    torch.save(model, pth)

    x = torch.randn(2, 4)
    preds, probs = inference.loadmodel_4test(str(tmp_path), "upstream.pth", x, device="cpu")
    assert preds.shape == (2,)
    assert probs.shape == (2,)


class TestCheckpointCache:
    """Module-level cache for loaded .pth checkpoints.

    PLM-E's predict_api calls run_ensemble in a (chunk, effector_type)
    nested loop — 5 types × ~17 chunks = ~85 calls per genome — and each
    call would re-unpickle ~30 small checkpoint files. The cache makes
    each (model_dir, model_name, device) triple a one-shot disk read.
    """

    def _stub(self, tmp_path, name="ckpt.pth"):
        m = SimpleMLP(input_dim=4, hidden_dim=8, hidden_layer=1, dropout=0.0).eval()
        torch.save(m, tmp_path / name)

    def _patch_counting_load(self, monkeypatch):
        count = {"n": 0}
        original = torch.load

        def counting_load(*args, **kwargs):
            count["n"] += 1
            return original(*args, **kwargs)

        monkeypatch.setattr(torch, "load", counting_load)
        return count

    def test_repeated_loads_hit_cache(self, tmp_path, monkeypatch):
        inference.clear_model_cache()
        self._stub(tmp_path)
        count = self._patch_counting_load(monkeypatch)

        x = torch.randn(2, 4)
        for _ in range(5):
            inference.loadmodel_4predict(str(tmp_path), "ckpt.pth", x, "cpu")
        assert count["n"] == 1

    def test_loadmodel_4predict_and_4test_share_cache(self, tmp_path, monkeypatch):
        inference.clear_model_cache()
        self._stub(tmp_path)
        count = self._patch_counting_load(monkeypatch)

        x = torch.randn(2, 4)
        inference.loadmodel_4predict(str(tmp_path), "ckpt.pth", x, "cpu")
        inference.loadmodel_4test(str(tmp_path), "ckpt.pth", x, "cpu")
        # Same (dir, name, device) → one disk read total.
        assert count["n"] == 1

    def test_different_devices_keyed_separately(self, tmp_path, monkeypatch):
        inference.clear_model_cache()
        self._stub(tmp_path)
        count = self._patch_counting_load(monkeypatch)

        # Hit the cache helper directly: we only need to verify the key
        # discriminates on device, not actually run inference (a "meta"
        # tensor can't, but it's a valid map_location).
        inference._load_model_cached(str(tmp_path), "ckpt.pth", "cpu")
        inference._load_model_cached(str(tmp_path), "ckpt.pth", "meta")
        assert count["n"] == 2

    def test_different_names_keyed_separately(self, tmp_path, monkeypatch):
        inference.clear_model_cache()
        self._stub(tmp_path, name="a.pth")
        self._stub(tmp_path, name="b.pth")
        count = self._patch_counting_load(monkeypatch)

        x = torch.randn(2, 4)
        inference.loadmodel_4predict(str(tmp_path), "a.pth", x, "cpu")
        inference.loadmodel_4predict(str(tmp_path), "b.pth", x, "cpu")
        assert count["n"] == 2

    def test_clear_model_cache_drops_entries(self, tmp_path, monkeypatch):
        inference.clear_model_cache()
        self._stub(tmp_path)
        count = self._patch_counting_load(monkeypatch)

        x = torch.randn(2, 4)
        inference.loadmodel_4predict(str(tmp_path), "ckpt.pth", x, "cpu")
        inference.clear_model_cache()
        inference.loadmodel_4predict(str(tmp_path), "ckpt.pth", x, "cpu")
        assert count["n"] == 2

    def test_loaded_model_is_in_eval_mode(self, tmp_path):
        # Upstream pre-cache code only called .eval() inside the
        # DataParallel branch -- non-DP ckpts could retain train mode,
        # leaving Dropout active during inference. The cache now eval()s
        # unconditionally; pin that.
        inference.clear_model_cache()
        m = SimpleMLP(input_dim=4, hidden_dim=8, hidden_layer=1, dropout=0.5)
        m.train()  # save in train mode to reproduce the silent-bug condition
        torch.save(m, tmp_path / "ckpt.pth")

        loaded = inference._load_model_cached(str(tmp_path), "ckpt.pth", "cpu")
        assert loaded.training is False
