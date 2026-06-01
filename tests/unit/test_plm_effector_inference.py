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
    monkeypatch.setattr(SimpleMLP, "__module__", "models")
    model = SimpleMLP(input_dim=4, hidden_dim=8, hidden_layer=1, dropout=0.0).eval()
    pth = tmp_path / "upstream.pth"
    torch.save(model, pth)

    x = torch.randn(2, 4)
    preds, probs = inference.loadmodel_4test(str(tmp_path), "upstream.pth", x, device="cpu")
    assert preds.shape == (2,)
    assert probs.shape == (2,)
