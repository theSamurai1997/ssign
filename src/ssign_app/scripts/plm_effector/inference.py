# Vendored from https://github.com/zhengdd0422/PLM-Effector's upstream
# `trainer.py` (CC-BY 3.0). Renamed to `inference.py` because the kept
# functions are all inference-time — there is no training code here.
# Upstream imports of `random`, `numpy`, and the large block of
# sklearn.metrics have been dropped; only `os` and `torch` remain --
# upstream's per-call `gc.collect()` + `torch.cuda.empty_cache()` were
# removed alongside the cache (loaded checkpoints must stay resident
# for the cache to be useful; see _MODEL_CACHE below).
import os
import sys

import torch
from torch.utils.data import DataLoader, TensorDataset

from . import models as _models

# Upstream's training script ran with `models.py` on the import path as a
# top-level module, so the pickled .pth files store class references as
# ``models.SimpleMLP`` etc. After vendoring, the same file lives at
# ``ssign_app.scripts.plm_effector.models``, so torch.load's unpickler
# can't resolve ``models`` and raises ModuleNotFoundError. Aliasing the
# vendored module into ``sys.modules`` makes the pickled qualified names
# resolve transparently. setdefault avoids clobbering a real ``models``
# package if one ever appears upstream.
sys.modules.setdefault("models", _models)


# Module-level cache of loaded checkpoints. PLM-E calls run_ensemble once
# per (effector_type, chunk) pair -- 85× per genome at chunk_size=256 --
# and each call loads ~30 small .pth files. With caching, each file is
# unpickled at most once per process. The cache holds the entire ensemble
# weight set (~750 MB at 5 effector types × 30 ckpts × ~5 MB each), which
# is small relative to the GPU's VRAM budget for PLM feature extraction.
# Keyed by (model_dir, model_name, str(device)) so the same checkpoint
# loaded on cpu vs cuda doesn't collide.
_MODEL_CACHE: "dict[tuple[str, str, str], torch.nn.Module]" = {}


def clear_model_cache() -> None:
    """Drop all cached checkpoint models. Test hygiene mostly; the process
    exit at the end of a PLM-E run also frees the cache."""
    _MODEL_CACHE.clear()


def _load_model_cached(model_dir: str, model_name: str, device) -> torch.nn.Module:
    """Return a checkpoint model, loading it the first time and reusing
    the in-memory copy on subsequent calls within the same process.

    FRAGILE: PyTorch 2.6 flipped ``torch.load``'s default from
    ``weights_only=False`` to ``weights_only=True``, which refuses to
    unpickle custom classes (e.g. models.SimpleMLP). PLM-Effector's
    upstream .pth files are whole-module saves (not state_dicts), so
    they require the legacy behaviour. The weights come from the
    upstream author URL via fetch_databases.sh — same trust boundary
    as every other .pth in the install.
    If this breaks: refactor upstream to save state_dicts instead, or
    use ``torch.serialization.add_safe_globals([SimpleMLP, ...])``.
    """
    key = (model_dir, model_name, str(device))
    cached = _MODEL_CACHE.get(key)
    if cached is not None:
        return cached
    model = torch.load(
        os.path.join(model_dir, model_name),
        map_location=device,
        weights_only=False,
    )
    if isinstance(model, torch.nn.DataParallel):
        model = model.module
    model.eval()
    _MODEL_CACHE[key] = model
    return model


def test_4predict_inbatch(model, features, device, batch_size=32):
    model.eval()
    dataset = TensorDataset(features)
    loader = DataLoader(dataset, batch_size=batch_size)

    all_preds = []
    all_probs = []

    with torch.no_grad():
        for batch in loader:
            x = batch[0].to(device)
            outputs = model(x).squeeze(dim=1)
            probs = torch.sigmoid(outputs)
            preds = (probs > 0.5).long()

            all_preds.append(preds.cpu())
            all_probs.append(probs.cpu())

    all_preds = torch.cat(all_preds).numpy()
    all_probs = torch.cat(all_probs).numpy()
    return all_preds, all_probs


def loadmodel_4predict(model_dir, model_name, x_test, device):
    model = _load_model_cached(model_dir, model_name, device)
    _, test_probs = test_4predict_inbatch(model, x_test, device, batch_size=128)
    return test_probs


def loadmodel_4test(model_dir, model_name, x_test, device):
    model = _load_model_cached(model_dir, model_name, device)
    test_preds, test_probs = test_4predict_inbatch(model, x_test, device, batch_size=128)
    return test_preds, test_probs
