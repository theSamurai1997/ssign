"""Top-level PLM-Effector prediction API.

Wraps `feature_extraction.extract_all_features` and `ensemble.run_ensemble`
into a single `predict()` call covering one or more effector types.
Writes one TSV per type matching the upstream output format plus an
explicit `passes_threshold` column so downstream ssign code does not
need to know each effector type's threshold.

Memory: features are extracted as per-chunk .npz files (default 256
proteins/chunk). The ensemble runs once per chunk per type, predictions
are accumulated, then concatenated and written as a TSV — so peak RAM
is bounded by one chunk's features + one chunk's predictions regardless
of input size.

Multi-type runs: PLM embeddings (ESM-1b / ESM-2 / ProtT5 / ProtBert) are
protein-intrinsic — the per-effector neural-net heads and XGBoost
stackers are the only type-specific bits. So when called with multiple
effector types, `predict` extracts features once and reuses them across
all types, saving ~75% wallclock on the K-12 tutorial.
"""

from __future__ import annotations

import logging
import os
import tempfile

import numpy as np

# torch / ensemble / feature_extraction are imported lazily inside the
# heavy functions so `write_predictions_tsv` (pure-numpy, unit-testable)
# can be imported without the ML deps installed.


logger = logging.getLogger(__name__)

_VALID_EFFECTOR_TYPES = ("T1SE", "T2SE", "T3SE", "T4SE", "T6SE")


def _resolve_device(device):
    """Accept None / "cuda" / "cpu" / torch.device and return a torch.device."""
    import torch

    if isinstance(device, torch.device):
        return device
    if device is None or device == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError(
                "PLM-Effector requires a CUDA GPU (default) but torch reports "
                "no CUDA device available. Set device='cpu' to force CPU "
                "inference (very slow, ~1-2 min per protein)."
            )
        return torch.device("cuda")
    if device == "cpu":
        return torch.device("cpu")
    raise ValueError(f"Unrecognised device: {device!r}")


def predict(
    proteins_fasta: str,
    weights_dir: str,
    effector_types,
    out_dir: str,
    device=None,
    batch_size: int = 5,
    chunk_size: int = 256,
    dtype: str = "bf16",
) -> dict:
    """Run PLM-Effector on a protein FASTA for one or more effector types.

    Writes `{out_dir}/{effector_type}.tsv` per type. Feature extraction
    runs once over the union of PLMs needed by all requested types.

    Args:
        proteins_fasta: input protein FASTA.
        weights_dir: directory containing `transformers_pretrained/<model>/`
            subdirs (HuggingFace PLMs) and `trained_models/<effector_type>_*.pth`
            + `<effector_type>_XGB_stackingmeta_model.json` per type.
        effector_types: a single type string or an iterable of types from
            T1SE / T2SE / T3SE / T4SE / T6SE.
        out_dir: directory to write per-type prediction TSVs into.
        device: None (auto-pick cuda), "cuda", "cpu", or a torch.device.
            CPU works but is very slow.
        batch_size: PLM forward-pass batch size. Drop to 1-2 if VRAM-constrained.
        chunk_size: proteins per feature-extraction chunk. Default 256 keeps
            peak host RAM ~10-17 GB per PLM (FP32).

    Returns:
        `{effector_type: (n_positive, out_path)}` per requested type.
    """
    if isinstance(effector_types, str):
        effector_types = [effector_types]
    effector_types = list(effector_types)
    bad = [t for t in effector_types if t not in _VALID_EFFECTOR_TYPES]
    if bad:
        raise ValueError(f"unrecognised effector_types {bad!r}; must be from {_VALID_EFFECTOR_TYPES}")
    if not effector_types:
        raise ValueError("effector_types is empty; nothing to predict")
    if not os.path.exists(proteins_fasta):
        raise FileNotFoundError(f"Input FASTA not found: {proteins_fasta}")
    if not os.path.isdir(weights_dir):
        raise FileNotFoundError(
            f"PLM-Effector weights directory not found: {weights_dir}\n"
            f"  Fetch with scripts/fetch_weights.sh or mirror manually from "
            f"http://www.mgc.ac.cn/PLM-Effector/downloads.html"
        )
    os.makedirs(out_dir, exist_ok=True)

    from .ensemble import run_ensemble
    from .feature_extraction import extract_all_features, iter_chunk_features, pretrained_types_for

    torch_device = _resolve_device(device)
    pretrained_types = pretrained_types_for(effector_types)
    logger.info(
        "PLM-Effector: extracting features for %s on %s (batch_size=%d, chunk_size=%d, dtype=%s, PLMs=%s)",
        ",".join(effector_types),
        torch_device,
        batch_size,
        chunk_size,
        dtype,
        ",".join(pretrained_types),
    )

    results = {t: {"ids": [], "stacked": [], "probs": [], "passes": []} for t in effector_types}
    with tempfile.TemporaryDirectory(prefix="plm_features_") as feature_dir:
        chunk_paths = extract_all_features(
            proteins_fasta=proteins_fasta,
            pretrained_types=pretrained_types,
            weights_dir=weights_dir,
            device=torch_device,
            feature_cache_dir=feature_dir,
            batch_size=batch_size,
            chunk_size=chunk_size,
            dtype=dtype,
        )
        n_chunks = len(next(iter(chunk_paths.values())))
        logger.info(
            "PLM-Effector: running ensemble for %d type(s) across %d chunk(s)",
            len(effector_types),
            n_chunks,
        )

        for chunk_features in iter_chunk_features(chunk_paths, delete_after_yield=True):
            for eff_type in effector_types:
                ids, stacked, final_probs, passes = run_ensemble(
                    features=chunk_features,
                    weights_dir=weights_dir,
                    effector_type=eff_type,
                    device=torch_device,
                )
                results[eff_type]["ids"].append(ids)
                results[eff_type]["stacked"].append(stacked)
                results[eff_type]["probs"].append(final_probs)
                results[eff_type]["passes"].append(passes)

    summary: dict = {}
    for eff_type, parts in results.items():
        seq_ids = np.concatenate(parts["ids"])
        stacked = np.concatenate(parts["stacked"], axis=0)
        final_probs = np.concatenate(parts["probs"])
        passes = np.concatenate(parts["passes"])
        n_base = stacked.shape[1]
        out_path = os.path.join(out_dir, f"{eff_type}.tsv")
        write_predictions_tsv(out_path, seq_ids, stacked, final_probs, passes, eff_type, n_base)
        n_positive = int(passes.sum())
        logger.info(
            "PLM-Effector: %s — wrote %d predictions (%d passing threshold) to %s",
            eff_type,
            len(seq_ids),
            n_positive,
            out_path,
        )
        summary[eff_type] = (n_positive, out_path)
    return summary


def write_predictions_tsv(
    out_path: str,
    seq_ids: np.ndarray,
    stacked: np.ndarray,
    final_probs: np.ndarray,
    passes: np.ndarray,
    effector_type: str,
    n_base_models: int,
) -> None:
    """Write predictions as TSV: `seq_id\\tmodel1..N\\tstacking\\tpasses_threshold\\teffector_type`.

    One row per input protein. Upstream's format only emitted rows above
    threshold; we emit all proteins and let downstream code filter so the
    output shape is predictable regardless of positive-hit count.
    """
    headers = (
        ["seq_id"] + [f"model{i + 1}" for i in range(n_base_models)] + ["stacking", "passes_threshold", "effector_type"]
    )
    with open(out_path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for i, sid in enumerate(seq_ids):
            sid_str = sid.decode() if isinstance(sid, bytes) else str(sid)
            # Strip leading '>' if present (upstream parser kept it)
            if sid_str.startswith(">"):
                sid_str = sid_str[1:].split()[0]
            row = (
                [sid_str]
                + [f"{stacked[i, j]:.4f}" for j in range(n_base_models)]
                + [
                    f"{final_probs[i]:.4f}",
                    "1" if passes[i] else "0",
                    effector_type,
                ]
            )
            f.write("\t".join(row) + "\n")
