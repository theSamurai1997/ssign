"""Feature extraction via pretrained protein language models.

Refactored from upstream's `transformers_pretrainedmodel-features_extract_4predict.py`
(https://github.com/zhengdd0422/PLM-Effector, CC-BY 3.0 — see package
__init__.py). Upstream shipped the logic as a script with argparse and a
file-based `.npz` intermediate protocol; this module exposes a callable
function that returns in-memory numpy arrays instead. The heavy lifting
(tokenise → batch forward pass) is unchanged.
"""

from __future__ import annotations

import gc
import logging
import os

import numpy as np
import torch

from .utils import (
    batch_extract_features,
    read_fasta_for_prediction,
    read_fasta_for_prediction_terminal,
)

logger = logging.getLogger(__name__)

_MODEL_SUBDIRS = {
    "ProtBert": "prot_bert",
    "ProtT5": "prot_t5_xl_uniref50",
    "esm1": "esm1b_t33_650M_UR50S",
    "esm2_t33": "esm2_t33_650M_UR50D",
}

_MODEL_MAX_LENGTH = {"ProtBert": 512, "ProtT5": 512, "esm1": 1024, "esm2_t33": 1024}
_MODEL_REAL_SEQUENCE_LEN = {
    "ProtBert": 510,
    "ProtT5": 511,
    "esm1": 1022,
    "esm2_t33": 1022,
}


def _load_model_and_tokenizer(pretrained_type: str, weights_dir: str, device):
    """Load a HuggingFace pretrained PLM from `weights_dir/transformers_pretrained/<subdir>/`.

    `low_cpu_mem_usage=True` halves the CPU-RAM peak during load by streaming
    weights into a meta-device-allocated model instead of reading the full
    state_dict into a temp buffer first. Bit-identical to the default load
    (no accuracy change). Requires `accelerate` (in ssign[extended]). Critical
    for ProtT5 (~12 GB FP32) on memory-constrained nodes.
    """
    from transformers import AutoModel, AutoTokenizer, T5EncoderModel, T5Tokenizer

    subdir = _MODEL_SUBDIRS[pretrained_type]
    model_path = os.path.join(weights_dir, "transformers_pretrained", subdir)
    if not os.path.isdir(model_path):
        raise FileNotFoundError(
            f"Pretrained weights directory not found: {model_path}\n"
            f"  Expected PLM-Effector's HuggingFace weights at "
            f"{weights_dir}/transformers_pretrained/{subdir}/"
        )

    if pretrained_type == "ProtBert":
        tokenizer = AutoTokenizer.from_pretrained(model_path, do_lower_case=False)
        model = AutoModel.from_pretrained(model_path, low_cpu_mem_usage=True)
    elif pretrained_type == "ProtT5":
        tokenizer = T5Tokenizer.from_pretrained(model_path, do_lower_case=False)
        model = T5EncoderModel.from_pretrained(model_path, low_cpu_mem_usage=True)
    else:
        tokenizer = AutoTokenizer.from_pretrained(model_path)
        model = AutoModel.from_pretrained(model_path, low_cpu_mem_usage=True)

    model = model.to(device)
    model.eval()
    return model, tokenizer


def extract_terminal_features(
    proteins_fasta: str,
    pretrained_type: str,
    weights_dir: str,
    device,
    batch_size: int = 5,
) -> dict:
    """Extract N- and C-terminal features for every protein in the FASTA.

    Returns a dict keyed by terminal (`"Nterminal"` / `"Cterminal"`) with
    values `{"embedding": ndarray, "attention_masks": ndarray, "seq_ids": ndarray}`
    matching the shape that the upstream `.npz` files used.
    """
    model, tokenizer = _load_model_and_tokenizer(pretrained_type, weights_dir, device)
    max_length = _MODEL_MAX_LENGTH[pretrained_type]
    real_seq_len = _MODEL_REAL_SEQUENCE_LEN[pretrained_type]

    out: dict = {}
    try:
        for terminal in ("Nterminal", "Cterminal"):
            if terminal == "Cterminal":
                ids, sequences = read_fasta_for_prediction_terminal(
                    proteins_fasta,
                    model_type=pretrained_type,
                    terminal=terminal,
                    maxlen=real_seq_len,
                )
                tokenizer.padding_side = "left"
            else:
                ids, sequences = read_fasta_for_prediction(proteins_fasta, model_type=pretrained_type)
                tokenizer.padding_side = "right"

            embeddings, attention_masks = batch_extract_features(
                sequences,
                pretrained_type,
                model,
                tokenizer,
                device,
                max_length=max_length,
                batch_size=batch_size,
            )

            out[terminal] = {
                "embedding": embeddings,
                "attention_masks": attention_masks,
                "seq_ids": np.array(ids),
            }

            if device.type == "cuda":
                torch.cuda.empty_cache()
            gc.collect()
    finally:
        del model, tokenizer
        if device.type == "cuda":
            torch.cuda.empty_cache()
        gc.collect()

    return out


def _save_features_npz(out_path: str, features: dict) -> None:
    """Write a single PLM's two-terminal features to a .npz file.

    Keys are flattened (n_emb/n_mask/n_ids + c_emb/c_mask/c_ids) so np.load
    can read them without `allow_pickle=True`.
    """
    np.savez(
        out_path,
        n_emb=features["Nterminal"]["embedding"],
        n_mask=features["Nterminal"]["attention_masks"],
        n_ids=features["Nterminal"]["seq_ids"],
        c_emb=features["Cterminal"]["embedding"],
        c_mask=features["Cterminal"]["attention_masks"],
        c_ids=features["Cterminal"]["seq_ids"],
    )


def _load_features_npz(in_path: str) -> dict:
    """Inverse of _save_features_npz."""
    data = np.load(in_path, allow_pickle=False)
    return {
        "Nterminal": {
            "embedding": data["n_emb"],
            "attention_masks": data["n_mask"],
            "seq_ids": data["n_ids"],
        },
        "Cterminal": {
            "embedding": data["c_emb"],
            "attention_masks": data["c_mask"],
            "seq_ids": data["c_ids"],
        },
    }


def _extract_one_plm_in_subprocess(
    proteins_fasta: str,
    pretrained_type: str,
    weights_dir: str,
    device,
    batch_size: int,
    out_npz: str,
) -> None:
    """Run extract_terminal_features for one PLM in a fresh Python process.

    PyTorch's caching allocator holds onto RAM and VRAM across `del model
    + empty_cache + gc.collect`, so stacking 3-4 PLMs in one process
    pushes peak memory above 32 GB on Imperial CX3's GPU nodes (T1SE on
    ESM-1b OOM-killed even with `low_cpu_mem_usage=True`). Spawning a
    subprocess per PLM bounds peak memory to one model — when the
    subprocess exits, the OS reclaims every page.
    """
    import subprocess
    import sys

    cmd = [
        sys.executable,
        "-m",
        "ssign_app.scripts.plm_effector.feature_extraction",
        "--proteins-fasta",
        proteins_fasta,
        "--pretrained-type",
        pretrained_type,
        "--weights-dir",
        weights_dir,
        "--device",
        str(device),
        "--batch-size",
        str(batch_size),
        "--output",
        out_npz,
    ]
    logger.info(f"PLM-Effector: spawning isolated subprocess for {pretrained_type}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # Surface stderr so the runner sees the real failure, not just exit -9.
        raise RuntimeError(
            f"PLM extraction subprocess for {pretrained_type} exited "
            f"with code {result.returncode}\n--- stderr ---\n{result.stderr[-2000:]}"
        )


def extract_all_features(
    proteins_fasta: str,
    effector_type: str,
    weights_dir: str,
    device,
    batch_size: int = 5,
    isolate_plms: bool = True,
) -> dict:
    """Run feature extraction for every PLM that `effector_type` needs.

    With `isolate_plms=True` (default), each PLM runs in a fresh Python
    subprocess so VRAM and host RAM are fully reclaimed between models.
    This is the only reliable way to keep peak memory bounded to a single
    PLM on memory-constrained nodes. Set `isolate_plms=False` to stack
    in-process (slightly faster on a node with abundant RAM).

    Returns a nested dict: `features[pretrained_type][terminal] = {"embedding", "attention_masks", "seq_ids"}`.
    """
    if effector_type == "T4SE":
        pretrained_types = ["esm1", "esm2_t33", "ProtBert", "ProtT5"]
    else:
        pretrained_types = ["esm1", "esm2_t33", "ProtT5"]

    features: dict = {}
    if isolate_plms:
        import tempfile

        with tempfile.TemporaryDirectory(prefix="plm_features_") as tmp:
            for pretrained_type in pretrained_types:
                npz_path = os.path.join(tmp, f"{pretrained_type}.npz")
                _extract_one_plm_in_subprocess(
                    proteins_fasta=proteins_fasta,
                    pretrained_type=pretrained_type,
                    weights_dir=weights_dir,
                    device=device,
                    batch_size=batch_size,
                    out_npz=npz_path,
                )
                features[pretrained_type] = _load_features_npz(npz_path)
                os.remove(npz_path)  # free disk during the loop
    else:
        for pretrained_type in pretrained_types:
            features[pretrained_type] = extract_terminal_features(
                proteins_fasta=proteins_fasta,
                pretrained_type=pretrained_type,
                weights_dir=weights_dir,
                device=device,
                batch_size=batch_size,
            )
    return features


def _cli_main() -> int:
    """Subprocess entry point used by _extract_one_plm_in_subprocess.

    Loads one PLM, extracts features for both terminals, writes to .npz,
    exits. The OS reclaims all memory at exit — the whole point of this
    indirection. Not for direct end-user use; called via
    `python -m ssign_app.scripts.plm_effector.feature_extraction ...`.
    """
    import argparse

    p = argparse.ArgumentParser(description=_cli_main.__doc__)
    p.add_argument("--proteins-fasta", required=True)
    p.add_argument("--pretrained-type", required=True, choices=list(_MODEL_SUBDIRS))
    p.add_argument("--weights-dir", required=True)
    p.add_argument("--device", required=True, help='torch device string, e.g. "cuda" or "cpu"')
    p.add_argument("--batch-size", type=int, default=5)
    p.add_argument("--output", required=True, help="Output .npz path")
    args = p.parse_args()

    device = torch.device(args.device)
    # Bound torch intra-op threads to the cgroup-allocated CPUs. Default
    # is host total, which on a shared 128-core node would launch 128
    # OMP threads per PLM forward pass and trash the 4-CPU allocation.
    try:
        from ssign_app.scripts.ssign_lib.resources import effective_cpu_count

        torch.set_num_threads(effective_cpu_count())
    except Exception as e:
        logger.warning(f"Could not set torch thread count: {e}")

    features = extract_terminal_features(
        proteins_fasta=args.proteins_fasta,
        pretrained_type=args.pretrained_type,
        weights_dir=args.weights_dir,
        device=device,
        batch_size=args.batch_size,
    )
    _save_features_npz(args.output, features)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(_cli_main())
