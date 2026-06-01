"""Prediction-time utilities for PLM-Effector.

Vendored subset of the upstream `utils.py` (Zheng et al. 2026,
https://github.com/zhengdd0422/PLM-Effector). See `__init__.py` for
licence and attribution. Training-only helpers (WeightedSampler,
compute_class_weights, custom_collate, load_model_without_dataparallel,
load_val_numpy, load_test_numpy, load_test_numpy_nopool,
process_Cterminal_sequence) have been dropped because they are not
needed for inference. Upstream `torch_geometric`, `scipy.sparse`, and
`Bio.PDB` imports were vestigial in the prediction path and have been
removed.
"""

from __future__ import annotations

import random
import re

import numpy as np

# torch is imported lazily inside the heavy helpers so the pure-Python
# FASTA parsers can be unit-tested on a minimal dev environment that
# doesn't have torch installed.


_BERT_LIKE = {"Bert", "BioBERT", "ProtBert", "ProtT5"}


def _normalise_sequence(seq: str, model_type: str) -> str:
    """Prepare a single protein sequence for a pretrained tokenizer.

    BERT-family tokenisers (ProtBert, ProtT5) expect space-separated
    residues; ESM tokenisers expect the raw string. All models want the
    rare residues U/Z/O/B replaced by X.
    """
    spaced = " ".join(seq) if model_type in _BERT_LIKE else seq
    return re.sub(r"[UZOB]", "X", spaced)


def read_fasta_for_prediction(fasta_path: str, model_type: str = "esm1"):
    """Return `(ids, sequences)` parsed from a FASTA file.

    IDs include the leading `>` as emitted by the upstream parser, which
    callers strip. Residues are normalised for the chosen pretrained
    tokenizer.
    """
    ids: list[str] = []
    sequences: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                ids.append(line)
            else:
                sequences.append(_normalise_sequence(line, model_type))
    return ids, sequences


def read_fasta_for_prediction_terminal(
    fasta_path: str,
    model_type: str = "esm1",
    terminal: str = "Cterminal",
    maxlen: int = 1022,
):
    """Like `read_fasta_for_prediction` but truncated to N- or C-terminal region.

    Long proteins get cut to the first `maxlen` residues (N-terminal) or the
    last `maxlen` residues (C-terminal). Short proteins pass through
    unchanged. PLM-Effector trains one model per terminal.
    """
    ids: list[str] = []
    sequences: list[str] = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                ids.append(line)
            else:
                truncated = line[-maxlen:] if terminal == "Cterminal" else line[:maxlen]
                sequences.append(_normalise_sequence(truncated, model_type))
    return ids, sequences


_DTYPE_ALIASES = {
    "fp32": None,
    "float32": None,
    "none": None,
    "": None,
    "bf16": "bfloat16",
    "bfloat16": "bfloat16",
    "fp16": "float16",
    "float16": "float16",
    "half": "float16",
}


def resolve_autocast_dtype(name):
    """Map a user-facing dtype name to the torch.dtype expected by autocast.

    Returns ``None`` to signal "no autocast" (default for fp32 / CPU paths).
    Accepts the same name as a torch.dtype, a string alias, or None.
    Unknown names raise ValueError so a typo at the CLI doesn't silently
    fall through to fp32.
    """
    import torch  # lazy

    if name is None:
        return None
    if isinstance(name, torch.dtype):
        return name
    key = str(name).strip().lower()
    if key not in _DTYPE_ALIASES:
        raise ValueError(
            f"Unknown PLM-Effector dtype: {name!r}. Expected one of: fp32 / bf16 / fp16 (or their long names)."
        )
    resolved = _DTYPE_ALIASES[key]
    if resolved is None:
        return None
    return getattr(torch, resolved)


def batch_extract_features(
    sequences,
    pretrained_type: str,
    model,
    tokenizer,
    device,
    max_length: int = 512,
    batch_size: int = 10,
    autocast_dtype=None,
):
    """Run a batched forward pass through a pretrained PLM.

    Returns `(features, attention_masks)` as numpy arrays shaped
    `(n_sequences, max_length, hidden_dim)` and
    `(n_sequences, max_length)` respectively.

    ``autocast_dtype`` lets the caller run the forward pass under
    ``torch.autocast`` (bfloat16 or float16). None (default) keeps the
    model in its loaded precision (fp32). bfloat16 is the safe default
    for inference on A40/A100/L40S because it preserves fp32's exponent
    range; fp16 trades extra speed on V100/T4 for narrower range and is
    occasionally seen to shift predictions on softmax-heavy paths.
    """
    import contextlib

    import torch  # Lazy — ssign's base install doesn't require torch

    features: list = []
    attention_masks: list = []

    use_autocast = autocast_dtype is not None and device.type == "cuda"

    def _amp_ctx():
        # Fresh autocast manager per batch — torch.autocast tracks state in
        # __enter__/__exit__, and re-using one instance across iterations is
        # technically allowed but unidiomatic enough to confuse a future reader.
        if use_autocast:
            return torch.autocast(device_type="cuda", dtype=autocast_dtype)
        return contextlib.nullcontext()

    for i in range(0, len(sequences), batch_size):
        batch = sequences[i : i + batch_size]
        inputs = tokenizer(
            batch,
            return_tensors="pt",
            padding="max_length",
            truncation=True,
            max_length=max_length,
        )
        inputs = {k: v.to(device) for k, v in inputs.items()}
        # ProtT5's forward() rejects extra tokenizer kwargs (e.g. token_type_ids)
        # so we pass only the two it accepts; the other PLMs are happy with **inputs.
        if pretrained_type == "ProtT5":
            fwd_kwargs = {"input_ids": inputs["input_ids"], "attention_mask": inputs["attention_mask"]}
        else:
            fwd_kwargs = inputs
        with torch.no_grad(), _amp_ctx():
            outputs = model(**fwd_kwargs)
        # numpy lacks a bfloat16 dtype, so always cast back to fp32 before
        # leaving the GPU. The cast is free relative to the forward pass.
        features.append(outputs.last_hidden_state.float().cpu())
        attention_masks.append(inputs["attention_mask"].cpu())
        del inputs, outputs
        if device.type == "cuda":
            torch.cuda.empty_cache()

    all_features = torch.cat(features, dim=0).numpy()
    all_attention_masks = torch.cat(attention_masks, dim=0).numpy()
    return all_features, all_attention_masks


def pool_features(features, attention_masks, pooling: str = "mean"):
    """Pool per-residue features into a single per-sequence vector.

    Supports `"mean"` (default) and `"max"` pooling, masking padding
    positions via the attention mask. Returns a `torch.Tensor`.
    """
    import torch  # Lazy — ssign's base install doesn't require torch

    if not torch.is_tensor(features):
        features = torch.tensor(features, dtype=torch.float32)
    if not torch.is_tensor(attention_masks):
        attention_masks = torch.tensor(attention_masks, dtype=torch.float32)

    attention_masks = attention_masks.unsqueeze(-1)

    if pooling == "mean":
        return (features * attention_masks).sum(dim=1) / attention_masks.sum(dim=1)
    if pooling == "max":
        return (features * attention_masks).max(dim=1).values
    raise ValueError(f"Unsupported pooling type {pooling!r}; use 'mean' or 'max'.")


def set_seed(seed: int) -> None:
    """Seed Python, NumPy, and PyTorch RNGs for reproducible inference."""
    import torch  # Lazy — ssign's base install doesn't require torch

    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
