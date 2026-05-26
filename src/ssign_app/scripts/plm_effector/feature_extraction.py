"""Feature extraction via pretrained protein language models.

Refactored from upstream's `transformers_pretrainedmodel-features_extract_4predict.py`
(https://github.com/zhengdd0422/PLM-Effector, CC-BY 3.0 — see package
__init__.py). Upstream shipped the logic as a script with argparse and a
file-based `.npz` intermediate protocol; this module exposes a callable
function that returns chunked numpy arrays on disk instead.

Chunking: a genome-scale FASTA (~4k proteins) accumulates ~22 GB of
per-protein embeddings per terminal per PLM at FP32, blowing a 32 GB
cgroup. We process proteins in chunks (default 256), write one .npz per
(PLM, chunk), and let the ensemble consume chunks one at a time.
Per-protein outputs are bit-identical to the all-at-once path because
padding is `max_length` (chunk size doesn't affect per-sequence tokens),
the PLMs have no BatchNorm and no cross-protein attention, and the
tokenizer is per-sequence.
"""

from __future__ import annotations

import gc
import glob
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

DEFAULT_CHUNK_SIZE = 256


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


def _read_both_terminals(proteins_fasta: str, pretrained_type: str):
    """Return `{"Nterminal": (ids, seqs), "Cterminal": (ids, seqs)}`."""
    real_seq_len = _MODEL_REAL_SEQUENCE_LEN[pretrained_type]
    n_ids, n_seqs = read_fasta_for_prediction(proteins_fasta, model_type=pretrained_type)
    c_ids, c_seqs = read_fasta_for_prediction_terminal(
        proteins_fasta,
        model_type=pretrained_type,
        terminal="Cterminal",
        maxlen=real_seq_len,
    )
    return {"Nterminal": (n_ids, n_seqs), "Cterminal": (c_ids, c_seqs)}


def iter_terminal_feature_chunks(
    proteins_fasta: str,
    pretrained_type: str,
    weights_dir: str,
    device,
    batch_size: int = 5,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
):
    """Yield per-chunk N+C terminal features for every protein in the FASTA.

    Each yield is a dict shaped like `extract_terminal_features`' return
    but covering only `chunk_size` proteins. The model is loaded once and
    reused across chunks (loading cost amortises). Padding is
    `max_length`, so per-protein outputs are bit-identical to processing
    all proteins at once.
    """
    # Read FASTA first so an empty input skips the model load (5-30 s + 2-4 GB).
    terminals = _read_both_terminals(proteins_fasta, pretrained_type)
    n_proteins = len(terminals["Nterminal"][0])
    if n_proteins == 0:
        return

    model, tokenizer = _load_model_and_tokenizer(pretrained_type, weights_dir, device)
    max_length = _MODEL_MAX_LENGTH[pretrained_type]

    try:
        for chunk_start in range(0, n_proteins, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_proteins)
            chunk_out: dict = {}

            for terminal, (ids, seqs) in terminals.items():
                # padding_side flips between terminals: C-terminal pads on
                # the left so the C-terminal residues sit at the right edge
                # (last attention positions), N-terminal pads right.
                tokenizer.padding_side = "left" if terminal == "Cterminal" else "right"

                chunk_seqs = seqs[chunk_start:chunk_end]
                chunk_ids = ids[chunk_start:chunk_end]

                embeddings, attention_masks = batch_extract_features(
                    chunk_seqs,
                    pretrained_type,
                    model,
                    tokenizer,
                    device,
                    max_length=max_length,
                    batch_size=batch_size,
                )

                chunk_out[terminal] = {
                    "embedding": embeddings,
                    "attention_masks": attention_masks,
                    "seq_ids": np.array(chunk_ids),
                }

            yield chunk_out

            if device.type == "cuda":
                torch.cuda.empty_cache()
            gc.collect()
    finally:
        del model, tokenizer
        if device.type == "cuda":
            torch.cuda.empty_cache()
        gc.collect()


def extract_terminal_features(
    proteins_fasta: str,
    pretrained_type: str,
    weights_dir: str,
    device,
    batch_size: int = 5,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict:
    """All-in-one terminal feature extraction (concatenates all chunks).

    Convenience wrapper around `iter_terminal_feature_chunks`. Allocates
    the full output in RAM — only safe for small inputs or test fixtures;
    production code should use the chunked path via `extract_all_features`.
    """
    chunks = list(
        iter_terminal_feature_chunks(
            proteins_fasta,
            pretrained_type,
            weights_dir,
            device,
            batch_size=batch_size,
            chunk_size=chunk_size,
        )
    )
    if not chunks:
        return {}
    if len(chunks) == 1:
        return chunks[0]
    out: dict = {}
    for terminal in ("Nterminal", "Cterminal"):
        out[terminal] = {
            "embedding": np.concatenate([c[terminal]["embedding"] for c in chunks], axis=0),
            "attention_masks": np.concatenate([c[terminal]["attention_masks"] for c in chunks], axis=0),
            "seq_ids": np.concatenate([c[terminal]["seq_ids"] for c in chunks], axis=0),
        }
    return out


def _save_features_npz(out_path: str, features: dict) -> None:
    """Write a single PLM's two-terminal features (one chunk's worth) to a .npz file.

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


def _chunk_path(out_dir: str, pretrained_type: str, chunk_idx: int) -> str:
    return os.path.join(out_dir, f"{pretrained_type}_chunk{chunk_idx:04d}.npz")


def _discover_chunk_paths(out_dir: str, pretrained_type: str) -> list[str]:
    return sorted(glob.glob(os.path.join(out_dir, f"{pretrained_type}_chunk*.npz")))


def iter_chunk_features(chunk_paths: dict, delete_after_yield: bool = True):
    """Yield per-chunk features dicts loaded from disk, optionally deleting after.

    Given the `{pretrained_type: [chunk_path_0, ...]}` mapping returned by
    `extract_all_features`, yields one `{pretrained_type: features_dict}`
    per chunk index. Each chunk's .npz files are deleted after yield to
    keep peak disk usage bounded.
    """
    if not chunk_paths:
        return
    n_chunks = len(next(iter(chunk_paths.values())))
    for chunk_idx in range(n_chunks):
        yield {pt: _load_features_npz(paths[chunk_idx]) for pt, paths in chunk_paths.items()}
        if delete_after_yield:
            for paths in chunk_paths.values():
                try:
                    os.remove(paths[chunk_idx])
                except OSError:
                    pass


def _extract_one_plm_in_subprocess(
    proteins_fasta: str,
    pretrained_type: str,
    weights_dir: str,
    device,
    batch_size: int,
    chunk_size: int,
    out_dir: str,
) -> list[str]:
    """Run iter_terminal_feature_chunks for one PLM in a fresh Python process.

    PyTorch's caching allocator holds onto RAM and VRAM across `del model
    + empty_cache + gc.collect`, so stacking 3-4 PLMs in one process
    pushes peak memory above 32 GB on Imperial CX3's GPU nodes. Spawning a
    subprocess per PLM bounds peak memory to one model — when the
    subprocess exits, the OS reclaims every page.

    Writes one `<pretrained_type>_chunkNNNN.npz` per chunk into `out_dir`.
    Returns the sorted list of chunk paths written.
    """
    import subprocess
    import sys
    from collections import deque

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
        "--chunk-size",
        str(chunk_size),
        "--output-dir",
        out_dir,
    ]
    logger.info(
        "PLM-Effector: spawning isolated subprocess for %s (chunk_size=%d)",
        pretrained_type,
        chunk_size,
    )

    # Stream stderr line-by-line to our stderr so the user sees per-chunk
    # progress instead of a multi-minute silence. Keep a rolling tail so we
    # can include the last few hundred lines in the failure message rather
    # than ask the user to scrollback.
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    tail: deque[str] = deque(maxlen=200)
    assert proc.stderr is not None
    for line in proc.stderr:
        sys.stderr.write(line)
        sys.stderr.flush()
        tail.append(line)
    returncode = proc.wait()
    if returncode != 0:
        raise RuntimeError(
            f"PLM extraction subprocess for {pretrained_type} exited "
            f"with code {returncode}\n--- stderr (tail) ---\n{''.join(tail)}"
        )

    chunk_paths = _discover_chunk_paths(out_dir, pretrained_type)
    if not chunk_paths:
        raise RuntimeError(
            f"PLM extraction subprocess for {pretrained_type} exited successfully but wrote no chunk files to {out_dir}"
        )
    return chunk_paths


def extract_all_features(
    proteins_fasta: str,
    effector_type: str,
    weights_dir: str,
    device,
    feature_cache_dir: str,
    batch_size: int = 5,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    isolate_plms: bool = True,
) -> dict:
    """Run feature extraction for every PLM that `effector_type` needs.

    With `isolate_plms=True` (default), each PLM runs in a fresh Python
    subprocess so VRAM and host RAM are fully reclaimed between models.
    Within each subprocess, proteins are processed in chunks of
    `chunk_size` so the all-protein embedding tensor never sits in RAM at
    once (~22 GB → ~1.3 GB per chunk for ESM-1b @ FP32).

    Args:
        feature_cache_dir: where to write per-chunk `<plm>_chunkNNNN.npz`
            files. Caller owns this directory; predict_api uses a tempdir.

    Returns:
        `{pretrained_type: [chunk_path_0, chunk_path_1, ...]}`. All PLMs
        produce the same number of chunks (proteins are split identically),
        and the chunk_idx ordering aligns across PLMs.
    """
    if effector_type == "T4SE":
        pretrained_types = ["esm1", "esm2_t33", "ProtBert", "ProtT5"]
    else:
        pretrained_types = ["esm1", "esm2_t33", "ProtT5"]

    chunk_paths: dict = {}
    for pretrained_type in pretrained_types:
        if isolate_plms:
            chunk_paths[pretrained_type] = _extract_one_plm_in_subprocess(
                proteins_fasta=proteins_fasta,
                pretrained_type=pretrained_type,
                weights_dir=weights_dir,
                device=device,
                batch_size=batch_size,
                chunk_size=chunk_size,
                out_dir=feature_cache_dir,
            )
        else:
            paths = []
            for chunk_idx, chunk_features in enumerate(
                iter_terminal_feature_chunks(
                    proteins_fasta=proteins_fasta,
                    pretrained_type=pretrained_type,
                    weights_dir=weights_dir,
                    device=device,
                    batch_size=batch_size,
                    chunk_size=chunk_size,
                )
            ):
                path = _chunk_path(feature_cache_dir, pretrained_type, chunk_idx)
                _save_features_npz(path, chunk_features)
                paths.append(path)
            chunk_paths[pretrained_type] = paths
    return chunk_paths


def _cli_main() -> int:
    """Subprocess entry point used by _extract_one_plm_in_subprocess.

    Loads one PLM, iterates over protein chunks, writes one .npz per
    chunk into --output-dir, exits. The OS reclaims all memory at exit
    — the whole point of this indirection. Not for direct end-user use;
    called via `python -m ssign_app.scripts.plm_effector.feature_extraction ...`.
    """
    import argparse

    p = argparse.ArgumentParser(description=_cli_main.__doc__)
    p.add_argument("--proteins-fasta", required=True)
    p.add_argument("--pretrained-type", required=True, choices=list(_MODEL_SUBDIRS))
    p.add_argument("--weights-dir", required=True)
    p.add_argument("--device", required=True, help='torch device string, e.g. "cuda" or "cpu"')
    p.add_argument("--batch-size", type=int, default=5)
    p.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    p.add_argument("--output-dir", required=True, help="Directory to write chunk .npz files into")
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

    os.makedirs(args.output_dir, exist_ok=True)
    n_written = 0
    for chunk_idx, chunk_features in enumerate(
        iter_terminal_feature_chunks(
            proteins_fasta=args.proteins_fasta,
            pretrained_type=args.pretrained_type,
            weights_dir=args.weights_dir,
            device=device,
            batch_size=args.batch_size,
            chunk_size=args.chunk_size,
        )
    ):
        path = _chunk_path(args.output_dir, args.pretrained_type, chunk_idx)
        _save_features_npz(path, chunk_features)
        n_written += 1
    logger.info(
        "%s: wrote %d chunk(s) (chunk_size=%d) to %s",
        args.pretrained_type,
        n_written,
        args.chunk_size,
        args.output_dir,
    )
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(_cli_main())
