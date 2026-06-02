#!/usr/bin/env python3
"""Sample N random non-SS-neighborhood proteins for enrichment background.

Reuses extract_neighborhood's helpers to identify the SS-component
neighborhood, then samples N proteins from the complement to build the
null background distribution. Used by --enrichment-stats: DLP/DSE run
on these null proteins in addition to the neighborhood, and the
resulting positive rates become p_DLP / p_DSE in the per-system
binomial test (enrichment_testing.py).

Determinism: seeded via --seed (default 42). Re-running with the same
seed, proteome, and ss_components yields the same null set.
"""

import argparse
import logging
import os as _os
import random
import sys as _sys

_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)

from extract_neighborhood import (  # noqa: E402
    get_neighborhood_proteins,
    load_gene_order,
    load_ss_components,
)
from ssign_lib.fasta_io import read_fasta, write_fasta  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def sample_null(all_ids, exclude, n, rng):
    """Pick up to n IDs from ``all_ids - exclude``. Returns a sorted list.

    Sorting the candidate pool before sampling makes the output
    deterministic given a seed -- dict iteration order isn't guaranteed
    to be stable across runs even with PYTHONHASHSEED, but a sorted list
    plus a seeded RNG is.
    """
    candidates = sorted(set(all_ids) - set(exclude))
    if not candidates:
        return []
    if n >= len(candidates):
        return candidates
    return sorted(rng.sample(candidates, n))


def main():
    parser = argparse.ArgumentParser(description="Sample non-SS-neighborhood proteins")
    parser.add_argument("--proteins", required=True, help="Full proteome FASTA")
    parser.add_argument("--gene-order", required=True)
    parser.add_argument("--ss-components", required=True)
    parser.add_argument("--window", type=int, default=3)
    parser.add_argument("--n", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-ids", required=True)
    args = parser.parse_args()

    proteome = read_fasta(args.proteins)
    if not proteome:
        logger.warning("Empty proteome -- writing empty null sample")
        open(args.out_fasta, "w").close()
        open(args.out_ids, "w").close()
        return

    ss_components = load_ss_components(args.ss_components)
    gene_order = load_gene_order(args.gene_order)
    neighborhood = get_neighborhood_proteins(gene_order, ss_components, args.window)

    rng = random.Random(args.seed)
    picked = sample_null(proteome.keys(), neighborhood, args.n, rng)

    write_fasta({pid: proteome[pid] for pid in picked}, args.out_fasta)
    with open(args.out_ids, "w") as f:
        for pid in picked:
            f.write(pid + "\n")

    pool = len(proteome) - len(neighborhood)
    logger.info(
        "Sampled %d null proteins (requested %d; pool %d non-neighborhood of %d total; %d in neighborhood; seed=%d)",
        len(picked),
        args.n,
        pool,
        len(proteome),
        len(neighborhood),
        args.seed,
    )


if __name__ == "__main__":
    main()
