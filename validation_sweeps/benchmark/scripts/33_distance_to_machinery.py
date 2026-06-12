#!/usr/bin/env python3
"""Phase 2: how far is each effector from its secretion machinery? (why proximity ±3 misses most.)

For every testable effector, the gene-order distance to the nearest anchored answer-key machinery
gene (from ceiling_per_effector). A strip/box per SS type on a log axis, with the ±3 proximity line,
shows directly that T2/T3/T4/T6 effectors sit a median of hundreds of genes from their apparatus
(genuine genome-dispersal), while T1SS clusters at distance 1 (operonic) with a few far exceptions.

Inputs : data/phase1/ceiling_per_effector.tsv
Output : data/phase2/figures/summary/04_distance_to_machinery.png
Run:     <repo>/.venv/bin/python scripts/33_distance_to_machinery.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

BENCH = Path(__file__).resolve().parents[1]
FIGDIR = BENCH / "data" / "phase2" / "figures" / "summary"
TYPES = ["T1SS", "T2SS", "T3SS", "T4SS", "T6SS"]
plt.rcParams.update(
    {
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.titlepad": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#444444",
    }
)


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    rows = [
        r
        for r in csv.DictReader(open(BENCH / "data" / "phase1" / "ceiling_per_effector.tsv"), delimiter="\t")
        if r.get("testable") == "yes"
    ]
    dists = {ss: [] for ss in TYPES}
    for r in rows:
        d = r.get("nearest_dist", "").strip()
        if d.isdigit():
            dists[r["ss_type"]].append(max(1, int(d)))  # clamp 0->1 for the log axis (adjacency)

    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    rng = np.random.RandomState(0)  # noqa: NPY002 — fixed seed, jitter only
    for i, ss in enumerate(TYPES):
        ys = np.array(dists[ss])
        if len(ys) == 0:
            continue
        x = i + (rng.rand(len(ys)) - 0.5) * 0.5
        ax.scatter(x, ys, s=14, alpha=0.45, color="#3F6F8C", edgecolors="none")
        med = float(np.median(ys))
        ax.plot([i - 0.32, i + 0.32], [med, med], color="#A93232", lw=2.5, zorder=5)
        ax.text(
            i,
            med * 1.25,
            f"median {med:.0f}",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="#A93232",
            fontweight="bold",
        )
    ax.axhline(3, color="#2E7D6F", lw=1.8, ls="--")
    ax.text(
        len(TYPES) - 0.5,
        3.4,
        "proximity reach (±3)",
        ha="right",
        va="bottom",
        color="#2E7D6F",
        fontsize=9,
        fontweight="bold",
    )
    ax.set_yscale("log")
    ax.set_xticks(range(len(TYPES)))
    ax.set_xticklabels([f"{ss}\n(n={len(dists[ss])})" for ss in TYPES])
    ax.set_ylabel("genes from effector to nearest machinery (log)")
    ax.set_title(
        "Why proximity ±3 misses most effectors: how far they sit from their apparatus\nT1SS clusters at adjacency; T2/T3/T4/T6 effectors are genome-dispersed (median 100s of genes)"
    )
    fig.tight_layout()
    fig.savefig(FIGDIR / "04_distance_to_machinery.png")
    plt.close(fig)
    print("wrote 04_distance_to_machinery.png")
    for ss in TYPES:
        ys = dists[ss]
        if ys:
            within3 = sum(1 for d in ys if d <= 3)
            print(f"  {ss}: n={len(ys)} median={np.median(ys):.0f} within±3={within3}")


if __name__ == "__main__":
    main()
