#!/usr/bin/env python3
"""
21_ceiling_figures.py  (Phase 1 task 5.5: ceiling figures)

Renders the Phase 1 ceiling result as four numbered figures:
  01  gold-set composition per SS type (reachable@3 / impossible@3 / untestable counts)
  02  ceiling vs proximity window (N = 3,5,7) per SS type  -- the headline
  03  effector -> nearest own-machinery distance ECDF per SS type, with N reference lines
  04  per-genome ceiling@5 for genomes with >=3 testable effectors (drill-down)

Input : data/phase1/ceiling_per_effector.tsv
Output: figures/01_*.png .. 04_*.png
Run:    .venv/bin/python scripts/21_ceiling_figures.py
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BENCH = Path(__file__).resolve().parents[1]
IN = BENCH / "data" / "phase1" / "ceiling_per_effector.tsv"
FIGDIR = BENCH / "figures"
NS = (3, 5, 7)
ORDER = ["T1SS", "T2SS", "T3SS", "T4SS", "T6SS"]

THEME = {
    "ss_colors": {
        "T1SS": "#3F8E8C",
        "T2SS": "#E0884B",
        "T3SS": "#6C8EAD",
        "T4SS": "#A93232",
        "T6SS": "#8E6CA3",
    },
    "reachable": "#3F8E8C",
    "impossible": "#C24A4A",
    "untestable": "#BFBFBF",
    "ref_line": "#444444",
}

plt.rcParams.update(
    {
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.titlepad": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#444444",
        "axes.labelcolor": "#222222",
        "xtick.color": "#444444",
        "ytick.color": "#444444",
    }
)


def load():
    with open(IN) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def dist_or_inf(r):
    """Gene distance to nearest own-machinery; +inf when off-replicon / not adjacent."""
    return int(r["nearest_dist"]) if r["nearest_dist"] != "" else np.inf


def fig01_composition(rows, out):
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    reach, imp, unt = [], [], []
    for ss in ORDER:
        rs = [r for r in rows if r["ss_type"] == ss]
        t = [r for r in rs if r["testable"] == "yes"]
        r3 = sum(r["reachable_n3"] == "true" for r in t)
        reach.append(r3)
        imp.append(len(t) - r3)
        unt.append(len(rs) - len(t))
    x = np.arange(len(ORDER))
    ax.bar(x, reach, color=THEME["reachable"], label="reachable @N=3")
    ax.bar(x, imp, bottom=reach, color=THEME["impossible"], label="impossible (testable, >3)")
    ax.bar(x, unt, bottom=np.array(reach) + np.array(imp), color=THEME["untestable"], label="untestable")
    for i, ss in enumerate(ORDER):
        tot = reach[i] + imp[i] + unt[i]
        ax.text(i, tot + 2, str(tot), ha="center", va="bottom", fontsize=9, color="#222")
    ax.set_xticks(x)
    ax.set_xticklabels(ORDER)
    ax.set_ylabel("verified effectors (count)")
    ax.set_title("Gold-set composition by reachability of the +/-3 proximity rule")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.set_ylim(top=max(reach[i] + imp[i] + unt[i] for i in range(len(ORDER))) * 1.18)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig02_ceiling_vs_window(rows, out):
    fig, ax = plt.subplots(figsize=(7, 4.4))
    for ss in ORDER:
        t = [r for r in rows if r["ss_type"] == ss and r["testable"] == "yes"]
        if not t:
            continue
        y = [100 * sum(r[f"reachable_n{n}"] == "true" for r in t) / len(t) for n in NS]
        ax.plot(NS, y, "-o", color=THEME["ss_colors"][ss], lw=2, ms=6, label=f"{ss} (n={len(t)})")
        ax.annotate(
            f"{y[-1]:.0f}%",
            (NS[-1], y[-1]),
            textcoords="offset points",
            xytext=(6, 0),
            va="center",
            fontsize=9,
            color=THEME["ss_colors"][ss],
        )
    ax.set_xticks(NS)
    ax.set_xlabel("proximity window N (genes either side of a machinery component)")
    ax.set_ylabel("ceiling: % of testable effectors reachable")
    ax.set_title("Ceiling vs proximity window, per secretion-system type")
    ax.set_ylim(-3, 100)
    ax.set_xlim(2.6, 7.8)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig03_distance_ecdf(rows, out):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for ss in ORDER:
        t = [r for r in rows if r["ss_type"] == ss and r["testable"] == "yes"]
        if not t:
            continue
        d = np.sort([dist_or_inf(r) for r in t])
        finite = d[np.isfinite(d)]
        if len(finite) == 0:
            continue
        # ECDF clipped at the largest finite distance; off-replicon points keep the curve <1
        xs = np.concatenate([[0.5], finite])
        ys = np.searchsorted(finite, xs, side="right") / len(d)
        ax.step(xs, 100 * ys, where="post", color=THEME["ss_colors"][ss], lw=2, label=f"{ss} (n={len(t)})")
    for k, n in enumerate(NS):
        ax.axvline(n, color=THEME["ref_line"], ls="--", lw=0.7, alpha=0.7)
        ax.text(n, 101 + 4 * (k % 2), f"N={n}", ha="center", va="bottom", fontsize=8, color=THEME["ref_line"])
    ax.set_xscale("log")
    ax.set_xlabel("gene-order distance to nearest own-instance machinery component (log)")
    ax.set_ylabel("% of testable effectors within distance")
    ax.set_title("How far effectors sit from their own machinery")
    ax.set_ylim(0, 108)
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig04_per_genome(rows, out):
    by = defaultdict(list)
    for r in rows:
        if r["testable"] == "yes":
            by[(r["ss_type"], r["refseq_genome"])].append(r)
    items = [(ss, g, rs) for (ss, g), rs in by.items() if len(rs) >= 3]
    items.sort(key=lambda t: (100 * sum(r["reachable_n5"] == "true" for r in t[2]) / len(t[2]), len(t[2])))
    labels, vals, cols = [], [], []
    for ss, g, rs in items:
        labels.append(f"{ss}  {g}  (n={len(rs)})")
        vals.append(100 * sum(r["reachable_n5"] == "true" for r in rs) / len(rs))
        cols.append(THEME["ss_colors"][ss])
    fig, ax = plt.subplots(figsize=(7.4, max(3.2, 0.34 * len(labels) + 1)))
    y = np.arange(len(labels))
    ax.barh(y, vals, color=cols)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("ceiling @N=5 (% of that genome's testable effectors reachable)")
    ax.set_title("Per-genome ceiling (genomes with >=3 testable effectors)")
    ax.set_xlim(0, 105)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    for stale in FIGDIR.glob("*.png"):
        if not stale.name[:2].isdigit():
            stale.unlink()
    rows = load()
    figs = [
        ("01_gold_set_composition.png", fig01_composition),
        ("02_ceiling_vs_window.png", fig02_ceiling_vs_window),
        ("03_distance_ecdf.png", fig03_distance_ecdf),
        ("04_per_genome_ceiling.png", fig04_per_genome),
    ]
    for name, fn in figs:
        fn(rows, FIGDIR / name)
    print("Figure index:")
    for name, _ in figs:
        print(f"  {name}")


if __name__ == "__main__":
    raise SystemExit(main())
