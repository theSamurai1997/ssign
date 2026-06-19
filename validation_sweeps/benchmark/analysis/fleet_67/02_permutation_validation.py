#!/usr/bin/env python3
"""Validate the binomial enrichment test against a circular-permutation null.

The production test is a one-sided binomial: it assumes each neighborhood
protein is independently secreted at the genome background rate. But secreted
genes cluster (operons, islands), so independence can inflate significance.
A circular-permutation null preserves that clustering: for a system whose
neighborhood holds k secreted proteins out of M, slide an M-gene window around
the genome's (circular) gene order and count secreted at every offset — the
observed k's rank in that distribution is a spatially-honest p-value.

This computes, per system (DLP and DSE), the permutation p-value and a binomial
p-value with an n=1000 *sampled* background (the user-facing default), BH-corrects
within genome (matching production), and compares against the binomial-with-exact
background calls the fleet actually made. Answers: does the binomial agree with
an assumption-light spatial test, and does n=1000 catch what permutation catches?

    .venv/bin/python 02_permutation_validation.py

Approximations (noted): the neighborhood is treated as one contiguous M-window
(true for contiguous systems); the sliding window spans all proteins incl. SS
components (the real neighborhood excludes them, but components are rarely
secreted); multi-contig genomes are ordered contig-then-position and wrapped.
"""

from __future__ import annotations

import csv
import os
import random
import re
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ssign_app.scripts.enrichment_testing import (
    bh_fdr,
    binom_pvalue,
    is_dlp_positive,
    is_dse_positive,
)

FLEET = "/tmp/ssign_fleet_67"
HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figures")
CONF = 0.8
SEED = 42
ALPHA = 0.05

THEME = {"both": "#3F8E8C", "binom_only": "#C44E52", "perm_only": "#7A5C9E", "neither": "#C9C9C9", "ref": "#444444"}
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
        "xtick.color": "#444444",
        "ytick.color": "#444444",
    }
)


def ordered_secreted_vectors(genome: str):
    """(dlp_vec, dse_vec) of secreted-positivity bools in circular gene order."""
    raw = os.path.join(FLEET, genome, "results", f"{genome}_results_raw.csv")
    recs = []
    with open(raw) as fh:
        for row in csv.DictReader(fh):
            try:
                start = int(float(row.get("start") or 0))
            except ValueError:
                start = 0
            recs.append((row.get("contig", ""), start, is_dlp_positive(row, CONF), is_dse_positive(row, CONF)))
    recs.sort(key=lambda r: (r[0], r[1]))
    dlp = np.array([r[2] for r in recs], dtype=int)
    dse = np.array([r[3] for r in recs], dtype=int)
    return dlp, dse


def perm_pvalue(vec: np.ndarray, M: int, k: int) -> float:
    """P(secreted-count in a random circular M-window >= k), +1 smoothed."""
    n = len(vec)
    if M <= 0 or M >= n or k <= 0:
        return 1.0
    vv = np.concatenate([vec, vec[:M]])
    cs = np.concatenate([[0], np.cumsum(vv)])
    counts = cs[M : M + n] - cs[:n]  # n circular windows of size M
    return float((np.sum(counts >= k) + 1) / (n + 1))


def parse_systems(genome: str):
    """System rows from the genome summary: (ss_type, tool, M, k, significant_binom_exact)."""
    summ = os.path.join(FLEET, genome, "results", f"{genome}_summary.txt")
    out = []
    if not os.path.exists(summ):
        return out
    with open(summ) as fh:
        for line in fh:
            if not re.match(r"\s*" + re.escape(genome) + r"\s+system\s", line):
                continue
            p = line.split()
            try:
                out.append(
                    {"ss_type": p[3], "tool": p[4], "M": int(p[5]), "k": int(p[6]), "sig_binom_exact": p[11] == "True"}
                )
            except (ValueError, IndexError):
                continue
    return out


def main():
    os.makedirs(FIGS, exist_ok=True)
    rng = random.Random(SEED)
    rows = []
    for g in sorted(os.listdir(FLEET)):
        systems = parse_systems(g)
        if not systems:
            continue
        dlp, dse = ordered_secreted_vectors(g)
        vecs = {"DLP": dlp, "DSE": dse}
        grows = []
        for s in systems:
            vec = vecs.get(s["tool"])
            if vec is None or len(vec) == 0:
                continue
            n = len(vec)
            p_perm = perm_pvalue(vec, s["M"], s["k"])
            # binomial with an n=1000 sampled background (default users get this)
            n_draw = min(1000, n)
            p_bg = vec[rng.sample(range(n), n_draw)].mean()
            p_n1000 = binom_pvalue(s["k"], s["M"], float(p_bg))
            grows.append({**s, "genome": g, "pvalue": p_perm, "p_n1000": p_n1000})
        # BH within genome (matches production), separately for each p source
        bh_fdr(grows, pvalue_key="pvalue", q_key="q_perm", sig_key="sig_perm")
        bh_fdr(grows, pvalue_key="p_n1000", q_key="q_n1000", sig_key="sig_n1000")
        rows.extend(grows)

    # --- agreement summary ---
    def counts(key_a, key_b):
        c = defaultdict(int)
        for r in rows:
            c[(r[key_a], r[key_b])] += 1
        return c

    n = len(rows)
    perm_sig = sum(r["sig_perm"] for r in rows)
    bex_sig = sum(r["sig_binom_exact"] for r in rows)
    bn1000_sig = sum(r["sig_n1000"] for r in rows)
    print(f"{n} per-system tests (DLP+DSE) across {len({r['genome'] for r in rows})} genomes\n")
    print(f"significant: permutation {perm_sig}  |  binomial-exact {bex_sig}  |  binomial-n1000 {bn1000_sig}\n")

    for ref, lab in (("sig_perm", "permutation"), ("sig_binom_exact", "binomial-exact")):
        c = counts("sig_n1000", ref)
        both = c[(True, True)]
        only_ref = c[(False, True)]
        only_n1000 = c[(True, False)]
        recall = both / (both + only_ref) if (both + only_ref) else float("nan")
        print(
            f"binomial-n1000 vs {lab}: both {both}, {lab}-only {only_ref}, n1000-only {only_n1000} "
            f"-> n1000 recovers {recall:.0%} of {lab}-significant systems"
        )
    # permutation vs binomial-exact: does the spatial null confirm the production calls?
    c = counts("sig_perm", "sig_binom_exact")
    print(
        f"\npermutation vs binomial-exact: both {c[(True, True)]}, binom-only {c[(False, True)]}, "
        f"perm-only {c[(True, False)]}, neither {c[(False, False)]}"
    )
    print(
        f"  -> of {bex_sig} production-significant systems, permutation confirms "
        f"{c[(True, True)]} ({c[(True, True)] / bex_sig:.0%})"
    )

    make_figure(rows)
    print("\nFigure: 05_permutation_vs_binomial.png")


def make_figure(rows):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2))

    # left: perm p vs binomial-n1000 p, colored by agreement
    def cat(r):
        if r["sig_perm"] and r["sig_n1000"]:
            return "both"
        if r["sig_n1000"]:
            return "binom_only"
        if r["sig_perm"]:
            return "perm_only"
        return "neither"

    for c in ("neither", "both", "binom_only", "perm_only"):
        pts = [r for r in rows if cat(r) == c]
        if not pts:
            continue
        a1.scatter(
            [max(r["p_n1000"], 1e-6) for r in pts],
            [max(r["pvalue"], 1e-6) for r in pts],
            s=18,
            color=THEME[c],
            alpha=0.7,
            label=f"{c.replace('_', '-')} ({len(pts)})",
        )
    a1.plot([1e-6, 1], [1e-6, 1], ls=":", color=THEME["ref"], lw=0.8)
    a1.axhline(ALPHA, ls="--", color=THEME["ref"], lw=0.6)
    a1.axvline(ALPHA, ls="--", color=THEME["ref"], lw=0.6)
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xlabel("binomial (n=1000) p-value")
    a1.set_ylabel("circular-permutation p-value")
    a1.set_title("Per-system: binomial vs spatial-permutation p-value")
    a1.legend(frameon=False, fontsize=8, loc="lower right")

    # right: significant-call counts per method
    perm = sum(r["sig_perm"] for r in rows)
    bex = sum(r["sig_binom_exact"] for r in rows)
    bn = sum(r["sig_n1000"] for r in rows)
    labels = ["permutation\n(spatial null)", "binomial\nexact bg", "binomial\nn=1000 bg"]
    vals = [perm, bex, bn]
    a2.bar(labels, vals, color=["#7A5C9E", "#3F8E8C", "#4C72B0"])
    for i, v in enumerate(vals):
        a2.text(i, v, str(v), ha="center", va="bottom", weight="bold")
    a2.set_ylabel("significant systems (q < 0.05)")
    a2.set_title("Significant systems by test method")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "05_permutation_vs_binomial.png"))
    plt.close(fig)


if __name__ == "__main__":
    main()
