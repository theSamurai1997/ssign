#!/usr/bin/env python3
"""Phase 1 enrichment-test validation on the PAO1 smoke run (job 3013556).

Two questions, both answered from the pulled whole-genome run output (no
re-prediction needed):

  A. Is the 200-protein null background representative? Recompute every
     DLP/DSE enrichment test with the EXACT whole-genome background rate
     (holding each scope's neighborhood hit-count k and size M fixed) and
     count how many significance calls flip.

  B. Why does PLM-Effector call ~25% of the proteome secreted? Characterize
     its max-probability distribution, per-type breakdown, and how the
     positive rate would fall under a stricter probability gate.

Reads the pipeline's own positivity rules + stats from
ssign_app.scripts.enrichment_testing so the recompute matches production.

    .venv/bin/python 01_analyze_background_and_plme.py
"""

from __future__ import annotations

import csv
import os
import re
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ssign_app.scripts.enrichment_testing import (
    bh_fdr,
    binom_pvalue,
    is_dlp_positive,
    is_dse_positive,
    is_plme_positive,
)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIGS = os.path.join(HERE, "figures")
RAW = os.path.join(DATA, "NC_002516.2_results_raw.csv")
SUMMARY = os.path.join(DATA, "NC_002516.2_summary.txt")
CONF = 0.8  # conf_threshold the run used (config.conf_threshold)

THEME = {
    "tool_colors": {"DLP": "#3F8E8C", "DSE": "#E0884B", "PLME": "#7A5C9E"},
    "sampled": "#C44E52",
    "exact": "#4C72B0",
    "neutral_bar": "#6C8EAD",
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


def load_raw() -> list[dict]:
    with open(RAW) as f:
        return list(csv.DictReader(f))


def exact_backgrounds(rows: list[dict]) -> dict[str, float]:
    """Whole-genome positive fraction per tool, using the production rules."""
    n = len(rows)
    return {
        "DLP": sum(1 for r in rows if is_dlp_positive(r, CONF)) / n,
        "DSE": sum(1 for r in rows if is_dse_positive(r, CONF)) / n,
        "PLME": sum(1 for r in rows if is_plme_positive(r)) / n,
    }


def parse_enrichment_table(path: str) -> list[dict]:
    """Pull the per-scope rows out of the run summary's ENRICHMENT block."""
    out = []
    with open(path) as f:
        for line in f:
            if not re.match(r"\s*NC_002516\.2\s+system", line):
                continue
            p = line.split()
            # sample scope_kind scope_id ss_type tool M k p_bg fold pvalue qvalue significant n_null
            out.append(
                {
                    "scope_id": p[2],
                    "ss_type": p[3],
                    "tool": p[4],
                    "M": int(p[5]),
                    "k": int(p[6]),
                    "p_bg_sampled": float(p[7]),
                    "pvalue_sampled": float(p[9]),
                    "significant_sampled": p[11] == "True",
                }
            )
    return out


def recompute_with_exact_bg(table: list[dict], bg: dict[str, float]) -> list[dict]:
    """Re-run the binomial + BH FDR with the exact whole-genome background.

    k and M are physical facts of the run (same predictions, same
    neighborhoods); only p_bg changes, so this isolates the background's
    effect on significance.
    """
    rows = []
    for t in table:
        p_true = bg[t["tool"]]
        rows.append(
            {
                **t,
                "p_bg_exact": p_true,
                "pvalue": binom_pvalue(t["k"], t["M"], p_true),
                "scope_id": t["scope_id"],
            }
        )
    bh_fdr(rows)  # adds qvalue + significant in place, BH across all (scope x tool)
    for r in rows:
        r["significant_exact"] = r.pop("significant")
        r["qvalue_exact"] = r.pop("qvalue")
    return rows


def plme_threshold_sweep(rows: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    """Positive rate if PLM-E required max_prob >= t, swept over t.

    Illustrative only: the production call is an OR over five per-type
    ensemble thresholds, not a single max-prob cut, but the max-prob
    distribution shows how many calls are marginal vs confident.
    """
    probs = np.array([float(r.get("plm_effector_max_prob") or 0) for r in rows if is_plme_positive(r)])
    ts = np.linspace(0.5, 1.0, 51)
    rate = np.array([(probs >= t).sum() / len(rows) for t in ts])
    return ts, rate, probs


def fig_background(bg: dict, recomp: list[dict]) -> str:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    tools = ["DLP", "DSE"]
    sampled = {t: next(r["p_bg_sampled"] for r in recomp if r["tool"] == t) for t in tools}
    x = np.arange(len(tools))
    w = 0.36
    ax1.bar(x - w / 2, [sampled[t] * 100 for t in tools], w, label="200-null sample", color=THEME["sampled"])
    ax1.bar(x + w / 2, [bg[t] * 100 for t in tools], w, label="whole-genome exact", color=THEME["exact"])
    for i, t in enumerate(tools):
        ax1.text(i - w / 2, sampled[t] * 100, f"{sampled[t] * 100:.2f}%", ha="center", va="bottom", fontsize=8)
        ax1.text(i + w / 2, bg[t] * 100, f"{bg[t] * 100:.2f}%", ha="center", va="bottom", fontsize=8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(tools)
    ax1.set_ylabel("background positive rate (%)")
    ax1.set_title("Background rate: 200-null sample vs whole genome")
    ax1.legend(frameon=False, fontsize=8)

    # significance count before/after
    sig_old = sum(1 for r in recomp if r["significant_sampled"])
    sig_new = sum(1 for r in recomp if r["significant_exact"])
    flips_lost = sum(1 for r in recomp if r["significant_sampled"] and not r["significant_exact"])
    flips_gain = sum(1 for r in recomp if not r["significant_sampled"] and r["significant_exact"])
    ax2.bar(["200-null", "whole-genome"], [sig_old, sig_new], color=[THEME["sampled"], THEME["exact"]])
    ax2.set_ylabel("significant (scope x tool) tests, q < 0.05")
    ax2.set_title(f"Significant calls: {sig_old} -> {sig_new}\n({flips_lost} lost, {flips_gain} gained with exact bg)")
    for i, v in enumerate([sig_old, sig_new]):
        ax2.text(i, v, str(v), ha="center", va="bottom")
    fig.tight_layout()
    out = os.path.join(FIGS, "01_background_estimate.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def fig_flips(recomp: list[dict]) -> str:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for r in recomp:
        po = max(r["pvalue_sampled"], 1e-7)
        pn = max(r["pvalue"], 1e-7)
        flipped = r["significant_sampled"] != r["significant_exact"]
        ax.scatter(
            po,
            pn,
            s=34,
            color=THEME["tool_colors"][r["tool"]],
            edgecolor="#A93232" if flipped else "none",
            linewidth=1.6 if flipped else 0,
            alpha=0.85,
        )
    lims = [1e-7, 1.5]
    ax.plot(lims, lims, ls=":", color=THEME["ref_line"], lw=0.8)
    ax.axhline(0.05, ls="--", color=THEME["ref_line"], lw=0.6, alpha=0.7)
    ax.axvline(0.05, ls="--", color=THEME["ref_line"], lw=0.6, alpha=0.7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lims)
    ax.set_ylim(*lims)
    ax.set_xlabel("pvalue (200-null background)")
    ax.set_ylabel("pvalue (whole-genome background)")
    ax.set_title("Per-scope p-value shift (red ring = significance flipped)")
    handles = [
        plt.Line2D([], [], marker="o", ls="", color=c, label=t) for t, c in THEME["tool_colors"].items() if t != "PLME"
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, title="tool")
    fig.tight_layout()
    out = os.path.join(FIGS, "02_significance_flips.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def fig_plme(rows: list[dict]) -> str:
    ts, rate, probs = plme_threshold_sweep(rows)
    n = len(rows)
    # plm_effector_type is a comma-joined list of the per-type ensembles that fired
    # for this protein; OR-count each type, and track multi-type (cross-reactive) calls.
    types = Counter()
    n_multi = 0
    for r in rows:
        if is_plme_positive(r):
            tys = [t.strip() for t in str(r.get("plm_effector_type") or "").split(",") if t.strip()]
            if len(tys) > 1:
                n_multi += 1
            for ty in tys:
                types[ty] += 1
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    # (a) max_prob distribution of positive calls
    axes[0].hist(probs, bins=30, color=THEME["tool_colors"]["PLME"], alpha=0.85)
    axes[0].set_xlabel("plm_effector_max_prob (positive calls)")
    axes[0].set_ylabel("proteins")
    axes[0].set_title(f"PLM-E positive calls: {len(probs)} / {n} = {len(probs) / n:.1%}")
    # (b) per-type
    if types:
        labs, vals = zip(*sorted(types.items(), key=lambda kv: -kv[1]))
        axes[1].bar(labs, vals, color=THEME["neutral_bar"])
        axes[1].set_ylabel("positive calls (type, OR-counted)")
        axes[1].set_title(f"PLM-E positives by type ({n_multi} multi-type)")
        for i, v in enumerate(vals):
            axes[1].text(i, v, str(v), ha="center", va="bottom", fontsize=8)
    # (c) threshold sweep
    axes[2].plot(ts, rate * 100, "-o", ms=3, color=THEME["tool_colors"]["PLME"])
    axes[2].axhline(1.65, ls="--", color=THEME["ref_line"], lw=0.7, alpha=0.8, label="DLP/DSE bg ~1.7%")
    axes[2].set_xlabel("max_prob threshold")
    axes[2].set_ylabel("genome positive rate (%)")
    axes[2].set_title("PLM-E positive rate vs stricter prob gate")
    axes[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    out = os.path.join(FIGS, "03_plme_genome_behaviour.png")
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> None:
    for stale in os.listdir(FIGS) if os.path.isdir(FIGS) else []:
        if stale.endswith(".png") and not stale[:2].isdigit():
            os.unlink(os.path.join(FIGS, stale))
    rows = load_raw()
    bg = exact_backgrounds(rows)
    table = parse_enrichment_table(SUMMARY)
    recomp = recompute_with_exact_bg(table, bg)

    n = len(rows)
    print(f"=== PAO1 (NC_002516.2), {n} proteins, conf={CONF} ===\n")
    print("Exact whole-genome background vs the 200-null estimate the run used:")
    for t in ("DLP", "DSE"):
        s = next(r["p_bg_sampled"] for r in recomp if r["tool"] == t)
        print(
            f"  {t:4s} exact={bg[t]:.4f} ({bg[t] * n:.0f}/{n})   sampled={s:.4f}   ratio exact/sampled={bg[t] / s:.1f}x"
        )
    print(f"  PLME exact={bg['PLME']:.4f} ({bg['PLME'] * n:.0f}/{n})   (not tested in this run)\n")

    sig_old = sum(1 for r in recomp if r["significant_sampled"])
    sig_new = sum(1 for r in recomp if r["significant_exact"])
    print(f"Significant (scope x tool) calls @ q<0.05:  200-null={sig_old}   whole-genome={sig_new}")
    print("Flips (significant under sampled bg, NOT under exact bg):")
    for r in recomp:
        if r["significant_sampled"] and not r["significant_exact"]:
            print(
                f"  {r['scope_id']:34s} {r['tool']:4s} k={r['k']}/M={r['M']}  "
                f"p {r['pvalue_sampled']:.4g}->{r['pvalue']:.4g}  q_exact={r['qvalue_exact']:.4g}"
            )
    gained = [r for r in recomp if not r["significant_sampled"] and r["significant_exact"]]
    if gained:
        print("Flips (newly significant under exact bg):")
        for r in gained:
            print(f"  {r['scope_id']:34s} {r['tool']:4s} k={r['k']}/M={r['M']}  q_exact={r['qvalue_exact']:.4g}")

    f1 = fig_background(bg, recomp)
    f2 = fig_flips(recomp)
    f3 = fig_plme(rows)
    print("\nFigure index:")
    print(f"  01  {os.path.basename(f1)}  — sampled vs exact background + significant-call counts")
    print(f"  02  {os.path.basename(f2)}  — per-scope p-value shift, flagged flips")
    print(f"  03  {os.path.basename(f3)}  — PLM-E max_prob dist, per-type, threshold sweep")


if __name__ == "__main__":
    main()
