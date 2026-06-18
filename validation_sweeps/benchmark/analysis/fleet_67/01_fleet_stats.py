#!/usr/bin/env python3
"""Fleet-wide stats over the 67-genome enrichment panel (job 3016930-35 + BX470251 rerun).

Three things Teo asked for (this script covers 1-2; calibration is 02_*):
  1. How the predictor scores + enrichment look overall, across all 67 genomes:
     - per-protein DLP extracellular / DSE score distributions (fleet-wide)
     - per-genome true background positive rate (the thing the enrichment test estimates)
     - enrichment results aggregated by SS type (significant calls, fold-enrichment)
  2. How well an n=1000 null sample estimates the true whole-genome background:
     the live runs used the EXACT background (whole-genome predictors), so we
     validate the n=1000 default by resampling from each genome's full predictions
     and measuring the estimate's error vs the exact value, across all 67 genomes.

    .venv/bin/python 01_fleet_stats.py
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

from ssign_app.scripts.enrichment_testing import is_dlp_positive, is_dse_positive

HERE = os.path.dirname(os.path.abspath(__file__))
FLEET = "/tmp/ssign_fleet_67"
FIGS = os.path.join(HERE, "figures")
CONF = 0.8
SEED = 42
RESAMPLE_REPS = 300
NULL_SIZES = [200, 1000]

THEME = {
    "dlp": "#3F8E8C",
    "dse": "#E0884B",
    "n200": "#C44E52",
    "n1000": "#4C72B0",
    "neutral": "#6C8EAD",
    "ref": "#444444",
    "hi": "#A93232",
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


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def load_genomes() -> list[dict]:
    """Per genome: id, n, dlp_probs[], dse_probs[], dlp_pos[], dse_pos[]."""
    out = []
    for g in sorted(os.listdir(FLEET)):
        raw = os.path.join(FLEET, g, "results", f"{g}_results_raw.csv")
        if not os.path.exists(raw):
            continue
        dlp_probs, dse_probs, dlp_pos, dse_pos = [], [], [], []
        with open(raw) as fh:
            for row in csv.DictReader(fh):
                dlp_probs.append(_f(row.get("dlp_extracellular_prob")))
                dse_probs.append(_f(row.get("dse_max_prob")))
                dlp_pos.append(is_dlp_positive(row, CONF))
                dse_pos.append(is_dse_positive(row, CONF))
        if not dlp_pos:
            continue
        out.append(
            {
                "id": g,
                "n": len(dlp_pos),
                "dlp_probs": np.array(dlp_probs),
                "dse_probs": np.array(dse_probs),
                "dlp_pos": np.array(dlp_pos, dtype=bool),
                "dse_pos": np.array(dse_pos, dtype=bool),
            }
        )
    return out


def parse_enrichment(genomes_ids: list[str]) -> list[dict]:
    rows = []
    for g in genomes_ids:
        summ = os.path.join(FLEET, g, "results", f"{g}_summary.txt")
        if not os.path.exists(summ):
            continue
        with open(summ) as fh:
            for line in fh:
                if not re.match(r"\s*" + re.escape(g) + r"\s+(system|broad_type)\s", line):
                    continue
                p = line.split()
                # sample kind scope_id ss_type tool M k p_bg fold pvalue qvalue significant n_null
                try:
                    rows.append(
                        {
                            "genome": g,
                            "kind": p[1],
                            "ss_type": p[3],
                            "tool": p[4],
                            "M": int(p[5]),
                            "k": int(p[6]),
                            "fold": float(p[8]) if p[8] not in ("", "None") else 0.0,
                            "qvalue": float(p[10]),
                            "significant": p[11] == "True",
                        }
                    )
                except (ValueError, IndexError):
                    continue
    return rows


# --- analysis 2: n=1000 resampling vs exact background ---


def resample_validation(genomes: list[dict]) -> dict:
    """Per genome with enough proteins: relative error of the sampled background
    estimate vs the exact whole-genome rate, for each null size."""
    rng = random.Random(SEED)
    res = {n: {"dlp": [], "dse": []} for n in NULL_SIZES}
    per_genome = []
    for gd in genomes:
        if gd["n"] < max(NULL_SIZES):
            continue
        idx = list(range(gd["n"]))
        row = {"id": gd["id"], "n": gd["n"]}
        for tool, posarr in (("dlp", gd["dlp_pos"]), ("dse", gd["dse_pos"])):
            p_true = posarr.mean()
            row[f"p_true_{tool}"] = p_true
            if p_true <= 0:
                continue
            for n in NULL_SIZES:
                errs = []
                for _ in range(RESAMPLE_REPS):
                    s = rng.sample(idx, n)
                    p_hat = posarr[s].mean()
                    errs.append(abs(p_hat - p_true) / p_true)
                med = float(np.median(errs))
                res[n][tool].append(med)
                row[f"relerr_{tool}_{n}"] = med
        per_genome.append(row)
    return {"agg": res, "per_genome": per_genome}


# --- figures ---


def fig_score_distributions(genomes, n_proteins):
    dlp = np.concatenate([g["dlp_probs"] for g in genomes])
    dse = np.concatenate([g["dse_probs"] for g in genomes])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.4))
    a1.hist(dlp, bins=50, color=THEME["dlp"], alpha=0.85)
    a1.axvline(CONF, ls="--", color=THEME["ref"], lw=0.8, label=f"threshold {CONF}")
    a1.set_yscale("log")
    a1.set_xlabel("DLP extracellular probability")
    a1.set_ylabel("proteins (log)")
    a1.set_title(f"DLP scores, all {n_proteins:,} proteins / 67 genomes\n{(dlp >= CONF).mean():.2%} ≥ {CONF}")
    a1.legend(frameon=False, fontsize=8)
    a2.hist(dse, bins=50, color=THEME["dse"], alpha=0.85)
    a2.axvline(CONF, ls="--", color=THEME["ref"], lw=0.8, label=f"threshold {CONF}")
    a2.set_yscale("log")
    a2.set_xlabel("DSE max secretion-type probability")
    a2.set_ylabel("proteins (log)")
    a2.set_title("DSE scores, all proteins\n(positivity also requires a real SS type)")
    a2.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "01_score_distributions.png"))
    plt.close(fig)


def fig_per_genome_background(genomes):
    gs = sorted(genomes, key=lambda g: g["n"])
    ns = [g["n"] for g in gs]
    p_dlp = [100 * g["dlp_pos"].mean() for g in gs]
    p_dse = [100 * g["dse_pos"].mean() for g in gs]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(ns, p_dlp, s=28, color=THEME["dlp"], label="DLP", alpha=0.85)
    ax.scatter(ns, p_dse, s=28, color=THEME["dse"], label="DSE", alpha=0.85)
    ax.set_xlabel("genome size (proteins)")
    ax.set_ylabel("true background positive rate (%)")
    ax.set_title("Per-genome background rate (the value the enrichment test must estimate)")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "02_per_genome_background.png"))
    plt.close(fig)


def fig_enrichment_by_sstype(enr):
    def broad(t):
        m = re.match(r"p?(T\d+)[a-z]*SS", t)
        return f"{m.group(1)}SS" if m else t

    by = defaultdict(lambda: {"DLP": [0, 0], "DSE": [0, 0]})  # bt -> tool -> [sig, total]
    for r in enr:
        if r["kind"] != "system":
            continue
        bt = broad(r["ss_type"])
        if r["tool"] in ("DLP", "DSE"):
            by[bt][r["tool"]][1] += 1
            if r["significant"]:
                by[bt][r["tool"]][0] += 1
    bts = sorted(by, key=lambda b: -(by[b]["DLP"][1] + by[b]["DSE"][1]))
    x = np.arange(len(bts))
    w = 0.38
    fig, ax = plt.subplots(figsize=(11, 4.8))
    for i, tool in enumerate(("DLP", "DSE")):
        sig = [by[b][tool][0] for b in bts]
        tot = [by[b][tool][1] for b in bts]
        ax.bar(
            x + (i - 0.5) * w,
            sig,
            w,
            color=THEME["dlp"] if tool == "DLP" else THEME["dse"],
            label=f"{tool} significant",
        )
        for xi, (s, t) in enumerate(zip(sig, tot)):
            ax.text(xi + (i - 0.5) * w, s, f"{s}/{t}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(bts, rotation=30, ha="right")
    ax.set_ylabel("significant systems (q<0.05)")
    ax.set_title("Enrichment significance by secretion-system type (per-system tests, 67 genomes)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "03_enrichment_by_sstype.png"))
    plt.close(fig)


def fig_nnull_validation(val):
    agg = val["agg"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, tool, title in ((a1, "dlp", "DLP"), (a2, "dse", "DSE")):
        data = [np.array(agg[n][tool]) * 100 for n in NULL_SIZES]
        bp = ax.boxplot(data, tick_labels=[f"n={n}" for n in NULL_SIZES], showfliers=False, patch_artist=True)
        for patch, c in zip(bp["boxes"], (THEME["n200"], THEME["n1000"])):
            patch.set_facecolor(c)
            patch.set_alpha(0.75)
        meds = [np.median(d) for d in data]
        for i, m in enumerate(meds, start=1):
            ax.text(i, m, f"{m:.0f}%", ha="center", va="bottom", fontsize=9, weight="bold")
        ax.set_ylabel("median |estimate − true| / true  (%)")
        ax.set_title(f"{title}: background-estimate error vs exact, across 67 genomes")
    fig.suptitle(
        "How well a size-n null sample estimates the true background (resampled)", y=1.02, fontsize=11, weight="bold"
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "04_nnull_validation.png"))
    plt.close(fig)


def main():
    os.makedirs(FIGS, exist_ok=True)
    genomes = load_genomes()
    total_prot = sum(g["n"] for g in genomes)
    print(
        f"loaded {len(genomes)} genomes, {total_prot:,} proteins total "
        f"(sizes {min(g['n'] for g in genomes)}–{max(g['n'] for g in genomes)})"
    )

    # --- 1. overall stats ---
    dlp_rate = np.mean([g["dlp_pos"].mean() for g in genomes])
    dse_rate = np.mean([g["dse_pos"].mean() for g in genomes])
    print(f"\nmean per-genome background: DLP {dlp_rate:.3%}, DSE {dse_rate:.3%}")
    print(
        f"per-genome DLP range {min(g['dlp_pos'].mean() for g in genomes):.2%}–{max(g['dlp_pos'].mean() for g in genomes):.2%}"
    )

    enr = parse_enrichment([g["id"] for g in genomes])
    sysrows = [r for r in enr if r["kind"] == "system"]
    n_sig = sum(1 for r in sysrows if r["significant"])
    print(
        f"\nenrichment: {len(sysrows)} (system × tool) tests across {len({r['genome'] for r in enr})} genomes; "
        f"{n_sig} significant (q<0.05)"
    )

    # --- 2. n=1000 validation ---
    val = resample_validation(genomes)
    print(f"\nn=1000 background-estimate validation ({len(val['per_genome'])} genomes ≥ {max(NULL_SIZES)} proteins):")
    for tool in ("dlp", "dse"):
        for n in NULL_SIZES:
            arr = np.array(val["agg"][n][tool]) * 100
            print(
                f"  {tool.upper()} n={n:<4d}: median rel-err {np.median(arr):5.1f}%  (IQR {np.percentile(arr, 25):.0f}–{np.percentile(arr, 75):.0f}%)"
            )

    fig_score_distributions(genomes, total_prot)
    fig_per_genome_background(genomes)
    fig_enrichment_by_sstype(enr)
    fig_nnull_validation(val)
    print("\nFigure index:")
    print("  01_score_distributions.png   — fleet-wide DLP/DSE per-protein score histograms")
    print("  02_per_genome_background.png — true background rate vs genome size")
    print("  03_enrichment_by_sstype.png  — significant enrichment calls per SS type")
    print("  04_nnull_validation.png      — n=200 vs n=1000 background-estimate error vs exact")


if __name__ == "__main__":
    main()
