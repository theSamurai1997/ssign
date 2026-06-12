#!/usr/bin/env python3
"""
27_phase2_figures.py  (Phase 2 task 6.6: figures for the full-panel actual-vs-ceiling result)

Reads the scored Phase-2 tables and renders the benchmark's headline figures:
  01  recall funnel        : testable -> findable@3/5/7 -> found  (both run tags)
  02  per-type recall      : ceiling@7 vs actual emitted, per SS type, both tags
  03  emission basis       : found effectors split own-type (legit) vs cross-type-only (accidental)
  04  discordant audit     : verdict partition of the 19 emitted-but-unreachable@7 effectors

Inputs : data/phase2/actual_vs_ceiling.<tag>.tsv   (per-type test/ceiling/actual)
         data/phase2/found_systems.<tag>.tsv        (per found effector, emission_basis)
Output : data/phase2/figures/01..04_*.png
Run:     .venv/bin/python scripts/27_phase2_figures.py
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BENCH = Path(__file__).resolve().parents[1]
P2 = BENCH / "data" / "phase2"
FIGDIR = P2 / "figures"
TAGS = {"panel_genbank_default": "T3SS excluded (as shipped)", "panel_genbank_t3ss": "T3SS included"}
SS_ORDER = ["T1SS", "T2SS", "T3SS", "T4SS", "T6SS"]

THEME = {
    "tag_colors": {"panel_genbank_default": "#3F8E8C", "panel_genbank_t3ss": "#E0884B"},
    "ceiling": "#6C8EAD",
    "found": "#A93232",
    "own": "#3F8E8C",
    "cross": "#A93232",
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


def load(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def vc(tag):
    return {r["ss_type"]: r for r in load(P2 / f"actual_vs_ceiling.{tag}.tsv")}


def fig01_funnel():
    """Testable -> findable@3/5/7 -> found, both tags. ceiling cols are run-independent."""
    a = vc("panel_genbank_default")["ALL"]
    stages = ["Testable\neffectors", "Findable\n@±3", "Findable\n@±5", "Findable\n@±7", "Found\nby ssign"]
    base = [int(a["n_testable"]), int(a["ceiling_n3"]), int(a["ceiling_n5"]), int(a["ceiling_n7"])]
    found = {t: int(vc(t)["ALL"]["actual_emitted"]) for t in TAGS}
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    x = np.arange(len(stages))
    # shared stages (testable + ceiling) as neutral bars; final "found" split by tag
    ax.bar(x[:4], base, color=THEME["neutral_bar"], width=0.6)
    w = 0.28
    for i, (t, lbl) in enumerate(TAGS.items()):
        ax.bar(x[4] + (i - 0.5) * w, found[t], width=w, color=THEME["tag_colors"][t], label=lbl)
    for xi, v in zip(x[:4], base):
        ax.text(xi, v + 4, str(v), ha="center", va="bottom", fontweight="bold")
    for i, t in enumerate(TAGS):
        ax.text(x[4] + (i - 0.5) * w, found[t] + 4, str(found[t]), ha="center", va="bottom", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(stages)
    ax.set_ylabel("effectors (testable set, n=499)")
    ax.set_title("Recall funnel: proximity could reach 127/499; ssign emits 39–51")
    ax.set_ylim(top=max(base) * 1.18)
    ax.legend(frameon=False, loc="upper right")
    _save(fig, "01_recall_funnel.png")


def fig02_per_type():
    fig, ax = plt.subplots(figsize=(9, 4.8))
    x = np.arange(len(SS_ORDER))
    w = 0.26
    d = vc("panel_genbank_default")
    ceil7 = [int(d[s]["ceiling_n7"]) for s in SS_ORDER]
    test = [int(d[s]["n_testable"]) for s in SS_ORDER]
    ax.bar(x - w, ceil7, width=w, color=THEME["ceiling"], label="findable @±7 (ceiling)")
    for i, t in enumerate(TAGS):
        vt = vc(t)
        act = [int(vt[s]["actual_emitted"]) for s in SS_ORDER]
        ax.bar(x + (i) * w, act, width=w, color=THEME["tag_colors"][t], label=f"found — {TAGS[t]}")
    for xi, c, n in zip(x, ceil7, test):
        ax.text(xi - w, c + 1, str(c), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\n(n={t})" for s, t in zip(SS_ORDER, test)])
    ax.set_ylabel("effectors")
    ax.set_title("Per-type recall vs proximity ceiling (T4SS: ceiling 10, found 0)")
    ax.legend(frameon=False, loc="upper left")
    _save(fig, "02_per_type_recall.png")


def fig03_emission_basis():
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True)
    for ax, t in zip(axes, TAGS):
        rows = load(P2 / f"found_systems.{t}.tsv")
        own = Counter(r["ss_type"] for r in rows if r["emission_basis"] == "own_type")
        cross = Counter(r["ss_type"] for r in rows if r["emission_basis"] == "cross_type_only")
        types = [s for s in SS_ORDER if own[s] + cross[s] > 0]
        x = np.arange(len(types))
        o = [own[s] for s in types]
        c = [cross[s] for s in types]
        ax.bar(x, o, color=THEME["own"], label="own-type system nearby (legit)")
        ax.bar(x, c, bottom=o, color=THEME["cross"], label="only a different-type system (accidental)")
        for xi, oi, ci in zip(x, o, c):
            if oi:
                ax.text(xi, oi / 2, str(oi), ha="center", va="center", color="white", fontsize=8)
            if ci:
                ax.text(xi, oi + ci / 2, str(ci), ha="center", va="center", color="white", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(types)
        n_own, n_cross = sum(o), sum(c)
        ax.set_title(f"{TAGS[t]}\nlegit {n_own} / accidental {n_cross}")
    axes[0].set_ylabel("found effectors")
    axes[1].legend(frameon=False, loc="upper right", fontsize=8)
    fig.suptitle(
        "Why each found effector was emitted: most via their own system type; only 6 are accidental cross-type",
        fontweight="bold",
        y=1.02,
    )
    _save(fig, "03_emission_basis.png")


def fig04_audit():
    # Mutually-exclusive partition of the 19 emitted-but-unreachable@7 effectors (discordant_audit.md)
    cats = {
        "Fully sound\n(resolves + correct)": ["TseM", "TseZ", "Tlde1A", "BipB", "BipC"],
        "Wrong / 404 DOI\n(biology sound)": ["celA", "plaA", "VirA", "CopN", "Tle1", "Tae4_Stm"],
        "Unidentifiable /\nunsupported row": ["EFF00142", "EFF00150", "TseA_T6SS1", "ChlaDub1"],
        "Misassigned\nSS type": ["BopA", "BopE"],
        "Duplicate row": ["Tle4", "TplE_alias_Tle4"],
    }
    labels = list(cats)
    counts = [len(v) for v in cats.values()]
    colors = ["#3F8E8C", "#C9A227", "#CC7A30", "#A93232", "#7A7A7A"]
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    y = np.arange(len(labels))[::-1]
    ax.barh(y, counts, color=colors)
    for yi, c in zip(y, counts):
        ax.text(c + 0.1, yi, str(c), va="center", fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("effectors (of 19 emitted-but-unreachable@±7 audited)")
    ax.set_title("Literature audit: only 5/19 rows fully sound; 14 carry a defect")
    ax.set_xlim(right=max(counts) + 1)
    _save(fig, "04_discordant_audit.png")


_INDEX = []


def _save(fig, name):
    fig.tight_layout()
    out = FIGDIR / name
    fig.savefig(out)
    plt.close(fig)
    _INDEX.append(name)


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    for stale in FIGDIR.glob("*.png"):
        if not stale.name[:2].isdigit():
            stale.unlink()
    fig01_funnel()
    fig02_per_type()
    fig03_emission_basis()
    fig04_audit()
    print("Figure index:")
    for n in _INDEX:
        print(f"  {n}")


if __name__ == "__main__":
    main()
