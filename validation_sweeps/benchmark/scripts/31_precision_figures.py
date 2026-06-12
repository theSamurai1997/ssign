#!/usr/bin/env python3
"""Phase 2 task 8.4: figures for the deterministic precision estimate.

Joins the two precision tiers (29 DB-confirmed floor + 30 annotation buckets) per emission and renders:
  01  precision bounds      : per SS type, floor (DB+gold confirmed) -> ceiling (1 - obvious-FP - apparatus)
  02  emission composition  : per SS type, the annotation-bucket stack (effector / hypothetical / other /
                              apparatus / housekeeping)
  03  overall picture       : the whole proximity emission set as one stack, with the DB floor and the
                              soft ceiling marked, showing the large unresolvable middle

Run on the as-shipped default tag (T3SS excluded). Inputs are run-independent given the run tables.
Inputs : data/phase2/emissions_{dbmatch,fpclass}.<tag>.tsv
Output : data/phase2/figures/precision/01..03_*.png
Run:     <repo>/.venv/bin/python scripts/31_precision_figures.py
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
P2 = BENCH / "data" / "phase2"
FIGDIR = P2 / "figures" / "precision"
TAG = "panel_genbank_default"
# Per-type figures cover these only (a type absent here is dropped from figs 01/02 but still counted
# in the all-proximity fig 03); keep in sync with the panel's nearby_ss_types.
SS_ORDER = ["T5aSS", "T6SSi", "T5bSS", "T1SS", "T5cSS", "T4aP", "T2SS"]
BUCKETS = ["effector", "hypothetical", "other", "apparatus", "housekeeping"]

THEME = {
    "effector": "#3F8E8C",
    "hypothetical": "#C9C2B6",
    "other": "#9AA8B5",
    "apparatus": "#6C8EAD",
    "housekeeping": "#A93232",
    "floor": "#234E70",
    "ceiling": "#3F8E8C",
    "ref": "#444444",
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
    }
)
_INDEX = []


def load():
    fp = {(r["unit_id"], r["locus_tag"]): r["fp_class"] for r in read_tsv(P2 / f"emissions_fpclass.{TAG}.tsv")}
    rows = []
    for r in read_tsv(P2 / f"emissions_dbmatch.{TAG}.tsv"):
        if r["substrate_source"] != "proximity":
            continue
        rows.append(
            {
                "ss": (r["nearby_ss_types"] or "?").split(";")[0].split(",")[0],
                "is_gold": r["is_gold"] == "yes",
                "confirmed": r["db_class"] == "confirmed",
                "fp": fp.get((r["unit_id"], r["locus_tag"]), "other"),
            }
        )
    return rows


def per_type(rows):
    """ss -> dict(n, floor, ceiling, bucket counts)."""
    out = {}
    for ss in SS_ORDER:
        sub = [r for r in rows if r["ss"] == ss]
        if not sub:
            continue
        n = len(sub)
        bc = Counter(r["fp"] for r in sub)
        hk_fp = sum(r["fp"] == "housekeeping" and not r["is_gold"] for r in sub)
        floor = sum(r["is_gold"] or r["confirmed"] for r in sub) / n
        ceiling = (n - hk_fp - bc["apparatus"]) / n
        out[ss] = {"n": n, "floor": floor, "ceiling": ceiling, "bc": bc}
    return out


def fig01_bounds(pt):
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    types = list(pt)
    y = np.arange(len(types))[::-1]
    for yi, ss in zip(y, types):
        f, c = pt[ss]["floor"], pt[ss]["ceiling"]
        ax.plot([f, c], [yi, yi], color="#BBBBBB", lw=6, solid_capstyle="round", zorder=1)
        ax.scatter([f], [yi], color=THEME["floor"], s=60, zorder=3)
        ax.scatter([c], [yi], color=THEME["ceiling"], s=60, zorder=3)
        ax.text(f, yi + 0.18, f"{f:.0%}", ha="center", va="bottom", fontsize=8, color=THEME["floor"])
        ax.text(c, yi + 0.18, f"{c:.0%}", ha="center", va="bottom", fontsize=8, color=THEME["ceiling"])
    ax.set_yticks(y)
    ax.set_yticklabels([f"{ss} (n={pt[ss]['n']})" for ss in types])
    ax.set_xlim(-0.02, 1.0)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    ax.set_xlabel("precision of proximity-called substrates")
    ax.scatter([], [], color=THEME["floor"], label="floor: DB/gold-confirmed (provable TP)")
    ax.scatter([], [], color=THEME["ceiling"], label="ceiling: not obviously non-secreted")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=8)
    ax.set_title("Precision is bounded, not pinned: a wide DB-floor → annotation-ceiling band")
    _save(fig, "01_precision_bounds.png")


def fig02_composition(pt):
    fig, ax = plt.subplots(figsize=(9, 4.8))
    types = list(pt)
    y = np.arange(len(types))[::-1]
    left = np.zeros(len(types))
    for b in BUCKETS:
        vals = np.array([pt[ss]["bc"][b] / pt[ss]["n"] for ss in types])
        ax.barh(y, vals, left=left, color=THEME[b], label=b)
        left += vals
    ax.set_yticks(y)
    ax.set_yticklabels([f"{ss} (n={pt[ss]['n']})" for ss in types])
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    ax.set_xlabel("share of proximity emissions")
    ax.legend(frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8)
    ax.set_title("Annotation composition: effector-like vs cytoplasmic-FP vs unresolvable middle")
    _save(fig, "02_emission_composition.png")


def fig03_overall(rows):
    n = len(rows)
    bc = Counter(r["fp"] for r in rows)
    hk_fp = sum(r["fp"] == "housekeeping" and not r["is_gold"] for r in rows)
    floor = sum(r["is_gold"] or r["confirmed"] for r in rows) / n
    ceiling = (n - hk_fp - bc["apparatus"]) / n
    fig, ax = plt.subplots(figsize=(9, 2.9))
    left = 0.0
    for b in BUCKETS:
        v = bc[b] / n
        ax.barh([0], [v], left=[left], color=THEME[b], label=f"{b} {v:.0%}")
        if v > 0.05:
            ax.text(left + v / 2, 0, f"{v:.0%}", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        left += v
    for x, lab, col in [
        (floor, f"DB floor {floor:.0%}", THEME["floor"]),
        (ceiling, f"ceiling {ceiling:.0%}", THEME["ref"]),
    ]:
        ax.axvline(x, color=col, lw=2, ls="--")
        ax.text(x, 0.62, lab, ha="center", va="bottom", fontsize=9, color=col, fontweight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, 0.9)
    ax.set_yticks([])
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0%}")
    ax.legend(frameon=False, ncol=5, loc="lower center", bbox_to_anchor=(0.5, -0.55), fontsize=8)
    ax.set_title(
        f"All {n} proximity substrate calls: ~{(bc['hypothetical'] + bc['other']) / n:.0%} unresolvable by DB or annotation"
    )
    _save(fig, "03_overall_precision.png")


def _save(fig, name):
    fig.tight_layout()
    fig.savefig(FIGDIR / name)
    plt.close(fig)
    _INDEX.append(name)


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    for stale in FIGDIR.glob("*.png"):
        if not stale.name[:2].isdigit():
            stale.unlink()
    rows = load()
    pt = per_type(rows)
    fig01_bounds(pt)
    fig02_composition(pt)
    fig03_overall(rows)
    print("Figure index:")
    for nme in _INDEX:
        print(f"  {nme}")


if __name__ == "__main__":
    main()
