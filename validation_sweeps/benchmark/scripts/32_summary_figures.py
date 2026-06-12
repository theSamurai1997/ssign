#!/usr/bin/env python3
"""Phase 2: the put-it-together summary figures (recall + the T3SS story + emission quality).

01  recall @ window 3 (as shipped)  : per SS type, of the literature-curated effectors, how many ssign
                                       FOUND, how many it COULD have (reachable within +/-3 of the
                                       machinery but missed), how many it could NEVER (machinery >3
                                       genes away = unreachable@3, or never put in front of ssign at
                                       all = non-testable: no genome, ORF absent, or no machinery anchor).
02  the T3SS story                   : T3SS effectors, default (excluded) vs included, showing both the
                                       exclusion AND that most are genome-dispersed (unreachable@3).
03  emission quality                 : of everything ssign emits (proximity calls), how many are
                                       reasonable / unresolvable / wrong, per assigned SS type.

Run on the as-shipped default tag; 02 also reads the T3SS-included tag. Inputs are run tables only.
Inputs : data/phase2/actual_per_effector.<tag>.tsv, emissions_{dbmatch,fpclass}.<tag>.tsv
Output : data/phase2/figures/summary/01..03_*.png
Run:     <repo>/.venv/bin/python scripts/32_summary_figures.py
"""

from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
P2 = BENCH / "data" / "phase2"
FIGDIR = P2 / "figures" / "summary"
TYPES = ["T1SS", "T2SS", "T3SS", "T4SS", "T6SS"]

C = {
    "found": "#2E7D6F",
    "reach_miss": "#E0A33B",
    "unreach": "#6C8EAD",
    "nontest": "#CDCDCD",
    "reasonable": "#2E7D6F",
    "unresolvable": "#C9C2B6",
    "wrong": "#A93232",
    "ref": "#333333",
}
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
_INDEX = []


def recall_tab(tag):
    """ss -> {total, found, reach_miss, unreach, nontest} at window 3."""
    a = read_tsv(P2 / f"actual_per_effector.{tag}.tsv")
    out = {}
    for ss in TYPES:
        rows = [r for r in a if r["ss_type"] == ss]
        found = sum(r["ssign_call"] == "emitted_secreted" for r in rows)
        nf = [r for r in rows if r["ssign_call"] != "emitted_secreted"]
        out[ss] = {
            "total": len(rows),
            "found": found,
            "reach_miss": sum(r["testable"] == "yes" and r["reachable_n3"] == "true" for r in nf),
            "unreach": sum(r["testable"] == "yes" and r["reachable_n3"] != "true" for r in nf),
            "nontest": sum(r["testable"] == "no" for r in nf),
        }
    return out


def fig01_recall(dft, t3t):
    # Testable effectors only (drop non-testable — they were never put in front of ssign). T3SS uses
    # the detection-enabled run; T1/T2/T4/T6 are identical across tags (only Flagellum/Tad/T3SS toggle).
    tab = {ss: (t3t if ss == "T3SS" else dft)[ss] for ss in TYPES}
    testable = {ss: tab[ss]["found"] + tab[ss]["reach_miss"] + tab[ss]["unreach"] for ss in TYPES}
    segs = [
        ("found", "found by ssign"),
        ("reach_miss", "reachable @±3, missed"),
        ("unreach", "unreachable @±3 (machinery >3 genes away)"),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    x = np.arange(len(TYPES))
    bottom = np.zeros(len(TYPES))
    for key, lab in segs:
        vals = np.array([tab[ss][key] for ss in TYPES])
        ax.bar(x, vals, bottom=bottom, color=C[key], label=lab, width=0.62)
        for xi, v, b in zip(x, vals, bottom):
            if v >= 10:  # tiny segments collide on short bars; the found/reachable annotation covers them
                ax.text(xi, b + v / 2, str(v), ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
        bottom += vals
    for xi, ss in zip(x, TYPES):
        reach = tab[ss]["found"] + tab[ss]["reach_miss"]
        ax.text(xi, testable[ss] + 3, f"found {tab[ss]['found']}/{reach}", ha="center", va="bottom", fontsize=8.5)
    labels = [f"{ss}*\n(n={testable[ss]})" if ss == "T3SS" else f"{ss}\n(n={testable[ss]})" for ss in TYPES]
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("testable effectors")
    ax.set_ylim(top=max(testable.values()) * 1.2)
    tot_f = sum(tab[ss]["found"] for ss in TYPES)
    tot_r = sum(tab[ss]["found"] + tab[ss]["reach_miss"] for ss in TYPES)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8.5)
    ax.set_title(
        f"Recall at proximity ±3: ssign found {tot_f} of {sum(testable.values())} testable effectors\n"
        f"(proximity could reach {tot_r}; the bar above each is found / reachable)"
    )
    ax.text(
        0.0,
        -0.26,
        "*T3SS shown with detection enabled (off by default).",
        transform=ax.transAxes,
        fontsize=7.5,
        color="#666",
    )
    _save(fig, "01_recall_window3.png")


def fig02_t3ss(dft, t3t):
    segs = [
        ("found", C["found"]),
        ("reach_miss", C["reach_miss"]),
        ("unreach", C["unreach"]),
        ("nontest", C["nontest"]),
    ]
    labels = ["default\n(T3SS excluded)", "T3SS included"]
    data = [dft["T3SS"], t3t["T3SS"]]
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    x = np.arange(2)
    bottom = np.zeros(2)
    for key, col in segs:
        vals = np.array([d[key] for d in data])
        ax.bar(x, vals, bottom=bottom, color=col, width=0.5)
        for xi, v, b in zip(x, vals, bottom):
            if v >= 5:
                ax.text(
                    xi,
                    b + v / 2,
                    str(v),
                    ha="center",
                    va="center",
                    color="white" if key != "nontest" else "#555",
                    fontsize=9,
                    fontweight="bold",
                )
        bottom += vals
    for xi, d in zip(x, data):
        ax.text(xi, d["total"] + 3, f"found {d['found']}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("T3SS effectors (n=237)")
    ax.set_ylim(top=260)
    ax.set_title("The T3SS story: excluded by default,\nand mostly out of proximity reach even when on")
    ax.text(
        0.5,
        -0.22,
        "The T3SS-on run works: MacSyFinder detects 30 injectisomes and found\nrises 3→15. But ~73% of T3SS effectors are genome-dispersed (unreachable\n@±3 of the injectisome), so proximity can't reach them. T3SS is OFF by\ndefault not because MacSyFinder can't find it — it can — but because\nDeepSecE floods T3SS with misclassified flagellar proteins.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8.2,
        color="#333",
    )
    _save(fig, "02_t3ss_story.png")


def emission_quality():
    db = read_tsv(P2 / "emissions_dbmatch.panel_genbank_default.tsv")
    fp = {
        (r["unit_id"], r["locus_tag"]): r["fp_class"]
        for r in read_tsv(P2 / "emissions_fpclass.panel_genbank_default.tsv")
    }

    def fam(s):
        s = (s or "?").split(";")[0].split(",")[0]
        if s.startswith("T6SS"):
            return "T6SS"
        if s.startswith("T5"):
            return "T5SS"
        return "T4SS" if s == "pT4SSt" else s

    cat = defaultdict(Counter)
    for r in db:
        if r["substrate_source"] != "proximity":
            continue
        cls = fp.get((r["unit_id"], r["locus_tag"]), "other")
        reasonable = r["db_class"] == "confirmed" or r["is_gold"] == "yes" or cls == "effector"
        bucket = "reasonable" if reasonable else ("wrong" if cls in ("housekeeping", "apparatus") else "unresolvable")
        cat[fam(r["nearby_ss_types"])][bucket] += 1
    return cat


def fig03_quality(cat):
    order = sorted((k for k in cat if sum(cat[k].values()) >= 20), key=lambda k: -sum(cat[k].values()))
    segs = [
        ("reasonable", "reasonable (DB/gold/effector annotation)"),
        ("unresolvable", "unresolvable (hypothetical / unclear)"),
        ("wrong", "wrong (cytoplasmic / machinery)"),
    ]
    fig, ax = plt.subplots(figsize=(9.5, 5.0))
    x = np.arange(len(order))
    bottom = np.zeros(len(order))
    for key, lab in segs:
        vals = np.array([cat[k][key] for k in order])
        ax.bar(x, vals, bottom=bottom, color=C[key], label=lab, width=0.62)
        for xi, v, b, k in zip(x, vals, bottom, order):
            if v / sum(cat[k].values()) >= 0.06:
                ax.text(
                    xi,
                    b + v / 2,
                    f"{v / sum(cat[k].values()):.0%}",
                    ha="center",
                    va="center",
                    color="white" if key != "unresolvable" else "#555",
                    fontsize=8.5,
                    fontweight="bold",
                )
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([f"{k}\n(n={sum(cat[k].values())})" for k in order])
    ax.set_ylabel("proximity substrate calls")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=8.2)
    tot = Counter()
    for k in cat:
        tot.update(cat[k])
    n = sum(tot.values())
    ax.set_title(
        f"Of everything ssign emits, how much is real? (proximity calls, by assigned type)\noverall: {tot['reasonable'] / n:.0%} reasonable · {tot['unresolvable'] / n:.0%} unresolvable · {tot['wrong'] / n:.0%} wrong  (n={n})"
    )
    _save(fig, "03_emission_quality.png")


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
    dft = recall_tab("panel_genbank_default")
    t3t = recall_tab("panel_genbank_t3ss")
    fig01_recall(dft, t3t)
    fig02_t3ss(dft, t3t)
    fig03_quality(emission_quality())
    print("Figure index:")
    for n in _INDEX:
        print(f"  {n}")


if __name__ == "__main__":
    main()
