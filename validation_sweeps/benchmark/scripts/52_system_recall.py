#!/usr/bin/env python3
"""System-level recall: per SS type, how many system INSTANCES had >=1 of their secreted proteins found.

The summary/01 recall figure counts secreted proteins (effectors). This variant counts SYSTEM INSTANCES:
a system is "found" if ssign emitted at least one of its (citation-verified) effectors. Same three-way
split as fig01, applied at the instance level:
  found       - the instance has >=1 effector emitted by ssign
  reach_miss  - 0 emitted, but >=1 effector is reachable @+/-3 (proximity could have reached it)
  unreach     - 0 emitted, every effector's machinery is >3 genes away

Uses the same cleaned answer key + T1SS staging fix as the other figures (clean_dataset). Effector ->
system instance comes from ceiling_per_effector.instance_id (joined by effector_locus, then uniprot).
T3SS uses the detection-enabled tag; the others use default (mirrors fig01).

Output: data/phase2/figures/summary/06_recall_systems.png
Run   : <repo>/.venv/bin/python scripts/52_system_recall.py
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import clean_dataset  # noqa: E402
from bench_io import read_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
P2 = BENCH / "data" / "phase2"
FIGDIR = P2 / "figures" / "summary"
CEILING = BENCH / "data" / "phase1" / "ceiling_per_effector.tsv"
TYPES = ["T1SS", "T2SS", "T3SS", "T4SS", "T6SS"]
C = {"found": "#2E7D6F", "reach_miss": "#E0A33B", "unreach": "#6C8EAD"}
plt.rcParams.update(
    {
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#444444",
    }
)


def _instance_lookup():
    cei = read_tsv(CEILING)
    by_loc = {r["effector_locus"]: r.get("instance_id", "") for r in cei if r.get("effector_locus", "").strip()}
    by_uni = {
        r["uniprot"]: r.get("instance_id", "") for r in cei if r.get("uniprot", "").strip() and r["uniprot"] != "-"
    }
    return by_loc, by_uni


def system_tab(tag, by_loc, by_uni, dropped_id):
    """ss -> {found, reach_miss, unreach} counting testable system instances."""
    rows = clean_dataset.load_clean_actual(P2 / f"actual_per_effector.{tag}.tsv")
    inst = defaultdict(lambda: {"found": False, "reach": False, "testable": False})
    for r in rows:
        if (r["gene"], r["uniprot"]) in dropped_id or (r["gene"], r.get("effector_locus", "")) in dropped_id:
            continue  # effector quarantined from the answer key
        iid = by_loc.get(r.get("effector_locus", "")) or by_uni.get(r.get("uniprot", ""))
        if not iid:
            continue
        v = inst[(r["ss_type"], iid)]
        if r["ssign_call"] == "emitted_secreted":
            v["found"] = True
        if r["testable"] == "yes":
            v["testable"] = True
            if r["reachable_n3"] == "true":
                v["reach"] = True
    out = {ss: {"found": 0, "reach_miss": 0, "unreach": 0} for ss in TYPES}
    for (ss, _iid), v in inst.items():
        if ss not in TYPES or not v["testable"]:
            continue
        if v["found"]:
            out[ss]["found"] += 1
        elif v["reach"]:
            out[ss]["reach_miss"] += 1
        else:
            out[ss]["unreach"] += 1
    return out


PLOT = ["T1SS", "T2SS", "T3SS", "T4SS", "T5SS", "T6SS"]
SEGS = [
    ("found", "found by ssign (>=1 secreted protein)"),
    ("reach_miss", "reachable, none found"),
    ("unreach", "unreachable (machinery >3 genes away)"),
]
# shorter wording for the poster key so the legend row is not too wide
SEG_LABELS_POSTER = {
    "found": "found by ssign",
    "reach_miss": "reachable, not found",
    "unreach": "unreachable (>3 genes)",
}


def render(tab, reachable, totals, tot_f, tot_r, outname, poster):
    """Stacked-bar system-recall figure. poster=True enlarges all type for a poster."""
    sz = (
        dict(fig=(14, 8.5), title=30, ylab=24, xtick=22, seg=20, ann=19, leg=17, barw=0.66)
        if poster
        else dict(fig=(9.8, 5.0), title=12, ylab=10, xtick=10, seg=9, ann=8.5, leg=8.5, barw=0.62)
    )
    fig, ax = plt.subplots(figsize=sz["fig"])
    x = np.arange(len(PLOT))
    bottom = np.zeros(len(PLOT))
    for key, lab in SEGS:
        vals = np.array([tab[ss][key] for ss in PLOT])
        ax.bar(x, vals, bottom=bottom, color=C[key], label=SEG_LABELS_POSTER[key] if poster else lab, width=sz["barw"])
        for xi, v, b in zip(x, vals, bottom):
            if v >= 2:
                ax.text(
                    xi,
                    b + v / 2,
                    str(v),
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=sz["seg"],
                    fontweight="bold",
                )
        bottom += vals
    for xi, ss in zip(x, PLOT):
        ax.text(
            xi,
            totals[ss] + 0.3,
            f"found {tab[ss]['found']}/{reachable[ss]}",
            ha="center",
            va="bottom",
            fontsize=sz["ann"],
        )
    ax.set_xticks(x)
    ax.set_xticklabels([f"{ss}\n(n={totals[ss]})" for ss in PLOT], fontsize=sz["xtick"])
    ax.tick_params(axis="y", labelsize=sz["xtick"])
    ax.set_ylabel("system instances", fontsize=sz["ylab"])
    ax.set_ylim(top=max(totals.values()) * 1.25)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=sz["leg"])
    ax.set_title(
        f"Ssign found {tot_f} out of {tot_r} reachable system instances", fontsize=sz["title"], fontweight="bold"
    )
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / outname, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outname}")


def main() -> int:
    by_loc, by_uni = _instance_lookup()
    dropped_id = clean_dataset.dropped_id()
    dft = system_tab("panel_genbank_default", by_loc, by_uni, dropped_id)
    t3t = system_tab("panel_genbank_t3ss", by_loc, by_uni, dropped_id)
    tab = {ss: (t3t if ss == "T3SS" else dft)[ss] for ss in TYPES}

    # T5SS is self-secreting -> not in the proximity tables; fold in the assembled T5SS recall (53).
    t5 = {"found": 0, "reach_miss": 0, "unreach": 0}
    for r in read_tsv(P2 / "t5ss_system_recall.tsv"):
        t5[r["status"]] += 1
    tab["T5SS"] = t5
    reachable = {ss: tab[ss]["found"] + tab[ss]["reach_miss"] for ss in PLOT}
    totals = {ss: sum(tab[ss].values()) for ss in PLOT}
    tot_f = sum(tab[ss]["found"] for ss in PLOT)
    tot_r = sum(reachable[ss] for ss in PLOT)

    render(tab, reachable, totals, tot_f, tot_r, "06_recall_systems.png", poster=False)
    render(tab, reachable, totals, tot_f, tot_r, "06_recall_systems_poster.png", poster=True)
    print(f"ssign found {tot_f}/{tot_r} reachable systems")
    for ss in PLOT:
        print(f"  {ss}: found {tab[ss]['found']}/{reachable[ss]} reachable (testable {totals[ss]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
