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


def system_tab(tag, by_loc, by_uni):
    """ss -> {found, reach_miss, unreach} counting testable system instances."""
    rows = clean_dataset.load_clean_actual(P2 / f"actual_per_effector.{tag}.tsv")
    inst = defaultdict(lambda: {"found": False, "reach": False, "testable": False})
    for r in rows:
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


def main() -> int:
    by_loc, by_uni = _instance_lookup()
    dft = system_tab("panel_genbank_default", by_loc, by_uni)
    t3t = system_tab("panel_genbank_t3ss", by_loc, by_uni)
    tab = {ss: (t3t if ss == "T3SS" else dft)[ss] for ss in TYPES}

    # T5SS is self-secreting -> not in the proximity tables; fold in the assembled T5SS recall (53).
    t5 = {"found": 0, "reach_miss": 0, "unreach": 0}
    for r in read_tsv(P2 / "t5ss_system_recall.tsv"):
        t5[r["status"]] += 1
    tab["T5SS"] = t5
    PLOT = ["T1SS", "T2SS", "T3SS", "T4SS", "T5SS", "T6SS"]
    reachable = {ss: tab[ss]["found"] + tab[ss]["reach_miss"] for ss in PLOT}
    totals = {ss: sum(tab[ss].values()) for ss in PLOT}

    segs = [
        ("found", "found by ssign (>=1 secreted protein)"),
        ("reach_miss", "reachable, none found"),
        ("unreach", "unreachable (machinery >3 genes away)"),
    ]
    fig, ax = plt.subplots(figsize=(9.8, 5.0))
    x = np.arange(len(PLOT))
    bottom = np.zeros(len(PLOT))
    for key, lab in segs:
        vals = np.array([tab[ss][key] for ss in PLOT])
        ax.bar(x, vals, bottom=bottom, color=C[key], label=lab, width=0.62)
        for xi, v, b in zip(x, vals, bottom):
            if v >= 2:
                ax.text(xi, b + v / 2, str(v), ha="center", va="center", color="white", fontsize=9, fontweight="bold")
        bottom += vals
    for xi, ss in zip(x, PLOT):
        ax.text(
            xi, totals[ss] + 0.3, f"found {tab[ss]['found']}/{reachable[ss]}", ha="center", va="bottom", fontsize=8.5
        )
    ax.set_xticks(x)
    ax.set_xticklabels([f"{ss}\n(n={totals[ss]})" for ss in PLOT])
    ax.set_ylabel("system instances")
    ax.set_ylim(top=max(totals.values()) * 1.25)
    tot_f = sum(tab[ss]["found"] for ss in PLOT)
    tot_r = sum(reachable[ss] for ss in PLOT)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8.5)
    ax.set_title(f"Ssign found {tot_f} out of {tot_r} reachable system instances")
    ax.text(
        0.0,
        -0.24,
        "T5SS self-secreting (T5a/c/d/e autotransporter = its own substrate; T5b TpsA secreted by an adjacent TpsB).",
        transform=ax.transAxes,
        fontsize=7.5,
        color="#666",
    )
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / "06_recall_systems.png")
    plt.close(fig)
    print(f"wrote 06_recall_systems.png — ssign found {tot_f}/{tot_r} reachable systems")
    for ss in PLOT:
        print(f"  {ss}: found {tab[ss]['found']}/{reachable[ss]} reachable (testable {totals[ss]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
