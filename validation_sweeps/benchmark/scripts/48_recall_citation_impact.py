#!/usr/bin/env python3
"""Recall figures after the full-table citation audit: with vs without the 121 unverifiable rows.

An effector's testable/reachable/found status is a property of its genome placement, not its citation,
so we do NOT re-run the placement pipeline. We attach each actual_per_effector row's citation_trust
(from positives_all + deepverify_removed), drop the audit-removed rows, and recompute recall over two
literature-effector populations:
  cleaned_with    = verified + unverifiable + fallback   (the kept 458)
  cleaned_without = verified + fallback                  (drop the 121 unverifiable -> 337)
pre_audit (all rows, the as-shipped denominator) is shown for context.

Recall categories mirror 32_summary_figures.recall_tab exactly (found / reachable@3-missed /
unreachable@3), testable rows only. T3SS uses the detection-enabled tag; the rest use default.

Outputs: data/phase2/figures/citation_impact/01_recall_with_unverifiable.png
         data/phase2/figures/citation_impact/02_recall_without_unverifiable.png
         data/phase2/figures/citation_impact/03_recall_impact.png
Run    : <repo>/.venv/bin/python scripts/48_recall_citation_impact.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
P2 = BENCH / "data" / "phase2"
DATASET = BENCH / "data" / "dataset"
FIGDIR = P2 / "figures" / "citation_impact"
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


def _real(u: str) -> str:
    u = (u or "").strip()
    return u if u and u != "-" else ""


def trust_lookups():
    """Per-effector citation class, keyed by uniprot / locus_tag / (gene, ss). Kept rows carry their
    tier; audit-removed rows carry 'removed'."""
    by_uni, by_locus, by_gene_ss = {}, {}, {}

    def add(rows, tier_fn):
        for r in rows:
            t = tier_fn(r)
            if u := _real(r.get("uniprot")):
                by_uni.setdefault(u, t)
            if lt := (r.get("locus_tag") or "").strip():
                by_locus.setdefault(lt, t)
            if g := (r.get("gene") or "").strip().lower():
                by_gene_ss.setdefault((g, (r.get("ss_type") or "").strip()), t)

    add(read_tsv(DATASET / "positives_all.tsv"), lambda r: r["citation_trust"])
    add(read_tsv(DATASET / "deepverify_removed.tsv"), lambda r: "removed")
    return by_uni, by_locus, by_gene_ss


def classify(tier: str) -> str:
    if tier.startswith("verified") or tier.startswith("fallback"):
        return "verified"
    if tier == "unverifiable":
        return "unverifiable"
    if tier == "removed":
        return "removed"
    return "unknown"  # 4 rows that don't join; kept in every population


def effector_class(row, lk) -> str:
    by_uni, by_locus, by_gene_ss = lk
    u, lt = _real(row.get("uniprot")), (row.get("effector_locus") or "").strip()
    g, ss = (row.get("gene") or "").strip().lower(), (row.get("ss_type") or "").strip()
    tier = by_uni.get(u) or by_locus.get(lt) or by_gene_ss.get((g, ss)) or "unknown"
    return classify(tier)


def recall_tab(tag, lk, keep):
    """ss -> {found, reach_miss, unreach} over rows whose citation class is in `keep`."""
    a = read_tsv(P2 / f"actual_per_effector.{tag}.tsv")
    out = {}
    for ss in TYPES:
        rows = [r for r in a if r["ss_type"] == ss and effector_class(r, lk) in keep]
        nf = [r for r in rows if r["ssign_call"] != "emitted_secreted"]
        out[ss] = {
            "found": sum(r["ssign_call"] == "emitted_secreted" for r in rows),
            "reach_miss": sum(r["testable"] == "yes" and r["reachable_n3"] == "true" for r in nf),
            "unreach": sum(r["testable"] == "yes" and r["reachable_n3"] != "true" for r in nf),
        }
    return out


def merged(dft, t3t):
    """T3SS from the detection-enabled tag, the rest from default (mirrors fig01)."""
    return {ss: (t3t if ss == "T3SS" else dft)[ss] for ss in TYPES}


def stacked(tab, title, fname):
    testable = {ss: sum(tab[ss].values()) for ss in TYPES}
    segs = [("found", "found by ssign"), ("reach_miss", "reachable @±3, missed"), ("unreach", "unreachable @±3")]
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    x = np.arange(len(TYPES))
    bottom = np.zeros(len(TYPES))
    for key, lab in segs:
        vals = np.array([tab[ss][key] for ss in TYPES])
        ax.bar(x, vals, bottom=bottom, color=C[key], label=lab, width=0.62)
        for xi, v, b in zip(x, vals, bottom):
            if v >= 8:
                ax.text(xi, b + v / 2, str(v), ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
        bottom += vals
    for xi, ss in zip(x, TYPES):
        reach = tab[ss]["found"] + tab[ss]["reach_miss"]
        ax.text(xi, testable[ss] + 2, f"found {tab[ss]['found']}/{reach}", ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{ss}*\n(n={testable[ss]})" if ss == "T3SS" else f"{ss}\n(n={testable[ss]})" for ss in TYPES])
    ax.set_ylabel("testable effectors")
    ax.set_ylim(top=max(testable.values()) * 1.25)
    tot_f = sum(tab[ss]["found"] for ss in TYPES)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8.5)
    ax.set_title(f"{title}\nssign found {tot_f} of {sum(testable.values())} testable effectors")
    ax.text(
        0.0,
        -0.24,
        "*T3SS shown with detection enabled (off by default).",
        transform=ax.transAxes,
        fontsize=7.5,
        color="#666",
    )
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / fname)
    plt.close(fig)
    print(f"wrote {fname}")


def impact(tabs, fname):
    """Grouped bars: recall % per ss type for each population."""
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    x = np.arange(len(TYPES))
    w = 0.26
    cols = {"pre_audit": "#BBBBBB", "cleaned_with": "#6C8EAD", "cleaned_without": "#2E7D6F"}
    for i, (name, tab) in enumerate(tabs.items()):
        pct, lab = [], []
        for ss in TYPES:
            testable = sum(tab[ss].values())
            r = 100 * tab[ss]["found"] / testable if testable else 0
            pct.append(r)
            lab.append(f"{tab[ss]['found']}/{testable}")
        bars = ax.bar(x + (i - 1) * w, pct, w, color=cols[name], label=name)
        for b, t in zip(bars, lab):
            ax.text(
                b.get_x() + b.get_width() / 2,
                b.get_height() + 1,
                t,
                ha="center",
                va="bottom",
                fontsize=6.5,
                rotation=90,
            )
    ax.set_xticks(x)
    ax.set_xticklabels(TYPES)
    ax.set_ylabel("recall (% of testable found)")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax.set_title("Recall per SS type before audit, vs after audit with / without the 121 unverifiable")
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / fname)
    plt.close(fig)
    print(f"wrote {fname}")


def main() -> int:
    lk = trust_lookups()
    dft, t3t = "panel_genbank_default", "panel_genbank_t3ss"
    pops = {
        "pre_audit": {"verified", "unverifiable", "removed", "unknown"},
        "cleaned_with": {"verified", "unverifiable", "unknown"},
        "cleaned_without": {"verified", "unknown"},
    }
    tabs = {name: merged(recall_tab(dft, lk, keep), recall_tab(t3t, lk, keep)) for name, keep in pops.items()}

    stacked(
        tabs["cleaned_with"],
        "Recall ±3, audited set incl. 121 unverifiable (458 rows)",
        "01_recall_with_unverifiable.png",
    )
    stacked(
        tabs["cleaned_without"],
        "Recall ±3, verified-only (121 unverifiable dropped, 337 rows)",
        "02_recall_without_unverifiable.png",
    )
    impact(tabs, "03_recall_impact.png")

    print("\nrecall (found / testable) per population:")
    hdr = f"{'ss':6s}" + "".join(f"{n:>20s}" for n in tabs)
    print(hdr)
    for ss in TYPES:
        line = f"{ss:6s}"
        for tab in tabs.values():
            t = sum(tab[ss].values())
            line += f"{tab[ss]['found']:>8d}/{t:<4d} ({100 * tab[ss]['found'] // t if t else 0:>3d}%)"
        print(line)
    line = f"{'ALL':6s}"
    for tab in tabs.values():
        f = sum(tab[ss]["found"] for ss in TYPES)
        t = sum(sum(tab[ss].values()) for ss in TYPES)
        line += f"{f:>8d}/{t:<4d} ({100 * f // t if t else 0:>3d}%)"
    print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
