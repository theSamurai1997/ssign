#!/usr/bin/env python3
"""Phase 2: true per-system PLM-E enrichment + null-size sweep on PAO1.

The job-3013556 run computed DLP/DSE enrichment only and wrote its
neighborhood files to job scratch (now gone), so per-system PLM-E
enrichment was never calculated. This reproduces the cheap detection
steps locally (~20 s macsyfinder, no GPU), reuses the run's whole-genome
predictions (from results_raw.csv), and runs the real enrichment test:

  - PLM-E gated at max_prob >= 0.8 (consistent with DLP/DSE), the rule
    Teo asked to make permanent.
  - null background swept over 200 / 1000 / all-non-neighborhood proteins,
    to show how the background estimate (and significance) converges.

This is the record of the analysis that justified dropping PLM-E from the
enrichment test. It was run against the PRE-change enrichment CLI (which still
accepted --plme); post-change, enrichment_testing.py no longer takes --plme, so
this driver is kept as documentation rather than a re-runnable script.

    .venv/bin/python 02_per_system_plme_enrichment.py
"""

from __future__ import annotations

import csv
import os
import random
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
SCRIPTS = os.path.join(REPO, "src", "ssign_app", "scripts")
PY = os.path.join(REPO, ".venv", "bin", "python")
GBFF = os.path.join(REPO, "validation_sweeps", "benchmark", "inputs_gb", "NC_002516.2.gbff")
RAW = os.path.join(HERE, "data", "NC_002516.2_results_raw.csv")
WORK = os.path.join(HERE, "phase2_work")
FIGS = os.path.join(HERE, "figures")
SAMPLE = "NC_002516.2"
CONF = 0.8
PLME_GATE = 0.8  # PLM-E max_prob gate (Teo: consistent with DLP/DSE)
NULL_SIZES = [200, 1000]  # plus "all" (whole non-neighborhood pool), added at runtime
SEED = 42

THEME = {
    "tool_colors": {"DLP": "#3F8E8C", "DSE": "#E0884B", "PLME": "#7A5C9E"},
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


def run(script: str, args: list[str], **kw) -> None:
    """Run an ssign pipeline script with the repo venv; raise on failure."""
    cmd = [PY, os.path.join(SCRIPTS, script), *args]
    r = subprocess.run(cmd, cwd=SCRIPTS, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        raise SystemExit(f"{script} failed (rc={r.returncode}):\n{r.stderr[-1500:]}")


def run_macsyfinder(proteins: str, out_dir: str) -> None:
    import shutil

    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    venv_bin = os.path.join(REPO, ".venv", "bin")
    # macsyfinder + its internal hmmsearch shim both live in the venv bin, which
    # isn't on the inherited PATH; prepend it so both resolve.
    env = {**os.environ, "PATH": venv_bin + os.pathsep + os.environ.get("PATH", "")}
    cmd = [
        os.path.join(venv_bin, "macsyfinder"),
        "--sequence-db",
        proteins,
        "--db-type",
        "ordered_replicon",
        "--models",
        "TXSScan",
        "all",
        "--out-dir",
        out_dir,
        "-w",
        "8",
        "--mute",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, env=env)
    if r.returncode != 0:
        raise SystemExit(f"macsyfinder failed (rc={r.returncode}):\n{r.stderr[-1500:]}")


def detect() -> dict:
    """Reproduce input-processing + detection + neighborhood; return file paths."""
    os.makedirs(WORK, exist_ok=True)
    f = {
        "proteins": os.path.join(WORK, "proteins.faa"),
        "gene_info": os.path.join(WORK, "gene_info.tsv"),
        "gene_order": os.path.join(WORK, "gene_order.tsv"),
        "msf_out": os.path.join(WORK, "macsyfinder_out"),
        "ss_components": os.path.join(WORK, "ss_components.tsv"),
        "valid_systems": os.path.join(WORK, "valid_systems.tsv"),
        "neigh_fasta": os.path.join(WORK, "neighborhood.faa"),
        "neigh_ids": os.path.join(WORK, "neighborhood_ids.tsv"),
    }
    print("  [1/5] extract proteins + gene info from GenBank ...")
    run(
        "extract_proteins.py",
        ["--input", GBFF, "--sample", SAMPLE, "--out-proteins", f["proteins"], "--out-gene-info", f["gene_info"]],
    )
    print("  [2/5] extract gene order ...")
    run("extract_gene_order.py", ["--gene-info", f["gene_info"], "--output", f["gene_order"]])
    print("  [3/5] macsyfinder (TXSScan, ordered_replicon) ...")
    run_macsyfinder(f["proteins"], f["msf_out"])
    print("  [4/5] validate systems -> ss_components (excluded: Flagellum,Tad; wholeness 0.8) ...")
    run(
        "validate_macsyfinder_systems.py",
        [
            "--msf-dir",
            f["msf_out"],
            "--gene-info",
            f["gene_info"],
            "--sample",
            SAMPLE,
            "--wholeness-threshold",
            "0.8",
            "--excluded-systems",
            "Flagellum,Tad",
            "--out-components",
            f["ss_components"],
            "--out-systems",
            f["valid_systems"],
        ],
    )
    print("  [5/5] extract neighborhood (window 3) ...")
    run(
        "extract_neighborhood.py",
        [
            "--gene-order",
            f["gene_order"],
            "--ss-components",
            f["ss_components"],
            "--proteins",
            f["proteins"],
            "--window",
            "3",
            "--output",
            f["neigh_fasta"],
            "--output-ids",
            f["neigh_ids"],
        ],
    )
    return f


def synth_prediction_tsvs(raw_rows: list[dict]) -> dict:
    """Write dlp/dse/plme TSVs from the run's whole-genome predictions.

    PLM-E passes_threshold is rewritten as (max_prob >= PLME_GATE) so the
    enrichment test's is_plme_positive applies the 0.8 confidence gate.
    """
    dlp = os.path.join(WORK, "dlp.tsv")
    dse = os.path.join(WORK, "dse.tsv")
    plme = os.path.join(WORK, "plme.tsv")
    with open(dlp, "w", newline="") as fd:
        w = csv.writer(fd, delimiter="\t")
        w.writerow(["locus_tag", "dlp_extracellular_prob"])
        for r in raw_rows:
            w.writerow([r["locus_tag"], r.get("dlp_extracellular_prob", "")])
    with open(dse, "w", newline="") as fd:
        w = csv.writer(fd, delimiter="\t")
        w.writerow(["locus_tag", "dse_ss_type", "dse_max_prob"])
        for r in raw_rows:
            w.writerow([r["locus_tag"], r.get("dse_ss_type", ""), r.get("dse_max_prob", "")])
    with open(plme, "w", newline="") as fd:
        w = csv.writer(fd, delimiter="\t")
        w.writerow(["locus_tag", "passes_threshold", "plm_effector_type", "plm_effector_max_prob"])
        for r in raw_rows:
            try:
                mp = float(r.get("plm_effector_max_prob") or 0)
            except ValueError:
                mp = 0.0
            gated = "1" if mp >= PLME_GATE else "0"
            w.writerow([r["locus_tag"], gated, r.get("plm_effector_type", ""), r.get("plm_effector_max_prob", "")])
    return {"dlp": dlp, "dse": dse, "plme": plme}


def null_pool(raw_rows: list[dict], f: dict) -> list[str]:
    """All locus_tags that are neither SS components nor in any neighborhood."""
    neigh = {line.strip() for line in open(f["neigh_ids"]) if line.strip()}
    comps = set()
    with open(f["ss_components"]) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("locus_tag"):
                comps.add(row["locus_tag"])
    excluded = neigh | comps
    return [r["locus_tag"] for r in raw_rows if r["locus_tag"] not in excluded]


def run_enrichment(f: dict, preds: dict, null_ids: list[str], tag: str) -> list[dict]:
    null_path = os.path.join(WORK, f"null_{tag}.txt")
    with open(null_path, "w") as fh:
        fh.write("\n".join(null_ids) + "\n")
    out = os.path.join(WORK, f"enrich_{tag}.tsv")
    run(
        "enrichment_testing.py",
        [
            "--ss-components",
            f["ss_components"],
            "--gene-order",
            f["gene_order"],
            "--dlp",
            preds["dlp"],
            "--dse",
            preds["dse"],
            "--plme",
            preds["plme"],
            "--null-ids",
            null_path,
            "--window",
            "3",
            "--conf-threshold",
            str(CONF),
            "--sample",
            SAMPLE,
            "--out",
            out,
        ],
    )
    with open(out) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> None:
    if not os.path.exists(GBFF):
        raise SystemExit(f"missing genome: {GBFF}")
    raw_rows = list(csv.DictReader(open(RAW)))
    print(f"PAO1 {SAMPLE}: {len(raw_rows)} proteins; reproducing detection ...")
    f = detect()
    preds = synth_prediction_tsvs(raw_rows)
    pool = null_pool(raw_rows, f)
    print(f"\nnon-neighborhood null pool: {len(pool)} proteins")

    rng = random.Random(SEED)
    sizes = [(str(n), rng.sample(pool, n)) for n in NULL_SIZES] + [("all", pool)]
    results = {tag: run_enrichment(f, preds, ids, tag) for tag, ids in sizes}

    # --- report: background + significance per tool per null size ---
    def pbg(rows, tool):
        return next((float(r["p_bg"]) for r in rows if r["tool"] == tool), float("nan"))

    def n_sig(rows, tool):
        return sum(1 for r in rows if r["tool"] == tool and r["significant"] == "True")

    print("\n=== background p_bg by null size (per tool) ===")
    print(f"  {'tool':5s} " + "  ".join(f"{tag:>9s}" for tag, _ in sizes))
    for tool in ("DLP", "DSE", "PLME"):
        print(f"  {tool:5s} " + "  ".join(f"{pbg(results[tag], tool):9.4f}" for tag, _ in sizes))
    print("\n=== significant (scope x tool) calls @ q<0.05 by null size ===")
    print(f"  {'tool':5s} " + "  ".join(f"{tag:>9s}" for tag, _ in sizes))
    for tool in ("DLP", "DSE", "PLME"):
        print(f"  {tool:5s} " + "  ".join(f"{n_sig(results[tag], tool):9d}" for tag, _ in sizes))

    # --- per-system PLME at the 1000-null background ---
    base = results["1000"]
    plme_rows = [r for r in base if r["tool"] == "PLME" and r["scope_kind"] == "system"]
    plme_rows.sort(key=lambda r: float(r["pvalue"]))
    print(f"\n=== per-system PLM-E (max_prob>=0.8 gate, 1000-null bg={pbg(base, 'PLME'):.3f}) ===")
    print(f"  {'scope_id':34s} {'M':>3s} {'k':>3s} {'fold':>6s} {'pvalue':>9s} {'qvalue':>9s}  sig")
    for r in plme_rows:
        print(
            f"  {r['scope_id']:34s} {r['M']:>3s} {r['k']:>3s} {str(r['fold_enrich']):>6s} "
            f"{float(r['pvalue']):9.4f} {float(r['qvalue']):9.4f}  {r['significant']}"
        )
    n_plme_sig = sum(1 for r in plme_rows if r["significant"] == "True")
    print(f"  -> {n_plme_sig}/{len(plme_rows)} systems PLM-E-enriched at q<0.05")

    make_figures(sizes, results)
    print("\nFigure index:")
    print("  04  04_null_size_sweep.png         — background + significant calls vs null size")
    print("  05  05_plme_per_system.png         — per-system PLM-E k/M and significance (0.8 gate)")


def make_figures(sizes, results) -> None:
    tags = [t for t, _ in sizes]
    xs = list(range(len(tags)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for tool in ("DLP", "DSE", "PLME"):
        ys = [next((float(r["p_bg"]) for r in results[t] if r["tool"] == tool), np.nan) * 100 for t in tags]
        ax1.plot(xs, ys, "-o", color=THEME["tool_colors"][tool], label=tool)
    ax1.set_xticks(xs)
    ax1.set_xticklabels(tags)
    ax1.set_xlabel("null sample size")
    ax1.set_ylabel("background positive rate (%)")
    ax1.set_title("Background estimate converges as null grows")
    ax1.legend(frameon=False, fontsize=8)
    for tool in ("DLP", "DSE", "PLME"):
        ys = [sum(1 for r in results[t] if r["tool"] == tool and r["significant"] == "True") for t in tags]
        ax2.plot(xs, ys, "-o", color=THEME["tool_colors"][tool], label=tool)
    ax2.set_xticks(xs)
    ax2.set_xticklabels(tags)
    ax2.set_xlabel("null sample size")
    ax2.set_ylabel("significant (scope x tool), q<0.05")
    ax2.set_title("Significant calls vs null size")
    ax2.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "04_null_size_sweep.png"))
    plt.close(fig)

    base = results["1000"]
    plme = [r for r in base if r["tool"] == "PLME" and r["scope_kind"] == "system"]
    plme.sort(key=lambda r: float(r["k"]) / float(r["M"]), reverse=True)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.32 * len(plme) + 1)))
    pbg = next((float(r["p_bg"]) for r in base if r["tool"] == "PLME"), np.nan)
    y = list(range(len(plme)))
    fracs = [float(r["k"]) / float(r["M"]) for r in plme]
    colors = ["#A93232" if r["significant"] == "True" else THEME["tool_colors"]["PLME"] for r in plme]
    ax.barh(y, fracs, color=colors)
    ax.axvline(pbg, ls="--", color=THEME["ref_line"], lw=0.8, label=f"PLM-E background {pbg:.2f}")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['ss_type']} ({r['k']}/{r['M']})" for r in plme], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("PLM-E-positive fraction of neighborhood (k/M)")
    ax.set_title("Per-system PLM-E density vs background (red = significant, q<0.05)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGS, "05_plme_per_system.png"))
    plt.close(fig)


if __name__ == "__main__":
    main()
