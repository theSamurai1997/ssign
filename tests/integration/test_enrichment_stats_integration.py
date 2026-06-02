"""End-to-end integration test for the --enrichment-stats path.

Drives the three new pieces (sample_null_proteins.py, enrichment_testing.py,
pool_enrichment_stats helper) as subprocesses against two synthetic
two-contig genomes. Predictions (DLP/DSE) are hand-synthesised because
the stats engine just consumes their TSV outputs -- the actual tool
binaries aren't part of what's being tested here.

Mirrors the test_quick_scripts_integration.py style: small fixtures, no
external binaries / databases, runs in well under a second.
"""

import csv
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = PROJECT_ROOT / "src" / "ssign_app" / "scripts"


def _run_script(script_name: str, args: list[str]) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(SCRIPTS_DIR / script_name)] + args
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60)


def _write_fasta(path: Path, seqs: dict):
    with open(path, "w") as f:
        for sid, seq in seqs.items():
            f.write(f">{sid}\n{seq}\n")


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict]):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _build_genome(root: Path, prefix: str, positives_neigh: set, positives_null: set):
    """Build proteome + gene_order + ss_components for one synthetic genome.

    Layout: 15 proteins, contig_A genes 0-9 + contig_B genes 0-4. Two T2SS
    components on contig_A at gene_index 5, 6 (sys_id <prefix>_T2SS_1).
    Returns dict of file paths.
    """
    root.mkdir(parents=True, exist_ok=True)
    proteome = {f"{prefix}_GA_{i:04d}": "M" + ("A" * 30) for i in range(10)}
    proteome.update({f"{prefix}_GB_{i:04d}": "M" + ("L" * 30) for i in range(5)})
    proteins = root / "proteins.faa"
    _write_fasta(proteins, proteome)

    gene_order = root / "gene_order.tsv"
    rows = []
    for i in range(10):
        rows.append({"contig": f"{prefix}_contig_A", "gene_index": str(i), "locus_tag": f"{prefix}_GA_{i:04d}"})
    for i in range(5):
        rows.append({"contig": f"{prefix}_contig_B", "gene_index": str(i), "locus_tag": f"{prefix}_GB_{i:04d}"})
    _write_tsv(gene_order, ["contig", "gene_index", "locus_tag"], rows)

    ss_components = root / "ss_components.tsv"
    _write_tsv(
        ss_components,
        ["locus_tag", "ss_type", "sys_id", "excluded"],
        [
            {"locus_tag": f"{prefix}_GA_0005", "ss_type": "T2SS", "sys_id": f"{prefix}_T2SS_1", "excluded": "False"},
            {"locus_tag": f"{prefix}_GA_0006", "ss_type": "T2SS", "sys_id": f"{prefix}_T2SS_1", "excluded": "False"},
        ],
    )

    # Hand-synthesise DLP + DSE outputs over the full proteome.
    dlp_rows = []
    dse_rows = []
    for sid in proteome:
        pos = sid in positives_neigh or sid in positives_null
        dlp_rows.append({"locus_tag": sid, "dlp_extracellular_prob": "0.95" if pos else "0.10"})
        if pos:
            dse_rows.append({"locus_tag": sid, "dse_ss_type": "T2SS", "dse_max_prob": "0.95"})
        else:
            dse_rows.append({"locus_tag": sid, "dse_ss_type": "Non-secreted", "dse_max_prob": "0.10"})
    dlp_tsv = root / "dlp.tsv"
    dse_tsv = root / "dse.tsv"
    _write_tsv(dlp_tsv, ["locus_tag", "dlp_extracellular_prob"], dlp_rows)
    _write_tsv(dse_tsv, ["locus_tag", "dse_ss_type", "dse_max_prob"], dse_rows)

    return {
        "proteins": proteins,
        "gene_order": gene_order,
        "ss_components": ss_components,
        "dlp": dlp_tsv,
        "dse": dse_tsv,
    }


class TestStatsEndToEnd:
    def test_two_genome_run_emits_per_genome_and_pooled(self, tmp_path):
        # Genome A: 2 positives in T2SS neighborhood, 0 in null pool.
        ga = _build_genome(
            tmp_path / "ga",
            prefix="A",
            positives_neigh={"A_GA_0003", "A_GA_0004"},
            positives_null=set(),
        )
        # Genome B: 2 positives in T2SS neighborhood, 0 in null pool.
        gb = _build_genome(
            tmp_path / "gb",
            prefix="B",
            positives_neigh={"B_GA_0003", "B_GA_0004"},
            positives_null=set(),
        )

        # ── 1. Sample the null pool for each genome ──
        for prefix, g in (("A", ga), ("B", gb)):
            null_fasta = tmp_path / f"{prefix}_null.faa"
            null_ids = tmp_path / f"{prefix}_null_ids.tsv"
            r = _run_script(
                "sample_null_proteins.py",
                [
                    "--proteins",
                    str(g["proteins"]),
                    "--gene-order",
                    str(g["gene_order"]),
                    "--ss-components",
                    str(g["ss_components"]),
                    "--window",
                    "3",
                    "--n",
                    "5",
                    "--seed",
                    "42",
                    "--out-fasta",
                    str(null_fasta),
                    "--out-ids",
                    str(null_ids),
                ],
            )
            assert r.returncode == 0, r.stderr
            g["null_ids"] = null_ids

        # ── 2. Run per-genome enrichment_testing ──
        per_genome_tsvs = []
        for prefix, g in (("A", ga), ("B", gb)):
            out = tmp_path / f"{prefix}_enrichment_stats.tsv"
            r = _run_script(
                "enrichment_testing.py",
                [
                    "--ss-components",
                    str(g["ss_components"]),
                    "--gene-order",
                    str(g["gene_order"]),
                    "--dlp",
                    str(g["dlp"]),
                    "--dse",
                    str(g["dse"]),
                    "--null-ids",
                    str(g["null_ids"]),
                    "--window",
                    "3",
                    "--conf-threshold",
                    "0.8",
                    "--sample",
                    f"genome_{prefix}",
                    "--out",
                    str(out),
                ],
            )
            assert r.returncode == 0, r.stderr
            per_genome_tsvs.append(out)

        # Each per-genome file: one T2SS system × 2 tools = 2 rows
        for tsv in per_genome_tsvs:
            rows = list(csv.DictReader(open(tsv), delimiter="\t"))
            assert len(rows) == 2
            assert {r["tool"] for r in rows} == {"DLP", "DSE"}
            assert {r["scope_kind"] for r in rows} == {"system"}

        # ── 3. Pool across the two genomes ──
        # Run via the public helper (the GUI's call site). The CLI alone
        # doesn't ship a pooled-mode subcommand -- Home.py invokes the
        # helper directly after _merge_genome_outputs.
        sys.path.insert(0, str(PROJECT_ROOT / "src"))
        from ssign_app.core.runner import pool_enrichment_stats

        pooled_path = tmp_path / "pooled_enrichment_stats.tsv"
        n_pooled = pool_enrichment_stats([str(p) for p in per_genome_tsvs], str(pooled_path))
        assert n_pooled == 2  # (T2SS, DLP) and (T2SS, DSE)

        pooled = list(csv.DictReader(open(pooled_path), delimiter="\t"))
        assert {r["scope_id"] for r in pooled} == {"T2SS"}
        # M_pool = 6 (genome A) + 6 (genome B) = 12; k_pool = 2 + 2 = 4
        for r in pooled:
            assert int(r["M"]) == 12
            assert int(r["k"]) == 4
            # Both genomes had 0 null positives → p_bg_pool = 0
            assert float(r["p_bg"]) == 0.0
            # Degenerate p_bg → pvalue = 1.0 (binom_pvalue returns 1 when p<=0)
            assert float(r["pvalue"]) == 1.0
