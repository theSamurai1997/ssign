"""Tests for enrichment_testing.py (A+ rewrite).

Replaces the old Fisher's-exact + dead-permutation tests. The new
module exposes pure functions for: broad-type collapse, DLP/DSE
positivity, BH FDR with monotone non-decreasing q, scope-level scoring,
and a CLI driver that wires everything together against the
neighborhood / null-sample TSVs.
"""

import os

import pytest
from _helpers import run_script_main, write_tsv
from enrichment_testing import (
    bh_fdr,
    binom_pvalue,
    broad_type,
    is_dlp_positive,
    is_dse_positive,
    score_scope,
)
from enrichment_testing import (
    main as enrichment_main,
)

# ---------------------------------------------------------------------------
# broad_type
# ---------------------------------------------------------------------------


class TestBroadType:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("T1SS", "T1SS"),
            ("T2SS", "T2SS"),
            ("T3SS", "T3SS"),
            ("T4SS", "T4SS"),
            ("pT4SSt", "T4SS"),
            ("T5aSS", "T5SS"),
            ("T5bSS", "T5SS"),
            ("T5cSS", "T5SS"),
            ("T6SSi", "T6SS"),
            ("T6SSii", "T6SS"),
            ("T6SSiii", "T6SS"),
            ("Flagellum", "Flagellum"),
            ("Tad", "Tad"),
            ("", ""),
        ],
    )
    def test_collapse(self, raw, expected):
        assert broad_type(raw) == expected


# ---------------------------------------------------------------------------
# Positivity rules
# ---------------------------------------------------------------------------


class TestDlpPositive:
    def test_above_threshold_is_positive(self):
        assert is_dlp_positive({"dlp_extracellular_prob": "0.9"}, conf=0.8)

    def test_at_threshold_is_positive(self):
        # >= threshold per proximity_analysis.py:153
        assert is_dlp_positive({"dlp_extracellular_prob": "0.8"}, conf=0.8)

    def test_below_threshold_is_negative(self):
        assert not is_dlp_positive({"dlp_extracellular_prob": "0.79"}, conf=0.8)

    def test_legacy_column_name(self):
        # cross_validate emits dlp_extracellular_prob, but the predictions
        # TSV may also have the legacy `extracellular_prob` field.
        assert is_dlp_positive({"extracellular_prob": "0.9"}, conf=0.8)

    def test_missing_column_is_negative(self):
        assert not is_dlp_positive({}, conf=0.8)

    def test_malformed_value_is_negative(self):
        assert not is_dlp_positive({"dlp_extracellular_prob": "n/a"}, conf=0.8)


class TestDsePositive:
    def test_secreted_above_threshold(self):
        assert is_dse_positive({"dse_ss_type": "T2SS", "dse_max_prob": "0.9"}, conf=0.8)

    def test_non_secreted_label_is_negative(self):
        assert not is_dse_positive({"dse_ss_type": "Non-secreted", "dse_max_prob": "0.99"}, conf=0.8)

    def test_t3ss_excluded(self):
        # DSE is unreliable on T3SS (per CLAUDE.md and proximity_analysis:162)
        assert not is_dse_positive({"dse_ss_type": "T3SS", "dse_max_prob": "0.99"}, conf=0.8)

    def test_empty_label_is_negative(self):
        assert not is_dse_positive({"dse_ss_type": "", "dse_max_prob": "0.99"}, conf=0.8)

    def test_below_threshold_is_negative(self):
        assert not is_dse_positive({"dse_ss_type": "T2SS", "dse_max_prob": "0.5"}, conf=0.8)


# ---------------------------------------------------------------------------
# Binomial test math
# ---------------------------------------------------------------------------


class TestBinomPvalue:
    def test_obvious_enrichment_gives_low_pvalue(self):
        # k=8 of M=10 with p=0.1 background → very low p-value
        p = binom_pvalue(8, 10, 0.1)
        assert p < 1e-4

    def test_at_expected_rate_gives_high_pvalue(self):
        # k=1 of M=10 with p=0.1 background → unsurprising
        p = binom_pvalue(1, 10, 0.1)
        assert p > 0.5

    def test_zero_n_returns_one(self):
        assert binom_pvalue(0, 0, 0.5) == 1.0

    def test_zero_p_returns_one(self):
        # Degenerate background — can't reject; surfaced as p=1.0
        assert binom_pvalue(5, 10, 0.0) == 1.0

    def test_p_one_returns_one(self):
        assert binom_pvalue(5, 10, 1.0) == 1.0


# ---------------------------------------------------------------------------
# BH FDR
# ---------------------------------------------------------------------------


class TestBhFdr:
    def test_assigns_qvalues_and_significant_flag(self):
        # 4 hypotheses, BH at alpha=0.05.
        # Sorted ascending: p = 0.001, 0.01, 0.04, 0.20
        # Raw q = 0.001*4/1, 0.01*4/2, 0.04*4/3, 0.20*4/4
        #       = 0.004,     0.02,     ~0.053,    0.20
        # Monotone fix from the back: q_4=0.2, q_3=0.053, q_2=0.02, q_1=0.004
        # Significant at alpha=0.05: rows with q < 0.05 → first two.
        rows = [
            {"scope_id": "a", "pvalue": 0.001},
            {"scope_id": "b", "pvalue": 0.20},
            {"scope_id": "c", "pvalue": 0.01},
            {"scope_id": "d", "pvalue": 0.04},
        ]
        bh_fdr(rows)
        by_id = {r["scope_id"]: r for r in rows}
        assert by_id["a"]["significant"] is True
        assert by_id["c"]["significant"] is True
        assert by_id["d"]["significant"] is False
        assert by_id["b"]["significant"] is False
        # Monotone non-decreasing q in ascending-p order
        ordered = sorted(rows, key=lambda r: r["pvalue"])
        qs = [r["qvalue"] for r in ordered]
        assert qs == sorted(qs)

    def test_qvalues_are_capped_at_one(self):
        rows = [{"scope_id": s, "pvalue": 0.9} for s in "abcd"]
        bh_fdr(rows)
        assert all(r["qvalue"] <= 1.0 for r in rows)

    def test_empty_input_no_op(self):
        rows = []
        bh_fdr(rows)
        assert rows == []


# ---------------------------------------------------------------------------
# score_scope
# ---------------------------------------------------------------------------


class TestScoreScope:
    def test_two_rows_one_per_tool(self):
        neigh = {"P1", "P2", "P3", "P4", "P5"}
        dlp = {p: {"dlp_extracellular_prob": "0.95"} for p in neigh}
        dse = {p: {"dse_ss_type": "T2SS", "dse_max_prob": "0.95"} for p in neigh}
        out = score_scope("sys_1", "T2SS", "system", neigh, dlp, dse, p_dlp=0.1, p_dse=0.1, conf=0.8)
        assert {r["tool"] for r in out} == {"DLP", "DSE"}
        assert all(r["M"] == 5 for r in out)
        assert all(r["k"] == 5 for r in out)

    def test_plme_never_emitted(self):
        # PLM-Effector is excluded from the enrichment test entirely: only DLP/DSE rows.
        neigh = {"P1", "P2", "P3", "P4"}
        dlp = {p: {"dlp_extracellular_prob": "0.95"} for p in neigh}
        dse = {p: {"dse_ss_type": "T2SS", "dse_max_prob": "0.95"} for p in neigh}
        out = score_scope("sys_1", "T2SS", "system", neigh, dlp, dse, p_dlp=0.1, p_dse=0.1, conf=0.8)
        assert {r["tool"] for r in out} == {"DLP", "DSE"}
        assert all(r["tool"] != "PLME" for r in out)

    def test_fold_enrich_empty_when_p_bg_zero(self):
        neigh = {"P1"}
        out = score_scope(
            "sys_1",
            "T2SS",
            "system",
            neigh,
            dlp={"P1": {"dlp_extracellular_prob": "0.95"}},
            dse={"P1": {"dse_ss_type": "T2SS", "dse_max_prob": "0.95"}},
            p_dlp=0.0,
            p_dse=0.0,
            conf=0.8,
        )
        assert all(r["fold_enrich"] == "" for r in out)

    def test_empty_neighborhood_returns_no_rows(self):
        out = score_scope("sys_1", "T2SS", "system", set(), {}, {}, p_dlp=0.1, p_dse=0.1, conf=0.8)
        assert out == []


# ---------------------------------------------------------------------------
# CLI driver — end-to-end with synthetic fixtures
# ---------------------------------------------------------------------------


def _write_pred_tsv(path, fieldnames, rows):
    return write_tsv(path, fieldnames, rows, delimiter="\t")


@pytest.fixture
def stats_fixture(tmp_dir, gene_order_tsv):
    """Two T2SS components on contig_A (GENE_0005, GENE_0006), +/-3 window.

    Neighborhood (excluding components): GENE_0002, GENE_0003, GENE_0004,
    GENE_0007, GENE_0008, GENE_0009 — six neighbors. Two of them (GENE_0003,
    GENE_0004) are heavily-positive in both tools so the binomial test
    should be significant against a p_bg of ~2/9 ≈ 0.22.
    """
    ss = os.path.join(tmp_dir, "ss_components.tsv")
    with open(ss, "w") as f:
        f.write("locus_tag\tss_type\tsys_id\texcluded\n")
        f.write("GENE_0005\tT2SS\tcontig_A_T2SS_1\tFalse\n")
        f.write("GENE_0006\tT2SS\tcontig_A_T2SS_1\tFalse\n")

    # Null sample: 9 proteins from contig_A_0..1 + all 5 contig_B + GENE_0009
    # (Actually contig_A_0/1 + GENEB_0..4 = 7. Plus we said 9 — drop the
    # constraint; what matters is the rate p_bg.)
    null_ids = os.path.join(tmp_dir, "null_ids.tsv")
    null_set = ["GENE_0000", "GENE_0001"] + [f"GENEB_{i:04d}" for i in range(5)]
    with open(null_ids, "w") as f:
        for nid in null_set:
            f.write(nid + "\n")

    # DLP / DSE predictions for the full set (neighborhood + null + components).
    # Make GENE_0003, GENE_0004 positive in the neighborhood; 1 of the 7
    # null proteins positive.
    dlp_rows = []
    dse_rows = []
    positives_neigh = {"GENE_0003", "GENE_0004"}
    positives_null = {"GENEB_0000"}  # 1/7 ≈ 0.14 background
    for i in range(10):
        L = f"GENE_{i:04d}"
        if L in positives_neigh:
            dlp_rows.append({"locus_tag": L, "dlp_extracellular_prob": "0.95"})
            dse_rows.append({"locus_tag": L, "dse_ss_type": "T2SS", "dse_max_prob": "0.95"})
        else:
            dlp_rows.append({"locus_tag": L, "dlp_extracellular_prob": "0.10"})
            dse_rows.append({"locus_tag": L, "dse_ss_type": "Non-secreted", "dse_max_prob": "0.10"})
    for i in range(5):
        L = f"GENEB_{i:04d}"
        if L in positives_null:
            dlp_rows.append({"locus_tag": L, "dlp_extracellular_prob": "0.95"})
            dse_rows.append({"locus_tag": L, "dse_ss_type": "T2SS", "dse_max_prob": "0.95"})
        else:
            dlp_rows.append({"locus_tag": L, "dlp_extracellular_prob": "0.10"})
            dse_rows.append({"locus_tag": L, "dse_ss_type": "Non-secreted", "dse_max_prob": "0.10"})

    dlp = _write_pred_tsv(
        os.path.join(tmp_dir, "dlp.tsv"),
        ["locus_tag", "dlp_extracellular_prob"],
        dlp_rows,
    )
    dse = _write_pred_tsv(
        os.path.join(tmp_dir, "dse.tsv"),
        ["locus_tag", "dse_ss_type", "dse_max_prob"],
        dse_rows,
    )
    return {"ss": ss, "null_ids": null_ids, "dlp": dlp, "dse": dse, "gene_order": gene_order_tsv}


class TestCliDriver:
    def _run(self, monkeypatch, tmp_dir, fx):
        out = os.path.join(tmp_dir, "stats.tsv")
        argv = [
            "enrichment_testing.py",
            "--ss-components",
            fx["ss"],
            "--gene-order",
            fx["gene_order"],
            "--dlp",
            fx["dlp"],
            "--dse",
            fx["dse"],
            "--null-ids",
            fx["null_ids"],
            "--window",
            "3",
            "--conf-threshold",
            "0.8",
            "--sample",
            "test",
            "--out",
            out,
        ]
        run_script_main(monkeypatch, enrichment_main, argv)
        return out

    def test_per_system_rows_emitted(self, monkeypatch, tmp_dir, stats_fixture):
        out = self._run(monkeypatch, tmp_dir, stats_fixture)
        import csv

        rows = list(csv.DictReader(open(out), delimiter="\t"))
        # One T2SS system × 2 tools = 2 rows (no broad-type aggregate since
        # there's only one system of this broad type)
        assert len(rows) == 2
        assert {r["tool"] for r in rows} == {"DLP", "DSE"}
        # Neighborhood excludes the two components themselves → M = 6
        assert all(r["M"] == "6" for r in rows)
        # GENE_0003 + GENE_0004 are positive → k = 2
        assert all(r["k"] == "2" for r in rows)
        # p_bg = 1/7 ≈ 0.142857
        for r in rows:
            assert abs(float(r["p_bg"]) - 1 / 7) < 1e-5

    def test_components_excluded_from_neighborhood(self, monkeypatch, tmp_dir, stats_fixture):
        out = self._run(monkeypatch, tmp_dir, stats_fixture)
        import csv

        rows = list(csv.DictReader(open(out), delimiter="\t"))
        # M=6 confirms components GENE_0005, GENE_0006 are excluded (without
        # exclusion M would be 8: indices 2..9 inclusive).
        assert all(r["M"] == "6" for r in rows)

    def test_no_components_writes_header_only(self, monkeypatch, tmp_dir, stats_fixture):
        # Empty ss_components.tsv (just header)
        empty_ss = os.path.join(tmp_dir, "empty_ss.tsv")
        with open(empty_ss, "w") as f:
            f.write("locus_tag\tss_type\tsys_id\texcluded\n")
        argv = [
            "enrichment_testing.py",
            "--ss-components",
            empty_ss,
            "--gene-order",
            stats_fixture["gene_order"],
            "--dlp",
            stats_fixture["dlp"],
            "--dse",
            stats_fixture["dse"],
            "--null-ids",
            stats_fixture["null_ids"],
            "--window",
            "3",
            "--conf-threshold",
            "0.8",
            "--sample",
            "test",
            "--out",
            os.path.join(tmp_dir, "stats.tsv"),
        ]
        run_script_main(monkeypatch, enrichment_main, argv)
        import csv

        rows = list(csv.DictReader(open(os.path.join(tmp_dir, "stats.tsv")), delimiter="\t"))
        assert rows == []
