"""Tests for sample_null_proteins.py.

Exercises the in-process sampler and the CLI driver. The 15-gene two-contig
fixture gives a well-defined complement (contig_A genes 0-1, all of
contig_B, none of the +/-3 window of GENE_0005/GENE_0006) so we can hand-
compute the candidate pool.
"""

import os
import random

import pytest
from _helpers import run_script_main
from extract_neighborhood import get_neighborhood_proteins, load_gene_order, load_ss_components
from sample_null_proteins import main as sample_main
from sample_null_proteins import sample_null

from ssign_app.scripts.ssign_lib.fasta_io import read_fasta, write_fasta

# ---------------------------------------------------------------------------
# sample_null — pure function
# ---------------------------------------------------------------------------


class TestSampleNull:
    def test_deterministic_for_same_seed(self):
        ids = [f"x{i}" for i in range(100)]
        a = sample_null(ids, exclude=set(), n=10, rng=random.Random(42))
        b = sample_null(ids, exclude=set(), n=10, rng=random.Random(42))
        assert a == b

    def test_different_seed_picks_different_ids(self):
        ids = [f"x{i}" for i in range(100)]
        a = sample_null(ids, exclude=set(), n=10, rng=random.Random(1))
        b = sample_null(ids, exclude=set(), n=10, rng=random.Random(2))
        assert a != b

    def test_never_includes_excluded(self):
        ids = [f"x{i}" for i in range(50)]
        exclude = {f"x{i}" for i in range(40)}
        picked = sample_null(ids, exclude=exclude, n=20, rng=random.Random(42))
        assert set(picked).isdisjoint(exclude)
        # Pool of 10 candidates, requested 20 → fallback to all 10
        assert len(picked) == 10

    def test_small_pool_returns_all(self):
        ids = [f"x{i}" for i in range(5)]
        picked = sample_null(ids, exclude=set(), n=100, rng=random.Random(42))
        assert set(picked) == set(ids)

    def test_empty_pool_returns_empty(self):
        picked = sample_null(["a", "b"], exclude={"a", "b"}, n=10, rng=random.Random(42))
        assert picked == []

    def test_zero_n_returns_empty(self):
        picked = sample_null(["a", "b", "c"], exclude=set(), n=0, rng=random.Random(42))
        assert picked == []


# ---------------------------------------------------------------------------
# CLI driver — end-to-end with the 15-gene fixture
# ---------------------------------------------------------------------------


@pytest.fixture
def proteome_fasta(tmp_dir):
    """15-protein FASTA matching the two_contig_genes fixture."""
    path = os.path.join(tmp_dir, "proteins.faa")
    seqs = {f"GENE_{i:04d}": "M" + ("A" * 50) for i in range(10)}
    seqs.update({f"GENEB_{i:04d}": "M" + ("L" * 50) for i in range(5)})
    write_fasta(seqs, path)
    return path


class TestCliDriver:
    def _run(self, monkeypatch, tmp_dir, proteome_fasta, gene_order_tsv, ss_components_tsv, n=5, seed=42):
        out_fasta = os.path.join(tmp_dir, "null.faa")
        out_ids = os.path.join(tmp_dir, "null_ids.tsv")
        argv = [
            "sample_null_proteins.py",
            "--proteins",
            proteome_fasta,
            "--gene-order",
            gene_order_tsv,
            "--ss-components",
            ss_components_tsv,
            "--window",
            "3",
            "--n",
            str(n),
            "--seed",
            str(seed),
            "--out-fasta",
            out_fasta,
            "--out-ids",
            out_ids,
        ]
        run_script_main(monkeypatch, sample_main, argv)
        return out_fasta, out_ids

    def test_picks_only_outside_neighborhood(
        self, monkeypatch, tmp_dir, proteome_fasta, gene_order_tsv, ss_components_tsv
    ):
        out_fasta, out_ids = self._run(monkeypatch, tmp_dir, proteome_fasta, gene_order_tsv, ss_components_tsv)
        # Compute the same neighborhood the script computes
        neighborhood = get_neighborhood_proteins(
            load_gene_order(gene_order_tsv),
            load_ss_components(ss_components_tsv),
            window=3,
        )
        picked_ids = [line.strip() for line in open(out_ids) if line.strip()]
        assert set(picked_ids).isdisjoint(neighborhood)
        # Both outputs are in sync (same ID set, same count)
        assert set(read_fasta(out_fasta).keys()) == set(picked_ids)

    def test_deterministic_across_invocations(
        self, monkeypatch, tmp_dir, proteome_fasta, gene_order_tsv, ss_components_tsv
    ):
        out_a, _ = self._run(monkeypatch, tmp_dir, proteome_fasta, gene_order_tsv, ss_components_tsv, n=5, seed=42)
        ids_a = sorted(read_fasta(out_a).keys())
        # Second run into a new tmp path
        new_dir = os.path.join(tmp_dir, "rerun")
        os.makedirs(new_dir)
        out_b = os.path.join(new_dir, "null.faa")
        out_ids_b = os.path.join(new_dir, "null_ids.tsv")
        argv = [
            "sample_null_proteins.py",
            "--proteins",
            proteome_fasta,
            "--gene-order",
            gene_order_tsv,
            "--ss-components",
            ss_components_tsv,
            "--window",
            "3",
            "--n",
            "5",
            "--seed",
            "42",
            "--out-fasta",
            out_b,
            "--out-ids",
            out_ids_b,
        ]
        run_script_main(monkeypatch, sample_main, argv)
        ids_b = sorted(read_fasta(out_b).keys())
        assert ids_a == ids_b

    def test_small_pool_returns_everything(
        self, monkeypatch, tmp_dir, proteome_fasta, gene_order_tsv, ss_components_tsv
    ):
        # Pool size: 15 total - 8 neighborhood (GENE_0002..GENE_0009) = 7.
        # Request 200 → fallback to all 7.
        _, out_ids = self._run(monkeypatch, tmp_dir, proteome_fasta, gene_order_tsv, ss_components_tsv, n=200)
        picked = {line.strip() for line in open(out_ids) if line.strip()}
        assert len(picked) == 7
        # The 7 are: GENE_0000, GENE_0001, GENEB_0000..GENEB_0004
        expected = {"GENE_0000", "GENE_0001"} | {f"GENEB_{i:04d}" for i in range(5)}
        assert picked == expected

    def test_empty_proteome_yields_empty_outputs(self, monkeypatch, tmp_dir, gene_order_tsv, ss_components_tsv):
        empty = os.path.join(tmp_dir, "empty.faa")
        open(empty, "w").close()
        out_fasta = os.path.join(tmp_dir, "null.faa")
        out_ids = os.path.join(tmp_dir, "null_ids.tsv")
        argv = [
            "sample_null_proteins.py",
            "--proteins",
            empty,
            "--gene-order",
            gene_order_tsv,
            "--ss-components",
            ss_components_tsv,
            "--window",
            "3",
            "--n",
            "10",
            "--seed",
            "42",
            "--out-fasta",
            out_fasta,
            "--out-ids",
            out_ids,
        ]
        run_script_main(monkeypatch, sample_main, argv)
        assert os.path.exists(out_fasta) and os.path.getsize(out_fasta) == 0
        assert os.path.exists(out_ids) and os.path.getsize(out_ids) == 0
