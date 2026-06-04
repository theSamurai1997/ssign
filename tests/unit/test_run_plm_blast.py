"""Unit tests for run_plm_blast.py — CSV parser.

Does not exercise the `plmblast.py` subprocess itself (that requires the
pLM-BLAST install plus a ~10 GB ECOD30 database). Tests target the
pure-Python helper that parses pLM-BLAST's CSV output.
"""

import os

from run_plm_blast import (
    _OUTPUT_FIELDNAMES,
    _reduce_to_top1,
    _to_output_row,
    _use_cuda_for_embedding,
    load_substrate_ids,
    parse_plmblast_csv,
    write_substrates_only_fasta,
)

# Representative pLM-BLAST CSV fixture. Real output columns verified on
# first integration run — swap in the actual captured output as a fixture
# once we run pLM-BLAST against the T1SS fixture on CX3.
_PLMBLAST_FIXTURE = (
    "qid,sid,sdesc,score,qstart,qend,tstart,tend\n"
    "GENE_00001,ecod_1a2bA1,Autotransporter beta-domain,0.912,10,145,5,140\n"
    "GENE_00001,ecod_3c4dB2,Hemolysin secretion ATP-binding,0.871,20,150,10,140\n"
    "GENE_00003,ecod_5e6fC1,TolC-like outer membrane channel,0.789,1,200,1,199\n"
)

# Pre-fix fixture (no sdesc column) — upstream emits this when the DB has
# no description column. Guards the empty-description fallback so we
# don't crash on DBs without metadata.
_PLMBLAST_FIXTURE_NO_SDESC = "qid,sid,score,qstart,qend,tstart,tend\nGENE_00001,ecod_1a2bA1,0.912,10,145,5,140\n"


class TestParsePlmblastCsv:
    def test_extracts_all_hits(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        assert len(entries) == 3

    def test_multiple_hits_per_protein_preserved(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        hits_for_gene1 = [e for e in entries if e["protein_id"] == "GENE_00001"]
        assert len(hits_for_gene1) == 2

    def test_captures_all_expected_fields(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        first = next(e for e in entries if e["protein_id"] == "GENE_00001")
        assert first["target_id"] == "ecod_1a2bA1"
        assert first["score"] == "0.912"
        assert first["qstart"] == "10"
        assert first["qend"] == "145"
        assert first["tstart"] == "5"
        assert first["tend"] == "140"

    def test_empty_qid_row_skipped(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(
                "qid,sid,score,qstart,qend,tstart,tend\n,ecod_empty,0.5,1,10,1,10\nGENE_X,ecod_real,0.9,1,100,1,100\n"
            )

        entries = parse_plmblast_csv(path)
        assert len(entries) == 1
        assert entries[0]["protein_id"] == "GENE_X"

    def test_empty_csv_returns_empty_list(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write("qid,sid,score,qstart,qend,tstart,tend\n")

        assert parse_plmblast_csv(path) == []

    def test_output_shape_is_stable(self, tmp_dir):
        """Every entry has the six expected fields for downstream consumers."""
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)

        entries = parse_plmblast_csv(path)
        required = {
            "protein_id",
            "target_id",
            "score",
            "qstart",
            "qend",
            "tstart",
            "tend",
        }
        for e in entries:
            assert required <= set(e.keys())


class TestParsePlmblastCsvSdesc:
    """The wrapper used to drop the upstream `sdesc` description column,
    leaving integrate_annotations.py with no `ecod_top1_description`
    value (task #80). These tests pin the new contract: pass through
    `sdesc` when present, fall back to empty when not."""

    def test_sdesc_captured_when_present(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)
        entries = parse_plmblast_csv(path)
        first = next(e for e in entries if e["protein_id"] == "GENE_00001")
        assert first["description"] == "Autotransporter beta-domain"

    def test_description_empty_when_sdesc_column_missing(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE_NO_SDESC)
        entries = parse_plmblast_csv(path)
        assert entries[0]["description"] == ""


class TestReduceToTop1:
    """Multi-hit-per-query reduction. Pre-fix, the wrapper emitted every
    hit, causing integrate_annotations to either explode the substrate
    row count or (more often) drop the file entirely. Top-1 keeps the
    highest-scoring hit per query."""

    def test_one_row_per_query(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)
        top1 = _reduce_to_top1(parse_plmblast_csv(path))
        ids = sorted(e["protein_id"] for e in top1)
        assert ids == ["GENE_00001", "GENE_00003"]

    def test_picks_highest_score(self, tmp_dir):
        path = os.path.join(tmp_dir, "plm_blast.csv")
        with open(path, "w") as f:
            f.write(_PLMBLAST_FIXTURE)
        top1 = _reduce_to_top1(parse_plmblast_csv(path))
        gene1 = next(e for e in top1 if e["protein_id"] == "GENE_00001")
        assert gene1["target_id"] == "ecod_1a2bA1"
        assert gene1["score"] == "0.912"

    def test_non_numeric_score_does_not_crash(self):
        entries = [
            {
                "protein_id": "G1",
                "target_id": "T1",
                "description": "",
                "score": "not-a-number",
                "qstart": "",
                "qend": "",
                "tstart": "",
                "tend": "",
            },
            {
                "protein_id": "G1",
                "target_id": "T2",
                "description": "",
                "score": "0.5",
                "qstart": "",
                "qend": "",
                "tstart": "",
                "tend": "",
            },
        ]
        top1 = _reduce_to_top1(entries)
        assert len(top1) == 1
        assert top1[0]["target_id"] == "T2"

    def test_empty_list_returns_empty(self):
        assert _reduce_to_top1([]) == []


class TestToOutputRow:
    """The wrapper's final-output column schema is the contract that
    integrate_annotations.TOOL_HIT_COLUMNS["pLM-BLAST"] relies on. This
    test pins the column names so a future rename here triggers a test
    failure before the silent-merge-skip bug returns."""

    def test_columns_match_consumer_contract(self):
        row = _to_output_row(
            {
                "protein_id": "GENE_X",
                "target_id": "ecod_99",
                "description": "Some domain",
                "score": "0.91",
                "qstart": "1",
                "qend": "100",
                "tstart": "1",
                "tend": "100",
            }
        )
        assert set(row.keys()) == set(_OUTPUT_FIELDNAMES)
        assert row["locus_tag"] == "GENE_X"
        assert row["ecod_top1_id"] == "ecod_99"
        assert row["ecod_top1_description"] == "Some domain"
        assert row["ecod_top1_score"] == "0.91"

    def test_output_fieldnames_match_integrate_annotations_contract(self):
        """The integrate_annotations consumer reads `ecod_top1_description`
        and joins on `locus_tag`. Pin both."""
        assert "locus_tag" in _OUTPUT_FIELDNAMES
        assert "ecod_top1_description" in _OUTPUT_FIELDNAMES


class TestSubstrateFiltering:
    """Runs before pLM-BLAST itself; verifies we only search substrates, not whole genome."""

    def test_load_substrate_ids_basic(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\tsample\n")
            f.write("G1\tsample1\n")
            f.write("G2\tsample1\n")
        assert load_substrate_ids(path) == {"G1", "G2"}

    def test_load_substrate_ids_skips_empty_rows(self, tmp_dir):
        path = os.path.join(tmp_dir, "substrates.tsv")
        with open(path, "w") as f:
            f.write("locus_tag\n")
            f.write("G1\n")
            f.write("\n")
            f.write("G2\n")
        assert load_substrate_ids(path) == {"G1", "G2"}

    def test_write_substrates_only_fasta_filters(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        with open(src, "w") as f:
            f.write(">G1\nMKT\n>G2\nMFV\n>G3\nMQK\n")
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G1", "G3"}, out)
        assert n == 2
        with open(out) as f:
            body = f.read()
        assert ">G1" in body and ">G3" in body
        assert ">G2" not in body

    def test_write_substrates_only_fasta_missing_ids_skipped(self, tmp_dir):
        src = os.path.join(tmp_dir, "proteins.faa")
        with open(src, "w") as f:
            f.write(">G1\nMKT\n")
        out = os.path.join(tmp_dir, "substrates.faa")
        n = write_substrates_only_fasta(src, {"G1", "NOT_IN_FASTA"}, out)
        assert n == 1


class TestUseCudaForEmbedding:
    """Regression coverage for the GPU auto-detect that controls --cuda on embeddings.py.

    Pre-fix, the wrapper never passed --cuda regardless of GPU availability,
    causing ProtT5 embedding to run on CPU on full-tier GPU nodes (~100x
    slower). These tests pin the detection contract so a future refactor
    doesn't silently regress it.
    """

    def test_force_cpu_env_var_overrides_detection(self, monkeypatch):
        monkeypatch.setenv("SSIGN_PLMBLAST_FORCE_CPU", "1")
        assert _use_cuda_for_embedding() is False

    def test_force_cpu_accepts_alternate_truthy_values(self, monkeypatch):
        for v in ("true", "TRUE", "yes", "Yes"):
            monkeypatch.setenv("SSIGN_PLMBLAST_FORCE_CPU", v)
            assert _use_cuda_for_embedding() is False, f"failed for value {v!r}"

    def test_force_cpu_falsy_value_does_not_override(self, monkeypatch):
        # An empty value (or "0", "false") should NOT force CPU — auto-detect wins.
        monkeypatch.setenv("SSIGN_PLMBLAST_FORCE_CPU", "0")
        # Result depends on whether torch + CUDA is available; we just verify
        # the env-var didn't short-circuit to False the way "1" does.
        try:
            import torch

            expected = bool(torch.cuda.is_available())
        except ImportError:
            expected = False
        assert _use_cuda_for_embedding() is expected

    def test_missing_torch_falls_back_to_cpu(self, monkeypatch):
        import builtins

        real_import = builtins.__import__

        def fake_import(name, *args, **kwargs):
            if name == "torch":
                raise ImportError("torch not installed (test stub)")
            return real_import(name, *args, **kwargs)

        monkeypatch.delenv("SSIGN_PLMBLAST_FORCE_CPU", raising=False)
        monkeypatch.setattr(builtins, "__import__", fake_import)
        assert _use_cuda_for_embedding() is False


class TestResolveEmbeddingsScriptErrors:
    """When `_resolve_plmblast_script` returns the bare-name fallback
    "plmblast.py", `_resolve_embeddings_script` used to compute
    `dirname(dirname(abspath("plmblast.py")))` against cwd and emit a
    "embeddings.py not found at <cwd>/.." message that looked like a
    path-config bug. Real cause: pLM-BLAST not installed. Test pins the
    clearer error message."""

    def test_bare_name_fallback_raises_install_error(self):
        import pytest
        from run_plm_blast import _resolve_embeddings_script

        with pytest.raises(RuntimeError) as exc:
            _resolve_embeddings_script("plmblast.py")
        msg = str(exc.value)
        assert "not installed" in msg
        assert "SSIGN_PLMBLAST_SCRIPT" in msg
        assert "git clone" in msg
        # Old misleading text must NOT appear.
        assert "embeddings.py not found next to" not in msg

    def test_real_path_with_missing_embeddings_raises_incomplete_install(self, tmp_path):
        # Simulate a partial clone: scripts/plmblast.py exists, but
        # embeddings.py one level up does not.
        import pytest
        from run_plm_blast import _resolve_embeddings_script

        scripts_dir = tmp_path / "pLM-BLAST" / "scripts"
        scripts_dir.mkdir(parents=True)
        fake_plmblast = scripts_dir / "plmblast.py"
        fake_plmblast.write_text("# stub")

        with pytest.raises(RuntimeError) as exc:
            _resolve_embeddings_script(str(fake_plmblast))
        assert "embeddings.py not found next to" in str(exc.value)

    def test_real_path_with_embeddings_returns_it(self, tmp_path):
        from run_plm_blast import _resolve_embeddings_script

        repo = tmp_path / "pLM-BLAST"
        (repo / "scripts").mkdir(parents=True)
        fake_plmblast = repo / "scripts" / "plmblast.py"
        fake_plmblast.write_text("# stub")
        fake_embed = repo / "embeddings.py"
        fake_embed.write_text("# stub")

        assert _resolve_embeddings_script(str(fake_plmblast)) == str(fake_embed)


class TestMainEndToEnd:
    """End-to-end test for run_plm_blast.main(): stub the upstream
    plmblast.py + embeddings.py subprocesses (so we don't need the real
    ~5 GB install or ECOD30 DB), feed in a synthetic raw CSV, and assert
    the wrapper's final TSV matches integrate_annotations.py's contract.

    This is the regression test for task #80 — guards against the
    silent-drop bug returning if anyone renames `_OUTPUT_FIELDNAMES`
    or skips the top-1 reduction.
    """

    def _setup_fake_plmblast(self, tmp_path, monkeypatch):
        """Create a fake pLM-BLAST tree (plmblast.py + embeddings.py) and
        a fake ECOD DB dir so the wrapper's preflight checks pass."""
        repo = tmp_path / "fake_plmblast"
        (repo / "scripts").mkdir(parents=True)
        (repo / "scripts" / "plmblast.py").write_text("# stub\n")
        (repo / "embeddings.py").write_text("# stub\n")
        monkeypatch.setenv("SSIGN_PLMBLAST_SCRIPT", str(repo / "scripts" / "plmblast.py"))

        ecod_db = tmp_path / "ecod_db"
        ecod_db.mkdir()
        return repo, ecod_db

    def _stub_subprocess(self, monkeypatch, raw_csv_content: str):
        """Stub subprocess.run so embeddings.py is a no-op (creates the
        expected .pt + .pt.csv files) and plmblast.py writes the raw
        CSV at the output path the wrapper requested."""
        import subprocess

        real_run = subprocess.run

        from pathlib import Path

        def fake_run(cmd, *args, **kwargs):
            # Embeddings call: [python, bootstrap, embed_script, "start",
            #                   embed_input, pt_path, "-embedder", "pt", ...]
            # plmblast call:   [python, plmblast.py, ecod_db,
            #                   query_emb_base, out_csv, "-cpc", ...]
            if "start" in cmd:
                pt_path = cmd[5]
                Path(pt_path).write_bytes(b"")
                Path(pt_path + ".csv").write_text("id\n")
                return subprocess.CompletedProcess(cmd, 0, "", "")
            # plmblast call — out_csv is positional arg 4.
            csv_out = cmd[4]
            Path(csv_out).write_text(raw_csv_content)
            return subprocess.CompletedProcess(cmd, 0, "", "")

        monkeypatch.setattr(subprocess, "run", fake_run)
        return real_run

    def test_output_columns_match_integrate_annotations_contract(self, tmp_path, monkeypatch):
        """The wrapper's TSV must carry locus_tag + ecod_top1_description
        so integrate_annotations can left-join on locus_tag and pick up
        the description for consensus voting."""
        import csv as _csv

        _repo, ecod_db = self._setup_fake_plmblast(tmp_path, monkeypatch)
        self._stub_subprocess(monkeypatch, _PLMBLAST_FIXTURE)

        # Inputs the wrapper expects.
        proteins = tmp_path / "proteins.faa"
        proteins.write_text(">GENE_00001\nMKT\n>GENE_00003\nMFV\n")
        substrates = tmp_path / "subs.tsv"
        substrates.write_text("locus_tag\nGENE_00001\nGENE_00003\n")
        out = tmp_path / "plm_blast.tsv"

        import sys as _sys

        monkeypatch.setattr(
            _sys,
            "argv",
            [
                "run_plm_blast.py",
                "--substrates",
                str(substrates),
                "--proteins",
                str(proteins),
                "--ecod-db",
                str(ecod_db),
                "--out",
                str(out),
            ],
        )
        # Bypass the GPU check (we're in a unit-test env).
        monkeypatch.setenv("SSIGN_PLMBLAST_FORCE_CPU", "1")

        from run_plm_blast import main

        rc = main()
        assert rc == 0
        assert out.exists()

        with open(out) as f:
            rows = list(_csv.DictReader(f, delimiter="\t"))

        # The contract: locus_tag + ecod_top1_description present, and
        # top-1 reduction applied (2 hits for GENE_00001 → 1 row).
        assert len(rows) == 2
        ids = sorted(r["locus_tag"] for r in rows)
        assert ids == ["GENE_00001", "GENE_00003"]
        gene1 = next(r for r in rows if r["locus_tag"] == "GENE_00001")
        assert gene1["ecod_top1_description"] == "Autotransporter beta-domain"
        assert gene1["ecod_top1_score"] == "0.912"
