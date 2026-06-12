#!/usr/bin/env python3
"""Phase 2 task 8.2: deterministic lower bound on emission precision via independent effector DBs.

A precision floor we can prove: of everything ssign emits, how many are homologous to an
experimentally-verified effector in a database we did NOT use as the recall gold set. Every such hit
is a true positive beyond the 51 gold effectors. This UNDERSTATES precision (it can only confirm
effectors already in a DB, never a genuinely novel one), so it is a floor, not the answer.

Reference set: SecReT4 (540 verified T4SS effectors) + SecReT6 (331 experimental T6SS effectors),
matched by protein homology (pyhmmer phmmer; no external aligner needed). An emission is:
  confirmed  - >=90% identity over >=80% of its length to a DB effector (E<1e-5): essentially the same
               protein -> a true positive.
  homologous - a significant hit (E<1e-3) but below the identity/coverage bar: suggestive only.
  none       - no significant DB hit.

Lower-bound precision = (gold + DB-confirmed) / emissions, reported for the proximity-called subset
(the real question) and split by SS type. T5SS-self emissions are reported separately.

Inputs : data/phase2/emissions.<tag>.tsv
         data/external_dbs/secret4/verified_effectors.fas, secret6/effector_exp_protein.fasta
Outputs: data/phase2/emissions_dbmatch.<tag>.tsv   (per emission: best DB hit, identity, class)
Run:     <repo>/.venv/bin/python scripts/29_db_confirmed.py --run-tag panel_genbank_default
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, SequenceFile, TextSequence

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
DBS = [
    ("SecReT4", BENCH / "data" / "external_dbs" / "secret4" / "verified_effectors.fas"),
    ("SecReT6", BENCH / "data" / "external_dbs" / "secret6" / "effector_exp_protein.fasta"),
]
ID_MIN, COV_MIN, E_CONF, E_HOM = 0.90, 0.80, 1e-5, 1e-3
ABC = Alphabet.amino()


def load_db():
    seqs = []
    for src, path in DBS:
        if not path.exists():
            sys.exit(f"missing DB: {path.relative_to(BENCH)}")
        with SequenceFile(str(path), digital=True, alphabet=ABC) as sf:
            for s in sf:
                s.name = f"{src}|{s.name.decode() if isinstance(s.name, bytes) else s.name}".encode()
                seqs.append(s)
    return DigitalSequenceBlock(ABC, seqs)


def pid_cov(aln, qlen):
    """(percent identity over aligned columns, fraction of the query covered)."""
    hs, ts = aln.hmm_sequence, aln.target_sequence
    cols = [(a, b) for a, b in zip(hs, ts) if a not in "-." and b not in "-."]
    if not cols:
        return 0.0, 0.0
    ident = sum(1 for a, b in cols if a.upper() == b.upper())
    return ident / len(cols), len(cols) / qlen if qlen else 0.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-tag", required=True)
    args = ap.parse_args()
    emissions = read_tsv(BENCH / "data" / "phase2" / f"emissions.{args.run_tag}.tsv")

    queries = []
    for i, e in enumerate(emissions):
        seq = (e.get("sequence") or "").strip().upper()
        if not seq:
            continue
        try:
            queries.append(TextSequence(name=str(i).encode(), sequence=seq).digitize(ABC))
        except Exception:
            pass  # non-standard residues that easel rejects -> treated as no-hit, logged below

    targets = load_db()
    best = {}  # query index -> (db_name, pid, cov, evalue)
    for hits in pyhmmer.hmmer.phmmer(queries, targets, cpus=0, E=E_HOM):
        qn = hits.query.name
        qi = int(qn.decode() if isinstance(qn, bytes) else qn)
        qlen = len((emissions[qi].get("sequence") or ""))
        top = next(iter(hits), None)
        if top is None:
            continue
        name = top.name.decode() if isinstance(top.name, bytes) else top.name
        pid, cov = pid_cov(top.best_domain.alignment, qlen)
        best[qi] = (name.split("|", 1)[0], round(pid, 3), round(cov, 3), top.evalue)

    rows = []
    for i, e in enumerate(emissions):
        db, pid, cov, ev = best.get(i, ("", 0.0, 0.0, ""))
        if db and pid >= ID_MIN and cov >= COV_MIN and ev != "" and ev < E_CONF:
            cls = "confirmed"
        elif db:
            cls = "homologous"
        else:
            cls = "none"
        rows.append(
            {
                **{
                    k: e[k]
                    for k in ("unit_id", "locus_tag", "substrate_source", "nearby_ss_types", "is_gold", "gold_gene")
                },
                "db_hit": db,
                "db_identity": pid,
                "db_coverage": cov,
                "db_evalue": ev,
                "db_class": cls,
            }
        )

    out = BENCH / "data" / "phase2" / f"emissions_dbmatch.{args.run_tag}.tsv"
    write_tsv(out, list(rows[0].keys()), rows)

    def is_tp(r):
        return r["is_gold"] == "yes" or r["db_class"] == "confirmed"

    print(f"wrote {out.relative_to(BENCH)}  ({len(rows)} emissions; {len(queries)} digitized)")
    for src in ("proximity", "T5SS-self"):
        sub = [r for r in rows if r["substrate_source"] == src]
        if not sub:
            continue
        tp = sum(is_tp(r) for r in sub)
        conf = sum(r["db_class"] == "confirmed" for r in sub)
        hom = sum(r["db_class"] == "homologous" for r in sub)
        gold = sum(r["is_gold"] == "yes" for r in sub)
        print(f"\n[{src}]  n={len(sub)}  gold={gold}  DB-confirmed={conf}  homologous={hom}")
        print(f"  precision FLOOR (gold|confirmed) = {tp}/{len(sub)} = {tp / len(sub):.1%}")
        per = defaultdict(lambda: [0, 0])
        for r in sub:
            k = (r["nearby_ss_types"] or "?").split(";")[0].split(",")[0]
            per[k][0] += is_tp(r)
            per[k][1] += 1
        for k in sorted(per, key=lambda k: -per[k][1]):
            t, n = per[k]
            print(f"    {k:8s} floor {t:3d}/{n:<4d} = {t / n:5.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
