#!/usr/bin/env python3
"""Split the full-table pass-1 sweep into a drop list + balanced agent deep-verify batches.

Per the scoping decision: rows whose deterministic verdict is FLAG_WRONG_TOPIC, FLAG_GENE_ABSENT, or
DOI_UNRESOLVED are dropped from the positive set (provenance too weak/broken to keep as a training
label). The survivors (CONSISTENT + INDETERMINATE) go to an agent deep-verify pass, grouped by the
sourcing paper so each agent reads a paper once and reports on every effector claimed from it.

Outputs:
  data/dataset/deepverify_dropped.tsv          - the 252 dropped rows + drop reason (audit trail)
  data/dataset/deepverify_input/batch_NN.json  - balanced agent input, one file per batch
Each batch JSON: {"batch": NN, "papers": [{doi, ss_type, effectors:[{gene,uniprot,locus_tag,
organism,sys_instance_id,verdict}]}]}. Bin-packed by effector count so batches are even work.

Run: python3 scripts/45_build_deepverify_input.py [n_batches]
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "data" / "dataset"
sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

KEEP = {"CONSISTENT", "INDETERMINATE"}
DROP = {"FLAG_WRONG_TOPIC", "FLAG_GENE_ABSENT", "DOI_UNRESOLVED", "FETCH_ERROR"}


def main(n_batches: int = 20) -> int:
    rows = read_tsv(DATASET / "citation_consistency_full.tsv")
    survivors = [r for r in rows if r["verdict"] in KEEP]
    dropped = [r for r in rows if r["verdict"] in DROP]

    drop_cols = ["gene", "uniprot", "locus_tag", "sys_instance_id", "ss_type", "organism", "sourcing_doi", "verdict"]
    write_tsv(DATASET / "deepverify_dropped.tsv", drop_cols, [{c: r[c] for c in drop_cols} for r in dropped])

    by_doi: dict[str, list[dict]] = defaultdict(list)
    for r in survivors:
        by_doi[r["sourcing_doi"]].append(r)

    papers = []
    for doi, rs in by_doi.items():
        papers.append(
            {
                "doi": doi,
                "ss_type": rs[0]["ss_type"],
                "effectors": [
                    {
                        "gene": r["gene"],
                        "uniprot": r["uniprot"],
                        "locus_tag": r["locus_tag"],
                        "organism": r["organism"],
                        "sys_instance_id": r["sys_instance_id"],
                        "verdict": r["verdict"],
                    }
                    for r in rs
                ],
            }
        )

    # Bin-pack largest-paper-first into the least-loaded batch -> even effector counts per agent.
    papers.sort(key=lambda p: -len(p["effectors"]))
    batches: list[list[dict]] = [[] for _ in range(n_batches)]
    load = [0] * n_batches
    for p in papers:
        i = load.index(min(load))
        batches[i].append(p)
        load[i] += len(p["effectors"])

    outdir = DATASET / "deepverify_input"
    outdir.mkdir(exist_ok=True)
    for old in outdir.glob("batch_*.json"):
        old.unlink()
    for i, b in enumerate(batches):
        if not b:
            continue
        (outdir / f"batch_{i:02d}.json").write_text(json.dumps({"batch": i, "papers": b}, indent=1))

    print(f"survivors: {len(survivors)} effectors across {len(by_doi)} papers")
    print(f"dropped:   {len(dropped)} -> data/dataset/deepverify_dropped.tsv")
    print(f"batches:   {sum(1 for b in batches if b)} files in data/dataset/deepverify_input/")
    print(f"  effectors/batch: min={min(load)} max={max(load)} (target ~{len(survivors) // n_batches})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 20))
