#!/usr/bin/env python3
"""Dataset group 2, tasks 2.1/2.2/2.4 (+ stages 2.3): instance-assign predicted effectors.

The model is instance-level: a positive is (protein, specific system instance), and the
instance supplies pair-features (gene-distance to that system's machinery). The benchmark
enumerated 90 literature-curated instances across 61 genomes (instances.tsv); the model
uses those, NOT a fresh MacSyFinder pass (same anti-circularity rule as the benchmark).

Per predicted survivor, count same-type benchmark instances in its genome (drift-tolerant
genome match via bench_index.accession_base):
  - exactly one  -> auto-assign (instance_source=auto). No nearest-machinery guessing.
  - two or more  -> ambiguous: stage for the 2.3 literature-audit agent (the corpus
                    sys_instance_id system-name label + primary_ref DOI are the leads).
  - genome not enumerated -> instance-unknown type-level positive (type_level=true,
                    pair-features will be null). Kept, never dropped.

Inputs:
  data/dataset/predicted_audited.tsv   (step 31)
  data/machinery/instances.tsv         (script 08)
Outputs (data/dataset/):
  predicted_instanced.tsv              all 325 + instance_id, instance_source, type_level
  predicted_instances_ambiguous.tsv    multi-instance rows staged for the 2.3 agent
"""

from __future__ import annotations

import importlib.util
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "dataset"
INSTANCES = ROOT / "data" / "machinery" / "instances.tsv"

_spec = importlib.util.spec_from_file_location("bi", Path(__file__).parent / "bench_index.py")
bi = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bi)
from bench_io import read_tsv, write_tsv  # noqa: E402  (scripts/ is on sys.path)


def main() -> int:
    pred = read_tsv(OUT / "predicted_audited.tsv")
    inst = read_tsv(INSTANCES)

    by_genome_type: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in inst:
        by_genome_type[(bi.accession_base(r["refseq_genome"]), r["ss_type"])].append(r)

    out_cols = list(pred[0].keys()) + ["instance_id", "instance_source", "type_level", "instance_candidates"]
    rows: list[dict] = []
    ambiguous: list[dict] = []

    for r in pred:
        r = dict(r)
        cands = by_genome_type.get((bi.accession_base(r["refseq_genome"]), r["ss_type"]), [])
        ids = [c["instance_id"] for c in cands]
        r["instance_candidates"] = ",".join(ids)
        if len(cands) == 1:
            r["instance_id"], r["instance_source"], r["type_level"] = ids[0], "auto", "no"
        elif len(cands) > 1:
            r["instance_id"], r["instance_source"], r["type_level"] = "", "", "no"  # pending 2.3
            ambiguous.append(r)
        else:
            r["instance_id"], r["instance_source"], r["type_level"] = "", "none", "yes"
        rows.append(r)

    write_tsv(OUT / "predicted_instanced.tsv", out_cols, rows)
    # Stage the ambiguous rows with just the columns the audit agent needs.
    amb_cols = [
        "gene",
        "uniprot",
        "locus_tag",
        "organism",
        "refseq_genome",
        "ss_type",
        "sys_instance_id",
        "primary_ref",
        "instance_candidates",
    ]
    write_tsv(OUT / "predicted_instances_ambiguous.tsv", amb_cols, ambiguous)

    src = Counter(r["instance_source"] for r in rows)
    tl = sum(1 for r in rows if r["type_level"] == "yes")
    print(f"predicted survivors: {len(rows)}")
    print(f"  auto-assigned (single same-type instance): {src.get('auto', 0)}")
    print(f"  ambiguous (>=2 instances, -> 2.3 audit):   {len(ambiguous)}")
    print(f"  type-level (genome not enumerated):        {tl}")
    if ambiguous:
        bt = Counter(r["ss_type"] for r in ambiguous)
        print(f"  ambiguous by type: {dict(sorted(bt.items()))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
