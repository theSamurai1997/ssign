#!/usr/bin/env python3
"""Dataset group 2, task 2.4: fold the literature-audit resolutions, emit positives_instanced.

Script 33 auto-assigned single-instance rows and left multi-instance rows pending. The 2.3
literature-audit agent resolved each ambiguous effector to a specific instance (with a
verbatim quote) or left it UNRESOLVED. This script applies those resolutions and emits the
final per-effector instanced table.

Guards: a RESOLVED instance_id is accepted only if it is one of that row's listed
instance_candidates (no out-of-set assignment). UNRESOLVED rows become instance-unknown
type-level positives (type_level=true) -- kept, never dropped.

Inputs:
  data/dataset/predicted_instanced.tsv          (step 33)
  data/dataset/instance_audit_results.json      (2.3 agent output)
Output:
  data/dataset/positives_instanced.tsv          (325 predicted, instance-resolved where possible)
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from bench_io import read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "dataset"


def main() -> int:
    rows = read_tsv(OUT / "predicted_instanced.tsv")
    audit_path = OUT / "instance_audit_results.json"
    if not audit_path.exists():
        raise SystemExit(
            f"missing {audit_path.name}: run the 2.3 literature-audit agent on "
            "predicted_instances_ambiguous.tsv first (it writes this file)."
        )
    audit = json.loads(audit_path.read_text())
    by_locus = {(a.get("locus_tag") or "").strip(): a for a in audit}

    applied, unresolved, out_of_set = 0, 0, 0
    for r in rows:
        if r.get("instance_source"):  # auto or none already set; only ambiguous pending are blank
            continue
        a = by_locus.get((r.get("locus_tag") or "").strip())
        cands = set((r.get("instance_candidates") or "").split(","))
        if a and a.get("resolution") == "RESOLVED" and a.get("instance_id", "") in cands:
            r["instance_id"] = a["instance_id"]
            r["instance_source"] = "literature"
            r["type_level"] = "no"
            r["instance_quote"] = a.get("quote", "")
            r["instance_source_doi"] = a.get("source_doi", "")
            applied += 1
        else:
            if a and a.get("resolution") == "RESOLVED":  # id not in candidate set -> reject
                out_of_set += 1
            r["instance_id"] = ""
            r["instance_source"] = "none"
            r["type_level"] = "yes"
            unresolved += 1

    header = list(rows[0].keys())
    for extra in ("instance_quote", "instance_source_doi"):
        if extra not in header:
            header.append(extra)
    write_tsv(OUT / "positives_instanced.tsv", header, rows)

    src = Counter(r["instance_source"] for r in rows)
    print(f"predicted rows: {len(rows)}")
    print(f"  auto:       {src.get('auto', 0)}")
    print(f"  literature: {src.get('literature', 0)}  (applied this run: {applied})")
    print(f"  type-level (none): {src.get('none', 0)}  (incl. {unresolved} ambiguous left unresolved)")
    if out_of_set:
        print(f"  WARNING: {out_of_set} agent assignments rejected (instance_id not in candidate set)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
