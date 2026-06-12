#!/usr/bin/env python3
"""Dataset group 1, task 1.1: extract the predicted-evidence corpus rows.

The secretion-classifier training set reuses the benchmark's audited corpus as its
positive backbone. The 582 *validated* rows already went through scripts 01-03 to
become the Phase-1 gold set. This script does the equivalent first step for the
*predicted*-evidence rows: load them, apply the SAME 2026-06-08 biology-error rules
and verification_status partitioning as script 01 (imported verbatim so there is one
source of truth), and stage the survivors for the network-verification step (31).

Nothing here is predicted-specific except the evidence_level filter; the partitioning
mirrors 01_build_gold_set so the predicted tier is held to the same bar as the gold set.

Outputs (to data/dataset/):
  30_predicted_raw.tsv             every predicted-evidence row, all SS types
  30_predicted_verified_clean.tsv  predicted AND VERIFIED, minus biology errors (no repair needed)
  30_predicted_repair_queue.tsv    predicted AND PARTIAL, minus biology errors (-> step 31 repair)
  30_predicted_dropped_status.tsv  predicted AND (FAIL | NEEDS_REVIEW): dropped, with reason
  30_predicted_biology_dropped.tsv apparatus/immunity/wrong-binding rows dropped regardless of status
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "source_corpus"
OUT = ROOT / "data" / "dataset"
OUT.mkdir(parents=True, exist_ok=True)

# Reuse script 01's frozen biology-drop rules + writer (single source of truth).
_spec = importlib.util.spec_from_file_location("build_gold_set", Path(__file__).parent / "01_build_gold_set.py")
_gold = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_gold)
biology_drop_reason = _gold.biology_drop_reason
write_tsv = _gold.write_tsv
from bench_io import by_type  # noqa: E402  (scripts/ is on sys.path)


def main() -> int:
    predicted: list[dict] = []
    header: list[str] = []
    for f in sorted(SRC.glob("T*_verified.tsv")):
        with f.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not header:
                header = list(reader.fieldnames)
            for r in reader:
                gene = (r.get("gene") or "").strip()
                if gene.startswith("#"):  # comment/section rows in the TSV
                    continue
                if (r.get("evidence_level") or "").strip() == "predicted":
                    predicted.append(r)

    write_tsv(OUT / "30_predicted_raw.tsv", header, predicted)

    verified, partial, dropped_status, biology_dropped = [], [], [], []
    for r in predicted:
        reason = biology_drop_reason(r.get("gene", ""), r.get("locus_tag", ""))
        if reason:
            biology_dropped.append({**r, "drop_reason": reason})
            continue
        status = (r.get("verification_status") or "").strip()
        if status == "VERIFIED":
            verified.append(r)
        elif status == "PARTIAL":
            partial.append(r)
        else:  # FAIL | NEEDS_REVIEW | (empty)
            dropped_status.append({**r, "drop_reason": f"verification_status={status or 'EMPTY'}"})

    write_tsv(OUT / "30_predicted_verified_clean.tsv", header, verified)
    write_tsv(OUT / "30_predicted_repair_queue.tsv", header, partial)
    write_tsv(OUT / "30_predicted_dropped_status.tsv", header + ["drop_reason"], dropped_status)
    write_tsv(OUT / "30_predicted_biology_dropped.tsv", header + ["drop_reason"], biology_dropped)

    print(f"predicted-evidence rows:        {len(predicted)}   {by_type(predicted)}")
    print(f"  VERIFIED clean (already kept):{len(verified):>4}   {by_type(verified)}")
    print(f"  repair queue (PARTIAL clean): {len(partial):>4}   {by_type(partial)}")
    print(f"  dropped (FAIL/NEEDS_REVIEW):  {len(dropped_status):>4}   {by_type(dropped_status)}")
    print(f"  biology-error drops:          {len(biology_dropped):>4}   {by_type(biology_dropped)}")
    if biology_dropped:
        print("\nbiology-dropped genes:")
        for r in biology_dropped:
            print(
                f"  [{(r.get('ss_type') or '').strip()}] {(r.get('gene') or '').strip():18} "
                f"{(r.get('locus_tag') or '').strip():14} {r['drop_reason']}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
