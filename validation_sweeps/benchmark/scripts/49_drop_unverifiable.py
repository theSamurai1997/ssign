#!/usr/bin/env python3
"""Drop the citation_trust=unverifiable rows from positives_all (Teo 2026-06-15: verified-only set).

After the full-table audit (47) positives_all = 458 = 330 verified_paper + 3 verified_external +
121 unverifiable + 4 fallback. The unverifiable rows are paywalled papers whose abstract did not name
the protein, no counter-evidence either way. Policy: keep only citation-verified rows. This removes
the 121, leaving 337. Reversible from the backup.

Backup : data/dataset/positives_all.pre_drop_unverifiable.tsv  (written once)
Log    : data/dataset/unverifiable_dropped.tsv
Run    : python3 scripts/49_drop_unverifiable.py
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "data" / "dataset"
sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

POS = DATASET / "positives_all.tsv"
BACKUP = DATASET / "positives_all.pre_drop_unverifiable.tsv"


def main() -> int:
    rows = read_tsv(POS)
    kept = [r for r in rows if r.get("citation_trust") != "unverifiable"]
    dropped = [r for r in rows if r.get("citation_trust") == "unverifiable"]
    if not dropped:
        print("no unverifiable rows present; nothing to drop")
        return 0
    if not BACKUP.exists():
        BACKUP.write_text(POS.read_text())
    write_tsv(POS, list(kept[0].keys()), kept)
    write_tsv(DATASET / "unverifiable_dropped.tsv", list(dropped[0].keys()), dropped)
    print(f"positives_all: {len(rows)} -> {len(kept)} kept, {len(dropped)} unverifiable dropped")
    print("dropped by ss_type:", dict(Counter(r["ss_type"] for r in dropped)))
    print("kept by ss_type:   ", dict(Counter(r["ss_type"] for r in kept)))
    print("kept by citation_trust:", dict(Counter(r["citation_trust"] for r in kept)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
