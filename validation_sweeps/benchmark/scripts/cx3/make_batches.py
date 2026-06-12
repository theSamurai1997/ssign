#!/usr/bin/env python3
"""
make_batches.py  (Phase 2 task 6.1b: split the panel into CX3 submit batches + pilot list)

CX3 gpu72 enforces a hidden ~2-concurrent-job-per-user cap (see memory cx3_gpu72_placement),
so the 67 panel genomes are run as a handful of serial-within-job batches submitted two at a
time, not 67 separate jobs. Writes:
  batches/pilot.txt      : the 4 T3SS-pilot units (run twice: default vs T3SS-included)
  batches/batch_NN.txt   : the full 67-unit panel, ~BATCH_SIZE units each
  batches/SUBMIT.sh      : the exact qsub lines (pilot first, then panel)

Run (locally, before rsyncing to CX3):
  .venv/bin/python scripts/cx3/make_batches.py
"""

from __future__ import annotations

import csv
from pathlib import Path

BENCH = Path(__file__).resolve().parents[2]
MANIFEST = BENCH / "data" / "phase2" / "panel_manifest.tsv"
BATCHDIR = BENCH / "scripts" / "cx3" / "batches"
BATCH_SIZE = 12

# T3SS-exclusion pilot: 3 T3SS-rich genomes spanning chromosomal + plasmid-borne injectisomes,
# plus PAO1 (flagella, NO injectisome T3SS) as the false-positive control that probes the
# DeepSecE-flagellar-misclassification concern behind ssign's stock T3SS exclusion.
PILOT = [
    "NC_003197.2",  # Salmonella Typhimurium LT2 (+pSLT) — dual SPI-1/SPI-2, 46 T3SS effectors
    "NC_004337.2",  # Shigella flexneri 2a (+pCP301) — plasmid-borne Mxi-Spa, tests plasmid input
    "NC_004578.1",  # Pseudomonas syringae pv. tomato DC3000 — Hrp T3SS, 25 effectors
    "NC_002516.2",  # P. aeruginosa PAO1 — flagella but no injectisome: T3SS false-positive control
]


def main():
    with open(MANIFEST) as fh:
        units = [r["unit_id"] for r in csv.DictReader(fh, delimiter="\t")]
    missing = [u for u in PILOT if u not in units]
    if missing:
        raise SystemExit(f"pilot units not in manifest: {missing}")

    BATCHDIR.mkdir(parents=True, exist_ok=True)
    (BATCHDIR / "pilot.txt").write_text("\n".join(PILOT) + "\n")

    batches = [units[i : i + BATCH_SIZE] for i in range(0, len(units), BATCH_SIZE)]
    names = []
    for n, batch in enumerate(batches, 1):
        name = f"batch_{n:02d}.txt"
        (BATCHDIR / name).write_text("\n".join(batch) + "\n")
        names.append(name)

    pbs = "$HOME/bench/run_benchmark_batch.pbs"
    lines = [
        "#!/bin/bash",
        "# Phase 2 submit plan. Defaults (INPUT_DIR, INPUT_DIR_GB, SSIGN_DB, SSIGN_VENV) come",
        "# from run_benchmark_batch.pbs; override with -v if your CX3 layout differs.",
        "# gpu72 enforces ~2 concurrent jobs/user -> submit in pairs, not all at once.",
        "set -eu",
        'B="$HOME/bench/batches"',
        "",
        "# 1) PILOT FIRST (4 genomes) — answers two questions before the full panel:",
        "#    input mode (GenBank vs Bakta) and T3SS exclusion. Submit ~2 at a time.",
        f"qsub -v BATCH_FILE=$B/pilot.txt,INPUT_MODE=genbank,RUN_TAG=pilot_genbank_default {pbs}",
        f"qsub -v BATCH_FILE=$B/pilot.txt,INPUT_MODE=fasta,RUN_TAG=pilot_fasta_default {pbs}",
        f"qsub -v BATCH_FILE=$B/pilot.txt,INPUT_MODE=genbank,INCLUDE_T3SS=1,RUN_TAG=pilot_genbank_t3ss {pbs}",
        "",
        "# 2) FULL PANEL — run only AFTER we review the pilot and lock input-mode + T3SS.",
        "#    Default below is GenBank (exact, no bridge); add a fasta pass if the pilot shows",
        "#    Bakta changes recall. Uncomment when ready.",
    ]
    for name in names:
        lines.append(f"# qsub -v BATCH_FILE=$B/{name},INPUT_MODE=genbank,RUN_TAG=panel_genbank {pbs}")
    (BATCHDIR / "SUBMIT.sh").write_text("\n".join(lines) + "\n")

    print(
        f"wrote {BATCHDIR.relative_to(BENCH)}/  ({len(names)} panel batches of <= {BATCH_SIZE} + pilot.txt + SUBMIT.sh)"
    )
    print(f"  pilot: {', '.join(PILOT)}")
    print(f"  panel: {len(units)} units across {len(names)} batches")


if __name__ == "__main__":
    raise SystemExit(main())
