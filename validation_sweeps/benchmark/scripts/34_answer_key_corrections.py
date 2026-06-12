#!/usr/bin/env python3
"""Phase 1/2 answer-key correction: 3 effectors whose machinery cannot be reliably anchored.

The unreachable-missed deep dive (UNREACHABLE_ANALYSIS.md) found two answer-key anchoring bugs:
  - TRP47, TRP32 (Ehrlichia chaffeensis): genuine T1SS substrates (the cited paper explicitly REJECTS
    T4SS for them), but our machinery resolver anchored them on VirB8 — a TYPE IV secretion gene. The
    real Hly T1SS apparatus is at a separate locus, is not labelled in the genome annotation, and ssign
    detects no secretion system in this genome at all.
  - frpC (Neisseria meningitidis): a real but genome-SCATTERED T1SS; the resolver product-matched a lone
    TolC 1340 genes away. There is no canonical T1SS apparatus adjacent to frpC.

In all three, we cannot reliably place the correct T1SS machinery, so "unreachable @±3" (which implies
a measured distance to real machinery) is the wrong label. The honest classification is the existing
`machinery_unanchored` reason -> non-testable (we can't fairly measure reachability). This does NOT change
ssign's found count (all three are not_emitted either way); it moves them from the unreachable bucket to
non-testable, which sharpens the T1SS recall picture (its remaining unreachables are the 2 genuine
biological exceptions, Serralysin + apxIIA).

Patches `testable=no`, `reason=machinery_unanchored`, and blanks reachable_n3/5/7 for the 3 effectors in
the Phase-1 ceiling table and both Phase-2 actual tables. Backs up each file once; logs every change.

Run: .venv/bin/python scripts/34_answer_key_corrections.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
TARGETS = {("frpC", "Q9JYV5"), ("TRP47", "Q2GHU2"), ("TRP32", "Q2GHT8")}
REASON = "machinery_unanchored"  # existing non-testable category
NOTE = "atypical/scattered T1SS; correct machinery not anchorable (was mis-anchored on VirB8/TolC)"
FILES = [
    BENCH / "data" / "phase1" / "ceiling_per_effector.tsv",
    BENCH / "data" / "phase2" / "actual_per_effector.panel_genbank_default.tsv",
    BENCH / "data" / "phase2" / "actual_per_effector.panel_genbank_t3ss.tsv",
]


def main() -> int:
    log = []
    for path in FILES:
        rows = read_tsv(path)
        header = list(rows[0].keys())
        reason_col = "reason" if "reason" in header else "ceiling_reason"
        for r in rows:
            if (r["gene"], r["uniprot"]) not in TARGETS:
                continue
            before = (r.get("testable"), r.get("reachable_n3"), r.get(reason_col))
            r["testable"] = "no"
            for c in ("reachable_n3", "reachable_n5", "reachable_n7"):
                if c in r:
                    r[c] = ""
            r[reason_col] = REASON
            if "testable_reason" in r:
                r["testable_reason"] = NOTE
            log.append((path.name, r["gene"], before, (r["testable"], r["reachable_n3"], r[reason_col])))
        backup = path.with_suffix(path.suffix + ".pre_anchor_fix")
        if not backup.exists():
            backup.write_text(path.read_text())
        write_tsv(path, header, rows)

    print(f"patched {len(TARGETS)} effectors across {len(FILES)} files ({len(log)} row edits):")
    for fn, gene, b, a in log:
        print(f"  {fn:42s} {gene:6s} testable {b[0]}->{a[0]}  reach_n3 {b[1]!r}->{a[1]!r}")
    print("backups: *.pre_anchor_fix")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
