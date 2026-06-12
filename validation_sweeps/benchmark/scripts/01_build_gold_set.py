#!/usr/bin/env python3
"""Phase 0a, steps 2.1-2.3: build the effector gold set up to the repair queue.

Reads the read-only source corpus (validated-evidence rows only), partitions by
the 2026-06-08 audit's verification_status, and removes biology errors (apparatus
components mislabeled as effectors, immunity/adaptor proteins) so they never enter
the citation-repair queue.

Network repair (UniProt rebuild, DOI re-citation) and the independent verification
pass are later steps (02_*, 03_*); this script is deterministic and offline.

Outputs (to data/gold_build/):
  01_validated_raw.tsv        every validated-evidence row, all SS types
  02_verified.tsv             validated AND VERIFIED (kept as-is, still subject to biology drop)
  02_partial.tsv              validated AND PARTIAL (repair candidates)
  02_dropped_status.tsv       validated AND (FAIL | NEEDS_REVIEW): dropped, with reason
  03_biology_dropped.tsv      apparatus/immunity rows dropped regardless of status, with reason
  03_repair_queue.tsv         PARTIAL minus biology errors -> feeds the network repair step
  03_verified_clean.tsv       VERIFIED minus biology errors -> already gold, no repair needed
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "source_corpus"
OUT = ROOT / "data" / "gold_build"
OUT.mkdir(parents=True, exist_ok=True)

SS_FILES = sorted(SRC.glob("T*_verified.tsv"))

# --- Biology-error rules, taken verbatim from AUDIT_2026-06-08.md ---------
# Apparatus components (Hcp/VgrG/PAAR) are the T6SS firing machinery, not cargo.
# Two "evolved VgrGs" carry C-terminal effector domains and ARE kept.
APPARATUS_GENE = re.compile(r"^(hcp|vgrg|paar)", re.IGNORECASE)
APPARATUS_KEEP_LOCI = {"PA0262", "VC_A0123"}  # VgrG2b, VgrG-3 (evolved VgrGs)

# Immunity / adaptor proteins named in the audit as non-effectors.
# Prefix match (not \b: underscore is a regex word-char, so "Tai_Atu4346" needs
# bare-prefix matching). Every match is logged so the drop list can be eyeballed.
IMMUNITY_GENE = re.compile(r"^(tsi|tdi|tldi|eagr|tai)", re.IGNORECASE)

# Explicit per-locus drops: real effectors but with a wrong instance/type/name
# binding the audit flagged "not rescuable". Keyed on locus_tag (never gene name,
# to avoid the prtA@XCV3671-vs-Dickeya collision). Teo's call 2026-06-10: drop.
EXPLICIT_DROP_LOCI = {
    "BAB1_1671": "audit: row is actually BspE not BspA (wrong name binding)",
    "RPR_RS04075": "audit: RARP-1 secreted by TolC/T1SS, not VirB/T4SS (wrong type)",
    "PA0820": "audit: TseT is an H2 effector at PA3907, not H1 at PA0820 (wrong instance+locus)",
    "BTH_II1883": "audit: TseM mislabeled T6SS-1 (chr I); actually a chr-II instance",
    "BTH_II1884": "audit: RbsB mislabeled T6SS-1 (chr I); actually a chr-II instance",
}


def biology_drop_reason(gene: str, locus: str) -> str | None:
    g = (gene or "").strip()
    locus = (locus or "").strip()
    if locus in EXPLICIT_DROP_LOCI:
        return EXPLICIT_DROP_LOCI[locus]
    if APPARATUS_GENE.match(g) and locus not in APPARATUS_KEEP_LOCI:
        return "apparatus-as-effector (T6SS structural Hcp/VgrG/PAAR)"
    if IMMUNITY_GENE.match(g):
        return "immunity/adaptor protein, not a secreted effector"
    return None


def write_tsv(path: Path, header: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in header})


def main() -> int:
    # 2.1 load validated-evidence rows from every SS table
    validated: list[dict] = []
    header: list[str] = []
    for f in SS_FILES:
        with f.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not header:
                header = list(reader.fieldnames)
            for r in reader:
                if (r.get("evidence_level") or "").strip() == "validated":
                    validated.append(r)
    write_tsv(OUT / "01_validated_raw.tsv", header, validated)

    # 2.2 partition by audit status
    verified = [r for r in validated if (r.get("verification_status") or "").strip() == "VERIFIED"]
    partial = [r for r in validated if (r.get("verification_status") or "").strip() == "PARTIAL"]
    dropped_status = [r for r in validated if (r.get("verification_status") or "").strip() in ("FAIL", "NEEDS_REVIEW")]
    for r in dropped_status:
        r["drop_reason"] = f"audit status {(r.get('verification_status') or '').strip()}"
    write_tsv(OUT / "02_verified.tsv", header, verified)
    write_tsv(OUT / "02_partial.tsv", header, partial)
    write_tsv(OUT / "02_dropped_status.tsv", header + ["drop_reason"], dropped_status)

    # 2.3 biology drop, applied to BOTH verified and partial (a structural protein
    #     is not an effector no matter how clean its citation is)
    biology_dropped: list[dict] = []
    verified_clean: list[dict] = []
    repair_queue: list[dict] = []
    for bucket, dest in ((verified, verified_clean), (partial, repair_queue)):
        for r in bucket:
            reason = biology_drop_reason(r.get("gene", ""), r.get("locus_tag", ""))
            if reason:
                r["drop_reason"] = reason
                biology_dropped.append(r)
            else:
                dest.append(r)
    write_tsv(OUT / "03_biology_dropped.tsv", header + ["drop_reason"], biology_dropped)
    write_tsv(OUT / "03_verified_clean.tsv", header, verified_clean)
    write_tsv(OUT / "03_repair_queue.tsv", header, repair_queue)

    # report
    def by_type(rows):
        c = {}
        for r in rows:
            c[(r.get("ss_type") or "?").strip()] = c.get((r.get("ss_type") or "?").strip(), 0) + 1
        return dict(sorted(c.items()))

    print(f"validated-evidence rows:        {len(validated)}")
    print(f"  VERIFIED:                     {len(verified)}   {by_type(verified)}")
    print(f"  PARTIAL:                      {len(partial)}   {by_type(partial)}")
    print(f"  dropped (FAIL/NEEDS_REVIEW):  {len(dropped_status)}   {by_type(dropped_status)}")
    print(f"biology-error drops:            {len(biology_dropped)}   {by_type(biology_dropped)}")
    print(
        f"  -> from VERIFIED:             {sum(1 for r in biology_dropped if (r.get('verification_status') or '').strip() == 'VERIFIED')}"
    )
    print(
        f"  -> from PARTIAL:              {sum(1 for r in biology_dropped if (r.get('verification_status') or '').strip() == 'PARTIAL')}"
    )
    print(f"VERIFIED clean (already gold):  {len(verified_clean)}   {by_type(verified_clean)}")
    print(f"repair queue (PARTIAL clean):   {len(repair_queue)}   {by_type(repair_queue)}")
    print("\nbiology-dropped genes:")
    for r in biology_dropped:
        print(
            f"  [{(r.get('ss_type') or '').strip()}] {(r.get('gene') or '').strip():18} {(r.get('locus_tag') or '').strip():14} {r['drop_reason']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
