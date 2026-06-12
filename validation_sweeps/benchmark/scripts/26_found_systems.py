#!/usr/bin/env python3
"""
26_found_systems.py  (Phase 2: how many ssign-DETECTED systems emitted the found effectors)

ssign's <sid>_results.csv carries a "# Secretion Systems (with secreted proteins)" chunk that
lists every MacSyFinder-detected system under a unique sys_id (e.g. NC_002516_proteins_T6SS_4),
its ss_type/wholeness, and the locus_tag of each component. We never need to re-run MacSyFinder:
the detected-systems table is right here in the output.

For each gold effector ssign EMITTED (ssign_call == emitted_secreted), find which detected
system(s) it sits within +/-3 genes of (ssign's own proximity rule), by gene-order distance from
the effector locus to any component locus of the system. Then count distinct systems, and split
by whether the detected system's type matches the effector's answer-key SS type (concordant) or
not (cross-type adjacency -> the protein was emitted via an unrelated nearby system).

Inputs : data/phase2/actual_per_effector.<run_tag>.tsv  (found effectors + ssign_locus)
         data/phase1/ceiling_per_effector.tsv            (effector -> refseq_genome)
         data/phase2/runs/<run_tag>/<unit>/results/<unit>_results.csv  (systems chunk)
Output : data/phase2/found_systems.<run_tag>.tsv  (one row per found effector -> attributed sys)
Run:     .venv/bin/python scripts/26_found_systems.py --run-tag panel_genbank_default
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import bench_index as bi  # noqa: E402
import bench_runout as bo  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
CEILING = BENCH / "data" / "phase1" / "ceiling_per_effector.tsv"
RUNS = BENCH / "data" / "phase2" / "runs"
WINDOW = 3  # ssign's proximity_window default (constants.py)


def _family(ss_type: str) -> str:
    """Collapse a TXSScan subtype to its coarse SS family so it compares to an effector's ss_type.
    T6SSi/T6SSii -> T6SS; T5aSS/T5bSS/T5cSS -> T5SS; T1SS/T2SS/T3SS/T4SS/T4aP/Tad kept as-is."""
    t = (ss_type or "").strip()
    if t.startswith("T6SS"):
        return "T6SS"
    if t.startswith("T5") and t.endswith("SS"):
        return "T5SS"
    return t


def parse_systems_chunk(results_csv: Path):
    """{sys_id: {"ss_type":..., "wholeness":..., "loci":[component locus_tags]}} from the
    '# Secretion Systems (with secreted proteins)' chunk."""
    with open(results_csv, newline="") as fh:
        lines = fh.read().split("\n")
    start = next((i for i, ln in enumerate(lines) if ln.strip().lower().startswith("# secretion systems (with")), None)
    if start is None:
        return {}
    end = next((i for i in range(start + 1, len(lines)) if lines[i].strip().startswith("#")), len(lines))
    block = [ln for ln in lines[start + 1 : end] if ln.strip()]
    if not block:
        return {}
    systems = {}
    for r in csv.DictReader(io.StringIO("\n".join(block))):
        sid = r.get("sys_id", "")
        if not sid:
            continue
        if r.get("record_type") == "system":
            systems.setdefault(sid, {"ss_type": r.get("ss_type", ""), "wholeness": r.get("wholeness", ""), "loci": []})
        elif r.get("record_type") == "component" and r.get("locus_tag"):
            systems.setdefault(sid, {"ss_type": r.get("ss_type", ""), "wholeness": "", "loci": []})
            systems[sid]["loci"].append(r["locus_tag"])
    return systems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-tag", required=True)
    args = ap.parse_args()
    run_root = RUNS / args.run_tag
    if not run_root.is_dir():
        sys.exit(f"no run dir: {run_root.relative_to(BENCH)}")

    idx = bi.load_from_tsv()
    genome_of = {
        (r["gene"], r["uniprot"], r["effector_locus"]): r.get("refseq_genome", "").strip() for r in bi.load_tsv(CEILING)
    }
    found = [
        r
        for r in bi.load_tsv(BENCH / "data" / "phase2" / f"actual_per_effector.{args.run_tag}.tsv")
        if r["ssign_call"] == bo.CALL_EMITTED
    ]

    sys_cache = {}  # unit -> systems dict

    def systems_for(unit):
        if unit not in sys_cache:
            rc = run_root / unit / "results" / f"{unit}_results.csv"
            sys_cache[unit] = parse_systems_chunk(rc) if rc.exists() else {}
        return sys_cache[unit]

    rows = []
    all_sids = set()
    concordant_sids, crosstype_sids = set(), set()
    for r in found:
        unit, eff_locus = r["unit_id"], r["ssign_locus"]
        genome = genome_of.get((r["gene"], r["uniprot"], r["effector_locus"]), "")
        eff_type = r["ss_type"]
        attributed = []  # (sys_id, sys_type, min_dist)
        for sid, s in systems_for(unit).items():
            dists = [d for lt in s["loci"] if (d := idx.gene_distance(genome, eff_locus, lt)) is not None]
            if dists and min(dists) <= WINDOW:
                attributed.append((sid, s["ss_type"], min(dists)))
        attributed.sort(key=lambda x: x[2])
        # ssign labels detected systems by TXSScan subtype (T6SSi, T5aSS, ...); effector ss_type is
        # the coarse family (T6SS). Compare families, else every T6SSi-via-T6SS match reads as cross.
        for sid, stype, _ in attributed:
            all_sids.add((unit, sid))
            (concordant_sids if _family(stype) == eff_type else crosstype_sids).add((unit, sid))
        atypes = [s[1] for s in attributed]
        if not attributed:
            basis = "no_system_in_window"
        elif any(_family(a) == eff_type for a in atypes):
            basis = "own_type"  # legit: a same-type system sits within +/-3 (right protein, right reason)
        else:
            basis = "cross_type_only"  # accidental: ONLY a different-type system is nearby
        rows.append(
            {
                "gene": r["gene"],
                "uniprot": r["uniprot"],
                "ss_type": eff_type,
                "unit_id": unit,
                "ssign_locus": eff_locus,
                "n_attributed_systems": len(attributed),
                "attributed_sys_ids": ";".join(s[0] for s in attributed),
                "attributed_types": ";".join(atypes),
                "min_dist": attributed[0][2] if attributed else "",
                "nearest_type": attributed[0][1] if attributed else "",
                "emission_basis": basis,
            }
        )

    out = BENCH / "data" / "phase2" / f"found_systems.{args.run_tag}.tsv"
    cols = list(rows[0].keys()) if rows else ["gene"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    own = sum(r["emission_basis"] == "own_type" for r in rows)
    cross = sum(r["emission_basis"] == "cross_type_only" for r in rows)
    no_sys = [r for r in rows if r["emission_basis"] == "no_system_in_window"]
    print(f"wrote {out.relative_to(BENCH)}  ({len(rows)} found effectors)")
    print(f"  distinct ssign-detected systems emitting >=1 found effector: {len(all_sids)}")
    print(f"    type-concordant systems (same SS type as effector) : {len(concordant_sids)}")
    print(f"    cross-type systems (different type nearby)          : {len(crosstype_sids)}")
    print("  WHY each found effector was emitted (emission basis):")
    print(f"    own_type        (a same-type system within +/-{WINDOW}, legit)        : {own}/{len(rows)}")
    print(f"    cross_type_only (ONLY a different-type system nearby, accidental): {cross}/{len(rows)}")
    print(
        f"    no_system_in_window (T5SS-self / non-proximity route)           : {len(no_sys)}"
        + (f"  ({', '.join(r['gene'] for r in no_sys)})" if no_sys else "")
    )


if __name__ == "__main__":
    raise SystemExit(main())
