#!/usr/bin/env python3
"""
19_effector_distance.py  (Phase 1 tasks 5.2 + 5.3: effector -> own-machinery distance + classify)

For each gold-set effector, measure the gene-order distance to the nearest component of ITS
OWN system instance, then classify reachable vs impossible at N = 3, 5, 7. This is the
ceiling: how many verified effectors the +/-N proximity rule COULD reach if ssign detected
the machinery perfectly. No ssign run is involved.

Two ceiling sources (column `ceiling_source`, set by script 17):
  - adjacency_scan (the 16 rescued T1SS effectors): the nearest-transporter distance was
    already measured by the adjacency scan (t1ss_ceiling.tsv); we read it straight off.
    The 3 unplaceable T1SS effectors are untestable (no genome).
  - answer_key (489 curated effectors + auto-assignable net-new): distance from the
    effector's CDS to the nearest anchored machinery locus of its own instance, both placed
    in gene order by the alias-aware index (bench_index).

Instance assignment for answer_key effectors:
  - curated rows: (ss_type, genome, sys_instance_id) -> instance_id (instances.tsv).
  - net-new external rows (no sys_instance_id): assigned ONLY when the genome carries exactly
    one same-type instance (unambiguous own-instance). Net-new effectors in multi-instance
    genomes, or in genomes with no curated same-type instance, are left UNTESTABLE
    (own-instance unknown) -- Teo's Checkpoint-A-consistent decision: no nearest-machinery
    guessing, which would circularly minimise the distance we are trying to measure.

Outcome per effector:
  - testable=yes  + nearest_dist (gene units) + reachable_n{3,5,7}. "Impossible at N" is simply
    testable and not reachable at that N (machinery known but further than N, or on a
    different replicon -> structurally unreachable by a +/-N window).
  - testable=no   + reason (no_genome / own_instance_unknown / no_instance_in_genome /
    machinery_unanchored / effector_locus_not_found).

Inputs : data/phase1/effector_gold_set_phase1.tsv, data/phase1/gene_order_index.tsv,
         data/machinery/instances.tsv, data/machinery/machinery_answer_key.tsv,
         data/t1ss_rescue/t1ss_ceiling.tsv
Output : data/phase1/ceiling_per_effector.tsv
Run:     .venv/bin/python scripts/19_effector_distance.py
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import bench_index as bi  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
P1 = BENCH / "data" / "phase1"
GOLD = P1 / "effector_gold_set_phase1.tsv"
INSTANCES = BENCH / "data" / "machinery" / "instances.tsv"
ANSWER_KEY = BENCH / "data" / "machinery" / "machinery_answer_key.tsv"
T1_CEILING = BENCH / "data" / "t1ss_rescue" / "t1ss_ceiling.tsv"
OUT = P1 / "ceiling_per_effector.tsv"

NS = (3, 5, 7)
FIELDS = [
    "gene",
    "uniprot",
    "ss_type",
    "refseq_genome",
    "instance_id",
    "ceiling_source",
    "effector_locus",
    "effector_ordinal",
    "effector_match",
    "nearest_locus",
    "nearest_gene",
    "nearest_tier",
    "nearest_dist",
    "reachable_n3",
    "reachable_n5",
    "reachable_n7",
    "testable",
    "reason",
    "placement_tier",
    "species_match",
]


load = bi.load_tsv
sv = bi.strip_version


def set_reachable(out, dist):
    """Set reachable_n3/5/7 flags from a numeric gene distance (or all-false if dist is None)."""
    for n in NS:
        out[f"reachable_n{n}"] = str(dist is not None and dist <= n).lower()


def main():
    idx = bi.load_from_tsv()
    gold = load(GOLD)
    instances = load(INSTANCES)
    answer = load(ANSWER_KEY)
    t1c = {r["uniprot"]: r for r in load(T1_CEILING)}

    # (ss_type, version-stripped genome, sys_instance_label) -> instance_id
    by_key = {
        (r["ss_type"], sv(r["refseq_genome"]), r["sys_instance_label"].strip()): r["instance_id"] for r in instances
    }
    # (ss_type, version-stripped genome) -> set of instance_ids (for unique-genome auto-assign)
    by_geno = defaultdict(set)
    for r in instances:
        by_geno[(r["ss_type"], sv(r["refseq_genome"]))].add(r["instance_id"])

    # instance_id -> anchored machinery loci [(locus_tag, gene, resolved_accession, tier)]
    mach = defaultdict(list)
    for r in answer:
        if r.get("locus_tag", "").strip():
            mach[r["instance_id"]].append(
                (r["locus_tag"].strip(), r.get("gene", ""), r.get("resolved_accession", ""), r.get("match_tier", ""))
            )

    def resolve_instance(r):
        ss, geno = r["ss_type"], sv(r["refseq_genome"])
        direct = by_key.get((ss, geno, r["sys_instance_id"].strip()))
        if direct:
            return direct, ""
        cand = by_geno.get((ss, geno), set())
        if len(cand) == 1:  # unambiguous own-instance
            return next(iter(cand)), ""
        if len(cand) == 0:
            return None, "no_instance_in_genome"
        return None, "own_instance_unknown"  # multi-instance, net-new: not guessed

    rows = []
    for r in gold:
        out = {f: "" for f in FIELDS}
        out.update(
            gene=r["gene"],
            uniprot=r["uniprot"],
            ss_type=r["ss_type"],
            refseq_genome=r["refseq_genome"],
            ceiling_source=r["ceiling_source"],
            placement_tier=r.get("placement_tier", ""),
            species_match=r.get("species_match", ""),
        )

        if r["ceiling_source"] == "adjacency_scan":
            out["instance_id"] = r["sys_instance_id"]
            if r.get("testable") == "no":
                out.update(testable="no", reason=r.get("testable_reason", "no_genome"))
                rows.append(out)
                continue
            cel = t1c[r["uniprot"]]
            out["effector_locus"] = cel["effector_locus"]
            out["nearest_locus"] = cel.get("abc_locus") or cel.get("mfp_locus") or ""
            out["nearest_gene"] = cel.get("abc_gene") or cel.get("mfp_gene") or ""
            out["nearest_tier"] = "adjacency"
            nd = cel["nearest_component_dist"]
            out["testable"] = "yes"
            if nd == "":  # NOT_FOUND -> machinery present in genome but
                out["nearest_dist"] = ""  # not in the operon (trans-secreted): impossible
                set_reachable(out, None)
                out["reason"] = "transporter_not_adjacent"
            else:
                d = int(nd)
                out["nearest_dist"] = d
                set_reachable(out, d)
            rows.append(out)
            continue

        # answer_key path
        if r["refseq_genome"].strip() in ("", "-"):
            out.update(testable="no", reason="no_genome")
            rows.append(out)
            continue
        iid, why = resolve_instance(r)
        out["instance_id"] = iid or ""
        if iid is None:
            out.update(testable="no", reason=why)
            rows.append(out)
            continue
        loci = mach.get(iid, [])
        if not loci:
            out.update(testable="no", reason="machinery_unanchored")
            rows.append(out)
            continue

        # locate the effector: by locus_tag, else by unique /gene symbol (the corpus tag
        # scheme can be absent from the RefSeq assembly, e.g. Yersinia YE_pYV#### on pYV).
        eff = idx.find(r["refseq_genome"], r["locus_tag"])
        match = "locus_tag"
        if eff is None:
            eff = idx.find_by_gene(r["refseq_genome"], r["gene"])
            match = "gene_symbol"
        if eff is None:
            out.update(testable="no", reason="effector_locus_not_found")
            rows.append(out)
            continue
        eff_rec, eff_ord, eff_cds = eff
        out["effector_locus"] = eff_cds["locus_tag"]
        out["effector_ordinal"] = eff_ord
        out["effector_match"] = match

        best = None  # (dist, locus, gene, tier)
        for lt, gname, racc, tier in loci:
            m = idx.find(r["refseq_genome"], lt)
            if m is None or m[0] != eff_rec:
                continue  # machinery unfindable or on a different replicon than the effector
            d = abs(eff_ord - m[1])
            if best is None or d < best[0]:
                best = (d, lt, gname, tier)
        out["testable"] = "yes"
        if best is None:  # anchored, but none on the effector's replicon
            out["nearest_dist"] = ""
            set_reachable(out, None)
            out["reason"] = "machinery_off_replicon"
        else:
            d, lt, gname, tier = best
            out.update(nearest_dist=d, nearest_locus=lt, nearest_gene=gname, nearest_tier=tier)
            set_reachable(out, d)
        rows.append(out)

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # report
    from collections import Counter

    testable = [r for r in rows if r["testable"] == "yes"]
    untest = [r for r in rows if r["testable"] == "no"]
    print(f"wrote {OUT.relative_to(BENCH)}  ({len(rows)} effectors)")
    bygene = [r for r in testable if r["effector_match"] == "gene_symbol"]
    print(f"  testable   : {len(testable)}  (located by gene-symbol fallback: {len(bygene)})")
    if bygene:
        print("    fallback placements:", ", ".join(f"{r['gene']}({r['ss_type']})" for r in bygene))
    print(f"  untestable : {len(untest)}  reasons {dict(Counter(r['reason'] for r in untest))}")
    print("\n  per SS type (testable / reachable@3 / @5 / @7):")
    for ss in ("T1SS", "T2SS", "T3SS", "T4SS", "T6SS"):
        t = [r for r in testable if r["ss_type"] == ss]
        if not t:
            continue
        r3 = sum(r["reachable_n3"] == "true" for r in t)
        r5 = sum(r["reachable_n5"] == "true" for r in t)
        r7 = sum(r["reachable_n7"] == "true" for r in t)
        u = sum(r["ss_type"] == ss for r in untest)
        print(
            f"    {ss}: testable={len(t):3d}  reachable {r3:3d}/{r5:3d}/{r7:3d}  "
            f"({100 * r3 / len(t):4.0f}%/{100 * r5 / len(t):4.0f}%/{100 * r7 / len(t):4.0f}%)  untestable={u}"
        )


if __name__ == "__main__":
    raise SystemExit(main())
