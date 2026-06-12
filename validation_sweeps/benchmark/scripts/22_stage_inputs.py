#!/usr/bin/env python3
"""
22_stage_inputs.py  (Phase 2 task 6.1: stage ssign input FASTAs for the panel)

ssign is run per assembly. The benchmark panel is the set of assemblies carrying a gold
effector (answer_key) plus the 15 representative genomes the orphan T1SS effectors were
placed into (adjacency_scan). Two things make the staging non-trivial:

  1. Multi-replicon assemblies. Some effectors and their OWN-instance machinery sit on
     different replicons of the same assembly: Brucella effectors on chr I but the virB
     T4SS on chr II; Yersinia pYV effectors + Ysc apparatus on the plasmid (NC_008791.1)
     alongside the chromosome; Vibrio chr I / chr II; Salmonella chromosome + pSLT. If we
     ran one replicon at a time, MacSyFinder could miss a system whose parts span replicons,
     and a chr-II-anchored effector's run would lack its chr-I machinery. So replicons that
     co-occur in an instance (effector genome + every machinery `resolved_accession`) are
     unioned into one run unit, and the input FASTA carries all of them.

  2. Accession-spelling drift. The corpus keys effectors to version-less / mixed-prefix
     accessions (`NC_002516` vs the cached `NC_002516.2`); we resolve each needed replicon
     to its cached file by version+prefix-stripped base.

The input FASTA is the assembly's NUCLEOTIDE sequence (one record per replicon): ssign
re-annotates it with Bakta from scratch, which is the real user path and what the Phase 2
Bakta->RefSeq bridge (task 6.3) assumes. We deliberately do NOT feed the RefSeq CDS
annotations (that would measure a cleaner-than-reality path).

The replicon set per unit is necessary-and-sufficient for recall: every gold effector's own
replicon and its own-instance machinery replicons are present, so no effector is missed for
lack of input. Bystander replicons that carry neither (and were never cached) are absent;
that cannot cause a false negative for these effectors and is recorded in the manifest.

Inputs : data/phase1/effector_gold_set_phase1.tsv
         data/machinery/instances.tsv, data/machinery/machinery_answer_key.tsv
         data/refseq_cache/*.gb
Outputs: inputs/<unit_id>.fasta                 (one per run unit; nucleotide, multi-record)
         data/phase2/panel_manifest.tsv         (unit_id, replicons, n_gold_effectors, kind)
         data/phase2/effector_unit_map.tsv       (uniprot -> unit_id, for task 6.4)
Run:     .venv/bin/python scripts/22_stage_inputs.py
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
import bench_index as bi  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
CACHE = BENCH / "data" / "refseq_cache"
GOLD = BENCH / "data" / "phase1" / "effector_gold_set_phase1.tsv"
INSTANCES = BENCH / "data" / "machinery" / "instances.tsv"
ANSWER_KEY = BENCH / "data" / "machinery" / "machinery_answer_key.tsv"
INPUTS = BENCH / "inputs"
P2 = BENCH / "data" / "phase2"
MANIFEST = P2 / "panel_manifest.tsv"
EFF_MAP = P2 / "effector_unit_map.tsv"


def norm(acc: str) -> str:
    """Version+prefix-stripped accession base (NC_002516.2 -> 002516), for drift-tolerant join."""
    acc = (acc or "").strip()
    acc = re.sub(r"^(NZ_|NC_|NT_|NW_)", "", acc)
    return acc.split(".")[0].lower()


load = bi.load_tsv


class UF:
    """Union-find over accession bases to group replicons into assemblies."""

    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def main():
    # cache base -> cached accession (prefer the longest/ versioned spelling)
    cache_by_base = defaultdict(list)
    for gb in CACHE.glob("*.gb"):
        cache_by_base[norm(gb.stem)].append(gb.stem)
    for base in cache_by_base:
        cache_by_base[base].sort(key=len, reverse=True)  # versioned/prefixed first

    gold = load(GOLD)
    instances = load(INSTANCES)
    mach = load(ANSWER_KEY)

    # per-instance replicon set: effector genome + every anchored machinery replicon
    mach_acc = defaultdict(set)
    for m in mach:
        acc = (m.get("resolved_accession") or "").strip()
        if acc:
            mach_acc[m["instance_id"]].add(acc)

    # collect the replicon set each normalized genome needs (effector genome + its machinery)
    unit_acc = defaultdict(set)  # genome-base -> {raw accessions}
    for i in instances:
        gb = norm(i["refseq_genome"])
        unit_acc[gb].add(i["refseq_genome"].strip())
        unit_acc[gb] |= mach_acc.get(i["instance_id"], set())
    for r in gold:
        if r["ceiling_source"] != "answer_key":
            continue
        g = r["refseq_genome"].strip()
        if g and g != "-":
            unit_acc[norm(g)].add(g)

    # union replicons that co-occur in any unit's set -> assemblies span accession bases
    uf = UF()
    for accs in unit_acc.values():
        bases = [norm(a) for a in accs]
        for b in bases[1:]:
            uf.union(bases[0], b)
    # adjacency placement genomes: each its own node (single replicon)
    placement = {
        r["refseq_genome"].strip()
        for r in gold
        if r["ceiling_source"] == "adjacency_scan" and r.get("testable") == "yes"
    }
    for a in placement:
        uf.find(norm(a))

    # assemble component -> raw accession set
    comp_acc = defaultdict(set)
    for accs in unit_acc.values():
        for a in accs:
            comp_acc[uf.find(norm(a))].add(a)
    for a in placement:
        comp_acc[uf.find(norm(a))].add(a)

    # resolve each component to a stable id (lowest cached accession) + cached files
    INPUTS.mkdir(exist_ok=True)
    P2.mkdir(parents=True, exist_ok=True)
    missing = []
    comp_files = {}  # comp root -> (unit_id, [cached stems])
    for root, accs in comp_acc.items():
        stems = set()
        for a in accs:
            hit = cache_by_base.get(norm(a))
            if not hit:
                missing.append(a)
                continue
            stems.add(hit[0])
        if not stems:
            continue
        unit_id = sorted(stems)[0]
        comp_files[root] = (unit_id, sorted(stems))

    # which assembly each gold effector belongs to (for task 6.4)
    base_to_unit = {}
    for root, (unit_id, stems) in comp_files.items():
        for s in stems:
            base_to_unit[norm(s)] = unit_id

    unit_stems = {unit_id: stems for unit_id, stems in comp_files.values()}

    # Two inputs per unit, so Phase 2 can compare input modes (Teo): a nucleotide FASTA
    # (ssign re-annotates with Bakta) and a GenBank with the RefSeq CDS preserved
    # (ssign --use-input-annotations -> locus_tags ARE the RefSeq tags the gold set is keyed
    # to, measured in the same gene order as the Phase 1 ceiling). The GenBank is the raw
    # cached records concatenated verbatim (no Biopython round-trip, so every qualifier and
    # old_locus_tag survives intact).
    GBFF = BENCH / "inputs_gb"
    GBFF.mkdir(exist_ok=True)
    n_records = {}
    for unit_id, stems in unit_stems.items():
        recs = []
        for stem in stems:
            recs.extend(SeqIO.parse(CACHE / f"{stem}.gb", "genbank"))
        SeqIO.write(recs, INPUTS / f"{unit_id}.fasta", "fasta")
        with open(GBFF / f"{unit_id}.gbff", "wb") as out:
            for stem in stems:
                out.write((CACHE / f"{stem}.gb").read_bytes())
        n_records[unit_id] = len(recs)

    # count gold effectors per unit + emit effector->unit map
    eff_rows = []
    per_unit_eff = defaultdict(int)
    per_unit_ak = defaultdict(int)  # answer_key effectors per unit (curated-machinery)
    for r in gold:
        g = r["refseq_genome"].strip()
        if not g or g == "-":
            continue
        unit = base_to_unit.get(norm(g), "")
        eff_rows.append(
            {
                "uniprot": r["uniprot"],
                "gene": r["gene"],
                "ss_type": r["ss_type"],
                "refseq_genome": g,
                "ceiling_source": r["ceiling_source"],
                "testable": r.get("testable", ""),
                "unit_id": unit,
            }
        )
        if unit:
            per_unit_eff[unit] += 1
            if r["ceiling_source"] == "answer_key":
                per_unit_ak[unit] += 1

    # a unit is a T1SS-rescue "placement" genome only if ALL its gold effectors are
    # adjacency_scan (no curated-machinery effector); union-merged genomes stay "panel".
    with open(MANIFEST, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["unit_id", "n_replicons", "replicons", "n_gold_effectors", "kind"])
        for unit_id in sorted(n_records):
            stems = sorted(unit_stems[unit_id])
            kind = "placement" if per_unit_ak[unit_id] == 0 else "panel"
            w.writerow([unit_id, n_records[unit_id], ",".join(stems), per_unit_eff[unit_id], kind])

    with open(EFF_MAP, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["uniprot", "gene", "ss_type", "refseq_genome", "ceiling_source", "testable", "unit_id"],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(eff_rows)

    unmapped = [r for r in eff_rows if not r["unit_id"]]
    print(f"wrote {len(n_records)} input FASTAs -> {INPUTS.relative_to(BENCH)}/")
    multi = {u: n for u, n in n_records.items() if n > 1}
    print(f"  multi-replicon units: {len(multi)} -> " + ", ".join(f"{u}({n})" for u, n in sorted(multi.items())))
    print(f"  total replicons packed: {sum(n_records.values())}")
    print(
        f"wrote {MANIFEST.relative_to(BENCH)} ({len(n_records)} units) + {EFF_MAP.relative_to(BENCH)} ({len(eff_rows)} effectors)"
    )
    if missing:
        print(f"  MISSING replicons (not in cache): {sorted(set(missing))}")
    if unmapped:
        print(f"  WARNING {len(unmapped)} effectors unmapped to any unit: {[r['uniprot'] for r in unmapped][:10]}")


if __name__ == "__main__":
    sys.exit(main())
