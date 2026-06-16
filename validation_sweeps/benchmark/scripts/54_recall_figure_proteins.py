#!/usr/bin/env python3
"""Per-protein list for EXACTLY the proteins represented in the system-recall figure (06).

The recall figure (52_system_recall.py) is an INSTANCE-level plot: each bar segment counts
system instances, and only TESTABLE instances are shown. So:
  - Not every gold-set effector is in the figure (effectors whose instance is untestable, or
    that map to no instance, are dropped).
  - Not every figure protein is in the gold set: T5SS (19) is self-secreting, has no proximity
    instance, and is NOT in effector_gold_set.tsv (which is T1/2/3/4/6 only). Those come from
    positives_all + the assembled T5SS recall table (53).

This script emits one row per protein that is represented in the figure, and ONLY those:
  T1/T2/T3/T4/T6 : every effector (from actual_per_effector, T3SS from the t3ss-detection tag,
                   others from default) whose instance_id is a TESTABLE instance shown in a bar.
                   Carries the full gold-set provenance + the ssign run result + the instance's
                   bar classification (found / reach_miss / unreach).
  T5SS           : the 19 genes in t5ss_system_recall.tsv, provenance from positives_all,
                   flagged in_gold_set=no, status from the T5SS recall table.

Output: data/phase2/recall_figure_proteins.tsv
Run   : <repo>/.venv/bin/python scripts/54_recall_figure_proteins.py
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import clean_dataset  # noqa: E402
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]
P1 = BENCH / "data" / "phase1"
P2 = BENCH / "data" / "phase2"
GOLD = BENCH / "data" / "gold_build" / "effector_gold_set.tsv"
POS = BENCH / "data" / "dataset" / "positives_all.tsv"
CEILING = P1 / "ceiling_per_effector.tsv"
OUT = P2 / "recall_figure_proteins.tsv"
# Answer-key audit (task #74-79): gene renames applied here; effector quarantines via clean_dataset.
CHANGELOG = P2 / "verification" / "audit_fix_changelog.tsv"

PROX_TYPES = ["T1SS", "T2SS", "T3SS", "T4SS", "T6SS"]  # T5SS handled separately

# Output schema: provenance block, then membership/status, then ssign run signals.
FIELDS = [
    # identity + provenance (gold set / positives_all)
    "ss_type",
    "subtype",
    "gene",
    "uniprot",
    "locus_tag",
    "organism",
    "refseq_genome",
    "sys_instance_id",
    "instance_id",
    "evidence_level",
    "primary_ref",
    "family",
    "length",
    "verification_status",
    "source",
    "gate",
    # figure membership
    "in_gold_set",
    "instance_status",
    "effector_testable",
    "reachable_n3",
    # ssign run result (blank for T5SS — not in the proximity panel table)
    "ssign_call",
    "match_method",
    "match_identity",
    "predicted_localization",
    "dlp_extracellular_prob",
    "dse_ss_type",
    "signalp_prediction",
    "plm_effector_secreted",
    "plm_effector_type",
    "n_tools_agreeing",
    "confidence_tier",
]


def instance_lookup():
    cei = read_tsv(CEILING)
    by_loc = {r["effector_locus"]: r["instance_id"] for r in cei if r.get("effector_locus", "").strip()}
    by_uni = {r["uniprot"]: r["instance_id"] for r in cei if r.get("uniprot", "").strip() and r["uniprot"] != "-"}
    return by_loc, by_uni


def classify_instances(rows, by_loc, by_uni):
    """Replicate 52_system_recall.system_tab: instance -> (testable, status)."""
    inst = defaultdict(lambda: {"found": False, "reach": False, "testable": False})
    eff_iid = {}
    for r in rows:
        iid = by_loc.get(r.get("effector_locus", "")) or by_uni.get(r.get("uniprot", ""))
        eff_iid[id(r)] = iid
        if not iid:
            continue
        v = inst[(r["ss_type"], iid)]
        if r["ssign_call"] == "emitted_secreted":
            v["found"] = True
        if r["testable"] == "yes":
            v["testable"] = True
            if r["reachable_n3"] == "true":
                v["reach"] = True
    status = {}
    for key, v in inst.items():
        if not v["testable"]:
            continue
        status[key] = "found" if v["found"] else ("reach_miss" if v["reach"] else "unreach")
    return status, eff_iid


def main() -> int:
    by_loc, by_uni = instance_lookup()
    gold_rows = read_tsv(GOLD)
    gold = {(r["gene"], r["uniprot"]): r for r in gold_rows}
    gold_by_gene = {r["gene"]: r for r in gold_rows}

    # Answer-key audit propagation: the actual_per_effector backbone carries the PRE-audit identity
    # (old gene/uniprot), so source provenance from the corrected gold row and drop quarantined ones.
    # dropped_id (shared with script 52) is keyed by genome identity, not (ss_type, gene) -- see helper.
    dropped_id = clean_dataset.dropped_id()
    gene_alias = {}  # old gene -> new gene (e.g. the MavF->SdbC rename)
    if CHANGELOG.exists():
        for c in read_tsv(CHANGELOG):
            if c["field"] == "gene" and c["old"] and c["new"]:
                gene_alias[c["old"]] = c["new"]

    def gold_provenance(gene, uniprot):
        """Corrected gold row for an actual-table effector (old uniprot may no longer key gold)."""
        g = gold.get((gene, uniprot)) or gold_by_gene.get(gene)
        if g is None and gene in gene_alias:
            g = gold_by_gene.get(gene_alias[gene])
        return g or {}

    # proximity types: T3SS from the t3ss-detection tag, the rest from default (mirrors the figure)
    dft = clean_dataset.load_clean_actual(P2 / "actual_per_effector.panel_genbank_default.tsv")
    t3t = clean_dataset.load_clean_actual(P2 / "actual_per_effector.panel_genbank_t3ss.tsv")
    status_d, iid_d = classify_instances(dft, by_loc, by_uni)
    status_t, iid_t = classify_instances(t3t, by_loc, by_uni)

    out_rows = []
    for rows, status, iidmap, want in (
        (dft, status_d, iid_d, set(PROX_TYPES) - {"T3SS"}),
        (t3t, status_t, iid_t, {"T3SS"}),
    ):
        for r in rows:
            ss = r["ss_type"]
            if ss not in want:
                continue
            if (r["gene"], r["uniprot"]) in dropped_id or (r["gene"], r.get("effector_locus", "")) in dropped_id:
                continue  # quarantined by the answer-key audit (genome-identity match)
            iid = iidmap[id(r)]
            inst_status = status.get((ss, iid))
            if inst_status is None:  # instance not shown in the figure (untestable / unmapped)
                continue
            g = gold_provenance(r["gene"], r["uniprot"])
            in_gold = bool(g)
            out_rows.append(
                {
                    "ss_type": ss,
                    "subtype": "",
                    "gene": g.get("gene", r["gene"]),
                    "uniprot": g.get("uniprot", r["uniprot"]) if in_gold else r["uniprot"],
                    "locus_tag": g.get("locus_tag", r.get("effector_locus", "")),
                    "organism": g.get("organism", ""),
                    "refseq_genome": g.get("refseq_genome", ""),
                    "sys_instance_id": g.get("sys_instance_id", ""),
                    "instance_id": iid,
                    "evidence_level": g.get("evidence_level", ""),
                    "primary_ref": g.get("primary_ref", ""),
                    "family": g.get("family", ""),
                    "length": g.get("length", ""),
                    "verification_status": g.get("verification_status", ""),
                    "source": g.get("source", ""),
                    "gate": g.get("gate", ""),
                    "in_gold_set": "yes" if in_gold else "no",
                    "instance_status": inst_status,
                    "effector_testable": r.get("testable", ""),
                    "reachable_n3": r.get("reachable_n3", ""),
                    "ssign_call": r.get("ssign_call", ""),
                    "match_method": r.get("match_method", ""),
                    "match_identity": r.get("match_identity", ""),
                    "predicted_localization": r.get("predicted_localization", ""),
                    "dlp_extracellular_prob": r.get("dlp_extracellular_prob", ""),
                    "dse_ss_type": r.get("dse_ss_type", ""),
                    "signalp_prediction": r.get("signalp_prediction", ""),
                    "plm_effector_secreted": r.get("plm_effector_secreted", ""),
                    "plm_effector_type": r.get("plm_effector_type", ""),
                    "n_tools_agreeing": r.get("n_tools_agreeing", ""),
                    "confidence_tier": r.get("confidence_tier", ""),
                }
            )

    # T5SS: 19 from the assembled recall table, provenance from positives_all (NOT in gold set)
    pos_t5 = {r["gene"]: r for r in read_tsv(POS) if r["ss_type"] == "T5SS"}
    for r in read_tsv(P2 / "t5ss_system_recall.tsv"):
        p = pos_t5.get(r["gene"], {})
        out_rows.append(
            {
                "ss_type": "T5SS",
                "subtype": r["subtype"],
                "gene": r["gene"],
                "uniprot": p.get("uniprot", ""),
                "locus_tag": p.get("locus_tag", ""),
                "organism": p.get("organism", ""),
                "refseq_genome": p.get("refseq_genome", ""),
                "sys_instance_id": p.get("sys_instance_id", ""),
                "instance_id": p.get("sys_instance_id", ""),
                "evidence_level": p.get("evidence_level", ""),
                "primary_ref": p.get("primary_ref", ""),
                "family": p.get("family", ""),
                "length": p.get("length", ""),
                "verification_status": p.get("verification_status", ""),
                "source": p.get("source", ""),
                "gate": "",
                "in_gold_set": "no",
                "instance_status": r["status"],
                "effector_testable": "yes",
                "reachable_n3": "true" if r["status"] != "unreach" else "false",
                "ssign_call": "",
                "match_method": "",
                "match_identity": "",
                "predicted_localization": "",
                "dlp_extracellular_prob": "",
                "dse_ss_type": "",
                "signalp_prediction": "",
                "plm_effector_secreted": "",
                "plm_effector_type": "",
                "n_tools_agreeing": "",
                "confidence_tier": "",
            }
        )

    # Collapse same-effector multi-instance duplicates: one gold protein assigned to >1 detected
    # system instance (e.g. EspA/EspZ in C. rodentium hit the real LEE *and* a spurious 2nd T3SS
    # instance, producing identical rows differing only in instance_id). Keep one row per genome
    # identity (ss_type, gene, locus_tag, organism), preferring a reachable 'found' instance over
    # 'unreach'. Only collapse rows with a real locus_tag; the collapse is logged, never silent.
    def _rank(s):
        return 0 if s == "found" else (2 if s == "unreach" else 1)

    best, singles, collapsed = {}, [], []
    for r in out_rows:
        if not r["locus_tag"]:
            singles.append(r)
            continue
        k = (r["ss_type"], r["gene"], r["locus_tag"], r["organism"])
        if k not in best:
            best[k] = r
        elif _rank(r["instance_status"]) < _rank(best[k]["instance_status"]):
            collapsed.append((best[k]["gene"], best[k]["instance_id"]))
            best[k] = r
        else:
            collapsed.append((r["gene"], r["instance_id"]))
    out_rows = singles + list(best.values())

    out_rows.sort(
        key=lambda r: (PROX_TYPES.index(r["ss_type"]) if r["ss_type"] in PROX_TYPES else 5, r["instance_id"], r["gene"])
    )
    write_tsv(OUT, FIELDS, out_rows)
    if collapsed:
        print(f"collapsed {len(collapsed)} duplicate instance-row(s): {collapsed}")

    from collections import Counter

    by_ss = Counter(r["ss_type"] for r in out_rows)
    not_gold = Counter(r["ss_type"] for r in out_rows if r["in_gold_set"] == "no")
    print(f"wrote {OUT.relative_to(BENCH)}: {len(out_rows)} proteins represented in the recall figure")
    for ss in PROX_TYPES + ["T5SS"]:
        print(f"  {ss}: {by_ss.get(ss, 0)} proteins" + (f"  ({not_gold[ss]} not in gold set)" if not_gold[ss] else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
