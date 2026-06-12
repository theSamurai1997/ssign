"""Dataset group 4, task 4.2: pair-features + system-features per positive.

A "pair-feature" is the gene-distance from an effector to the nearest machinery locus of its
ASSIGNED system instance (|ordinal difference| on the same replicon), plus the proximity-window
reachability flags the benchmark uses (within +/-3/5/7 genes). System-features are the SS type
and the instance's machinery component count.

Three branches, by how the instance distance is obtained:
  - validated  -> already computed by the benchmark; pulled from ceiling_per_effector.tsv
                  (nearest_dist / nearest_tier / nearest_locus / reachable_n3,n5,n7). pair_source=ceiling.
  - predicted  -> instance is a machinery-answer-key id; compute the distance here from the
    (instanced)   gene-order index to each of that instance's machinery loci, take the min.
                  pair_source=computed, or =unreachable if no machinery locus is on the same
                  resolved replicon (e.g. genome absent from the index).
  - type-level / self_secreted / no instance -> no instance to measure against; pair-features
                  null. pair_source = none (instance-unknown) or self (T5a/c/d/e autotransporter).

This is the only group-4 task that does NOT need the Phase-2 ssign run output: it is a function
of the labels + the benchmark's own gene-order index and machinery answer key. The per-protein
tool signals (4.1) and the PU unlabeled set (4.4) are still gated on the runs.

Inputs:  data/dataset/positives_all.tsv, data/phase1/ceiling_per_effector.tsv,
         data/machinery/machinery_answer_key.tsv, data/phase1/gene_order_index.tsv
Output:  data/dataset/pair_features.tsv
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from bench_index import load_from_tsv
from bench_io import read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "dataset"
WINDOWS = (3, 5, 7)

OUT_COLS = [
    "protein_id",
    "gene",
    "uniprot",
    "locus_tag",
    "refseq_genome",
    "ss_type",
    "evidence_tier",
    "instance_source",
    "type_level",
    "sys_instance_id",
    "n_machinery",
    "nearest_dist",
    "nearest_tier",
    "nearest_locus",
    "reachable_n3",
    "reachable_n5",
    "reachable_n7",
    "pair_source",
]


def _protein_id(r: dict) -> str:
    """Stable per-protein key (also the reference id the ESM embedding step 4.3 will cache on)."""
    return (r.get("uniprot") or "").strip() or (r.get("locus_tag") or "").strip() or f"{r['gene']}|{r['refseq_genome']}"


def _null_pair(row: dict, source: str) -> dict:
    row.update(
        nearest_dist="",
        nearest_tier="",
        nearest_locus="",
        reachable_n3="",
        reachable_n5="",
        reachable_n7="",
        pair_source=source,
    )
    return row


def _set_pair(row: dict, dist: int, tier: str, locus: str, source: str) -> dict:
    row["nearest_dist"], row["nearest_tier"], row["nearest_locus"] = dist, tier, locus
    for w in WINDOWS:
        row[f"reachable_n{w}"] = "true" if dist <= w else "false"
    row["pair_source"] = source
    return row


def main() -> int:
    pos = read_tsv(OUT / "positives_all.tsv")
    ceiling = read_tsv(ROOT / "data" / "phase1" / "ceiling_per_effector.tsv")
    answer_key = read_tsv(ROOT / "data" / "machinery" / "machinery_answer_key.tsv")

    ceil_by = {(r["gene"], r["uniprot"], r["refseq_genome"]): r for r in ceiling}
    machinery = defaultdict(list)  # instance_id -> [machinery locus_tag, ...]
    for r in answer_key:
        if r.get("locus_tag", "").strip():
            machinery[r["instance_id"]].append(r["locus_tag"].strip())
    idx = load_from_tsv()

    rows = []
    for p in pos:
        row = {
            "protein_id": _protein_id(p),
            "gene": p["gene"],
            "uniprot": p.get("uniprot", ""),
            "locus_tag": p.get("locus_tag", ""),
            "refseq_genome": p["refseq_genome"],
            "ss_type": p["ss_type"],
            "evidence_tier": p["evidence_tier"],
            "instance_source": p["instance_source"],
            "type_level": p["type_level"],
            "sys_instance_id": p["sys_instance_id"],
            "n_machinery": len(machinery.get(p["sys_instance_id"], [])),
        }
        inst = p["sys_instance_id"].strip()

        if p["instance_source"] == "self":  # T5a/c/d/e autotransporter is its own system
            rows.append(_null_pair(row, "self"))
        elif p["type_level"] == "yes" or not inst:  # instance-unknown -> no pair to measure
            rows.append(_null_pair(row, "none"))
        elif p["evidence_tier"] == "validated":  # benchmark already computed it
            c = ceil_by.get((p["gene"], p["uniprot"], p["refseq_genome"]))
            if c and c["nearest_dist"].strip():
                rows.append(_set_pair(row, int(c["nearest_dist"]), c["nearest_tier"], c["nearest_locus"], "ceiling"))
            else:
                rows.append(_null_pair(row, "unreachable"))
        else:  # predicted with an answer-key instance: compute distance to its machinery
            best = None
            for mloc in machinery.get(inst, []):
                d = idx.gene_distance(p["refseq_genome"], p["locus_tag"], mloc)
                if d is not None and (best is None or d < best[0]):
                    best = (d, mloc)
            rows.append(_set_pair(row, best[0], "", best[1], "computed") if best else _null_pair(row, "unreachable"))

    write_tsv(OUT / "pair_features.tsv", OUT_COLS, rows)

    src = Counter(r["pair_source"] for r in rows)
    have = sum(1 for r in rows if r["nearest_dist"] != "")
    print(f"positives: {len(rows)}   pair-feature present: {have}   null: {len(rows) - have}")
    print(f"  pair_source: {dict(sorted(src.items()))}")
    for tier in ("validated", "predicted"):
        sub = [r for r in rows if r["evidence_tier"] == tier]
        h = sum(1 for r in sub if r["nearest_dist"] != "")
        print(
            f"  {tier}: {h}/{len(sub)} with distance   {dict(sorted(Counter(r['pair_source'] for r in sub).items()))}"
        )
    reach = {w: sum(1 for r in rows if r.get(f"reachable_n{w}") == "true") for w in WINDOWS}
    print(f"  reachable within +/- {dict(reach)}  (of {have} placed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
