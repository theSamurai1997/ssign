#!/usr/bin/env python3
"""Merge the agent deep-verify batch outputs into one per-row trust verdict over the 673 survivors.

The agent pass wrote one JSON per batch under data/dataset/pass2_full_raw/. This script aligns each
batch's results back to its input effectors (positionally, validating the gene at each position, since
homolog rows like apxIA x27 share gene+blank-instance and differ only by organism), normalizes the
two old-schema batches, fills deterministic fallbacks for the handful the agents could not cover
(Coxiella select-agent paper blocked by a content filter; 2 rows the agents dropped), assigns a final
trust tier per row, and cross-tabs trust against evidence_tier and ss_type.

It does NOT mutate positives_all.tsv -- that is the apply step (47), gated on the keep/drop policy for
the REFUTED classes (a predicted-tier homolog cited to a single-organism family paper is an expected
provenance pattern, not necessarily a wrong label; a 'validated' row refuted wrong_system is a real defect).

Inputs : data/dataset/deepverify_input/batch_*.json   (agent inputs; carry organism)
         data/dataset/pass2_full_raw/batch_*.json      (agent outputs; status/reason/source/quote)
         data/dataset/citation_consistency_full.tsv    (pass-1 verdict, for fallback)
         data/dataset/positives_all.tsv                (evidence_tier)
Output : data/dataset/deepverify_results_full.tsv      (one row per survivor + trust_tier)
Run    : python3 scripts/46_merge_deepverify.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "data" / "dataset"
INDIR = DATASET / "deepverify_input"
OUTDIR = DATASET / "pass2_full_raw"
sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

# Batches that were split into children; their own output never existed -> skip, use the children.
SUPERSEDED = {"batch_18.json", "batch_18c.json"}


def infer_reason(note: str) -> str:
    """Old-schema NOT_SUPPORTED rows (batches 00/04) carry no reason field; recover it from the note."""
    n = (note or "").lower()
    if "organism" in n or "strain" in n or "homolog" in n:
        return "wrong_organism"
    if "system" in n or "t3ss" in n or "t6ss" in n or "t4ss" in n:
        return "wrong_system"
    if "protein" in n or "different protein" in n:
        return "wrong_protein"
    return "unstructured"


def trust_tier(status: str, reason: str, source: str) -> str:
    if status == "SUPPORTED":
        return "verified_paper" if source == "cited_paper" else "verified_external"
    if status == "INACCESSIBLE":
        return "unverifiable"
    if status in ("REFUTED", "NOT_SUPPORTED"):
        return "refuted_" + (reason or "unstructured")
    return "unknown"


def main() -> int:
    # pass-1 verdict by (gene, inst, organism) for the deterministic fallback
    p1 = {
        (r["gene"], r["sys_instance_id"], r["organism"]): r["verdict"]
        for r in read_tsv(DATASET / "citation_consistency_full.tsv")
    }
    # evidence_tier by the same key
    etier = {
        ((p.get("gene") or "").strip(), (p.get("sys_instance_id") or "").strip(), (p.get("organism") or "").strip()): (
            p.get("evidence_tier") or ""
        ).strip()
        for p in read_tsv(DATASET / "positives_all.tsv")
    }

    batch_files = sorted(INDIR.glob("batch_*.json"))
    rows: list[dict] = []
    for inf in batch_files:
        if inf.name in SUPERSEDED:
            continue
        inp = json.loads(inf.read_text())
        effectors = [(p["doi"], p["ss_type"], e) for p in inp["papers"] for e in p["effectors"]]
        outf = OUTDIR / inf.name
        results = json.loads(outf.read_text())["results"] if outf.exists() else []

        aligned = len(results) == len(effectors) and all(
            (results[i].get("gene") == effectors[i][2]["gene"]) for i in range(len(results))
        )
        for i, (doi, ss, e) in enumerate(effectors):
            key = (e["gene"], e["sys_instance_id"], e["organism"])
            res = results[i] if (aligned and i < len(results)) else _match(results, e)
            if res is None:  # agent dropped this row or whole batch blocked -> deterministic fallback
                rows.append(_fallback_row(doi, ss, e, p1, etier))
                continue
            status = res.get("status", "")
            reason = res.get("reason", "") or (infer_reason(res.get("note", "")) if status == "NOT_SUPPORTED" else "")
            rows.append(
                {
                    "gene": e["gene"],
                    "uniprot": e["uniprot"],
                    "locus_tag": e["locus_tag"],
                    "sys_instance_id": e["sys_instance_id"],
                    "ss_type": ss,
                    "organism": e["organism"],
                    "evidence_tier": etier.get(key, ""),
                    "sourcing_doi": doi,
                    "status": status,
                    "reason": reason,
                    "source": res.get("source", ""),
                    "trust_tier": trust_tier(status, reason, res.get("source", "")),
                    "quote": (res.get("quote", "") or "").replace("\t", " ").replace("\n", " ")[:300],
                    "note": (res.get("note", "") or "").replace("\t", " ").replace("\n", " ")[:200],
                }
            )

    if not rows:
        print("no survivor rows found; nothing to merge")
        return 1
    write_tsv(DATASET / "deepverify_results_full.tsv", list(rows[0].keys()), rows)

    tiers = Counter(r["trust_tier"] for r in rows)
    print(f"merged {len(rows)} survivor rows from {len(batch_files) - len(SUPERSEDED)} batches")
    print("\ntrust tiers:")
    for t, c in tiers.most_common():
        print(f"  {c:4d}  {t}")
    print("\ntrust x evidence_tier:")
    for tier in ("validated", "predicted"):
        sub = Counter(r["trust_tier"] for r in rows if r["evidence_tier"] == tier)
        print(f"  {tier}: {dict(sub)}")
    print("\nverified (paper+external) by ss_type:")
    for ss in sorted({r["ss_type"] for r in rows}):
        sub = [r for r in rows if r["ss_type"] == ss]
        v = sum(1 for r in sub if r["trust_tier"].startswith("verified"))
        print(f"  {ss}: {v}/{len(sub)} verified")
    return 0


def _match(results: list[dict], e: dict):
    """Fallback when positional alignment fails: match by (gene, instance)."""
    for r in results:
        if r.get("gene") == e["gene"] and r.get("sys_instance_id") == e["sys_instance_id"]:
            return r
    return None


def _fallback_row(doi: str, ss: str, e: dict, p1: dict, etier: dict) -> dict:
    key = (e["gene"], e["sys_instance_id"], e["organism"])
    pass1 = p1.get(key, "")
    return {
        "gene": e["gene"],
        "uniprot": e["uniprot"],
        "locus_tag": e["locus_tag"],
        "sys_instance_id": e["sys_instance_id"],
        "ss_type": ss,
        "organism": e["organism"],
        "evidence_tier": etier.get(key, ""),
        "sourcing_doi": doi,
        "status": "BLOCKED",
        "reason": "agent_unavailable",
        "source": "none",
        "trust_tier": "fallback_" + pass1.lower(),
        "quote": "",
        "note": "deep-verify unavailable (content filter / dropped); pass-1 verdict " + pass1,
    }


if __name__ == "__main__":
    raise SystemExit(main())
