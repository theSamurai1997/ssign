#!/usr/bin/env python3
"""Phase 0a, steps 2b.3-2b.5: reconcile external DBs against the corpus gold set.

Inputs:
  data/gold_build/gold_set_corpus.tsv            corpus gold set (523)
  data/external_dbs/secret6/secret6_records.jsonl SecReT6 experimental T6SS effectors
  data/external_dbs/secret4/secret4_mapped.json   SecReT4 verified T4SS effectors (locus-mapped, gated)
  data/external_dbs/citation_overrides.tsv       manual 2b.4 citation pass: source_id -> DOI/gene
                                                 (KEEP) or HOLD, for entries whose DB exposed no DOI

Dedup key: locus_tag is the universal cross-DB key (present for ~all external
entries and every gold row); UniProt accession is a secondary key where the
external record exposes a real accession (not a UniParc id or gene-based entry name).

Net-new entries are split by a detection gate (Teo 2026-06-10, "gate on detected
systems"):
  ACTIVE  system type ssign detects (T6SS; T4SS type-IVA VirB/D4) -> folded in now,
          subject to the standard verification bar.
  GATED   T4SS type-IVB (Dot/Icm): genome-dispersed, ssign detection unconfirmed ->
          recorded and held, folded only if Phase 2 shows ssign detects the system.

Standard bar (2b.4) for an ACTIVE net-new entry to be folded: a locus_tag (so Phase 1
can position it) AND a citation (DOI or PMID). Entries lacking a citation are kept
in a NEEDS_CITATION queue, not silently added.

Outputs (data/external_dbs/):
  reconciliation_report.tsv   per-DB tallies + coverage gaps
  net_new_active.tsv          folded ACTIVE net-new (source-tagged)
  net_new_gated.tsv           held GATED net-new (Dot/Icm)
  net_new_needs_citation.tsv  ACTIVE net-new missing a citation
  net_new_apparatus_dropped.tsv  bare Hcp/VgrG/PAAR + immunity dropped (Phase 0a biology rule)
And (data/gold_build/):
  effector_gold_set.tsv       corpus + folded ACTIVE net-new, each tagged source/gate
"""

from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "data" / "gold_build"
EXT = ROOT / "data" / "external_dbs"

ACC_RE = re.compile(r"^([A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}|[OPQ][0-9][A-Z0-9]{3}[0-9])$")

# Apparatus/immunity drop, mirroring scripts/01_build_gold_set.py (keep in sync). The
# corpus dropped structural Hcp/VgrG/PAAR and immunity proteins as non-cargo; net-new
# must apply the same rule or the fold silently re-introduces what 01 removed. For net-new
# the "evolved" effectors (Hcp/VgrG/PAAR fused to a C-terminal toxin) are flagged by a real
# effector subtype, so we only drop a bare apparatus name when its subtype is empty/'-'.
APPARATUS_GENE = re.compile(r"^(hcp|vgrg|paar)", re.IGNORECASE)
IMMUNITY_GENE = re.compile(r"^(tsi|tdi|tldi|eagr|tai)", re.IGNORECASE)


def biology_drop(gene: str, subtype: str) -> str:
    """'' to keep, else a drop reason. Bare apparatus = apparatus name + no toxin subtype."""
    g = (gene or "").strip()
    st = (subtype or "").strip()
    if APPARATUS_GENE.match(g) and st in ("", "-"):
        return "apparatus-as-effector (bare T6SS structural Hcp/VgrG/PAAR)"
    if IMMUNITY_GENE.match(g):
        return "immunity/adaptor protein, not a secreted effector"
    return ""


def load_overrides() -> dict:
    """source_id -> {doi, pmid, gene, verdict, note} from the manual citation pass (2b.4)."""
    p = EXT / "citation_overrides.tsv"
    if not p.exists():
        return {}
    return {r["source_id"]: r for r in csv.DictReader(p.open(), delimiter="\t")}


def norm_loc(s: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", (s or "").upper())


def loc_prefix(s: str) -> str:
    """Locus namespace = leading alpha run of the normalised locus (genome-specific).

    'BAB1_1058' -> 'BAB', 'lpg2526' -> 'LPG', 'PA1844' -> 'PA'. Two loci sharing a
    prefix come from the same genome's annotation, so a prefix match means the entry
    sits in a genome the corpus already curates (in-panel)."""
    m = re.match(r"([A-Z]+)", norm_loc(s))
    return m.group(1) if m else ""


def norm_genome(s: str) -> str:
    """Genome accession without version, with the RefSeq 'NZ_' wrapper stripped so a
    RefSeq accession (NZ_HG326223) matches the INSDC accession it wraps (HG326223)."""
    g = (s or "").split(".")[0].strip().upper()
    return g[3:] if g.startswith("NZ_") else g


def norm_acc(uni: str) -> str:
    """Extract a UniProtKB accession from a raw 'Uniprot ID' field, else ''.

    Handles plain accessions (Q9I2Q1), TrEMBL entry names (Q9HYC2_PSEAE -> Q9HYC2),
    and rejects UniParc ids (UPI...) and gene-based names (VGRG3_VIBCH)."""
    u = (uni or "").strip().upper()
    if not u or u == "-" or u.startswith("UPI"):
        return ""
    cand = u.split("_", 1)[0]
    return cand if ACC_RE.match(cand) else ""


def read_gold() -> list[dict]:
    return list(csv.DictReader((GOLD / "gold_set_corpus.tsv").open(), delimiter="\t"))


def load_secret6() -> list[dict]:
    out = []
    for line in (EXT / "secret6" / "secret6_records.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        out.append(
            {
                "source_db": "SecReT6",
                "source_id": r.get("source_id", ""),
                "ss_type": "T6SS",
                "gate": "ACTIVE",
                "subtype": r.get("effector_type", ""),
                "gene": r.get("gene", ""),
                "organism": "",
                "locus_tag": r.get("locus_tag", ""),
                "uniprot_raw": r.get("uniprot", ""),
                "refseq_genome": r.get("refseq_genome", ""),
                "doi": r.get("doi", ""),
                "pmid": r.get("pmid", ""),
            }
        )
    return out


def load_secret4() -> list[dict]:
    out = []
    for r in json.loads((EXT / "secret4" / "secret4_mapped.json").read_text()):
        out.append(
            {
                "source_db": "SecReT4",
                "source_id": r.get("source_id", ""),
                "ss_type": "T4SS",
                "gate": r.get("gate", "ACTIVE"),
                "subtype": r.get("t4ss_subtype", ""),
                "gene": "",
                "organism": r.get("ncbi_organism", "") or r.get("organism", ""),
                "locus_tag": r.get("locus_tag", ""),
                "uniprot_raw": "",
                "refseq_genome": "",
                "doi": r.get("doi", ""),
                "pmid": "",
            }
        )
    return out


def main() -> int:
    gold = read_gold()
    gold_loc = {norm_loc(r["locus_tag"]) for r in gold if r["locus_tag"].strip()}
    gold_acc = {r["uniprot"].strip().upper() for r in gold if r["uniprot"].strip() and r["uniprot"].strip() != "-"}
    # in-panel matching is per SS type, so a short locus prefix can't collide across types
    gold_prefix: dict[str, set] = {}
    gold_genome: dict[str, set] = {}
    for r in gold:
        st = r["ss_type"]
        if r["locus_tag"].strip():
            gold_prefix.setdefault(st, set()).add(loc_prefix(r["locus_tag"]))
        if r["refseq_genome"].strip():
            gold_genome.setdefault(st, set()).add(norm_genome(r["refseq_genome"]))
    for st in gold_prefix:
        gold_prefix[st].discard("")

    external = load_secret6() + load_secret4()
    overrides = load_overrides()

    # classify each external record
    for e in external:
        st = e["ss_type"]
        nl = norm_loc(e["locus_tag"])
        acc = norm_acc(e["uniprot_raw"])
        e["acc"] = acc
        # manual citation pass (2b.4): attach a verified DOI + gene name to an entry whose
        # source DB exposed none, so it can clear the cite bar; HOLD verdicts attach nothing.
        ov = overrides.get(e["source_id"])
        if ov and (ov.get("verdict") or "").strip() == "KEEP":
            e["doi"] = e["doi"] or ov.get("doi", "")
            e["pmid"] = e["pmid"] or ov.get("pmid", "")
            e["gene"] = e["gene"] or ov.get("gene", "")
        is_dup = (nl and nl in gold_loc) or (acc and acc in gold_acc)
        has_cite = bool(e["doi"]) or bool(e["pmid"])
        has_loc = bool(nl)
        e["drop_reason"] = biology_drop(e["gene"], e["subtype"])
        # in-panel = the entry sits in a genome the corpus already curates. Exact genome
        # accession is authoritative; the locus-prefix fallback is only used when the entry
        # carries no genome (prefix alone over-matches across strains, e.g. PA14 vs PAO1).
        eg = norm_genome(e["refseq_genome"])
        if eg:
            in_panel = eg in gold_genome.get(st, set())
        else:
            in_panel = loc_prefix(e["locus_tag"]) in gold_prefix.get(st, set())
        e["panel_status"] = "IN_PANEL" if in_panel else "PANEL_EXPANSION"
        if is_dup:
            e["class"] = "DUPLICATE"
        elif e["drop_reason"]:
            e["class"] = "DROP_APPARATUS"
        elif e["gate"] == "GATED":
            e["class"] = "NET_NEW_GATED"
        elif has_loc and has_cite:
            e["class"] = "NET_NEW_ACTIVE"
        elif has_loc and not has_cite:
            e["class"] = "NEEDS_CITATION"
        else:
            e["class"] = "REJECT_UNPLACEABLE"

    # dedup within the external set itself (same locus appearing twice) for net-new tables
    def dedup_by_locus(rows: list[dict]) -> list[dict]:
        seen, out = set(), []
        for r in sorted(rows, key=lambda x: (x["source_db"], x["source_id"])):
            k = (r["ss_type"], norm_loc(r["locus_tag"]))
            if k in seen:
                continue
            seen.add(k)
            out.append(r)
        return out

    active = dedup_by_locus([e for e in external if e["class"] == "NET_NEW_ACTIVE"])
    gated = dedup_by_locus([e for e in external if e["class"] == "NET_NEW_GATED"])
    needcite = dedup_by_locus([e for e in external if e["class"] == "NEEDS_CITATION"])
    apparatus = dedup_by_locus([e for e in external if e["class"] == "DROP_APPARATUS"])

    # ---- write net-new tables ----
    cols = [
        "source_db",
        "source_id",
        "ss_type",
        "gate",
        "panel_status",
        "subtype",
        "gene",
        "organism",
        "locus_tag",
        "acc",
        "refseq_genome",
        "doi",
        "pmid",
    ]

    def write(path: Path, rows: list[dict], header: list[str]) -> None:
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in header})

    # in-panel active net-new fold in now; panel-expansion is held for the Checkpoint-A
    # panel decision (folding it would silently expand the genome panel).
    active_in = [e for e in active if e["panel_status"] == "IN_PANEL"]
    active_exp = [e for e in active if e["panel_status"] == "PANEL_EXPANSION"]

    write(EXT / "net_new_active.tsv", active, cols)
    write(EXT / "net_new_panel_expansion.tsv", active_exp, cols)
    write(EXT / "net_new_gated.tsv", gated, cols)
    write(EXT / "net_new_needs_citation.tsv", needcite, cols)
    write(EXT / "net_new_apparatus_dropped.tsv", apparatus, cols + ["drop_reason"])

    # ---- effector_gold_set.tsv = corpus + folded IN_PANEL ACTIVE net-new ----
    GHEAD = list(gold[0].keys())
    final = [dict(r, source="corpus", gate="ACTIVE") for r in gold]
    for e in active_in:
        row = {k: "" for k in GHEAD}
        row.update(
            {
                "gene": e["gene"] or e["source_id"],
                "uniprot": e["acc"],
                "locus_tag": e["locus_tag"],
                "organism": e["organism"],
                "refseq_genome": e["refseq_genome"],
                "ss_type": e["ss_type"],
                "evidence_level": "validated",
                "primary_ref": e["doi"] or (f"PMID:{e['pmid']}" if e["pmid"] else ""),
                "verification_status": "EXTERNAL_DB",
                "verification_notes": f"{e['source_db']} net-new ({e['subtype']})",
            }
        )
        row["source"] = e["source_db"]
        row["gate"] = "ACTIVE"
        final.append(row)
    write(GOLD / "effector_gold_set.tsv", final, GHEAD + ["source", "gate"])

    # ---- reconciliation_report.tsv ----
    rep_rows = []
    for db in ("SecReT4", "SecReT6"):
        sub = [e for e in external if e["source_db"] == db]
        c = Counter(e["class"] for e in sub)
        rep_rows.append(
            {
                "source_db": db,
                "total_parsed": len(sub),
                "duplicate": c.get("DUPLICATE", 0),
                "active_in_panel": len([e for e in active_in if e["source_db"] == db]),
                "active_panel_expansion": len([e for e in active_exp if e["source_db"] == db]),
                "net_new_gated": len([e for e in gated if e["source_db"] == db]),
                "needs_citation": len([e for e in needcite if e["source_db"] == db]),
                "apparatus_dropped": len([e for e in apparatus if e["source_db"] == db]),
                "reject_unplaceable": c.get("REJECT_UNPLACEABLE", 0),
            }
        )
    # coverage gaps as note-only pseudo-rows (blank counts)
    rep_rows.append(
        {"source_db": "EffectiveDB", "note": "excluded: prediction-only, no experimental per-effector citations"}
    )
    rep_rows.append(
        {"source_db": "BastionHub(T1/T2)", "note": "DOWN; only curated T1SS/T2SS source -> T1/T2 net-new pending"}
    )
    write(
        EXT / "reconciliation_report.tsv",
        rep_rows,
        [
            "source_db",
            "total_parsed",
            "duplicate",
            "active_in_panel",
            "active_panel_expansion",
            "net_new_gated",
            "needs_citation",
            "apparatus_dropped",
            "reject_unplaceable",
            "note",
        ],
    )

    # ---- console summary ----
    print("=== reconciliation ===")
    for r in rep_rows[:2]:
        print(
            f"  {r['source_db']}: parsed {r['total_parsed']} | dup {r['duplicate']} | "
            f"active(in-panel {r['active_in_panel']}, expansion {r['active_panel_expansion']}) | "
            f"GATED {r['net_new_gated']} | needs-cite {r['needs_citation']} | reject {r['reject_unplaceable']}"
        )
    print(
        f"\nfolded IN_PANEL active net-new: {len(active_in)}  (by type: {dict(Counter(e['ss_type'] for e in active_in))})"
    )
    print(f"held PANEL_EXPANSION active:    {len(active_exp)}  (new genomes -> Checkpoint A panel decision)")
    print(f"held GATED net-new:             {len(gated)}  (Dot/Icm, pending Phase 2)")
    print(f"needs-citation queue:           {len(needcite)}")
    print(
        f"apparatus dropped (Phase 0a rule): {len(apparatus)}  (by type: {dict(Counter(e['ss_type'] for e in apparatus))})"
    )
    print(f"\neffector_gold_set.tsv: {len(final)} rows ({len(gold)} corpus + {len(active_in)} folded)")
    print(f"  by ss_type: {dict(Counter(r['ss_type'] for r in final))}")
    print("\ncoverage gaps: EffectiveDB (prediction-only, excluded); BastionHub down (T1/T2SS net-new pending)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
