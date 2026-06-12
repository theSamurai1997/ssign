#!/usr/bin/env python3
"""
10_resolve_machinery.py  (Phase 0b, task 3.5)

Map each literature-named machinery gene (from data/machinery_raw/<id>.json) to a
RefSeq locus_tag + coordinates in the genome of its system instance.

Design rules (locked):
  - RefSeq is used here ONLY as a coordinate / ID lookup, never to discover machinery.
  - Matching is alias-aware: a paper gene name is matched, case-insensitively, against
    the /gene, /gene_synonym, /old_locus_tag and /locus_tag qualifiers of every gene/CDS
    feature in the organism.
  - Replicon-wide: matching ranges over ALL cached replicons that share the instance
    genome's ORGANISM annotation, not just the instance's own accession. (Brucella virB
    machinery is on chromosome II while its effectors are on chromosome I; Phase 1 then
    scores those effectors cross-replicon = impossible.)
  - Two match tiers, kept separate so the 3.6 verification pass can scrutinise them:
      tier 1  alias   : exact (case-insensitive) hit on a name qualifier  -> high confidence
      tier 2  product : whole-word hit of the gene token inside /product  -> needs review
    Genes matching neither are flagged match_tier=none (the "non-resolving" set).

Inputs:
  data/machinery_raw/*.json        90 curated instances
  data/machinery/instances.tsv     instance_id -> refseq_genome (version-less accession)
  data/refseq_cache/*.gb           61 cached GenBank replicons (gbwithparts)

Output:
  data/machinery/machinery_resolved.tsv   one row per (instance, machinery gene)

Run:
  python scripts/10_resolve_machinery.py
"""

from __future__ import annotations

import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

BENCH = Path(__file__).resolve().parents[1]
CACHE_DIR = BENCH / "data" / "refseq_cache"
RAW_DIR = BENCH / "data" / "machinery_raw"
INSTANCES_TSV = BENCH / "data" / "machinery" / "instances.tsv"
OUT_TSV = BENCH / "data" / "machinery" / "machinery_resolved.tsv"

NAME_QUALS = ("gene", "gene_synonym", "old_locus_tag", "locus_tag")

FIELDS = [
    "instance_id",
    "ss_type",
    "organism",
    "instance_accession",
    "status",
    "gene",
    "family",
    "role",
    "match_tier",
    "match_field",
    "matched_name",
    "locus_tag",
    "resolved_accession",
    "start",
    "end",
    "strand",
    "cross_replicon",
    "n_candidates",
    "candidates",
    "doi",
    "pmid",
]


def strip_version(acc: str) -> str:
    return acc.split(".")[0].strip()


def load_genomes():
    """Parse every cached replicon once.

    Returns:
      acc2org:   version-stripped accession -> organism string (from the replicon
                 that the instance accession points at)
      org2feats: organism -> list of feature dicts across all its replicons
                 {accession, locus_tag, start, end, strand, aliases:set(lowercased), product}
    """
    acc2org = {}
    org2feats = defaultdict(list)
    for gb in sorted(CACHE_DIR.glob("*.gb")):
        rec = SeqIO.read(str(gb), "genbank")
        organism = rec.annotations.get("organism", rec.id)
        # register both the versioned id and the version-stripped form so a
        # version-less instances.tsv accession still resolves
        acc2org[strip_version(rec.id)] = organism
        acc2org[strip_version(rec.name)] = organism

        # merge gene + CDS features that share a locus_tag into one record
        by_locus = {}
        for feat in rec.features:
            if feat.type not in ("gene", "CDS"):
                continue
            q = feat.qualifiers
            lt = q.get("locus_tag", [None])[0]
            key = lt or f"{rec.id}:{int(feat.location.start)}"
            entry = by_locus.setdefault(
                key,
                {
                    "accession": rec.id,
                    "locus_tag": lt,
                    "start": int(feat.location.start) + 1,  # 1-based, inclusive
                    "end": int(feat.location.end),
                    "strand": feat.location.strand,
                    "aliases": set(),
                    "product": "",
                },
            )
            for ql in NAME_QUALS:
                for val in q.get(ql, []):
                    entry["aliases"].add(val.lower())
            prod = q.get("product", [""])[0]
            if prod and not entry["product"]:
                entry["product"] = prod
        org2feats[organism].extend(by_locus.values())
    return acc2org, org2feats


def build_alias_index(feats):
    """alias(lowercased) -> list of feature dicts (one organism's features)."""
    idx = defaultdict(list)
    for f in feats:
        for a in f["aliases"]:
            idx[a].append(f)
    return idx


# whole-word token match inside a /product string, case-insensitive
def product_hits(token, feats):
    pat = re.compile(rf"(?<![A-Za-z0-9]){re.escape(token)}(?![A-Za-z0-9])", re.IGNORECASE)
    return [f for f in feats if f["product"] and pat.search(f["product"])]


# A paper gene carries two usable names: the native symbol (g['gene'], e.g. EpsC)
# and the unified family symbol (g['family'], e.g. GspC). RefSeq may annotate either,
# so resolution tries both. The family field can be descriptive ("SctW/translocon
# (tip)"), so we keep only its leading symbol-like token.
def candidate_names(g):
    names = []
    gene = g.get("gene", "").strip()
    if gene:
        names.append(gene)
    fam = g.get("family", "").strip()
    m = re.match(r"^([A-Za-z][A-Za-z0-9]{1,7})\b", fam)
    if m and m.group(1).lower() not in (n.lower() for n in names):
        names.append(m.group(1))
    return names


def main():
    inst2acc = {}
    with open(INSTANCES_TSV) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            inst2acc[row["instance_id"]] = strip_version(row["refseq_genome"])

    acc2org, org2feats = load_genomes()
    alias_idx_by_org = {org: build_alias_index(feats) for org, feats in org2feats.items()}

    rows = []
    per_inst = defaultdict(
        lambda: {"genes": 0, "alias": 0, "product": 0, "none": 0, "cross": 0, "status": "", "ss_type": "", "org": ""}
    )

    for jf in sorted(RAW_DIR.glob("*.json")):
        with open(jf) as fh:
            d = json.load(fh)
        iid = d["instance_id"]
        inst_acc = inst2acc.get(iid)
        organism = acc2org.get(inst_acc)
        feats = org2feats.get(organism, [])
        idx = alias_idx_by_org.get(organism, {})
        per_inst[iid].update(status=d["status"], ss_type=d["ss_type"], org=organism or "?")

        for g in d["machinery_genes"]:
            name = g["gene"].strip()
            names = candidate_names(g)
            per_inst[iid]["genes"] += 1

            # tier 1: alias hit on any candidate name; tier 2: /product word-match
            cands, tier, field, matched = [], "none", "", ""
            for nm in names:
                hit = idx.get(nm.lower(), [])
                if hit:
                    cands, tier, field, matched = hit, "alias", "name_qualifier", nm
                    break
            if not cands:
                for nm in names:
                    hit = product_hits(nm, feats)
                    if hit:
                        cands, tier, field, matched = hit, "product", "product", nm
                        break

            # dedup, then order: same-replicon first, then by coordinate
            seen, uniq = set(), []
            for c in cands:
                k = (c["accession"], c["locus_tag"])
                if k not in seen:
                    seen.add(k)
                    uniq.append(c)
            uniq.sort(key=lambda c: (strip_version(c["accession"]) != inst_acc, c["start"]))

            # base fields shared by resolved and non-resolving rows
            row = {f: "" for f in FIELDS}
            row.update(
                instance_id=iid,
                ss_type=d["ss_type"],
                organism=organism or "",
                instance_accession=inst_acc or "",
                status=d["status"],
                gene=name,
                family=g.get("family", ""),
                role=g.get("role", ""),
                match_tier=tier,
                match_field=field,
                matched_name=matched,
                n_candidates=len(uniq),
                doi=g.get("doi", ""),
                pmid=g.get("pmid", ""),
            )

            if uniq:
                per_inst[iid][tier] += 1
                best = uniq[0]
                cross = strip_version(best["accession"]) != inst_acc
                if cross:
                    per_inst[iid]["cross"] += 1
                row.update(
                    locus_tag=best["locus_tag"] or "",
                    resolved_accession=best["accession"],
                    start=best["start"],
                    end=best["end"],
                    strand=best["strand"],
                    cross_replicon=str(cross).lower(),
                    candidates=";".join(f"{c['accession']}:{c['locus_tag']}" for c in uniq[1:]),
                )
            else:
                per_inst[iid]["none"] += 1
            rows.append(row)

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # ---- report ----
    tot = len(rows)
    n_alias = sum(1 for r in rows if r["match_tier"] == "alias")
    n_prod = sum(1 for r in rows if r["match_tier"] == "product")
    n_none = sum(1 for r in rows if r["match_tier"] == "none")
    n_cross = sum(1 for r in rows if r["cross_replicon"] == "true")
    n_ambig = sum(1 for r in rows if int(r["n_candidates"] or 0) > 1)
    zero_inst = [iid for iid, s in per_inst.items() if s["alias"] + s["product"] == 0]

    print(f"\nwrote {OUT_TSV.relative_to(BENCH)}")
    print(f"machinery genes total : {tot}")
    print(f"  alias-resolved      : {n_alias} ({n_alias / tot:.0%})")
    print(f"  product-resolved    : {n_prod} ({n_prod / tot:.0%})  [needs 3.6 review]")
    print(f"  resolved (any tier) : {n_alias + n_prod} ({(n_alias + n_prod) / tot:.0%})")
    print(f"  NON-RESOLVING       : {n_none} ({n_none / tot:.0%})")
    print(f"  cross-replicon hits : {n_cross}")
    print(f"  ambiguous (>1 cand) : {n_ambig}")
    print(f"\ninstances total       : {len(per_inst)}")
    print(f"  with >=1 gene placed: {len(per_inst) - len(zero_inst)}")
    print(f"  with ZERO placed    : {len(zero_inst)}  <- need attention")
    for iid in sorted(zero_inst):
        s = per_inst[iid]
        print(f"      {iid:10s} {s['ss_type']:5s} {s['status']:8s} {s['genes']:3d} genes  {s['org']}")

    # per-SS-type resolve rate
    print("\nper SS type (resolved / total genes):")
    by_type = defaultdict(lambda: [0, 0])
    for r in rows:
        by_type[r["ss_type"]][1] += 1
        if r["match_tier"] != "none":
            by_type[r["ss_type"]][0] += 1
    for t in sorted(by_type):
        got, tt = by_type[t]
        print(f"  {t:6s} {got:4d}/{tt:<4d} ({got / tt:.0%})")


if __name__ == "__main__":
    sys.exit(main())
