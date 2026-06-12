#!/usr/bin/env python3
"""
15_t1ss_adjacency.py  (T1SS rescue: verify operon adjacency -> ceiling)

For each rescued T1SS effector placed in a genome (script 14), confirm the thing that
makes T1SS special: the effector is encoded next to its own secretion machinery (an
ABC transporter of the HlyB family + an HlyD-family membrane-fusion protein, the TolC
channel being recruited from elsewhere). We do NOT assume this -- we read the gene order
around the placed effector and look for the transporter signature within +/-7 genes.

Why this is enough (and not circular): the HlyB/HlyD families ARE the literature
definition of T1SS machinery, and we read them off the genome annotation, not from
MacSyFinder. So:
  - adjacency confirmed  -> the effector is reachable by the +/-N proximity rule
                            (ceiling = reachable), by construction not by ssign.
  - adjacency NOT found   -> flagged exception (e.g. a trans-encoded transporter);
                            that effector is a genuine 'far/impossible' case, not an
                            ssign failure.
Any later ssign MISS of a confirmed-adjacent effector is therefore a true ssign false
negative, not a structural impossibility. This makes T1SS the cleanest recall test.

Classification reads each neighbour's /gene + /product:
  ABC : type I secretion ATPase/permease, HlyB family, or a known native ATPase gene
        (hlyB prtD aprD cyaB lktB apxIB hasD rsaD lipB ...)
  MFP : membrane-fusion protein, HlyD family, or native adaptor gene
        (hlyD prtE aprE cyaD lktD apxID hasE rsaE ...)
  OMF : TolC-family / type I secretion outer-membrane protein (recorded, not required)

Input : data/t1ss_rescue/t1ss_rescued.tsv  (placed rows: placement_nuc + start/stop)
Genomes fetched via scripts/09_fetch_refseq.py fetch() into data/refseq_cache/
Output: data/t1ss_rescue/t1ss_ceiling.tsv

Run:
  .venv/bin/python scripts/15_t1ss_adjacency.py
"""

from __future__ import annotations

import csv
import importlib
import re
import sys
from pathlib import Path

from Bio import SeqIO

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).parent))
fetch_mod = importlib.import_module("09_fetch_refseq")  # reuse fetch()

RESCUE = BENCH / "data" / "t1ss_rescue"
IN_TSV = RESCUE / "t1ss_rescued.tsv"
CACHE = BENCH / "data" / "refseq_cache"
OUT_TSV = RESCUE / "t1ss_ceiling.tsv"

WINDOW = 7  # scan +/-7 genes; classify reachable at 3/5/7

# native gene symbols of the canonical T1SS transporter triad, by role
ABC_GENES = {"hlyb", "prtd", "aprd", "cyab", "lktb", "apxib", "hasd", "rsad", "lipb", "cvab", "tlid", "siad"}
MFP_GENES = {"hlyd", "prte", "apre", "cyad", "lktd", "apxid", "hase", "rsae", "lipc", "cvaa", "tlic"}
OMF_GENES = {"tolc", "prtf", "aprf", "cyae", "rsafa", "rsafb", "hasf"}

ABC_PROD = re.compile(
    r"type i(?:\s|-)?secretion.*(atpase|permease|abc|atp-binding)"
    r"|hlyb"
    r"|abc.{0,25}(atp-binding|atpase).*secret"
    r"|secretion.*abc.{0,25}(atp|permease)",
    re.IGNORECASE,
)
MFP_PROD = re.compile(
    r"membrane[ -]fusion"
    r"|hlyd"
    r"|type i(?:\s|-)?secretion.*membrane",
    re.IGNORECASE,
)
OMF_PROD = re.compile(
    r"\btolc\b"
    r"|type i(?:\s|-)?secretion.*outer[ -]membrane"
    r"|outer[ -]membrane (channel|protein tolc)",
    re.IGNORECASE,
)


def classify(gene: str, product: str) -> set[str]:
    g = (gene or "").lower()
    roles = set()
    if g in ABC_GENES or ABC_PROD.search(product or ""):
        roles.add("ABC")
    if g in MFP_GENES or MFP_PROD.search(product or ""):
        roles.add("MFP")
    if g in OMF_GENES or OMF_PROD.search(product or ""):
        roles.add("OMF")
    return roles


def ordered_cds(rec):
    """Ordered list of CDS on a replicon: (start, end, locus_tag, gene, product, roles)."""
    feats = []
    for f in rec.features:
        if f.type != "CDS":
            continue
        q = f.qualifiers
        lt = q.get("locus_tag", [""])[0]
        gene = q.get("gene", [""])[0]
        product = q.get("product", [""])[0]
        feats.append(
            {
                "start": int(f.location.start) + 1,
                "end": int(f.location.end),
                "locus_tag": lt,
                "gene": gene,
                "product": product,
                "roles": classify(gene, product),
            }
        )
    feats.sort(key=lambda d: d["start"])
    return feats


def find_effector_idx(feats, start, stop):
    """Index of the CDS overlapping the placement coordinates.

    Requires the overlap to cover >=50% of the placement span, so a coordinate that
    falls in an intergenic gap and only clips a neighbour is rejected (returns None)
    rather than silently anchoring the scan on the wrong CDS. The placement coords are
    the matched protein's own coords on THIS replicon, so a true hit overlaps ~100%.
    """
    span = max(1, stop - start)
    best, best_ov = None, 0
    for i, f in enumerate(feats):
        ov = min(f["end"], stop) - max(f["start"], start)
        if ov > best_ov:
            best, best_ov = i, ov
    if best is None or best_ov < 0.5 * span:
        return None
    return best


FIELDS = [
    "gene",
    "uniprot",
    "placement_nuc",
    "tier",
    "species_match",
    "effector_locus",
    "effector_start",
    "abc_locus",
    "abc_gene",
    "abc_dist",
    "mfp_locus",
    "mfp_gene",
    "mfp_dist",
    "omf_locus",
    "omf_dist",
    "nearest_component_dist",
    "reachable_n3",
    "reachable_n5",
    "reachable_n7",
    "adjacency",
    "note",
]


def scan_adjacency(feats, ei):
    """Scan +/-WINDOW around effector index `ei` for the T1SS transporter signature.

    Single source of truth for the adjacency rule, reused by script 16. Returns
    (status, fields) where fields populates the adjacency-related FIELDS columns.
    status: CONFIRMED (ABC+MFP both nearby) / PARTIAL (one component) / NOT_FOUND.
    """
    eff = feats[ei]
    comp = {"ABC": None, "MFP": None, "OMF": None}  # role -> (dist, feat)
    dump = []
    for d in range(-WINDOW, WINDOW + 1):
        j = ei + d
        if j < 0 or j >= len(feats) or j == ei:
            continue
        f = feats[j]
        if f["roles"]:
            dump.append(f"{d:+d}:{'/'.join(sorted(f['roles']))}:{f['gene'] or f['locus_tag']}")
            for role in f["roles"]:
                if comp[role] is None or abs(d) < comp[role][0]:
                    comp[role] = (abs(d), f)

    fields = {"effector_locus": eff["locus_tag"], "effector_start": eff["start"]}
    for role, key in (("ABC", "abc"), ("MFP", "mfp"), ("OMF", "omf")):
        if comp[role]:
            dist, f = comp[role]
            fields[f"{key}_dist"] = dist
            fields[f"{key}_locus"] = f["locus_tag"]
            if key != "omf":
                fields[f"{key}_gene"] = f["gene"]

    dists = [comp[role][0] for role in ("ABC", "MFP", "OMF") if comp[role]]
    nearest = min(dists) if dists else None
    fields["nearest_component_dist"] = nearest if nearest is not None else ""
    for n in (3, 5, 7):
        fields[f"reachable_n{n}"] = str(nearest is not None and nearest <= n).lower()
    status = "CONFIRMED" if (comp["ABC"] and comp["MFP"]) else ("PARTIAL" if dists else "NOT_FOUND")
    fields["adjacency"] = status
    fields["note"] = " ".join(dump)
    return status, fields


def adjacency_for(nuc, start, stop):
    """Fetch a replicon, locate the effector at coords, return (status, fields)."""
    if not fetch_mod.fetch(nuc):
        return "FETCH_FAILED", {"adjacency": "FETCH_FAILED", "note": "could not fetch placement genome"}
    rec = SeqIO.read(str(CACHE / f"{nuc}.gb"), "genbank")
    feats = ordered_cds(rec)
    ei = find_effector_idx(feats, int(start), int(stop))
    if ei is None:
        return "EFFECTOR_NOT_FOUND", {
            "adjacency": "EFFECTOR_NOT_FOUND",
            "note": f"no CDS covering >=50% of {start}-{stop} on {nuc}",
        }
    return scan_adjacency(feats, ei)


def main():
    placed = [r for r in csv.DictReader(open(IN_TSV), delimiter="\t") if r["final_status"].startswith("placed")]
    print(f"verifying operon adjacency for {len(placed)} placed T1SS effectors\n", file=sys.stderr)

    rows = []
    for r in placed:
        nuc = r["placement_nuc"].strip()
        out = {f: "" for f in FIELDS}
        out.update(
            gene=r["gene"], uniprot=r["uniprot"], placement_nuc=nuc, tier=r["tier"], species_match=r["species_match"]
        )
        status, fields = adjacency_for(nuc, r["start"], r["stop"])
        out.update(fields)
        rows.append(out)
        print(
            f"  {r['gene']:10s} {nuc:20s} eff={out.get('effector_locus') or '?':16s} "
            f"{out['adjacency']:10s} nearest={out['nearest_component_dist']} [{out['note']}]",
            file=sys.stderr,
        )

    with open(OUT_TSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    conf = sum(1 for r in rows if r["adjacency"] == "CONFIRMED")
    part = sum(1 for r in rows if r["adjacency"] == "PARTIAL")
    none = sum(1 for r in rows if r["adjacency"] in ("NOT_FOUND", "EFFECTOR_NOT_FOUND", "FETCH_FAILED"))
    r3 = sum(1 for r in rows if r["reachable_n3"] == "true")
    r5 = sum(1 for r in rows if r["reachable_n5"] == "true")
    r7 = sum(1 for r in rows if r["reachable_n7"] == "true")
    print(f"\nwrote {OUT_TSV.relative_to(BENCH)}")
    print(f"  placed effectors        : {len(rows)}")
    print(f"  adjacency CONFIRMED (ABC+MFP): {conf}")
    print(f"  PARTIAL (one component)      : {part}")
    print(f"  NOT_FOUND / error            : {none}")
    print(f"  reachable @N=3 / 5 / 7       : {r3} / {r5} / {r7}")
    if part or none:
        print("\n  exceptions to inspect:")
        for r in rows:
            if r["adjacency"] != "CONFIRMED":
                print(
                    f"    {r['gene']:10s} {r['uniprot']:8s} {r['placement_nuc']:20s} {r['adjacency']:16s} {r['note']}"
                )


if __name__ == "__main__":
    raise SystemExit(main())
