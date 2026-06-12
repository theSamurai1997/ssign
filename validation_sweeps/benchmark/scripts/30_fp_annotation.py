#!/usr/bin/env python3
"""Phase 2 task 8.3: bound emission precision from above by obvious false positives (annotation).

The DB floor (8.2) only covers T4SS/T6SS. To bound the other direction, and the SS types with no
independent DB, classify each emission by its functional annotation. A protein annotated as a
transcriptional regulator, a chemotaxis protein, or a ribosomal subunit, sitting next to a secretion
system, is almost certainly a proximity bystander, NOT a secreted effector. The clearly-cytoplasmic
fraction is therefore an obvious-false-positive rate; one minus it (plus machinery) is a soft
precision ceiling. This is a heuristic on annotation text, not ground truth, so buckets are reported
transparently, not collapsed to a single number.

Buckets (first matching rule wins, order matters):
  apparatus    - the secretion machinery itself emitted as a substrate (ShlB/FhaC TpsB pore, TssA,
                 Hcp/VgrG/PAAR, phage-tail, pilin/fimbrial). A real error, but a different class.
  housekeeping - clearly cytoplasmic / non-secreted core function (transcriptional regulator,
                 chemotaxis, cell division, ribosome/translation/replication, central metabolism,
                 membrane transporters). The obvious-FP bucket.
  effector     - annotation consistent with a secreted effector (autotransporter, protease/peptidase,
                 toxin, hemolysin, adhesin, lipase, nuclease, two-partner secretion, RTX, HasA...).
  hypothetical - hypothetical / uncharacterized / DUF: genuinely unknown (could be novel effector).
  other        - anything unmatched.

Inputs : data/phase2/emissions.<tag>.tsv
Outputs: data/phase2/emissions_fpclass.<tag>.tsv
Run:     .venv/bin/python scripts/30_fp_annotation.py --run-tag panel_genbank_default
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from bench_io import read_tsv, write_tsv  # noqa: E402

BENCH = Path(__file__).resolve().parents[1]

# Order matters: apparatus before effector so the ShlB/FhaC TpsB pore (annotated "hemolysin
# secretion/activation protein") is read as machinery, not as a hemolysin effector.
RULES = [
    (
        "apparatus",
        (
            "secretion/activation",
            "shlb",
            "fhac",
            "hecb",
            "tpsb",
            "type vi secretion system protein",
            "secretion system-associated",
            "secretion system accessory",
            "contractile sheath",
            "baseplate",
            "tssa",
            "tssb",
            "tssc",
            "tssj",
            "tssk",
            "tssl",
            "tssm",
            "tagh",
            "tagj",
            "hcp",
            "vgrg",
            "paar",
            "phage tail",
            "fimbrial",
            "pilin",
            "pilus",
            "secretin",
            "type iv pili",
            "flagell",
        ),
    ),
    (
        "housekeeping",
        (
            "cytoplasmic protein",
            "transcriptional regulator",
            "lysr",
            "helix-turn-helix",
            "methyl-accepting chemotaxis",
            "chemotaxis",
            "cell division",
            "mfs transporter",
            "abc transporter",
            "permease",
            "ribosomal",
            "ribosome",
            "elongation factor",
            "trna",
            "aminoacyl",
            "dna polymerase",
            "rna polymerase",
            "gyrase",
            "helicase",
            "chaperon",
            "groel",
            "dnak",
            "atp synthase",
            "nadh",
            "dehydrogenase",
            "shikimate",
            "dehydroquinate",
            "ftsz",
            "transketolase",
            "oxidoreductase",
            "isomerase",
            "two-component",
            "sensor histidine",
            "response regulator",
            "h-ns",
            "histone",
            "transposase",
            "integrase",
            "recombinase",
            "penicillin-binding",
            "reductase",
            "carboxylase",
            "decarboxylase",
            "hydratase",
            "aldolase",
            "mutase",
            "pyridoxal phosphate",
            "tonb-dependent receptor",
        ),
    ),
    (
        "effector",
        (
            "autotransporter",
            "toxin",
            "effector",
            "hemolysin",
            "haemolysin",
            "hemagglutinin",
            "adhesin",
            "yada",
            "protease",
            "peptidase",
            "serralysin",
            "metalloprotease",
            "lipase",
            "phospholipase",
            "esterase",
            "nuclease",
            "intimin",
            "invasin",
            "rtx",
            "ig-like",
            "ankyrin",
            "leucine-rich",
            "amidase",
            "transglutaminase",
            "deamidase",
            "adp-ribosyl",
            "heme acquisition",
            "hasa",
            "two-partner secretion",
            "exported protein",
            "espr",
            "deubiquitinase",
        ),
    ),
    ("hypothetical", ("hypothetical", "uncharacterized", "duf", "domain of unknown", "unknown function")),
]


def classify(text: str) -> str:
    # First match wins, housekeeping before effector ON PURPOSE: a protein whose annotation carries
    # both a housekeeping and an effector token (rare) is scored the conservative way (counts as a
    # likely FP, lowering the ceiling) rather than inflating apparent precision.
    t = text.lower()
    for bucket, kws in RULES:
        if any(k in t for k in kws):
            return bucket
    return "other"


def annot(e: dict) -> str:
    for c in ("detailed_annotation", "gbff_annotation", "product", "broad_annotation"):
        if e.get(c):
            return e[c]
    return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-tag", required=True)
    args = ap.parse_args()
    emissions = read_tsv(BENCH / "data" / "phase2" / f"emissions.{args.run_tag}.tsv")

    rows = []
    for e in emissions:
        rows.append(
            {
                "unit_id": e["unit_id"],
                "locus_tag": e["locus_tag"],
                "substrate_source": e["substrate_source"],
                "nearby_ss_types": e["nearby_ss_types"],
                "is_gold": e["is_gold"],
                "annotation": annot(e),
                "fp_class": classify(annot(e)),
            }
        )

    out = BENCH / "data" / "phase2" / f"emissions_fpclass.{args.run_tag}.tsv"
    write_tsv(out, list(rows[0].keys()), rows)

    order = ["effector", "hypothetical", "apparatus", "housekeeping", "other"]
    print(f"wrote {out.relative_to(BENCH)}  ({len(rows)} emissions)")
    for src in ("proximity", "T5SS-self"):
        sub = [r for r in rows if r["substrate_source"] == src]
        if not sub:
            continue
        c = Counter(r["fp_class"] for r in sub)
        n = len(sub)
        # Obvious-FP counts housekeeping among UNLABELLED rows: a known-gold effector that happens to
        # carry a housekeeping annotation (genome-annotation noise) is a true positive, not an FP.
        hk_fp = sum(r["fp_class"] == "housekeeping" and r["is_gold"] == "no" for r in sub)
        print(f"\n[{src}]  n={n}")
        for b in order:
            print(f"  {b:13s} {c[b]:4d}  {c[b] / n:5.1%}")
        print(
            f"  -> obvious-FP (housekeeping, unlabelled) = {hk_fp}/{n} = {hk_fp / n:.1%}  |  "
            f"soft precision ceiling (non-FP, non-apparatus) = {(n - hk_fp - c['apparatus']) / n:.1%}"
        )
        # gold rows should almost never land in housekeeping; a sanity flag
        gold_hk = [r for r in sub if r["is_gold"] == "yes" and r["fp_class"] == "housekeeping"]
        if gold_hk:
            print(
                f"  !! {len(gold_hk)} GOLD effectors classified housekeeping (rule miss): {[r['locus_tag'] for r in gold_hk]}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
