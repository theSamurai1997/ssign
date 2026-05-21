#!/usr/bin/env python3
"""Handle T5SS (Type V Secretion System) self-secreting autotransporters.

T5aSS proteins ARE their own substrates (self-secreting). Each MacSyFinder
T5aSS component is then classified by a geometric Pfam-domain filter that
matches the published method (Reid et al., split_t5a_domains.py 2026):

    passenger_length = (barrel_start - LINKER_LENGTH) - (sp_end + 1) + 1

where ``barrel_start`` is the envelope start of the PF03797 (autotransporter
β-barrel) hit on the protein, ``sp_end`` is the SignalP cleavage-site position,
and ``LINKER_LENGTH`` (30 aa) is the alpha-helical linker between passenger
and barrel. The HMM (PF03797 + the porin PF13505 for OMP discrimination) is
bundled with the package at ``src/ssign_app/data/t5aSS/``.

Classification (``domain_group`` column of the domains TSV):
- Classical AT          — PF03797 + passenger >= 100 aa     → KEEP as substrate
- Minimal passenger     — PF03797 + passenger 1-99 aa       → KEEP as substrate
- Barrel-only           — PF03797 with passenger == 0 aa    → DROP (pseudogene fragment)
- OMP/Porin (no AT barrel) — PF13505 only, no PF03797       → DROP (non-AT outer-membrane protein)
- Unclassified-AT       — MacSyFinder called T5aSS but neither HMM hit       → KEEP (lenient)
"""

import argparse
import csv
import logging
import os as _os
import sys as _sys
from importlib.resources import files

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

_scripts_dir = _os.path.dirname(_os.path.abspath(__file__))
if _scripts_dir not in _sys.path:
    _sys.path.insert(0, _scripts_dir)

from ssign_lib.constants import LINKER_LENGTH, MIN_PASSENGER_LENGTH  # noqa: E402

BUNDLED_HMMS = {"PF03797": "PF03797.hmm", "PF13505": "PF13505.hmm"}
DROP_GROUPS = {"Barrel-only", "OMP/Porin (no AT barrel)"}


def _hmm_path(filename: str) -> str:
    """Resolve a bundled HMM path via importlib.resources."""
    return str(files("ssign_app.data.t5aSS") / filename)


def scan_bundled_pfams(proteins_fasta: str) -> dict[str, dict[str, tuple[int, int]]]:
    """Run pyhmmer for each bundled HMM. Returns {locus_tag: {pfam: (env_from, env_to)}}.

    Uses each HMM's gathering threshold (GA cutoff). Both PF03797 and PF13505
    have GA set in their headers (26 and 28.7 bits respectively).

    FRAGILE: pyhmmer.hmmsearch API. If pyhmmer changes the call signature or
    domain-attribute names, this breaks. Tested with pyhmmer 0.10.x. If this
    breaks: check src/ssign_app/shims/hmmsearch.py for the equivalent call
    pattern — it uses the same API and will fail in the same way.
    """
    import pyhmmer
    from pyhmmer.easel import Alphabet, SequenceFile
    from pyhmmer.plan7 import HMMFile

    alphabet = Alphabet.amino()
    with SequenceFile(proteins_fasta, digital=True, alphabet=alphabet) as sf:
        targets = sf.read_block()

    hits: dict[str, dict[str, tuple[int, int]]] = {}
    for pfam_id, hmm_filename in BUNDLED_HMMS.items():
        with HMMFile(_hmm_path(hmm_filename)) as hf:
            hmm = next(iter(hf))
        top_hits = next(pyhmmer.hmmsearch(hmm, targets, bit_cutoffs="gathering"))
        for hit in top_hits:
            name = hit.name.decode()
            if not hit.domains:
                continue
            best = min(hit.domains, key=lambda d: d.i_evalue)
            hits.setdefault(name, {})[pfam_id] = (best.env_from, best.env_to)
    return hits


def classify_t5a(
    pfam_hits: dict[str, tuple[int, int]],
    sp_end: int | None,
) -> tuple[str, int]:
    """Apply geometric classifier. Returns (domain_group, passenger_length).

    When SignalP gave no cleavage site, sp_end falls back to 1 so a real
    autotransporter is kept rather than mis-dropped on a SignalP miss.
    """
    has_barrel = "PF03797" in pfam_hits
    has_porin = "PF13505" in pfam_hits

    if not has_barrel and has_porin:
        return "OMP/Porin (no AT barrel)", 0
    if not has_barrel:
        return "Unclassified-AT", 0

    barrel_start, _barrel_end = pfam_hits["PF03797"]
    effective_sp_end = sp_end if (sp_end is not None and sp_end > 0) else 1
    passenger_length = max(0, (barrel_start - LINKER_LENGTH) - (effective_sp_end + 1) + 1)

    if passenger_length >= MIN_PASSENGER_LENGTH:
        return "Classical AT", passenger_length
    if passenger_length >= 1:
        return "Minimal passenger", passenger_length
    return "Barrel-only", 0


def _parse_sp_end(raw: str) -> int | None:
    """Parse SignalP CS position (formats vary: '22', '22-23', '22.0', '')."""
    if not raw:
        return None
    token = raw.split("-")[0].strip()
    try:
        return int(float(token))
    except (ValueError, TypeError):
        return None


def main():
    parser = argparse.ArgumentParser(description="Handle T5SS autotransporters")
    parser.add_argument("--ss-components", required=True)
    parser.add_argument("--predictions", required=True)
    parser.add_argument("--proteins", required=True, help="Proteins FASTA for Pfam scan")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out-substrates", required=True)
    parser.add_argument("--out-domains", required=True)
    args = parser.parse_args()

    t5ss_components = []
    with open(args.ss_components) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("ss_type", "").startswith("T5"):
                t5ss_components.append(row)

    predictions = {}
    with open(args.predictions) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            predictions[row["locus_tag"]] = row

    pfam_hits = scan_bundled_pfams(args.proteins) if t5ss_components else {}

    substrates = []
    classifications = []
    for comp in t5ss_components:
        locus = comp["locus_tag"]
        ss_type = comp.get("ss_type", "T5aSS")
        pred = predictions.get(locus, {})
        sp_end = _parse_sp_end(pred.get("signalp_cs_position", ""))

        # Geometric filter only applies to T5aSS (T5bSS/T5cSS use their own
        # MacSyFinder profiles and don't fragment the same way).
        if ss_type == "T5aSS":
            domain_group, passenger_length = classify_t5a(pfam_hits.get(locus, {}), sp_end)
        else:
            domain_group, passenger_length = f"{ss_type}-component", 0

        classifications.append(
            {
                "locus_tag": locus,
                "sample_id": args.sample,
                "ss_type": ss_type,
                "domain_group": domain_group,
                "passenger_length": passenger_length,
                "sp_end": sp_end if sp_end is not None else "",
                "barrel_start": pfam_hits.get(locus, {}).get("PF03797", ("", ""))[0],
            }
        )

        if domain_group in DROP_GROUPS:
            continue

        try:
            dlp_prob = float(pred.get("dlp_extracellular_prob", pred.get("extracellular_prob", 0)))
        except (ValueError, TypeError):
            dlp_prob = 0.0

        substrates.append(
            {
                "locus_tag": locus,
                "sample_id": args.sample,
                "tool": "T5SS-self",
                "nearby_ss_types": ss_type,
                "dlp_extracellular_prob": dlp_prob,
                "predicted_localization": pred.get("predicted_localization", ""),
                "dlp_max_localization": pred.get("dlp_max_localization", ""),
                "dlp_max_probability": pred.get("dlp_max_probability", ""),
                "dse_ss_type": pred.get("dse_ss_type", ""),
                "dse_max_prob": pred.get("dse_max_prob", ""),
                "signalp_prediction": pred.get("signalp_prediction", ""),
                "signalp_probability": pred.get("signalp_probability", ""),
                "signalp_cs_position": pred.get("signalp_cs_position", ""),
                "product": pred.get("product", ""),
            }
        )

    sub_fields = [
        "locus_tag",
        "sample_id",
        "tool",
        "nearby_ss_types",
        "dlp_extracellular_prob",
        "predicted_localization",
        "dlp_max_localization",
        "dlp_max_probability",
        "dse_ss_type",
        "dse_max_prob",
        "signalp_prediction",
        "signalp_probability",
        "signalp_cs_position",
        "product",
    ]
    with open(args.out_substrates, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=sub_fields, delimiter="\t")
        writer.writeheader()
        for s in substrates:
            writer.writerow(s)

    domain_fields = [
        "locus_tag",
        "sample_id",
        "ss_type",
        "domain_group",
        "passenger_length",
        "sp_end",
        "barrel_start",
    ]
    with open(args.out_domains, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=domain_fields, delimiter="\t")
        writer.writeheader()
        for c in classifications:
            writer.writerow(c)

    n_dropped = sum(1 for c in classifications if c["domain_group"] in DROP_GROUPS)
    logger.info(
        "Found %d T5SS self-substrates in %s (dropped %d as barrel-only/OMP)",
        len(substrates),
        args.sample,
        n_dropped,
    )


if __name__ == "__main__":
    main()
