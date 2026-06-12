#!/usr/bin/env python3
"""Dataset group 6, task 6.1: deterministic citation-consistency sweep (ssign-found scope).

The benchmark's 19-row discordant audit found that effector-sourcing DOIs in the answer key
are frequently wrong (resolve to an unrelated paper) or 404, even on rows marked VERIFIED. That
is a training-label provenance defect, not a recall problem. This script measures it
deterministically (no agents) over the ssign-FOUND gold effectors so the result is reproducible.

For each of the 51 `emitted_secreted` rows in the Phase-2 actual-call table:
  1. join to positives_all.tsv for the sourcing DOI (primary_ref, else instance_source_doi) +
     organism (the actual-call table carries neither);
  2. confirm the DOI resolves (DOI.org Handle API, reusing 03_finalize_gold_set.doi_resolves);
  3. fetch CrossRef title + abstract + container-title (cached) and ask the only question CrossRef
     can answer deterministically: does the row's own gene name or its genus appear in the paper?

Verdict per row:
  CONSISTENT       - gene or genus appears in the resolved paper's title/abstract/journal.
  FLAG_WRONG_TOPIC - abstract present, neither gene/genus AND no SS-type phrase appears -> the DOI
                     points at a paper in an unrelated field (the strong wrong-DOI signal:
                     Toxoplasma invasion for an E. coli effector, soil ecology for a Legionella one).
  FLAG_GENE_ABSENT - abstract present, the SS-type topic IS there but neither gene nor genus is named
                     -> a plausibly-correct machinery paper that just doesn't name this protein in its
                     abstract; weaker, still queue for pass-2 to confirm.
  DOI_UNRESOLVED   - DOI.org does not resolve the handle (404 / never registered), or the row had no DOI.
  FETCH_ERROR      - DOI resolves but CrossRef was unreachable after retries (e.g. 429); re-run to settle
                     (kept distinct so a throttled run is never silently miscounted as INDETERMINATE).
  INDETERMINATE    - resolves but CrossRef gives no abstract (publisher never deposited one) or the
                     DOI is not in CrossRef at all; a title-only check can't refute -> queue for the
                     pass-2 agent re-audit (task 6.2), do NOT call it a mismatch.

The 19 rows already hand-audited in the benchmark (discordant_audit.md) all fall inside this set,
so their manual verdict is merged in as `prior_audit` to measure deterministic-vs-manual agreement.

Inputs : data/phase2/actual_per_effector.panel_genbank_t3ss.tsv  (the 51 found effectors)
         data/dataset/positives_all.tsv                          (sourcing DOI + organism)
Outputs: data/dataset/citation_consistency_found.tsv             (one row per found effector)
         data/dataset/.doi_cache.json        (DOI.org resolution cache; shared with 31)
         data/dataset/.crossref_cache.json   (CrossRef content cache; re-runs are offline)
Run:     .venv/bin/python scripts/41_citation_consistency.py
"""

from __future__ import annotations

import html
import importlib.util
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "data" / "dataset"
PHASE2 = ROOT / "data" / "phase2"
CEILING = ROOT / "data" / "phase1" / "ceiling_per_effector.tsv"
DOI_CACHE = DATASET / ".doi_cache.json"
CROSSREF_CACHE = DATASET / ".crossref_cache.json"
# Distinct from 03's UA on purpose: the `mailto:` token opts CrossRef's polite pool in (03 only hits
# DOI.org, which has no polite pool). Don't "unify" the two into one string.
UA = "ssign-benchmark-citation-check/1.0 (mailto:teoreid@gmail.com)"

# Reuse the benchmark's DOI.org resolver verbatim (single source of truth for resolution).
_spec = importlib.util.spec_from_file_location("finalize_gold_set", Path(__file__).parent / "03_finalize_gold_set.py")
_fz = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_fz)
import bench_index as bi  # noqa: E402  (scripts/ on sys.path)
from bench_io import by_type, norm_doi, read_tsv, write_tsv  # noqa: E402  (scripts/ on sys.path)


def _real_uni(u: str) -> str:
    """A usable UniProt accession, treating the literal '-' placeholder (common in the found set)
    as absent. Joining on '-' silently collapses every placeholder row onto one positives entry."""
    u = (u or "").strip()
    return u if u and u != "-" else ""


def resolve_positive(found_row: dict, cbyloc: dict, by_inst: dict, by_uni: dict, by_gsg: dict):
    """Bridge a found effector to its positives_all row (-> sourcing DOI + organism).

    Keyed via the gold-set system instance, NOT raw UniProt: ~1/3 of found effectors carry
    uniprot='-', so a UniProt-first join collapses them. found -> ceiling by (effector_locus, gene)
    yields the instance_id, which keys positives uniquely. UniProt / (gene,ss,genome) are fallbacks.
    Returns (positives_row | {}, join_method)."""
    c = cbyloc.get((found_row["effector_locus"], found_row["gene"]), {})
    inst = (c.get("instance_id") or "").strip()
    ss = found_row["ss_type"]
    g = found_row["gene"].lower()
    uni = _real_uni(found_row.get("uniprot")) or _real_uni(c.get("uniprot"))
    gb = bi.accession_base((c.get("refseq_genome") or "").strip()) if c.get("refseq_genome") else ""
    if inst and (inst, g, ss) in by_inst:
        return by_inst[(inst, g, ss)], "instance"
    if uni and uni in by_uni:
        return by_uni[uni], "uniprot"
    if gb and (g, ss, gb) in by_gsg:
        return by_gsg[(g, ss, gb)], "gene_ss_genome"
    return {}, "UNRESOLVED"


# SS-type -> phrases a paper would use for it; supporting signal only (gene/genus are decisive).
SS_PHRASES = {
    "T1SS": ["t1ss", "type i secretion", "type 1 secretion", "abc transporter", "atp-binding cassette"],
    "T2SS": ["t2ss", "type ii secretion", "type 2 secretion"],
    "T3SS": ["t3ss", "type iii secretion", "type 3 secretion", "injectisome"],
    "T4SS": ["t4ss", "type iv secretion", "type 4 secretion"],
    "T5SS": ["t5ss", "type v secretion", "autotransporter", "two-partner secretion"],
    "T6SS": ["t6ss", "type vi secretion", "type 6 secretion"],
}

# Mutually-exclusive partition of the 19 hand-audited rows (benchmark discordant_audit.md / fig04).
PRIOR_AUDIT = {
    **{g: "sound" for g in ["TseM", "TseZ", "Tlde1A", "BipB", "BipC"]},
    **{g: "wrong_or_404_doi" for g in ["celA", "plaA", "VirA", "CopN", "Tle1", "Tae4_Stm"]},
    **{g: "unidentifiable_row" for g in ["EFF00142", "EFF00150", "TseA_T6SS1", "ChlaDub1"]},
    **{g: "misassigned_ss_type" for g in ["BopA", "BopE"]},
    **{g: "duplicate_row" for g in ["Tle4", "TplE_alias_Tle4"]},
}


def crossref_content(doi: str, cache: dict) -> dict | None:
    """CrossRef /works/{doi} -> {title, abstract, container, in_crossref}; None on hard failure.
    in_crossref=False means a 404 from CrossRef (DOI may still resolve via a non-CrossRef agency)."""
    key = norm_doi(doi)
    if not key:
        return None
    if key in cache:
        return cache[key]
    url = "https://api.crossref.org/works/" + urllib.parse.quote(key, safe="")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    out: dict | None = None
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                msg = json.load(resp).get("message", {})
            out = {
                "in_crossref": True,
                "title": " ".join(msg.get("title", []) or []),
                "abstract": msg.get("abstract", "") or "",
                "container": " ".join(msg.get("container-title", []) or []),
            }
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                out = {"in_crossref": False, "title": "", "abstract": "", "container": ""}
                break
            time.sleep(1 + attempt)
        except Exception:
            time.sleep(1 + attempt)
    if out is not None:
        cache[key] = out
    return out


def _token_hit(token: str, haystack: str) -> bool:
    """Word-boundary, case-insensitive match so short gene names ('celA', 'Tle1') don't match
    inside longer words. JATS abstract markup is stripped to plain text by the caller."""
    token = token.strip()
    if len(token) < 3:  # 1-2 char tokens are too noisy to assert on
        return False
    return re.search(r"(?<![A-Za-z0-9])" + re.escape(token.lower()) + r"(?![A-Za-z0-9])", haystack) is not None


def crossref_haystack(content: dict | None) -> str:
    """Lowercased, tag-and-entity-stripped title+abstract+container for token matching. Shared by the
    pass-1 check and the pass-2 verifier so both strip JATS markup AND HTML entities identically."""
    if not content:
        return ""
    raw = " ".join([content["title"], content["abstract"], content["container"]])
    return html.unescape(re.sub(r"<[^>]+>", " ", raw)).lower()


def effector_hit(gene: str, genus: str, ss: str, haystack: str) -> bool:
    """True if the effector's gene, its genus, or one of its SS-type phrases appears in the paper text."""
    return bool(
        (gene and _token_hit(gene, haystack))
        or (genus and _token_hit(genus, haystack))
        or any(ph in haystack for ph in SS_PHRASES.get(ss, []))
    )


def main() -> int:
    found = [
        r
        for r in read_tsv(PHASE2 / "actual_per_effector.panel_genbank_t3ss.tsv")
        if r["ssign_call"] == "emitted_secreted"
    ]
    pos = read_tsv(DATASET / "positives_all.tsv")
    cbyloc = {(c["effector_locus"].strip(), c["gene"].strip()): c for c in read_tsv(CEILING)}
    by_inst, by_uni, by_gsg = {}, {}, {}
    for p in pos:
        g = (p.get("gene") or "").strip().lower()
        ss = (p.get("ss_type") or "").strip()
        inst = (p.get("sys_instance_id") or "").strip()
        u = _real_uni(p.get("uniprot"))
        gb = bi.accession_base((p.get("refseq_genome") or "").strip()) if p.get("refseq_genome") else ""
        if inst:
            by_inst.setdefault((inst, g, ss), p)
        if u:
            by_uni.setdefault(u, p)
        if gb:
            by_gsg.setdefault((g, ss, gb), p)

    doi_cache: dict = json.loads(DOI_CACHE.read_text()) if DOI_CACHE.exists() else {}
    cr_cache: dict = json.loads(CROSSREF_CACHE.read_text()) if CROSSREF_CACHE.exists() else {}

    rows: list[dict] = []
    for r in found:
        gene = (r.get("gene") or "").strip()
        uni = _real_uni(r.get("uniprot"))
        ss = r["ss_type"]
        p, join_method = resolve_positive(r, cbyloc, by_inst, by_uni, by_gsg)
        organism = (p.get("organism") or "").strip()
        genus = organism.split()[0] if organism else ""
        doi = norm_doi(p.get("primary_ref", "")) or norm_doi(p.get("instance_source_doi", ""))

        # Only pause for actual network calls; a fully-cached re-run shouldn't sleep per row.
        needs_net = bool(doi) and (doi not in doi_cache or doi not in cr_cache)
        resolves = _fz.doi_resolves(doi, doi_cache) if doi else False
        content = crossref_content(doi, cr_cache) if doi else None
        if needs_net:
            time.sleep(0.05)

        in_crossref = bool(content and content["in_crossref"])
        has_abstract = bool(content and content["abstract"])
        haystack = crossref_haystack(content)
        gene_hit = bool(gene) and _token_hit(gene, haystack)
        genus_hit = bool(genus) and _token_hit(genus, haystack)
        ss_hit = any(ph in haystack for ph in SS_PHRASES.get(ss, []))

        # doi_resolves() is False for an empty DOI too, so this one arm covers "no DOI on the row".
        if not resolves:
            verdict = "DOI_UNRESOLVED"
        elif doi and content is None:
            verdict = "FETCH_ERROR"  # CrossRef unreachable after retries (e.g. 429) -> re-run, not a real call
        elif gene_hit or genus_hit:
            verdict = "CONSISTENT"
        elif not in_crossref or not has_abstract:
            verdict = "INDETERMINATE"  # title-only -> can't refute deterministically; pass-2 agent
        elif ss_hit:
            verdict = "FLAG_GENE_ABSENT"  # right SS topic, gene/genus not named -> imprecise, check
        else:
            verdict = "FLAG_WRONG_TOPIC"  # abstract present, unrelated field -> wrong-paper DOI

        rows.append(
            {
                "gene": gene,
                "uniprot": uni,
                "ss_type": ss,
                "organism": organism,
                "unit_id": r["unit_id"],
                "join_method": join_method,
                "pos_sys_instance_id": p.get("sys_instance_id", ""),
                "pos_uniprot": p.get("uniprot", ""),
                "pos_locus_tag": p.get("locus_tag", ""),
                "sourcing_doi": doi,
                "doi_resolves": "yes" if resolves else "no",
                "in_crossref": "yes" if in_crossref else "no",
                "has_abstract": "yes" if has_abstract else "no",
                "gene_in_paper": "yes" if gene_hit else "no",
                "genus_in_paper": "yes" if genus_hit else "no",
                "ss_phrase_in_paper": "yes" if ss_hit else "no",
                "verdict": verdict,
                "prior_audit": PRIOR_AUDIT.get(gene, ""),
                "crossref_title": content["title"] if content else "",
            }
        )

    DOI_CACHE.write_text(json.dumps(doi_cache, indent=0))
    CROSSREF_CACHE.write_text(json.dumps(cr_cache, indent=0))
    if not rows:
        print("no emitted_secreted rows found; nothing to check")
        return 0
    cols = list(rows[0].keys())
    write_tsv(DATASET / "citation_consistency_found.tsv", cols, rows)

    verdicts = Counter(r["verdict"] for r in rows)
    flag_types = ("FLAG_WRONG_TOPIC", "FLAG_GENE_ABSENT", "DOI_UNRESOLVED", "FETCH_ERROR")
    flagged = sorted((r for r in rows if r["verdict"] in flag_types), key=lambda r: r["verdict"])
    print(f"wrote data/dataset/citation_consistency_found.tsv  ({len(rows)} ssign-found effectors)")
    print("verdicts:", dict(verdicts), f"  ({by_type(rows)})")
    print(f"\nflagged for pass-2 ({len(flagged)}):")
    for r in flagged:
        print(f"  {r['verdict']:16s} {r['ss_type']:5s} {r['gene']:18s} {r['sourcing_doi']}")
        print(f"                   title: {r['crossref_title'][:88]}")

    # Deterministic-vs-manual agreement on the 19 already-audited rows.
    audited = [r for r in rows if r["prior_audit"]]
    print(f"\ndeterministic vs prior manual audit ({len(audited)} overlapping rows):")
    for r in audited:
        print(f"  {r['gene']:18s} prior={r['prior_audit']:20s} -> deterministic={r['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
