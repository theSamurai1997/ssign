# Deep-verify contract (read this in full before judging anything)

You are auditing training-label provenance for a bacterial secretion-system effector dataset. For
each claimed effector, you decide whether its CITED paper actually supports that the named protein is
a secreted substrate/effector of the claimed secretion system in the named organism.

## ANTI-HALLUCINATION CONTRACT (absolute)
- Never fabricate a quote, DOI, gene, or locus tag.
- A verbatim quote must be copied EXACTLY (character for character) from text you actually retrieved.
- Inventing a quote, or "reconstructing" one from memory, is the worst possible failure. When unsure, INACCESSIBLE.

## How to resolve a paper
1. `convert_article_ids` (id_type="doi") -> PMID/PMCID.
2. If a PMCID exists -> `get_full_text_article` (best: full text).
3. Else `get_article_metadata` for the abstract, or WebFetch `https://doi.org/<doi>`.
4. Record the resolved title. If nothing resolves at all -> every effector from it is INACCESSIBLE.

## The three statuses (the REFUTED vs INACCESSIBLE line is the whole point)
- **SUPPORTED** — the retrieved text shows the named protein is a secreted substrate/effector of the
  claimed system in (or consistent with) the named organism. Give a verbatim quote (<=240 chars,
  copied exactly) that names the protein and its secretion/effector role.
- **REFUTED** — you read ENOUGH of the paper to establish POSITIVE counter-evidence. Requires a `reason`:
  - `wrong_organism` — the paper's experimental subject is a different organism and it makes no claim
    about this organism's protein (e.g. a nomenclature/review paper on species A cited as the source
    for a homolog in unrelated species B). This is a real provenance defect — use it.
  - `wrong_system` — the paper is about a different secretion system than the claimed `ss_type`.
  - `wrong_protein` — the paper is about a different protein than the named gene.
  - `no_effector_evidence` — the paper discusses this exact protein but gives no secretion/effector evidence.
- **INACCESSIBLE** — you could NOT retrieve enough text to judge. quote "". This includes:
  paywalled/403 full text where the gene is simply not in the abstract; no PMC and no abstract; DOI
  resolves to a page you cannot read. **Absence of the gene from text you could not fully retrieve is
  INACCESSIBLE, never REFUTED.** REFUTED always needs positive evidence the paper is about something else.

`source` = "cited_paper" when the quote/judgment is from the cited DOI; "external" if you could only
confirm via UniProt/another reference (put that URL/DOI in `note`); "none" otherwise.

## Output
Write the output path you were given as JSON:
`{"batch":N,"results":[{"doi","resolved_title","gene","sys_instance_id","status","reason","source","quote","note"}]}`
Exactly one result object per (paper, effector) in your input batch. `reason` is "" unless status is REFUTED.
Keep notes terse. Final chat message: one line `SUPPORTED=x REFUTED=y INACCESSIBLE=z`. The file is the deliverable.
