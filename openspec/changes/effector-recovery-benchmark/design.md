## Context

ssign predicts secretion-system substrates only for proteins within +/-N genes of a detected SS component (default N=3). This is conservative by design and almost certainly misses effectors encoded far from their machinery. No published number says how big that miss is, per SS type, against experimentally-verified effectors. This change builds that number.

A prior attempt (`validation_sweeps/analysis/05_*`-`11_*`) failed two ways: it scored ssign against an incomplete reference (only the best-studied system per genome, so e.g. "T1SS unreachable" was an artifact of zero T1SS effectors in the list), and it used MacSyFinder-derived or hand-typed apparatus positions, which either makes the test circular (MacSyFinder is ssign's core) or introduces transcription errors. That work is discarded.

What already exists and is reused: a literature-curated effector corpus in the `secretion-classifier` repo (`data/verified/`, 949 rows, T1-T6SS) with per-row gene, uniprot, locus_tag, organism, refseq_genome, ss_type, sys_instance_id, evidence_level, primary_ref, and a 2026-06-08 audit assigning each row VERIFIED / PARTIAL / FAIL. The corpus is effectors only; it contains no apparatus genes.

## Goals / Non-Goals

**Goals:**
- A per-SS-type benchmark giving three numbers: ceiling (reachable by proximity), impossible (structurally out of reach), actual (ssign-recovered).
- An effector gold set where every retained row is experimentally validated and citation-verifiable.
- A machinery answer key built independently of ssign, every apparatus gene traceable to a paper quote and a confirmed RefSeq locus_tag.
- A reproducible, documented analysis section under `validation_sweeps/` that states its own strengths and limitations.

**Non-Goals:**
- Fixing the ssign whole-genome-flag substrate-reduction bug (noted, filed separately, not fixed here).
- Using the `predicted`-evidence corpus rows (left untouched for a later task).
- Re-curating external-database entries from their primary literature: we verify each net-new addition (evidence + resolvable citation + locus_tag resolves) but trust the source DB's own curation rather than re-reading every paper.
- A T5SS ceiling/actual number (T5SS is out of the main benchmark; only a preliminary side study).
- Any change to ssign runtime code or dependencies.

## Decisions

**D1. Reuse the curated corpus instead of fresh literature mining.** The 949-row corpus already pairs effectors to systems with citations, and the audit already separated trustworthy biology from broken citations. Rebuilding from zero would discard that. Alternative (fresh mining) rejected: slower, and re-derives work already audited.

**D2. Scope to `validated`, repair recoverable `PARTIAL`, keep only `validated`+`VERIFIED`.** Teo wants experimentally verified effectors. Of 587 validated rows, 352 are already VERIFIED, 223 are PARTIAL (biology sound, citation/uniprot broken, repairable per the audit's documented recipe), 12 are FAIL/NEEDS_REVIEW. Repair re-cites DOIs and rebuilds uniprot from the trustworthy locus_tag; it never re-decides biology. Alternative (VERIFIED-only, skip repair) rejected: would gut T1SS/T4SS/T6SS, which are mostly PARTIAL due to bad DOIs, and bias the benchmark toward T3SS (already 351/355 verified).

**D2b. Reconcile against external curated effector databases.** The corpus reflects what its curators found, so its coverage is a selection bias the benchmark would inherit. Cross-checking against independently-curated databases bounds that bias. Sources: SecReT4 (T4SS, ~1,884 effectors), SecReT6 (T6SS), SecretEPDB and EffectiveDB (manually-curated T3/T4/T6), and BastionHub, which is the only resource with real T1SS (~180) and T2SS (~80) depth. Procedure: pull each, map its identifiers (UniProt/RefSeq/locus_tag) onto our genomes, take net-new effectors, and verify each addition to the same bar as the corpus (experimental evidence + a resolvable citation). Alternative (corpus-only) rejected by the user in favour of the stronger unbiased claim. Constraint: T5SS, T1SS-beyond-BastionHub, and any DB-absent type stay literature-only. Access caveat: BastionHub and SecretEPDB are Monash-hosted and currently returning 503; the durable fallback is the BastionHub NAR-2021 supplementary tables or a Wayback snapshot, so the phase is not blocked on the live site.

**D3. Machinery from literature only (Option B); RefSeq is coordinate lookup, not discovery.** Three candidate sources were tested empirically:
- MacSyFinder: rejected outright, it is the system under test (circular).
- UniProt GO terms: rejected, coverage is near-zero (dotA, a core Dot/Icm gene, has no GO terms; only 1 Legionella protein carries a T4SS GO term against a ~27-gene apparatus).
- NCBI RefSeq `/product` names: they do carry apparatus identity ("type IVB secretion system protein IcmX") and are pipeline-independent, but Teo chose not to let any annotation pipeline even suggest candidate genes. So the literature names the machinery; RefSeq only converts a confirmed gene name to coordinates.

**D4. One curation agent per system instance, not per SS type.** The answer-key unit is the system instance (98 of them, keyed by genome+type+sys_instance_id). One agent per type would own ~32 instances across ~19 genomes, too much context, shallow work, hallucination risk. Per-instance keeps each task bounded and verifiable. Type knowledge is supplied as a short per-type "apparatus brief" handed to each instance agent, so domain expertise without breadth overload.

**D5. Strict four-status anti-hallucination schema.** Each instance agent returns exactly one of COMPLETE / PARTIAL / REFERENCE_ONLY / NONE_KNOWN. Every gene needs a verbatim quote from its naming paper. Agents record gene names, never locus_tags. A separated-paper rule tells agents that the effector's citation often is not the apparatus paper, so follow citations to the founding paper. This makes every claim checkable and gives the agent honest "I couldn't get it" exits instead of forcing fabrication.

**D6. Two verification gates mirror the effector audit.** After extraction: a scripted resolve step maps each paper-named gene to a RefSeq locus_tag + coordinates in the stated genome (flagging non-resolving names), then an independent verification pass re-checks DOI-resolves, paper-names-the-gene, locus_tag-correct. Same shape the original audit used on effectors.

**D7. Distance to the effector's own instance.** Reachability is gene-order distance to the nearest component of the effector's own system instance (matched by sys_instance_id), same replicon. Genomes have multiple same-type systems (PAO1 has three T6SS); counting distance to any same-type component elsewhere would inflate the ceiling. Sweep N in {3,5,7} to show how the ceiling moves with the window.

**D8. T5SS is excluded from the main benchmark and handled as a side study (Option C).** The corpus has no T5SS rows and no curated database tracks T5SS, because its product stays attached to the cell surface rather than being released: T5aSS (monomeric) and T5cSS (trimeric) passengers are fused to their own barrel (same gene), and only T5bSS (two-partner) releases a separate exoprotein (TpsA) encoded adjacent to its TpsB transporter. Measuring a proximity ceiling for a system whose cargo is its own gene is vacuous, so T5SS is left out of the headline ceiling/actual numbers. Separately, ssign flags predicted secreted proteins near T5aSS clusters; that observation is logged as an exploratory side study (chase it in the literature, report preliminary) with no claim of ground truth, since nothing curated exists to benchmark it against. Alternatives rejected: treating T5SS as an in-benchmark "100% by construction" type (misleading, the number carries no information) and full T5SS curation (no source data).

**D9. Bakta-to-RefSeq coordinate bridge for Phase 2.** ssign re-annotates with Bakta and drops original locus_tags, so actual-recall matching uses reciprocal coordinate overlap on the same contig, >=90% of the longer feature. Unbridged proteins are recorded as such, never force-matched.

**D10. Two gated checkpoints.** After Phase 0 the gold set, the answer key, and the proposed genome panel are shown to Teo before any Phase 1 analysis or Phase 2 ssign run. The panel is chosen data-first (genomes with enough verified effectors), not pre-decided.

## Risks / Trade-offs

- **Literature-curation hallucination** → verbatim-quote requirement, four-status honest-exit schema, independent verification pass, and scripted locus_tag resolution that flags rather than guesses.
- **Corpus selection bias** (only well-studied systems are in the corpus) → mitigated by reconciling against external curated DBs (D2b); residual bias (all DBs over-represent model pathogens) stated as a limitation in the docs.
- **External-DB access outage** (BastionHub and SecretEPDB return 503; BastionHub is the sole T1SS/T2SS source) → fall back to the BastionHub NAR-2021 supplementary tables or a Wayback snapshot, and continue chasing the maintainers' reply; T1SS/T2SS depth degrades gracefully to corpus-only if it stays unrecoverable.
- **Cross-DB identifier mapping and dedup** (each DB uses different IDs and may double-count an effector we already have) → map every external entry to a UniProt/locus_tag key, dedup against the corpus on that key, and only verify/keep net-new entries.
- **Thin T1SS after the validated bar** (only 26 validated T1SS rows, spread across ~168 genomes) → T1SS ceiling may rest on few effectors; report the n and caveat it rather than overclaiming.
- **locus_tag format drift** (RefSeq primary `LPG_RS*` vs legacy `lpg*`; corpus uses legacy) → resolve via old_locus_tag/gene alias maps when locating effectors and machinery.
- **Bridge edge cases** (heavy fragmentation, split genes, multi-mapping) → >=90% reciprocal threshold and an explicit unbridged bucket; counts of unbridged reported, not hidden.
- **Plasmid-borne effectors absent from ssign input** (Yersinia's 7 T3SS effectors live on pCD1, which was not in the benchmarked input) → detect replicon-not-in-input and mark as not-in-input, distinct from "missed".
- **Repair does not recover every PARTIAL row** → rows that cannot reach VERIFIED are dropped, with reasons; final gold-set count is an output, not a target.
- **Genomes not yet cached** (only 6 of ~65 genomes have RefSeq downloaded) → a fetch step pulls the rest before resolution.

## Migration Plan

1. Phase 0a: build verified effector gold set (validated scope, repair, drop, finalize).
2. Phase 0b: build machinery answer key (per-instance literature curation, resolve, verify).
3. Checkpoint: present gold set + answer key + proposed genome panel to Teo for approval.
4. Phase 1: ceiling/impossible analysis at N=3,5,7; figures + docs section.
5. Checkpoint: panel approved before any ssign run.
6. Phase 2: ssign runs on the panel; bridge; actual-vs-ceiling; finalize docs.
7. Remove the old `analysis/05_*`-`11_*` scripts and figures `09`/`10`/`11`. Rollback for any phase is `git revert`; nothing here touches ssign runtime, so there is no production rollback surface.

## Open Questions

- **BastionHub recovery:** whether live access returns (maintainers contacted) or we proceed from the NAR-2021 supplement; decided at Phase 0a execution time based on availability.
- **Phase 2 compute:** which CX3 queue/recipe and whether the panel size fits the per-user job cap; resolve when the panel is set.
- **Final gold-set size per SS type:** an output of the repair pass; affects which types have enough effectors to report a stable ceiling.
