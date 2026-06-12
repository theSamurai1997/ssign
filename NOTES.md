# Deferred work

Tracks items skipped during tasks. One bullet per item: what, why, trigger to revisit.

## effector-recovery-benchmark — CITATION-INTEGRITY: NEXT TASK (2026-06-12)

**CONFIRMED real** (not an agent hallucination): 3/3 spot-checks via CrossRef API (deterministic)
match the audit agents exactly — celA/plaA DOI = a soil-16S-ecology paper, VirA DOI = a mouse
tooth-development paper, Tle4 DOI = a protein-design paper. Genes are REAL (real loci/genomes);
it's the sourcing-DOI metadata that's wrong, + 2 genuine SS-type mislabels (BopA/BopE = Bsa T3SS
not T6SS). Our own `doi_resolves` col is blank for all 19 and `verification_status`=VERIFIED, so
"VERIFIED" never meant "DOI points to the right paper."

**CAVEAT — don't extrapolate:** the 19 audited are the most-suspicious slice (emitted-but-machinery-
far), NOT a random sample. 14/19 defective there says nothing reliable about the other ~480.
Prevalence across the full set is unknown → must be measured.

**Precision blind spot:** ssign emitted 1,933 (default) / 2,321 (t3ss) secreted proteins panel-wide;
only 39/51 are gold effectors → **1,894 / 2,270 emissions are unvalidated** (novel TP or FP, no
ground truth). Benchmark measures recall only, says nothing about precision.

**AGREED NEXT STEP (Teo, 2026-06-12):** deterministic citation-consistency sweep (CrossRef
title+abstract vs the row's gene/organism/SS-type; flag mismatches; NO agents in pass 1), then
targeted re-audit of flagged rows. **SCOPE: only the ssign-FOUND effectors for now** (the ~51
emitted_secreted gold effectors in `actual_per_effector.panel_genbank_t3ss.tsv`), to save time —
not the full 582. Items 3 (triage wrong-DOI-vs-mislabel by training impact) + 4 (precision estimate
for the 1.9k unlabeled emissions) deferred to discussion when we get there.

**PASS-1 DONE (2026-06-12)** — secretion-classifier-dataset task 6.1, `scripts/41_citation_consistency.py`
→ `data/dataset/citation_consistency_found.tsv`. **A join bug was caught and fixed mid-task** (the first
reported numbers were partly wrong): joining found→positives on raw UniProt collapsed the 17 found
effectors with `uniprot='-'` (a truthy dash) onto one arbitrary positives row (12 shared DOI nature11433)
→ garbage organism/DOI for those rows. Fixed by bridging found→`ceiling_per_effector` by (effector_locus,
gene) → gold `instance_id`, then keying positives by (instance_id, gene, ss). All 51 now resolve via
instance, 35 distinct DOIs (`join_method` column records how).

**Corrected result (51 found):** **21 CONSISTENT** (14 by gene name, 7 genus-only = weaker),
**8 FLAG_WRONG_TOPIC** (DOI in an unrelated field — celA/plaA→soil ecology, Map×2→Toxoplasma,
VirA→tooth development, BipB→colicins, aprA→CS1 pilin, and **EspZ→the T6SS-discovery PNAS paper**),
**3 FLAG_GENE_ABSENT** (prtB/prtA1/aprA cite a real Prt/Lip ABC-exporter T1SS paper that just doesn't
name the protease in its abstract — probably OK, confirm in pass-2), **3 DOI_UNRESOLVED** (YspE/CopN/
ChlaDub1 — classic ASM/JBC/Science DOIs failing the DOI.org handle check; likely legacy false-negatives
or typos), **16 INDETERMINATE** (CrossRef record exists but no abstract → can't refute deterministically).
Net: the deterministic method adjudicates ~30/51; 16 indeterminate + 14 flagged = 30 go to pass-2.
Validated vs the 19 prior manual verdicts: re-flagged celA/plaA/VirA independently AND caught **BipB's
wrong DOI the manual audit missed** (it was "sound" only because its agent was AUP-blocked; the DOI is a
colicin paper). Also: **aprA has two instances with two different DOIs** (one wrong CS1-pilin, one OK-ish
lipase) → the defect is per-row metadata, not per-protein. Caveat: ssign-found ~10% slice; full-582
prevalence still unmeasured.

**PASS-2 DONE (task 6.2)** — `pass2_results.tsv` + `pass2_results.md`; raw returns `pass2_raw/batch_*.json`;
verify `scripts/42_pass2_verify.py`. 30 non-CONSISTENT rows re-audited: 5 literature agents +
batch C (YspE/VirA/BipB/CopN/ChlaDub1) done in-session via PubMed (agent AUP-blocked on B. pseudomallei).
Every returned DOI deterministically re-verified (registered + on CrossRef + gene/genus present):
**26 rows now carry a verified-real on-topic primary DOI** (21 RESOLVED corrected, 3 CONFIRMED,
2 MISASSIGNED), **0 fabricated/unverified DOIs** — the agents agreed with the deterministic check on all
26, so they were reliable here.

**Dataset edits this produced (for the training set, NOT yet applied):**
- **2 MISASSIGNED ss_type:** BopA + BopE (B. thailandensis BTH_II rows) are **Bsa T3SS** effectors, not
  T6SS (name-collision; BopA=Cullinane 2008, BopE=Stevens 2003). Fix ss_type or drop. These are 2 of the
  benchmark's 6 "accidental cross-type" emissions.
- **4 NOT_FOUND (drop/down-tier):** ChlaDub1 (no T3SS-secretion evidence, DUB-only per Misaghi 2006);
  P.entomophila aprA (homology-only, no ortholog-specific T1SS paper); EFF00136 + EFF00150 (BTH DUF3274
  proteins not named in their cited paper).
- **1 duplicate:** Tle4 (idx15) = TplE_alias_Tle4 (idx24) = PA1510; BOTH also had wrong DOIs (protein-
  design / lung-regeneration papers) → correct = Jiang 2016 Cell Rep 10.1016/j.celrep.2016.07.012. De-dupe.
- EFF00142's pass-1 "unidentifiable" tag was WRONG: Russell 2012 names BTH_I2691 as a T6SS-1 substrate → CONFIRMED.
- EspZ's DOI was the T6SS-discovery PNAS paper → corrected to Kanack 2005 (10.1128/IAI.73.7.4327-4337.2005).

**ITEM 3 DONE (2026-06-12, task 6.3)** — applied to `positives_all.tsv` via `scripts/43_apply_citation_corrections.py`
(post-build overlay, **930→925 rows**). Triaged by training impact: **22 pure-provenance DOI fixes** (label
unchanged), **2 label-changing** (BopA/BopE T6SS→T3SS, re-labelled to type-level T3SS + DOI fixed), **5 removed**
(4 unsupported + TplE_alias_Tle4 duplicate → `positives_removed_citation.tsv`). Backup =
`positives_all.pre_citation_audit.tsv`; per-field log = `citation_corrections_log.tsv`. Verified surgical
(ss_type counts moved exactly T1−1/T3+1/T6−5; other instances of EspZ/Tle4 untouched), idempotent, and
byte-identical on re-run from the backup.

**Load-bearing caveats:**
- **Overlay, not in the build chain:** re-running `37_fold_t5ss.py` regenerates the uncorrected 930-row table,
  so `43` MUST be re-applied after any `37` rebuild. Documented in `data/dataset/README.md` "Citation-audit overlay".
- **Downstream re-run:** `40_pair_features.py` (and the eventual feature matrix) must be rebuilt — 2 ss_type
  changes + 5 dropped rows. pair_features.tsv is now stale.
- **Out of scope (real, deferred):** the audit only touched the ssign-FOUND rows. Other non-found instances may
  carry the SAME wrong DOIs — confirmed: a second EspZ instance (T3SS_24) still cites the T6SS-discovery PNAS
  paper. Fixing those is the full-corpus citation sweep (the original "prevalence unknown" item), still deferred.

**ITEM 4 DONE (2026-06-12, effector-recovery-benchmark §8)** — deterministic precision estimate, scripts
28→29→30→31, `data/phase2/PRECISION.md` + `figures/precision/01..03`. Teo chose deterministic tiers only
(DB-confirmed floor + obvious-FP), agent-sampled adjudication deferred. Denominator = the 1,572 proximity
substrate calls (default); 361 T5SS-self assessed separately (correct-by-construction; 68% effector /
0.3% housekeeping → T5SS detection sound). **Result: proximity precision is a wide band, ~3% provable
floor (DB-confirmed homology to SecReT4/6; T6SS the only well-covered type at 8.3%) → ~75% soft ceiling
(not-obviously-non-secreted), with ~19% clearly-FP housekeeping, ~6% machinery, ~11% annotation-effector,
and ~64% (hypothetical+other) unresolvable by DB or annotation.** Pairs with recall (8–10% emitted of
testable): the proximity rule is both permissive AND low-recall → the case for the classifier, which can
adjudicate the unresolvable middle. Tier-1 used pyhmmer phmmer (no external aligner; the repo venv is at
`../../.venv`, NOT `benchmark/.venv` — stdlib scripts silently fell back to system python3).

**Caveats baked into PRECISION.md:** no negative ground truth (ceiling overstates); annotation-based tier
inherits genome-annotation errors (NleB mis-annotated "IS3 transposase"); DB floor covers only T4/T6SS.

**STILL OPEN:** the full-corpus citation sweep (point 41→42→43 at the unaudited rows; defect confirmed
beyond the found slice), and optionally a stratified agent-adjudicated precision sample to place a point
estimate inside the 3–75% band. Both deferred to discussion with Teo.

### Earlier deferred repair items (from the 19-row audit; full table in `data/phase2/discordant_audit.md`):
- **Wrong/nonexistent sourcing DOI (≥9):** celA, plaA, VirA, ChlaDub1(404), CopN(404), Tle4, TplE, Tle1, Tae4_Stm. Re-source before trusting.
- **Genuine misassignment:** BopA, BopE labelled T6SS but are Bsa **T3SS** effectors (Bop name-collision) → reassign. ChlaDub1 T3SS call unsupported (inclusion-membrane DUB) → review.
- **Unidentifiable/unsupported rows:** EFF00142, EFF00150 (opaque DB placeholders), TseA_T6SS1 (cited paper never names it) → drop or repair.
- **Duplicate row:** Tle4 == TplE_alias_Tle4 (PA1510/Q9I3K2) → de-dupe.
- Same defect class likely across the wider corpus → run citation audit on ALL positives, not just these 19.
- BipB/BipC agents blocked by AUP filter on B. pseudomallei content; checked manually (sound Bsa T3SS translocon). If a future agent audit is needed, run those two by hand.

## effector-recovery-benchmark — FULL PANEL DONE (2026-06-12)

12 jobs (2988054-2988065) finished clean on CX3 gpu72: both tags 67/67 genomes, 0 empty/failed, 196,170 protein rows each. Synced to `data/phase2/runs/panel_genbank_{default,t3ss}/`; scored with scripts 24 + 25. Outputs `data/phase2/actual_{per_effector,vs_ceiling}.panel_genbank_{default,t3ss}.tsv`.

**Exact ssign-detected systems behind the found effectors** (`scripts/26_found_systems.py`, `found_systems.*.tsv`): default 39 found effectors via **32** distinct detected systems; t3ss 51 via **39**. Emission basis (why emitted): default own-type(legit) **33/39**, cross-type-only(accidental) **6/39**; t3ss own-type **45/51**, cross-type-only **6/51**. So the large majority of found effectors sit next to their OWN system type; only **6** are accidental cross-type emissions (the celA/plaA-via-T4aP, BopA/BopE-via-T5cSS, VirA/ChlaDub1-via-T5aSS set, same as the audit's cross-type rows). NB: required a family-normalization fix in script 26 (ssign labels detected systems by TXSScan subtype `T6SSi`, not the coarse `T6SS`); an earlier uncorrected run wrongly showed ~half accidental.

**Full-panel recall (emitted of testable, testable=499):**

| SS | testable | ceiling@7 | actual default | actual t3ss |
|----|---|---|---|---|
| T1SS | 25 | 80% | 60% | 60% |
| T2SS | 77 | 4% | 5% | 5% |
| T3SS | 227 | 33% | 1% | 7% |
| T4SS | 98 | 10% | 0% | 0% |
| T6SS | 72 | 26% | 24% | 24% |
| ALL | 499 | 25% | **8%** | **10%** |

Headline: even within +/-7 genes (ceiling 25%), as-shipped ssign emits 8% (T3SS-excluded) / 10% (T3SS-included). The ceiling-to-actual gap (25% -> 8-10%) is the core motivation for the secretion-classifier. T4SS actual=0% is the starkest gap (ceiling 10%, emits nothing). T3SS default=1% is by-design (excluded); the t3ss pass lifts it to 7%, still far under the 33% ceiling. 16-19 effectors emitted but ceiling-unreachable@7 = ssign's detected machinery disagrees with the literature answer key (mostly T6SSi self-adjacency, a few T5aSS/T4aP); benign, listed in script-25 output.

Denominator honesty: 14 no_run (effectors with no genome accession, all non-testable) + 75 not_in_input (only 6 testable; ssign's ORF set didn't carry those 6 loci = forced misses already counted against recall). Nothing inflated.

**Next (now unblocked):** 6.6/6.7 figures + docs; dataset group-4 4.1 feature join (the `actual_per_effector.*.tsv` tables ARE the per-protein tool signals). Benchmark-side bridge/SUBMIT/script changes still UNCOMMITTED (not yet approved).

## effector-recovery-benchmark — Phase 2 pilot results + FASTA bridge bug (2026-06-11)

4-genome pilot (NC_002516.2/003197.2/004337.2/004578.1) ran on CX3, 3 tags x 4 genomes, all 20/20 steps. Synced to `data/phase2/runs/`; scored with scripts 24 + 25. Outputs: `data/phase2/actual_{per_effector,vs_ceiling}.<tag>.tsv`.

- **In-panel pilot numbers (genbank_default):** of 118 testable effectors in the 4 genomes, ceiling-reachable@7 = 30 (25%), ssign emitted 4 = **3% of testable / 13% of proximity-reachable@7**. T3SS-enabled tag: 6 emitted = 5% / 20%. The gap is large and real (supports the classifier; doc 00 decision tree "gap >20% -> proceed"). NB script 25's printed table dilutes "actual" across all 499 effectors incl. 445 no_run (not-in-this-4-genome-pilot) -> shows 1%; the in-panel 3-5% is the honest pilot number. **Full panel run needed for real figures (6.6/6.7).**
- **FASTA-mode bridge — FIXED 2026-06-11 (coordinate bridge).** Was 0/137 matched: Bakta renames locus_tags (`NFOBNJ_00001`) and the seq fallback read a `sequence` column ssign's results_raw never emits. Fix: `bench_runout.RunOutput.find_by_coord` + `by_coord` index match on (contig_base, strand, 3'-stop) since Bakta calls the same ORF at the same RefSeq coordinate; script 24 looks up effector coords from the gene-order index. Result: **0→116 matched, and the emitted set is identical to genbank mode** (4 effectors, exact), confirming correctness. Remaining 21 unmatched are effectors with no locus_tag at all (unbridgeable). If ssign ever emits an aa `sequence` column, the seq bridge reactivates automatically (kept as a later tier).
- **3 emitted-but-ceiling-unreachable@7** (VirA/Q7BU69 T3SS nearby=T5aSS; Tle4 + TplE_alias_Tle4 T6SS nearby=T6SSi): ssign emitted these off a *different* nearby system than the literature answer-key assigns. This is the benchmark task-6b "ssign machinery != literature instance" signal; capture in the 6b write-up, not a bug.

## CX3 environment

- **signalp6 + deeplocpro PATH on fresh nodes** (task #22). Now auto-discovered via `_find_in_conda_envs` in `core/runner.py`. Binaries execute directly without `conda activate`, which works for pip-installed Python entry points because the shebang pins the env's Python. **Open risk:** if either tool's torch build needs system CUDA libs (rather than its own bundled CUDA wheels), `LD_LIBRARY_PATH` won't include `<env>/lib/` and the subprocess will crash with `libcudnn.so.X not found`. Mitigation if it ever fires: prepend `<env>/lib` to `LD_LIBRARY_PATH` in `run_signalp.py` / `run_deeplocpro.py` when the binary path lives under a conda env. Trigger to revisit: a CUDA-related ImportError from either wrapper.

## PLM-Effector performance

- **Ensemble model checkpoint cache** (task #16). `run_ensemble` now called 85× instead of 5×; re-reads ~150 MB of small files per call. Page cache makes real cost ~150 MB total, so lower priority than the raw number suggests. Revisit: if PLM-E step wallclock is still the long pole after `0445d94` cross-type caching is validated.
- **FP16/BF16 + `--batch-size 32`**. Could drop PLM-E ~74m → ~12-15m on whole-genome runs. Needs validation that FP16 doesn't shift predictions. Now low priority: PLM-E runs on the SS neighborhood (~128 proteins) by default, so absolute wallclock is small. Only revisit if `--plme-whole-genome` becomes a common workflow.

## Disk sizes (measured 2026-06-03 on CX3, scripts/audit_disk_sizes.py)

- **base 2 GB / extended 140 GB / full 1.3 TB**. BLAST nr is the long pole at 802 GB; users without nr-cross-genome BLASTp save 1.2 TB by staying on extended. UniRef30 dominates HH-suite at 261 GB. Several install.md per-tool estimates were ~50% off (EggNOG 25→47, IPS 24→35, PLM-E 18→26, HH-suite 55→340). Updated in docs/how-to/install.md.

## Torch.load safety

- **Migrate `run_deepsece.py` to `weights_only=True` + `add_safe_globals`.** DeepSecE's checkpoint is a state_dict (already loaded via `model.load_state_dict(...)`), not a whole-module pickle like PLM-Effector's, so it's a safe candidate for the stricter loader PyTorch 2.6 introduced. PLM-E itself can't migrate (whole-module saves; would need an upstream refactor). Revisit: any time a new ssign dep bumps torch and `weights_only=False` triggers a deprecation warning.

## Annotation parallel-group scheduling

- **Wave-scheduling or finish-and-redistribute** (deferred 2026-06-03). After fixing the 3-4x oversubscription bug (commit `82cece9`), each annotation tool now gets `effective_cpu_count / N` cores in the parallel group. Remaining inefficiency: when one tool (e.g. IPS at 36m) finishes before the others (EggNOG, pLM-BLAST), the cores it released sit idle because the surviving tools were started with a fixed thread count and can't add workers mid-run. Three options considered: (1) two-wave scheduling (fast tools then slow tools), (2) restart-survivors-with-more-threads (throws away expensive embed/diamond work), (3) leave as-is. Going with (3) for now. **Trigger to revisit:** once we have step_timings.csv from a few post-fix runs across different genome sizes — if there's a consistent 10-20%+ wallclock floor that's clearly tail-tool dominated, the wave-scheduling experiment becomes worth it.

## Benchmark: T1SS effector rescue (effector-recovery-benchmark, 2026-06-11)

The corpus left 19 of 28 validated T1SS effectors without a genome (old UniProt entries cite 1980s-90s EMBL gene clones, not assemblies), so they were untestable for the proximity ceiling. Per Teo's decision (representative-species placement; floor set to >=90% id after he asked to include the 93% prtA), `scripts/12_rescue_t1ss_ipg.py` (IPG, identical protein) + `13_rescue_t1ss_blast.py` (remote blastp, >=90% id / >=90% cov) + `14_finalize_t1ss_rescue.py` placed **16/19** into RefSeq genomes (9 IPG-identical, 7 BLAST; 13 exact-species, 3 genus-only). Output: `data/t1ss_rescue/t1ss_rescued.tsv`. (prtA Q07295 = 93.1% in same-species D. chrysanthemi SR64; the rest 98-100%.)

- **3 still unplaced, with reasons (genuinely unplaceable, not a decision):**
  - `prtG` (Q07162, Dickeya chrysanthemi): best genome match only **61%** id. Genuinely divergent.
  - `prtA` (P82115, Photorhabdus sp. Az29) + `lktA` (P55123, "Pasteurella haemolytica-like sp. 5943B"): "no significant similarity" because those exact strains have **no RefSeq genome**; the search was correctly species-restricted. *A genus-level fallback (flagged genus_only) could rescue the conserved leukotoxin; the Photorhabdus one is a partial UniProt entry and riskier.* Trigger to revisit: if Teo later wants >=18/19, add a genus-level BLAST fallback.
- **Adjacency verified (script 15 + 16, 2026-06-11):** instead of assuming ceiling=100%, we read the gene order around each placed effector for the T1SS transporter (ABC/HlyB-family + HlyD-family MFP; literature families off RefSeq annotation, ssign-independent). **14/16 CONFIRMED adjacent** (transporter 1-3 genes away -> reachable at N=3). **2 genuine exceptions:** apxIIA (apxIICA operon lacks a transporter; trans-secreted by ApxIB/D from the distant apxI operon) + serralysin (Lip/LipBCD transporter at a separate locus, confirmed across 13 complete Serratia genomes). These 2 are real "impossible" cases. `16_t1ss_replace_fragmented.py` fixed hlyA (was on a 3 kb single-CDS WGS contig that truncated hlyCABD -> re-placed into complete genome NZ_CP031766.1, operon intact). Tables: `data/t1ss_rescue/t1ss_ceiling.tsv`.
- **Phase 1 integration still TODO:** fold the 16 placements into `effector_gold_set.tsv` (add `placement_tier` ipg_identical/representative_strain + `species_match` columns), add the placement genomes to the Phase 2 panel, and record the T1SS ceiling as 14 reachable + 2 impossible. Placement genomes already cached in `refseq_cache/`. Trigger: start of Phase 1.
- 3 placements are sister-species (genus_only): `prtA1` P.luminescens->P.akhurstii, `prtB`/`prtC` D.chrysanthemi->D.dadantii. Sequence is 100% identical (IPG), so defensible, but flagged.

## Benchmark: Phase 1 ceiling (effector-recovery-benchmark, 2026-06-11)

Phase 1 complete: ceiling = % of testable verified effectors the +/-N proximity rule could reach. **T1 80% / T2 1-4% / T3 21-33% / T4 6-10% / T6 22-26%; 499 testable, 83 untestable.** Scripts 17-21 + `scripts/bench_index.py`; outputs in `data/phase1/`, figures `figures/01-04`, writeup `docs/phase1_ceiling.md`. Deferred / judgment-call items:

- **33 net-new SecReT6 effectors dropped to untestable** (Teo's Checkpoint-A-consistent call): 25 in multi-T6SS genomes + 8 in Serratia SMDB11 with no curated T6SS instance. We do NOT guess the nearest instance (circular). T6SS ceiling rests on the 72 rigorously-assignable effectors. **Trigger:** if Phase 2 makes T6SS coverage worth expanding, literature-assign each net-new effector to its specific T6SS (SecReT6 citations), then re-run 19.
- **pulA (T2SS_16) untestable.** Its instance was dropped at Checkpoint A (X12831 = 1.8 kb EMBL pul fragment, not a genome); not backfilled. The founding T2SS effector is therefore not in the ceiling. **Trigger:** if wanted, place pulA into a complete K. pneumoniae genome + curate the gsp machinery there (real Phase 0b work).
- **26 effector_locus_not_found** = corpus locus_tag scheme absent from the available RefSeq assembly (different assembly/annotation), no unique gene symbol to bridge. Genuinely unplaceable in cache. **Trigger:** coordinate-based placement via UniProt->genome if completeness matters.
- **Gene-symbol fallback (9 effectors, incl. canonical Yersinia Yops on pYV).** Located by *unique* /gene symbol when the tag scheme was missing; flagged `effector_match=gene_symbol` in `ceiling_per_effector.tsv`. Uniqueness-gated (no paralog ambiguity), but a weaker identifier than locus_tag. Documented in `docs/phase1_ceiling.md` strengths.
- **off-replicon machinery -> impossible.** Effectors whose own-instance machinery anchored only to a different replicon are counted impossible (structurally unreachable by a +/-N window), not untestable. Reason `machinery_off_replicon` in the per-effector table.

## ssign output: emit genomic coordinates (base pipeline)

- **`results.csv` and `results_raw.csv` carry no genomic coordinates** (only `locus_tag`, `sequence`, per-tool signals). The intermediate `gene_info` step already has `contig,start,end,strand` (`extract_proteins.py:384`) but `_build_master_csv` / `_build_raw_csv` (`core/runner.py`) drop them. **Add `contig,start,end,strand` to the raw CSV at least** (left-join from gene_info, which is already the raw base). **Why it matters:** any coordinate-based downstream use (operon context, the effector-recovery benchmark's Bakta->RefSeq bridge, the secretion-classifier model's positional features) currently has to recover coordinates by re-parsing the input or matching on sequence. The benchmark Phase 2 worked around it with locus_tag + protein-sequence identity matching, but emitting coordinates removes the workaround and is generally useful. **Trigger:** next base-pipeline pass; small change, gene_info is already the join base in `_build_raw_csv`. (Flagged 2026-06-11 during effector-recovery Phase 2.)

## Statistics

- **(Resolved 2026-06-02)** The broken permutation + biased Fisher path was replaced by the A+ rewrite: opt-in `--enrichment-stats` flag, null sample of N=200 random non-SS-neighborhood proteins per genome, scipy binomial test per real SS system + per broad type, BH FDR. Multi-genome runs also emit a pooled view. See `enrichment_testing.py` + `pool_enrichment_stats` in `core/runner.py`.

## Shared TSV/parsing helpers (#75a simplify follow-up, resolved 2026-06-05)

- **Tolerant int parsing — resolved.** `ssign_lib/parsing.py:parse_int_or_none(value, allow_range=False)` is the single source. `t5_passenger` imports it directly; `t5ss_handler._parse_sp_end` is now a one-liner wrapping it with `allow_range=True` so call sites still read as `_parse_sp_end(...)`.
- **"Load TSV → dict by key" — resolved.** `ssign_lib/tsv_io.py:load_tsv_by_key(path, key_columns, missing_ok=True)` is the single source. `cross_validate_predictions._load_tsv_by_locus` passes the tolerant fallback chain `("locus_tag", "protein_id", "seq_id")`; `enrichment_testing.load_predictions_keyed` keeps strict-locus_tag-only via `key_columns=("locus_tag",)`. `t5_passenger.load_t5_classifications` was left as-is — it does substantial per-field coercion on top of the TSV read, so the shared helper would only save the open-and-iterate boilerplate.

## secretion-classifier-dataset — deferred re-citation (group 1, 2026-06-11)

- **25 predicted positives kept with UNRESOLVED citations** (`evidence_tier=predicted`, `citation_status=UNRESOLVED` in `data/dataset/predicted_audited.tsv`). The audit kept them (protein/locus/instance label is sound; the broken DOI is metadata, not a label, same bar as the validated gold set which also never drops on citation alone). 8 distinct broken DOIs: lktA/mbxA/movA RTX leukotoxins (11), prtA/prtC proteases (4), Ehrlichia T4SS EBP/ECH (3, IAI.00513-13), Shigella OspD2/3 (2, science.1175302), Anaplasma T4SS APH (2, cmi.13405), Chlamydia CT813 (1), V. vulnificus rtxA (1, JB.187.10.3392), Dickeya pelZ (1, annurev-micro-102215). 36 of the original 61 were already recovered (Apx ×32 → Frey 1993 10.1099/00221287-139-8-1723; ApxIVA ×2 → Schaller 1999; V. cholerae rtxA ×2 → Boardman 2004) in `doi_recite.jsonl`. **Trigger:** before any public release that cites the predicted tier, run a per-paper re-citation for these 8 (one paper each, all well-characterized families); add verified DOIs to `doi_recite.jsonl` and re-run scripts 31→32. Held back because applying a shaky umbrella citation to benchmark provenance is worse than an honest flag. The agent re-citation route hit an API policy block on toxin literature, so this was done inline; the remaining 8 need careful manual lookups.

## secretion-classifier-dataset — predicted instance audit flags (group 2, 2026-06-11)

The 2.3 literature-audit agent (27 RESOLVED / 17 UNRESOLVED of 44 ambiguous) flagged 3 corpus data-quality issues to fix upstream (assignment stands at gene/family level; not blocking):
- **C. rodentium ICC168 NleG locus_tags mislabeled.** Corpus has ROD_31501/ROD_21621/ROD_25881/ROD_15971; NCBI annotates those as exported-protein/hisA/prophage/hypothetical. Real ICC168 NleG tags are ROD_16511 / ROD_48891 / ROD_40971. **Trigger:** correct these 4 locus_tags in `T3SS_verified.tsv` before feature-join (else bench_runout won't find them in run output).
- **VopX / VPA1374 name-locus mismatch.** UniProt annotates VPA1374 as "Uncharacterized"; "VopX" is a V. cholerae effector. Left UNRESOLVED. **Trigger:** verify the V. parahaemolyticus locus for VopX (or drop if a mislabel).
- **EC042_4675 (Tle3_Sci2) unresolvable.** Locus returns zero NCBI/UniProt hits and sits ~100 genes outside both E. coli 042 T6SS clusters; documented EAEC 042 Tle is Tle1 (Sci-1). Left UNRESOLVED. **Trigger:** re-source the correct locus from the bioRxiv ref (10.1101/2025.02.11.637775) or drop.

## ssign extended tier is GPU-gated (PLM-Effector cannot fall back to CPU via runner)

- **Finding (2026-06-11, while diagnosing CX3 pilot queue):** `run_plm_effector.py` defaults `--device cuda` and `plm_effector/predict_api.py:_resolve_device` **raises RuntimeError if CUDA is absent** (only `device='cpu'` explicitly avoids it). `core/runner.py:_step_plm_effector` (~line 2377) builds the PLM-E argv WITHOUT a `--device` flag, so it always uses the cuda default. Net: **extended/full tier crashes at the PLM-Effector step on any GPU-free node** (CPU-only HPC nodes, Amine's Mac, laptops). DeepSecE already auto-selects CPU (`run_deepsece.py:_select_device`), so PLM-E is the lone hard GPU gate.
- **Impact:** can't schedule benchmark Phase-2 runs on the abundant CPU queues (must wait for gpu72 GPU). Also contradicts the "runs anywhere / zero-maintenance" longevity pitch and breaks Mac/laptop extended runs.
- **Fix:** add a `plm_effector_device` config knob + env (mirror DeepSecE's `_select_device` auto-detect: cuda if available else cpu) and pass `--device` through in `_step_plm_effector`. Keep cuda the default when a GPU exists (CPU is ~1-2 min/protein per the code's own warning, but the +/-3 neighborhood keeps protein counts small). **Trigger:** next base-pipeline pass, or sooner if the CX3 GPU queue stays multi-hour and we want to run the panel on CPU nodes.

## secretion-classifier-dataset — T5SS sourcing gaps (group 3, 2026-06-11)

T5SS sourced by subtype (scripts 36/37): 29 proposed → **23 placed + verified** (DOI.org-resolved + locus in genome), folded into `positives_all.tsv` (930 total). Subtypes: T5a 7, T5b 11, T5c 3, T5d 1, T5e 1. self_secreted=true for the 12 autotransporters (T5a/c/d/e), false for the 11 T5b TpsA substrates.

Deferred:
- **Autotransporter sourcing was cut short** by a wifi/socket drop after 18 entries (agent was targeting ~25-40). The brief's anchor list has more uncovered T5a/T5c examples (SPATEs Pic/Sat/Vat, IgA proteases, more TAAs: Hia/UspA2/BadA/SadA). **Trigger:** if more T5SS depth is wanted, re-run one autotransporter sourcing agent (it appends to `t5ss_raw/`, then 36→37 re-fold). Low urgency — T5SS is a supplementary/partly-exploratory tier.
- **6 unplaceable** (`t5ss_unplaceable.tsv`): hap (HI_0248 absent from NC_000907.1 — likely wrong tag), tsh (APECO1_O1CoBM73 tag mismatch), hia/badA/uspA1/invA (agent gave no usable locus_tag / no genome accession). **Trigger:** re-source these specific loci if T5SS coverage matters; each needs one targeted NCBI lookup.
- **GCF→nucleotide resolution is NCBI-flaky.** `36.replicons_for` uses assembly→elink→nuccore but NCBI eutils dropped connections (EOF) repeatedly; the 5 Haemophilus T5b were rescued by hardcoding GCF→NC in `tps.json` (NC_000907.1/NC_002940.2/NC_017452.1). If re-sourcing, prefer instructing agents to emit nucleotide (NC_/NZ_) accessions, not GCF.

## secretion-classifier-dataset — model-handoff deferrals (group 5, 2026-06-11)

`positives_all.tsv` (930) is the label-complete deliverable; `data/dataset/README.md` documents schema + conventions. Items below belong to later changes, not this one:

- **validated:predicted weight ratio is unset.** `evidence_tier` is recorded but the loss weight is a model-training decision. **Trigger:** the model-training change (nnPU loss); start from validated=1.0 and sweep predicted ∈ {0.3-0.7}.
- **Type-level positives (298 rows, `type_level=yes`) have no instance, so pair-features are null.** Whether they feed the instance model or only a separate type-level head is undecided. **Trigger:** model architecture design; if instance-only, hold the 298 out or route them to the type-level head.
- **T5SS subtype depth is a first pass** (23 placed: T5a 7 / T5b 11 / T5c 3 / T5d 1 / T5e 1). Depth decision (source more vs leave as supplementary) waits on whether the self-secreted tier helps. **Trigger:** after the first model eval shows whether self_secreted rows move the needle; if yes, re-run autotransporter sourcing (see group-3 note).
- **T5a-neighbor DLP/DSE observation stays exploratory.** The "autotransporter neighbours look secreted" signal lives in benchmark task 6b as a side-study, never a training label. **Trigger:** none for this dataset; revisit only if 6b finds it predictive.
- **Feature side (group 4) is blocked on CX3 Phase-2 panel runs.** `training_dataset.tsv` (positives + DLP/DSE/SignalP/PLM-E + ESM ref + gene-distance + PU unlabeled set) joins `results_raw` via `bench_runout` once runs land. **Trigger:** Phase-2 pilots flip Q→R on CX3 gpu72; then rsync runs back and run group-4 scripts (4.1-4.5, not yet written). The 4.5 assembler's output schema is now **pinned**: it must emit the columns in `secretion-classifier/secretion_classifier/schema.py` (REQUIRED_COLUMNS), which the model's `TrainingDataset` loader validates against. 4.2 pair-features (`pair_features.tsv`) already supplies the pair/system columns.
- **Model-prep core is built** (2026-06-11, separate repo `reidmat/secretion-classifier`, commits f641d32..3b04d54): data-independent loader/losses/prior/model/splits/metrics/sweep, 33 tests green. Only the trainer + ESM extraction + `training_dataset.tsv` remain, all gated on the same CX3 runs. See memory [[project_secretion_classifier_model]].
- **Label sentinels are a bare-string cross-file contract.** `evidence_tier`, `instance_source`, `type_level`, `self_secreted` values are bare literals shared across scripts 36/37 and (soon) group 4 (e.g. `self_secreted == "true"` in 37 depends on 36 writing `str(...).lower()`). A typo mislabels rows silently. **Trigger:** when writing the group-4 feature scripts, hoist these into a small shared label-constants block (e.g. in `bench_io.py`) and reuse across 36/37/group-4 rather than re-typing literals.
