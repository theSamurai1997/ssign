## 1. Predicted-effector audit (unblocked)

- [x] 1.1 Extract the 347 `predicted`-evidence rows from the source corpus into a working table; log per-type counts. (`30_extract_predicted.py`: 347 → 171 VERIFIED + 154 PARTIAL kept, 17 status-drop, 5 biology-drop.)
- [x] 1.2 Verify each predicted row to the gold-set bar: DOI resolves (Crossref, reuse the script 02/03 + `.doi_cache` pattern) + UniProt/locus_tag cross-check; emit `predicted_audited.tsv` + a provenance table with each row's fate. (`31_verify_predicted.py`: 325 kept, 0 unplaceable; 61 broken-DOI rows re-cited via agent under verbatim-quote+Crossref contract.)
- [x] 1.3 Combine surviving predicted with the validated gold set into a tiered positive table with an `evidence_tier` column (validated|predicted); report retained/excluded counts per type. (`32_tiered_positives.py`: 907 positives = 582 validated + 325 predicted.)

## 2. Predicted instance assignment (unblocked)

- [x] 2.1 For each verified predicted effector, count same-type instances in its genome (reuse `instances.tsv` + `bench_index` drift-tolerant genome match). (`33`: 139 in enumerated genomes, 186 not.)
- [x] 2.2 Auto-assign the single-same-type-instance cases; mark `instance_source=auto`. (95 auto-assigned.)
- [x] 2.3 Literature-audit agent for multi-instance genomes: read the sourcing DOI, attempt to resolve the specific instance with a verbatim quote; mark `instance_source=literature` on success. Anti-hallucination contract; no nearest-machinery guessing. (44 ambiguous → 27 RESOLVED w/ verbatim quote, 17 UNRESOLVED; all ids in-candidate-set; 3 corpus tag flags in NOTES.)
- [x] 2.4 Emit unresolved rows as instance-unknown type-level positives (`type_level=yes`, pair-features null); write `positives_instanced.tsv` + provenance. (`35`: 95 auto + 27 literature + 203 type-level = 325.)

## 3. T5SS sourcing (unblocked)

- [x] 3.1 Write per-subtype apparatus/effector briefs (T5aSS, T5bSS, T5cSS; note T5dSS/T5eSS) with seed Pfam families (PF03797, PF03895, TpsA/TpsB-POTRA) and the label convention (self_secreted true for a/c, false for b). (`data/dataset/t5ss_sourcing_brief.md`.)
- [x] 3.2 Agent sourcing (one per subtype) under the anti-hallucination contract: verbatim quote + resolvable DOI + real locus_tag; emit `t5ss_raw/<subtype>.json`. (2 agents: `autotransporters.json` 18 self-secreted T5a/c/d/e, `tps.json` 11 T5b. Autotransporter agent cut short by a wifi/socket drop after 18 — incremental writes saved them.)
- [x] 3.3 Resolve each sourced gene to a RefSeq locus_tag + coordinates (reuse the machinery-resolver approach, script 10 pattern); fetch any missing genomes into the cache. (`36`: GCF-assembly→nucleotide-replicon resolution + 09_fetch_refseq; old_locus_tag-aware.)
- [x] 3.4 Independent verification pass (mirror benchmark 3.6): DOI registered + locus_tag exists in genome + gene named in quote; emit `t5ss_effectors.tsv` with `ss_type`, subtype, `self_secreted`, anchor status. (`36`: 23/29 placed, all 23 verified, DOI.org-checked; 6 unplaceable flagged.)
- [x] 3.5 Fold the placed T5SS rows into the positive table (gold-set schema); record unplaceable ones separately. (`37`: `positives_all.tsv` = 930 = 907 + 23 T5SS; `t5ss_unplaceable.tsv` holds the 6.)

## 4. Feature matrix assembly (gated on Phase 2 panel runs)

- [ ] 4.1 Per (protein, instance) positive, locate the protein in ssign run output via `bench_runout` (locus_tag / sequence / >=90% identity bridge) and pull DLP/DSE/SignalP/PLM-E signals; report feature-unavailable rows with reasons.
- [x] 4.2 Attach pair-features (gene-distance to the assigned system's nearest machinery locus, from the gene-order index) and system-features (SS type, component count); null pair-features for instance-unknown positives. (`40_pair_features.py` → `pair_features.tsv`, 547/930 placed: 456 ceiling + 91 computed; 261 type-level + 12 self + 110 unreachable null. Run-independent. Also canonicalized validated `sys_instance_id` from the ceiling so the join keys one id space.)
- [ ] 4.3 Add the ESM embedding reference per protein (embedding computed/cached separately; store a stable id/path, not the vector inline).
- [ ] 4.4 Build the PU unlabeled candidate set: non-effector proteins passing the EffectorP-style secretion filter (DLP/DSE/SignalP/TM), minus high-identity-to-positive homologs (MMseqs2/`bench_index`); label unlabeled, not negative.
- [ ] 4.5 Emit `training_dataset.tsv` (positives + unlabeled, with evidence_tier, labels, features, provenance) + a dataset card; report pending-run genome count.

## 5. Wrap-up

- [x] 5.1 Dataset README: schema, label conventions (incl. self_secreted, type_level), evidence tiers, PU set semantics, and the run-dependency. (`data/dataset/README.md`.)
- [x] 5.2 Record deferred items (model training, weight ratio, T5SS depth, T5a-neighbor exploratory link to benchmark 6b) in NOTES/backlog with triggers. (NOTES "model-handoff deferrals" block.)
- [x] 5.3 Simplify pass over the new scripts; verify reproducibility from the benchmark outputs. (4-agent review; fixed `37._key` non-uniqueness + assert, `36` esearch field-tag/uniqueness guard, double-`find` walrus. **Found+fixed a real bug:** `37` was folding the pre-instance predicted branch, dropping the group-2 instance assignment + drifting the canonical id column; now unions validated+predicted-instanced+T5SS with one `sys_instance_id`. Offline chain 32→33→35→37 reproduces `positives_all.tsv` byte-identical.)

## 6. Citation-consistency sweep (ssign-found scope)

Motivated by the benchmark's 19-row discordant audit (`../effector-recovery-benchmark`, RESULTS.md):
the answer key's effector-sourcing DOIs are frequently wrong (resolve to unrelated papers) or 404,
even on rows marked VERIFIED. That is a training-label provenance defect. Scope this first pass to
the ~51 ssign-FOUND gold effectors (Teo 2026-06-12, save time; not the full 582). All 19 already-
audited rows fall inside this set, so the deterministic check is validated against the manual one.

- [x] 6.1 Pass-1 deterministic check (`scripts/41_citation_consistency.py` → `citation_consistency_found.tsv`).
  **Join-bug fix (important):** the first cut joined found→positives on raw UniProt, but ~1/3 of found
  effectors carry `uniprot='-'`, and `'-'` is truthy, so they all collapsed onto one arbitrary positives
  row (12 shared DOI `nature11433`) → garbage organism/DOI/verdict for 17 rows. Fixed: bridge
  found→`ceiling_per_effector` by (effector_locus, gene) to get the gold `instance_id`, then key positives
  by (instance_id, gene, ss) (UniProt / gene+ss+genome as fallbacks). All 51 now resolve via instance,
  **35 distinct DOIs** (was badly collapsed), recorded in a `join_method` column. Corrected CrossRef
  content check: **21 CONSISTENT** (14 by gene name, 7 genus-only), **8 FLAG_WRONG_TOPIC** (DOI in an
  unrelated field: celA/plaA→soil ecology, Map×2→Toxoplasma, VirA→tooth development, BipB→colicins,
  aprA→CS1 pilin, **EspZ→the T6SS-discovery PNAS paper**), **3 FLAG_GENE_ABSENT** (prtB/prtA1/aprA cite a
  real Prt/Lip ABC-exporter paper that just doesn't name the protease in its abstract), **3 DOI_UNRESOLVED**
  (YspE/CopN/ChlaDub1 — classic-publisher DOIs that fail the DOI.org handle check; may be legacy
  false-negatives or typos, pass-2 settles), **16 INDETERMINATE** (CrossRef record exists but no abstract
  → can't refute deterministically). Validated against the 19 prior manual verdicts: re-flagged
  celA/plaA/VirA independently AND caught **BipB's wrong DOI the manual pass missed** (its agent was
  AUP-blocked; only biology, not the DOI, had been checked). Verdict enum encodes WRONG_TOPIC vs
  GENE_ABSENT vs FETCH_ERROR so the table is self-interpreting.
- [x] 6.2 Targeted re-audit of the 30 non-CONSISTENT rows (`pass2_input.json` → `pass2_raw/batch_*.json`
  → `42_pass2_verify.py` → `pass2_results.tsv` + `pass2_results.md`). 5 literature agents under the
  verbatim-quote contract + batch C (YspE/VirA/BipB/CopN/ChlaDub1) done in-session via PubMed (agent
  AUP-blocked on B. pseudomallei). Every returned DOI deterministically re-verified (registered +
  on CrossRef + gene/genus present): **26 rows now carry a verified-real, on-topic primary DOI** (21
  RESOLVED to a corrected DOI, 3 CONFIRMED, 2 MISASSIGNED), **0 UNVERIFIED/fabricated DOIs**. Dataset
  edits flagged: **2 MISASSIGNED ss_type** (BopA/BopE T6SS→T3SS, Bsa effectors, name-collision);
  **4 NOT_FOUND drop/down-tier** (ChlaDub1 no T3SS-secretion evidence; P.entomophila aprA homology-only;
  EFF00136/EFF00150 DUF3274 not in cited paper); **1 duplicate** (Tle4 = TplE_alias_Tle4 = PA1510, both
  also had wrong DOIs → Jiang 2016). EFF00142's prior "unidentifiable" tag was wrong (Russell 2012 names
  BTH_I2691 → CONFIRMED). The agents agreed with the deterministic re-check on all 26 sourced rows.
- [x] 6.3 Apply the audit to the training labels (`43_apply_citation_corrections.py`, post-build overlay
  on `positives_all.tsv`, 930→925). Triage by training impact: **22 pure-provenance fixes** (corrected DOI
  only, label unchanged); **2 label-changing** (BopA/BopE T6SS→T3SS, re-labelled to type-level T3SS
  positives + DOI fixed); **5 removed** (4 unsupported + 1 duplicate, quarantined to
  `positives_removed_citation.tsv`). Backup `positives_all.pre_citation_audit.tsv`; full change log
  `citation_corrections_log.tsv`. Idempotent + reproducible (re-run from backup is byte-identical).
  Overlay must re-run after any `37` rebuild; `40_pair_features.py` needs re-run (2 ss_type changes + 5
  drops). **Out of scope (noted):** other non-found instances may still carry the same wrong DOIs (e.g.
  a second EspZ instance still cites the T6SS-discovery paper) — that's the full-corpus sweep, deferred.
