# Gap-validation experiment

Quantify how many true substrates ssign's proximity filter misses,
per SS type, on well-characterized bacterial genomes.

## Hypothesis

If ssign's default proximity filter (±3 genes around each MacSyFinder
component) misses a substantial fraction of literature-validated
substrates, then a sequence-aware classifier that uses tool outputs
as features could recover them.

Null result: proximity-filtered ssign recovers ≥95% of literature-
validated substrates per genome per SS type. Project ends.

## Validation genomes

Six well-characterized genomes selected from the data audit
(`02_data_audit.md`, Table 2):

| Genome | SS focus | Catalogued | Validated | RefSeq |
|---|---|---|---|---|
| *Legionella pneumophila* Philadelphia 1 | T4SS Dot/Icm | ~330 | ~150-200 | NC_002942 |
| *Coxiella burnetii* RSA 493 (Nine Mile) | T4SS Dot/Icm | ~150 | ~80 | NC_002971 |
| *Salmonella* Typhimurium LT2 | T3SS-1 / T3SS-2 | ~40 | ~35 | NC_003197 |
| *Yersinia pestis* CO92 | T3SS Yop | 7 | 7 | NC_003143 + NC_003131 |
| *Pseudomonas aeruginosa* PAO1 | T3 + T6 | ~3 T3, ~25 T6 | ~15 | NC_002516 |
| *Vibrio cholerae* N16961 | T2 + T6 | ~10 T2, ~6 T6 | ~12 | NC_002505 + NC_002506 |

The substrate ground-truth tables are being curated by a background
agent (output pending). When delivered they live in
`data/ground_truth/<genome>.tsv` with columns: substrate_name,
uniprot_id, locus_tag, ss_type, evidence_level, primary_reference.

## Two conditions per genome

A. **Proximity-filtered** (ssign default)
B. **Whole-genome**: `--dlp-whole-genome --dse-whole-genome --sp-whole-genome --plme-whole-genome`

12 jobs total (6 genomes × 2 configs).

## Scoring

For each genome, for each SS type:

1. Take the literature substrate list.
2. Map ssign-predicted substrates back to UniProt / locus tag.
3. Compute:
   - `recall_proximity = |proximity_hits ∩ ground_truth| / |ground_truth|`
   - `recall_whole_genome = |whole_genome_hits ∩ ground_truth| / |ground_truth|`
   - `gap = recall_whole_genome - recall_proximity`
   - `fp_blowup = |whole_genome_hits| / |proximity_hits|` (substrate-level)

4. Stratify by evidence level: validated-only vs all-catalogued.

## Decision rule

- Average `gap` over six genomes, T1-T6 pooled.
- **<5%** → null result. Document, stop, save the multimodal-classifier
  project for someone else.
- **5-20%** → look at *which* substrates are missed. If they are
  niche-pathogen-specific (e.g., orphan effectors in Legionella), maybe
  not worth chasing. If they're well-known classes ssign systematically
  misses, proceed.
- **>20%** → strong motivation. Proceed to POC classifier.

## Compute

12 ssign runs at extended tier, 32c/64gb/RTX6000. ~45 min each = ~9
GPU-h total, fits in one CX3 evening's worth of throughput.

## Output

`research/secretion_classifier/data/gap_results.tsv` — one row per
(genome, ss_type, condition) with all metrics. Notebook at
`scripts/score_gap.py`.
