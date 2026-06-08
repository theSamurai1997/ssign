# Data audit and sourcing pipeline

Findings from the 2026-06-08 data audit and dedup-pipeline research
agents.

## Available curated substrate data (verified 2026)

Raw positive counts per SS type, before dedup:

| Source | T1SE | T2SE | T3SE | T4SE | T6SE | Notes |
|---|---|---|---|---|---|---|
| **PLM-Effector** training+test | 181 | 80 | 571 | 585 | 327 | CD-HIT 90%; anonymized FASTA IDs (no species) |
| **DeepSecE** training | 128 | 68 | 406 | 507 | 232 | CD-HIT 60%; anonymized IDs |
| **SecReT4** verified | 0 | 0 | 0 | 540 | 0 | Gene names + species kept |
| **SecReT6** experimental | 0 | 0 | 0 | 0 | 331 | Gene names + species kept |
| **NR union (exact-seq)** | **181** | **80** | **571** | **705** | **610** | After cross-DB exact dedup |

**Hard finding: PLM-Effector is a 98-100% exact-sequence superset of
DeepSecE.** Treat PLM-E as the canonical positive set. SecReT4 adds
120 genuinely new T4SE; SecReT6 adds 283 new T6SE. Pooled NR positive
total: ~2,150.

After 30%-identity MMseqs2 clustering (the recommended dedup
threshold), expect 600-1,200 NR positives total. T1SE and T2SE remain
the data-limiting classes.

## BastionHub

BastionHub (Monash University, NAR 2020 49:D651) was an aggregator
database that combined curated substrates from earlier sources. As of
2026-06, the site returns 503 (likely retired). Wayback Machine has
the homepage HTML but not the FASTA downloads.

**Why this matters:** PLM-Effector's training set was partly built
from BastionHub records. Without the original, we can't perfectly
replicate their data composition. **In practice the impact is small**:
the same underlying substrates flow forward into SecReT4/6, DeepSecE,
and PLM-E, which are all still accessible. We can build a comparable
training set from those four sources without BastionHub.

**Action if we need it later:** email the Monash ERC group to request
a snapshot. Pre-archive the data we have on Zenodo immediately so we
don't repeat the loss.

## Re-BLASTing for source-genome metadata

The PLM-E and DeepSecE FASTAs use anonymized integer IDs (`>T4SE_28`).
No species or accession info. To attach genomic-context features
(MacSyFinder hit, DLP score, SignalP class) we need to know each
protein's source genome.

**Workflow** (from data-sourcing agent):
```
mmseqs easy-search all_anonymized_positives.fasta swissprot \
    hits.m8 tmp \
    --min-seq-id 0.95 -c 0.9 --cov-mode 0 -s 7.5 \
    --format-output "query,target,pident,evalue,theader" \
    --threads 32
```

- ~1500 anonymized proteins vs SwissProt: ~10 min on GPU node.
- Expected hit rate >=95% identity: 60-80% (known effectors are
  over-represented in SwissProt).
- Fallback chain for misses: UniRef90, then TrEMBL, then NCBI RefSeq.
- **Do not fall back below 80% identity**: below that the "source
  genome" is a guess and will pollute genome-context features.

**Decision**: re-BLAST is needed for the main training set (we need
genome context for those ~1500 proteins). NOT needed for the
sequence-only POC (where we just need the FASTA and the per-type
label). At <10 min compute it's free; do it once and cache the
mapping.

## Recommended dedup pipeline

Tool: **MMseqs2 v17+** (not CD-HIT; CD-HIT loses sensitivity below 50%
identity).

**Important nuance**: in the multimodal architecture, the same
substrate sequence from two different source genomes is NOT a
duplicate training example, because its tabular features (MacSyFinder
hits in the surrounding genes, DeepLocPro score in its own proteome,
genomic-context features) are different. Keep cross-genome
occurrences as separate training examples. They contribute real new
information.

The only deduplication step that runs at the protein level:

1. Pool all positives per SS type, prefix IDs with source DB.
2. Exact-dedup at 100% identity AND same source genome only (catches
   trivial within-DB duplicates from format issues; not cross-genome).

Splitting is handled at the cluster level by GraphPart so that
cross-genome orthologs of the same protein land in the SAME fold
(train OR test, never split). This preserves test-set integrity
without dropping training data:

3. Cluster all positives at **30% identity, 80% bidirectional
   coverage** (`--min-seq-id 0.3 -c 0.8 --cov-mode 0`). This produces
   clusters but does NOT drop sequences.
4. Pass the clusters + (pos/neg) labels to **GraphPart** (Teufel et
   al. NARGAB 2023) which assigns each cluster to train/val/test such
   that no cross-partition pair exceeds 30% identity. Every member of
   a cluster goes to the same fold.
5. Also do a leave-one-genus-out eval (e.g., hold out *Legionellales*
   for T4SE) as a complementary test of cross-genus generalization.

Report per-fold MCC + PR-AUC, not accuracy.

## Negative set construction

Two sources:

1. **DeepSecE's curated non-effectors**: 1,577 proteins, available in
   the TXSE-Dataset.tar.gz download.
2. **Genome-context negatives**: every protein in a validation genome
   that does NOT match (>=50% id, >=70% cov) any pooled positive.

**Critical filter** (EffectorP 3.0 recipe, applied per protein):
keep only candidates where:
- DeepLocPro extracellular probability > 0.3, OR
- DeepSecE per-type score > 0.3, OR
- SignalP signal-peptide positive, OR
- has TM helix

This drops ribosomal / DNA-pol / metabolism proteins (~85% of any
genome) and leaves ~5k informative negatives across 5-10 validation
genomes. The classifier then has to actually learn the discriminative
signal, not just "ribosomal != effector."

See `04_negative_sets.md` for the full PU-learning recipe.

## Full pipeline runtime

On CX3 v1_gpu72 (32 cores, L40S or RTX 6000):
- Clustering: ~5 min total across all SS types
- SwissProt provenance search: ~10 min
- GraphPart needle splitting (slowest step, O(N^2)): 30-60 min for
  T4SE (largest, ~700 reps after clustering)

**Whole pipeline: under 2 hours on one node.** Build the dataset once,
cache the result, train any number of models off it.

## Reproducibility

Pin tool versions in `scripts/build_dataset.sh`:
- mmseqs2 17.b804f
- graph-part 0.1.4
- seqkit 2.8.2
- emboss 6.6.0
- biopython 1.84

Mirror downloaded raw FASTAs to Zenodo with a SHA-256 manifest so the
SecReT / DeepSecE / PLM-E URLs going down later doesn't kill
reproducibility.
