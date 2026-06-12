# Phase 1: the proximity ceiling

How many experimentally-verified effectors *could* ssign's "secreted protein lives within
±N genes of a secretion-system component" rule ever recover, if detection were perfect? This
is the ceiling. It needs no ssign run: it depends only on where each effector sits in gene
order relative to its own system's machinery. Phase 2 measures what ssign actually recovers
and compares against this ceiling.

## Result

Ceiling = % of **testable** effectors whose nearest own-instance machinery component is within
N genes. Denominator excludes untestable effectors (see below); it is never padded with them.

| SS type | gold | testable | untestable | ceiling N=3 | N=5 | N=7 |
|---------|-----:|---------:|-----------:|------------:|----:|----:|
| T1SS    |  28  |   25     |   3        |   80%       | 80% | 80% |
| T2SS    |  83  |   77     |   6        |    1%       |  3% |  4% |
| T3SS    | 237  |  227     |  10        |   21%       | 29% | 33% |
| T4SS    | 111  |   98     |  13        |    6%       |  8% | 10% |
| T6SS    | 123  |   72     |  51        |   22%       | 24% | 26% |
| **ALL** | 582  |  499     |  83        |   18%       | 23% | 25% |

Figures: `figures/01_gold_set_composition.png` (overview), `02_ceiling_vs_window.png`
(headline), `03_distance_ecdf.png` (the distance distribution behind the ceiling),
`04_per_genome_ceiling.png` (drill-down). Tables: `data/phase1/ceiling_summary.tsv`,
`ceiling_by_genome.tsv`, `ceiling_per_effector.tsv`.

## What the numbers mean

The proximity rule is a strong constraint for some systems and almost useless for others, and
that is the point of measuring the ceiling rather than assuming one.

- **T1SS ~80%.** T1SS is operonic: the effector is encoded beside its own ABC transporter +
  membrane-fusion-protein. We confirmed adjacency directly from gene order (the rescued-T1SS
  adjacency scan, §below). The two impossible cases are real biology, not noise: apxIIA is
  *trans*-secreted by a transporter at a distant operon, and serralysin's Lip transporter sits
  at a separate locus.
- **T2SS ~1-4%.** T2SS substrates (secreted hydrolases, toxins, e.g. plcH, plcB) are scattered
  across the genome, thousands of genes from the gsp/Xcp apparatus. Only effectors of the
  accessory Hxc system (lapA, 2 genes away) and a few others fall in window. The proximity rule
  fundamentally cannot recover most T2SS effectors; this is the clearest negative result.
- **T4SS ~6-10%.** Likewise scattered: Legionella/Coxiella-style effector repertoires are
  distributed genome-wide, not clustered at the apparatus.
- **T3SS ~21-33%** and **T6SS ~22-26%.** Intermediate. A clustered minority (LEE-encoded T3SS
  effectors, T6SS auxiliary-cluster effectors near Hcp/VgrG/PAAR) sits in window; the rest are
  dispersed.

## Method

1. **Gene-order index** (`bench_index.py`, `18_build_gene_order.py`). Every cached RefSeq
   genome is parsed into per-replicon CDS lists ordered by coordinate; each CDS gets a 0-based
   ordinal. Distance between two loci = `|ordinal difference|` on the same replicon, which is
   exactly the rule's "±N genes". Matching is drift-tolerant: locus_tags are folded
   (underscore-insensitive, so corpus `ECs4550` = RefSeq `ECs_4550`), and genome accessions are
   resolved across version and RefSeq-prefix drift (`NC_002516` = `NC_002516.2`, `HG326223` =
   `NZ_HG326223.1`).

2. **Distance to own machinery** (`19_effector_distance.py`). Each effector is matched to its
   system instance (`sys_instance_id`), and its distance is taken to the nearest *anchored*
   machinery component of **that same instance**, from the literature-derived answer key
   (Phase 0b). Machinery positions are read from the genome annotation, **never from
   MacSyFinder** (ssign's own detector), so the ceiling is ssign-independent and not circular.

3. **The 16 rescued T1SS effectors** (Phase 0a, §`tasks.md` 4.3) carry their distance from the
   adjacency scan instead: their literature machinery is unknown, so we read the HlyB/HlyD-family
   transporter straight off the gene order around the placed effector. Confirmed adjacent (1-3
   genes) for 14, genuinely non-adjacent for 2.

4. **Classify** reachable (≤N) vs impossible (>N, or machinery on a different replicon) at
   N = 3, 5, 7. Aggregate per type and per genome (`20_aggregate_ceiling.py`).

## Testable vs untestable (the honest denominator)

An effector is **untestable** when we cannot measure its distance, and untestable effectors are
reported but excluded from the ceiling fraction so it is never inflated. 83 of 582 are
untestable:

- **own_instance_unknown (25)** + **no_instance_in_genome (8)**: net-new external-DB T6SS
  effectors whose own instance can't be pinned (multi-T6SS genomes, or a genome with no curated
  T6SS). Per the Checkpoint-A decision we do **not** guess the nearest instance (that would
  circularly minimise the distance). This is why T6SS has a large untestable share.
- **no_genome (11)** + 3 unplaceable T1SS: effector has no usable genome (corpus placeholder, or
  no RefSeq assembly for the strain).
- **effector_locus_not_found (26)**: the corpus locus_tag scheme is absent from the available
  assembly (different annotation scheme/assembly), and no unique gene symbol bridges it.
- **machinery_unanchored (10)**: the effector's instance is curated but its machinery genes did
  not resolve to loci (e.g. PAO1 H3-T6SS annotated as hypothetical).

## Strengths

- **ssign-independent.** Machinery positions come from literature + genome annotation, not from
  ssign's detector, so the ceiling is a fair external bound.
- **Identity-based placement.** Effectors are placed by stable locus_tag, or by a *unique*
  gene symbol where the tag scheme is missing (9 effectors, incl. the canonical Yersinia Yops on
  pYV; flagged `effector_match=gene_symbol` in `ceiling_per_effector.tsv`). No fuzzy matching.
- **Conservative untestable handling.** Anything we cannot place rigorously is dropped from the
  denominator, not assumed reachable.

## Limitations

- **Corpus bias.** The gold set is T3SS-heavy (237) and T1SS-thin (28, 16 of them rescued into
  representative same-species genomes). Per-type ceilings are robust; a single "average across
  effectors" number is dominated by T3SS.
- **T6SS untestable share is large** (51/123) because of the net-new effectors held at
  Checkpoint A. The T6SS ceiling rests on the 72 rigorously-assignable effectors.
- **Ceiling, not recall.** This says nothing about whether ssign detects the machinery or emits
  the protein. That is Phase 2; actual recall is bounded above by these numbers.
