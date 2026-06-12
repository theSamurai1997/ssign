# Effector-recovery benchmark — precision estimate (full panel, 2026-06-12)

The recall result (`RESULTS.md`) asks how many *known* effectors ssign emits. Precision asks the
opposite: of everything ssign emits as a substrate, how much is real? ssign emits **1,933** substrate
calls (as-shipped default, T3SS excluded) / 2,321 (T3SS included) across the 67 panel genomes; only
51 are gold effectors. There is **no ground-truth negative set**, and whether a novel protein is a
real secreted effector is often unprovable, so precision is reported as a **deterministic range, not a
single number** (agent-sampled adjudication, which could narrow it, was deferred).

Reproduce: `28_emissions.py` → `29_db_confirmed.py` → `30_fp_annotation.py` → `31_precision_figures.py`,
each `--run-tag panel_genbank_default` (and `panel_genbank_t3ss`). Figures in `figures/precision/01..03`.

## Scope: the proximity calls are the question

| substrate_source | n (default) | what it is |
|---|---|---|
| `proximity` | 1,572 | called a substrate because it sits within ±3 genes of a detected system — **the precision question** |
| `T5SS-self` | 361 | the autotransporter passenger domain, its own substrate — correct by construction |

T5SS-self is a detection sanity check, not a precision unknown: 68% are annotated as autotransporters/
effectors and **0.3%** as cytoplasmic, consistent with the T5SS detection being sound. Everything below
is the 1,572 proximity calls.

## Two deterministic bounds (Figure 01, 03)

**Floor, what we can prove (tier 1).** Homology (pyhmmer phmmer, ≥90% id / ≥80% cov) to an
experimentally-verified effector in a database we did *not* use as the gold set: SecReT4 (540) +
SecReT6 (331). Every hit is a true positive. **Floor = 3.2% overall** (gold + 14 novel DB-confirmed).
This is a hard underestimate: it can only confirm effectors already in a DB. It is also only meaningful
for the SS types those DBs cover, **T6SS** reaches **8.3%**; T1/T2/T5SS read ~0 because *no independent
DB covers them*, not because they are low-precision.

**Ceiling, what is not obviously wrong (tier 2).** Annotation buckets over all 1,572:

| bucket | share | meaning |
|---|---|---|
| effector-like | 11.1% | autotransporter, protease, toxin, adhesin, HasA, two-partner... |
| hypothetical | 28.5% | hypothetical / uncharacterized / DUF — genuinely unknown |
| other | 35.2% | a long tail of varied annotations — mostly ambiguous |
| apparatus | 5.7% | the secretion machinery itself emitted as a substrate (TpsB pore, Tss, Hcp/VgrG) |
| **housekeeping** | **19.4%** | transcriptional regulators, chemotaxis, cell division, ribosome, central metabolism → **obvious false positives** |

So **≥19% are clearly wrong** (a ribosomal protein or LysR regulator is not a secreted effector), and
the **soft ceiling is ~75%** (everything not obviously non-secreted and not machinery).

## The headline: a wide band, dominated by the unresolvable middle

Precision of ssign's proximity substrate calls sits between a **~3% provable floor and a ~75% soft
ceiling**, and **~64% (hypothetical + other) is unresolvable by either database or annotation**. Per
type (Figure 01): floor / ceiling = T1SS 6/84, T6SS 8/70, T2SS 2/79, T5bSS 0/81, T4aP 3/59.

This is the precision complement of the recall gap. Recall: ssign emits only 8–10% of testable known
effectors. Precision: of what it does emit, ~1/5 is clearly wrong and only ~1/10 is annotation-
confirmable as an effector. The proximity rule is both **permissive** (large unresolvable middle, ~1/5
clear FPs) and **low-recall**, and the large ambiguous middle is exactly what a learned classifier
could adjudicate where proximity and annotation cannot. That is the case for the secretion-classifier.

## Caveats

- **No negative ground truth.** The ceiling counts "not obviously non-secreted" as potential TP; many
  of the hypothetical/other middle are likely bystanders, so true precision is well below the ceiling.
- **Annotation-based.** Tier 2 inherits genome-annotation errors (e.g. the gold effector NleB is
  mis-annotated "IS3 transposase" by Bakta and would read as housekeeping; gold rows are excluded from
  the FP count for exactly this reason). The keyword buckets are a heuristic, ~35% stays "other".
- **DB floor covers only T4SS/T6SS.** SecReT4/6 give no independent anchor for T1/T2/T5SS.
- **Deferred:** a stratified agent-adjudicated sample (~100 emissions) would place a point estimate
  with a CI inside this band; not run here.
