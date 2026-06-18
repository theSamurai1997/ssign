# enrichment-stats Specification

## Purpose
TBD - created by archiving change enrichment-background-and-plme-default-off. Update Purpose after archive.
## Requirements
### Requirement: Enrichment background sample size
The enrichment test SHALL default to a null sample of 1000 non-neighborhood proteins when estimating the genome background positive rate by sampling.

#### Scenario: Default null size
- **WHEN** an enrichment-enabled run estimates the background by sampling and no null size is overridden
- **THEN** the null sample SHALL contain 1000 proteins (subject to the available non-neighborhood pool size)

#### Scenario: Pool smaller than the default
- **WHEN** the non-neighborhood protein pool has fewer than 1000 members
- **THEN** the null sample SHALL use the entire available pool without error

### Requirement: Exact background when predictions are whole-genome
When per-protein predictions exist for the whole proteome, the enrichment test SHALL estimate the background from ALL non-neighborhood, non-component proteins rather than a random subsample.

#### Scenario: Whole-genome predictor run
- **WHEN** the DLP/DSE predictions were produced for every protein in the genome
- **THEN** the background rate SHALL be computed over the full set of non-neighborhood, non-component proteins
- **AND** no random subsampling SHALL be applied

#### Scenario: Neighborhood-only predictor run
- **WHEN** predictions exist only for the neighborhood and the sampled null proteins
- **THEN** the background SHALL be estimated from the null sample of the configured size

### Requirement: Enrichment predictor set excludes PLM-Effector
The enrichment test SHALL test DLP and DSE only and SHALL NOT include PLM-Effector as an enrichment predictor.

#### Scenario: Enrichment output columns
- **WHEN** the enrichment test runs on any genome
- **THEN** the output SHALL contain rows for tools DLP and DSE only
- **AND** SHALL contain no PLME rows even if a PLM-Effector output file is present

