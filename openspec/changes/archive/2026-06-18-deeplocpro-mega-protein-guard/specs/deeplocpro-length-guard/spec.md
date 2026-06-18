## ADDED Requirements

### Requirement: Over-length proteins do not crash DeepLocPro
The DeepLocPro step SHALL NOT pass sequences longer than a configured maximum to the DeepLocPro model, so that a single over-length protein cannot fail the step (and the run).

#### Scenario: A genome contains a mega-protein
- **WHEN** the input contains a protein longer than the configured maximum length (default 5000 aa)
- **THEN** that protein SHALL be withheld from the DeepLocPro model invocation
- **AND** DeepLocPro SHALL run to completion on the remaining sequences
- **AND** the DeepLocPro step SHALL be reported as succeeded

#### Scenario: All proteins within the limit
- **WHEN** every input protein is at or below the maximum length
- **THEN** all proteins SHALL be passed to DeepLocPro unchanged

### Requirement: Skipped proteins are surfaced, not silently dropped
Each protein withheld for exceeding the maximum length SHALL appear in the DeepLocPro output marked as not predicted, and SHALL be logged.

#### Scenario: Output row for a skipped protein
- **WHEN** a protein is withheld for being over-length
- **THEN** the output SHALL contain a row for that protein with a not-predicted localization and zero localization probabilities
- **AND** a warning SHALL be logged identifying the skipped protein(s) and their length(s)

#### Scenario: Downstream treats a skipped protein as non-secreted
- **WHEN** a skipped protein has no extracellular probability in the DeepLocPro output
- **THEN** cross-validation SHALL treat it as not extracellular (non-secreted by DeepLocPro), without error

### Requirement: Configurable maximum length
The maximum sequence length SHALL be a single configurable value with a default of 5000 aa.

#### Scenario: Override the default
- **WHEN** the maximum length is overridden (constant or environment)
- **THEN** the partition SHALL use the overridden value
