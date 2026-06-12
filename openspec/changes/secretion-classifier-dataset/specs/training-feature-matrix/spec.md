## ADDED Requirements

### Requirement: Positive examples join labels to per-protein and pair features
The system SHALL emit one positive row per (effector protein, assigned system instance), joining the per-protein tool signals from ssign run output (DeepLocPro, DeepSecE, SignalP, PLM-Effector, and an ESM embedding reference) with pair-features (gene-distance from the protein to its assigned system's machinery) and system-features (SS type, component count). Instance-unknown positives SHALL be emitted with pair-features null.

#### Scenario: Instance-resolved positive
- **WHEN** an effector is assigned to a specific instance and its genome has ssign run output
- **THEN** a positive row is written with protein features, pair-features (gene-distance), system-features, evidence_tier, and label=1

#### Scenario: Effector protein not located in the run output
- **WHEN** an effector cannot be matched to any protein in the ssign run output
- **THEN** the row is recorded as feature-unavailable with a reason, not silently dropped

### Requirement: A PU candidate/unlabeled set is built with a biological hard-negative filter
The system SHALL construct the unlabeled candidate set from non-effector proteins, restricted to plausibly-secreted candidates (EffectorP-style: extracellular/secretion signal from DeepLocPro, DeepSecE, SignalP, or a TM helix), and SHALL remove candidates with high sequence identity to any positive. The output SHALL be suitable for positive-unlabeled (nnPU) training, not treated as confirmed negatives.

#### Scenario: Candidate passes the secretion filter
- **WHEN** a non-effector protein shows an extracellular/secretion signal and is not highly similar to a positive
- **THEN** it is included in the unlabeled candidate set, labelled unlabeled (not negative)

#### Scenario: Candidate is a likely positive by homology
- **WHEN** a candidate is highly sequence-similar to a known positive
- **THEN** it is excluded from the unlabeled set to avoid mislabeling

### Requirement: Feature assembly is gated on the availability of run output
The system SHALL clearly separate the label-side outputs (auditable without ssign runs) from the feature matrix (which requires the Phase 2 panel ssign runs), and SHALL report which genomes lack run output rather than emitting partial feature rows as if complete.

#### Scenario: Run output missing for a genome
- **WHEN** a positive's genome has no ssign run output yet
- **THEN** its label-side row exists but its feature row is marked pending-run, and the count of pending genomes is reported
