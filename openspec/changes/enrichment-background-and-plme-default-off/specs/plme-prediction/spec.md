## ADDED Requirements

### Requirement: PLM-Effector off by default
ssign SHALL NOT run PLM-Effector by default. PLM-Effector SHALL run only when the user explicitly enables it.

#### Scenario: Default run
- **WHEN** a user runs ssign without specifying any PLM-Effector option
- **THEN** PLM-Effector SHALL NOT execute
- **AND** the pipeline SHALL complete using DLP, DSE, and SignalP predictions

#### Scenario: Explicit opt-in
- **WHEN** a user explicitly enables PLM-Effector
- **THEN** PLM-Effector SHALL execute as part of the run

### Requirement: PLM-Effector remains installable
PLM-Effector SHALL remain installable as part of the extended tier. This change SHALL NOT remove its package, weights, or dependency declarations.

#### Scenario: Extended-tier install
- **WHEN** the extended tier is installed
- **THEN** PLM-Effector and its dependencies SHALL be available for opt-in use and for classifier-feature generation

### Requirement: PLM-Effector positivity confidence gate
When PLM-Effector is used as a binary positive/negative call, a protein SHALL count as PLM-Effector-positive only if its maximum per-type probability is at least 0.8.

#### Scenario: Below the gate
- **WHEN** a protein's PLM-Effector max probability is below 0.8
- **THEN** it SHALL NOT be counted as PLM-Effector-positive, even if a native per-type threshold was passed

#### Scenario: At or above the gate
- **WHEN** a protein's PLM-Effector max probability is at least 0.8
- **THEN** it SHALL be counted as PLM-Effector-positive

#### Scenario: Consistency with other predictors
- **WHEN** PLM-Effector positivity is evaluated
- **THEN** the 0.8 confidence gate SHALL match the DLP and DSE positivity convention
