## ADDED Requirements

### Requirement: Unambiguous predicted effectors are auto-assigned to their instance
The system SHALL assign a verified predicted effector to a specific system instance when its genome contains exactly one detected/curated system of that SS type. No nearest-machinery guessing SHALL be used.

#### Scenario: Genome has a single same-type system
- **WHEN** a predicted effector's genome has exactly one system instance of its SS type
- **THEN** the effector is assigned to that instance and marked instance-resolved (auto)

### Requirement: Ambiguous predicted effectors get a literature-audit resolution attempt
The system SHALL, for predicted effectors in genomes with two or more same-type instances, attempt to resolve the specific instance from the sourcing literature (the primary-reference DOI), recording a verbatim supporting quote when resolved.

#### Scenario: Literature names the specific system
- **WHEN** the sourcing paper attributes the effector to a named/locatable system among the genome's instances
- **THEN** the effector is assigned to that instance with a verbatim quote and DOI recorded (instance-resolved: literature)

#### Scenario: Literature does not disambiguate
- **WHEN** the sourcing paper does not identify which same-type instance the effector belongs to
- **THEN** the effector is retained as an instance-unknown type-level positive, flagged accordingly

### Requirement: Instance-unknown positives are usable but distinguished
The system SHALL emit instance-unknown positives in a way that preserves their SS type and protein features while marking that pair-features requiring a specific instance (e.g. gene-distance to machinery) are unavailable.

#### Scenario: Type-level positive carried into the dataset
- **WHEN** an instance-unknown positive is written to the dataset
- **THEN** it carries its SS type and protein features, with instance-dependent pair-features left null and a type-level flag set
