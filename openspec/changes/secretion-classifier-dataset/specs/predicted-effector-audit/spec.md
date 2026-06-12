## ADDED Requirements

### Requirement: Predicted rows are verified to the gold-set standard
The system SHALL apply the same verification used for the validated gold set to every `predicted`-evidence corpus row: the primary-reference DOI MUST resolve, and the UniProt accession / locus_tag MUST cross-check against the cited record. Rows that fail verification SHALL be excluded with a recorded reason.

#### Scenario: Predicted row with a resolving DOI and consistent identifiers
- **WHEN** a predicted row's DOI resolves and its UniProt/locus_tag are consistent with the cited record
- **THEN** the row is retained and recorded as verified-predicted

#### Scenario: Predicted row with a broken citation or mismatched identifier
- **WHEN** a predicted row's DOI does not resolve or its identifiers contradict the cited record
- **THEN** the row is excluded and the failure reason is logged in a provenance table

### Requirement: Evidence tier is recorded for training weight
The system SHALL tag every retained positive with an `evidence_tier` (validated vs predicted) so downstream training can weight validated examples more heavily than predicted ones. No predicted row SHALL be dropped solely for being predicted once it passes verification.

#### Scenario: Tiered positive set emitted
- **WHEN** the audited validated and predicted rows are combined into the positive set
- **THEN** each row carries an evidence_tier field and validated/predicted counts are reported
