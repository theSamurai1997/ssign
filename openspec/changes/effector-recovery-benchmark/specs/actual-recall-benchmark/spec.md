## ADDED Requirements

### Requirement: Genome panel run through ssign

The benchmark SHALL run ssign on a genome panel chosen data-first from the gold set (genomes carrying enough verified effectors to be informative), with the panel presented to the user before the runs.

#### Scenario: Panel approved before runs
- **WHEN** the gold set and machinery answer key are complete
- **THEN** the proposed genome panel is shown to the user for approval before any ssign run starts

### Requirement: Bakta-to-RefSeq coordinate bridge

Because ssign re-annotates with Bakta and does not preserve original locus_tags, each ssign protein SHALL be bridged back to a RefSeq locus_tag by reciprocal coordinate overlap on the same contig, accepting a match at >=90% of the longer feature length.

#### Scenario: Confident bridge
- **WHEN** a Bakta CDS and a RefSeq CDS on the same contig overlap by >=90% of the longer feature
- **THEN** the Bakta locus_tag is bridged to that RefSeq locus_tag

#### Scenario: No confident bridge
- **WHEN** no RefSeq CDS meets the overlap threshold
- **THEN** the protein is recorded as unbridged rather than force-matched

### Requirement: Actual versus ceiling comparison

For each SS type, the benchmark SHALL count the verified effectors ssign actually emits as secreted and report actual recall against the ceiling. Actual recovery SHALL never exceed the ceiling for the same window.

#### Scenario: Actual within ceiling
- **WHEN** actual recall is computed for a SS type at window N
- **THEN** it is less than or equal to the ceiling for that type at window N

### Requirement: Documented input edge cases

Known input edge cases SHALL be detected and documented rather than silently producing zeros: plasmid-borne effectors require the plasmid to be present in ssign's input, and ssign's whole-genome flags are known to reduce substrate counts.

#### Scenario: Plasmid-borne effectors flagged
- **WHEN** an effector lies on a replicon absent from ssign's input (e.g. a virulence plasmid)
- **THEN** the result marks those effectors as not-in-input rather than counting them as missed predictions

#### Scenario: Whole-genome flag anomaly noted
- **WHEN** a whole-genome-flag run yields fewer substrates than the default run for the same genome
- **THEN** the discrepancy is recorded as a known ssign behavior, not treated as ground-truth recall
