## ADDED Requirements

### Requirement: T5SS effectors are sourced by subtype under the anti-hallucination contract
The system SHALL source T5SS examples (absent from the corpus) per subtype, each example requiring a verbatim quote from the characterizing paper, a resolvable DOI, and a real locus_tag. Locus_tags SHALL NEVER be invented; an example without a resolvable locus_tag is recorded as unplaceable, not fabricated.

#### Scenario: T5SS example with full provenance
- **WHEN** an agent proposes a T5SS example with a verbatim quote, resolvable DOI, and a locus_tag that exists in the cited genome
- **THEN** the example is accepted into the T5SS set with its subtype and provenance

#### Scenario: T5SS example missing a verifiable locus_tag
- **WHEN** a proposed example cannot be tied to a real locus_tag in a genome
- **THEN** it is recorded as unplaceable with its citation, and excluded from the placed set

### Requirement: T5SS subtype determines the label convention
The system SHALL distinguish T5bSS (two-partner secretion) from T5aSS/T5cSS (autotransporters). T5bSS TpsA proteins SHALL be labelled as normal (protein, instance) substrates. T5aSS/T5cSS autotransporters SHALL be labelled as self-secreted positives with `self_secreted=true`, representing the secreted protein itself rather than a separate substrate.

#### Scenario: Two-partner secretion substrate
- **WHEN** a sourced example is a T5bSS TpsA substrate with its TpsB pore
- **THEN** it is recorded as a (protein, instance) positive with self_secreted=false

#### Scenario: Autotransporter
- **WHEN** a sourced example is a T5aSS or T5cSS autotransporter
- **THEN** it is recorded as a self-secreted positive (self_secreted=true) tied to its own detected system
