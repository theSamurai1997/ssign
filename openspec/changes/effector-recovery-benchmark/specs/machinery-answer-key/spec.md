## ADDED Requirements

### Requirement: Pipeline-independent machinery source

Apparatus gene membership SHALL be determined from primary literature only. MacSyFinder, TXSScan, ssign, or any annotation pipeline SHALL NOT be used to decide which genes are machinery, because they are the system under test and would make the benchmark circular.

#### Scenario: MacSyFinder excluded from answer key
- **WHEN** determining which genes are apparatus components of a system instance
- **THEN** no membership decision originates from MacSyFinder or ssign output

### Requirement: Per-instance curation with bounded scope

The answer key SHALL be built one system instance at a time (keyed by `refseq_genome` + `ss_type` + `sys_instance_id`), one curation agent per instance, each agent given a type-specific apparatus brief listing the canonical machinery gene families for that SS type.

#### Scenario: One instance per agent
- **WHEN** curating machinery for a genome that contains multiple system instances
- **THEN** each instance is curated by its own agent task with its own output record

### Requirement: Strict anti-hallucination output schema

Each curation agent SHALL return exactly one status per instance: `COMPLETE`, `PARTIAL`, `REFERENCE_ONLY`, or `NONE_KNOWN`. Every machinery gene reported SHALL carry a verbatim quote from the paper that names it. Agents SHALL record gene names only and SHALL NOT invent or guess locus_tags.

#### Scenario: Gene without quote rejected
- **WHEN** an agent lists a machinery gene with no verbatim supporting quote from a named paper
- **THEN** that gene is not accepted into the answer key

#### Scenario: Reference-only outcome
- **WHEN** an agent locates the paper that characterizes the apparatus but cannot extract the gene set (paywall or unreadable supplementary)
- **THEN** it returns `REFERENCE_ONLY` with the DOI(s) and no fabricated gene list

#### Scenario: No fabricated locus_tags
- **WHEN** a paper names a machinery gene only by gene name
- **THEN** the agent records the gene name and leaves locus_tag resolution to the scripted resolve step

### Requirement: Separated-paper handling

When the effector's corpus citation does not itself characterize the apparatus, the agent SHALL follow citations or search for the founding paper that does, and report that paper's DOI as the source.

#### Scenario: Apparatus defined in a different paper
- **WHEN** the effector's `primary_ref` paper does not list the machinery genes
- **THEN** the agent identifies and cites the paper that does, rather than returning `NONE_KNOWN`

### Requirement: Scripted coordinate resolution

A scripted step SHALL map each paper-named machinery gene to a RefSeq `locus_tag` and genome coordinates within that specific `refseq_genome`, using RefSeq strictly as a coordinate lookup. Gene names that do not resolve SHALL be flagged, not guessed.

#### Scenario: Name resolves to coordinates
- **WHEN** a confirmed gene name matches a CDS in the stated genome
- **THEN** its locus_tag, contig, and start/end coordinates are recorded

#### Scenario: Name fails to resolve
- **WHEN** a confirmed gene name has no matching CDS in the stated genome
- **THEN** it is flagged for manual review and not assigned coordinates

### Requirement: Independent verification pass

A verification pass independent of the curation agent SHALL confirm, per accepted gene, that the cited DOI resolves, that the paper actually names the gene as a machinery component, and that the resolved locus_tag is correct for the genome.

#### Scenario: Verification confirms a gene
- **WHEN** the DOI resolves, the paper names the gene, and the locus_tag checks out
- **THEN** the gene is marked verified in the answer key

#### Scenario: Verification rejects a gene
- **WHEN** any of those three checks fails
- **THEN** the gene is removed or flagged, with the failing check recorded
