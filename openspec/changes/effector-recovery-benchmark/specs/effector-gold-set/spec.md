## ADDED Requirements

### Requirement: Experimentally-validated scope

The gold set SHALL contain only effectors whose corpus `evidence_level` is `validated` (experimental evidence of secretion). Rows with `evidence_level` of `predicted` SHALL be excluded from the gold set and left untouched in the source corpus.

#### Scenario: Predicted row excluded
- **WHEN** a corpus row has `evidence_level = predicted`
- **THEN** it does not appear in the gold set, and its row in the source corpus is unchanged

#### Scenario: Validated row admitted for verification
- **WHEN** a corpus row has `evidence_level = validated`
- **THEN** it enters the verification pipeline and is eligible for the gold set if it passes

### Requirement: Citation and identifier repair without changing biology

Recoverable `PARTIAL` rows SHALL be repaired by re-citing the correct primary reference and rebuilding the `uniprot` accession from the trustworthy `locus_tag` via the UniProt REST API. Repair SHALL NOT change which protein, organism, or system a row asserts.

#### Scenario: Broken DOI re-cited
- **WHEN** a row's `primary_ref` DOI fails to resolve or resolves to an unrelated paper
- **THEN** the row is re-cited to the correct paper for that protein, and the gene/organism/ss_type fields are left unchanged

#### Scenario: UniProt accession rebuilt from locus_tag
- **WHEN** a row's `uniprot` accession maps to the wrong protein or organism
- **THEN** a new accession is derived from the row's `locus_tag` and recorded, with the locus_tag treated as the source of truth

### Requirement: Removal of non-effector and unverifiable rows

Rows that are biology errors (apparatus components such as Hcp/VgrG/PAAR mislabeled as effectors, immunity or adaptor proteins, eukaryotic-host artifacts, wrong-system or wrong-instance bindings) and rows that remain `FAIL` after review SHALL be removed from the gold set, with the removal reason recorded.

#### Scenario: Apparatus-as-effector dropped
- **WHEN** a row names a structural T6SS component (Hcp, VgrG, or PAAR) without an effector domain
- **THEN** the row is removed from the gold set and the reason is recorded

#### Scenario: Host artifact dropped
- **WHEN** a row's organism cannot host the asserted secretion system (e.g. a eukaryotic-host protein)
- **THEN** the row is removed with the reason recorded

### Requirement: Verification-only use of reference databases

UniProt and RefSeq SHALL be used only to verify and locate proteins (stable ID, RefSeq cross-reference, genome coordinates, annotation cross-check), never as a source to decide that a protein is an effector of a given secretion system.

#### Scenario: No discovery from databases
- **WHEN** building the gold set
- **THEN** no effector-to-system assignment originates from a UniProt or RefSeq query; every such assignment traces to a literature citation

### Requirement: External curated-database reconciliation

The gold set SHALL be reconciled against externally-curated effector databases (SecReT4, SecReT6, SecretEPDB, EffectiveDB, and BastionHub for T1SS/T2SS) to reduce the source corpus's selection bias. Each external entry SHALL be mapped to a UniProt/locus_tag key and deduplicated against the corpus; only net-new entries are added, and each addition SHALL pass the same admission bar (experimental support + a resolvable citation + a locus_tag that resolves in its genome). Where a database is unreachable, a durable fallback (published supplementary tables or an archived snapshot) SHALL be used, and any database that cannot be obtained SHALL be recorded as a coverage gap rather than silently skipped.

#### Scenario: Net-new external effector added
- **WHEN** an external database lists an experimentally-supported effector that is absent from the corpus and resolves to a CDS in a covered genome
- **THEN** it is verified to the standard bar and added to the gold set, tagged with its source database

#### Scenario: Duplicate external entry ignored
- **WHEN** an external entry maps to a UniProt/locus_tag key already present in the corpus
- **THEN** it is not added again, and the existing row is retained

#### Scenario: Unreachable database recorded as a gap
- **WHEN** a database cannot be retrieved live or via fallback
- **THEN** its absence is recorded as an explicit coverage gap in the gold-set provenance, not omitted silently

### Requirement: Final gold set composition

The final gold set SHALL contain only rows that are both `validated` and `VERIFIED`, and every retained row SHALL carry a resolvable primary-reference DOI and a `locus_tag` that resolves to a real CDS in its `refseq_genome`.

#### Scenario: Row admitted
- **WHEN** a row is `validated`, passes verification to `VERIFIED`, has a resolvable DOI, and its locus_tag resolves in the stated genome
- **THEN** it is included in the final gold set

#### Scenario: Row rejected
- **WHEN** a row fails any of those conditions
- **THEN** it is excluded, and the exclusion reason is recorded
