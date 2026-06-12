## ADDED Requirements

### Requirement: Distance measured to the effector's own system instance

For each verified effector, reachability SHALL be measured as the gene-order distance to the nearest machinery component of its own system instance (matched by `sys_instance_id`) on the same replicon. Distance to components of a different instance or a same-type system elsewhere in the genome SHALL NOT count.

#### Scenario: Same-instance component counts
- **WHEN** an effector and a machinery component of its own system instance lie on the same contig
- **THEN** their gene-order separation is a candidate distance for that effector

#### Scenario: Other-instance component ignored
- **WHEN** the only nearby component belongs to a different system instance than the effector
- **THEN** that component does not make the effector reachable

### Requirement: Reachability classification across windows

Each effector SHALL be classified reachable or structurally-impossible at each proximity window N in {3, 5, 7}: reachable if its nearest own-instance component is within N genes, impossible otherwise.

#### Scenario: Reachable at a window
- **WHEN** an effector's nearest own-instance component is 4 genes away
- **THEN** it is impossible at N=3 and reachable at N=5 and N=7

#### Scenario: Impossible at all tested windows
- **WHEN** an effector's nearest own-instance component is 40 genes away
- **THEN** it is impossible at N=3, 5, and 7

### Requirement: Reporting granularity

Results SHALL be reported as the reachable fraction (the ceiling) and its complement (impossible) broken down per SS type and per genome, at each window.

#### Scenario: Per-type ceiling reported
- **WHEN** the analysis completes for a SS type
- **THEN** it reports the reachable fraction at N=3, 5, and 7 for that type

### Requirement: Computed without the system under test

The ceiling analysis SHALL use only the verified gold set, the machinery answer key, and RefSeq gene order. It SHALL NOT invoke ssign or MacSyFinder.

#### Scenario: No pipeline call
- **WHEN** computing the ceiling
- **THEN** no ssign or MacSyFinder process is run and no pipeline output is read

### Requirement: T5SS excluded from the main benchmark

T5SS SHALL be excluded from the headline ceiling/actual numbers, because its product stays cell-surface-attached rather than released and no curated database tracks it, so a proximity ceiling carries no information. The benchmark's reported SS types are T1, T2, T3, T4, and T6. The separate observation that ssign flags predicted secreted proteins near T5aSS clusters SHALL be handled as a clearly-labelled preliminary side study, with no ground-truth claim.

#### Scenario: T5SS absent from headline results
- **WHEN** ceiling and actual results are reported per SS type
- **THEN** T5SS is not among them, and its exclusion and rationale are documented

#### Scenario: T5aSS-neighbor observation reported as preliminary
- **WHEN** the T5aSS-neighbor side study is presented
- **THEN** it is labelled exploratory and preliminary, and makes no recall or ceiling claim

