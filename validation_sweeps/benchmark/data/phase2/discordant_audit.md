# Discordant-effector literature audit (2026-06-12)

The 19 effectors that ssign **emitted but the answer-key machinery sits >7 genes away**
(`actual_per_effector.panel_genbank_t3ss.tsv`, `emitted_secreted` & `testable` & `reachable_n7 != true`).
One literature agent per effector-sourcing DOI re-read the paper under an anti-hallucination
contract (verbatim quote + resolvable DOI required). BipB/BipC agents were blocked by the AUP
filter on *B. pseudomallei* content and were checked manually.

## Headline: the answer key has a citation-integrity problem

Most of these rows are **not** biological misassignments, the effector→SS-type call is usually
right, but the **sourcing DOIs are frequently wrong or nonexistent**, and a few rows are
genuinely mislabelled or unidentifiable. This was surfaced by auditing only the 19 discordant
rows; the same defects likely exist in the wider corpus and warrant a full citation pass.

| # | effector | uniprot | our type/inst | DOI status | verdict | ssign emitted via |
|---|---|---|---|---|---|---|
| 1 | celA | Q5ZU89 | T2SS_05 | **WRONG** (resolves to a soil-ecology 16S paper) | biology plausibly T2SS (Cianciotto lit), but cite is wrong | T4aP (coincidental) |
| 2 | plaA | Q5ZX07 | T2SS_05 | **WRONG** (same soil paper) | as celA | T4aP (coincidental) |
| 3 | VirA | Q7BU69 | T3SS_20 | **WRONG** (resolves to a mouse tooth-development paper) | VirA IS a real T3SS effector (ext. confirmed) | T5aSS (coincidental) |
| 4 | ChlaDub1 | - | T3SS_01 | **404 / unregistered** | **questionable**: Cdu1 is a secreted inclusion-membrane DUB; no source ties it to T3SS | T5aSS (likely a Pmp autotransporter; coincidental) |
| 5 | CopN | - | T3SS_01 | **404 / unregistered** | CopN IS a real T3SS substrate (plug/gatekeeper); fix DOI (likely jbc.M115.670232) | T3SS (type-concordant) |
| 6 | BipB | A0A8A4DZ70 | T3SS_21 | not re-checked (agent blocked) | established Bsa T3SS translocon; assignment sound | T3SS (type-concordant) |
| 7 | BipC | A0A8A4DTS9 | T3SS_21 | not re-checked (agent blocked) | established Bsa T3SS translocon/effector; assignment sound | T3SS (type-concordant) |
| 8 | BopA | - | T6SS_09 | resolves | **MISASSIGNED**: BopA is a Bsa **T3SS** effector; 2024 paper does not even contain "BopA" | T5cSS (coincidental) |
| 9 | BopE | Q2T6X8 | T6SS_09 | resolves | **MISASSIGNED**: BopE is a Bsa **T3SS** effector (Stevens 2003); 2024 paper uses it only as a control panel | T5cSS (coincidental) |
| 10 | TseM | - | T6SS_09 | resolves | CONFIRMED T6SS (paper: **T6SS-4**, BTH_II1883); anchor-distance artifact | T6SS (concordant) |
| 11 | TseZ | Q2T421 | T6SS_09 | resolves | CONFIRMED T6SS (**T6SS4**); anchor-distance artifact | T6SS (concordant) |
| 12 | Tle4 | - | T6SS_06 | **WRONG** (resolves to a Baker protein-design paper) | CONFIRMED H2-T6SS; **duplicate of #13** | T6SS (concordant) |
| 13 | TplE_alias_Tle4 | Q9I3K2 | T6SS_06 | **WRONG** (resolves to a lung-regeneration paper) | CONFIRMED H2-T6SS (=PA1510); **dup of #12**, de-dupe | T6SS (concordant) |
| 14 | EFF00142 | - | T6SS_10 | resolves | **UNIDENTIFIABLE** DB-placeholder id; not in cited paper | T6SS (concordant) |
| 15 | EFF00150 | - | T6SS_10 | resolves | **UNIDENTIFIABLE** DB-placeholder id; not in cited paper | T6SS (concordant) |
| 16 | Tle1 | - | T6SS_10 | **MIS-CITED** (cited 2012 amidase paper has no Tle1; Tle1 is from 2013 Nature nature12074) | biology real (T6SS phospholipase), cite wrong | T6SS (concordant) |
| 17 | TseA_T6SS1 | Q2SV43 | T6SS_10 | resolves | **UNSUPPORTED**: cited paper never names it; UniProt = VasX-domain BTH_I2691; no substrate characterized | T6SS (concordant) |
| 18 | Tae4_Stm | - | T6SS_17 | **WRONG** (resolves to a Drosophila-immunity paper; correct = chom.2012.04.007) | CONFIRMED T6SS (SPI-6 amidase effector) | T6SS (concordant) |
| 19 | Tlde1A | - | T6SS_17 | resolves | CONFIRMED T6SS (SPI-6 L,D-transpeptidase effector) | T6SS (concordant) |

## Tally

- **Genuine biological misassignment (wrong SS type):** BopA, BopE (both T6SS→Bsa T3SS, "Bop"
  name-collision). ChlaDub1 questionable (no T3SS evidence; it's an inclusion-membrane DUB).
- **Unidentifiable / unsupported rows (drop or repair):** EFF00142, EFF00150, TseA_T6SS1.
- **Duplicate row:** Tle4 == TplE_alias_Tle4 (same protein PA1510/Q9I3K2).
- **Wrong or nonexistent sourcing DOI (≥9):** celA, plaA, VirA, ChlaDub1(404), CopN(404), Tle4,
  TplE, Tle1, Tae4_Stm. (Effector biology often still correct, the metadata is wrong.)
- **Sound assignment, discordance is only answer-key anchor distance:** TseM, TseZ, Tlde1A, CopN,
  BipB, BipC (and the confirmed T6SS rows above).

## Implications

1. **Recall numbers are unaffected** (recall is measured against machinery gene-position, not the
   effector's sourcing DOI). The headline 8%/10% stands.
2. **ssign's emission "for the wrong reason":** celA, plaA, VirA, ChlaDub1, BopA, BopE were emitted
   only because an unrelated nearby system (T4aP/T5aSS/T5cSS) sat within ±3 genes, not their own
   system. This is the cross-type-adjacency cost quantified in `found_systems.*.tsv`.
3. **Dataset action (deferred):** run the citation-integrity audit across the full corpus (not just
   these 19); drop/repair the 3 unidentifiable rows; reassign BopA/BopE to T3SS; de-dupe Tle4/TplE;
   re-source the ≥9 wrong DOIs. Tracked in NOTES.md.
