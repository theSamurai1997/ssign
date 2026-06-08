# Ground-truth substrate lists for gap-validation experiment

Curated by background research agent on 2026-06-08. Six well-characterized
bacterial genomes with literature-validated SS effector inventories.

## Summary

| Genome | SS targets | Catalogued | Validated subset (in TSV) | Reference (most recent) |
|---|---|---|---|---|
| L. pneumophila Phila 1 | T4SS | ~368 | 43 | doi:10.1038/s44320-024-00076-z (2024) |
| C. burnetii RSA 493 | T4SS | ~150 | 41 | doi:10.1146/annurev-micro-020518-115904 (2019) |
| S. Typhimurium LT2 | T3SS-1 + T3SS-2 | ~44 | 44 | doi:10.1016/j.chom.2017.07.009 (2017) |
| Y. pestis CO92 | T3SS Ysc | 7 | 7 | doi:10.1016/j.micpath.2016.02.013 (2016) |
| P. aeruginosa PAO1 | T3SS + T6SSx3 | ~24 | 24 | doi:10.3389/fcimb.2016.00061 (2016) |
| V. cholerae N16961 | T2SS + T6SS | ~28 | 27 | doi:10.1074/jbc.M110.211078 (2011) |
| **Total** |  | **~621** | **186** |  |

## Caveats from the agent

1. **Legionella**: the 368-entry catalog (Karlowicz 2024) collapses
   ortholog mappings from many strains onto Phila 1 lpg numbers. When
   scoring recall against NC_002942, drop any entry whose lpg tag is
   absent from the strain's own GFF.
2. **Coxiella**: ~50 of the literature-cited "effectors" are
   computational only (PmrA-regulon ML, Tier-2 Carey 2011 hits never
   BlaM-confirmed). Tag those as "predicted" rather than "validated".
3. **UniProt accessions marked `-`**: no strain-specific reviewed
   UniProt entry exists. Match by locus tag instead.
4. **Salmonella LT2**: lacks SopE phage. Do not score against SopE.
5. **P. aeruginosa PAO1**: lacks ExoU. Do not score against ExoU.

## Recall scoring (recommended)

For each genome, compute recall at two cuts:
- **Validated core**: the entries in `<genome>.tsv` with
  `evidence=validated`. This is the must-detect set.
- **All-catalogued**: union of validated + predicted from the broader
  literature. Inflates the denominator but matches how some published
  benchmarks are reported.

Report side-by-side in `gap_results.tsv`.
