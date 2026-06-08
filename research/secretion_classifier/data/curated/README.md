# Curated instance-labeled effector tables

Five per-SS-type TSVs of literature-curated effector inventories,
produced by background research agents on 2026-06-08. Each row is one
known substrate of a specific secretion-system instance, with
organism + locus tag + sys_instance_id where applicable.

## Schema

`gene  uniprot  locus_tag  organism  refseq_genome  ss_type  sys_instance_id  evidence_level  primary_ref`

(T1SS adds three extra columns: `family`, `length`, `proteome` — useful
for downstream filtering.)

- `evidence_level`: "validated" (experimental translocation demonstrated
  in vitro or in vivo, e.g. BlaM / CyaA / TEM reporter assay) or
  "predicted" (TrEMBL-annotated homolog of a validated effector, or
  bioinformatic-only assignment without translocation confirmation).
- `sys_instance_id`: cluster identifier where multiple SS instances of
  the same type exist in one genome (e.g. PAO1 T6SS H1/H2/H3, EHEC
  Sakai LEE + 14 prophage IDs). Blank when single-instance (most
  T2SS / single-cluster T4SS organisms).
- `uniprot`: `-` if no reviewed UniProt entry exists for the gene;
  match by locus_tag against the RefSeq genome instead.

## Counts (2026-06-08)

| File | Total | Validated | Predicted | Notes |
|---|---|---|---|---|
| `T1SS_curated.tsv` | 212 | 26 | 186 | Predicted are TrEMBL homologs of validated RTX-toxin / metalloprotease / adhesin families |
| `T2SS_curated.tsv` | 101 | 91 | 10 | V. cholerae Eps (23) + Legionella Lsp (19) + Dickeya Out (19) + PAO1 Xcp (15) anchor the set |
| `T3SS_curated.tsv` | 355 | 237 | 118 | Multi-cluster splits captured: EHEC Sakai LEE + 14 prophages, Salmonella SPI-1/SPI-2, Yersinia Ysc + Ysa, V. parahaemolyticus T3SS-1/T3SS-2 |
| `T4SS_curated.tsv` | 131 | 114 | 17 | ~104 net-new beyond SecReT4. Legionella post-2019 Karlowicz atlas + Brucella + Bartonella VirB-Bep + Anaplasma/Ehrlichia/Wolbachia/Rickettsia |
| `T6SS_curated.tsv` | 150 | 119 | 16 | PAO1 H1/H2/H3 + Burkholderia thailandensis T6SS-1 to T6SS-5 + V. cholerae aux1/2/3 + Serratia + Salmonella SPI-6 |
| **Total** | **949** | **587** | **347** | |

Combined with the pre-existing SecReT4-derived 504 T4 cluster labels
and the validation-genome ground truth (186 entries in
`../ground_truth/`), total instance-labeled positives available for
Tier-2 training: **~1,400-1,500** (after dedup for overlaps between
curated/SecReT4/ground_truth on T4 and T6).

## Provenance

Each agent mined a mix of UniProt REST API queries, published review
papers with machine-readable supplementary tables, per-organism
specialized databases (pathogens3d.org for Legionella,
pseudomonas-syringae.org for P. syringae Hops, etc.), VFDB, and
original characterization papers for locus tag mapping. See agent
inline summaries in task records for the exact source list per SS type.

Primary refs in the TSVs are DOIs of the original characterization
paper.

## Dedup before training

These tables are NOT cross-deduplicated against:

- The pre-existing 504 SecReT4-derived T4 cluster labels (de-anon
  recovery from PLM-E + DSE)
- The 50 SecReT6 coordinate-joined T6 entries
- The `../ground_truth/*.tsv` for the 6 validation genomes

Before training, run an exact-sequence dedup pass at 100% identity
across the union to catch literal overlaps, then a 30%-identity
MMseqs2 clustering for held-out split integrity. See
`02_data_audit.md` for the full pipeline recipe.

## Caveats

- T1SS validated set is small (26) because experimental translocation
  has only been confirmed for canonical lab strains. The 186 predicted
  entries are mostly TrEMBL homologs — introduce some label noise.
- T6SS multi-cluster coverage in Burkholderia thailandensis is mostly
  T6SS-1 and T6SS-5; T6SS-2/3/4 are listed as apparatus-only because
  effectors aren't characterized in the literature.
- Some `locus_tag` cells need a follow-up tBLASTn pass against the
  source RefSeq genome (flagged `-`) before training.
