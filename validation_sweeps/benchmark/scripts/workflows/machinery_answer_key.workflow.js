export const meta = {
  name: 'machinery-answer-key',
  description: 'Literature-only machinery (apparatus) curation, one agent per secretion-system instance',
  phases: [
    { title: 'T1SS' }, { title: 'T2SS' }, { title: 'T3SS' },
    { title: 'T4SS' }, { title: 'T6SS' },
  ],
}

const INSTANCES = [
  {
    "instance_id": "T3SS_08",
    "ss_type": "T3SS",
    "refseq_genome": "NC_003131.1",
    "organism": "Yersinia pestis CO92",
    "sys_instance_label": "Ysc-Yop",
    "n_effectors": "7",
    "effector_loci": "YPCD1.06,YPCD1.19c,YPCD1.20,YPCD1.26c,YPCD1.67c,YPCD1.71c,YPCD1.72c"
  },
  {
    "instance_id": "T3SS_21",
    "ss_type": "T3SS",
    "refseq_genome": "NC_006351.1",
    "organism": "Burkholderia pseudomallei K96243",
    "sys_instance_label": "Bsa",
    "n_effectors": "9",
    "effector_loci": "BPSS1516,BPSS1524,BPSS1525,BPSS1526,BPSS1527,BPSS1528,BPSS1529,BPSS1531,BPSS1532"
  },
  {
    "instance_id": "T3SS_27",
    "ss_type": "T3SS",
    "refseq_genome": "NC_011601.1",
    "organism": "Escherichia coli O127:H6 E2348/69 (EPEC)",
    "sys_instance_label": "prophage_PP6",
    "n_effectors": "3",
    "effector_loci": "E2348C_1707,E2348C_1710,E2348C_1721"
  },
  {
    "instance_id": "T3SS_29",
    "ss_type": "T3SS",
    "refseq_genome": "NC_013716.1",
    "organism": "Citrobacter rodentium ICC168",
    "sys_instance_label": "prophage_CRP1",
    "n_effectors": "2",
    "effector_loci": "ROD_12141,ROD_12161"
  },
  {
    "instance_id": "T3SS_30",
    "ss_type": "T3SS",
    "refseq_genome": "NC_013716.1",
    "organism": "Citrobacter rodentium ICC168",
    "sys_instance_label": "prophage_CRP2",
    "n_effectors": "3",
    "effector_loci": "ROD_06241,ROD_06251,ROD_07221"
  },
  {
    "instance_id": "T3SS_31",
    "ss_type": "T3SS",
    "refseq_genome": "NC_013716.1",
    "organism": "Citrobacter rodentium ICC168",
    "sys_instance_label": "prophage_CRP3",
    "n_effectors": "3",
    "effector_loci": "ROD_29641,ROD_29651,ROD_29661"
  },
  {
    "instance_id": "T3SS_32",
    "ss_type": "T3SS",
    "refseq_genome": "NC_013716.1",
    "organism": "Citrobacter rodentium ICC168",
    "sys_instance_label": "prophage_CRP4",
    "n_effectors": "1",
    "effector_loci": "ROD_25861"
  },
  {
    "instance_id": "T4SS_01",
    "ss_type": "T4SS",
    "refseq_genome": "NC_000915.1",
    "organism": "Helicobacter pylori 26695",
    "sys_instance_label": "Cag",
    "n_effectors": "1",
    "effector_loci": "HP_0547"
  },
  {
    "instance_id": "T4SS_02",
    "ss_type": "T4SS",
    "refseq_genome": "NC_002942.5",
    "organism": "Legionella pneumophila Philadelphia 1",
    "sys_instance_label": "DotIcm",
    "n_effectors": "56",
    "effector_loci": "lpg0130,lpg0160,lpg0246,lpg0376,lpg0390,lpg0437,lpg0589,lpg0621,lpg0696,lpg0758,lpg0945,lpg0962,lpg1368,lpg1488,lpg1621,lpg1683,lpg1701,lpg1851,lpg1888,lpg1924,lpg1950,lpg1959,lpg1962,lpg2098,lpg2144,lpg2147,lpg2148,lpg2155,lpg2240,lpg2248,lpg2311,lpg2391,lpg2410,lpg2422,lpg2452,lpg2455,lpg2456,lpg2465,lpg2504,lpg2510,lpg2522,lpg2526,lpg2529,lpg2552,lpg2580,lpg2603,lpg2718,lpg2813,lpg2815,lpg2818,lpg2829,lpg2831,lpg2862,lpg2887,lpg2965,lpg3022"
  },
  {
    "instance_id": "T4SS_03",
    "ss_type": "T4SS",
    "refseq_genome": "NC_002971.4",
    "organism": "Coxiella burnetii RSA493",
    "sys_instance_label": "DotIcm",
    "n_effectors": "4",
    "effector_loci": "CBU0041,CBU0626,CBU0937,CBU1217"
  },
  {
    "instance_id": "T4SS_04",
    "ss_type": "T4SS",
    "refseq_genome": "NC_002978.6",
    "organism": "Wolbachia pipientis wMel",
    "sys_instance_label": "VirB",
    "n_effectors": "3",
    "effector_loci": "E0495_02270,E0495_02275,E0495_06075"
  },
  {
    "instance_id": "T4SS_06",
    "ss_type": "T4SS",
    "refseq_genome": "NC_003317.1",
    "organism": "Brucella melitensis 16M",
    "sys_instance_label": "VirB",
    "n_effectors": "1",
    "effector_loci": "BMEI1674"
  },
  {
    "instance_id": "T4SS_07",
    "ss_type": "T4SS",
    "refseq_genome": "NC_005956.1",
    "organism": "Bartonella henselae Houston-1",
    "sys_instance_label": "VirB-Bep",
    "n_effectors": "7",
    "effector_loci": "BH13370,BH13390,BH13400,BH13410,BH13420,BH13430,BH13440"
  },
  {
    "instance_id": "T4SS_08",
    "ss_type": "T4SS",
    "refseq_genome": "NC_007618.1",
    "organism": "Brucella abortus 2308",
    "sys_instance_label": "VirB",
    "n_effectors": "13",
    "effector_loci": "BAB1_0123,BAB1_0244,BAB1_0275,BAB1_0279,BAB1_0466,BAB1_0782,BAB1_1043,BAB1_1090,BAB1_1232,BAB1_1279,BAB1_1305,BAB1_1325,BAB1_1492"
  },
  {
    "instance_id": "T4SS_09",
    "ss_type": "T4SS",
    "refseq_genome": "NC_007624.1",
    "organism": "Brucella abortus 2308",
    "sys_instance_label": "VirB",
    "n_effectors": "3",
    "effector_loci": "BAB2_0123,BAB2_0466,BAB2_0541"
  },
  {
    "instance_id": "T4SS_10",
    "ss_type": "T4SS",
    "refseq_genome": "NC_007797.1",
    "organism": "Anaplasma phagocytophilum HZ",
    "sys_instance_label": "VirB",
    "n_effectors": "4",
    "effector_loci": "APH_0233,APH_0740,APH_0859,APH_1387"
  },
  {
    "instance_id": "T4SS_11",
    "ss_type": "T4SS",
    "refseq_genome": "NC_007799.1",
    "organism": "Ehrlichia chaffeensis Arkansas",
    "sys_instance_label": "VirB",
    "n_effectors": "4",
    "effector_loci": "ECH_0039,ECH_0822,ECH_0825,ECH_1027"
  },
  {
    "instance_id": "T4SS_12",
    "ss_type": "T4SS",
    "refseq_genome": "NC_010263.3",
    "organism": "Rickettsia rickettsii Iowa",
    "sys_instance_label": "VirB",
    "n_effectors": "1",
    "effector_loci": "RrIowa_0604"
  },
  {
    "instance_id": "T4SS_13",
    "ss_type": "T4SS",
    "refseq_genome": "NZ_CP019789.1",
    "organism": "Bartonella schoenbuchensis R1",
    "sys_instance_label": "VbhB-VbhT",
    "n_effectors": "1",
    "effector_loci": "B11C_100026"
  },
  {
    "instance_id": "T6SS_01",
    "ss_type": "T6SS",
    "refseq_genome": "NC_002505.1",
    "organism": "Vibrio cholerae N16961",
    "sys_instance_label": "Aux1",
    "n_effectors": "1",
    "effector_loci": "VC_1418"
  },
  {
    "instance_id": "T6SS_02",
    "ss_type": "T6SS",
    "refseq_genome": "NC_002506.1",
    "organism": "Vibrio cholerae N16961",
    "sys_instance_label": "Aux2",
    "n_effectors": "1",
    "effector_loci": "VC_A0020"
  },
  {
    "instance_id": "T6SS_03",
    "ss_type": "T6SS",
    "refseq_genome": "NC_002506.1",
    "organism": "Vibrio cholerae N16961",
    "sys_instance_label": "Aux3",
    "n_effectors": "1",
    "effector_loci": "VC_A0285"
  },
  {
    "instance_id": "T6SS_04",
    "ss_type": "T6SS",
    "refseq_genome": "NC_002506.1",
    "organism": "Vibrio cholerae N16961",
    "sys_instance_label": "Large-T6SS",
    "n_effectors": "2",
    "effector_loci": "VC_A0117,VC_A0123"
  },
  {
    "instance_id": "T6SS_06",
    "ss_type": "T6SS",
    "refseq_genome": "NC_002516.2",
    "organism": "Pseudomonas aeruginosa PAO1",
    "sys_instance_label": "H2-T6SS",
    "n_effectors": "8",
    "effector_loci": "PA0259,PA0260,PA0262,PA1510,PA1510,PA3290,PA3487,PA3905"
  },
  {
    "instance_id": "T6SS_07",
    "ss_type": "T6SS",
    "refseq_genome": "NC_002516.2",
    "organism": "Pseudomonas aeruginosa PAO1",
    "sys_instance_label": "H3-T6SS",
    "n_effectors": "3",
    "effector_loci": "PA2374,PA4922,PA5089"
  },
  {
    "instance_id": "T6SS_08",
    "ss_type": "T6SS",
    "refseq_genome": "NC_003063.2",
    "organism": "Agrobacterium tumefaciens C58",
    "sys_instance_label": "T6SS",
    "n_effectors": "3",
    "effector_loci": "Atu3640,Atu4347,Atu4350"
  },
  {
    "instance_id": "T6SS_09",
    "ss_type": "T6SS",
    "refseq_genome": "NC_007650.1",
    "organism": "Burkholderia thailandensis E264",
    "sys_instance_label": "T6SS-5",
    "n_effectors": "5",
    "effector_loci": "BTH_II0865,BTH_II0866,BTH_II0874,BTH_II0875,BTH_II0876"
  },
  {
    "instance_id": "T6SS_10",
    "ss_type": "T6SS",
    "refseq_genome": "NC_007651.1",
    "organism": "Burkholderia thailandensis E264",
    "sys_instance_label": "T6SS-1",
    "n_effectors": "5",
    "effector_loci": "BTH_I0070,BTH_I2691,BTH_I2692,BTH_I2693,BTH_I2694"
  },
  {
    "instance_id": "T6SS_11",
    "ss_type": "T6SS",
    "refseq_genome": "NC_009085.1",
    "organism": "Acinetobacter baumannii ATCC 17978",
    "sys_instance_label": "T6SS",
    "n_effectors": "5",
    "effector_loci": "A1S_1296,ACX60_00605,ACX60_11695,ACX60_15365,ACX60_17660"
  },
  {
    "instance_id": "T6SS_12",
    "ss_type": "T6SS",
    "refseq_genome": "NC_010465.1",
    "organism": "Yersinia pseudotuberculosis YPIII",
    "sys_instance_label": "T6SS-3",
    "n_effectors": "1",
    "effector_loci": "YPK_0952"
  },
  {
    "instance_id": "T6SS_13",
    "ss_type": "T6SS",
    "refseq_genome": "NC_010465.1",
    "organism": "Yersinia pseudotuberculosis YPIII",
    "sys_instance_label": "T6SS-4",
    "n_effectors": "2",
    "effector_loci": "YPK_3548,YPK_3549"
  },
  {
    "instance_id": "T6SS_14",
    "ss_type": "T6SS",
    "refseq_genome": "NC_011000.1",
    "organism": "Burkholderia cenocepacia J2315",
    "sys_instance_label": "T6SS-6_BC",
    "n_effectors": "1",
    "effector_loci": "BCAL0824"
  },
  {
    "instance_id": "T6SS_15",
    "ss_type": "T6SS",
    "refseq_genome": "NC_013508.1",
    "organism": "Edwardsiella piscicida EIB202",
    "sys_instance_label": "EVP",
    "n_effectors": "4",
    "effector_loci": "ETAE_0972,ETAE_0975,ETAE_0976,ETAE_0977"
  },
  {
    "instance_id": "T6SS_16",
    "ss_type": "T6SS",
    "refseq_genome": "NC_013508.1",
    "organism": "Edwardsiella piscicida EIB202",
    "sys_instance_label": "MGE",
    "n_effectors": "1",
    "effector_loci": "ETAE_2037"
  },
  {
    "instance_id": "T6SS_17",
    "ss_type": "T6SS",
    "refseq_genome": "NC_016856.1",
    "organism": "Salmonella enterica Typhimurium 14028",
    "sys_instance_label": "SPI-6",
    "n_effectors": "3",
    "effector_loci": "STM14_0334,STM14_0336,STM14_0337"
  },
  {
    "instance_id": "T6SS_18",
    "ss_type": "T6SS",
    "refseq_genome": "NC_017626.1",
    "organism": "Escherichia coli 042",
    "sys_instance_label": "Sci-1",
    "n_effectors": "2",
    "effector_loci": "EC042_4528,EC042_4533"
  },
  {
    "instance_id": "T6SS_19",
    "ss_type": "T6SS",
    "refseq_genome": "NC_017626.1",
    "organism": "Escherichia coli 042",
    "sys_instance_label": "Sci-2",
    "n_effectors": "2",
    "effector_loci": "EC042_4660,EC042_4661"
  },
  {
    "instance_id": "T6SS_20",
    "ss_type": "T6SS",
    "refseq_genome": "NZ_CP012003",
    "organism": "Acinetobacter baumannii AB307-0294",
    "sys_instance_label": "T6SS",
    "n_effectors": "1",
    "effector_loci": "A1S_1292"
  },
  {
    "instance_id": "T6SS_21",
    "ss_type": "T6SS",
    "refseq_genome": "NZ_HG326223.1",
    "organism": "Serratia marcescens Db10",
    "sys_instance_label": "T6SS",
    "n_effectors": "10",
    "effector_loci": "-,-,SMDB11_2270,SMDB11_2278,SMDB11_2279,SMDB11_2280,SMDB11_2342,SMDB11_2343,SMDB11_2369,SMDB11_4259"
  }
]

function buildPrompt(it) {
  return `You are curating the MACHINERY (apparatus) gene list for ONE secretion-system instance, from the PUBLISHED LITERATURE ONLY. This is the ground-truth answer key for a benchmark that tests a proximity-based effector predictor, so it must be ssign-independent and free of hallucination.

## The instance
- instance_id: ${it.instance_id}
- ss_type: ${it.ss_type}
- organism: ${it.organism}
- refseq_genome (for your reference only, DO NOT mine it): ${it.refseq_genome}
- sys_instance_label: ${it.sys_instance_label || '(none given)'}
- this instance's curated effectors (context, NOT machinery): ${it.effector_loci || '(none)'} (${it.n_effectors} effectors)

## Step 1 - read your type brief
Read the file scripts/apparatus_briefs/${it.ss_type}.md . It lists the canonical machinery gene families and the system-specific synonyms a paper may use. Use it as your recognition vocabulary. Also read data/machinery_raw/T2SS_01.json to see the exact output format and quality bar (a COMPLETE pilot).

## Hard rules (anti-circularity, anti-hallucination)
1. LITERATURE ONLY. Find the apparatus genes from published papers (use the PubMed MCP tools - search_articles, get_article_metadata, get_full_text_article, convert_article_ids - and WebSearch / WebFetch for full text). Do NOT use MacSyFinder, Pfam/HMM scans, or any automated detector: MacSyFinder is the tool under test, so using it is circular. Do NOT open or mine the RefSeq GenBank to discover genes; coordinate/locus resolution is a separate scripted step done later.
2. Record gene NAMES ONLY (e.g. "VirB4", "HrcC", "EpsD", "ClpV", "VgrG1"). NEVER invent or guess a locus_tag. If a paper gives a locus_tag you may put it in the note, but the gene field is the name.
3. Every machinery gene needs a VERBATIM quote from a paper that names that gene as a component of THIS system in (or transferable to) this organism, plus a resolvable DOI and the PMID. The quote must actually contain the gene name (or unambiguously the component). If you cannot find a naming quote, do not list the gene.
4. Separated-paper rule: the apparatus is often described in a different paper than the effector. Follow citations to the founding/characterisation paper for the machinery. The apparatus paper may be on a model strain of the same system family (e.g. the archetypal VirB/D4, Ysc, Sct system); that is acceptable as long as the paper clearly names the components of this system type and you note any strain transfer.
5. ${it.ss_type === 'T6SS' ? 'T6SS SPECIAL: Hcp, VgrG and PAAR ARE machinery here (secreted tube/spike/tip) and belong in the answer key, even though they are excluded from the effector set. Count ClpV (TssH). Do NOT list delivered toxin effectors (Tse/Tae/Tle/Tde/Tme/Rhs-CT) or immunity proteins.' : 'Do NOT list the secreted substrate/effector, escort chaperones, or transcriptional regulators as machinery.'}
6. Do NOT pad. Only genes you can defend with a quote. Quality over completeness.

## Status (pick exactly one)
- COMPLETE: the core apparatus of this system is documented gene-by-gene in the literature with quotes.
- PARTIAL: some core components documented with quotes, but clear gaps remain.
- REFERENCE_ONLY: papers describe/operate this system but do not name the individual apparatus genes (so you can cite the system but list few/no genes).
- NONE_KNOWN: no literature found that names apparatus genes for this instance.

## Output
Write a file at data/machinery_raw/${it.instance_id}.json with EXACTLY this schema (valid JSON, no trailing commas):
{
  "instance_id": "${it.instance_id}",
  "ss_type": "${it.ss_type}",
  "organism": "${it.organism}",
  "refseq_genome": "${it.refseq_genome}",
  "sys_instance_label": "${it.sys_instance_label}",
  "status": "COMPLETE|PARTIAL|REFERENCE_ONLY|NONE_KNOWN",
  "machinery_genes": [
    {"gene":"", "family":"", "role":"", "doi":"", "pmid":"", "quote":"", "note":""}
  ],
  "papers": [
    {"doi":"", "pmid":"", "first_author":"", "year":"", "role":"apparatus|review|effector"}
  ],
  "notes": "what the system is, which strain the apparatus paper used if transferred, what you deliberately excluded, and any low-confidence calls"
}

Verify each DOI resolves and each PMID is real before writing (a get_article_metadata call on the PMID confirms both). After writing the file, return ONE line: "${it.instance_id} <status> <N machinery_genes>". Your returned text is data, not a message to a human.`
}

const results = await parallel(
  INSTANCES.map((it) => () =>
    agent(buildPrompt(it), { label: it.instance_id, phase: it.ss_type })
  )
)

const summary = results.map((r, i) => `${INSTANCES[i].instance_id}: ${r ? String(r).trim().slice(0, 120) : 'NULL (agent failed)'}`)
log(`Curation finished. ${results.filter(Boolean).length}/${INSTANCES.length} returned.`)
return { count: INSTANCES.length, returned: results.filter(Boolean).length, summary }
