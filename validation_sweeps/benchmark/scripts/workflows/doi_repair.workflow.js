export const meta = {
  name: 'doi-repair',
  description: 'Repair broken effector citation DOIs via literature agents with verbatim-quote verification',
  phases: [{ title: 'Repair', detail: 'one literature agent per broken-DOI effector row' }],
}

// 54 PARTIAL effector rows whose stored primary-reference DOI is broken/wrong.
// Generated from data/gold_build/04_repair_residual.tsv (doi_status=BROKEN_RESIDUAL),
// minus the 5 pilot rows already done.
const ROWS = [
{
"gene": "apxIIA",
"uniprot": "P15377",
"organism": "Actinobacillus pleuropneumoniae (Haemophilus pleuropneumoniae)",
"ss_type": "T1SS",
"locus_tag": "",
"refseq_genome": "",
"doi_old": "10.1128/iai.61.5.2056-2063.1993",
"notes": "(a) BAD DOI: DOI returns 404 on Crossref and journal site; not findable in PubMed. Frey 1993 Apx paper is in Gene 123:51 (10.1016/0378-1119(93)90539-e) or Mol Microbiol.; (b) family (ApxIIA (A. pleuropneumoniae RTX-II)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate ide"
},
{
"gene": "apxIIIA",
"uniprot": "P55131",
"organism": "Actinobacillus pleuropneumoniae (Haemophilus pleuropneumoniae)",
"ss_type": "T1SS",
"locus_tag": "",
"refseq_genome": "",
"doi_old": "10.1128/iai.61.5.2056-2063.1993",
"notes": "(a) BAD DOI: DOI returns 404 on Crossref and journal site; not findable in PubMed. Frey 1993 Apx paper is in Gene 123:51 (10.1016/0378-1119(93)90539-e) or Mol Microbiol.; (b) family (ApxIIIA (A. pleuropneumoniae RTX-III)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate i"
},
{
"gene": "frpC",
"uniprot": "Q9JYV5",
"organism": "Neisseria meningitidis serogroup B (strain ATCC BAA-335 / MC58)",
"ss_type": "T1SS",
"locus_tag": "NMB1415",
"refseq_genome": "AE002098",
"doi_old": "10.3390/toxins9100327",
"notes": "(a) WRONG DOI: DOI resolves to Menezes 2017 cyanobacterial toxins in Portuguese recreational waters, NOT FrpC. Wrong reference.; (b) family (FrpC (Neisseria iron-regulated RTX)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate identity is independently supported; (c) locu"
},
{
"gene": "hasA",
"uniprot": "Q54450",
"organism": "Serratia marcescens",
"ss_type": "T1SS",
"locus_tag": "",
"refseq_genome": "",
"doi_old": "10.1002/j.1460-2075.1994.tb06846.x",
"notes": "(a) WRONG DOI: DOI resolves to 'Note added' placeholder, not a research paper. Letoffe HasA 1994 EMBO J 13:5804 has a different DOI.; (b) family (HasA (heme-acquisition hemophore)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate identity is independently supported; (c) l"
},
{
"gene": "ltxA",
"uniprot": "P16462",
"organism": "Aggregatibacter actinomycetemcomitans (Actinobacillus actinomycetemcomitans) (Haemophilus actinomycetemcomitans)",
"ss_type": "T1SS",
"locus_tag": "",
"refseq_genome": "",
"doi_old": "10.1016/S0021-9258(19)84909-X",
"notes": "(a) BAD DOI: JBC DOI returns 404 on doi.org. Real Lally 1989 leukotoxin cloning paper is in Biochem Biophys Res Commun (PMID 2647082).; (b) family (LktA (Mannheimia leukotoxin)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate identity is independently supported; (c) locu"
},
{
"gene": "lktA",
"uniprot": "P55117",
"organism": "Bibersteinia trehalosi (Pasteurella trehalosi)",
"ss_type": "T1SS",
"locus_tag": "",
"refseq_genome": "",
"doi_old": "10.1016/S0021-9258(19)84909-X",
"notes": "(a) BAD DOI: JBC DOI returns 404 on doi.org. Real Lally 1989 leukotoxin cloning paper is in Biochem Biophys Res Commun (PMID 2647082).; (b) family (LktA (Mannheimia leukotoxin)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate identity is independently supported; (c) locu"
},
{
"gene": "lktA",
"uniprot": "Q9ETX2",
"organism": "Mannheimia glucosida",
"ss_type": "T1SS",
"locus_tag": "",
"refseq_genome": "",
"doi_old": "10.1016/S0021-9258(19)84909-X",
"notes": "(a) BAD DOI: JBC DOI returns 404 on doi.org. Real Lally 1989 leukotoxin cloning paper is in Biochem Biophys Res Commun (PMID 2647082).; (b) family (LktA (Mannheimia leukotoxin)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate identity is independently supported; (c) locu"
},
{
"gene": "lktA",
"uniprot": "P0C082",
"organism": "Mannheimia haemolytica (Pasteurella haemolytica)",
"ss_type": "T1SS",
"locus_tag": "",
"refseq_genome": "",
"doi_old": "10.1016/S0021-9258(19)84909-X",
"notes": "(a) BAD DOI: JBC DOI returns 404 on doi.org. Real Lally 1989 leukotoxin cloning paper is in Biochem Biophys Res Commun (PMID 2647082).; (b) family (LktA (Mannheimia leukotoxin)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate identity is independently supported; (c) locu"
},
{
"gene": "lktA",
"uniprot": "P55123",
"organism": "Pasteurella haemolytica-like sp. (strain 5943B)",
"ss_type": "T1SS",
"locus_tag": "",
"refseq_genome": "",
"doi_old": "10.1016/S0021-9258(19)84909-X",
"notes": "(a) BAD DOI: JBC DOI returns 404 on doi.org. Real Lally 1989 leukotoxin cloning paper is in Biochem Biophys Res Commun (PMID 2647082).; (b) family (LktA (Mannheimia leukotoxin)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate identity is independently supported; (c) locu"
},
{
"gene": "rtxA",
"uniprot": "A1YKW7",
"organism": "Kingella kingae",
"ss_type": "T1SS",
"locus_tag": "",
"refseq_genome": "",
"doi_old": "10.1128/JB.186.20.6985-6991.2004",
"notes": "(a) BAD DOI: DOI returns 404 on doi.org/Crossref. Cannot verify; Kehl-Fie/St Geme Kingella RtxA paper is at PMID 15466052 or similar (different DOI).; (b) family (RtxA (Kingella kingae)) is a validated T1SS substrate family per Linhartova 2010 (10.1111/j.1574-6976.2010.00231.x) and other reviews; cited DOI is wrong but family-level substrate identity is independently supported;"
},
{
"gene": "pelA",
"uniprot": "P0C1A3",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_03370",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "pelB",
"uniprot": "E0SJQ6",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_04193",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "pelC",
"uniprot": "E0SJQ7",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_04192",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "pelD",
"uniprot": "E0SAZ3",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_03372",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "pelE",
"uniprot": "P0C1A5",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_03371",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "pelI",
"uniprot": "E0SM75",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_00058",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "pelL",
"uniprot": "P0C1A7",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_02794",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "pemA",
"uniprot": "P0C1A9",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_03374",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "pemB",
"uniprot": "Q47474",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_03435",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "paeY",
"uniprot": "E0SAZ4",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_03373",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "paeX",
"uniprot": "E0SDD1",
"organism": "Dickeya dadantii 3937",
"ss_type": "T2SS",
"locus_tag": "Dda3937_03363",
"refseq_genome": "NC_014500.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "aerA",
"uniprot": "P09167",
"organism": "Aeromonas hydrophila ATCC 7966",
"ss_type": "T2SS",
"locus_tag": "aerA",
"refseq_genome": "NC_008570.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "lip",
"uniprot": "P40600",
"organism": "Aeromonas hydrophila ATCC 7966",
"ss_type": "T2SS",
"locus_tag": "lip",
"refseq_genome": "NC_008570.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "amyA",
"uniprot": "P41131",
"organism": "Aeromonas hydrophila ATCC 7966",
"ss_type": "T2SS",
"locus_tag": "amyA",
"refseq_genome": "NC_008570.1",
"doi_old": "doi:10.1146/annurev-micro-102215-095506",
"notes": "DOI:DOI does not resolve; not present in Annu Rev Microbiol volume 70 (2016) which uses prefix annurev-micro-102215. Citation is broken/fabricated."
},
{
"gene": "prtA",
"uniprot": "Q3BPB1",
"organism": "Xanthomonas euvesicatoria 85-10",
"ss_type": "T2SS",
"locus_tag": "XCV3671",
"refseq_genome": "NC_007508.1",
"doi_old": "doi:10.1128/jb.00322-15",
"notes": "UniProt:Q3BFK4 returned HTTP 404 (entry obsolete/missing). Cannot verify."
},
{
"gene": "stmPr1",
"uniprot": "Q93IQ4",
"organism": "Stenotrophomonas maltophilia K279a",
"ss_type": "T2SS",
"locus_tag": "Smlt0861",
"refseq_genome": "NC_010943.1",
"doi_old": "doi:10.1128/IAI.00388-15",
"notes": "DOI:DOI does not resolve. The actual DuMont/Karaba/Cianciotto 2015 StmPr1/StmPr2 T2SS paper is 10.1128/IAI.00672-15. Citation is wrong."
},
{
"gene": "zmpA",
"uniprot": "Q7WSN3",
"organism": "Burkholderia cenocepacia J2315",
"ss_type": "T2SS",
"locus_tag": "BCAL2284",
"refseq_genome": "NC_011000.1",
"doi_old": "doi:10.1099/mic.0.27353-0",
"notes": "DOI:DOI does not resolve. The Kooi/Corbett/Sokol 2005 ZmpA paper is 10.1128/JB.187.13.4421-4429.2005 (J Bacteriol). Citation is wrong."
},
{
"gene": "GALLS",
"uniprot": "-",
"organism": "Agrobacterium rhizogenes",
"ss_type": "T4SS",
"locus_tag": "-",
"refseq_genome": "-",
"doi_old": "10.1094/MPMI-19-1145",
"notes": "DOI=NOT_FOUND (DOI 404 at crossref). Biology=PARTIAL: GALLS is A. rhizogenes VirB-Ri functional equivalent of VirE2 (Hodges 2006 MPMI-19-1145). No UniProt/locus_tag. Instance correct."
},
{
"gene": "CagF",
"uniprot": "O25276",
"organism": "Helicobacter pylori 26695",
"ss_type": "T4SS",
"locus_tag": "HP_0543",
"refseq_genome": "NC_000915.1",
"doi_old": "10.1074/jbc.M501633200",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=PARTIAL: CagF is a Cag PAI chaperone (HP_0543) for CagA, but uniprot O25257 actually maps to HP_0520 (Cag1) in UniProt, not HP_0543. uniprot/locus_tag mismatch. Note: CagF is a T4SS chaperone, not strictly a translocated substrate."
},
{
"gene": "TcpB_Btp1",
"uniprot": "Q8YF53",
"organism": "Brucella melitensis 16M",
"ss_type": "T4SS",
"locus_tag": "BMEI1674",
"refseq_genome": "NC_003317.1",
"doi_old": "10.1128/IAI.00750-12",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: TcpB/Btp1 Q8YF53/BMEI1674 confirmed in UniProt as B. melitensis virulence factor with TIR domain. Brucella has only VirB. | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "TcpB_Btp1_Ba",
"uniprot": "Q2YPC4",
"organism": "Brucella abortus 2308",
"ss_type": "T4SS",
"locus_tag": "BAB1_0279",
"refseq_genome": "NC_007618.1",
"doi_old": "10.1128/IAI.00750-12",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: TcpB/Btp1 Q2YPC4/BAB1_0279 confirmed in UniProt as B. abortus paralog; VirB substrate. | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "BtpB",
"uniprot": "Q2YNA0",
"organism": "Brucella abortus 2308",
"ss_type": "T4SS",
"locus_tag": "BAB1_0782",
"refseq_genome": "NC_007618.1",
"doi_old": "10.1128/IAI.00750-12",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: BtpB BAB1_0782 is the second TIR-domain effector of Brucella (Salcedo, Smith). VirB substrate. | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "TRP47",
"uniprot": "Q2GHU2",
"organism": "Ehrlichia chaffeensis Arkansas",
"ss_type": "T4SS",
"locus_tag": "ECH_0166",
"refseq_genome": "NC_007799.1",
"doi_old": "10.1128/IAI.05225-11",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: TRP47 Q2GHU2/ECH_0166 confirmed in UniProt as gp47; secreted via T4SS (Luo 2008). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "TRP32",
"uniprot": "Q2GHT8",
"organism": "Ehrlichia chaffeensis Arkansas",
"ss_type": "T4SS",
"locus_tag": "ECH_0170",
"refseq_genome": "NC_007799.1",
"doi_old": "10.1128/IAI.00153-10",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: TRP32 Q2GHT8/ECH_0170 confirmed in UniProt as Variable_length_PCR_target protein. Documented Ehrlichia T4SS substrate (Luo, Wakeel). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "Etf-1",
"uniprot": "Q2GG12",
"organism": "Ehrlichia chaffeensis Arkansas",
"ss_type": "T4SS",
"locus_tag": "ECH_0825",
"refseq_genome": "NC_007799.1",
"doi_old": "10.1073/pnas.1218674110",
"notes": "DOI=NOT_FOUND (DOI 404 (Etf-1 PNAS expected)). Biology=GOOD: Etf-1 Q2GG12/ECH_0825 confirmed in UniProt; Etf-1 is Ehrlichia translocated factor-1, VirB substrate (Liu 2012). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "Sca4",
"uniprot": "Q52658",
"organism": "Rickettsia conorii Malish 7",
"ss_type": "T4SS",
"locus_tag": "RC0667",
"refseq_genome": "NC_003103.1",
"doi_old": "10.1126/science.abm9836",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: Sca4 Q52658/RC0667 confirmed in UniProt as R. conorii surface antigen. Lamason 2018 + Aistleitner 2020 documented Sca4 as Rickettsia VirB substrate. | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "SdjA",
"uniprot": "Q5ZTK6",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg2155",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1038/s41589-020-0700-0",
"notes": "DOI=NOT_FOUND (DOI 404 (SdjA Nat Chem Biol expected)). Biology=GOOD: SdjA lpg2155: Wan 2019 Nat Chem Biol identified SdjA as PR-toxin DotIcm substrate. | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "SidH",
"uniprot": "Q5ZRQ1",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg2829",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1126/science.aab4151",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: SidH lpg2829: large DotIcm substrate (Luo 2007). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "SidI",
"uniprot": "Q5ZSL3",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg2504",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1073/pnas.1810122115",
"notes": "DOI=NOT_FOUND (DOI 404 (Lem4/SidI PNAS expected)). Biology=GOOD: SidI lpg2504: DotIcm substrate inhibiting host translation (Shen 2009). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "Lem4",
"uniprot": "Q5ZX05",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg0928",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1073/pnas.1810122115",
"notes": "DOI=NOT_FOUND (DOI 404 (Lem4/SidI PNAS expected)). Biology=GOOD: Lem4 lpg0928: DotIcm substrate (HAD-superfamily phosphatase, Quaile 2018). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "RidL",
"uniprot": "Q5ZT54",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg2311",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1083/jcb.201406048",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: RidL lpg2311: retromer-binding DotIcm substrate (Finsel 2013, Yao 2018). (DOI in row 10.1083/jcb.201406048 is 404 but biology validated.) | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "MitF",
"uniprot": "Q5ZRR2",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg2818",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1126/science.aab1567",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: MitF lpg2818: DotIcm substrate (Escoll 2017). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "Lgt1",
"uniprot": "Q5ZVS2",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg1368",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1126/science.1158333",
"notes": "DOI=NOT_FOUND (DOI 404; not a real Science DOI). Biology=GOOD: Lgt1 lpg1368: DotIcm substrate (glucosyltransferase, Belyi 2008). (Bad DOI in row.) | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "Lgt2",
"uniprot": "Q5ZRL9",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg2862",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1126/science.1158333",
"notes": "DOI=NOT_FOUND (DOI 404; not a real Science DOI). Biology=GOOD: Lgt2 lpg2862: DotIcm substrate (Belyi 2008). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "Lgt3",
"uniprot": "Q5ZVF2",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg1488",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1126/science.1158333",
"notes": "DOI=NOT_FOUND (DOI 404; not a real Science DOI). Biology=GOOD: Lgt3 lpg1488: DotIcm substrate (Belyi 2008). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "MavQ",
"uniprot": "Q5ZRR7",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg2813",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1083/jcb.201906133",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: MavQ lpg2813: DotIcm substrate, PI 3-kinase regulator (Hsieh 2021). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "RavD",
"uniprot": "Q5ZZ51",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg0160",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1074/jbc.M113.502120",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: RavD lpg0160: DotIcm substrate (Luo 2013). DUB activity (Pike 2019). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "LotA",
"uniprot": "Q5ZTB4",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg2248",
"refseq_genome": "NC_002942.5",
"doi_old": "10.15252/embj.2019103258",
"notes": "DOI=NOT_FOUND (DOI 404 at crossref (LotA EMBO J expected)). Biology=GOOD: LotA lpg2248: Legionella OTU deubiquitinase (Kubori/Nagai 2018 EMBO J). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "LotB",
"uniprot": "Q5ZV21",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg1621",
"refseq_genome": "NC_002942.5",
"doi_old": "10.15252/embj.2020107656",
"notes": "DOI=NOT_FOUND (DOI 404 at crossref (LotB/LotC EMBO J expected)). Biology=GOOD: LotB lpg1621: Legionella OTU DUB (Shin 2020 EMBO J). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "LotC",
"uniprot": "Q5ZSI8",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg2529",
"refseq_genome": "NC_002942.5",
"doi_old": "10.15252/embj.2020107656",
"notes": "DOI=NOT_FOUND (DOI 404 at crossref (LotB/LotC EMBO J expected)). Biology=GOOD: LotC lpg2529: Legionella OTU DUB (Shin 2020). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "SidP",
"uniprot": "Q5ZZ81",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg0130",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1074/jbc.M111.291849",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: SidP lpg0130: PI phosphatase DotIcm substrate (Toulabi 2013). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "Lem3",
"uniprot": "Q5ZXN5",
"organism": "Legionella pneumophila Philadelphia 1",
"ss_type": "T4SS",
"locus_tag": "lpg0696",
"refseq_genome": "NC_002942.5",
"doi_old": "10.1038/nature10337",
"notes": "DOI=WRONG (Resolves to Lgr5/R-spondin, NOT Lem3). Biology=GOOD: Lem3 lpg0696: DotIcm substrate, de-AMPylase paralog (Tan 2011). (Bad DOI in row.) | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "CinF",
"uniprot": "Q83FA5",
"organism": "Coxiella burnetii RSA493",
"ss_type": "T4SS",
"locus_tag": "CBU0041",
"refseq_genome": "NC_002971.4",
"doi_old": "10.1073/pnas.2102337118",
"notes": "DOI=NOT_FOUND (DOI 404 (CinF Coxiella PNAS expected)). Biology=GOOD: CinF CBU0041: Coxiella Dot/Icm effector (Carey 2011). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
},
{
"gene": "CBU0626",
"uniprot": "Q83DR9",
"organism": "Coxiella burnetii RSA493",
"ss_type": "T4SS",
"locus_tag": "CBU0626",
"refseq_genome": "NC_002971.4",
"doi_old": "10.1128/IAI.00558-19",
"notes": "DOI=NOT_FOUND (DOI 404). Biology=GOOD: CBU0626: Coxiella Dot/Icm substrate (Larson/Heinzen). | (a) FAILS: cited DOI does not point to the right paper; (b)/(c)/(d) PASS."
}
]

const TYPE_HINT = {
  T1SS: 'Type I secretion system (RTX exporter: ABC transporter + membrane-fusion protein + TolC). The substrate is the secreted RTX toxin/protein itself.',
  T2SS: 'Type II secretion system (the Gsp/Out main terminal branch). The substrate is a folded exoenzyme/toxin secreted across the outer membrane.',
  T4SS: 'Type IV secretion system (e.g. VirB/VirD4, Dot/Icm). The substrate is a protein effector translocated into a host or recipient cell.',
}

const SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['gene','uniprot','ss_type','status','corrected_doi','doi_resolves','is_review','primary_doi_if_review','paper_title','paper_year','paper_journal','verbatim_quote','quote_source','confidence','reasoning','searched'],
  properties: {
    gene: { type: 'string' },
    uniprot: { type: 'string' },
    ss_type: { type: 'string' },
    status: { type: 'string', enum: ['RESOLVED', 'CONTRADICTED', 'NOT_FOUND'] },
    corrected_doi: { type: ['string', 'null'] },
    doi_resolves: { type: 'boolean' },
    is_review: { type: 'boolean' },
    primary_doi_if_review: { type: ['string', 'null'] },
    paper_title: { type: ['string', 'null'] },
    paper_year: { type: ['string', 'null'] },
    paper_journal: { type: ['string', 'null'] },
    verbatim_quote: { type: ['string', 'null'] },
    quote_source: { type: ['string', 'null'] },
    confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
    reasoning: { type: 'string' },
    searched: { type: 'string' },
  },
}

function promptFor(r) {
  return [
    'You are verifying and repairing a literature citation for ONE experimentally-validated bacterial secretion-system effector (a SECRETED SUBSTRATE / cargo, NOT machinery). This feeds a benchmark; correctness and ZERO hallucination matter more than completeness.',
    '',
    'PROTEIN:',
    '- gene: ' + r.gene,
    '- UniProt: ' + (r.uniprot || '(none)'),
    '- organism: ' + r.organism,
    '- secretion system type: ' + r.ss_type + ' -- ' + (TYPE_HINT[r.ss_type] || ''),
    '- locus_tag: ' + (r.locus_tag || '(none)') + '  genome: ' + (r.refseq_genome || '(none)'),
    '',
    'PROBLEM: the stored primary-reference DOI is broken/wrong: ' + r.doi_old,
    'Audit note (hints, may be partly wrong): ' + r.notes,
    '',
    'TASK: Find the DOI of a peer-reviewed paper that EXPERIMENTALLY establishes that THIS specific protein is secreted by / translocated by / is a substrate of the named secretion system in this organism.',
    '',
    'RULES (anti-hallucination, strict):',
    '1. Prefer the primary experimental paper that demonstrated secretion/translocation. The audit note often names the author+year (e.g. "Belyi 2008", "Luo 2007"); use that to find the real DOI. A review is acceptable ONLY if it explicitly states the protein is experimentally-verified cargo; then give BOTH the review DOI and the primary DOI (in primary_doi_if_review).',
    '2. You MUST confirm the corrected DOI actually resolves. Check Crossref (https://api.crossref.org/works/{doi}) or https://doi.org/{doi} or PubMed. NEVER return a DOI you have not confirmed resolves to the paper you describe.',
    '3. Provide a VERBATIM quote (copied exactly) from that paper abstract or body stating the protein is secreted/translocated/a substrate/effector of this system. If you cannot find such a sentence, you have NOT verified secretion.',
    '4. If the experimental literature CONTRADICTS the secretion-system label (the protein is shown NOT to be secreted by this system, or is secreted by a different route), set status="CONTRADICTED" and give the verbatim contradicting quote + its DOI.',
    '5. NEVER fabricate a DOI, title, author, year, or quote. If you cannot confidently find a qualifying paper, set status="NOT_FOUND" and list what you searched.',
    '6. Tools: use the PubMed MCP tools (search_articles, lookup_article_by_citation, get_article_metadata, get_full_text_article) and WebSearch / WebFetch. Crossref REST is open and needs no auth.',
    '',
    'Echo gene/uniprot/ss_type back exactly. status is one of RESOLVED / CONTRADICTED / NOT_FOUND. Use null for any unknown string field.',
  ].join('\n')
}

phase('Repair')
log('Repairing ' + ROWS.length + ' broken-DOI effector citations')
const results = await parallel(ROWS.map((r) => () =>
  agent(promptFor(r), { label: 'doi:' + r.gene, phase: 'Repair', schema: SCHEMA })
))
const records = results.filter(Boolean)
const byStatus = {}
for (const x of records) byStatus[x.status] = (byStatus[x.status] || 0) + 1
log('done: ' + JSON.stringify(byStatus))
return { total: ROWS.length, returned: records.length, byStatus, records }
