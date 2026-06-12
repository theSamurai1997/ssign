# Why ssign misses effectors it "could never have found" (unreachable @±3)

The recall figure shows 386 testable effectors that ssign missed AND that sit >3 genes from their
curated machinery (out of proximity reach). Are those misses *expected* from the known genomic biology
of each system, or do they hide pipeline / answer-key problems? Short answer: **almost all are expected
biology** — effectors of T2/T3/T4SS are genuinely genome-dispersed, so a ±3 proximity rule structurally
cannot reach them. T1SS is the operonic exception (median distance = 1 gene), and its few misses are
real biological exceptions, not failures. Two genuine answer-key issues surfaced (below).

Method: per-effector gene-order distance to the nearest anchored machinery gene (`ceiling_per_effector`),
plus a per-system literature review (one agent per system, primary papers cited inline). Distances:
`figures/summary/04_distance_to_machinery.png`.

| system | n testable | median genes to machinery | within ±3 | expected to be reachable? |
|---|---|---|---|---|
| **T1SS** | 23 | **1** | 20 | **yes — operonic** (RTX toxin beside its ABC transporter) |
| T2SS | 64 | 302 | 1 | no — substrates scattered genome-wide |
| T3SS | 213 | 45 | 48 | partly — local complement + dispersed majority |
| T4SS | 84 | 203 | 6 | no — effectors dispersed, portable C-terminal signal |
| T6SS | 72 | 232 | 16 | mixed — see anchor caveat |

## T1SS: operonic as expected; the 5 misses are real exceptions (not bugs)

T1SS effectors sit a **median of 1 gene** from their transporter — your intuition is right. Only 5 of
25 are unreachable, and the literature confirms every one is genuine non-operonic biology:

- **Serralysin** (*Serratia marcescens*): secreted by the **LipBCD generalist exporter at a separate
  locus** (shared with LipA, SlaA). The corpus's own DOI (10.1128/jb.179.15.4754) states LipBCD secretes
  both LipA and PrtA. Genuine.
- **ApxIIA** (*Actinobacillus pleuropneumoniae*): textbook **secretion in trans** — the apxIICA operon
  has no transporter; ApxIB/D from the separate apxI operon export it (10.1099/00221287-139-8-1723).
- **FrpC** (*Neisseria meningitidis*): a **functional but genome-scattered** T1SS — the cited paper is
  literally titled "Unusual genetic organization of a functional type I protein secretion system"
  (10.1128/IAI.73.9.5554) and reports hlyB, hlyD/tolC scattered and unlinked to frpC. Only frpD (an OM
  *anchoring* lipoprotein, not a transporter) is adjacent. **Re: the "ssign should find the TolC" point**
  — ssign detected *no* T1SS in this genome and there is no TolC/HlyB/HlyD near frpC; the apparatus is
  real but dispersed, so there was nothing adjacent to detect. The answer key's "TolC at 1340 genes" is a
  spurious product-match.
- **TRP47 / TRP32** (*Ehrlichia chaffeensis*): genuine T1SS (Hly apparatus at a separate locus, TRP genes
  scattered) — and the literature **explicitly tested and rejected T4SS** for these exact proteins
  (10.3389/fcimb.2011.00022). **Answer-key bug:** our machinery resolver anchored them on **VirB8 (a T4SS
  gene)**. The reachability verdict (unreachable) is still correct, but the anchor is wrong.

**Verdict: 5/5 genuine biology.** T1SS is far less reliably operonic than the canonical HlyA-HlyBD picture:
shared exporters, in-trans secretion, and scattered apparatus all break adjacency.

## T2SS / T3SS / T4SS: dispersal is the expected biology

- **T2SS — proximity captures almost none, as expected.** The apparatus is one compact 12-16-gene operon,
  but substrates are a functionally diverse set (toxins, proteases, lipases, cellulases) recruited
  *post-translationally by folded-state recognition in the periplasm* — no genetic linkage. Legionella's
  Lsp substrates (PlaA, CelA, ProA…) are scattered chromosome-wide. The one clustered case (Klebsiella
  pullulanase pulA, embedded in the pul operon) is the textbook exception. (Korotkov & Sandkvist; Douzi/
  Filloux/Voulhoux 10.1098/rstb.2011.0204.)
- **T3SS — few reachable, by horizontal-acquisition biology.** Effectors are acquired piecemeal on
  prophages, plasmids and separate islands and coordinated by *shared regulation* (HrpL, SsrB, HilA), not
  adjacency (Brown & Finlay 10.4161/mge.1.2.16733). Chromosomal systems (Salmonella SPI, EPEC/EHEC LEE,
  P. syringae Hrp, Chlamydia) keep a small local complement and disperse the rest; only the *plasmid*
  systems (Yersinia pYV, Shigella) co-localize apparatus + effectors — which is why T3SS's median (45) is
  lower than T2/T4.
- **T4SS — almost none reachable; effectors fully decoupled.** Substrates are recognized by a *portable
  C-terminal signal* read by the coupling protein, so location is unconstrained: Legionella/Coxiella
  Dot/Icm translocate ~300 effectors scattered genome-wide (Schroeder 10.3389/fcimb.2017.00528; Carey
  10.1371/journal.ppat.1002056). Exceptions where proximity *would* work: Bartonella bep cluster,
  Helicobacter cag (n=1), the small Agrobacterium vir set.

For all three, the high unreachable fraction is the **correct characterization of the biology**, not a
ssign or benchmark failure. Proximity is simply the wrong model for systems whose effectors aren't
co-encoded with the machine.

## T6SS: intermediate — and the benchmark anchor undercounts it

T6SS splits. Many "specialized/cargo" effectors are encoded **immediately beside the vgrG / paar / hcp**
spike gene that carries them, with the cognate **immunity gene directly downstream** (a strong,
T6SS-specific signal; Whitney 10.1074/jbc.M113.488320; vgrG-linkage as an effector marker
10.1073/pnas.1406110111). **But** those spike genes are often orphans / auxiliary modules far from the
*core* tss cluster (TssM/ClpV). Our answer key anchors on the core cluster, so an effector one gene from
its own orphan vgrG can read as ">200 genes from machinery."

**Actionable caveat:** the T6SS ceiling (and "unreachable@3" count) is likely an *under*estimate — a
vgrG/paar/hcp anchor + an effector–immunity-pair test would recover the spike-loaded class. ssign itself
detects vgrG as a T6SS component, so some of the 44 "unreachable" T6SS effectors may in fact be reachable
once anchored correctly. This is the one system where the proximity rule's apparent failure is partly a
benchmark-anchor choice, not biology.

## Follow-ups this surfaced

1. **Answer-key fix:** TRP47/TRP32 machinery anchored on T4SS VirB8 — re-anchor on the Hly T1SS (or mark
   machinery_unanchored). Verdict (unreachable) unchanged, but the anchor gene is wrong.
2. **T6SS re-anchoring:** add a nearest-vgrG/paar/hcp distance + immunity-pair signal alongside the
   core-cluster anchor; recompute the T6SS ceiling. Likely raises T6SS reachable@3.
3. **frpC** answer-key TolC anchor (1340 genes, product-tier) is spurious — drop or mark unanchored.

## Bottom line

The "couldn't have found" effectors are overwhelmingly **expected biology**: T1SS is operonic (and
ssign does well there), while T2/T3/T4SS effectors are genuinely genome-dispersed and proximity cannot
reach them by design. This is the strongest argument that a *learned, sequence/context model* (not a
fixed proximity window) is needed to recover effectors of the dispersed systems — and that a single
global proximity rule is the wrong abstraction across systems whose effector-genomics differ this much.
