# T5SS effector-sourcing brief (secretion-classifier dataset, group 3)

T5SS (autotransporter / two-partner secretion) is absent from the benchmark corpus
because it has no separable apparatus to enumerate, the translocator is part of (or a
dedicated partner of) the secreted protein itself. For the classifier dataset we still
want T5SS positives, sourced **by subtype** under the same anti-hallucination contract
as the machinery answer key.

## The label convention (load-bearing)

T5SS splits into "self-secreted" autotransporters and a two-partner exception:

| Subtype | What is secreted | `self_secreted` | Why |
|---|---|---|---|
| **T5aSS** classical autotransporter | the **passenger** domain of the same polypeptide | **true** | passenger + β-barrel translocator are one gene; the protein secretes itself |
| **T5bSS** two-partner secretion (TPS) | **TpsA** exoprotein, a *separate* gene from its **TpsB** pore | **false** | TpsA is a genuine (protein, instance) substrate, co-localized with its TpsB transporter, treat like any other substrate |
| **T5cSS** trimeric autotransporter (TAA) | the trimeric passenger of the same polypeptide | **true** | self-secreted, like T5a but trimeric (YadA-like) |
| **T5dSS** (rare) | patatin-like passenger of the same polypeptide | **true** | hybrid autotransporter (passenger + AT β-domain + POTRA); self-secreted |
| **T5eSS** inverted autotransporter | passenger (C-terminal) of the same polypeptide | **true** | inverted domain order, still self-secreted |

Why this matters for the model: T5b TpsA rows teach normal substrate selection (a
distinct protein near its transporter). T5a/c/d/e rows are **self-secreted positives**,
they teach "is this autotransporter really a secreted protein," i.e. a confidence
score on MacSyFinder's autotransporter call, not substrate-vs-machinery discrimination.
Keep them flagged so the model/eval can down-weight or hold them out.

## Seed Pfam families (recognition vocabulary for the agents)

| Subtype | Translocator domain (Pfam) | Notes |
|---|---|---|
| T5aSS | **PF03797** (Autotransporter β-domain, "Surface_Ag_2") | C-terminal ~275-aa 12-stranded β-barrel; passenger N-terminal |
| T5bSS | TpsB: **PF01103** (Omp85/Bac_surface_Ag) + **POTRA** domains; TpsA: **PF05860** (Haemagg_act, TPS N-terminal secretion domain) | two genes, usually adjacent (*tpsBA* / *fhaC-fhaB*) |
| T5cSS | **PF03895** (YadA_anchor / trimeric AT C-terminal) | trimeric; head-stalk-anchor architecture |
| T5dSS | PF03797-like AT β-domain fused to a periplasmic POTRA + patatin (PF01734) passenger | PlpD is the type example |
| T5eSS | intimin/invasin β-domain (inverted: barrel N-terminal, passenger C-terminal) | Big_2 / intimin domains |

## Canonical examples per subtype (anchors the agents can start from; NOT a substitute for sourcing)

- **T5aSS (self_secreted=true):** SPATEs, EspP/EspC (E. coli), Hbp/Tsh, Pic/Sat/Vat,
  IgA1 protease (Neisseria), AIDA-I, Antigen 43 (Ag43), BrkA + Pertactin (Bordetella),
  IcsA/VirG (Shigella), App/Hap (Haemophilus).
- **T5bSS (self_secreted=false; the TpsA is the substrate):** FhaB/FHA + FhaC
  (Bordetella), HMW1 + HMW1B, HMW2 (Haemophilus), ShlA + ShlB (Serratia), HpmA + HpmB
  (Proteus), HxuA + HxuB (Haemophilus), CdiA + CdiB (E. coli / Burkholderia CDI),
  LspA1/A2 (Haemophilus ducreyi).
- **T5cSS (self_secreted=true):** YadA (Yersinia), Hia / Hsf (Haemophilus), UspA1 /
  UspA2 (Moraxella), BadA (Bartonella henselae), Vomp (B. quintana), SadA (Salmonella),
  NadA (Neisseria meningitidis), EibA-F (E. coli), DsrA (Haemophilus ducreyi).
- **T5dSS (self_secreted=true):** PlpD (Pseudomonas aeruginosa).
- **T5eSS (self_secreted=true):** Intimin/EaeA (EHEC/EPEC), Invasin/InvA (Yersinia).

## Anti-hallucination contract (same as benchmark task 3.4)

Each sourced example requires, verbatim in the agent output:
1. a **verbatim quote** from the characterizing paper naming this protein as a secreted
   T5SS substrate/autotransporter of the stated subtype;
2. a **resolvable DOI** (will be Crossref/DOI.org-verified downstream, agent word not trusted);
3. a **real locus_tag** that exists in a named RefSeq genome. Locus_tags are NEVER
   invented; an example without a resolvable locus_tag is recorded as **unplaceable**,
   not fabricated, and excluded from the placed set.

Output schema (per example, JSON): `gene, organism, refseq_genome, locus_tag, ss_type="T5SS",
subtype (T5aSS|T5bSS|T5cSS|T5dSS|T5eSS), self_secreted (true|false), uniprot (optional),
primary_ref (DOI), quote, note`. For T5bSS also give the partner TpsB gene/locus in `note`.

Depth: T5a/T5b/T5c are well-characterized, source broadly. T5d/T5e are rare, report
whatever has a clean primary reference; do not pad.
