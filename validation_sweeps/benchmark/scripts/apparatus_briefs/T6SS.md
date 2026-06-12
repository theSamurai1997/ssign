# T6SS apparatus brief (Type VI secretion system)

**Mechanism (plain):** an inverted contractile phage tail anchored in the membrane. A
sheath contracts and fires a spike-tipped tube out of the cell (and into a neighbouring
cell), carrying toxin effectors. ~13 core genes (TssA-TssM) plus accessories.

**IMPORTANT for this benchmark:** Hcp, VgrG and PAAR are the secreted tube/spike/tip and
ARE machinery (they belong in the answer key). They were excluded from the *effector* gold
set as "apparatus-as-effector"; that exclusion does not apply here. Cargo toxins fused to
an evolved Hcp/VgrG/PAAR, or delivered as separate Tse/Tae/Tle/Tde/Tme proteins, are
effectors, NOT machinery.

| Unified | Role | Common synonyms |
|---|---|---|
| TssA | tail assembly / cap | TssA, ImpA, SciA |
| TssB | sheath (inner) | TssB, VipA, ImpB |
| TssC | sheath (outer) | TssC, VipB, ImpC |
| TssD | Hcp: hexameric inner tube | Hcp, ImpC?; "Hcp1/Hcp2" |
| TssE | baseplate (gp25-like) | TssE |
| TssF | baseplate | TssF, ImpG, SciC |
| TssG | baseplate | TssG, ImpH |
| TssH | ClpV AAA+ ATPase (recycles sheath) | TssH, ClpV, SciG |
| TssI | VgrG: trimeric spike (tube tip) | VgrG, "VgrG1/2/3" |
| TssJ | OM lipoprotein (membrane complex) | TssJ, SciN, Lip |
| TssK | baseplate–membrane link | TssK, ImpJ, SciP |
| TssL | IM (membrane complex) | TssL, ImpK, SciP/DotU-like |
| TssM | IM (membrane complex, large) | TssM, IcmF, ImpL |
| PAAR | spike-sharpening tip | PAAR-repeat protein |
| accessory | regulation/adaptors | TagA-TagN, TagF, Fha/TagH, PpkA/PppA (threonine phosphorylation), TssM-associated |

**Notes for curation:**
- The corpus `sys_instance_label` is often "T6SS" or a representative locus. Many genomes
  carry one main T6SS locus; some have 2-3 (count only the instance the effectors belong to
  if the paper distinguishes them).
- Count Hcp/VgrG/PAAR as machinery. Do NOT count the delivered toxin effectors
  (Tse/Tae/Tle/Tde/Tme/Tre/Rhs-CT, amidases, nucleases, etc.) or their immunity proteins.
- ClpV (TssH) is essential machinery (sheath recycling), count it.
