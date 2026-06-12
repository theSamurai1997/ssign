# T4SS apparatus brief (Type IV secretion system)

**Mechanism (plain):** a versatile machine that translocates proteins (and sometimes DNA)
into host or recipient cells. Two architecture families: **type-IVA** (the VirB/D4
archetype, ~12 genes) and **type-IVB** (Legionella/Coxiella Dot/Icm, ~27 genes). The
translocated effectors are NOT machinery.

## Type-IVA core (VirB/D4 archetype)

| Family | Role | Common synonyms (Agrobacterium / Brucella / Bordetella Ptl / Helicobacter Cag / Bartonella) |
|---|---|---|
| VirB1 | lytic transglycosylase (cell-wall remodelling) | VirB1; Ptl? |
| VirB2 | major pilin | VirB2; CagC |
| VirB3 | IM, pilus assembly | VirB3; PtlB |
| VirB4 | IM ATPase (energiser) | VirB4; PtlC; CagE; TrwK |
| VirB5 | minor pilin / tip adhesin | VirB5; CagL |
| VirB6 | polytopic IM channel | VirB6; PtlD; CagW |
| VirB7 | OM lipoprotein (core complex) | VirB7; CagT |
| VirB8 | IM assembly factor | VirB8; PtlE; CagV |
| VirB9 | OM core-complex (with VirB7/10) | VirB9; PtlF; CagX |
| VirB10 | IM-OM channel ring | VirB10; PtlG; CagY |
| VirB11 | cytoplasmic ATPase | VirB11; PtlH; HP0525/CagAlpha; TrwD |
| VirD4 | type-IV coupling protein (substrate receptor ATPase) | VirD4; PtlI?; Cagbeta/HP0524; TrwB |

Bordetella pertussis toxin secretion uses the **Ptl** (PtlA-I) operon; *Helicobacter pylori*
Cag uses **Cag** genes (cagT/V/W/X/Y + cagE + HP0524/0525). Brucella, Bartonella, Anaplasma,
Ehrlichia, Wolbachia, Rickettsia use **virB1-11 + virD4** directly.

## Type-IVB core (Dot/Icm, Legionella / Coxiella)

A distinct ~27-gene system: **dotA-dotU, icmB-icmX** (e.g. DotA, DotB ATPase, DotC, DotD,
DotF, DotG/IcmE, DotH/IcmK, DotO/IcmB ATPase, DotU, IcmF, IcmS, IcmW, IcmR, IcmQ, IcmX...).
Effectors are genome-dispersed (proximity rule does not apply); in this benchmark the
Dot/Icm effector flood is **GATED** pending Phase 2 detection.

**Notes for curation:**
- The corpus `sys_instance_label` is often "VirB" or a representative locus; identify the
  family (VirB/D4 vs Dot/Icm vs Cag vs Ptl) from the organism + paper.
- Do NOT list translocated effectors (VirE2/VirE3/VirF/VirD5 are effectors, not apparatus;
  CagA is an effector, the Cag genes above are apparatus).
- *virD2* (relaxase) and *virC1/2* are T-DNA processing, not the translocation channel;
  note but do not count as channel apparatus.
