# T1SS apparatus brief (Type I secretion system)

**Mechanism (plain):** a one-step pump that takes an unfolded substrate straight from
the cytoplasm to outside the cell, skipping the periplasm. Built from just three
proteins spanning both membranes. The substrate (e.g. an RTX toxin or adhesin) carries
a C-terminal secretion signal and is NOT part of the machinery.

**Core machinery (3 components):**

| Family | Role | Common synonyms / examples |
|---|---|---|
| ABC transporter | inner-membrane ATPase that energises export and forms the IM channel | HlyB, PrtD, LipB, AprD, CyaB, TolC-associated ABC; Pfam ABC_membrane + ABC_tran |
| Membrane fusion protein (MFP) / adaptor | periplasm-spanning adaptor, IM-anchored, bridges ABC to the OMP | HlyD, PrtE, AprE, LipC, CvaA; "HlyD family" / RND adaptor |
| Outer-membrane protein (OMP) | TolC-family channel forming the exit duct | TolC, PrtF, AprF; "TolC family" outer membrane efflux |

**Notes for curation:**
- The three genes are usually adjacent to the substrate gene (e.g. *hlyCABD*), but TolC
  is frequently encoded elsewhere and shared with efflux systems, count it as machinery
  only if the paper assigns it to this T1SS.
- *hlyC* (acyltransferase that activates HlyA) is a maturation enzyme, NOT a transport
  apparatus component; note it but do not count it as a core export component.
- Do NOT list the secreted substrate (HlyA, RTX toxin, protease, lipase, adhesin).
