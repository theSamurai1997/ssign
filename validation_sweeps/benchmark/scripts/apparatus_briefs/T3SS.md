# T3SS apparatus brief (Type III secretion system / injectisome)

**Mechanism (plain):** a molecular syringe that injects effector proteins directly from
the bacterial cytoplasm into a host cell through a hollow needle. ~20+ core genes. The
injected effectors are NOT machinery; neither are the chaperones that escort them.

**Unified Sct nomenclature** (Hueck/Wagner) with the major family synonyms. Plant-pathogen
(Hrp/Hrc) and animal-pathogen (Ysc/Inv/Mxi/Ssa/Esc) systems are homologous.

| Unified | Role | Synonyms (Salmonella SPI1 / Yersinia Ysc / EPEC LEE / Shigella / plant Hrp) |
|---|---|---|
| SctC | outer-membrane secretin ring | InvG / YscC / EscC / MxiD / HrcC |
| SctJ | inner ring (lipoprotein) | PrgK / YscJ / EscJ / MxiJ / HrcJ |
| SctD | inner MS-ring | PrgH / YscD / EscD / MxiG / HrpQ |
| SctR | export gate | SpaP / YscR / EscR / Spa24 / HrcR |
| SctS | export gate | SpaQ / YscS / EscS / Spa9 / HrcS |
| SctT | export gate | SpaR / YscT / EscT / Spa29 / HrcT |
| SctU | export-gate autoprotease (switch) | SpaS / YscU / EscU / Spa40 / HrcU |
| SctV | major export-apparatus ring | InvA / YscV / EscV / MxiA / HrcV |
| SctN | ATPase | InvC / YscN / EscN / Spa47 / HrcN |
| SctO | ATPase stalk | InvI / YscO / EscO / Spa13 / HrpO |
| SctL | ATPase stator / C-ring link | OrgB / YscL / EscL / MxiN / HrpE |
| SctQ | cytoplasmic C-ring / sorting platform | SpaO / YscQ / EscQ / Spa33 / HrcQ |
| SctK | C-ring base | OrgA / YscK / EscK / MxiK |
| SctI | inner rod | PrgJ / YscI / EscI / MxiI / HrpB |
| SctF | needle subunit | PrgI / YscF / EscF / MxiH / HrpA (Hrp pilus) |
| SctW | gatekeeper | InvE / YopN+TyeA / SepL / MxiC / HrpJ |
| tip / translocon | host-membrane pore (translocators) | SipB/C/D / YopB/D+LcrV / EspA/B/D / IpaB/C/D / HrpF,HrpK |
| SctP | needle-length ruler | InvJ / YscP / EscP / Spa32 / HrpP |
| pilotin | secretin assembly | InvH / YscW / — / MxiM |

**Family labels in the corpus `sys_instance_label`:** Bsa, Hrp, Hrp1, LEE, SPI-1 (Inv/Spa/Prg),
SPI-2 (Ssa), Ysc, Mxi/Spa, Psc, Esc. Treat the whole named locus as the apparatus.

**Notes for curation:**
- Flagellar export genes (Flh/Fli/Flg) are a homologous but SEPARATE system, do not count
  them as injectisome machinery even though SctU/SctV resemble FlhB/FlhA.
- Do NOT list secreted effectors, translocated chaperones (e.g. SycE, CesT, IpgC), or
  transcriptional regulators (e.g. HrpL, InvF, ExsA) as apparatus.
- The translocon proteins (SipB/C, IpaB/C, EspB/D, YopB/D) are apparatus (they form the
  host-membrane pore), not effectors, despite being secreted.
