# T2SS apparatus brief (Type II secretion system)

**Mechanism (plain):** a two-step system. Substrates first cross the inner membrane
folded-up via Sec or Tat, then the T2SS pushes them through the outer membrane using a
short piston-like "pseudopilus". The substrate (e.g. a toxin, lipase, protease) is NOT
machinery.

**Core machinery (~12-15 genes, named Gsp C-O/S; each system uses a 3-letter prefix):**

| Unified | Role | Common system prefixes/synonyms |
|---|---|---|
| GspD | secretin: the outer-membrane channel (12-15mer ring) | XcpQ, EpsD, OutD, PulD, XpsD, HxcQ, LspD, ExeD, EtpD, GspD |
| GspC | inner-membrane platform, links to secretin | XcpP, EpsC, OutC, PulC |
| GspE | cytoplasmic secretion ATPase | XcpR, EpsE, OutE, PulE, XpsE |
| GspF | integral IM platform component | XcpS, EpsF, OutF, PulF |
| GspL | IM platform, binds ATPase | XcpY, EpsL, OutL, PulL |
| GspM | IM platform | XcpZ, EpsM, OutM, PulM |
| GspG | major pseudopilin (builds the piston) | XcpT, EpsG, OutG, PulG, XpsG |
| GspH, GspI, GspJ, GspK | minor pseudopilins (pseudopilus tip/initiation) | XcpU-X, EpsH-K, OutH-K |
| GspO | prepilin peptidase (matures pseudopilins; often shared with T4P) | XcpA/PilD, VcpD, PulO |
| GspN, GspS | accessory (not in all systems) | EpsN, PulS (pilotin) |
| GspA, GspB | accessory IM (some systems, e.g. Exe/Vibrio) | ExeA, ExeB |

**Naming:** the corpus `sys_instance_label` already gives the prefix and gene range,
e.g. "Xcp (xcpP-Z)", "Eps (epsC-N)", "Pul (pulC-O)", "Out (outC-S)", "Lsp (lspD-K)",
"Hxc (hxcQ-Z)", "Xps (xpsE-N)", "Txc (txcQ-Z)", "Yts1 (yts1A-N)". Treat the whole named
operon (C through O/S) as the apparatus.

**Notes for curation:**
- The prepilin peptidase (GspO/PilD) is often encoded with the T4 pilus system and shared;
  count it if the paper ties it to this T2SS.
- Do NOT list the secreted substrate, nor the upstream Sec/Tat translocon (that is general
  IM transport, not the T2SS apparatus).
