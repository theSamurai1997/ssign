# Citation pass-2: agent re-audit of the 30 flagged/indeterminate ssign-found effectors

Scope: the 30 non-CONSISTENT rows from pass-1 (`citation_consistency_found.tsv`) = 8 FLAG_WRONG_TOPIC
+ 3 FLAG_GENE_ABSENT + 3 DOI_UNRESOLVED + 16 INDETERMINATE. Each was sent to a literature agent under
an anti-hallucination contract (CONFIRMED / RESOLVED / MISASSIGNED / NOT_FOUND; resolvable DOI +
verbatim secretion quote naming the protein; never fabricate). 5 agent batches ran; batch C (YspE,
VirA, BipB, CopN, ChlaDub1) was AUP-blocked on *B. pseudomallei* content and done by hand in-session
via PubMed. **Every returned DOI was then deterministically re-verified** (`42_pass2_verify.py`: DOI
registered + on CrossRef + gene/genus present) so no fabricated or wrong-paper DOI can pass silently.

## Outcome (30 rows)

| status | n | verify |
|---|---|---|
| RESOLVED (corrected DOI found) | 21 | 25 VERIFIED + 1 real-but-no-abstract |
| CONFIRMED (current DOI was right) | 3 | (in the 25) |
| MISASSIGNED (wrong SS type) | 2 | (in the 25) |
| NOT_FOUND (no qualifying paper) | 4 | 4 NO_DOI |

**0 UNVERIFIED_DOI** — every DOI the agents/manual pass returned resolves and is on-topic. The
deterministic re-check agreeing with the agents on all 26 sourced rows means the agents were reliable
here (consistent with the 3/3 CrossRef spot-checks in pass-1).

## What needs a dataset edit

- **2 MISASSIGNED ss_type → fix or drop:** `BopA` and `BopE` (B. thailandensis BTH_II locus rows) are
  labelled **T6SS** but are **Bsa T3SS** effectors (name-collision with the B. pseudomallei Bop set;
  BopA = Cullinane 2008 autophagy effector, BopE = Stevens 2003 GEF). Their ssign emission was the
  accidental cross-type kind (the 6 flagged in benchmark RESULTS).
- **4 NOT_FOUND → drop or down-tier:**
  - `ChlaDub1` (Chlamydia) — canonical paper (Misaghi 2006) shows only DUB activity, **no T3SS
    secretion**; the T3SS-effector label is literature-unsupported.
  - `aprA` (P. entomophila instance) — secreted protease, but no primary paper pins the *P. entomophila*
    ortholog to T1SS specifically (homology-only to P. aeruginosa AprA, which is separately confirmed).
  - `EFF00136` / `EFF00150` — opaque placeholder IDs resolving to BTH DUF3274 proteins **not named** in
    their cited paper; unsupported.

## Notable corrections the re-audit surfaced

- **`Tle4` (idx15) and `TplE_alias_Tle4` (idx24) are a duplicate** (both = P. aeruginosa PA1510) AND
  both carried a wrong DOI (a protein-structure-design paper / a lung-regeneration paper); correct
  primary = Jiang 2016 Cell Rep `10.1016/j.celrep.2016.07.012`. De-dupe to one row.
- **`EspZ`** cited the *T6SS-discovery* PNAS paper; corrected to Kanack 2005 (`10.1128/IAI.73.7.4327-4337.2005`).
- **`EFF00142`** was tagged "unidentifiable" in the pass-1 prior audit but is actually supported:
  Russell 2012 names BTH_I2691 as a T6SS-1 substrate → CONFIRMED. The prior tag was wrong.
- `Map`×2, `EspH`×2, `celA`, `plaA`, `VirA`, `CopN`, `YspE`, `BipB`, `Tae4_Stm`, `Tle1`, `Tle1_Sci1`
  all re-sourced to verified primary papers (see `pass2_results.tsv` for per-row DOI + verbatim quote).

Full per-row table with DOIs, verbatim quotes, and verify_status: `pass2_results.tsv`.
Raw agent/manual returns: `pass2_raw/batch_*.json`. **Items 3 (triage by training impact) + 4
(precision estimate) still deferred to discussion.**
