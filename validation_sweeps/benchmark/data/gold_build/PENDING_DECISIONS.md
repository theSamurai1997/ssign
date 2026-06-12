# Pending decisions, Phase 0a gold set

## D-A: 5 "real effector, wrong metadata" rows  [RESOLVED 2026-06-10: DROP all 5]

Teo's call: drop all 5 (conservative, matches audit). Encoded as `EXPLICIT_DROP_LOCI`
in `scripts/01_build_gold_set.py`; they appear in `03_biology_dropped.tsv`.


These survived the pattern-based biology drop and sit in `03_repair_queue.tsv`.
Unlike the apparatus/immunity drops, each IS a real secreted effector; only its
instance/type/name binding is wrong. The audit ("Biology misclassifications, NOT
rescuable") said drop them, but it also recorded the correct value for each, so
fixing is possible. Instance binding matters here because the ceiling measures
distance to the effector's OWN system instance.

| gene | locus (table) | type | audit finding | fix-per-audit |
|---|---|---|---|---|
| BspA | BAB1_1671 | T4SS | row is actually BspE, not BspA (Myeni 2013) | rename gene -> BspE, keep |
| RARP-1 | RPR_RS04075 | T4SS | secreted by TolC/T1SS, not VirB/T4SS (Kaur 2012) | reclassify to T1SS, or drop from T4SS |
| TseT_H1 | PA0820 | T6SS | real effector but H2 @ PA3907, not H1 @ PA0820 (Burkinshaw 2018) | re-bind to H2 + change locus to PA3907 |
| TseM_Mn_binding | BTH_II1883 | T6SS | real effector, mislabeled T6SS-1 (chr I); actually chr II instance | re-bind to the chr-II instance |
| RbsB_metal_scav | BTH_II1884 | T6SS | real effector, mislabeled T6SS-1 (chr I); actually chr II instance | re-bind to the chr-II instance |

Default if undecided: **drop** (conservative, matches the audit's stated action).
Teo's call: drop all 5, or fix-per-audit (keeps 5 real effectors, costs a small
manual correction each). TseT_H1 fix changes the genome position, not just a label.

## Note: name-vs-locus matching

`prtA @ XCV3671` (Xanthomonas) is NOT a biology error; the audit's prt drop was
Dickeya-specific (`Dda3937_*`, already FAIL). Caught only because we checked by
locus, not gene name. Any future drop rule must key on locus_tag, not gene name.
