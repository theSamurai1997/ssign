# Why ssign misses effectors it COULD have found (reachable @±3 but not emitted)

These are ssign's true false negatives: the literature machinery is within ±3 genes (so proximity
*could* reach them), yet ssign did not emit them. 62 effectors (T1SS 5, T2SS 0, T3SS 40, T4SS 6,
T6SS 11; T3SS counted with detection enabled). The diagnosis is uniform and important.

## 50 of 62 (81%) are DETECTION failures, not prediction failures

ssign runs the secreted-protein predictors (DLP/DeepSecE/SignalP/PLM-Effector) **only on the ±3
neighborhood of a ssign-DETECTED secretion system**. So if MacSyFinder/TXSScan doesn't detect the
system, the effector is never handed to the predictors at all, every tool signal is blank, and it
can't be emitted. That is exactly what we see (`35_false_negatives.py`, fig `05`):

| system | reachable-missed | detection miss (no system detected) | processed but not emitted |
|---|---|---|---|
| T1SS | 5 | **5** | 0 |
| T3SS | 40 | **37** | 3 |
| T4SS | 6 | **6** | 0 |
| T6SS | 11 | 2 | **9** |

So for T1/T3/T4SS the misses are almost entirely the **detection step failing** — the literature
machinery is adjacent (that's why they're reachable), but ssign's MacSyFinder didn't call a system,
so the effector was never evaluated. T6SS is the exception (below).

## T1SS: all 5 are detection failures on textbook RTX toxins

The 5 T1SS misses are **HlyA, ApxIA, LtxA, LktA, ZapA** — canonical RTX toxins. ssign detected **zero
secretion systems** in all 5 genomes. The cleanest case is **HlyA** (E. coli α-hemolysin, the
textbook T1SS substrate) in the *complete* genome NZ_CP031766.1: the full operon is intact and
annotated — **hlyC – hlyA – hlyB – hlyD** (the HlyB/HlyD transporter is 1-2 genes from the toxin) —
yet ssign called no T1SS.

**Likely cause:** TXSScan's T1SS model requires the outer-membrane channel **TolC** alongside the
ABC transporter (HlyB) + MFP (HlyD). But TolC is a shared housekeeping gene encoded *elsewhere* in the
genome, not in the toxin operon, so MacSyFinder cannot assemble a complete-enough T1SS (wholeness
≥0.8) from the operon alone → no system → effector never evaluated. This is a systematic ssign T1SS
detection limitation for RTX toxins, and it explains why the operonic system (where proximity *should*
shine) still loses its most famous substrates. (Worth confirming against the TXSScan T1SS model
definition; flagged as a ssign-side follow-up.)

## T6SS is the informative exception: prediction/filter rejection

For T6SS, 9 of 11 misses are **processed_not_emitted** — ssign DID detect a T6SS and ran the predictors
on the protein, but cross-validation / substrate filtering did not emit it. These are the genuinely
interesting prediction-side misses (the protein was a candidate but rejected), as opposed to the
detection-side misses dominating the other systems. Worth a targeted drill-down on these 9.

## Implication

ssign's recall is **bottlenecked by secretion-system DETECTION** (MacSyFinder/TXSScan sensitivity),
not by the secreted-protein predictors — the predictors never even see most missed effectors. Two
takeaways:

1. **ssign-side fix:** T1SS detection needs to tolerate a non-co-localized TolC (or accept a
   HlyB+HlyD operon as a T1SS). This alone would recover the RTX-toxin class. T3SS/T4SS detection
   misses deserve the same scrutiny.
2. **Classifier argument:** a learned, per-protein model that is NOT gated on system detection +
   proximity would recover these — HlyA and friends are unmistakable secreted toxins that any sequence
   model flags. The detection+proximity gating is precisely what loses them. This is the recall-side
   complement to the precision argument.

Per-effector detail: `actual_per_effector.*.tsv` (filter `reachable_n3=true`, `ssign_call!=emitted_secreted`).
