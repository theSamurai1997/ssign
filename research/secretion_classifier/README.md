# Secretion-system substrate classifier (research branch)

A side project exploring whether a multimodal classifier (frozen PLM
embedding + tabular features from MacSyFinder / DeepLocPro / SignalP /
DeepSecE) can do better than ssign's proximity-based substrate filtering.

This directory is for planning, scoping, data audits, and prototype
scripts. None of it is on the ssign user path yet. If something here
becomes shippable, it moves to `src/ssign_app/`.

## Status

- 2026-06-08: Project scoped. Gap-validation experiment designed.
  Background research agents running on PLM backbone choice, PU-learning
  recipe, data sourcing.

## Layout

```
research/secretion_classifier/
├── README.md                  (this file)
├── 00_overview.md             motivation + decision tree
├── 01_gap_validation.md       the recall-quantification experiment
├── 02_data_audit.md           what training data actually exists in 2026
├── 03_architecture.md         multimodal classifier design (open)
├── 04_negative_sets.md        PU learning + class imbalance (open)
├── 05_compute_plan.md         GPU-hour and disk budget (open)
├── data/
│   └── ground_truth/          per-genome curated substrate lists
├── scripts/                   prototype scripts
└── notes/                     working notes, agent outputs
```

## How to read this

Start with `00_overview.md`. The gap-validation experiment in
`01_gap_validation.md` is the first concrete action item. Architecture
choices in `03_*` and beyond are not committed yet.
