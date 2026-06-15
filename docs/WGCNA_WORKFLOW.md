# WGCNA Workflow

`pipeline.yml` is the active registry. This document explains the interpretation
layers around the WGCNA scripts without changing the registered run order.

## Inputs

Core WGCNA construction consumes dataset-scoped mapped protein contrast files:

```text
data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv
```

Downstream WGCNA interpretation may also consume:

- `config/marker_panels/wgcna_reference_marker_sets.csv`
- `results/tables/03_qc_exploration/05_empirical_roi_marker_discovery/empirical_roi_marker_sets.csv`
- `data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv`
- `data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv`
- `results/tables/04_differential_expression_enrichment/neuropil_reference_annotation/<dataset>/`
- `data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx`

Microglia is interpreted as microglia-enriched ROI/local microenvironment
proteomics. It is not purified microglia, and downstream annotation should not
turn neuropil overlap or local microenvironment signal into purified immune
activation.

## Script Layers

### 01_WGCNA.r

Role: network/module construction.

This script builds the WGCNA network, detects modules, calculates eigengenes,
creates stable module definitions, clusters supermodules, and writes module and
supermodule construction outputs. It is the only registered WGCNA stage that
recomputes core WGCNA state.

Run:

```bash
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset <dataset>
```

Inspect:

- `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/`
- `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/`

Safe to rerun: no for routine interpretation updates. Rerun intentionally when
input matrices, WGCNA settings, or module construction choices change.

### 03_score_module_activity.R

Role: module scoring.

This script scores existing module definitions against dataset matrices and
metadata. It supports source-scoped definitions such as `wgcna`, `overlap`, or
`custom`. It does not rebuild WGCNA modules.

Run:

```bash
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset <dataset>
```

Inspect:

- `results/tables/06_modules_WGCNA/module_score/<dataset>/`

Safe to rerun: yes.

### 05_module_supermodule_group_effects.r

Role: group-effect modelling.

This script models module and supermodule effects between CON, RES, and SUS,
including spatial-unit and spatial-adjusted/global summaries where estimable. It
uses existing module/supermodule state and records model warnings/fallbacks.

Run:

```bash
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset <dataset>
```

Inspect:

- `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv`

Safe to rerun: yes.

### 06_annotate_module_microenvironment.r

Role: biological annotation.

This script annotates existing modules and supermodules with marker, empirical
ROI, neuropil-reference, microenvironment, and label rationale evidence. For
microglia, labels should remain conservative ROI/local microenvironment labels
unless immune/phagolysosomal/complement/inflammatory evidence is explicit.

Run:

```bash
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset <dataset>
```

Inspect:

- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv`
- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv`
- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_display_label_audit.csv`
- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_targeted_signature_overlap_details.csv`

Safe to rerun: yes.

For microglia, targeted signature evidence is split into auditable classes:
microglia-enriched empirical, microglia-enriched reference-supported, curated
microglia-relevant program, neuropil-shared, and ambiguous. Curated programs
are not claim-ready by themselves; single generic mitochondrial overlaps are
cautionary mitochondrial/oxidative-stress annotations rather than purified
microglia specificity. The module annotation table keeps legacy targeted-overlap
columns and adds explicit driver/caution columns such as
`targeted_signature_primary_driver`, `targeted_signature_driver_class`,
`targeted_signature_driver_signature`, `targeted_signature_driver_padj`,
`targeted_signature_driver_NES`, and
`targeted_signature_driver_overlap_proteins`, plus
`n_unique_targeted_signatures`, `n_unique_targeted_overlap_proteins`, and
`curated_program_overlap_warning`.

### 07_wgcna_interpretable_summary.r

Role: manuscript-facing summary.

This script joins group effects, annotations, and optional DE/GSEA overlap into
readable tables, source data, and plots. It should preserve both broad plotting
labels and full annotation labels so figures stay compact without hiding the
underlying rationale.

Run:

```bash
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset <dataset>
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset all
```

Inspect:

- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_interpretable_summary.xlsx`
- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_group_effects_interpretable.csv`
- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_label_audit.csv`
- `results/source_data/06_modules_WGCNA/interpretable_summary/<dataset>/`

Safe to rerun: yes.

## Interpretation Checklist

- Use `01_WGCNA.r` outputs to describe module construction, not final biological
  claims.
- Use `05_module_supermodule_group_effects.r` for primary module/supermodule
  group-effect evidence.
- Use `06_annotate_module_microenvironment.r` and
  `07_wgcna_interpretable_summary.r` to explain biological context and figure
  labels.
- For microglia ROI results, separate targeted microglia signature overlap,
  microglia-associated ROI signal, shared microglia-neuropil ROI signal, and
  genuine immune/phagolysosomal/complement/inflammatory evidence.
- Do not subtract neuropil evidence or reinterpret ROI/local microenvironment
  labels as purified microglia regulation.
