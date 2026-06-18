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
- `config/marker_panels/microenvironment_marker_panels.csv`
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
Supplemental microenvironment panels are versioned outside code in
`config/marker_panels/microenvironment_marker_panels.csv`; each marker records
source type, reference, allowed use, claim role, and caution notes. The run
manifest records the config path and hash.

Run:

```bash
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset <dataset>
```

Inspect:

- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_biological_annotation.csv`
- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_biological_annotation.csv`
- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_supermodule_display_label_audit.csv`
- `results/tables/06_modules_WGCNA/module_annotation/<dataset>/WGCNA_module_targeted_signature_overlap_details.csv`
- `results/reviewer_audit/wgcna_microenvironment_threshold_sensitivity.csv`
- `results/reviewer_audit/wgcna_label_confidence_audit.csv`
- `results/reviewer_audit/wgcna_annotation_source_audit.csv`

Safe to rerun: yes.

Microenvironment annotation is threshold-aware. The primary marker fraction
threshold remains `0.10`, while reviewer sensitivity is reported at `0.05`,
`0.10`, and `0.20`. Threshold-unstable, low-confidence, or insufficient
annotations are retained as annotation/context only and should not be promoted
to biological claims without passing the separate model, animal, robustness,
and evidence-independence gates.

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

### 08_microglia_neuropil_independence.R

Role: microglia-neuropil independence sensitivity.

This microglia-only script fits paired baseline and neuron-neuropil-adjusted
models after matching microglia ROI samples to neuron-neuropil covariates by
`AnimalID + Region`. Predeclared adjustment families live in
`config/microglia_neuropil_independence.yml` and are reported separately from
the exploratory strongest-Spearman match.

Run:

```bash
Rscript 06_modules_WGCNA/08_microglia_neuropil_independence.R --dataset microglia
```

Inspect:

- `results/tables/06_modules_WGCNA/microglia_neuropil_independence/microglia/microglia_neuropil_independence_effects.csv`
- `results/tables/06_modules_WGCNA/microglia_neuropil_independence/microglia/microglia_module_neuropil_independence_classification.csv`
- `results/reviewer_audit/microglia_neuropil_independence_claim_gate.csv`
- `results/reviewer_audit/microglia_neuropil_covariate_selection_audit.csv`

Only `predeclared_primary` and `predeclared_secondary` rows can make a
neuropil-independence gate eligible, and only when the baseline microglia
contrast has a claim-relevant primary effect first. By default this means
`primary_effect_status = FDR_pass` under the configured threshold in
`config/microglia_neuropil_independence.yml`; nominal-only effects are
diagnostic unless the config explicitly changes that policy. Null or
non-significant primary effects are classified as diagnostic/inconclusive
no-primary-effect rows and cannot enable microglia-specific wording.

The `exploratory_best_spearman` row is retained to diagnose whether a
best-matched neuron-neuropil covariate explains the signal, but because it is
selected using endpoint correlation it is not claim-enabling.

Independence classifications are conservative: `neuropil_independent`,
`partially_neuropil_adjusted`, `neuropil_sensitive`,
`inconclusive_low_power`, `inconclusive_missing_match`,
`diagnostic_no_primary_effect`, `inconclusive_no_primary_effect`,
`mixed_or_covariate_sensitive`, and `exploratory_only`. Near-zero baseline
effects are protected by `effect_before_near_zero` and
`percent_attenuation_reliable`; unstable attenuation does not support claim
eligibility. `neuropil_sensitive` blocks or downgrades microglia-specific
interpretation. Inconclusive/no-primary rows remain diagnostic/contextual.
Microglia ROI or local-microenvironment wording may remain suggestive only when
it stays explicitly ROI/local-microenvironment and no required independence
gate fails.

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
- Use predeclared neuropil-adjustment rows, not exploratory best-Spearman rows,
  when evaluating whether a microglia-specific or cell-intrinsic phrase is
  claim-gate eligible.
- Do not subtract neuropil evidence or reinterpret ROI/local microenvironment
  labels as purified microglia regulation.
# WGCNA Claim-Gate Inference Notes

Primary WGCNA inference for group effects is the module or supermodule eigengene model exported by `06_modules_WGCNA/05_module_supermodule_group_effects.r`. The claim-grade columns in `module_group_effects.csv` and `supermodule_group_effects.csv` record the model family, formula, emmeans status, rank/singularity diagnostics, animal random-effect use, and biological replicate unit used for each row.

Fallback tests are diagnostic only. If emmeans fails or a two-group t-test substitute is emitted, the numerical estimate is retained for review, but `primary_model_stable = FALSE`, `claim_allowed_model = FALSE`, and `model_downgrade_reason` includes `diagnostic_only_model_fallback`.

The biological replicate unit is explicit. `animal_level_status` is one of `animal_level`, `repeated_sample_mixed_model`, `sample_level_or_unclear`, `insufficient_animals`, or `missing_animal_id`; `pseudoreplication_guard` records the corresponding reviewer gate status.

Robustness for WGCNA claims is summarized in `results/reviewer_audit/wgcna_robustness_claim_gate.csv`. The biological claims table may allow WGCNA group-effect claims only when model, animal-level, QC, and robustness gates all pass under the strict claim gate.

Reviewer-facing WGCNA claim rows separate model fit from statistical strength. `model_fit_status` describes whether the model fit itself was claim-grade, while `statistical_evidence_status` distinguishes FDR-passing rows from nominal-only or FDR-failing rows. `primary_model_status` is kept for compatibility and should be read as a gate shorthand, not a diagnosis by itself.

`WGCNA_module_group_effect` and `WGCNA_supermodule_group_effect` claim rows are the direct outputs of the group-effect models. `WGCNA_module` claim rows derived from `WGCNA_module_evidence_rank.csv` are ranked evidence summaries for review and context; unless their model, robustness, animal-level, and evidence-independence support is complete, they should remain suggestive/contextual rather than primary manuscript evidence.

Blocked WGCNA rows are retained for audit but must not be phrased as evidence-supported manuscript claims. `blocked_claim_wording_audit.csv` checks this wording contract, and `wgcna_claim_source_audit.csv` summarizes WGCNA claim source, model-fit status, statistical-evidence status, robustness gate, and claim-use class.
