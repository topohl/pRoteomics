# WGCNA Workflow

`pipeline.yml` is the active registry. This document explains the interpretation
layers around the WGCNA scripts without changing the registered run order.

## Inputs

Core WGCNA construction stages dataset-scoped imputed expression workbooks with
the original protein-group annotation columns:

```text
data/processed/01_preprocessing/impute/*pgmatrix_imputed_<dataset>*.xlsx
```

Expression columns are keyed by canonical `ProteinGroupID`. The complete member
set is retained in `WGCNA_canonical_feature_table.csv` and
`WGCNA_protein_group_member_bridge.csv`. `ProteinID` is a deprecated downstream
alias of `ProteinGroupID`; `RepresentativeUniProt` and `FeatureDisplayLabel` are
display fields only. Quantitatively valid multi-gene, partially mapped and
unresolved groups remain in the network but are excluded from primary gene-level
GO annotation. Confirmed mixed-species or contaminant groups are excluded with an
explicit reason.

Identity construction is shared with Phase 1A. An existing `ProteinGroupID` or
opaque stable source feature ID is preserved/preferred. Membership fields such as
`Protein.Group` are canonicalized and matched through the sorted mapped-member set,
with a canonical source-membership discriminator only when it carries additional
stable feature information (for example an isoform-qualified group). This makes
member order irrelevant and lets entry-name contrast rows and accession-based
expression rows converge when their mapped membership is equivalent. Rows that
remain indistinguishable after these fields are applied fail validation; row order
and filenames are never used to repair identity.

Legacy module consumers may still read `UniProt`/`GeneSymbol`; those columns now
mirror representative display annotations and must not be treated as quantitative
keys. Such consumers remain compatibility-only until their own protein-group-aware
migration. The canonical handoff key is `ProteinGroupID`, with `MemberUniProts`
and `GeneSymbols` carrying complete mapped membership.

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

Supermodules are eigengene meta-modules/co-varying module blocks, constructed
with average linkage on `1 - signed Pearson module-eigengene correlation`.
They are not inferred from shared proteins: WGCNA modules partition proteins,
so hub/protein overlap is audit-only. The authoritative key is
`dataset + SupermoduleID`; display and biological labels are metadata and are
never used to join, group, deduplicate, or identify a cluster.

Supermodule GO support is restricted to `ModuleProteinSetType == "all"` and a
member module supports a term only at `p.adjust <= 0.05`. High confidence
requires support in every member module; medium requires at least two and at
least half. Other multi-module clusters remain mixed/unresolved, while
singletons are explicitly labelled `singleton`. Singleton structural status
does not lower the reviewed biological-label confidence inherited from the
member module. The full cut-height sensitivity grid is
`0.25, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65`.

Run:

```bash
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset <dataset>
```

Inspect:

- `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/`
- `results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/supermodules/`
- `.../supermodules/wgcna_supermodule_GO_term_support_audit.csv`
- `.../supermodules/wgcna_supermodule_biological_coherence.csv`
- `.../supermodules/wgcna_supermodule_validation_summary.csv`
- `.../supermodules/supermodule_clustering_sensitivity.csv`

Safe to rerun: no for routine interpretation updates. Rerun intentionally when
input matrices, WGCNA settings, or module construction choices change.
States predating the `protein_group_id_v1` feature-key contract, or states whose
ordered feature fingerprint differs, are rejected and require a full WGCNA rerun.

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

This Phase 2 script requires a publishable Phase 1 identity contract and uses
it as the sole module/supermodule identity and membership authority. It builds
an exact one-to-one bridge from canonical modules to frozen-state eigengene
columns, reconstructs canonical supermodule endpoints with the tested
`make_supermodule_eigengenes()` definition, and audits sample, AnimalID,
endpoint, PCA-orientation, model, contrast, and FDR provenance.

Neuron soma and microglia use Region as the spatial unit. Neuron neuropil uses
the observed Region-Layer unit; no generic layer main effect is fitted.
Repeated samples require an AnimalID random intercept. Failed or diagnostic
attempts never appear in the primary statistical tables, and Stage 05 does not
silently downgrade mixed models or emmeans contrasts to `lm` or t-tests.

Run:

```bash
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset <dataset> --level both
```

Inspect:

- `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_composition.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_endpoint_provenance.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_model_validation.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_sample_inclusion_audit.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_legacy_output_staleness_audit.csv`

Canonical publication is atomic and requires `--level both` so the
dataset-wide FDR family cannot change with a partial run.
`FDR_within_dataset_level` is BH within dataset and level;
`FDR_dataset_all_levels` is BH across both levels within a dataset.
`FDR_global` is the exact compatibility alias of the latter and is not a
cross-dataset adjustment.

Phase 2 intentionally leaves legacy Stage 05 figures, brackets, heatmaps,
selected interpretations, biological labels, and Stage 06-13 products
unchanged. Downstream migration is a separate phase.

Safe to rerun: yes, after the identity contract is publishable.

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

Automatic proposals remain centralized in `R/wgcna_labeling_utils.R`. For the
current microglia network, Stage 07 additionally loads
`config/wgcna_labels/microglia.csv`, validates exactly 13 modules and 9
supermodules against current Stage 01 membership, and makes reviewed labels
authoritative by `dataset + level + entity_id`. Biological process,
subcellular context, and ROI context remain separate. Automatic proposals and
GO evidence remain provenance and are fingerprinted. `final_plot_label` is an
exact compatibility alias for `canonical_plot_label`.

Run:

```bash
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset <dataset>
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset all
```

Inspect:

- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_interpretable_summary.xlsx`
- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_group_effects_interpretable.csv`
- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_label_audit.csv`
- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_label_candidates.csv`
- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_final_label_lookup.csv`
- `results/source_data/06_modules_WGCNA/interpretable_summary/<dataset>/`

Safe to rerun: yes.

### 08_wgcna_publication_figures.R

Role: final reviewed microglia all-supermodule publication layer.

This script reads the Stage 07 canonical lookup and reuses Stage 05 raw
eigengenes, estimates, confidence intervals, FDR values, model diagnostics,
sample/animal counts, and PC1 loadings. It does not fit models or alter WGCNA.
It produces architecture, 3 x 3 global eigengenes, SUS-RES forest, SUS-RES
spatial heatmap, multi-supermodule loadings (SM01, SM03, SM09 only), and the
supplementary all-contrast matrix as SVG/PDF plus one source CSV per figure.

Run:

```powershell
Rscript .\06_modules_WGCNA\08_wgcna_publication_figures.R --dataset microglia
```

Safe to rerun: yes.

### 09_microglia_neuropil_independence.R

Role: microglia-neuropil independence sensitivity.

This microglia-only script fits paired baseline and neuron-neuropil-adjusted
models after matching microglia ROI samples to neuron-neuropil covariates by
`AnimalID + Region`. Predeclared adjustment families live in
`config/microglia_neuropil_independence.yml` and are reported separately from
the exploratory strongest-Spearman match.

Run:

```bash
Rscript 06_modules_WGCNA/09_microglia_neuropil_independence.R --dataset microglia
```

Inspect:

- `results/tables/06_modules_WGCNA/microglia_neuropil_independence/microglia/microglia_neuropil_independence_effects.csv`
- `results/tables/06_modules_WGCNA/microglia_neuropil_independence/microglia/microglia_module_neuropil_independence_classification.csv`
- `results/reviewer_audit/microglia_neuropil_independence_claim_gate.csv`
- `results/reviewer_audit/microglia_neuropil_covariate_selection_audit.csv`
- `results/reviewer_audit/microglia_neuropil_independence_endpoint_scope_audit.csv`

Outputs identify module, targeted-signature, supermodule, or unavailable
endpoint scope explicitly. The current direct endpoints are module eigengenes
and targeted-signature scores. Module diagnostics can support module claims,
but cannot certify supermodule claims; absent a direct supermodule endpoint,
the latter report `no_direct_supermodule_independence_test`.

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

## Optional Microglia Readiness And Claim Handoff

`12_microglia_wgcna_nature_readiness_audit.R` is additive: it keeps the
historical primary network and memberships unchanged, recognizes AnimalID as
the repeated-measure unit, and treats standard WGCNA permutation preservation
as an unblocked ROI-row descriptive diagnostic rather than a claim gate. The
strict within-animal nonspatial sensitivity is diagnostic only; it cannot
downgrade a module’s primary architecture classification. Region-organized
covariance may be meaningful biology rather than technical failure.

`12b_finalize_microglia_wgcna_nature_readiness_audit.R` performs integrity and
report packaging only. `13_wgcna_claim_readiness.R` is the non-circular
manuscript handoff. Stages 05, 06, 07 and 12 do not consume Stage 13. The
reviewed registry is authoritative for final microglia labels; automatic GO and
marker names are candidate/provenance evidence only. Singleton compatibility
identities inherit the one member’s biological label and confidence, while
biological-process and ROI-context confidence remain distinct.

Stage 13 retains 22 technical rows for stable downstream joins. The six
singleton supermodule rows are `compatibility_alias` identities targeting
their one current module and are never independently manuscript-claimable.
SM01, SM03 and SM09 are the only higher-order-block claim identities. Stage 13
selects exactly `SUS - RES` / `spatial_adjusted_global` /
`global_spatial_adjusted` and carries the unchanged Stage 05 estimate, SE,
95% CI, p-value, within-dataset-level FDR, global FDR, model, replicate-unit,
source-file and source-key provenance. Claim status uses `FDR_global`.

Configured future-network defaults and historical selected values are separate
parameter concepts. The intended defaults remain neuron neuropil `0.55`,
neuron soma `0.35`, and microglia `0.45`. The frozen current microglia network
was generated with an explicit historical `0.40` override, so its saved
memberships and provenance remain associated with `0.40`. Stage 12 resolves
that selected value from the saved state, Stage 01 manifest, and clustering
sensitivity record rather than substituting the future `0.45` default.

```powershell
Rscript .\06_modules_WGCNA\12_microglia_wgcna_nature_readiness_audit.R --animal-bootstrap 500
Rscript .\06_modules_WGCNA\08b_microglia_wgcna_readiness_publication_figures.R --dataset microglia
Rscript .\06_modules_WGCNA\13_wgcna_claim_readiness.R --dataset microglia
```
