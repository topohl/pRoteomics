# WGCNA Modules And Group Effects

This folder builds WGCNA modules/supermodules and answers the primary biological
question:

**Which WGCNA modules and supermodules differ between CON, RES, and SUS across
microglia, neuron_soma, and neuron_neuropil?**

The primary inference layer is:

```powershell
Rscript analysis/05_wgcna/test_module_phenotypes.R --dataset <dataset> --level both
```

`build_wgcna_modules.R` builds networks and exports QC/descriptive screens. Its
module-trait and condition/eigengene heatmaps are useful for exploration, but
they are not the final group-effect inference.

## Representative GO comparison heatmaps

`render_module_go_heatmaps.R` turns the Stage 01 module GO results
into three complementary views: a broad representative module heatmap, a
broad representative supermodule heatmap, and a focused manuscript-style
supermodule dot matrix. The default is Biological Process (`--ontology BP`);
use `--ontology all` to generate each view separately for BP, MF, and CC. The
focused view selects up to three FDR-supported terms per supermodule by default
(`--focused-terms-per-supermodule 3`) and compares their deduplicated union
across every supermodule. Selection walks the existing evidence rank in order
and conservatively skips an exact GO ancestor/descendant only when its
significant-row contributing-gene Jaccard overlap is at least 0.50, or any term
pair with near-identical contributing genes (Jaccard at least 0.80). The
candidate-level focused selection audit records selected, redundancy-skipped,
and display-limit terms; no nonsignificant terms are added.
Module cells are capped `-log10(BH FDR)` from the `all` module-protein ORA
results, and white means the term is not FDR-significant. Supermodule cells are
the mean module score across their member modules, so a term shared by more
members is stronger. Representative supermodule terms are selected from the
full module GO result before aggregation, preventing the module panel's
top-term cutoff from hiding a recurrent theme. This deliberately summarizes
existing module evidence and does **not** perform or imply a new pooled
supermodule ORA test. In the focused panel, dot colour is the mean member-module
capped `-log10(BH FDR)` and dot size is the fraction of member modules with BH
FDR <= 0.05; absent dots mean no member module supports that term at the cutoff.
The contributing genes are the slash-delimited Entrez IDs already retained in
Stage 01 `geneID`, unioned only across FDR-significant member-module rows. If
that field or installed GO hierarchy is unavailable, the focused audit records
the deterministic hierarchy-only or evidence-rank-only fallback mode.

After the three dataset-focused sources have been generated consistently, run
`render_module_go_heatmaps.R --dataset all --ontology BP` to create the
manuscript-facing coordinated figure under the same output family in `all/`.
It vertically aligns Neuropil, Soma, and Microglia-enriched ROI panels while
retaining each dataset's own GO rows and SM columns. Colour and size use common,
unnormalized scales and one shared legend. Exact GO `TermID` recurrence across
datasets is recorded descriptively; it is not a cross-dataset enrichment test,
meta-analysis, combined p-value, FDR, or convergence claim.

The figures use a Nature-style double-column width (7.2 inches/183 mm), vector
SVG/PDF export, compact module IDs, and module columns grouped under their
data-driven supermodule. The accompanying chart-contract CSV records the
question, visual structure, palette, and displayed measure for each panel.

Outputs are under
`results/{tables,figures}/06_modules_WGCNA/01b_module_supermodule_GO_heatmaps/<dataset>/`.

## Recommended Run Order

```powershell
Rscript analysis/05_wgcna/build_wgcna_modules.R --dataset <dataset>
Rscript analysis/05_wgcna/render_module_go_heatmaps.R --dataset <dataset>
Rscript analysis/05_wgcna/build_module_identity_contract.R --dataset <dataset>
Rscript analysis/05_wgcna/test_module_phenotypes.R --dataset <dataset> --level both
Rscript analysis/05_wgcna/compare_module_enrichment_overlap.R --dataset <dataset>
Rscript analysis/05_wgcna/annotate_module_microenvironment.R --dataset <dataset>
Rscript analysis/05_wgcna/summarize_module_interpretation.R --dataset <dataset>
Rscript analysis/05_wgcna/render_module_figures.R --dataset microglia
Rscript analysis/05_wgcna/test_microglia_neuropil_independence.R --dataset microglia
Rscript analysis/05_wgcna/summarize_module_complex_architecture.R --dataset <dataset>
Rscript analysis/05_wgcna/audit_module_robustness.R --dataset <dataset>
Rscript analysis/05_wgcna/audit_microglia_module_claims.R --animal-bootstrap 500
Rscript analysis/05_wgcna/render_microglia_module_figures.R --dataset microglia
Rscript analysis/05_wgcna/audit_module_claim_readiness.R --dataset microglia
```

Run the final cross-dataset summary after all datasets are complete:

```powershell
Rscript analysis/05_wgcna/summarize_module_interpretation.R --dataset all
```

## Canonical Identity Contract

`build_module_identity_contract.R` is a read-only Phase 1 publication step. It
does not recompute WGCNA or modify Stage 01 outputs. It publishes current
module and supermodule identity under:

`results/tables/06_modules_WGCNA/identity_contract/<dataset>/`

Module identity comes from the frozen state and the Stage 01 downstream module
definitions. Only the exact identifier forms `ME#RRGGBB`, `#RRGGBB`, and
`WGCNA_#RRGGBB` may be normalized to `WGCNA_#RRGGBB`; labels, row order,
eigengene order, and approximate matching are prohibited.

For neuronal datasets, supermodule identity comes from the active
provenance-selected `wgcna_supermodule_eigengene_clusters.csv` only when its
cut height, output-manifest hash, sensitivity metadata, module coverage, and
one-to-one membership all agree. For microglia, the verified current
`wgcna_module_supermodule_annotation.csv` remains authoritative under the same
coverage checks. Exact member-module composition, not a repeated `SMxx`
string, defines supermodule identity.

Stage 05-13, publication-score, and circular-atlas files are compatibility
audit inputs only. The identity-contract script never repairs them or uses
their memberships or biological labels as authority. A required validation
failure writes diagnostic hashes, validation, compatibility, and status
outputs, then refuses to publish entity or membership contracts.

## Dataset Spatial Units

`microglia` and `neuron_soma` use `region`.

`neuron_neuropil` uses `region_layer`.

For microglia, `config/wgcna_labels/microglia.csv` is the authoritative reviewed
biological-label source. Automatic GO/marker labels remain candidate/provenance
evidence. Singleton supermodule compatibility identities inherit their member
module label and confidence; singleton status is structural metadata, not a
biological-confidence penalty. Stage 12 is optional manuscript-readiness audit;
its strict nonspatial sensitivity is diagnostic only. Stage 12b packages the
audit and does not make scientific inferences. Stage 13 is the canonical,
non-circular claim-readiness handoff for manuscript/global consumers.
Stage 13 retains all stable technical identities, but singleton supermodule
IDs are compatibility aliases for their one member module and cannot form
separate manuscript claims. The atomic Phase 3 source migration is complete:
Stage 07 now provides the sole claim-facing inferential handoff, and downstream
consumers use its tier-specific FDR family, claim gate, and exact Stage 05
source key.

`test_module_phenotypes.R` is the Phase 2 quantitative boundary.
It requires a publishable Phase 1 identity contract and treats that contract as
the only supermodule-membership authority. The frozen-state module eigengenes
are joined to canonical modules through an exact, audited identifier bridge.
For neuronal modules, the only permitted normalization is `ME#RRGGBB` or
`#RRGGBB` to `WGCNA_#RRGGBB`; microglia uses frozen-state `WGCNA_mNN`
metadata.

Stage 05 records `SpatialUnitType` and uses `spatial_unit` for the tested
region or region-layer. Neuron soma and microglia use Region; neuron neuropil
uses the observed Region-Layer unit and never fits a generic layer main effect.
The canonical scopes are `within_spatial_unit`, `spatial_adjusted_global`, and
`stress_by_spatial_interaction`.

## Interpretation Hierarchy

Use evidence in this order:

1. Primary adjusted WGCNA eigengene and supermodule models from `test_module_phenotypes.R`
2. Secondary module/program score robustness and behavior coupling from `score_module_activity.R`
3. DE/GSEA overlap support from `compare_module_enrichment_overlap.R`
4. Descriptive module-trait and condition heatmaps from `build_wgcna_modules.R`

Stage 05 first averages technical source rows within each hemisphere, then
gives the one or two observed hemispheres equal weight within each
animal-spatial-unit without imputation. It exports both hemisphere provenance
and a canonical `aggregated_row_sha256` content hash. Under the
`sha256_utf8_lf_v1` hash contract, fields are encoded as UTF-8, separated by
literal LF bytes, terminated by one final LF byte, and hashed without BOM or
native text-mode translation. The prespecified primary endpoint is
`SUS - RES` at `spatial_adjusted_global` /
`global_spatial_adjusted`, fitted with ML as
`eigengene ~ StressGroup + SpatialUnit + (1 | AnimalID)`. Modules and
higher-order multimodule supermodules have separate primary BH families.
Contextual global `RES - CON` and `SUS - CON` rows have contrast-specific
secondary BH families and never enter the primary family. Local contrasts use
exploratory localization families. Spatial heterogeneity is represented by one
nested-ML likelihood-ratio omnibus test per independent endpoint; conditional
contrasts are separate exploratory follow-ups.

`model_diagnostic_scope` makes the interpretation of top-level diagnostics
explicit. Ordinary named contrasts use `single_fitted_model`. Interaction
omnibus rows use `composite_reduced_full`: their validity and stability are
deterministic composite properties, their top-level raw single-fit variance
fields are typed `NA`, and all numerical diagnostics remain in `reduced_*` and
`full_*`. Conditional interaction follow-ups use `full_interaction_model`, so
their top-level diagnostics describe the full interaction model.

A converged, full-rank, finite single-`AnimalID`-intercept fit whose variance
ratio (`random_intercept_variance / residual_variance`) is at or below `1e-4`
remains valid with
`model_stability_status = boundary_random_intercept_zero`,
`primary_model_stable = FALSE`, and an explicit warning. The independently
recorded `lme4::isSingular(..., tol = 1e-4)` result is diagnostic only;
disagreement sets `diagnostic_review_required = TRUE` without replacing the
ratio classification or changing inference eligibility. Rank deficiency,
optimizer or non-boundary convergence failure, non-estimability, nonfinite
results, sample-contract failure, and singular complex random structures are
invalid. No t-test fallback is produced.

`FDR_primary_global`, `FDR_secondary_global`,
`FDR_interaction_omnibus`, and `FDR_local_exploratory` encode the prespecified
biological families. `FDR_conservative_all_tests` is a reviewer sensitivity
across independent primary and secondary tests. Deprecated `FDR_global`, where
present, is an exact alias of that conservative sensitivity and is not the
claim gate. Singleton supermodule compatibility rows inherit the member-module
statistics and diagnostics but all FDR fields are `NA`.

Stage 05 establishes model validity, statistical support, and hypothesis
independence only. Every row has
`manuscript_claim_ready = not_assessed_stage05`; the status output separately
records `stage05_output_status = phase2b_statistical_outputs_complete` and
`publication_status = not_assessed_stage05`.

Main outputs:

- `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_endpoint_provenance.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_model_validation.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_sample_inclusion_audit.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_animal_spatial_unit_values.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_hemisphere_values.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_interaction_conditional_followup.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_sensitivity.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_left_right_concordance.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_downstream_consumer_migration_audit.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_contract_status.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_legacy_output_staleness_audit.csv`

Phase 2B does not regenerate Stage 05 figures, label outputs, marker-trait
correlations, selected interpretations, or Stage 06-13 products. Those
existing files are preserved and enumerated by hash as stale auxiliary outputs
until an authorized output refresh. The generated contract status predates the
completed Phase 3 source migration and deliberately reports
`downstream_compatible = FALSE`, `downstream_contract_status =
phase3_migration_required`, and `should_block_execution = TRUE`. This is an
advisory generated artifact, not an indication that claim-facing source
consumers still select legacy FDRs.

## Supermodule Annotation

`build_wgcna_modules.R` keeps data-driven eigengene clustering and sensitivity outputs.
These are eigengene meta-modules (co-varying module blocks) constructed by
average linkage on `1 - signed Pearson module-eigengene correlation`; they are
not protein-overlap clusters. WGCNA modules partition proteins, so protein or
hub overlap is retained only as a partition-integrity audit and never supports
supermodule coherence or confidence.

Every computational operation uses `dataset + SupermoduleID`; biological and
display labels are metadata. GO naming uses `ModuleProteinSetType == "all"`
and counts a member module only when `p.adjust <= 0.05`. High confidence
requires the same term in every member module. Medium confidence requires at
least two and at least half. Other multi-module clusters are conservatively
mixed/unresolved; singletons are explicitly `singleton`. That structural status
does not lower the reviewed biological-label confidence inherited from the
member module.
Future-network dataset defaults are `0.55` for neuron neuropil, `0.35` for
neuron soma, and `0.45` for microglia. These defaults are distinct from the
selected value recorded for a saved network. The frozen current microglia
network was generated using a historical explicit `0.40` override; its
memberships remain tied to `0.40`. The sensitivity grid is
`0.25, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65`.

Manual labels absent from the active dataset are retained for audit but marked
with `present_in_dataset = FALSE`, `annotation_scope =
manual_absent_from_dataset`, and `manual_label_status =
manual_label_absent_from_dataset`.

## Secondary Module Scores

`score_module_activity.R` is a secondary module/program scoring and
behavior-coupling layer. It preserves mapping trace, coverage QC, replicate QC,
and behavior-coupling exports. It records `PROTEOMICS_MODULE_DEFINITION_SOURCE`
or the dataset fallback in `module_score_run_metadata.csv`; primary WGCNA
eigengene group effects still come from `test_module_phenotypes.R`.
When the score source is `wgcna`, the script also exports secondary
supermodule eigengene score tables and supermodule directional robustness plots.

Two module-score source modes are useful for `neuron_neuropil`:

- `overlap`: curated biological programs from recurrent overlap proteins. This
  remains the default fallback for neuron neuropil module scoring.
- `wgcna`: data-driven co-expression modules from `build_wgcna_modules.R`. The default
  pipeline also runs this as an additional neuron-neuropil score pass so the
  curated-program plot and WGCNA module-score effect-size plot are both present.

`microglia` and `neuron_soma` default to `wgcna`. Score tables preserve
`ModuleID` as the stable technical join key and add `ModuleDisplayLabel` for
readable plotting/export labels.

Example explicit score-source runs:

```powershell
$env:PROTEOMICS_MODULE_DEFINITION_SOURCE = "wgcna"
Rscript analysis/05_wgcna/score_module_activity.R --dataset microglia
Rscript analysis/05_wgcna/score_module_activity.R --dataset neuron_soma

$env:PROTEOMICS_MODULE_DEFINITION_SOURCE = "overlap"
Rscript analysis/05_wgcna/score_module_activity.R --dataset neuron_neuropil

$env:PROTEOMICS_MODULE_DEFINITION_SOURCE = "wgcna"
Rscript analysis/05_wgcna/score_module_activity.R --dataset neuron_neuropil
Remove-Item Env:\PROTEOMICS_MODULE_DEFINITION_SOURCE
```

## Example PowerShell Commands

```powershell
Rscript analysis/05_wgcna/build_wgcna_modules.R --dataset microglia
Rscript analysis/05_wgcna/test_module_phenotypes.R --dataset microglia --level both
Rscript analysis/05_wgcna/compare_module_enrichment_overlap.R --dataset microglia
Rscript analysis/05_wgcna/annotate_module_microenvironment.R --dataset microglia
Rscript analysis/05_wgcna/summarize_module_interpretation.R --dataset microglia
```

```powershell
Rscript analysis/05_wgcna/build_wgcna_modules.R --dataset neuron_soma
Rscript analysis/05_wgcna/test_module_phenotypes.R --dataset neuron_soma --level both
Rscript analysis/05_wgcna/compare_module_enrichment_overlap.R --dataset neuron_soma
Rscript analysis/05_wgcna/annotate_module_microenvironment.R --dataset neuron_soma
Rscript analysis/05_wgcna/summarize_module_interpretation.R --dataset neuron_soma
```

```powershell
Rscript analysis/05_wgcna/build_wgcna_modules.R --dataset neuron_neuropil
Rscript analysis/05_wgcna/test_module_phenotypes.R --dataset neuron_neuropil --level both
Rscript analysis/05_wgcna/compare_module_enrichment_overlap.R --dataset neuron_neuropil
Rscript analysis/05_wgcna/annotate_module_microenvironment.R --dataset neuron_neuropil
Rscript analysis/05_wgcna/summarize_module_interpretation.R --dataset neuron_neuropil
```

Note: `build_module_spatial_networks.R` is a legacy optional helper. The canonical
spatial-network stage lives in `07_spatial_networks/`.
