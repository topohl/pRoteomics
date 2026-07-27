# WGCNA Modules And Group Effects

This folder builds WGCNA modules/supermodules and answers the primary biological
question:

**Which WGCNA modules and supermodules differ between CON, RES, and SUS across
microglia, neuron_soma, and neuron_neuropil?**

The primary inference layer is:

```powershell
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset <dataset> --level both
```

`01_WGCNA.r` builds networks and exports QC/descriptive screens. Its
module-trait and condition/eigengene heatmaps are useful for exploration, but
they are not the final group-effect inference.

## Recommended Run Order

```powershell
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset <dataset>
Rscript 06_modules_WGCNA/00_wgcna_identity_contract.R --dataset <dataset>
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset <dataset> --level both
Rscript 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r --dataset <dataset>
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset <dataset>
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset <dataset>
Rscript 06_modules_WGCNA/08_wgcna_publication_figures.R --dataset microglia
Rscript 06_modules_WGCNA/09_microglia_neuropil_independence.R --dataset microglia
Rscript 06_modules_WGCNA/10_module_complex_architecture.r --dataset <dataset>
Rscript 06_modules_WGCNA/11_module_robustness_sensitivity.r --dataset <dataset>
Rscript 06_modules_WGCNA/12_microglia_wgcna_nature_readiness_audit.R --animal-bootstrap 500
Rscript 06_modules_WGCNA/08b_microglia_wgcna_readiness_publication_figures.R --dataset microglia
Rscript 06_modules_WGCNA/13_wgcna_claim_readiness.R --dataset microglia
```

Run the final cross-dataset summary after all datasets are complete:

```powershell
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset all
```

## Canonical Identity Contract

`00_wgcna_identity_contract.R` is a read-only Phase 1 publication step. It
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
separate manuscript claims. Phase 2B changes the Stage 05 statistical contract,
so Stage 13 and every other downstream consumer remain on their prior semantics
until the atomic Phase 3 migration. Do not run those consumers against Phase 2B
outputs.

`05_module_supermodule_group_effects.r` is the Phase 2 quantitative boundary.
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

1. Primary adjusted WGCNA eigengene and supermodule models from `05_module_supermodule_group_effects.r`
2. Secondary module/program score robustness and behavior coupling from `03_score_module_activity.R`
3. DE/GSEA overlap support from `04_wgcna_de_gsea_overlap.r`
4. Descriptive module-trait and condition heatmaps from `01_WGCNA.r`

Stage 05 averages available left/right endpoint values within each
animal-spatial-unit before inference. The prespecified primary endpoint is
`SUS - RES` at `spatial_adjusted_global` /
`global_spatial_adjusted`, fitted with ML as
`eigengene ~ StressGroup + SpatialUnit + (1 | AnimalID)`. Modules and
higher-order multimodule supermodules have separate primary BH families.
Contextual global `RES - CON` and `SUS - CON` rows have contrast-specific
secondary BH families and never enter the primary family. Local contrasts use
exploratory localization families. Spatial heterogeneity is represented by one
nested-ML likelihood-ratio omnibus test per independent endpoint; conditional
contrasts are separate exploratory follow-ups.

A converged, full-rank, finite single-`AnimalID`-intercept fit whose random
intercept variance is at the boundary remains valid with
`model_stability_status = boundary_random_intercept_zero`,
`primary_model_stable = FALSE`, and an explicit warning. Rank deficiency,
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
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_interaction_conditional_followup.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_sensitivity.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_left_right_concordance.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_downstream_consumer_migration_audit.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_contract_status.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/WGCNA_group_effect_legacy_output_staleness_audit.csv`

Phase 2B does not regenerate Stage 05 figures, label outputs, marker-trait
correlations, selected interpretations, or Stage 06-13 products. Those
existing files are preserved and enumerated by hash as stale auxiliary outputs
until the Phase 3 consumer migration. The contract status deliberately reports
`downstream_compatible = FALSE`, `downstream_contract_status =
phase3_migration_required`, and `should_block_execution = TRUE`. This is an
advisory contract, not an enforced runtime block.

## Supermodule Annotation

`01_WGCNA.r` keeps data-driven eigengene clustering and sensitivity outputs.
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

`03_score_module_activity.R` is a secondary module/program scoring and
behavior-coupling layer. It preserves mapping trace, coverage QC, replicate QC,
and behavior-coupling exports. It records `PROTEOMICS_MODULE_DEFINITION_SOURCE`
or the dataset fallback in `module_score_run_metadata.csv`; primary WGCNA
eigengene group effects still come from `05_module_supermodule_group_effects.r`.
When the score source is `wgcna`, the script also exports secondary
supermodule eigengene score tables and supermodule directional robustness plots.

Two module-score source modes are useful for `neuron_neuropil`:

- `overlap`: curated biological programs from recurrent overlap proteins. This
  remains the default fallback for neuron neuropil module scoring.
- `wgcna`: data-driven co-expression modules from `01_WGCNA.r`. The default
  pipeline also runs this as an additional neuron-neuropil score pass so the
  curated-program plot and WGCNA module-score effect-size plot are both present.

`microglia` and `neuron_soma` default to `wgcna`. Score tables preserve
`ModuleID` as the stable technical join key and add `ModuleDisplayLabel` for
readable plotting/export labels.

Example explicit score-source runs:

```powershell
$env:PROTEOMICS_MODULE_DEFINITION_SOURCE = "wgcna"
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset microglia
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset neuron_soma

$env:PROTEOMICS_MODULE_DEFINITION_SOURCE = "overlap"
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset neuron_neuropil

$env:PROTEOMICS_MODULE_DEFINITION_SOURCE = "wgcna"
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset neuron_neuropil
Remove-Item Env:\PROTEOMICS_MODULE_DEFINITION_SOURCE
```

## Example PowerShell Commands

```powershell
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset microglia
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset microglia --level both
Rscript 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r --dataset microglia
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset microglia
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset microglia
```

```powershell
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset neuron_soma
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset neuron_soma --level both
Rscript 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r --dataset neuron_soma
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset neuron_soma
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset neuron_soma
```

```powershell
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset neuron_neuropil
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset neuron_neuropil --level both
Rscript 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r --dataset neuron_neuropil
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset neuron_neuropil
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset neuron_neuropil
```

Note: `06_module_spatial_networks.r` is a legacy optional helper. The canonical
spatial-network stage lives in `07_spatial_networks/`.
