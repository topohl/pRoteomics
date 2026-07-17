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
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset <dataset> --level both
Rscript 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r --dataset <dataset>
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset <dataset>
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset <dataset>
Rscript 06_modules_WGCNA/10_module_complex_architecture.r --dataset <dataset>
Rscript 06_modules_WGCNA/11_module_robustness_sensitivity.r --dataset <dataset>
```

Run the final cross-dataset summary after all datasets are complete:

```powershell
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset all
```

## Dataset Spatial Units

`microglia` and `neuron_soma` use `region`.

`neuron_neuropil` uses `region_layer`.

`05_module_supermodule_group_effects.r` records this in `SpatialUnitType` and
uses `spatial_unit` for the tested region or region-layer.

## Interpretation Hierarchy

Use evidence in this order:

1. Primary adjusted WGCNA eigengene and supermodule models from `05_module_supermodule_group_effects.r`
2. Secondary module/program score robustness and behavior coupling from `03_score_module_activity.R`
3. DE/GSEA overlap support from `04_wgcna_de_gsea_overlap.r`
4. Descriptive module-trait and condition heatmaps from `01_WGCNA.r`

The group-effect tables include ranked module and supermodule rows with model
metadata, FDRs, and `evidence_status`:
`robust_FDR`, `suggestive_FDR10`, `nominal_only`, `model_unstable`, or
`not_supported`.

Main outputs:

- `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv`
- `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv`
- `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_group_effects_interpretable.csv`
- `results/tables/06_modules_WGCNA/interpretable_summary/all/WGCNA_cross_dataset_supermodule_program_summary.csv`

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
mixed/unresolved; singletons are explicitly `singleton` with low confidence.
The sensitivity grid is `0.25, 0.35, 0.45, 0.50, 0.55, 0.65`.

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
