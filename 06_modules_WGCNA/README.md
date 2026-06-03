# WGCNA Modules And Downstream Interpretation

`01_WGCNA.r` remains the primary network-construction script. The downstream
scripts added here test group effects, add biological annotation, and create
human-readable summaries without changing the primary WGCNA protein universe.

Run per dataset:

```bash
Rscript 06_modules_WGCNA/01_WGCNA.r --dataset <dataset>
Rscript 03_qc_exploration/04b_import_reference_marker_sources.r
Rscript 03_qc_exploration/05_empirical_roi_marker_discovery.r
Rscript 03_qc_exploration/06_wgcna_marker_trait_export.r --dataset <dataset>
Rscript 04_differential_expression_enrichment/04_neuropil_contamination_annotation.r --dataset microglia
Rscript 06_modules_WGCNA/05_module_supermodule_group_effects.r --dataset <dataset>
Rscript 06_modules_WGCNA/06_annotate_module_microenvironment.r --dataset <dataset>
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset <dataset>
Rscript 06_modules_WGCNA/07_wgcna_interpretable_summary.r --dataset all
```

Valid datasets are `neuron_neuropil`, `neuron_soma`, and `microglia`.

For microglia-enriched ROI interpretation, also run the neuropil reference
annotation when DE/GSEA manifests exist:

```bash
Rscript 04_differential_expression_enrichment/04_neuropil_contamination_annotation.r --dataset microglia
```

Important constraints:

- Microglia ROI samples are not treated as purified microglia.
- Neuropil intensities and logFC values are not subtracted.
- Neuropil-overlap proteins are not removed from primary microglia WGCNA.
- Marker and neuropil evidence is used only for annotation, reporting, plotting,
  and sensitivity interpretation.
- Marker sets are loaded from `config/marker_panels/wgcna_reference_marker_sets.csv`,
  then empirical ROI sets if present, then legacy hard-coded panels only as fallback.
- Reference marker import is offline by default; live downloads require
  `PROTEOMICS_REFERENCE_MARKERS_ALLOW_DOWNLOAD=true`.

Main answers:

- Changed modules: `results/tables/06_modules_WGCNA/group_effects/<dataset>/module_group_effects.csv`
- Changed supermodules: `results/tables/06_modules_WGCNA/group_effects/<dataset>/supermodule_group_effects.csv`
- Interpretable supermodule summary: `results/tables/06_modules_WGCNA/interpretable_summary/<dataset>/WGCNA_supermodule_group_effects_interpretable.csv`
- Microglia ROI support class: `results/tables/06_modules_WGCNA/module_annotation/microglia/WGCNA_module_biological_annotation.csv` and `WGCNA_supermodule_biological_annotation.csv`

Group-effect outputs use explicit effect scopes:

- `within_spatial_unit`: `eigengene ~ StressGroup + Sex + Batch`
- `spatial_adjusted_global`: `eigengene ~ StressGroup + SpatialLabel + Sex + Batch`
- `stress_by_spatial_interaction`: `eigengene ~ StressGroup * SpatialLabel + Sex + Batch`

When repeated `AnimalID` values are present and `lmerTest` is installed, the
group-effect script adds `(1 | AnimalID)` and records `model_type =
lmerTest_lmer`; otherwise it uses `lm` and records the fallback. Marker-trait
correlations are written as annotation-only sensitivity outputs under
`results/tables/06_modules_WGCNA/group_effects/<dataset>/`.

Microglia ROI module classes are conservative:
`microglia_supported`, `microglia_state_or_activation_supported`,
`shared_microenvironment`, `neuropil_sensitive`,
`other_cellular_or_vascular_sensitive`, or `ambiguous`.

Supermodule naming is auditable. `01_WGCNA.r` exports data-driven IDs/labels,
curated labels, final labels, label source/confidence/rationale, manual-review
flags, and clustering sensitivity across cut heights 0.25-0.45.

Optional existing downstream:

```bash
Rscript 06_modules_WGCNA/03_score_module_activity.R --dataset <dataset>
Rscript 06_modules_WGCNA/04_wgcna_de_gsea_overlap.r --dataset <dataset>
Rscript 06_modules_WGCNA/06_module_spatial_networks.r --dataset <dataset> --module-definition-source wgcna
```

Note: `06_module_spatial_networks.r` already used the `06_` prefix before the
new annotation script was added. It was not renamed to avoid breaking existing
calls.
