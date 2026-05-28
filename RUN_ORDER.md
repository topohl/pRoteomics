# Recommended run order

This file gives a conservative execution order for the pRoteomics workflows. Individual projects may skip modules that are not relevant.

Canonical machine-readable outputs now belong under `data/processed/<module>/`. Publication tables, figures, source data, logs and reports belong under `results/tables`, `results/figures`, `results/source_data`, `results/logs` and `results/reports`.

## 1. Preprocessing and metadata harmonization

Dataset-scoped runs can be launched by setting `PROTEOMICS_DATASET` once:

```powershell
$env:PROTEOMICS_DATASET = "microglia"
Rscript 01_preprocessing/03_gct_extractR.r --dry-run
Rscript 02_id_mapping/01_MapThatProt_batch.r --dry-run
Rscript 04_differential_expression_enrichment/01_clusterProfiler.r --dry-run
Rscript 04_differential_expression_enrichment/02_compareGO.r --dry-run
```

Or with the launcher:

```powershell
Rscript run_dataset_pipeline.R --dataset microglia --dry-run
Rscript run_dataset_pipeline.R --dataset microglia
```

Valid dataset families are `neuron_neuropil`, `neuron_soma`, and `microglia`. The shared `R/dataset_config.R` helper resolves `PROTEOMICS_DATASET` first, with backward-compatible fallback to `PROTEOMICS_COMPARISON` and `PROTEOMICS_GCT_COMPARISON`.

```text
01_preprocessing/
```

Expected stable outputs:

```text
data/metadata/sample_metadata_clean.tsv
data/processed/01_preprocessing/protein_matrix_raw.tsv
data/processed/01_preprocessing/protein_matrix_normalized.tsv
data/processed/01_preprocessing/protein_matrix_imputed.tsv
```

The exact matrix names may differ by project. For PRIDE preparation, the most important file is `sample_metadata_clean.tsv`, including a `raw_file_name` column.

## 2. UniProt and identifier mapping

```text
02_id_mapping/
```

Expected outputs:

```text
data/processed/02_id_mapping/protein_id_mapping.tsv
data/processed/02_id_mapping/mapped_protein_matrix.tsv
data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/*.csv
```

The current canonical contrast handoff is:

```powershell
$env:PROTEOMICS_DATASET = "microglia"
Rscript 01_preprocessing/03_gct_extractR.r --dry-run
Rscript 02_id_mapping/01_MapThatProt_batch.r --dry-run
```

`03_gct_extractR.r` writes split contrast CSVs to:

```text
data/processed/01_preprocessing/gct_extractR/<comparison>/forward/
data/processed/01_preprocessing/gct_extractR/<comparison>/reverse/
```

Run the GCT extraction and ID mapping once per biological dataset family, for example `neuron_neuropil`, `neuron_soma`, and `microglia`. `01_MapThatProt_batch.r` defaults to `PROTEOMICS_DATASET=neuron_neuropil` and `PROTEOMICS_MAP_DIRECTION=forward`; set `PROTEOMICS_DATASET` for each family. It writes clusterProfiler-ready files to:

```text
data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/
```

Set `PROTEOMICS_MAP_DIRECTION=reverse` only when intentionally producing reverse contrasts.

## 3. QC and exploratory analysis

```text
03_qc_exploration/
```

Typical outputs include PCA plots, variance summaries, protein/peptide QC, and sample-level QC reports.

## 4. Differential expression and enrichment

```text
04_differential_expression_enrichment/
```

Typical outputs include differential abundance tables, clusterProfiler/GSEA outputs, compareGO outputs, and publication-style enrichment figures.

The manifest-backed order is:

```r
source("04_differential_expression_enrichment/01_clusterProfiler.r")
source("04_differential_expression_enrichment/02_compareGO.r")
```

Dry-run validation without launching analyses:

```bash
Rscript 04_differential_expression_enrichment/01_clusterProfiler.r --dry-run
Rscript 04_differential_expression_enrichment/02_compareGO.r --dry-run
```

Set `analysis.dataset` in `config/clusterProfiler_config.yml` and `dataset` in `config/compareGO_config.yml` before each dataset family run. Leave mapped paths empty in the configs to use the dataset-aware defaults. `01_clusterProfiler.r` reads mapped contrast CSVs and writes:

```text
data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv
results/source_data/04_differential_expression_enrichment/clusterProfiler/<dataset>/
results/figures/04_differential_expression_enrichment/clusterProfiler/<dataset>/
results/tables/04_differential_expression_enrichment/clusterProfiler/<dataset>/
results/logs/04_differential_expression_enrichment/clusterProfiler/<dataset>/
results/reports/04_differential_expression_enrichment/clusterProfiler/<dataset>/
```

`02_compareGO.r` requires the manifest `dataset` column, filters manifest rows by `dataset`, and consumes only manifest-selected clusterProfiler outputs plus the mapped/log2FC files listed in that manifest. It writes:

```text
data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv
results/tables/04_differential_expression_enrichment/compareGO/<dataset>/
results/figures/04_differential_expression_enrichment/compareGO/<dataset>/
results/source_data/04_differential_expression_enrichment/compareGO/<dataset>/
results/logs/04_differential_expression_enrichment/compareGO/<dataset>/
results/reports/04_differential_expression_enrichment/compareGO/<dataset>/
```

## 5. Cell-type enrichment

```text
05_celltype_enrichment_EWCE/
```

Typical outputs include EWCE result tables, cell-type enrichment plots, and measured-proteome-aware enrichment summaries.

Phase 3 canonicalized `01_EWCE_E9.r` so it reads the processed proteomics matrix and sample metadata through `R/paths.R` and writes to:

```text
data/processed/05_celltype_enrichment_EWCE/EWCE_E9/
results/tables/05_celltype_enrichment_EWCE/EWCE_E9/
results/figures/05_celltype_enrichment_EWCE/EWCE_E9/
results/source_data/05_celltype_enrichment_EWCE/EWCE_E9/
results/logs/05_celltype_enrichment_EWCE/EWCE_E9/
results/reports/05_celltype_enrichment_EWCE/EWCE_E9/
```

Dry-run:

```bash
Rscript 05_celltype_enrichment_EWCE/01_EWCE_E9.r --dry-run
```

## 6. Module and WGCNA analyses

```text
06_modules_WGCNA/
```

Typical outputs include module assignments, module scores, module preservation statistics, and trait correlations.

Phase 3 canonicalized the safer downstream/helper scripts:

```bash
Rscript 06_modules_WGCNA/02_module_spatial_networks.r --dry-run
Rscript 06_modules_WGCNA/03_overlap_modules.r --dry-run
Rscript 06_modules_WGCNA/91_module_score.r --dry-run
```

`91_module_score.r` is the canonical active module-score script. By default it consumes overlap/GSEA-derived modules. To score WGCNA modules from `01_WGCNA.r`, run with `PROTEOMICS_MODULE_DEFINITION_SOURCE=WGCNA`; the script then consumes `WGCNA_modules_long.xlsx` and, when available, `wgcna_final_model_state.rds` for eigengene scores. `90_module_score_v0.0.1.r` is retained only as an older reference and should not be used in the active run order.

`01_WGCNA.r` uses `R/paths.R` and writes an input manifest/hash table plus run manifest, but it still consumes precombined variancePartition-style matrices (`male.data.xlsx` and `sample_info.xlsx`). Generate those upstream or set `PROTEOMICS_WGCNA_EXPR_XLSX` and `PROTEOMICS_WGCNA_META_XLSX` explicitly. It exports stable color-based WGCNA module definitions under `results/tables/06_modules_WGCNA/01_WGCNA/modules/`.

## 7. Spatial network analyses

```text
07_spatial_networks/
```

Typical outputs include region/layer network edge tables, differential networks, and bootstrap stability summaries.

Phase 4 canonicalized the spatial-network producer and Phase 3 canonicalized downstream network scripts:

```bash
Rscript 07_spatial_networks/01_network_spatial_relations.r --dry-run
Rscript 07_spatial_networks/02_differential_networks.r --dry-run
Rscript 07_spatial_networks/03_bootstrap_network_stability.r --dry-run
Rscript 07_spatial_networks/04_bootstrap_differential_network_stability.r --dry-run
Rscript 07_spatial_networks/05_bootstrap_differential_network_figures.r --dry-run
Rscript 07_spatial_networks/06_chord_diagram.r --dry-run
```

These scripts expect the canonical spatial object:

```text
data/processed/07_spatial_networks/network_spatial_relations/network_spatial_relations_objects.rds
```

`01_network_spatial_relations.r` now writes that canonical object while preserving the downstream RDS structure expected by the other network and behavior scripts.

## 8. Behavior and physiology coupling

```text
08_behavior_physio_coupling/
```

Typical outputs include correlation tables and figures linking proteomic modules or networks with behavioral/physiological phenotypes.

Phase 3 canonicalized `02_network_behavior_coupling.r`. It expects behavior inputs under:

```text
data/external/behavior/
```

Dry-run:

```bash
Rscript 08_behavior_physio_coupling/02_network_behavior_coupling.r --dry-run
```

## 9. PRIDE / ProteomeXchange and journal export

```r
source("09_export_pride_journal/01_make_pride_manifest.R")
source("09_export_pride_journal/02_make_sample_metadata.R")
source("09_export_pride_journal/03_make_supplementary_tables.R")
source("09_export_pride_journal/04_validate_pride_submission.R")
source("09_export_pride_journal/05_make_methods_summary.R")
```

Dry-run validation:

```bash
Rscript 09_export_pride_journal/01_make_pride_manifest.R --dry-run
Rscript 09_export_pride_journal/02_make_sample_metadata.R --dry-run
Rscript 09_export_pride_journal/03_make_supplementary_tables.R --dry-run
Rscript 09_export_pride_journal/04_validate_pride_submission.R --dry-run
Rscript 09_export_pride_journal/05_make_methods_summary.R --dry-run
```

Before running this step, copy local large files into the gitignored `pride_submission/` folder as needed:

```text
pride_submission/metadata/
pride_submission/processed_data/
pride_submission/supplementary_tables/
pride_submission/methods/
pride_submission/manifests/
pride_submission/validation/
```

This step generates:

```text
pride_submission/manifests/pride_file_manifest.tsv
pride_submission/metadata/sample_metadata.tsv
pride_submission/metadata/sdrf_like_metadata.tsv
pride_submission/validation/validation_report.tsv
pride_submission/methods/methods_summary.md
```

Review `validation_report.tsv` and resolve all `FAIL` entries before upload.
