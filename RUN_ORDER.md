# Recommended run order

This file gives a conservative execution order for the pRoteomics workflows. Individual projects may skip modules that are not relevant.

Canonical machine-readable outputs now belong under `data/processed/<module>/`. Publication tables, figures, source data, logs and reports belong under `results/tables`, `results/figures`, `results/source_data`, `results/logs` and `results/reports`.

## 1. Preprocessing and metadata harmonization

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
data/processed/02_id_mapping/mapped/forward/per_file/*.csv
```

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

`01_clusterProfiler.r` reads mapped contrast CSVs and writes:

```text
data/processed/04_differential_expression_enrichment/clusterProfiler/clusterProfiler_manifest.csv
results/source_data/04_differential_expression_enrichment/clusterProfiler/
results/figures/04_differential_expression_enrichment/clusterProfiler/
results/tables/04_differential_expression_enrichment/clusterProfiler/
results/logs/04_differential_expression_enrichment/clusterProfiler/
results/reports/04_differential_expression_enrichment/clusterProfiler/
```

`02_compareGO.r` consumes only manifest-selected clusterProfiler outputs and the mapped/log2FC files listed in that manifest. It writes:

```text
data/processed/04_differential_expression_enrichment/compareGO/compareGO_input_manifest.csv
results/tables/04_differential_expression_enrichment/compareGO/
results/figures/04_differential_expression_enrichment/compareGO/
results/source_data/04_differential_expression_enrichment/compareGO/
results/logs/04_differential_expression_enrichment/compareGO/
results/reports/04_differential_expression_enrichment/compareGO/
```

## 5. Cell-type enrichment

```text
05_celltype_enrichment_EWCE/
```

Typical outputs include EWCE result tables, cell-type enrichment plots, and measured-proteome-aware enrichment summaries.

## 6. Module and WGCNA analyses

```text
06_modules_WGCNA/
```

Typical outputs include module assignments, module scores, module preservation statistics, and trait correlations.

## 7. Spatial network analyses

```text
07_spatial_networks/
```

Typical outputs include region/layer network edge tables, differential networks, and bootstrap stability summaries.

## 8. Behavior and physiology coupling

```text
08_behavior_physio_coupling/
```

Typical outputs include correlation tables and figures linking proteomic modules or networks with behavioral/physiological phenotypes.

## 9. PRIDE / ProteomeXchange and journal export

```r
source("09_export_pride_journal/01_make_pride_manifest.R")
source("09_export_pride_journal/02_make_sample_metadata.R")
source("09_export_pride_journal/03_make_supplementary_tables.R")
source("09_export_pride_journal/04_validate_pride_submission.R")
source("09_export_pride_journal/05_make_methods_summary.R")
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
