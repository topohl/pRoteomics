# Recommended run order

This file gives a conservative execution order for the pRoteomics workflows. Individual projects may skip modules that are not relevant.

## 1. Preprocessing and metadata harmonization

```text
01_preprocessing/
```

Expected stable outputs:

```text
results/analysis_ready/sample_metadata_clean.tsv
results/analysis_ready/protein_matrix_raw.tsv
results/analysis_ready/protein_matrix_normalized.tsv
results/analysis_ready/protein_matrix_imputed.tsv
```

The exact matrix names may differ by project. For PRIDE preparation, the most important file is `sample_metadata_clean.tsv`, including a `raw_file_name` column.

## 2. UniProt and identifier mapping

```text
02_id_mapping/
```

Expected outputs:

```text
results/analysis_ready/protein_id_mapping.tsv
results/analysis_ready/mapped_protein_matrix.tsv
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

## 9. PRIDE / ProteomeXchange package preparation

```r
source("09_pride_submission/00_make_pride_package.R")
```

Before running this step, copy local large files into the gitignored `PRIDE_package/` folder:

```text
PRIDE_package/01_raw/
PRIDE_package/02_search_results/
PRIDE_package/03_processed_quantification/
PRIDE_package/04_downstream_analysis/
```

This step generates:

```text
PRIDE_package/00_metadata/sdrf_proteomics.tsv
PRIDE_package/00_metadata/pride_file_manifest.tsv
PRIDE_package/00_metadata/checksum.md5
PRIDE_package/00_metadata/validation_report.tsv
```

Review `validation_report.tsv` and resolve all `FAIL` entries before upload.
