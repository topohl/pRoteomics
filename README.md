![Logo](https://github.com/user-attachments/files/28182224/bf2400e5f1d00757e1bd57e3abdc823cadab95b39e3705765fb82bdde967186b.tiff)

# pRoteomics

pRoteomics is a repository of R-based workflows for spatial and systems-level proteomics analyses.

The repository has evolved from a primarily clusterProfiler/GSEA workflow into a broader analysis framework integrating:

- spatial proteomics
- differential enrichment analysis
- EWCE cell-type enrichment
- WGCNA and module analyses
- spatial region/layer network analysis
- bootstrap-based network stability analyses
- behavior and physiology coupling
- PRIDE/ProteomeXchange package preparation support

The current focus is high-resolution spatial proteomics across hippocampal region/layer structures combined with systems-level network and behavioral analyses.

---

# Repository structure

```text
01_preprocessing/
02_id_mapping/
03_qc_exploration/
04_differential_expression_enrichment/
05_celltype_enrichment_EWCE/
06_modules_WGCNA/
07_spatial_networks/
08_behavior_physio_coupling/
09_pride_submission/
09_export_pride_journal/
90_testing/
99_deprecated/
```

Canonical data and result folders:

```text
data/raw/
data/metadata/
data/external/
data/processed/<module>/
results/figures/<module>/
results/tables/<module>/
results/source_data/<module>/
results/logs/<module>/
results/reports/<module>/
pride_submission/
```

Use `R/paths.R` for repo-relative paths. Local overrides should go in `config/*.local.yml` or environment variables such as `PROTEOMICS_PROJECT_ROOT`; do not commit machine-specific paths.

---

# Pipeline overview

```text
raw proteomics matrices + metadata
→ preprocessing / metadata harmonization
→ UniProt and ID mapping
→ QC and exploratory structure analysis
→ differential enrichment and GSEA
→ EWCE cell-type enrichment
→ WGCNA and module analyses
→ spatial network analyses
→ bootstrap stability analyses
→ behavior and physiology coupling
→ PRIDE / ProteomeXchange package metadata generation
```

Recommended manuscript workflow:

```bash
for ds in neuron_neuropil neuron_soma microglia; do
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage all --dry-run
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage core
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage qc
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage enrichment
  Rscript run_dataset_pipeline.R --dataset "$ds" --stage modules
done
Rscript 09_export_pride_journal/06_make_biological_claims_table.R
```

Resolve dry-run failures before launching full analyses. The claim table is an
evidence index for manuscript drafting; it uses `NA` rather than inventing
missing statistics.

---

# Main analysis modules

## 01_preprocessing

Purpose:
- metadata formatting
- matrix harmonization
- imputation
- merged metadata generation

Representative scripts:
- `01_impute.r`
- `03_gct_extractR.r`
- `04_format_metadata.r`

---

## 02_id_mapping

Purpose:
- UniProt mapping
- clusterProfiler-compatible ID conversion
- WGCNA-compatible identifier harmonization

Representative scripts:
- `01_MapThatProt.r`
- `02_MapThatProt_batch.r`

The canonical contrast handoff is now `01_preprocessing/03_gct_extractR.r` to `02_id_mapping/01_MapThatProt_batch.r`, producing clusterProfiler-ready mapped files at `data/processed/02_id_mapping/mapped/<dataset>/forward/per_file/`. Run this per dataset family, such as `neuron_neuropil`, `neuron_soma`, or `microglia`.

---

## 03_qc_exploration

Purpose:
- PCA
- variance partitioning
- rank abundance analysis
- protein/peptide QC

Representative scripts:
- `00_dataset_qc_report.r`
- `05_pca_confounding_qc.r`
- `06_variance_partitioning.r`

`00_dataset_qc_report.r` is the canonical one-stop dataset QC report. It reads
the dataset-resolved matrix and metadata, then writes missingness, imputation
footprint, sample/protein counts, PCA, metadata structure, abundance
distribution, and outlier tables/figures under
`results/*/03_qc_exploration/00_dataset_qc_report/<dataset>/`.

---

## 04_differential_expression_enrichment

Purpose:
- GO enrichment
- GSEA
- pathway comparison
- publication-style enrichment figures

Representative scripts:
- `01_clusterProfiler.r`
- `02_compareGO.r`
- `03_biological_program_summary.r`
- `04_compare_pathways.r`

`01_clusterProfiler.r` now writes a manifest at:

```text
data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv
```

`02_compareGO.r` consumes that manifest instead of recursively discovering arbitrary CSVs. It filters by the required `dataset` column so neuron neuropil, neuron soma, and microglia runs do not overwrite or mix with each other. This preserves ontology, comparison, route category, route unit, simplification state, plot-used status, input hashes and config hashes across the clusterProfiler to compareGO handoff.

Dry-run checks are available:

```bash
Rscript 04_differential_expression_enrichment/01_clusterProfiler.r --dry-run
Rscript 04_differential_expression_enrichment/02_compareGO.r --dry-run
Rscript 04_differential_expression_enrichment/03_biological_program_summary.r --dataset neuron_neuropil --dry-run
```

`03_biological_program_summary.r` conservatively groups redundant GO/GSEA terms
into seven manuscript programs: RNA/RNP processing, ribosome/translation,
mitochondria/OXPHOS/metabolism, proteostasis/ubiquitin/protein folding,
synapse/vesicle organization, cytoskeleton/motility, and
development/patterning. It writes tidy, wide, heatmap-ready, and SVG heatmap
outputs for each dataset.

---

## 05_celltype_enrichment_EWCE

Purpose:
- EWCE analyses
- spatial cell-type interpretation
- measured-proteome-aware enrichment workflows

Representative scripts:
- `01_EWCE_E9.r`

`01_EWCE_E9.r` now uses canonical module folders and supports `--dry-run`.

---

## 06_modules_WGCNA

Purpose:
- WGCNA
- module scoring
- module preservation
- overlap-based module generation

Representative scripts:
- `01_WGCNA.r`
- `05_wgcna_de_gsea_overlap.r`
- `91_module_score.r`

`01_WGCNA.r` is the canonical dataset-aware WGCNA entrypoint. It supports
`--dataset neuron_neuropil`, `--dataset neuron_soma`, `--dataset microglia`,
and `--dry-run`; stages WGCNA inputs from canonical processed matrices and
metadata when available; and writes dataset-scoped module contracts, logs,
state, tables, figures, and manifests under
`data/processed/06_modules_WGCNA/01_WGCNA/<dataset>/` and
`results/*/06_modules_WGCNA/01_WGCNA/<dataset>/`. Manually staged
variancePartition-style workbooks are legacy fallback inputs only via
`PROTEOMICS_WGCNA_EXPR_XLSX` and `PROTEOMICS_WGCNA_META_XLSX`.

---

## 07_spatial_networks

Purpose:
- anatomical region/layer relationship networks
- differential network analysis
- bootstrap network validation
- chord/network visualization

Representative scripts:
- `01_network_spatial_relations.r`
- `02_differential_networks.r`
- `03_bootstrap_network_stability.r`

Phase 4 canonicalized `01_network_spatial_relations.r`, and Phase 3 canonicalized `02_differential_networks.r`, `03_bootstrap_network_stability.r`, `04_bootstrap_differential_network_stability.r`, `05_bootstrap_differential_network_figures.r`, and `06_chord_diagram.r`. The producer writes the spatial object at `data/processed/07_spatial_networks/network_spatial_relations/network_spatial_relations_objects.rds`, which downstream network and behavior scripts consume.

---

## 08_behavior_physio_coupling

Purpose:
- connect proteomics and network structure with behavior and physiology
- movement/stress score integration
- systems-level phenotype coupling

Representative scripts:
- `01_correlate_proteomics_with_behavior.r`
- `02_network_behavior_coupling.r`

`02_network_behavior_coupling.r` now reads behavior files from `data/external/behavior/`, writes canonical outputs under `results/*/08_behavior_physio_coupling/network_behavior_coupling/`, and supports `--dry-run`.

---

## 09_pride_submission

Purpose:
- create a local PRIDE/ProteomeXchange package skeleton
- generate SDRF-style sample-to-file metadata
- generate a file manifest
- calculate MD5 checksums
- validate whether raw files, search outputs, processed results, and metadata are traceable

Main command:

```r
source("09_pride_submission/00_make_pride_package.R")
```

This creates or updates local files under the gitignored `PRIDE_package/` folder:

```text
PRIDE_package/00_metadata/sdrf_proteomics.tsv
PRIDE_package/00_metadata/pride_file_manifest.tsv
PRIDE_package/00_metadata/checksum.md5
PRIDE_package/00_metadata/validation_report.tsv
```

Large raw/vendor mass spectrometry files should be uploaded to PRIDE, not committed to GitHub.

Detailed instructions are in:

```text
09_pride_submission/README_PRIDE.md
```

## 09_export_pride_journal

Purpose:
- generate PRIDE-ready manifests
- stage sample metadata and SDRF-like metadata
- stage publication supplementary tables
- validate deposition readiness
- write methods/provenance summaries

Main scripts:

```text
09_export_pride_journal/01_make_pride_manifest.R
09_export_pride_journal/02_make_sample_metadata.R
09_export_pride_journal/03_make_supplementary_tables.R
09_export_pride_journal/04_validate_pride_submission.R
09_export_pride_journal/05_make_methods_summary.R
```

Each export script supports `--dry-run` to validate expected folders and inputs without staging or writing analysis products.

---

# Running the pipeline

A recommended execution order is documented in:

`RUN_ORDER.md`

The repository also contains:

- `90_testing/` for exploratory or developmental workflows
- `99_deprecated/` for archived legacy scripts retained for reproducibility

---

# System requirements

Recommended:

- R >= 4.2
- RStudio or VS Code
- macOS/Linux preferred for large workflows

Core package ecosystem includes:

- tidyverse
- clusterProfiler
- enrichplot
- WGCNA
- limma
- EWCE
- ggplot2
- openxlsx
- igraph
- ggraph
- patchwork
- pheatmap
- ComplexHeatmap
- lme4/lmerTest
- mgcv

Additional dependencies vary between workflows.

---

# Outputs

Typical outputs include:

- publication-ready SVG/PDF figures
- enrichment tables
- network edge tables
- bootstrap stability summaries
- module score matrices
- EWCE outputs
- source-data tables
- QC reports
- PRIDE package metadata, manifests, checksums, and validation reports

---

# Reproducibility philosophy

The repository intentionally preserves:

- exploratory workflows
- older analysis versions
- testing scripts
- alternative implementations

This is done to maintain reproducibility of intermediate biological findings and figure-generation pipelines.

Older scripts are preferentially moved to `99_deprecated/` rather than deleted.

---

# Author

Tobias Pohl
