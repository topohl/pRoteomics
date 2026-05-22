<img width="759" alt="Screenshot 2025-04-18 at 14 47 44" src="https://github.com/user-attachments/assets/be03dd69-eea1-4f5f-bdef-a05ba77392fe" />

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
90_testing/
99_deprecated/
```

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

---

## 03_qc_exploration

Purpose:
- PCA
- variance partitioning
- rank abundance analysis
- protein/peptide QC

Representative scripts:
- `03_pcaPlot.r`
- `05_pcaPlot_v3.r`
- `07_varPart.r`

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
- `03_compare_pathways.r`

---

## 05_celltype_enrichment_EWCE

Purpose:
- EWCE analyses
- spatial cell-type interpretation
- measured-proteome-aware enrichment workflows

Representative scripts:
- `01_EWCE_E9.r`

---

## 06_modules_WGCNA

Purpose:
- WGCNA
- module scoring
- module preservation
- overlap-based module generation

Representative scripts:
- `01_WGCNA.r`
- `02_WGCNAtraitpreservation.r`
- `03_module_spatial_networks.r`

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

---

## 08_behavior_physio_coupling

Purpose:
- connect proteomics and network structure with behavior and physiology
- movement/stress score integration
- systems-level phenotype coupling

Representative scripts:
- `01_correlate_proteomics_with_behavior.r`
- `02_network_behavior_coupling.r`

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

