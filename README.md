# Gene Set Enrichment Analysis (GSEA) with clusterProfiler

This repository provides an end-to-end example for performing Gene Set Enrichment Analysis (GSEA) using the R package **clusterProfiler** along with several complementary packages. The workflow demonstrates how to install and load the necessary packages, prepare gene expression data, and run both Gene Ontology (GO) and KEGG pathway enrichment analyses. Multiple visualization outputs such as dotplots, enrichment maps, network plots, ridgeplots, GSEA plots, PubMed trend plots, and heatmaps are generated throughout the analysis.

Detailed documentation for clusterProfiler can be found here:  
[clusterProfiler Documentation](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)

## Table of Contents

- [Overview](#overview)
- [Installation and Requirements](#installation-and-requirements)
- [Usage](#usage)
- [File Structure](#file-structure)
- [Citation](#citation)
- [License](#license)

## Overview

- **Purpose:**  
  This project provides a detailed workflow for performing GSEA on gene expression data using **clusterProfiler**. It covers both GO and KEGG analyses along with comprehensive visualizations.

- **Key Features:**  
  - **Automatic Package Installation:** Automatically installs required Bioconductor and CRAN packages.
  - **Data Preparation:** Loads gene expression data from a CSV file, renames columns, and generates a sorted gene list.
  - **GO Enrichment Analysis:** Uses the `gseGO` function to perform GO enrichment and generates various plots including dotplots, enrichment maps, network plots, and ridgeplots.
  - **KEGG Pathway Enrichment Analysis:** Converts gene symbols to ENTREZ IDs, performs KEGG GSEA using `gseKEGG`, and integrates with Pathview for pathway diagrams.
  - **Additional Visualizations:** Includes a PubMed trend plot and a GO enrichment heatmap.
  - **Result Saving:** Uses custom functions to save output plots in the specified results directory.

## Installation and Requirements

- **Prerequisites:**
  - R (and optionally RStudio) installed on your system.
  - Internet connectivity to download packages and annotation/pathway data.

- **Required Packages:**  
  The analysis utilizes the following R packages:
  - **Bioconductor Packages:** `clusterProfiler`, `pathview`, `enrichplot`, `DOSE`, `org.Mm.eg.db`
  - **CRAN Packages:** `ggplot2`, `ggnewscale`, `cowplot`, `ggridges`, `europepmc`, `ggpubr`, `ggrepel`, `ggsci`, `ggthemes`, `ggExtra`, `ggforce`, `ggalluvial`, `lattice`, `latticeExtra`

- **Installation Instructions:**  
  The script is designed to automatically check for and install any missing packages on runtime. Simply run the R script after adjusting the working directory and file paths to match your environment.

## Usage

1. **Clone the Repository:**  
   Clone or download this repository to your local machine.

2. **Update the Working Directory:**  
   Open the R script (e.g., `script.R`) and modify the variable `working_dir` so that it points to your local working directory. Ensure your CSV file is placed in the `/Datasets/` subdirectory.

3. **Run the Analysis:**  
   The script performs the following steps:
   - **Package Setup:**  
     Checks for and installs all required packages using a helper function (`install_and_load`).
   - **Data Preparation:**  
     Loads a CSV file with gene expression data, renames the first column to `gene_symbol`, and creates a sorted gene list based on `log2fc` values.
   - **GO Enrichment:**  
     Executes GSEA using `gseGO` with predefined parameters and generates several plots (dotplot, enrichment map, network plot, ridgeplot, and a GSEA plot).
   - **KEGG Pathway Enrichment:**  
     Converts gene symbols to ENTREZ IDs, performs KEGG GSEA using `gseKEGG`, and creates similar visualizations. Pathway diagrams are generated using Pathview.
   - **Additional Visualizations:**  
     Generates a PubMed trend plot for selected enriched terms and a GO enrichment heatmap.
   - **Saving Outputs:**  
     Uses custom functions (e.g., `save_plot`) to store generated plots in the `/Results/` directory.

4. **Adjust and Extend:**  
   Feel free to modify analysis parameters (e.g., p-value cutoffs, category limits) or integrate additional plots as needed.

## File Structure

├── README.md&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;// This file: provides an overview and instructions. <br/>
├── clusterProfiler.R&nbsp;&nbsp;&nbsp;&nbsp;// The R script containing the GSEA workflow.<br/>
├── Datasets/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;// Directory with sample CSV files containing gene expression data. <br/>
└── Results/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;// Directory where output plots and results will be saved.<br/>

## Citation

If you use this workflow or any part of this analysis in your research, please cite **clusterProfiler** as follows:

**LG Wang, Y Han, QY He. _clusterProfiler: an R package for comparing biological themes among gene clusters._ OMICS: A Journal of Integrative Biology, 2012, 16(5):284-287.**  
DOI: [10.1089/omi.2011.0118](http://dx.doi.org/10.1089/omi.2011.0118)
