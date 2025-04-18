# pRoteomics

pRoteomics is a repository of R scripts designed for advanced proteomics analyses. The two primary scripts included focus on processing gene data for clusterProfiler analysis and performing comprehensive Gene Set Enrichment Analysis (GSEA) along with pathway visualization.

## Overview

This repository contains two main workflows:

- **UniProt ID Mapping for ClusterProfiler Analysis**  
  This script processes gene data containing UniProtKB IDs (e.g., "GENE_MOUSE"), maps them to corresponding UniProt Accessions using an official UniProt mapping file, and saves the mapped results for downstream clusterProfiler analyses.

- **Gene Set Enrichment Analysis (GSEA) Workflow using clusterProfiler**  
  This script executes GSEA and KEGG pathway enrichment on mouse datasets. It includes gene list sorting based on log₂ fold change, various functional and publication-quality plot generations, and KEGG pathway visualizations using the Pathview package.

## System Requirements

- **R** (version ≥ 4.0 recommended)
- Required R packages:
  - **Data Manipulation and I/O:** `dplyr`, `stringr`, `tidyr`, `purrr`, `readr`, `pacman`
  - **Enrichment and Plotting:** `clusterProfiler`, `pathview`, `enrichplot`, `DOSE`, `ggplot2`, `cowplot`, `ggridges`, `ggpubr`, `ggrepel`, `ggsci`, `ggthemes`, `ggExtra`, `ggforce`, `ggalluvial`, `lattice`, `latticeExtra`
  - **Bioconductor and OrgDB:** `BiocManager`, `org.Mm.eg.db`
  - **Additional Visualization:** `ggplotify`, `svglite`

- **Input Files:**
  - A gene data CSV file (named based on cell types, e.g., `TESTmcherryG1_mcherryG2.csv` or `mcherryG1_mcherryG2.csv`) located in the `Datasets` folder.
  - The UniProt mapping file (`MOUSE_10090_idmapping.dat`) downloaded and placed in the `Datasets` folder.

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/pRoteomics.git
   cd pRoteomics
   ```

2. **Prepare the Environment:**

   - Open the project in RStudio or run the scripts using the command line via Rscript.
   - Ensure the Datasets folder contains the required gene data CSV file and the UniProt mapping file.

3. **Dependency Management:**

   The scripts include code to automatically check for and install any missing packages. Ensure you have an active internet connection for package installation.

## Usage

### UniProt ID Mapping Script

This script maps gene symbols (e.g., "GENE_MOUSE") from your dataset to their corresponding UniProt Accessions.

- **Workflow Steps:**
  - Setup environment and load necessary libraries.
  - Define project directories and input file paths.
  - Validate the existence of input gene data and UniProt mapping file.
  - Preprocess gene identifiers to extract and clean primary IDs.
  - Parse the UniProt mapping file to extract a mapping from UniProtKB-ID to Accession.
  - Merge the gene data with the mapping information.
  - Save the final mapped output as a CSV file for further downstream analysis.

- **Run the Script:**

   ```bash
   Rscript path/to/UniProt_ID_Mapping.R
   ```

### Gene Set Enrichment Analysis (GSEA) Workflow Script

This script performs multiple enrichment analyses and creates high-quality visualizations for the analysis of mouse gene datasets.

- **Workflow Steps:**
  - Install and load all required packages (both CRAN and Bioconductor).
  - Define directories and load the gene dataset.
  - Prepare and sort gene lists based on log₂ fold change values.
  - Execute Gene Ontology (GO) GSEA across all ontologies.
  - Generate various plots such as dotplots, enrichment maps, network plots, ridgeplots, and GSEA curves.
  - Conduct ORA on top regulated genes and generate corresponding visualizations.
  - Convert gene IDs from UniProt to ENTREZID for KEGG analysis.
  - Perform KEGG GSEA and generate KEGG pathway visualizations using the Pathview package.
  - Optionally, execute additional enrichment analysis with EnrichGO and produce heatmaps.

- **Run the Script:**

   ```bash
   Rscript path/to/GSEA_Workflow.R
   ```

## Outputs

- **CSV Files:**
  
  Mapped gene data and any saved enrichment results are written to a Results directory.

- **Visualizations:**
  
  Publication-quality plots such as dotplots, network plots, ridgeplots, and KEGG pathway images are saved as SVG files within the Results folder.

## Contributing

Contributions, bug reports, and feature requests are welcome. Please open an issue or submit a pull request to help improve the repository.

## Author

Tobias Pohl