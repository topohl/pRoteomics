#' @title Comparative GO Enrichment Analysis Using clusterProfiler
#'
#' @description This script performs a comprehensive analysis to compare the enrichment of Gene Ontology (GO) terms across multiple gene clusters. The workflow includes:
#' \itemize{
#'   \item Installing and loading necessary packages if they are not already installed.
#'   \item Reading gene lists from CSV files (each representing a cell type or condition) in a specified directory.
#'   \item Converting gene symbols to Entrez IDs using the bitr function from clusterProfiler for accurate enrichment analysis.
#'   \item Executing a comparative GO enrichment (Biological Process) analysis across the gene clusters using the compareCluster function.
#'   \item Visualizing the enrichment results with a dotplot enhanced by ggplot2 customizations.
#' }
#'
#' @details The analysis utilizes the clusterProfiler package for enrichment testing and the enrichplot package for visualizations. Additionally, it leverages organism-specific annotations from the org.Mm.eg.db package (Mus musculus) to ensure accurate mapping of gene identifiers.
#'
#' @note Ensure that the file paths (e.g., for gene list CSV files and project directories) are correctly specified according to your local environment. The script may require adjustments if your gene identifiers are in different formats.
#'
#' @author
#' Tobias Pohl
#'

# Install required packages if not already installed
required_pkgs <- c("clusterProfiler", "org.Mm.eg.db", "enrichplot")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

if (!require("pacman")) install.packages("pacman")
pacman::p_load(clusterProfiler, org.Mm.eg.db, enrichplot, ggplot2, dplyr, tidyr, ggrepel)

# Example: Read in gene lists from files (one list per cell type/condition)
# Each file should contain a single column of gene symbols or Entrez IDs
file_list <- list.files(path = "/Users/tobiaspohl/Documents/clusterProfiler/Results/core_enrichment", pattern = "*.csv", full.names = TRUE)
gene_lists <- lapply(file_list, function(f) scan(f, what = "", quiet = TRUE))
names(gene_lists) <- tools::file_path_sans_ext(basename(file_list))

# If your genes are in SYMBOL, convert to ENTREZID (recommended for clusterProfiler)
gene_lists_entrez <- lapply(gene_lists, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
})

# Comparative GO enrichment (Biological Process)
ck <- compareCluster(
  geneCluster = gene_lists_entrez,
  fun = "enrichGO",
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# View results
head(as.data.frame(ck))

# Visualization: Dotplot for comparison
dotplot(ck, showCategory = 10) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))