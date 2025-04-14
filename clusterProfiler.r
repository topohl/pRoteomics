#' Gene Set Enrichment Analysis (GSEA) Workflow using clusterProfiler
#'
#' This script performs Gene Set Enrichment Analysis (GSEA), KEGG pathway enrichment,
#' and visualization of functional gene sets using the `clusterProfiler` ecosystem.
#' It supports mouse datasets and creates various publication-quality plots.
#'
#' @details
#' The analysis includes:
#' - Setup of required packages and environment
#' - Loading and sorting of gene lists based on log2 fold change
#' - Gene Ontology (GO) GSEA
#' - KEGG GSEA (symbol to ENTREZID conversion)
#' - Functional plots: dotplots, network plots, ridgeplots, GSEA plots
#' - PubMed trend analysis for top terms
#' - GO enrichment analysis and heatmaps
#'
#' @references
#' clusterProfiler vignette: \url{https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html}
#'
#' @author
#' Tobias Pohl

# ----------------------------------------------------
# Install and load packages
# ----------------------------------------------------

# Check and install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define required packages
required_packages <- c(
  "clusterProfiler", "pathview", "enrichplot", "DOSE", "ggplot2", "ggnewscale",
  "cowplot", "ggridges", "europepmc", "ggpubr", "ggrepel", "ggsci", "ggthemes",
  "ggExtra", "ggforce", "ggalluvial", "lattice", "latticeExtra", "BiocManager",
  "org.Mm.eg.db", "ggplotify", "svglite"
)

# Install missing Bioconductor packages
bioc_packages <- c(
  "clusterProfiler", "pathview", "enrichplot", "DOSE", "org.Mm.eg.db"
)

missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc)
}

install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
}
invisible(lapply(required_packages, install_and_load))

# ----------------------------------------------------
# Define working directory, cell types, and data paths
# ----------------------------------------------------

# Set directories
working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp3_Nlgn3_development/LaserDissProteomics/GSEA"
results_dir <- file.path(working_dir, "Results")
setwd(working_dir)

# Define cell types
cell_types <- c("neuropil", "microglia")
file_name <- paste0(paste(cell_types, collapse = "_"), ".csv")
data_path <- file.path(working_dir, "Datasets", file_name)

# ----------------------------------------------------
# Load and prepare data
# ----------------------------------------------------

# Load and prepare gene data
df <- read.csv(data_path, header = TRUE)
colnames(df)[1] <- "gene_symbol"

original_gene_list <- df$log2fc
names(original_gene_list) <- df$gene_symbol
gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)

# ----------------------------------------------------
# Perform Gene Set Enrichment Analysis (GSEA)
# ----------------------------------------------------

# Set annotation package and load it
organism <- "org.Mm.eg.db" 
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# Gene Set Enrichment Analysis (GSEA)
gse <- gseGO(
  geneList = gene_list, ont = "ALL", keyType = "SYMBOL",
  minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, verbose = TRUE,
  OrgDb = organism, pAdjustMethod = "none"
)

# Helper to save plots
save_plot <- function(plot, filename) {
  ggsave(file.path(results_dir, filename), plot, units = "cm", dpi = 300)
}

# Dotplot
require(DOSE)
dot_title <- paste("GSEA of", paste(cell_types, collapse = " over "))
p1 <- clusterProfiler::dotplot(gse, showCategory = 10, split = ".sign") +
  facet_grid(. ~ .sign) +
  labs(title = dot_title)
save_plot(p1, paste("GSEAdotplot_", paste(cell_types, collapse = "_"), ".svg"))

# ----------------------------------------------------
# Enrichment Map, Network Plot, and Ridgeplot of GSEA
# ----------------------------------------------------

# Enrichment Map and Network Plot
emapplot(pairwise_termsim(gse), showCategory = 10)
save_plot(emapplot(pairwise_termsim(gse), showCategory = 10),
          paste("GSEAemap_", paste(cell_types, collapse = "_"), ".svg"))

cnetplot(gse, categorySize = "pvalue", foldChange = gene_list)
save_plot(cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
          paste("GSEAcnet_", paste(cell_types, collapse = "_"), ".svg"))

# Ridgeplot: Visualizing enrichment distribution
ridgeplot_gse <- ridgeplot(gse) +
  labs(x = "Enrichment Distribution", title = "GSEA Ridgeplot") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
save_plot(ridgeplot_gse, paste("GSEA_Ridgeplot_", paste(cell_types, collapse = "_"), ".svg"))

# GSEA Plot: Highlighting specific gene sets
gseaplot_gse <- gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1)
save_plot(gseaplot_gse, paste("GSEA_Plot_", paste(cell_types, collapse = "_"), ".svg"))

# PubMed Trend Plot: Analyzing publication trends for top enriched terms
top_terms <- head(gse@result$Description, 3)
pmcplot_gse <- pmcplot(top_terms, 2010:2018, proportion = FALSE) +
  labs(title = "Publication Trends for Top Enriched Terms")
save_plot(pmcplot_gse, paste("GSEA_PubMed_Trends_", paste(cell_types, collapse = "_"), ".svg"))

# ----------------------------------------------------
# KEGG Gene Set Enrichment Analysis
# ----------------------------------------------------

# KEGG GSEA Analysis
# Map SYMBOL to ENTREZID
ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

# Remove duplicated SYMBOLs
dedup_ids <- ids[!duplicated(ids$SYMBOL), ]

# Merge df with ENTREZID mapping
df2 <- merge(df, dedup_ids, by.x = "gene_symbol", by.y = "SYMBOL")

# Now use ENTREZID as names
kegg_gene_list <- df2$log2fc
names(kegg_gene_list) <- df2$ENTREZID

# Remove NAs and sort
kegg_gene_list <- sort(na.omit(kegg_gene_list), decreasing = TRUE)

kegg_organism <- "mmu"
kk2 <- gseKEGG(
  geneList = kegg_gene_list, organism = "mmu",
  minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "none",
  keyType = "ncbi-geneid"
)

# KEGG Dotplot and Network Plot
kegg_dot_title <- paste("KEGG GSEA Enriched Pathways of", paste(cell_types, collapse = " over "))
p2 <- clusterProfiler::dotplot(kk2, showCategory = 10, title = kegg_dot_title, split = ".sign") +
  facet_grid(. ~ .sign)
save_plot(p2, paste("KEGGpathway_", paste(cell_types, collapse = "_"), ".svg"))

emapplot(pairwise_termsim(kk2), showCategory = 10)
cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list)
ridgeplot(kk2) + labs(x = "Enrichment distribution")
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

# KEGG Pathway Visualizations using Pathview
if (!requireNamespace("pathview", quietly = TRUE)) install.packages("pathview")
library(pathview)

pathview(gene.data = kegg_gene_list, pathway.id = "mmu00030", species = kegg_organism)
pathview(gene.data = kegg_gene_list, pathway.id = "mmu00030", species = kegg_organism, kegg.native = FALSE)

# ----------------------------------------------------
# EnrichGO Analysis and Heatmap
# ----------------------------------------------------

# EnrichGO Analysis and Heatmap
go_enrich <- enrichGO(
  gene = names(gene_list), universe = names(gene_list), OrgDb = organism, 
  keyType = 'SYMBOL', readable = TRUE, ont = "ALL", pvalueCutoff = 1, 
  qvalueCutoff = 1, pAdjustMethod = "none", minGSSize = 3, maxGSSize = 800
)

# Heatmap Plot
p3 <- heatplot(go_enrich, foldChange = gene_list, showCategory = 5)
save_plot(p4, paste("GOheatmap_", paste(cell_types, collapse = "_"), ".svg"))
cowplot::plot_grid(p1, p3, ncol = 1, labels = LETTERS[1:2])
