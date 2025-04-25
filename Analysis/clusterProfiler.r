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

# Define cell types
cell_types <- c("cfos3", "cfos4")

# Set directories
# working_dir <- "/Users/tobiaspohl/Documents/clusterProfiler"
working_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
results_dir <- file.path(working_dir, "Results", paste(cell_types, collapse = "_"))
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}
setwd(working_dir)

file_name <- paste0(paste(cell_types, collapse = "_"), ".csv")
data_path <- file.path(working_dir, "Datasets", "mapped", file_name)

# set organism
organism <- "org.Mm.eg.db"

# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------

# Helper to save plots
save_plot <- function(plot, filename) {
  ggsave(file.path(results_dir, filename), plot, units = "cm", dpi = 300)
}

plot_dot <- function(dataset, cell_types, results_dir = results_dir) {
  dot_title <- paste("GSEA of", paste(cell_types, collapse = " over "))
  p <- clusterProfiler::dotplot(dataset, showCategory = 10, split = ".sign") +
    facet_wrap(~ .sign, nrow = 1) +
    labs(title = dot_title, x = "Gene Ratio", y = "Gene Set") +
    scale_fill_viridis_c(option = "cividis") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 12, face = "plain"),
      panel.grid.major = element_line(color = "grey85", size = 0.3),
      panel.grid.minor = element_blank()
    )
  analysis_type <- deparse(substitute(dataset))
  plot_filename <- paste0(analysis_type, "_", paste(cell_types, collapse = "_"), ".svg")
  save_plot(p, plot_filename)
  return(p)
}

# ----------------------------------------------------
# Load and prepare data
# ----------------------------------------------------

# Load and prepare gene data
df <- read.csv(data_path, header = TRUE)
colnames(df)[1] <- "gene_symbol"
original_gene_list <- df$log2fc
names(original_gene_list) <- df$gene_symbol
gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)
# remove duplicates
gene_list <- gene_list[!duplicated(names(gene_list))]

top_df <- read.csv(data_path, header = TRUE)
colnames(top_df)[1] <- "gene_symbol"
top_gene_list <- top_df$log2fc
names(top_gene_list) <- top_df$gene_symbol
top_gene_list <- sort(na.omit(top_gene_list), decreasing = TRUE)
# Select top N genes by absolute fold change
top_genes <- names(top_gene_list)[abs(top_gene_list) > 1]  # Adjust threshold as needed
top_genes <- sort(top_gene_list[top_genes], decreasing = TRUE)
top_genes <- names(top_genes)

# ----------------------------------------------------
# Perform Gene Set Enrichment Analysis (GSEA)
# ----------------------------------------------------

# Gene Set Enrichment Analysis (GSEA)
gse <- gseGO(
  geneList = gene_list, ont = "CC", keyType = "UNIPROT",
  minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, verbose = TRUE,
  OrgDb = organism, pAdjustMethod = "BH"
)

plot_dot(gse, cell_types)

# save core_enrichment results
core_enrichment <- gse@result$core_enrichment
# Save the GSEA results
core_dir <- file.path(results_dir, "core_enrichment")
if (!dir.exists(core_dir)) {
  dir.create(core_dir, recursive = TRUE)
}
write.csv(gse@result, file = file.path(core_dir, paste("coreEnrichment_", paste(cell_types, collapse = "_"), ".csv")))

# ----------------------------------------------------
# Focus on the top regulated genes
# ----------------------------------------------------

# Set annotation package and load it
#organism <- "org.Mm.eg.db" 
#BiocManager::install(organism, character.only = TRUE)
#library(organism, character.only = TRUE)

# Perform enrichment analysis (ORA)
#library(org.Mm.eg.db)
ora <- enrichGO(gene = top_genes, ont = "CC", keyType = "UNIPROT",
             minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1,
             OrgDb = organism, pAdjustMethod = "none")

# Dotplot for ORA
p7 <- clusterProfiler::dotplot(ora, showCategory = 10) +
  labs(title = "ORA of Top Regulated Genes") +
  scale_x_continuous(limits = c(0, 1))  # Scale Gene Ratio from 0 to 1
save_plot(p7, paste("ORA_dotplot_", paste(cell_types, collapse = "_"), ".svg"))
# Save the ORA results
#write.csv(ora@result, file = file.path(results_dir, paste("ORA_results_", paste(cell_types, collapse = "_"), ".csv")))

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
pmcplot_gse <- pmcplot(top_terms, 2010:2025, proportion = FALSE) +
  labs(title = "Publication Trends for Top Enriched Terms")
save_plot(pmcplot_gse, paste("GSEA_PubMed_Trends_", paste(cell_types, collapse = "_"), ".svg"))

# ----------------------------------------------------
# KEGG Gene Set Enrichment Analysis
# ----------------------------------------------------

# KEGG GSEA Analysis
# Map SYMBOL to ENTREZID
ids <- bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

# Remove duplicated UNIPROT IDs
dedup_ids <- ids[!duplicated(ids$UNIPROT), ]

# Merge df with ENTREZID mapping
df2 <- merge(df, dedup_ids, by.x = "gene_symbol", by.y = "UNIPROT")

# Now use ENTREZID as names
kegg_gene_list <- df2$log2fc
names(kegg_gene_list) <- df2$ENTREZID
# Remove duplicates
kegg_gene_list <- kegg_gene_list[!duplicated(names(kegg_gene_list))]

# Remove NAs and sort
kegg_gene_list <- sort(na.omit(kegg_gene_list), decreasing = TRUE)

kegg_organism <- "mmu"
kk2 <- gseKEGG(
  geneList = kegg_gene_list, organism = "mmu",
  minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
  keyType = "ncbi-geneid"
)

# KEGG Dotplot and Network Plot
#kegg_dot_title <- paste("KEGG GSEA Enriched Pathways of", paste(cell_types, collapse = " over "))
#p2 <- clusterProfiler::dotplot(kk2, showCategory = 10, title = kegg_dot_title, split = ".sign") +
#  facet_grid(. ~ .sign)
#save_plot(p2, paste("KEGGpathway_", paste(cell_types, collapse = "_"), ".svg"))

plot_dot(kk2, cell_types)

emapplot(pairwise_termsim(kk2), showCategory = 10)
cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list)
ridgeplot(kk2) + labs(x = "Enrichment distribution")
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

# KEGG Pathway Visualizations using Pathview
if (!requireNamespace("pathview", quietly = TRUE)) install.packages("pathview")
library(pathview)

# KEGG Pathway Visualizations using Pathview

# Define additional pathway IDs for visualization.
path_ids <- c(
  "mmu04110", "mmu04115", "mmu04114", "mmu04113", "mmu04112", "mmu04111",
  "mmu04116", "mmu04117", "mmu04118", "mmu04119", "mmu04720", "mmu04721",
  "mmu04722", "mmu04725", "mmu04726", "mmu04727", "mmu04724", "mmu04080",
  "mmu00030"
)

# Create a subdirectory for pathview results and get its absolute path.
pathview_dir <- normalizePath(file.path(results_dir, "pathview"), winslash = "/", mustWork = FALSE)
if (!dir.exists(pathview_dir)) {
  dir.create(pathview_dir, recursive = TRUE)
}

oldwd <- getwd()
setwd(pathview_dir)
# Generate pathview plots for each pathway in the list using absolute paths.
lapply(path_ids, function(pid) {
  pathview(
    gene.data = kegg_gene_list, 
    pathway.id = pid, 
    species = kegg_organism, 
  )
})
setwd(oldwd)
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
#p3 <- heatplot(go_enrich, foldChange = gene_list, showCategory = 5)
#save_plot(p4, paste("GOheatmap_", paste(cell_types, collapse = "_"), ".svg"))
#cowplot::plot_grid(p1, p3, ncol = 1, labels = LETTERS[1:2])
