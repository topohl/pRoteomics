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

# Modularized package installation and loading functions

# Function to ensure BiocManager is installed
checkBiocManager <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
}

# Function to install missing Bioconductor packages
installBioC <- function(bioc_packages) {
  missing <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    BiocManager::install(missing)
  }
}

# Function to install and load a CRAN package
installAndLoadCRAN <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

# Main function to setup all required packages
setupPackages <- function() {
  checkBiocManager()

  required_packages <- c(
    "clusterProfiler", "pathview", "enrichplot", "DOSE", "ggplot2", "ggnewscale",
    "cowplot", "ggridges", "europepmc", "ggpubr", "ggrepel", "ggsci", "ggthemes",
    "ggExtra", "ggforce", "ggalluvial", "lattice", "latticeExtra", "BiocManager",
    "org.Mm.eg.db", "ggplotify", "svglite"
  )

  bioc_packages <- c("clusterProfiler", "pathview", "enrichplot", "DOSE", "org.Mm.eg.db")

  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}

# Call the main setup function to install and load all necessary packages
setupPackages()

# ----------------------------------------------------
# Define working directory, cell types, and data paths
# ----------------------------------------------------

# Define cell types
cell_types <- c("mcherry2", "mcherry4")

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

# Set comparison ontology
ont <- "BP"  # Cellular Component
# ont <- "BP"  # Biological Process
# ont <- "MF"  # Molecular Function
# ont <- "ALL"  # All Ontologies

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

# Save the GSEA results in the primary directory using ont as the subfolder
core_dir <- file.path(results_dir, "core_enrichment", ont)
if (!dir.exists(core_dir)) {
  dir.create(core_dir, recursive = TRUE)
}
write.csv(gse@result,
          file = file.path(core_dir, paste0("coreEnrichment_", ont, "_", paste(cell_types, collapse = "_"), ".csv")))

# Additionally, save the GSEA results in the specified Datasets directory using ont as the subfolder
additional_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment", ont)
if (!dir.exists(additional_dir)) {
  dir.create(additional_dir, recursive = TRUE)
}
write.csv(gse@result,
          file = file.path(additional_dir, paste0("coreEnrichment_", ont, "_", paste(cell_types, collapse = "_"), ".csv")))

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
  "mmu00030",
  "mmu04151",  # P!3K-Akt signaling pathway
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
    low = "#6698CC",
    mid = "white",
    high = "#F08C21",
    file.type = "svg"
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


# ----------------------------------------------------
# KEGG GSEA Analysis with Predefined UniProt IDs
# ----------------------------------------------------

# Define your UniProt gene list
#selected_uniprot <- c("Q8K3R3", "Q8R071", "P11798", "P28652", "P20444", "P68404", "P05132", "P68181")
#selected_uniprot <- c(
#  "Q9JL06", "Q9EQ31", "Q149R9", "P61793", "Q6NS65", "Q9QY96", "Q8JZL2", "P56942",
#  "Q8BZ39", "Q5H8A1", "P12023", "P27784", "P51670", "P10107", "O88536", "P11859",
#  "Q9QXK8", "O55040", "Q6S9I3", "P32299", "Q61125", "P62880", "P29387", "Q61011",
#  "P62881", "P62874", "Q80SZ7", "Q9DAS9", "Q61016", "P63078", "Q61012", "Q61017",
#  "P63216", "Q9CXP8", "Q9JMF3", "P50153", "P61953", "P63213", "O08899", "O08850",
#  "Q9CX84", "Q9JL25", "Q8K443", "P97428", "V9GXQ3", "Q9QZB0", "Q99PG4", "Q9DC04",
#  "P21278", "P21279", "P30678", "P30677", "P63085", "P28867", "P20444", "Q99MK8",
#  "A3KGF7", "Q9Z1B3", "Q91UZ1", "P51432", "P35991", "Q9CWR0", "Q62245", "Q0KL02",
#  "A2CG49", "P97717", "P97714", "P97718", "P30987", "Q5U431", "P55095", "Q61606",
#  "Q8BZP8", "P0C0P8", "P35455", "O08908", "Q64143", "P26450", "P42337", "O08849",
#  "A0A5F8MPV0", "Q9Z329", "P70227", "P11881", "Q99LR1", "Q8R2Y0", "Q6WQJ1", "Q91WC9",
#  "O35678", "P16054", "P23298", "Q02111", "Q61143", "Q9QZC1", "Q9WVC5", "O88673",
#  "D3YXJ0", "Q6P5E8", "Q6NS52", "Q91WG7", "Q80UP3", "E9PUQ8", "D3YWQ0", "Q9R1C6",
#  "F6UUZ3", "Q8VCK6", "Q62361", "P21761", "P41539", "P30549", "Q920H4", "Q9ERZ3",
#  "P12657", "P55099", "P47937", "P56481", "O08786", "P09240", "Q3UVX5", "P97772",
#  "P70174", "Q8K4Z6", "P70259", "O54798", "O54799", "P21729", "Q8R1I2", "Q9CR53",
#  "Q924H0", "E9Q468", "Q9WVA8", "Q9QXU7", "Q14A28", "Q9JKL1", "Q8K458", "Q9QXZ9",
#  "P47993", "Q9R0M1", "P43117", "Q8BFU7", "Q99P50", "Q9EQX0", "Q62035", "O88855",
#  "Q9JJL9", "Q8CE23", "G3UWA8", "Q920A1", "P30548", "Q99JA4", "P49650", "O55241",
#  "P58308", "Q8BMC0", "P97926", "P35454", "P35383", "P19221", "O08675", "O88634",
#  "P30558", "P55086", "P29754", "Q9ERK9", "Q91V45", "Q6Y4S4", "Q7TMA4", "Q8BFQ3",
#  "Q9Z282", "Q8BUD0", "Q61038", "P35375", "Q8BLG2", "P48757", "P58307", "Q9WU02",
#  "Q62463", "Q3UFD7", "P13562", "Q01776", "Q76JU9", "Q765I1", "Q9QZQ3", "Q8VIH9",
#  "P48302", "Q61614", "P22387", "P22389", "P48299", "P70310", "O88319", "Q9D3P9",
#  "P35363", "Q02152", "P34968", "Q9WUT3", "P18653", "P18654", "Q01147-1", "Q60631",
#  "Q06186", "Q01279", "Q61411", "P32883", "P08556", "P28862", "Q9WVS8", "Q63844"
#)

# Define the selected UniProt IDs for the custom GSEA analysis
selected_uniprot <- c(
  "P21279", "P21278", "Q9Z1B3", "P51432", "P11881", "P68404",
  "P63318", "P0DP26", "P0DP27", "P11798", "P28652", "Q61411",
  "Q99N57", "P31938", "P63085", "Q63844", "Q8BWG8", "Q91YI4", "V9GXQ9"
)

# Convert to ENTREZIDs
selected_entrez <- bitr(
  selected_uniprot,
  fromType = "UNIPROT",
  toType = "ENTREZID",
  OrgDb = "org.Mm.eg.db"
)

# Extract fold changes for these genes from your dataframe
# Assuming 'df' contains 'gene_symbol' = UNIPROT and 'log2fc'
selected_df <- merge(df, selected_entrez, by.x = "gene_symbol", by.y = "UNIPROT")
selected_kegg_list <- selected_df$log2fc
names(selected_kegg_list) <- selected_df$ENTREZID

# Remove duplicates and sort
selected_kegg_list <- selected_kegg_list[!duplicated(names(selected_kegg_list))]
selected_kegg_list <- sort(na.omit(selected_kegg_list), decreasing = TRUE)

# Run GSEA on selected gene set
gsea_kegg_selected <- gseKEGG(
  geneList = selected_kegg_list,
  organism = "mmu",
  minGSSize = 3,
  maxGSSize = 800,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  keyType = "ncbi-geneid",
  verbose = TRUE
)

# Plot and save
plot_dot(gsea_kegg_selected, cell_types)
write.csv(
  gsea_kegg_selected@result,
  file = file.path(results_dir, paste0("KEGG_GSEA_Predefined_UniProt_", paste(cell_types, collapse = "_"), ".csv"))
)

# Optional plots
emapplot(pairwise_termsim(gsea_kegg_selected), showCategory = 10)
cnetplot(gsea_kegg_selected, categorySize = "pvalue", foldChange = selected_kegg_list)
ridgeplot(gsea_kegg_selected)
gseaplot(gsea_kegg_selected, by = "all", title = gsea_kegg_selected$Description[1], geneSetID = 1)



# ----------------------------------------------------
# Custom GSEA: NK3R-signalling pathway (UniProt IDs)
# ----------------------------------------------------

# Define your custom gene set as TERM2GENE data.frame
#nk3r_genes <- c("Q8K3R3", "Q8R071", "P11798", "P28652", "P20444", "P68404", "P05132", "P68181")

#nk3r_genes <- c(
#  "Q9JL06", "Q9EQ31", "Q149R9", "P61793", "Q6NS65", "Q9QY96", "Q8JZL2", "P56942",
#  "Q8BZ39", "Q5H8A1", "P12023", "P27784", "P51670", "P10107", "O88536", "P11859",
#  "Q9QXK8", "O55040", "Q6S9I3", "P32299", "Q61125", "P62880", "P29387", "Q61011",
#  "P62881", "P62874", "Q80SZ7", "Q9DAS9", "Q61016", "P63078", "Q61012", "Q61017",
#  "P63216", "Q9CXP8", "Q9JMF3", "P50153", "P61953", "P63213", "O08899", "O08850",
#  "Q9CX84", "Q9JL25", "Q8K443", "P97428", "V9GXQ3", "Q9QZB0", "Q99PG4", "Q9DC04",
#  "P21278", "P21279", "P30678", "P30677", "P63085", "P28867", "P20444", "Q99MK8",
#  "A3KGF7", "Q9Z1B3", "Q91UZ1", "P51432", "P35991", "Q9CWR0", "Q62245", "Q0KL02",
#  "A2CG49", "P97717", "P97714", "P97718", "P30987", "Q5U431", "P55095", "Q61606",
#  "Q8BZP8", "P0C0P8", "P35455", "O08908", "Q64143", "P26450", "P42337", "O08849",
#  "A0A5F8MPV0", "Q9Z329", "P70227", "P11881", "Q99LR1", "Q8R2Y0", "Q6WQJ1", "Q91WC9",
#  "O35678", "P16054", "P23298", "Q02111", "Q61143", "Q9QZC1", "Q9WVC5", "O88673",
#  "D3YXJ0", "Q6P5E8", "Q6NS52", "Q91WG7", "Q80UP3", "E9PUQ8", "D3YWQ0", "Q9R1C6",
#  "F6UUZ3", "Q8VCK6", "Q62361", "P21761", "P41539", "P30549", "Q920H4", "Q9ERZ3",
#  "P12657", "P55099", "P47937", "P56481", "O08786", "P09240", "Q3UVX5", "P97772",
#  "P70174", "Q8K4Z6", "P70259", "O54798", "O54799", "P21729", "Q8R1I2", "Q9CR53",
#  "Q924H0", "E9Q468", "Q9WVA8", "Q9QXU7", "Q14A28", "Q9JKL1", "Q8K458", "Q9QXZ9",
#  "P47993", "Q9R0M1", "P43117", "Q8BFU7", "Q99P50", "Q9EQX0", "Q62035", "O88855",
#  "Q9JJL9", "Q8CE23", "G3UWA8", "Q920A1", "P30548", "Q99JA4", "P49650", "O55241",
#  "P58308", "Q8BMC0", "P97926", "P35454", "P35383", "P19221", "O08675", "O88634",
#  "P30558", "P55086", "P29754", "Q9ERK9", "Q91V45", "Q6Y4S4", "Q7TMA4", "Q8BFQ3",
#  "Q9Z282", "Q8BUD0", "Q61038", "P35375", "Q8BLG2", "P48757", "P58307", "Q9WU02",
#  "Q62463", "Q3UFD7", "P13562", "Q01776", "Q76JU9", "Q765I1", "Q9QZQ3", "Q8VIH9",
#  "P48302", "Q61614", "P22387", "P22389", "P48299", "P70310", "O88319", "Q9D3P9",
#  "P35363", "Q02152", "P34968", "Q9WUT3", "P18653", "P18654", "Q01147-1", "Q60631",
#  "Q06186", "Q01279", "Q61411", "P32883", "P08556", "P28862", "Q9WVS8", "Q63844"
#)

#nk3r_genes <- c(
#  "P21279", "P21278", "Q9Z1B3", "P51432", "P11881", "P68404",
#  "P63318", "P0DP26", "P0DP27", "P11798", "P28652", "Q61411",
#  "Q99N57", "P31938", "P63085", "Q63844", "Q8BWG8", "Q91YI4", 
#  "V9GXQ9")

# nk3r core genes
nk3r_genes <- c(
  "P21279", "P21278", "P51432", "P11881", "P63318",
  "P68404", "P0DP26", "P0DP27", "P0DP28", "P11798",
  "P28652", "P47937", "P47713"
)

term2gene_nk3r <- data.frame(
  term = rep("NK3R-signalling", length(nk3r_genes)),
  gene = nk3r_genes
)

# Prepare gene list: log2FCs named by UniProt ID
custom_gene_list <- df$log2fc
names(custom_gene_list) <- df$gene_symbol
custom_gene_list <- sort(na.omit(custom_gene_list), decreasing = TRUE)
custom_gene_list <- custom_gene_list[!duplicated(names(custom_gene_list))]

# Run GSEA with your custom term2gene set
gsea_nk3r <- clusterProfiler::GSEA(
  geneList = custom_gene_list,
  TERM2GENE = term2gene_nk3r,
  pvalueCutoff = 1,
  minGSSize = 1,
  maxGSSize = 500,
  verbose = TRUE
)

# Save and visualize
plot_dot(gsea_nk3r, cell_types)
openxlsx::write.xlsx(gsea_nk3r@result,
           file = file.path(results_dir, paste0("Custom_GSEA_NK3R_", paste(cell_types, collapse = "_"), ".xlsx")))

# Additional plot
gseaplot(gsea_nk3r, by = "all", title = "NK3R-signalling", geneSetID = 1)


# ----------------------------------------------------
# Celltype scoring
# ----------------------------------------------------

# load celltypes csv file
celltype_file <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler",
                           "Datasets/celltypes.xlsx")
# read the celltype file
library(openxlsx)
celltype_df <- read.xlsx(celltype_file, sheet = 1)

uniprot_mapping_file_path <- file.path(working_dir, "Datasets", "MOUSE_10090_idmapping.dat")

cat("Loading UniProt mapping file from:", uniprot_mapping_file_path, "\n")
uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Extract Gene_Name mappings
gene_names <- uniprot_df %>%
  filter(V2 == "Gene_Name") %>%
  distinct(V1, .keep_all = TRUE) %>%
  rename(
    UniProtID = V1,
    Gene_Name = V3
  ) %>%
  select(UniProtID, Gene_Name)

# Extract UniProtKB-ID mappings
uniprot_ids <- uniprot_df %>%
  filter(V2 == "UniProtKB-ID") %>%
  distinct(V1, .keep_all = TRUE) %>%
  rename(
    UniProtID = V1,
    UniProtKB_ID = V3
  ) %>%
  select(UniProtID, UniProtKB_ID)

# Combine both mappings into one gene_map
gene_map <- gene_names %>%
  left_join(uniprot_ids, by = "UniProtID")

# load gct file
read_gct <- function(file_path) {
  # Read the file, skipping first two lines (version and size info)
  gct_data <- read.delim(file_path, skip = 2, header = FALSE, check.names = FALSE)
  
  # Rename "Name" to "Gene" if needed
  if ("Name" %in% colnames(gct_data)) {
    colnames(gct_data)[1] <- "Gene"
  }
  
  return(gct_data)
}

gct_file <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/neha_spatial_group.gct"
gct_data <- read_gct(gct_file)

head(gct_data)

# === 2. Extract metadata ===
metadata <- gct_data[1:4, ]  # Keep all columns including the first column (IDs)
colnames(metadata) <- paste0("V", seq_len(ncol(metadata)))  # Temporary colnames

sample_id      <- as.character(unlist(gct_data[1, -1]))
celltype       <- as.character(unlist(gct_data[2, -1]))
groups          <- as.character(unlist(gct_data[3, -1]))
celltype_group <- as.character(unlist(gct_data[4, -1]))

metadata_df <- data.frame(
  ColName         = paste0("V", 2:(ncol(gct_data))),  # for joining later
  Sample          = sample_id,
  Celltype        = celltype,
  Group           = groups,
  Celltype_Group  = celltype_group,
  stringsAsFactors = FALSE
)

# === 3. Extract expression data ===
expr_data <- gct_data[-c(1:4), ]
colnames(expr_data)[1] <- "Protein_ID"

# Restore sample names as colnames
colnames(expr_data)[-1] <- sample_id

# === 4. Reshape to long format ===
expr_long <- expr_data %>%
  pivot_longer(
    cols = -Protein_ID,
    names_to = "Sample",
    values_to = "Expression"
  )

# === 5. Join with metadata ===
expr_annotated <- expr_long %>%
  left_join(metadata_df %>% select(-ColName), by = "Sample") %>%
  mutate(Expression = as.numeric(Expression)) %>%
  mutate(Protein_ID = sub(";.*", "", Protein_ID))  # Keep only the first entry before ";" 

# Join based on Protein_ID matching UniProtKB_ID
expr_annotated <- expr_annotated %>%
  left_join(gene_map, by = c("Protein_ID" = "UniProtKB_ID"))

# Add final label column (use Gene_Name if available, otherwise fall back to Protein_ID)
expr_annotated <- expr_annotated %>%
  relocate(Gene_Name, .before = 1) %>%
  relocate(UniProtID, .after = Gene_Name) %>%
  relocate(Expression, .after = last_col())

# === Done! ===
glimpse(expr_annotated)




# === 1. Reshape marker list ===
celltype_long <- celltype_df %>%
  pivot_longer(cols = everything(), names_to = "Celltype_Class", values_to = "Gene_Label") %>%
  filter(!is.na(Gene_Label) & Gene_Label != "")

# === 2. Filter expr_annotated_named for matching genes ===
marker_expr <- expr_annotated_named %>%
  filter(Gene_Label %in% celltype_long$Gene_Label)

# === 3. Add Celltype_Class info to each expression entry ===
marker_expr <- marker_expr %>%
  left_join(celltype_long, by = "Gene_Label")

# === 4. Compute average expression for each Celltype_Class per Celltype_Group ===
celltype_scores <- marker_expr %>%
  group_by(Celltype_Group, Celltype_Class) %>%
  summarise(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  )

# === 5. (Optional) Pivot wider for a table format ===
celltype_score_table <- celltype_scores %>%
  pivot_wider(names_from = Celltype_Class, values_from = Mean_Expression)

# create graph of celltype scores
library(ggplot2)
celltype_score_plot <- ggplot(celltype_scores, aes(x = Celltype_Group, y = Mean_Expression, fill = Celltype_Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Celltype Scores by Group", x = "Celltype Group", y = "Mean Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot in working directory
ggsave(file.path(working_dir, "celltype_scores_by_group.svg"), celltype_score_plot, units = "cm", dpi = 300)