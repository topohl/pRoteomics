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
    "org.Mm.eg.db", "ggplotify", "svglite", "tidyr", "dplyr", "pheatmap", "proxy",
    "tibble", "stringr"
  )

  bioc_packages <- c("clusterProfiler", "pathview", "enrichplot", "DOSE", "org.Mm.eg.db")

  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}

# Call the main setup function to install and load all necessary packages
setupPackages()

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
# Define working directory, cell types, and data paths
# ----------------------------------------------------

# define folder that contains the .csv files with the differential expression results
file_direction <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/neuron-phenotypeWithinUnit"

# get all csv files (full paths)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

# Loop: for each file set data_path and derive cell_types from the filename (without extension)
# Note: this loop must enclose the remainder of the script so the analysis runs per file.
for (data_path in file_list) {
  file_name <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # derive cell types by splitting the filename on underscores
  # e.g. "CA1_slmres_CA1_slmcon.csv" -> c("CA1","slmres","CA1","slmcon")
  # keep the full vector so downstream code can choose how to use it
  cell_types <- unlist(strsplit(name_no_ext, "_"))

  # a compact tag for naming output files / folders
  file_tag <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)

  message("Processing file: ", file_name, "  -> cell_types: ", paste(cell_types, collapse = ", "))

  # set data_path and cell_types for the rest of the script to use
  # data_path (full path) is already set; cell_types and file_tag are available to downstream code

  # --- BEGIN analysis for this file ---
  # The rest of your script should run here for the current data_path.
  # IMPORTANT: make sure there's a matching closing `}` at the very end of the script
  # to close this for-loop so subsequent iterations run.

# Define cell types
cell_types <- c("CA1_slmres", "CA1_slmcon")

# Set directories
# working_dir <- "/Users/tobiaspohl/Documents/clusterProfiler"
working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
results_dir <- file.path(working_dir, "Results", paste(cell_types, collapse = "_"))
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}
setwd(working_dir)

file_name <- paste0(paste(cell_types, collapse = "_"), ".csv")
data_path <- file.path(working_dir, "Datasets", "mapped", "neuron-phenotypeWithinUnit", file_name)

# set organism
organism <- "org.Mm.eg.db"

# Set comparison ontology
ont <- "BP"  # Biological Process
# ont <- "CC"  # Cellular Component
# ont <- "BP"  # Biological Process
# ont <- "MF"  # Molecular Function
# ont <- "ALL"  # All Ontologies

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
  geneList = gene_list, ont = ont, keyType = "UNIPROT",
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
additional_dir <- file.path(results_dir, "core_enrichment", ont)
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
  "mmu04151"  # P!3K-Akt signaling pathway
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
                           "Datasets/celltypes_long.xlsx")
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

grouping_factor <- "Celltype"

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
  filter(!grepl("Background", Celltype)) %>%
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
marker_expr <- expr_annotated %>%
  filter(Gene_Name %in% celltype_long$Gene_Label)

# === 3. Add Celltype_Class info to each expression entry ===
marker_expr <- marker_expr %>%
  left_join(celltype_long, by = c("Gene_Name" = "Gene_Label"))

# === 4. Compute average expression for each Celltype_Class per grouping_factor ===
celltype_scores <- marker_expr %>%
  group_by(!!sym(grouping_factor), Celltype_Class) %>%
  summarise(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  )

# === 5. (Optional) Pivot wider for a table format ===
celltype_score_table <- celltype_scores %>%
  pivot_wider(names_from = Celltype_Class, values_from = Mean_Expression)

# create graph of celltype scores
library(ggplot2)
celltype_score_plot <- ggplot(celltype_scores, aes(x = !!sym(grouping_factor), y = Mean_Expression, fill = Celltype_Class)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Celltype Scores by Group",
    x = grouping_factor,  # Beschriftung bleibt dynamisch
    y = "Mean Expression",
    fill = "Cell Class"
  ) +
  scale_fill_manual(values = c("#f36d07", "#455A64", "#c7c7c7")) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_line(color = "black", size = 1),
    legend.position = c(0.9, 0.9),
    legend.background = element_rect(color = "black", fill = NA),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 16)
  )

# Save the plot in working directory
ggsave(file.path(working_dir, "celltype_scores_by_group.svg"), celltype_score_plot, units = "cm", dpi = 300)

# --- 1. Gene-Klasse-Mapping vorbereiten
gene_class_df <- celltype_long %>%
  distinct(Gene_Label, Celltype_Class) %>%
  rename(Gene_Name = Gene_Label)

annotation_row <- data.frame(
  Celltype_Class = gene_class_df$Celltype_Class
)
rownames(annotation_row) <- gene_class_df$Gene_Name

# --- 2. Mittelwert pro Gen und grouping_factor berechnen
marker_matrix <- marker_expr %>%
  group_by(Gene_Name, !!sym(grouping_factor)) %>%
  summarise(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  { print(grouping_factor); print(colnames(.)); . } %>%
  pivot_wider(
    names_from = !!sym(grouping_factor),
    values_from = Mean_Expression
  ) %>%
  column_to_rownames("Gene_Name") %>%
  as.matrix()

# --- 3. Clean matrix: remove rows/columns with all NA or insufficient data
marker_matrix_clean <- marker_matrix[rowSums(is.na(marker_matrix)) < ncol(marker_matrix), ]
marker_matrix_clean <- marker_matrix_clean[, colSums(is.na(marker_matrix_clean)) < nrow(marker_matrix_clean)]
marker_matrix_clean <- marker_matrix_clean[rowSums(!is.na(marker_matrix_clean)) >= 2, ]
marker_matrix_clean <- marker_matrix_clean[, colSums(!is.na(marker_matrix_clean)) >= 2]

# --- 4. Optional: Entferne einzelne problematische Gene
marker_matrix_clean <- marker_matrix_clean[!grepl("^Oasl2\\s*$", rownames(marker_matrix_clean)), ]

# --- 5. Ersetze NaN/Inf durch NA
marker_matrix_clean[is.nan(marker_matrix_clean)] <- NA
marker_matrix_clean[is.infinite(marker_matrix_clean)] <- NA

# --- 6. annotation_row nach cleaned matrix filtern
annotation_row <- annotation_row[rownames(marker_matrix_clean), , drop = FALSE]

# --- 7. Farbschema für Annotationen
annotation_colors <- list(
  Celltype_Class = c(
    "gaba"   = "#f36d07",
    "vglut1" = "#455A64",
    "vglut2" = "#c7c7c7"
  )
)

# --- 8. Heatmap speichern
pheatmap(
  marker_matrix_clean,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
  na_col = "grey",
  main = paste("Marker Expression per", grouping_factor),
  fontsize = 12,
  fontsize_row = 10,
  fontsize_col = 12,
  cellheight = 12,
  cellwidth = 12,
  border_color = NA,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  filename = file.path(working_dir, paste0("marker_expression_heatmap_", grouping_factor, ".pdf"))
)


# Z-transform expression per gene across all samples
expr_scaled <- marker_expr %>%
  group_by(Gene_Name) %>%
  mutate(Z_Expression = scale(Expression)) %>%
  ungroup()

# Then average z-scores per grouping_factor × Celltype_Class
celltype_z_scores <- expr_scaled %>%
  group_by(grouping_factor, Celltype_Class) %>%
  summarise(
    Mean_Z = mean(Z_Expression, na.rm = TRUE),
    .groups = "drop"
  )

heatmap_plot <- ggplot(celltype_z_scores, aes(x = reorder(Celltype_Class, Mean_Z, FUN = median), y = grouping_factor, fill = Mean_Z)) +
  geom_tile(color = "grey80", size = 0.2, width = 0.5) +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0, name = "Z-Score",
    guide = guide_colorbar(barwidth = 0.8, barheight = 20)
  ) +
  labs(
    title = "Celltype Marker Signature by Group",
    x = "Marker Class",
    y = "Sample Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave(file.path(working_dir, "celltype_scores_heatmap.svg"), heatmap_plot, width = 16, height = 9, units = "cm", dpi = 300)



# 1. Sum z-scores by grouping_factor and Celltype_Class
proportional_scores <- expr_scaled %>%
  group_by(grouping_factor, Celltype_Class) %>%
  summarise(
    Mean_Z = mean(Z_Expression, na.rm = TRUE)
    .groups = "drop"
  )

# 2. Calculate proportions per grouping_factor
proportional_scores <- proportional_scores %>%
  group_by(grouping_factor) %>%
  mutate(
    Proportion = Total_Z / sum(Total_Z, na.rm = TRUE) * 100
  ) %>%
  ungroup()

  prop_barplot <- ggplot(proportional_scores, aes(x = grouping_factor, y = Proportion, fill = Celltype_Class)) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(
    title = "Proportional Celltype Signature per Group",
    x = "Celltype Group",
    y = "Proportion of Total Marker Signal (%)",
    fill = "Celltype Class"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "#f0f0f0")
  )

ggsave(file.path(working_dir, "celltype_relative_proportions.svg"), prop_barplot, width = 16, height = 10, units = "cm", dpi = 300)





# ----------------------------------------------------
# Define input directory and enumerate CSVs
# ----------------------------------------------------
file_direction <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/neuron-phenotypeWithinUnit"

file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

# ----------------------------------------------------
# Process each file
# ----------------------------------------------------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)
  cell_types  <- unlist(strsplit(name_no_ext, "_"))
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  # derive cell_types from the actual file name by splitting at the SECOND underscore
  name_no_ext <- tools::file_path_sans_ext(file_name)
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L)) {
    # no underscore -> single entry (keep as character vector)
    cell_types <- c(name_no_ext)
  } else if (length(us) == 1L) {
    # only one underscore -> split into two parts
    cell_types <- strsplit(name_no_ext, "_", fixed = TRUE)[[1]]
  } else {
    # two-or-more underscores -> split at the second underscore, preserve rest of the name
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    cell_types <- c(left, right)
  }
  cell_types <- trimws(cell_types)
  message("Processing file: ", file_name, " -> cell_types: ", paste(cell_types, collapse = ", "))

  # ----------------------------------------------------
  # Set working/results directories PER FILE
  # ----------------------------------------------------
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  results_dir <- file.path(working_dir, "Results", file_tag)
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  setwd(working_dir)

  # Set organism and ontology
  organism <- "org.Mm.eg.db"
  ont <- "BP"  # change to "CC", "MF", or "ALL" if desired

  # ----------------------------------------------------
  # Load and prepare data (uses data_path from the loop)
  # ----------------------------------------------------
  df <- read.csv(data_path, header = TRUE)
  colnames(df)[1] <- "gene_symbol"

  original_gene_list <- df$log2fc
  names(original_gene_list) <- df$gene_symbol
  gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  top_df <- df
  colnames(top_df)[1] <- "gene_symbol"
  top_gene_list <- top_df$log2fc
  names(top_gene_list) <- top_df$gene_symbol
  top_gene_list <- sort(na.omit(top_gene_list), decreasing = TRUE)
  top_genes <- names(top_gene_list)[abs(top_gene_list) > 1]
  top_genes <- sort(top_gene_list[top_genes], decreasing = TRUE)
  top_genes <- names(top_genes)

  # ----------------------------------------------------
  # GSEA GO
  # ----------------------------------------------------
  gse <- gseGO(
    geneList = gene_list, ont = "BP", keyType = "UNIPROT",
    minGSSize = 10, maxGSSize = 800, pvalueCutoff = 0.9, verbose = TRUE,
    OrgDb = organism, pAdjustMethod = "BH"
  )

  #gse <- pairwise_termsim(gse)

  plot_dot(gse, cell_types)

  core_dir <- file.path(results_dir, "core_enrichment", ont)
  if (!dir.exists(core_dir)) dir.create(core_dir, recursive = TRUE)
  write.csv(gse@result,
            file = file.path(core_dir, paste0("coreEnrichment_", ont, "_", file_tag, ".csv")))
  


  message("Simplifying GSEA results (cutoff = 0.70) ...")
  gse2 <- tryCatch(
    simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
    error = function(e) {
      message("  Simplify failed: ", conditionMessage(e))
      NULL
    }
  )

  out_file <- file.path(core_dir, paste0("coreEnrichment_simplify_", ont, "_", file_tag, ".csv"))
  if (!is.null(gse2)) {
    write.csv(gse2@result, file = out_file)
    message("  Simplified GSEA results saved to: ", out_file, " (", nrow(gse2@result), " terms)")
  } else {
    message("  No simplified results to save.")
  }

  # check for .sign column
  
  plot_dot(gse2, cell_types)

  # ----------------------------------------------------
  # ORA on top genes
  # ----------------------------------------------------
  ora <- enrichGO(
    gene = top_genes, ont = "BP", keyType = "UNIPROT",
    minGSSize = 10, maxGSSize = 800, pvalueCutoff = 0.9,
    OrgDb = organism, pAdjustMethod = "none"
  )

  plot_dot(ora, cell_types)
  write.csv(ora@result,
            file = file.path(core_dir, paste0("ORA_topGenes_", ont, "_", file_tag, ".csv")))
  
  #ora <- pairwise_termsim(ora)
  # simplify ora results
  message("Simplifying ORA results (cutoff = 0.70) ...")
  ora2 <- tryCatch(
    simplify(ora, cutoff = 0.7, by = "p.adjust", select_fun = min),
    error = function(e) {
      message("  Simplify failed: ", conditionMessage(e))
      NULL
    }
  )
  out_file <- file.path(core_dir, paste0("ORA_topGenes_simplify_", ont, "_", file_tag, ".csv"))
  if (!is.null(ora2)) {
    write.csv(ora2@result, file = out_file)
    message("  Simplified ORA results saved to: ", out_file, " (", nrow(ora2@result), " terms)")
  } else {
    message("  No simplified ORA results to save.")
  }

  plot_dot(ora2, cell_types)

  p7 <- clusterProfiler::dotplot(ora, showCategory = 10) +
    labs(title = "ORA of Top Regulated Genes") +
    scale_x_continuous(limits = c(0, 1))
  save_plot(p7, paste0("ORA_dotplot_", file_tag, ".svg"))

  # ----------------------------------------------------
  # Enrichment map / cnet / ridge / gseaplot
  # ----------------------------------------------------
  save_plot(emapplot(pairwise_termsim(gse), showCategory = 10),
            paste0("GSEAemap_", file_tag, ".svg"))
  save_plot(cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
            paste0("GSEAcnet_", file_tag, ".svg"))

  # also plot simplified version
  save_plot(emapplot(pairwise_termsim(gse2), showCategory = 10),
            paste0("GSEAemap_simplified_", file_tag, ".svg"))
  save_plot(cnetplot(gse2, categorySize = "pvalue", foldChange = gene_list),
            paste0("GSEAcnet_simplified_", file_tag, ".svg"))

  ridgeplot_gse <- ridgeplot(gse) +
    labs(x = "Enrichment Distribution", title = "GSEA Ridgeplot") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  save_plot(ridgeplot_gse, paste0("GSEA_Ridgeplot_", file_tag, ".svg"))

  if (nrow(gse@result) > 0) {
    gseaplot_gse <- gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1)
    save_plot(gseaplot_gse, paste0("GSEA_Plot_", file_tag, ".svg"))
  }

  # PubMed trends
  if (nrow(gse@result) >= 1) {
    top_terms <- head(gse@result$Description, 3)
    pmcplot_gse <- pmcplot(top_terms, 2010:2025, proportion = FALSE) +
      labs(title = "Publication Trends for Top Enriched Terms")
    save_plot(pmcplot_gse, paste0("GSEA_PubMed_Trends_", file_tag, ".svg"))
  }

  # ----------------------------------------------------
  # KEGG GSEA
  # ----------------------------------------------------
  ids <- bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  dedup_ids <- ids[!duplicated(ids$UNIPROT), ]
  df2 <- merge(df, dedup_ids, by.x = "gene_symbol", by.y = "UNIPROT")

  kegg_gene_list <- df2$log2fc
  names(kegg_gene_list) <- df2$ENTREZID
  kegg_gene_list <- kegg_gene_list[!duplicated(names(kegg_gene_list))]
  kegg_gene_list <- sort(na.omit(kegg_gene_list), decreasing = TRUE)

  kegg_organism <- "mmu"
  kk2 <- gseKEGG(
    geneList = kegg_gene_list, organism = "mmu",
    minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
    keyType = "ncbi-geneid"
  )

  plot_dot(kk2, cell_types)

  #kk2 <- pairwise_termsim(kk2)

  # simplify kegg results
  #message("Simplifying KEGG GSEA results (cutoff = 0.70) ...")
  #kk3 <- tryCatch(
  #  simplify(kk2, cutoff = 0.7, by = "p.adjust", select_fun = min),
  #  error = function(e) {
  #    message("  Simplify failed: ", conditionMessage(e))
  #    NULL
  #  }
  #)
  #out_file <- file.path(core_dir, paste0("coreEnrichment_KEGG_simplify_", file_tag, ".csv"))
  #if (!is.null(kk3)) {
  #  write.csv(kk2@result, file = out_file)
  #  message("  Simplified KEGG GSEA results saved to: ", out_file, " (", nrow(kk2@result), " terms)")
  #} else {
  #  message("  No simplified KEGG results to save.")
  #}

  #plot_dot(kk3, cell_types)

  # Optional plots (save if desired)
  save_plot(emapplot(pairwise_termsim(kk2), showCategory = 10),
            paste0("KEGGemap_", file_tag, ".svg"))
  save_plot(cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list),
            paste0("KEGGcnet_", file_tag, ".svg"))
  save_plot(ridgeplot(kk2) + labs(x = "Enrichment distribution"),
            paste0("KEGG_ridge_", file_tag, ".svg"))
  if (nrow(kk2@result) > 0) {
    save_plot(gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1),
              paste0("KEGG_gseaplot_", file_tag, ".svg"))
  }

  # Pathview
  if (!requireNamespace("pathview", quietly = TRUE)) install.packages("pathview")
  library(pathview)
  path_ids <- c("mmu04110","mmu04115","mmu04114","mmu04113","mmu04112","mmu04111",
                "mmu04116","mmu04117","mmu04118","mmu04119","mmu04720","mmu04721",
                "mmu04722","mmu04725","mmu04726","mmu04727","mmu04724","mmu04080",
                "mmu00030","mmu04151")
  pathview_dir <- normalizePath(file.path(results_dir, "pathview"), winslash = "/", mustWork = FALSE)
  if (!dir.exists(pathview_dir)) dir.create(pathview_dir, recursive = TRUE)
  oldwd <- getwd(); setwd(pathview_dir)
  lapply(path_ids, function(pid) {
    pathview(
      gene.data = kegg_gene_list,
      pathway.id = pid,
      species = kegg_organism,
      low = "#6698CC", mid = "white", high = "#F08C21",
      file.type = "svg"
    )
  })
  setwd(oldwd)

  # ----------------------------------------------------
  # GO enrich (ALL) heatmap scaffold (kept, save if needed)
  # ----------------------------------------------------
  go_enrich <- enrichGO(
    gene = names(gene_list), universe = names(gene_list), OrgDb = organism,
    keyType = 'SYMBOL', readable = TRUE, ont = "ALL", pvalueCutoff = 1,
    qvalueCutoff = 1, pAdjustMethod = "none", minGSSize = 3, maxGSSize = 800
  )

  # ----------------------------------------------------
  # KEGG GSEA on predefined UniProt IDs (kept)
  # ----------------------------------------------------
  #selected_uniprot <- c(
  #  "P21279","P21278","Q9Z1B3","P51432","P11881","P68404",
  #  "P63318","P0DP26","P0DP27","P11798","P28652","Q61411",
  #  "Q99N57","P31938","P63085","Q63844","Q8BWG8","Q91YI4","V9GXQ9"
  #)
  #selected_entrez <- bitr(selected_uniprot, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  #selected_df <- merge(df, selected_entrez, by.x = "gene_symbol", by.y = "UNIPROT")
  #selected_kegg_list <- selected_df$log2fc; names(selected_kegg_list) <- selected_df$ENTREZID
  #selected_kegg_list <- sort(unique(na.omit(selected_kegg_list)), decreasing = TRUE)
  #gsea_kegg_selected <- gseKEGG(
  #  geneList = selected_kegg_list, organism = "mmu",
  #  minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1,
  #  pAdjustMethod = "BH", keyType = "ncbi-geneid", verbose = TRUE
  #)
  #plot_dot(gsea_kegg_selected, cell_types)
  #write.csv(gsea_kegg_selected@result,
  #          file = file.path(results_dir, paste0("KEGG_GSEA_Predefined_UniProt_", file_tag, ".csv")))

  # ----------------------------------------------------
  # Custom GSEA (NK3R-signalling)
  # ----------------------------------------------------
  #nk3r_genes <- c("P21279","P21278","P51432","P11881","P63318",
  #                "P68404","P0DP26","P0DP27","P0DP28","P11798",
  #                "P28652","P47937","P47713")
  #term2gene_nk3r <- data.frame(term = rep("NK3R-signalling", length(nk3r_genes)),
  #                             gene = nk3r_genes)

  #custom_gene_list <- df$log2fc; names(custom_gene_list) <- df$gene_symbol
  #custom_gene_list <- sort(na.omit(custom_gene_list), decreasing = TRUE)
  #custom_gene_list <- custom_gene_list[!duplicated(names(custom_gene_list))]

  #gsea_nk3r <- clusterProfiler::GSEA(
  #  geneList = custom_gene_list, TERM2GENE = term2gene_nk3r,
  #  pvalueCutoff = 1, minGSSize = 1, maxGSSize = 500, verbose = TRUE
  #)

  #plot_dot(gsea_nk3r, cell_types)
  #openxlsx::write.xlsx(gsea_nk3r@result,
  #  file = file.path(results_dir, paste0("Custom_GSEA_NK3R_", file_tag, ".xlsx"))
  #)
  #if (nrow(gsea_nk3r@result) > 0) {
  #  save_plot(gseaplot(gsea_nk3r, by = "all", title = "NK3R-signalling", geneSetID = 1),
  #            paste0("Custom_GSEA_NK3R_plot_", file_tag, ".svg"))
  #}

  # ----------------------------------------------------
  # Celltype scoring block (saved under per-file results)
  # ----------------------------------------------------
  #celltype_file <- file.path("S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler",
  #                           "Datasets/celltypes_long.xlsx")
  #library(openxlsx)
  #celltype_df <- read.xlsx(celltype_file, sheet = 1)

  #uniprot_mapping_file_path <- file.path(working_dir, "Datasets", "MOUSE_10090_idmapping.dat")
  #cat("Loading UniProt mapping file from:", uniprot_mapping_file_path, "\n")
  #uniprot_df <- read.delim(uniprot_mapping_file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

 # gene_names <- uniprot_df %>%
 #   dplyr::filter(V2 == "Gene_Name") %>%
 #   dplyr::distinct(V1, .keep_all = TRUE) %>%
 #   dplyr::rename(UniProtID = V1, Gene_Name = V3) %>%
 #   dplyr::select(UniProtID, Gene_Name)

 # uniprot_ids <- uniprot_df %>%
 #   dplyr::filter(V2 == "UniProtKB-ID") %>%
 #   dplyr::distinct(V1, .keep_all = TRUE) %>%
 #   dplyr::rename(UniProtID = V1, UniProtKB_ID = V3) %>%
 #   dplyr::select(UniProtID, UniProtKB_ID)

  #gene_map <- gene_names %>% dplyr::left_join(uniprot_ids, by = "UniProtID")

  #read_gct <- function(file_path) {
  #  gct_data <- read.delim(file_path, skip = 2, header = FALSE, check.names = FALSE)
  #  if ("Name" %in% colnames(gct_data)) colnames(gct_data)[1] <- "Gene"
  #  gct_data
  #}

  #gct_file <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/neha_spatial_group.gct"
  #gct_data <- read_gct(gct_file)

  #grouping_factor <- "Celltype"

  #metadata <- gct_data[1:4, ]
  #colnames(metadata) <- paste0("V", seq_len(ncol(metadata)))

  #sample_id       <- as.character(unlist(gct_data[1, -1]))
  #celltype        <- as.character(unlist(gct_data[2, -1]))
  #groups          <- as.character(unlist(gct_data[3, -1]))
  #celltype_group  <- as.character(unlist(gct_data[4, -1]))

  #metadata_df <- data.frame(
  #  ColName = paste0("V", 2:(ncol(gct_data))),
  #  Sample = sample_id,
  #  Celltype = celltype,
  #  Group = groups,
  #  Celltype_Group = celltype_group,
  #  stringsAsFactors = FALSE
  #)

  #expr_data <- gct_data[-c(1:4), ]
  #colnames(expr_data)[1] <- "Protein_ID"
  #colnames(expr_data)[-1] <- sample_id

  #expr_long <- tidyr::pivot_longer(
  #  expr_data, cols = -Protein_ID, names_to = "Sample", values_to = "Expression"
  #)

  #expr_annotated <- expr_long %>%
  #  dplyr::left_join(metadata_df %>% dplyr::select(-ColName), by = "Sample") %>%
  #  dplyr::mutate(Expression = as.numeric(Expression)) %>%
  #  dplyr::mutate(Protein_ID = sub(";.*", "", Protein_ID)) %>%
  #  dplyr::left_join(gene_map, by = c("Protein_ID" = "UniProtKB_ID")) %>%
  #  dplyr::filter(!grepl("Background", Celltype)) %>%
  #  dplyr::relocate(Gene_Name, .before = 1) %>%
  #  dplyr::relocate(UniProtID, .after = Gene_Name) %>%
  #  dplyr::relocate(Expression, .after = dplyr::last_col())

  #celltype_long <- celltype_df %>%
  #  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Celltype_Class", values_to = "Gene_Label") %>%
  #  dplyr::filter(!is.na(Gene_Label) & Gene_Label != "")

  #marker_expr <- expr_annotated %>%
  #  dplyr::filter(Gene_Name %in% celltype_long$Gene_Label) %>%
  #  dplyr::left_join(celltype_long, by = c("Gene_Name" = "Gene_Label"))

  #celltype_scores <- marker_expr %>%
  #  dplyr::group_by(!!rlang::sym(grouping_factor), Celltype_Class) %>%
  #  dplyr::summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

  #celltype_score_plot <- ggplot2::ggplot(celltype_scores,
  #  ggplot2::aes(x = !!rlang::sym(grouping_factor), y = Mean_Expression, fill = Celltype_Class)) +
  #  ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  #  ggplot2::labs(title = "Celltype Scores by Group", x = grouping_factor, y = "Mean Expression", fill = "Cell Class") +
  #  ggplot2::scale_fill_manual(values = c("#f36d07", "#455A64", "#c7c7c7")) +
  #  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
  #  ggplot2::theme_bw(base_size = 16) +
  #  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 20),
  #                 axis.title = ggplot2::element_text(face = "bold", size = 18),
  #                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, color = "black", size = 16),
  #                 axis.text.y = ggplot2::element_text(color = "black", size = 16),
  #                 panel.grid.major = ggplot2::element_blank(),
  #                 panel.grid.minor = ggplot2::element_blank(),
  #                 panel.border = ggplot2::element_blank(),
  #                 axis.ticks = ggplot2::element_line(color = "black", size = 1),
  #                 legend.position = c(0.9, 0.9),
  #                 legend.background = ggplot2::element_rect(color = "black", fill = NA),
  #                 legend.title = ggplot2::element_text(face = "bold", size = 16),
  #                 legend.text = ggplot2::element_text(size = 16))

  #ggplot2::ggsave(file.path(results_dir, paste0("celltype_scores_by_group_", file_tag, ".svg")),
  #                celltype_score_plot, units = "cm", dpi = 300)

  #gene_class_df <- celltype_long %>% dplyr::distinct(Gene_Label, Celltype_Class) %>% dplyr::rename(Gene_Name = Gene_Label)
  #annotation_row <- data.frame(Celltype_Class = gene_class_df$Celltype_Class)
  #rownames(annotation_row) <- gene_class_df$Gene_Name

  #marker_matrix <- marker_expr %>%
  #  dplyr::group_by(Gene_Name, !!rlang::sym(grouping_factor)) %>%
  #  dplyr::summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  #  tidyr::pivot_wider(names_from = !!rlang::sym(grouping_factor), values_from = Mean_Expression) %>%
  #  tibble::column_to_rownames("Gene_Name") %>%
  #  as.matrix()

  #marker_matrix_clean <- marker_matrix[rowSums(is.na(marker_matrix)) < ncol(marker_matrix), ]
  #marker_matrix_clean <- marker_matrix_clean[, colSums(is.na(marker_matrix_clean)) < nrow(marker_matrix_clean)]
  #marker_matrix_clean <- marker_matrix_clean[rowSums(!is.na(marker_matrix_clean)) >= 2, ]
  #marker_matrix_clean <- marker_matrix_clean[, colSums(!is.na(marker_matrix_clean)) >= 2]
  #marker_matrix_clean <- marker_matrix_clean[!grepl("^Oasl2\\s*$", rownames(marker_matrix_clean)), ]
  #marker_matrix_clean[is.nan(marker_matrix_clean)] <- NA
  #marker_matrix_clean[is.infinite(marker_matrix_clean)] <- NA

  #annotation_row <- annotation_row[rownames(marker_matrix_clean), , drop = FALSE]
  #annotation_colors <- list(
  #  Celltype_Class = c(gaba = "#f36d07", vglut1 = "#455A64", vglut2 = "#c7c7c7")
  #)

  #pheatmap::pheatmap(
  #  marker_matrix_clean, cluster_rows = TRUE, cluster_cols = FALSE,
  #  color = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
  #  na_col = "grey", main = paste("Marker Expression per", grouping_factor),
  #  fontsize = 12, fontsize_row = 10, fontsize_col = 12,
  #  cellheight = 12, cellwidth = 12, border_color = NA,
  #  annotation_row = annotation_row, annotation_colors = annotation_colors,
  #  filename = file.path(results_dir, paste0("marker_expression_heatmap_", grouping_factor, "_", file_tag, ".pdf"))
  #)

  #expr_scaled <- marker_expr %>%
  #  dplyr::group_by(Gene_Name) %>%
  #  dplyr::mutate(Z_Expression = scale(Expression)) %>%
  #  dplyr::ungroup()

  #celltype_z_scores <- expr_scaled %>%
  #  dplyr::group_by(!!rlang::sym(grouping_factor), Celltype_Class) %>%
  #  dplyr::summarise(Mean_Z = mean(Z_Expression, na.rm = TRUE), .groups = "drop")

  #heatmap_plot <- ggplot2::ggplot(
  #  celltype_z_scores,
  #  ggplot2::aes(x = reorder(Celltype_Class, Mean_Z, FUN = median), y = !!rlang::sym(grouping_factor), fill = Mean_Z)
  #) +
  #  ggplot2::geom_tile(color = "grey80", size = 0.2, width = 0.5) +
  #  ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
  #                                midpoint = 0, name = "Z-Score",
  #                                guide = ggplot2::guide_colorbar(barwidth = 0.8, barheight = 20)) +
  #  ggplot2::labs(title = "Celltype Marker Signature by Group", x = "Marker Class", y = "Sample Group") +
  #  ggplot2::theme_minimal(base_size = 12) +
  #  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, margin = margin(b = 10)),
  #                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
  #                 axis.text.y = ggplot2::element_text(size = 10),
  #                 axis.title = ggplot2::element_text(size = 12),
  #                 panel.grid = ggplot2::element_blank(),
  #                 legend.title = ggplot2::element_text(size = 12),
  #                 legend.text = ggplot2::element_text(size = 10))

  #ggplot2::ggsave(file.path(results_dir, paste0("celltype_scores_heatmap_", file_tag, ".svg")),
  #                heatmap_plot, width = 16, height = 9, units = "cm", dpi = 300)

  #proportional_scores <- expr_scaled %>%
  #  dplyr::group_by(!!rlang::sym(grouping_factor), Celltype_Class) %>%
  #  dplyr::summarise(Mean_Z = mean(Z_Expression, na.rm = TRUE), .groups = "drop") %>%
  #  dplyr::group_by(!!rlang::sym(grouping_factor)) %>%
  #  dplyr::mutate(Proportion = Mean_Z / sum(Mean_Z, na.rm = TRUE) * 100) %>%
  #  dplyr::ungroup()

  #prop_barplot <- ggplot2::ggplot(proportional_scores,
  #  ggplot2::aes(x = !!rlang::sym(grouping_factor), y = Proportion, fill = Celltype_Class)) +
  #  ggplot2::geom_bar(stat = "identity", width = 0.8) +
  #  ggplot2::labs(
  #    title = "Proportional Celltype Signature per Group",
  #    x = "Celltype Group", y = "Proportion of Total Marker Signal (%)", fill = "Celltype Class"
  #  ) +
  #  ggplot2::scale_fill_brewer(palette = "Set2") +
  #  ggplot2::theme_minimal(base_size = 14) +
  #  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
  #                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
  #                 panel.grid = ggplot2::element_line(color = "#f0f0f0"))

  #ggplot2::ggsave(file.path(results_dir, paste0("celltype_relative_proportions_", file_tag, ".svg")),
  #                prop_barplot, width = 16, height = 10, units = "cm", dpi = 300)

  # ---- END single-file processing ----
}





























# ---- Parallel backend setup ----
library(foreach)
library(doParallel)

n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
cl <- parallel::makeCluster(n_cores)        # Windows-safe cluster backend
doParallel::registerDoParallel(cl)

# Ensure required packages are installed/loaded on master BEFORE spawning workers
setupPackages()

# Define constants & helpers outside the loop so they can be exported
working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
file_direction <- file.path(working_dir, "Datasets", "mapped", "neuron-phenotypeWithinUnit")

save_plot <- function(plot, filename) {
  ggsave(file.path(results_dir, filename), plot, units = "cm", dpi = 300)
}
plot_dot <- function(dataset, cell_types, results_dir) {
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
  ggsave(file.path(results_dir, plot_filename), p, units = "cm", dpi = 300)
  p
}

# Enumerate input files
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

# ---- Parallel loop over files ----
res_list <- foreach(
  data_path = file_list,
  .packages = c(
    "clusterProfiler","pathview","enrichplot","DOSE","ggplot2","ggnewscale",
    "cowplot","ggridges","europepmc","ggpubr","ggrepel","ggsci","ggthemes",
    "ggExtra","ggforce","ggalluvial","lattice","latticeExtra","org.Mm.eg.db",
    "ggplotify","svglite","tidyr","dplyr","pheatmap","proxy","tibble"
  ),
  .export = c("working_dir","save_plot","plot_dot"),
  .errorhandling = "pass",
  .inorder = FALSE
) %dopar% {

  # Per-task libraries that sometimes need explicit attach on workers
  suppressPackageStartupMessages({
    library(clusterProfiler); library(enrichplot); library(ggplot2)
    library(DOSE); library(dplyr); library(tidyr); library(tibble)
  })

  # Derive names
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # cell_types parsing: split at second underscore if present
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[21]]
  if (identical(us[21], -1L)) {
    cell_types <- c(name_no_ext)
  } else if (length(us) == 1L) {
    cell_types <- strsplit(name_no_ext, "_", fixed = TRUE)[[21]]
  } else {
    pos2 <- us[22]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    cell_types <- c(left, right)
  }
  cell_types <- trimws(cell_types)
  file_tag   <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)

  # Results dir for this file (avoid setwd in parallel)
  results_dir <- file.path(working_dir, "Results", file_tag)
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  organism <- "org.Mm.eg.db"
  ont <- "BP"

  # Load data
  df <- read.csv(data_path, header = TRUE)
  colnames(df)[21] <- "gene_symbol"   # UniProt accession column

  # Ranked list by UniProt accession
  original_gene_list <- df$log2fc
  names(original_gene_list) <- df$gene_symbol
  gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  # Top genes for ORA
  top_gene_list <- sort(na.omit(df$log2fc), decreasing = TRUE)
  names(top_gene_list) <- df$gene_symbol
  top_genes <- names(top_gene_list)[abs(top_gene_list) > 1]
  top_genes <- names(sort(top_gene_list[top_genes], decreasing = TRUE))

  # ---- GO GSEA with UniProt keyType ----
  gse <- tryCatch({
    gseGO(
      geneList = gene_list, ont = "CC", keyType = "UNIPROT",
      minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, verbose = FALSE,
      OrgDb = organism, pAdjustMethod = "BH"
    )
  }, error = function(e) NULL)

  if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
    plot_dot(gse, cell_types, results_dir)
    core_dir <- file.path(results_dir, "core_enrichment", ont)
    if (!dir.exists(core_dir)) dir.create(core_dir, recursive = TRUE)
    write.csv(gse@result, file = file.path(core_dir, paste0("coreEnrichment_", ont, "_", file_tag, ".csv")))
    # Selected plots
    p_emap <- emapplot(pairwise_termsim(gse), showCategory = 10)
    ggsave(file.path(results_dir, paste0("GSEAemap_", file_tag, ".svg")), p_emap, units = "cm", dpi = 300)
    p_cnet <- cnetplot(gse, categorySize = "pvalue", foldChange = gene_list)
    ggsave(file.path(results_dir, paste0("GSEAcnet_", file_tag, ".svg")), p_cnet, units = "cm", dpi = 300)
    p_ridge <- ridgeplot(gse) + labs(x = "Enrichment Distribution", title = "GSEA Ridgeplot") +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    ggsave(file.path(results_dir, paste0("GSEA_Ridgeplot_", file_tag, ".svg")), p_ridge, units = "cm", dpi = 300)
    if (nrow(gse@result) > 0) {
      p_gsea <- gseaplot(gse, by = "all", title = gse@result$Description[21], geneSetID = 1)
      ggsave(file.path(results_dir, paste0("GSEA_Plot_", file_tag, ".svg")), p_gsea, units = "cm", dpi = 300)
      top_terms <- head(gse@result$Description, 3)
      p_pmc <- pmcplot(top_terms, 2010:2025, proportion = FALSE) + labs(title = "Publication Trends for Top Enriched Terms")
      ggsave(file.path(results_dir, paste0("GSEA_PubMed_Trends_", file_tag, ".svg")), p_pmc, units = "cm", dpi = 300)
    }
  }

  # ---- ORA on top genes (UniProt) ----
  ora <- tryCatch({
    enrichGO(
      gene = top_genes, ont = "CC", keyType = "UNIPROT",
      minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1,
      OrgDb = organism, pAdjustMethod = "none"
    )
  }, error = function(e) NULL)
  if (!is.null(ora) && nrow(as.data.frame(ora)) > 0) {
    p7 <- clusterProfiler::dotplot(ora, showCategory = 10) +
      labs(title = "ORA of Top Regulated Genes") +
      scale_x_continuous(limits = c(0, 1))
    ggsave(file.path(results_dir, paste0("ORA_dotplot_", file_tag, ".svg")), p7, units = "cm", dpi = 300)
  }

  # ---- KEGG GSEA (map UniProt -> ENTREZID) ----
  ids <- tryCatch({
    bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  }, error = function(e) NULL)

  kk2 <- NULL
  if (!is.null(ids) && nrow(ids) > 0) {
    ids <- ids[!is.na(ids$ENTREZID), ]
    ids <- ids[!duplicated(ids$UNIPROT), ]
    if (nrow(ids) > 0) {
      df2 <- merge(df, ids, by.x = "gene_symbol", by.y = "UNIPROT")
      kegg_gene_list <- df2$log2fc
      names(kegg_gene_list) <- df2$ENTREZID
      kegg_gene_list <- kegg_gene_list[!duplicated(names(kegg_gene_list))]
      kegg_gene_list <- sort(na.omit(kegg_gene_list), decreasing = TRUE)
      if (length(kegg_gene_list) > 0) {
        kk2 <- tryCatch({
          gseKEGG(
            geneList = kegg_gene_list, organism = "mmu",
            minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
            keyType = "ncbi-geneid"
          )
        }, error = function(e) NULL)
      }
    }
  }

  if (!is.null(kk2) && nrow(as.data.frame(kk2)) > 0) {
    plot_dot(kk2, cell_types, results_dir)
    p_emap2 <- emapplot(pairwise_termsim(kk2), showCategory = 10)
    ggsave(file.path(results_dir, paste0("KEGGemap_", file_tag, ".svg")), p_emap2, units = "cm", dpi = 300)
    p_cnet2 <- cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list)
    ggsave(file.path(results_dir, paste0("KEGGcnet_", file_tag, ".svg")), p_cnet2, units = "cm", dpi = 300)
    p_ridge2 <- ridgeplot(kk2) + labs(x = "Enrichment distribution")
    ggsave(file.path(results_dir, paste0("KEGG_ridge_", file_tag, ".svg")), p_ridge2, units = "cm", dpi = 300)
    if (nrow(kk2@result) > 0) {
      p_gsea2 <- gseaplot(kk2, by = "all", title = kk2$Description[21], geneSetID = 1)
      ggsave(file.path(results_dir, paste0("KEGG_gseaplot_", file_tag, ".svg")), p_gsea2, units = "cm", dpi = 300)
    }

    # Pathview block (avoid setwd in parallel; use out.suffix to differentiate)
    if (!requireNamespace("pathview", quietly = TRUE)) install.packages("pathview")
    pv_dir <- normalizePath(file.path(results_dir, "pathview"), winslash = "/", mustWork = FALSE)
    if (!dir.exists(pv_dir)) dir.create(pv_dir, recursive = TRUE)
    path_ids <- c("mmu04110","mmu04115","mmu04114","mmu04113","mmu04112","mmu04111",
                  "mmu04116","mmu04117","mmu04118","mmu04119","mmu04720","mmu04721",
                  "mmu04722","mmu04725","mmu04726","mmu04727","mmu04724","mmu04080",
                  "mmu00030","mmu04151")
    invisible(lapply(path_ids, function(pid) {
      pathview::pathview(
        gene.data  = kegg_gene_list,
        pathway.id = pid,
        species    = "mmu",
        out.suffix = file_tag,         # ensure unique file names per worker
        kegg.native = TRUE,
        low = "#6698CC", mid = "white", high = "#F08C21",
        file.type = "svg",
        kegg.dir = pv_dir,             # output directory
        out.dir  = pv_dir
      )
    }))
  }

  # Return a minimal status object for each file
  list(file = file_name,
       go_terms = if (!is.null(gse)) nrow(as.data.frame(gse)) else 0L,
       kegg_terms = if (!is.null(kk2)) nrow(as.data.frame(kk2)) else 0L)
}

# Stop cluster when done
parallel::stopCluster(cl)




























# =========================
# Parallel GSEA pipeline (Option A, final)
# =========================

# ---- Install & load packages on master ----
setupPackages <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  cran_pkgs <- c(
    "ggplot2","ggnewscale","cowplot","ggridges","europepmc","ggpubr","ggrepel",
    "ggsci","ggthemes","ggExtra","ggforce","ggalluvial","lattice","latticeExtra",
    "ggplotify","svglite","tidyr","dplyr","pheatmap","proxy","tibble",
    "foreach","doParallel"
  )
  missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_cran)) install.packages(missing_cran, dependencies = TRUE)

  bioc_pkgs <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","pathview")
  missing_bioc <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

  invisible(lapply(c(cran_pkgs, bioc_pkgs), function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))
}
setupPackages()

# ---- Setup parallel backend ----
library(foreach)
library(doParallel)
n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# ---- Constants & helpers ----
working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
file_direction <- file.path(working_dir, "Datasets", "mapped", "neuron-phenotypeWithinUnit")

save_plot <- function(plot, filename, results_dir) {
  ggplot2::ggsave(file.path(results_dir, filename), plot, units = "cm", dpi = 300)
}

plot_dot <- function(dataset, cell_types, results_dir) {
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
  ggplot2::ggsave(file.path(results_dir, plot_filename), p, units = "cm", dpi = 300)
  p
}

# ---- Enumerate input files ----
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

# ---- Parallel loop ----
res_list <- foreach(
  data_path = file_list,
  .packages = c(
    "clusterProfiler","pathview","enrichplot","DOSE","ggplot2","ggnewscale",
    "cowplot","ggridges","europepmc","ggpubr","ggrepel","ggsci","ggthemes",
    "ggExtra","ggforce","ggalluvial","lattice","latticeExtra","org.Mm.eg.db",
    "ggplotify","svglite","tidyr","dplyr","pheatmap","proxy","tibble"
  ),
  .export = c("working_dir","save_plot","plot_dot"),
  .errorhandling = "pass",
  .inorder = FALSE
) %dopar% {

  # Per-worker package attach to open independent SQLite connections
  suppressPackageStartupMessages({
    library(clusterProfiler); library(enrichplot); library(ggplot2)
    library(DOSE); library(dplyr); library(tidyr); library(tibble)
    library(org.Mm.eg.db)
  })

  # ---- Names & outputs ----
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Split at second underscore using proper gregexpr indexing
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[21]]
  if (identical(us[21], -1L)) {
    cell_types <- c(name_no_ext)
  } else if (length(us) == 1L) {
    cell_types <- strsplit(name_no_ext, "_", fixed = TRUE)[[21]]
  } else {
    pos2 <- us[22]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    cell_types <- c(left, right)
  }
  cell_types <- trimws(cell_types)
  file_tag   <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)

  results_dir <- file.path(working_dir, "Results", file_tag)
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  organism <- "org.Mm.eg.db"
  ont <- "BP"

  # ---- Load input ----
  df <- read.csv(data_path, header = TRUE, check.names = FALSE)
  # Rename first column to UniProt "gene_symbol"
  colnames(df)[21] <- "gene_symbol"

  # ---- Ranked list (UniProt) ----
  original_gene_list <- df$log2fc
  names(original_gene_list) <- df$gene_symbol
  gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  # ---- Top genes for ORA ----
  top_gene_list <- sort(na.omit(df$log2fc), decreasing = TRUE)
  names(top_gene_list) <- df$gene_symbol
  top_genes <- names(top_gene_list)[abs(top_gene_list) > 1]
  top_genes <- names(sort(top_gene_list[top_genes], decreasing = TRUE))

  # ---- GO GSEA (UniProt) ----
  gse <- tryCatch({
    gseGO(
      geneList = gene_list, ont = "CC", keyType = "UNIPROT",
      minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, verbose = FALSE,
      OrgDb = organism, pAdjustMethod = "BH"
    )
  }, error = function(e) NULL)

  if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
    plot_dot(gse, cell_types, results_dir)
    core_dir <- file.path(results_dir, "core_enrichment", ont)
    if (!dir.exists(core_dir)) dir.create(core_dir, recursive = TRUE)
    write.csv(gse@result, file = file.path(core_dir, paste0("coreEnrichment_", ont, "_", file_tag, ".csv")))

    p_emap <- enrichplot::emapplot(pairwise_termsim(gse), showCategory = 10)
    ggplot2::ggsave(file.path(results_dir, paste0("GSEAemap_", file_tag, ".svg")), p_emap, units = "cm", dpi = 300)

    p_cnet <- enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list)
    ggplot2::ggsave(file.path(results_dir, paste0("GSEAcnet_", file_tag, ".svg")), p_cnet, units = "cm", dpi = 300)

    p_ridge <- enrichplot::ridgeplot(gse) +
      labs(x = "Enrichment Distribution", title = "GSEA Ridgeplot") +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    ggplot2::ggsave(file.path(results_dir, paste0("GSEA_Ridgeplot_", file_tag, ".svg")), p_ridge, units = "cm", dpi = 300)

    if (nrow(gse@result) > 0) {
      p_gsea <- enrichplot::gseaplot(gse, by = "all", title = gse@result$Description[21], geneSetID = 1)
      ggplot2::ggsave(file.path(results_dir, paste0("GSEA_Plot_", file_tag, ".svg")), p_gsea, units = "cm", dpi = 300)

      top_terms <- head(gse@result$Description, 3)
      p_pmc <- europepmc::pmcplot(top_terms, 2010:2025, proportion = FALSE) + labs(title = "Publication Trends for Top Enriched Terms")
      ggplot2::ggsave(file.path(results_dir, paste0("GSEA_PubMed_Trends_", file_tag, ".svg")), p_pmc, units = "cm", dpi = 300)
    }
  }

  # ---- ORA (UniProt) ----
  ora <- tryCatch({
    enrichGO(
      gene = top_genes, ont = "CC", keyType = "UNIPROT",
      minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1,
      OrgDb = organism, pAdjustMethod = "none"
    )
  }, error = function(e) NULL)

  if (!is.null(ora) && nrow(as.data.frame(ora)) > 0) {
    p7 <- clusterProfiler::dotplot(ora, showCategory = 10) +
      labs(title = "ORA of Top Regulated Genes") +
      scale_x_continuous(limits = c(0, 1))
    ggplot2::ggsave(file.path(results_dir, paste0("ORA_dotplot_", file_tag, ".svg")), p7, units = "cm", dpi = 300)
  }

  # ---- KEGG GSEA (UniProt -> ENTREZID) ----
  kk2 <- NULL
  kegg_gene_list <- NULL

  ids <- tryCatch({
    bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  }, error = function(e) NULL)

  if (!is.null(ids) && nrow(ids) > 0) {
    ids <- ids[!is.na(ids$ENTREZID), ]
    ids <- ids[!duplicated(ids$UNIPROT), ]
    if (nrow(ids) > 0) {
      df2 <- merge(df, ids, by.x = "gene_symbol", by.y = "UNIPROT")

      kegg_gene_list <- df2$log2fc
      names(kegg_gene_list) <- df2$ENTREZID
      kegg_gene_list <- kegg_gene_list[!duplicated(names(kegg_gene_list))]
      kegg_gene_list <- sort(na.omit(kegg_gene_list), decreasing = TRUE)

      if (length(kegg_gene_list) > 0) {
        kk2 <- tryCatch({
          gseKEGG(
            geneList = kegg_gene_list, organism = "mmu",
            minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
            keyType = "ncbi-geneid"
          )
        }, error = function(e) NULL)
      }
    }
  }

  if (!is.null(kk2) && nrow(as.data.frame(kk2)) > 0) {
    plot_dot(kk2, cell_types, results_dir)

    p_emap2 <- enrichplot::emapplot(pairwise_termsim(kk2), showCategory = 10)
    ggplot2::ggsave(file.path(results_dir, paste0("KEGGemap_", file_tag, ".svg")), p_emap2, units = "cm", dpi = 300)

    p_cnet2 <- enrichplot::cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list)
    ggplot2::ggsave(file.path(results_dir, paste0("KEGGcnet_", file_tag, ".svg")), p_cnet2, units = "cm", dpi = 300)

    p_ridge2 <- enrichplot::ridgeplot(kk2) + labs(x = "Enrichment distribution")
    ggplot2::ggsave(file.path(results_dir, paste0("KEGG_ridge_", file_tag, ".svg")), p_ridge2, units = "cm", dpi = 300)

    if (nrow(kk2@result) > 0) {
      p_gsea2 <- enrichplot::gseaplot(kk2, by = "all", title = kk2$Description[21], geneSetID = 1)
      ggplot2::ggsave(file.path(results_dir, paste0("KEGG_gseaplot_", file_tag, ".svg")), p_gsea2, units = "cm", dpi = 300)
    }

    # ---- Pathview: stagger & guard to avoid DB-lock contention ----
    if (!requireNamespace("pathview", quietly = TRUE)) {
      install.packages("pathview")
      suppressPackageStartupMessages(library(pathview))
    } else {
      suppressPackageStartupMessages(library(pathview))
    }

    # Stagger worker starts a bit
    Sys.sleep(runif(1, 0, 1.0))

    pv_dir <- normalizePath(file.path(results_dir, "pathview"), winslash = "/", mustWork = FALSE)
    if (!dir.exists(pv_dir)) dir.create(pv_dir, recursive = TRUE)

    path_ids <- c(
      "mmu04110","mmu04115","mmu04114","mmu04113","mmu04112","mmu04111",
      "mmu04116","mmu04117","mmu04118","mmu04119","mmu04720","mmu04721",
      "mmu04722","mmu04725","mmu04726","mmu04727","mmu04724","mmu04080",
      "mmu00030","mmu04151"
    )

    invisible(lapply(path_ids, function(pid) {
      try({
        pathview::pathview(
          gene.data   = kegg_gene_list,
          pathway.id  = pid,
          species     = "mmu",
          out.suffix  = file_tag,
          kegg.native = TRUE,
          low = "#6698CC", mid = "white", high = "#F08C21",
          file.type   = "svg",
          kegg.dir    = pv_dir,
          out.dir     = pv_dir
        )
      }, silent = TRUE)
    }))
  }

  # ---- Return per-file status ----
  list(file = file_name,
       go_terms = if (!is.null(gse)) nrow(as.data.frame(gse)) else 0L,
       kegg_terms = if (!is.null(kk2)) nrow(as.data.frame(kk2)) else 0L)
}

# ---- Teardown parallel backend ----
parallel::stopCluster(cl)

