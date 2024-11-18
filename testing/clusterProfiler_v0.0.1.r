# Gene Set Enrichment Analysis (GSEA) for Neuron Types
# Documentation: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

# This script performs Gene Set Enrichment Analysis (GSEA) for neuron types using the clusterProfiler package.
# Steps included:
# 1. Install and load necessary Bioconductor and R packages.
# 2. Load the dataset and filter for neuron-type-specific markers.
# 3. Prepare the gene list for enrichment analysis.
# 4. Perform Gene Ontology (GO) enrichment analysis.
# 5. Visualize GO enrichment results using dotplot, ridgeplot, and heatplot.
# 6. Perform KEGG pathway analysis.
# 7. Visualize KEGG pathway analysis results using dotplot, emapplot, and cnetplot.
# 8. Visualize specific KEGG pathways using the pathview package.
# 9. Save the generated plots to the results directory.

# Note: Ensure that the working directory and dataset paths are correctly set before running the script.

# Install Bioconductor and required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "pathview", "enrichplot", "DOSE"))

# Additional R packages
required_packages <- c("ggplot2", "ggnewscale", "cowplot", "ggridges", 
                       "ggpubr", "ggrepel", "ggsci", "ggthemes", "ggExtra", 
                       "ggforce", "ggalluvial", "lattice", "latticeExtra")
install_and_load <- function(package) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
        library(package, character.only = TRUE)
    }
}
lapply(required_packages, install_and_load)

# Load Bioconductor annotation package for the organism (mouse example)
organism <- "org.Mm.eg.db"  # For mouse
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# Set directories
working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp3_Nlgn3_development/LaserDissProteomics/GSEA"
results_dir <- file.path(working_dir, "Results")
setwd(working_dir)

# Load dataset
data_path <- file.path(working_dir, "Datasets", "your_dataset.csv")
df <- read.csv(data_path, header = TRUE)
colnames(df)[1] <- "gene_symbol"

# Focus on neuron-type-specific markers
neuro_markers <- c("GAD1", "GAD2", "SLC32A1", "SLC17A7", "SLC17A6", "GRIA1", "GRIN1")
filtered_df <- df[df$gene_symbol %in% neuro_markers, ]

# Prepare gene list for enrichment analysis
filtered_df <- filtered_df[!is.na(filtered_df$log2fc), ]
gene_list <- sort(filtered_df$log2fc, decreasing = TRUE)
names(gene_list) <- filtered_df$gene_symbol

# Perform Gene Ontology (GO) Enrichment Analysis
go_enrich <- enrichGO(
    gene = names(gene_list), OrgDb = organism, keyType = "SYMBOL",
    ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE
)

# Visualize GO Enrichment results
library(enrichplot)
dotplot(go_enrich, showCategory = 10) + labs(title = "GO Enrichment for Neuron Types")
ridgeplot(go_enrich) + labs(x = "Enrichment distribution")
heatplot(go_enrich, foldChange = gene_list, showCategory = 5)

# Perform KEGG Pathway Analysis
library(clusterProfiler)
kegg_organism <- "mmu"  # Mouse KEGG prefix
kegg_results <- gseKEGG(
    geneList = gene_list, organism = kegg_organism, nPerm = 10000,
    minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, pAdjustMethod = "none"
)

# Visualize KEGG Pathway Analysis results
dotplot(kegg_results, showCategory = 10, title = "KEGG GSEA for Neurotransmitter Pathways")
emapplot(kegg_results)
cnetplot(kegg_results, categorySize = "pvalue", foldChange = gene_list)

# Visualize specific KEGG Pathways using Pathview
library(pathview)
pathview(
    gene.data = gene_list, pathway.id = "mmu04727",  # GABAergic synapse
    species = kegg_organism
)
pathview(
    gene.data = gene_list, pathway.id = "mmu04724",  # Glutamatergic synapse
    species = kegg_organism
)

# Save the generated plots to the results directory
save_plot <- function(plot, filename) {
    ggsave(file.path(results_dir, filename), plot, units = "cm", dpi = 300)
}
save_plot(dotplot(go_enrich, showCategory = 10), "GO_Enrichment.png")
save_plot(dotplot(kegg_results, showCategory = 10), "KEGG_Enrichment.png")

# Save Ridgeplot of KEGG Pathways
ridgeplot(kegg_results) + labs(x = "Enrichment distribution")
save_plot(ridgeplot(kegg_results), "KEGG_Ridgeplot.png")

# Additional Synapse Pathways (GABAergic, Glutamatergic)
pathview(
    gene.data = gene_list, pathway.id = "mmu04727",  # GABAergic synapse
    species = kegg_organism, kegg.native = TRUE
)
pathview(
    gene.data = gene_list, pathway.id = "mmu04724",  # Glutamatergic synapse
    species = kegg_organism, kegg.native = TRUE
)