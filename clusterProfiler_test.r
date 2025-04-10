# Introduction: Gene Set Enrichment Analysis (GSEA) using clusterProfiler
# Documentation: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

# Check and install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "pathview", "enrichplot", "DOSE"))

# Load libraries with installation if needed
required_packages <- c("clusterProfiler", "enrichplot", "ggplot2", "ggnewscale", 
                       "cowplot", "ggridges", "europepmc", "ggpubr", "ggrepel", 
                       "ggsci", "ggthemes", "ggExtra", "ggforce", "ggalluvial", 
                       "lattice", "latticeExtra", "DOSE")

install_and_load <- function(package) {
    if (!require(package, character.only = TRUE)) {
        install.packages(package)
        library(package, character.only = TRUE)
    }
}

lapply(required_packages, install_and_load)

# Set annotation package and load it
organism <- "org.Mm.eg.db" 
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# Set directories
#working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp3_Nlgn3_development/LaserDissProteomics/GSEA"
working_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/"
results_dir <- file.path(working_dir, "Results")
setwd(working_dir)

# Define cell types
group <- "g1"
cell_types <- c("neuron", "cfos")
file_name <- paste(group, paste(cell_types, collapse = "_"), sep = "_")
data_path <- file.path(working_dir, "Datasets", paste0(file_name, ".csv"))
#data_path <- file.path(working_dir, "Datasets", paste(group, paste(cell_types, collapse = "_"), sep = "_"), ".csv")

# Load and prepare gene data
df <- read.csv(data_path, header = TRUE)
colnames(df)[1] <- "gene_symbol"
original_gene_list <- df$log2fc
names(original_gene_list) <- df$gene_symbol
gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)

# Gene Set Enrichment Analysis (GSEA)
gse <- gseGO(
    geneList = gene_list, ont = "ALL", keyType = "SYMBOL", nPerm = 10000,
    minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, verbose = TRUE,
    OrgDb = organism, pAdjustMethod = "none"
)

# Plotting function for reuse
save_plot <- function(plot, filename) {
    ggsave(file.path(results_dir, filename), plot, units = "cm", dpi = 300)
}

# Dotplot
dot_title <- paste("GSEA of", paste(cell_types, collapse = " over "))
p1 <- dotplot(gse, showCategory = 10, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = dot_title)
save_plot(p1, paste("GSEAdotplot_", paste(cell_types, collapse = "_"), ".png"))

# Enrichment Map and Network Plot
emapplot(pairwise_termsim(gse), showCategory = 10)
cnetplot(gse, categorySize = "pvalue", foldChange = gene_list)

# Ridgeplot
ridgeplot(gse) + labs(x = "Enrichment distribution")

# GSEA Plot
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

# PubMed Trend Plot
terms <- gse$Description[1:3]
pmcplot(terms, 2010:2018, proportion = FALSE)

# KEGG GSEA Analysis
ids <- bitr(names(original_gene_list), fromType = "SYMBOL", toType = "UNIPROT", OrgDb = organism)
dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]), ]
df2 <- df[df$gene_symbol %in% dedup_ids$SYMBOL, ]
df2$gene_symbol <- dedup_ids$UNIPROT
kegg_gene_list <- sort(na.omit(df2$log2fc), decreasing = TRUE)
names(kegg_gene_list) <- df2$gene_symbol

kegg_organism <- "mmu"
kk2 <- gseKEGG(
    geneList = kegg_gene_list, organism = kegg_organism, nPerm = 10000,
    minGSSize = 3, maxGSSize = 800, pvalueCutoff = 0.05, pAdjustMethod = "none",
    keyType = "uniprot"
)

# KEGG Dotplot and Network Plot
kegg_dot_title <- paste("KEGG GSEA Enriched Pathways of", paste(cell_types, collapse = " over "))
p2 <- dotplot(kk2, showCategory = 10, title = kegg_dot_title, split = ".sign") +
    facet_grid(. ~ .sign)
save_plot(p2, paste("KEGGpathway_", paste(cell_types, collapse = "_"), ".png"))

emapplot(kk2)
cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list)
ridgeplot(kk2) + labs(x = "Enrichment distribution")
gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

# KEGG Pathway Visualizations using Pathview
if (!requireNamespace("pathview", quietly = TRUE)) install.packages("pathview")
library(pathview)

pathview(gene.data = kegg_gene_list, pathway.id = "mmu00030", species = kegg_organism)
pathview(gene.data = kegg_gene_list, pathway.id = "mmu00030", species = kegg_organism, kegg.native = FALSE)

# EnrichGO Analysis and Heatmap
go_enrich <- enrichGO(
    gene = names(gene_list), universe = names(gene_list), OrgDb = organism, 
    keyType = 'SYMBOL', readable = TRUE, ont = "CC", pvalueCutoff = 0.05, 
    qvalueCutoff = 0.10
)

# Heatmap Plot
p3 <- heatplot(go_enrich, foldChange = gene_list, showCategory = 5)
cowplot::plot_grid(p1, p3, ncol = 1, labels = LETTERS[1:2])
