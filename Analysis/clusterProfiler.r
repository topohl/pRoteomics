#' ============================================================
#' HIGH-THROUGHPUT Gene Set Enrichment Analysis (GSEA) Workflow
#' WITH PARALLEL PROCESSING
#' ============================================================
#' 
#' This script performs comprehensive GSEA, KEGG, ORA, custom pathway analysis
#' for MULTIPLE cell type comparisons in parallel using doParallel.
#' Celltype scoring runs separately at the end.
#'
#' @author Tobias Pohl
#' ============================================================

#' ============================================================
#' HIGH-THROUGHPUT GSEA Workflow with FUTURE (BETTER)
#' ============================================================

# ----------------------------------------------------
# 1. PACKAGE SETUP
# ----------------------------------------------------
checkBiocManager <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
}

installBioC <- function(bioc_packages) {
  missing <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) BiocManager::install(missing)
}

installAndLoadCRAN <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
}

setupPackages <- function() {
  checkBiocManager()
  required_packages <- c("clusterProfiler", "pathview", "enrichplot", "DOSE", "ggplot2", "ggnewscale",
                         "cowplot", "ggridges", "europepmc", "ggpubr", "ggrepel", "ggsci", "ggthemes",
                         "ggExtra", "ggforce", "ggalluvial", "lattice", "latticeExtra", "BiocManager",
                         "org.Mm.eg.db", "ggplotify", "svglite", "tidyr", "dplyr", "pheatmap", "proxy",
                         "tibble", "openxlsx", "future", "future.apply")
  bioc_packages <- c("clusterProfiler", "pathview", "enrichplot", "DOSE", "org.Mm.eg.db")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}

setupPackages()

# ----------------------------------------------------
# 2. DIRECTORY ORGANIZATION FUNCTIONS
# ----------------------------------------------------
create_analysis_dirs <- function(base_dir, comparison_name, ontology) {
  dirs <- list(
    results = file.path(base_dir, "Results", comparison_name),
    go_ont = file.path(base_dir, "Results", comparison_name, "GO", ontology),
    kegg = file.path(base_dir, "Results", comparison_name, "KEGG"),
    ora = file.path(base_dir, "Results", comparison_name, "ORA"),
    custom = file.path(base_dir, "Results", comparison_name, "Custom"),
    pathview = file.path(base_dir, "Results", comparison_name, "pathview"),
    plots_go = file.path(base_dir, "Plots", comparison_name, paste0("GO_", ontology)),
    plots_kegg = file.path(base_dir, "Plots", comparison_name, "KEGG"),
    plots_ora = file.path(base_dir, "Plots", comparison_name, "ORA"),
    plots_custom = file.path(base_dir, "Plots", comparison_name, "Custom"),
    core_enrich = file.path(base_dir, "Datasets/core_enrichment", ontology)
  )
  lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))
  return(dirs)
}

save_plot_organized <- function(plot, filename, directory) {
  ggsave(file.path(directory, filename), plot, units = "cm", dpi = 300)
}

read_gct <- function(file_path) {
  gct_data <- read.delim(file_path, skip = 2, header = FALSE, check.names = FALSE)
  if ("Name" %in% colnames(gct_data)) {
    colnames(gct_data)[1] <- "Gene"
  }
  return(gct_data)
}

# ----------------------------------------------------
# 3. DATA INPUT/COMPARISONS
# ----------------------------------------------------
mapped_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/mapped"
comparison_files <- list.files(mapped_dir, pattern = "\\.csv$", full.names = FALSE)
comparison_list <- lapply(comparison_files, function(f) {
  parts <- strsplit(sub("\\.csv$", "", f), "_")[[1]]
  if (length(parts) == 2 && all(nzchar(parts))) parts else NULL
})
comparison_list <- Filter(Negate(is.null), comparison_list)

if (length(comparison_list) == 0) {
  warning("No valid comparison files found in: ", mapped_dir)
}

working_base <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
working_dir <- working_base
mapped_data_base <- file.path(working_base, "Datasets/mapped")
organism <- "org.Mm.eg.db"
ont <- "BP"

nk3r_genes <- c("P21279", "P21278", "P51432", "P11881", "P63318", "P68404", "P0DP26", "P0DP27", "P0DP28", "P11798", "P28652", "P47937", "P47713")
selected_uniprot <- c("P21279", "P21278", "Q9Z1B3", "P51432", "P11881", "P68404", "P63318", "P0DP26", "P0DP27", "P11798", "P28652", "Q61411", "Q99N57", "P31938", "P63085", "Q63844", "Q8BWG8", "Q91YI4", "V9GXQ9")
path_ids <- c("mmu04110", "mmu04115", "mmu04114", "mmu04113", "mmu04112", "mmu04111", "mmu04116", "mmu04117", "mmu04118", "mmu04119", "mmu04720", "mmu04721", "mmu04722", "mmu04725", "mmu04726", "mmu04727", "mmu04724", "mmu04080", "mmu00030", "mmu04151")

# ----------------------------------------------------
# 4. SETUP PARALLEL PROCESSING WITH FUTURE
# ----------------------------------------------------
library(future)
library(future.apply)

# Set up multisession plan
n_cores <- availableCores() - 1
plan(multisession, workers = n_cores)

cat("Setting up parallel processing with", n_cores, "cores using future\n")

# ----------------------------------------------------
# 5. ANALYSIS FUNCTION
# ----------------------------------------------------
analyze_comparison <- function(cell_types, working_base, mapped_data_base, organism, ont, 
                               nk3r_genes, selected_uniprot, path_ids) {
  
  # Load required libraries in worker
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(pathview)
    library(enrichplot)
    library(DOSE)
    library(ggplot2)
    library(tidyr)
    library(dplyr)
    library(openxlsx)
    library(org.Mm.eg.db)
  })
  
  comparison_name <- paste(cell_types, collapse = "_")
  
  tryCatch({
    # Setup directories
    dirs <- create_analysis_dirs(working_base, comparison_name, ont)
    
    # Define data file path
    file_name <- paste0(comparison_name, ".csv")
    data_path <- file.path(mapped_data_base, file_name)
    
    if (!file.exists(data_path)) {
      return(list(status = "FAILED", comparison = comparison_name, error = "File not found"))
    }
    
    # Load and prepare data
    df <- read.csv(data_path, header = TRUE, stringsAsFactors = FALSE)
    if (!"log2fc" %in% colnames(df)) {
      if ("logFC" %in% colnames(df)) {
        colnames(df)[colnames(df) == "logFC"] <- "log2fc"
      } else {
        return(list(status = "FAILED", comparison = comparison_name, error = "No log2fc column"))
      }
    }
    
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
    # GSEA (GO)
    # ----------------------------------------------------
    gse <- gseGO(geneList = gene_list, ont = ont, keyType = "UNIPROT", 
                 minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, 
                 verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH")
    
    gse_plot <- clusterProfiler::dotplot(gse, showCategory = 10, split = ".sign") +
      facet_wrap(~ .sign, nrow = 1) +
      labs(title = paste("GSEA", ont, "of", comparison_name), x = "Gene Ratio", y = "Gene Set") +
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
    
    save_plot_organized(gse_plot, paste0("GSEA_", ont, "_dotplot.svg"), dirs$plots_go)
    write.csv(gse@result, file = file.path(dirs$core_enrich, paste0(comparison_name, ".csv")), row.names = FALSE)
    write.csv(gse@result, file = file.path(dirs$go_ont, paste0("GSEA_", ont, "_results.csv")), row.names = FALSE)
    
    # GSEA Visualization Plots
    save_plot_organized(emapplot(pairwise_termsim(gse), showCategory = 10), paste0("GSEA_", ont, "_emap.svg"), dirs$plots_go)
    save_plot_organized(cnetplot(gse, categorySize = "pvalue", foldChange = gene_list), paste0("GSEA_", ont, "_cnet.svg"), dirs$plots_go)
    
    ridgeplot_gse <- ridgeplot(gse) + labs(x = "Enrichment Distribution", title = paste("GSEA", ont, "Ridgeplot")) + theme_minimal()
    save_plot_organized(ridgeplot_gse, paste0("GSEA_", ont, "_ridge.svg"), dirs$plots_go)
    
    gseaplot_gse <- gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1)
    save_plot_organized(gseaplot_gse, paste0("GSEA_", ont, "_plot.svg"), dirs$plots_go)
    
    top_terms <- head(gse@result$Description, 3)
    pmcplot_gse <- pmcplot(top_terms, 2010:2025, proportion = FALSE) + labs(title = paste("Publication Trends -", ont))
    save_plot_organized(pmcplot_gse, paste0("GSEA_", ont, "_pubmed.svg"), dirs$plots_go)
    
    # ORA
    ora <- enrichGO(gene = top_genes, ont = ont, keyType = "UNIPROT", 
                    minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, 
                    OrgDb = organism, pAdjustMethod = "BH")
    
    ora_plot <- clusterProfiler::dotplot(ora, showCategory = 10) + 
      labs(title = paste("ORA", ont, "- Top Regulated Genes")) + 
      scale_x_continuous(limits = c(0, 1))
    save_plot_organized(ora_plot, paste0("ORA_", ont, "_dotplot.svg"), dirs$plots_ora)
    write.csv(ora@result, file = file.path(dirs$ora, paste0("ORA_", ont, "_results.csv")), row.names = FALSE)
    
    # KEGG GSEA
    ids <- bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = organism)
    dedup_ids <- ids[!duplicated(ids$UNIPROT), ]
    df2 <- merge(df, dedup_ids, by.x = "gene_symbol", by.y = "UNIPROT")
    
    kegg_gene_list <- df2$log2fc
    names(kegg_gene_list) <- df2$ENTREZID
    kegg_gene_list <- kegg_gene_list[!duplicated(names(kegg_gene_list))]
    kegg_gene_list <- sort(na.omit(kegg_gene_list), decreasing = TRUE)
    
    kk2 <- gseKEGG(geneList = kegg_gene_list, organism = "mmu", 
                   minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, 
                   pAdjustMethod = "BH", keyType = "ncbi-geneid", verbose = FALSE)
    
    kegg_dot <- clusterProfiler::dotplot(kk2, showCategory = 10, split = ".sign") + 
      facet_wrap(~ .sign, nrow = 1) + 
      labs(title = "KEGG GSEA") + 
      theme_minimal()
    save_plot_organized(kegg_dot, "KEGG_dotplot.svg", dirs$plots_kegg)
    save_plot_organized(emapplot(pairwise_termsim(kk2), showCategory = 10), "KEGG_emap.svg", dirs$plots_kegg)
    save_plot_organized(cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list), "KEGG_cnet.svg", dirs$plots_kegg)
    save_plot_organized(ridgeplot(kk2), "KEGG_ridge.svg", dirs$plots_kegg)
    save_plot_organized(gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1), "KEGG_plot.svg", dirs$plots_kegg)
    write.csv(kk2@result, file = file.path(dirs$kegg, "KEGG_GSEA_results.csv"), row.names = FALSE)
    
    # Pathview
    pathview_dir <- normalizePath(dirs$pathview, winslash = "/", mustWork = FALSE)
    oldwd <- getwd()
    setwd(pathview_dir)
    
    lapply(path_ids, function(pid) {
      pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu", 
               low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg")
    })
    
    setwd(oldwd)
    
    # EnrichGO Analysis
    go_enrich <- enrichGO(gene = names(gene_list), universe = names(gene_list), 
                          OrgDb = organism, keyType = 'UNIPROT', readable = TRUE, 
                          ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, 
                          pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 800)
    
    p4 <- clusterProfiler::dotplot(go_enrich, showCategory = 20, split = "ONTOLOGY") +
      facet_grid(ONTOLOGY ~ ., scales = "free_y") +
      labs(title = "GO Enrichment Dotplot", x = "Gene Ratio", y = "GO Term", color = "p.adjust", size = "Count") +
      scale_color_viridis_c(option = "magma", direction = -1) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    save_plot_organized(p4, paste0("GOenrich_dotplot.svg"), dirs$plots_go)
    write.csv(go_enrich@result, file = file.path(dirs$go_ont, paste0("enrichGO_ALL_results.csv")), row.names = FALSE)
    
    # KEGG GSEA with Predefined UniProt IDs
    selected_entrez <- bitr(selected_uniprot, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = organism)
    selected_df <- merge(df, selected_entrez, by.x = "gene_symbol", by.y = "UNIPROT")
    selected_kegg_list <- selected_df$log2fc
    names(selected_kegg_list) <- selected_df$ENTREZID
    selected_kegg_list <- selected_kegg_list[!duplicated(names(selected_kegg_list))]
    selected_kegg_list <- sort(na.omit(selected_kegg_list), decreasing = TRUE)
    
    gsea_kegg_selected <- gseKEGG(geneList = selected_kegg_list, organism = "mmu", 
                                   minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, 
                                   pAdjustMethod = "BH", keyType = "ncbi-geneid", verbose = FALSE)
    
    kegg_selected_dot <- clusterProfiler::dotplot(gsea_kegg_selected, showCategory = 10, split = ".sign") + 
      facet_wrap(~ .sign, nrow = 1) + 
      labs(title = "KEGG GSEA (Predefined)") + 
      theme_minimal()
    save_plot_organized(kegg_selected_dot, "KEGG_Predefined_dotplot.svg", dirs$plots_kegg)
    write.csv(gsea_kegg_selected@result, file = file.path(dirs$kegg, "KEGG_GSEA_Predefined_UniProt.csv"), row.names = FALSE)
    
    save_plot_organized(emapplot(pairwise_termsim(gsea_kegg_selected), showCategory = 10), "KEGG_Predefined_emap.svg", dirs$plots_kegg)
    save_plot_organized(cnetplot(gsea_kegg_selected, categorySize = "pvalue", foldChange = selected_kegg_list), "KEGG_Predefined_cnet.svg", dirs$plots_kegg)
    save_plot_organized(ridgeplot(gsea_kegg_selected), "KEGG_Predefined_ridge.svg", dirs$plots_kegg)
    save_plot_organized(gseaplot(gsea_kegg_selected, by = "all", title = gsea_kegg_selected$Description[1], geneSetID = 1), "KEGG_Predefined_plot.svg", dirs$plots_kegg)
    
    # Custom GSEA: NK3R-signalling
    term2gene_nk3r <- data.frame(term = rep("NK3R-signalling", length(nk3r_genes)), gene = nk3r_genes)
    custom_gene_list <- df$log2fc
    names(custom_gene_list) <- df$gene_symbol
    custom_gene_list <- sort(na.omit(custom_gene_list), decreasing = TRUE)
    custom_gene_list <- custom_gene_list[!duplicated(names(custom_gene_list))]
    
    gsea_nk3r <- clusterProfiler::GSEA(geneList = custom_gene_list, TERM2GENE = term2gene_nk3r, 
                                       pvalueCutoff = 1, minGSSize = 1, maxGSSize = 500, verbose = FALSE)
    
    nk3r_dot <- clusterProfiler::dotplot(gsea_nk3r, showCategory = 10, split = ".sign") + 
      facet_wrap(~ .sign, nrow = 1) + 
      labs(title = "NK3R-signalling GSEA") + 
      theme_minimal()
    save_plot_organized(nk3r_dot, "NK3R_dotplot.svg", dirs$plots_custom)
    
    nk3r_plot <- gseaplot(gsea_nk3r, by = "all", title = "NK3R-signalling", geneSetID = 1)
    save_plot_organized(nk3r_plot, "NK3R_gsea_plot.svg", dirs$plots_custom)
    openxlsx::write.xlsx(gsea_nk3r@result, file = file.path(dirs$custom, "NK3R_GSEA_results.xlsx"))
    
    return(list(status = "SUCCESS", comparison = comparison_name))
    
  }, error = function(e) {
    return(list(status = "ERROR", comparison = comparison_name, error = conditionMessage(e)))
  })
}

# ----------------------------------------------------
# 6. RUN PARALLEL ANALYSIS WITH FUTURE
# ----------------------------------------------------
cat("\n==============================================\n")
cat("STARTING PARALLEL GSEA ANALYSIS\n")
cat("==============================================\n\n")

results <- future_lapply(comparison_list, function(cell_types) {
  analyze_comparison(cell_types, working_base, mapped_data_base, organism, ont, 
                    nk3r_genes, selected_uniprot, path_ids)
}, future.seed = TRUE)

# Reset to sequential processing
plan(sequential)

# Print summary
cat("\n==============================================\n")
cat("PARALLEL ANALYSIS SUMMARY\n")
cat("==============================================\n\n")

for (result in results) {
  if (result$status == "SUCCESS") {
    cat("✓", result$comparison, "- COMPLETED\n")
  } else {
    cat("✗", result$comparison, "- FAILED:", result$error, "\n")
  }
}

cat("\n==============================================\n")
cat("ALL COMPARISONS COMPLETED!\n")
cat("==============================================\n\n")


# ----------------------------------------------------
# 5. CELLTYPE SCORING ANALYSIS (Sequential)
# ----------------------------------------------------
cat("\n==============================================\n")
cat("STARTING CELLTYPE SCORING ANALYSIS\n")
cat("==============================================\n\n")

grouping_factor <- "Celltype"

# Load celltypes markers
celltype_file <- file.path(working_base, "Datasets/celltypes_long.xlsx")
celltype_df <- read.xlsx(celltype_file, sheet = 1)

# Load UniProt mapping
uniprot_mapping_file_path <- file.path(working_dir, "Datasets", "MOUSE_10090_idmapping.dat")
cat("Loading UniProt mapping file from:", uniprot_mapping_file_path, "\n")
uniprot_df <- read.delim(uniprot_mapping_file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract Gene_Name mappings
gene_names <- uniprot_df %>%
  filter(V2 == "Gene_Name") %>%
  distinct(V1, .keep_all = TRUE) %>%
  rename(UniProtID = V1, Gene_Name = V3) %>%
  select(UniProtID, Gene_Name)

# Extract UniProtKB-ID mappings
uniprot_ids <- uniprot_df %>%
  filter(V2 == "UniProtKB-ID") %>%
  distinct(V1, .keep_all = TRUE) %>%
  rename(UniProtID = V1, UniProtKB_ID = V3) %>%
  select(UniProtID, UniProtKB_ID)

# Combine both mappings
gene_map <- gene_names %>%
  left_join(uniprot_ids, by = "UniProtID")

# Load GCT file
gct_file <- file.path(working_base, "Datasets/gct/data/pg.matrix_filtered_pcaAdjusted.gct")
gct_data <- read_gct(gct_file)

# Extract metadata
fetch_meta <- function(label) {
  idx <- which(gct_data[, 1] == label)
  if (length(idx) == 0) stop("Metadata label not found in GCT: ", label)
  as.character(unlist(gct_data[idx, -1]))
}

sample_id      <- fetch_meta("sampleNumber")
celltype       <- fetch_meta("celltype")
groups         <- fetch_meta("ExpGroup")
celltype_group <- fetch_meta("celltype_group")

stopifnot(length(sample_id) == length(celltype), length(groups) == length(celltype_group), length(sample_id) == length(groups))

metadata_df <- data.frame(
  ColName         = paste0("V", 2:(ncol(gct_data))),
  Sample          = sample_id,
  Celltype        = celltype,
  Group           = groups,
  Celltype_Group  = celltype_group,
  stringsAsFactors = FALSE
)

# Extract expression data
expr_data <- gct_data[-c(1:4), ]
colnames(expr_data)[1] <- "Protein_ID"
colnames(expr_data)[-1] <- sample_id

# Reshape to long format
expr_long <- expr_data %>%
  pivot_longer(cols = -Protein_ID, names_to = "Sample", values_to = "Expression")

# Join with metadata
expr_annotated <- expr_long %>%
  left_join(metadata_df %>% select(-ColName), by = "Sample") %>%
  mutate(Expression = as.numeric(Expression)) %>%
  mutate(Protein_ID = sub(";.*", "", Protein_ID))

# Join with gene map
expr_annotated <- expr_annotated %>%
  left_join(gene_map, by = c("Protein_ID" = "UniProtKB_ID"))

# Add final label column
expr_annotated <- expr_annotated %>%
  filter(!grepl("Background", Celltype)) %>%
  relocate(Gene_Name, .before = 1) %>%
  relocate(UniProtID, .after = Gene_Name) %>%
  relocate(Expression, .after = last_col())

# Reshape marker list
celltype_long <- celltype_df %>%
  pivot_longer(cols = everything(), names_to = "Celltype_Class", values_to = "Gene_Label") %>%
  filter(!is.na(Gene_Label) & Gene_Label != "")

# Filter for matching genes
marker_expr <- expr_annotated %>%
  filter(Gene_Name %in% celltype_long$Gene_Label)

# Add Celltype_Class info
marker_expr <- marker_expr %>%
  left_join(celltype_long, by = c("Gene_Name" = "Gene_Label"))

# Compute average expression
celltype_scores <- marker_expr %>%
  group_by(!!sym(grouping_factor), Celltype_Class) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

# Pivot wider for table format
celltype_score_table <- celltype_scores %>%
  pivot_wider(names_from = Celltype_Class, values_from = Mean_Expression)

# Create bar plot of celltype scores
celltype_score_plot <- ggplot(celltype_scores, aes(x = !!sym(grouping_factor), y = Mean_Expression, fill = Celltype_Class)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Celltype Scores by Group", x = grouping_factor, y = "Mean Expression", fill = "Cell Class") +
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

ggsave(file.path(working_dir, "celltype_scores_by_group.svg"), celltype_score_plot, units = "cm", dpi = 300)

# Create pheatmap of marker expression
gene_class_df <- celltype_long %>%
  distinct(Gene_Label, Celltype_Class) %>%
  rename(Gene_Name = Gene_Label)

# Keep track of which class each gene belongs to
gene_class_expanded <- celltype_long %>%
  rename(Gene_Name = Gene_Label) %>%
  group_by(Gene_Name) %>%
  mutate(
    Gene_Name_Unique = if(n() > 1) {
      paste0(Gene_Name, "_", Celltype_Class)  # Add class suffix for duplicates
    } else {
      Gene_Name
    }
  ) %>%
  ungroup()

annotation_row <- gene_class_expanded %>%
  distinct(Gene_Name_Unique, Celltype_Class) %>%
  column_to_rownames("Gene_Name_Unique")

# Calculate mean expression per gene and grouping_factor
marker_expr_unique <- marker_expr %>%
  left_join(
    gene_class_expanded %>% select(Gene_Name, Celltype_Class, Gene_Name_Unique),
    by = c("Gene_Name", "Celltype_Class")
  )

marker_matrix <- marker_expr_unique %>%
  group_by(Gene_Name_Unique, !!sym(grouping_factor)) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = !!sym(grouping_factor), values_from = Mean_Expression) %>%
  column_to_rownames("Gene_Name_Unique") %>%
  as.matrix()

# Rest of cleaning code stays the same
marker_matrix_clean <- marker_matrix[rowSums(is.na(marker_matrix)) < ncol(marker_matrix), ]
marker_matrix_clean <- marker_matrix_clean[, colSums(is.na(marker_matrix_clean)) < nrow(marker_matrix_clean)]
marker_matrix_clean <- marker_matrix_clean[rowSums(!is.na(marker_matrix_clean)) >= 2, ]
marker_matrix_clean <- marker_matrix_clean[, colSums(!is.na(marker_matrix_clean)) >= 2]
marker_matrix_clean <- marker_matrix_clean[!grepl("^Oasl2\\s*$", rownames(marker_matrix_clean)), ]
marker_matrix_clean[is.nan(marker_matrix_clean)] <- NA
marker_matrix_clean[is.infinite(marker_matrix_clean)] <- NA

# Filter annotation_row to match cleaned matrix
annotation_row <- annotation_row[rownames(marker_matrix_clean), , drop = FALSE]

# Color scheme for annotations
annotation_colors <- list(Celltype_Class = c("gaba" = "#f36d07", "vglut1" = "#455A64", "vglut2" = "#c7c7c7"))

# Save pheatmap
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

# Z-score heatmap
expr_scaled <- marker_expr %>%
  group_by(Gene_Name) %>%
  mutate(Z_Expression = as.numeric(scale(Expression))) %>%
  ungroup()

celltype_z_scores <- expr_scaled %>%
  group_by(!!sym(grouping_factor), Celltype_Class) %>%
  summarise(Mean_Z = mean(Z_Expression, na.rm = TRUE), .groups = "drop")

heatmap_plot <- ggplot(celltype_z_scores, aes(x = reorder(Celltype_Class, Mean_Z, FUN = median), y = .data[[grouping_factor]], fill = Mean_Z)) +
  geom_tile(color = "grey80", size = 0.2, width = 0.5) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "Z-Score", guide = guide_colorbar(barwidth = 0.8, barheight = 20)) +
  labs(title = "Celltype Marker Signature by Group", x = "Marker Class", y = "Sample Group") +
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

cat("\n==============================================\n")
cat("ENTIRE PIPELINE COMPLETED!\n")
cat("==============================================\n")
