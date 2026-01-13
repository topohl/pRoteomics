#' ============================================================
#' HIGH-THROUGHPUT Gene Set Enrichment Analysis (GSEA) Workflow
#' WITH PARALLEL PROCESSING AND GO TERM SIMPLIFICATION
#' ============================================================
#' 
#' This script performs comprehensive GSEA, KEGG, ORA, custom pathway analysis
#' for MULTIPLE cell type comparisons in parallel using future.
#' Includes GO term redundancy removal via simplify().
#' Celltype scoring runs separately at the end.
#'
#' @author Tobias Pohl
#' ============================================================

# ----------------------------------------------------
# 0. SIMPLIFICATION SETTINGS
# ----------------------------------------------------
# *** MAIN CONTROL: Set whether to perform simplification at all ***
PERFORM_SIMPLIFICATION <- FALSE  # Set to FALSE to skip simplification entirely

# Control which version to use for plots (only matters if PERFORM_SIMPLIFICATION = TRUE)
USE_SIMPLIFIED_FOR_PLOTS <- FALSE  # TRUE = simplified plots, FALSE = full results

# Semantic similarity cutoff for simplify() function (only used if simplification is enabled)
# Range: 0.4 (strict, removes more) to 0.9 (lenient, keeps more)
# Default: 0.7 (moderate redundancy removal)
SIMPLIFY_CUTOFF <- 0.7

# Note: 
# - If PERFORM_SIMPLIFICATION = FALSE, only full results are computed and saved
# - If PERFORM_SIMPLIFICATION = TRUE, both versions are computed and saved

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
                         "tibble", "openxlsx", "future", "future.apply", "GOSemSim")
  bioc_packages <- c("clusterProfiler", "pathview", "enrichplot", "DOSE", "org.Mm.eg.db", "GOSemSim")
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
ont <- "CC"

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
    library(GOSemSim)
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
# GSEA (GO) WITH OPTIONAL SIMPLIFICATION
# ----------------------------------------------------
gse <- gseGO(geneList = gene_list, ont = ont, keyType = "UNIPROT", 
             minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, 
             verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH")


# Perform simplification ONLY if requested
gse_simplified <- NULL  # Initialize as NULL

if (PERFORM_SIMPLIFICATION && !is.null(gse) && nrow(gse@result) > 0) {
  gse_temp <- tryCatch({
    simplify(gse, cutoff = SIMPLIFY_CUTOFF, by = "p.adjust", select_fun = min)
  }, error = function(e) {
    warning("Simplify failed for ", comparison_name, ": ", e$message)
    NULL  # Return NULL on failure
  })

  # Validate S4 class
  if (!is.null(gse_temp) && (is(gse_temp, "gseaResult") || is(gse_temp, "enrichResult"))) {
    gse_simplified <- gse_temp
  }
}


# Determine which version to use for plotting
if (PERFORM_SIMPLIFICATION && !is.null(gse_simplified)) {
  gse_for_plot <- if(USE_SIMPLIFIED_FOR_PLOTS) gse_simplified else gse
  plot_suffix <- if(USE_SIMPLIFIED_FOR_PLOTS) "_simplified" else "_full"
} else {
  # No simplification performed - use full results only
  gse_for_plot <- gse
  plot_suffix <- "_full"
}


# Generate plots using selected version
if (!is.null(gse_for_plot) && nrow(gse_for_plot@result) > 0) {
  plot_label <- if(PERFORM_SIMPLIFICATION && USE_SIMPLIFIED_FOR_PLOTS && !is.null(gse_simplified)) {
    "(Simplified)"
  } else {
    "(Full)"
  }

  gse_plot <- clusterProfiler::dotplot(gse_for_plot, showCategory = 10, split = ".sign") +
    facet_wrap(~ .sign, nrow = 1) +
    labs(title = paste("GSEA", ont, "of", comparison_name, plot_label), 
         x = "Gene Ratio", y = "Gene Set") +
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

  save_plot_organized(gse_plot, paste0("GSEA_", ont, "_dotplot", plot_suffix, ".svg"), dirs$plots_go)

  # Other plots
  save_plot_organized(emapplot(pairwise_termsim(gse_for_plot), showCategory = 10), 
                     paste0("GSEA_", ont, "_emap", plot_suffix, ".svg"), dirs$plots_go)
  save_plot_organized(cnetplot(gse_for_plot, categorySize = "pvalue", foldChange = gene_list), 
                     paste0("GSEA_", ont, "_cnet", plot_suffix, ".svg"), dirs$plots_go)

  ridgeplot_gse <- ridgeplot(gse_for_plot) + 
    labs(x = "Enrichment Distribution", 
         title = paste("GSEA", ont, "Ridgeplot", plot_label)) + 
    theme_minimal()
  save_plot_organized(ridgeplot_gse, paste0("GSEA_", ont, "_ridge", plot_suffix, ".svg"), dirs$plots_go)

  gseaplot_gse <- gseaplot(gse_for_plot, by = "all", 
                           title = gse_for_plot@result$Description[1], geneSetID = 1)
  save_plot_organized(gseaplot_gse, paste0("GSEA_", ont, "_plot", plot_suffix, ".svg"), dirs$plots_go)

  top_terms <- head(gse_for_plot@result$Description, 3)
  pmcplot_gse <- pmcplot(top_terms, 2010:2025, proportion = FALSE) + 
    labs(title = paste("Publication Trends -", ont, plot_label))
  save_plot_organized(pmcplot_gse, paste0("GSEA_", ont, "_pubmed", plot_suffix, ".svg"), dirs$plots_go)
}


# Save results based on whether simplification was performed
write.csv(gse@result, 
          file = file.path(dirs$go_ont, paste0("GSEA_", ont, "_results_full.csv")), 
          row.names = FALSE)

if (PERFORM_SIMPLIFICATION && !is.null(gse_simplified)) {
  # Save simplified version
  write.csv(gse_simplified@result, 
            file = file.path(dirs$go_ont, paste0("GSEA_", ont, "_results_simplified.csv")), 
            row.names = FALSE)
}

# Save the version used for plots to core_enrich
gse_for_export <- if(PERFORM_SIMPLIFICATION && USE_SIMPLIFIED_FOR_PLOTS && !is.null(gse_simplified)) {
  gse_simplified
} else {
  gse
}

write.csv(gse_for_export@result, 
          file = file.path(dirs$core_enrich, paste0(comparison_name, plot_suffix, ".csv")), 
          row.names = FALSE)

# ----------------------------------------------------
# ORA WITH OPTIONAL SIMPLIFICATION
# ----------------------------------------------------
ora <- enrichGO(gene = top_genes, ont = ont, keyType = "UNIPROT", 
                minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, 
                OrgDb = organism, pAdjustMethod = "BH")


# Simplify ONLY if requested
ora_simplified <- NULL

if (PERFORM_SIMPLIFICATION && nrow(ora@result) > 0) {
  ora_simplified <- tryCatch({
    simplify(ora, cutoff = SIMPLIFY_CUTOFF, by = "p.adjust", select_fun = min)
  }, error = function(e) {
    warning("ORA simplify failed for ", comparison_name, ": ", e$message)
    NULL
  })

  # Validate S4 class
  if (!is.null(ora_simplified) && !is(ora_simplified, "enrichResult")) {
    ora_simplified <- NULL
  }
}


# Determine version for plotting
if (PERFORM_SIMPLIFICATION && !is.null(ora_simplified)) {
  ora_for_plot <- if(USE_SIMPLIFIED_FOR_PLOTS) ora_simplified else ora
} else {
  ora_for_plot <- ora
}


if (nrow(ora_for_plot@result) > 0) {
  plot_label <- if(PERFORM_SIMPLIFICATION && USE_SIMPLIFIED_FOR_PLOTS && !is.null(ora_simplified)) {
    "(Simplified)"
  } else {
    "(Full)"
  }

  ora_plot <- clusterProfiler::dotplot(ora_for_plot, showCategory = 10) + 
    labs(title = paste("ORA", ont, "- Top Regulated Genes", plot_label)) + 
    scale_x_continuous(limits = c(0, 1))
  save_plot_organized(ora_plot, paste0("ORA_", ont, "_dotplot", plot_suffix, ".svg"), dirs$plots_ora)
}


# Save results
write.csv(ora@result, 
          file = file.path(dirs$ora, paste0("ORA_", ont, "_results_full.csv")), 
          row.names = FALSE)

if (PERFORM_SIMPLIFICATION && !is.null(ora_simplified)) {
  write.csv(ora_simplified@result, 
            file = file.path(dirs$ora, paste0("ORA_", ont, "_results_simplified.csv")), 
            row.names = FALSE)
}
    
    # ----------------------------------------------------
    # KEGG GSEA (No simplification - KEGG specific)
    # ----------------------------------------------------
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
    
    if (nrow(kk2@result) > 0) {
      save_plot_organized(gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1), "KEGG_plot.svg", dirs$plots_kegg)
    }
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
    
# ----------------------------------------------------
# EnrichGO Analysis (ALL ontologies) WITH OPTIONAL SIMPLIFICATION
# ----------------------------------------------------
go_enrich <- enrichGO(gene = names(gene_list), universe = names(gene_list), 
                      OrgDb = organism, keyType = 'UNIPROT', readable = TRUE, 
                      ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1, 
                      pAdjustMethod = "BH", minGSSize = 3, maxGSSize = 800)


# Initialize simplified version as NULL
go_enrich_simplified <- NULL


# Process simplification ONLY if requested
if (PERFORM_SIMPLIFICATION) {
  simplified_results_list <- list()

  for (ont_type in c("BP", "CC", "MF")) {
    # Run enrichGO for this specific ontology
    go_single <- tryCatch({
      enrichGO(gene = names(gene_list), 
               universe = names(gene_list), 
               OrgDb = organism, 
               keyType = 'UNIPROT', 
               readable = TRUE, 
               ont = ont_type,
               pvalueCutoff = 1, 
               qvalueCutoff = 1, 
               pAdjustMethod = "BH", 
               minGSSize = 3, 
               maxGSSize = 800)
    }, error = function(e) {
      warning("enrichGO failed for ", ont_type, ": ", e$message)
      NULL
    })

    # Simplify if results exist
    if (!is.null(go_single) && nrow(go_single@result) > 0) {
      go_simplified <- tryCatch({
        simplify(go_single, cutoff = SIMPLIFY_CUTOFF, by = "p.adjust", select_fun = min)
      }, error = function(e) {
        warning("Simplify failed for ", ont_type, ": ", e$message)
        go_single
      })

      # Store results as data.frame
      if (is(go_simplified, "enrichResult")) {
        simplified_results_list[[ont_type]] <- go_simplified@result
      } else if (is.data.frame(go_simplified)) {
        simplified_results_list[[ont_type]] <- go_simplified
      } else {
        simplified_results_list[[ont_type]] <- go_single@result
      }
    }
  }

  # Combine simplified results
  if (length(simplified_results_list) > 0) {
    go_enrich_simplified <- go_enrich
    go_enrich_simplified@result <- do.call(rbind, simplified_results_list)
  }
}


# Choose version for plotting
if (PERFORM_SIMPLIFICATION && !is.null(go_enrich_simplified)) {
  go_enrich_for_plot <- if(USE_SIMPLIFIED_FOR_PLOTS) go_enrich_simplified else go_enrich
} else {
  go_enrich_for_plot <- go_enrich
}


# Plot with selected version
if (nrow(go_enrich_for_plot@result) > 0) {
  plot_label <- if(PERFORM_SIMPLIFICATION && USE_SIMPLIFIED_FOR_PLOTS && !is.null(go_enrich_simplified)) {
    "(Simplified)"
  } else {
    "(Full)"
  }

  p4 <- clusterProfiler::dotplot(go_enrich_for_plot, showCategory = 20, split = "ONTOLOGY") +
    facet_grid(ONTOLOGY ~ ., scales = "free_y") +
    labs(title = paste("GO Enrichment Dotplot", plot_label), 
         x = "Gene Ratio", y = "GO Term", color = "p.adjust", size = "Count") +
    scale_color_viridis_c(option = "magma", direction = -1) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot_organized(p4, paste0("GOenrich_dotplot", plot_suffix, ".svg"), dirs$plots_go)
}


# Save results
write.csv(go_enrich@result, 
          file = file.path(dirs$go_ont, paste0("enrichGO_ALL_results_full.csv")), 
          row.names = FALSE)

if (PERFORM_SIMPLIFICATION && !is.null(go_enrich_simplified)) {
  write.csv(go_enrich_simplified@result, 
            file = file.path(dirs$go_ont, paste0("enrichGO_ALL_results_simplified.csv")), 
            row.names = FALSE)
}
    
    # ----------------------------------------------------
    # KEGG GSEA with Predefined UniProt IDs
    # ----------------------------------------------------
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
    
    if (nrow(gsea_kegg_selected@result) > 0) {
      save_plot_organized(gseaplot(gsea_kegg_selected, by = "all", title = gsea_kegg_selected$Description[1], geneSetID = 1), "KEGG_Predefined_plot.svg", dirs$plots_kegg)
    }
    
    # ----------------------------------------------------
    # Custom GSEA: NK3R-signalling (No simplification needed)
    # ----------------------------------------------------
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
    
    if (nrow(gsea_nk3r@result) > 0) {
      nk3r_plot <- gseaplot(gsea_nk3r, by = "all", title = "NK3R-signalling", geneSetID = 1)
      save_plot_organized(nk3r_plot, "NK3R_gsea_plot.svg", dirs$plots_custom)
    }
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
# 7. CELLTYPE SCORING ANALYSIS (Sequential)
# ----------------------------------------------------
cat("\n==============================================\n")
cat("STARTING CELLTYPE SCORING ANALYSIS\n")
cat("==============================================\n\n")

# Define and create output directories for Celltype Scoring
celltype_results_base <- file.path(working_base, "Results", "Celltype_Scoring")
celltype_dirs <- list(
  base = celltype_results_base,
  plots = file.path(celltype_results_base, "Plots"),
  tables = file.path(celltype_results_base, "Tables"),
  heatmaps = file.path(celltype_results_base, "Heatmaps")
)

lapply(celltype_dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE))

grouping_factor <- "Celltype"

# Load celltypes markers
celltype_file <- file.path(working_base, "Datasets/celltypes_long.xlsx")
celltype_df <- read.xlsx(celltype_file, sheet = 1)

# Load UniProt mapping
uniprot_mapping_file_path <- file.path(working_dir, "Datasets", "MOUSE_10090_idmapping.dat")
cat("Loading UniProt mapping file from:", uniprot_mapping_file_path, "\n")
if (!file.exists(uniprot_mapping_file_path)) {
  stop("UniProt mapping file not found at: ", uniprot_mapping_file_path, "\nPlease verify the file path.")
}
uniprot_df <- read.delim(uniprot_mapping_file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract Gene_Name mappings
gene_names <- uniprot_df %>%
  dplyr::filter(V2 == "Gene_Name") %>%
  dplyr::distinct(V1, .keep_all = TRUE) %>%
  dplyr::rename(UniProtID = V1, Gene_Name = V3) %>%
  dplyr::select(UniProtID, Gene_Name)

# Extract UniProtKB-ID mappings - FIXED
uniprot_ids <- uniprot_df %>%
  dplyr::filter(V2 == "UniProtKB-ID") %>%
  dplyr::distinct(V1, .keep_all = TRUE) %>%
  dplyr::rename(UniProtID = V1, UniProtKB_ID = V3) %>%
  dplyr::select(UniProtID, UniProtKB_ID)

# Combine both mappings
gene_map <- gene_names %>%
  left_join(uniprot_ids, by = "UniProtID")

# Load GCT file
gct_file <- file.path(working_base, "Datasets/gct/data/pg.matrix_filtered_pcaAdjusted_unnormalized.gct")
if (!file.exists(gct_file)) {
  stop("GCT file not found at: ", gct_file, "\nPlease verify the file path.")
}
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
  left_join(metadata_df %>% dplyr::select(-ColName), by = "Sample") %>%
  mutate(Expression = suppressWarnings(as.numeric(Expression))) %>%
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

# ----------------------------------------------------
# PREPARE ANNOTATIONS (Static)
# ----------------------------------------------------
# Keep track of which class each gene belongs to
gene_class_expanded <- celltype_long %>%
  dplyr::rename(Gene_Name = Gene_Label) %>%
  group_by(Gene_Name) %>%
  mutate(
    Gene_Name_Unique = Gene_Name,  # Version without suffix
    Gene_Name_WithClass = if(n() > 1) {
      paste0(Gene_Name, "_", Celltype_Class)  # Add class suffix for duplicates
    } else {
      Gene_Name
    }
  ) %>%
  ungroup()

# Create multi-class annotation for genes in multiple classes
annotation_row_unique_base <- gene_class_expanded %>%
  group_by(Gene_Name_Unique) %>%
  summarise(Celltype_Class = paste(sort(unique(Celltype_Class)), collapse = "+"), .groups = "drop") %>%
  column_to_rownames("Gene_Name_Unique")

# Expanded color scheme for multi-class genes
unique_classes <- unique(annotation_row_unique_base$Celltype_Class)
base_colors <- c("gaba" = "#f36d07", "vglut1" = "#455A64", "vglut2" = "#c7c7c7")
extra_colors <- c("#00BCD4", "#8BC34A", "#FFC107", "#FF5722")

annotation_colors <- list(
  Celltype_Class = c(
    base_colors,
    setNames(extra_colors[1:(length(unique_classes) - length(base_colors))], 
             setdiff(unique_classes, names(base_colors)))
  )
)

# Create annotation for withclass version
annotation_row_withclass_base <- gene_class_expanded %>%
  distinct(Gene_Name_WithClass, Celltype_Class) %>%
  column_to_rownames("Gene_Name_WithClass")

# Color scheme for single-class annotations
annotation_colors_withclass <- list(Celltype_Class = c("gaba" = "#f36d07", "vglut1" = "#455A64", "vglut2" = "#c7c7c7"))

# ----------------------------------------------------
# DEFINE ANALYSIS FUNCTION
# ----------------------------------------------------
run_celltype_analysis <- function(current_data, suffix) {
  cat("  Generating plots for:", suffix, "\n")
  
  # Compute average expression
  celltype_scores <- current_data %>%
    group_by(!!sym(grouping_factor), Celltype_Class) %>%
    summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")
  
  # Save scores to Excel
  write.xlsx(celltype_scores, file = file.path(celltype_dirs$tables, paste0("celltype_scores", suffix, ".xlsx")))

  # Create bar plot of celltype scores
  celltype_score_plot <- ggplot(celltype_scores, aes(x = !!sym(grouping_factor), y = Mean_Expression, fill = Celltype_Class)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    labs(title = paste("Celltype Scores", suffix), x = grouping_factor, y = "Mean Expression", fill = "Cell Class") +
    scale_fill_manual(values = c("gaba" = "#f36d07", "vglut1" = "#455A64", "vglut2" = "#c7c7c7")) +
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
  
  ggsave(file.path(celltype_dirs$plots, paste0("celltype_scores", suffix, ".svg")), celltype_score_plot, units = "cm", dpi = 300)
  
  # Calculate mean expression per gene and grouping_factor
  marker_expr_unique <- current_data %>%
    left_join(
      gene_class_expanded %>% dplyr::select(Gene_Name, Celltype_Class, Gene_Name_Unique, Gene_Name_WithClass),
      by = c("Gene_Name", "Celltype_Class")
    )
  
  # ===== HEATMAP 1: Gene_Name_Unique (no suffix) =====
  marker_matrix_unique <- marker_expr_unique %>%
    group_by(Gene_Name_Unique, !!sym(grouping_factor)) %>%
    summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(grouping_factor), values_from = Mean_Expression) %>%
    column_to_rownames("Gene_Name_Unique") %>%
    as.matrix()
  
  # Clean matrix (unique version)
  marker_matrix_unique_clean <- marker_matrix_unique[rowSums(is.na(marker_matrix_unique)) < ncol(marker_matrix_unique), , drop=FALSE]
  if(ncol(marker_matrix_unique_clean) > 0 && nrow(marker_matrix_unique_clean) > 0) {
    marker_matrix_unique_clean <- marker_matrix_unique_clean[, colSums(is.na(marker_matrix_unique_clean)) < nrow(marker_matrix_unique_clean), drop=FALSE]
    # Ensure we have enough data for heatmap
    if(nrow(marker_matrix_unique_clean) >= 2 && ncol(marker_matrix_unique_clean) >= 2) {
        marker_matrix_unique_clean <- marker_matrix_unique_clean[rowSums(!is.na(marker_matrix_unique_clean)) >= 1, , drop=FALSE]
        marker_matrix_unique_clean <- marker_matrix_unique_clean[!grepl("^Oasl2\\s*$", rownames(marker_matrix_unique_clean)), , drop=FALSE]
        marker_matrix_unique_clean[is.nan(marker_matrix_unique_clean)] <- NA
        marker_matrix_unique_clean[is.infinite(marker_matrix_unique_clean)] <- NA
        
        # Filter annotation_row to match cleaned matrix
        annotation_row_unique_clean <- annotation_row_unique_base[rownames(marker_matrix_unique_clean), , drop = FALSE]
        
        # Save pheatmap (unique version)
        pheatmap(
          marker_matrix_unique_clean,
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          color = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
          na_col = "grey",
          main = paste("Marker Expression", suffix, "(Unique Names)"),
          fontsize = 12,
          fontsize_row = 10,
          fontsize_col = 12,
          cellheight = 12,
          cellwidth = 12,
          border_color = NA,
          annotation_row = annotation_row_unique_clean,
          annotation_colors = annotation_colors,
          filename = file.path(celltype_dirs$heatmaps, paste0("marker_expression_heatmap", suffix, "_unique.pdf"))
        )
        
        # ===== HEATMAP 1b: Sorted by Celltype_Class =====
        annotation_row_sorted <- annotation_row_unique_clean %>% arrange(Celltype_Class)
        marker_matrix_sorted <- marker_matrix_unique_clean[rownames(annotation_row_sorted), , drop = FALSE]
        
        pheatmap(
          marker_matrix_sorted,
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          color = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
          na_col = "grey",
          main = paste("Marker Expression", suffix, "(Sorted by Class)"),
          fontsize = 12,
          fontsize_row = 10,
          fontsize_col = 12,
          cellheight = 12,
          cellwidth = 12,
          border_color = NA,
          annotation_row = annotation_row_sorted,
          annotation_colors = annotation_colors,
          filename = file.path(celltype_dirs$heatmaps, paste0("marker_expression_heatmap", suffix, "_sorted.pdf"))
        )
    }
  }

  # ===== HEATMAP 2: Gene_Name_WithClass (with suffix) =====
  marker_matrix_withclass <- marker_expr_unique %>%
    group_by(Gene_Name_WithClass, !!sym(grouping_factor)) %>%
    summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(grouping_factor), values_from = Mean_Expression) %>%
    column_to_rownames("Gene_Name_WithClass") %>%
    as.matrix()
  
  # Clean matrix (withclass version)
  marker_matrix_withclass_clean <- marker_matrix_withclass[rowSums(is.na(marker_matrix_withclass)) < ncol(marker_matrix_withclass), , drop=FALSE]
  if(ncol(marker_matrix_withclass_clean) > 0 && nrow(marker_matrix_withclass_clean) > 0) {
      marker_matrix_withclass_clean <- marker_matrix_withclass_clean[, colSums(is.na(marker_matrix_withclass_clean)) < nrow(marker_matrix_withclass_clean), drop=FALSE]
      
      if(nrow(marker_matrix_withclass_clean) >= 2 && ncol(marker_matrix_withclass_clean) >= 2) {
          marker_matrix_withclass_clean <- marker_matrix_withclass_clean[rowSums(!is.na(marker_matrix_withclass_clean)) >= 1, , drop=FALSE]
          marker_matrix_withclass_clean <- marker_matrix_withclass_clean[!grepl("^Oasl2\\s*$", rownames(marker_matrix_withclass_clean)), , drop=FALSE]
          marker_matrix_withclass_clean[is.nan(marker_matrix_withclass_clean)] <- NA
          marker_matrix_withclass_clean[is.infinite(marker_matrix_withclass_clean)] <- NA
          
          # Filter annotation_row to match cleaned matrix
          annotation_row_withclass <- annotation_row_withclass_base[rownames(marker_matrix_withclass_clean), , drop = FALSE]
          
          # Save pheatmap (withclass version)
          pheatmap(
            marker_matrix_withclass_clean,
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            color = colorRampPalette(c("#6698CC", "white", "#F08C21"))(100),
            na_col = "grey",
            main = paste("Marker Expression", suffix, "(With Class Suffix)"),
            fontsize = 12,
            fontsize_row = 10,
            fontsize_col = 12,
            cellheight = 12,
            cellwidth = 12,
            border_color = NA,
            annotation_row = annotation_row_withclass,
            annotation_colors = annotation_colors_withclass,
            filename = file.path(celltype_dirs$heatmaps, paste0("marker_expression_heatmap", suffix, "_withclass.pdf"))
          )
      }
  }
  
  # Z-score heatmap
  expr_scaled <- current_data %>%
    group_by(Gene_Name) %>%
    mutate(Z_Expression = as.numeric(scale(Expression))) %>%
    ungroup()
  
  celltype_z_scores <- expr_scaled %>%
    group_by(!!sym(grouping_factor), Celltype_Class) %>%
    summarise(Mean_Z = mean(Z_Expression, na.rm = TRUE), .groups = "drop")
  
  heatmap_plot <- ggplot(celltype_z_scores, aes(x = reorder(Celltype_Class, Mean_Z, FUN = median), y = .data[[grouping_factor]], fill = Mean_Z)) +
    geom_tile(color = "grey80", size = 0.2, width = 0.5) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "Z-Score", guide = guide_colorbar(barwidth = 0.8, barheight = 20)) +
    labs(title = paste("Celltype Marker Signature", suffix), x = "Marker Class", y = "Sample Group") +
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
  
  ggsave(file.path(celltype_dirs$heatmaps, paste0("celltype_scores_heatmap", suffix, ".svg")), heatmap_plot, width = 16, height = 9, units = "cm", dpi = 300)
  
  # ===== NEW: Z-score heatmap sorted by Celltype_Class =====
  heatmap_plot_sorted <- ggplot(celltype_z_scores, aes(x = Celltype_Class, y = .data[[grouping_factor]], fill = Mean_Z)) +
    geom_tile(color = "grey80", size = 0.2, width = 0.5) +
    scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, name = "Z-Score", guide = guide_colorbar(barwidth = 0.8, barheight = 20)) +
    labs(title = paste("Celltype Marker Signature", suffix, "(Sorted)"), x = "Marker Class", y = "Sample Group") +
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
  
  ggsave(file.path(celltype_dirs$heatmaps, paste0("celltype_scores_heatmap_sorted", suffix, ".svg")), heatmap_plot_sorted, width = 16, height = 9, units = "cm", dpi = 300)
}

# ----------------------------------------------------
# EXECUTE ANALYSIS
# ----------------------------------------------------

# 1. Run on ALL data (Combined)
cat("Processing ALL data combined...\n")
run_celltype_analysis(marker_expr, "_All")

# 2. Run separately for each ExpGroup
unique_groups <- unique(marker_expr$Group)
cat("Processing separate groups:", paste(unique_groups, collapse=", "), "\n")

for (grp in unique_groups) {
  # Sanitize group name for filename
  safe_grp <- gsub("[^a-zA-Z0-9]", "_", grp)
  
  # Filter data for this group
  group_data <- marker_expr %>% filter(Group == grp)
  
  # Run analysis
  run_celltype_analysis(group_data, paste0("_", safe_grp))
}

# ----------------------------------------------------
# 8. COMPARISON BETWEEN EXP GROUPS
# ----------------------------------------------------
cat("\nProcessing Comparison between ExpGroups...\n")

# Calculate mean expression per Celltype_Class and Group
group_comparison <- marker_expr %>%
  group_by(Group, Celltype_Class) %>%
  summarise(Mean_Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

# Save comparison data to Excel
write.xlsx(group_comparison, file = file.path(celltype_dirs$tables, "celltype_scores_comparison_ExpGroups.xlsx"))

# Create comparison plot
comparison_plot <- ggplot(group_comparison, aes(x = Group, y = Mean_Expression, fill = Celltype_Class)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  labs(title = "Celltype Scores Comparison by ExpGroup", x = "Experimental Group", y = "Mean Expression", fill = "Cell Class") +
  scale_fill_manual(values = c("gaba" = "#f36d07", "vglut1" = "#455A64", "vglut2" = "#c7c7c7")) +
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
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 16)
  )

ggsave(file.path(celltype_dirs$plots, "celltype_scores_comparison_ExpGroups.svg"), comparison_plot, units = "cm", dpi = 300, width = 25, height = 15)

# ----------------------------------------------------
# 9. ADVANCED INTERACTION VISUALIZATION
# ----------------------------------------------------
cat("\nGenerating Advanced Interaction Plots...\n")

# 1. Balloon Plot (Dot Plot) - FANCY VERSION
# Visualizes Mean Expression (Size) and Z-Score (Color)
# Grouping by grouping_factor (Celltype), Group (ExpGroup), and Celltype_Class
interaction_stats <- marker_expr %>%
  group_by(!!sym(grouping_factor), Group, Celltype_Class) %>%
  summarise(
    Mean_Expr = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Celltype_Class) %>%
  mutate(
    Z_Score = scale(Mean_Expr),
    Rel_Expr = Mean_Expr / max(Mean_Expr)
  ) %>%
  ungroup()

balloon_plot <- ggplot(interaction_stats, aes(x = !!sym(grouping_factor), y = Celltype_Class)) +
  # Add a subtle background tile to define the grid
  geom_tile(fill = "white", color = "grey95", size = 0.5) +
  # The bubbles with a nicer stroke and alpha
  geom_point(aes(size = Mean_Expr, fill = Z_Score), shape = 21, color = "grey20", stroke = 1, alpha = 0.9) +
  # Facet by Experimental Group
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  # Improved color scale (High contrast diverging)
  scale_fill_gradient2(low = "#313695", mid = "white", high = "#a50026", midpoint = 0, 
                       name = "Z-Score\n(Row-scaled)") +
  # Adjusted size range
  scale_size_continuous(range = c(4, 14), name = "Mean\nExpression") +
  # Clean labels
  labs(
    title = paste("Interaction:", grouping_factor, "vs Marker Class by ExpGroup"),
    subtitle = "Balloon Plot: Size=Expression, Color=Z-score. Faceted by Experimental Group.",
    x = grouping_factor,
    y = "Marker Class"
  ) +
  # Polished theme
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12, margin = margin(b = 20)),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12, face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold", size = 12),
    axis.title = element_text(face = "bold", color = "grey30"),
    panel.grid.major = element_line(color = "grey92", linetype = "dotted"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 10),
    legend.background = element_rect(fill = "grey98", color = NA),
    legend.key = element_blank()
  ) +
  # Ensure size legend dots are neutral color
  guides(size = guide_legend(override.aes = list(fill = "grey50")))

ggsave(file.path(celltype_dirs$plots, "interaction_balloon_plot.svg"), balloon_plot, 
       width = 28, height = 16, units = "cm", dpi = 300)

# 2. Fancy Distribution Plot (Violin + Boxplot + Jitter)
# Shows the full distribution of data points behind the means
# X-axis: grouping_factor, Fill: Group, Facet: Celltype_Class
fancy_dist_plot <- ggplot(marker_expr, aes(x = !!sym(grouping_factor), y = Expression, fill = Group)) +
  geom_violin(alpha = 0.4, color = NA, trim = FALSE, scale = "width", position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15, color = "grey20", outlier.shape = NA, alpha = 0.8, position = position_dodge(width = 0.8)) +
  geom_point(aes(color = Group), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 0.5, alpha = 0.4, show.legend = FALSE) +
  facet_wrap(~ Celltype_Class, scales = "free_y", ncol = 1) +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) + 
  scale_color_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  labs(
    title = paste("Expression Distribution by", grouping_factor, ", Group and Marker Class"),
    subtitle = "Violin plots showing density with overlaid boxplots. Colored by Experimental Group.",
    x = grouping_factor,
    y = "Protein Expression",
    fill = "Exp Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    strip.background = element_rect(fill = "#2c3e50"),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(celltype_dirs$plots, "interaction_distribution_fancy.svg"), fancy_dist_plot, 
       width = 24, height = 28, units = "cm", dpi = 300)

# ----------------------------------------------------
# 10. DIRECT DIFFERENCE & PROFILE VISUALIZATION
# ----------------------------------------------------
cat("\nGenerating Difference and Profile Plots...\n")

# A. DEVIATION PLOT (Log2FC vs Average)
# This shows how each ExpGroup differs from the average expression of that marker in that celltype
# 1. Calculate Consensus Mean (Average across all groups for each Celltype/Class)
consensus_means <- marker_expr %>%
  group_by(!!sym(grouping_factor), Celltype_Class) %>%
  summarise(Global_Mean = mean(Expression, na.rm = TRUE), .groups = "drop")

# 2. Join with Group Means and Calculate Log2FC
deviation_data <- interaction_stats %>%
  left_join(consensus_means, by = c(grouping_factor, "Celltype_Class")) %>%
  mutate(
    Log2FC_vs_Mean = log2(Mean_Expr / Global_Mean),
    # Handle potential infinite values if mean is 0
    Log2FC_vs_Mean = ifelse(is.infinite(Log2FC_vs_Mean), 0, Log2FC_vs_Mean)
  )

# 3. Plot Deviation
deviation_plot <- ggplot(deviation_data, aes(x = Celltype_Class, y = Log2FC_vs_Mean, fill = Group)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", size = 0.2) +
  facet_wrap(as.formula(paste("~", grouping_factor)), scales = "fixed") +
  scale_fill_viridis_d(option = "D", begin = 0.1, end = 0.9) +
  labs(
    title = "Deviation from Consensus Profile",
    subtitle = "Log2 Fold Change of each Group vs the Average of all Groups",
    x = "Marker Class",
    y = "Log2 Fold Change (vs Average)",
    fill = "Exp Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    axis.text.x = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(file.path(celltype_dirs$plots, "interaction_deviation_lollipop.svg"), deviation_plot, 
       width = 24, height = 18, units = "cm", dpi = 300)


# B. RADAR PROFILE PLOT (Polar Coordinates)
# This shows the "shape" of the marker identity (GABA vs VGLUT1 vs VGLUT2)
# Since we have 3 axes, this forms a triangle profile for each group
radar_plot <- ggplot(interaction_stats, aes(x = Celltype_Class, y = Mean_Expr, group = Group, color = Group, fill = Group)) +
  geom_polygon(alpha = 0.2, size = 1) +
  geom_point(size = 2) +
  coord_polar() +
  facet_wrap(as.formula(paste("~", grouping_factor)), ncol = 2) +
  scale_color_viridis_d(option = "D", begin = 0.1, end = 0.9) +
  scale_fill_viridis_d(option = "D", begin = 0.1, end = 0.9) +
  labs(
    title = "Marker Identity Profiles (Radar Chart)",
    subtitle = "Shape indicates the balance between GABA, VGLUT1, and VGLUT2 markers",
    x = NULL,
    y = "Mean Expression",
    color = "Exp Group",
    fill = "Exp Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    axis.text.y = element_blank(), # Hide radial axis labels to reduce clutter
    axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", size = 11),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.major = element_line(color = "grey85", linetype = "dashed"),
    legend.position = "bottom"
  )

ggsave(file.path(celltype_dirs$plots, "interaction_radar_profile.svg"), radar_plot, 
       width = 22, height = 22, units = "cm", dpi = 300)

# C. DIFFERENCE HEATMAP
#' Create and Save Consensus Deviation Heatmap
#'
#' This section generates a high-density heatmap visualizing the Log2 Fold Changes (Log2FC)
#' of specific markers relative to the global mean across different experimental groups.
#'
#' @description
#' The plot uses `geom_tile` to display deviations:
#' *   **X-axis:** Marker Class (`Celltype_Class`).
#' *   **Y-axis:** Experimental Group (`Group`).
#' *   **Faceting:** The plot is faceted horizontally by the dynamic `grouping_factor`.
#' *   **Color Scale:** A diverging gradient is used:
#'     *   **Blue (#313695):** Negative deviation (lower than mean).
#'     *   **White:** No deviation (equal to mean).
#'     *   **Red (#a50026):** Positive deviation (higher than mean).
#'
#' @note **Why are some tiles grey?**
#' In `ggplot2`, missing values (`NA`) in the variable mapped to the `fill` aesthetic
#' (here, `Log2FC_vs_Mean`) are rendered in grey by default. If tiles appear grey,
#' it indicates that there is no calculated Log2 Fold Change data for that specific
#' combination of `Celltype_Class` and `Group`. This often happens if the protein/marker
#' was not detected in that specific group or if the calculation resulted in an undefined value.
#'
#' @output Saves the plot as "interaction_deviation_heatmap.svg" in the defined plots directory.
# A compact, high-density view of deviations
diff_heatmap <- ggplot(deviation_data, aes(x = Celltype_Class, y = Group, fill = Log2FC_vs_Mean)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(data = subset(deviation_data, is.na(Log2FC_vs_Mean)), 
            aes(label = "NA"), color = "grey50", size = 3) +
  facet_grid(as.formula(paste("~", grouping_factor)), scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#313695", mid = "white", high = "#a50026", midpoint = 0,
                       name = "Log2FC\n(vs Mean)", na.value = "grey99") +
  labs(
    title = "Consensus Deviation Heatmap",
    subtitle = "Heatmap of Log2 Fold Changes relative to the global mean across groups",
    x = "Marker Class",
    y = "Experimental Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", angle = 90, hjust = 0.5),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(celltype_dirs$plots, "interaction_deviation_heatmap.svg"), diff_heatmap,
       width = 20, height = 10, units = "cm", dpi = 300)

      # ----------------------------------------------------
      # 11. DIRECT COMPARISON: ExpGroup 2 vs ExpGroup 4
      # ----------------------------------------------------
      cat("\nGenerating Direct Comparison: ExpGroup 2 vs ExpGroup 4...\n")

      # Filter for the specific groups
      target_groups <- c("ExpGroup 2", "ExpGroup 4")
      direct_comp_data <- interaction_stats %>%
        filter(Group %in% target_groups) %>%
        select(!!sym(grouping_factor), Celltype_Class, Group, Mean_Expr) %>%
        pivot_wider(names_from = Group, values_from = Mean_Expr)

      # Check if both groups exist in the data
      if (all(target_groups %in% colnames(direct_comp_data))) {
        
        # Calculate Log2 Fold Change (Group 2 / Group 4)
        # Note: Adjust direction (2 vs 4 or 4 vs 2) as needed. Here: log2(Group2 / Group4)
        direct_comp_data <- direct_comp_data %>%
          mutate(
            Log2FC_2vs4 = log2(`ExpGroup 2` / `ExpGroup 4`),
            # Handle infinite/NaN if expression is 0
            Log2FC_2vs4 = ifelse(is.infinite(Log2FC_2vs4) | is.nan(Log2FC_2vs4), 0, Log2FC_2vs4)
          )
        
        # Save table
        write.xlsx(direct_comp_data, file = file.path(celltype_dirs$tables, "direct_comparison_ExpGroup2_vs_ExpGroup4.xlsx"))
        
        # Plot
        direct_comp_plot <- ggplot(direct_comp_data, aes(x = Celltype_Class, y = Log2FC_2vs4, fill = Log2FC_2vs4 > 0)) +
          geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
          geom_col(color = "black", size = 0.2, width = 0.7) +
          facet_wrap(as.formula(paste("~", grouping_factor)), scales = "fixed") +
          scale_fill_manual(values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"), 
                            labels = c("TRUE" = "Higher in ExpGroup 2", "FALSE" = "Higher in ExpGroup 4"),
                            name = "Direction") +
          labs(
            title = "Direct Comparison: ExpGroup 2 vs ExpGroup 4",
            subtitle = "Positive values indicate higher expression in ExpGroup 2",
            x = "Marker Class",
            y = "Log2 Fold Change (ExpGroup 2 / ExpGroup 4)"
          ) +
          theme_minimal(base_size = 14) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11),
            axis.text.x = element_text(face = "bold"),
            strip.background = element_rect(fill = "grey90", color = NA),
            strip.text = element_text(face = "bold"),
            panel.grid.major.x = element_blank(),
            legend.position = "bottom"
          )
        
        ggsave(file.path(celltype_dirs$plots, "direct_comparison_ExpGroup2_vs_ExpGroup4.svg"), 
               direct_comp_plot, width = 20, height = 15, units = "cm", dpi = 300)
        
      } else {
        cat("Warning: Could not perform direct comparison. One or both groups ('ExpGroup 2', 'ExpGroup 4') not found in data.\n")
        cat("Available groups:", paste(unique(interaction_stats$Group), collapse = ", "), "\n")
      }

cat("\n==============================================\n")
cat("ENTIRE PIPELINE COMPLETED!\n")
cat("==============================================\n")
