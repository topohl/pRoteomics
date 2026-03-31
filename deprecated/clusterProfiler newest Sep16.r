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

# Helper to save plots (use named args for ggsave)
save_plot <- function(plot, filename) {
  ggsave(filename = file.path(results_dir, filename),
         plot = plot, units = "cm", dpi = 300)
}

# Unified split dotplot builder with consistent aesthetics
# Recomputes .sign from NES for gseaResult if missing;
# returns a consistently styled split plot or NULL if not applicable
build_split_dotplot <- function(gsea_obj, cell_types, title_prefix = "GSEA of") {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))

  # Only meaningful for GSEA objects that provide NES
  if (inherits(gsea_obj, "gseaResult")) {
    if (!".sign" %in% colnames(df)) {
      if (!"NES" %in% colnames(df)) return(invisible(NULL))
      df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
      # attach to result so enrichplot::fortify sees it
      gsea_obj@result$.sign <- df$.sign
    }
    # Ensure both facets exist; if only one, it will still be styled identically
    p <- enrichplot::dotplot(
      gsea_obj,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      split = ".sign",
      font.size = 12,
      label_format = 30
    ) +
      ggplot2::facet_wrap(~ .sign, nrow = 1) +
      ggplot2::labs(
        title = paste(title_prefix, paste(cell_types, collapse = " over ")),
        x = "Gene Ratio", y = "Gene Set"
      ) +
      ggplot2::scale_fill_viridis_c(option = "cividis") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 10),
        strip.text = ggplot2::element_text(size = 12, face = "plain"),
        panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = ggplot2::element_blank()
      )
    return(p)
  }

  # Non-GSEA results have no NES/.sign; return NULL to avoid inconsistent style
  invisible(NULL)
}

# Generic dotplot that does not force split; useful for ORA or single-sign cases
plot_dot <- function(dataset, cell_types, results_dir = results_dir, basename = NULL) {
  if (is.null(dataset)) { message("plot_dot: dataset is NULL, skipping."); return(invisible(NULL)) }
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) { message("plot_dot: empty results, skipping."); return(invisible(NULL)) }

  dot_title <- paste("GSEA of", paste(cell_types, collapse = " over "))
  p <- enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                           showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = dot_title, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::scale_fill_viridis_c(option = "cividis") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 12, face = "plain"),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (!is.null(basename)) {
    save_plot(p, paste0(basename, ".svg"))
  } else {
    analysis_type <- deparse(substitute(dataset))
    save_plot(p, paste0(analysis_type, "_", paste(cell_types, collapse = "_"), ".svg"))
  }
  return(p)
}

# ----------------------------------------------------
# Define input directory and enumerate CSVs
# ----------------------------------------------------
file_direction <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/neuron-regionLayerBaseline"

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
    geneList = gene_list, ont = ont, keyType = "UNIPROT",
    minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, verbose = TRUE,
    OrgDb = organism, pAdjustMethod = "BH"
  )

  # show gse structure
  str(gse)

  plot_dot(gse, cell_types)

  # Unsimplified split with unified styling
  p_gse_split <- build_split_dotplot(gse, cell_types)
  save_plot(p_gse_split, paste0("gseGO_", ont, "_", paste(cell_types, collapse = "_"), "_dotplot_split.svg"))


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

  # Simplified split with unified styling (rebuilds .sign from NES)
  p_gse2_split <- build_split_dotplot(gse2, cell_types)
  if (!is.null(p_gse2_split)) {
  save_plot(p_gse2_split, paste0("gseGO_", ont, "_", paste(cell_types, collapse = "_"), "_simplified_split.svg"))
  }

  # show gse2 structure
  str(gse2)

  plot_dot(gse2, cell_types)

  # ----------------------------------------------------
  # ORA on top genes
  # ----------------------------------------------------
  ora <- enrichGO(
    gene = top_genes, ont = "BP", keyType = "UNIPROT",
    minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
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
}























############################################################ var 1
# clusterProfiler + WGCNA modules
# Robust ID handling + improved folder structure
############################################################ var 1

## ---------- Package setup (unchanged from before) ----------
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
  required_packages <- c(
    "clusterProfiler", "pathview", "enrichplot", "DOSE", "ggplot2", "ggnewscale",
    "cowplot", "ggridges", "europepmc", "ggpubr", "ggrepel", "ggsci", "ggthemes",
    "ggExtra", "ggforce", "ggalluvial", "lattice", "latticeExtra", "BiocManager",
    "org.Mm.eg.db", "ggplotify", "svglite", "tidyr", "dplyr", "pheatmap", "proxy",
    "tibble", "stringr", "tools", "AnnotationDbi", "igraph", "readr", "forcats", "vroom"
  )
  bioc_packages <- c("clusterProfiler", "pathview", "enrichplot", "DOSE", "org.Mm.eg.db")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Helpers ----------
# results_dir is set within the loop; save_plot writes to subfolders
save_plot <- function(plot, filename, subdir = "", width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  out_dir <- file.path(results_dir, subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(out_dir, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}

build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::scale_fill_viridis_c(option = "cividis") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::scale_fill_viridis_c(option = "cividis") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

## ---------- Load WGCNA modules ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/tables_modules"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv))
modules_df <- modules_df |> dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "")
modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(as.character(x)))

module_sizes <- tibble::tibble(Module = names(modules_list),
                               Size = vapply(modules_list, length, integer(1))) |>
  dplyr::arrange(dplyr::desc(Size))

## ---------- Fisher DEG enrichment helper ----------
fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
  deg_set <- intersect(unique(deg_genes), universe_genes)
  bg_set  <- setdiff(universe_genes, deg_set)
  out <- lapply(names(modules_list), function(m) {
    mod_genes <- intersect(modules_list[[m]], universe_genes)
    a <- length(intersect(mod_genes, deg_set))
    b <- length(setdiff(mod_genes, deg_set))
    c <- length(setdiff(deg_set, mod_genes))
    d <- length(setdiff(bg_set, mod_genes))
    mat <- matrix(c(a,b,c,d), nrow=2)
    ft <- fisher.test(mat, alternative = "greater")
    data.frame(
      Module = m, InModule_DE = a, InModule_NotDE = b,
      NotInModule_DE = c, NotInModule_NotDE = d,
      OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  res <- dplyr::bind_rows(out)
  res$Padj <- p.adjust(res$Pvalue, method = "BH")
  dplyr::arrange(res, Padj, Pvalue)
}

## ---------- ID helpers ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Input comparisons ----------
file_direction <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/neuron-regionLayerBaseline"
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # context label from filename
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L)) {
    context <- name_no_ext
  } else if (length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories (new structure) -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- file.path(working_dir, "Results", file_tag)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  dir_meta <- file.path(results_dir, "00_meta");      dir.create(dir_meta, showWarnings = FALSE)
  dir_glob <- file.path(results_dir, "10_global");    dir.create(dir_glob, showWarnings = FALSE)
  dir_mod  <- file.path(results_dir, "20_modules");   dir.create(dir_mod, showWarnings = FALSE)
  dir_kegg <- file.path(results_dir, "90_pathview");  dir.create(dir_kegg, showWarnings = FALSE)

  # per-ontology subfolders under global and modules
  gl_GO_BP <- file.path(dir_glob, "GO_BP"); dir.create(gl_GO_BP, showWarnings = FALSE, recursive = TRUE)
  gl_GO_MF <- file.path(dir_glob, "GO_MF"); dir.create(gl_GO_MF, showWarnings = FALSE, recursive = TRUE)
  gl_GO_CC <- file.path(dir_glob, "GO_CC"); dir.create(gl_GO_CC, showWarnings = FALSE, recursive = TRUE)
  gl_KEGG  <- file.path(dir_glob, "KEGG");  dir.create(gl_KEGG,  showWarnings = FALSE, recursive = TRUE)

  mod_GO_BP <- file.path(dir_mod, "GO_BP"); dir.create(mod_GO_BP, showWarnings = FALSE, recursive = TRUE)
  mod_GO_MF <- file.path(dir_mod, "GO_MF"); dir.create(mod_GO_MF, showWarnings = FALSE, recursive = TRUE)
  mod_GO_CC <- file.path(dir_mod, "GO_CC"); dir.create(mod_GO_CC, showWarnings = FALSE, recursive = TRUE)

  # write module summary once per comparison
  readr::write_csv(module_sizes, file.path(dir_mod, "module_list.csv"))

  # save session info and parameters
  capture.output(sessionInfo(), file = file.path(dir_meta, "sessionInfo.txt"))
  writeLines(c(
    paste0("data_path: ", data_path),
    paste0("context: ", context),
    paste0("timestamp: ", Sys.time())
  ), con = file.path(dir_meta, "params.txt"))

  ## ----- Data load with robust ID hygiene -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  # expect first column = UNIPROT id; standardize names
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))

  # remove rows with missing id or log2fc
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id))

  # resolve duplicates: keep the entry with max absolute log2fc per id
  df <- df |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # Build ranked vector (names length == values length guaranteed)
  gene_list_tbl <- df |>
    dplyr::select(id, log2fc) |>
    dplyr::arrange(dplyr::desc(log2fc))
  gene_list <- gene_list_tbl$log2fc
  names(gene_list) <- gene_list_tbl$id

  universe_genes <- names(gene_list)

  # DEGs for ORA (adjust thresholds as needed)
  deg_tbl <- df |>
    dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  organism <- "org.Mm.eg.db"
  # change ontologies here if needed
  onts <- c("BP", "MF", "CC")

  ## ----- GLOBAL: GO GSEA and ORA -----
  for (ont in onts) {
    # global GSEA
    gse <- tryCatch(clusterProfiler::gseGO(
      geneList = gene_list, ont = ont, keyType = "UNIPROT",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, verbose = FALSE,
      OrgDb = organism, pAdjustMethod = "BH"
    ), error = function(e) NULL)

    subdir_glob <- switch(ont,
                          "BP" = file.path("10_global","GO_BP"),
                          "MF" = file.path("10_global","GO_MF"),
                          "CC" = file.path("10_global","GO_CC"))

    if (!is.null(gse)) {
      readr::write_csv(gse@result, file.path(results_dir, subdir_glob,
                                             paste0("gseGO_", ont, "_", file_tag, ".csv")))
      save_plot(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)),
                paste0("gseGO_", ont, "_split_", file_tag, ".svg"), subdir = subdir_glob)
      save_plot(enrichplot::emapplot(enrichplot::pairwise_termsim(gse), showCategory = 10),
                paste0("gseGO_", ont, "_emap_", file_tag, ".svg"), subdir = subdir_glob)
      save_plot(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"), subdir = subdir_glob)
      save_plot(enrichplot::ridgeplot(gse) + ggplot2::labs(x = "Enrichment Distribution"),
                paste0("gseGO_", ont, "_ridge_", file_tag, ".svg"), subdir = subdir_glob)
      if (nrow(gse@result) > 0) {
        save_plot(enrichplot::gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1),
                  paste0("gseGO_", ont, "_gseaplot_", file_tag, ".svg"), subdir = subdir_glob)
      }
      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) NULL)
      if (!is.null(gse2)) {
        readr::write_csv(gse2@result, file.path(results_dir, subdir_glob,
                                                paste0("gseGO_", ont, "_simplified_", file_tag, ".csv")))
        save_plot(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)),
                  paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"), subdir = subdir_glob)
      }
    }

    # global ORA
    if (length(top_genes) >= 10) {
      ora <- tryCatch(clusterProfiler::enrichGO(
        gene = top_genes, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        OrgDb = organism, pAdjustMethod = "none"
      ), error = function(e) NULL)
      if (!is.null(ora)) {
        readr::write_csv(ora@result, file.path(results_dir, subdir_glob,
                                               paste0("ORA_", ont, "_", file_tag, ".csv")))
        save_plot(plot_dot(ora, paste0("ORA (", ont, ") — ", context)),
                  paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"), subdir = subdir_glob)
        save_plot(enrichplot::emapplot(enrichplot::pairwise_termsim(ora), showCategory = 10),
                  paste0("ORA_", ont, "_emap_", file_tag, ".svg"), subdir = subdir_glob)
        save_plot(enrichplot::cnetplot(ora, showCategory = 5),
                  paste0("ORA_", ont, "_cnet_", file_tag, ".svg"), subdir = subdir_glob)
      }
    }
  } # end ont loop

  ## ----- GLOBAL: KEGG GSEA -----
  ids <- clusterProfiler::bitr(universe_genes, fromType = "UNIPROT", toType = "ENTREZID",
                               OrgDb = "org.Mm.eg.db")
  ids <- ids[!duplicated(ids$UNIPROT), ]
  df2 <- dplyr::inner_join(df, ids, by = c("id" = "UNIPROT"))
  kegg_gene_list <- df2$log2fc
  names(kegg_gene_list) <- df2$ENTREZID
  kegg_gene_list <- sort(kegg_gene_list[!is.na(kegg_gene_list)], decreasing = TRUE)
  kk2 <- tryCatch(clusterProfiler::gseKEGG(
    geneList = kegg_gene_list, organism = "mmu",
    minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
    keyType = "ncbi-geneid"
  ), error = function(e) NULL)
  if (!is.null(kk2)) {
    readr::write_csv(kk2@result, file.path(gl_KEGG, paste0("gseKEGG_", file_tag, ".csv")))
    save_plot(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                  showCategory = 10, label_format = 30),
              paste0("gseKEGG_dotplot_", file_tag, ".svg"), subdir = "10_global/KEGG")
    save_plot(enrichplot::emapplot(enrichplot::pairwise_termsim(kk2), showCategory = 10),
              paste0("gseKEGG_emap_", file_tag, ".svg"), subdir = "10_global/KEGG")
    save_plot(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
              paste0("gseKEGG_cnet_", file_tag, ".svg"), subdir = "10_global/KEGG")
    save_plot(enrichplot::ridgeplot(kk2) + ggplot2::labs(x = "Enrichment distribution"),
              paste0("gseKEGG_ridge_", file_tag, ".svg"), subdir = "10_global/KEGG")
    if (nrow(kk2@result) > 0) {
      save_plot(enrichplot::gseaplot(kk2, by = "all", title = kk2@result$Description[1], geneSetID = 1),
                paste0("gseKEGG_gseaplot_", file_tag, ".svg"), subdir = "10_global/KEGG")
    }
  }
  # Pathview
  path_ids <- c("mmu04110","mmu04115","mmu04114","mmu04113","mmu04112","mmu04111",
                "mmu04116","mmu04117","mmu04118","mmu04119","mmu04720","mmu04721",
                "mmu04722","mmu04725","mmu04726","mmu04727","mmu04724","mmu04080",
                "mmu00030","mmu04151")
  oldwd <- getwd(); setwd(dir_kegg)
  lapply(path_ids, function(pid) {
    pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                       low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg")
  })
  setwd(oldwd)

  ## ----- MODULE-AWARE: Fisher + per-module GSEA/ORA -----
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  readr::write_csv(fisher_tbl, file.path(dir_mod, "fisher_deg", paste0("fisher_", file_tag, ".csv")))
  # quick barplot of top modules
  dir.create(file.path(dir_mod, "fisher_deg"), showWarnings = FALSE)
  if (nrow(fisher_tbl) > 0) {
    topF <- head(fisher_tbl[order(fisher_tbl$Padj, fisher_tbl$Pvalue), ], 20)
    pbar <- ggplot2::ggplot(topF, ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                               y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Module DEG enrichment — ", context),
                    x = "Module", y = "-log10(BH-adjusted p)")
    save_plot(pbar, paste0("fisher_bar_", file_tag, ".svg"), subdir = "20_modules/fisher_deg", width = 24, height = 20)
  }

  sel_modules <- fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)]
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) next

    # choose ontology subdir
    for (ont in onts) {
      m_base <- file.path("20_modules", paste0("GO_", ont), paste0("Module_", m))
      dir.create(file.path(results_dir, m_base), recursive = TRUE, showWarnings = FALSE)

      # per-module ranked list
      m_gene_list <- gene_list[names(gene_list) %in% m_genes]
      if (length(m_gene_list) >= 20) {
        m_gse <- tryCatch(clusterProfiler::gseGO(
          geneList = sort(m_gene_list, decreasing = TRUE),
          ont = ont, keyType = "UNIPROT",
          minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) NULL)
        if (!is.null(m_gse) && nrow(m_gse@result) > 0) {
          readr::write_csv(m_gse@result, file.path(results_dir, m_base,
                                                   paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv")))
          save_plot(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                    paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"),
                    subdir = m_base)
          save_plot(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust",
                                        showCategory = 10, label_format = 30),
                    paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"),
                    subdir = m_base)
        }
      }

      # per-module ORA on DEGs
      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 10) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none"
        ), error = function(e) NULL)
        if (!is.null(m_ora) && nrow(m_ora@result) > 0) {
          readr::write_csv(m_ora@result, file.path(results_dir, m_base,
                                                   paste0("ORA_", ont, "_", m, "_", file_tag, ".csv")))
          save_plot(plot_dot(m_ora, paste0("Module ", m, " — ORA (", ont, ") — ", context)),
                    paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"),
                    subdir = m_base)
          save_plot(enrichplot::emapplot(enrichplot::pairwise_termsim(m_ora), showCategory = 10),
                    paste0("ORA_", ont, "_", m, "_emap_", file_tag, ".svg"),
                    subdir = m_base)
          save_plot(enrichplot::cnetplot(m_ora, showCategory = 5),
                    paste0("ORA_", ont, "_", m, "_cnet_", file_tag, ".svg"),
                    subdir = m_base)
        }
      }
    } # end ont loop
  } # end modules loop
} # end comparisons loop






















############################################################
# clusterProfiler + WGCNA modules (mouse)
# EntryName->Accession mapping for modules, fgsea GSEA, KEGG validation
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","UniProt.ws", "vroom"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST","UniProt.ws")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Plot helpers ----------
save_plot <- function(plot, filename, subdir = "", width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  out_dir <- file.path(results_dir, subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(out_dir, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

## ---------- Load WGCNA modules (EntryName format, e.g., 5NT1A_MOUSE) ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/tables_modules"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "")

modules_list_entry <- split(modules_df$Gene, modules_df$Module)
modules_list_entry <- lapply(modules_list_entry, function(x) unique(as.character(x)))

module_sizes_entry <- tibble::tibble(Module = names(modules_list_entry),
                                     Size_EntryNames = vapply(modules_list_entry, length, integer(1)))

## ---------- Convert EntryName -> Accession ----------
convert_entry_names_to_accessions <- function(entry_vec, taxon_id = 10090) {
  entry_vec <- unique(entry_vec[!is.na(entry_vec) & nzchar(entry_vec)])
  if (!length(entry_vec)) return(character(0))

  accs <- NULL
  # Try UniProt.ws
  try({
    up <- UniProt.ws::UniProt.ws(taxId = taxon_id)
    res <- UniProt.ws::select(up,
                              keys = entry_vec,
                              columns = c("accession","id"),
                              keytype = "UniProtKB")
    res <- res[!is.na(res$accession) & !is.na(res$id), ]
    res <- res[!duplicated(res$id), ]  # keep first accession per entry
    map <- res$accession
    names(map) <- res$id
    accs <- unname(map[entry_vec])
  }, silent = TRUE)

  # Fallback: UniProt REST search if needed
  if (is.null(accs) || all(is.na(accs))) {
    q <- paste(paste0('(', paste(sprintf('id:%s', entry_vec), collapse = ' OR '), ')'),
               sprintf('AND organism_id:%d', taxon_id))
    url <- "https://rest.uniprot.org/uniprotkb/search?fields=accession,id&format=tsv&query="
    tsv <- try(readr::read_tsv(paste0(url, utils::URLencode(q, reserved = TRUE)),
                               show_col_types = FALSE), silent = TRUE)
    if (!inherits(tsv, "try-error") && nrow(tsv)) {
      map <- tsv$Entry
      names(map) <- tsv$`Entry Name`
      accs <- unname(map[entry_vec])
    }
  }
  accs[is.na(accs)] <- NA_character_
  accs
}

modules_list <- lapply(modules_list_entry, function(x) convert_entry_names_to_accessions(x, taxon_id = 10090))
modules_list <- lapply(modules_list, function(v) unique(v[!is.na(v) & nzchar(v)]))

module_sizes_acc <- tibble::tibble(Module = names(modules_list),
                                   Size_Accessions = vapply(modules_list, length, integer(1)))
module_sizes <- dplyr::left_join(module_sizes_entry, module_sizes_acc, by = "Module")

## ---------- Fisher enrichment helper ----------
fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
  deg_set <- intersect(unique(deg_genes), universe_genes)
  bg_set  <- setdiff(universe_genes, deg_set)
  out <- lapply(names(modules_list), function(m) {
    mod_genes <- intersect(modules_list[[m]], universe_genes)
    a <- length(intersect(mod_genes, deg_set))
    b <- length(setdiff(mod_genes, deg_set))
    c <- length(setdiff(deg_set, mod_genes))
    d <- length(setdiff(bg_set, mod_genes))
    mat <- matrix(c(a,b,c,d), nrow=2)
    ft <- fisher.test(mat, alternative = "greater")
    data.frame(
      Module = m, InModule_DE = a, InModule_NotDE = b,
      NotInModule_DE = c, NotInModule_NotDE = d,
      OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  res <- dplyr::bind_rows(out)
  res$Padj <- p.adjust(res$Pvalue, method = "BH")
  dplyr::arrange(res, Padj, Pvalue)
}

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Inputs ----------
file_direction <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/neuron-regionLayerBaseline"
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- file.path(working_dir, "Results", file_tag)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  dir_meta <- file.path(results_dir, "00_meta");      dir.create(dir_meta, recursive = TRUE, showWarnings = FALSE)
  dir_glob <- file.path(results_dir, "10_global");    dir.create(dir_glob, recursive = TRUE, showWarnings = FALSE)
  dir_mod  <- file.path(results_dir, "20_modules");   dir.create(dir_mod,  recursive = TRUE, showWarnings = FALSE)
  dir_kegg <- file.path(results_dir, "90_pathview");  dir.create(dir_kegg, recursive = TRUE, showWarnings = FALSE)

  gl_GO_BP <- file.path(dir_glob, "GO_BP"); dir.create(gl_GO_BP, recursive = TRUE, showWarnings = FALSE)
  gl_GO_MF <- file.path(dir_glob, "GO_MF"); dir.create(gl_GO_MF, recursive = TRUE, showWarnings = FALSE)
  gl_GO_CC <- file.path(dir_glob, "GO_CC"); dir.create(gl_GO_CC, recursive = TRUE, showWarnings = FALSE)
  gl_KEGG  <- file.path(dir_glob, "KEGG");  dir.create(gl_KEGG,  recursive = TRUE, showWarnings = FALSE)

  mod_GO_BP <- file.path(dir_mod, "GO_BP"); dir.create(mod_GO_BP, recursive = TRUE, showWarnings = FALSE)
  mod_GO_MF <- file.path(dir_mod, "GO_MF"); dir.create(mod_GO_MF, recursive = TRUE, showWarnings = FALSE)
  mod_GO_CC <- file.path(dir_mod, "GO_CC"); dir.create(mod_GO_CC, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir_mod, "fisher_deg"), recursive = TRUE, showWarnings = FALSE)

  capture.output(sessionInfo(), file = file.path(dir_meta, "sessionInfo.txt"))
  readr::write_csv(module_sizes, file.path(dir_mod, "module_list_map_stats.csv"))
  writeLines(c(
    paste0("data_path: ", data_path),
    paste0("context: ", context),
    paste0("timestamp: ", Sys.time())
  ), con = file.path(dir_meta, "params.txt"))

  ## ----- Load comparison; robust ID hygiene (assume UNIPROT accessions in col1) -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  gene_list_tbl <- df |>
    dplyr::select(id, log2fc) |>
    dplyr::arrange(dplyr::desc(log2fc))
  gene_list <- gene_list_tbl$log2fc
  names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)

  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  organism <- "org.Mm.eg.db"
  onts <- c("BP") # extend to c("BP","MF","CC") if desired

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    subdir_glob <- switch(ont,
      "BP" = file.path("10_global","GO_BP"),
      "MF" = file.path("10_global","GO_MF"),
      "CC" = file.path("10_global","GO_CC")
    )

    gse <- tryCatch(clusterProfiler::gseGO(
      geneList = gene_list, ont = ont, keyType = "UNIPROT",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
      by = "fgsea", nPermSimple = 10000, eps = 0,
      verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
    ), error = function(e) NULL)
    if (!is.null(gse)) {
      readr::write_csv(gse@result, file.path(results_dir, subdir_glob, paste0("gseGO_", ont, "_", file_tag, ".csv")))
      save_plot(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)),
                paste0("gseGO_", ont, "_split_", file_tag, ".svg"), subdir = subdir_glob)
      save_plot(enrichplot::emapplot(enrichplot::pairwise_termsim(gse), showCategory = 10),
                paste0("gseGO_", ont, "_emap_", file_tag, ".svg"), subdir = subdir_glob)
      save_plot(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"), subdir = subdir_glob)
      save_plot(enrichplot::ridgeplot(gse) + ggplot2::labs(x = "Enrichment Distribution"),
                paste0("gseGO_", ont, "_ridge_", file_tag, ".svg"), subdir = subdir_glob)
      if (nrow(gse@result) > 0) {
        save_plot(enrichplot::gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1),
                  paste0("gseGO_", ont, "_gseaplot_", file_tag, ".svg"), subdir = subdir_glob)
      }
      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) NULL)
      if (!is.null(gse2)) {
        readr::write_csv(gse2@result, file.path(results_dir, subdir_glob, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv")))
        save_plot(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)),
                  paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"), subdir = subdir_glob)
      }
    }

    if (length(top_genes) >= 10) {
      ora <- tryCatch(clusterProfiler::enrichGO(
        gene = top_genes, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
      ), error = function(e) NULL)
      if (!is.null(ora)) {
        readr::write_csv(ora@result, file.path(results_dir, subdir_glob, paste0("ORA_", ont, "_", file_tag, ".csv")))
        save_plot(plot_dot(ora, paste0("ORA (", ont, ") — ", context)),
                  paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"), subdir = subdir_glob)
        save_plot(enrichplot::emapplot(enrichplot::pairwise_termsim(ora), showCategory = 10),
                  paste0("ORA_", ont, "_emap_", file_tag, ".svg"), subdir = subdir_glob)
        save_plot(enrichplot::cnetplot(ora, showCategory = 5),
                  paste0("ORA_", ont, "_cnet_", file_tag, ".svg"), subdir = subdir_glob)
      }
    }
  }

  ## ----- GLOBAL: KEGG (mapIds to 1:1) -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  keep <- !is.na(entrez_map)
  kegg_gene_list <- gene_list[keep]
  names(kegg_gene_list) <- unname(entrez_map[keep])
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

  kk2 <- tryCatch(clusterProfiler::gseKEGG(
    geneList = kegg_gene_list, organism = "mmu",
    minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
    keyType = "ncbi-geneid"
  ), error = function(e) NULL)
  if (!is.null(kk2)) {
    readr::write_csv(kk2@result, file.path(gl_KEGG, paste0("gseKEGG_", file_tag, ".csv")))
    save_plot(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust", showCategory = 10, label_format = 30),
              paste0("gseKEGG_dotplot_", file_tag, ".svg"), subdir = "10_global/KEGG")
    save_plot(enrichplot::emapplot(enrichplot::pairwise_termsim(kk2), showCategory = 10),
              paste0("gseKEGG_emap_", file_tag, ".svg"), subdir = "10_global/KEGG")
    save_plot(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
              paste0("gseKEGG_cnet_", file_tag, ".svg"), subdir = "10_global/KEGG")
    save_plot(enrichplot::ridgeplot(kk2) + ggplot2::labs(x = "Enrichment distribution"),
              paste0("gseKEGG_ridge_", file_tag, ".svg"), subdir = "10_global/KEGG")
    if (nrow(kk2@result) > 0) {
      save_plot(enrichplot::gseaplot(kk2, by = "all", title = kk2@result$Description[1], geneSetID = 1),
                paste0("gseKEGG_gseaplot_", file_tag, ".svg"), subdir = "10_global/KEGG")
    }
  }

  ## ----- PATHVIEW with validated KEGG IDs -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  oldwd <- getwd(); setwd(dir_kegg)
  invisible(lapply(pv_ids, function(pid) {
    try(pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"), silent = TRUE)
  }))
  setwd(oldwd)

  ## ----- MODULE-AWARE: Fisher + per-module GO -----
  # Overlap diagnostics per module
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  readr::write_csv(overlap_counts, file.path(dir_mod, "module_overlap_counts.csv"))

  # Fisher enrichment
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  readr::write_csv(fisher_tbl, file.path(results_dir, "20_modules", "fisher_deg", paste0("fisher_", file_tag, ".csv")))

  sel_modules <- fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)]
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) next

    for (ont in onts) {
      m_base <- file.path("20_modules", paste0("GO_", ont), paste0("Module_", m))
      dir.create(file.path(results_dir, m_base), recursive = TRUE, showWarnings = FALSE)

      # Per-module ranked list
      m_gene_list <- gene_list[names(gene_list) %in% m_genes]

      # Per-module GSEA (lower minGSSize for module-restricted)
      if (length(m_gene_list) >= 10) {
        m_gse <- tryCatch(clusterProfiler::gseGO(
          geneList = sort(m_gene_list, decreasing = TRUE),
          ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) NULL)
        if (!is.null(m_gse) && nrow(m_gse@result) > 0) {
          readr::write_csv(m_gse@result, file.path(results_dir, m_base, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv")))
          save_plot(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                    paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"), subdir = m_base)
          save_plot(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust", showCategory = 10, label_format = 30),
                    paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"), subdir = m_base)
        } else {
          readr::write_csv(tibble::tibble(note = "no GSEA terms"), file.path(results_dir, m_base, paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv")))
        }
      }

      # Per-module ORA (universe explicit; lower minGSSize)
      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) NULL)
        if (!is.null(m_ora) && nrow(m_ora@result) > 0) {
          readr::write_csv(m_ora@result, file.path(results_dir, m_base, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv")))
          save_plot(plot_dot(m_ora, paste0("Module ", m, " — ORA (", ont, ") — ", context)),
                    paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"), subdir = m_base)
          save_plot(enrichplot::emapplot(enrichplot::pairwise_termsim(m_ora), showCategory = 10),
                    paste0("ORA_", ont, "_", m, "_emap_", file_tag, ".svg"), subdir = m_base)
          save_plot(enrichplot::cnetplot(m_ora, showCategory = 5),
                    paste0("ORA_", ont, "_", m, "_cnet_", file_tag, ".svg"), subdir = m_base)
        } else {
          readr::write_csv(tibble::tibble(note = "no ORA terms"), file.path(results_dir, m_base, paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv")))
        }
      } else {
        readr::write_csv(tibble::tibble(note = "insufficient DEGs for ORA"),
                         file.path(results_dir, m_base, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv")))
      }
    } # ont
  } # module
} # comparison



























############################################################
# clusterProfiler + WGCNA modules (mouse)
# EntryName->Accession mapping, diagnostics, fgsea GSEA, KEGG validation
# Expanded subfolders: Tables vs Plots for global and per-module outputs
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","UniProt.ws", "vroom"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST","UniProt.ws")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Directory helpers ----------
mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP   = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                  P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF   = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                  P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC   = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                  P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                  P = mk(results_dir, "10_global","KEGG","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots"))
    ),
    pathview = mk(results_dir, "90_pathview")
  )
}
save_table <- function(df, dirpath, fname) {
  readr::write_csv(df, file.path(dirpath, fname))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}

## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text  = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || !nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

## ---------- Load WGCNA modules (EntryName format) ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "")

modules_list_entry <- split(modules_df$Gene, modules_df$Module)
modules_list_entry <- lapply(modules_list_entry, function(x) unique(as.character(x)))

module_sizes_entry <- tibble::tibble(Module = names(modules_list_entry),
                                     N_entry_names = vapply(modules_list_entry, length, integer(1)))

## ---------- Convert EntryName -> Accession ----------
convert_entry_names_to_accessions <- function(entry_vec, taxon_id = 10090) {
  entry_vec <- unique(entry_vec[!is.na(entry_vec) & nzchar(entry_vec)])
  if (!length(entry_vec)) return(character(0))

  accs <- NULL
  # Try UniProt.ws first
  try({
    up <- UniProt.ws::UniProt.ws(taxId = taxon_id)
    res <- UniProt.ws::select(up,
                              keys = entry_vec,
                              columns = c("accession","id"),
                              keytype = "UniProtKB")
    res <- res[!is.na(res$accession) & !is.na(res$id), ]
    res <- res[!duplicated(res$id), ]  # one accession per entryname
    map <- res$accession
    names(map) <- res$id
    accs <- unname(map[entry_vec])
  }, silent = TRUE)

  # Fallback: UniProt REST search
  if (is.null(accs) || all(is.na(accs))) {
    q <- paste(paste0('(', paste(sprintf('id:%s', entry_vec), collapse = ' OR '), ')'),
               sprintf('AND organism_id:%d', taxon_id))
    url <- "https://rest.uniprot.org/uniprotkb/search?fields=accession,id&format=tsv&query="
    tsv <- try(readr::read_tsv(paste0(url, utils::URLencode(q, reserved = TRUE)),
                               show_col_types = FALSE), silent = TRUE)
    if (!inherits(tsv, "try-error") && nrow(tsv)) {
      map <- tsv$Entry
      names(map) <- tsv$`Entry Name`
      accs <- unname(map[entry_vec])
    }
  }
  accs[is.na(accs)] <- NA_character_
  accs
}
modules_list <- lapply(modules_list_entry, function(x) convert_entry_names_to_accessions(x, taxon_id = 10090))
modules_list <- lapply(modules_list, function(v) unique(v[!is.na(v) & nzchar(v)]))

## ---------- Mapping diagnostics ----------
map_stats_per_module <- tibble::tibble(Module = character(),
                                       N_entry_names = integer(),
                                       N_mapped_accessions = integer(),
                                       N_unmapped = integer(),
                                       Pct_mapped = numeric())
unmapped_tbl <- tibble::tibble(Module = character(), EntryName = character(), Accession = character())

map_stats_per_module <- tibble::tibble(
  Module = names(modules_list_entry),
  N_entry_names = vapply(modules_list_entry, length, integer(1)),
  N_mapped_accessions = vapply(modules_list, length, integer(1))
) |>
  dplyr::mutate(N_unmapped = N_entry_names - N_mapped_accessions,
                Pct_mapped = round(100 * N_mapped_accessions / pmax(1, N_entry_names), 1))

entry_to_acc <- tryCatch({
  lapply(names(modules_list_entry), function(m) {
    e <- modules_list_entry[[m]]
    a <- convert_entry_names_to_accessions(e, taxon_id = 10090)
    tibble::tibble(Module = m, EntryName = e, Accession = a)
  }) |> dplyr::bind_rows()
}, error = function(e) tibble::tibble(Module = character(), EntryName = character(), Accession = character()))
unmapped_tbl <- entry_to_acc |>
  dplyr::filter(is.na(Accession) | !nzchar(Accession)) |>
  dplyr::arrange(Module, EntryName)

total_entries <- sum(map_stats_per_module$N_entry_names)
total_mapped  <- sum(map_stats_per_module$N_mapped_accessions)
total_unmapped <- total_entries - total_mapped
pct <- if (total_entries > 0) 100 * total_mapped / total_entries else 0
message(sprintf("[UniProt mapping] EntryNames: %d, mapped: %d (%.1f%%), unmapped: %d",
                total_entries, total_mapped, pct, total_unmapped))

## ---------- Fisher enrichment helper ----------
fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
  deg_set <- intersect(unique(deg_genes), universe_genes)
  bg_set  <- setdiff(universe_genes, deg_set)
  out <- lapply(names(modules_list), function(m) {
    mod_genes <- intersect(modules_list[[m]], universe_genes)
    a <- length(intersect(mod_genes, deg_set))
    b <- length(setdiff(mod_genes, deg_set))
    c <- length(setdiff(deg_set, mod_genes))
    d <- length(setdiff(bg_set, mod_genes))
    mat <- matrix(c(a,b,c,d), nrow=2)
    ft <- fisher.test(mat, alternative = "greater")
    data.frame(
      Module = m, InModule_DE = a, InModule_NotDE = b,
      NotInModule_DE = c, NotInModule_NotDE = d,
      OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
  res <- dplyr::bind_rows(out)
  res$Padj <- p.adjust(res$Pvalue, method = "BH")
  dplyr::arrange(res, Padj, Pvalue)
}

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Inputs ----------
comparison <- "neuron-regionLayerBaseline"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Save mapping diagnostics (once per comparison for archiving)
  capture.output(sessionInfo(), file = file.path(paths$meta, "sessionInfo.txt"))
  save_table(map_stats_per_module, paths$modules$root, "module_entryname_to_accession_stats.csv")
  if (nrow(unmapped_tbl) > 0) {
    save_table(unmapped_tbl, paths$modules$root, "module_unmapped_entrynames.csv")
  } else {
    save_table(tibble::tibble(note = "All EntryNames mapped to accessions"), paths$modules$root, "module_unmapped_entrynames.csv")
  }
  writeLines(c(
    paste0("data_path: ", data_path),
    paste0("context: ", context),
    paste0("timestamp: ", Sys.time())
  ), con = file.path(paths$meta, "params.txt"))

  ## ----- Load comparison (assume UniProt accessions in col1) -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  gene_list_tbl <- df |>
    dplyr::select(id, log2fc) |>
    dplyr::arrange(dplyr::desc(log2fc))
  gene_list <- gene_list_tbl$log2fc
  names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)

  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  organism <- "org.Mm.eg.db"
  onts <- c("BP") # extend to c("BP","MF","CC") if needed

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    gse <- tryCatch(clusterProfiler::gseGO(
      geneList = gene_list, ont = ont, keyType = "UNIPROT",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
      by = "fgsea", nPermSimple = 10000, eps = 0,
      verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
    ), error = function(e) NULL)
    if (!is.null(gse)) {
      save_table(gse@result, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))
      save_plot_dir(enrichplot::emapplot(enrichplot::pairwise_termsim(gse), showCategory = 10), dstP,
                    paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
      save_plot_dir(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list), dstP,
                    paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))
      save_plot_dir(enrichplot::ridgeplot(gse) + ggplot2::labs(x = "Enrichment Distribution"), dstP,
                    paste0("gseGO_", ont, "_ridge_", file_tag, ".svg"))
      if (nrow(gse@result) > 0) {
        save_plot_dir(enrichplot::gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1), dstP,
                      paste0("gseGO_", ont, "_gseaplot_", file_tag, ".svg"))
      }
      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) NULL)
      if (!is.null(gse2)) {
        save_table(gse2@result, dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
      }
    }

    if (length(top_genes) >= 10) {
      ora <- tryCatch(clusterProfiler::enrichGO(
        gene = top_genes, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
      ), error = function(e) NULL)
      if (!is.null(ora)) {
        save_table(ora@result, dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        save_plot_dir(enrichplot::emapplot(enrichplot::pairwise_termsim(ora), showCategory = 10), dstP,
                      paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
        save_plot_dir(enrichplot::cnetplot(ora, showCategory = 5), dstP,
                      paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
      }
    }
  }

  ## ----- GLOBAL: KEGG (mapIds to 1:1) -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  keep <- !is.na(entrez_map)
  kegg_gene_list <- gene_list[keep]
  names(kegg_gene_list) <- unname(entrez_map[keep])
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

  kk2 <- tryCatch(clusterProfiler::gseKEGG(
    geneList = kegg_gene_list, organism = "mmu",
    minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
    keyType = "ncbi-geneid"
  ), error = function(e) NULL)
  if (!is.null(kk2)) {
    save_table(kk2@result, paths$glob$KEGG$T,
               paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = 10, label_format = 30),
                  paths$glob$KEGG$P, 
                  paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    save_plot_dir(enrichplot::emapplot(enrichplot::pairwise_termsim(kk2),
                                       showCategory = 10),
                  paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
    save_plot_dir(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                  paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))
    save_plot_dir(enrichplot::ridgeplot(kk2) +
                    ggplot2::labs(x = "Enrichment distribution"),
                  paths$glob$KEGG$P, paste0("gseKEGG_ridge_", file_tag, ".svg"))
    if (nrow(kk2@result) > 0) {
      save_plot_dir(enrichplot::gseaplot(kk2, by = "all",
                                         title = kk2@result$Description[1],
                                         geneSetID = 1),
                    paths$glob$KEGG$P,
                    paste0("gseKEGG_gseaplot_", file_tag, ".svg"))
    }
  }

  ## ----- PATHVIEW with validated KEGG IDs + status log -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && nrow(kk2@result) > 0) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  if (length(kegg_gene_list) > 0 && all(!is.na(as.numeric(kegg_gene_list)))) {
    writeLines(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
  } else {
    writeLines("No valid ENTREZ-mapped genes for Pathview (kegg_gene_list empty).",
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
  }

  ## ----- MODULE-AWARE: overlaps, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")

  # UpSet-style overlap (top modules by DEG overlap) with guards
  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), ]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids
  save_table(mat, paths$modules$overlap$T, paste0("upset_input_", file_tag, ".csv"))

  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))
  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      library(ComplexUpset)
      p_up <- ggplot(mat, aes(x = NULL)) +
        ComplexUpset::upset(
          intersect = colnames(mat)[colnames(mat) != "id"],
          base_annotations = list(
            'Intersection size' = ComplexUpset::intersection_size()
          ),
          set_sizes = ComplexUpset::upset_set_size(), # replaces set_size(...)
          min_size = 5
        )
      save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      library(UpSetR)
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
              paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }


  # Barplot of overlap counts and percentages
  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, ]
  if (nrow(ov2) > 0) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2$OverlapDEG,20))*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  }

  # Fisher enrichment of DEGs in modules
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))

  # Per-module analyses
  sel_modules <- fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)]
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) next

    for (ont in onts) {
      # Per-module nested folders
      m_root     <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea   <- mk(m_root, "Tables", "gsea")
      m_P_gsea   <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull<- mk(m_root, "Tables", "ora_full")
      m_P_orafull<- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]

      # Per-module GSEA
      if (length(m_gene_list) >= 10) {
        m_gse <- tryCatch(clusterProfiler::gseGO(
          geneList = sort(m_gene_list, decreasing = TRUE),
          ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) NULL)
        if (!is.null(m_gse) && nrow(m_gse@result) > 0) {
          save_table(m_gse@result, m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          save_plot_dir(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust", showCategory = 10, label_format = 30),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      }

      # ORA on DEG-intersection
      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) NULL)
        if (!is.null(m_ora) && nrow(m_ora@result) > 0) {
          save_table(m_ora@result, m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
          save_plot_dir(enrichplot::emapplot(enrichplot::pairwise_termsim(m_ora), showCategory = 10),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_emap_", file_tag, ".svg"))
          save_plot_dir(enrichplot::cnetplot(m_ora, showCategory = 5),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_cnet_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "insufficient DEGs for ORA"), m_T_oradeg,
                   paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on full module gene set (additional)
      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) NULL)
        if (!is.null(m_ora_all) && nrow(m_ora_all@result) > 0) {
          save_table(m_ora_all@result, m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      }
    } # ont
  } # module
} # comparison






























############################################################
# clusterProfiler + WGCNA modules (mouse)
# IDs: UniProt accessions already in input (no mapping step)
# Adds: robust guards, diagnostics, run logs, README, tie ranks
# Outputs are split into Tables vs Plots with explicit status files
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","vroom", "UniProt.ws", "BiocParallel",
    "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges", "SummarizedExperiment",
    "Biobase", "GO.db", "UpSetR", "ComplexUpset"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Directory + I/O helpers ----------
mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP   = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                  P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF   = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                  P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC   = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                  P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                  P = mk(results_dir, "10_global","KEGG","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots"))
    ),
    pathview = mk(results_dir, "90_pathview")
  )
}
save_table <- function(df, dirpath, fname) {
  readr::write_csv(df, file.path(dirpath, fname))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}
# Simple run logger per comparison file
log_line <- function(paths, tag, text) {
  lf <- file.path(paths$meta, paste0("runlog_", tag, ".txt"))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text), file = lf, append = TRUE)
}
# Stable ranking helper for tied p-values
stable_rank <- function(df) {
  if (is.null(df) || nrow(df)==0) return(df)
  # accept common p-value columns
  if (!("pvalue" %in% names(df)) && ("p.adjust" %in% names(df))) df$pvalue <- df$p.adjust
  if (!("p.adjust" %in% names(df)) && ("pvalue" %in% names(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df$rank_p <- rank(df$pvalue, ties.method = "min")
  df$rank_padj <- rank(df$p.adjust, ties.method = "min")
  df[order(df$p.adjust, df$pvalue, df$rank_padj, df$rank_p, decreasing = FALSE), ]
}

# Check if enrichment result has terms (min_n = 1 by default)
has_terms <- function(x, min_n = 1) {
  tryCatch({
    df <- as.data.frame(x)
    !is.null(df) && nrow(df) >= min_n
  }, error = function(e) FALSE)
}

## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text  = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

## ---------- Load WGCNA modules (already UniProt in Gene) ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "") |>
  dplyr::mutate(Gene = as.character(Gene))

# Module gene sets (UniProt accessions)
modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(x))
module_sizes <- tibble::tibble(
  Module = names(modules_list),
  N_accessions = vapply(modules_list, length, integer(1))
)

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Inputs ----------
comparison <- "neuron-regionLayerBaseline"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Meta, README, module sizes
  capture.output(sessionInfo(), file = file.path(paths$meta, "sessionInfo.txt"))
  save_table(module_sizes, paths$modules$root, "module_sizes.csv")
  writeLines(c(
    paste0("data_path: ", data_path),
    paste0("context: ", context),
    paste0("timestamp: ", Sys.time()),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (via mapIds)",
    "DEG threshold for ORA: |log2fc| > 1"
  ), con = file.path(paths$meta, "params.txt"))

  readme_lines <- c(
    "# Analysis README",
    paste0("Comparison: ", comparison),
    paste0("Context: ", context),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (mapIds from org.Mm.eg.db, multiVals='first')",
    "DEG threshold for ORA: |log2fc| > 1",
    "Ontologies default: BP",
    "Folders:",
    "  10_global/GO_*: GSEA and ORA across full universe",
    "  10_global/KEGG: gseKEGG results and mapping diagnostics",
    "  20_modules/: overlap, Fisher, per-module GO (GSEA/ORA)",
    "  90_pathview/: KEGG pathway renderings and status logs",
    "Check 00_meta/runlog_*.txt for stepwise notes."
  )
  writeLines(readme_lines, con = file.path(results_dir, "README.txt"))
  log_line(paths, file_tag, "Initialized results folders and README.")

  ## ----- Load comparison (assume UniProt accessions in col1) -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  gene_list_tbl <- df |>
    dplyr::select(id, log2fc) |>
    dplyr::arrange(dplyr::desc(log2fc))
  gene_list <- gene_list_tbl$log2fc
  names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)

  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  # Per-file diagnostics
  save_table(tibble::tibble(
    n_rows = nrow(df),
    n_universe = length(universe_genes),
    n_deg_abs1 = length(top_genes),
    deg_threshold = 1
  ), paths$meta, paste0("summary_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("Universe=%d; DEGs(|log2fc|>1)=%d", length(universe_genes), length(top_genes)))

  organism <- "org.Mm.eg.db"
  onts <- c("BP") # extend if needed

   ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    gse <- tryCatch(clusterProfiler::gseGO(
      geneList = gene_list, ont = ont, keyType = "UNIPROT",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
      by = "fgsea", nPermSimple = 10000, eps = 0,
      verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
    ), error = function(e) { log_line(paths, file_tag, paste("gseGO error:", e$message)); NULL })

    if (!is.null(gse)) {
      res_stable <- stable_rank(gse@result)
      save_table(res_stable, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))
      save_plot_dir(enrichplot::emapplot(enrichplot::pairwise_termsim(gse), showCategory = 10), dstP,
                    paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
      save_plot_dir(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list), dstP,
                    paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))
      save_plot_dir(enrichplot::ridgeplot(gse) + ggplot2::labs(x = "Enrichment Distribution"), dstP,
                    paste0("gseGO_", ont, "_ridge_", file_tag, ".svg"))
      if (nrow(gse@result) > 0) {
        save_plot_dir(enrichplot::gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1), dstP,
                      paste0("gseGO_", ont, "_gseaplot_", file_tag, ".svg"))
      }
      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) { log_line(paths, file_tag, paste("simplify error:", e$message)); NULL })
      if (!is.null(gse2)) {
        save_table(stable_rank(gse2@result), dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
      }
    } else {
      save_table(tibble::tibble(note = "gseGO returned NULL"), dstT, paste0("gseGO_", ont, "_", file_tag, "_empty.csv"))
    }

    if (length(top_genes) >= 10) {
      ora <- tryCatch(clusterProfiler::enrichGO(
        gene = top_genes, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
      ), error = function(e) { log_line(paths, file_tag, paste("ORA error:", e$message)); NULL })
      if (!is.null(ora)) {
        save_table(stable_rank(ora@result), dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        save_plot_dir(enrichplot::emapplot(enrichplot::pairwise_termsim(ora), showCategory = 10), dstP,
                      paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
        save_plot_dir(enrichplot::cnetplot(ora, showCategory = 5), dstP,
                      paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
      } else {
        save_table(tibble::tibble(note = "ORA returned NULL"), dstT, paste0("ORA_", ont, "_", file_tag, "_empty.csv"))
      }
    } else {
      save_table(tibble::tibble(note = "Global ORA skipped: <10 DEGs"),
                 dstT, paste0("ORA_", ont, "_", file_tag, "_skipped.csv"))
      log_line(paths, file_tag, paste("Global ORA skipped for", ont, "due to insufficient DEGs."))
    }
  }

  ## ----- GLOBAL: KEGG (mapIds to 1:1 ENTREZ) -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  keep <- !is.na(entrez_map)
  kegg_gene_list <- gene_list[keep]
  names(kegg_gene_list) <- unname(entrez_map[keep])
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

  # KEGG mapping diagnostics
  kmap_tbl <- tibble::tibble(
    id = universe_genes,
    ENTREZID = unname(entrez_map[universe_genes]),
    in_kegg = !is.na(entrez_map[universe_genes])
  )
  save_table(kmap_tbl, paths$glob$KEGG$T, paste0("kegg_mapping_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("KEGG mapIds: %d/%d mapped (%.1f%%)",
                                    sum(kmap_tbl$in_kegg), nrow(kmap_tbl), 100*mean(kmap_tbl$in_kegg)))

  kk2 <- NULL
  if (sum(kmap_tbl$in_kegg) >= 10) {
    kk2 <- tryCatch(clusterProfiler::gseKEGG(
      geneList = kegg_gene_list, organism = "mmu",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
      keyType = "ncbi-geneid"
    ), error = function(e) { log_line(paths, file_tag, paste("gseKEGG error:", e$message)); NULL })
  } else {
    log_line(paths, file_tag, "gseKEGG skipped: <10 mapped ENTREZIDs.")
  }

  if (!is.null(kk2) && nrow(kk2@result) > 0) {
    save_table(stable_rank(kk2@result), paths$glob$KEGG$T,
               paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = 10, label_format = 30),
                  paths$glob$KEGG$P, paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    save_plot_dir(enrichplot::emapplot(enrichplot::pairwise_termsim(kk2), showCategory = 10),
                  paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
    save_plot_dir(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                  paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))
    save_plot_dir(enrichplot::ridgeplot(kk2) + ggplot2::labs(x = "Enrichment distribution"),
                  paths$glob$KEGG$P, paste0("gseKEGG_ridge_", file_tag, ".svg"))
    if (nrow(kk2@result) > 0) {
      save_plot_dir(enrichplot::gseaplot(kk2, by = "all",
                                         title = kk2@result$Description[1],
                                         geneSetID = 1),
                    paths$glob$KEGG$P, paste0("gseKEGG_gseaplot_", file_tag, ".svg"))
    }
  } else {
    save_table(tibble::tibble(note = "No KEGG enrichment results (gseKEGG null or empty)."),
               paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, "_empty.csv"))
    log_line(paths, file_tag, "gseKEGG produced no results.")
  }

  ## ----- PATHVIEW with validated KEGG IDs + status log -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && nrow(kk2@result) > 0) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  save_table(tibble::tibble(pv_ids = pv_ids), paths$pathview, paste0("pathview_ids_", file_tag, ".csv"))

  if (length(kegg_gene_list) > 0 && all(is.finite(as.numeric(kegg_gene_list)))) {
    writeLines(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
    log_line(paths, file_tag, sprintf("Pathview attempted for %d pathways.", length(pv_ids)))
  } else {
    writeLines(c(
      "No valid ENTREZ-mapped genes for Pathview or gene list non-numeric.",
      paste0("length(kegg_gene_list)=", length(kegg_gene_list))
    ), con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Pathview skipped: invalid kegg_gene_list.")
  }

  ## ----- MODULE-AWARE: overlaps, UpSet, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")

  # UpSet-style overlap (top modules by DEG overlap) with guards
  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), ]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids
  save_table(mat, paths$modules$overlap$T, paste0("upset_input_", file_tag, ".csv"))

  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))
  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      library(ComplexUpset)
      p_up <- ggplot(mat, aes(x = NULL)) +
        ComplexUpset::upset(
          intersect = colnames(mat)[colnames(mat) != "id"],
          base_annotations = list(
            'Intersection size' = ComplexUpset::intersection_size(),
            'Set size' = ComplexUpset::set_size()
          ),
          min_size = 5
        )
      save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      library(UpSetR)
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
               paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }

  # Barplot of overlap counts and percentages
  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, ]
  if (nrow(ov2) > 0) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2$OverlapDEG,20))*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No non-zero overlaps to plot.", con = file.path(paths$modules$overlap$P, paste0("overlap_status_", file_tag, ".txt")))
  }

  # Fisher enrichment of DEGs in modules + plot
  fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
    deg_set <- intersect(unique(deg_genes), universe_genes)
    bg_set  <- setdiff(universe_genes, deg_set)
    out <- lapply(names(modules_list), function(m) {
      mod_genes <- intersect(modules_list[[m]], universe_genes)
      a <- length(intersect(mod_genes, deg_set))
      b <- length(setdiff(mod_genes, deg_set))
      c <- length(setdiff(deg_set, mod_genes))
      d <- length(setdiff(bg_set, mod_genes))
      mat <- matrix(c(a,b,c,d), nrow=2)
      ft <- fisher.test(mat, alternative = "greater")
      data.frame(
        Module = m, InModule_DE = a, InModule_NotDE = b,
        NotInModule_DE = c, NotInModule_NotDE = d,
        OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    res <- dplyr::bind_rows(out)
    res$Padj <- p.adjust(res$Pvalue, method = "BH")
    dplyr::arrange(res, Padj, Pvalue)
  }
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))

  if (nrow(fisher_tbl) > 0) {
    ft_plot <- ggplot2::ggplot(head(fisher_tbl[order(fisher_tbl$Padj, fisher_tbl$Pvalue), ], 20),
                               ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                            y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#6a51a3") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Fisher enrichment (DEGs in modules) — ", context),
                    x = "Module", y = "-log10(Padj)")
    save_plot_dir(ft_plot, paths$modules$fisher$P, paste0("fisher_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No Fisher rows to plot.", con = file.path(paths$modules$fisher$P, paste0("fisher_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Fisher plot skipped: empty table.")
  }

    # Per-module analyses
  sel_modules <- fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)]
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) {
      log_line(paths, file_tag, paste("Module", m, "skipped: <10 genes in universe."))
      next
    }

    for (ont in onts) {
      # Per-module nested folders
      m_root     <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea   <- mk(m_root, "Tables", "gsea")
      m_P_gsea   <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull<- mk(m_root, "Tables", "ora_full")
      m_P_orafull<- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]

      # Per-module GSEA
      if (length(m_gene_list) >= 10) {
        m_gse <- tryCatch(clusterProfiler::gseGO(
          geneList = sort(m_gene_list, decreasing = TRUE),
          ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "gseGO error:", e$message)); NULL })
        if (!is.null(m_gse) && nrow(m_gse@result) > 0) {
          save_table(stable_rank(m_gse@result), m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          save_plot_dir(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust", showCategory = 10, label_format = 30),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "Module GSEA skipped: <10 genes."), m_T_gsea,
                   paste0("gseGO_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on DEG-intersection
      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA DEG error:", e$message)); NULL })
        if (!is.null(m_ora) && nrow(m_ora@result) > 0) {
          save_table(stable_rank(m_ora@result), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
          save_plot_dir(enrichplot::emapplot(enrichplot::pairwise_termsim(m_ora), showCategory = 10),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_emap_", file_tag, ".svg"))
          save_plot_dir(enrichplot::cnetplot(m_ora, showCategory = 5),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_cnet_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(
          note = "insufficient DEGs for ORA",
          n_module_genes_in_universe = length(m_genes),
          n_deg_in_module = length(intersect(top_genes, m_genes)),
          deg_threshold = 1
        ), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on full module gene set
      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA full error:", e$message)); NULL })
        if (!is.null(m_ora_all) && nrow(m_ora_all@result) > 0) {
          save_table(stable_rank(m_ora_all@result), m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "ORA full module skipped: <5 genes."),
                   m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }
    } # ont
  } # module
} # for each comparison file






























############################################################
# clusterProfiler + WGCNA modules (mouse)
# IDs: UniProt accessions already in input (no mapping step)
# Robust guards for empty results, diagnostics, run logs, README
# ORA/GSEA tables gain stable rank columns to handle padj ties
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","vroom", "UniProt.ws", "BiocParallel",
    "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges", "SummarizedExperiment",
    "Biobase", "GO.db", "UpSetR", "ComplexUpset"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Directory + I/O helpers ----------
mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP   = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                  P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF   = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                  P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC   = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                  P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                  P = mk(results_dir, "10_global","KEGG","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots"))
    ),
    pathview = mk(results_dir, "90_pathview")
  )
}
save_table <- function(df, dirpath, fname) {
  readr::write_csv(df, file.path(dirpath, fname))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}
# Simple run logger per comparison file
log_line <- function(paths, tag, text) {
  lf <- file.path(paths$meta, paste0("runlog_", tag, ".txt"))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text), file = lf, append = TRUE)
}
# Stable ranking helper for tied p-values
stable_rank <- function(df) {
  if (is.null(df) || nrow(df)==0) return(df)
  if (!("pvalue" %in% names(df)) && ("p.adjust" %in% names(df))) df$pvalue <- df$p.adjust
  if (!("p.adjust" %in% names(df)) && ("pvalue" %in% names(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df$rank_p <- rank(df$pvalue, ties.method = "min")
  df$rank_padj <- rank(df$p.adjust, ties.method = "min")
  df[order(df$p.adjust, df$pvalue, df$rank_padj, df$rank_p, decreasing = FALSE), ]
}
# Result presence check
has_terms <- function(x, min_n = 1) {
  tryCatch({
    df <- as.data.frame(x)
    !is.null(df) && nrow(df) >= min_n
  }, error = function(e) FALSE)
}

## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text  = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

## ---------- Load WGCNA modules (already UniProt in Gene) ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "") |>
  dplyr::mutate(Gene = as.character(Gene))

# Module gene sets (UniProt accessions)
modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(x))
module_sizes <- tibble::tibble(
  Module = names(modules_list),
  N_accessions = vapply(modules_list, length, integer(1))
)

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Inputs ----------
#comparison <- "neuron-regionLayerBaseline"
comparison <- "neuron-phenotypeWithinUnit"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Meta, README, module sizes
  capture.output(sessionInfo(), file = file.path(paths$meta, "sessionInfo.txt"))
  save_table(module_sizes, paths$modules$root, "module_sizes.csv")
  writeLines(c(
    paste0("data_path: ", data_path),
    paste0("context: ", context),
    paste0("timestamp: ", Sys.time()),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (via mapIds)",
    "DEG threshold for ORA: |log2fc| > 1"
  ), con = file.path(paths$meta, "params.txt"))

  readme_lines <- c(
    "# Analysis README",
    paste0("Comparison: ", comparison),
    paste0("Context: ", context),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (mapIds from org.Mm.eg.db, multiVals='first')",
    "DEG threshold for ORA: |log2fc| > 1",
    "Ontologies default: BP",
    "Folders:",
    "  10_global/GO_*: GSEA and ORA across full universe",
    "  10_global/KEGG: gseKEGG results and mapping diagnostics",
    "  20_modules/: overlap, Fisher, per-module GO (GSEA/ORA)",
    "  90_pathview/: KEGG pathway renderings and status logs",
    "Check 00_meta/runlog_*.txt for stepwise notes."
  )
  writeLines(readme_lines, con = file.path(results_dir, "README.txt"))
  log_line(paths, file_tag, "Initialized results folders and README.")

  ## ----- Load comparison (assume UniProt accessions in col1) -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  gene_list_tbl <- df |>
    dplyr::select(id, log2fc) |>
    dplyr::arrange(dplyr::desc(log2fc))
  gene_list <- gene_list_tbl$log2fc
  names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)

  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  # Per-file diagnostics
  save_table(tibble::tibble(
    n_rows = nrow(df),
    n_universe = length(universe_genes),
    n_deg_abs1 = length(top_genes),
    deg_threshold = 1
  ), paths$meta, paste0("summary_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("Universe=%d; DEGs(|log2fc|>1)=%d", length(universe_genes), length(top_genes)))

  organism <- "org.Mm.eg.db"
  onts <- c("BP", "MF", "CC") # extend if needed

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    gse <- tryCatch(clusterProfiler::gseGO(
      geneList = gene_list, ont = ont, keyType = "UNIPROT",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
      by = "fgsea", nPermSimple = 10000, eps = 0,
      verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
    ), error = function(e) { log_line(paths, file_tag, paste("gseGO error:", e$message)); NULL })

    if (!is.null(gse) && has_terms(gse, 1)) {
      res_stable <- stable_rank(gse@result)
      save_table(res_stable, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))

      if (nrow(gse@result) >= 2) {
        ts_gse <- tryCatch(enrichplot::pairwise_termsim(gse), error = function(e) NULL)
        if (!is.null(ts_gse) && "compareClusterResult" %in% slotNames(ts_gse) &&
            nrow(ts_gse@compareClusterResult) > 0) {
          save_plot_dir(enrichplot::emapplot(ts_gse, showCategory = 10), dstP,
                        paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
        } else {
          writeLines("emapplot skipped: no enriched term after similarity calc.",
                     file.path(dstP, paste0("gseGO_", ont, "_emap_", file_tag, "_skipped.txt")))
        }
      }

      p_cnet <- tryCatch(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                         error = function(e) NULL)
      if (!is.null(p_cnet)) save_plot_dir(p_cnet, dstP, paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))

      rp <- tryCatch(enrichplot::ridgeplot(gse) + ggplot2::labs(x = "Enrichment Distribution"),
                     error = function(e) NULL)
      if (!is.null(rp)) save_plot_dir(rp, dstP, paste0("gseGO_", ont, "_ridge_", file_tag, ".svg"))

      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) { log_line(paths, file_tag, paste("simplify error:", e$message)); NULL })
      if (!is.null(gse2) && has_terms(gse2, 1)) {
        save_table(stable_rank(gse2@result), dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
      }
    } else {
      save_table(tibble::tibble(note = "gseGO empty or NULL"), dstT, paste0("gseGO_", ont, "_", file_tag, "_empty.csv"))
      writeLines("All plots skipped: no enriched term.", file.path(dstP, paste0("gseGO_", ont, "_", file_tag, "_skipped.txt")))
    }

    if (length(top_genes) >= 10) {
      ora <- tryCatch(clusterProfiler::enrichGO(
        gene = top_genes, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
      ), error = function(e) { log_line(paths, file_tag, paste("ORA error:", e$message)); NULL })
      if (!is.null(ora) && has_terms(ora, 1)) {
        save_table(stable_rank(ora@result), dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        if (nrow(ora@result) >= 2) {
          ts_ora <- tryCatch(enrichplot::pairwise_termsim(ora), error = function(e) NULL)
          if (!is.null(ts_ora) && "compareClusterResult" %in% slotNames(ts_ora) &&
              nrow(ts_ora@compareClusterResult) > 0) {
            save_plot_dir(enrichplot::emapplot(ts_ora, showCategory = 10), dstP,
                          paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
          } else {
            writeLines("emapplot skipped: no enriched term after similarity calc.",
                       file.path(dstP, paste0("ORA_", ont, "_emap_", file_tag, "_skipped.txt")))
          }
        }
        p_cnet2 <- tryCatch(enrichplot::cnetplot(ora, showCategory = 5), error = function(e) NULL)
        if (!is.null(p_cnet2)) save_plot_dir(p_cnet2, dstP, paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
      } else {
        save_table(tibble::tibble(note = "ORA returned NULL or empty"), dstT, paste0("ORA_", ont, "_", file_tag, "_empty.csv"))
      }
    } else {
      save_table(tibble::tibble(note = "Global ORA skipped: <10 DEGs"),
                 dstT, paste0("ORA_", ont, "_", file_tag, "_skipped.csv"))
      log_line(paths, file_tag, paste("Global ORA skipped for", ont, "due to insufficient DEGs."))
    }
  }

  ## ----- GLOBAL: KEGG (mapIds to 1:1 ENTREZ) -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  keep <- !is.na(entrez_map)
  kegg_gene_list <- gene_list[keep]
  names(kegg_gene_list) <- unname(entrez_map[keep])
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

  # KEGG mapping diagnostics
  kmap_tbl <- tibble::tibble(
    id = universe_genes,
    ENTREZID = unname(entrez_map[universe_genes]),
    in_kegg = !is.na(entrez_map[universe_genes])
  )
  save_table(kmap_tbl, paths$glob$KEGG$T, paste0("kegg_mapping_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("KEGG mapIds: %d/%d mapped (%.1f%%)",
                                    sum(kmap_tbl$in_kegg), nrow(kmap_tbl), 100*mean(kmap_tbl$in_kegg)))

  kk2 <- NULL
  if (sum(kmap_tbl$in_kegg) >= 10) {
    kk2 <- tryCatch(clusterProfiler::gseKEGG(
      geneList = kegg_gene_list, organism = "mmu",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
      keyType = "ncbi-geneid"
    ), error = function(e) { log_line(paths, file_tag, paste("gseKEGG error:", e$message)); NULL })
  } else {
    log_line(paths, file_tag, "gseKEGG skipped: <10 mapped ENTREZIDs.")
  }

  if (!is.null(kk2) && has_terms(kk2, 1)) {
    save_table(stable_rank(kk2@result), paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = 10, label_format = 30),
                  paths$glob$KEGG$P, paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    if (nrow(kk2@result) >= 2) {
      ts_kk <- tryCatch(enrichplot::pairwise_termsim(kk2), error = function(e) NULL)
      if (!is.null(ts_kk) && "compareClusterResult" %in% slotNames(ts_kk) &&
          nrow(ts_kk@compareClusterResult) > 0) {
        save_plot_dir(enrichplot::emapplot(ts_kk, showCategory = 10),
                      paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
      } else {
        writeLines("emapplot skipped: no enriched term after similarity calc.",
                   file.path(paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, "_skipped.txt")))
      }
    }
    p_cnet3 <- tryCatch(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                        error = function(e) NULL)
    if (!is.null(p_cnet3)) save_plot_dir(p_cnet3, paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))
    rp2 <- tryCatch(enrichplot::ridgeplot(kk2) + ggplot2::labs(x = "Enrichment distribution"),
                    error = function(e) NULL)
    if (!is.null(rp2)) save_plot_dir(rp2, paths$glob$KEGG$P, paste0("gseKEGG_ridge_", file_tag, ".svg"))
    if (nrow(kk2@result) > 0) {
      gp <- tryCatch(enrichplot::gseaplot(kk2, by = "all",
                                          title = kk2@result$Description[1], geneSetID = 1),
                     error = function(e) NULL)
      if (!is.null(gp)) save_plot_dir(gp, paths$glob$KEGG$P, paste0("gseKEGG_gseaplot_", file_tag, ".svg"))
    }
  } else {
    save_table(tibble::tibble(note = "No KEGG enrichment results (gseKEGG null or empty)."),
               paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, "_empty.csv"))
    log_line(paths, file_tag, "gseKEGG produced no results.")
  }

  ## ----- PATHVIEW with validated KEGG IDs + status log -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && has_terms(kk2, 1)) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  save_table(tibble::tibble(pv_ids = pv_ids), paths$pathview, paste0("pathview_ids_", file_tag, ".csv"))

  if (length(kegg_gene_list) > 0 && all(is.finite(as.numeric(kegg_gene_list)))) {
    writeLines(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
    log_line(paths, file_tag, sprintf("Pathview attempted for %d pathways.", length(pv_ids)))
  } else {
    writeLines(c(
      "No valid ENTREZ-mapped genes for Pathview or gene list non-numeric.",
      paste0("length(kegg_gene_list)=", length(kegg_gene_list))
    ), con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Pathview skipped: invalid kegg_gene_list.")
  }

  ## ----- MODULE-AWARE: overlaps, UpSet, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")

  # UpSet-style overlap (top modules by DEG overlap) with guards
  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), ]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids
  save_table(mat, paths$modules$overlap$T, paste0("upset_input_", file_tag, ".csv"))

  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))
  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      library(ComplexUpset)
      p_up <- ggplot(mat, aes(x = NULL)) +
        ComplexUpset::upset(
          intersect = colnames(mat)[colnames(mat) != "id"],
          base_annotations = list(
            'Intersection size' = ComplexUpset::intersection_size(),
            'Set size' = ComplexUpset::set_size()
          ),
          min_size = 5
        )
      save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      library(UpSetR)
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
               paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }

  # Barplot of overlap counts and percentages
  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, ]
  if (nrow(ov2) > 0) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2, 20)$OverlapDEG)*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No non-zero overlaps to plot.", con = file.path(paths$modules$overlap$P, paste0("overlap_status_", file_tag, ".txt")))
  }

  # Fisher enrichment of DEGs in modules + plot
  fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
    deg_set <- intersect(unique(deg_genes), universe_genes)
    bg_set  <- setdiff(universe_genes, deg_set)
    out <- lapply(names(modules_list), function(m) {
      mod_genes <- intersect(modules_list[[m]], universe_genes)
      a <- length(intersect(mod_genes, deg_set))
      b <- length(setdiff(mod_genes, deg_set))
      c <- length(setdiff(deg_set, mod_genes))
      d <- length(setdiff(bg_set, mod_genes))
      mat <- matrix(c(a,b,c,d), nrow=2)
      ft <- fisher.test(mat, alternative = "greater")
      data.frame(
        Module = m, InModule_DE = a, InModule_NotDE = b,
        NotInModule_DE = c, NotInModule_NotDE = d,
        OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    res <- dplyr::bind_rows(out)
    res$Padj <- p.adjust(res$Pvalue, method = "BH")
    dplyr::arrange(res, Padj, Pvalue)
  }
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))

  if (nrow(fisher_tbl) > 0) {
    ft_plot <- ggplot2::ggplot(head(fisher_tbl[order(fisher_tbl$Padj, fisher_tbl$Pvalue), ], 20),
                               ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                            y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#6a51a3") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Fisher enrichment (DEGs in modules) — ", context),
                    x = "Module", y = "-log10(Padj)")
    save_plot_dir(ft_plot, paths$modules$fisher$P, paste0("fisher_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No Fisher rows to plot.", con = file.path(paths$modules$fisher$P, paste0("fisher_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Fisher plot skipped: empty table.")
  }

  ## ----- Per-module analyses -----
  sel_modules <- fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)]
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) {
      log_line(paths, file_tag, paste("Module", m, "skipped: <10 genes in universe."))
      next
    }

    for (ont in onts) {
      # Per-module nested folders
      m_root     <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea   <- mk(m_root, "Tables", "gsea")
      m_P_gsea   <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull<- mk(m_root, "Tables", "ora_full")
      m_P_orafull<- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]

      # Per-module GSEA
      if (length(m_gene_list) >= 10) {
        m_gse <- tryCatch(clusterProfiler::gseGO(
          geneList = sort(m_gene_list, decreasing = TRUE),
          ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "gseGO error:", e$message)); NULL })
        if (!is.null(m_gse) && has_terms(m_gse, 1)) {
          save_table(stable_rank(m_gse@result), m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          p_dot <- tryCatch(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust",
                                                showCategory = 10, label_format = 30), error = function(e) NULL)
          if (!is.null(p_dot)) save_plot_dir(p_dot, m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "Module GSEA skipped: <10 genes."), m_T_gsea,
                   paste0("gseGO_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on DEG-intersection
      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA DEG error:", e$message)); NULL })
        if (!is.null(m_ora) && has_terms(m_ora, 1)) {
          save_table(stable_rank(m_ora@result), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
          if (nrow(m_ora@result) >= 2) {
            ts_mora <- tryCatch(enrichplot::pairwise_termsim(m_ora), error = function(e) NULL)
            if (!is.null(ts_mora) && "compareClusterResult" %in% slotNames(ts_mora) &&
                nrow(ts_mora@compareClusterResult) > 0) {
              save_plot_dir(enrichplot::emapplot(ts_mora, showCategory = 10),
                            m_P_oradeg, paste0("ORA_", ont, "_", m, "_emap_", file_tag, ".svg"))
            } else {
              writeLines("emapplot skipped: no enriched term after similarity calc.",
                         file.path(m_P_oradeg, paste0("ORA_", ont, "_", m, "_emap_", file_tag, "_skipped.txt")))
            }
          }
          p_cnet_mora <- tryCatch(enrichplot::cnetplot(m_ora, showCategory = 5), error = function(e) NULL)
          if (!is.null(p_cnet_mora)) save_plot_dir(p_cnet_mora, m_P_oradeg, paste0("ORA_", ont, "_", m, "_cnet_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(
          note = "insufficient DEGs for ORA",
          n_module_genes_in_universe = length(m_genes),
          n_deg_in_module = length(intersect(top_genes, m_genes)),
          deg_threshold = 1
        ), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on full module gene set
      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA full error:", e$message)); NULL })
        if (!is.null(m_ora_all) && has_terms(m_ora_all, 1)) {
          save_table(stable_rank(m_ora_all@result), m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "ORA full module skipped: <5 genes."),
                   m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }
    } # ont
  } # module
} # for each comparison file






































############################################################
# clusterProfiler + WGCNA modules (mouse)
# IDs: UniProt accessions already in input (no mapping step)
# Deterministic tie-breaking for fgsea (order by id), no value jitter
# ComplexUpset set sizes panel disabled for version-agnostic plots
# Robust guards, diagnostics, run logs, README
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","vroom", "UniProt.ws", "BiocParallel",
    "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges", "SummarizedExperiment",
    "Biobase", "GO.db", "UpSetR", "ComplexUpset"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Directory + I/O helpers ----------
mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP   = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                  P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF   = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                  P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC   = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                  P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                  P = mk(results_dir, "10_global","KEGG","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots"))
    ),
    pathview = mk(results_dir, "90_pathview")
  )
}
save_table <- function(df, dirpath, fname) {
  readr::write_csv(df, file.path(dirpath, fname))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}
log_line <- function(paths, tag, text) {
  lf <- file.path(paths$meta, paste0("runlog_", tag, ".txt"))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text), file = lf, append = TRUE)
}
stable_rank <- function(df) {
  if (is.null(df) || nrow(df)==0) return(df)
  if (!("pvalue" %in% names(df)) && ("p.adjust" %in% names(df))) df$pvalue <- df$p.adjust
  if (!("p.adjust" %in% names(df)) && ("pvalue" %in% names(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df$rank_p <- rank(df$pvalue, ties.method = "min")
  df$rank_padj <- rank(df$p.adjust, ties.method = "min")
  df[order(df$p.adjust, df$pvalue, df$rank_padj, df$rank_p, decreasing = FALSE), ]
}
has_terms <- function(x, min_n = 1) {
  tryCatch({
    df <- as.data.frame(x)
    !is.null(df) && nrow(df) >= min_n
  }, error = function(e) FALSE)
}
muffle_fgsea_ties <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("There are ties in the preranked stats", conditionMessage(m))) {
      invokeRestart("muffleMessage")
    }
  })
}

## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text  = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

## ---------- Load WGCNA modules (already UniProt in Gene) ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "") |>
  dplyr::mutate(Gene = as.character(Gene))

modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(x))
module_sizes <- tibble::tibble(
  Module = names(modules_list),
  N_accessions = vapply(modules_list, length, integer(1))
)

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Inputs ----------
#comparison <- "neuron-regionLayerBaseline"
comparison <- "neuron-phenotypeWithinUnit"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Meta, README, module sizes
  capture.output(sessionInfo(), file = file.path(paths$meta, "sessionInfo.txt"))
  save_table(module_sizes, paths$modules$root, "module_sizes.csv")
  writeLines(c(
    paste0("data_path: ", data_path),
    paste0("context: ", context),
    paste0("timestamp: ", Sys.time()),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (via mapIds)",
    "DEG threshold for ORA: |log2fc| > 1",
    "GSEA ties: deterministically broken by sorting ties by id (values unchanged)"
  ), con = file.path(paths$meta, "params.txt"))

  readme_lines <- c(
    "# Analysis README",
    paste0("Comparison: ", comparison),
    paste0("Context: ", context),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (mapIds from org.Mm.eg.db, multiVals='first')",
    "DEG threshold for ORA: |log2fc| > 1",
    "Ontologies default: BP",
    "ComplexUpset: set sizes panel disabled for portability",
    "GSEA ties: sorted by id for deterministic order; values are NOT jittered",
    "Folders:",
    "  10_global/GO_*: GSEA and ORA across full universe",
    "  10_global/KEGG: gseKEGG results and mapping diagnostics",
    "  20_modules/: overlap, Fisher, per-module GO (GSEA/ORA)",
    "  90_pathview/: KEGG pathway renderings and status logs",
    "Check 00_meta/runlog_*.txt for stepwise notes."
  )
  writeLines(readme_lines, con = file.path(results_dir, "README.txt"))
  log_line(paths, file_tag, "Initialized results folders and README.")

  ## ----- Load comparison (assume UniProt accessions in col1) -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # Deterministic tie-breaking by id, values unchanged
  gene_list_tbl <- df |>
    dplyr::select(id, log2fc) |>
    dplyr::arrange(dplyr::desc(log2fc), id)
  gene_list <- gene_list_tbl$log2fc
  names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)

  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  # Per-file diagnostics
  pct_tied <- 100 * mean(duplicated(gene_list))
  save_table(tibble::tibble(
    n_rows = nrow(df),
    n_universe = length(universe_genes),
    n_deg_abs1 = length(top_genes),
    deg_threshold = 1,
    pct_tied_stats = pct_tied
  ), paths$meta, paste0("summary_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("Universe=%d; DEGs(|log2fc|>1)=%d; tied stats=%.1f%%",
                                    length(universe_genes), length(top_genes), pct_tied))

  organism <- "org.Mm.eg.db"
  onts <- c("BP") # extend if needed

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
      geneList = gene_list, ont = ont, keyType = "UNIPROT",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
      by = "fgsea", nPermSimple = 10000, eps = 0,
      verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
    ), error = function(e) { log_line(paths, file_tag, paste("gseGO error:", e$message)); NULL }) )

    if (!is.null(gse) && has_terms(gse, 1)) {
      res_stable <- stable_rank(gse@result)
      res_stable$note <- "fgsea ties handled by deterministic id sort (values unchanged)"
      save_table(res_stable, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))

      if (nrow(gse@result) >= 2) {
        ts_gse <- tryCatch(enrichplot::pairwise_termsim(gse), error = function(e) NULL)
        if (!is.null(ts_gse) && "compareClusterResult" %in% slotNames(ts_gse) &&
            nrow(ts_gse@compareClusterResult) > 0) {
          save_plot_dir(enrichplot::emapplot(ts_gse, showCategory = 10), dstP,
                        paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
        } else {
          writeLines("emapplot skipped: no enriched term after similarity calc.",
                     file.path(dstP, paste0("gseGO_", ont, "_emap_", file_tag, "_skipped.txt")))
        }
      }

      p_cnet <- tryCatch(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                         error = function(e) NULL)
      if (!is.null(p_cnet)) save_plot_dir(p_cnet, dstP, paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))

      rp <- tryCatch(enrichplot::ridgeplot(gse) + ggplot2::labs(x = "Enrichment Distribution"),
                     error = function(e) NULL)
      if (!is.null(rp)) save_plot_dir(rp, dstP, paste0("gseGO_", ont, "_ridge_", file_tag, ".svg"))

      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) { log_line(paths, file_tag, paste("simplify error:", e$message)); NULL })
      if (!is.null(gse2) && has_terms(gse2, 1)) {
        save_table(stable_rank(gse2@result), dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
      }
    } else {
      save_table(tibble::tibble(note = "gseGO empty or NULL"), dstT, paste0("gseGO_", ont, "_", file_tag, "_empty.csv"))
      writeLines("All plots skipped: no enriched term.", file.path(dstP, paste0("gseGO_", ont, "_", file_tag, "_skipped.txt")))
    }

    if (length(top_genes) >= 10) {
      ora <- tryCatch(clusterProfiler::enrichGO(
        gene = top_genes, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
      ), error = function(e) { log_line(paths, file_tag, paste("ORA error:", e$message)); NULL })
      if (!is.null(ora) && has_terms(ora, 1)) {
        save_table(stable_rank(ora@result), dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        if (nrow(ora@result) >= 2) {
          ts_ora <- tryCatch(enrichplot::pairwise_termsim(ora), error = function(e) NULL)
          if (!is.null(ts_ora) && "compareClusterResult" %in% slotNames(ts_ora) &&
              nrow(ts_ora@compareClusterResult) > 0) {
            save_plot_dir(enrichplot::emapplot(ts_ora, showCategory = 10), dstP,
                          paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
          } else {
            writeLines("emapplot skipped: no enriched term after similarity calc.",
                       file.path(dstP, paste0("ORA_", ont, "_emap_", file_tag, "_skipped.txt")))
          }
        }
        p_cnet2 <- tryCatch(enrichplot::cnetplot(ora, showCategory = 5), error = function(e) NULL)
        if (!is.null(p_cnet2)) save_plot_dir(p_cnet2, dstP, paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
      } else {
        save_table(tibble::tibble(note = "ORA returned NULL or empty"), dstT, paste0("ORA_", ont, "_", file_tag, "_empty.csv"))
      }
    } else {
      save_table(tibble::tibble(note = "Global ORA skipped: <10 DEGs"),
                 dstT, paste0("ORA_", ont, "_", file_tag, "_skipped.csv"))
      log_line(paths, file_tag, paste("Global ORA skipped for", ont, "due to insufficient DEGs."))
    }
  }

  ## ----- GLOBAL: KEGG (mapIds to 1:1 ENTREZ) -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  keep <- !is.na(entrez_map)
  kegg_gene_list <- gene_list[keep]
  names(kegg_gene_list) <- unname(entrez_map[keep])
  kegg_gene_list <- kegg_gene_list[order(kegg_gene_list, decreasing = TRUE, method = "radix")]

  # KEGG mapping diagnostics
  kmap_tbl <- tibble::tibble(
    id = universe_genes,
    ENTREZID = unname(entrez_map[universe_genes]),
    in_kegg = !is.na(entrez_map[universe_genes])
  )
  save_table(kmap_tbl, paths$glob$KEGG$T, paste0("kegg_mapping_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("KEGG mapIds: %d/%d mapped (%.1f%%)",
                                    sum(kmap_tbl$in_kegg), nrow(kmap_tbl), 100*mean(kmap_tbl$in_kegg)))

  kk2 <- NULL
  if (sum(kmap_tbl$in_kegg) >= 10) {
    kk2 <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseKEGG(
      geneList = kegg_gene_list, organism = "mmu",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
      keyType = "ncbi-geneid"
    ), error = function(e) { log_line(paths, file_tag, paste("gseKEGG error:", e$message)); NULL }) )
  } else {
    log_line(paths, file_tag, "gseKEGG skipped: <10 mapped ENTREZIDs.")
  }

  if (!is.null(kk2) && has_terms(kk2, 1)) {
    save_table(stable_rank(kk2@result), paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = 10, label_format = 30),
                  paths$glob$KEGG$P, paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    if (nrow(kk2@result) >= 2) {
      ts_kk <- tryCatch(enrichplot::pairwise_termsim(kk2), error = function(e) NULL)
      if (!is.null(ts_kk) && "compareClusterResult" %in% slotNames(ts_kk) &&
          nrow(ts_kk@compareClusterResult) > 0) {
        save_plot_dir(enrichplot::emapplot(ts_kk, showCategory = 10),
                      paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
      } else {
        writeLines("emapplot skipped: no enriched term after similarity calc.",
                   file.path(paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, "_skipped.txt")))
      }
    }
    p_cnet3 <- tryCatch(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                        error = function(e) NULL)
    if (!is.null(p_cnet3)) save_plot_dir(p_cnet3, paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))
    rp2 <- tryCatch(enrichplot::ridgeplot(kk2) + ggplot2::labs(x = "Enrichment distribution"),
                    error = function(e) NULL)
    if (!is.null(rp2)) save_plot_dir(rp2, paths$glob$KEGG$P, paste0("gseKEGG_ridge_", file_tag, ".svg"))
    if (nrow(kk2@result) > 0) {
      gp <- tryCatch(enrichplot::gseaplot(kk2, by = "all",
                                          title = kk2@result$Description[1], geneSetID = 1),
                     error = function(e) NULL)
      if (!is.null(gp)) save_plot_dir(gp, paths$glob$KEGG$P, paste0("gseKEGG_gseaplot_", file_tag, ".svg"))
    }
  } else {
    save_table(tibble::tibble(note = "No KEGG enrichment results (gseKEGG null or empty)."),
               paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, "_empty.csv"))
    log_line(paths, file_tag, "gseKEGG produced no results.")
  }

  ## ----- PATHVIEW with validated KEGG IDs + status log -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && has_terms(kk2, 1)) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  save_table(tibble::tibble(pv_ids = pv_ids), paths$pathview, paste0("pathview_ids_", file_tag, ".csv"))

  if (length(kegg_gene_list) > 0 && all(is.finite(as.numeric(kegg_gene_list)))) {
    writeLines(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
    log_line(paths, file_tag, sprintf("Pathview attempted for %d pathways.", length(pv_ids)))
  } else {
    writeLines(c(
      "No valid ENTREZ-mapped genes for Pathview or gene list non-numeric.",
      paste0("length(kegg_gene_list)=", length(kegg_gene_list))
    ), con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Pathview skipped: invalid kegg_gene_list.")
  }

  ## ----- MODULE-AWARE: overlaps, UpSet, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")

  # UpSet-style overlap (top modules by DEG overlap) with version-agnostic settings
  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), ]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids
  save_table(mat, paths$modules$overlap$T, paste0("upset_input_", file_tag, ".csv"))

  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))
  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      library(ComplexUpset)
      p_up <- ggplot(mat, aes(x = NULL)) +
        ComplexUpset::upset(
          intersect = colnames(mat)[colnames(mat) != "id"],
          base_annotations = list(
            'Intersection size' = ComplexUpset::intersection_size()
          ),
          set_sizes = FALSE,   # disable for portability across versions
          min_size = 5
        )
      save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      library(UpSetR)
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
               paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }

  # Barplot of overlap counts and percentages
  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, ]
  if (nrow(ov2) > 0) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2, 20)$OverlapDEG)*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No non-zero overlaps to plot.", con = file.path(paths$modules$overlap$P, paste0("overlap_status_", file_tag, ".txt")))
  }

  # Fisher enrichment of DEGs in modules + plot
  fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
    deg_set <- intersect(unique(deg_genes), universe_genes)
    bg_set  <- setdiff(universe_genes, deg_set)
    out <- lapply(names(modules_list), function(m) {
      mod_genes <- intersect(modules_list[[m]], universe_genes)
      a <- length(intersect(mod_genes, deg_set))
      b <- length(setdiff(mod_genes, deg_set))
      c <- length(setdiff(deg_set, mod_genes))
      d <- length(setdiff(bg_set, mod_genes))
      mat <- matrix(c(a,b,c,d), nrow=2)
      ft <- fisher.test(mat, alternative = "greater")
      data.frame(
        Module = m, InModule_DE = a, InModule_NotDE = b,
        NotInModule_DE = c, NotInModule_NotDE = d,
        OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    res <- dplyr::bind_rows(out)
    res$Padj <- p.adjust(res$Pvalue, method = "BH")
    dplyr::arrange(res, Padj, Pvalue)
  }
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))

  if (nrow(fisher_tbl) > 0) {
    ft_plot <- ggplot2::ggplot(head(fisher_tbl[order(fisher_tbl$Padj, fisher_tbl$Pvalue), ], 20),
                               ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                            y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#6a51a3") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Fisher enrichment (DEGs in modules) — ", context),
                    x = "Module", y = "-log10(Padj)")
    save_plot_dir(ft_plot, paths$modules$fisher$P, paste0("fisher_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No Fisher rows to plot.", con = file.path(paths$modules$fisher$P, paste0("fisher_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Fisher plot skipped: empty table.")
  }

  ## ----- Per-module analyses -----
  sel_modules <- fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)]
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) {
      log_line(paths, file_tag, paste("Module", m, "skipped: <10 genes in universe."))
      next
    }

    for (ont in onts) {
      # Per-module nested folders
      m_root     <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea   <- mk(m_root, "Tables", "gsea")
      m_P_gsea   <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull<- mk(m_root, "Tables", "ora_full")
      m_P_orafull<- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]
      # Order by value desc, break ties by id (names) deterministically
      m_gene_list <- m_gene_list[order(m_gene_list, decreasing = TRUE, method = "radix")]
      if (length(m_gene_list) >= 10) {
        m_gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
          geneList = m_gene_list,
          ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "gseGO error:", e$message)); NULL }) )
        if (!is.null(m_gse) && has_terms(m_gse, 1)) {
          resm <- stable_rank(m_gse@result); resm$note <- "fgsea ties handled by deterministic id sort"
          save_table(resm, m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          p_dot <- tryCatch(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust",
                                                showCategory = 10, label_format = 30), error = function(e) NULL)
          if (!is.null(p_dot)) save_plot_dir(p_dot, m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "Module GSEA skipped: <10 genes."), m_T_gsea,
                   paste0("gseGO_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on DEG-intersection
      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA DEG error:", e$message)); NULL })
        if (!is.null(m_ora) && has_terms(m_ora, 1)) {
          save_table(stable_rank(m_ora@result), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
          if (nrow(m_ora@result) >= 2) {
            ts_mora <- tryCatch(enrichplot::pairwise_termsim(m_ora), error = function(e) NULL)
            if (!is.null(ts_mora) && "compareClusterResult" %in% slotNames(ts_mora) &&
                nrow(ts_mora@compareClusterResult) > 0) {
              save_plot_dir(enrichplot::emapplot(ts_mora, showCategory = 10),
                            m_P_oradeg, paste0("ORA_", ont, "_", m, "_emap_", file_tag, ".svg"))
            } else {
              writeLines("emapplot skipped: no enriched term after similarity calc.",
                         file.path(m_P_oradeg, paste0("ORA_", ont, "_", m, "_emap_", file_tag, "_skipped.txt")))
            }
          }
          p_cnet_mora <- tryCatch(enrichplot::cnetplot(m_ora, showCategory = 5), error = function(e) NULL)
          if (!is.null(p_cnet_mora)) save_plot_dir(p_cnet_mora, m_P_oradeg, paste0("ORA_", ont, "_", m, "_cnet_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
        } else {
        save_table(tibble::tibble(
          note = "insufficient DEGs for ORA",
          n_module_genes_in_universe = length(m_genes),
          n_deg_in_module = length(intersect(top_genes, m_genes)),
          deg_threshold = 1
        ), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on full module gene set
      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA full error:", e$message)); NULL })
        if (!is.null(m_ora_all) && has_terms(m_ora_all, 1)) {
          save_table(stable_rank(m_ora_all@result), m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "ORA full module skipped: <5 genes."),
                   m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }
    } # ont
  } # module
} # for each comparison file
































############################################################
# clusterProfiler + WGCNA modules (mouse)
# IDs: UniProt accessions already in input (no mapping step)
# Deterministic tie-breaking for fgsea (order by id), no value jitter
# ComplexUpset set sizes panel disabled for version-agnostic plots
# Robust guards, diagnostics, run logs, README
############################################################


## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","vroom", "UniProt.ws", "BiocParallel",
    "BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges", "SummarizedExperiment",
    "Biobase", "GO.db", "UpSetR", "ComplexUpset"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()


## ---------- Directory + I/O helpers ----------
mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP   = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                  P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF   = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                  P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC   = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                  P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                  P = mk(results_dir, "10_global","KEGG","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots"))
    ),
    pathview = mk(results_dir, "90_pathview")
  )
}
save_table <- function(df, dirpath, fname) {
  readr::write_csv(df, file.path(dirpath, fname))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}
log_line <- function(paths, tag, text) {
  lf <- file.path(paths$meta, paste0("runlog_", tag, ".txt"))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text), file = lf, append = TRUE)
}
stable_rank <- function(df) {
  if (is.null(df) || nrow(df)==0) return(df)
  if (!("pvalue" %in% names(df)) && ("p.adjust" %in% names(df))) df$pvalue <- df$p.adjust
  if (!("p.adjust" %in% names(df)) && ("pvalue" %in% names(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df$rank_p <- rank(df$pvalue, ties.method = "min")
  df$rank_padj <- rank(df$p.adjust, ties.method = "min")
  df[order(df$p.adjust, df$pvalue, df$rank_padj, df$rank_p, decreasing = FALSE), ]
}
has_terms <- function(x, min_n = 1) {
  tryCatch({
    df <- as.data.frame(x)
    !is.null(df) && nrow(df) >= min_n
  }, error = function(e) FALSE)
}
muffle_fgsea_ties <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("There are ties in the preranked stats", conditionMessage(m))) {
      invokeRestart("muffleMessage")
    }
  })
}


## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text  = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}


## ---------- Load WGCNA modules (already UniProt in Gene) ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "") |>
  dplyr::mutate(Gene = as.character(Gene))

modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(x))
module_sizes <- tibble::tibble(
  Module = names(modules_list),
  N_accessions = vapply(modules_list, length, integer(1))
)


## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db


## ---------- Inputs ----------
#comparison <- "neuron-regionLayerBaseline"
comparison <- "neuron-phenotypeWithinUnit"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)


## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))


## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Meta, README, module sizes
  capture.output(sessionInfo(), file = file.path(paths$meta, "sessionInfo.txt"))
  save_table(module_sizes, paths$modules$root, "module_sizes.csv")
  writeLines(c(
    paste0("data_path: ", data_path),
    paste0("context: ", context),
    paste0("timestamp: ", Sys.time()),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (via mapIds)",
    "DEG threshold for ORA: |log2fc| > 1",
    "GSEA ties: deterministically broken by sorting ties by id (values unchanged)"
  ), con = file.path(paths$meta, "params.txt"))

  readme_lines <- c(
    "# Analysis README",
    paste0("Comparison: ", comparison),
    paste0("Context: ", context),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (mapIds from org.Mm.eg.db, multiVals='first')",
    "DEG threshold for ORA: |log2fc| > 1",
    "Ontologies default: BP",
    "ComplexUpset: set sizes panel disabled for portability",
    "GSEA ties: sorted by id for deterministic order; values are NOT jittered",
    "Folders:",
    "  10_global/GO_*: GSEA and ORA across full universe",
    "  10_global/KEGG: gseKEGG results and mapping diagnostics",
    "  20_modules/: overlap, Fisher, per-module GO (GSEA/ORA)",
    "  90_pathview/: KEGG pathway renderings and status logs",
    "Check 00_meta/runlog_*.txt for stepwise notes."
  )
  writeLines(readme_lines, con = file.path(results_dir, "README.txt"))
  log_line(paths, file_tag, "Initialized results folders and README.")

  ## ----- Load comparison (assume UniProt accessions in col1) -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # Deterministic tie-breaking by id, values unchanged
  gene_list_tbl <- df |>
    dplyr::select(id, log2fc) |>
    dplyr::arrange(dplyr::desc(log2fc), id)
  gene_list <- gene_list_tbl$log2fc
  names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)

  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  # Per-file diagnostics
  pct_tied <- 100 * mean(duplicated(gene_list))
  save_table(tibble::tibble(
    n_rows = nrow(df),
    n_universe = length(universe_genes),
    n_deg_abs1 = length(top_genes),
    deg_threshold = 1,
    pct_tied_stats = pct_tied
  ), paths$meta, paste0("summary_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("Universe=%d; DEGs(|log2fc|>1)=%d; tied stats=%.1f%%",
                                    length(universe_genes), length(top_genes), pct_tied))

  organism <- "org.Mm.eg.db"
  onts <- c("BP") # extend if needed

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
      geneList = gene_list, ont = ont, keyType = "UNIPROT",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
      by = "fgsea", nPermSimple = 10000, eps = 0,
      verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
    ), error = function(e) { log_line(paths, file_tag, paste("gseGO error:", e$message)); NULL }) )

    if (!is.null(gse) && has_terms(gse, 1)) {
      res_stable <- stable_rank(gse@result)
      res_stable$note <- "fgsea ties handled by deterministic id sort (values unchanged)"
      save_table(res_stable, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))

      if (nrow(gse@result) >= 2) {
        ts_gse <- tryCatch(enrichplot::pairwise_termsim(gse), error = function(e) NULL)
        if (!is.null(ts_gse) && "compareClusterResult" %in% slotNames(ts_gse) &&
            nrow(ts_gse@compareClusterResult) > 0) {
          save_plot_dir(enrichplot::emapplot(ts_gse, showCategory = 10), dstP,
                        paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
        } else {
          writeLines("emapplot skipped: no enriched term after similarity calc.",
                     file.path(dstP, paste0("gseGO_", ont, "_emap_", file_tag, "_skipped.txt")))
        }
      }

      p_cnet <- tryCatch(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                         error = function(e) NULL)
      if (!is.null(p_cnet)) save_plot_dir(p_cnet, dstP, paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))

      rp <- tryCatch(enrichplot::ridgeplot(gse) + ggplot2::labs(x = "Enrichment Distribution"),
                     error = function(e) NULL)
      if (!is.null(rp)) save_plot_dir(rp, dstP, paste0("gseGO_", ont, "_ridge_", file_tag, ".svg"))

      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) { log_line(paths, file_tag, paste("simplify error:", e$message)); NULL })
      if (!is.null(gse2) && has_terms(gse2, 1)) {
        save_table(stable_rank(gse2@result), dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
      }
    } else {
      save_table(tibble::tibble(note = "gseGO empty or NULL"), dstT, paste0("gseGO_", ont, "_", file_tag, "_empty.csv"))
      writeLines("All plots skipped: no enriched term.", file.path(dstP, paste0("gseGO_", ont, "_", file_tag, "_skipped.txt")))
    }

    if (length(top_genes) >= 10) {
      ora <- tryCatch(clusterProfiler::enrichGO(
        gene = top_genes, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
      ), error = function(e) { log_line(paths, file_tag, paste("ORA error:", e$message)); NULL })
      if (!is.null(ora) && has_terms(ora, 1)) {
        save_table(stable_rank(ora@result), dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        if (nrow(ora@result) >= 2) {
          ts_ora <- tryCatch(enrichplot::pairwise_termsim(ora), error = function(e) NULL)
          if (!is.null(ts_ora) && "compareClusterResult" %in% slotNames(ts_ora) &&
              nrow(ts_ora@compareClusterResult) > 0) {
            save_plot_dir(enrichplot::emapplot(ts_ora, showCategory = 10), dstP,
                          paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
          } else {
            writeLines("emapplot skipped: no enriched term after similarity calc.",
                       file.path(dstP, paste0("ORA_", ont, "_emap_", file_tag, "_skipped.txt")))
          }
        }
        p_cnet2 <- tryCatch(enrichplot::cnetplot(ora, showCategory = 5), error = function(e) NULL)
        if (!is.null(p_cnet2)) save_plot_dir(p_cnet2, dstP, paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
      } else {
        save_table(tibble::tibble(note = "ORA returned NULL or empty"), dstT, paste0("ORA_", ont, "_", file_tag, "_empty.csv"))
      }
    } else {
      save_table(tibble::tibble(note = "Global ORA skipped: <10 DEGs"),
                 dstT, paste0("ORA_", ont, "_", file_tag, "_skipped.csv"))
      log_line(paths, file_tag, paste("Global ORA skipped for", ont, "due to insufficient DEGs."))
    }
  }

  ## ----- GLOBAL: KEGG (mapIds to 1:1 ENTREZ) -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  keep <- !is.na(entrez_map)
  kegg_gene_list <- gene_list[keep]
  names(kegg_gene_list) <- unname(entrez_map[keep])
  kegg_gene_list <- kegg_gene_list[order(kegg_gene_list, decreasing = TRUE, method = "radix")]

  # KEGG mapping diagnostics
  kmap_tbl <- tibble::tibble(
    id = universe_genes,
    ENTREZID = unname(entrez_map[universe_genes]),
    in_kegg = !is.na(entrez_map[universe_genes])
  )
  save_table(kmap_tbl, paths$glob$KEGG$T, paste0("kegg_mapping_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("KEGG mapIds: %d/%d mapped (%.1f%%)",
                                    sum(kmap_tbl$in_kegg), nrow(kmap_tbl), 100*mean(kmap_tbl$in_kegg)))

  kk2 <- NULL
  if (sum(kmap_tbl$in_kegg) >= 10) {
    kk2 <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseKEGG(
      geneList = kegg_gene_list, organism = "mmu",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
      keyType = "ncbi-geneid"
    ), error = function(e) { log_line(paths, file_tag, paste("gseKEGG error:", e$message)); NULL }) )
  } else {
    log_line(paths, file_tag, "gseKEGG skipped: <10 mapped ENTREZIDs.")
  }

  if (!is.null(kk2) && has_terms(kk2, 1)) {
    save_table(stable_rank(kk2@result), paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = 10, label_format = 30),
                  paths$glob$KEGG$P, paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    if (nrow(kk2@result) >= 2) {
      ts_kk <- tryCatch(enrichplot::pairwise_termsim(kk2), error = function(e) NULL)
      if (!is.null(ts_kk) && "compareClusterResult" %in% slotNames(ts_kk) &&
          nrow(ts_kk@compareClusterResult) > 0) {
        save_plot_dir(enrichplot::emapplot(ts_kk, showCategory = 10),
                      paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
      } else {
        writeLines("emapplot skipped: no enriched term after similarity calc.",
                   file.path(paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, "_skipped.txt")))
      }
    }
    p_cnet3 <- tryCatch(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                        error = function(e) NULL)
    if (!is.null(p_cnet3)) save_plot_dir(p_cnet3, paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))
    rp2 <- tryCatch(enrichplot::ridgeplot(kk2) + ggplot2::labs(x = "Enrichment distribution"),
                    error = function(e) NULL)
    if (!is.null(rp2)) save_plot_dir(rp2, paths$glob$KEGG$P, paste0("gseKEGG_ridge_", file_tag, ".svg"))
    if (nrow(kk2@result) > 0) {
      gp <- tryCatch(enrichplot::gseaplot(kk2, by = "all",
                                          title = kk2@result$Description[1], geneSetID = 1),
                     error = function(e) NULL)
      if (!is.null(gp)) save_plot_dir(gp, paths$glob$KEGG$P, paste0("gseKEGG_gseaplot_", file_tag, ".svg"))
    }
  } else {
    save_table(tibble::tibble(note = "No KEGG enrichment results (gseKEGG null or empty)."),
               paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, "_empty.csv"))
    log_line(paths, file_tag, "gseKEGG produced no results.")
  }

  ## ----- PATHVIEW with validated KEGG IDs + status log -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && has_terms(kk2, 1)) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  save_table(tibble::tibble(pv_ids = pv_ids), paths$pathview, paste0("pathview_ids_", file_tag, ".csv"))

  if (length(kegg_gene_list) > 0 && all(is.finite(as.numeric(kegg_gene_list)))) {
    writeLines(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
    log_line(paths, file_tag, sprintf("Pathview attempted for %d pathways.", length(pv_ids)))
  } else {
    writeLines(c(
      "No valid ENTREZ-mapped genes for Pathview or gene list non-numeric.",
      paste0("length(kegg_gene_list)=", length(kegg_gene_list))
    ), con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Pathview skipped: invalid kegg_gene_list.")
  }

  ## ----- MODULE-AWARE: overlaps, UpSet, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")

  # UpSet-style overlap (top modules by DEG overlap) with version-agnostic settings
  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), ]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))

  # Indicator matrix: one row per id, one col per set; keep id as separate column
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids

  # Compute basic check for non-zero overlaps
  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))

  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      # ComplexUpset expects data first; convert indicators to logicals
      library(ComplexUpset)
      set_cols <- setdiff(colnames(mat), "id")
      mat[set_cols] <- lapply(mat[set_cols], function(x) as.logical(x))
      p_up <- ComplexUpset::upset(
        data = mat,
        intersect = set_cols,
        base_annotations = list('Intersection size' = ComplexUpset::intersection_size()),
        set_sizes = FALSE,
        min_size = 5
      )
      save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      # UpSetR fallback
      library(UpSetR)
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
    } else {
      save_table(tibble::tibble(note="No UpSet library available (ComplexUpset/UpSetR)."),
                 paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped_no_pkg.csv"))
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
               paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }

  # Barplot of overlap counts and percentages
  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, ]
  if (nrow(ov2) > 0) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2, 20)$OverlapDEG)*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No non-zero overlaps to plot.", con = file.path(paths$modules$overlap$P, paste0("overlap_status_", file_tag, ".txt")))
  }

  # Fisher enrichment of DEGs in modules + plot
  fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
    deg_set <- intersect(unique(deg_genes), universe_genes)
    bg_set  <- setdiff(universe_genes, deg_set)
    out <- lapply(names(modules_list), function(m) {
      mod_genes <- intersect(modules_list[[m]], universe_genes)
      a <- length(intersect(mod_genes, deg_set))
      b <- length(setdiff(mod_genes, deg_set))
      c <- length(setdiff(deg_set, mod_genes))
      d <- length(setdiff(bg_set, mod_genes))
      mat <- matrix(c(a,b,c,d), nrow=2)
      ft <- fisher.test(mat, alternative = "greater")
      data.frame(
        Module = m, InModule_DE = a, InModule_NotDE = b,
        NotInModule_DE = c, NotInModule_NotDE = d,
        OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    res <- dplyr::bind_rows(out)
    res$Padj <- p.adjust(res$Pvalue, method = "BH")
    dplyr::arrange(res, Padj, Pvalue)
  }
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))

  if (nrow(fisher_tbl) > 0) {
    ft_plot <- ggplot2::ggplot(head(fisher_tbl[order(fisher_tbl$Padj, fisher_tbl$Pvalue), ], 20),
                               ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                            y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#6a51a3") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Fisher enrichment (DEGs in modules) — ", context),
                    x = "Module", y = "-log10(Padj)")
    save_plot_dir(ft_plot, paths$modules$fisher$P, paste0("fisher_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No Fisher rows to plot.", con = file.path(paths$modules$fisher$P, paste0("fisher_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Fisher plot skipped: empty table.")
  }

  ## ----- Per-module analyses -----
  sel_modules <- fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)]
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) {
      log_line(paths, file_tag, paste("Module", m, "skipped: <10 genes in universe."))
      next
    }

    for (ont in onts) {
      # Per-module nested folders
      m_root     <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea   <- mk(m_root, "Tables", "gsea")
      m_P_gsea   <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull<- mk(m_root, "Tables", "ora_full")
      m_P_orafull<- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]
      # Order by value desc, break ties by id (names) deterministically
      m_gene_list <- m_gene_list[order(m_gene_list, decreasing = TRUE, method = "radix")]
      if (length(m_gene_list) >= 10) {
        m_gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
          geneList = m_gene_list,
          ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "gseGO error:", e$message)); NULL }) )
        if (!is.null(m_gse) && has_terms(m_gse, 1)) {
          resm <- stable_rank(m_gse@result); resm$note <- "fgsea ties handled by deterministic id sort"
          save_table(resm, m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          p_dot <- tryCatch(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust",
                                                showCategory = 10, label_format = 30), error = function(e) NULL)
          if (!is.null(p_dot)) save_plot_dir(p_dot, m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "Module GSEA skipped: <10 genes."), m_T_gsea,
                   paste0("gseGO_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on DEG-intersection
      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA DEG error:", e$message)); NULL })
        if (!is.null(m_ora) && has_terms(m_ora, 1)) {
          save_table(stable_rank(m_ora@result), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
          if (nrow(m_ora@result) >= 2) {
            ts_mora <- tryCatch(enrichplot::pairwise_termsim(m_ora), error = function(e) NULL)
            if (!is.null(ts_mora) && "compareClusterResult" %in% slotNames(ts_mora) &&
                nrow(ts_mora@compareClusterResult) > 0) {
              save_plot_dir(enrichplot::emapplot(ts_mora, showCategory = 10),
                            m_P_oradeg, paste0("ORA_", ont, "_", m, "_emap_", file_tag, ".svg"))
            } else {
              writeLines("emapplot skipped: no enriched term after similarity calc.",
                         file.path(m_P_oradeg, paste0("ORA_", ont, "_", m, "_emap_", file_tag, "_skipped.txt")))
            }
          }
          p_cnet_mora <- tryCatch(enrichplot::cnetplot(m_ora, showCategory = 5), error = function(e) NULL)
          if (!is.null(p_cnet_mora)) save_plot_dir(p_cnet_mora, m_P_oradeg, paste0("ORA_", ont, "_", m, "_cnet_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(
          note = "insufficient DEGs for ORA",
          n_module_genes_in_universe = length(m_genes),
          n_deg_in_module = length(intersect(top_genes, m_genes)),
          deg_threshold = 1
        ), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on full module gene set
      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA full error:", e$message)); NULL })
        if (!is.null(m_ora_all) && has_terms(m_ora_all, 1)) {
          save_table(stable_rank(m_ora_all@result), m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "ORA full module skipped: <5 genes."),
                   m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }
    } # ont
  } # module
} # for each comparison file





























############################################################
# clusterProfiler + WGCNA modules (mouse)
# IDs: UniProt accessions already in input (no mapping step)
# Deterministic tie-breaking for fgsea (order by id), no value jitter
# ComplexUpset set sizes panel disabled for version-agnostic plots
# Robust guards, diagnostics, run logs, README
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","vroom","UniProt.ws","BiocParallel",
    "BiocGenerics","S4Vectors","IRanges","GenomeInfoDb","GenomicRanges","SummarizedExperiment",
    "Biobase","GO.db","UpSetR","ComplexUpset"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Directory + I/O helpers ----------
mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP   = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                  P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF   = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                  P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC   = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                  P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                  P = mk(results_dir, "10_global","KEGG","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots"))
    ),
    pathview = mk(results_dir, "90_pathview")
  )
}
save_table <- function(df, dirpath, fname) {
  readr::write_csv(df, file.path(dirpath, fname))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}
log_line <- function(paths, tag, text) {
  lf <- file.path(paths$meta, paste0("runlog_", tag, ".txt"))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text), file = lf, append = TRUE)
}
stable_rank <- function(df) {
  if (is.null(df) || nrow(df)==0) return(df)
  if (!("pvalue" %in% names(df)) && ("p.adjust" %in% names(df))) df$pvalue <- df$p.adjust
  if (!("p.adjust" %in% names(df)) && ("pvalue" %in% names(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df$rank_p <- rank(df$pvalue, ties.method = "min")
  df$rank_padj <- rank(df$p.adjust, ties.method = "min")
  df[order(df$p.adjust, df$pvalue, df$rank_padj, df$rank_p, decreasing = FALSE), ]
}
has_terms <- function(x, min_n = 1) {
  tryCatch({
    df <- as.data.frame(x)
    !is.null(df) && nrow(df) >= min_n
  }, error = function(e) FALSE)
}
muffle_fgsea_ties <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("There are ties in the preranked stats", conditionMessage(m))) {
      invokeRestart("muffleMessage")
    }
  })
}

## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text  = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

## ---------- Load WGCNA modules (already UniProt in Gene) ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "") |>
  dplyr::mutate(Gene = as.character(Gene))

modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(x))
module_sizes <- tibble::tibble(
  Module = names(modules_list),
  N_accessions = vapply(modules_list, length, integer(1))
)

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Inputs ----------
#comparison <- "neuron-regionLayerBaseline"
comparison <- "neuron-phenotypeWithinUnit"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Meta, README, module sizes
  capture.output(sessionInfo(), file = file.path(paths$meta, "sessionInfo.txt"))
  save_table(module_sizes, paths$modules$root, "module_sizes.csv")
  writeLines(c(
    paste0("data_path: ", data_path),
    paste0("context: ", context),
    paste0("timestamp: ", Sys.time()),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (via mapIds)",
    "DEG threshold for ORA: |log2fc| > 1",
    "GSEA ties: deterministically broken by sorting ties by id (values unchanged)"
  ), con = file.path(paths$meta, "params.txt"))

  readme_lines <- c(
    "# Analysis README",
    paste0("Comparison: ", comparison),
    paste0("Context: ", context),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (mapIds from org.Mm.eg.db, multiVals='first')",
    "DEG threshold for ORA: |log2fc| > 1",
    "Ontologies default: BP",
    "ComplexUpset: set sizes panel disabled for portability",
    "GSEA ties: sorted by id for deterministic order; values are NOT jittered",
    "Folders:",
    "  10_global/GO_*: GSEA and ORA across full universe",
    "  10_global/KEGG: gseKEGG results and mapping diagnostics",
    "  20_modules/: overlap, Fisher, per-module GO (GSEA/ORA)",
    "  90_pathview/: KEGG pathway renderings and status logs",
    "Check 00_meta/runlog_*.txt for stepwise notes."
  )
  writeLines(readme_lines, con = file.path(results_dir, "README.txt"))
  log_line(paths, file_tag, "Initialized results folders and README.")

  ## ----- Load comparison (assume UniProt accessions in col1) -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # Deterministic tie-breaking by id, values unchanged
  gene_list_tbl <- df |>
    dplyr::select(id, log2fc) |>
    dplyr::arrange(dplyr::desc(log2fc), id)
  gene_list <- gene_list_tbl$log2fc
  names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)

  # Extended logging for universe
  pct_tied <- 100 * mean(duplicated(gene_list))
  log_line(paths, file_tag, sprintf("Universe genes=%d; DEGs(|log2fc|>1)=%d; pct_tied=%.2f%%",
                                    length(universe_genes), sum(abs(df$log2fc)>1), pct_tied))
  log_line(paths, file_tag, sprintf("Top gene_list head: %s",
    paste(utils::head(sprintf("%s:%.3f", names(gene_list), gene_list)), collapse=", ")))

  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  # Per-file diagnostics table
  save_table(tibble::tibble(
    n_rows = nrow(df),
    n_universe = length(universe_genes),
    n_deg_abs1 = length(top_genes),
    deg_threshold = 1,
    pct_tied_stats = pct_tied
  ), paths$meta, paste0("summary_", file_tag, ".csv"))

  organism <- "org.Mm.eg.db"
  onts <- c("BP")

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    log_line(paths, file_tag, sprintf("Starting gseGO ont=%s; n=%d; anyDuplicated(names)=%d",
                                      ont, length(gene_list), anyDuplicated(names(gene_list))))

    gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
      geneList = gene_list, ont = ont, keyType = "UNIPROT",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
      by = "fgsea", nPermSimple = 10000, eps = 0,
      verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
    ), error = function(e) { log_line(paths, file_tag, paste("gseGO error:", e$message)); NULL }) )

    if (!is.null(gse) && has_terms(gse, 1)) {
      res_stable <- stable_rank(gse@result)
      res_stable$note <- "fgsea ties handled by deterministic id sort (values unchanged)"
      save_table(res_stable, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))

      log_line(paths, file_tag, sprintf("gseGO terms=%d; min p.adjust=%.3g",
        nrow(gse@result), min(gse@result$p.adjust, na.rm=TRUE)))

      if (nrow(gse@result) >= 2) {
        ts_gse <- tryCatch(enrichplot::pairwise_termsim(gse), error = function(e) NULL)
        if (!is.null(ts_gse) && "compareClusterResult" %in% slotNames(ts_gse) &&
            nrow(ts_gse@compareClusterResult) > 0) {
          save_plot_dir(enrichplot::emapplot(ts_gse, showCategory = 10), dstP,
                        paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
        } else {
          writeLines("emapplot skipped: no enriched term after similarity calc.",
                     file.path(dstP, paste0("gseGO_", ont, "_emap_", file_tag, "_skipped.txt")))
        }
      }

      p_cnet <- tryCatch(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                         error = function(e) NULL)
      if (!is.null(p_cnet)) save_plot_dir(p_cnet, dstP, paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))

      rp <- tryCatch(enrichplot::ridgeplot(gse) + ggplot2::labs(x = "Enrichment Distribution"),
                     error = function(e) NULL)
      if (!is.null(rp)) save_plot_dir(rp, dstP, paste0("gseGO_", ont, "_ridge_", file_tag, ".svg"))

      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) { log_line(paths, file_tag, paste("simplify error:", e$message)); NULL })
      if (!is.null(gse2) && has_terms(gse2, 1)) {
        save_table(stable_rank(gse2@result), dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
        log_line(paths, file_tag, sprintf("simplify: %d -> %d terms",
          nrow(gse@result), nrow(gse2@result)))
      }
    } else {
      save_table(tibble::tibble(note = "gseGO empty or NULL"), dstT, paste0("gseGO_", ont, "_", file_tag, "_empty.csv"))
      writeLines("All plots skipped: no enriched term.", file.path(dstP, paste0("gseGO_", ont, "_", file_tag, "_skipped.txt")))
    }

    if (length(top_genes) >= 10) {
      ora <- tryCatch(clusterProfiler::enrichGO(
        gene = top_genes, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
      ), error = function(e) { log_line(paths, file_tag, paste("ORA error:", e$message)); NULL })
      if (!is.null(ora) && has_terms(ora, 1)) {
        save_table(stable_rank(ora@result), dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        if (nrow(ora@result) >= 2) {
          ts_ora <- tryCatch(enrichplot::pairwise_termsim(ora), error = function(e) NULL)
          if (!is.null(ts_ora) && "compareClusterResult" %in% slotNames(ts_ora) &&
              nrow(ts_ora@compareClusterResult) > 0) {
            save_plot_dir(enrichplot::emapplot(ts_ora, showCategory = 10), dstP,
                          paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
          } else {
            writeLines("emapplot skipped: no enriched term after similarity calc.",
                       file.path(dstP, paste0("ORA_", ont, "_emap_", file_tag, "_skipped.txt")))
          }
        }
        p_cnet2 <- tryCatch(enrichplot::cnetplot(ora, showCategory = 5), error = function(e) NULL)
        if (!is.null(p_cnet2)) save_plot_dir(p_cnet2, dstP, paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
        log_line(paths, file_tag, sprintf("ORA terms=%d; min pvalue=%.3g",
          nrow(ora@result), min(ora@result$pvalue, na.rm=TRUE)))
      } else {
        save_table(tibble::tibble(note = "ORA returned NULL or empty"), dstT, paste0("ORA_", ont, "_", file_tag, "_empty.csv"))
      }
    } else {
      save_table(tibble::tibble(note = "Global ORA skipped: <10 DEGs"),
                 dstT, paste0("ORA_", ont, "_", file_tag, "_skipped.csv"))
      log_line(paths, file_tag, paste("Global ORA skipped for", ont, "due to insufficient DEGs."))
    }
  }

  ## ----- GLOBAL: KEGG (UNIPROT -> ENTREZ; de-dup; detailed log) -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  n_univ <- length(universe_genes)
  n_mapped <- sum(!is.na(entrez_map[universe_genes]))
  n_unmapped <- n_univ - n_mapped
  log_line(paths, file_tag, sprintf("KEGG mapIds: %d/%d mapped (%.1f%%), %d unmapped",
                                    n_mapped, n_univ, 100*n_mapped/n_univ, n_unmapped))

  # KEGG mapping diagnostics table
  kmap_tbl <- tibble::tibble(
    id = universe_genes,
    ENTREZID = unname(entrez_map[universe_genes]),
    in_kegg = !is.na(entrez_map[universe_genes])
  )
  save_table(kmap_tbl, paths$glob$KEGG$T, paste0("kegg_mapping_", file_tag, ".csv"))

  # Build ranked stats named by ENTREZ, then de-duplicate names
  kegg_gene_list <- gene_list[!is.na(entrez_map)]
  names(kegg_gene_list) <- unname(entrez_map[!is.na(entrez_map)])

  dup_mask <- duplicated(names(kegg_gene_list))
  n_dups <- sum(dup_mask)
  if (n_dups > 0) {
    dup_ids <- unique(names(kegg_gene_list)[dup_mask])
    log_line(paths, file_tag, sprintf("gseKEGG: duplicate ENTREZ names=%d; examples: %s",
      n_dups, paste(utils::head(dup_ids, 10), collapse=", ")))
  }

  # Strategy: keep first occurrence per ENTREZ (optionally aggregate)
  keep_idx <- !duplicated(names(kegg_gene_list))
  removed_dups <- sum(!keep_idx)
  kegg_gene_list <- kegg_gene_list[keep_idx]

  # Optional aggregation alternative (commented):
  # kegg_gene_list <- tapply(kegg_gene_list, names(kegg_gene_list), mean)

  # Sort decreasing for fgsea
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE, method = "radix")

  log_line(paths, file_tag, sprintf(
    "gseKEGG input: n=%d; anyDuplicated(names)=%d; pct ties=%.2f%%; removed_dups=%d",
    length(kegg_gene_list),
    anyDuplicated(names(kegg_gene_list)),
    100*mean(duplicated(kegg_gene_list)),
    removed_dups
  ))

  kk2 <- NULL
  if (length(kegg_gene_list) >= 10) {
    kk2 <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseKEGG(
      geneList = kegg_gene_list, organism = "mmu",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
      keyType = "ncbi-geneid"
    ), error = function(e) { log_line(paths, file_tag, paste("gseKEGG error:", e$message)); NULL }) )
  } else {
    log_line(paths, file_tag, "gseKEGG skipped: <10 unique ENTREZ in ranked list after de-dup.")
  }

  if (!is.null(kk2) && has_terms(kk2, 1)) {
    save_table(stable_rank(kk2@result), paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = 10, label_format = 30),
                  paths$glob$KEGG$P, paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    if (nrow(kk2@result) >= 2) {
      ts_kk <- tryCatch(enrichplot::pairwise_termsim(kk2), error = function(e) NULL)
      if (!is.null(ts_kk) && "compareClusterResult" %in% slotNames(ts_kk) &&
          nrow(ts_kk@compareClusterResult) > 0) {
        save_plot_dir(enrichplot::emapplot(ts_kk, showCategory = 10),
                      paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
      } else {
        writeLines("emapplot skipped: no enriched term after similarity calc.",
                   file.path(paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, "_skipped.txt")))
      }
    }
    p_cnet3 <- tryCatch(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                        error = function(e) NULL)
    if (!is.null(p_cnet3)) save_plot_dir(p_cnet3, paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))
    rp2 <- tryCatch(enrichplot::ridgeplot(kk2) + ggplot2::labs(x = "Enrichment distribution"),
                    error = function(e) NULL)
    if (!is.null(rp2)) save_plot_dir(rp2, paths$glob$KEGG$P, paste0("gseKEGG_ridge_", file_tag, ".svg"))
    if (nrow(kk2@result) > 0) {
      gp <- tryCatch(enrichplot::gseaplot(kk2, by = "all",
                                          title = kk2@result$Description[1], geneSetID = 1),
                     error = function(e) NULL)
      if (!is.null(gp)) save_plot_dir(gp, paths$glob$KEGG$P, paste0("gseKEGG_gseaplot_", file_tag, ".svg"))
    }
  } else {
    save_table(tibble::tibble(note = "No KEGG enrichment results (gseKEGG null or empty)."),
               paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, "_empty.csv"))
    log_line(paths, file_tag, "gseKEGG produced no results.")
  }

  ## ----- PATHVIEW with validated KEGG IDs + status log -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && has_terms(kk2, 1)) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  save_table(tibble::tibble(pv_ids = pv_ids), paths$pathview, paste0("pathview_ids_", file_tag, ".csv"))

  if (length(kegg_gene_list) > 0 && all(is.finite(as.numeric(kegg_gene_list)))) {
    writeLines(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
    log_line(paths, file_tag, sprintf("Pathview attempted for %d pathways.", length(pv_ids)))
  } else {
    writeLines(c(
      "No valid ENTREZ-mapped genes for Pathview or gene list non-numeric.",
      paste0("length(kegg_gene_list)=", length(kegg_gene_list))
    ), con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Pathview skipped: invalid kegg_gene_list.")
  }

  ## ----- MODULE-AWARE: overlaps, UpSet, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")
  log_line(paths, file_tag, sprintf("Modules=%d; with DEG overlap>0: %d",
                                    length(modules_list), sum(overlap_counts$OverlapDEG>0)))

  # UpSet-style overlap (top modules by DEG overlap) with version-agnostic settings
  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), ]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids

  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))
  log_line(paths, file_tag, sprintf("UpSet candidates: sets=%d; all_ids=%d; nonzero_sets=%d",
                                    length(sets), length(all_ids), overlap_nonzero))

  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      library(ComplexUpset)
      set_cols <- setdiff(colnames(mat), "id")
      mat[set_cols] <- lapply(mat[set_cols], function(x) as.logical(x))
      p_up <- ComplexUpset::upset(
        data = mat,
        intersect = set_cols,
        base_annotations = list('Intersection size' = ComplexUpset::intersection_size()),
        set_sizes = FALSE,
        min_size = 5
      )
      save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
      log_line(paths, file_tag, "ComplexUpset plot saved.")
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      library(UpSetR)
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
      log_line(paths, file_tag, "UpSetR plot saved.")
    } else {
      save_table(tibble::tibble(note="No UpSet library available (ComplexUpset/UpSetR)."),
                 paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped_no_pkg.csv"))
      log_line(paths, file_tag, "UpSet skipped: no package available.")
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
               paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }

  # Barplot of overlap counts and percentages
  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, ]
  if (nrow(ov2) > 0) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2, 20)$OverlapDEG)*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No non-zero overlaps to plot.", con = file.path(paths$modules$overlap$P, paste0("overlap_status_", file_tag, ".txt")))
  }

  # Fisher enrichment of DEGs in modules + plot
  fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
    deg_set <- intersect(unique(deg_genes), universe_genes)
    bg_set  <- setdiff(universe_genes, deg_set)
    out <- lapply(names(modules_list), function(m) {
      mod_genes <- intersect(modules_list[[m]], universe_genes)
      a <- length(intersect(mod_genes, deg_set))
      b <- length(setdiff(mod_genes, deg_set))
      c <- length(setdiff(deg_set, mod_genes))
      d <- length(setdiff(bg_set, mod_genes))
      mm <- matrix(c(a,b,c,d), nrow=2)
      ft <- fisher.test(mm, alternative = "greater")
      data.frame(
        Module = m, InModule_DE = a, InModule_NotDE = b,
        NotInModule_DE = c, NotInModule_NotDE = d,
        OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    res <- dplyr::bind_rows(out)
    res$Padj <- p.adjust(res$Pvalue, method = "BH")
    dplyr::arrange(res, Padj, Pvalue)
  }
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("Fisher tested modules=%d; significant (Padj<0.05)=%d",
    nrow(fisher_tbl), sum(fisher_tbl$Padj<0.05, na.rm=TRUE)))

  if (nrow(fisher_tbl) > 0) {
    ft_plot <- ggplot2::ggplot(head(fisher_tbl[order(fisher_tbl$Padj, fisher_tbl$Pvalue), ], 20),
                               ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                            y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#6a51a3") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Fisher enrichment (DEGs in modules) — ", context),
                    x = "Module", y = "-log10(Padj)")
    save_plot_dir(ft_plot, paths$modules$fisher$P, paste0("fisher_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No Fisher rows to plot.", con = file.path(paths$modules$fisher$P, paste0("fisher_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Fisher plot skipped: empty table.")
  }

  ## ----- Per-module analyses -----
  sel_modules <- fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)]
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    log_line(paths, file_tag, sprintf("Module %s: n_in_universe=%d; n_deg=%d",
      m, length(m_genes), length(intersect(m_genes, top_genes))))
    if (length(m_genes) < 10) {
      log_line(paths, file_tag, paste("Module", m, "skipped: <10 genes in universe."))
      next
    }

    for (ont in onts) {
      # Per-module nested folders
      m_root     <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea   <- mk(m_root, "Tables", "gsea")
      m_P_gsea   <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull<- mk(m_root, "Tables", "ora_full")
      m_P_orafull<- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]
      m_gene_list <- m_gene_list[order(m_gene_list, decreasing = TRUE, method = "radix")]
      log_line(paths, file_tag, sprintf("Module %s gseGO ont=%s; n=%d; anyDuplicated(names)=%d",
        m, ont, length(m_gene_list), anyDuplicated(names(m_gene_list))))

      if (length(m_gene_list) >= 10) {
        m_gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
          geneList = m_gene_list,
          ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "gseGO error:", e$message)); NULL }) )
        if (!is.null(m_gse) && has_terms(m_gse, 1)) {
          resm <- stable_rank(m_gse@result); resm$note <- "fgsea ties handled by deterministic id sort"
          save_table(resm, m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          p_dot <- tryCatch(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust",
                                                showCategory = 10, label_format = 30), error = function(e) NULL)
          if (!is.null(p_dot)) save_plot_dir(p_dot, m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
          log_line(paths, file_tag, sprintf("Module %s gseGO terms=%d; min p.adjust=%.3g",
            m, nrow(m_gse@result), min(m_gse@result$p.adjust, na.rm=TRUE)))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "Module GSEA skipped: <10 genes."), m_T_gsea,
                   paste0("gseGO_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on DEG-intersection
      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA DEG error:", e$message)); NULL })
        if (!is.null(m_ora) && has_terms(m_ora, 1)) {
          save_table(stable_rank(m_ora@result), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
          log_line(paths, file_tag, sprintf("Module %s ORA(DEG) terms=%d; min pvalue=%.3g",
            m, nrow(m_ora@result), min(m_ora@result$pvalue, na.rm=TRUE)))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(
          note = "insufficient DEGs for ORA",
          n_module_genes_in_universe = length(m_genes),
          n_deg_in_module = length(intersect(top_genes, m_genes)),
          deg_threshold = 1
        ), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on full module gene set
      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA full error:", e$message)); NULL })
        if (!is.null(m_ora_all) && has_terms(m_ora_all, 1)) {
          save_table(stable_rank(m_ora_all@result), m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
          log_line(paths, file_tag, sprintf("Module %s ORA(full) terms=%d; min p.adjust=%.3g",
            m, nrow(m_ora_all@result), min(m_ora_all@result$p.adjust, na.rm=TRUE)))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "ORA full module skipped: <5 genes."),
                   m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }
    } # ont
  } # module
} # for each comparison file


































############################################################
# clusterProfiler + WGCNA modules (mouse) — extended figures
# - Executive summary export (key counts, top terms, paths)
# - Redundancy-reduced enrichment themes (network clusters)
# - Leading-edge gene heatmaps for top GSEA terms
# - NES vs -log10(FDR) landscape plot
# - Optional Reactome ORA/GSEA (ReactomePA)
# - Rich runlog: provenance, parameters, timings, previews
# Notes:
#   Best-practice figure types: dot/bar, emap, cnet, heatmaps,
#   enrichment landscape, redundancy clustering, module panels.
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","vroom","UniProt.ws","BiocParallel",
    "BiocGenerics","S4Vectors","IRanges","GenomeInfoDb","GenomicRanges","SummarizedExperiment",
    "Biobase","GO.db","UpSetR","ComplexUpset","igraph","ggraph"
  )
  # Optional Reactome; will be guarded
  optional_bioc <- c("ReactomePA")
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST","ReactomePA")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Directory + I/O helpers ----------
mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP   = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                  P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF   = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                  P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC   = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                  P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                  P = mk(results_dir, "10_global","KEGG","Plots")),
      Reactome = list(T = mk(results_dir, "10_global","Reactome","Tables"),
                      P = mk(results_dir, "10_global","Reactome","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots")),
      themes  = mk(results_dir, "20_modules","themes_networks")
    ),
    pathview = mk(results_dir, "90_pathview"),
    summary = mk(results_dir, "99_summary")
  )
}
save_table <- function(df, dirpath, fname) {
  readr::write_csv(df, file.path(dirpath, fname))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}
log_line <- function(paths, tag, text) {
  lf <- file.path(paths$meta, paste0("runlog_", tag, ".txt"))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text), file = lf, append = TRUE)
}
stable_rank <- function(df) {
  if (is.null(df) || nrow(df)==0) return(df)
  if (!("pvalue" %in% names(df)) && ("p.adjust" %in% names(df))) df$pvalue <- df$p.adjust
  if (!("p.adjust" %in% names(df)) && ("pvalue" %in% names(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df$rank_p <- rank(df$pvalue, ties.method = "min")
  df$rank_padj <- rank(df$p.adjust, ties.method = "min")
  df[order(df$p.adjust, df$pvalue, df$rank_padj, df$rank_p, decreasing = FALSE), ]
}
has_terms <- function(x, min_n = 1) {
  tryCatch({
    df <- as.data.frame(x)
    !is.null(df) && nrow(df) >= min_n
  }, error = function(e) FALSE)
}
muffle_fgsea_ties <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("There are ties in the preranked stats", conditionMessage(m))) {
      invokeRestart("muffleMessage")
    }
  })
}

# Timing helper
with_timing <- function(expr, paths, tag, label) {
  t0 <- Sys.time()
  res <- NULL; err <- NULL
  tryCatch({ res <- eval.parent(substitute(expr)) },
           error = function(e) { err <<- e })
  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  log_line(paths, tag, sprintf("%s elapsed=%.2fs %s",
                               label, elapsed,
                               if (!is.null(err)) paste("| error:", err$message) else "" ))
  if (!is.null(err)) stop(err)
  res
}

## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text  = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Redundancy reduction: build term-sim graph and cluster
build_theme_graph <- function(enrich_obj, top_n = 30, sim_cut = 0.35) {
  if (is.null(enrich_obj) || !has_terms(enrich_obj, 2)) return(NULL)
  ts <- tryCatch(enrichplot::pairwise_termsim(enrich_obj), error = function(e) NULL)
  if (is.null(ts) || !"similarity" %in% slotNames(ts)) return(NULL)
  df <- as.data.frame(ts@compareClusterResult)
  if (!nrow(df)) return(NULL)
  df <- df[order(df$p.adjust, df$pvalue), ]
  df <- head(df, top_n)
  sim <- ts@termsim
  sim <- sim[rownames(sim) %in% df$ID, colnames(sim) %in% df$ID, drop = FALSE]
  edges <- as.data.frame(as.table(sim), stringsAsFactors = FALSE)
  edges <- edges[edges$Freq >= sim_cut & edges$Var1 != edges$Var2, ]
  if (!nrow(edges)) return(NULL)
  g <- igraph::graph_from_data_frame(edges[, c("Var1","Var2")], directed = FALSE)
  igraph::V(g)$label <- df$Description[match(igraph::V(g)$name, df$ID)]
  igraph::V(g)$padj <- df$p.adjust[match(igraph::V(g)$name, df$ID)]
  igraph::V(g)$size <- -log10(pmax(1e-300, igraph::V(g)$padj))
  g
}

plot_theme_graph <- function(g, title_txt = "Enrichment themes") {
  if (is.null(g)) return(NULL)
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.2) +
    ggraph::geom_node_point(ggplot2::aes(size = size, color = size)) +
    ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE, size = 3) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title_txt, color = "-log10(FDR)", size = "-log10(FDR)")
}

# GSEA landscape scatter: NES vs -log10(FDR)
plot_gsea_landscape <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj) || !has_terms(gsea_obj, 1)) return(NULL)
  df <- as.data.frame(gsea_obj@result)
  if (!nrow(df) || !"NES" %in% names(df)) return(NULL)
  df$mlogFDR <- -log10(pmax(1e-300, df$p.adjust))
  ggplot2::ggplot(df, ggplot2::aes(x = NES, y = mlogFDR)) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(title = title_txt, x = "NES", y = "-log10(FDR)") +
    ggplot2::theme_minimal()
}

# Leading-edge heatmap prep: build a compact heatmap matrix from top terms
prep_leading_edge_matrix <- function(gsea_obj, gene_list, top_terms = 5) {
  if (is.null(gsea_obj) || !has_terms(gsea_obj, 1)) return(NULL)
  res <- as.data.frame(gsea_obj@result)
  res <- res[order(res$p.adjust, -abs(res$NES)), ]
  res <- head(res, top_terms)
  if (!nrow(res)) return(NULL)
  # geneID is /-separated genes; for UNIPROT keyType, cnetplot uses internal mapping, but leading edge genes are not always provided explicitly.
  # Use core enrichment column if available in result (clusterProfiler stores leading edge as 'core_enrichment')
  if (!"core_enrichment" %in% colnames(res)) return(NULL)
  le_genes <- unique(unlist(strsplit(paste(res$core_enrichment, collapse = "/"), "/")))
  le_genes <- intersect(le_genes, names(gene_list))
  if (!length(le_genes)) return(NULL)
  mat <- matrix(gene_list[le_genes], nrow = length(le_genes), dimnames = list(le_genes, "log2fc"))
  mat
}

plot_leading_edge_heatmap <- function(mat, title_txt) {
  if (is.null(mat)) return(NULL)
  df <- data.frame(Gene = rownames(mat), log2fc = mat[,1], row.names = NULL, stringsAsFactors = FALSE)
  df$Gene <- factor(df$Gene, levels = df$Gene[order(df$log2fc, decreasing = TRUE)])
  ggplot2::ggplot(df, ggplot2::aes(x = "Leading edge", y = Gene, fill = log2fc)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    ggplot2::labs(title = title_txt, x = NULL, y = "Genes (leading edge)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
}

## ---------- Load WGCNA modules ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "") |>
  dplyr::mutate(Gene = as.character(Gene))

modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(x))
module_sizes <- tibble::tibble(
  Module = names(modules_list),
  N_accessions = vapply(modules_list, length, integer(1))
)

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Inputs ----------
comparison <- "neuron-phenotypeWithinUnit"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Provenance snapshot
  capture.output(sessionInfo(), file = file.path(paths$meta, "sessionInfo.txt"))
  writeLines(capture.output(Sys.getenv()), con = file.path(paths$meta, "env.txt"))
  writeLines(paste("Start:", Sys.time()), con = file.path(paths$meta, "run_start.txt"))
  save_table(module_sizes, paths$modules$root, "module_sizes.csv")

  # Parameters README
  writeLines(c(
    "# Analysis README",
    paste0("Comparison: ", comparison),
    paste0("Context: ", context),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (mapIds from org.Mm.eg.db, multiVals='first')",
    "DEG threshold for ORA: |log2fc| > 1",
    "Ontologies default: BP",
    "ComplexUpset: set sizes panel disabled for portability",
    "GSEA ties: sorted by id for deterministic order; values are NOT jittered",
    "Folders:",
    "  10_global/GO_*: GSEA and ORA across full universe",
    "  10_global/KEGG: gseKEGG results and mapping diagnostics",
    "  10_global/Reactome: optional Reactome enrichment",
    "  20_modules/: overlap, Fisher, per-module GO (GSEA/ORA), themes",
    "  90_pathview/: KEGG pathway renderings and status logs",
    "  99_summary/: executive summary",
    "Check 00_meta/runlog_*.txt for stepwise notes."
  ), con = file.path(results_dir, "README.txt"))
  log_line(paths, file_tag, "Initialized results folders and README.")

  ## ----- Load comparison -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # Rank list
  gene_list_tbl <- df |>
    dplyr::select(id, log2fc) |>
    dplyr::arrange(dplyr::desc(log2fc), id)
  gene_list <- gene_list_tbl$log2fc
  names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)
  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id
  pct_tied <- 100 * mean(duplicated(gene_list))
  log_line(paths, file_tag, sprintf("Universe=%d; DEGs(|log2fc|>1)=%d; pct_tied=%.2f%%; preview=%s",
    length(universe_genes), length(top_genes), pct_tied,
    paste(utils::head(sprintf("%s:%.3f", names(gene_list), gene_list)), collapse=", ")))

  organism <- "org.Mm.eg.db"
  onts <- c("BP")

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    gse <- with_timing({
      muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
        geneList = gene_list, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        by = "fgsea", nPermSimple = 10000, eps = 0,
        verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
      ), error = function(e) { log_line(paths, file_tag, paste("gseGO error:", e$message)); NULL }) )
    }, paths, file_tag, sprintf("gseGO(%s)", ont))

    if (!is.null(gse) && has_terms(gse, 1)) {
      res_stable <- stable_rank(gse@result)
      res_stable$note <- "fgsea ties handled by deterministic id sort (values unchanged)"
      save_table(res_stable, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))
      save_plot_dir(plot_gsea_landscape(gse, paste0("GSEA landscape (", ont, ") — ", context)),
                    dstP, paste0("gseGO_", ont, "_landscape_", file_tag, ".svg"))

      # Themes graph (redundancy reduction)
      g_theme <- build_theme_graph(gse, top_n = 30, sim_cut = 0.35)
      save_plot_dir(plot_theme_graph(g_theme, paste0("GO themes (", ont, ") — ", context)),
                    paths$modules$themes, paste0("GO_", ont, "_themes_", file_tag, ".svg"))

      # Leading-edge compact heatmap
      le_mat <- prep_leading_edge_matrix(gse, gene_list, top_terms = 5)
      save_plot_dir(plot_leading_edge_heatmap(le_mat, paste0("Leading edge genes (", ont, ") — ", context)),
                    dstP, paste0("gseGO_", ont, "_leading_edge_", file_tag, ".svg"))

      # Enrichmentmap-style visuals
      if (nrow(gse@result) >= 2) {
        ts_gse <- tryCatch(enrichplot::pairwise_termsim(gse), error = function(e) NULL)
        if (!is.null(ts_gse) && "compareClusterResult" %in% slotNames(ts_gse) &&
            nrow(ts_gse@compareClusterResult) > 0) {
          save_plot_dir(enrichplot::emapplot(ts_gse, showCategory = 10), dstP,
                        paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
        } else {
          writeLines("emapplot skipped: no enriched term after similarity calc.",
                     file.path(dstP, paste0("gseGO_", ont, "_emap_", file_tag, "_skipped.txt")))
        }
      }
      p_cnet <- tryCatch(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                         error = function(e) NULL)
      if (!is.null(p_cnet)) save_plot_dir(p_cnet, dstP, paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))

      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) { log_line(paths, file_tag, paste("simplify error:", e$message)); NULL })
      if (!is.null(gse2) && has_terms(gse2, 1)) {
        save_table(stable_rank(gse2@result), dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
        log_line(paths, file_tag, sprintf("simplify: %d -> %d terms",
          nrow(gse@result), nrow(gse2@result)))
      }
    } else {
      save_table(tibble::tibble(note = "gseGO empty or NULL"), dstT, paste0("gseGO_", ont, "_", file_tag, "_empty.csv"))
      writeLines("All plots skipped: no enriched term.", file.path(dstP, paste0("gseGO_", ont, "_", file_tag, "_skipped.txt")))
    }

    if (length(top_genes) >= 10) {
      ora <- with_timing({
        tryCatch(clusterProfiler::enrichGO(
          gene = top_genes, ont = ont, keyType = "UNIPROT",
          minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("ORA error:", e$message)); NULL })
      }, paths, file_tag, sprintf("ORA(%s)", ont))
      if (!is.null(ora) && has_terms(ora, 1)) {
        save_table(stable_rank(ora@result), dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        if (nrow(ora@result) >= 2) {
          ts_ora <- tryCatch(enrichplot::pairwise_termsim(ora), error = function(e) NULL)
          if (!is.null(ts_ora) && "compareClusterResult" %in% slotNames(ts_ora) &&
              nrow(ts_ora@compareClusterResult) > 0) {
            save_plot_dir(enrichplot::emapplot(ts_ora, showCategory = 10), dstP,
                          paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
          } else {
            writeLines("emapplot skipped: no enriched term after similarity calc.",
                       file.path(dstP, paste0("ORA_", ont, "_emap_", file_tag, "_skipped.txt")))
          }
        }
        p_cnet2 <- tryCatch(enrichplot::cnetplot(ora, showCategory = 5), error = function(e) NULL)
        if (!is.null(p_cnet2)) save_plot_dir(p_cnet2, dstP, paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
      } else {
        save_table(tibble::tibble(note = "ORA returned NULL or empty"), dstT, paste0("ORA_", ont, "_", file_tag, "_empty.csv"))
      }
    } else {
      save_table(tibble::tibble(note = "Global ORA skipped: <10 DEGs"),
                 dstT, paste0("ORA_", ont, "_", file_tag, "_skipped.csv"))
      log_line(paths, file_tag, paste("Global ORA skipped for", ont, "due to insufficient DEGs."))
    }
  }

  ## ----- GLOBAL: Reactome (optional) -----
  have_reactome <- requireNamespace("ReactomePA", quietly = TRUE)
  if (have_reactome) {
    # Map UniProt -> Entrez for Reactome
    entrez_map_re <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                           column = "ENTREZID", multiVals = "first")
    re_universe <- unname(entrez_map_re[!is.na(entrez_map_re)])
    re_deg <- unname(entrez_map_re[match(top_genes, names(entrez_map_re))])
    re_deg <- re_deg[!is.na(re_deg)]
    # ORA
    re_ora <- tryCatch(ReactomePA::enrichPathway(gene = re_deg, organism = "mouse",
                                                 universe = re_universe, pvalueCutoff = 1,
                                                 pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 800),
                       error = function(e) NULL)
    if (!is.null(re_ora) && has_terms(re_ora, 1)) {
      save_table(stable_rank(re_ora@result), paths$glob$Reactome$T, paste0("Reactome_ORA_", file_tag, ".csv"))
      save_plot_dir(plot_dot(re_ora, paste0("Reactome ORA — ", context)),
                    paths$glob$Reactome$P, paste0("Reactome_ORA_dotplot_", file_tag, ".svg"))
    }
    # GSEA
    re_gene_list <- gene_list[!is.na(entrez_map_re)]
    names(re_gene_list) <- unname(entrez_map_re[!is.na(entrez_map_re)])
    # de-dup
    re_gene_list <- re_gene_list[!duplicated(names(re_gene_list))]
    re_gene_list <- sort(re_gene_list, decreasing = TRUE, method = "radix")
    re_gse <- tryCatch(ReactomePA::gsePathway(geneList = re_gene_list, organism = "mouse",
                                              minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
                                              pAdjustMethod = "BH", eps = 0, nPermSimple = 10000),
                       error = function(e) NULL)
    if (!is.null(re_gse) && has_terms(re_gse, 1)) {
      save_table(stable_rank(re_gse@result), paths$glob$Reactome$T, paste0("Reactome_GSEA_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(re_gse, paste0("Reactome GSEA — ", context)),
                    paths$glob$Reactome$P, paste0("Reactome_GSEA_split_", file_tag, ".svg"))
      save_plot_dir(plot_gsea_landscape(re_gse, paste0("Reactome GSEA landscape — ", context)),
                    paths$glob$Reactome$P, paste0("Reactome_GSEA_landscape_", file_tag, ".svg"))
    }
  } else {
    log_line(paths, file_tag, "ReactomePA not available; Reactome enrichment skipped.")
  }

  ## ----- GLOBAL: KEGG (UNIPROT -> ENTREZ; de-dup; detailed log) -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  n_univ <- length(universe_genes)
  n_mapped <- sum(!is.na(entrez_map[universe_genes]))
  n_unmapped <- n_univ - n_mapped
  log_line(paths, file_tag, sprintf("KEGG mapIds: %d/%d mapped (%.1f%%), %d unmapped",
                                    n_mapped, n_univ, 100*n_mapped/n_univ, n_unmapped))

  # Diagnostics table
  kmap_tbl <- tibble::tibble(
    id = universe_genes,
    ENTREZID = unname(entrez_map[universe_genes]),
    in_kegg = !is.na(entrez_map[universe_genes])
  )
  save_table(kmap_tbl, paths$glob$KEGG$T, paste0("kegg_mapping_", file_tag, ".csv"))

  # Ranked stats named by ENTREZ, de-dup
  kegg_gene_list <- gene_list[!is.na(entrez_map)]
  names(kegg_gene_list) <- unname(entrez_map[!is.na(entrez_map)])
  dup_mask <- duplicated(names(kegg_gene_list))
  n_dups <- sum(dup_mask)
  if (n_dups > 0) {
    dup_ids <- unique(names(kegg_gene_list)[dup_mask])
    log_line(paths, file_tag, sprintf("gseKEGG: duplicate ENTREZ names=%d; examples: %s",
      n_dups, paste(utils::head(dup_ids, 10), collapse=", ")))
  }
  keep_idx <- !duplicated(names(kegg_gene_list))
  removed_dups <- sum(!keep_idx)
  kegg_gene_list <- kegg_gene_list[keep_idx]
  kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE, method = "radix")
  log_line(paths, file_tag, sprintf(
    "gseKEGG input: n=%d; anyDuplicated(names)=%d; pct ties=%.2f%%; removed_dups=%d",
    length(kegg_gene_list),
    anyDuplicated(names(kegg_gene_list)),
    100*mean(duplicated(kegg_gene_list)),
    removed_dups
  ))

  kk2 <- NULL
  if (length(kegg_gene_list) >= 10) {
    kk2 <- with_timing({
      muffle_fgsea_ties( tryCatch(clusterProfiler::gseKEGG(
        geneList = kegg_gene_list, organism = "mmu",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
        keyType = "ncbi-geneid"
      ), error = function(e) { log_line(paths, file_tag, paste("gseKEGG error:", e$message)); NULL }) )
    }, paths, file_tag, "gseKEGG")
  } else {
    log_line(paths, file_tag, "gseKEGG skipped: <10 unique ENTREZ in ranked list after de-dup.")
  }

  if (!is.null(kk2) && has_terms(kk2, 1)) {
    save_table(stable_rank(kk2@result), paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = 10, label_format = 30),
                  paths$glob$KEGG$P, paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    if (nrow(kk2@result) >= 2) {
      ts_kk <- tryCatch(enrichplot::pairwise_termsim(kk2), error = function(e) NULL)
      if (!is.null(ts_kk) && "compareClusterResult" %in% slotNames(ts_kk) &&
          nrow(ts_kk@compareClusterResult) > 0) {
        save_plot_dir(enrichplot::emapplot(ts_kk, showCategory = 10),
                      paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
      } else {
        writeLines("emapplot skipped: no enriched term after similarity calc.",
                   file.path(paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, "_skipped.txt")))
      }
    }
    save_plot_dir(plot_gsea_landscape(kk2, paste0("KEGG GSEA landscape — ", context)),
                  paths$glob$KEGG$P, paste0("gseKEGG_landscape_", file_tag, ".svg"))
    p_cnet3 <- tryCatch(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                        error = function(e) NULL)
    if (!is.null(p_cnet3)) save_plot_dir(p_cnet3, paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))
    rp2 <- tryCatch(enrichplot::ridgeplot(kk2) + ggplot2::labs(x = "Enrichment distribution"),
                    error = function(e) NULL)
    if (!is.null(rp2)) save_plot_dir(rp2, paths$glob$KEGG$P, paste0("gseKEGG_ridge_", file_tag, ".svg"))
  } else {
    save_table(tibble::tibble(note = "No KEGG enrichment results (gseKEGG null or empty)."),
               paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, "_empty.csv"))
    log_line(paths, file_tag, "gseKEGG produced no results.")
  }

  ## ----- PATHVIEW with validated KEGG IDs + status log -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && has_terms(kk2, 1)) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  save_table(tibble::tibble(pv_ids = pv_ids), paths$pathview, paste0("pathview_ids_", file_tag, ".csv"))

  if (length(kegg_gene_list) > 0 && all(is.finite(as.numeric(kegg_gene_list)))) {
    writeLines(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
    log_line(paths, file_tag, sprintf("Pathview attempted for %d pathways.", length(pv_ids)))
  } else {
    writeLines(c(
      "No valid ENTREZ-mapped genes for Pathview or gene list non-numeric.",
      paste0("length(kegg_gene_list)=", length(kegg_gene_list))
    ), con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Pathview skipped: invalid kegg_gene_list.")
  }

  ## ----- MODULE-AWARE: overlaps, UpSet, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")
  log_line(paths, file_tag, sprintf("Modules=%d; with DEG overlap>0: %d",
                                    length(modules_list), sum(overlap_counts$OverlapDEG>0)))

  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), ]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids

  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))
  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      library(ComplexUpset)
      set_cols <- setdiff(colnames(mat), "id")
      mat[set_cols] <- lapply(mat[set_cols], function(x) as.logical(x))
      p_up <- ComplexUpset::upset(
        data = mat,
        intersect = set_cols,
        base_annotations = list('Intersection size' = ComplexUpset::intersection_size()),
        set_sizes = FALSE,
        min_size = 5
      )
      save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      library(UpSetR)
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
    } else {
      save_table(tibble::tibble(note="No UpSet library available (ComplexUpset/UpSetR)."),
                 paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped_no_pkg.csv"))
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
               paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }

  # Overlap bar
  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, ]
  if (nrow(ov2) > 0) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2, 20)$OverlapDEG)*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  }

  # Fisher enrichment of DEGs in modules + plot
  fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
    deg_set <- intersect(unique(deg_genes), universe_genes)
    bg_set  <- setdiff(universe_genes, deg_set)
    out <- lapply(names(modules_list), function(m) {
      mod_genes <- intersect(modules_list[[m]], universe_genes)
      a <- length(intersect(mod_genes, deg_set))
      b <- length(setdiff(mod_genes, deg_set))
      c <- length(setdiff(deg_set, mod_genes))
      d <- length(setdiff(bg_set, mod_genes))
      mm <- matrix(c(a,b,c,d), nrow=2)
      ft <- fisher.test(mm, alternative = "greater")
      data.frame(
        Module = m, InModule_DE = a, InModule_NotDE = b,
        NotInModule_DE = c, NotInModule_NotDE = d,
        OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    res <- dplyr::bind_rows(out)
    res$Padj <- p.adjust(res$Pvalue, method = "BH")
    dplyr::arrange(res, Padj, Pvalue)
  }
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))
  if (nrow(fisher_tbl) > 0) {
    ft_plot <- ggplot2::ggplot(head(fisher_tbl[order(fisher_tbl$Padj, fisher_tbl$Pvalue), ], 20),
                               ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                            y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#6a51a3") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Fisher enrichment (DEGs in modules) — ", context),
                    x = "Module", y = "-log10(Padj)")
    save_plot_dir(ft_plot, paths$modules$fisher$P, paste0("fisher_bar_", file_tag, ".svg"), 24, 20)
  }

  ## ----- Per-module analyses -----
  sel_modules <- fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)]
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) {
      log_line(paths, file_tag, paste("Module", m, "skipped: <10 genes in universe."))
      next
    }

    for (ont in onts) {
      # Per-module nested folders
      m_root     <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea   <- mk(m_root, "Tables", "gsea")
      m_P_gsea   <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull<- mk(m_root, "Tables", "ora_full")
      m_P_orafull<- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]
      m_gene_list <- m_gene_list[order(m_gene_list, decreasing = TRUE, method = "radix")]
      if (length(m_gene_list) >= 10) {
        m_gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
          geneList = m_gene_list,
          ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "gseGO error:", e$message)); NULL }) )
        if (!is.null(m_gse) && has_terms(m_gse, 1)) {
          resm <- stable_rank(m_gse@result); resm$note <- "fgsea ties handled by deterministic id sort"
          save_table(resm, m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          p_dot <- tryCatch(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust",
                                                showCategory = 10, label_format = 30), error = function(e) NULL)
          if (!is.null(p_dot)) save_plot_dir(p_dot, m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "Module GSEA skipped: <10 genes."), m_T_gsea,
                   paste0("gseGO_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on DEG-intersection
      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA DEG error:", e$message)); NULL })
        if (!is.null(m_ora) && has_terms(m_ora, 1)) {
          save_table(stable_rank(m_ora@result), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(
          note = "insufficient DEGs for ORA",
          n_module_genes_in_universe = length(m_genes),
          n_deg_in_module = length(intersect(top_genes, m_genes)),
          deg_threshold = 1
        ), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      # ORA on full module gene set
      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA full error:", e$message)); NULL })
        if (!is.null(m_ora_all) && has_terms(m_ora_all, 1)) {
          save_table(stable_rank(m_ora_all@result), m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "ORA full module skipped: <5 genes."),
                   m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }
    } # ont
  } # module

  ## ----- EXECUTIVE SUMMARY (99_summary) -----
  # Collect top highlights
  summary_items <- list()
  # Universe/DEG
  summary_items$universe <- data.frame(
    metric = c("universe_size","deg_abs_gt1","pct_tied"),
    value  = c(length(universe_genes), length(top_genes), round(pct_tied,2))
  )
  # Top GO/KEGG/Reactome if available
  add_top <- function(path, pattern, key) {
    fp <- list.files(path, pattern = pattern, full.names = TRUE)
    if (length(fp)) {
      df <- tryCatch(readr::read_csv(fp[1], show_col_types = FALSE), error = function(e) NULL)
      if (!is.null(df) && nrow(df)>0) {
        df <- df[order(df$p.adjust %||% df$pvalue), ]
        summary_items[[key]] <- head(df[, c("ID","Description","p.adjust","pvalue")], 5)
      }
    }
  }
  `%||%` <- function(a,b) if (!is.null(a)) a else b
  add_top(paths$glob$BP$T, paste0("^gseGO_BP_", file_tag, "\\.csv$"), "top_go_gsea")
  add_top(paths$glob$KEGG$T, paste0("^gseKEGG_", file_tag, "\\.csv$"), "top_kegg_gsea")
  add_top(paths$glob$Reactome$T, paste0("^Reactome_GSEA_", file_tag, "\\.csv$"), "top_reactome_gsea")
  add_top(paths$modules$fisher$T, paste0("^fisher_", file_tag, "\\.csv$"), "top_modules_fisher")

  # Write a concise markdown summary
  md <- c(
    "# Executive Summary",
    paste0("- Context: ", context),
    paste0("- Universe size: ", summary_items$universe$value[summary_items$universe$metric=="universe_size"]),
    paste0("- DEGs |log2fc|>1: ", summary_items$universe$value[summary_items$universe$metric=="deg_abs_gt1"]),
    paste0("- pct tied stats: ", summary_items$universe$value[summary_items$universe$metric=="pct_tied"], "%"),
    "",
    "## Top GO GSEA (BP)",
    if (!is.null(summary_items$top_go_gsea)) paste0("* ", apply(summary_items$top_go_gsea, 1, function(r) sprintf("%s — %s (FDR=%.2g)", r[["ID"]], r[["Description"]], as.numeric(r[["p.adjust"]])))) else "No terms",
    "",
    "## Top KEGG GSEA",
    if (!is.null(summary_items$top_kegg_gsea)) paste0("* ", apply(summary_items$top_kegg_gsea, 1, function(r) sprintf("%s — %s (FDR=%.2g)", r[["ID"]], r[["Description"]], as.numeric(r[["p.adjust"]])))) else "No terms",
    "",
    "## Top Reactome GSEA",
    if (!is.null(summary_items$top_reactome_gsea)) paste0("* ", apply(summary_items$top_reactome_gsea, 1, function(r) sprintf("%s — %s (FDR=%.2g)", r[["ID"]], r[["Description"]], as.numeric(r[["p.adjust"]])))) else "Skipped/None",
    "",
    "## Top Modules by Fisher (Padj)",
    if (!is.null(summary_items$top_modules_fisher)) paste0("* ", apply(summary_items$top_modules_fisher, 1, function(r) sprintf("%s — OR=%.2f (Padj=%.2g)", r[["Module"]], as.numeric(r[["OddsRatio"]]), as.numeric(r[["Padj"]])))) else "No modules"
  )
  md <- unlist(md)
  writeLines(md, con = file.path(paths$summary, paste0("SUMMARY_", file_tag, ".md")))

  writeLines(paste("End:", Sys.time()), con = file.path(paths$meta, "run_end.txt"))
  log_line(paths, file_tag, "Completed all analyses.")
}





































############################################################
# clusterProfiler + WGCNA modules (mouse) — extended figures
# Adds: robust executive summary with guards against NULLs,
# defensive checks before order()/sort()/cluster calls.
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","vroom","UniProt.ws","BiocParallel",
    "BiocGenerics","S4Vectors","IRanges","GenomeInfoDb","GenomicRanges","SummarizedExperiment",
    "Biobase","GO.db","UpSetR","ComplexUpset","igraph","ggraph"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Directory + I/O helpers ----------
mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP   = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                  P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF   = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                  P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC   = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                  P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                  P = mk(results_dir, "10_global","KEGG","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots"))
    ),
    pathview = mk(results_dir, "90_pathview"),
    summary = mk(results_dir, "99_summary")
  )
}
save_table <- function(df, dirpath, fname) {
  readr::write_csv(df, file.path(dirpath, fname))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}
log_line <- function(paths, tag, text) {
  lf <- file.path(paths$meta, paste0("runlog_", tag, ".txt"))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text), file = lf, append = TRUE)
}
stable_rank <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  if (!("pvalue" %in% names(df)) && ("p.adjust" %in% names(df))) df$pvalue <- df$p.adjust
  if (!("p.adjust" %in% names(df)) && ("pvalue" %in% names(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df$rank_p <- rank(df$pvalue, ties.method = "min")
  df$rank_padj <- rank(df$p.adjust, ties.method = "min")
  df[order(df$p.adjust, df$pvalue, df$rank_padj, df$rank_p, decreasing = FALSE), ]
}
has_terms <- function(x, min_n = 1) {
  tryCatch({
    df <- as.data.frame(x)
    !is.null(df) && nrow(df) >= min_n
  }, error = function(e) FALSE)
}
muffle_fgsea_ties <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("There are ties in the preranked stats", conditionMessage(m))) {
      invokeRestart("muffleMessage")
    }
  })
}

# Defensive helpers
is_nonempty_num <- function(x) is.numeric(x) && length(x) > 0
is_nonempty_char <- function(x) is.character(x) && length(x) > 0

with_timing <- function(expr, paths, tag, label) {
  t0 <- Sys.time()
  res <- NULL; err <- NULL
  tryCatch({ res <- eval.parent(substitute(expr)) },
           error = function(e) { err <<- e })
  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  log_line(paths, tag, sprintf("%s elapsed=%.2fs %s",
                               label, elapsed,
                               if (!is.null(err)) paste("| error:", err$message) else "" ))
  if (!is.null(err)) stop(err)
  res
}

## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text  = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}

plot_gsea_landscape <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj) || !has_terms(gsea_obj, 1)) return(NULL)
  df <- as.data.frame(gsea_obj@result)
  if (!nrow(df) || !"NES" %in% names(df)) return(NULL)
  df$mlogFDR <- -log10(pmax(1e-300, df$p.adjust))
  ggplot2::ggplot(df, ggplot2::aes(x = NES, y = mlogFDR)) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(title = title_txt, x = "NES", y = "-log10(FDR)") +
    ggplot2::theme_minimal()
}

prep_leading_edge_matrix <- function(gsea_obj, gene_list, top_terms = 5) {
  if (is.null(gsea_obj) || !has_terms(gsea_obj, 1)) return(NULL)
  res <- as.data.frame(gsea_obj@result)
  if (!nrow(res)) return(NULL)
  # guard order() inputs
  if (!("p.adjust" %in% names(res)) || !("NES" %in% names(res))) return(NULL)
  res <- res[order(res$p.adjust, -abs(res$NES)), , drop = FALSE]
  res <- head(res, top_terms)
  if (!nrow(res)) return(NULL)
  if (!"core_enrichment" %in% colnames(res)) return(NULL)
  le_genes <- unique(unlist(strsplit(paste(res$core_enrichment, collapse = "/"), "/")))
  le_genes <- intersect(le_genes, names(gene_list))
  if (!length(le_genes)) return(NULL)
  mat <- matrix(gene_list[le_genes], nrow = length(le_genes), dimnames = list(le_genes, "log2fc"))
  mat
}

plot_leading_edge_heatmap <- function(mat, title_txt) {
  if (is.null(mat)) return(NULL)
  df <- data.frame(Gene = rownames(mat), log2fc = mat[,1], row.names = NULL, stringsAsFactors = FALSE)
  if (!nrow(df)) return(NULL)
  df$Gene <- factor(df$Gene, levels = df$Gene[order(df$log2fc, decreasing = TRUE)])
  ggplot2::ggplot(df, ggplot2::aes(x = "Leading edge", y = Gene, fill = log2fc)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    ggplot2::labs(title = title_txt, x = NULL, y = "Genes (leading edge)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
}

## ---------- Load WGCNA modules ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "") |>
  dplyr::mutate(Gene = as.character(Gene))

modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(x))
module_sizes <- tibble::tibble(
  Module = names(modules_list),
  N_accessions = vapply(modules_list, length, integer(1))
)

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Inputs ----------
comparison <- "neuron-phenotypeWithinUnit"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Provenance snapshot
  capture.output(sessionInfo(), file = file.path(paths$meta, "sessionInfo.txt"))
  writeLines(paste("Start:", Sys.time()), con = file.path(paths$meta, "run_start.txt"))
  save_table(module_sizes, paths$modules$root, "module_sizes.csv")

  # Parameters README
  writeLines(c(
    "# Analysis README",
    paste0("Comparison: ", comparison),
    paste0("Context: ", context),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (mapIds from org.Mm.eg.db, multiVals='first')",
    "DEG threshold for ORA: |log2fc| > 1",
    "Ontologies default: BP",
    "ComplexUpset: set sizes panel disabled for portability",
    "GSEA ties: sorted by id for deterministic order; values are NOT jittered",
    "Folders:",
    "  10_global/GO_*: GSEA and ORA across full universe",
    "  10_global/KEGG: gseKEGG results and mapping diagnostics",
    "  20_modules/: overlap, Fisher, per-module GO (GSEA/ORA)",
    "  90_pathview/: KEGG pathway renderings and status logs",
    "  99_summary/: executive summary"
  ), con = file.path(results_dir, "README.txt"))
  log_line(paths, file_tag, "Initialized results folders and README.")

  ## ----- Load comparison -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  gene_list_tbl <- df |> dplyr::select(id, log2fc) |> dplyr::arrange(dplyr::desc(log2fc), id)
  gene_list <- gene_list_tbl$log2fc; names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)
  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  pct_tied <- 100 * mean(duplicated(gene_list))
  save_table(tibble::tibble(
    n_rows = nrow(df),
    n_universe = length(universe_genes),
    n_deg_abs1 = length(top_genes),
    deg_threshold = 1,
    pct_tied_stats = pct_tied
  ), paths$meta, paste0("summary_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("Universe=%d; DEGs(|log2fc|>1)=%d; tied stats=%.2f%%",
                                    length(universe_genes), length(top_genes), pct_tied))

  organism <- "org.Mm.eg.db"
  onts <- c("BP")

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    gse <- with_timing({
      muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
        geneList = gene_list, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        by = "fgsea", nPermSimple = 10000, eps = 0,
        verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
      ), error = function(e) { log_line(paths, file_tag, paste("gseGO error:", e$message)); NULL }) )
    }, paths, file_tag, sprintf("gseGO(%s)", ont))

    if (!is.null(gse) && has_terms(gse, 1)) {
      res_stable <- stable_rank(gse@result)
      res_stable$note <- "fgsea ties handled by deterministic id sort (values unchanged)"
      save_table(res_stable, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))
      save_plot_dir(plot_gsea_landscape(gse, paste0("GSEA landscape (", ont, ") — ", context)),
                    dstP, paste0("gseGO_", ont, "_landscape_", file_tag, ".svg"))

      le_mat <- prep_leading_edge_matrix(gse, gene_list, top_terms = 5)
      save_plot_dir(plot_leading_edge_heatmap(le_mat, paste0("Leading edge genes (", ont, ") — ", context)),
                    dstP, paste0("gseGO_", ont, "_leading_edge_", file_tag, ".svg"))

      if (nrow(gse@result) >= 2) {
        ts_gse <- tryCatch(enrichplot::pairwise_termsim(gse), error = function(e) NULL)
        if (!is.null(ts_gse) && "compareClusterResult" %in% slotNames(ts_gse) &&
            nrow(ts_gse@compareClusterResult) > 0) {
          save_plot_dir(enrichplot::emapplot(ts_gse, showCategory = 10), dstP,
                        paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
        } else {
          writeLines("emapplot skipped: no enriched term after similarity calc.",
                     file.path(dstP, paste0("gseGO_", ont, "_emap_", file_tag, "_skipped.txt")))
        }
      }

      p_cnet <- tryCatch(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                         error = function(e) NULL)
      if (!is.null(p_cnet)) save_plot_dir(p_cnet, dstP, paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))

      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) { log_line(paths, file_tag, paste("simplify error:", e$message)); NULL })
      if (!is.null(gse2) && has_terms(gse2, 1)) {
        save_table(stable_rank(gse2@result), dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
      }
    } else {
      save_table(tibble::tibble(note = "gseGO empty or NULL"), dstT, paste0("gseGO_", ont, "_", file_tag, "_empty.csv"))
      writeLines("All plots skipped: no enriched term.", file.path(dstP, paste0("gseGO_", ont, "_", file_tag, "_skipped.txt")))
    }

    if (length(top_genes) >= 10) {
      ora <- with_timing({
        tryCatch(clusterProfiler::enrichGO(
          gene = top_genes, ont = ont, keyType = "UNIPROT",
          minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("ORA error:", e$message)); NULL })
      }, paths, file_tag, sprintf("ORA(%s)", ont))
      if (!is.null(ora) && has_terms(ora, 1)) {
        save_table(stable_rank(ora@result), dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        if (nrow(ora@result) >= 2) {
          ts_ora <- tryCatch(enrichplot::pairwise_termsim(ora), error = function(e) NULL)
          if (!is.null(ts_ora) && "compareClusterResult" %in% slotNames(ts_ora) &&
              nrow(ts_ora@compareClusterResult) > 0) {
            save_plot_dir(enrichplot::emapplot(ts_ora, showCategory = 10), dstP,
                          paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
          } else {
            writeLines("emapplot skipped: no enriched term after similarity calc.",
                       file.path(dstP, paste0("ORA_", ont, "_emap_", file_tag, "_skipped.txt")))
          }
        }
        p_cnet2 <- tryCatch(enrichplot::cnetplot(ora, showCategory = 5), error = function(e) NULL)
        if (!is.null(p_cnet2)) save_plot_dir(p_cnet2, dstP, paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
      } else {
        save_table(tibble::tibble(note = "ORA returned NULL or empty"), dstT, paste0("ORA_", ont, "_", file_tag, "_empty.csv"))
      }
    } else {
      save_table(tibble::tibble(note = "Global ORA skipped: <10 DEGs"),
                 dstT, paste0("ORA_", ont, "_", file_tag, "_skipped.csv"))
      log_line(paths, file_tag, paste("Global ORA skipped for", ont, "due to insufficient DEGs."))
    }
  }

  ## ----- GLOBAL: KEGG with de-dup -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  kmap_tbl <- tibble::tibble(
    id = universe_genes,
    ENTREZID = unname(entrez_map[universe_genes]),
    in_kegg = !is.na(entrez_map[universe_genes])
  )
  save_table(kmap_tbl, paths$glob$KEGG$T, paste0("kegg_mapping_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("KEGG mapIds: %d/%d mapped (%.1f%%)",
                                    sum(kmap_tbl$in_kegg), nrow(kmap_tbl), 100*mean(kmap_tbl$in_kegg)))

  kegg_gene_list <- gene_list[!is.na(entrez_map)]
  names(kegg_gene_list) <- unname(entrez_map[!is.na(entrez_map)])
  keep_idx <- !duplicated(names(kegg_gene_list))
  removed_dups <- sum(!keep_idx)
  kegg_gene_list <- kegg_gene_list[keep_idx]
  kegg_gene_list <- kegg_gene_list[order(kegg_gene_list, decreasing = TRUE, method = "radix")]
  log_line(paths, file_tag, sprintf("gseKEGG input: n=%d; removed_dups=%d; anyDuplicated=%d; pct ties=%.2f%%",
                                    length(kegg_gene_list), removed_dups, anyDuplicated(names(kegg_gene_list)),
                                    100*mean(duplicated(kegg_gene_list))))

  kk2 <- NULL
  if (length(kegg_gene_list) >= 10) {
    kk2 <- with_timing({
      muffle_fgsea_ties( tryCatch(clusterProfiler::gseKEGG(
        geneList = kegg_gene_list, organism = "mmu",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
        keyType = "ncbi-geneid"
      ), error = function(e) { log_line(paths, file_tag, paste("gseKEGG error:", e$message)); NULL }) )
    }, paths, file_tag, "gseKEGG")
  } else {
    log_line(paths, file_tag, "gseKEGG skipped: <10 unique ENTREZ in ranked list after de-dup.")
  }

  if (!is.null(kk2) && has_terms(kk2, 1)) {
    save_table(stable_rank(kk2@result), paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = 10, label_format = 30),
                  paths$glob$KEGG$P, paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    if (nrow(kk2@result) >= 2) {
      ts_kk <- tryCatch(enrichplot::pairwise_termsim(kk2), error = function(e) NULL)
      if (!is.null(ts_kk) && "compareClusterResult" %in% slotNames(ts_kk) &&
          nrow(ts_kk@compareClusterResult) > 0) {
        save_plot_dir(enrichplot::emapplot(ts_kk, showCategory = 10),
                      paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
      } else {
        writeLines("emapplot skipped: no enriched term after similarity calc.",
                   file.path(paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, "_skipped.txt")))
      }
    }
    save_plot_dir(plot_gsea_landscape(kk2, paste0("KEGG GSEA landscape — ", context)),
                  paths$glob$KEGG$P, paste0("gseKEGG_landscape_", file_tag, ".svg"))
    p_cnet3 <- tryCatch(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                        error = function(e) NULL)
    if (!is.null(p_cnet3)) save_plot_dir(p_cnet3, paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))
  } else {
    save_table(tibble::tibble(note = "No KEGG enrichment results (gseKEGG null or empty)."),
               paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, "_empty.csv"))
    log_line(paths, file_tag, "gseKEGG produced no results.")
  }

  ## ----- PATHVIEW -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && has_terms(kk2, 1)) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  save_table(tibble::tibble(pv_ids = pv_ids), paths$pathview, paste0("pathview_ids_", file_tag, ".csv"))

  if (length(kegg_gene_list) > 0 && all(is.finite(as.numeric(kegg_gene_list)))) {
    writeLines(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
    log_line(paths, file_tag, sprintf("Pathview attempted for %d pathways.", length(pv_ids)))
  } else {
    writeLines(c(
      "No valid ENTREZ-mapped genes for Pathview or gene list non-numeric.",
      paste0("length(kegg_gene_list)=", length(kegg_gene_list))
    ), con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Pathview skipped: invalid kegg_gene_list.")
  }

  ## ----- MODULE-AWARE: overlaps, UpSet, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")

  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), , drop = FALSE]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids

  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))
  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      library(ComplexUpset)
      set_cols <- setdiff(colnames(mat), "id")
      mat[set_cols] <- lapply(mat[set_cols], function(x) as.logical(x))
      p_up <- ComplexUpset::upset(
        data = mat,
        intersect = set_cols,
        base_annotations = list('Intersection size' = ComplexUpset::intersection_size()),
        set_sizes = FALSE,
        min_size = 5
      )
      save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      library(UpSetR)
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
    } else {
      save_table(tibble::tibble(note="No UpSet library available (ComplexUpset/UpSetR)."),
                 paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped_no_pkg.csv"))
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
               paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }

  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, , drop = FALSE]
  if (nrow(ov2) > 0 && is_nonempty_num(ov2$OverlapDEG)) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2, 20)$OverlapDEG)*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No non-zero overlaps to plot.", con = file.path(paths$modules$overlap$P, paste0("overlap_status_", file_tag, ".txt")))
  }

  fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
    deg_set <- intersect(unique(deg_genes), universe_genes)
    bg_set  <- setdiff(universe_genes, deg_set)
    out <- lapply(names(modules_list), function(m) {
      mod_genes <- intersect(modules_list[[m]], universe_genes)
      a <- length(intersect(mod_genes, deg_set))
      b <- length(setdiff(mod_genes, deg_set))
      c <- length(setdiff(deg_set, mod_genes))
      d <- length(setdiff(bg_set, mod_genes))
      mm <- matrix(c(a,b,c,d), nrow=2)
      ft <- fisher.test(mm, alternative = "greater")
      data.frame(
        Module = m, InModule_DE = a, InModule_NotDE = b,
        NotInModule_DE = c, NotInModule_NotDE = d,
        OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    res <- dplyr::bind_rows(out)
    res$Padj <- p.adjust(res$Pvalue, method = "BH")
    dplyr::arrange(res, Padj, Pvalue)
  }
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))

  if (nrow(fisher_tbl) > 0 && is_nonempty_num(fisher_tbl$Padj)) {
    ord_idx <- order(fisher_tbl$Padj, fisher_tbl$Pvalue)
    ft_plot <- ggplot2::ggplot(head(fisher_tbl[ord_idx, , drop = FALSE], 20),
                               ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                            y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#6a51a3") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Fisher enrichment (DEGs in modules) — ", context),
                    x = "Module", y = "-log10(Padj)")
    save_plot_dir(ft_plot, paths$modules$fisher$P, paste0("fisher_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No Fisher rows to plot.", con = file.path(paths$modules$fisher$P, paste0("fisher_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Fisher plot skipped: empty table.")
  }

  ## ----- Per-module analyses -----
  sel_modules <- if (nrow(fisher_tbl) > 0) fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)] else character(0)
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) {
      log_line(paths, file_tag, paste("Module", m, "skipped: <10 genes in universe."))
      next
    }

    for (ont in onts) {
      m_root     <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea   <- mk(m_root, "Tables", "gsea")
      m_P_gsea   <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull<- mk(m_root, "Tables", "ora_full")
      m_P_orafull<- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]
      m_gene_list <- m_gene_list[order(m_gene_list, decreasing = TRUE, method = "radix")]
      if (length(m_gene_list) >= 10) {
        m_gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
          geneList = m_gene_list, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "gseGO error:", e$message)); NULL }) )
        if (!is.null(m_gse) && has_terms(m_gse, 1)) {
          resm <- stable_rank(m_gse@result); resm$note <- "fgsea ties handled by deterministic id sort"
          save_table(resm, m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          p_dot <- tryCatch(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust",
                                                showCategory = 10, label_format = 30), error = function(e) NULL)
          if (!is.null(p_dot)) save_plot_dir(p_dot, m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "Module GSEA skipped: <10 genes."), m_T_gsea,
                   paste0("gseGO_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA DEG error:", e$message)); NULL })
        if (!is.null(m_ora) && has_terms(m_ora, 1)) {
          save_table(stable_rank(m_ora@result), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(
          note = "insufficient DEGs for ORA",
          n_module_genes_in_universe = length(m_genes),
          n_deg_in_module = length(intersect(top_genes, m_genes)),
          deg_threshold = 1
        ), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA full error:", e$message)); NULL })
        if (!is.null(m_ora_all) && has_terms(m_ora_all, 1)) {
          save_table(stable_rank(m_ora_all@result), m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "ORA full module skipped: <5 genes."),
                   m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }
    } # ont
  } # module

  ## ----- EXECUTIVE SUMMARY (robust) -----
  summary_items <- list()
  summary_items$universe <- data.frame(
    metric = c("universe_size","deg_abs_gt1","pct_tied"),
    value  = c(length(universe_genes), length(top_genes), round(pct_tied,2))
  )
  # Top GO/KEGG tables if they exist
  get_first_csv <- function(path, pattern) {
    fps <- list.files(path, pattern = pattern, full.names = TRUE)
    if (length(fps)) tryCatch(readr::read_csv(fps[1], show_col_types = FALSE), error = function(e) NULL) else NULL
  }
  summary_items$top_go_gsea <- get_first_csv(paths$glob$BP$T, paste0("^gseGO_BP_", file_tag, "\\.csv$"))
  summary_items$top_kegg_gsea <- get_first_csv(paths$glob$KEGG$T, paste0("^gseKEGG_", file_tag, "\\.csv$"))
  summary_items$top_modules_fisher <- {
    ft <- tryCatch(readr::read_csv(file.path(paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv")), show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(ft)) ft[order(ft$Padj, ft$Pvalue), , drop = FALSE] else NULL
  }

  safe_univ_val <- function(x, metric, default = NA_character_) {
    if (is.null(x) || !nrow(x) || !"metric" %in% names(x) || !"value" %in% names(x)) return(default)
    v <- x$value[x$metric == metric]
    if (length(v) == 0) return(default)
    as.character(v[[1]])
  }
  fmt_top <- function(df, fmt) {
    if (is.null(df) || !nrow(df)) return(character(0))
    need <- intersect(c("ID","Description","p.adjust"), names(df))
    if (length(need) < 3) return(character(0))
    df <- df[order(df$p.adjust), , drop = FALSE]
    apply(head(df, 5), 1, function(r) sprintf(fmt, r[["ID"]], r[["Description"]], as.numeric(r[["p.adjust"]])))
  }

  md_lines <- list(
    "# Executive Summary",
    paste0("- Context: ", context),
    paste0("- Universe size: ", safe_univ_val(summary_items$universe, "universe_size", "NA")),
    paste0("- DEGs |log2fc|>1: ", safe_univ_val(summary_items$universe, "deg_abs_gt1", "NA")),
    paste0("- pct tied stats: ", safe_univ_val(summary_items$universe, "pct_tied", "NA"), "%"),
    "",
    "## Top GO GSEA (BP)"
  )
  go_lines <- fmt_top(summary_items$top_go_gsea, "%s — %s (FDR=%.2g)")
  md_lines <- c(md_lines, if (length(go_lines)) paste0("* ", go_lines) else "No terms", "", "## Top KEGG GSEA")
  kegg_lines <- fmt_top(summary_items$top_kegg_gsea, "%s — %s (FDR=%.2g)")
  md_lines <- c(md_lines, if (length(kegg_lines)) paste0("* ", kegg_lines) else "No terms", "", "## Top Reactome GSEA")
  # Reactome optional: skip here to keep script self-contained, or add if collected
  md_lines <- c(md_lines, "Skipped/None", "", "## Top Modules by Fisher (Padj)")
  if (!is.null(summary_items$top_modules_fisher) && nrow(summary_items$top_modules_fisher)) {
    mod_tbl <- head(summary_items$top_modules_fisher, 5)
    if (all(c("Module","OddsRatio","Padj") %in% names(mod_tbl))) {
      mod_lines <- apply(mod_tbl, 1, function(r) sprintf("%s — OR=%.2f (Padj=%.2g)", r[["Module"]], as.numeric(r[["OddsRatio"]]), as.numeric(r[["Padj"]])))
      md_lines <- c(md_lines, paste0("* ", mod_lines))
    } else {
      md_lines <- c(md_lines, "No modules (columns missing).")
    }
  } else {
    md_lines <- c(md_lines, "No modules")
  }
  md <- unlist(md_lines, use.names = FALSE)
  writeLines(md, con = file.path(paths$summary, paste0("SUMMARY_", file_tag, ".md")))

  writeLines(paste("End:", Sys.time()), con = file.path(paths$meta, "run_end.txt"))
  log_line(paths, file_tag, "Completed all analyses.")
}




































############################################################
# clusterProfiler + WGCNA modules (mouse) — extended figures
# Adds: robust executive summary with guards against NULLs,
# defensive checks before order()/sort()/cluster calls.
# Restored: Reactome (ORA+GSEA), themes_networks, richer KEGG
# diagnostics incl. ridgeplot, env.txt provenance, README + summary hooks.
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","vroom","UniProt.ws","BiocParallel",
    "BiocGenerics","S4Vectors","IRanges","GenomeInfoDb","GenomicRanges","SummarizedExperiment",
    "Biobase","GO.db","UpSetR","ComplexUpset","igraph","ggraph"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST","ReactomePA")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Directory + I/O helpers ----------
mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP    = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                   P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF    = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                   P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC    = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                   P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG  = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                   P = mk(results_dir, "10_global","KEGG","Plots")),
      Reactome = list(T = mk(results_dir, "10_global","Reactome","Tables"),
                      P = mk(results_dir, "10_global","Reactome","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots")),
      themes  = mk(results_dir, "20_modules","themes_networks")
    ),
    pathview = mk(results_dir, "90_pathview"),
    summary  = mk(results_dir, "99_summary")
  )
}
save_table <- function(df, dirpath, fname) {
  readr::write_csv(df, file.path(dirpath, fname))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height)
}
log_line <- function(paths, tag, text) {
  lf <- file.path(paths$meta, paste0("runlog_", tag, ".txt"))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text), file = lf, append = TRUE)
}
stable_rank <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  if (!("pvalue" %in% names(df)) && ("p.adjust" %in% names(df))) df$pvalue <- df$p.adjust
  if (!("p.adjust" %in% names(df)) && ("pvalue" %in% names(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df$rank_p <- rank(df$pvalue, ties.method = "min")
  df$rank_padj <- rank(df$p.adjust, ties.method = "min")
  df[order(df$p.adjust, df$pvalue, df$rank_padj, df$rank_p, decreasing = FALSE), ]
}
has_terms <- function(x, min_n = 1) {
  tryCatch({
    df <- as.data.frame(x)
    !is.null(df) && nrow(df) >= min_n
  }, error = function(e) FALSE)
}
muffle_fgsea_ties <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("There are ties in the preranked stats", conditionMessage(m))) {
      invokeRestart("muffleMessage")
    }
  })
}

# Defensive helpers
is_nonempty_num <- function(x) is.numeric(x) && length(x) > 0
is_nonempty_char <- function(x) is.character(x) && length(x) > 0

with_timing <- function(expr, paths, tag, label) {
  t0 <- Sys.time()
  res <- NULL; err <- NULL
  tryCatch({ res <- eval.parent(substitute(expr)) },
           error = function(e) { err <<- e })
  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  log_line(paths, tag, sprintf("%s elapsed=%.2fs %s",
                               label, elapsed,
                               if (!is.null(err)) paste("| error:", err$message) else "" ))
  if (!is.null(err)) stop(err)
  res
}

## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = 10, split = ".sign", font.size = 12, label_format = 30
  ) +
    ggplot2::facet_wrap(~ .sign, nrow = 1) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      strip.text  = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = 10, font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
}
plot_gsea_landscape <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj) || !has_terms(gsea_obj, 1)) return(NULL)
  df <- as.data.frame(gsea_obj@result)
  if (!nrow(df) || !"NES" %in% names(df)) return(NULL)
  df$mlogFDR <- -log10(pmax(1e-300, df$p.adjust))
  ggplot2::ggplot(df, ggplot2::aes(x = NES, y = mlogFDR)) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(title = title_txt, x = "NES", y = "-log10(FDR)") +
    ggplot2::theme_minimal()
}
prep_leading_edge_matrix <- function(gsea_obj, gene_list, top_terms = 5) {
  if (is.null(gsea_obj) || !has_terms(gsea_obj, 1)) return(NULL)
  res <- as.data.frame(gsea_obj@result)
  if (!nrow(res)) return(NULL)
  if (!("p.adjust" %in% names(res)) || !("NES" %in% names(res))) return(NULL)
  res <- res[order(res$p.adjust, -abs(res$NES)), , drop = FALSE]
  res <- head(res, top_terms)
  if (!nrow(res)) return(NULL)
  if (!"core_enrichment" %in% colnames(res)) return(NULL)
  le_genes <- unique(unlist(strsplit(paste(res$core_enrichment, collapse = "/"), "/")))
  le_genes <- intersect(le_genes, names(gene_list))
  if (!length(le_genes)) return(NULL)
  mat <- matrix(gene_list[le_genes], nrow = length(le_genes), dimnames = list(le_genes, "log2fc"))
  mat
}
plot_leading_edge_heatmap <- function(mat, title_txt) {
  if (is.null(mat)) return(NULL)
  df <- data.frame(Gene = rownames(mat), log2fc = mat[,1], row.names = NULL, stringsAsFactors = FALSE)
  if (!nrow(df)) return(NULL)
  df$Gene <- factor(df$Gene, levels = df$Gene[order(df$log2fc, decreasing = TRUE)])
  ggplot2::ggplot(df, ggplot2::aes(x = "Leading edge", y = Gene, fill = log2fc)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    ggplot2::labs(title = title_txt, x = NULL, y = "Genes (leading edge)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
}

## ---------- Theme graph helpers (restored) ----------
build_theme_graph <- function(enrich_obj, top_n = 30, sim_cut = 0.35) {
  if (is.null(enrich_obj) || !has_terms(enrich_obj, 2)) return(NULL)
  ts <- tryCatch(enrichplot::pairwise_termsim(enrich_obj), error = function(e) NULL)
  if (is.null(ts) || !"compareClusterResult" %in% slotNames(ts)) return(NULL)
  dfr <- tryCatch(as.data.frame(ts@compareClusterResult), error = function(e) NULL)
  if (is.null(dfr) || !nrow(dfr)) return(NULL)
  dfr <- dfr[order(dfr$p.adjust, dfr$pvalue), , drop = FALSE]
  ids <- head(dfr$ID, top_n)
  sim <- ts@termsim[rownames(ts@termsim) %in% ids, colnames(ts@termsim) %in% ids, drop = FALSE]
  if (is.null(sim) || !nrow(sim) || !ncol(sim)) return(NULL)
  ed <- as.data.frame(as.table(sim), stringsAsFactors = FALSE)
  ed <- ed[ed$Freq >= sim_cut & ed$Var1 != ed$Var2, , drop = FALSE]
  if (!nrow(ed)) return(NULL)
  g <- igraph::graph_from_data_frame(ed, directed = FALSE)
  igraph::V(g)$label <- dfr$Description[match(igraph::V(g)$name, dfr$ID)]
  igraph::V(g)$size  <- -log10(pmax(1e-300, dfr$p.adjust[match(igraph::V(g)$name, dfr$ID)]))
  g
}
plot_theme_graph <- function(g, title_txt = "Enrichment themes") {
  if (is.null(g)) return(NULL)
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.25) +
    ggraph::geom_node_point(ggplot2::aes(size = size, color = size)) +
    ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE, size = 3) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title_txt, color = "-log10(FDR)", size = "-log10(FDR)")
}

## ---------- Load WGCNA modules ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "") |>
  dplyr::mutate(Gene = as.character(Gene))

modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(x))
module_sizes <- tibble::tibble(
  Module = names(modules_list),
  N_accessions = vapply(modules_list, length, integer(1))
)

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Inputs ----------
comparison <- "neuron-phenotypeWithinUnit"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Provenance snapshot
  capture.output(sessionInfo(), file = file.path(paths$meta, "sessionInfo.txt"))
  writeLines(capture.output(Sys.getenv()), con = file.path(paths$meta, "env.txt"))
  writeLines(paste("Start:", Sys.time()), con = file.path(paths$meta, "run_start.txt"))
  save_table(module_sizes, paths$modules$root, "module_sizes.csv")

  # Parameters README
  writeLines(c(
    "# Analysis README",
    paste0("Comparison: ", comparison),
    paste0("Context: ", context),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (mapIds from org.Mm.eg.db, multiVals='first')",
    "DEG threshold for ORA: |log2fc| > 1",
    "Ontologies default: BP",
    "ComplexUpset: set sizes panel disabled for portability",
    "GSEA ties: sorted by id for deterministic order; values are NOT jittered",
    "Folders:",
    "  10_global/GO_*: GSEA and ORA across full universe",
    "  10_global/KEGG: gseKEGG results and mapping diagnostics",
    "  10_global/Reactome: optional Reactome enrichment",
    "  20_modules/: overlap, Fisher, per-module GO (GSEA/ORA)",
    "  20_modules/themes_networks: similarity-based term theme graphs",
    "  90_pathview/: KEGG pathway renderings and status logs",
    "  99_summary/: executive summary"
  ), con = file.path(results_dir, "README.txt"))
  log_line(paths, file_tag, "Initialized results folders and README.")

  ## ----- Load comparison -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  gene_list_tbl <- df |> dplyr::select(id, log2fc) |> dplyr::arrange(dplyr::desc(log2fc), id)
  gene_list <- gene_list_tbl$log2fc; names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)
  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  pct_tied <- 100 * mean(duplicated(gene_list))
  save_table(tibble::tibble(
    n_rows = nrow(df),
    n_universe = length(universe_genes),
    n_deg_abs1 = length(top_genes),
    deg_threshold = 1,
    pct_tied_stats = pct_tied
  ), paths$meta, paste0("summary_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("Universe=%d; DEGs(|log2fc|>1)=%d; tied stats=%.2f%%",
                                    length(universe_genes), length(top_genes), pct_tied))

  organism <- "org.Mm.eg.db"
  onts <- c("BP")

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    gse <- with_timing({
      muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
        geneList = gene_list, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        by = "fgsea", nPermSimple = 10000, eps = 0,
        verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
      ), error = function(e) { log_line(paths, file_tag, paste("gseGO error:", e$message)); NULL }) )
    }, paths, file_tag, sprintf("gseGO(%s)", ont))

    if (!is.null(gse) && has_terms(gse, 1)) {
      res_stable <- stable_rank(gse@result)
      res_stable$note <- "fgsea ties handled by deterministic id sort (values unchanged)"
      save_table(res_stable, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))
      save_plot_dir(plot_gsea_landscape(gse, paste0("GSEA landscape (", ont, ") — ", context)),
                    dstP, paste0("gseGO_", ont, "_landscape_", file_tag, ".svg"))

      le_mat <- prep_leading_edge_matrix(gse, gene_list, top_terms = 5)
      save_plot_dir(plot_leading_edge_heatmap(le_mat, paste0("Leading edge genes (", ont, ") — ", context)),
                    dstP, paste0("gseGO_", ont, "_leading_edge_", file_tag, ".svg"))

      if (nrow(gse@result) >= 2) {
        ts_gse <- tryCatch(enrichplot::pairwise_termsim(gse), error = function(e) NULL)
        if (!is.null(ts_gse) && "compareClusterResult" %in% slotNames(ts_gse) &&
            nrow(ts_gse@compareClusterResult) > 0) {
          save_plot_dir(enrichplot::emapplot(ts_gse, showCategory = 10), dstP,
                        paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
        } else {
          writeLines("emapplot skipped: no enriched term after similarity calc.",
                     file.path(dstP, paste0("gseGO_", ont, "_emap_", file_tag, "_skipped.txt")))
        }
      }

      p_cnet <- tryCatch(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                         error = function(e) NULL)
      if (!is.null(p_cnet)) save_plot_dir(p_cnet, dstP, paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))

      # GO themes (restored)
      g_theme <- build_theme_graph(gse, top_n = 30, sim_cut = 0.35)
      save_plot_dir(plot_theme_graph(g_theme, paste0("GO themes (", ont, ") — ", context)),
                    paths$modules$themes, paste0("GO_", ont, "_themes_", file_tag, ".svg"))

      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) { log_line(paths, file_tag, paste("simplify error:", e$message)); NULL })
      if (!is.null(gse2) && has_terms(gse2, 1)) {
        log_line(paths, file_tag, sprintf("GO %s simplify: %d -> %d terms",
                                          ont, nrow(gse@result), nrow(gse2@result)))
        save_table(stable_rank(gse2@result), dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
      }
    } else {
      save_table(tibble::tibble(note = "gseGO empty or NULL"), dstT, paste0("gseGO_", ont, "_", file_tag, "_empty.csv"))
      writeLines("All plots skipped: no enriched term.", file.path(dstP, paste0("gseGO_", ont, "_", file_tag, "_skipped.txt")))
    }

    if (length(top_genes) >= 10) {
      ora <- with_timing({
        tryCatch(clusterProfiler::enrichGO(
          gene = top_genes, ont = ont, keyType = "UNIPROT",
          minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("ORA error:", e$message)); NULL })
      }, paths, file_tag, sprintf("ORA(%s)", ont))
      if (!is.null(ora) && has_terms(ora, 1)) {
        save_table(stable_rank(ora@result), dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        if (nrow(ora@result) >= 2) {
          ts_ora <- tryCatch(enrichplot::pairwise_termsim(ora), error = function(e) NULL)
          if (!is.null(ts_ora) && "compareClusterResult" %in% slotNames(ts_ora) &&
              nrow(ts_ora@compareClusterResult) > 0) {
            save_plot_dir(enrichplot::emapplot(ts_ora, showCategory = 10), dstP,
                          paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
          } else {
            writeLines("emapplot skipped: no enriched term after similarity calc.",
                       file.path(dstP, paste0("ORA_", ont, "_emap_", file_tag, "_skipped.txt")))
          }
        }
        p_cnet2 <- tryCatch(enrichplot::cnetplot(ora, showCategory = 5), error = function(e) NULL)
        if (!is.null(p_cnet2)) save_plot_dir(p_cnet2, dstP, paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
      } else {
        save_table(tibble::tibble(note = "ORA returned NULL or empty"), dstT, paste0("ORA_", ont, "_", file_tag, "_empty.csv"))
      }
    } else {
      save_table(tibble::tibble(note = "Global ORA skipped: <10 DEGs"),
                 dstT, paste0("ORA_", ont, "_", file_tag, "_skipped.csv"))
      log_line(paths, file_tag, paste("Global ORA skipped for", ont, "due to insufficient DEGs."))
    }
  }

  ## ----- GLOBAL: KEGG with de-dup -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  kmap_tbl <- tibble::tibble(
    id = universe_genes,
    ENTREZID = unname(entrez_map[universe_genes]),
    in_kegg = !is.na(entrez_map[universe_genes])
  )
  save_table(kmap_tbl, paths$glob$KEGG$T, paste0("kegg_mapping_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("KEGG mapIds: %d/%d mapped (%.1f%%)",
                                    sum(kmap_tbl$in_kegg), nrow(kmap_tbl), 100*mean(kmap_tbl$in_kegg)))

  kegg_gene_list <- gene_list[!is.na(entrez_map)]
  names(kegg_gene_list) <- unname(entrez_map[!is.na(entrez_map)])
  keep_idx <- !duplicated(names(kegg_gene_list))
  removed_dups <- sum(!keep_idx)
  kegg_gene_list <- kegg_gene_list[keep_idx]
  kegg_gene_list <- kegg_gene_list[order(kegg_gene_list, decreasing = TRUE, method = "radix")]

  # Extra diagnostics (restored)
  nmapped <- sum(!is.na(entrez_map))
  nunmapped <- sum(is.na(entrez_map))
  dupe_names <- names(kegg_gene_list)[duplicated(names(kegg_gene_list))]
  log_line(paths, file_tag, sprintf("KEGG mapIds diagnostics: mapped=%d, unmapped=%d, dup_ENTREZ_examples=%s",
                                    nmapped, nunmapped, paste(head(unique(dupe_names), 5), collapse=",")))
  log_line(paths, file_tag, sprintf("gseKEGG input: n=%d; removed_dups=%d; anyDuplicated=%d; pct ties=%.2f%%",
                                    length(kegg_gene_list), removed_dups, anyDuplicated(names(kegg_gene_list)),
                                    100*mean(duplicated(kegg_gene_list))))

  kk2 <- NULL
  if (length(kegg_gene_list) >= 10) {
    kk2 <- with_timing({
      muffle_fgsea_ties( tryCatch(clusterProfiler::gseKEGG(
        geneList = kegg_gene_list, organism = "mmu",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
        keyType = "ncbi-geneid"
      ), error = function(e) { log_line(paths, file_tag, paste("gseKEGG error:", e$message)); NULL }) )
    }, paths, file_tag, "gseKEGG")
  } else {
    log_line(paths, file_tag, "gseKEGG skipped: <10 unique ENTREZ in ranked list after de-dup.")
  }

  if (!is.null(kk2) && has_terms(kk2, 1)) {
    save_table(stable_rank(kk2@result), paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = 10, label_format = 30),
                  paths$glob$KEGG$P, paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    if (nrow(kk2@result) >= 2) {
      ts_kk <- tryCatch(enrichplot::pairwise_termsim(kk2), error = function(e) NULL)
      if (!is.null(ts_kk) && "compareClusterResult" %in% slotNames(ts_kk) &&
          nrow(ts_kk@compareClusterResult) > 0) {
        save_plot_dir(enrichplot::emapplot(ts_kk, showCategory = 10),
                      paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
      } else {
        writeLines("emapplot skipped: no enriched term after similarity calc.",
                   file.path(paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, "_skipped.txt")))
      }
    }
    save_plot_dir(plot_gsea_landscape(kk2, paste0("KEGG GSEA landscape — ", context)),
                  paths$glob$KEGG$P, paste0("gseKEGG_landscape_", file_tag, ".svg"))
    p_cnet3 <- tryCatch(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                        error = function(e) NULL)
    if (!is.null(p_cnet3)) save_plot_dir(p_cnet3, paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))

    # KEGG ridgeplot (optional restored)
    kkdf <- as.data.frame(kk2@result)
    if (nrow(kkdf) && "p.adjust" %in% names(kkdf)) {
      kkdf$mlogFDR <- -log10(pmax(1e-300, kkdf$p.adjust))
      rp <- ggplot2::ggplot(kkdf, ggplot2::aes(x = mlogFDR, y = "KEGG", fill = "KEGG")) +
        ggridges::geom_density_ridges(alpha = 0.6, color = NA, scale = 1) +
        ggplot2::theme_minimal() + ggplot2::labs(title = paste0("KEGG enrichment density — ", context),
                                                 x = "-log10(FDR)", y = NULL) +
        ggplot2::theme(legend.position = "none")
      save_plot_dir(rp, paths$glob$KEGG$P, paste0("gseKEGG_ridge_", file_tag, ".svg"), 20, 12)
    }

    # KEGG themes (restored)
    k_theme <- tryCatch(build_theme_graph(kk2, 30, 0.35), error = function(e) NULL)
    save_plot_dir(plot_theme_graph(k_theme, paste0("KEGG themes — ", context)),
                  paths$modules$themes, paste0("KEGG_themes_", file_tag, ".svg"))
  } else {
    save_table(tibble::tibble(note = "No KEGG enrichment results (gseKEGG null or empty)."),
               paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, "_empty.csv"))
    log_line(paths, file_tag, "gseKEGG produced no results.")
  }

  ## ----- GLOBAL: Reactome (optional, restored) -----
  have_reactome <- requireNamespace("ReactomePA", quietly = TRUE)
  if (have_reactome) {
    entrez_map_re <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                           column = "ENTREZID", multiVals = "first")
    re_universe <- unname(entrez_map_re[!is.na(entrez_map_re)])

    # ORA (DEG intersect)
    re_deg <- unname(entrez_map_re[match(top_genes, names(entrez_map_re))])
    re_deg <- re_deg[!is.na(re_deg)]
    re_ora <- tryCatch(ReactomePA::enrichPathway(
      gene = re_deg, organism = "mouse", universe = re_universe,
      pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 800
    ), error = function(e) NULL)
    if (!is.null(re_ora) && has_terms(re_ora, 1)) {
      save_table(stable_rank(re_ora@result), paths$glob$Reactome$T, paste0("ReactomeORA_", file_tag, ".csv"))
      save_plot_dir(plot_dot(re_ora, paste0("Reactome ORA — ", context)),
                    paths$glob$Reactome$P, paste0("ReactomeORA_dotplot_", file_tag, ".svg"))
    } else {
      save_table(tibble::tibble(note = "Reactome ORA null/empty"),
                 paths$glob$Reactome$T, paste0("ReactomeORA_", file_tag, "_empty.csv"))
    }

    # GSEA
    re_gene_list <- gene_list[!is.na(entrez_map_re)]
    names(re_gene_list) <- unname(entrez_map_re[!is.na(entrez_map_re)])
    keep_idx_re <- !duplicated(names(re_gene_list))
    re_gene_list <- re_gene_list[keep_idx_re]
    re_gene_list <- re_gene_list[order(re_gene_list, decreasing = TRUE, method = "radix")]

    re_gse <- tryCatch(ReactomePA::gsePathway(
      geneList = re_gene_list, organism = "mouse",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
      pAdjustMethod = "BH", eps = 0, nPermSimple = 10000
    ), error = function(e) NULL)

    if (!is.null(re_gse) && has_terms(re_gse, 1)) {
      save_table(stable_rank(re_gse@result), paths$glob$Reactome$T, paste0("ReactomeGSEA_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(re_gse, paste0("Reactome GSEA — ", context)),
                    paths$glob$Reactome$P, paste0("ReactomeGSEA_split_", file_tag, ".svg"))
      save_plot_dir(plot_gsea_landscape(re_gse, paste0("Reactome GSEA landscape — ", context)),
                    paths$glob$Reactome$P, paste0("ReactomeGSEA_landscape_", file_tag, ".svg"))

      # Reactome themes (restored)
      r_theme <- build_theme_graph(re_gse, 30, 0.35)
      save_plot_dir(plot_theme_graph(r_theme, paste0("Reactome themes — ", context)),
                    paths$modules$themes, paste0("Reactome_themes_", file_tag, ".svg"))
    } else {
      save_table(tibble::tibble(note = "Reactome GSEA null/empty"),
                 paths$glob$Reactome$T, paste0("ReactomeGSEA_", file_tag, "_empty.csv"))
    }
  } else {
    log_line(paths, file_tag, "ReactomePA not available; Reactome enrichment skipped.")
  }

  ## ----- PATHVIEW -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && has_terms(kk2, 1)) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  save_table(tibble::tibble(pv_ids = pv_ids), paths$pathview, paste0("pathview_ids_", file_tag, ".csv"))

  if (length(kegg_gene_list) > 0 && all(is.finite(as.numeric(kegg_gene_list)))) {
    writeLines(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
               con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
    log_line(paths, file_tag, sprintf("Pathview attempted for %d pathways.", length(pv_ids)))
  } else {
    writeLines(c(
      "No valid ENTREZ-mapped genes for Pathview or gene list non-numeric.",
      paste0("length(kegg_gene_list)=", length(kegg_gene_list))
    ), con = file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Pathview skipped: invalid kegg_gene_list.")
  }

  ## ----- MODULE-AWARE: overlaps, UpSet, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")

  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), , drop = FALSE]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids

  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))
  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      library(ComplexUpset)
      set_cols <- setdiff(colnames(mat), "id")
      mat[set_cols] <- lapply(mat[set_cols], function(x) as.logical(x))
      p_up <- ComplexUpset::upset(
        data = mat,
        intersect = set_cols,
        base_annotations = list('Intersection size' = ComplexUpset::intersection_size()),
        set_sizes = FALSE,
        min_size = 5
      )
      save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      library(UpSetR)
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
    } else {
      save_table(tibble::tibble(note="No UpSet library available (ComplexUpset/UpSetR)."),
                 paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped_no_pkg.csv"))
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
               paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }

  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, , drop = FALSE]
  if (nrow(ov2) > 0 && is_nonempty_num(ov2$OverlapDEG)) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2, 20)$OverlapDEG)*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No non-zero overlaps to plot.", con = file.path(paths$modules$overlap$P, paste0("overlap_status_", file_tag, ".txt")))
  }

  fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
    deg_set <- intersect(unique(deg_genes), universe_genes)
    bg_set  <- setdiff(universe_genes, deg_set)
    out <- lapply(names(modules_list), function(m) {
      mod_genes <- intersect(modules_list[[m]], universe_genes)
      a <- length(intersect(mod_genes, deg_set))
      b <- length(setdiff(mod_genes, deg_set))
      c <- length(setdiff(deg_set, mod_genes))
      d <- length(setdiff(bg_set, mod_genes))
      mm <- matrix(c(a,b,c,d), nrow=2)
      ft <- fisher.test(mm, alternative = "greater")
      data.frame(
        Module = m, InModule_DE = a, InModule_NotDE = b,
        NotInModule_DE = c, NotInModule_NotDE = d,
        OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    res <- dplyr::bind_rows(out)
    res$Padj <- p.adjust(res$Pvalue, method = "BH")
    dplyr::arrange(res, Padj, Pvalue)
  }
  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))

  if (nrow(fisher_tbl) > 0 && is_nonempty_num(fisher_tbl$Padj)) {
    ord_idx <- order(fisher_tbl$Padj, fisher_tbl$Pvalue)
    ft_plot <- ggplot2::ggplot(head(fisher_tbl[ord_idx, , drop = FALSE], 20),
                               ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                            y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#6a51a3") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Fisher enrichment (DEGs in modules) — ", context),
                    x = "Module", y = "-log10(Padj)")
    save_plot_dir(ft_plot, paths$modules$fisher$P, paste0("fisher_bar_", file_tag, ".svg"), 24, 20)
  } else {
    writeLines("No Fisher rows to plot.", con = file.path(paths$modules$fisher$P, paste0("fisher_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Fisher plot skipped: empty table.")
  }

  ## ----- Per-module analyses -----
  sel_modules <- if (nrow(fisher_tbl) > 0) fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)] else character(0)
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) {
      log_line(paths, file_tag, paste("Module", m, "skipped: <10 genes in universe."))
      next
    }

    for (ont in onts) {
      m_root     <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea   <- mk(m_root, "Tables", "gsea")
      m_P_gsea   <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull<- mk(m_root, "Tables", "ora_full")
      m_P_orafull<- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]
      m_gene_list <- m_gene_list[order(m_gene_list, decreasing = TRUE, method = "radix")]
      if (length(m_gene_list) >= 10) {
        m_gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
          geneList = m_gene_list, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "gseGO error:", e$message)); NULL }) )
        if (!is.null(m_gse) && has_terms(m_gse, 1)) {
          resm <- stable_rank(m_gse@result); resm$note <- "fgsea ties handled by deterministic id sort"
          save_table(resm, m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          p_dot <- tryCatch(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust",
                                                showCategory = 10, label_format = 30), error = function(e) NULL)
          if (!is.null(p_dot)) save_plot_dir(p_dot, m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "Module GSEA skipped: <10 genes."), m_T_gsea,
                   paste0("gseGO_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA DEG error:", e$message)); NULL })
        if (!is.null(m_ora) && has_terms(m_ora, 1)) {
          save_table(stable_rank(m_ora@result), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(
          note = "insufficient DEGs for ORA",
          n_module_genes_in_universe = length(m_genes),
          n_deg_in_module = length(intersect(top_genes, m_genes)),
          deg_threshold = 1
        ), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA full error:", e$message)); NULL })
        if (!is.null(m_ora_all) && has_terms(m_ora_all, 1)) {
          save_table(stable_rank(m_ora_all@result), m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "ORA full module skipped: <5 genes."),
                   m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }
    } # ont
  } # module

  ## ----- EXECUTIVE SUMMARY (robust) -----
  summary_items <- list()
  summary_items$universe <- data.frame(
    metric = c("universe_size","deg_abs_gt1","pct_tied"),
    value  = c(length(universe_genes), length(top_genes), round(pct_tied,2))
  )
  # Top GO/KEGG/Reactome tables if they exist
  get_first_csv <- function(path, pattern) {
    fps <- list.files(path, pattern = pattern, full.names = TRUE)
    if (length(fps)) tryCatch(readr::read_csv(fps[1], show_col_types = FALSE), error = function(e) NULL) else NULL
  }
  summary_items$top_go_gsea <- get_first_csv(paths$glob$BP$T, paste0("^gseGO_BP_", file_tag, "\\.csv$"))
  summary_items$top_kegg_gsea <- get_first_csv(paths$glob$KEGG$T, paste0("^gseKEGG_", file_tag, "\\.csv$"))
  summary_items$top_reactome_gsea <- get_first_csv(paths$glob$Reactome$T, paste0("^ReactomeGSEA_", file_tag, "\\.csv$"))
  summary_items$top_modules_fisher <- {
    ft <- tryCatch(readr::read_csv(file.path(paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv")), show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(ft)) ft[order(ft$Padj, ft$Pvalue), , drop = FALSE] else NULL
  }

  safe_univ_val <- function(x, metric, default = NA_character_) {
    if (is.null(x) || !nrow(x) || !"metric" %in% names(x) || !"value" %in% names(x)) return(default)
    v <- x$value[x$metric == metric]
    if (length(v) == 0) return(default)
    as.character(v[[1]])
  }
  fmt_top <- function(df, fmt) {
    if (is.null(df) || !nrow(df)) return(character(0))
    need <- intersect(c("ID","Description","p.adjust"), names(df))
    if (length(need) < 3) return(character(0))
    df <- df[order(df$p.adjust), , drop = FALSE]
    apply(head(df, 5), 1, function(r) sprintf(fmt, r[["ID"]], r[["Description"]], as.numeric(r[["p.adjust"]])))
  }

  md_lines <- list(
    "# Executive Summary",
    paste0("- Context: ", context),
    paste0("- Universe size: ", safe_univ_val(summary_items$universe, "universe_size", "NA")),
    paste0("- DEGs |log2fc|>1: ", safe_univ_val(summary_items$universe, "deg_abs_gt1", "NA")),
    paste0("- pct tied stats: ", safe_univ_val(summary_items$universe, "pct_tied", "NA"), "%"),
    "",
    "## Top GO GSEA (BP)"
  )
  go_lines <- fmt_top(summary_items$top_go_gsea, "%s — %s (FDR=%.2g)")
  md_lines <- c(md_lines, if (length(go_lines)) paste0("* ", go_lines) else "No terms", "",
                "## Top KEGG GSEA")
  kegg_lines <- fmt_top(summary_items$top_kegg_gsea, "%s — %s (FDR=%.2g)")
  md_lines <- c(md_lines, if (length(kegg_lines)) paste0("* ", kegg_lines) else "No terms", "",
                "## Top Reactome GSEA")
  reactome_lines <- fmt_top(summary_items$top_reactome_gsea, "%s — %s (FDR=%.2g)")
  md_lines <- c(md_lines, if (length(reactome_lines)) paste0("* ", reactome_lines) else "Skipped/None", "",
                "## Top Modules by Fisher (Padj)")
  if (!is.null(summary_items$top_modules_fisher) && nrow(summary_items$top_modules_fisher)) {
    mod_tbl <- head(summary_items$top_modules_fisher, 5)
    if (all(c("Module","OddsRatio","Padj") %in% names(mod_tbl))) {
      mod_lines <- apply(mod_tbl, 1, function(r) sprintf("%s — OR=%.2f (Padj=%.2g)",
                                                         r[["Module"]], as.numeric(r[["OddsRatio"]]), as.numeric(r[["Padj"]])))
      md_lines <- c(md_lines, paste0("* ", mod_lines))
    } else {
      md_lines <- c(md_lines, "No modules (columns missing).")
    }
  } else {
    md_lines <- c(md_lines, "No modules")
  }
  md <- unlist(md_lines, use.names = FALSE)
  writeLines(md, con = file.path(paths$summary, paste0("SUMMARY_", file_tag, ".md")))

  writeLines(paste("End:", Sys.time()), con = file.path(paths$meta, "run_end.txt"))
  log_line(paths, file_tag, "Completed all analyses.")
}































############################################################
# clusterProfiler + WGCNA modules (mouse) — extended figures
# Fixed: reliable directory creation, robust SUMMARY writer,
# deterministic summary file selection, non-empty themes logs,
# ggsave(create.dir=TRUE) for ggplot2 >= 3.5.0.
############################################################

## ---------- Package setup ----------
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
  required_packages <- c(
    "clusterProfiler","enrichplot","DOSE","fgsea","ggplot2","ggridges","svglite",
    "org.Mm.eg.db","AnnotationDbi","tibble","dplyr","tidyr","stringr","readr",
    "purrr","forcats","tools","pathview","KEGGREST","vroom","UniProt.ws","BiocParallel",
    "BiocGenerics","S4Vectors","IRanges","GenomeInfoDb","GenomicRanges","SummarizedExperiment",
    "Biobase","GO.db","UpSetR","ComplexUpset","igraph","ggraph"
  )
  bioc_packages <- c("clusterProfiler","enrichplot","DOSE","org.Mm.eg.db","fgsea","KEGGREST","ReactomePA")
  installBioC(bioc_packages)
  invisible(lapply(required_packages, installAndLoadCRAN))
}
setupPackages()

## ---------- Directory + I/O helpers ----------
mk <- function(...) {
  p <- file.path(...)
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}
ensure_parent <- function(fp) {
  d <- dirname(fp)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  fp
}
subdir_paths <- function(results_dir, file_tag) {
  list(
    meta = mk(results_dir, "00_meta"),
    glob = list(
      BP     = list(T = mk(results_dir, "10_global","GO_BP","Tables"),
                    P = mk(results_dir, "10_global","GO_BP","Plots")),
      MF     = list(T = mk(results_dir, "10_global","GO_MF","Tables"),
                    P = mk(results_dir, "10_global","GO_MF","Plots")),
      CC     = list(T = mk(results_dir, "10_global","GO_CC","Tables"),
                    P = mk(results_dir, "10_global","GO_CC","Plots")),
      KEGG   = list(T = mk(results_dir, "10_global","KEGG","Tables"),
                    P = mk(results_dir, "10_global","KEGG","Plots")),
      Reactome = list(T = mk(results_dir, "10_global","Reactome","Tables"),
                      P = mk(results_dir, "10_global","Reactome","Plots"))
    ),
    modules = list(
      root    = mk(results_dir, "20_modules"),
      overlap = list(T = mk(results_dir, "20_modules","overlap_plots","Tables"),
                     P = mk(results_dir, "20_modules","overlap_plots","Plots")),
      fisher  = list(T = mk(results_dir, "20_modules","fisher_deg","Tables"),
                     P = mk(results_dir, "20_modules","fisher_deg","Plots")),
      themes  = mk(results_dir, "20_modules","themes_networks")
    ),
    pathview = mk(results_dir, "90_pathview"),
    summary  = mk(results_dir, "99_summary")
  )
}
save_table <- function(df, dirpath, fname) {
  if (is.null(df)) return(invisible(NULL))
  if (!dir.exists(dirpath)) dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(df, ensure_parent(file.path(dirpath, fname)))
}
save_plot_dir <- function(plot, dirpath, filename, width = 20, height = 16) {
  if (is.null(plot)) return(invisible(NULL))
  if (!dir.exists(dirpath)) dir.create(dirpath, recursive = TRUE, showWarnings = FALSE)
  # Final normalization for ggplot2 >= 3.5.0
  plot <- plot + ggplot2::theme(
    axis.ticks.length         = grid::unit(2, "pt"),
    axis.ticks.length.x       = grid::unit(2, "pt"),
    axis.ticks.length.y       = grid::unit(2, "pt"),
    axis.minor.ticks.length   = grid::unit(1, "pt"),
    axis.minor.ticks.length.x = grid::unit(1, "pt"),
    axis.minor.ticks.length.y = grid::unit(1, "pt")
  )
  ggplot2::ggsave(filename = file.path(dirpath, filename), plot = plot,
                  units = "cm", dpi = 300, width = width, height = height,
                  create.dir = TRUE)
}
write_lines_safe <- function(lines, fp) {
  ensure_parent(fp)
  writeLines(lines, con = fp)
}
log_line <- function(paths, tag, text) {
  lf <- ensure_parent(file.path(paths$meta, paste0("runlog_", tag, ".txt")))
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text), file = lf, append = TRUE)
}
stable_rank <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  if (!("pvalue" %in% names(df)) && ("p.adjust" %in% names(df))) df$pvalue <- df$p.adjust
  if (!("p.adjust" %in% names(df)) && ("pvalue" %in% names(df))) df$p.adjust <- p.adjust(df$pvalue, method = "BH")
  df$rank_p <- rank(df$pvalue, ties.method = "min")
  df$rank_padj <- rank(df$p.adjust, ties.method = "min")
  df[order(df$p.adjust, df$pvalue, df$rank_padj, df$rank_p, decreasing = FALSE), ]
}
has_terms <- function(x, min_n = 1) {
  tryCatch({
    df <- as.data.frame(x)
    !is.null(df) && nrow(df) >= min_n
  }, error = function(e) FALSE)
}
muffle_fgsea_ties <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("There are ties in the preranked stats", conditionMessage(m))) {
      invokeRestart("muffleMessage")
    }
  })
}
# Defensive helpers
is_nonempty_num <- function(x) is.numeric(x) && length(x) > 0
is_nonempty_char <- function(x) is.character(x) && length(x) > 0

with_timing <- function(expr, paths, tag, label) {
  t0 <- Sys.time()
  res <- NULL; err <- NULL
  tryCatch({ res <- eval.parent(substitute(expr)) },
           error = function(e) { err <<- e })
  t1 <- Sys.time()
  elapsed <- as.numeric(difftime(t1, t0, units = "secs"))
  log_line(paths, tag, sprintf("%s elapsed=%.2fs %s",
                               label, elapsed,
                               if (!is.null(err)) paste("| error:", err$message) else "" ))
  if (!is.null(err)) stop(err)
  res
}

## ---------- Plot helpers ----------
build_split_dotplot <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj)) return(invisible(NULL))
  df <- tryCatch(as.data.frame(gsea_obj), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(invisible(NULL))
  if (!inherits(gsea_obj, "gseaResult")) return(invisible(NULL))
  if (!".sign" %in% colnames(df)) {
    if (!"NES" %in% colnames(df)) return(invisible(NULL))
    df$.sign <- ifelse(df$NES > 0, "upregulated", "downregulated")
    gsea_obj@result$.sign <- df$.sign
  }
  enrichplot::dotplot(
    gsea_obj, x = "GeneRatio", color = "p.adjust",
    showCategory = min(10, nrow(df)), split = ".sign", font.size = 12, label_format = 30
    ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 10),
    strip.text  = ggplot2::element_text(size = 12),
    panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    axis.ticks.length = grid::unit(2, "pt")
  )
}
plot_dot <- function(dataset, title_txt) {
  if (is.null(dataset)) return(invisible(NULL))
  resdf <- tryCatch(as.data.frame(dataset), error = function(e) NULL)
  if (is.null(resdf) || nrow(resdf) == 0) return(invisible(NULL))
  enrichplot::dotplot(dataset, x = "GeneRatio", color = "p.adjust",
                      showCategory = min(10, nrow(resdf)), font.size = 12, label_format = 30) +
    ggplot2::labs(title = title_txt, x = "Gene Ratio", y = "Gene Set") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
    plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 10),
    panel.grid.major = ggplot2::element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    axis.ticks.length = grid::unit(2, "pt")
  )
}
plot_gsea_landscape <- function(gsea_obj, title_txt) {
  if (is.null(gsea_obj) || !has_terms(gsea_obj, 1)) return(NULL)
  df <- as.data.frame(gsea_obj@result)
  if (!nrow(df) || !"NES" %in% names(df)) return(NULL)
  df$mlogFDR <- -log10(pmax(1e-300, df$p.adjust))
  ggplot2::ggplot(df, ggplot2::aes(x = NES, y = mlogFDR)) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(title = title_txt, x = "NES", y = "-log10(FDR)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.ticks.length = grid::unit(2, "pt"))

}
prep_leading_edge_matrix <- function(gsea_obj, gene_list, top_terms = 5) {
  if (is.null(gsea_obj) || !has_terms(gsea_obj, 1)) return(NULL)
  res <- as.data.frame(gsea_obj@result)
  if (!nrow(res)) return(NULL)
  if (!("p.adjust" %in% names(res)) || !("NES" %in% names(res))) return(NULL)
  res <- res[order(res$p.adjust, -abs(res$NES)), , drop = FALSE]
  res <- head(res, top_terms)
  if (!nrow(res)) return(NULL)
  if (!"core_enrichment" %in% colnames(res)) return(NULL)
  le_genes <- unique(unlist(strsplit(paste(res$core_enrichment, collapse = "/"), "/")))
  le_genes <- intersect(le_genes, names(gene_list))
  if (!length(le_genes)) return(NULL)
  mat <- matrix(gene_list[le_genes], nrow = length(le_genes), dimnames = list(le_genes, "log2fc"))
  mat
}
plot_leading_edge_heatmap <- function(mat, title_txt) {
  if (is.null(mat)) return(NULL)
  df <- data.frame(Gene = rownames(mat), log2fc = mat[,1], row.names = NULL, stringsAsFactors = FALSE)
  if (!nrow(df)) return(NULL)
  df$Gene <- factor(df$Gene, levels = df$Gene[order(df$log2fc, decreasing = TRUE)])
  ggplot2::ggplot(df, ggplot2::aes(x = "Leading edge", y = Gene, fill = log2fc)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
    ggplot2::labs(title = title_txt, x = NULL, y = "Genes (leading edge)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(2, "pt")
    )

}

## ---------- Theme graph helpers (restored, more permissive) ----------
build_theme_graph <- function(enrich_obj, top_n = 50, sim_cut = 0.2) {
  if (is.null(enrich_obj) || !has_terms(enrich_obj, 2)) return(NULL)
  ts <- tryCatch(enrichplot::pairwise_termsim(enrich_obj), error = function(e) NULL)
  if (is.null(ts) || !"compareClusterResult" %in% slotNames(ts)) return(NULL)
  dfr <- tryCatch(as.data.frame(ts@compareClusterResult), error = function(e) NULL)
  if (is.null(dfr) || !nrow(dfr)) return(NULL)
  dfr <- dfr[order(dfr$p.adjust, dfr$pvalue), , drop = FALSE]
  ids <- head(dfr$ID, top_n)
  sim <- ts@termsim[rownames(ts@termsim) %in% ids, colnames(ts@termsim) %in% ids, drop = FALSE]
  if (is.null(sim) || !nrow(sim) || !ncol(sim)) return(NULL)
  ed <- as.data.frame(as.table(sim), stringsAsFactors = FALSE)
  ed <- ed[ed$Freq >= sim_cut & ed$Var1 != ed$Var2, , drop = FALSE]
  if (!nrow(ed)) return(NULL)
  g <- igraph::graph_from_data_frame(ed, directed = FALSE)
  igraph::V(g)$label <- dfr$Description[match(igraph::V(g)$name, dfr$ID)]
  igraph::V(g)$size  <- -log10(pmax(1e-300, dfr$p.adjust[match(igraph::V(g)$name, dfr$ID)]))
  g
}
plot_theme_graph <- function(g, title_txt = "Enrichment themes") {
  if (is.null(g)) return(NULL)
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(alpha = 0.25) +
    ggraph::geom_node_point(ggplot2::aes(size = size, color = size)) +
    ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE, size = 3) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme_void() +
    ggplot2::theme(axis.ticks.length = grid::unit(2, "pt")) +
    ggplot2::labs(title = title_txt, color = "-log10(FDR)", size = "-log10(FDR)")
}

## ---------- Load WGCNA modules ----------
modules_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output/wgcna_modules_mapped"
module_files <- list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE)
stopifnot("No WGCNA module CSVs found" = length(module_files) > 0)

read_module_csv <- function(fp) {
  df <- read.csv(fp, header = TRUE, check.names = FALSE)
  colnames(df)[1:4] <- c("Gene","Module","ModuleMembership","GeneSignificance")
  df
}
modules_df <- dplyr::bind_rows(lapply(module_files, read_module_csv)) |>
  dplyr::filter(!is.na(Module), Module != "", !is.na(Gene), Gene != "") |>
  dplyr::mutate(Gene = as.character(Gene))

modules_list <- split(modules_df$Gene, modules_df$Module)
modules_list <- lapply(modules_list, function(x) unique(x))
module_sizes <- tibble::tibble(
  Module = names(modules_list),
  N_accessions = vapply(modules_list, length, integer(1))
)

## ---------- Annotation DB ----------
library(org.Mm.eg.db)
orgdb <- org.Mm.eg.db

## ---------- Global parameters ----------
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))  # for unit()

# Ensure axis.ticks.length are grid::unit objects for ggplot2 >= 3.5
ggplot2::theme_set(ggplot2::theme_get() + ggplot2::theme(
  axis.ticks.length   = grid::unit(2, "pt"),
  axis.ticks.length.x = grid::unit(2, "pt"),
  axis.ticks.length.y = grid::unit(2, "pt")
))

## ---------- Inputs ----------
comparison <- "neuron-phenotypeWithinUnit"
file_direction <- paste0("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/mapped/", comparison)
file_list <- list.files(file_direction, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_direction)

## ---------- Valid KEGG pathway IDs (mmu) ----------
kegg_df <- KEGGREST::keggList("pathway", "mmu") |> tibble::enframe(name = NULL, value = "raw")
kegg_df <- tidyr::separate(kegg_df, raw, into = c("kegg_id","name"), sep = "\\s+", extra = "merge", fill = "right")
valid_kegg_ids <- unique(sub("^path:", "", kegg_df$kegg_id))

## ---------- Main loop ----------
for (data_path in file_list) {
  file_name    <- basename(data_path)
  name_no_ext  <- tools::file_path_sans_ext(file_name)

  # Context label
  us <- gregexpr("_", name_no_ext, fixed = TRUE)[[1]]
  if (identical(us[1], -1L) || length(us) == 1L) {
    context <- name_no_ext
  } else {
    pos2 <- us[2]
    left  <- substr(name_no_ext, 1, pos2 - 1)
    right <- substr(name_no_ext, pos2 + 1, nchar(name_no_ext))
    context <- paste0(left, " vs ", right)
  }

  ## ----- Directories -----
  working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
  setwd(working_dir)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  results_dir <- mk(working_dir, "Results", comparison, file_tag)
  paths <- subdir_paths(results_dir, file_tag)

  # Provenance snapshot
  capture.output(sessionInfo(), file = ensure_parent(file.path(paths$meta, "sessionInfo.txt")))
  write_lines_safe(capture.output(Sys.getenv()), file.path(paths$meta, "env.txt"))
  write_lines_safe(paste("Start:", Sys.time()), file.path(paths$meta, "run_start.txt"))
  save_table(module_sizes, paths$modules$root, "module_sizes.csv")

  # Parameters README
  write_lines_safe(c(
    "# Analysis README",
    paste0("Comparison: ", comparison),
    paste0("Context: ", context),
    "Identifiers: UniProt accessions",
    "GO keyType: UNIPROT; KEGG keyType: ncbi-geneid (mapIds from org.Mm.eg.db, multiVals='first')",
    "DEG threshold for ORA: |log2fc| > 1",
    "Ontologies default: BP",
    "ComplexUpset: set sizes panel disabled for portability",
    "GSEA ties: sorted by id for deterministic order; values are NOT jittered",
    "Folders:",
    "  10_global/GO_*: GSEA and ORA across full universe",
    "  10_global/KEGG: gseKEGG results and mapping diagnostics",
    "  10_global/Reactome: optional Reactome enrichment",
    "  20_modules/: overlap, Fisher, per-module GO (GSEA/ORA)",
    "  20_modules/themes_networks: similarity-based term theme graphs",
    "  90_pathview/: KEGG pathway renderings and status logs",
    "  99_summary/: executive summary"
  ), file.path(results_dir, "README.txt"))
  log_line(paths, file_tag, "Initialized results folders and README.")

  # Whole-iteration error capture to guarantee summary write
  iteration_error <- NULL
  on.exit({
    # EXECUTIVE SUMMARY (robust, always attempt)
    summary_items <- list()
    summary_items$universe <- data.frame(
      metric = c("universe_size","deg_abs_gt1","pct_tied"),
      value  = c(length(universe_genes), length(top_genes), round(pct_tied,2))
    )

    # Deterministic CSV pickers
    pick_first <- function(dir, patterns_by_priority) {
      if (!dir.exists(dir)) return(NULL)
      fps <- unlist(lapply(patterns_by_priority, function(p) {
        list.files(dir, pattern = p, full.names = TRUE)
      }), use.names = FALSE)
      fps <- unique(fps)
      if (length(fps) == 0) return(NULL)
      tryCatch(readr::read_csv(fps[1], show_col_types = FALSE), error = function(e) NULL)
    }
    summary_items$top_go_gsea <- pick_first(paths$glob$BP$T,
      c(paste0("^gseGO_BP_simplified_", file_tag, "\\.csv$"),
        paste0("^gseGO_BP_", file_tag, "\\.csv$")) )
    summary_items$top_kegg_gsea <- pick_first(paths$glob$KEGG$T,
      c(paste0("^gseKEGG_", file_tag, "\\.csv$")) )
    summary_items$top_reactome_gsea <- pick_first(paths$glob$Reactome$T,
      c(paste0("^ReactomeGSEA_", file_tag, "\\.csv$")) )
    summary_items$top_modules_fisher <- {
      fp <- file.path(paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))
      if (file.exists(fp)) {
        ft <- tryCatch(readr::read_csv(fp, show_col_types = FALSE), error = function(e) NULL)
        if (!is.null(ft)) ft[order(ft$Padj, ft$Pvalue), , drop = FALSE] else NULL
      } else NULL
    }

    safe_univ_val <- function(x, metric, default = NA_character_) {
      if (is.null(x) || !nrow(x) || !"metric" %in% names(x) || !"value" %in% names(x)) return(default)
      v <- x$value[x$metric == metric]
      if (length(v) == 0) return(default)
      as.character(v[[1]])
    }
    fmt_top <- function(df, fmt) {
      if (is.null(df) || !nrow(df)) return(character(0))
      need <- intersect(c("ID","Description","p.adjust"), names(df))
      if (length(need) < 3) return(character(0))
      df <- df[order(df$p.adjust), , drop = FALSE]
      apply(head(df, 5), 1, function(r) sprintf(fmt, r[["ID"]], r[["Description"]], as.numeric(r[["p.adjust"]])))
    }

    md_lines <- list(
      "# Executive Summary",
      paste0("- Context: ", context),
      paste0("- Universe size: ", safe_univ_val(summary_items$universe, "universe_size", "NA")),
      paste0("- DEGs |log2fc|>1: ", safe_univ_val(summary_items$universe, "deg_abs_gt1", "NA")),
      paste0("- pct tied stats: ", safe_univ_val(summary_items$universe, "pct_tied", "NA"), "%"),
      ""
    )
    go_lines <- fmt_top(summary_items$top_go_gsea, "%s — %s (FDR=%.2g)")
    md_lines <- c(md_lines, "## Top GO GSEA (BP)",
                  if (length(go_lines)) paste0("* ", go_lines) else "No terms", "")
    kegg_lines <- fmt_top(summary_items$top_kegg_gsea, "%s — %s (FDR=%.2g)")
    md_lines <- c(md_lines, "## Top KEGG GSEA",
                  if (length(kegg_lines)) paste0("* ", kegg_lines) else "No terms", "")
    reactome_lines <- fmt_top(summary_items$top_reactome_gsea, "%s — %s (FDR=%.2g)")
    md_lines <- c(md_lines, "## Top Reactome GSEA",
                  if (length(reactome_lines)) paste0("* ", reactome_lines) else "Skipped/None", "")
    md_lines <- c(md_lines, "## Top Modules by Fisher (Padj)")
    if (!is.null(summary_items$top_modules_fisher) && nrow(summary_items$top_modules_fisher)) {
      mod_tbl <- head(summary_items$top_modules_fisher, 5)
      if (all(c("Module","OddsRatio","Padj") %in% names(mod_tbl))) {
        mod_lines <- apply(mod_tbl, 1, function(r) sprintf("%s — OR=%.2f (Padj=%.2g)",
                                                           r[["Module"]], as.numeric(r[["OddsRatio"]]), as.numeric(r[["Padj"]])))
        md_lines <- c(md_lines, paste0("* ", mod_lines))
      } else {
        md_lines <- c(md_lines, "No modules (columns missing).")
      }
    } else {
      md_lines <- c(md_lines, "No modules")
    }
    if (!is.null(iteration_error)) {
      md_lines <- c(md_lines, "", "## Notes", paste0("* Error encountered earlier: ", iteration_error))
    }
    md <- unlist(md_lines, use.names = FALSE)
    write_lines_safe(md, file.path(paths$summary, paste0("SUMMARY_", file_tag, ".md")))
    write_lines_safe(paste("End:", Sys.time()), file.path(paths$meta, "run_end.txt"))
    log_line(paths, file_tag, "Completed all analyses (summary written).")
  }, add = TRUE)

  ## ----- Load comparison -----
  df_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  colnames(df_raw)[1] <- "id"
  stopifnot("log2fc column missing" = "log2fc" %in% colnames(df_raw))
  df <- df_raw |>
    dplyr::filter(!is.na(id), id != "", !is.na(log2fc)) |>
    dplyr::mutate(id = as.character(id)) |>
    dplyr::group_by(id) |>
    dplyr::slice_max(order_by = abs(log2fc), n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  gene_list_tbl <- df |> dplyr::select(id, log2fc) |> dplyr::arrange(dplyr::desc(log2fc), id)
  gene_list <- gene_list_tbl$log2fc; names(gene_list) <- gene_list_tbl$id
  universe_genes <- names(gene_list)
  deg_tbl <- df |> dplyr::filter(abs(log2fc) > 1)
  top_genes <- deg_tbl$id

  pct_tied <- 100 * mean(duplicated(gene_list))
  save_table(tibble::tibble(
    n_rows = nrow(df),
    n_universe = length(universe_genes),
    n_deg_abs1 = length(top_genes),
    deg_threshold = 1,
    pct_tied_stats = pct_tied
  ), paths$meta, paste0("summary_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("Universe=%d; DEGs(|log2fc|>1)=%d; tied stats=%.2f%%",
                                    length(universe_genes), length(top_genes), pct_tied))

  organism <- "org.Mm.eg.db"
  onts <- c("BP", "MF", "CC")

  ## ----- GLOBAL: GO -----
  for (ont in onts) {
    dstT <- switch(ont, "BP" = paths$glob$BP$T, "MF" = paths$glob$MF$T, "CC" = paths$glob$CC$T)
    dstP <- switch(ont, "BP" = paths$glob$BP$P, "MF" = paths$glob$MF$P, "CC" = paths$glob$CC$P)

    gse <- with_timing({
      muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
        geneList = gene_list, ont = ont, keyType = "UNIPROT",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
        by = "fgsea", nPermSimple = 10000, eps = 0,
        verbose = TRUE, OrgDb = organism, pAdjustMethod = "BH"
      ), error = function(e) { log_line(paths, file_tag, paste("gseGO error:", e$message)); NULL }) )
    }, paths, file_tag, sprintf("gseGO(%s)", ont))

    if (!is.null(gse) && has_terms(gse, 1)) {
      res_stable <- stable_rank(gse@result)
      res_stable$note <- "fgsea ties handled by deterministic id sort (values unchanged)"
      save_table(res_stable, dstT, paste0("gseGO_", ont, "_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(gse, paste0("GSEA (", ont, ") — ", context)), dstP,
                    paste0("gseGO_", ont, "_split_", file_tag, ".svg"))
      save_plot_dir(plot_gsea_landscape(gse, paste0("GSEA landscape (", ont, ") — ", context)),
                    dstP, paste0("gseGO_", ont, "_landscape_", file_tag, ".svg"))

      le_mat <- prep_leading_edge_matrix(gse, gene_list, top_terms = 5)
      save_plot_dir(plot_leading_edge_heatmap(le_mat, paste0("Leading edge genes (", ont, ") — ", context)),
                    dstP, paste0("gseGO_", ont, "_leading_edge_", file_tag, ".svg"))

      if (nrow(gse@result) >= 2) {
        ts_gse <- tryCatch(enrichplot::pairwise_termsim(gse), error = function(e) NULL)
        if (!is.null(ts_gse) && "compareClusterResult" %in% slotNames(ts_gse) &&
            nrow(ts_gse@compareClusterResult) > 0) {
          save_plot_dir(enrichplot::emapplot(ts_gse, showCategory = 10), dstP,
                        paste0("gseGO_", ont, "_emap_", file_tag, ".svg"))
        } else {
          write_lines_safe("emapplot skipped: no enriched term after similarity calc.",
                           file.path(dstP, paste0("gseGO_", ont, "_emap_", file_tag, "_skipped.txt")))
        }
      }

      p_cnet <- tryCatch(enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list),
                         error = function(e) NULL)
      if (!is.null(p_cnet)) save_plot_dir(p_cnet, dstP, paste0("gseGO_", ont, "_cnet_", file_tag, ".svg"))

      # GO themes (restored, robust)
      g_theme <- tryCatch(build_theme_graph(gse, top_n = 50, sim_cut = 0.2), error = function(e) NULL)
      if (!is.null(g_theme)) {
        save_plot_dir(plot_theme_graph(g_theme, paste0("GO themes (", ont, ") — ", context)),
                      paths$modules$themes, paste0("GO_", ont, "_themes_", file_tag, ".svg"))
      } else {
        write_lines_safe("Theme graph skipped: insufficient terms or edges after similarity threshold.",
                         file.path(paths$modules$themes, paste0("GO_", ont, "_themes_", file_tag, "_skipped.txt")))
      }

      gse2 <- tryCatch(clusterProfiler::simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min),
                       error = function(e) { log_line(paths, file_tag, paste("simplify error:", e$message)); NULL })
      if (!is.null(gse2) && has_terms(gse2, 1)) {
        log_line(paths, file_tag, sprintf("GO %s simplify: %d -> %d terms",
                                          ont, nrow(gse@result), nrow(gse2@result)))
        save_table(stable_rank(gse2@result), dstT, paste0("gseGO_", ont, "_simplified_", file_tag, ".csv"))
        save_plot_dir(build_split_dotplot(gse2, paste0("GSEA simplified (", ont, ") — ", context)), dstP,
                      paste0("gseGO_", ont, "_simplified_split_", file_tag, ".svg"))
      }
    } else {
      save_table(tibble::tibble(note = "gseGO empty or NULL"), dstT, paste0("gseGO_", ont, "_", file_tag, "_empty.csv"))
      write_lines_safe("All plots skipped: no enriched term.",
                       file.path(dstP, paste0("gseGO_", ont, "_", file_tag, "_skipped.txt")))
    }

    # ORA
    if (length(top_genes) >= 10) {
      ora <- with_timing({
        tryCatch(clusterProfiler::enrichGO(
          gene = top_genes, ont = ont, keyType = "UNIPROT",
          minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("ORA error:", e$message)); NULL })
      }, paths, file_tag, sprintf("ORA(%s)", ont))
      if (!is.null(ora) && has_terms(ora, 1)) {
        save_table(stable_rank(ora@result), dstT, paste0("ORA_", ont, "_", file_tag, ".csv"))
        save_plot_dir(plot_dot(ora, paste0("ORA (", ont, ") — ", context)), dstP,
                      paste0("ORA_", ont, "_dotplot_", file_tag, ".svg"))
        if (nrow(ora@result) >= 2) {
          ts_ora <- tryCatch(enrichplot::pairwise_termsim(ora), error = function(e) NULL)
          if (!is.null(ts_ora) && "compareClusterResult" %in% slotNames(ts_ora) &&
              nrow(ts_ora@compareClusterResult) > 0) {
            save_plot_dir(enrichplot::emapplot(ts_ora, showCategory = 10), dstP,
                          paste0("ORA_", ont, "_emap_", file_tag, ".svg"))
          } else {
            write_lines_safe("emapplot skipped: no enriched term after similarity calc.",
                             file.path(dstP, paste0("ORA_", ont, "_emap_", file_tag, "_skipped.txt")))
          }
        }
        p_cnet2 <- tryCatch(enrichplot::cnetplot(ora, showCategory = 5), error = function(e) NULL)
        if (!is.null(p_cnet2)) save_plot_dir(p_cnet2, dstP, paste0("ORA_", ont, "_cnet_", file_tag, ".svg"))
      } else {
        save_table(tibble::tibble(note = "ORA returned NULL or empty"), dstT, paste0("ORA_", ont, "_", file_tag, "_empty.csv"))
      }
    } else {
      save_table(tibble::tibble(note = "Global ORA skipped: <10 DEGs"),
                 dstT, paste0("ORA_", ont, "_", file_tag, "_skipped.csv"))
      log_line(paths, file_tag, paste("Global ORA skipped for", ont, "due to insufficient DEGs."))
    }
  }

  ## ----- GLOBAL: KEGG with de-dup -----
  entrez_map <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                      column = "ENTREZID", multiVals = "first")
  kmap_tbl <- tibble::tibble(
    id = universe_genes,
    ENTREZID = unname(entrez_map[universe_genes]),
    in_kegg = !is.na(entrez_map[universe_genes])
  )
  save_table(kmap_tbl, paths$glob$KEGG$T, paste0("kegg_mapping_", file_tag, ".csv"))
  log_line(paths, file_tag, sprintf("KEGG mapIds: %d/%d mapped (%.1f%%)",
                                    sum(kmap_tbl$in_kegg), nrow(kmap_tbl), 100*mean(kmap_tbl$in_kegg)))

  kegg_gene_list <- gene_list[!is.na(entrez_map)]
  names(kegg_gene_list) <- unname(entrez_map[!is.na(entrez_map)])
  keep_idx <- !duplicated(names(kegg_gene_list))
  removed_dups <- sum(!keep_idx)
  kegg_gene_list <- kegg_gene_list[keep_idx]
  kegg_gene_list <- kegg_gene_list[order(kegg_gene_list, decreasing = TRUE, method = "radix")]

  # Extra diagnostics
  nmapped <- sum(!is.na(entrez_map))
  nunmapped <- sum(is.na(entrez_map))
  dupe_names <- names(kegg_gene_list)[duplicated(names(kegg_gene_list))]
  log_line(paths, file_tag, sprintf("KEGG mapIds diagnostics: mapped=%d, unmapped=%d, dup_ENTREZ_examples=%s",
                                    nmapped, nunmapped, paste(head(unique(dupe_names), 5), collapse=",")))
  log_line(paths, file_tag, sprintf("gseKEGG input: n=%d; removed_dups=%d; anyDuplicated=%d; pct ties=%.2f%%",
                                    length(kegg_gene_list), removed_dups, anyDuplicated(names(kegg_gene_list)),
                                    100*mean(duplicated(kegg_gene_list))))

  kk2 <- NULL
  if (length(kegg_gene_list) >= 10) {
    kk2 <- with_timing({
      muffle_fgsea_ties( tryCatch(clusterProfiler::gseKEGG(
        geneList = kegg_gene_list, organism = "mmu",
        minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, pAdjustMethod = "BH",
        keyType = "ncbi-geneid"
      ), error = function(e) { log_line(paths, file_tag, paste("gseKEGG error:", e$message)); NULL }) )
    }, paths, file_tag, "gseKEGG")
  } else {
    log_line(paths, file_tag, "gseKEGG skipped: <10 unique ENTREZ in ranked list after de-dup.")
  }

  if (!is.null(kk2) && has_terms(kk2, 1)) {
    save_table(stable_rank(kk2@result), paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, ".csv"))
    save_plot_dir(enrichplot::dotplot(kk2, x = "GeneRatio", color = "p.adjust",
                                      showCategory = min(10, nrow(kk2@result)), label_format = 30),
                  paths$glob$KEGG$P, paste0("gseKEGG_dotplot_", file_tag, ".svg"))
    if (nrow(kk2@result) >= 2) {
      ts_kk <- tryCatch(enrichplot::pairwise_termsim(kk2), error = function(e) NULL)
      if (!is.null(ts_kk) && "compareClusterResult" %in% slotNames(ts_kk) &&
          nrow(ts_kk@compareClusterResult) > 0) {
        save_plot_dir(enrichplot::emapplot(ts_kk, showCategory = 10),
                      paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, ".svg"))
      } else {
        write_lines_safe("emapplot skipped: no enriched term after similarity calc.",
                         file.path(paths$glob$KEGG$P, paste0("gseKEGG_emap_", file_tag, "_skipped.txt")))
      }
    }
    save_plot_dir(plot_gsea_landscape(kk2, paste0("KEGG GSEA landscape — ", context)),
                  paths$glob$KEGG$P, paste0("gseKEGG_landscape_", file_tag, ".svg"))
    p_cnet3 <- tryCatch(enrichplot::cnetplot(kk2, categorySize = "pvalue"),
                        error = function(e) NULL)
    if (!is.null(p_cnet3)) save_plot_dir(p_cnet3, paths$glob$KEGG$P, paste0("gseKEGG_cnet_", file_tag, ".svg"))

    # KEGG ridgeplot
    kkdf <- as.data.frame(kk2@result)
    if (nrow(kkdf) && "p.adjust" %in% names(kkdf)) {
      kkdf$mlogFDR <- -log10(pmax(1e-300, kkdf$p.adjust))
      rp <- ggplot2::ggplot(kkdf, ggplot2::aes(x = mlogFDR, y = "KEGG", fill = "KEGG")) +
        ggridges::geom_density_ridges(alpha = 0.6, color = NA, scale = 1) +
        ggplot2::theme_minimal() + ggplot2::labs(title = paste0("KEGG enrichment density — ", context),
                                                 x = "-log10(FDR)", y = NULL) +
        ggplot2::theme(legend.position = "none")
      save_plot_dir(rp, paths$glob$KEGG$P, paste0("gseKEGG_ridge_", file_tag, ".svg"), 20, 12)
    }

    # KEGG themes
    k_theme <- tryCatch(build_theme_graph(kk2, 50, 0.2), error = function(e) NULL)
    if (!is.null(k_theme)) {
      save_plot_dir(plot_theme_graph(k_theme, paste0("KEGG themes — ", context)),
                    paths$modules$themes, paste0("KEGG_themes_", file_tag, ".svg"))
    } else {
      write_lines_safe("Theme graph skipped: insufficient terms or edges after similarity threshold.",
                       file.path(paths$modules$themes, paste0("KEGG_themes_", file_tag, "_skipped.txt")))
    }
  } else {
    save_table(tibble::tibble(note = "No KEGG enrichment results (gseKEGG null or empty)."),
               paths$glob$KEGG$T, paste0("gseKEGG_", file_tag, "_empty.csv"))
    log_line(paths, file_tag, "gseKEGG produced no results.")
  }

  ## ----- GLOBAL: Reactome (optional) -----
  have_reactome <- requireNamespace("ReactomePA", quietly = TRUE)
  if (have_reactome) {
    entrez_map_re <- AnnotationDbi::mapIds(orgdb, keys = universe_genes, keytype = "UNIPROT",
                                           column = "ENTREZID", multiVals = "first")
    re_universe <- unname(entrez_map_re[!is.na(entrez_map_re)])

    # ORA (DEG intersect)
    re_deg <- unname(entrez_map_re[match(top_genes, names(entrez_map_re))])
    re_deg <- re_deg[!is.na(re_deg)]
    re_ora <- tryCatch(ReactomePA::enrichPathway(
      gene = re_deg, organism = "mouse", universe = re_universe,
      pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 800
    ), error = function(e) NULL)
    if (!is.null(re_ora) && has_terms(re_ora, 1)) {
      save_table(stable_rank(re_ora@result), paths$glob$Reactome$T, paste0("ReactomeORA_", file_tag, ".csv"))
      save_plot_dir(plot_dot(re_ora, paste0("Reactome ORA — ", context)),
                    paths$glob$Reactome$P, paste0("ReactomeORA_dotplot_", file_tag, ".svg"))
    } else {
      save_table(tibble::tibble(note = "Reactome ORA null/empty"),
                 paths$glob$Reactome$T, paste0("ReactomeORA_", file_tag, "_empty.csv"))
    }

    # GSEA
    re_gene_list <- gene_list[!is.na(entrez_map_re)]
    names(re_gene_list) <- unname(entrez_map_re[!is.na(entrez_map_re)])
    keep_idx_re <- !duplicated(names(re_gene_list))
    re_gene_list <- re_gene_list[keep_idx_re]
    re_gene_list <- re_gene_list[order(re_gene_list, decreasing = TRUE, method = "radix")]

    re_gse <- tryCatch(ReactomePA::gsePathway(
      geneList = re_gene_list, organism = "mouse",
      minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1,
      pAdjustMethod = "BH", eps = 0, nPermSimple = 10000
    ), error = function(e) NULL)

    if (!is.null(re_gse) && has_terms(re_gse, 1)) {
      save_table(stable_rank(re_gse@result), paths$glob$Reactome$T, paste0("ReactomeGSEA_", file_tag, ".csv"))
      save_plot_dir(build_split_dotplot(re_gse, paste0("Reactome GSEA — ", context)),
                    paths$glob$Reactome$P, paste0("ReactomeGSEA_split_", file_tag, ".svg"))
      save_plot_dir(plot_gsea_landscape(re_gse, paste0("Reactome GSEA landscape — ", context)),
                    paths$glob$Reactome$P, paste0("ReactomeGSEA_landscape_", file_tag, ".svg"))

      # Reactome themes
      r_theme <- build_theme_graph(re_gse, 50, 0.2)
      if (!is.null(r_theme)) {
        save_plot_dir(plot_theme_graph(r_theme, paste0("Reactome themes — ", context)),
                      paths$modules$themes, paste0("Reactome_themes_", file_tag, ".svg"))
      } else {
        write_lines_safe("Theme graph skipped: insufficient terms or edges after similarity threshold.",
                         file.path(paths$modules$themes, paste0("Reactome_themes_", file_tag, "_skipped.txt")))
      }
    } else {
      save_table(tibble::tibble(note = "Reactome GSEA null/empty"),
                 paths$glob$Reactome$T, paste0("ReactomeGSEA_", file_tag, "_empty.csv"))
    }
  } else {
    log_line(paths, file_tag, "ReactomePA not available; Reactome enrichment skipped.")
  }

  ## ----- PATHVIEW -----
  pv_ids <- c("mmu04110","mmu04115","mmu04114","mmu04720","mmu04721","mmu04722",
              "mmu04725","mmu04726","mmu04727","mmu04724","mmu04080","mmu00030","mmu04151")
  if (!is.null(kk2) && has_terms(kk2, 1)) {
    pv_ids <- unique(c(head(kk2@result$ID, 10), pv_ids))
  }
  pv_ids <- intersect(pv_ids, valid_kegg_ids)
  save_table(tibble::tibble(pv_ids = pv_ids), paths$pathview, paste0("pathview_ids_", file_tag, ".csv"))

  if (length(kegg_gene_list) > 0 && all(is.finite(as.numeric(kegg_gene_list)))) {
    write_lines_safe(sprintf("n_entrez_in_gene_list=%d", length(kegg_gene_list)),
                     file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    oldwd <- getwd(); setwd(paths$pathview)
    invisible(lapply(pv_ids, function(pid) {
      try(
        pathview::pathview(gene.data = kegg_gene_list, pathway.id = pid, species = "mmu",
                           low = "#6698CC", mid = "white", high = "#F08C21", file.type = "svg"),
        silent = TRUE
      )
    }))
    setwd(oldwd)
    log_line(paths, file_tag, sprintf("Pathview attempted for %d pathways.", length(pv_ids)))
  } else {
    write_lines_safe(c(
      "No valid ENTREZ-mapped genes for Pathview or gene list non-numeric.",
      paste0("length(kegg_gene_list)=", length(kegg_gene_list))
    ), file.path(paths$pathview, paste0("pathview_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Pathview skipped: invalid kegg_gene_list.")
  }

  ## ----- MODULE-AWARE: overlaps, UpSet, Fisher, per-module GO -----
  overlap_counts <- tibble::tibble(
    Module = names(modules_list),
    OverlapRanked = vapply(modules_list, function(g) length(intersect(g, universe_genes)), integer(1)),
    OverlapDEG    = vapply(modules_list, function(g) length(intersect(g, top_genes)), integer(1))
  )
  save_table(overlap_counts, paths$modules$overlap$T, "module_overlap_counts.csv")

  ov <- overlap_counts[order(-overlap_counts$OverlapDEG), , drop = FALSE]
  top_mods <- head(ov$Module[ov$OverlapDEG > 0], 10)
  sets <- list(DEGs = top_genes)
  for (m in top_mods) sets[[paste0("Module_", m)]] <- modules_list[[m]]
  all_ids <- unique(unlist(sets, use.names = FALSE))
  mat <- lapply(sets, function(s) as.integer(all_ids %in% s)) |> as.data.frame()
  colnames(mat) <- names(sets); mat$id <- all_ids

  overlap_nonzero <- sum(vapply(sets, function(s) length(intersect(s, all_ids))>0, logical(1)))
  if (overlap_nonzero >= 2 && sum(colSums(as.matrix(mat[, names(sets), drop = FALSE]))) > 0) {
    if (requireNamespace("ComplexUpset", quietly = TRUE)) {
      library(ComplexUpset)
      set_cols <- setdiff(colnames(mat), "id")
      if (length(set_cols) >= 2) {
        mat[set_cols] <- lapply(mat[set_cols], function(x) as.logical(x))
        p_up <- ComplexUpset::upset(
          data = mat,
          intersect = set_cols,
          base_annotations = list('Intersection size' = ComplexUpset::intersection_size()),
          set_sizes = FALSE,
          min_size = 5
        )
        save_plot_dir(p_up, paths$modules$overlap$P, paste0("UpSet_topModules_", file_tag, ".svg"), 28, 20)
      } else {
        save_table(tibble::tibble(note="Insufficient sets (need >=2 logical sets)."),
                   paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped_too_few_sets.csv"))
      }
    } else if (requireNamespace("UpSetR", quietly = TRUE)) {
      library(UpSetR)
      ensure_parent(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")))
      grDevices::pdf(file.path(paths$modules$overlap$P, paste0("UpSetR_topModules_", file_tag, ".pdf")), width = 12, height = 8)
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), nintersects = 30)
      grDevices::dev.off()
    } else {
      save_table(tibble::tibble(note="No UpSet library available (ComplexUpset/UpSetR)."),
                 paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped_no_pkg.csv"))
    }
  } else {
    save_table(tibble::tibble(note="Insufficient non-zero sets for UpSet; plot skipped."),
               paths$modules$overlap$T, paste0("UpSet_", file_tag, "_skipped.csv"))
    log_line(paths, file_tag, "UpSet skipped due to zero overlaps.")
  }

  ov2 <- ov[ov$OverlapRanked > 0 | ov$OverlapDEG > 0, , drop = FALSE]
  if (nrow(ov2) > 0 && is_nonempty_num(ov2$OverlapDEG)) {
    ov2$Pct_in_module_DEG <- round(100 * ov2$OverlapDEG / pmax(1, ov2$OverlapRanked), 1)
    pbar <- ggplot2::ggplot(head(ov2, 20),
                            ggplot2::aes(x = forcats::fct_reorder(Module, OverlapDEG), y = OverlapDEG)) +
      ggplot2::geom_col(fill = "#4575b4") +
      ggplot2::coord_flip() +
      ggplot2::geom_text(ggplot2::aes(label = paste0(OverlapDEG, " (", Pct_in_module_DEG, "%)")),
                         hjust = -0.1, size = 3) +
      ggplot2::ylim(0, max(head(ov2, 20)$OverlapDEG)*1.15 + 1) +
      ggplot2::labs(title = paste0("DEG overlap with modules — ", context),
                    x = "Module", y = "DEG overlap (count)")
    save_plot_dir(pbar, paths$modules$overlap$P, paste0("Module_DEG_overlap_bar_", file_tag, ".svg"), 24, 20)
  } else {
    write_lines_safe("No non-zero overlaps to plot.",
                     file.path(paths$modules$overlap$P, paste0("overlap_status_", file_tag, ".txt")))
  }

  fisher_deg_enrichment <- function(deg_genes, universe_genes, modules_list) {
    deg_set <- intersect(unique(deg_genes), universe_genes)
    bg_set  <- setdiff(universe_genes, deg_set)
    out <- lapply(names(modules_list), function(m) {
      mod_genes <- intersect(modules_list[[m]], universe_genes)
      a <- length(intersect(mod_genes, deg_set))
      b <- length(setdiff(mod_genes, deg_set))
      c <- length(setdiff(deg_set, mod_genes))
      d <- length(setdiff(bg_set, mod_genes))
      mm <- matrix(c(a,b,c,d), nrow=2)
      ft <- fisher.test(mm, alternative = "greater")
      data.frame(
        Module = m, InModule_DE = a, InModule_NotDE = b,
        NotInModule_DE = c, NotInModule_NotDE = d,
        ModuleSize = length(mod_genes), DEG_in_module_pct = ifelse((a+b)>0, 100*a/(a+b), NA_real_),
        OddsRatio = unname(ft$estimate), Pvalue = ft$p.value,
        stringsAsFactors = FALSE
      )
    })
    res <- dplyr::bind_rows(out)
    res$Padj <- p.adjust(res$Pvalue, method = "BH")
    dplyr::arrange(res, Padj, Pvalue)
  }

  fisher_tbl <- fisher_deg_enrichment(deg_genes = top_genes,
                                      universe_genes = universe_genes,
                                      modules_list = modules_list)
  save_table(fisher_tbl, paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv"))

  if (nrow(fisher_tbl) > 0 && is_nonempty_num(fisher_tbl$Padj)) {
    ord_idx <- order(fisher_tbl$Padj, fisher_tbl$Pvalue)
    ft_plot <- ggplot2::ggplot(head(fisher_tbl[ord_idx, , drop = FALSE], 20),
                               ggplot2::aes(x = forcats::fct_reorder(Module, -log10(Padj)),
                                            y = -log10(Padj))) +
      ggplot2::geom_col(fill = "#6a51a3") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = paste0("Fisher enrichment (DEGs in modules) — ", context),
                    x = "Module", y = "-log10(Padj)")
    save_plot_dir(ft_plot, paths$modules$fisher$P, paste0("fisher_bar_", file_tag, ".svg"), 24, 20)
  } else {
    write_lines_safe("No Fisher rows to plot.",
                     file.path(paths$modules$fisher$P, paste0("fisher_status_", file_tag, ".txt")))
    log_line(paths, file_tag, "Fisher plot skipped: empty table.")
  }

  ## ----- Per-module analyses -----
  sel_modules <- if (nrow(fisher_tbl) > 0) fisher_tbl$Module[order(fisher_tbl$Padj, fisher_tbl$Pvalue)] else character(0)
  for (m in sel_modules) {
    m_genes <- intersect(modules_list[[m]], universe_genes)
    if (length(m_genes) < 10) {
      log_line(paths, file_tag, paste("Module", m, "skipped: <10 genes in universe."))
      next
    }

    for (ont in onts) {
      m_root      <- mk(results_dir, "20_modules", paste0("GO_", ont), paste0("Module_", m))
      m_T_gsea    <- mk(m_root, "Tables", "gsea")
      m_P_gsea    <- mk(m_root, "Plots",  "gsea")
      m_T_oradeg  <- mk(m_root, "Tables", "ora_deg")
      m_P_oradeg  <- mk(m_root, "Plots",  "ora_deg")
      m_T_orafull <- mk(m_root, "Tables", "ora_full")
      m_P_orafull <- mk(m_root, "Plots",  "ora_full")

      m_gene_list <- gene_list[names(gene_list) %in% m_genes]
      m_gene_list <- m_gene_list[order(m_gene_list, decreasing = TRUE, method = "radix")]
      if (length(m_gene_list) >= 10) {
        m_gse <- muffle_fgsea_ties( tryCatch(clusterProfiler::gseGO(
          geneList = m_gene_list, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          by = "fgsea", nPermSimple = 10000, eps = 0,
          verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH"
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "gseGO error:", e$message)); NULL }) )
        if (!is.null(m_gse) && has_terms(m_gse, 1)) {
          resm <- stable_rank(m_gse@result); resm$note <- "fgsea ties handled by deterministic id sort"
          save_table(resm, m_T_gsea, paste0("gseGO_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(build_split_dotplot(m_gse, paste0("Module ", m, " — GSEA (", ont, ") — ", context)),
                        m_P_gsea, paste0("gseGO_", ont, "_", m, "_split_", file_tag, ".svg"))
          p_dot <- tryCatch(enrichplot::dotplot(m_gse, x = "GeneRatio", color = "p.adjust",
                                                showCategory = min(10, nrow(m_gse@result)), label_format = 30), error = function(e) NULL)
          if (!is.null(p_dot)) save_plot_dir(p_dot, m_P_gsea, paste0("gseGO_", ont, "_", m, "_dot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no GSEA terms"), m_T_gsea,
                     paste0("gseGO_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "Module GSEA skipped: <10 genes."), m_T_gsea,
                   paste0("gseGO_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      m_top <- intersect(top_genes, m_genes)
      if (length(m_top) >= 5) {
        m_ora <- tryCatch(clusterProfiler::enrichGO(
          gene = m_top, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "none", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA DEG error:", e$message)); NULL })
        if (!is.null(m_ora) && has_terms(m_ora, 1)) {
          save_table(stable_rank(m_ora@result), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora, paste0("Module ", m, " — ORA (DEG intersect, ", ont, ") — ", context)),
                        m_P_oradeg, paste0("ORA_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (DEG intersect)"), m_T_oradeg,
                     paste0("ORA_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(
          note = "insufficient DEGs for ORA",
          n_module_genes_in_universe = length(m_genes),
          n_deg_in_module = length(intersect(top_genes, m_genes)),
          deg_threshold = 1
        ), m_T_oradeg, paste0("ORA_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }

      m_all <- m_genes
      if (length(m_all) >= 5) {
        m_ora_all <- tryCatch(clusterProfiler::enrichGO(
          gene = m_all, ont = ont, keyType = "UNIPROT",
          minGSSize = 5, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = organism, pAdjustMethod = "BH", universe = universe_genes
        ), error = function(e) { log_line(paths, file_tag, paste("Module", m, "ORA full error:", e$message)); NULL })
        if (!is.null(m_ora_all) && has_terms(m_ora_all, 1)) {
          save_table(stable_rank(m_ora_all@result), m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, ".csv"))
          save_plot_dir(plot_dot(m_ora_all, paste0("Module ", m, " — ORA full set (", ont, ") — ", context)),
                        m_P_orafull, paste0("ORA_full_", ont, "_", m, "_dotplot_", file_tag, ".svg"))
        } else {
          save_table(tibble::tibble(note = "no ORA terms (full module)"), m_T_orafull,
                     paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_empty.csv"))
        }
      } else {
        save_table(tibble::tibble(note = "ORA full module skipped: <5 genes."),
                   m_T_orafull, paste0("ORA_fullModule_", ont, "_", m, "_", file_tag, "_skipped.csv"))
      }
    } # ont
  } # module
    ## ----- EXECUTIVE SUMMARY (robust) -----
  summary_items <- list()
  summary_items$universe <- data.frame(
    metric = c("universe_size","deg_abs_gt1","pct_tied"),
    value  = c(length(universe_genes), length(top_genes), round(pct_tied,2))
  )
  # Top GO/KEGG/Reactome tables if they exist
  get_first_csv <- function(path, pattern) {
    fps <- list.files(path, pattern = pattern, full.names = TRUE)
    if (length(fps)) tryCatch(readr::read_csv(fps[1], show_col_types = FALSE), error = function(e) NULL) else NULL
  }
  summary_items$top_go_gsea <- get_first_csv(paths$glob$BP$T, paste0("^gseGO_BP_", file_tag, "\\.csv$"))
  summary_items$top_kegg_gsea <- get_first_csv(paths$glob$KEGG$T, paste0("^gseKEGG_", file_tag, "\\.csv$"))
  summary_items$top_reactome_gsea <- get_first_csv(paths$glob$Reactome$T, paste0("^ReactomeGSEA_", file_tag, "\\.csv$"))
  summary_items$top_modules_fisher <- {
    ft <- tryCatch(readr::read_csv(file.path(paths$modules$fisher$T, paste0("fisher_", file_tag, ".csv")), show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(ft)) ft[order(ft$Padj, ft$Pvalue), , drop = FALSE] else NULL
  }

  safe_univ_val <- function(x, metric, default = NA_character_) {
    if (is.null(x) || !nrow(x) || !"metric" %in% names(x) || !"value" %in% names(x)) return(default)
    v <- x$value[x$metric == metric]
    if (length(v) == 0) return(default)
    as.character(v[[1]])
  }
  fmt_top <- function(df, fmt) {
    if (is.null(df) || !nrow(df)) return(character(0))
    need <- intersect(c("ID","Description","p.adjust"), names(df))
    if (length(need) < 3) return(character(0))
    df <- df[order(df$p.adjust), , drop = FALSE]
    apply(head(df, 5), 1, function(r) sprintf(fmt, r[["ID"]], r[["Description"]], as.numeric(r[["p.adjust"]])))
  }

  md_lines <- list(
    "# Executive Summary",
    paste0("- Context: ", context),
    paste0("- Universe size: ", safe_univ_val(summary_items$universe, "universe_size", "NA")),
    paste0("- DEGs |log2fc|>1: ", safe_univ_val(summary_items$universe, "deg_abs_gt1", "NA")),
    paste0("- pct tied stats: ", safe_univ_val(summary_items$universe, "pct_tied", "NA"), "%"),
    "",
    "## Top GO GSEA (BP)"
  )
  go_lines <- fmt_top(summary_items$top_go_gsea, "%s — %s (FDR=%.2g)")
  md_lines <- c(md_lines, if (length(go_lines)) paste0("* ", go_lines) else "No terms", "",
                "## Top KEGG GSEA")
  kegg_lines <- fmt_top(summary_items$top_kegg_gsea, "%s — %s (FDR=%.2g)")
  md_lines <- c(md_lines, if (length(kegg_lines)) paste0("* ", kegg_lines) else "No terms", "",
                "## Top Reactome GSEA")
  reactome_lines <- fmt_top(summary_items$top_reactome_gsea, "%s — %s (FDR=%.2g)")
  md_lines <- c(md_lines, if (length(reactome_lines)) paste0("* ", reactome_lines) else "Skipped/None", "",
                "## Top Modules by Fisher (Padj)")
  if (!is.null(summary_items$top_modules_fisher) && nrow(summary_items$top_modules_fisher)) {
    mod_tbl <- head(summary_items$top_modules_fisher, 5)
    if (all(c("Module","OddsRatio","Padj") %in% names(mod_tbl))) {
      mod_lines <- apply(mod_tbl, 1, function(r) sprintf("%s — OR=%.2f (Padj=%.2g)",
                                                         r[["Module"]], as.numeric(r[["OddsRatio"]]), as.numeric(r[["Padj"]])))
      md_lines <- c(md_lines, paste0("* ", mod_lines))
    } else {
      md_lines <- c(md_lines, "No modules (columns missing).")
    }
  } else {
    md_lines <- c(md_lines, "No modules")
  }
  md <- unlist(md_lines, use.names = FALSE)
  writeLines(md, con = file.path(paths$summary, paste0("SUMMARY_", file_tag, ".md")))

  writeLines(paste("End:", Sys.time()), con = file.path(paths$meta, "run_end.txt"))
  log_line(paths, file_tag, "Completed all analyses.")
} # end for data_path






