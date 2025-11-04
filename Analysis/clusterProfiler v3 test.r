# =========================
# Phase 0: Install & load
# =========================
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

# =========================
# Phase 1: Compute (sequential, OrgDb used here)
# =========================
split_at_second_underscore <- function(x) {
  us <- gregexpr("_", x, fixed = TRUE)[[1]]
  if (identical(us[1], -1L)) return(c(x))
  if (length(us) == 1L) return(strsplit(x, "_", fixed = TRUE)[[1]])
  pos2 <- us[2]; c(substr(x, 1, pos2 - 1), substr(x, pos2 + 1, nchar(x)))
}

save_plot <- function(plot, filename, results_dir) {
  ggplot2::ggsave(file.path(results_dir, filename), plot, units = "cm", dpi = 300)
}

plot_dot <- function(dataset, cell_types, results_dir, prefix = NULL) {
  dot_title <- paste("GSEA of", paste(cell_types, collapse = " over "))
  p <- clusterProfiler::dotplot(dataset, showCategory = 10, split = ".sign") +
    facet_wrap(~ .sign, nrow = 1) +
    labs(title = dot_title, x = "Gene Ratio", y = "Gene Set") +
    scale_fill_viridis_c(option = "cividis") +
    theme_minimal(base_size = 12)
  fname <- if (is.null(prefix)) {
    paste0(deparse(substitute(dataset)), "_", paste(cell_types, collapse = "_"), ".svg")
  } else {
    paste0(prefix, "_", paste(cell_types, collapse = "_"), ".svg")
  }
  ggplot2::ggsave(file.path(results_dir, fname), p, units = "cm", dpi = 300)
  p
}

working_dir  <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
file_dir     <- file.path(working_dir, "Datasets", "mapped", "neuron-phenotypeWithinUnit")
results_root <- file.path(working_dir, "Results")
dir.create(results_root, showWarnings = FALSE, recursive = TRUE)

file_list <- list.files(file_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found in: ", file_dir)

organism <- "org.Mm.eg.db"
ont <- "BP"

for (data_path in file_list) {
  file_name   <- basename(data_path)
  name_no_ext <- tools::file_path_sans_ext(file_name)
  file_tag    <- gsub("[^A-Za-z0-9_]", "_", name_no_ext)
  cell_types  <- trimws(split_at_second_underscore(name_no_ext))

  results_dir <- file.path(results_root, file_tag)
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

  df <- read.csv(data_path, header = TRUE, check.names = FALSE)
  if (!"log2fc" %in% names(df)) stop("Missing 'log2fc' in ", file_name)
  colnames(df)[1] <- "gene_symbol"

  original_gene_list <- df$log2fc
  names(original_gene_list) <- sub("-\\d+$", "", df$gene_symbol)
  gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]

  # Use original_gene_list (which already carries names) so names() length matches values
  top_gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)
  top_gene_list <- top_gene_list[!duplicated(names(top_gene_list))]
  top_genes <- names(top_gene_list)[abs(top_gene_list) > 1]
  top_genes <- names(sort(top_gene_list[top_genes], decreasing = TRUE))

  gse <- tryCatch({
    gseGO(
      geneList = gene_list, ont = ont, keyType = "UNIPROT",
      minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1, verbose = TRUE,
      OrgDb = organism, pAdjustMethod = "BH"
    )
  }, error = function(e) NULL)

  if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
    core_dir <- file.path(results_dir, "core_enrichment", ont)
    dir.create(core_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(gse@result, file = file.path(core_dir, paste0("coreEnrichment_", ont, "_", file_tag, ".csv")))
  }

  ids <- tryCatch({
    bitr(names(original_gene_list), fromType = "UNIPROT", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
  }, error = function(e) NULL)

  kk2 <- NULL; kegg_gene_list <- NULL
  if (!is.null(ids) && nrow(ids) > 0) {
    dedup_ids <- ids[!duplicated(ids$UNIPROT), ]
    df2 <- merge(df, dedup_ids, by.x = "gene_symbol", by.y = "UNIPROT")
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

  saveRDS(list(
    file_tag        = file_tag,
    cell_types      = cell_types,
    results_dir     = results_dir,
    gene_list       = gene_list,
    top_genes       = top_genes,
    gse             = gse,
    kegg_gene_list  = kegg_gene_list,
    kk2             = kk2
  ), file = file.path(results_dir, "compute_artifacts.rds"))
}

# =========================
# Phase 2: Rendering (%dopar%, no OrgDb usage)
# =========================
library(foreach); library(doParallel)
n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
render_log <- file.path(results_root, "render_log.txt")
cl <- parallel::makeCluster(n_cores, outfile = render_log)
doParallel::registerDoParallel(cl)

render_pkgs <- c("ggplot2","enrichplot","clusterProfiler","pathview","ggnewscale",
                 "cowplot","ggridges","ggpubr","ggrepel","ggsci","ggthemes",
                 "ggExtra","ggforce","ggplotify","svglite","tidyr","dplyr","tibble")
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages(lapply(render_pkgs, function(p) library(p, character.only = TRUE)))
}))

rds_files <- list.files(results_root, pattern = "compute_artifacts.rds$", recursive = TRUE, full.names = TRUE)
if (length(rds_files) == 0) stop("No compute_artifacts.rds files found; Phase 1 produced none.")

invisible(foreach(
  rds = rds_files,
  .packages = render_pkgs,
  .errorhandling = "pass",
  .inorder = FALSE
) %dopar% {
  obj <- readRDS(rds)
  file_tag        <- obj$file_tag
  cell_types      <- obj$cell_types
  results_dir     <- obj$results_dir
  gene_list       <- obj$gene_list
  top_genes       <- obj$top_genes
  gse             <- obj$gse
  kegg_gene_list  <- obj$kegg_gene_list
  kk2             <- obj$kk2

  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
    plot_dot(gse, cell_types, results_dir)
    p_emap <- enrichplot::emapplot(pairwise_termsim(gse), showCategory = 10)
    save_plot(p_emap, paste0("GSEAemap_", file_tag, ".svg"), results_dir)

    p_cnet <- enrichplot::cnetplot(gse, categorySize = "pvalue", foldChange = gene_list)
    save_plot(p_cnet, paste0("GSEAcnet_", file_tag, ".svg"), results_dir)

    p_ridge <- enrichplot::ridgeplot(gse) +
      labs(x = "Enrichment Distribution", title = "GSEA Ridgeplot") +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    save_plot(p_ridge, paste0("GSEA_Ridgeplot_", file_tag, ".svg"), results_dir)

    if (nrow(gse@result) > 0) {
      p_gsea <- enrichplot::gseaplot(gse, by = "all", title = gse@result$Description[1], geneSetID = 1)
      save_plot(p_gsea, paste0("GSEA_Plot_", file_tag, ".svg"), results_dir)

      top_terms <- head(gse@result$Description, 3)
      p_pmc <- europepmc::pmcplot(top_terms, 2010:2025, proportion = FALSE) +
        labs(title = "Publication Trends for Top Enriched Terms")
      save_plot(p_pmc, paste0("GSEA_PubMed_Trends_", file_tag, ".svg"), results_dir)
    }

    if (!is.null(top_genes) && length(top_genes) > 0) {
      ora <- tryCatch({
        enrichGO(
          gene = top_genes, ont = "CC", keyType = "UNIPROT",
          minGSSize = 3, maxGSSize = 800, pvalueCutoff = 1,
          OrgDb = "org.Mm.eg.db", pAdjustMethod = "none"
        )
      }, error = function(e) NULL)
      if (!is.null(ora) && nrow(as.data.frame(ora)) > 0) {
        p7 <- clusterProfiler::dotplot(ora, showCategory = 10) +
          labs(title = "ORA of Top Regulated Genes") +
          scale_x_continuous(limits = c(0, 1))
        save_plot(p7, paste0("ORA_dotplot_", file_tag, ".svg"), results_dir)
      }
    }
  }

  if (!is.null(kk2) && nrow(as.data.frame(kk2)) > 0) {
    plot_dot(kk2, cell_types, results_dir)
    p_emap2 <- enrichplot::emapplot(pairwise_termsim(kk2), showCategory = 10)
    save_plot(p_emap2, paste0("KEGGemap_", file_tag, ".svg"), results_dir)

    p_cnet2 <- enrichplot::cnetplot(kk2, categorySize = "pvalue", foldChange = gene_list)
    save_plot(p_cnet2, paste0("KEGGcnet_", file_tag, ".svg"), results_dir)

    p_ridge2 <- enrichplot::ridgeplot(kk2) + labs(x = "Enrichment distribution")
    save_plot(p_ridge2, paste0("KEGG_ridge_", file_tag, ".svg"), results_dir)

    if (nrow(kk2@result) > 0) {
      p_gsea2 <- enrichplot::gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)
      save_plot(p_gsea2, paste0("KEGG_gseaplot_", file_tag, ".svg"), results_dir)
    }

    pathview_dir <- normalizePath(file.path(results_dir, "pathview"), winslash = "/", mustWork = FALSE)
    if (!dir.exists(pathview_dir)) dir.create(pathview_dir, recursive = TRUE)
    path_ids <- c("mmu04110","mmu04115","mmu04114","mmu04113","mmu04112","mmu04111",
                  "mmu04116","mmu04117","mmu04118","mmu04119","mmu04720","mmu04721",
                  "mmu04722","mmu04725","mmu04726","mmu04727","mmu04724","mmu04080",
                  "mmu00030","mmu04151")
    oldwd <- getwd()
    setwd(pathview_dir)
    invisible(lapply(path_ids, function(pid) {
      try({
        pathview::pathview(
          gene.data  = kegg_gene_list,
          pathway.id = pid,
          species    = "mmu",
          low  = "#6698CC", mid = "white", high = "#F08C21",
          file.type = "svg"
        )
      }, silent = TRUE)
    }))
    setwd(oldwd)
  }

  TRUE
})  # end foreach
    # end invisible

parallel::stopCluster(cl)
