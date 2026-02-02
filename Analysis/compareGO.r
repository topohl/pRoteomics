#' compareGO.r - Comparative Gene Ontology Enrichment Analysis and Visualization
#'
#' @description
#' This script performs comparative GO enrichment analysis across multiple experiments.
#' It:
#'   - Reads multiple CSV files with enrichment data from a specified directory.
#'     Filenames (without extensions) are used as labels for each comparison.
#'   - Combines all data into a single data frame, adding a "Comparison" column.
#'   - Selects top enriched terms per comparison based on absolute NES (Normalized Enrichment Score).
#'   - Filters the combined dataset to include only these top terms for consistent visualization.
#'   - Reorders comparisons in plots based on the maximum absolute NES.
#'   - Adds a significance label ("✱") for gene sets with adjusted p-value < 0.05.
#'
#' @details
#' The script generates:
#'   1. A heatmap showing differential enrichment across comparisons.
#'   2. A dot plot visualizing NES (color) and significance (-log10(p.adjust)) for each term.
#' Additional outputs:
#'   - Core gene lists per term and comparison.
#'   - A binary matrix of core gene presence/absence for Jaccard similarity calculation.
#'   - Expanded core enrichment heatmaps for individual terms.
#'
#' @section File Inputs:
#'   - CSV files from:
#'     "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/core_enrichment"
#'     (pattern "*.csv").
#'
#' @section Outputs:
#'   - Heatmaps and dot plots of enrichment profiles.
#'   - Core gene tables (CSV).
#'   - SVG and PNG files for all plots.
#'
#' @note
#'   Ensure all file paths and required packages are set up before running.
#'
#' @author
#'   Tobias Pohl

# -----------------------------------------------------
# Load Required Libraries
# -----------------------------------------------------

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
# Ensure 'rlang' version >= 1.1.7 is installed before running this script.
if (!requireNamespace("rlang", quietly=TRUE) || packageVersion("rlang") < "1.1.7") {
    stop("Please install 'rlang' version >= 1.1.7 manually before running this script.")
}
if (!requireNamespace("simplifyEnrichment", quietly=TRUE)) {
    BiocManager::install("simplifyEnrichment", force = TRUE)
}
library(simplifyEnrichment)

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, stringr, ggpubr, ggthemes, dplyr, tidyr, purrr,
  readr, pheatmap, tibble, tidyverse, RColorBrewer, writexl, scales, ggrepel)

# -----------------------------------------------------
# Define Theme and Helper Functions
# -----------------------------------------------------
#' Custom ggplot2 theme for consistent plots
theme_custom <- function(base_size = 8, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(face = "bold", size = base_size + 1),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(0.8, "lines"),
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

#' Dynamically calculate plot dimensions based on data content
calc_dims <- function(df_plot) {
  n_cols <- length(unique(as.character(df_plot$Comparison)))
  n_rows <- length(unique(as.character(df_plot$Description)))
  w <- max(4, 4 + (n_cols * 0.4)) 
  h <- max(6, 2 + (n_rows * 0.3))
  return(list(w = w, h = h))
}


#' Simplify GO terms from multiple enrichment result lists
simplifyGOFromMultipleLists_custom <- function (lt, go_id_column = NULL, padj_column = NULL, padj_cutoff = 0.01,
    filter = function(x) any(x < padj_cutoff), default = 1, ont = NULL,
    db = "org.Hs.eg.db", measure = "Sim_XGraSM_2013", heatmap_param = list(NULL),
    show_barplot = TRUE, method = "binary_cut", control = list(),
    min_term = NULL, verbose = TRUE, column_title = NULL, ...)
{
    n = length(lt)
    if (is.data.frame(lt[[1]])) {
        if (is.null(go_id_column)) {
            go_id_column = which(sapply(lt[[1]], function(x) all(grepl("^GO:\\d+$",
                x))))[1]
            if (length(go_id_column) == 0) {
                if (!is.null(rownames(lt[[1]]))) {
                  go_id_column = rownames
                  if (is.null(rownames(lt[[1]]))) {
                    stop_wrap("Cannot find the GO ID column in the data frames. Please explicitly set argument `go_id_column`.")
                  }
                  if (verbose) {
                    message_wrap("Use row names of the data frame as `go_id_column`.")
                  }
                }
                else {
                  stop_wrap("Cannot find the GO ID column in the data frames. Please explicitly set argument `go_id_column`.")
                }
            }
            else {
                if (verbose) {
                  message_wrap(qq("Use column '@{colnames(lt[[1]])[go_id_column]}' as `go_id_column`."))
                }
            }
        }
        if (is.null(padj_column)) {
            cn = colnames(lt[[1]])
            ind = test_padj_column(cn)
            if (length(ind)) {
                padj_column = ind
                if (verbose) {
                  message_wrap(qq("Use column '@{colnames(lt[[1]])[padj_column]}' as `padj_column`."))
                }
            }
            else {
                stop_wrap("Cannot find the column the contains adjusted p-values in the data frames. Please explicitly set argument `padj_column`.")
            }
        }
        lt = lapply(lt, function(x) {
            if (is.function(go_id_column)) {
                structure(x[, padj_column], names = go_id_column(x))
            }
            else {
                structure(x[, padj_column], names = x[, go_id_column])
            }
        })
        return(simplifyGOFromMultipleLists(lt, padj_cutoff = padj_cutoff,
            filter = filter, default = default, ont = ont, db = db,
            measure = measure, heatmap_param = heatmap_param,
            show_barplot = show_barplot, method = method, control = control,
            min_term = min_term, verbose = verbose, column_title = column_title,
            ...))
    }
    else if (is.character(lt[[1]])) {
        lt = lapply(lt, function(x) structure(rep(1, length(x)),
            names = x))
        return(simplifyGOFromMultipleLists(lt, default = 0, filter = function(x) TRUE,
            ont = ont, db = db, measure = measure, show_barplot = show_barplot,
            method = method, heatmap_param = list(transform = function(x) x,
                breaks = c(0, 1), col = c("transparent", "#2f97ffbe"),
                name = "", labels = c("not available", "available")),
            control = control, min_term = min_term, verbose = verbose,
            column_title = column_title, ...))
    }
    heatmap_param2 = list(transform = NULL, breaks = NULL, col = NULL,
        labels = NULL, name = "padj")
    for (nm in names(heatmap_param)) {
        heatmap_param2[[nm]] = heatmap_param[[nm]]
    }
    transform = heatmap_param2$transform
    if (is.null(transform))
        transform = function(x) -log10(x)
    breaks = heatmap_param2$breaks
    col = heatmap_param2$col
    labels = heatmap_param2$labels
    name = heatmap_param2$name
    if (is.null(name))
        name = ""
    if (is.null(breaks) && is.null(col)) {
        digit = ceiling(-log10(padj_cutoff))
        base = padj_cutoff * 10^digit
        breaks = c(1, padj_cutoff, base * 10^(-digit * 2))
        col = c("green", "white", "red")
        labels = gt_render(c("1", qq("@{base}x10<sup>-@{digit}</sup>"),
            qq("@{base}x10<sup>-@{digit*2}</sup>")))
    }
    else if (!is.null(breaks) && !is.null(col)) {
        if (length(breaks) != length(col)) {
            stop_wrap("Length of `breaks` must be the same as the length of `col`.")
        }
    }
    all_go_id = unique(unlist(lapply(lt, names)))
    if (!all(grepl("^GO:\\d+$", all_go_id))) {
        stop_wrap("Only GO ID is allowed.")
    }
    m = matrix(default, nrow = length(all_go_id), ncol = n)
    rownames(m) = all_go_id
    colnames(m) = names(lt)
    if (is.null(colnames))
        colnames = paste0("Group", 1:n)
    for (i in 1:n) {
        m[names(lt[[i]]), i] = lt[[i]]
    }
    l = apply(m, 1, function(x) {
        if (all(is.na(x))) {
            FALSE
        }
        else {
            l = filter(x[!is.na(x)])
            if (length(l) == 1) {
                return(l)
            }
            else {
                return(any(l))
            }
        }
    })
    m = m[l, , drop = FALSE]
    m = t(apply(m, 1, transform))
    if (verbose)
        message(qq("@{nrow(m)}/@{length(all_go_id)} GO IDs left for clustering."))
    if (length(unique(m[!is.na(m)])) <= 2) {
        col = structure(col, names = breaks)
    }
    else {
        if (is.null(breaks) && is.null(col)) {
            col = NULL
        }
        else if (!is.null(breaks) && !is.null(col)) {
            if (length(breaks) != length(col)) {
                stop_wrap("Length of `breaks` and `col` should be the same.")
            }
            col = colorRamp2(transform(breaks), col)
        }
        else {
            stop_wrap("Arguments `breaks` and `col` should be set at the same time.")
        }
    }
    all_go_id = rownames(m)
    sim_mat = GO_similarity(all_go_id, ont = ont, db = db, measure = measure)
    all_go_id = rownames(sim_mat)
    heatmap_legend_param = list()
    heatmap_legend_param$at = transform(breaks)
    heatmap_legend_param$labels = if (is.null(labels))
        breaks
    else labels
    heatmap_legend_param$title = name
    mm = m[all_go_id, , drop = FALSE]
    if (show_barplot) {
        draw_ht = function(align_to) {
            s = sapply(align_to, function(index) max(apply(mm[index,
                ], 2, function(x) sum(x >= transform(padj_cutoff)))))
            max = max(s)
            by = diff(grid.pretty(c(0, max)))[1]
            Heatmap(mm, col = col, name = if (name == "") 
                NULL
            else name, show_row_names = FALSE, cluster_columns = FALSE,
                border = "black", heatmap_legend_param = heatmap_legend_param,
                width = unit(0.5, "cm") * n, use_raster = TRUE,
                left_annotation = rowAnnotation(empty = anno_block(width = unit(1.2,
                  "cm"), panel_fun = function(index) grid.text(qq("Number of significant GO terms in each cluster (padj < @{padj_cutoff})"),
                  unit(0, "npc"), 0.5, just = "top", rot = 90,
                  gp = gpar(fontsize = 10))), bar = anno_link(align_to = align_to,
                  side = "left", gap = unit(3, "mm"), link_gp = gpar(fill = "#DDDDDD",
                    col = "#AAAAAA"), internal_line = FALSE,
                  panel_fun = function(index) {
                    v = apply(mm[index, ], 2, function(x) sum(x >=
                      transform(padj_cutoff)))
                    grid.text(v[2])
                    pushViewport(viewport())
                    grid.rect(gp = gpar(fill = "#DDDDDD", col = "#DDDDDD"))
                    grid.lines(c(1, 0, 0, 1), c(0, 0, 1, 1),
                      gp = gpar(col = "#AAAAAA"), default.units = "npc")
                    pushViewport(viewport(xscale = c(0.5, length(v) +
                      0.5), yscale = c(0, max(v)), height = unit(1,
                      "npc") - unit(2, "mm")))
                    grid.rect(seq_along(v), 0, width = 0.6, height = unit(v,
                      "native"), default.units = "native", just = "bottom",
                      gp = gpar(fill = "#444444", col = "#444444"))
                    if (length(index)/nrow(mm) > 0.05) {
                      grid.yaxis(at = seq(0, max(v), by = by),
                        gp = gpar(col = "#444444", cex = 0.6))
                    }
                    popViewport()
                    popViewport()
                  }, size = s/sum(s) * (unit(1, "npc") - unit(3,
                    "mm") * (length(align_to) - 1) - unit(2,
                    "mm") * length(align_to)) + unit(2, "mm"))),
                post_fun = function(ht) {
                  decorate_annotation("bar", {
                    nc = ncol(mm)
                    grid.text(colnames(mm), (seq_len(nc) - 0.5)/nc *
                      (unit(1, "npc") - unit(5, "mm")), y = -ht_opt$COLUMN_ANNO_PADDING,
                      default.units = "npc", just = "right",
                      rot = 90)
                  })
                })
        }
    }
    else {
        draw_ht = Heatmap(mm, col = col, name = if (name == "")
            NULL
        else name, show_row_names = FALSE, cluster_columns = FALSE,
            border = "black", heatmap_legend_param = heatmap_legend_param,
            width = unit(0.5, "cm") * n, use_raster = TRUE)
    }
    if (is.null(min_term))
        min_term = round(nrow(sim_mat) * 0.02)
    if (is.null(column_title))
        column_title = qq("@{length(all_go_id)} GO terms clustered by '@{method}'")
    simplifyGO(sim_mat, ht_list = draw_ht, method = method, verbose = verbose,
        min_term = min_term, control = control, column_title = column_title,
        ...)
}

# -----------------------------------------------------
# Set Analysis Parameters and Directory Structure
# -----------------------------------------------------

# Gene Ontology domain (MF, BP, or CC)
ont <- "BP"  # Biological Process

# Ensemble profiling method
ensemble_profiling <- "learning_signature"

# Experimental condition
condition <- "effects_inhibition_memory_ensemble"

# Base project path for all input/output
base_project_path <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"

# Set working directory
setwd(base_project_path)

uniprot_mapping_file_path <- file.path(
  base_project_path,
  "Datasets",
  "MOUSE_10090_idmapping.dat"
)

# Load clusterProfiler and organism database
if (!require("clusterProfiler")) BiocManager::install("clusterProfiler")
if (!require("org.Mm.eg.db")) BiocManager::install("org.Mm.eg.db")
library(clusterProfiler)
library(org.Mm.eg.db)

# Read UniProt mapping file (UniprotID <-> Gene_Name)
uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  filter(V2 == "Gene_Name") %>%
  dplyr::select(UniprotID = V1, Gene_Name = V3)

# -----------------------------------------------------
# Import Enrichment Data Files
# -----------------------------------------------------

# List all enrichment CSV files for the specified condition
input_dir <- file.path(base_project_path, "Datasets", "core_enrichment", ont, ensemble_profiling, condition)
file_paths <- list.files(
  path = input_dir,
  pattern = "*.csv",
  full.names = TRUE
)
names(file_paths) <- basename(file_paths) %>% str_remove(".csv")

# -----------------------------------------------------
# Define Output Directories
# -----------------------------------------------------

main_output_dir <- file.path(base_project_path, "Results", "compareGO", ont, ensemble_profiling, condition)
dir.create(main_output_dir, showWarnings = FALSE, recursive = TRUE)

subdirs <- list(
  tables = file.path(main_output_dir, "01_Tables_and_Data"),
  plots_main = file.path(main_output_dir, "02_Main_Plots"),
  core_enrichment_plots = file.path(main_output_dir, "03_Core_Enrichment_Heatmaps"),
  gene_lists = file.path(main_output_dir, "04_Gene_Lists"),
  volcanoes = file.path(main_output_dir, "05_Volcano_Plots"),
  sig_proteins = file.path(main_output_dir, "06_Significant_Proteins"),
  go_enrichment = file.path(main_output_dir, "07_Regulated_Protein_GO"),
  gene_centric = file.path(main_output_dir, "08_Gene_Centric_Analysis")
)
lapply(subdirs, dir.create, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------
# Read and Combine Enrichment Data
# -----------------------------------------------------

enrichment_list <- lapply(file_paths, read.csv)
names(enrichment_list) <- names(file_paths)

combined_df <- bind_rows(
  lapply(names(enrichment_list), function(name) {
    df <- enrichment_list[[name]]
    df$Comparison <- name
    return(df)
  })
)

# Remove whitespace from comparison names to avoid duplicates
combined_df$Comparison <- str_trim(combined_df$Comparison)

# Save combined enrichment data to Excel
writexl::write_xlsx(
  combined_df, 
  path = file.path(subdirs$tables, paste0("Combined_Enrichment_Data.xlsx"))
)

# Extract GO IDs for each gene from core_enrichment column
gene_go_ids <- combined_df %>%
  dplyr::select(GO_ID = ID, core_enrichment) %>%
  mutate(core_enrichment = strsplit(core_enrichment, "/")) %>%
  unnest(core_enrichment) %>%
  dplyr::rename(Gene = core_enrichment) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(GO_IDs = paste(unique(GO_ID), collapse = "; "), .groups = "drop")

# -----------------------------------------------------
# Select Top GO Terms for Comparison
# -----------------------------------------------------

# Parameters for top term selection
significant_only <- TRUE      
top10_terms <- TRUE           
top10_per_comp <- TRUE        
top5_up_down_per_comp <- TRUE 

# Helper: Remove redundant terms based on gene overlap (Jaccard index)
check_redundancy <- function(terms_df, combined_df, threshold = 0.7) {
  get_genes_for_term <- function(desc, df) {
    row_data <- df %>% 
      filter(Description == desc) %>% 
      arrange(desc(abs(NES))) %>% 
      dplyr::slice(1)
    if(nrow(row_data) == 0) return(character(0))
    unlist(strsplit(as.character(row_data$core_enrichment), "/"))
  }
  unique_candidates <- unique(terms_df$Description)
  term_gene_map <- lapply(unique_candidates, function(term) {
    sort(unique(get_genes_for_term(term, combined_df)))
  })
  names(term_gene_map) <- unique_candidates
  kept_terms <- character(0)
  for (candidate in unique_candidates) {
    candidate_genes <- term_gene_map[[candidate]]
    is_redundant <- FALSE
    for (accepted in kept_terms) {
      accepted_genes <- term_gene_map[[accepted]]
      intersect_len <- length(intersect(candidate_genes, accepted_genes))
      union_len <- length(union(candidate_genes, accepted_genes))
      jaccard <- ifelse(union_len > 0, intersect_len / union_len, 0)
      if (jaccard > threshold) {
        is_redundant <- TRUE
        break
      }
    }
    if (!is_redundant) {
      kept_terms <- c(kept_terms, candidate)
    }
  }
  return(kept_terms)
}

# Store 5-up-5-down selection for separate plotting
top5_up_down_df <- NULL

# Select top terms per comparison (non-redundant, up/down split)
if (top5_up_down_per_comp) {
  n_per_direction <- 10  # 10 up + 10 down = 20 in total per comparison
  if (significant_only) {
    df_sig <- combined_df %>% filter(p.adjust < 0.05)
  } else {
    df_sig <- combined_df
  }
  pick_top_n_non_redundant <- function(sub_df, n_target, direction = "up") {
    if (direction == "up") {
      sorted_cand <- sub_df %>% arrange(desc(NES))
    } else {
      sorted_cand <- sub_df %>% arrange(NES) 
    }
    candidates <- sorted_cand$Description
    term_signatures <- list()
    for(desc in unique(candidates)) {
      row_data <- sub_df %>% filter(Description == desc) %>% dplyr::slice(1)
      term_signatures[[desc]] <- sort(unique(unlist(strsplit(as.character(row_data$core_enrichment), "/"))))
    }
    kept <- character(0)
    for (term in candidates) {
      if (length(kept) >= n_target) break
      if (term %in% kept) next 
      curr_genes <- term_signatures[[term]]
      is_redundant <- FALSE
      for (k in kept) {
        k_genes <- term_signatures[[k]]
        intersect_len <- length(intersect(curr_genes, k_genes))
        union_len <- length(union(curr_genes, k_genes))
        jac <- ifelse(union_len > 0, intersect_len / union_len, 0)
        if (jac > 0.7) {
          is_redundant <- TRUE
          break
        }
      }
      if (!is_redundant) {
        kept <- c(kept, term)
      }
    }
    sub_df %>% filter(Description %in% kept) %>% 
      group_by(Description) %>% 
      filter(row_number() == 1) %>% 
      ungroup()
  }

  unique_comps <- unique(df_sig$Comparison)
  results_list <- list()

  for (comp in unique_comps) {
    message(paste("Processing redundancy check for comparison:", comp))
    comp_df <- df_sig %>% filter(Comparison == comp)
    up_candidates <- comp_df %>% filter(NES > 0)
    top_up <- pick_top_n_non_redundant(up_candidates, n_per_direction, "up")
    down_candidates <- comp_df %>% filter(NES < 0)
    top_down <- pick_top_n_non_redundant(down_candidates, n_per_direction, "down")
    results_list[[comp]] <- bind_rows(top_up, top_down)
  }
  
  top_df_per_comp <- bind_rows(results_list) %>%
    distinct(Comparison, ID, .keep_all = TRUE)
  top5_up_down_df <- top_df_per_comp
} else if (top10_per_comp) {
  if (significant_only) {
    top_df_per_comp <- combined_df %>%
      filter(p.adjust < 0.05) %>%
      group_by(Comparison) %>%
      arrange(desc(abs(NES))) %>%
      dplyr::slice_head(n = 10) %>%
      ungroup()
  } else {
    top_df_per_comp <- combined_df %>%
      group_by(Comparison) %>%
      arrange(desc(abs(NES))) %>%
      dplyr::slice_head(n = 10) %>%
      ungroup()
  }
} else {
  if (significant_only) {
    top_df_per_comp <- combined_df %>%
      filter(p.adjust < 0.05)
  } else {
    top_df_per_comp <- combined_df
  }
}
# -----------------------------------------------------
# Select Overall Top Terms (Standard and Refined)
# -----------------------------------------------------
# For each comparison, select top 5 up and top 5 down terms (by NES, non-redundant).
# "Standard": just top/bottom 5 by NES per comparison.
# "Refined": non-redundant (by Jaccard), still top/bottom 5 per comparison.

target_n_terms <- 5

# Helper: get genes for a term
get_term_genes <- function(desc, df) {
  row_data <- df %>% 
    filter(Description == desc) %>% 
    arrange(desc(abs(NES))) %>% 
    dplyr::slice(1)
  genes <- unlist(strsplit(as.character(row_data$core_enrichment), "/"))
  return(sort(unique(genes)))
}

# Non-redundant top-N selection by Jaccard
pick_top_n_non_redundant <- function(sub_df, n_target, direction = "up") {
  if (direction == "up") {
    sorted_cand <- sub_df %>% arrange(desc(NES))
  } else {
    sorted_cand <- sub_df %>% arrange(NES)
  }
  candidates <- sorted_cand$Description
  term_signatures <- list()
  for(desc in unique(candidates)) {
    row_data <- sub_df %>% filter(Description == desc) %>% dplyr::slice(1)
    term_signatures[[desc]] <- sort(unique(unlist(strsplit(as.character(row_data$core_enrichment), "/"))))
  }
  kept <- character(0)
  for (term in candidates) {
    if (length(kept) >= n_target) break
    if (term %in% kept) next 
    curr_genes <- term_signatures[[term]]
    is_redundant <- FALSE
    for (k in kept) {
      k_genes <- term_signatures[[k]]
      intersect_len <- length(intersect(curr_genes, k_genes))
      union_len <- length(union(curr_genes, k_genes))
      jac <- ifelse(union_len > 0, intersect_len / union_len, 0)
      if (jac > 0.7) {
        is_redundant <- TRUE
        break
      }
    }
    if (!is_redundant) {
      kept <- c(kept, term)
    }
  }
  sub_df %>% filter(Description %in% kept) %>% 
    group_by(Description) %>% 
    filter(row_number() == 1) %>% 
    ungroup()
}

# --- STANDARD: top/bottom 5 by NES per comparison (no redundancy check) ---
if (significant_only) {
  candidate_pool_df <- combined_df %>%
    filter(p.adjust < 0.05) %>%
    arrange(desc(abs(NES)))
} else {
  candidate_pool_df <- combined_df %>%
    arrange(desc(abs(NES)))
}

unique_comps <- unique(candidate_pool_df$Comparison)
results_list_up_std <- list()
results_list_down_std <- list()
results_list_up_ref <- list()
results_list_down_ref <- list()

for (comp in unique_comps) {
  comp_df <- candidate_pool_df %>% filter(Comparison == comp)
  # Standard: just top/bottom 5 by NES
  up_std <- comp_df %>% filter(NES > 0) %>% arrange(desc(NES)) %>% dplyr::slice_head(n = target_n_terms)
  down_std <- comp_df %>% filter(NES < 0) %>% arrange(NES) %>% dplyr::slice_head(n = target_n_terms)
  results_list_up_std[[comp]] <- up_std
  results_list_down_std[[comp]] <- down_std
  # Refined: non-redundant top/bottom 5 by NES
  up_ref <- pick_top_n_non_redundant(comp_df %>% filter(NES > 0), target_n_terms, "up")
  down_ref <- pick_top_n_non_redundant(comp_df %>% filter(NES < 0), target_n_terms, "down")
  results_list_up_ref[[comp]] <- up_ref
  results_list_down_ref[[comp]] <- down_ref
}

top_df_up_std <- bind_rows(results_list_up_std)
top_df_down_std <- bind_rows(results_list_down_std)
top_df_up_ref <- bind_rows(results_list_up_ref)
top_df_down_ref <- bind_rows(results_list_down_ref)

top_terms_up_std <- unique(top_df_up_std$Description)
top_terms_down_std <- unique(top_df_down_std$Description)
top_terms_up_ref <- unique(top_df_up_ref$Description)
top_terms_down_ref <- unique(top_df_down_ref$Description)

top_terms_standard <- unique(c(top_terms_up_std, top_terms_down_std))
top_terms_refined <- unique(c(top_terms_up_ref, top_terms_down_ref))

# Set default top_terms for downstream plotting (use refined)
top_terms <- top_terms_refined

# -----------------------------------------------------
# Plot Comparison: Standard vs Refined Selection
# -----------------------------------------------------

# Custom color palette: upregulation = #faa51a, downregulation = #4c87c6, 0 = white
custom_palette <- colorRampPalette(c("#4c87c6", "white", "#faa51a"))

plot_data_standard <- combined_df %>%
  filter(Description %in% top_terms_standard) %>%
  mutate(Selection_Type = "Standard")
plot_data_refined <- combined_df %>%
  filter(Description %in% top_terms_refined) %>%
  mutate(Selection_Type = "Refined")
comparison_plot_data <- bind_rows(plot_data_standard, plot_data_refined) %>%
  mutate(Comparison = factor(Comparison)) 

max_abs_nes <- max(abs(comparison_plot_data$NES), na.rm = TRUE)

comp_plot <- ggplot(comparison_plot_data, aes(
  x = Comparison,
  y = reorder(Description, NES, FUN = mean),
  color = NES,
  size = -log10(p.adjust)
)) +
  geom_point() +
  facet_grid(Selection_Type ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradientn(
    colours = custom_palette(100),
    name = "NES",
    limits = c(-max_abs_nes, max_abs_nes)
  ) +
  scale_size_continuous(range = c(2, 6)) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) + 
  labs(
    title = "Selection method comparison",
    y = "GO Term",
    x = "Comparison"
  ) +
  theme_custom(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"))

n_cols_comp <- length(unique(as.character(comparison_plot_data$Comparison)))
dims_comp <- list()
dims_comp$w <- max(4, 6 + (n_cols_comp * 0.6)) 
dims_comp$h <- calc_dims(comparison_plot_data)$h

ggsave(
  filename = file.path(subdirs$plots_main, "Comparison_Selection_Methods_Dotplot.svg"),
  plot = comp_plot,
  width = dims_comp$w, 
  height = dims_comp$h
)

# -----------------------------------------------------
# Prepare Data for Heatmap Visualization
# -----------------------------------------------------

lookup_df <- combined_df %>%
  filter(Description %in% top_terms)
comparison_order <- lookup_df %>%
  group_by(Comparison) %>%
  summarize(max_abs_NES = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_abs_NES)) %>%
  pull(Comparison)
lookup_df <- lookup_df %>%
  mutate(Comparison = factor(Comparison, levels = comparison_order))

# -----------------------------------------------------
# Generate Enrichment Heatmap
# -----------------------------------------------------

lookup_df <- lookup_df %>%
  mutate(sig_label = ifelse(p.adjust < 0.05, "*", ""))

heatmap_data <- lookup_df %>%
  dplyr::select(Description, Comparison, NES) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = NES)

heatmap_data <- as.data.frame(heatmap_data)
rownames(heatmap_data) <- heatmap_data$Description
heatmap_data <- heatmap_data[, -which(names(heatmap_data) == "Description")]
heatmap_data <- as.matrix(heatmap_data)
heatmap_data[is.na(heatmap_data)] <- 0 
heatmap_data_export <- tibble::rownames_to_column(as.data.frame(heatmap_data), var = "RowNames")
writexl::write_xlsx(heatmap_data_export, path = file.path(subdirs$tables, "Matrix_Heatmap_NES.xlsx"))
heatmap_labels <- lookup_df %>%
  dplyr::select(Description, Comparison, sig_label) %>%
  tidyr::pivot_wider(names_from = Comparison, values_from = sig_label) %>%
  tibble::column_to_rownames("Description")
plot_height <- max(6, nrow(heatmap_data) * 0.3)
max_abs_val <- max(abs(heatmap_data), na.rm=TRUE)
my_breaks <- seq(-max_abs_val, max_abs_val, length.out = 100)

# Custom color palette: downregulation = #4c87c6, 0 = white, upregulation = #faa51a
my_colors <- colorRampPalette(c("#4c87c6", "white", "#faa51a"))(100)

heatmap_plot <- pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = heatmap_labels,
  number_color = "black",
  color = my_colors,
  breaks = my_breaks,
  fontsize = 8,
  fontsize_number = 7,
  border_color = "white",
  cellwidth = 20,
  cellheight = 10,
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 45,
  silent = TRUE,
  legend = TRUE
)
heatmap_row_order <- rownames(heatmap_data)[heatmap_plot$tree_row$order]
heatmap_col_order <- colnames(heatmap_data)[heatmap_plot$tree_col$order]
output_file <- file.path(subdirs$plots_main, "Heatmap_Enrichment_Comparisons.svg")
svg(output_file, width = 6, height = plot_height, family = "sans")
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# -----------------------------------------------------
# Generate Dot Plots for Enrichment Results
# -----------------------------------------------------

# Helper: Order dotplot axes by clustering
order_dotplot <- function(df, desc_col = "Description", comp_col = "Comparison", val_col = "NES") {
  mat <- df %>%
    dplyr::select(all_of(c(desc_col, comp_col, val_col))) %>%
    filter(!is.na(!!sym(val_col))) %>%
    pivot_wider(names_from = all_of(comp_col), values_from = all_of(val_col), values_fill = 0) %>%
    column_to_rownames(desc_col) %>%
    as.matrix()
  hc_col <- hclust(dist(t(mat)))
  col_order <- colnames(mat)[hc_col$order]
  dom_comp <- apply(mat, 1, function(x) colnames(mat)[which.max(abs(x))])
  dom_val  <- apply(mat, 1, function(x) {
    val <- x[which.max(abs(x))]
    return(val)
  })
  ordering_df <- data.frame(
    Term = rownames(mat),
    Comp = dom_comp,
    Val = dom_val
  )
  ordering_df$Comp <- factor(ordering_df$Comp, levels = col_order)
  ordering_df <- ordering_df %>%
    mutate(
      Comp_Index = as.integer(Comp),
      Direction = ifelse(Val > 0, 1, 0),
      Abs_Val = abs(Val)
    ) %>%
    arrange(desc(Comp_Index), Direction, Abs_Val)
  row_order <- ordering_df$Term
  return(list(row_order = row_order, col_order = col_order))
}

# Dotplot: Top terms per comparison
df_per_comp_sub <- combined_df %>%
  filter(Description %in% unique(top_df_per_comp$Description))
ordering_res <- order_dotplot(df_per_comp_sub)
full_data_for_plot <- df_per_comp_sub %>%
  mutate(
    Comparison = factor(Comparison, levels = ordering_res$col_order),
    Description = factor(Description, levels = ordering_res$row_order)
  )

dotplot <- ggplot(full_data_for_plot, aes(
  x = Comparison,
  y = Description, 
  color = NES,
  size = -log10(p.adjust)
)) +
  geom_point(alpha = 1) +
  scale_color_gradientn(
    colours = custom_palette(100),
    name = "NES",
    limits = c(-max(abs(full_data_for_plot$NES), na.rm=TRUE), max(abs(full_data_for_plot$NES), na.rm=TRUE))
  ) +
  scale_size_continuous(
    name = expression(-log[10](italic(P)[adj])),
    range = c(1.5, 5),
    limits = c(0, NA)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), drop = FALSE) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(
    title = "Top terms per comparison",
    x = NULL,
    y = NULL
  ) +
  theme_custom(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color="black"),
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  )

# Dotplot: Overall top terms (refined)
df_top_sub <- combined_df %>%
  filter(Description %in% top_terms)
ordering_res_top <- order_dotplot(df_top_sub)
top_terms_all_data <- df_top_sub %>%
  mutate(
    Comparison = factor(Comparison, levels = ordering_res$col_order),
    Description = factor(Description, levels = ordering_res$row_order)
  )

dotplot_top <- ggplot(top_terms_all_data, aes(
  x = Comparison,
  y = Description, 
  color = NES,
  size = -log10(p.adjust)
)) +
  geom_point(alpha = 1) +
  scale_color_gradientn(
    colours = custom_palette(100),
    name = "NES",
    limits = c(-max(abs(top_terms_all_data$NES), na.rm=TRUE), max(abs(top_terms_all_data$NES), na.rm=TRUE))
  ) +
  scale_size_continuous(
    name = expression(-log[10](italic(P)[adj])),
    range = c(1.5, 5),
    limits = c(0, NA)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), drop = FALSE) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(
    title = "Top terms (overall)",
    x = NULL,
    y = NULL
  ) +
  theme_custom(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color="black"),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2)
  )

dims_dotplot <- calc_dims(full_data_for_plot)
dims_dotplot_top <- calc_dims(top_terms_all_data)
output_dotplot <- file.path(subdirs$plots_main, "Dotplot_Enrichment_TopGenes_PerComp.svg")
output_dotplot_top <- file.path(subdirs$plots_main, "Dotplot_Enrichment_TopGenes_Overall.svg")
ggsave(output_dotplot, plot = dotplot, width = dims_dotplot$w, height = dims_dotplot$h, dpi = 300)
ggsave(output_dotplot_top, plot = dotplot_top, width = dims_dotplot_top$w, height = dims_dotplot_top$h, dpi = 300)

# Dotplot: Top 5 up/down terms per comparison (if selected)
if (!is.null(top5_up_down_df) && nrow(top5_up_down_df) > 0) {
  selected_descriptions <- unique(top5_up_down_df$Description)
  df_split_sub <- combined_df %>%
    filter(Description %in% selected_descriptions)
  ordering_res_split <- order_dotplot(df_split_sub)
  plot_data_split <- df_split_sub %>%
    mutate(
      Comparison = factor(Comparison, levels = ordering_res_split$col_order),
      Description = factor(Description, levels = ordering_res_split$row_order)
    )
  max_abs_nes_split <- max(abs(plot_data_split$NES), na.rm = TRUE)
  dotplot_split <- ggplot(plot_data_split, aes(
    x = Comparison,
    y = Description, 
    color = NES,
    size = -log10(p.adjust)
  )) +
    geom_point(alpha = 1) +
    scale_color_gradientn(
      colours = custom_palette(100),
      name = "NES",
      limits = c(-max_abs_nes_split, max_abs_nes_split)
    ) +
    scale_size_continuous(
      name = expression(-log[10](italic(P)[adj])),
      range = c(1.5, 5),
      limits = c(0, NA)
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.1, 0.1)), drop = FALSE) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) +
    labs(
      title = "Top 5 up/down terms",
      x = NULL,
      y = NULL
    ) +
    theme_custom(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color="black"),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2)
    )
  dims_split <- calc_dims(plot_data_split)
  output_dotplot_split <- file.path(subdirs$plots_main, "Dotplot_Enrichment_Top5_Up_Down_PerComp.svg")
  ggsave(output_dotplot_split, plot = dotplot_split, width = dims_split$w, height = dims_split$h, dpi = 300)
}

# -----------------------------------------------------
# Plot Enrichment for Selected Genes
# -----------------------------------------------------

gene_list <- c("P55099", "Q6NXX1")
long_df <- combined_df %>%
  mutate(core_gene = strsplit(as.character(core_enrichment), "/|;|,|\\s+")) %>%
  unnest(core_gene)
filtered_df <- long_df %>%
  filter(core_gene %in% gene_list) %>%
  mutate(Description = factor(Description, levels = unique(Description)))
nes_min <- min(filtered_df$NES, na.rm = TRUE)
nes_max <- max(filtered_df$NES, na.rm = TRUE)
max_abs_nes <- max(abs(nes_min), abs(nes_max))
plot_selected_genes <- ggplot(filtered_df, aes(
  x = Comparison,
  y = Description,
  color = NES,
  size = -log10(p.adjust),
  shape = core_gene
)) +
  geom_point(alpha = 1) +
  scale_color_gradientn(
    colours = custom_palette(100),
    limits = c(-max_abs_nes, max_abs_nes),
    name = "NES"
  ) +
  scale_size_continuous(
    name = expression(-log[10](italic(P)[adj])),
    range = c(2, 6)
  ) +
  labs(
    title = "Enrichment for selected genes",
    x = NULL,
    y = NULL,
    shape = "Gene"
  ) +
  theme_custom(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
output_plot_file <- file.path(subdirs$plots_main, "Plot_Selected_Genes.svg")
ggsave(output_plot_file, plot = plot_selected_genes, width = 8, height = 6, dpi = 300)

# -----------------------------------------------------
# Gene-Centric Enrichment Summary Plot
# -----------------------------------------------------

gene_enrichment_summary <- long_df %>%
  filter(core_gene %in% gene_list) %>%
  group_by(core_gene, Comparison) %>%
  summarise(
    count = n(),
    mean_NES = mean(NES, na.rm = TRUE),
    max_padj = min(p.adjust, na.rm = TRUE),
    .groups = "drop"
  )

plot_gene_summary <- ggplot(gene_enrichment_summary, aes(
  x = Comparison,
  y = core_gene,
  size = count,
  fill = mean_NES
)) +
  geom_point(shape = 21, stroke = 0.2, color="black") +
  scale_size_continuous(
    name = "Occurrences",
    range = c(2, 8)
  ) +
  scale_fill_gradientn(
    colours = rev(brewer.pal(n = 7, name = "RdBu")),
    name = "Mean NES",
    limits = c(-max_abs_nes, max_abs_nes)
  ) +
  labs(
    title = "Gene-centric enrichment",
    x = NULL,
    y = "Gene"
  ) +
  theme_custom(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

gene_summary_file <- file.path(subdirs$gene_centric, "Gene_Centric_Summary.svg")
ggsave(gene_summary_file, plot = plot_gene_summary, width = 5, height = 4, dpi = 300)

# -----------------------------------------------------
# Extract and Export Core Genes per Comparison
# -----------------------------------------------------

core_gene_sets <- lapply(names(enrichment_list), function(name) {
  df <- enrichment_list[[name]]
  df_long <- df %>%
    dplyr::select(Description, core_enrichment) %>%
    mutate(core_enrichment = str_split(core_enrichment, "/")) %>%
    unnest(core_enrichment) %>%
    mutate(Comparison = name)
  return(df_long)
})

names(core_gene_sets) <- names(enrichment_list)
core_genes_df <- bind_rows(core_gene_sets)
write.csv(core_genes_df,
          file = file.path(subdirs$gene_lists, "Core_Genes_All_Comparisons.csv"),
          row.names = FALSE)
core_genes_df %>%
  group_by(Comparison, Description) %>%
  summarise(Genes = list(unique(core_enrichment)), .groups = "drop") %>%
  rowwise() %>%
  mutate(file_name = file.path(
    subdirs$gene_lists,
    paste0("CoreGene_", Comparison, "_", substr(make.names(Description), 1, 50), ".csv")
  )) %>%
  pwalk(~ write.csv(data.frame(Gene = ..3), file = ..4, row.names = FALSE))

# -----------------------------------------------------
# Compute Jaccard Similarity of Core Genes Across Comparisons
# -----------------------------------------------------

core_gene_sets_aggregated <- core_genes_df %>%
  group_by(Comparison) %>%
  summarise(Genes = list(unique(core_enrichment)))
all_genes <- unique(unlist(core_gene_sets_aggregated$Genes))
binary_matrix <- sapply(core_gene_sets_aggregated$Genes, function(gene_list) all_genes %in% gene_list)
rownames(binary_matrix) <- all_genes
colnames(binary_matrix) <- core_gene_sets_aggregated$Comparison
jaccard_similarity <- function(x, y) {
  intersect = sum(x & y)
  union = sum(x | y)
  if (union == 0) return(NA)
  return(intersect / union)
}
jaccard_matrix <- outer(1:ncol(binary_matrix), 1:ncol(binary_matrix), Vectorize(function(i, j) {
  jaccard_similarity(binary_matrix[, i], binary_matrix[, j])
}))
rownames(jaccard_matrix) <- colnames(binary_matrix)
colnames(jaccard_matrix) <- colnames(binary_matrix)

# Choose color for numbers based on background
# If similarity > 0.5 (darker blue), use white, else black
number_col_mat <- matrix(ifelse(jaccard_matrix > 0.5, "white", "black"),
                         nrow = nrow(jaccard_matrix), ncol = ncol(jaccard_matrix),
                         dimnames = dimnames(jaccard_matrix))

jaccard_heatmap <- pheatmap(
  jaccard_matrix,
  main = "Jaccard similarity",
  color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
  display_numbers = TRUE,
  number_color = number_col_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 14,
  fontsize_number = 12, 
  border_color = "white"
)
jaccard_sim_matrix <- file.path(subdirs$plots_main, "Heatmap_Jaccard_Similarity.svg")
svg(jaccard_sim_matrix, width = 7, height = 6, family = "sans")
grid::grid.draw(jaccard_heatmap$gtable)
dev.off()

# -----------------------------------------------------
# Generate Heatmap of Core Enrichment Genes
# -----------------------------------------------------

core_long_df <- bind_rows(
  lapply(names(enrichment_list), function(name) {
    df <- enrichment_list[[name]]
    df %>%
      dplyr::select(Description, NES, core_enrichment) %>%
      mutate(core_enrichment = strsplit(core_enrichment, "/")) %>%
      tidyr::unnest(core_enrichment) %>%
      dplyr::rename(Gene = core_enrichment) %>%
      mutate(Comparison = name)
  })
)

heatmap_df <- core_long_df %>%
  group_by(Gene, Comparison) %>%
  summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Comparison, values_from = NES)
heatmap_matrix <- heatmap_df %>%
  column_to_rownames("Gene") %>%
  as.matrix()

heatmap_matrix[is.na(heatmap_matrix)] <- 0
nes_min <- min(heatmap_matrix, na.rm = TRUE)
nes_max <- max(heatmap_matrix, na.rm = TRUE)
max_abs_nes <- max(abs(nes_min), abs(nes_max))
breaks <- seq(-max_abs_nes, max_abs_nes, length.out = 101)
# Use the custom_palette defined earlier for color
heatmap_plot <- pheatmap(
  heatmap_matrix,
  color = custom_palette(100),
  breaks = breaks,
  main = "Core enrichment heatmap",
  fontsize = 12,
  fontsize_row = 1, # effectively hides row labels (protein names)
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA
)
output_file <- file.path(subdirs$plots_main, "Heatmap_Overall_Core_Enrichment.svg")
svg(output_file, width = 8, height = 8, family = "sans")
grid::grid.draw(heatmap_plot$gtable)
dev.off()

# Identify top 25 up/down genes by NES for supplementary tables/plots
long_heatmap_df <- heatmap_df %>%
  pivot_longer(-Gene, names_to = "Comparison", values_to = "NES")
top25_up <- long_heatmap_df %>%
  group_by(Gene) %>%
  summarize(max_nes = max(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_nes)) %>%
  dplyr::slice_head(n = 25) %>%
  mutate(Direction = "Up")
top25_down <- long_heatmap_df %>%
  group_by(Gene) %>%
  summarize(min_nes = min(NES, na.rm = TRUE), .groups = "drop") %>%
  arrange(min_nes) %>%
  dplyr::slice_head(n = 25) %>%
  mutate(Direction = "Down")
top_bottom_genes <- bind_rows(top25_up, top25_down)
gene_go <- core_long_df %>%
  dplyr::select(Gene, Description) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(GO_Terms = paste(unique(Description), collapse = "; "), .groups = "drop")
top_bottom_genes <- top_bottom_genes %>%
  left_join(gene_go, by = "Gene")
write_xlsx(
  top_bottom_genes,
  file.path(subdirs$gene_lists, "TopBottom_GeneList.xlsx")
)

# Generate heatmap for top 25 up/down genes
top_bottom_matrix <- heatmap_matrix[rownames(heatmap_matrix) %in% top_bottom_genes$Gene, , drop = FALSE]
gene_order <- heatmap_df %>%
  filter(Gene %in% rownames(top_bottom_matrix)) %>%
  pivot_longer(-Gene, names_to = "Comparison", values_to = "NES") %>%
  group_by(Gene) %>%
  summarize(order_metric = max(abs(NES), na.rm = TRUE)) %>%
  arrange(desc(order_metric)) %>%
  pull(Gene)

top_bottom_matrix <- top_bottom_matrix[gene_order, , drop = FALSE]
uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  filter(V2 == "Gene_Name") %>%
  dplyr::select(UniprotID = V1, Gene_Name = V3)
mapped_names <- left_join(
  tibble(Gene = rownames(top_bottom_matrix)),
  uniprot_df,
  by = c("Gene" = "UniprotID")
)

rownames(top_bottom_matrix) <- ifelse(
  is.na(mapped_names$Gene_Name),
  mapped_names$Gene,
  mapped_names$Gene_Name
)

nes_min <- min(top_bottom_matrix, na.rm = TRUE)
nes_max <- max(top_bottom_matrix, na.rm = TRUE)
max_abs_nes <- max(abs(nes_min), abs(nes_max))
breaks <- seq(-max_abs_nes, max_abs_nes, length.out = 101)
top_bottom_plot <- pheatmap(
  top_bottom_matrix,
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
  breaks = breaks,
  main = "Core enrichment overview",
  fontsize = 8,
  fontsize_row = 6,
  fontsize_col = 8,
  fontsize_number = 6,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  border_color = "white",
  border_width = 0.5,
  cellwidth = 10,
  cellheight = 7,
  angle_col = 45,
  treeheight_row = 15,
  treeheight_col = 15,
  legend = TRUE,
  legend_breaks = c(round(-max_abs_nes,1), 0, round(max_abs_nes,1)),
  show_rownames = TRUE,
  show_colnames = TRUE
)
svg(file.path(subdirs$plots_main, "Heatmap_TopBottom_Enrichment.svg"), width = 4, height = 8, family="sans")
grid::grid.draw(top_bottom_plot$gtable)
dev.off()

# -----------------------------------------------------
# Generate Per-Comparison Heatmaps and Excel Files
# -----------------------------------------------------

uniprot_df <- read.delim(
  uniprot_mapping_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

uniprot_subset <- uniprot_df %>%
  filter(V2 %in% c("Gene_Name", "UniProtKB-ID")) %>%
  pivot_wider(names_from = V2, values_from = V3, values_fn = list) %>%
  unnest(cols = everything()) %>%
  distinct(V1, .keep_all = TRUE) %>%
  dplyr::rename(UniprotID = V1)
gene_synonyms <- uniprot_df %>%
  filter(V2 == "Gene_Synonym") %>%
  distinct(V1, .keep_all = TRUE) %>%
  dplyr::rename(UniprotID = V1, Gene_Synonym = V3)
gene_descriptions <- core_long_df %>%
  dplyr::select(Gene, Description) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarize(Description = paste(unique(Description), collapse = "; "), .groups = "drop")
log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)

log2fc_df_list <- lapply(log2fc_files, function(f) {
  df <- read_csv(f, col_types = cols())
  df$Comparison <- tools::file_path_sans_ext(basename(f))
  df <- df %>%
    dplyr::rename(
      padj = adj.P.Val,
      log2fc = logFC
    ) %>%
    as.data.frame()
  return(df)
})

log2fc_long_df <- bind_rows(log2fc_df_list)
unique_comparisons <- unique(log2fc_long_df$Comparison)
for (comp in unique_comparisons) {
  full_comp_df <- log2fc_long_df %>%
    filter(Comparison == comp, !is.na(log2fc), !is.na(padj))
  # Volcano plot for each comparison
  volcano_df <- full_comp_df %>%
    filter(!str_detect(gene_symbol, "_")) %>%
    left_join(uniprot_subset %>% dplyr::select(UniprotID, Gene_Name),
              by = c("gene_symbol" = "UniprotID")) %>%
    mutate(Significance = case_when(
      padj < 0.05 & log2fc > 0 ~ "up",
      padj < 0.05 & log2fc < 0 ~ "down",
      TRUE ~ "n.s."
    ))
  top5_up <- volcano_df %>%
    filter(Significance == "up") %>%
    arrange(desc(log2fc)) %>%
    dplyr::slice_head(n = 5)
  top5_down <- volcano_df %>%
    filter(Significance == "down") %>%
    arrange(log2fc) %>%
    dplyr::slice_head(n = 5)
  top_genes_to_label <- bind_rows(top5_up, top5_down)
  max_x_val <- max(abs(volcano_df$log2fc), 0, na.rm = TRUE)
  x_limit <- ceiling(max_x_val) + 0.5
  max_y_val <- max(abs(-log10(volcano_df$padj)), 0, na.rm = TRUE)
  y_limit <- ceiling(max_y_val) + 1
  volcano_plot <- ggplot(volcano_df, aes(x = log2fc, y = -log10(padj), label = Gene_Name)) +
    geom_point(aes(color = Significance), shape = 16, alpha = 0.6, size = 6) + # larger dots, 0.6 alpha
    ggrepel::geom_text_repel(
      data = top_genes_to_label,
      aes(label = Gene_Name),
      size = 6, # larger font
      min.segment.length = 0,
      max.overlaps = Inf,
      box.padding = 0.7
    ) +
    scale_color_manual(
      values = c("up" = "#f36e21", "down" = "#465b65", "n.s." = "#BABABA"),
      breaks = c("up", "down", "n.s.")
    ) +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    labs(
      title = paste("Volcano plot:", comp),
      x = expression(log[2] ~ Fold ~ Change),
      y = expression(-log[10](italic(P)[adj]))
    ) +
    theme_custom(base_size = 16) + # larger base font
    theme(
      legend.position = "none"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30", linewidth = 0.3)
  volcano_file <- file.path(subdirs$volcanoes, paste0("Volcano_", comp, ".svg"))
  message(paste("Saving volcano plot to:", volcano_file))
  ggsave(
    filename = volcano_file,
    plot = volcano_plot,
    width = 5,
    height = 5
  )
  comp_df <- full_comp_df %>% 
    filter(padj < 0.05) %>%
    ungroup()
  if (nrow(comp_df) == 0) {
    message(paste("  No significant genes found for", comp, "- skipping heatmap and Excel."))
    next
  }
  comp_df <- comp_df %>%
    mutate(Direction = case_when(
      log2fc > 0 ~ "Up",
      log2fc < 0 ~ "Down",
      TRUE ~ "Neutral"
    ))
  top_bottom_genes <- bind_rows(
    comp_df %>%
      filter(Direction == "Up") %>%
      arrange(desc(log2fc)) %>%
      dplyr::slice_head(n = 25),
    comp_df %>%
      filter(Direction == "Down") %>%
      arrange(log2fc) %>%
      dplyr::slice_head(n = 25)
  ) %>%
    distinct(gene_symbol, Comparison, .keep_all = TRUE)
  mapped_top_bottom <- top_bottom_genes %>%
    left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
    left_join(gene_synonyms, by = c("gene_symbol" = "UniprotID")) %>%
    left_join(gene_descriptions, by = c("gene_symbol" = "Gene")) %>%
    left_join(gene_go_ids, by = c("gene_symbol" = "Gene")) %>%
    mutate(
      Gene_Name = as.character(Gene_Name),
      `UniProtKB-ID` = as.character(`UniProtKB-ID`),
      Gene_Synonym = as.character(Gene_Synonym),
      Gene_Label = ifelse(is.na(Gene_Name) | Gene_Name == "", gene_symbol, Gene_Name),
      Gene_Label = make.unique(as.character(Gene_Label)),
      abs_log2fc = abs(log2fc),
      UniProt_Link = paste0("https://www.uniprot.org/uniprot/", gene_symbol)
    )
  mapped_top_bottom <- mapped_top_bottom %>%
    group_by(Direction) %>%
    arrange(Direction, desc(abs_log2fc)) %>%
    mutate(Rank = row_number()) %>%
    ungroup()
  heatmap_matrix <- mapped_top_bottom %>%
    dplyr::select(Gene_Label, log2fc) %>%
    tibble::column_to_rownames("Gene_Label") %>%
    as.matrix()
  if (nrow(heatmap_matrix) >= 2) {
    max_abs <- max(abs(heatmap_matrix), na.rm = TRUE)
    breaks <- seq(-max_abs, max_abs, length.out = 101)
    heatmap_plot <- pheatmap(
      heatmap_matrix,
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
      breaks = breaks,
      main = NA,
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize_row = 8,
      show_colnames = FALSE,
      border_color = "white",
      silent = TRUE,
      legend = FALSE
    )
    heatmap_file <- file.path(subdirs$plots_main, paste0("Heatmap_Log2FC_", comp, ".svg"))
    message(paste("  Saving heatmap to:", heatmap_file))
    svg(heatmap_file, width = 2.5, height = 5, family="sans")
    grid::grid.draw(heatmap_plot$gtable)
    dev.off()
  }
  excel_out <- mapped_top_bottom %>%
    dplyr::select(
      Gene_Name,
      `UniProtKB-ID`,
      Uniprot_Accession = gene_symbol,
      Gene_Synonym,
      Comparison,
      log2fc,
      abs_log2fc,
      padj,
      Direction,
      Rank,
      GO_IDs,
      Description,
      UniProt_Link
    ) %>%
    mutate(across(everything(), as.character)) %>%
    distinct()
  excel_file <- file.path(subdirs$sig_proteins, paste0("SigProteins_TopBottom_", comp, ".xlsx"))
  write_xlsx(excel_out, excel_file)
  all_sig_proteins <- comp_df %>%
    left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
    left_join(gene_synonyms, by = c("gene_symbol" = "UniprotID")) %>%
    left_join(gene_descriptions, by = c("gene_symbol" = "Gene")) %>%
    left_join(gene_go_ids, by = c("gene_symbol" = "Gene")) %>%
    mutate(
      Gene_Name = as.character(Gene_Name),
      `UniProtKB-ID` = as.character(`UniProtKB-ID`),
      Gene_Synonym = as.character(Gene_Synonym),
      abs_log2fc = abs(log2fc),
      UniProt_Link = paste0("https://www.uniprot.org/uniprot/", gene_symbol)
    )
  upregulated <- all_sig_proteins %>%
    filter(Direction == "Up") %>%
    arrange(desc(log2fc)) %>%
    mutate(Rank = row_number()) %>%
    dplyr::select(
      Gene_Name,
      `UniProtKB-ID`,
      Uniprot_Accession = gene_symbol,
      Gene_Synonym,
      Comparison,
      log2fc,
      abs_log2fc,
      padj,
      Direction,
      Rank,
      GO_IDs,
      Description,
      UniProt_Link
    ) %>%
    mutate(across(everything(), as.character)) %>%
    distinct()
  if (nrow(upregulated) > 0) {
    up_file <- file.path(subdirs$sig_proteins, paste0("SigProteins_Upregulated_", comp, ".xlsx"))
    write_xlsx(upregulated, up_file)
  }
  downregulated <- all_sig_proteins %>%
    filter(Direction == "Down") %>%
    arrange(log2fc) %>%
    mutate(Rank = row_number()) %>%
    dplyr::select(
      Gene_Name,
      `UniProtKB-ID`,
      Uniprot_Accession = gene_symbol,
      Gene_Synonym,
      Comparison,
      log2fc,
      abs_log2fc,
      padj,
      Direction,
      Rank,
      GO_IDs,
      Description,
      UniProt_Link
    ) %>%
    mutate(across(everything(), as.character)) %>%
    distinct()
  if (nrow(downregulated) > 0) {
    down_file <- file.path(subdirs$sig_proteins, paste0("SigProteins_Downregulated_", comp, ".xlsx"))
    write_xlsx(downregulated, down_file)
  }
}

# -----------------------------------------------------
# Run GO Enrichment on Up/Downregulated Proteins
# -----------------------------------------------------

enrichment_summary_log <- data.frame(
  Comparison = character(),
  Direction = character(),
  Input_Genes_Count = integer(),
  Enriched_Terms_Count = integer(),
  Status = character(),
  Saved_CSV = character(),
  Saved_Plot = character(),
  Timestamp = character(),
  stringsAsFactors = FALSE
)
comparison_files <- list.files(
  path = subdirs$sig_proteins,
  pattern = "SigProteins_(Upregulated|Downregulated)_.*\\.xlsx$",
  full.names = TRUE
)
if (length(comparison_files) > 0) {
  # Extract comparison names robustly from Upregulated or Downregulated files
  comparisons <- unique(gsub("SigProteins_(Upregulated|Downregulated)_|\\.xlsx", "", 
                             basename(comparison_files)))
  message("Found ", length(comparisons), " comparisons to process: ", paste(comparisons, collapse = ", "))

  for (comp in comparisons) {
    message(paste("Processing GO enrichment for comparison:", comp))
    
    # --- Process Upregulated Genes ---
    file_up <- file.path(subdirs$sig_proteins, paste0("SigProteins_Upregulated_", comp, ".xlsx"))
    
    if (file.exists(file_up)) {
      df_up <- readxl::read_xlsx(file_up)
      up_genes <- df_up %>% filter(!is.na(Uniprot_Accession)) %>% pull(Uniprot_Accession) %>% unique()
      
      # Use UNIPROT IDs directly
      if (length(up_genes) > 0) {
        ego_up <- tryCatch(
          enrichGO(
            gene = up_genes,
            OrgDb = org.Mm.eg.db,
            keyType = "UNIPROT",    # Changed from ENTREZID to UNIPROT
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            minGSSize = 10
          ),
          error = function(e) NULL
        )
        
        if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
          out_csv <- file.path(subdirs$go_enrichment, paste0("GO_", comp, "_Upregulated_BP.csv"))
          write.csv(as.data.frame(ego_up), out_csv, row.names = FALSE)
          
          ego_df <- as.data.frame(ego_up) %>%
            dplyr::mutate(Count = as.numeric(Count)) %>%
            dplyr::filter(Count >= 10) %>%
            arrange(p.adjust) %>%
            dplyr::slice_head(n = 15) %>%
            mutate(Description = stringr::str_wrap(Description, width = 50)) %>%
            mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
            arrange(GeneRatio) %>%
            mutate(Description = factor(Description, levels = unique(Description)))
          
          if (nrow(ego_df) > 0) {
            lollipop_palette <- colorRampPalette(c("#4c87c6", "#faa51a"))
            p_lollipop <- ggplot(ego_df, aes(x = GeneRatio, y = Description)) +
              geom_segment(aes(x = 0, xend = GeneRatio, y = Description, yend = Description), color = "#bdbdbd", linewidth = 1) +
              geom_point(aes(size = Count, color = p.adjust), alpha = 0.9) +
              scale_color_gradientn(colours = lollipop_palette(100), name = expression(italic(P)[adj])) +
              scale_size_continuous(range = c(2, 7), name = "Gene Count") +
              scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)), name = "Gene Ratio") +
              labs(
                title = paste0("Upregulated GO - ", comp),
                y = "GO Term"
              ) +
              theme_custom(base_size = 10) +
              theme(
                axis.text.y = element_text(size = 9),
                axis.text.x = element_text(size = 9),
                legend.position = "right"
              )
            
            out_plot <- file.path(subdirs$go_enrichment, paste0("GO_Lollipop_", comp, "_Upregulated_BP.svg"))
            ggsave(out_plot, plot = p_lollipop, width = 7, height = 5)
          }
        }
      }
    } else {
      message("  No Upregulated file found for: ", comp)
    }
    
    # --- Process Downregulated Genes ---
    file_down <- file.path(subdirs$sig_proteins, paste0("SigProteins_Downregulated_", comp, ".xlsx"))
    
    if (file.exists(file_down)) {
      df_down <- readxl::read_xlsx(file_down)
      down_genes <- df_down %>% filter(!is.na(Uniprot_Accession)) %>% pull(Uniprot_Accession) %>% unique()
      
      # Use UNIPROT IDs directly
      if (length(down_genes) > 0) {
        ego_down <- tryCatch(
          enrichGO(
            gene = down_genes,
            OrgDb = org.Mm.eg.db,
            keyType = "UNIPROT",    # Changed from ENTREZID to UNIPROT
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2,
            minGSSize = 10
          ),
          error = function(e) NULL
        )
        
        if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
          out_csv <- file.path(subdirs$go_enrichment, paste0("GO_", comp, "_Downregulated_BP.csv"))
          write.csv(as.data.frame(ego_down), out_csv, row.names = FALSE)
          
          ego_df <- as.data.frame(ego_down) %>%
            dplyr::mutate(Count = as.numeric(Count)) %>%
            dplyr::filter(Count >= 10) %>%
            arrange(p.adjust) %>%
            dplyr::slice_head(n = 15) %>%
            mutate(Description = stringr::str_wrap(Description, width = 50)) %>%
            mutate(GeneRatio = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
            arrange(GeneRatio) %>%
            mutate(Description = factor(Description, levels = unique(Description)))
          
          if (nrow(ego_df) > 0) {
            lollipop_palette <- colorRampPalette(c("#4c87c6", "#faa51a"))
            p_lollipop <- ggplot(ego_df, aes(x = GeneRatio, y = Description)) +
              geom_segment(aes(x = 0, xend = GeneRatio, y = Description, yend = Description), color = "#bdbdbd", linewidth = 1) +
              geom_point(aes(size = Count, color = p.adjust), alpha = 0.9) +
              scale_color_gradientn(colours = lollipop_palette(100), name = expression(italic(P)[adj])) +
              scale_size_continuous(range = c(2, 7), name = "Gene Count") +
              scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)), name = "Gene Ratio") +
              labs(
                title = paste0("Downregulated GO - ", comp),
                y = "GO Term"
              ) +
              theme_custom(base_size = 10) +
              theme(
                axis.text.y = element_text(size = 9),
                axis.text.x = element_text(size = 9),
                legend.position = "right"
              )
            
            out_plot <- file.path(subdirs$go_enrichment, paste0("GO_Lollipop_", comp, "_Downregulated_BP.svg"))
            ggsave(out_plot, plot = p_lollipop, width = 7, height = 5)
          }
        }
      }
    } else {
      message("  No Downregulated file found for: ", comp)
    }
  }
}

# -----------------------------------------------------
# Generate Individual Core Enrichment Heatmaps
# -----------------------------------------------------

log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)
names(log2fc_files) <- basename(log2fc_files) %>%
  str_remove(".csv") %>%
  str_trim()
log2fc_list <- lapply(names(log2fc_files), function(comp) {
  df <- read.csv(log2fc_files[comp])
  df$Comparison <- comp
  df
})
log2fc_df <- bind_rows(log2fc_list)
log2fc_df <- log2fc_df %>%
  dplyr::rename(
    log2fc = logFC,
    pvalue = P.Value,
    padj = adj.P.Val
  )
core_long_df %>%
  group_by(Description) %>%
  group_split() %>%
  walk(function(df_term) {
    description <- unique(df_term$Description)
    df_term_log2fc <- df_term %>%
      left_join(log2fc_df, by = c("Gene" = "gene_symbol", "Comparison"))
    matrix_df <- df_term_log2fc %>%
      dplyr::select(Gene, Comparison, log2fc) %>%
      group_by(Gene, Comparison) %>%
      summarize(log2fc = mean(log2fc, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Comparison, values_from = log2fc) %>%
      column_to_rownames("Gene")
    colnames(matrix_df) <- str_trim(colnames(matrix_df))
    heatmap_matrix <- as.matrix(matrix_df)
    heatmap_matrix[is.na(heatmap_matrix)] <- 0
    safe_desc <- substr(make.names(description), 1, 30)
    filename <- file.path(subdirs$core_enrichment_plots, paste0("Heatmap_Core_", safe_desc, ".svg"))
    cluster_rows_option <- nrow(heatmap_matrix) > 1
    cluster_cols_option <- ncol(heatmap_matrix) > 1
    max_val <- max(abs(heatmap_matrix), na.rm = TRUE)
    breaks <- seq(-max_val, max_val, length.out = 101)
    plot_height <- max(4, nrow(heatmap_matrix) * 0.15 + 1)
    svg(filename, width = 5, height = plot_height, family = "sans")
    pheatmap(
      heatmap_matrix,
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
      breaks = breaks,
      main = substr(description, 1, 40),
      fontsize = 8,
      fontsize_row = 7,
      fontsize_col = 8,
      cluster_rows = cluster_rows_option,
      cluster_cols = cluster_cols_option,
      border_color = "white",
      angle_col = 45,
      cellwidth = 12,
      cellheight = 10,
      legend = TRUE
    )
    dev.off()
  })

# -----------------------------------------------------
# Gene-Centric Log2 Fold-Change (log2FC) Analysis
# -----------------------------------------------------

gene_list <- sort(unique(c(
  "P21279", "P21278", "P51432", "P11881", "P63318",
  "P68404", "P0DP26", "P0DP27", "P0DP28", "P11798",
  "P28652", "P47937", "P47713"
)))
log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)
log2fc_list <- lapply(log2fc_files, function(file) {
  read_csv(file) %>% 
    mutate(
      Comparison = tools::file_path_sans_ext(basename(file))
    )
})
log2fc_long <- bind_rows(log2fc_list) %>%
  dplyr::rename(
    log2fc = logFC,
    pvalue = P.Value,
    padj = adj.P.Val
  )
log2fc_filtered <- log2fc_long %>%
  filter(gene_symbol %in% gene_list)
gene_fc_summary <- log2fc_filtered %>%
  group_by(gene_symbol, Comparison) %>%
  summarise(
    count = n(),
    mean_log2fc = mean(log2fc, na.rm = TRUE),
    min_padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  )
max_abs_fc <- max(abs(gene_fc_summary$mean_log2fc), na.rm = TRUE)
plot_gene_fc <- ggplot(gene_fc_summary, aes(
  x = Comparison,
  y = gene_symbol,
  size = count,
  fill = mean_log2fc
)) +
  geom_point(shape = 21, stroke = 0.2, color="black") +
  scale_size_continuous(guide = "none", range = c(2, 6)) +  
  scale_fill_gradientn(
    colours = custom_palette(100),
    name = "log2 FC",
    limits = c(-max_abs_fc, max_abs_fc)
  ) +
  labs(
    title = "Gene-centric log2 FC",
    x = NULL,
    y = "Gene"
  ) +
  theme_custom(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
fc_plot_file <- file.path(subdirs$gene_centric, paste0("GeneCentric_Log2FC_Bubble_", condition, ".svg"))
ggsave(fc_plot_file, plot = plot_gene_fc, width = 4, height = 5, dpi = 300)
summary_file <- file.path(subdirs$gene_centric, paste0("GeneCentric_Log2FC_Summary.xlsx"))
writexl::write_xlsx(gene_fc_summary, path = summary_file)

# -----------------------------------------------------
# Generate Supplementary Tables (Merged)
# -----------------------------------------------------

supp_enrichment <- combined_df %>%
  dplyr::select(Comparison, ID, Description, setSize, 
                pvalue, p.adjust, qvalue, NES, core_enrichment, rank, leading_edge)
supp_sig_proteins <- log2fc_long %>%
  filter(padj < 0.05) %>%
  left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
  left_join(gene_synonyms, by = c("gene_symbol" = "UniprotID")) %>% 
  left_join(gene_descriptions, by = c("gene_symbol" = "Gene")) %>%
  mutate(Direction = ifelse(log2fc > 0, "Up", "Down")) %>%
  dplyr::select(Comparison, UniprotID = gene_symbol, Gene_Name, Gene_Synonym, Description,
                log2fc, pvalue, padj, Direction) %>%
  arrange(Comparison, padj)
supp_gene_centric <- gene_fc_summary %>%
  left_join(uniprot_subset, by = c("gene_symbol" = "UniprotID")) %>%
  dplyr::rename(Gene_Name = Gene_Name)
supp_nes_matrix <- lookup_df %>%
  dplyr::select(Description, Comparison, NES) %>%
  pivot_wider(names_from = Comparison, values_from = NES)
supp_binary_matrix <- as.data.frame(binary_matrix) %>%
  tibble::rownames_to_column(var = "Gene")
supp_list <- list(
  "GO_Enrichment_Results" = supp_enrichment,
  "Significant_Proteins" = supp_sig_proteins,
  "Gene_Centric_Log2FC" = supp_gene_centric,
  "NES_Matrix_TopTerms" = supp_nes_matrix,
  "Core_Genes_Binary_Matrix" = supp_binary_matrix
)
supp_output_path <- file.path(subdirs$tables, "Supplementary_Data.xlsx")
writexl::write_xlsx(supp_list, path = supp_output_path)

# perform simplifyGOFromMultipleLists()
# -----------------------------------------------------
# Run simplifyGOFromMultipleLists() for each comparison (up/downregulated proteins)
# -----------------------------------------------------

# Prepare GO enrichment results as a list of data frames for each comparison (up/downregulated proteins)
log2fc_files <- list.files(
  path = file.path(base_project_path, "Datasets", "mapped", ensemble_profiling, condition),
  pattern = "*.csv",
  full.names = TRUE
)
go_results_list <- list()
for (f in log2fc_files) {
  comp_name <- tools::file_path_sans_ext(basename(f))
  df <- readr::read_csv(f, show_col_types = FALSE)
  # Upregulated
  up_genes <- df %>% filter(logFC > 0) %>% pull(gene_symbol) %>% unique()
  if (length(up_genes) > 0) {
    ego_up <- tryCatch(
      enrichGO(
        gene = up_genes,
        OrgDb = org.Mm.eg.db,
        keyType = "UNIPROT",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      ),
      error = function(e) NULL
    )
    if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
      ego_up_df <- as.data.frame(ego_up)
      # Only filter after collecting, and keep at least a few terms if possible
      if (nrow(ego_up_df) > 0) {
        filtered <- ego_up_df %>% dplyr::filter(p.adjust < 0.05)
        if (nrow(filtered) < 3 && nrow(ego_up_df) >= 3) {
          # If too few significant, use top 3 by p.adjust
          filtered <- ego_up_df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = 3)
        }
        if (nrow(filtered) > 0) {
          go_results_list[[paste0(comp_name, "_Up")]] <- filtered
        }
      }
    }
  }
  # Downregulated
  down_genes <- df %>% filter(logFC < 0) %>% pull(gene_symbol) %>% unique()
  if (length(down_genes) > 0) {
    ego_down <- tryCatch(
      enrichGO(
        gene = down_genes,
        OrgDb = org.Mm.eg.db,
        keyType = "UNIPROT",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      ),
      error = function(e) NULL
    )
    if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
      ego_down_df <- as.data.frame(ego_down)
      if (nrow(ego_down_df) > 0) {
        filtered <- ego_down_df %>% dplyr::filter(p.adjust < 0.05)
        if (nrow(filtered) < 3 && nrow(ego_down_df) >= 3) {
          filtered <- ego_down_df %>% dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = 3)
        }
        if (nrow(filtered) > 0) {
          go_results_list[[paste0(comp_name, "_Down")]] <- filtered
        }
      }
    }
  }
}
# Remove empty or too small lists
go_results_list <- go_results_list[sapply(go_results_list, nrow) > 0]

# --- simplifyGOFromMultipleLists: one plot for all comparisons ---
if (length(go_results_list) >= 2) {
  simplifygo_plot_file <- file.path(subdirs$plots_main, "simplifyGO_MultipleLists.svg")
  svg(simplifygo_plot_file, width = 10, height = 8, family = "sans")
  simplifyGOFromMultipleLists_custom(
    go_results_list,
    go_id_column = "ID",
    padj_column = "p.adjust",
    padj_cutoff = 0.01,
    ont = "BP",
    db = "org.Mm.eg.db",
    measure = "Sim_XGraSM_2013",
    method = "binary_cut",
    verbose = TRUE,
    column_title = "simplifyGO: All comparisons",
    # Custom color: white (high p-value), orange (low p-value)
    heatmap_param = list(
      col = c("white", "#faa51a"),
      breaks = c(1, 0.0001),
      name = "padj"
    )
  )
  dev.off()
  # Save cluster assignments
  all_go <- purrr::imap_dfr(go_results_list, ~dplyr::mutate(.x, Source = .y))
  go_ids <- unique(all_go$ID)
  mat <- simplifyEnrichment::GO_similarity(go_ids, ont = "BP")
  clust <- simplifyEnrichment::cluster_terms(mat)
  # Only proceed if clust is not NULL, not empty, and matches the number of GO IDs
  if (!is.null(clust) && length(clust) > 0 && length(clust) == length(go_ids)) {
    clust_df <- data.frame(
      ID = names(clust),
      Cluster = as.integer(clust),
      stringsAsFactors = FALSE
    )
    go_annot <- all_go %>%
      dplyr::select(ID, Source) %>%
      distinct() %>%
      group_by(ID) %>%
      summarise(Source = paste(Source, collapse = "; "), .groups = "drop")
    desc_map <- all_go %>% dplyr::select(ID, Description) %>% distinct()
    out_df <- clust_df %>%
      left_join(go_annot, by = "ID") %>%
      left_join(desc_map, by = "ID")
    simplifygo_table_file <- file.path(subdirs$tables, "simplifyGO_MultipleLists_Results.xlsx")
    writexl::write_xlsx(out_df, simplifygo_table_file)
  } else {
    message("No clusters were found by simplifyEnrichment::cluster_terms or cluster size mismatch; skipping cluster assignment export.")
  }
}

message("Analysis complete. Files saved to: ", main_output_dir)





