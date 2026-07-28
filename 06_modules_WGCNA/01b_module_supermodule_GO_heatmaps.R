#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/01b_module_supermodule_GO_heatmaps.R
# Stage: modules_downstream
# Scope: dataset_specific
# Consumes: required Stage 01 module GO enrichment and supermodule membership.
# Produces: representative module and supermodule GO comparison heatmaps.
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia with --dataset.
# Notes: Uses all-module ORA results. Supermodule cells summarize member-module
#        evidence and are not pooled enrichment tests.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_go_comparison_utils.R"))
source(repo_path("R", "validation_utils.R"))

run <- wgcna_cli()
DATASET <- run$dataset
args <- run$args
arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[[1]] == length(args)) return(default)
  args[[hit[[1]] + 1L]]
}
ontology_arg <- toupper(arg_value("--ontology", "BP"))
ONTOLOGIES <- if (identical(ontology_arg, "ALL")) c("BP", "MF", "CC") else unique(trimws(strsplit(ontology_arg, ",", fixed = TRUE)[[1]]))
if (length(setdiff(ONTOLOGIES, c("BP", "MF", "CC")))) stop("--ontology must be BP, MF, CC, a comma-separated subset, or all.", call. = FALSE)
FDR_CUTOFF <- suppressWarnings(as.numeric(arg_value("--fdr", "0.05")))
if (!is.finite(FDR_CUTOFF) || FDR_CUTOFF <= 0 || FDR_CUTOFF >= 1) stop("--fdr must be between 0 and 1.", call. = FALSE)
SCORE_CAP <- suppressWarnings(as.numeric(arg_value("--score-cap", "6")))
if (!is.finite(SCORE_CAP) || SCORE_CAP <= 0) stop("--score-cap must be positive.", call. = FALSE)
TOP_TERMS_PER_MODULE <- suppressWarnings(as.integer(arg_value("--top-terms-per-module", "2")))
TOP_TERMS_PER_SUPERMODULE <- suppressWarnings(as.integer(arg_value("--top-terms-per-supermodule", "3")))
MAX_TERMS <- suppressWarnings(as.integer(arg_value("--max-terms", "24")))
if (any(!is.finite(c(TOP_TERMS_PER_MODULE, TOP_TERMS_PER_SUPERMODULE, MAX_TERMS))) || any(c(TOP_TERMS_PER_MODULE, TOP_TERMS_PER_SUPERMODULE, MAX_TERMS) < 1L)) {
  stop("Term-count parameters must be positive integers.", call. = FALSE)
}

PATHS <- wgcna_downstream_paths("01b_module_supermodule_GO_heatmaps", DATASET)
MODULES_DIR <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", DATASET, "modules")
SUPERMODULES_DIR <- path_results("tables", "06_modules_WGCNA", "01_WGCNA", DATASET, "supermodules")
GO_FILE <- file.path(MODULES_DIR, "WGCNA_module_GO_enrichment_long.csv")
MODULE_MAP_FILE <- file.path(MODULES_DIR, "module_name_map.csv")
MEMBERSHIP_FILE <- file.path(
  SUPERMODULES_DIR,
  if (identical(DATASET, "microglia")) "wgcna_module_supermodule_annotation.csv" else "wgcna_supermodule_eigengene_clusters.csv"
)

if (run$dry_run) {
  invisible(lapply(unlist(PATHS), dir_create))
  dry_run_line("Script", "06_modules_WGCNA/01b_module_supermodule_GO_heatmaps.R")
  dry_run_line("Dataset", DATASET)
  dry_run_line("Ontologies", paste(ONTOLOGIES, collapse = ", "), "INFO")
  dry_run_line("Module GO enrichment", GO_FILE, if (file.exists(GO_FILE)) "PASS" else "FAIL")
  dry_run_line("Module label map", MODULE_MAP_FILE, if (file.exists(MODULE_MAP_FILE)) "PASS" else "FAIL")
  dry_run_line("Authoritative Stage 01 supermodule membership", MEMBERSHIP_FILE, if (file.exists(MEMBERSHIP_FILE)) "PASS" else "FAIL")
  dry_run_line("Output figures", PATHS$figures, "INFO")
  quit(status = if (all(file.exists(c(GO_FILE, MODULE_MAP_FILE, MEMBERSHIP_FILE)))) 0L else 1L, save = "no")
}

packages <- c("dplyr", "readr", "ggplot2", "scales")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))

read_required <- function(path, label) {
  x <- safe_read_csv(path)
  if (is.null(x)) stop("Missing or unreadable ", label, ": ", path, call. = FALSE)
  as.data.frame(x, stringsAsFactors = FALSE)
}

go <- read_required(GO_FILE, "module GO enrichment")
module_universe <- read_required(MODULE_MAP_FILE, "module label map")
membership <- read_required(MEMBERSHIP_FILE, "Stage 01 supermodule membership")
validate_wgcna_module_go_enrichment(go)
if (!"ModuleID" %in% names(module_universe)) stop("Module label map is missing ModuleID.", call. = FALSE)
if (length(setdiff(c("ModuleID", "SupermoduleID"), names(membership)))) stop("Supermodule membership is missing ModuleID and/or SupermoduleID.", call. = FALSE)

module_universe <- module_universe |>
  dplyr::mutate(ModuleID = as.character(.data$ModuleID)) |>
  dplyr::distinct(.data$ModuleID, .keep_all = TRUE)
module_label_col <- first_present_col(module_universe, c("ModuleLabel_Final", "final_label", "ModuleLabel_GO_BP"))
module_universe$ModuleDisplayLabel <- if (!is.na(module_label_col)) paste0(module_universe$ModuleID, " · ", as.character(module_universe[[module_label_col]])) else module_universe$ModuleID
module_universe$ModulePlotLabel <- sub("^WGCNA_", "", module_universe$ModuleID)
membership <- membership |>
  dplyr::mutate(ModuleID = as.character(.data$ModuleID), SupermoduleID = as.character(.data$SupermoduleID)) |>
  dplyr::filter(!is.na(.data$ModuleID), nzchar(.data$ModuleID), !is.na(.data$SupermoduleID), nzchar(.data$SupermoduleID)) |>
  dplyr::distinct(.data$ModuleID, .keep_all = TRUE)
if (nrow(membership) != nrow(module_universe) || !setequal(membership$ModuleID, module_universe$ModuleID)) {
  stop("Stage 01 supermodule membership must cover every module in module_name_map.csv exactly once.", call. = FALSE)
}
super_label_col <- first_present_col(membership, c("Supermodule_DisplayLabel", "Supermodule_FinalLabel", "Supermodule_LongLabel"))
membership$SupermoduleDisplayLabel <- if (!is.na(super_label_col)) as.character(membership[[super_label_col]]) else membership$SupermoduleID
blank_super_label <- is.na(membership$SupermoduleDisplayLabel) | !nzchar(membership$SupermoduleDisplayLabel)
membership$SupermoduleDisplayLabel[blank_super_label] <- membership$SupermoduleID[blank_super_label]
module_universe <- dplyr::left_join(
  module_universe,
  membership[c("ModuleID", "SupermoduleID", "SupermoduleDisplayLabel")],
  by = "ModuleID", relationship = "one-to-one"
) |>
  dplyr::arrange(.data$SupermoduleID, .data$ModuleID)

# Nature-style chart contract: read module biology in the detailed panel,
# then read recurring biology across supermodules in the companion panel.
NATURE_DOUBLE_COLUMN_IN <- 7.2
chart_contract <- data.frame(
  panel = c("module", "supermodule"),
  analytical_question = c(
    "Which representative GO terms annotate each co-regulated WGCNA module?",
    "Which GO themes recur across the member modules of each supermodule?"
  ),
  chart_family = "heatmap",
  displayed_measure = c(
    "capped -log10(BH FDR) from module ORA",
    "mean capped -log10(BH FDR) across member modules"
  ),
  visual_structure = c("modules grouped by data-driven supermodule", "one column per data-driven supermodule"),
  palette = "sequential purple, white = not FDR-significant",
  output_width_in = NATURE_DOUBLE_COLUMN_IN,
  stringsAsFactors = FALSE
)
chart_contract_file <- file.path(PATHS$tables, "WGCNA_GO_heatmap_chart_contract.csv")
write_csv_safe2(chart_contract, chart_contract_file)

plot_heatmap <- function(x, entity_col, entity_label_col, title, subtitle, output_stub, facet_col = NULL) {
  if (!nrow(x)) return(invisible(NULL))
  term_levels <- unique(x[order(x$selection_rank, x$Description), "Description"])
  entity_levels <- unique(x[order(x[[entity_col]]), entity_label_col])
  x$Description <- factor(x$Description, levels = rev(term_levels))
  x[[entity_label_col]] <- factor(x[[entity_label_col]], levels = entity_levels)
  if (!is.null(facet_col) && facet_col %in% names(x)) {
    x[[facet_col]] <- factor(x[[facet_col]], levels = unique(x[order(x[[facet_col]]), facet_col]))
  }
  width <- NATURE_DOUBLE_COLUMN_IN
  height <- max(4.2, 1.55 + 0.23 * length(term_levels))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data[[entity_label_col]], y = .data$Description, fill = .data$EnrichmentScore)) +
    ggplot2::geom_tile(colour = "#FFFFFF", linewidth = 0.3) +
    ggplot2::scale_fill_gradientn(
      colours = c("#FFFFFF", "#EFEDF5", "#BDBDD9", "#807DBA", "#54278F"),
      limits = c(0, SCORE_CAP), oob = scales::squish,
      name = expression("-log"[10] * "(BH FDR)")
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 1, size = 6.5),
      axis.text.y = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7, colour = "#444444"),
      strip.text.x = ggplot2::element_text(face = "bold", size = 6.5),
      strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = NA),
      legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 6.5),
      plot.margin = ggplot2::margin(4, 5, 4, 4, unit = "pt")
    )
  if (!is.null(facet_col) && facet_col %in% names(x)) {
    p <- p + ggplot2::facet_grid(stats::as.formula(paste("~", facet_col)), scales = "free_x", space = "free_x")
  }
  for (extension in c("svg", "pdf", "png")) {
    ggplot2::ggsave(file.path(PATHS$figures, paste0(output_stub, ".", extension)), p, width = width, height = height, units = "in", limitsize = FALSE, dpi = 300)
  }
  invisible(p)
}

status_rows <- list()
outputs <- chart_contract_file
for (ontology in ONTOLOGIES) {
  full_module_result <- make_full_module_go_term_matrix(go, module_universe, ontology = ontology, fdr_cutoff = FDR_CUTOFF, score_cap = SCORE_CAP)
  module_result <- make_module_go_heatmap_data(go, module_universe, ontology = ontology, fdr_cutoff = FDR_CUTOFF, score_cap = SCORE_CAP, top_terms_per_module = TOP_TERMS_PER_MODULE, max_terms = MAX_TERMS)
  module_matrix <- module_result$matrix
  module_selection <- module_result$selection
  module_matrix_file <- file.path(PATHS$tables, paste0("WGCNA_module_GO_heatmap_source_", ontology, ".csv"))
  module_selection_file <- file.path(PATHS$tables, paste0("WGCNA_module_GO_representative_terms_", ontology, ".csv"))
  write_csv_safe2(module_matrix, module_matrix_file)
  write_csv_safe2(module_selection, module_selection_file)
  outputs <- c(outputs, module_matrix_file, module_selection_file)
  if (nrow(module_matrix)) {
    module_matrix <- dplyr::left_join(
      module_matrix,
      module_universe[c("ModuleID", "ModuleDisplayLabel", "ModulePlotLabel", "SupermoduleID", "SupermoduleDisplayLabel")],
      by = "ModuleID"
    )
    plot_heatmap(module_matrix, "ModuleID", "ModulePlotLabel",
      paste0(DATASET, " WGCNA modules: representative GO-", ontology, " terms"),
      paste0("Modules are grouped by data-driven supermodule. Cell colour: capped -log10(BH FDR); white = not FDR <= ", FDR_CUTOFF, "."),
      paste0("WGCNA_module_GO_heatmap_", ontology), facet_col = "SupermoduleID")
  }

  super_result <- make_supermodule_go_heatmap_data(full_module_result$matrix, membership, ontology = ontology, fdr_cutoff = FDR_CUTOFF, score_cap = SCORE_CAP, top_terms_per_supermodule = TOP_TERMS_PER_SUPERMODULE, max_terms = MAX_TERMS)
  super_matrix <- super_result$matrix
  super_selection <- super_result$selection
  super_matrix_file <- file.path(PATHS$tables, paste0("WGCNA_supermodule_GO_heatmap_source_", ontology, ".csv"))
  super_selection_file <- file.path(PATHS$tables, paste0("WGCNA_supermodule_GO_representative_terms_", ontology, ".csv"))
  write_csv_safe2(super_matrix, super_matrix_file)
  write_csv_safe2(super_selection, super_selection_file)
  outputs <- c(outputs, super_matrix_file, super_selection_file)
  if (nrow(super_matrix)) {
    plot_heatmap(super_matrix, "SupermoduleID", "SupermoduleID",
      paste0(DATASET, " WGCNA supermodules: representative GO-", ontology, " terms"),
      paste0("Supermodule IDs map to biological labels in the source table. Cell colour: recurrence-weighted mean member-module -log10(BH FDR); white = no member FDR <= ", FDR_CUTOFF, ". Not a pooled ORA test."),
      paste0("WGCNA_supermodule_GO_heatmap_", ontology))
  }
  status_rows[[length(status_rows) + 1L]] <- data.frame(
    dataset = DATASET, ontology = ontology, module_status = module_result$status, supermodule_status = super_result$status,
    n_module_terms = nrow(module_selection), n_supermodule_terms = nrow(super_selection),
    fdr_cutoff = FDR_CUTOFF, score_cap = SCORE_CAP,
    supermodule_measure = "recurrence_weighted_member_module_evidence_not_pooled_ORA; terms_selected_from_full_module_GO_table",
    membership_source = MEMBERSHIP_FILE, stringsAsFactors = FALSE
  )
}
status <- do.call(rbind, status_rows)
status_file <- file.path(PATHS$tables, "WGCNA_GO_heatmap_status.csv")
write_csv_safe2(status, status_file)
outputs <- c(outputs, status_file)
write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"), inputs = c(GO_FILE, MODULE_MAP_FILE, MEMBERSHIP_FILE), outputs = outputs,
  parameters = list(dataset = DATASET, ontologies = ONTOLOGIES, fdr_cutoff = FDR_CUTOFF, score_cap = SCORE_CAP,
                    top_terms_per_module = TOP_TERMS_PER_MODULE, top_terms_per_supermodule = TOP_TERMS_PER_SUPERMODULE, max_terms = MAX_TERMS),
  notes = c("Module cells are capped -log10 BH-FDR from all-protein module ORA.", "Supermodule terms are selected from the full module GO table before member-module evidence is summarized; they are not pooled ORA tests.")
)
message("Wrote representative WGCNA GO heatmaps to: ", PATHS$figures)
