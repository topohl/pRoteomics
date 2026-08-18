#!/usr/bin/env Rscript
# ================================================================
# Script: 06_modules_WGCNA/01b_module_supermodule_GO_heatmaps.R
# Stage: modules_downstream
# Scope: dataset_specific plus focused all-datasets compositor.
# Consumes: required Stage 01 module GO enrichment and supermodule membership.
# Produces: representative module and supermodule GO comparison heatmaps.
# Dataset behavior: runs broad/focused outputs per dataset; --dataset all
#                   combines existing focused sources without rebuilding WGCNA.
# Notes: Uses all-module ORA results. Supermodule cells summarize member-module
#        evidence and are not pooled enrichment tests.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "plotting_nature.R"))
source(repo_path("R", "wgcna_downstream_utils.R"))
source(repo_path("R", "wgcna_go_comparison_utils.R"))
source(repo_path("R", "validation_utils.R"))

run <- wgcna_cli(allow_all = TRUE)
DATASET <- run$dataset
IS_COMBINED <- identical(DATASET, "all")
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
FOCUSED_TERMS_PER_SUPERMODULE <- suppressWarnings(as.integer(arg_value("--focused-terms-per-supermodule", "3")))
FOCUSED_HIERARCHY_GENE_JACCARD <- 0.50
FOCUSED_NEAR_IDENTICAL_GENE_JACCARD <- 0.80
MAX_TERMS <- suppressWarnings(as.integer(arg_value("--max-terms", "24")))
if (any(!is.finite(c(TOP_TERMS_PER_MODULE, TOP_TERMS_PER_SUPERMODULE, FOCUSED_TERMS_PER_SUPERMODULE, MAX_TERMS))) ||
    any(c(TOP_TERMS_PER_MODULE, TOP_TERMS_PER_SUPERMODULE, FOCUSED_TERMS_PER_SUPERMODULE, MAX_TERMS) < 1L)) {
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
FOCUSED_DATASETS <- wgcna_focused_dataset_order()
focused_dataset_table_file <- function(dataset, stem, ontology) {
  path_results(
    "tables", "06_modules_WGCNA", "01b_module_supermodule_GO_heatmaps", dataset,
    paste0(stem, "_", ontology, ".csv")
  )
}

if (run$dry_run) {
  invisible(lapply(unlist(PATHS), dir_create))
  dry_run_line("Script", "06_modules_WGCNA/01b_module_supermodule_GO_heatmaps.R")
  dry_run_line("Dataset", DATASET)
  dry_run_line("Ontologies", paste(ONTOLOGIES, collapse = ", "), "INFO")
  dry_run_line("Focused terms per supermodule", FOCUSED_TERMS_PER_SUPERMODULE, "INFO")
  if (IS_COMBINED) {
    combined_inputs <- unlist(lapply(ONTOLOGIES, function(ontology) {
      unlist(lapply(FOCUSED_DATASETS, function(dataset) c(
        focused_dataset_table_file(dataset, "WGCNA_supermodule_GO_focused_source", ontology),
        focused_dataset_table_file(dataset, "WGCNA_supermodule_GO_focused_selected_terms", ontology)
      )))
    }))
    for (input in combined_inputs) {
      dry_run_line("Dataset-focused source", input, if (file.exists(input)) "PASS" else "FAIL")
    }
    dry_run_line("Combined output figures", PATHS$figures, "INFO")
    quit(status = if (all(file.exists(combined_inputs))) 0L else 1L, save = "no")
  }
  dry_run_line("Module GO enrichment", GO_FILE, if (file.exists(GO_FILE)) "PASS" else "FAIL")
  dry_run_line("Module label map", MODULE_MAP_FILE, if (file.exists(MODULE_MAP_FILE)) "PASS" else "FAIL")
  dry_run_line("Authoritative Stage 01 supermodule membership", MEMBERSHIP_FILE, if (file.exists(MEMBERSHIP_FILE)) "PASS" else "FAIL")
  dry_run_line("Output figures", PATHS$figures, "INFO")
  quit(status = if (all(file.exists(c(GO_FILE, MODULE_MAP_FILE, MEMBERSHIP_FILE)))) 0L else 1L, save = "no")
}

packages <- c("dplyr", "readr", "ggplot2", "scales", "svglite", if (IS_COMBINED) "patchwork")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))

read_required <- function(path, label) {
  x <- safe_read_csv(path)
  if (is.null(x)) stop("Missing or unreadable ", label, ": ", path, call. = FALSE)
  as.data.frame(x, stringsAsFactors = FALSE)
}

NATURE_DOUBLE_COLUMN_IN <- mm_to_in(nature_dimensions_mm()[["double_column"]])
COMBINED_FOCUSED_WIDTH_MM <- nature_dimensions_mm()[["double_column"]]
MANUSCRIPT_TEXT <- nature_manuscript_text_sizes_pt()

save_go_figure <- function(plot, filename, width, height, units) {
  extension <- tolower(tools::file_ext(filename))
  device <- switch(
    extension,
    svg = svglite::svglite,
    pdf = grDevices::cairo_pdf,
    png = if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png",
    stop("Unsupported GO figure extension: ", extension, call. = FALSE)
  )
  ggplot2::ggsave(
    filename, plot, width = width, height = height, units = units,
    device = device, limitsize = FALSE, dpi = 300, bg = "white"
  )
}

build_combined_focused_panel <- function(x, colour_limits, size_limits) {
  term_rows <- x[!duplicated(x$TermID), c(
    "TermID", "Description", "row_order", "cross_dataset_recurrence_exact_TermID"
  ), drop = FALSE]
  term_rows <- term_rows[order(term_rows$row_order, term_rows$TermID), , drop = FALSE]
  term_rows$GOPlotLabel <- term_rows$Description
  term_labels <- stats::setNames(
    vapply(term_rows$GOPlotLabel, function(one) paste(strwrap(one, width = 64L), collapse = "\n"), character(1)),
    term_rows$TermID
  )
  x$FocusedTermKey <- factor(x$TermID, levels = rev(term_rows$TermID))
  supermodule_order <- wgcna_canonical_supermodule_order(x$SupermoduleID)
  supermodule_labels <- unique(x[c("SupermoduleID", "SupermodulePlotLabel")])
  supermodule_labels <- supermodule_labels[match(supermodule_order, supermodule_labels$SupermoduleID), , drop = FALSE]
  x$SupermodulePlotLabel <- factor(x$SupermodulePlotLabel, levels = supermodule_labels$SupermodulePlotLabel)
  supported <- x[x$display_supported_dot, , drop = FALSE]

  ggplot2::ggplot(x, ggplot2::aes(x = .data$SupermodulePlotLabel, y = .data$FocusedTermKey)) +
    ggplot2::geom_point(
      data = supported,
      ggplot2::aes(
        colour = .data$mean_member_module_enrichment_score,
        size = .data$fraction_member_modules_FDR_significant
      ),
      shape = 16
    ) +
    ggplot2::scale_colour_gradientn(
      colours = nature_palette("support"),
      limits = colour_limits, oob = scales::squish,
      name = "Mean member-module -log10(BH FDR)",
      guide = ggplot2::guide_colourbar(
        title.position = "top", title.hjust = 0.5,
        barwidth = grid::unit(1.4, "in"), barheight = grid::unit(0.10, "in"),
        order = 1
      )
    ) +
    ggplot2::scale_size_continuous(
      limits = size_limits, range = c(1.5, 4.2),
      breaks = c(0.25, 0.5, 0.75, 1), labels = scales::label_number(accuracy = 0.01),
      name = "Fraction of member modules with BH FDR <= 0.05",
      guide = ggplot2::guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, order = 2)
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(labels = term_labels, drop = FALSE) +
    ggplot2::labs(title = unique(x$dataset_display_label), x = NULL, y = NULL) +
    theme_nature_base(base_size = MANUSCRIPT_TEXT[["normal"]], base_family = "Arial") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = MANUSCRIPT_TEXT[["dense"]], colour = "#333333"),
      axis.text.y = ggplot2::element_text(size = MANUSCRIPT_TEXT[["dense"]], colour = "#222222", lineheight = 0.94),
      plot.title = ggplot2::element_text(face = "bold", size = MANUSCRIPT_TEXT[["title"]], margin = ggplot2::margin(b = 2, unit = "pt")),
      legend.title = ggplot2::element_text(size = MANUSCRIPT_TEXT[["normal"]]),
      legend.text = ggplot2::element_text(size = MANUSCRIPT_TEXT[["normal"]]),
      plot.margin = ggplot2::margin(2, 3, 2, 2, unit = "pt")
    )
}

if (IS_COMBINED) {
  combined_outputs <- character(0)
  combined_inputs <- character(0)
  combined_status_rows <- list()
  combined_contract_rows <- list()
  for (ontology in ONTOLOGIES) {
    source_files <- stats::setNames(vapply(FOCUSED_DATASETS, function(dataset) {
      focused_dataset_table_file(dataset, "WGCNA_supermodule_GO_focused_source", ontology)
    }, character(1)), FOCUSED_DATASETS)
    selection_files <- stats::setNames(vapply(FOCUSED_DATASETS, function(dataset) {
      focused_dataset_table_file(dataset, "WGCNA_supermodule_GO_focused_selected_terms", ontology)
    }, character(1)), FOCUSED_DATASETS)
    source_tables <- lapply(names(source_files), function(dataset) {
      read_required(source_files[[dataset]], paste0(dataset, " focused source"))
    })
    names(source_tables) <- names(source_files)
    selection_tables <- lapply(names(selection_files), function(dataset) {
      read_required(selection_files[[dataset]], paste0(dataset, " focused selected terms"))
    })
    names(selection_tables) <- names(selection_files)
    expected_selection_token <- paste0("max_", FOCUSED_TERMS_PER_SUPERMODULE, "_per_supermodule")
    if (any(!vapply(selection_tables, function(x) {
      nrow(x) && all(grepl(expected_selection_token, x$selection_basis, fixed = TRUE))
    }, logical(1)))) {
      stop(
        "Combined focused inputs were not all generated with --focused-terms-per-supermodule ",
        FOCUSED_TERMS_PER_SUPERMODULE, ". Regenerate dataset-specific focused sources first.",
        call. = FALSE
      )
    }
    combined <- combine_focused_supermodule_go_sources(source_tables, selection_tables)
    if (!isTRUE(all.equal(combined$colour_limits, c(0, SCORE_CAP)))) {
      stop("Combined focused source ScoreCap does not match --score-cap; regenerate dataset-specific sources consistently.", call. = FALSE)
    }
    if (!isTRUE(all.equal(combined$fdr_cutoff, FDR_CUTOFF))) {
      stop("Combined focused source FDRCutoff does not match --fdr; regenerate dataset-specific sources consistently.", call. = FALSE)
    }

    combined_source_file <- file.path(PATHS$tables, paste0("WGCNA_supermodule_GO_focused_source_all_datasets_", ontology, ".csv"))
    combined_selection_file <- file.path(PATHS$tables, paste0("WGCNA_supermodule_GO_focused_selected_terms_all_datasets_", ontology, ".csv"))
    write_csv_safe2(combined$matrix, combined_source_file)
    write_csv_safe2(combined$selection, combined_selection_file)

    panels <- lapply(FOCUSED_DATASETS, function(dataset) {
      build_combined_focused_panel(
        combined$matrix[combined$matrix$dataset == dataset, , drop = FALSE],
        colour_limits = combined$colour_limits,
        size_limits = combined$size_limits
      )
    })
    n_terms <- vapply(FOCUSED_DATASETS, function(dataset) {
      length(unique(combined$matrix$TermID[combined$matrix$dataset == dataset]))
    }, integer(1))
    panel_heights <- 0.68 + 0.12 * n_terms
    combined_height_mm <- min(
      nature_dimensions_mm()[["maximum_height"]],
      max(105, 48 + 2.34 * sum(n_terms))
    )
    combined_plot <- patchwork::wrap_plots(
      panels, ncol = 1, heights = panel_heights, guides = "collect"
    ) +
      patchwork::plot_annotation(
        title = paste0("Focused WGCNA supermodule GO-", ontology, " comparison"),
        subtitle = "FDR-supported member-module GO terms after conservative redundancy pruning; shared evidence scales.",
        caption = paste0(
          "No dot = no member module with BH FDR <= ", FDR_CUTOFF,
          ". Member-module summaries; not pooled supermodule enrichment tests.  \u2020 Singleton SM; GO summary reflects its constituent module."
        ),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = MANUSCRIPT_TEXT[["title"]]),
          plot.subtitle = ggplot2::element_text(size = MANUSCRIPT_TEXT[["normal"]], colour = "#444444"),
          plot.caption = ggplot2::element_text(size = MANUSCRIPT_TEXT[["dense"]], colour = "#555555", hjust = 0)
        )
      ) &
      ggplot2::theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "center"
      )

    output_stub <- paste0("WGCNA_supermodule_GO_focused_comparison_all_datasets_", ontology)
    figure_files <- file.path(PATHS$figures, paste0(output_stub, ".", c("svg", "pdf", "png")))
    for (figure_file in figure_files) {
      save_go_figure(
        combined_plot, figure_file,
        width = COMBINED_FOCUSED_WIDTH_MM, height = combined_height_mm,
        units = "mm"
      )
    }

    combined_inputs <- c(combined_inputs, unname(source_files), unname(selection_files))
    combined_outputs <- c(combined_outputs, combined_source_file, combined_selection_file, figure_files)
    recurring_terms <- unique(combined$selection[
      combined$selection$n_datasets_selected >= 2L,
      c("TermID", "Description", "n_datasets_selected", "datasets_selected"),
      drop = FALSE
    ])
    combined_status_rows[[length(combined_status_rows) + 1L]] <- data.frame(
      dataset = "all", ontology = ontology, combined_focused_panel_generated = length(figure_files) == 3L,
      focused_terms_per_supermodule_requested = FOCUSED_TERMS_PER_SUPERMODULE,
      n_datasets = length(FOCUSED_DATASETS), n_source_rows = nrow(combined$matrix),
      n_selected_term_rows = nrow(combined$selection),
      n_exact_terms_selected_in_multiple_datasets = nrow(recurring_terms),
      combined_colour_measure = "mean_member_module_enrichment_score",
      combined_colour_scale = paste(combined$colour_limits, collapse = ";"),
      combined_size_measure = "fraction_member_modules_FDR_significant",
      combined_size_scale = paste(combined$size_limits, collapse = ";"),
      redundancy_modes = paste(sort(unique(combined$selection$redundancy_mode)), collapse = ";"),
      hierarchy_gene_jaccard_threshold = FOCUSED_HIERARCHY_GENE_JACCARD,
      near_identical_gene_jaccard_threshold = FOCUSED_NEAR_IDENTICAL_GENE_JACCARD,
      cross_dataset_recurrence = "descriptive exact TermID selection only; no cross-dataset inference",
      combined_figure_width_in = mm_to_in(COMBINED_FOCUSED_WIDTH_MM),
      combined_figure_height_in = mm_to_in(combined_height_mm),
      stringsAsFactors = FALSE
    )
    combined_contract_rows[[length(combined_contract_rows) + 1L]] <- data.frame(
      panel = "supermodule_focused_all_datasets",
      ontology = ontology,
      analytical_question = "Which representative biological processes distinguish WGCNA supermodule annotation profiles within and across compartments?",
      chart_family = "vertically stacked coordinated dot matrices",
      displayed_measure = "shared colour = mean member-module capped -log10(BH FDR); shared size = fraction of member modules with BH FDR <= 0.05",
      visual_structure = "dataset-specific term rows and canonical SM columns; one shared legend; no global GO-term union",
      palette = "shared sequential purple",
      selection_contract = paste0(
        "up to ", FOCUSED_TERMS_PER_SUPERMODULE,
        " FDR-supported terms per SM after sequential GO hierarchy/gene-overlap pruning; hierarchy+Jaccard >= ",
        FOCUSED_HIERARCHY_GENE_JACCARD, " or Jaccard >= ", FOCUSED_NEAR_IDENTICAL_GENE_JACCARD
      ),
      inference_scope = "descriptive member-module summaries and exact TermID recurrence; no pooled or cross-dataset inference",
      output_width_in = mm_to_in(COMBINED_FOCUSED_WIDTH_MM),
      output_height_in = mm_to_in(combined_height_mm),
      stringsAsFactors = FALSE
    )
  }
  combined_status_file <- file.path(PATHS$tables, "WGCNA_GO_heatmap_status.csv")
  combined_contract_file <- file.path(PATHS$tables, "WGCNA_GO_heatmap_chart_contract.csv")
  write_csv_safe2(do.call(rbind, combined_status_rows), combined_status_file)
  write_csv_safe2(do.call(rbind, combined_contract_rows), combined_contract_file)
  combined_outputs <- c(combined_outputs, combined_status_file, combined_contract_file)
  write_run_manifest(
    file.path(PATHS$logs, "run_manifest.yml"),
    inputs = unique(combined_inputs), outputs = combined_outputs,
    parameters = list(
      dataset = "all", ontologies = ONTOLOGIES, fdr_cutoff = FDR_CUTOFF,
      score_cap = SCORE_CAP, focused_terms_per_supermodule = FOCUSED_TERMS_PER_SUPERMODULE,
      focused_hierarchy_gene_jaccard_threshold = FOCUSED_HIERARCHY_GENE_JACCARD,
      focused_near_identical_gene_jaccard_threshold = FOCUSED_NEAR_IDENTICAL_GENE_JACCARD,
      source_datasets = FOCUSED_DATASETS
    ),
    notes = c(
      "Combines existing dataset-specific focused sources as coordinated panels; does not rebuild WGCNA or GO enrichment.",
      "Colour and size use shared unnormalized member-module summary scales across datasets.",
      "Focused representatives are display-only redundancy-pruned terms; dataset-specific audit tables retain every supported candidate.",
      "Cross-dataset recurrence is descriptive exact TermID identity only; no p-value, FDR, meta-analysis, or convergence claim is computed."
    )
  )
  message("Wrote combined focused WGCNA GO comparisons to: ", PATHS$figures)
  quit(status = 0L, save = "no")
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
chart_contract <- data.frame(
  panel = c("module", "supermodule", "supermodule_focused_comparison"),
  analytical_question = c(
    "Which representative GO terms annotate each co-regulated WGCNA module?",
    "Which GO themes recur across the member modules of each supermodule?",
    "Which representative biological processes distinguish the GO annotation profiles of WGCNA supermodules?"
  ),
  chart_family = c("heatmap", "heatmap", "dot matrix"),
  displayed_measure = c(
    "capped -log10(BH FDR) from module ORA",
    "mean capped -log10(BH FDR) across member modules",
    "colour = mean member-module capped -log10(BH FDR); size = fraction of member modules with BH FDR <= 0.05"
  ),
  visual_structure = c(
    "modules grouped by data-driven supermodule",
    "one column per data-driven supermodule",
    "union of FDR-supported representative terms compared across every data-driven supermodule"
  ),
  palette = "sequential purple, white = not FDR-significant",
  selection_contract = c(
    "broad representative selection; no focused redundancy pruning",
    "broad representative selection; no focused redundancy pruning",
    paste0(
      "up to ", FOCUSED_TERMS_PER_SUPERMODULE,
      " FDR-supported terms per SM after sequential GO hierarchy/gene-overlap pruning; hierarchy+Jaccard >= ",
      FOCUSED_HIERARCHY_GENE_JACCARD, " or Jaccard >= ", FOCUSED_NEAR_IDENTICAL_GENE_JACCARD
    )
  ),
  output_width_in = NATURE_DOUBLE_COLUMN_IN,
  inference_scope = c(
    "module ORA",
    "member-module evidence; not pooled supermodule ORA",
    "member-module evidence; not pooled supermodule ORA"
  ),
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
    save_go_figure(
      p, file.path(PATHS$figures, paste0(output_stub, ".", extension)),
      width = width, height = height, units = "in"
    )
  }
  invisible(p)
}

wrap_go_description <- function(x, width = 43L) {
  vapply(as.character(x), function(one) paste(strwrap(one, width = width), collapse = "\n"), character(1))
}

plot_focused_supermodule_comparison <- function(x, ontology, output_stub) {
  if (!nrow(x)) return(invisible(NULL))
  supermodule_order <- wgcna_canonical_supermodule_order(x$SupermoduleID)
  term_rows <- x[!duplicated(x$TermID), c("TermID", "Description", "row_order"), drop = FALSE]
  term_rows <- term_rows[order(term_rows$row_order), , drop = FALSE]
  x$FocusedTermKey <- factor(x$TermID, levels = rev(term_rows$TermID))
  supermodule_labels <- unique(x[c("SupermoduleID", "SupermodulePlotLabel")])
  supermodule_labels <- supermodule_labels[match(supermodule_order, supermodule_labels$SupermoduleID), , drop = FALSE]
  x$SupermodulePlotLabel <- factor(x$SupermodulePlotLabel, levels = supermodule_labels$SupermodulePlotLabel)
  term_labels <- stats::setNames(wrap_go_description(term_rows$Description), term_rows$TermID)
  supported <- x[x$display_supported_dot, , drop = FALSE]
  singleton_note <- if (any(x$is_singleton_supermodule)) "  \u2020 Singleton SM; GO summary reflects its constituent module." else ""
  is_figure3_main <- identical(DATASET, "neuron_neuropil") && identical(ontology, "BP")
  wrapped_term_lines <- sum(vapply(term_labels, function(one) length(strsplit(one, "\n", fixed = TRUE)[[1]]), integer(1)))
  height <- if (is_figure3_main) max(94, 28 + 4.8 * wrapped_term_lines) else max(4.2, 1.9 + 0.31 * nrow(term_rows))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data$SupermodulePlotLabel, y = .data$FocusedTermKey)) +
    ggplot2::geom_point(
      data = supported,
      ggplot2::aes(
        colour = .data$mean_member_module_enrichment_score,
        size = .data$fraction_member_modules_FDR_significant
      ),
      shape = 16
    ) +
    ggplot2::scale_colour_gradientn(
      colours = nature_palette("support"),
      limits = c(0, SCORE_CAP), oob = scales::squish,
      name = if (is_figure3_main) "Mean module\n-log10(BH FDR)" else "Mean member-module\n-log10(BH FDR)"
    ) +
    ggplot2::scale_size_continuous(
      limits = c(0, 1), range = if (is_figure3_main) c(2.0, 4.8) else c(1.4, 5.0),
      breaks = c(0.25, 0.5, 0.75, 1), labels = scales::label_number(accuracy = 0.01),
      name = if (is_figure3_main) "Supported module fraction" else "Fraction of member modules\nwith BH FDR <= 0.05"
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(labels = term_labels) +
    ggplot2::labs(
      title = if (is_figure3_main) NULL else paste0(DATASET, " WGCNA supermodules: focused GO-", ontology, " comparison"),
      subtitle = if (is_figure3_main) NULL else "Up to three FDR-supported member-module processes after conservative redundancy pruning.",
      caption = if (is_figure3_main) NULL else paste0(
        "No dot = no member module with BH FDR <= ", FDR_CUTOFF,
        ". Member-module summary; not a pooled supermodule enrichment test.", singleton_note
      ),
      x = NULL, y = NULL
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_colourbar(
        order = 1, title.position = "top",
        barwidth = grid::unit(if (is_figure3_main) 20 else 25, "mm"),
        barheight = grid::unit(2.1, "mm")
      ),
      size = ggplot2::guide_legend(order = 2, title.position = "top", nrow = 1)
    ) +
    (if (is_figure3_main) {
      theme_nature_manuscript_panel(base_size = MANUSCRIPT_TEXT[["normal"]], base_family = "Arial", axes = FALSE, publication_legible = TRUE)
    } else {
      theme_nature_base(base_size = 7, base_family = "Arial")
    }) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = if (is_figure3_main) MANUSCRIPT_TEXT[["normal"]] else 6.5),
      axis.text.y = ggplot2::element_text(size = if (is_figure3_main) MANUSCRIPT_TEXT[["dense"]] else 6.5, colour = "#222222", lineheight = 1.02),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 7, colour = "#444444"),
      plot.caption = ggplot2::element_text(size = if (is_figure3_main) MANUSCRIPT_TEXT[["dense"]] else 6.5, colour = "#555555", hjust = 0),
      legend.position = if (is_figure3_main) "bottom" else "right",
      legend.box = if (is_figure3_main) "vertical" else NULL,
      legend.title = ggplot2::element_text(size = if (is_figure3_main) MANUSCRIPT_TEXT[["normal"]] else 6.5),
      legend.text = ggplot2::element_text(size = if (is_figure3_main) MANUSCRIPT_TEXT[["normal"]] else 6.5),
      plot.margin = ggplot2::margin(2, 2, 2, 2, unit = "pt")
    )
  figure_files <- file.path(PATHS$figures, paste0(output_stub, ".", c("svg", "pdf", "png")))
  for (figure_file in figure_files) {
    save_go_figure(
      p, figure_file,
      width = if (is_figure3_main) 110 else NATURE_DOUBLE_COLUMN_IN,
      height = height,
      units = if (is_figure3_main) "mm" else "in"
    )
  }
  invisible(figure_files)
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

  focused_full_module_result <- make_full_module_go_term_matrix(
    go, module_universe, ontology = ontology, fdr_cutoff = FDR_CUTOFF,
    score_cap = SCORE_CAP, retain_gene_membership = TRUE
  )
  hierarchy <- wgcna_go_hierarchy(unique(focused_full_module_result$matrix$TermID), ontology = ontology)
  focused_result <- make_focused_supermodule_go_comparison(
    focused_full_module_result$matrix, membership, ontology = ontology,
    fdr_cutoff = FDR_CUTOFF, score_cap = SCORE_CAP,
    focused_terms_per_supermodule = FOCUSED_TERMS_PER_SUPERMODULE,
    go_hierarchy = hierarchy,
    hierarchy_gene_jaccard_threshold = FOCUSED_HIERARCHY_GENE_JACCARD,
    near_identical_gene_jaccard_threshold = FOCUSED_NEAR_IDENTICAL_GENE_JACCARD
  )
  focused_matrix <- focused_result$matrix
  focused_selection <- focused_result$selection
  focused_audit <- focused_result$audit
  if (nrow(focused_matrix)) focused_matrix$dataset <- DATASET else focused_matrix$dataset <- character(0)
  if (nrow(focused_selection)) focused_selection$dataset <- DATASET else focused_selection$dataset <- character(0)
  if (nrow(focused_audit)) focused_audit$dataset <- DATASET else focused_audit$dataset <- character(0)
  focused_matrix <- focused_matrix[c("dataset", setdiff(names(focused_matrix), "dataset"))]
  focused_selection <- focused_selection[c("dataset", setdiff(names(focused_selection), "dataset"))]
  focused_audit <- focused_audit[c("dataset", setdiff(names(focused_audit), "dataset"))]
  focused_matrix_file <- file.path(PATHS$tables, paste0("WGCNA_supermodule_GO_focused_source_", ontology, ".csv"))
  focused_selection_file <- file.path(PATHS$tables, paste0("WGCNA_supermodule_GO_focused_selected_terms_", ontology, ".csv"))
  focused_audit_file <- file.path(PATHS$tables, paste0("WGCNA_supermodule_GO_focused_selection_audit_", ontology, ".csv"))
  write_csv_safe2(focused_matrix, focused_matrix_file)
  write_csv_safe2(focused_selection, focused_selection_file)
  write_csv_safe2(focused_audit, focused_audit_file)
  outputs <- c(outputs, focused_matrix_file, focused_selection_file, focused_audit_file)
  focused_figure_files <- character(0)
  if (nrow(focused_matrix)) {
    focused_figure_files <- plot_focused_supermodule_comparison(
      focused_matrix, ontology,
      paste0("WGCNA_supermodule_GO_focused_comparison_", ontology)
    )
    outputs <- c(outputs, focused_figure_files)
  }
  supermodule_ids <- wgcna_canonical_supermodule_order(membership$SupermoduleID)
  focused_selectors <- unique(as.character(focused_selection$selected_for_supermodule))
  focused_selectors <- focused_selectors[!is.na(focused_selectors) & nzchar(focused_selectors)]
  selection_status_counts <- table(factor(
    focused_audit$selection_status,
    levels = c(
      "selected", "skipped_redundant_GO_ancestor_descendant",
      "skipped_near_identical_gene_set", "not_reached_after_three_selected"
    )
  ))
  status_rows[[length(status_rows) + 1L]] <- data.frame(
    dataset = DATASET, ontology = ontology, module_status = module_result$status, supermodule_status = super_result$status,
    n_module_terms = nrow(module_selection), n_supermodule_terms = nrow(super_selection),
    fdr_cutoff = FDR_CUTOFF, score_cap = SCORE_CAP,
    supermodule_measure = "recurrence_weighted_member_module_evidence_not_pooled_ORA; terms_selected_from_full_module_GO_table",
    membership_source = MEMBERSHIP_FILE,
    focused_status = focused_result$status,
    focused_panel_generated = length(focused_figure_files) == 3L,
    focused_terms_per_supermodule_requested = FOCUSED_TERMS_PER_SUPERMODULE,
    n_unique_focused_terms = length(unique(focused_selection$TermID)),
    n_supermodules_total = length(supermodule_ids),
    n_supermodules_with_FDR_supported_representative_term = length(intersect(supermodule_ids, focused_selectors)),
    n_supermodules_without_FDR_supported_representative_term = length(setdiff(supermodule_ids, focused_selectors)),
    supermodules_without_FDR_supported_representative_term = paste(setdiff(supermodule_ids, focused_selectors), collapse = ";"),
    focused_measure = "member-module evidence; not pooled supermodule ORA",
    focused_redundancy_mode = focused_result$redundancy_mode,
    focused_hierarchy_source = hierarchy$source,
    focused_hierarchy_version = hierarchy$version,
    focused_hierarchy_status = if (isTRUE(hierarchy$available)) "available" else hierarchy$reason,
    focused_gene_membership_column = if ("geneID" %in% names(go)) "geneID" else NA_character_,
    focused_gene_membership_delimiter = if ("geneID" %in% names(go)) "/" else NA_character_,
    hierarchy_gene_jaccard_threshold = FOCUSED_HIERARCHY_GENE_JACCARD,
    near_identical_gene_jaccard_threshold = FOCUSED_NEAR_IDENTICAL_GENE_JACCARD,
    n_FDR_supported_candidates = nrow(focused_audit),
    n_focused_selected = unname(selection_status_counts[["selected"]]),
    n_skipped_redundant_GO_ancestor_descendant = unname(selection_status_counts[["skipped_redundant_GO_ancestor_descendant"]]),
    n_skipped_near_identical_gene_set = unname(selection_status_counts[["skipped_near_identical_gene_set"]]),
    n_excluded_after_three_selected = unname(selection_status_counts[["not_reached_after_three_selected"]]),
    stringsAsFactors = FALSE
  )
}
status <- do.call(rbind, status_rows)
status_file <- file.path(PATHS$tables, "WGCNA_GO_heatmap_status.csv")
write_csv_safe2(status, status_file)
outputs <- c(outputs, status_file)
write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"), inputs = c(GO_FILE, MODULE_MAP_FILE, MEMBERSHIP_FILE), outputs = outputs,
  parameters = list(dataset = DATASET, ontologies = ONTOLOGIES, fdr_cutoff = FDR_CUTOFF, score_cap = SCORE_CAP,
                    top_terms_per_module = TOP_TERMS_PER_MODULE, top_terms_per_supermodule = TOP_TERMS_PER_SUPERMODULE,
                    focused_terms_per_supermodule = FOCUSED_TERMS_PER_SUPERMODULE,
                    focused_hierarchy_gene_jaccard_threshold = FOCUSED_HIERARCHY_GENE_JACCARD,
                    focused_near_identical_gene_jaccard_threshold = FOCUSED_NEAR_IDENTICAL_GENE_JACCARD,
                    max_terms = MAX_TERMS),
  notes = c(
    "Module cells are capped -log10 BH-FDR from all-protein module ORA.",
    "Supermodule terms are selected from the full module GO table before member-module evidence is summarized; they are not pooled ORA tests.",
    "Focused representatives are selected sequentially from evidence-ranked FDR-supported candidates using exact GO hierarchy and significant-row gene overlap when available; every candidate is retained in the selection audit.",
    "Focused dot colour is mean member-module capped -log10 BH FDR; dot size is the fraction of member modules with BH FDR <= the configured cutoff."
  )
)
message("Wrote representative WGCNA GO heatmaps and focused comparisons to: ", PATHS$figures)
