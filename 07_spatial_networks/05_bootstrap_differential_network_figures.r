# ================================================================
# Script: 07_spatial_networks/05_bootstrap_differential_network_figures.r
# Stage: networks
# Scope: dataset_specific
# Consumes: required results/tables/07_spatial_networks/bootstrap_differential_network_stability/01_Tables/bootstrap_differential_edge_stability_summary.csv; results/tables/07_spatial_networks/bootstrap_differential_network_stability/01_Tables/bootstrap_differential_edge_values_long.csv; optional results/tables/07_spatial_networks/bootstrap_differential_network_stability/01_Tables/candidate_edge_differential_stability_summary.csv.
# Produces: results/figures/07_spatial_networks/bootstrap_differential_network_figures/.
# Dataset behavior: runs for neuron_neuropil,neuron_soma according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Figure layer after bootstrap differential network stability.
# ================================================================

# ================================================================
# Consumes:
#   - bootstrap differential network tables from canonical module output
# Produces:
#   - publication-ready figures and source tables in canonical folders
# File contract:
#   - docs/active_script_io_audit.tsv object 07_spatial_networks/05_bootstrap_differential_network_figures.r
# ================================================================
# Publication-style figures for bootstrap differential spatial networks
# ================================================================
# Purpose:
#   Generate publication-ready figures from outputs created by:
#   Analysis/bootstrap_differential_network_stability.r
#
# Main figure outputs:
#   1) Forest plot of stable differential edges with bootstrap intervals
#   2) Mean DeltaR / differential-frequency heatmap
#   3) Bootstrap DeltaR distributions for candidate edges
#   4) Minimal stable differential network panels
#
# Run after bootstrap_differential_network_stability.r has completed.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
MODULE_ID <- "07_spatial_networks"
SUBSTEP_ID <- "bootstrap_differential_network_figures"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)
NETWORK_DATASET <- current_dataset()
assert_dataset_capability(NETWORK_DATASET, "layer", analysis = "bootstrap differential spatial network figure export")

required_pkgs <- c(
  "dplyr", "tidyr", "stringr", "purrr", "tibble", "ggplot2",
  "readr", "ggraph", "igraph", "svglite", "scales", "forcats"
)
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), ". Install them explicitly before running this script.", call. = FALSE)
params <- list(
  bootstrap_dir = path_results("tables", "07_spatial_networks", "bootstrap_differential_network_stability"),

  stable_frequency_threshold = 0.70,
  top_n_edges_per_comparison = 8,

  candidate_edges = tibble::tribble(
    ~Source, ~Target,
    "CA1_slm", "CA2_slm",
    "CA1_slm", "CA3_sr",
    "CA2_slm", "CA3_sr",
    "CA1_sr",  "DG_mo",
    "CA1_so",  "DG_mo",
    "CA1_sr",  "CA3_sr",
    "CA2_slm", "DG_po",
    "CA1_so",  "CA1_sr"
  ),

  palette = list(
    lower_in_B = "#457B9D",
    higher_in_B = "#E63946",
    neutral = "#C6C3BB",
    text = "black"
  )
)

if (is_dry_run()) {
  dry_run_line("Script", "07_spatial_networks/05_bootstrap_differential_network_figures.r")
  dry_run_line("Bootstrap table directory", params$bootstrap_dir, if (dir.exists(params$bootstrap_dir)) "PASS" else "FAIL")
  dry_run_line("Output folders", paste(unlist(CANONICAL_PATHS), collapse = "; "))
  quit(status = if (dir.exists(params$bootstrap_dir)) 0 else 1, save = "no")
}
if (!dir.exists(params$bootstrap_dir)) stop("Bootstrap table directory not found: ", params$bootstrap_dir, call. = FALSE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

message2 <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)

safe_name <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

make_dirs <- function(base_dir) {
  dirs <- list(
    figures = CANONICAL_PATHS$figures,
    tables = CANONICAL_PATHS$source_data
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

theme_publication_boot <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "sans", colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}

normalize_edge <- function(df) {
  df %>%
    dplyr::mutate(
      Source0 = pmin(.data$Source, .data$Target),
      Target0 = pmax(.data$Source, .data$Target),
      Source = .data$Source0,
      Target = .data$Target0
    ) %>%
    dplyr::select(-"Source0", -"Target0")
}

add_edge_labels <- function(df) {
  df %>%
    dplyr::mutate(
      Edge = paste(.data$Source, .data$Target, sep = " - "),
      DirectionLabel = dplyr::case_when(
        .data$MeanDeltaR < 0 ~ "lower in B",
        .data$MeanDeltaR > 0 ~ "higher in B",
        TRUE ~ "no change"
      ),
      DirectionLabel = factor(.data$DirectionLabel, levels = c("lower in B", "higher in B", "no change"))
    )
}

load_bootstrap_outputs <- function(params) {
  tables_dir <- file.path(params$bootstrap_dir, "01_Tables")

  edge_summary_file <- file.path(tables_dir, "bootstrap_differential_edge_stability_summary.csv")
  boot_long_file <- file.path(tables_dir, "bootstrap_differential_edge_values_long.csv")
  candidate_file <- file.path(tables_dir, "candidate_edge_differential_stability_summary.csv")

  if (!file.exists(edge_summary_file)) stop("Missing edge summary file: ", edge_summary_file)
  if (!file.exists(boot_long_file)) stop("Missing bootstrap long file: ", boot_long_file)

  edge_summary <- readr::read_csv(edge_summary_file, show_col_types = FALSE) %>%
    normalize_edge() %>%
    add_edge_labels()

  boot_long <- readr::read_csv(boot_long_file, show_col_types = FALSE) %>%
    normalize_edge() %>%
    dplyr::mutate(Edge = paste(.data$Source, .data$Target, sep = " - "))

  candidate_summary <- if (file.exists(candidate_file)) {
    readr::read_csv(candidate_file, show_col_types = FALSE) %>%
      normalize_edge() %>%
      add_edge_labels()
  } else {
    cand <- params$candidate_edges %>% normalize_edge()
    edge_summary %>% dplyr::inner_join(cand, by = c("Source", "Target"))
  }

  list(edge_summary = edge_summary, boot_long = boot_long, candidate_summary = candidate_summary)
}

select_stable_edges <- function(edge_summary, params) {
  edge_summary %>%
    dplyr::filter(.data$DifferentialFrequency >= params$stable_frequency_threshold) %>%
    dplyr::group_by(.data$Comparison) %>%
    dplyr::slice_max(order_by = .data$MeanAbsDeltaR, n = params$top_n_edges_per_comparison, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Edge = forcats::fct_reorder(.data$Edge, .data$MeanDeltaR)
    )
}

plot_forest_stable_edges <- function(edge_summary, outfile, params) {
  dat <- select_stable_edges(edge_summary, params)
  if (nrow(dat) == 0) {
    warning("No stable edges available for forest plot.")
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$MeanDeltaR, y = .data$Edge)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$DeltaR_Q025, xmax = .data$DeltaR_Q975, colour = .data$DirectionLabel),
      height = 0,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$DirectionLabel, size = .data$DifferentialFrequency), shape = 21, colour = "black", stroke = 0.25) +
    ggplot2::facet_wrap(~ Comparison, scales = "free_y", ncol = 1) +
    ggplot2::scale_colour_manual(values = c("lower in B" = params$palette$lower_in_B, "higher in B" = params$palette$higher_in_B, "no change" = params$palette$neutral), name = NULL) +
    ggplot2::scale_fill_manual(values = c("lower in B" = params$palette$lower_in_B, "higher in B" = params$palette$higher_in_B, "no change" = params$palette$neutral), name = NULL) +
    ggplot2::scale_size_continuous(range = c(1.5, 3.2), name = "Bootstrap\nfrequency") +
    ggplot2::labs(x = expression(Delta * "Spearman " * rho), y = NULL, title = "Stable spatial proteomic rewiring") +
    theme_publication_boot(base_size = 7) +
    ggplot2::theme(legend.position = "right")

  ggplot2::ggsave(outfile, p, width = 95, height = 125, units = "mm", device = svglite::svglite)
  invisible(p)
}

plot_delta_frequency_heatmap <- function(edge_summary, outfile, params) {
  dat <- select_stable_edges(edge_summary, params) %>%
    dplyr::mutate(
      Label = sprintf("%.0f%%", 100 * .data$DifferentialFrequency),
      Edge = forcats::fct_reorder(.data$Edge, .data$MeanDeltaR)
    )

  if (nrow(dat) == 0) {
    warning("No stable edges available for heatmap.")
    return(invisible(NULL))
  }

  max_abs <- max(abs(dat$MeanDeltaR), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$Comparison, y = .data$Edge, fill = .data$MeanDeltaR)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Label), size = 2.0, colour = "black") +
    ggplot2::scale_fill_gradient2(
      low = params$palette$lower_in_B,
      mid = "white",
      high = params$palette$higher_in_B,
      midpoint = 0,
      limits = c(-max_abs, max_abs),
      name = expression(Delta * rho)
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = "Bootstrap-supported differential edges") +
    theme_publication_boot(base_size = 7) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  ggplot2::ggsave(outfile, p, width = 100, height = 90, units = "mm", device = svglite::svglite)
  invisible(p)
}

plot_candidate_delta_distributions <- function(boot_long, candidate_summary, outfile, params) {
  cand <- params$candidate_edges %>% normalize_edge() %>% dplyr::mutate(Keep = TRUE)

  dat <- boot_long %>%
    dplyr::inner_join(cand, by = c("Source", "Target")) %>%
    dplyr::mutate(
      Edge = paste(.data$Source, .data$Target, sep = " - "),
      DirectionLabel = dplyr::case_when(
        .data$DeltaR < 0 ~ "lower in B",
        .data$DeltaR > 0 ~ "higher in B",
        TRUE ~ "no change"
      ),
      DirectionLabel = factor(.data$DirectionLabel, levels = c("lower in B", "higher in B", "no change"))
    )

  if (nrow(dat) == 0) {
    warning("No candidate bootstrap distributions available.")
    return(invisible(NULL))
  }

  summary_dat <- candidate_summary %>%
    dplyr::mutate(Edge = paste(.data$Source, .data$Target, sep = " - "))

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$DeltaR, fill = .data$DirectionLabel)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey45") +
    ggplot2::geom_density(alpha = 0.55, linewidth = 0.25, colour = "black", adjust = 1.0) +
    ggplot2::geom_point(
      data = summary_dat,
      ggplot2::aes(x = .data$MeanDeltaR, y = 0),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.8,
      fill = "white",
      colour = "black",
      stroke = 0.25
    ) +
    ggplot2::facet_grid(Edge ~ Comparison, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("lower in B" = params$palette$lower_in_B, "higher in B" = params$palette$higher_in_B, "no change" = params$palette$neutral), name = NULL) +
    ggplot2::labs(x = expression(Delta * "Spearman " * rho), y = "Bootstrap density", title = "Candidate edge bootstrap distributions") +
    theme_publication_boot(base_size = 7) +
    ggplot2::theme(legend.position = "top")

  ggplot2::ggsave(outfile, p, width = 180, height = 145, units = "mm", device = svglite::svglite)
  invisible(p)
}

plot_minimal_networks <- function(edge_summary, out_dir, params) {
  dat_all <- edge_summary %>%
    dplyr::filter(.data$DifferentialFrequency >= params$stable_frequency_threshold)

  if (nrow(dat_all) == 0) {
    warning("No stable edges available for minimal network plots.")
    return(invisible(NULL))
  }

  for (cmp in unique(dat_all$Comparison)) {
    edges <- dat_all %>%
      dplyr::filter(.data$Comparison == cmp) %>%
      dplyr::mutate(
        DirectionLabel = dplyr::case_when(
          .data$MeanDeltaR < 0 ~ "lower in B",
          .data$MeanDeltaR > 0 ~ "higher in B",
          TRUE ~ "no change"
        )
      )

    if (nrow(edges) == 0) next

    nodes <- unique(c(edges$Source, edges$Target)) %>%
      tibble::tibble(name = .) %>%
      dplyr::mutate(
        Region = stringr::str_extract(.data$name, "CA1|CA2|CA3|DG"),
        Layer = stringr::str_replace(.data$name, "^(CA1|CA2|CA3|DG)_", "")
      )

    g <- igraph::graph_from_data_frame(
      edges %>%
        dplyr::transmute(
          from = .data$Source,
          to = .data$Target,
          DifferentialFrequency = .data$DifferentialFrequency,
          MeanDeltaR = .data$MeanDeltaR,
          DirectionLabel = .data$DirectionLabel
        ),
      directed = FALSE,
      vertices = nodes
    )

    p <- ggraph::ggraph(g, layout = "fr") +
      ggraph::geom_edge_link(
        ggplot2::aes(width = .data$DifferentialFrequency, alpha = .data$DifferentialFrequency, colour = .data$DirectionLabel),
        show.legend = TRUE
      ) +
      ggraph::geom_node_point(shape = 21, fill = "white", colour = "black", size = 3.5, stroke = 0.35) +
      ggraph::geom_node_text(ggplot2::aes(label = .data$name), repel = TRUE, size = 2.4) +
      ggraph::scale_edge_colour_manual(values = c("lower in B" = params$palette$lower_in_B, "higher in B" = params$palette$higher_in_B, "no change" = params$palette$neutral), name = NULL) +
      ggraph::scale_edge_width(range = c(0.25, 1.8), guide = "none") +
      ggraph::scale_edge_alpha(range = c(0.35, 0.95), guide = "none") +
      ggplot2::labs(title = paste0("Stable differential network: ", cmp)) +
      ggplot2::theme_void(base_size = 7) +
      ggplot2::theme(
        text = ggplot2::element_text(family = "sans", colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 8),
        legend.position = "bottom"
      )

    ggplot2::ggsave(
      file.path(out_dir, paste0(safe_name(cmp), "_minimal_stable_differential_network.svg")),
      p,
      width = 90,
      height = 85,
      units = "mm",
      device = svglite::svglite
    )
  }

  invisible(NULL)
}

# -------------------------------
# Main
# -------------------------------
message2("Starting Publication-style bootstrap differential network figure generation")
dirs <- make_dirs(params$bootstrap_dir)
write_session_info(file.path(CANONICAL_PATHS$logs, "sessionInfo.txt"))
loaded <- load_bootstrap_outputs(params)

edge_summary <- loaded$edge_summary
boot_long <- loaded$boot_long
candidate_summary <- loaded$candidate_summary

plot_forest_stable_edges(
  edge_summary,
  file.path(dirs$figures, "stable_differential_edges_forest.svg"),
  params
)

plot_delta_frequency_heatmap(
  edge_summary,
  file.path(dirs$figures, "stable_differential_edges_delta_frequency_heatmap.svg"),
  params
)

plot_candidate_delta_distributions(
  boot_long,
  candidate_summary,
  file.path(dirs$figures, "candidate_edge_bootstrap_delta_distributions.svg"),
  params
)

plot_minimal_networks(edge_summary, dirs$figures, params)

# Export a compact figure-ready table for slide/manuscript labels.
figure_table <- select_stable_edges(edge_summary, params) %>%
  dplyr::transmute(
    Comparison,
    Edge,
    Mean_R_A,
    Mean_R_B,
    MeanDeltaR,
    DeltaR_Q025,
    DeltaR_Q975,
    DifferentialFrequency,
    Prob_DeltaR_LessThan0,
    Prob_DeltaR_GreaterThan0,
    Interpretation = dplyr::case_when(
      MeanDeltaR < 0 ~ "weaker/lower in group B",
      MeanDeltaR > 0 ~ "stronger/higher in group B",
      TRUE ~ "no directional change"
    )
  )

readr::write_csv(figure_table, file.path(dirs$tables, "publication_ready_stable_differential_edges.csv"))
write_run_manifest(
  file.path(CANONICAL_PATHS$logs, "run_manifest.yml"),
  inputs = list(bootstrap_dir = params$bootstrap_dir),
  outputs = list(figures = dirs$figures, tables = dirs$tables, source_data = dirs$source_data),
  parameters = params,
  notes = "Figure script records stable frequency threshold and candidate edges used for plotting."
)

message2("Finished figure generation")
message2("Figure directory: ", dirs$figures)
