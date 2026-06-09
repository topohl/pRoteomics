# ================================================================
# Script: 07_spatial_networks/04_bootstrap_differential_network_stability.r
# Stage: networks
# Scope: dataset_specific
# Consumes: required results/tables/07_spatial_networks/differential_networks/<dataset>/; optional none.
# Produces: results/tables/07_spatial_networks/bootstrap_differential_network_stability/01_Tables/bootstrap_differential_edge_stability_summary.csv.
# Dataset behavior: runs for neuron_neuropil,neuron_soma according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Bootstrap differential network stability.
# ================================================================

# ================================================================
# Consumes:
#   - group-specific spatial network RDS from canonical spatial-network output
# Produces:
#   - bootstrap differential network tables, figures, network files and RDS cache in canonical folders
# File contract:
#   - docs/active_script_io_audit.tsv object 07_spatial_networks/04_bootstrap_differential_network_stability.r
# ================================================================
# Bootstrap stability analysis for differential spatial networks
# ================================================================
# Purpose:
#   Test whether group differences in region/layer proteomic similarity
#   networks are stable under bootstrap resampling.
#
# Required first step:
#   Run Analysis/network_spatial_relations.r with run_group_specific_networks = TRUE.
#
# Key design:
#   This script resamples samples WITHIN each ExpGroup x RegionLayer stratum,
#   recomputes each region/layer mean expression profile, recomputes Spearman
#   spatial similarity networks, and then computes DeltaR = R_B - R_A for each
#   group comparison.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "spatial_network_utils.R"))
MODULE_ID <- "07_spatial_networks"
SUBSTEP_ID <- "bootstrap_differential_network_stability"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)
NETWORK_DATASET <- current_dataset()
assert_dataset_capability(NETWORK_DATASET, "region", analysis = "bootstrap differential spatial network stability analysis")
spatial_unit <- if (NETWORK_DATASET == "neuron_neuropil") "region_layer" else "region"

required_pkgs <- c(
  "dplyr", "tidyr", "stringr", "purrr", "tibble", "ggplot2",
  "pheatmap", "igraph", "ggraph", "openxlsx", "svglite"
)
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing required R package(s): ", paste(missing, collapse = ", "), ". Install them explicitly before running this script.", call. = FALSE)
resolve_spatial_rds <- function() {
  override <- Sys.getenv("PROTEOMICS_SPATIAL_NETWORK_OBJECT", unset = "")
  if (nzchar(override)) return(normalizePath(override, winslash = "/", mustWork = FALSE))
  scoped <- path_processed("07_spatial_networks", "network_spatial_relations", NETWORK_DATASET, spatial_unit, "network_spatial_relations_objects.rds")
  if (file.exists(scoped)) return(scoped)
  path_processed("07_spatial_networks", "network_spatial_relations", "network_spatial_relations_objects.rds")
}

params <- list(
  spatial_rds = resolve_spatial_rds(),
  output_dir = CANONICAL_PATHS$reports,

  group_order = c("CON", "RES", "SUS"),
  comparisons = list(c("CON", "RES"), c("CON", "SUS"), c("RES", "SUS")),

  n_boot = 500,
  seed = 1234,

  # Used to call an edge present in a single group during a bootstrap.
  edge_abs_r_threshold = 0.50,

  # Used to call a differential edge in a bootstrap.
  delta_abs_r_threshold = 0.10,

  # Used for consensus differential tables/plots.
  differential_frequency_threshold = 0.70,

  candidate_edges = tibble::tribble(
    ~Source, ~Target,
    "CA1_slm", "CA2_slm",
    "CA1_slm", "CA3_sr",
    "CA2_slm", "CA3_sr",
    "CA1_so",  "CA1_sr"
  )
)

if (is_dry_run()) {
  dry_run_line("Script", "07_spatial_networks/04_bootstrap_differential_network_stability.r")
  dry_run_line("Dataset", NETWORK_DATASET)
  dry_run_line("Spatial unit", spatial_unit)
  dry_run_line("Spatial RDS", params$spatial_rds, if (file.exists(params$spatial_rds)) "PASS" else "FAIL")
  dry_run_line("Output folders", paste(unlist(CANONICAL_PATHS), collapse = "; "))
  quit(status = if (file.exists(params$spatial_rds)) 0 else 1, save = "no")
}
if (!file.exists(params$spatial_rds)) stop("spatial_rds not found: ", params$spatial_rds, call. = FALSE)
invisible(lapply(required_pkgs, library, character.only = TRUE))

message2 <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
safe_name <- function(x) stringr::str_replace_all(as.character(x), "[^A-Za-z0-9_\\-]+", "_")

make_dirs <- function(base_dir) {
  dirs <- list(
    base = base_dir,
    tables = CANONICAL_PATHS$tables,
    figures = CANONICAL_PATHS$figures,
    networks = file.path(CANONICAL_PATHS$processed, "network_files"),
    logs = CANONICAL_PATHS$logs
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
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

clean_unit_matrix <- function(unit_matrix) {
  unit_matrix <- as.matrix(unit_matrix)
  storage.mode(unit_matrix) <- "numeric"
  unit_matrix[!is.finite(unit_matrix)] <- NA_real_

  keep_complete_enough <- rowSums(!is.na(unit_matrix)) >= 2
  unit_matrix <- unit_matrix[keep_complete_enough, , drop = FALSE]

  if (nrow(unit_matrix) < 3 || ncol(unit_matrix) < 3) return(NULL)

  row_sd <- apply(unit_matrix, 1, stats::sd, na.rm = TRUE)
  keep_variable <- is.finite(row_sd) & row_sd > 0
  unit_matrix <- unit_matrix[keep_variable, , drop = FALSE]

  if (nrow(unit_matrix) < 3 || ncol(unit_matrix) < 3) return(NULL)
  unit_matrix
}

make_unit_matrix_stratified_boot <- function(expr, sample_md, group) {
  expr <- as.matrix(expr)
  storage.mode(expr) <- "numeric"
  expr[!is.finite(expr)] <- NA_real_

  md_group <- sample_md %>%
    dplyr::filter(
      .data$ExpGroup == .env$group,
      !is.na(.data$RegionLayer),
      .data$SampleColumn %in% colnames(expr)
    ) %>%
    dplyr::distinct(.data$SampleColumn, .keep_all = TRUE)

  if (nrow(md_group) < 4) return(NULL)

  region_layers <- sort(unique(as.character(md_group$RegionLayer)))
  unit_profiles <- list()

  for (rl in region_layers) {
    cols <- md_group$SampleColumn[as.character(md_group$RegionLayer) == rl]
    cols <- cols[cols %in% colnames(expr)]
    if (length(cols) < 2) next

    boot_cols <- sample(cols, size = length(cols), replace = TRUE)
    prof <- rowMeans(expr[, boot_cols, drop = FALSE], na.rm = TRUE)
    prof[!is.finite(prof)] <- NA_real_
    unit_profiles[[rl]] <- prof
  }

  if (length(unit_profiles) < 3) return(NULL)

  unit_matrix <- do.call(cbind, unit_profiles)
  rownames(unit_matrix) <- rownames(expr)
  colnames(unit_matrix) <- names(unit_profiles)

  clean_unit_matrix(unit_matrix)
}

cor_edges <- function(unit_matrix, group) {
  unit_matrix <- clean_unit_matrix(unit_matrix)
  if (is.null(unit_matrix)) return(tibble::tibble())

  cor_mat <- suppressWarnings(cor(unit_matrix, method = "spearman", use = "pairwise.complete.obs"))
  cor_mat[!is.finite(cor_mat)] <- NA_real_
  units <- colnames(cor_mat)

  if (length(units) < 2) return(tibble::tibble())

  pairs <- utils::combn(units, 2, simplify = FALSE)
  edges <- purrr::map_dfr(pairs, function(pair) {
    r <- cor_mat[pair[1], pair[2]]
    tibble::tibble(
      Source = pair[1],
      Target = pair[2],
      SpearmanR = as.numeric(r)
    )
  })

  edges %>%
    dplyr::filter(!is.na(.data$SpearmanR)) %>%
    dplyr::mutate(
      Group = group,
      AbsSpearmanR = abs(.data$SpearmanR)
    )
}

bootstrap_group_edges <- function(expr, sample_md, group) {
  unit_matrix <- make_unit_matrix_stratified_boot(expr, sample_md, group)
  if (is.null(unit_matrix)) return(tibble::tibble())
  cor_edges(unit_matrix, group)
}

compare_edge_sets <- function(edges_a, edges_b, group_a, group_b, params) {
  ea <- edges_a %>%
    dplyr::transmute(Source = .data$Source, Target = .data$Target, R_A = .data$SpearmanR, AbsR_A = .data$AbsSpearmanR)
  eb <- edges_b %>%
    dplyr::transmute(Source = .data$Source, Target = .data$Target, R_B = .data$SpearmanR, AbsR_B = .data$AbsSpearmanR)

  dplyr::full_join(ea, eb, by = c("Source", "Target")) %>%
    dplyr::mutate(
      R_A = dplyr::coalesce(.data$R_A, 0),
      R_B = dplyr::coalesce(.data$R_B, 0),
      AbsR_A = dplyr::coalesce(.data$AbsR_A, 0),
      AbsR_B = dplyr::coalesce(.data$AbsR_B, 0),
      Group_A = group_a,
      Group_B = group_b,
      Comparison = paste0(group_a, "_vs_", group_b),
      DeltaR = .data$R_B - .data$R_A,
      AbsDeltaR = abs(.data$DeltaR),
      Present_A = .data$AbsR_A >= params$edge_abs_r_threshold,
      Present_B = .data$AbsR_B >= params$edge_abs_r_threshold,
      Differential = .data$AbsDeltaR >= params$delta_abs_r_threshold,
      Direction = dplyr::case_when(
        .data$DeltaR > 0 ~ "higher_in_B",
        .data$DeltaR < 0 ~ "lower_in_B",
        TRUE ~ "no_change"
      ),
      EdgeClass = dplyr::case_when(
        .data$Present_A & !.data$Present_B ~ "lost_in_B",
        !.data$Present_A & .data$Present_B ~ "gained_in_B",
        .data$Present_A & .data$Present_B & .data$DeltaR >= params$delta_abs_r_threshold ~ "strengthened_in_B",
        .data$Present_A & .data$Present_B & .data$DeltaR <= -params$delta_abs_r_threshold ~ "weakened_in_B",
        .data$Present_A & .data$Present_B ~ "preserved",
        TRUE ~ "weak_in_both"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(.data$AbsDeltaR))
}

bootstrap_iteration <- function(expr, sample_md, params, i) {
  group_edges <- purrr::map(params$group_order, function(grp) {
    bootstrap_group_edges(expr, sample_md, grp)
  })
  names(group_edges) <- params$group_order

  purrr::map_dfr(params$comparisons, function(cmp) {
    group_a <- cmp[1]
    group_b <- cmp[2]

    edges_a <- group_edges[[group_a]]
    edges_b <- group_edges[[group_b]]

    if (is.null(edges_a) || is.null(edges_b) || nrow(edges_a) == 0 || nrow(edges_b) == 0) {
      return(tibble::tibble())
    }

    compare_edge_sets(edges_a, edges_b, group_a, group_b, params)
  }) %>%
    dplyr::mutate(Bootstrap = i)
}

summarise_differential_bootstrap <- function(boot_tbl, params) {
  boot_tbl %>%
    dplyr::group_by(.data$Comparison, .data$Group_A, .data$Group_B, .data$Source, .data$Target) %>%
    dplyr::summarise(
      Mean_R_A = mean(.data$R_A, na.rm = TRUE),
      Mean_R_B = mean(.data$R_B, na.rm = TRUE),
      MeanDeltaR = mean(.data$DeltaR, na.rm = TRUE),
      SDDeltaR = sd(.data$DeltaR, na.rm = TRUE),
      MedianDeltaR = median(.data$DeltaR, na.rm = TRUE),
      DeltaR_Q025 = as.numeric(stats::quantile(.data$DeltaR, 0.025, na.rm = TRUE)),
      DeltaR_Q975 = as.numeric(stats::quantile(.data$DeltaR, 0.975, na.rm = TRUE)),
      MeanAbsDeltaR = mean(.data$AbsDeltaR, na.rm = TRUE),
      DifferentialFrequency = mean(.data$AbsDeltaR >= params$delta_abs_r_threshold, na.rm = TRUE),
      LowerInB_Frequency = mean(.data$DeltaR <= -params$delta_abs_r_threshold, na.rm = TRUE),
      HigherInB_Frequency = mean(.data$DeltaR >= params$delta_abs_r_threshold, na.rm = TRUE),
      LostInB_Frequency = mean(.data$EdgeClass == "lost_in_B", na.rm = TRUE),
      GainedInB_Frequency = mean(.data$EdgeClass == "gained_in_B", na.rm = TRUE),
      Preserved_Frequency = mean(.data$EdgeClass == "preserved", na.rm = TRUE),
      Prob_DeltaR_LessThan0 = mean(.data$DeltaR < 0, na.rm = TRUE),
      Prob_DeltaR_GreaterThan0 = mean(.data$DeltaR > 0, na.rm = TRUE),
      NBoot = dplyr::n(),
      StableDifferentialEdge = .data$DifferentialFrequency >= params$differential_frequency_threshold,
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$Comparison, dplyr::desc(.data$DifferentialFrequency), dplyr::desc(.data$MeanAbsDeltaR))
}

extract_candidate_edges <- function(summary_tbl, candidate_edges) {
  cand <- candidate_edges %>% normalize_edge()

  summary_tbl %>%
    normalize_edge() %>%
    dplyr::inner_join(cand, by = c("Source", "Target")) %>%
    dplyr::arrange(.data$Source, .data$Target, .data$Comparison)
}

plot_delta_heatmap <- function(summary_tbl, comparison, outfile) {
  dat <- summary_tbl %>% dplyr::filter(.data$Comparison == comparison)
  if (nrow(dat) == 0) return(invisible(NULL))

  units <- sort(unique(c(dat$Source, dat$Target)))
  mat <- matrix(NA_real_, nrow = length(units), ncol = length(units), dimnames = list(units, units))

  for (i in seq_len(nrow(dat))) {
    r <- dat[i, ]
    mat[r$Source, r$Target] <- r$MeanDeltaR
    mat[r$Target, r$Source] <- r$MeanDeltaR
  }
  diag(mat) <- 0

  max_abs <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

  pheatmap::pheatmap(
    mat,
    color = colorRampPalette(c("#457B9D", "white", "#E63946"))(101),
    breaks = seq(-max_abs, max_abs, length.out = 102),
    border_color = NA,
    fontsize = 7,
    main = paste0("Bootstrap mean DeltaR: ", comparison),
    filename = outfile,
    width = 6,
    height = 5
  )
}

plot_stable_differential_network <- function(summary_tbl, comparison, outfile, params) {
  edges <- summary_tbl %>%
    dplyr::filter(
      .data$Comparison == comparison,
      .data$DifferentialFrequency >= params$differential_frequency_threshold
    )

  if (nrow(edges) == 0) {
    warning("No stable differential edges passed threshold for ", comparison)
    return(invisible(NULL))
  }

  nodes <- unique(c(edges$Source, edges$Target)) %>%
    tibble::tibble(name = .) %>%
    dplyr::distinct(.data$name, .keep_all = TRUE)

  g <- igraph::graph_from_data_frame(
    edges %>%
      dplyr::transmute(
        from = .data$Source,
        to = .data$Target,
        DifferentialFrequency = .data$DifferentialFrequency,
        MeanDeltaR = .data$MeanDeltaR,
        Direction = dplyr::if_else(.data$MeanDeltaR >= 0, "higher_in_B", "lower_in_B")
      ),
    directed = FALSE,
    vertices = nodes
  )

  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(aes(width = DifferentialFrequency, alpha = DifferentialFrequency, linetype = Direction), colour = "grey35") +
    ggraph::geom_node_point(shape = 21, fill = "white", colour = "black", size = 4, stroke = 0.4) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 2.8) +
    ggraph::scale_edge_width(range = c(0.2, 1.8), guide = "none") +
    ggraph::scale_edge_alpha(range = c(0.25, 0.9), guide = "none") +
    ggplot2::labs(title = paste0("Stable differential network: ", comparison)) +
    ggplot2::theme_void(base_size = 8) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  ggplot2::ggsave(outfile, p, width = 150, height = 120, units = "mm", device = svglite::svglite)
  invisible(p)
}

write_graphml <- function(summary_tbl, comparison, outfile, params) {
  edges <- summary_tbl %>%
    dplyr::filter(
      .data$Comparison == comparison,
      .data$DifferentialFrequency >= params$differential_frequency_threshold
    )

  if (nrow(edges) == 0) return(invisible(NULL))

  nodes <- unique(c(edges$Source, edges$Target)) %>%
    tibble::tibble(name = .) %>%
    dplyr::distinct(.data$name, .keep_all = TRUE)

  g <- igraph::graph_from_data_frame(
    edges %>% dplyr::transmute(
      from = .data$Source,
      to = .data$Target,
      DifferentialFrequency = .data$DifferentialFrequency,
      MeanDeltaR = .data$MeanDeltaR,
      SDDeltaR = .data$SDDeltaR,
      LowerInB_Frequency = .data$LowerInB_Frequency,
      HigherInB_Frequency = .data$HigherInB_Frequency
    ),
    directed = FALSE,
    vertices = nodes
  )
  igraph::write_graph(g, outfile, format = "graphml")
}

# -------------------------------
# Main
# -------------------------------
set.seed(params$seed)
dirs <- make_dirs(params$output_dir)
write_session_info(file.path(dirs$logs, "sessionInfo.txt"))
log_file <- file.path(dirs$logs, paste0("bootstrap_differential_network_stability_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message2("Starting bootstrap differential network stability analysis")
message2("Bootstrap iterations: ", params$n_boot)
message2("Edge threshold |r| >= ", params$edge_abs_r_threshold)
message2("Delta threshold |DeltaR| >= ", params$delta_abs_r_threshold)

if (!file.exists(params$spatial_rds)) stop("spatial_rds not found. Run network_spatial_relations.r first: ", params$spatial_rds)
obj <- readRDS(params$spatial_rds)
expr <- obj$expression_matrix
sample_md <- obj$sample_metadata

if (!all(c("SampleColumn", "RegionLayer", "ExpGroup") %in% names(sample_md))) {
  stop("sample_metadata must contain SampleColumn, RegionLayer, and ExpGroup.")
}

message2("Available group counts:")
print(sample_md %>% dplyr::count(.data$ExpGroup))

message2("Samples matching expression matrix by group:")
print(sample_md %>% dplyr::mutate(InExpressionMatrix = .data$SampleColumn %in% colnames(expr)) %>% dplyr::count(.data$ExpGroup, .data$InExpressionMatrix))

message2("Region/layer counts by group:")
print(sample_md %>% dplyr::filter(.data$SampleColumn %in% colnames(expr)) %>% dplyr::count(.data$ExpGroup, .data$RegionLayer))

message2("One-pass edge sanity check:")
sanity <- purrr::map_dfr(params$group_order, function(grp) {
  unit_matrix <- make_unit_matrix_stratified_boot(expr, sample_md, grp)
  e <- if (is.null(unit_matrix)) tibble::tibble() else cor_edges(unit_matrix, grp)
  tibble::tibble(
    Group = grp,
    UnitRows = if (is.null(unit_matrix)) 0L else nrow(unit_matrix),
    UnitCols = if (is.null(unit_matrix)) 0L else ncol(unit_matrix),
    NEdges = nrow(e)
  )
})
print(sanity)
if (any(sanity$NEdges == 0)) stop("At least one group produced zero edges in the sanity check.")

boot_list <- vector("list", params$n_boot)
for (i in seq_len(params$n_boot)) {
  if (i %% 10 == 0) message2("Bootstrap iteration ", i, " / ", params$n_boot)
  boot_list[[i]] <- bootstrap_iteration(expr, sample_md, params, i)
}

boot_tbl <- dplyr::bind_rows(boot_list)
if (nrow(boot_tbl) == 0) stop("No bootstrap differential edge results were produced.")

utils::write.csv(boot_tbl, file.path(dirs$tables, "bootstrap_differential_edge_values_long.csv"), row.names = FALSE)

edge_summary <- summarise_differential_bootstrap(boot_tbl, params)
utils::write.csv(edge_summary, file.path(dirs$tables, "bootstrap_differential_edge_stability_summary.csv"), row.names = FALSE)

edge_validation <- edge_summary %>%
  dplyr::mutate(
    edge_id = make_edge_id(.data$Source, .data$Target),
    dataset = NETWORK_DATASET,
    group = .data$Comparison,
    observed_r = .data$MeanDeltaR,
    bootstrap_median = .data$MedianDeltaR,
    bootstrap_ci_low = .data$DeltaR_Q025,
    bootstrap_ci_high = .data$DeltaR_Q975,
    permutation_p = pmin(.data$Prob_DeltaR_LessThan0, .data$Prob_DeltaR_GreaterThan0, na.rm = TRUE) * 2,
    permutation_p = pmin(.data$permutation_p, 1),
    fdr = stats::p.adjust(.data$permutation_p, method = "BH"),
    n_animals = dplyr::n_distinct(sample_md$AnimalID),
    n_proteins = nrow(expr),
    stability_score = .data$DifferentialFrequency,
    interpretation_strength = vapply(seq_len(dplyr::n()), function(i) {
      network_interpretation_strength(fdr = fdr[[i]], stability_score = stability_score[[i]], n_animals = n_animals[[i]])
    }, character(1))
  ) %>%
  dplyr::rename(source = Source, target = Target) %>%
  dplyr::select(
    "edge_id", "source", "target", "dataset", "group",
    "observed_r", "bootstrap_median", "bootstrap_ci_low", "bootstrap_ci_high",
    "permutation_p", "fdr", "n_animals", "n_proteins", "stability_score",
    "interpretation_strength", dplyr::everything()
  )
utils::write.csv(edge_validation, file.path(dirs$tables, "edge_bootstrap_permutation_validation_summary.csv"), row.names = FALSE)

candidate_summary <- extract_candidate_edges(edge_summary, params$candidate_edges)
utils::write.csv(candidate_summary, file.path(dirs$tables, "candidate_edge_differential_stability_summary.csv"), row.names = FALSE)

stable_edges <- edge_summary %>%
  dplyr::filter(.data$DifferentialFrequency >= params$differential_frequency_threshold)
utils::write.csv(stable_edges, file.path(dirs$tables, "stable_differential_edges.csv"), row.names = FALSE)

for (cmp in unique(edge_summary$Comparison)) {
  cmp_safe <- safe_name(cmp)
  plot_delta_heatmap(edge_summary, cmp, file.path(dirs$figures, paste0(cmp_safe, "_bootstrap_mean_deltaR_heatmap.pdf")))
  plot_stable_differential_network(edge_summary, cmp, file.path(dirs$figures, paste0(cmp_safe, "_stable_differential_network.svg")), params)

  cmp_edges <- stable_edges %>% dplyr::filter(.data$Comparison == cmp)
  cmp_nodes <- unique(c(cmp_edges$Source, cmp_edges$Target)) %>% tibble::tibble(name = .)
  utils::write.csv(cmp_edges, file.path(dirs$networks, paste0(cmp_safe, "_stable_differential_edges.csv")), row.names = FALSE)
  utils::write.csv(cmp_nodes, file.path(dirs$networks, paste0(cmp_safe, "_stable_differential_nodes.csv")), row.names = FALSE)
  write_graphml(edge_summary, cmp, file.path(dirs$networks, paste0(cmp_safe, "_stable_differential.graphml")), params)
}

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Edge_Summary")
openxlsx::writeData(wb, "Edge_Summary", edge_summary)
openxlsx::addWorksheet(wb, "Candidate_Edges")
openxlsx::writeData(wb, "Candidate_Edges", candidate_summary)
openxlsx::addWorksheet(wb, "Stable_Differential")
openxlsx::writeData(wb, "Stable_Differential", stable_edges)
openxlsx::saveWorkbook(wb, file.path(dirs$tables, "bootstrap_differential_network_stability_summary.xlsx"), overwrite = TRUE)

saveRDS(
  list(
    params = params,
    bootstrap_long = boot_tbl,
    edge_summary = edge_summary,
    candidate_summary = candidate_summary,
    stable_edges = stable_edges,
    sessionInfo = sessionInfo()
  ),
  file.path(dirs$logs, "bootstrap_differential_network_stability_objects.rds")
)
write_run_manifest(
  file.path(dirs$logs, "run_manifest.yml"),
  inputs = list(spatial_rds = params$spatial_rds),
  outputs = list(
    bootstrap_long = file.path(dirs$tables, "bootstrap_differential_edge_values_long.csv"),
    edge_summary = file.path(dirs$tables, "bootstrap_differential_edge_stability_summary.csv"),
    edge_validation = file.path(dirs$tables, "edge_bootstrap_permutation_validation_summary.csv"),
    stable_edges = file.path(dirs$tables, "stable_differential_edges.csv"),
    networks = dirs$networks
  ),
  parameters = params,
  notes = "Bootstrap seed, edge threshold, delta threshold and differential frequency threshold are recorded in parameters."
)

message2("Finished bootstrap differential network stability analysis")
message2("Interpret DifferentialFrequency as how often |DeltaR| exceeds the threshold across bootstrap resamples.")
