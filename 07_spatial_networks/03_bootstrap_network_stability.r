# ================================================================
# Consumes:
#   - spatial network RDS from data/processed/07_spatial_networks/network_spatial_relations/
# Produces:
#   - bootstrap stability tables, figures, network files and RDS cache in canonical folders
# File contract:
#   - docs/active_script_io_audit.tsv object 07_spatial_networks/03_bootstrap_network_stability.r
# ================================================================
# Bootstrap stability analysis for spatial proteomics networks
# ================================================================
# Purpose:
#   Estimate robustness of region/layer network edges by repeated bootstrap
#   resampling of samples. This provides edge stability metrics and consensus
#   networks to avoid overinterpreting unstable single edges.
#
# Required first step:
#   Run Analysis/network_spatial_relations.r
#
# Outputs:
#   - edge stability frequencies
#   - bootstrap distributions of edge weights
#   - consensus networks
#   - node stability summaries
#   - stability heatmaps and network plots
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
MODULE_ID <- "07_spatial_networks"
SUBSTEP_ID <- "bootstrap_network_stability"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

params <- list(
  spatial_rds = path_processed("07_spatial_networks", "network_spatial_relations", "network_spatial_relations_objects.rds"),
  output_dir = CANONICAL_PATHS$reports,
  n_boot = 250,
  edge_abs_r_threshold = 0.50,
  consensus_frequency_threshold = 0.70,
  seed = 1234
)

if (is_dry_run()) {
  dry_run_line("Script", "07_spatial_networks/03_bootstrap_network_stability.r")
  dry_run_line("Spatial RDS", params$spatial_rds, if (file.exists(params$spatial_rds)) "PASS" else "FAIL")
  dry_run_line("Output folders", paste(unlist(CANONICAL_PATHS), collapse = "; "))
  quit(status = if (file.exists(params$spatial_rds)) 0 else 1, save = "no")
}
if (!file.exists(params$spatial_rds)) stop("spatial_rds not found: ", params$spatial_rds, call. = FALSE)

required_pkgs <- c("dplyr", "tidyr", "stringr", "purrr", "tibble", "ggplot2", "pheatmap", "igraph", "ggraph", "openxlsx", "svglite")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")
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

aggregate_region_layer_expression <- function(expr, sample_md) {
  sample_order <- match(colnames(expr), sample_md$SampleColumn)
  md <- sample_md[sample_order, ]

  long <- as.data.frame(expr) %>%
    tibble::rownames_to_column("Protein") %>%
    tidyr::pivot_longer(-Protein, names_to = "SampleColumn", values_to = "Expression") %>%
    left_join(md %>% dplyr::select(SampleColumn, RegionLayer), by = "SampleColumn")

  agg <- long %>%
    filter(!is.na(RegionLayer)) %>%
    group_by(Protein, RegionLayer) %>%
    summarise(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = RegionLayer, values_from = MeanExpression) %>%
    tibble::column_to_rownames("Protein") %>%
    as.matrix()

  agg[rowSums(!is.na(agg)) > 1, , drop = FALSE]
}

cor_edges <- function(unit_matrix) {
  cor_mat <- suppressWarnings(cor(unit_matrix, method = "spearman", use = "pairwise.complete.obs"))
  units <- colnames(cor_mat)
  expand.grid(Source = units, Target = units, stringsAsFactors = FALSE) %>%
    as_tibble() %>%
    filter(Source < Target) %>%
    rowwise() %>%
    mutate(SpearmanR = cor_mat[Source, Target]) %>%
    ungroup() %>%
    mutate(AbsSpearmanR = abs(SpearmanR))
}

bootstrap_iteration <- function(expr, sample_md, threshold) {
  samp <- sample(sample_md$SampleColumn, size = nrow(sample_md), replace = TRUE)

  md_boot <- tibble::tibble(
    SampleColumn = samp,
    BootstrapSampleColumn = paste0(samp, "__boot", seq_along(samp))
  ) %>%
    dplyr::left_join(sample_md, by = "SampleColumn") %>%
    dplyr::mutate(SampleColumn = .data$BootstrapSampleColumn) %>%
    dplyr::select(-.data$BootstrapSampleColumn)

  expr_boot <- expr[, samp, drop = FALSE]
  colnames(expr_boot) <- md_boot$SampleColumn

  unit_matrix <- aggregate_region_layer_expression(expr_boot, md_boot)
  edges <- cor_edges(unit_matrix)

  edges %>%
    dplyr::mutate(Present = .data$AbsSpearmanR >= threshold)
}

summarise_bootstrap_edges <- function(edge_boot_tbl, threshold) {
  edge_boot_tbl %>%
    group_by(Source, Target) %>%
    summarise(
      MeanR = mean(SpearmanR, na.rm = TRUE),
      SDR = sd(SpearmanR, na.rm = TRUE),
      MedianR = median(SpearmanR, na.rm = TRUE),
      MinR = min(SpearmanR, na.rm = TRUE),
      MaxR = max(SpearmanR, na.rm = TRUE),
      StabilityFrequency = mean(Present, na.rm = TRUE),
      NBoot = n(),
      StableEdge = StabilityFrequency >= threshold,
      .groups = "drop"
    ) %>%
    arrange(desc(StabilityFrequency), desc(abs(MeanR)))
}

plot_stability_heatmap <- function(edge_summary, outfile) {
  units <- sort(unique(c(edge_summary$Source, edge_summary$Target)))
  mat <- matrix(NA_real_, nrow = length(units), ncol = length(units), dimnames = list(units, units))

  for (i in seq_len(nrow(edge_summary))) {
    r <- edge_summary[i, ]
    mat[r$Source, r$Target] <- r$StabilityFrequency
    mat[r$Target, r$Source] <- r$StabilityFrequency
  }
  diag(mat) <- 1

  pheatmap::pheatmap(
    mat,
    color = colorRampPalette(c("white", "#457B9D", "#E63946"))(101),
    border_color = NA,
    fontsize = 7,
    main = "Bootstrap edge stability frequency",
    filename = outfile,
    width = 6,
    height = 5
  )
}

plot_consensus_network <- function(edge_summary, outfile, threshold) {
  edges <- edge_summary %>% filter(StabilityFrequency >= threshold)
  if (nrow(edges) == 0) {
    warning("No consensus edges passed threshold.")
    return(invisible(NULL))
  }

  nodes <- unique(c(edges$Source, edges$Target)) %>% tibble(name = .)
  g <- igraph::graph_from_data_frame(edges %>% transmute(from = Source, to = Target, StabilityFrequency, MeanR), directed = FALSE, vertices = nodes)

  p <- ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(aes(width = StabilityFrequency, alpha = StabilityFrequency), colour = "grey35") +
    ggraph::geom_node_point(shape = 21, fill = "white", colour = "black", size = 4, stroke = 0.4) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 2.8) +
    scale_edge_width(range = c(0.2, 1.8), guide = "none") +
    scale_edge_alpha(range = c(0.25, 0.9), guide = "none") +
    labs(title = paste0("Consensus network (frequency >= ", threshold, ")")) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(outfile, p, width = 150, height = 120, units = "mm", device = svglite::svglite)
  invisible(p)
}

node_stability_summary <- function(edge_summary) {
  bind_rows(
    edge_summary %>% transmute(Node = Source, StabilityFrequency, MeanR),
    edge_summary %>% transmute(Node = Target, StabilityFrequency, MeanR)
  ) %>%
    group_by(Node) %>%
    summarise(
      NEdges = n(),
      MeanStability = mean(StabilityFrequency, na.rm = TRUE),
      MedianStability = median(StabilityFrequency, na.rm = TRUE),
      MaxStability = max(StabilityFrequency, na.rm = TRUE),
      MeanAbsR = mean(abs(MeanR), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(MeanStability), desc(MaxStability))
}

write_graphml <- function(edge_summary, outfile, threshold) {
  edges <- edge_summary %>% filter(StabilityFrequency >= threshold)
  if (nrow(edges) == 0) return(invisible(NULL))
  nodes <- unique(c(edges$Source, edges$Target)) %>% tibble(name = .)
  g <- igraph::graph_from_data_frame(edges %>% transmute(from = Source, to = Target, StabilityFrequency, MeanR), directed = FALSE, vertices = nodes)
  igraph::write_graph(g, outfile, format = "graphml")
}

# -------------------------------
# Main
# -------------------------------
set.seed(params$seed)
dirs <- make_dirs(params$output_dir)
write_session_info(file.path(dirs$logs, "sessionInfo.txt"))
log_file <- file.path(dirs$logs, paste0("bootstrap_network_stability_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message2("Starting bootstrap network stability analysis")
if (!file.exists(params$spatial_rds)) stop("spatial_rds not found. Run network_spatial_relations.r first.")
obj <- readRDS(params$spatial_rds)
expr <- obj$expression_matrix
sample_md <- obj$sample_metadata

message2("Bootstrap iterations: ", params$n_boot)
message2("Edge threshold |r| >= ", params$edge_abs_r_threshold)

boot_list <- vector("list", params$n_boot)
for (i in seq_len(params$n_boot)) {
  if (i %% 10 == 0) message2("Bootstrap iteration ", i, " / ", params$n_boot)
  boot_list[[i]] <- bootstrap_iteration(expr, sample_md, params$edge_abs_r_threshold) %>% mutate(Bootstrap = i)
}

boot_tbl <- bind_rows(boot_list)
utils::write.csv(boot_tbl, file.path(dirs$tables, "bootstrap_edge_values_long.csv"), row.names = FALSE)

edge_summary <- summarise_bootstrap_edges(boot_tbl, params$consensus_frequency_threshold)
utils::write.csv(edge_summary, file.path(dirs$tables, "bootstrap_edge_stability_summary.csv"), row.names = FALSE)

node_summary <- node_stability_summary(edge_summary)
utils::write.csv(node_summary, file.path(dirs$tables, "bootstrap_node_stability_summary.csv"), row.names = FALSE)

consensus_edges <- edge_summary %>% filter(StabilityFrequency >= params$consensus_frequency_threshold)
utils::write.csv(consensus_edges, file.path(dirs$tables, "consensus_edges.csv"), row.names = FALSE)

plot_stability_heatmap(edge_summary, file.path(dirs$figures, "bootstrap_edge_stability_heatmap.pdf"))
plot_consensus_network(edge_summary, file.path(dirs$figures, "consensus_network.svg"), params$consensus_frequency_threshold)

utils::write.csv(consensus_edges, file.path(dirs$networks, "consensus_network_edges.csv"), row.names = FALSE)
utils::write.csv(unique(c(consensus_edges$Source, consensus_edges$Target)) %>% tibble(name = .), file.path(dirs$networks, "consensus_network_nodes.csv"), row.names = FALSE)
write_graphml(edge_summary, file.path(dirs$networks, "consensus_network.graphml"), params$consensus_frequency_threshold)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Edge_Summary")
openxlsx::writeData(wb, "Edge_Summary", edge_summary)
openxlsx::addWorksheet(wb, "Node_Summary")
openxlsx::writeData(wb, "Node_Summary", node_summary)
openxlsx::addWorksheet(wb, "Consensus_Edges")
openxlsx::writeData(wb, "Consensus_Edges", consensus_edges)
openxlsx::saveWorkbook(wb, file.path(dirs$tables, "bootstrap_network_stability_summary.xlsx"), overwrite = TRUE)

saveRDS(
  list(params = params, bootstrap_long = boot_tbl, edge_summary = edge_summary, node_summary = node_summary, consensus_edges = consensus_edges, sessionInfo = sessionInfo()),
  file.path(dirs$logs, "bootstrap_network_stability_objects.rds")
)

message2("Finished bootstrap network stability analysis")
message2("Consensus edges are substantially more reliable than single-run edges.")
