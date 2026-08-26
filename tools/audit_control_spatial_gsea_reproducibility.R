#!/usr/bin/env Rscript

# Read-only targeted comparison of one existing Figure-2 GO-BP result against
# the deterministic production-strength replay of the identical ranked list.

source(file.path("R", "paths.R"))
source(repo_path("R", "protein_group_enrichment_utils.R"))
source(repo_path("R", "clusterprofiler_reproducibility.R"))
source(repo_path("R", "control_spatial_identity_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  hit <- match(flag, args)
  if (!is.na(hit) && hit < length(args)) args[[hit + 1L]] else default
}
dataset <- arg_value("--dataset", "neuron_soma")
contrast <- arg_value("--contrast", "CA1_vs_mean_other_soma_regions")

if (!requireNamespace("readr", quietly = TRUE) ||
    !requireNamespace("clusterProfiler", quietly = TRUE) ||
    !requireNamespace("org.Mm.eg.db", quietly = TRUE) ||
    !requireNamespace("yaml", quietly = TRUE)) {
  stop(
    "readr, clusterProfiler, org.Mm.eg.db, and yaml are required.",
    call. = FALSE
  )
}

table_root <- path_results(
  "tables", "04_differential_expression_enrichment",
  "control_spatial_identity_validation", "global"
)
protein_path <- file.path(table_root, "anatomical_protein_contrasts.csv")
go_path <- file.path(table_root, "control_anatomical_go_bp_gsea.csv")
if (!file.exists(protein_path) || !file.exists(go_path)) {
  stop("Existing Figure-2 analytical tables are required.", call. = FALSE)
}

protein <- readr::read_csv(
  protein_path, show_col_types = FALSE, progress = FALSE
)
protein <- protein[
  protein$dataset == dataset & protein$contrast == contrast, , drop = FALSE
]
if (!nrow(protein)) stop("Requested protein contrast was not found.", call. = FALSE)
ranked <- build_enrichment_gene_inputs(protein, strict = TRUE)$ranked

cfg <- validate_clusterprofiler_gsea_config(yaml::read_yaml(
  repo_path("config", "clusterProfiler_config.yml")
))
identity <- control_spatial_gsea_reproducibility(
  dataset, contrast, "gseGO_BP",
  cfg$analysis$gsea_seed_base, cfg$analysis$n_perm_simple
)
replayed <- as.data.frame(run_seeded_clusterprofiler_gsea(
  clusterProfiler::gseGO,
  gsea_seed = identity$gsea_seed,
  n_perm_simple = identity$n_perm_simple,
  geneList = ranked, OrgDb = org.Mm.eg.db::org.Mm.eg.db,
  keyType = "SYMBOL", ont = "BP",
  pvalueCutoff = 1, verbose = FALSE
))
existing <- readr::read_csv(go_path, show_col_types = FALSE, progress = FALSE)
existing <- existing[
  existing$dataset == dataset & existing$contrast == contrast, , drop = FALSE
]

matched <- match(existing$ID, replayed$ID)
common <- !is.na(matched)
numeric_difference <- function(column) {
  if (!any(common)) return(NA_real_)
  max(abs(
    as.numeric(existing[[column]][common]) -
      as.numeric(replayed[[column]][matched[common]])
  ), na.rm = TRUE)
}
selected_ids <- function(data) {
  data <- data[
    is.finite(data$NES) & data$NES > 0 &
      is.finite(data$p.adjust) & data$p.adjust < 0.05,
    , drop = FALSE
  ]
  data <- data[order(data$p.adjust, -data$NES, data$ID), , drop = FALSE]
  utils::head(as.character(data$ID), 2L)
}

cat("dataset=", dataset, "\n", sep = "")
cat("contrast=", contrast, "\n", sep = "")
cat("gsea_seed=", identity$gsea_seed, "\n", sep = "")
cat("n_perm_simple=", identity$n_perm_simple, "\n", sep = "")
cat("existing_terms=", nrow(existing), "\n", sep = "")
cat("replayed_terms=", nrow(replayed), "\n", sep = "")
cat("identical_GO_ID_set=", setequal(existing$ID, replayed$ID), "\n", sep = "")
cat("common_GO_IDs=", sum(common), "\n", sep = "")
cat("max_abs_NES_difference=", numeric_difference("NES"), "\n", sep = "")
cat("max_abs_nominal_p_difference=", numeric_difference("pvalue"), "\n", sep = "")
cat("max_abs_adjusted_p_difference=", numeric_difference("p.adjust"), "\n", sep = "")
cat(
  "identical_common_leading_edges=",
  identical(
    as.character(existing$core_enrichment[common]),
    as.character(replayed$core_enrichment[matched[common]])
  ), "\n", sep = ""
)
existing_threshold <- sort(as.character(existing$ID[
  is.finite(existing$NES) & existing$NES > 0 &
    is.finite(existing$p.adjust) & existing$p.adjust < 0.05
]))
replayed_threshold <- sort(as.character(replayed$ID[
  is.finite(replayed$NES) & replayed$NES > 0 &
    is.finite(replayed$p.adjust) & replayed$p.adjust < 0.05
]))
cat(
  "identical_display_threshold_IDs=",
  identical(existing_threshold, replayed_threshold), "\n", sep = ""
)
cat("existing_display_threshold_count=", length(existing_threshold), "\n", sep = "")
cat("replayed_display_threshold_count=", length(replayed_threshold), "\n", sep = "")
cat(
  "threshold_IDs_removed_by_replay=",
  paste(setdiff(existing_threshold, replayed_threshold), collapse = ";"),
  "\n", sep = ""
)
cat(
  "threshold_IDs_added_by_replay=",
  paste(setdiff(replayed_threshold, existing_threshold), collapse = ";"),
  "\n", sep = ""
)
existing_top <- selected_ids(existing)
replayed_top <- selected_ids(replayed)
describe_ids <- function(data, ids) {
  idx <- match(ids, data$ID)
  paste(paste0(ids, "=", data$Description[idx]), collapse = ";")
}
cat("existing_top_two=", describe_ids(existing, existing_top), "\n", sep = "")
cat("replayed_top_two=", describe_ids(replayed, replayed_top), "\n", sep = "")
