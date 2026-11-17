#!/usr/bin/env Rscript
# ================================================================
# Script: analysis/differential_abundance/summarize_top_term_drivers.R
# Stage: enrichment
# Scope: dataset_specific
# Consumes: required results/differential_abundance/compare_go_enrichment/<dataset>/models/compareGO_input_manifest.csv; data/processed/04_differential_expression_enrichment/compareGO/<dataset>/compareGO_input_manifest.csv; config/compareGO_config.yml; optional none declared in pipeline.yml
# Produces: results/differential_abundance/summarize_top_term_drivers/<dataset>/tables/top_term_driver_genes.csv; results/differential_abundance/summarize_top_term_drivers/<dataset>/tables/selected_top_terms.csv; results/differential_abundance/summarize_top_term_drivers/<dataset>/manifests/run_manifest.yml
# Dataset behavior: runs for neuron_neuropil,neuron_soma,microglia according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Top GO terms per comparison and the genes that drive them, from the canonical compareGO bundle (GSEA_GO, BP, phenotype_within_unit).
# ================================================================
#
# Restores, on canonical inputs, the top-term selection and the top-gene
# driver table (section 10) of the archived compareGO tail
# (archive/04_differential_expression_enrichment/legacy/02_compareGO_superseded_tail.r,
# provenance only; nothing here sources or parses it). Rules are in
# R/enrichment/top_term_driver_utils.R:
#   - best instance (largest |NES|) of each term in each comparison;
#   - candidate pool p.adjust < 0.05 (significant_only);
#   - per comparison the top target_n_terms terms up and down by NES with
#     setSize >= min_set_size, dropping a term whose leading edge has Jaccard >
#     redundancy_threshold with a term already kept;
#   - leading edge = term-gene provenance rows with core_enrichment_member and
#     gene_level_claim_allowed TRUE;
#   - for the selected terms, over all comparisons: Freq (comparisons in which
#     the gene is in the term's leading edge), Mean_NES (mean |NES|),
#     Comparisons, n_significant_comparisons; top 100 genes per term.
#
# Usage: Rscript analysis/differential_abundance/summarize_top_term_drivers.R --dataset <dataset> [--dry-run]

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "enrichment_io.R"))
source(repo_path("R", "differential_abundance_paths.R"))
source(repo_path("R", "top_term_driver_utils.R"))

ANALYSIS_ID <- "summarize_top_term_drivers"
SCRIPT_ID <- "analysis/differential_abundance/summarize_top_term_drivers.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
DATASET <- current_dataset_from_cli()
PATHS <- differential_abundance_dirs(ANALYSIS_ID, scope = DATASET, create = FALSE)

config_candidates <- c(repo_path("config", "compareGO_config.local.yml"),
                       repo_path("config", "compareGO_config.yml"))
config_file <- config_candidates[file.exists(config_candidates)][1]
config <- if (!is.na(config_file) && requireNamespace("yaml", quietly = TRUE)) yaml::read_yaml(config_file) else list()
PARAMS <- top_term_driver_parameters(config)

manifest_file <- canonical_comparego_manifest_path(DATASET)
driver_file <- file.path(PATHS$tables, "top_term_driver_genes.csv")
selection_file <- file.path(PATHS$tables, "selected_top_terms.csv")
run_manifest_file <- file.path(PATHS$manifests, "run_manifest.yml")

if (is_dry_run()) {
  if (file.exists(manifest_file)) {
    read_canonical_comparego_manifest(manifest_file, DATASET, require_files = FALSE)
  }
  dry_run_line("Script", SCRIPT_ID)
  dry_run_line("Dataset", DATASET)
  dry_run_line("compareGO manifest", manifest_file, if (file.exists(manifest_file)) "PASS" else "FAIL")
  dry_run_line("Selection config", if (is.na(config_file)) "defaults" else config_file, "INFO")
  dry_run_line("Parameters", paste(names(PARAMS), vapply(PARAMS, as.character, character(1)),
                                   sep = "=", collapse = "; "), "INFO")
  dry_run_line("Top-term driver genes", driver_file)
  dry_run_line("Selected top terms", selection_file)
  dry_run_line("Run manifest", run_manifest_file)
  quit(status = if (file.exists(manifest_file)) 0 else 1, save = "no")
}

required_pkgs <- c("dplyr", "readr", "yaml")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}
if (!file.exists(manifest_file)) {
  stop("Canonical compareGO manifest not found for dataset ", DATASET, ": ", manifest_file, call. = FALSE)
}
dir_create(PATHS$tables)
dir_create(PATHS$manifests)

write_csv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  readr::write_csv(x, tmp, na = "")
  if (!file.exists(tmp)) stop("Incomplete write: ", path, call. = FALSE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not publish completed output: ", path, call. = FALSE)
  invisible(path)
}

# ---- canonical inputs ------------------------------------------------------

message("Reading the canonical compareGO bundle for ", DATASET, " ...")
bundle <- read_canonical_comparego_bundle(
  manifest_file, DATASET,
  provenance_columns = c("member_accessions", "core_enrichment_member")
)
admitted <- bundle$manifest[
  bundle$manifest$route_category == TOP_TERM_DRIVER_ROUTE &
    bundle$manifest$result_type == TOP_TERM_DRIVER_RESULT_TYPE &
    toupper(bundle$manifest$ontology) == TOP_TERM_DRIVER_ONTOLOGY &
    bundle$manifest$comparego_analysis_status == "included", , drop = FALSE]
if (!nrow(admitted)) {
  stop("Canonical compareGO manifest has no admitted ", TOP_TERM_DRIVER_RESULT_TYPE, "/",
    TOP_TERM_DRIVER_ONTOLOGY, " ", TOP_TERM_DRIVER_ROUTE, " rows for ", DATASET, ".", call. = FALSE)
}
route_units <- unique(admitted[c("comparison", "route_unit")])
if (anyDuplicated(route_units$comparison)) {
  stop("Canonical compareGO manifest has more than one route unit for a comparison in ", DATASET, ".", call. = FALSE)
}
comparisons <- sort(unique(as.character(admitted$comparison)), method = "radix")

terms <- bundle$terms[bundle$terms$comparison %in% comparisons &
  bundle$terms$result_type == TOP_TERM_DRIVER_RESULT_TYPE &
  toupper(bundle$terms$ontology) == TOP_TERM_DRIVER_ONTOLOGY, , drop = FALSE]
provenance <- bundle$provenance[bundle$provenance$comparison %in% comparisons, , drop = FALSE]
n_provenance_rows <- nrow(bundle$provenance)
bundle$provenance <- NULL
core_rows <- top_term_core_rows(provenance)
rm(provenance)
invisible(gc())

# ---- selection -------------------------------------------------------------

term_table <- top_term_table(terms)
best <- top_term_best_instances(term_table)
leading_edges <- top_term_leading_edges(core_rows)
selected <- top_term_selection(best, leading_edges, PARAMS)
top_terms <- unique(selected$Description)
if (!length(top_terms)) {
  message("No term passes the selection for ", DATASET, "; writing empty tables.")
}
drivers <- top_term_driver_genes(term_table, core_rows, top_terms, PARAMS)

selected_out <- data.frame(
  dataset = rep(DATASET, nrow(selected)),
  comparison = selected$comparison,
  route_unit = route_units$route_unit[match(selected$comparison, route_units$comparison)],
  direction = selected$direction,
  selection_rank = selected$selection_rank,
  term_id = selected$term_id,
  Description = selected$Description,
  NES = selected$NES,
  p.adjust = selected$p.adjust,
  setSize = selected$setSize,
  n_leading_edge_genes = selected$n_leading_edge_genes,
  max_jaccard_with_kept = selected$max_jaccard_with_kept,
  stringsAsFactors = FALSE
)
selected_out <- selected_out[order(selected_out$comparison, selected_out$direction != "up",
  selected_out$selection_rank, method = "radix"), , drop = FALSE]

write_csv_atomic(drivers, driver_file)
write_csv_atomic(selected_out, selection_file)

n_significant <- length(unique(best$comparison[!is.na(best$p.adjust) & best$p.adjust < PARAMS$fdr_threshold]))
counts <- list(
  comparisons = length(comparisons),
  comparisons_with_significant_term = n_significant,
  selected_rows = nrow(selected_out),
  top_terms = length(top_terms),
  driver_rows = nrow(drivers),
  driver_genes = length(unique(drivers$Gene)),
  provenance_rows_read = n_provenance_rows,
  claim_safe_core_rows = nrow(core_rows)
)
message(paste(names(counts), unlist(counts), sep = ": ", collapse = "\n"))

write_run_manifest(
  run_manifest_file,
  inputs = list(
    manifest = manifest_file,
    term_comparison = bundle$term_source,
    term_gene_provenance = bundle$provenance_source,
    analysis_status = bundle$status_source,
    config = if (is.na(config_file)) NULL else config_file
  ),
  outputs = list(
    top_term_driver_genes = driver_file,
    selected_top_terms = selection_file
  ),
  parameters = c(
    list(dataset = DATASET, result_type = TOP_TERM_DRIVER_RESULT_TYPE,
         ontology = TOP_TERM_DRIVER_ONTOLOGY, route_category = TOP_TERM_DRIVER_ROUTE),
    PARAMS,
    list(counts = counts)
  ),
  notes = paste(
    "Top-term selection and top-gene driver table of the archived compareGO tail,",
    "reimplemented on the canonical compareGO bundle. Leading edge = claim-safe core",
    "(core_enrichment_member and gene_level_claim_allowed). Freq and Mean_NES count all",
    "comparisons in which the selected term has a finite NES; n_significant_comparisons",
    "counts those with the term at p.adjust < 0.05."
  )
)
message("Top-term drivers written to: ", PATHS$tables)
