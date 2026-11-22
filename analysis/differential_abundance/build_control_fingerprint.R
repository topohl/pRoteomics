#!/usr/bin/env Rscript
# ================================================================
# Script: analysis/differential_abundance/build_control_fingerprint.R
# Stage: enrichment
# Scope: global
# Consumes: required results/differential_abundance/validate_control_spatial_identity/global/tables/anatomical_protein_contrasts.csv; results/differential_abundance/validate_control_spatial_identity/global/tables/kaulich_signature_mapping.csv; config/marker_panels/wgcna_reference_marker_sets.csv; +9 more; optional results/spatial_validation/build_protein_spatial_atlas/global/tables/protein_baseline_spatial_profile.csv; +2 more
# Produces: results/differential_abundance/build_control_fingerprint/global/tables/con_baseline_profile_long.csv; results/differential_abundance/build_control_fingerprint/global/tables/figure2_spatial_fingerprint_scores.csv; results/differential_abundance/build_control_fingerprint/global/tables/figure2_spatial_fingerprint_proteins.csv; +2 more
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Control-only spatial fingerprint of Extended Data Fig. 6c/e (three CON animals). Descriptive; no model, no test.
# ================================================================
#
# Kept for Extended Data Fig. 6c/e (Tobias, 2026-11-22): the panels show only
# the three CON animals, who are the same with eight animals as with nine. Ported
# from figures/spatial_v6_fingerprint_selection.R, retired from Exp9_manuscript
# in 8e22dd9, and the CON baseline of audits/spatial/spatial_v6_baseline_profile.R.
# Output file names are unchanged so the figure readers keep working.
#
#   con_baseline_profile_long.csv  CON mean per protein and spatial
#     unit (level-2, hemispheres combined) for every measured protein, with its
#     within-protein z-score; same code path as the protein spatial atlas.
#   figure2_spatial_fingerprint_scores.csv    form B (ED6c): mean z of the
#     members of each external signature (Kaulich 2025, cell-class panels).
#   figure2_spatial_fingerprint_proteins.csv  form A (ED6e): the top two
#     proteins of the enriched side of each CON-only anatomical contrast,
#     shown across every spatial unit.
#   figure2_spatial_fingerprint_selection.csv where every row comes from.
#
# Changes from the retired producer: gene symbols come from the canonical
# feature table (official_gene_symbol) instead of the WGCNA shortlist, and the
# run stops if a selected protein is missing from the proteins table (H1F3 was;
# see R/spatial/control_fingerprint_utils.R).
#
# Usage: Rscript analysis/differential_abundance/build_control_fingerprint.R [--dry-run]

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "qc_exploration_utils.R"))
source(repo_path("R", "protigy_input_utils.R"))
source(repo_path("R", "spatial_systems_data_utils.R"))
source(repo_path("R", "differential_abundance_paths.R"))
source(repo_path("R", "control_fingerprint_utils.R"))

ANALYSIS_ID <- "build_control_fingerprint"
SCRIPT_ID <- "analysis/differential_abundance/build_control_fingerprint.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
PATHS <- differential_abundance_dirs(ANALYSIS_ID, scope = "global", create = FALSE)
DATASETS <- valid_datasets()

args <- commandArgs(trailingOnly = TRUE)
dataset_arg <- if ("--dataset" %in% args) args[which(args == "--dataset")[1] + 1L] else "global"
if (!tolower(dataset_arg) %in% c("global", "all")) {
  stop("This is a global script; use --dataset global or all.", call. = FALSE)
}

validation_table <- function(name) {
  canonical_result_path("differential_abundance", "validate_control_spatial_identity",
                        "global", "tables", name)
}
inputs <- list(
  anatomical_contrasts = validation_table("anatomical_protein_contrasts.csv"),
  kaulich_mapping = validation_table("kaulich_signature_mapping.csv"),
  reference_panels = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"),
  spatial_order = repo_path("config", "manuscript_spatial_order.yml")
)
atlas_file <- canonical_result_path("spatial_validation", "build_protein_spatial_atlas",
                                    "global", "tables", "protein_baseline_spatial_profile.csv")
outputs <- list(
  baseline = file.path(PATHS$tables, "con_baseline_profile_long.csv"),
  scores = file.path(PATHS$tables, "figure2_spatial_fingerprint_scores.csv"),
  proteins = file.path(PATHS$tables, "figure2_spatial_fingerprint_proteins.csv"),
  selection = file.path(PATHS$tables, "figure2_spatial_fingerprint_selection.csv")
)
run_manifest_file <- file.path(PATHS$manifests, "run_manifest.yml")

# The dry run writes nothing, not even the input-resolution ledger.
if (is_dry_run()) {
  dry_run_line("Script", SCRIPT_ID)
  ok <- TRUE
  for (nm in names(inputs)) {
    hit <- file.exists(inputs[[nm]])
    ok <- ok && hit
    dry_run_line(nm, inputs[[nm]], if (hit) "PASS" else "FAIL")
  }
  for (ds in DATASETS) {
    matrices <- Sys.glob(path_processed("01_preprocessing", "impute",
      sprintf("*_pgmatrix_imputed_%s_*samples_missing70pct.xlsx", ds)))
    metadata <- preprocessing_module_score_metadata(ds)
    ok <- ok && length(matrices) > 0L && file.exists(metadata)
    dry_run_line(paste0("expression_", ds),
                 if (length(matrices)) paste(basename(matrices), collapse = ", ") else "none",
                 if (length(matrices)) "PASS" else "FAIL")
    dry_run_line(paste0("metadata_", ds), metadata, if (file.exists(metadata)) "PASS" else "FAIL")
  }
  dry_run_line("atlas_baseline (optional check)", atlas_file,
               if (file.exists(atlas_file)) "PASS" else "WARN")
  for (nm in names(outputs)) dry_run_line(paste0("would write ", nm), outputs[[nm]])
  dry_run_line("would write run manifest", run_manifest_file)
  quit(status = if (ok) 0L else 1L, save = "no")
}

for (pkg in c("readr", "yaml", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required.", call. = FALSE)
}
missing_inputs <- names(inputs)[!file.exists(unlist(inputs))]
if (length(missing_inputs)) {
  stop("missing_required_input: ", paste(missing_inputs, "=", unlist(inputs[missing_inputs]),
                                         collapse = "; "), call. = FALSE)
}
read_table <- function(path) {
  as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
}
spatial_order <- cfp_spatial_order(inputs$spatial_order)

# ---- CON baseline over the full protein universe ---------------------------

baseline_rows <- list()
symbol_rows <- list()
expression_files <- character()
metadata_files <- character()
for (ds in DATASETS) {
  message("Building the CON baseline for ", ds)
  ds_inputs <- resolve_dataset_inputs(ds, purpose = "wgcna", script = SCRIPT_ID,
                                      stage = "differential_expression_enrichment")
  md <- preprocessing_module_score_metadata(ds)
  if (!file.exists(md)) stop("Canonical merged metadata not found: ", md, call. = FALSE)
  canonical <- qc_load_canonical_expression(ds_inputs$expression_file, md,
                                            dataset = ds, strict = TRUE)
  cfp_require_columns(canonical$feature_table, c("ProteinGroupID", "official_gene_symbol"),
                      paste("Canonical feature table of", ds))
  lv <- sps_build_spatial_levels(ds, canonical = canonical)
  baseline_rows[[ds]] <- cfp_con_profile(lv$level2$mat, lv$level2$meta, ds, spatial_order)
  symbol_rows[[ds]] <- data.frame(
    dataset = ds,
    ProteinGroupID = as.character(canonical$feature_table$ProteinGroupID),
    official_gene_symbol = as.character(canonical$feature_table$official_gene_symbol),
    stringsAsFactors = FALSE)
  expression_files[[ds]] <- ds_inputs$expression_file
  metadata_files[[ds]] <- md
}
base <- do.call(rbind, baseline_rows)
rownames(base) <- NULL
symbols <- do.call(rbind, symbol_rows)
base$GeneSymbol <- symbols$official_gene_symbol[match(
  paste(base$dataset, base$ProteinGroupID), paste(symbols$dataset, symbols$ProteinGroupID))]
base$gene <- cfp_gene_key(base$GeneSymbol)
if (any(base$n_con_animals != 3L)) {
  stop("The control fingerprint expects three CON animals in every dataset.", call. = FALSE)
}

atlas_check <- if (file.exists(atlas_file)) {
  cfp_reconcile_atlas(base, read_table(atlas_file), spatial_order)
} else {
  list(shared_cells = 0L, max_abs_difference = NA_real_)
}

# ---- form B: external signature scores (ED6c) -----------------------------

signatures <- cfp_signature_sets(read_table(inputs$kaulich_mapping),
                                 read_table(inputs$reference_panels))
scores <- cfp_annotate_units(cfp_score_signatures(base, signatures), spatial_order)

# ---- form A: named proteins from the CON-only contrasts (ED6e) ------------

pick <- cfp_select_contrast_proteins(read_table(inputs$anatomical_contrasts))
proteins <- cfp_annotate_units(cfp_fingerprint_proteins(base, pick), spatial_order)
cfp_check_selection_in_proteins(pick, proteins)
selection <- cfp_selection_table(pick, signatures)

# ---- write ------------------------------------------------------------------

base$value_definition <- paste0(
  "con_mean_log2 = mean over CON animals of the canonical bilateral ",
  "animal-level (level-2) imputed log2 abundance in that spatial unit; ",
  "con_z = within-protein standardisation of con_mean_log2 across that ",
  "protein's spatial units. Descriptive only: no model, no test, no p-value.")
base$selection_status <- "phenotype_blind_full_universe"
base$gene_symbol_source <- "official_gene_symbol of the canonical feature table"
base$gene <- NULL

dir_create(PATHS$tables)
dir_create(PATHS$manifests)
write_csv_safe(base, outputs$baseline)
write_csv_safe(scores, outputs$scores)
write_csv_safe(proteins, outputs$proteins)
write_csv_safe(selection, outputs$selection)

counts <- list(
  baseline_rows = nrow(base),
  proteins_per_dataset = as.list(tapply(base$ProteinGroupID, base$dataset,
                                        function(z) length(unique(z)))),
  proteins_without_symbol = length(unique(base$ProteinGroupID[is.na(base$GeneSymbol)])),
  signatures = length(unique(scores$signature)),
  scored_cells = sum(!is.na(scores$score)),
  contrasts = length(unique(pick$contrast)),
  selected_genes = length(unique(pick$gene)),
  atlas_shared_cells = atlas_check$shared_cells,
  atlas_max_abs_difference = atlas_check$max_abs_difference
)
write_run_manifest(
  run_manifest_file,
  inputs = c(inputs, list(atlas_baseline = if (file.exists(atlas_file)) atlas_file else NULL),
             as.list(stats::setNames(expression_files, paste0("expression_", names(expression_files)))),
             as.list(stats::setNames(metadata_files, paste0("metadata_", names(metadata_files))))),
  outputs = outputs,
  parameters = list(
    groups = "CON only",
    min_signature_members = CFP_MIN_SIGNATURE_MEMBERS,
    top_proteins_per_contrast = CFP_TOP_PROTEINS_PER_CONTRAST,
    reference_sets = as.list(CFP_REFERENCE_SETS),
    gene_symbol_source = "official_gene_symbol of the canonical feature table",
    counts = counts
  ),
  notes = paste(
    "Control-only spatial fingerprint for Extended Data Fig. 6c/e, kept 2026-11-22.",
    "Rows come from external signatures, external cell-class panels or the",
    "prespecified CON-only anatomical contrasts; no stress group enters selection",
    "or values. Ported from the retired figures/spatial_v6_fingerprint_selection.R."
  )
)

cat("\n===== control spatial fingerprint =====\n")
cat("baseline rows:", nrow(base), "\n")
cat("signatures:", counts$signatures, " scored cells:", counts$scored_cells, "\n")
cat("selected genes by contrast:\n")
for (cc in sort(unique(pick$contrast))) {
  cat(sprintf("  %-36s %s\n", cc, paste(pick$gene[pick$contrast == cc], collapse = ", ")))
}
if (atlas_check$shared_cells > 0L) {
  cat("atlas check: ", atlas_check$shared_cells, " shared cells, max abs difference ",
      format(atlas_check$max_abs_difference, scientific = TRUE), "\n", sep = "")
} else {
  cat("atlas check: skipped (no atlas baseline yet)\n")
}
message("Control fingerprint written to: ", PATHS$tables)
