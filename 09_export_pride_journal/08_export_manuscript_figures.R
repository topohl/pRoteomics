#!/usr/bin/env Rscript

# ================================================================
# Script: 09_export_pride_journal/08_export_manuscript_figures.R
# Stage: export
# Scope: global
# Consumes: required results/figures/; optional none.
# Produces: results/manuscript/figure_1/; results/manuscript/figure_2/; results/manuscript/extended_data/; +1 more.
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Manuscript figure export.
# ================================================================

# Collect final manuscript figures from module outputs (journal output; not PRIDE-required).

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "validation_utils.R"))
source(repo_path("R", "export_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
dry_run <- is_dry_run()

manuscript_dirs <- path_results("manuscript", c("figure_1", "figure_2", "extended_data"))

# Explicit root list, never an unrestricted recursive results/figures scan: the
# EWCE root is deliberately the canonical branch only, and an open scan would
# also pull comparison/repro trees. manuscript_panels holds the Figure-3 panels
# and was previously omitted, so no Figure-3 panel could ever be exported.
candidate_roots <- c(
  path_results("figures", "03_qc_exploration"),
  path_results("figures", "04_differential_expression_enrichment"),
  canonical_ewce_figure_root(),
  path_results("figures", "06_modules_WGCNA"),
  path_results("figures", "07_spatial_networks"),
  path_results("figures", "08_behavior_physio_coupling"),
  path_results("figures", "08_biological_interpretation"),
  path_results("figures", "10_biological_integration"),
  path_results("figures", "manuscript_panels")
)

# Selection is read-only, so it is computed before the dry-run guard and
# reported by it. Nothing above this point may mutate the filesystem:
# --dry-run must be side-effect free.
candidates <- unlist(lapply(candidate_roots[dir.exists(candidate_roots)], list.files, pattern = "\\.(svg|pdf|png)$", recursive = TRUE, full.names = TRUE), use.names = FALSE)
candidates <- candidates[is_exportable_result_path(candidates)]
candidates <- drop_legacy_dataset_suffixed_aliases(candidates)
candidates <- drop_orphan_figure_families(candidates)

if (isTRUE(dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/08_export_manuscript_figures.R")
  dry_run_line("Candidate figure roots", paste(candidate_roots, collapse = "; "))
  dry_run_line("Output root", path_results("manuscript"))
  dry_run_line("Selected files", length(candidates))
  quit(status = 0, save = "no")
}

invisible(lapply(manuscript_dirs, dir_create))

figures_root <- path_results("figures")
candidate_rel <- relative_to(candidates, figures_root)

# Built column-wise. relative_to(), file.info() and file.exists() are all
# vectorised, so this issues a handful of filesystem round trips instead of
# several per file, which matters on a network share.
figure_audit_table <- function(paths, targets, rel_paths) {
  ext <- tolower(tools::file_ext(paths))
  sibling_source <- file.path(
    path_results("source_data"),
    sub("\\.(svg|pdf|png)$", ".csv", rel_paths, ignore.case = TRUE)
  )
  data.frame(
    source_file = paths,
    target_file = targets,
    file_ext = ext,
    file_size_bytes = as.numeric(file.info(paths)$size),
    editable_vector = ext %in% c("svg", "pdf"),
    has_png_companion = file.exists(sub("\\.[^.]+$", ".png", paths)),
    source_data_candidate = sibling_source,
    source_data_candidate_exists = file.exists(sibling_source),
    publication_readiness_note = ifelse(
      ext == "svg", "editable_vector_check_fonts_labels_and_panel_size",
      ifelse(ext == "pdf", "vector_or_pdf_check_fonts_labels_and_panel_size",
             "raster_check_resolution_and_source_data")
    ),
    stringsAsFactors = FALSE
  )
}

manifest <- data.frame(
  source_file = candidates,
  target_file = manuscript_figure_target_paths(
    candidate_rel, path_results("manuscript", "extended_data")
  ),
  stringsAsFactors = FALSE
)
# Fail closed before touching the payload, then copy and verify. Nothing below
# runs unless every single file landed, so the manifests and the run manifest
# can never describe a partial export.
assert_unique_export_targets(manifest$target_file, manifest$source_file)
over_budget <- manifest$target_file[nchar(manifest$target_file) > manuscript_figure_path_budget()]
if (length(over_budget)) {
  stop(
    length(over_budget), " target path(s) exceed the ",
    manuscript_figure_path_budget(), "-character budget; first: ",
    basename(over_budget[[1]]), call. = FALSE
  )
}
copy_export_targets(manifest$source_file, manifest$target_file)

figure_audit <- if (nrow(manifest)) {
  figure_audit_table(manifest$source_file, manifest$target_file, candidate_rel)
} else {
  data.frame(
    source_file = character(), target_file = character(), file_ext = character(),
    file_size_bytes = numeric(), editable_vector = logical(), has_png_companion = logical(),
    source_data_candidate = character(), source_data_candidate_exists = logical(),
    publication_readiness_note = character(), stringsAsFactors = FALSE
  )
}

manifest_path <- path_results("manuscript", "figure_export_manifest.csv")
audit_path <- path_results("manuscript", "figure_publication_audit.csv")
utils::write.csv(manifest, manifest_path, row.names = FALSE)
utils::write.csv(figure_audit, audit_path, row.names = FALSE)
write_run_manifest(
  path_results("logs", "09_export_pride_journal", "manuscript_figures", "run_manifest.yml"),
  inputs = list(figures = candidates),
  outputs = list(manifest = manifest_path, figure_audit = audit_path, manuscript_root = path_results("manuscript")),
  notes = "Collect-only manuscript figure export; no analyses are recomputed. Target filenames preserve source-relative context to avoid basename collisions."
)
message("Manuscript figure export manifest written: ", manifest_path)
