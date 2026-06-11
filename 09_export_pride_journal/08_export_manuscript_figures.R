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

args <- commandArgs(trailingOnly = TRUE)
dry_run <- is_dry_run()

manuscript_dirs <- path_results("manuscript", c("figure_1", "figure_2", "extended_data"))
invisible(lapply(manuscript_dirs, dir_create))

candidate_roots <- c(
  path_results("figures", "03_qc_exploration"),
  path_results("figures", "04_differential_expression_enrichment"),
  path_results("figures", "06_modules_WGCNA"),
  path_results("figures", "07_spatial_networks"),
  path_results("figures", "08_behavior_physio_coupling")
)

if (isTRUE(dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/08_export_manuscript_figures.R")
  dry_run_line("Candidate figure roots", paste(candidate_roots, collapse = "; "))
  dry_run_line("Output root", path_results("manuscript"))
  quit(status = 0, save = "no")
}

candidates <- unlist(lapply(candidate_roots[dir.exists(candidate_roots)], list.files, pattern = "\\.(svg|pdf|png)$", recursive = TRUE, full.names = TRUE), use.names = FALSE)

figure_target_name <- function(path) {
  rel <- relative_to(path, path_results("figures"))
  stem <- tools::file_path_sans_ext(rel)
  paste0(safe_filename(stem), ".", tolower(tools::file_ext(path)))
}

figure_audit_row <- function(path, target) {
  ext <- tolower(tools::file_ext(path))
  info <- file.info(path)
  sibling_source <- file.path(path_results("source_data"), sub("\\.(svg|pdf|png)$", ".csv", relative_to(path, path_results("figures")), ignore.case = TRUE))
  data.frame(
    source_file = path,
    target_file = target,
    file_ext = ext,
    file_size_bytes = if (nrow(info)) as.numeric(info$size) else NA_real_,
    editable_vector = ext %in% c("svg", "pdf"),
    has_png_companion = file.exists(sub("\\.[^.]+$", ".png", path)),
    source_data_candidate = sibling_source,
    source_data_candidate_exists = file.exists(sibling_source),
    publication_readiness_note = if (ext == "svg") {
      "editable_vector_check_fonts_labels_and_panel_size"
    } else if (ext == "pdf") {
      "vector_or_pdf_check_fonts_labels_and_panel_size"
    } else {
      "raster_check_resolution_and_source_data"
    },
    stringsAsFactors = FALSE
  )
}

manifest <- data.frame(
  source_file = candidates,
  target_file = file.path(path_results("manuscript", "extended_data"), vapply(candidates, figure_target_name, character(1))),
  stringsAsFactors = FALSE
)
if (nrow(manifest)) {
  file.copy(manifest$source_file, manifest$target_file, overwrite = TRUE)
}

figure_audit <- if (nrow(manifest)) {
  do.call(rbind, Map(figure_audit_row, manifest$source_file, manifest$target_file))
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
