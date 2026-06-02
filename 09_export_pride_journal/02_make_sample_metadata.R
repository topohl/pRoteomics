#!/usr/bin/env Rscript
# Creates SDRF-like sample metadata from pg_matrix-era metadata (not raw-MS acquisition design).
# Fails clearly when mandatory fields cannot be resolved.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "export_helpers.R"))

ensure_pride_dirs()
config <- load_export_config()
cli <- export_cli_args()
datasets <- resolve_export_datasets(cli$dataset, config)

read_metadata_workbook <- function() {
  meta_paths <- metadata_input_paths(config)
  xlsx <- meta_paths[grepl("\\.xlsx$", meta_paths, ignore.case = TRUE)]
  if (length(xlsx) && requireNamespace("readxl", quietly = TRUE)) {
    return(as.data.frame(readxl::read_excel(xlsx[[1]])))
  }
  meta <- read_sample_metadata()
  if (nrow(meta)) return(as.data.frame(meta))
  data.frame()
}

standardize_meta_names <- function(df) {
  names(df) <- standardize_names(names(df))
  df
}

meta_raw <- read_metadata_workbook()
if (!nrow(meta_raw)) {
  if (isTRUE(cli$dry_run)) {
    dry_run_line("Script", "09_export_pride_journal/02_make_sample_metadata.R")
    dry_run_line("Sample metadata rows", 0, "FAIL")
    quit(status = 1, save = "no")
  }
  stop(
    "No sample metadata found. Expected one of: ",
    paste(metadata_input_paths(config), collapse = ", "),
    call. = FALSE
  )
}

meta <- standardize_meta_names(meta_raw)

resolve_col <- function(df, candidates, required = TRUE, label = paste(candidates, collapse = "/")) {
  hit <- intersect(standardize_names(candidates), names(df))
  if (!length(hit)) {
    if (isTRUE(required)) {
      stop("Mandatory metadata column missing: ", label, call. = FALSE)
    }
    return(rep(NA_character_, nrow(df)))
  }
  as.character(df[[hit[[1]]]])
}

if (!"dataset" %in% names(meta)) {
  if ("celltype_layer" %in% names(meta) || "celltype" %in% names(meta)) {
    ct <- resolve_col(meta, c("celltype_layer", "celltype"), required = FALSE)
    meta$dataset <- ifelse(
      grepl("neuropil", ct, ignore.case = TRUE), "neuron_neuropil",
      ifelse(grepl("soma", ct, ignore.case = TRUE), "neuron_soma",
        ifelse(grepl("microgl", ct, ignore.case = TRUE), "microglia", NA_character_)
      )
    )
  }
}
if (!"dataset" %in% names(meta) || all(is.na(meta$dataset))) {
  if (length(datasets) == 1L) {
    meta$dataset <- datasets[[1]]
  } else {
    stop("Cannot resolve dataset for sample metadata. Add a dataset column or run with --dataset <one>.", call. = FALSE)
  }
}

meta$dataset <- normalize_dataset(meta$dataset)
meta <- meta[meta$dataset %in% datasets, , drop = FALSE]
if (!nrow(meta)) {
  stop("No metadata rows for requested dataset(s): ", paste(datasets, collapse = ", "), call. = FALSE)
}

sample_id <- resolve_col(meta, c("sample_id", "sample", "samplecolumn"))
if (any(is.na(sample_id) | !nzchar(sample_id))) {
  stop("Mandatory sample_id values are missing or empty.", call. = FALSE)
}

build_export_table <- function(meta, dataset) {
  sub <- meta[meta$dataset == dataset, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  interpretation <- config$interpretation_scope[[dataset]] %||% dataset_interpretation(dataset)
  sid <- resolve_col(sub, c("sample_id", "sample"))
  pg_col <- resolve_col(sub, c("shortmicrogliame", "sample_number", "samplecolumn", "sample_id"), required = FALSE)
  pg_col <- ifelse(is.na(pg_col) | !nzchar(pg_col), sid, pg_col)
  raw_name <- resolve_col(sub, c("raw_file_name", "rawfile", "file", "filename"), required = FALSE)
  if (all(is.na(raw_name) | !nzchar(raw_name))) {
    raw_name <- ifelse(grepl("\\.(d|raw|wiff|mzml)$", sid, ignore.case = TRUE), basename(sid), NA_character_)
  }
  data.frame(
    dataset = dataset,
    sample_id = sid,
    animal_id = resolve_col(sub, c("animal_id", "mouse_id", "subject_id", "animalid"), required = FALSE),
    group = resolve_col(sub, c("group", "condition", "phenotype", "treatment"), required = FALSE),
    condition = resolve_col(sub, c("condition", "group", "phenotype"), required = FALSE),
    sex = resolve_col(sub, c("sex"), required = FALSE),
    batch = resolve_col(sub, c("batch", "run_batch", "ms_batch"), required = FALSE),
    region = resolve_col(sub, c("region"), required = FALSE),
    layer = resolve_col(sub, c("layer"), required = FALSE),
    roi_cell_compartment = resolve_col(sub, c("celltype_roi", "celltype_layer", "celltype", "roi"), required = FALSE),
    pg_matrix_source_column = pg_col,
    biological_replicate = resolve_col(sub, c("biological_replicate", "bio_rep", "replicategroup", "animal_id", "mouse_id", "animalid"), required = FALSE),
    technical_replicate = resolve_col(sub, c("technical_replicate", "tech_rep", "side", "replicate_side", "sample_number"), required = FALSE),
    interpretation_scope = interpretation,
    raw_file_name = raw_name,
    stringsAsFactors = FALSE
  )
}

export_meta <- do.call(rbind, lapply(datasets, function(ds) build_export_table(meta, ds)))

sdrf <- data.frame(
  source_name = export_meta$sample_id,
  characteristics_organism = "Mus musculus",
  characteristics_sex = export_meta$sex,
  characteristics_phenotype = export_meta$condition,
  characteristics_region = export_meta$region,
  characteristics_layer = export_meta$layer,
  characteristics_cell_compartment = export_meta$roi_cell_compartment,
  characteristics_dataset = export_meta$dataset,
  characteristics_interpretation_scope = export_meta$interpretation_scope,
  assay_name = export_meta$pg_matrix_source_column,
  comment_pg_matrix_column = export_meta$pg_matrix_source_column,
  comment_data_file = export_meta$raw_file_name,
  stringsAsFactors = FALSE
)

out_clean <- pride_submission_dir("metadata", "sample_metadata.tsv")
out_sdrf <- pride_submission_dir("metadata", "sdrf_like_metadata.tsv")

if (isTRUE(cli$dry_run)) {
  dry_run_line("Script", "09_export_pride_journal/02_make_sample_metadata.R")
  dry_run_line("Datasets", paste(datasets, collapse = ", "))
  dry_run_line("Sample metadata rows", nrow(export_meta), if (nrow(export_meta) > 0) "PASS" else "FAIL")
  dry_run_line("Metadata target", out_clean)
  dry_run_line("SDRF-like target", out_sdrf)
  quit(status = if (nrow(export_meta) > 0) 0 else 1, save = "no")
}

write_tsv(export_meta, out_clean)
write_tsv(sdrf, out_sdrf)
message("Wrote sample metadata: ", out_clean)
message("Wrote SDRF-like metadata (pg_matrix-era): ", out_sdrf)
