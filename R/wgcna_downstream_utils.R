# Shared helpers for WGCNA downstream interpretation scripts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "module_contracts.R"))

WGCNA_ROI_NOTE <- "microglia-enriched ROI/local microenvironment; annotation only, not purity correction."

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

wgcna_cli <- function(default_dataset = "neuron_neuropil", allow_all = FALSE) {
  args <- commandArgs(trailingOnly = TRUE)
  value_after <- function(flag, default = "") {
    hit <- which(args == flag)
    if (!length(hit) || hit[[1]] == length(args)) return(default)
    args[[hit[[1]] + 1L]]
  }
  dataset <- value_after("--dataset", Sys.getenv("PROTEOMICS_DATASET", unset = default_dataset))
  if (allow_all && identical(tolower(dataset), "all")) {
    dataset <- "all"
  } else {
    dataset <- validate_dataset(dataset, source = "--dataset")
    Sys.setenv(PROTEOMICS_DATASET = dataset)
  }
  list(
    args = args,
    dataset = dataset,
    level = value_after("--level", "both"),
    dry_run = is_dry_run()
  )
}

wgcna_downstream_paths <- function(substep, dataset) {
  create_module_dirs("06_modules_WGCNA", file.path(substep, dataset))
}

safe_read_csv <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NULL)
  if (requireNamespace("readr", quietly = TRUE)) {
    return(tryCatch(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL))
  }
  tryCatch(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}

write_csv_safe2 <- function(x, path) {
  dir_create(dirname(path))
  if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(x, path, na = "") else utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

write_table_and_source <- function(x, table_dir, source_dir, filename) {
  out <- list(
    table = write_csv_safe2(x, file.path(table_dir, filename)),
    source = write_csv_safe2(x, file.path(source_dir, filename))
  )
  invisible(out)
}

first_present_col <- function(df, candidates) {
  if (is.null(df)) return(NA_character_)
  norm <- function(x) tolower(gsub("[^a-z0-9]", "", x))
  hit <- match(norm(candidates), norm(names(df)))
  hit <- hit[!is.na(hit)]
  if (!length(hit)) return(NA_character_)
  names(df)[hit[[1]]]
}

normalize_gene_token <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("_MOUSE$", "", x, ignore.case = TRUE)
  gsub("[^A-Z0-9]", "", x)
}

split_tokens <- function(x) {
  x <- as.character(x)
  out <- unlist(strsplit(x, "[/;,|[:space:]]+"), use.names = FALSE)
  out <- toupper(trimws(out))
  unique(out[nzchar(out) & !is.na(out)])
}

wgcna_marker_sets <- function() {
  list(
    microglia = c("Aif1", "Tmem119", "P2ry12", "Cx3cr1", "Csf1r", "C1qa", "C1qb", "C1qc", "Hexb", "Tyrobp", "Trem2", "Mertk", "Itgam", "Spi1"),
    neuronal_synaptic_neuropil = c("Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2", "Snap25", "Syn1", "Syp", "Dlg4", "Camk2a", "Camk2b", "Map2", "Nefl", "Nefm", "Vamp2"),
    neuropil_synaptic_neuronal = c("Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2", "Snap25", "Syn1", "Syp", "Dlg4", "Camk2a", "Camk2b", "Map2", "Nefl", "Nefm", "Vamp2"),
    nuclear_soma = c("H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3", "Matr3", "Srsf3", "Ddx39b"),
    astrocyte = c("Gfap", "Aqp4", "Aldh1l1", "Slc1a2", "Slc1a3", "Aldoc", "Glul", "Gja1"),
    oligodendrocyte_myelin = c("Mbp", "Mog", "Plp1", "Cnp", "Mag", "Mobp", "Cldn11"),
    endothelial_pericyte_vascular = c("Pecam1", "Cldn5", "Kdr", "Flt1", "Rgs5", "Pdgfrb", "Vtn", "Acta2"),
    mitochondrial_oxphos = c("Ndufs1", "Ndufa9", "Sdha", "Uqcrc2", "Cox4i1", "Atp5f1a", "Atp5f1b"),
    ribosomal_translation = c("Rpl3", "Rpl4", "Rpl5", "Rps3", "Rps6", "Eef1a1", "Eef2"),
    rnp_rna_processing = c("Hnrnpa2b1", "Hnrnpc", "Sfpq", "Snrnp70", "Ddx5", "Ddx17", "Pabpc1")
  )
}

standardize_wgcna_metadata <- function(meta, dataset) {
  meta <- as.data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE)
  sample_col <- first_present_col(meta, c("Sample", "sample", "SampleID", "sample_id", "row.names"))
  if (is.na(sample_col)) meta$Sample <- rownames(meta) else meta$Sample <- as.character(meta[[sample_col]])
  for (target in c("Region", "Layer", "Sex", "Batch")) {
    col <- first_present_col(meta, c(target, tolower(target), if (target == "Batch") c("plate", "run", "batch_id") else character()))
    meta[[target]] <- if (!is.na(col)) as.character(meta[[col]]) else NA_character_
  }
  group_col <- first_present_col(meta, c("StressGroup", "ExpGroup", "condition", "Group", "group", "group2"))
  meta$StressGroup <- if (!is.na(group_col)) toupper(as.character(meta[[group_col]])) else NA_character_
  meta$StressGroup <- dplyr::case_when(
    grepl("^CON|^CTRL|CONTROL", meta$StressGroup, ignore.case = TRUE) ~ "CON",
    grepl("^RES", meta$StressGroup, ignore.case = TRUE) ~ "RES",
    grepl("^SUS", meta$StressGroup, ignore.case = TRUE) ~ "SUS",
    TRUE ~ meta$StressGroup
  )
  if (!"ExpGroup" %in% names(meta)) meta$ExpGroup <- meta$StressGroup
  meta$RegionLayer <- ifelse(!is.na(meta$Region) & nzchar(meta$Region) & !is.na(meta$Layer) & nzchar(meta$Layer), paste(meta$Region, meta$Layer, sep = "_"), NA_character_)
  spatial_col <- if (dataset == "neuron_neuropil" && any(!is.na(meta$RegionLayer))) "RegionLayer" else "Region"
  meta$SpatialUnit <- spatial_col
  meta$SpatialLabel <- as.character(meta[[spatial_col]])
  meta
}

resolve_wgcna_files <- function(dataset) {
  list(
    state = path_processed("06_modules_WGCNA", "01_WGCNA", dataset, "wgcna_final_model_state.rds"),
    definitions = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_definitions_for_downstream.csv"),
    module_summary = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_summary.csv"),
    go = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules", "WGCNA_module_GO_enrichment_long.csv"),
    supermodule_annotation = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "supermodules", "wgcna_module_supermodule_annotation.csv"),
    supermodule_summary = path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "supermodules", "wgcna_supermodule_summary.csv"),
    marker_traits = path_results("tables", "03_qc_exploration", "06_wgcna_marker_trait_export", dataset, "wgcna_marker_traits_by_sample.csv"),
    neuropil_annotation = path_results("tables", "04_differential_expression_enrichment", "neuropil_reference_annotation", "microglia", "microglia_neuropil_annotation_latest.csv")
  )
}

load_wgcna_state <- function(path) {
  if (!file.exists(path)) stop("Missing WGCNA final state: ", path, call. = FALSE)
  readRDS(path)
}

extract_module_eigengenes <- function(state) {
  MEs <- state$mergedMEs %||% state$MEs %||% state$moduleEigengenes
  if (is.null(MEs)) stop("WGCNA state does not contain mergedMEs/module eigengenes.", call. = FALSE)
  MEs <- as.data.frame(MEs, check.names = FALSE)
  if (!"Sample" %in% names(MEs)) MEs <- tibble::rownames_to_column(MEs, "Sample")
  MEs
}

module_col_to_id <- function(x) {
  out <- sub("^ME", "", as.character(x))
  out
}

make_supermodule_eigengenes <- function(module_eigengenes, super_map) {
  me_cols <- setdiff(names(module_eigengenes), "Sample")
  super_map <- super_map[super_map$module_eigengene %in% me_cols & !is.na(super_map$SupermoduleID), , drop = FALSE]
  rows <- list(Sample = module_eigengenes$Sample)
  comp <- list()
  for (sid in unique(super_map$SupermoduleID)) {
    members <- unique(super_map$module_eigengene[super_map$SupermoduleID == sid])
    vals <- module_eigengenes[, members, drop = FALSE]
    vals <- vals[, vapply(vals, function(z) stats::var(as.numeric(z), na.rm = TRUE) > 0, logical(1)), drop = FALSE]
    if (!ncol(vals)) next
    if (ncol(vals) == 1L) {
      score <- as.numeric(vals[[1]])
    } else {
      pc <- stats::prcomp(vals, center = TRUE, scale. = TRUE)$x[, 1L]
      mean_vec <- rowMeans(vals, na.rm = TRUE)
      score <- if (stats::cor(pc, mean_vec, use = "pairwise.complete.obs") < 0) -pc else pc
    }
    col <- paste0("SM__", safe_filename(sid))
    rows[[col]] <- as.numeric(score)
    comp[[length(comp) + 1L]] <- data.frame(
      supermodule_id = sid,
      supermodule_eigengene = col,
      n_member_modules = length(members),
      member_modules = paste(members, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  list(eigengenes = as.data.frame(rows, check.names = FALSE), composition = dplyr::bind_rows(comp))
}

empty_group_effects <- function(dataset, level, reason) {
  data.frame(
    dataset = dataset, level = level, module_id = NA_character_, supermodule_id = NA_character_,
    module_label = NA_character_, supermodule_label = NA_character_, spatial_unit = NA_character_,
    contrast = NA_character_, estimate = NA_real_, SE = NA_real_, statistic = NA_real_,
    p_value = NA_real_, FDR_within_dataset_level = NA_real_, FDR_global = NA_real_,
    direction = NA_character_, n_samples = 0L, formula_used = NA_character_,
    dropped_covariates = NA_character_, rank_deficient_model = NA, model_warning = reason,
    stringsAsFactors = FALSE
  )
}

required_group_effect_columns <- c(
  "dataset", "level", "module_id", "supermodule_id", "module_label", "supermodule_label",
  "spatial_unit", "contrast", "estimate", "SE", "statistic", "p_value",
  "FDR_within_dataset_level", "FDR_global", "direction", "n_samples",
  "formula_used", "dropped_covariates", "rank_deficient_model", "model_warning"
)

required_module_annotation_columns <- c(
  "dataset", "ModuleID", "ModuleColor", "n_proteins", "microenvironment_class", "interpretation_note"
)

required_interpretable_columns <- c(
  "dataset", "level", "contrast", "estimate", "p_value", "FDR_global", "interpretation_sentence"
)
