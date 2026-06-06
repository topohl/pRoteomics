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

shorten_supermodule_label <- function(x, max_chars = 45) {
  vapply(as.character(x), function(z) {
    z <- trimws(z)
    if (is.na(z) || !nzchar(z)) return("Unresolved / mixed")
    z <- gsub("\\s+([0-9]+)\\s+[Mm]odules?$", "", z)
    z <- gsub("\\b[Mm]odules?\\b", "", z)
    z <- gsub("\\s*\\([^)]*modules?[^)]*\\)", "", z, ignore.case = TRUE)
    z <- gsub("\\b(process|regulation|pathway|of|the|cellular|biological|positive|negative)\\b", "", z, ignore.case = TRUE)
    z <- gsub("\\s+", " ", z)
    parts <- trimws(unlist(strsplit(z, "\\s*[;/|]\\s*", perl = TRUE), use.names = FALSE))
    parts <- parts[nzchar(parts)]
    if (length(parts)) z <- paste(utils::head(parts, 2), collapse = " / ")
    z <- stringr::str_squish(z)
    if (!nzchar(z)) z <- "Unresolved / mixed"
    if (nchar(z) > max_chars) {
      words <- unlist(strsplit(z, "\\s+"), use.names = FALSE)
      keep <- character()
      for (word in words) {
        cand <- paste(c(keep, word), collapse = " ")
        if (nchar(cand) > max_chars) break
        keep <- c(keep, word)
      }
      z <- if (length(keep)) paste(keep, collapse = " ") else substr(z, 1, max_chars)
    }
    z
  }, character(1))
}

macroprogram_display <- function(x) {
  vapply(as.character(x), function(z) {
    z0 <- tolower(trimws(z %||% ""))
    if (!nzchar(z0)) return("Unresolved / mixed")
    if (grepl("ecm|adhesion|basement membrane|collagen|laminin|integrin", z0)) return("Perivascular ECM / adhesion")
    if (grepl("mitochondr|respiratory|oxidative|\\batp\\b|\\btca\\b|acetyl-coa", z0)) return("Mitochondrial metabolism")
    if (grepl("\\brna\\b|ribosome|translation|splice|\\brnp\\b|ncrna", z0)) return("RNA / translation")
    if (grepl("synapse|vesicle|postsynaptic|actin|cytoskeleton", z0)) return("Synaptic / cytoskeletal")
    if (grepl("microglia|phagolysosomal|immune", z0)) return("Microglia state")
    if (grepl("neuropil|neuronal", z0)) return("Neuropil / neuronal")
    if (grepl("\\bbbb\\b|endothelial|pericyte|vascular", z0)) return("Vascular / BBB")
    "Unresolved / mixed"
  }, character(1))
}

supermodule_microenvironment_label <- function(cls, dataset = current_dataset()) {
  vapply(as.character(cls), function(z) {
    z <- tolower(trimws(z %||% ""))
    if (!nzchar(z)) return(NA_character_)
    if (z %in% c("vascular_basement_membrane_ecm", "vascular/ecm")) return("Perivascular ECM")
    if (z %in% c("vascular_bbb_mural", "vascular")) return("Vascular / BBB")
    if (z %in% c("neuropil_sensitive", "neuropil")) return("Neuropil reference overlap")
    if (z %in% c("astrocyte_or_endfoot_sensitive", "astrocyte")) return("Astrocyte / endfoot")
    if (z %in% c("oligodendrocyte_or_myelin_sensitive", "oligodendrocyte/myelin")) return("Oligodendrocyte / myelin")
    if (z %in% c("ambiguous_or_mixed", "shared_microenvironment", "mixed")) return("Mixed / unresolved")
    if (identical(as.character(dataset), "microglia") && z %in% c("microglia_supported", "microglia_state_or_activation_supported")) return("Microglia-associated ROI")
    NA_character_
  }, character(1))
}

compose_supermodule_display_label <- function(supermodule_id, short_label) {
  id <- as.character(supermodule_id)
  id[is.na(id) | !nzchar(id)] <- "SM??"
  label <- shorten_supermodule_label(short_label, max_chars = 30)
  label[is.na(label) | !nzchar(label)] <- "Mixed / unresolved"
  paste0(id, " \u00b7 ", label)
}

classify_supermodule_label_confidence <- function(n_modules, go_class = "unresolved",
                                                  has_coherent_hubs = FALSE,
                                                  microenvironment_class = NA_character_,
                                                  high_unmapped_fraction = FALSE) {
  n_modules <- suppressWarnings(as.integer(n_modules))
  singleton <- is.na(n_modules) | n_modules <= 1L
  go_class <- as.character(go_class %||% "unresolved")
  micro_label <- supermodule_microenvironment_label(microenvironment_class)
  micro_supported <- !is.na(micro_label) & nzchar(micro_label) & !micro_label %in% "Mixed / unresolved"
  mixed_micro <- !is.na(micro_label) & micro_label == "Mixed / unresolved"
  go_supported <- go_class %in% c("GO_supported", "data_driven_GO")
  suggestive_go <- go_class %in% c("suggestive_GO", "manual_only")
  if (singleton || isTRUE(high_unmapped_fraction) || mixed_micro) return("low")
  if ((go_supported || micro_supported) && isTRUE(has_coherent_hubs)) return("high")
  if ((go_supported && micro_supported) || (suggestive_go && isTRUE(has_coherent_hubs)) || (micro_supported && isTRUE(has_coherent_hubs))) return("medium")
  if (isTRUE(has_coherent_hubs) || suggestive_go || go_supported || micro_supported) return("low")
  "unresolved"
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
    microglia = c("Aif1", "Tmem119", "P2ry12", "Cx3cr1", "Csf1r", "Hexb",
  "Fcrls", "Olfml3", "Sall1", "Siglech", "Gpr34", "P2ry13",
  "Tgfbr1", "Mertk", "Spi1", "C1qa", "C1qb", "C1qc", "Tyrobp", "Trem2", "Apoe",
  "Lpl", "Ctsb", "Ctsd", "Ctsz", "Lgals3", "Itgam",
  "Cd68", "Fcgr3", "Clec7a"),
    neuronal_synaptic_neuropil = c("Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2",
  "Snap25", "Snap91", "Syn1", "Syn2", "Syp", "Syt1",
  "Vamp2", "Vamp1", "Dlg4", "Dlg3", "Shank1", "Shank2",
  "Homer1", "Camk2a", "Camk2b", "Gria1", "Gria2",
  "Grin1", "Grin2a", "Map2", "Tubb3", "Nefl", "Nefm",
  "Nefh", "Rbfox3"),
    neuropil_synaptic_neuronal = c("Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2",
  "Snap25", "Snap91", "Syn1", "Syn2", "Syp", "Syt1",
  "Vamp2", "Vamp1", "Dlg4", "Dlg3", "Shank1", "Shank2",
  "Homer1", "Camk2a", "Camk2b", "Gria1", "Gria2",
  "Grin1", "Grin2a", "Map2", "Tubb3", "Nefl", "Nefm",
  "Nefh", "Rbfox3"),
    nuclear_soma = c("H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3",
  "Matr3", "Srsf3", "Ddx39b", "Lmna", "Lmnb1", "Lmnb2",
  "Ncl", "Npm1", "Nono", "Sfpq", "Hnrnpm", "Hnrnpu",
  "Rbm39", "Top2b", "Hist1h1c"),
    astrocyte = c("Gfap", "Aqp4", "Aldh1l1", "Slc1a2", "Slc1a3",
  "Aldoc", "Glul", "Gja1", "S100b", "Fabp7",
  "Clu", "Sparcl1", "Fgfr3", "Pla2g7", "Mlc1"),
    oligodendrocyte_myelin = c("Mbp", "Mog", "Plp1", "Cnp", "Mag", "Mobp", "Cldn11",
  "Myrf", "Olig1", "Olig2", "Sox10", "Mal", "Fa2h",
  "Ugt8a", "Opalin", "Ermn", "Tspan2"),
    endothelial_pericyte_vascular = c("Pecam1", "Cldn5", "Kdr", "Flt1", "Tek", "Klf2",
  "Klf4", "Slco1a4", "Abcb1a", "Plvap", "Vwf", "Emcn"),
    mitochondrial_oxphos = c("Ndufs1", "Ndufs2", "Ndufv1", "Ndufa9", "Ndufb8",
  "Sdha", "Sdhb", "Uqcrc1", "Uqcrc2", "Cyc1",
  "Cox4i1", "Cox5a", "Cox6c", "Atp5f1a", "Atp5f1b",
  "Atp5f1c", "Atp5mc1", "Atp5mc2", "Atp5pb"),
    ribosomal_translation = c("Rpl3", "Rpl4", "Rpl5", "Rpl7", "Rpl10", "Rpl13",
  "Rps3", "Rps6", "Rps8", "Rps14", "Rps18",
  "Eef1a1", "Eef1a2", "Eef2", "Eif3a", "Eif4a1",
  "Eif4g1", "Eif5a"),
    rnp_rna_processing = c("Hnrnpa2b1", "Hnrnpc", "Hnrnpm", "Hnrnpu",
  "Sfpq", "Nono", "Snrnp70", "Snrpa", "Snrpb",
  "Ddx5", "Ddx17", "Ddx39b", "Pabpc1", "Pabpc4",
  "Rbm39", "Srsf1", "Srsf3", "Srsf7")
  )
}

wgcna_registry_required_columns <- c(
  "marker_set", "cell_class", "cell_state", "gene_symbol", "source_type",
  "source_name", "source_reference", "selection_rule", "confidence", "use_for", "notes"
)

read_wgcna_marker_registry <- function(path = Sys.getenv("PROTEOMICS_WGCNA_MARKER_REGISTRY_FILE", unset = "")) {
  if (!nzchar(path)) path <- repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  registry <- safe_read_csv(path)
  if (is.null(registry) || !nrow(registry)) return(NULL)
  missing <- setdiff(wgcna_registry_required_columns, names(registry))
  if (length(missing)) stop("Marker registry is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  registry$gene_symbol <- as.character(registry$gene_symbol)
  registry$gene_token <- normalize_gene_token(registry$gene_symbol)
  registry <- registry[nzchar(registry$gene_token) & !is.na(registry$gene_token), , drop = FALSE]
  attr(registry, "marker_registry_file") <- path
  attr(registry, "marker_registry_version") <- paste(unique(registry$source_name), collapse = ";")
  registry
}

marker_registry_to_sets <- function(registry) {
  if (is.null(registry) || !nrow(registry)) return(list())
  split(as.character(registry$gene_symbol), as.character(registry$marker_set)) |>
    lapply(function(x) unique(x[nzchar(normalize_gene_token(x))]))
}

read_empirical_roi_marker_sets <- function(path = Sys.getenv("PROTEOMICS_WGCNA_EMPIRICAL_MARKER_FILE", unset = "")) {
  if (!nzchar(path)) path <- path_results("tables", "03_qc_exploration", "05_empirical_roi_marker_discovery", "empirical_roi_marker_sets.csv")
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  empirical <- safe_read_csv(path)
  if (is.null(empirical) || !nrow(empirical)) return(NULL)
  required <- c("marker_set", "GeneSymbol")
  missing <- setdiff(required, names(empirical))
  if (length(missing)) stop("Empirical marker file is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  empirical$gene_symbol <- as.character(empirical$GeneSymbol)
  empirical$gene_token <- normalize_gene_token(empirical$gene_symbol)
  empirical <- empirical[nzchar(empirical$gene_token) & !is.na(empirical$gene_token), , drop = FALSE]
  attr(empirical, "empirical_marker_file") <- path
  attr(empirical, "empirical_marker_set_version") <- paste(unique(empirical$marker_source %||% "empirical_roi_marker_sets"), collapse = ";")
  empirical
}

load_wgcna_marker_sets <- function(include_empirical = TRUE, include_legacy_aliases = TRUE, quiet = FALSE) {
  registry <- read_wgcna_marker_registry()
  empirical <- if (isTRUE(include_empirical)) read_empirical_roi_marker_sets() else NULL
  sets <- list()
  metadata <- data.frame(marker_set = character(), marker_source = character(), source_file = character(), stringsAsFactors = FALSE)

  if (!is.null(registry)) {
    ref_sets <- marker_registry_to_sets(registry)
    sets <- c(sets, ref_sets)
    metadata <- rbind(metadata, data.frame(
      marker_set = names(ref_sets),
      marker_source = "reference_registry",
      source_file = attr(registry, "marker_registry_file") %||% NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.null(empirical)) {
    emp_sets <- split(as.character(empirical$GeneSymbol), as.character(empirical$marker_set)) |>
      lapply(function(x) unique(x[nzchar(normalize_gene_token(x))]))
    sets <- c(sets, emp_sets)
    metadata <- rbind(metadata, data.frame(
      marker_set = names(emp_sets),
      marker_source = "empirical_roi_marker_sets",
      source_file = attr(empirical, "empirical_marker_file") %||% NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  if (!length(sets)) {
    if (!isTRUE(quiet)) warning("Falling back to legacy hard-coded WGCNA marker panels; run 03_qc_exploration/04b_import_reference_marker_sources.r to create the registry.", call. = FALSE)
    sets <- wgcna_marker_sets()
    metadata <- data.frame(marker_set = names(sets), marker_source = "legacy_hardcoded_fallback", source_file = NA_character_, stringsAsFactors = FALSE)
  } else if (isTRUE(include_legacy_aliases)) {
    legacy <- wgcna_marker_sets()
    alias_map <- c(
      microglia = "canonical_microglia_homeostatic",
      neuronal_synaptic_neuropil = "canonical_neuronal_synaptic_neuropil",
      neuropil_synaptic_neuronal = "canonical_neuronal_synaptic_neuropil",
      nuclear_soma = "canonical_neuronal_soma_nuclear",
      astrocyte = "canonical_astrocyte",
      oligodendrocyte_myelin = "canonical_oligodendrocyte_myelin",
      endothelial_pericyte_vascular = "canonical_endothelial_vascular",
      mitochondrial_oxphos = "canonical_mitochondrial_oxphos",
      ribosomal_translation = "canonical_ribosomal_translation",
      rnp_rna_processing = "canonical_rnp_rna_processing"
    )
    for (alias in names(alias_map)) {
      target <- unname(alias_map[[alias]])
      if (!alias %in% names(sets)) sets[[alias]] <- sets[[target]] %||% legacy[[alias]]
    }
  }

  sets <- sets[!duplicated(names(sets))]
  attr(sets, "marker_source_metadata") <- metadata
  attr(sets, "marker_registry_version") <- if (!is.null(registry)) attr(registry, "marker_registry_version") else NA_character_
  attr(sets, "empirical_marker_set_version") <- if (!is.null(empirical)) attr(empirical, "empirical_marker_set_version") else NA_character_
  sets
}

standardize_wgcna_metadata <- function(meta, dataset) {
  meta <- as.data.frame(meta, check.names = FALSE, stringsAsFactors = FALSE)
  sample_col <- first_present_col(meta, c("Sample", "sample", "SampleID", "sample_id", "row.names"))
  if (is.na(sample_col)) meta$Sample <- rownames(meta) else meta$Sample <- as.character(meta[[sample_col]])
  animal_col <- first_present_col(meta, c("AnimalID", "Animal", "MouseID", "mouse_id", "animal_id", "subject", "donor"))
  meta$AnimalID <- if (!is.na(animal_col)) as.character(meta[[animal_col]]) else NA_character_
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
    dataset = dataset, level = level, endpoint_id = NA_character_, endpoint_label = NA_character_,
    module_id = NA_character_, supermodule_id = NA_character_,
    module_label = NA_character_, supermodule_label = NA_character_, spatial_unit = NA_character_,
    effect_scope = NA_character_, SpatialUnitType = NA_character_, model_type = NA_character_,
    has_repeated_animals = NA, n_animals = NA_integer_,
    contrast = NA_character_, estimate = NA_real_, SE = NA_real_, statistic = NA_real_,
    p_value = NA_real_, FDR_within_dataset_level = NA_real_, FDR_global = NA_real_,
    evidence_status = "not_supported",
    direction = NA_character_, n_samples = 0L, formula_requested = NA_character_, formula_used = NA_character_,
    dropped_covariates = NA_character_, rank_deficient_model = NA, model_warning = reason,
    stringsAsFactors = FALSE
  )
}

required_group_effect_columns <- c(
  "dataset", "level", "endpoint_id", "endpoint_label", "module_id", "supermodule_id", "module_label", "supermodule_label",
  "spatial_unit", "effect_scope", "SpatialUnitType", "model_type", "has_repeated_animals",
  "n_animals", "contrast", "estimate", "SE", "statistic", "p_value",
  "FDR_within_dataset_level", "FDR_global", "evidence_status", "direction", "n_samples",
  "formula_requested", "formula_used", "dropped_covariates", "rank_deficient_model", "model_warning"
)

required_module_annotation_columns <- c(
  "dataset", "ModuleID", "ModuleColor", "n_proteins", "microenvironment_class", "interpretation_note"
)

required_interpretable_columns <- c(
  "dataset", "level", "contrast", "estimate", "p_value", "FDR_global", "interpretation_sentence"
)
