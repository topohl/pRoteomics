#!/usr/bin/env Rscript
#
# Import cached/package-accessible reference marker evidence into a normalized WGCNA marker registry.
#
# This script builds config/marker_panels/wgcna_reference_marker_sets.csv from:
#   1. conservative canonical fallback panels,
#   2. optional cached external cell-type marker sources listed in reference_marker_sources.yml,
#   3. optional automatically cached GO/MGI and SynGO evidence for compartment-fidelity markers,
#   4. a generated source-controlled compartment fidelity marker table:
#      config/marker_panels/compartment_fidelity_marker_sets.csv
#
# The compartment-fidelity layer is designed for 04c/04d QC:
#   Soma markers, Neuropil markers, Microglia/PVM markers.
# It is QC/interpreter support only; it is not purity estimation and must not be
# derived from CON/RES/SUS differential abundance results.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "wgcna_downstream_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
dry_run <- is_dry_run()
SUBSTEP_ID <- "reference_marker_import"
PATHS <- create_module_dirs("03_qc_exploration", file.path("reference_marker_import"))
manifest_file <- Sys.getenv(
  "PROTEOMICS_REFERENCE_MARKER_SOURCES_YML",
  unset = path_external("reference_markers", "reference_marker_sources.yml")
)
allow_download <- tolower(Sys.getenv("PROTEOMICS_REFERENCE_MARKERS_ALLOW_DOWNLOAD", unset = "false")) %in% c("1", "true", "yes", "y")
top_n <- suppressWarnings(as.integer(Sys.getenv("PROTEOMICS_REFERENCE_MARKER_TOP_N_PER_CLASS", unset = "100")))
if (!is.finite(top_n) || top_n < 1L) top_n <- 100L
min_specificity <- suppressWarnings(as.numeric(Sys.getenv("PROTEOMICS_REFERENCE_MARKER_MIN_SPECIFICITY", unset = NA_character_)))

# External-preferred default: use Allen, GO/MGI, and SynGO for fidelity scoring
# when available. Manual curated seeds are fallback-only unless the policy is
# overridden.
fidelity_source_policy <- tolower(Sys.getenv("PROTEOMICS_FIDELITY_SOURCE_POLICY", unset = "external_preferred"))
if (!fidelity_source_policy %in% c("external_preferred", "manual_plus_external", "external_only", "manual_only")) {
  fidelity_source_policy <- "external_preferred"
}
auto_include_syngo <- tolower(Sys.getenv("PROTEOMICS_FIDELITY_AUTO_INCLUDE_SYNGo", unset = "true")) %in% c("1", "true", "yes", "y")
max_syngo_per_subpanel <- suppressWarnings(as.integer(Sys.getenv("PROTEOMICS_FIDELITY_MAX_SYNGo_PER_SUBPANEL", unset = "50")))
if (!is.finite(max_syngo_per_subpanel) || max_syngo_per_subpanel < 1L) max_syngo_per_subpanel <- 50L
auto_include_go_mgi <- tolower(Sys.getenv("PROTEOMICS_FIDELITY_AUTO_INCLUDE_GO_MGI", unset = "true")) %in% c("1", "true", "yes", "y")
max_go_mgi_per_subpanel <- suppressWarnings(as.integer(Sys.getenv("PROTEOMICS_FIDELITY_MAX_GO_MGI_PER_SUBPANEL", unset = "30")))
if (!is.finite(max_go_mgi_per_subpanel) || max_go_mgi_per_subpanel < 1L) max_go_mgi_per_subpanel <- 30L
go_mgi_include_evidence <- Sys.getenv("PROTEOMICS_FIDELITY_GO_MGI_EVIDENCE_CODES", unset = "EXP,IDA,IMP,IPI,IGI")
go_mgi_include_subpanels <- Sys.getenv(
  "PROTEOMICS_FIDELITY_GO_MGI_SUBPANELS",
  unset = "active_zone,postsynaptic_density,synaptic_vesicle,dendritic_spine,nuclear_speck,nuclear_matrix,spliceosome,rnp,chromatin"
)

# URLs are stable for GO/MGI current releases. For SynGO, prefer a manifest/env URL
# because public download URLs may change.
mgi_gaf_url <- Sys.getenv("PROTEOMICS_MGI_GAF_URL", unset = "https://current.geneontology.org/annotations/mgi.gaf.gz")
go_basic_obo_url <- Sys.getenv("PROTEOMICS_GO_BASIC_OBO_URL", unset = "https://current.geneontology.org/ontology/go-basic.obo")
mgi_homology_url <- Sys.getenv("PROTEOMICS_MGI_HOMOLOGY_URL", unset = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt")
syngo_url_env <- Sys.getenv("PROTEOMICS_SYNGo_URL", unset = "")
syngo_file_name_env <- Sys.getenv("PROTEOMICS_SYNGo_FILE_NAME", unset = "annotations.xlsx")
syngo_genes_file_name_env <- Sys.getenv("PROTEOMICS_SYNGo_GENES_FILE_NAME", unset = "genes.xlsx")
syngo_ontologies_file_name_env <- Sys.getenv("PROTEOMICS_SYNGo_ONTOLOGIES_FILE_NAME", unset = "ontologies.xlsx")

if (dry_run) {
  invisible(lapply(unlist(PATHS), dir_create))
  dry_run_line("Script", "03_qc_exploration/04b_import_reference_marker_sources.r")
  dry_run_line("Manifest", manifest_file, if (file.exists(manifest_file)) "PASS" else "FAIL")
  dry_run_line("Live download allowed", allow_download, "INFO")
  dry_run_line("GO/MGI GAF cache", path_external("reference_markers", "go_mgi", "raw", "mgi.gaf.gz"), "INFO")
  dry_run_line("GO OBO cache", path_external("reference_markers", "go_mgi", "raw", "go-basic.obo"), "INFO")
  dry_run_line("MGI human-mouse homology cache", path_external("reference_markers", "go_mgi", "raw", "HOM_MouseHumanSequence.rpt"), "INFO")
  dry_run_line("SynGO cache folder", path_external("reference_markers", "syngo", "raw"), "INFO")
  dry_run_line("Compartment fidelity marker table", repo_path("config", "marker_panels", "compartment_fidelity_marker_sets.csv"), "INFO")
  dry_run_line("Output registry", repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"), "INFO")
  quit(status = if (file.exists(manifest_file)) 0L else 1L, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "svglite", "readr", "yaml", "readxl")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

if (!file.exists(manifest_file)) stop("Missing marker source manifest: ", manifest_file, call. = FALSE)
manifest <- yaml::read_yaml(manifest_file)
sources <- manifest$sources %||% list()

# -----------------------------------------------------------------------------
# Conservative canonical fallback panels
# -----------------------------------------------------------------------------

canonical_marker_panels <- list(
  canonical_microglia_homeostatic = c("P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Hexb", "Fcrls", "Sall1", "Siglech", "Gpr34", "Mertk", "Aif1"),
  canonical_microglia_phagolysosomal_state = c("Tyrobp", "Trem2", "Apoe", "Lpl", "Cst7", "Ctsb", "Ctsd", "Lgals3", "Itgax", "Axl", "C1qa", "C1qb", "C1qc"),
  canonical_neuronal_synaptic_neuropil = c("Snap25", "Syp", "Syn1", "Vamp2", "Stx1a", "Stxbp1", "Dlg4", "Camk2a", "Camk2b", "Map2", "Nefl", "Nefm", "Rbfox3", "Grin1", "Gria1"),
  canonical_neuronal_soma_nuclear = c("Rbfox3", "Map2", "Tubb3", "H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3", "Matr3", "Srsf3", "Ddx39b"),
  canonical_astrocyte = c("Aqp4", "Gfap", "Aldh1l1", "Slc1a2", "Slc1a3", "Glul", "Aldoc", "Gja1", "S100b"),
  canonical_oligodendrocyte_myelin = c("Mbp", "Plp1", "Mog", "Cnp", "Mag", "Mobp", "Cldn11", "Myrf", "Olig1", "Olig2"),
  canonical_opc = c("Pdgfra", "Cspg4", "Vcan", "Sox10", "Olig1", "Olig2"),
  canonical_endothelial_vascular = c("Cldn5", "Pecam1", "Kdr", "Flt1", "Slco1a4"),
  canonical_pericyte_vascular = c("Rgs5", "Pdgfrb", "Vtn", "Acta2"),
  canonical_peripheral_myeloid_caution = c("Lyz2", "Cd74", "H2-Ab1", "Fcgr1", "Ccr2", "Ly6c2", "Itgam"),
  canonical_mitochondrial_oxphos = c("Ndufs1", "Ndufa9", "Sdha", "Uqcrc2", "Cox4i1", "Atp5f1a", "Atp5f1b"),
  canonical_ribosomal_translation = c("Rpl3", "Rpl4", "Rpl5", "Rps3", "Rps6", "Eef1a1", "Eef2"),
  canonical_rnp_rna_processing = c("Hnrnpa2b1", "Hnrnpc", "Sfpq", "Snrnp70", "Ddx5", "Ddx17", "Pabpc1")
)

panel_meta <- data.frame(
  marker_set = names(canonical_marker_panels),
  cell_class = c("microglia", "microglia", "neuron", "neuron", "astrocyte", "oligodendrocyte", "opc", "endothelial", "pericyte", "peripheral_myeloid", "mitochondrial", "ribosomal", "rnp"),
  cell_state = c("homeostatic_identity", "phagolysosomal_complement_state", "synaptic_neuropil", "soma_nuclear", "canonical", "myelin", "canonical", "vascular", "vascular", "caution", "oxphos", "translation", "rna_processing"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

is_regular_file <- function(path) {
  isTRUE(file.exists(path)) && !isTRUE(file.info(path)$isdir)
}

report_reference_marker_status <- function(source_name, status, detail = NULL) {
  msg <- paste0("[reference markers] ", source_name, ": ", status)
  if (!is.null(detail) && nzchar(as.character(detail))) {
    msg <- paste0(msg, " (", detail, ")")
  }
  message(msg)
}

read_any_table <- function(path) {
  if (!is_regular_file(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(tryCatch(readRDS(path), error = function(e) NULL))
  if (ext == "xlsx" || ext == "xls") {
    if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
    return(tryCatch(as.data.frame(readxl::read_excel(path), check.names = FALSE), error = function(e) NULL))
  }
  if (ext == "tsv" || ext == "txt") return(tryCatch(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL))
  if (ext == "csv") return(tryCatch(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL))
  NULL
}

collapse_unique <- function(x) paste(sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))]))), collapse = ";")

parse_env_list <- function(x) {
  x <- unlist(strsplit(as.character(x %||% ""), ",", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  x[nzchar(x)]
}

file_hash <- function(path) {
  if (!is_regular_file(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  unname(tools::md5sum(path))
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

resolve_local_path <- function(path) {
  path <- as.character(path %||% "")
  if (!nzchar(path)) return(getwd())
  # Avoid repo_path() here: some path helpers in this repo are strict about
  # file-vs-directory semantics. Marker-source local_dir entries are directories.
  if (grepl("^[A-Za-z]:[/\\]", path) || startsWith(path, "/") || startsWith(path, "~")) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(getwd(), path), winslash = "/", mustWork = FALSE)
}

download_if_missing <- function(url, dest, allow_download = FALSE, source_name = "unknown") {
  safe_dir_create(dirname(dest))
  if (is_regular_file(dest)) {
    report_reference_marker_status(source_name, "cached", dest)
    return(data.frame(
      source_name = source_name,
      url = url,
      file = dest,
      status = "cached",
      downloaded_at = NA_character_,
      hash = file_hash(dest),
      stringsAsFactors = FALSE
    ))
  }
  if (dir.exists(dest)) {
    report_reference_marker_status(source_name, "cache path is a directory; waiting for a file", dest)
  }
  if (!allow_download || is.na(url) || !nzchar(url)) {
    report_reference_marker_status(
      source_name,
      if (allow_download) "missing optional file; no download URL configured" else "missing optional file; download disabled",
      dest
    )
    return(data.frame(
      source_name = source_name,
      url = url,
      file = dest,
      status = "missing_download_disabled",
      downloaded_at = NA_character_,
      hash = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  tmp <- tempfile(fileext = paste0(".", tools::file_ext(dest)))
  ok <- tryCatch({
    utils::download.file(url, tmp, mode = "wb", quiet = FALSE)
    if (!file.exists(tmp) || file.info(tmp)$size <= 0) stop("Downloaded file is empty")
    file.rename(tmp, dest)
  }, error = function(e) {
    warning("Download failed for ", source_name, ": ", conditionMessage(e), call. = FALSE)
    FALSE
  })
  final_status <- if (isTRUE(ok) && is_regular_file(dest)) "downloaded" else "download_failed"
  report_reference_marker_status(source_name, final_status, dest)
  data.frame(
    source_name = source_name,
    url = url,
    file = dest,
    status = final_status,
    downloaded_at = if (identical(final_status, "downloaded")) as.character(Sys.time()) else NA_character_,
    hash = if (is_regular_file(dest)) file_hash(dest) else NA_character_,
    stringsAsFactors = FALSE
  )
}

mouse_symbol_guess <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x <- gsub("_MOUSE$", "", x, ignore.case = TRUE)
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_character_)
    # Already mouse-like symbols such as Snap25 or H2-Ab1 are left unchanged.
    if (grepl("[a-z]", s)) return(s)
    paste0(substr(s, 1, 1), tolower(substr(s, 2, nchar(s))))
  }, character(1))
}

make_registry_row <- function(marker_set, fidelity_marker_class, fidelity_subpanel, cell_class, cell_state,
                              genes, source_type, source_name, source_reference, source_term_or_category,
                              selection_rule, confidence, use_for, include_in_fidelity_score = TRUE,
                              notes = "", evidence_level = NA_character_) {
  data.frame(
    marker_set = marker_set,
    fidelity_marker_class = fidelity_marker_class,
    fidelity_subpanel = fidelity_subpanel,
    cell_class = cell_class,
    cell_state = cell_state,
    gene_symbol = unique(as.character(genes)),
    source_type = source_type,
    source_name = source_name,
    source_reference = source_reference,
    source_term_or_category = source_term_or_category,
    evidence_level = evidence_level,
    selection_rule = selection_rule,
    confidence = confidence,
    use_for = use_for,
    include_in_fidelity_score = as.character(include_in_fidelity_score),
    notes = notes,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol))
}

# Empty-template helper. Several external sources are optional. Returning a bare
# data.frame() is unsafe because downstream dplyr code expects the fidelity
# metadata columns to exist even when no rows were imported.
empty_fidelity_candidate_table <- function() {
  data.frame(
    marker_set = character(),
    fidelity_marker_class = character(),
    fidelity_subpanel = character(),
    cell_class = character(),
    cell_state = character(),
    gene_symbol = character(),
    source_type = character(),
    source_name = character(),
    source_reference = character(),
    source_term_or_category = character(),
    evidence_level = character(),
    selection_rule = character(),
    confidence = character(),
    use_for = character(),
    include_in_fidelity_score = character(),
    notes = character(),
    stringsAsFactors = FALSE
  )
}

ensure_fidelity_candidate_cols <- function(df) {
  template <- empty_fidelity_candidate_table()
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(template)
  missing <- setdiff(names(template), names(df))
  for (col in missing) df[[col]] <- template[[col]][NA_integer_]
  df[, names(template), drop = FALSE]
}

# -----------------------------------------------------------------------------
# Existing generic external marker import
# -----------------------------------------------------------------------------

coerce_marker_table <- function(df, source) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(NULL)
  gene_col <- first_present_col(df, c("gene_symbol", "GeneSymbol", "gene", "Gene", "symbol", "marker", "Marker"))
  class_col <- first_present_col(df, c("cell_class", "cell_type", "celltype", "CellType", "panel_name", "class", "cluster", "reference_cell_label"))
  score_col <- first_present_col(df, c("marker_rank_score", "specificity_score", "specificity_ratio", "specificity", "auc", "score", "logFC", "avg_log2FC", "rank"))
  mean_target_col <- first_present_col(df, c("mean_target", "mean_expr", "mean_expression", "avg_expr", "avg_expression"))
  mean_other_col <- first_present_col(df, c("mean_other", "mean_other_mean_expr", "max_other_mean_expr"))
  pct_target_col <- first_present_col(df, c("pct_target", "frac_clusters_expr_gt_0", "pct_expr", "pct.1"))
  pct_other_col <- first_present_col(df, c("pct_other", "pct.2"))
  evidence_col <- first_present_col(df, c("evidence_type", "selection_rule"))
  use_for_col <- first_present_col(df, c("use_for"))
  reference_col <- first_present_col(df, c("reference_source", "source_reference"))
  if (is.na(gene_col) || is.na(class_col)) return(NULL)
  out <- data.frame(
    source_name = source$source_name %||% "unknown",
    source_version = source$source_version %||% "cached_local",
    source_type = source$source_type %||% "local_file",
    gene_symbol = as.character(df[[gene_col]]),
    cell_class = as.character(df[[class_col]]),
    cell_state = NA_character_,
    reference_cell_label = as.character(df[[class_col]]),
    mean_target = if (!is.na(mean_target_col)) suppressWarnings(as.numeric(df[[mean_target_col]])) else NA_real_,
    mean_other = if (!is.na(mean_other_col)) suppressWarnings(as.numeric(df[[mean_other_col]])) else NA_real_,
    logFC_target_vs_other = if (!is.na(score_col) && grepl("log", score_col, ignore.case = TRUE)) suppressWarnings(as.numeric(df[[score_col]])) else NA_real_,
    specificity_score = if (!is.na(score_col)) suppressWarnings(as.numeric(df[[score_col]])) else NA_real_,
    pct_target = if (!is.na(pct_target_col)) suppressWarnings(as.numeric(df[[pct_target_col]])) else NA_real_,
    pct_other = if (!is.na(pct_other_col)) suppressWarnings(as.numeric(df[[pct_other_col]])) else NA_real_,
    rank_within_source = NA_integer_,
    selection_rule = if (!is.na(evidence_col)) as.character(df[[evidence_col]]) else "cached_marker_table",
    selected = TRUE,
    confidence = "external_cached",
    notes = source$notes %||% "",
    source_reference = if (!is.na(reference_col)) as.character(df[[reference_col]]) else source$source_reference %||% source$source_name %||% "unknown",
    use_for = if (!is.na(use_for_col)) as.character(df[[use_for_col]]) else "annotation_reporting_sensitivity_interpretation",
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(normalize_gene_token(out$gene_symbol)), , drop = FALSE]
  out <- out |>
    dplyr::group_by(.data$source_name, .data$cell_class) |>
    dplyr::arrange(dplyr::desc(dplyr::coalesce(.data$specificity_score, abs(.data$logFC_target_vs_other), 0)), .by_group = TRUE) |>
    dplyr::mutate(rank_within_source = dplyr::row_number()) |>
    dplyr::ungroup()
  if (is.finite(min_specificity)) out <- out |> dplyr::filter(is.na(.data$specificity_score) | .data$specificity_score >= min_specificity)
  out |> dplyr::filter(.data$rank_within_source <= top_n)
}

# -----------------------------------------------------------------------------
# GO/MGI parsing and candidate extraction
# -----------------------------------------------------------------------------

read_go_obo_terms <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  term_starts <- which(lines == "[Term]")
  if (!length(term_starts)) return(NULL)
  term_ends <- c(term_starts[-1] - 1L, length(lines))

  parse_one <- function(block) {
    get_first <- function(prefix) {
      x <- block[startsWith(block, prefix)]
      if (length(x)) sub(prefix, "", x[[1]], fixed = TRUE) else NA_character_
    }
    id <- get_first("id: ")
    name <- get_first("name: ")
    namespace <- get_first("namespace: ")
    obsolete <- any(block == "is_obsolete: true")
    is_a <- sub(" !.*$", "", sub("is_a: ", "", block[startsWith(block, "is_a: ")], fixed = TRUE))
    part_of <- block[grepl("^relationship: part_of ", block)]
    part_of <- sub(" !.*$", "", sub("relationship: part_of ", "", part_of))
    data.frame(
      go_id = id,
      name = name,
      namespace = namespace,
      is_obsolete = obsolete,
      parents = paste(c(is_a, part_of), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(lapply(seq_along(term_starts), function(i) {
    parse_one(lines[term_starts[[i]]:term_ends[[i]]])
  })) |>
    dplyr::filter(!is.na(.data$go_id), nzchar(.data$go_id))
}

go_descendants <- function(terms, root_names) {
  if (is.null(terms) || !nrow(terms)) return(character())
  roots <- terms$go_id[tolower(terms$name) %in% tolower(root_names)]
  if (!length(roots)) return(character())

  parent_map <- strsplit(terms$parents, ";", fixed = TRUE)
  names(parent_map) <- terms$go_id
  parent_map <- lapply(parent_map, function(x) x[nzchar(x) & !is.na(x)])

  children <- lapply(roots, function(x) character())
  names(children) <- roots
  all_ids <- terms$go_id
  current <- roots
  seen <- roots
  repeat {
    kids <- all_ids[vapply(parent_map, function(parents) any(parents %in% current), logical(1))]
    kids <- setdiff(kids, seen)
    if (!length(kids)) break
    seen <- unique(c(seen, kids))
    current <- kids
  }
  seen
}

read_mgi_gaf <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(
    readr::read_tsv(path, comment = "!", col_names = FALSE, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) < 9L) return(NULL)
  names(df)[seq_len(min(ncol(df), 17L))] <- c(
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB_Reference",
    "Evidence_Code", "With_From", "Aspect", "DB_Object_Name", "DB_Object_Synonym",
    "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension",
    "Gene_Product_Form_ID"
  )[seq_len(min(ncol(df), 17L))]
  df |>
    dplyr::filter(.data$Aspect == "C") |>
    dplyr::transmute(
      gene_symbol = as.character(.data$DB_Object_Symbol),
      go_id = as.character(.data$GO_ID),
      evidence_code = as.character(.data$Evidence_Code),
      source_reference = as.character(.data$DB_Reference),
      assigned_by = as.character(.data$Assigned_By),
      annotation_date = as.character(.data$Date)
    ) |>
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol), !is.na(.data$go_id), nzchar(.data$go_id))
}

extract_go_mgi_candidates <- function(gaf_file, obo_file) {
  gaf <- read_mgi_gaf(gaf_file)
  terms <- read_go_obo_terms(obo_file)
  if (is.null(gaf) || is.null(terms) || !nrow(gaf) || !nrow(terms)) return(data.frame())

  term_map <- terms |> dplyr::select(go_id, go_name = name)

  neuropil_terms <- list(
    fidelity_neuropil_synapse = c("synapse"),
    fidelity_neuropil_presynapse = c("presynapse"),
    fidelity_neuropil_postsynapse = c("postsynapse"),
    fidelity_neuropil_synaptic_vesicle = c("synaptic vesicle"),
    fidelity_neuropil_active_zone = c("presynaptic active zone"),
    fidelity_neuropil_postsynaptic_density = c("postsynaptic density"),
    fidelity_neuropil_dendritic_spine = c("dendritic spine"),
    fidelity_neuropil_neuron_projection = c("axon", "dendrite", "neuron projection")
  )

  soma_terms <- list(
    fidelity_soma_nucleus = c("nucleus"),
    fidelity_soma_nucleoplasm = c("nucleoplasm"),
    fidelity_soma_chromatin = c("chromatin"),
    fidelity_soma_nuclear_matrix = c("nuclear matrix"),
    fidelity_soma_nuclear_speck = c("nuclear speck"),
    fidelity_soma_spliceosome = c("spliceosomal complex"),
    fidelity_soma_rnp = c("ribonucleoprotein complex")
  )

  make_candidate <- function(marker_set, root_names, marker_class, subpanel, cell_state) {
    ids <- go_descendants(terms, root_names)
    if (!length(ids)) return(data.frame())
    gaf |>
      dplyr::filter(.data$go_id %in% ids) |>
      dplyr::left_join(term_map, by = "go_id") |>
      dplyr::distinct(.data$gene_symbol, .data$go_id, .data$evidence_code, .keep_all = TRUE) |>
      dplyr::transmute(
        marker_set = marker_set,
        fidelity_marker_class = marker_class,
        fidelity_subpanel = subpanel,
        cell_class = "neuron",
        cell_state = cell_state,
        gene_symbol = mouse_symbol_guess(.data$gene_symbol),
        source_type = "go_mgi_cached_download",
        source_name = "GO_MGI",
        source_reference = paste0("mgi.gaf.gz; go-basic.obo; root_terms=", paste(root_names, collapse = ";")),
        source_term_or_category = paste0(.data$go_id, ":", .data$go_name),
        evidence_level = .data$evidence_code,
        selection_rule = "GO/MGI cellular component descendant candidate; not auto-included in final score unless curated seed or explicitly enabled",
        confidence = "external_go_candidate",
        use_for = "compartment_fidelity_candidate_audit",
        include_in_fidelity_score = "FALSE",
        notes = "GO/MGI candidate only. Broad terms require manual review before final fidelity scoring.",
        stringsAsFactors = FALSE
      )
  }

  dplyr::bind_rows(
    lapply(names(neuropil_terms), function(ms) {
      make_candidate(ms, neuropil_terms[[ms]], "Neuropil markers", sub("^fidelity_neuropil_", "", ms), "synaptic_neuropil_candidate")
    }),
    lapply(names(soma_terms), function(ms) {
      make_candidate(ms, soma_terms[[ms]], "Soma markers", sub("^fidelity_soma_", "", ms), "soma_nuclear_candidate")
    })
  ) |>
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol)) |>
    dplyr::distinct(.data$marker_set, .data$gene_symbol, .data$source_term_or_category, .keep_all = TRUE)
}

# -----------------------------------------------------------------------------
# Human-to-mouse homology support for HGNC-based sources such as SynGO
# -----------------------------------------------------------------------------

read_mgi_homology <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)

  names(df) <- gsub("\\s+", "_", names(df))
  class_col <- first_present_col(df, c("DB_Class_Key", "HomoloGene_ID", "Homology_Class", "class_key"))
  org_col <- first_present_col(df, c("Common_Organism_Name", "Organism", "organism"))
  symbol_col <- first_present_col(df, c("Symbol", "symbol"))
  hgnc_col <- first_present_col(df, c("HGNC_ID", "hgnc_id"))

  if (is.na(class_col) || is.na(org_col) || is.na(symbol_col)) return(NULL)

  tab <- df |>
    dplyr::transmute(
      homology_class = as.character(.data[[class_col]]),
      organism = as.character(.data[[org_col]]),
      symbol = as.character(.data[[symbol_col]]),
      hgnc_id = if (!is.na(hgnc_col)) as.character(.data[[hgnc_col]]) else NA_character_
    ) |>
    dplyr::filter(!is.na(.data$homology_class), nzchar(.data$homology_class), !is.na(.data$symbol), nzchar(.data$symbol))

  mouse <- tab |>
    dplyr::filter(grepl("mouse|mus musculus", .data$organism, ignore.case = TRUE)) |>
    dplyr::group_by(.data$homology_class) |>
    dplyr::summarise(mouse_symbol = dplyr::first(.data$symbol), .groups = "drop")

  human <- tab |>
    dplyr::filter(grepl("human|homo sapiens", .data$organism, ignore.case = TRUE)) |>
    dplyr::group_by(.data$homology_class) |>
    dplyr::summarise(
      hgnc_symbol = dplyr::first(.data$symbol),
      hgnc_id = dplyr::first(.data$hgnc_id[!is.na(.data$hgnc_id) & nzchar(.data$hgnc_id)]),
      .groups = "drop"
    )

  human |>
    dplyr::left_join(mouse, by = "homology_class") |>
    dplyr::filter(!is.na(.data$mouse_symbol), nzchar(.data$mouse_symbol)) |>
    dplyr::mutate(
      hgnc_symbol_key = normalize_gene_token(.data$hgnc_symbol),
      hgnc_id_key = toupper(gsub("[^A-Za-z0-9]", "", as.character(.data$hgnc_id)))
    ) |>
    dplyr::distinct(.data$hgnc_symbol_key, .data$hgnc_id_key, .keep_all = TRUE)
}

map_hgnc_to_mouse <- function(hgnc_symbol, hgnc_id = NA_character_, homology_map = NULL) {
  hgnc_symbol <- as.character(hgnc_symbol)
  hgnc_id <- as.character(hgnc_id)
  if (is.null(homology_map) || !nrow(homology_map)) return(mouse_symbol_guess(hgnc_symbol))

  symbol_key <- normalize_gene_token(hgnc_symbol)
  id_key <- toupper(gsub("[^A-Za-z0-9]", "", hgnc_id))
  out <- mouse_symbol_guess(hgnc_symbol)

  idx_symbol <- match(symbol_key, homology_map$hgnc_symbol_key)
  idx_id <- match(id_key, homology_map$hgnc_id_key)
  mapped <- ifelse(!is.na(idx_symbol), homology_map$mouse_symbol[idx_symbol], NA_character_)
  mapped <- ifelse(is.na(mapped) & !is.na(idx_id), homology_map$mouse_symbol[idx_id], mapped)
  out[!is.na(mapped) & nzchar(mapped)] <- mapped[!is.na(mapped) & nzchar(mapped)]
  out
}

# -----------------------------------------------------------------------------
# SynGO parsing and candidate extraction
# -----------------------------------------------------------------------------

find_source <- function(source_names) {
  for (src in sources) {
    if ((src$source_name %||% "") %in% source_names) return(src)
  }
  NULL
}

classify_syngo_subpanel <- function(go_name, go_domain = NA_character_) {
  x_l <- tolower(as.character(go_name))
  domain_l <- tolower(as.character(go_domain))
  dplyr::case_when(
    !is.na(domain_l) & domain_l == "cc" & grepl("synaptic vesicle|vesicle", x_l) ~ "presynaptic_vesicle",
    !is.na(domain_l) & domain_l == "cc" & grepl("active zone", x_l) ~ "active_zone",
    !is.na(domain_l) & domain_l == "cc" & grepl("postsynaptic density|post-synaptic density|psd", x_l) ~ "postsynaptic_density",
    !is.na(domain_l) & domain_l == "cc" & grepl("postsynap", x_l) ~ "postsynapse",
    !is.na(domain_l) & domain_l == "cc" & grepl("presynap", x_l) ~ "presynapse",
    !is.na(domain_l) & domain_l == "cc" & grepl("synapse|synaptic membrane|synaptic specialization", x_l) ~ "synapse_other",
    !is.na(domain_l) & domain_l == "bp" & grepl("synaptic vesicle", x_l) ~ "bp_synaptic_vesicle_process",
    !is.na(domain_l) & domain_l == "bp" & grepl("neurotransmitter|synaptic transmission|synaptic signaling|synapse organization|presynaptic|postsynaptic", x_l) ~ "bp_synaptic_process",
    is.na(domain_l) & grepl("synaptic vesicle|vesicle", x_l) ~ "presynaptic_vesicle",
    is.na(domain_l) & grepl("active zone", x_l) ~ "active_zone",
    is.na(domain_l) & grepl("postsynaptic density|post-synaptic density|psd", x_l) ~ "postsynaptic_density",
    is.na(domain_l) & grepl("postsynap", x_l) ~ "postsynapse",
    is.na(domain_l) & grepl("presynap", x_l) ~ "presynapse",
    is.na(domain_l) & grepl("synaptic", x_l) ~ "synaptic_process",
    TRUE ~ NA_character_
  )
}

syngo_source_strength <- function(df) {
  # SynGO annotations.xlsx provides boolean/evidence columns. Prefer rows with
  # traceable experimental evidence and CC localization for strict neuropil QC.
  cols <- names(df)
  ev_cols <- cols[grepl("^evidence_", cols, ignore.case = TRUE)]
  if (!length(ev_cols)) return(rep("unknown", nrow(df)))

  ev <- df[, ev_cols, drop = FALSE]
  ev_text <- apply(ev, 1, function(x) paste(as.character(x[!is.na(x)]), collapse = ";"))
  dplyr::case_when(
    grepl("mic:|ephys:|ophys:|biochem:|target:antibody|target:tag", ev_text, ignore.case = TRUE) ~ "traceable_experimental",
    grepl("TRUE", ev_text, ignore.case = TRUE) ~ "nontraceable_author_statement",
    nzchar(ev_text) ~ "evidence_reported",
    TRUE ~ "unknown"
  )
}

find_existing_file <- function(dir, candidates) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (x in candidates) {
    p <- file.path(dir, x)
    if (is_regular_file(p)) return(p)
  }
  if (!dir.exists(dir)) return(file.path(dir, candidates[[1]]))
  hits <- unlist(lapply(candidates, function(x) list.files(dir, pattern = paste0("^", gsub("\\.", "\\\\.", x), "$"), full.names = TRUE, ignore.case = TRUE)), use.names = FALSE)
  hits <- hits[vapply(hits, is_regular_file, logical(1))]
  if (length(hits)) hits[[1]] else file.path(dir, candidates[[1]])
}

extract_syngo_candidates <- function(annotation_file, genes_file = NA_character_, ontologies_file = NA_character_,
                                     homology_map = NULL, source_reference = "SynGO cached annotation table") {
  if (!is_regular_file(annotation_file)) return(empty_fidelity_candidate_table())

  ann <- read_any_table(annotation_file)
  if (is.null(ann) || !nrow(ann)) return(empty_fidelity_candidate_table())

  # SynGO 1.3 complete data export expected columns include:
  # id, uniprot_id, hgnc_id, hgnc_symbol, pubmed_id, go_id, go_name, go_domain,
  # evidence_biological_system, evidence_protein_targeting,
  # evidence_experiment_assay, evidence_nontracable_author_statement.
  gene_col <- first_present_col(ann, c("hgnc_symbol", "HGNC symbol", "HGNC", "gene_symbol", "GeneSymbol", "gene", "Gene", "symbol"))
  hgnc_id_col <- first_present_col(ann, c("hgnc_id", "HGNC ID", "HGNC_ID"))
  go_name_col <- first_present_col(ann, c("go_name", "GO name", "GO_name", "go_term", "GO_term", "term", "Term", "description", "Description"))
  go_id_col <- first_present_col(ann, c("go_id", "GO_ID", "GOid", "goid", "GO"))
  go_domain_col <- first_present_col(ann, c("go_domain", "GO domain", "domain", "Domain"))
  pubmed_col <- first_present_col(ann, c("pubmed_id", "PubMed", "PMID", "pmid"))
  uniprot_col <- first_present_col(ann, c("uniprot_id", "UniProt", "uniprot"))

  if (is.na(gene_col) || is.na(go_name_col)) return(empty_fidelity_candidate_table())

  ann$source_strength <- syngo_source_strength(ann)
  ann$mouse_symbol <- map_hgnc_to_mouse(
    ann[[gene_col]],
    if (!is.na(hgnc_id_col)) ann[[hgnc_id_col]] else NA_character_,
    homology_map = homology_map
  )

  out <- ann |>
    dplyr::transmute(
      hgnc_symbol = as.character(.data[[gene_col]]),
      hgnc_id = if (!is.na(hgnc_id_col)) as.character(.data[[hgnc_id_col]]) else NA_character_,
      gene_symbol = as.character(.data$mouse_symbol),
      go_id = if (!is.na(go_id_col)) as.character(.data[[go_id_col]]) else NA_character_,
      go_name = as.character(.data[[go_name_col]]),
      go_domain = if (!is.na(go_domain_col)) as.character(.data[[go_domain_col]]) else NA_character_,
      pubmed_id = if (!is.na(pubmed_col)) as.character(.data[[pubmed_col]]) else NA_character_,
      uniprot_id = if (!is.na(uniprot_col)) as.character(.data[[uniprot_col]]) else NA_character_,
      evidence_level = as.character(.data$source_strength)
    ) |>
    dplyr::mutate(
      fidelity_subpanel = classify_syngo_subpanel(.data$go_name, .data$go_domain),
      strict_cc = tolower(.data$go_domain) == "cc" & .data$fidelity_subpanel %in% c("presynaptic_vesicle", "active_zone", "postsynaptic_density", "postsynapse", "presynapse", "synapse_other"),
      supportive_bp = tolower(.data$go_domain) == "bp" & .data$fidelity_subpanel %in% c("bp_synaptic_vesicle_process", "bp_synaptic_process")
    ) |>
    dplyr::filter(!is.na(.data$gene_symbol), nzchar(.data$gene_symbol), !is.na(.data$fidelity_subpanel), nzchar(.data$fidelity_subpanel)) |>
    dplyr::group_by(.data$gene_symbol, .data$fidelity_subpanel) |>
    dplyr::arrange(dplyr::desc(.data$strict_cc), dplyr::desc(.data$evidence_level == "traceable_experimental"), .by_group = TRUE) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$fidelity_subpanel) |>
    dplyr::arrange(dplyr::desc(.data$strict_cc), dplyr::desc(.data$evidence_level == "traceable_experimental"), .data$gene_symbol, .by_group = TRUE) |>
    dplyr::mutate(rank_within_subpanel = dplyr::row_number()) |>
    dplyr::ungroup()

  if (!nrow(out)) return(empty_fidelity_candidate_table())

  out |>
    dplyr::transmute(
      marker_set = paste0("fidelity_neuropil_syngo_", .data$fidelity_subpanel),
      fidelity_marker_class = "Neuropil markers",
      fidelity_subpanel = paste0("syngo_", .data$fidelity_subpanel),
      cell_class = "neuron",
      cell_state = dplyr::case_when(
        .data$strict_cc ~ "synaptic_neuropil_localization",
        .data$supportive_bp ~ "synaptic_neuropil_process_supportive",
        TRUE ~ "synaptic_neuropil_candidate"
      ),
      gene_symbol = .data$gene_symbol,
      source_type = "syngo_v1_3_cached_export",
      source_name = "SynGO",
      source_reference = paste0(source_reference, "; annotations=", basename(annotation_file),
                                ifelse(file.exists(genes_file), paste0("; genes=", basename(genes_file)), ""),
                                ifelse(file.exists(ontologies_file), paste0("; ontologies=", basename(ontologies_file)), "")),
      source_term_or_category = paste0(.data$go_domain, ":", .data$go_id, ":", .data$go_name,
                                       ifelse(!is.na(.data$pubmed_id) & nzchar(.data$pubmed_id), paste0(";PMID=", .data$pubmed_id), ""),
                                       ifelse(!is.na(.data$uniprot_id) & nzchar(.data$uniprot_id), paste0(";UniProt=", .data$uniprot_id), ""),
                                       ifelse(!is.na(.data$hgnc_symbol) & nzchar(.data$hgnc_symbol), paste0(";HGNC=", .data$hgnc_symbol), "")),
      evidence_level = .data$evidence_level,
      selection_rule = paste0("SynGO 1.3 annotation export; strict_cc=", .data$strict_cc,
                              "; supportive_bp=", .data$supportive_bp,
                              "; auto_include=", auto_include_syngo,
                              "; max_per_subpanel=", max_syngo_per_subpanel),
      confidence = dplyr::case_when(
        .data$strict_cc & .data$evidence_level == "traceable_experimental" ~ "external_syngo_strict_cc_experimental",
        .data$strict_cc ~ "external_syngo_strict_cc",
        .data$supportive_bp ~ "external_syngo_supportive_bp",
        TRUE ~ "external_syngo_candidate"
      ),
      use_for = "compartment_fidelity_candidate_audit",
      include_in_fidelity_score = as.character(auto_include_syngo & .data$strict_cc & .data$rank_within_subpanel <= max_syngo_per_subpanel),
      notes = "SynGO 1.3 candidate. Strict CC annotations are the preferred neuropil source; BP terms are supportive/audit-only unless manually curated.",
      stringsAsFactors = FALSE
    ) |>
    dplyr::distinct(.data$marker_set, .data$gene_symbol, .data$source_term_or_category, .keep_all = TRUE)
}

# -----------------------------------------------------------------------------
# Curated fidelity seed table
# -----------------------------------------------------------------------------

curated_fidelity_seed_registry <- function() {
  dplyr::bind_rows(
    make_registry_row(
      marker_set = "canonical_neuronal_synaptic_neuropil",
      fidelity_marker_class = "Neuropil markers",
      fidelity_subpanel = "presynaptic_vesicle_release",
      cell_class = "neuron",
      cell_state = "synaptic_neuropil",
      genes = c("Snap25", "Syp", "Syn1", "Syn2", "Vamp2", "Stx1a", "Stxbp1", "Syt1", "Sv2a", "Rab3a"),
      source_type = "curated_compartment_fidelity_seed",
      source_name = "manual_synaptic_seed",
      source_reference = "SynGO/GO/MGI/literature-supported synaptic vesicle and presynaptic release marker seed",
      source_term_or_category = "presynaptic vesicle/release",
      selection_rule = "curated conservative seed; not CON/RES/SUS-derived",
      confidence = "curated_conservative",
      use_for = "compartment_fidelity",
      include_in_fidelity_score = TRUE,
      notes = "Neuropil fidelity seed: presynaptic/synaptic vesicle/release proteins."
    ),
    make_registry_row(
      marker_set = "canonical_neuronal_synaptic_neuropil",
      fidelity_marker_class = "Neuropil markers",
      fidelity_subpanel = "postsynaptic_density",
      cell_class = "neuron",
      cell_state = "synaptic_neuropil",
      genes = c("Dlg4", "Dlg3", "Homer1", "Shank1", "Shank2", "Shank3", "Grin1", "Grin2a", "Grin2b", "Gria1", "Gria2", "Camk2a", "Camk2b"),
      source_type = "curated_compartment_fidelity_seed",
      source_name = "manual_psd_seed",
      source_reference = "SynGO/GO/MGI/literature-supported postsynaptic-density marker seed",
      source_term_or_category = "postsynaptic density/excitatory synapse",
      selection_rule = "curated conservative seed; not CON/RES/SUS-derived",
      confidence = "curated_conservative",
      use_for = "compartment_fidelity",
      include_in_fidelity_score = TRUE,
      notes = "Neuropil fidelity seed: postsynaptic-density and excitatory synapse proteins."
    ),
    make_registry_row(
      marker_set = "canonical_neuronal_synaptic_neuropil",
      fidelity_marker_class = "Neuropil markers",
      fidelity_subpanel = "neurite_cytoskeleton",
      cell_class = "neuron",
      cell_state = "synaptic_neuropil",
      genes = c("Nefl", "Nefm", "Nefh", "Gap43", "Stmn2", "Dpysl2", "Map2"),
      source_type = "curated_compartment_fidelity_seed",
      source_name = "manual_neurite_seed",
      source_reference = "GO/MGI/literature-supported neurite/axon/dendrite marker seed",
      source_term_or_category = "neurite/axon/dendrite cytoskeleton",
      selection_rule = "curated conservative seed; not CON/RES/SUS-derived",
      confidence = "curated_conservative",
      use_for = "compartment_fidelity",
      include_in_fidelity_score = TRUE,
      notes = "Neuropil supportive seed: neurite/axon/dendrite cytoskeleton."
    ),
    make_registry_row(
      marker_set = "canonical_neuronal_soma_nuclear",
      fidelity_marker_class = "Soma markers",
      fidelity_subpanel = "soma_nuclear_strict",
      cell_class = "neuron",
      cell_state = "soma_nuclear",
      genes = c("Rbfox3", "Matr3", "Srsf3", "Ddx39b", "Nono", "Sfpq", "Hnrnpa2b1", "Hnrnpc"),
      source_type = "curated_compartment_fidelity_seed",
      source_name = "manual_soma_nuclear_seed",
      source_reference = "GO/MGI/literature-supported nuclear/RNA-processing soma marker seed",
      source_term_or_category = "nuclear/RNA-processing/soma-associated",
      selection_rule = "curated conservative seed; not CON/RES/SUS-derived",
      confidence = "curated_conservative",
      use_for = "compartment_fidelity",
      include_in_fidelity_score = TRUE,
      notes = "Soma fidelity seed: nuclear/RNA-processing markers. Not formal soma purity estimation."
    ),
    make_registry_row(
      marker_set = "canonical_neuronal_soma_nuclear",
      fidelity_marker_class = "Soma markers",
      fidelity_subpanel = "chromatin_histone_supportive",
      cell_class = "neuron",
      cell_state = "soma_nuclear",
      genes = c("H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3", "H1f0"),
      source_type = "curated_compartment_fidelity_seed",
      source_name = "manual_chromatin_seed",
      source_reference = "GO/MGI/literature-supported chromatin/histone nuclear marker seed",
      source_term_or_category = "chromatin/histone/nuclear",
      selection_rule = "curated supportive seed; not CON/RES/SUS-derived",
      confidence = "curated_supportive",
      use_for = "compartment_fidelity",
      include_in_fidelity_score = TRUE,
      notes = "Soma supportive seed: chromatin/histone nuclear markers. Interpret cautiously because not neuron-specific."
    ),
    make_registry_row(
      marker_set = "canonical_neuronal_soma_nuclear",
      fidelity_marker_class = "Soma markers",
      fidelity_subpanel = "neuronal_soma_supportive",
      cell_class = "neuron",
      cell_state = "soma_nuclear_supportive",
      genes = c("Tubb3", "Map2"),
      source_type = "curated_compartment_fidelity_seed",
      source_name = "manual_neuronal_soma_supportive_seed",
      source_reference = "curated neuronal soma/cytoskeletal supportive markers",
      source_term_or_category = "neuronal soma supportive",
      selection_rule = "curated supportive seed; not CON/RES/SUS-derived",
      confidence = "curated_supportive",
      use_for = "compartment_fidelity",
      include_in_fidelity_score = TRUE,
      notes = "Supportive neuronal markers; not strict soma/nuclear markers."
    ),
    make_registry_row(
      marker_set = "canonical_microglia_homeostatic",
      fidelity_marker_class = "Microglia/PVM markers",
      fidelity_subpanel = "homeostatic_microglia",
      cell_class = "microglia",
      cell_state = "homeostatic_identity",
      genes = c("P2ry12", "Tmem119", "Cx3cr1", "Csf1r", "Hexb", "Fcrls", "Sall1", "Siglech", "Gpr34", "Mertk", "Aif1"),
      source_type = "curated_compartment_fidelity_seed",
      source_name = "manual_microglia_homeostatic_seed",
      source_reference = "curated homeostatic microglia marker seed; compatible with Allen-derived microglia/PVM reference panels",
      source_term_or_category = "homeostatic microglia",
      selection_rule = "curated conservative seed; not CON/RES/SUS-derived",
      confidence = "curated_conservative",
      use_for = "compartment_fidelity",
      include_in_fidelity_score = TRUE,
      notes = "Microglia identity markers. Proteomics detectability may be limited."
    ),
    make_registry_row(
      marker_set = "canonical_microglia_phagolysosomal_state",
      fidelity_marker_class = "Microglia/PVM markers",
      fidelity_subpanel = "phagolysosomal_complement_microglia",
      cell_class = "microglia",
      cell_state = "phagolysosomal_complement_state",
      genes = c("Tyrobp", "Trem2", "Apoe", "Lpl", "Cst7", "Ctsb", "Ctsd", "Lgals3", "Itgax", "Axl", "C1qa", "C1qb", "C1qc"),
      source_type = "curated_compartment_fidelity_seed",
      source_name = "manual_microglia_phagolysosomal_seed",
      source_reference = "curated phagolysosomal/complement microglia marker seed",
      source_term_or_category = "phagolysosomal/complement microglia",
      selection_rule = "curated conservative seed; not CON/RES/SUS-derived",
      confidence = "curated_conservative",
      use_for = "compartment_fidelity",
      include_in_fidelity_score = TRUE,
      notes = "Microglia state markers; not interpreted as pure identity markers."
    )
  )
}

# -----------------------------------------------------------------------------
# Download/cache GO/MGI and SynGO sources, then process candidates
# -----------------------------------------------------------------------------

download_logs <- list()

# GO/MGI cache.
go_mgi_raw <- path_external("reference_markers", "go_mgi", "raw")
go_mgi_processed <- path_external("reference_markers", "go_mgi", "processed")
safe_dir_create(go_mgi_raw)
safe_dir_create(go_mgi_processed)
mgi_gaf_file <- file.path(go_mgi_raw, "mgi.gaf.gz")
go_obo_file <- file.path(go_mgi_raw, "go-basic.obo")
mgi_homology_file <- file.path(go_mgi_raw, "HOM_MouseHumanSequence.rpt")

download_logs[["mgi_gaf"]] <- download_if_missing(mgi_gaf_url, mgi_gaf_file, allow_download, "GO_MGI_mgi_gaf")
download_logs[["go_obo"]] <- download_if_missing(go_basic_obo_url, go_obo_file, allow_download, "GO_basic_obo")
download_logs[["mgi_homology"]] <- download_if_missing(mgi_homology_url, mgi_homology_file, allow_download, "MGI_human_mouse_homology")

homology_map <- read_mgi_homology(mgi_homology_file)
if (!is.null(homology_map) && nrow(homology_map)) {
  write_csv_safe2(homology_map, file.path(go_mgi_processed, "mgi_human_mouse_homology_map.csv"))
}

go_mgi_candidates <- ensure_fidelity_candidate_cols(extract_go_mgi_candidates(mgi_gaf_file, go_obo_file))
if (nrow(go_mgi_candidates)) {
  write_csv_safe2(go_mgi_candidates, file.path(go_mgi_processed, "go_mgi_compartment_marker_candidates.csv"))
}

# SynGO cache. Prefer manifest source with source_name syngo/SynGO if present;
# otherwise use env vars. SynGO 1.3 complete export is expected as annotations.xlsx,
# genes.xlsx, and ontologies.xlsx. annotations.xlsx is the primary evidence table.
syngo_src <- find_source(c("syngo", "SynGO", "SYNGO"))
syngo_raw <- if (!is.null(syngo_src)) resolve_local_path(syngo_src$local_dir %||% path_external("reference_markers", "syngo", "raw")) else path_external("reference_markers", "syngo", "raw")
syngo_processed <- path_external("reference_markers", "syngo", "processed")
safe_dir_create(syngo_raw)
safe_dir_create(syngo_processed)

syngo_url <- syngo_src$url %||% syngo_url_env
syngo_file_name <- syngo_src$file_name %||% syngo_file_name_env
syngo_genes_file_name <- syngo_src$genes_file_name %||% syngo_genes_file_name_env
syngo_ontologies_file_name <- syngo_src$ontologies_file_name %||% syngo_ontologies_file_name_env

syngo_file <- find_existing_file(syngo_raw, c(syngo_file_name, "annotations.xlsx", "syngo_annotations.xlsx", "SynGO_annotations.xlsx", "annotations.tsv", "annotations.csv"))
syngo_genes_file <- find_existing_file(syngo_raw, c(syngo_genes_file_name, "genes.xlsx", "syngo_genes.xlsx", "SynGO_genes.xlsx"))
syngo_ontologies_file <- find_existing_file(syngo_raw, c(syngo_ontologies_file_name, "ontologies.xlsx", "syngo_ontologies.xlsx", "SynGO_ontologies.xlsx"))

if (nzchar(syngo_url) || is_regular_file(syngo_file)) {
  download_logs[["syngo_annotations"]] <- download_if_missing(syngo_url, syngo_file, allow_download, "SynGO_annotations")
} else {
  report_reference_marker_status("SynGO_annotations", "missing optional file; no download URL configured", syngo_file)
}
# genes.xlsx and ontologies.xlsx are metadata/support files. If they are not present,
# annotations.xlsx can still be parsed.
if (!is_regular_file(syngo_genes_file)) {
  report_reference_marker_status("SynGO_genes", "missing optional support file", syngo_genes_file)
  syngo_genes_file <- NA_character_
}
if (!is_regular_file(syngo_ontologies_file)) {
  report_reference_marker_status("SynGO_ontologies", "missing optional support file", syngo_ontologies_file)
  syngo_ontologies_file <- NA_character_
}

syngo_candidates <- ensure_fidelity_candidate_cols(
  extract_syngo_candidates(
    syngo_file,
    genes_file = syngo_genes_file,
    ontologies_file = syngo_ontologies_file,
    homology_map = homology_map,
    source_reference = if (nzchar(syngo_url)) syngo_url else "SynGO 1.3 complete data export cached locally"
  )
)
if (nrow(syngo_candidates)) {
  write_csv_safe2(syngo_candidates, file.path(syngo_processed, "syngo_compartment_marker_candidates.csv"))
}

# Consolidated download/cache manifest. This must exist even when no optional
# downloadable source was requested, because the source-status and run-manifest
# sections below append it if present.
download_log <- dplyr::bind_rows(download_logs)
if (is.null(download_log) || !is.data.frame(download_log)) {
  download_log <- data.frame(
    source_name = character(),
    url = character(),
    file = character(),
    status = character(),
    downloaded_at = character(),
    hash = character(),
    stringsAsFactors = FALSE
  )
}
if (nrow(download_log)) {
  write_csv_safe2(download_log, file.path(PATHS$tables, "reference_marker_download_manifest.csv"))
  write_csv_safe2(download_log, path_external("reference_markers", "reference_marker_download_manifest.csv"))
}

# External candidate registries. Manual fallback is decided after generic
# sources are imported, so Allen microglia/PVM can satisfy that class.
fidelity_seed <- curated_fidelity_seed_registry()
syngo_include_registry <- ensure_fidelity_candidate_cols(syngo_candidates) |>
  dplyr::filter(.data$include_in_fidelity_score %in% "TRUE") |>
  dplyr::mutate(
    # Use the existing recognized marker panel for compatibility with current 04c,
    # while preserving the subpanel/provenance columns.
    marker_set = "canonical_neuronal_synaptic_neuropil",
    use_for = "compartment_fidelity"
  )

go_mgi_evidence_keep <- parse_env_list(go_mgi_include_evidence)
go_mgi_subpanel_keep <- parse_env_list(go_mgi_include_subpanels)
go_mgi_include_registry <- ensure_fidelity_candidate_cols(go_mgi_candidates) |>
  dplyr::filter(
    auto_include_go_mgi,
    .data$evidence_level %in% go_mgi_evidence_keep,
    .data$fidelity_subpanel %in% go_mgi_subpanel_keep
  ) |>
  dplyr::mutate(
    evidence_priority = match(.data$evidence_level, go_mgi_evidence_keep),
    gene_token = normalize_gene_token(.data$gene_symbol)
  ) |>
  dplyr::filter(nzchar(.data$gene_token)) |>
  dplyr::arrange(.data$fidelity_marker_class, .data$fidelity_subpanel, .data$evidence_priority, .data$source_term_or_category, .data$gene_symbol) |>
  dplyr::distinct(.data$fidelity_marker_class, .data$fidelity_subpanel, .data$gene_token, .keep_all = TRUE) |>
  dplyr::group_by(.data$fidelity_marker_class, .data$fidelity_subpanel) |>
  dplyr::slice_head(n = max_go_mgi_per_subpanel) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    include_in_fidelity_score = "TRUE",
    use_for = "compartment_fidelity",
    selection_rule = paste0(
      "GO/MGI cellular component descendant candidate auto-included; evidence in ",
      paste(go_mgi_evidence_keep, collapse = "/"),
      "; max_per_subpanel=", max_go_mgi_per_subpanel
    ),
    confidence = paste0(.data$confidence, "_auto_included")
  ) |>
  dplyr::select(-"evidence_priority", -"gene_token")

go_mgi_include_keys <- go_mgi_include_registry |>
  dplyr::transmute(
    marker_set,
    gene_token = normalize_gene_token(.data$gene_symbol),
    source_term_or_category
  ) |>
  dplyr::distinct()

go_mgi_candidates_flagged <- ensure_fidelity_candidate_cols(go_mgi_candidates) |>
  dplyr::mutate(gene_token = normalize_gene_token(.data$gene_symbol)) |>
  dplyr::left_join(
    go_mgi_include_keys |> dplyr::mutate(go_mgi_auto_included = TRUE),
    by = c("marker_set", "gene_token", "source_term_or_category")
  ) |>
  dplyr::mutate(
    include_in_fidelity_score = ifelse(.data$go_mgi_auto_included %in% TRUE, "TRUE", .data$include_in_fidelity_score)
  ) |>
  dplyr::select(-"gene_token", -"go_mgi_auto_included")

# -----------------------------------------------------------------------------
# Generic external cell-type/reference marker import from manifest
# -----------------------------------------------------------------------------

source_status <- list()
external_evidence <- list()
for (src in sources) {
  src_name <- src$source_name %||% "unknown"
  src_dir <- resolve_local_path(src$local_dir %||% "")
  status <- "missing_optional"
  n_rows <- 0L
  if (identical(src$source_type, "package") && requireNamespace(src$package %||% "", quietly = TRUE)) {
    status <- "package_available_not_auto_extracted"
  }
  # SynGO is handled above as a compartment-fidelity source. Still allow generic
  # parsing if a compatible marker table is present, but do not double-count status.
  patterns <- unlist(src$file_patterns %||% c("*.csv", "*.tsv", "*.xlsx", "*.rds"), use.names = FALSE)
  files <- character()
  if (dir.exists(src_dir)) {
    regex <- paste0("(", paste(gsub("\\*", ".*", gsub("\\.", "\\\\.", patterns)), collapse = "|"), ")$")
    files <- list.files(src_dir, pattern = regex, full.names = TRUE, ignore.case = TRUE)
    files <- files[vapply(files, is_regular_file, logical(1))]
  }
  if (length(files)) {
    parsed <- lapply(files, function(f) coerce_marker_table(read_any_table(f), src))
    parsed <- parsed[!vapply(parsed, is.null, logical(1))]
    if (length(parsed)) {
      tab <- dplyr::bind_rows(parsed)
      external_evidence[[src_name]] <- tab
      n_rows <- nrow(tab)
      status <- "imported_cached_file"
    }
  }
  if (allow_download && identical(status, "missing_optional")) status <- "download_allowed_but_not_implemented"
  if ((src_name %in% c("syngo", "SynGO", "SYNGO")) && is_regular_file(syngo_file)) {
    status <- paste0(status, ";syngo_compartment_processed")
  }
  report_reference_marker_status(src_name, status, src_dir)
  source_status[[src_name]] <- data.frame(source_name = src_name, status = status, n_rows = n_rows, local_dir = src_dir, stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------
# Build final marker registry
# -----------------------------------------------------------------------------

fallback_registry <- dplyr::bind_rows(lapply(names(canonical_marker_panels), function(ms) {
  meta <- panel_meta[panel_meta$marker_set == ms, , drop = FALSE]
  data.frame(
    marker_set = ms,
    fidelity_marker_class = NA_character_,
    fidelity_subpanel = NA_character_,
    cell_class = meta$cell_class,
    cell_state = meta$cell_state,
    gene_symbol = canonical_marker_panels[[ms]],
    source_type = "curated_fallback",
    source_name = "proteomics_curated_fallback",
    source_reference = "03_qc_exploration/04b_import_reference_marker_sources.r curated conservative fallback",
    source_term_or_category = NA_character_,
    evidence_level = NA_character_,
    selection_rule = "conservative_curated_fallback; annotation_only",
    confidence = "curated_conservative",
    use_for = "annotation_reporting_sensitivity_interpretation",
    include_in_fidelity_score = "FALSE",
    notes = "Marker evidence only; not purity correction, not WGCNA filtering, not CON/RES/SUS-derived.",
    stringsAsFactors = FALSE
  )
}))

evidence <- dplyr::bind_rows(external_evidence)
if (!nrow(evidence)) {
  evidence <- fallback_registry |>
    dplyr::transmute(
      source_name,
      source_version = as.character(Sys.Date()),
      source_type,
      gene_symbol,
      cell_class,
      cell_state,
      reference_cell_label = cell_class,
      mean_target = NA_real_,
      mean_other = NA_real_,
      logFC_target_vs_other = NA_real_,
      specificity_score = NA_real_,
      pct_target = NA_real_,
      pct_other = NA_real_,
      rank_within_source = dplyr::row_number(),
      selection_rule,
      selected = TRUE,
      confidence,
      notes
    )
} else {
  evidence <- evidence |>
    dplyr::mutate(selected = .data$rank_within_source <= top_n)
}

registry_external <- if (nrow(dplyr::bind_rows(external_evidence))) {
  dplyr::bind_rows(external_evidence) |>
    dplyr::filter(.data$selected %in% TRUE) |>
    dplyr::transmute(
      marker_set = paste0("reference_", safe_filename(tolower(.data$cell_class))),
      fidelity_marker_class = NA_character_,
      fidelity_subpanel = NA_character_,
      cell_class,
      cell_state,
      gene_symbol,
      source_type,
      source_name,
      source_reference,
      source_term_or_category = NA_character_,
      evidence_level = NA_character_,
      selection_rule,
      confidence,
      use_for,
      include_in_fidelity_score = "FALSE",
      notes
    )
} else {
  data.frame()
}

allen_fidelity_registry <- if (nrow(registry_external)) {
  registry_external |>
    dplyr::filter(.data$marker_set == "reference_microglia_pvm") |>
    dplyr::mutate(
      fidelity_marker_class = "Microglia/PVM markers",
      fidelity_subpanel = "allen_microglia_pvm",
      use_for = "compartment_fidelity",
      include_in_fidelity_score = "TRUE",
      selection_rule = paste0(.data$selection_rule, "; external_preferred_fidelity_source")
    )
} else {
  empty_fidelity_candidate_table()
}

external_fidelity_registry <- if (identical(fidelity_source_policy, "manual_only")) {
  empty_fidelity_candidate_table()
} else {
  dplyr::bind_rows(syngo_include_registry, go_mgi_include_registry, allen_fidelity_registry)
}

external_fidelity_classes_available <- external_fidelity_registry |>
  dplyr::filter(!is.na(.data$fidelity_marker_class), nzchar(.data$fidelity_marker_class)) |>
  dplyr::pull(.data$fidelity_marker_class) |>
  unique()

manual_fidelity_registry <- dplyr::case_when(
  fidelity_source_policy == "manual_only" ~ "all",
  fidelity_source_policy == "manual_plus_external" ~ "all",
  fidelity_source_policy == "external_only" ~ "none",
  TRUE ~ "missing_classes_only"
)
manual_fidelity_registry <- if (identical(manual_fidelity_registry, "all")) {
  fidelity_seed
} else if (identical(manual_fidelity_registry, "missing_classes_only")) {
  fidelity_seed |>
    dplyr::filter(!.data$fidelity_marker_class %in% external_fidelity_classes_available)
} else {
  empty_fidelity_candidate_table()
}

manual_include_keys <- manual_fidelity_registry |>
  dplyr::transmute(
    marker_set,
    gene_token = normalize_gene_token(.data$gene_symbol),
    source_name,
    manual_included = TRUE
  ) |>
  dplyr::distinct()

fidelity_seed_flagged <- fidelity_seed |>
  dplyr::mutate(gene_token = normalize_gene_token(.data$gene_symbol)) |>
  dplyr::left_join(manual_include_keys, by = c("marker_set", "gene_token", "source_name")) |>
  dplyr::mutate(
    include_in_fidelity_score = ifelse(.data$manual_included %in% TRUE, "TRUE", "FALSE"),
    use_for = ifelse(.data$manual_included %in% TRUE, .data$use_for, "compartment_fidelity_fallback_not_used"),
    notes = ifelse(.data$manual_included %in% TRUE, .data$notes, paste0(.data$notes, " Manual fallback not used because external fidelity source is available."))
  ) |>
  dplyr::select(-"gene_token", -"manual_included")

fidelity_candidates <- dplyr::bind_rows(fidelity_seed_flagged, go_mgi_candidates_flagged, ensure_fidelity_candidate_cols(syngo_candidates)) |>
  dplyr::mutate(gene_token = normalize_gene_token(.data$gene_symbol)) |>
  dplyr::filter(nzchar(.data$gene_token)) |>
  dplyr::distinct(.data$marker_set, .data$gene_token, .data$source_name, .data$source_term_or_category, .keep_all = TRUE) |>
  dplyr::select(-"gene_token")

fidelity_registry <- dplyr::bind_rows(external_fidelity_registry, manual_fidelity_registry) |>
  dplyr::mutate(gene_token = normalize_gene_token(.data$gene_symbol)) |>
  dplyr::filter(nzchar(.data$gene_token)) |>
  dplyr::distinct(.data$marker_set, .data$gene_token, .keep_all = TRUE) |>
  dplyr::select(-"gene_token")

fidelity_class_source_policy <- dplyr::bind_rows(
  external_fidelity_registry |>
    dplyr::count(.data$fidelity_marker_class, .data$source_name, .data$source_type, name = "n_markers") |>
    dplyr::mutate(role = "external_scored"),
  fidelity_seed_flagged |>
    dplyr::count(.data$fidelity_marker_class, .data$source_name, .data$source_type, .data$include_in_fidelity_score, name = "n_markers") |>
    dplyr::mutate(role = ifelse(.data$include_in_fidelity_score %in% "TRUE", "manual_fallback_scored", "manual_fallback_available_not_used")) |>
    dplyr::select(-"include_in_fidelity_score")
) |>
  dplyr::mutate(policy = fidelity_source_policy) |>
  dplyr::arrange(.data$fidelity_marker_class, .data$role, .data$source_name)

write_table_and_source(fidelity_candidates, PATHS$tables, PATHS$source_data, "compartment_fidelity_marker_candidates.csv")
write_table_and_source(fidelity_registry, PATHS$tables, PATHS$source_data, "compartment_fidelity_marker_sets.csv")
write_csv_safe2(fidelity_registry, repo_path("config", "marker_panels", "compartment_fidelity_marker_sets.csv"))
write_csv_safe2(fidelity_class_source_policy, file.path(PATHS$tables, "reference_marker_fidelity_source_policy.csv"))

message("[reference markers] fidelity source policy: ", fidelity_source_policy)
message("[reference markers] external scored sources: ", paste(sort(unique(external_fidelity_registry$source_name)), collapse = ", "))
if (nrow(manual_fidelity_registry)) {
  message("[reference markers] manual fallback scored for classes: ", paste(sort(unique(manual_fidelity_registry$fidelity_marker_class)), collapse = ", "))
} else {
  message("[reference markers] manual fallback scored for classes: none")
}

fallback_fidelity_class <- function(marker_set) {
  marker_set_l <- tolower(as.character(marker_set))
  dplyr::case_when(
    marker_set_l %in% c("canonical_neuronal_soma_nuclear", "nuclear_soma") ~ "Soma markers",
    marker_set_l %in% c("canonical_neuronal_synaptic_neuropil", "neuronal_synaptic_neuropil") ~ "Neuropil markers",
    marker_set_l %in% c("canonical_microglia_homeostatic", "canonical_microglia_phagolysosomal_state", "reference_microglia_pvm", "microglia_homeostatic", "microglia_phagolysosomal") ~ "Microglia/PVM markers",
    TRUE ~ NA_character_
  )
}

fallback_registry <- fallback_registry |>
  dplyr::mutate(fallback_fidelity_marker_class = fallback_fidelity_class(.data$marker_set)) |>
  dplyr::filter(
    fidelity_source_policy %in% c("manual_only", "manual_plus_external") |
      is.na(.data$fallback_fidelity_marker_class) |
      !.data$fallback_fidelity_marker_class %in% external_fidelity_classes_available
  ) |>
  dplyr::select(-"fallback_fidelity_marker_class")

# Bind order matters: fidelity_registry comes first so overlapping genes in the
# recognized canonical soma/neuropil/microglia marker sets retain explicit
# fidelity metadata rather than generic fallback metadata.
registry <- dplyr::bind_rows(fidelity_registry, fallback_registry, registry_external) |>
  dplyr::mutate(gene_token = normalize_gene_token(.data$gene_symbol)) |>
  dplyr::filter(nzchar(.data$gene_token)) |>
  dplyr::distinct(.data$marker_set, .data$gene_token, .keep_all = TRUE) |>
  dplyr::select(-"gene_token")

ranked <- registry |>
  dplyr::group_by(.data$marker_set) |>
  dplyr::mutate(rank = dplyr::row_number(), n_markers = dplyr::n()) |>
  dplyr::ungroup()

# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------

write_table_and_source(evidence, PATHS$tables, PATHS$source_data, "reference_marker_evidence_long.csv")
write_table_and_source(ranked, PATHS$tables, PATHS$source_data, "reference_marker_sets_ranked.csv")
write_csv_safe2(registry, repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"))

status_tbl <- dplyr::bind_rows(source_status)
if (nrow(download_log)) {
  status_tbl <- dplyr::bind_rows(
    status_tbl,
    download_log |>
      dplyr::transmute(source_name, status, n_rows = NA_integer_, local_dir = dirname(.data$file))
  )
}

if (nrow(status_tbl)) {
  write_csv_safe2(status_tbl, file.path(PATHS$tables, "reference_marker_source_status.csv"))
  p_status <- ggplot2::ggplot(status_tbl, ggplot2::aes(x = .data$source_name, y = dplyr::coalesce(.data$n_rows, 0L), fill = .data$status)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Imported evidence rows", fill = "Status") +
    ggplot2::theme_classic(base_size = 8)
  ggplot2::ggsave(file.path(PATHS$figures, "reference_marker_source_summary.svg"), p_status, width = 150, height = 90, units = "mm", device = svglite::svglite)
}

candidate_source_summary <- fidelity_candidates |>
  dplyr::group_by(.data$source_name, .data$source_type, .data$fidelity_marker_class, .data$evidence_level, .data$include_in_fidelity_score) |>
  dplyr::summarise(
    n_candidate_rows = dplyr::n(),
    n_unique_gene_symbols = dplyr::n_distinct(.data$gene_symbol),
    n_marker_sets = dplyr::n_distinct(.data$marker_set),
    marker_sets = collapse_unique(.data$marker_set),
    source_terms_or_categories = collapse_unique(.data$source_term_or_category),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$source_name, .data$fidelity_marker_class, .data$include_in_fidelity_score)
if (nrow(candidate_source_summary)) {
  write_csv_safe2(candidate_source_summary, file.path(PATHS$tables, "compartment_fidelity_candidate_source_summary.csv"))
}

counts <- registry |> dplyr::count(.data$cell_class, name = "n_markers")
p_counts <- ggplot2::ggplot(counts, ggplot2::aes(x = .data$cell_class, y = .data$n_markers)) +
  ggplot2::geom_col(fill = "#2F6F73") +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Registry markers") +
  ggplot2::theme_classic(base_size = 8)
ggplot2::ggsave(file.path(PATHS$figures, "reference_marker_cellclass_counts.svg"), p_counts, width = 120, height = 90, units = "mm", device = svglite::svglite)

fid_counts <- fidelity_registry |>
  dplyr::count(.data$fidelity_marker_class, .data$fidelity_subpanel, name = "n_markers")
if (nrow(fid_counts)) {
  p_fid <- ggplot2::ggplot(fid_counts, ggplot2::aes(x = .data$fidelity_subpanel, y = .data$n_markers, fill = .data$fidelity_marker_class)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Fidelity markers", fill = "Fidelity class") +
    ggplot2::theme_classic(base_size = 8)
  ggplot2::ggsave(file.path(PATHS$figures, "compartment_fidelity_marker_counts.svg"), p_fid, width = 150, height = 95, units = "mm", device = svglite::svglite)
}

notes <- c(
  "# Reference marker import and compartment fidelity marker registry", "",
  "This step builds a normalized marker registry for downstream QC/interpreter support.", "",
  "Key outputs:",
  "- config/marker_panels/wgcna_reference_marker_sets.csv: combined marker registry used by 04c.",
  "- config/marker_panels/compartment_fidelity_marker_sets.csv: conservative soma/neuropil/microglia-PVM fidelity markers.",
  "- compartment_fidelity_marker_candidates.csv: broader GO/MGI and SynGO candidate audit table.", "",
  "- reference_marker_source_status.csv: cache/import status for each manifest/download source.",
  "- reference_marker_fidelity_source_policy.csv: straight source-policy readout showing external scored sources and manual fallback use.",
  "- compartment_fidelity_candidate_source_summary.csv: GO/MGI, SynGO, and curated candidate counts by source and inclusion flag.", "",
  "Interpretation constraints:",
  "- Fidelity markers are not formal purity/deconvolution markers.",
  "- Neuropil markers are based on synaptic/PSD/vesicle/neurite evidence rather than broad neuronal cell-type references.",
  "- Soma markers are labelled soma/nuclear and include nuclear/RNA-processing/chromatin/supportive neuronal proteins.",
  "- Microglia markers are labelled Microglia/PVM when PVM-compatible reference panels are used downstream.",
  "- Default fidelity source policy is external_preferred: Allen, GO/MGI, and SynGO are used when available; manual curated seeds are scored only for missing fidelity classes.",
  "- GO/MGI broad candidates are auto-included by default using selected evidence codes/subpanels; set PROTEOMICS_FIDELITY_AUTO_INCLUDE_GO_MGI=false to keep them audit-only.",
  "- SynGO 1.3 complete export is parsed from annotations.xlsx as the primary evidence table; genes.xlsx and ontologies.xlsx are treated as support/provenance files when present.",
  "- SynGO candidates are auto-included by default, and only strict CC localization rows are eligible by default.",
  "- Marker selection is source/provenance based and not CON/RES/SUS-derived.", "",
  paste0("Live download allowed: ", allow_download),
  paste0("Fidelity source policy: ", fidelity_source_policy),
  paste0("SynGO auto-include: ", auto_include_syngo),
  paste0("GO/MGI auto-include: ", auto_include_go_mgi),
  paste0("GO/MGI auto-include evidence codes: ", paste(go_mgi_evidence_keep, collapse = ", ")),
  paste0("GO/MGI auto-include subpanels: ", paste(go_mgi_subpanel_keep, collapse = ", ")),
  paste0("Fidelity markers in final table: ", nrow(fidelity_registry)),
  paste0("GO/MGI markers auto-included in final table: ", nrow(go_mgi_include_registry)),
  paste0("GO/MGI candidates: ", nrow(go_mgi_candidates)),
  paste0("SynGO candidates: ", nrow(syngo_candidates))
)
writeLines(notes, file.path(PATHS$reports, "reference_marker_import_interpretation_notes.md"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(
    manifest = manifest_file,
    source_dirs = vapply(sources, function(x) resolve_local_path(x$local_dir %||% ""), character(1)),
    mgi_gaf = mgi_gaf_file,
    go_basic_obo = go_obo_file,
    mgi_homology = mgi_homology_file,
    syngo_annotations = syngo_file,
    syngo_genes = syngo_genes_file,
    syngo_ontologies = syngo_ontologies_file
  ),
  outputs = list(
    evidence = file.path(PATHS$tables, "reference_marker_evidence_long.csv"),
    ranked = file.path(PATHS$tables, "reference_marker_sets_ranked.csv"),
    registry = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"),
    fidelity_markers = repo_path("config", "marker_panels", "compartment_fidelity_marker_sets.csv"),
    fidelity_candidates = file.path(PATHS$tables, "compartment_fidelity_marker_candidates.csv"),
    fidelity_candidate_source_summary = file.path(PATHS$tables, "compartment_fidelity_candidate_source_summary.csv"),
    download_manifest = file.path(PATHS$tables, "reference_marker_download_manifest.csv"),
    source_status = file.path(PATHS$tables, "reference_marker_source_status.csv"),
    fidelity_source_policy = file.path(PATHS$tables, "reference_marker_fidelity_source_policy.csv"),
    figures = PATHS$figures
  ),
  parameters = list(
    allow_download = allow_download,
    fidelity_source_policy = fidelity_source_policy,
    top_n_per_class = top_n,
    min_specificity = min_specificity,
    auto_include_syngo = auto_include_syngo,
    max_syngo_per_subpanel = max_syngo_per_subpanel,
    auto_include_go_mgi = auto_include_go_mgi,
    max_go_mgi_per_subpanel = max_go_mgi_per_subpanel,
    go_mgi_include_evidence = go_mgi_evidence_keep,
    go_mgi_include_subpanels = go_mgi_subpanel_keep,
    mgi_gaf_url = mgi_gaf_url,
    go_basic_obo_url = go_basic_obo_url,
    mgi_homology_url = mgi_homology_url,
    syngo_url = syngo_url,
    syngo_file_name = syngo_file_name,
    syngo_genes_file_name = syngo_genes_file_name,
    syngo_ontologies_file_name = syngo_ontologies_file_name
  ),
  notes = "Reference marker evidence is cached/versioned/optional and used for annotation/reporting/sensitivity only. Compartment fidelity markers are explicit QC/interpreter markers, not purity/deconvolution estimates."
)

message("Reference marker import complete. Registry: ", repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv"))
message("Compartment fidelity markers: ", repo_path("config", "marker_panels", "compartment_fidelity_marker_sets.csv"))
