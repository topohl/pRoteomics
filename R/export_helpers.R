# Shared helpers for 09_export_pride_journal (manifest-driven, pg_matrix-onward).

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("validate_dataset", mode = "function")) {
  source(repo_path("R", "dataset_config.R"))
}
if (!exists("read_sample_metadata", mode = "function")) {
  source(repo_path("R", "pride_helpers.R"))
}
if (!exists("resolve_dataset_inputs", mode = "function")) {
  source(repo_path("R", "dataset_inputs.R"))
}

export_config_path <- function() {
  repo_path("09_export_pride_journal", "config", "export_config.yml")
}

load_export_config <- function(path = export_config_path()) {
  if (!file.exists(path)) {
    stop("Export config not found: ", path, call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read export config.", call. = FALSE)
  }
  cfg <- yaml::read_yaml(path)
  cfg$config_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  cfg
}

export_cli_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  value_after <- function(flag, default = "") {
    hit <- which(args == flag)
    if (!length(hit) || hit[[1]] == length(args)) return(default)
    args[[hit[[1]] + 1]]
  }
  has <- function(flag) flag %in% args
  list(
    dataset = value_after("--dataset", default = Sys.getenv("PROTEOMICS_DATASET", unset = "all")),
    export_level = value_after("--export-level", default = "pg_matrix_onward"),
    include_derived_results = has("--include-derived-results"),
    recursive = has("--recursive"),
    dry_run = has("--dry-run") || is_dry_run(),
    skip_validation = has("--skip-validation"),
    skip_supplementary = has("--skip-supplementary"),
    skip_manuscript = has("--skip-manuscript"),
    skip_claims = has("--skip-claims")
  )
}

resolve_export_datasets <- function(dataset_arg, config) {
  if (is.null(dataset_arg) || !nzchar(dataset_arg) || identical(tolower(dataset_arg), "all")) {
    return(config$datasets)
  }
  validate_dataset(dataset_arg, source = "--dataset")
}

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("openssl", quietly = TRUE)) {
    hash <- openssl::sha256(file(path))
    return(paste0(format(hash), collapse = ""))
  }
  unname(tools::md5sum(path))
}

list_files_nonrecursive <- function(dir, pattern = NULL) {
  if (!dir.exists(dir)) return(character())
  list.files(dir, pattern = pattern, full.names = TRUE, recursive = FALSE, all.files = FALSE, no.. = TRUE)
}

list_files_shallow <- function(dir, pattern = NULL, max_depth = 2L) {
  if (!dir.exists(dir)) return(character())
  if (max_depth <= 0L) return(character())
  hits <- list_files_nonrecursive(dir, pattern)
  if (max_depth > 1L) {
    subdirs <- list.dirs(dir, full.names = TRUE, recursive = FALSE)
    subdirs <- subdirs[basename(subdirs) != basename(dir)]
    deeper <- unlist(lapply(subdirs, list_files_shallow, pattern = pattern, max_depth = max_depth - 1L), use.names = FALSE)
    hits <- unique(c(hits, deeper))
  }
  hits[file.exists(hits) & !file.info(hits)$isdir]
}

glob_paths <- function(patterns, root = repo_root()) {
  patterns <- unique(as.character(patterns))
  patterns <- patterns[nzchar(patterns)]
  if (!length(patterns)) return(character())
  hits <- unlist(lapply(patterns, function(p) {
    full <- if (grepl("^([A-Za-z]:|/|~)", p)) normalizePath(path.expand(p), winslash = "/", mustWork = FALSE) else repo_path(p)
    if (file.exists(full)) return(full)
    Sys.glob(full)
  }), use.names = FALSE)
  hits <- unique(normalizePath(hits[hits != "" & file.exists(hits)], winslash = "/", mustWork = FALSE))
  hits[!file.info(hits)$isdir]
}

normalize_export_path <- function(path) {
  gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = FALSE))
}

is_noncanonical_ewce_export_path <- function(path) {
  normalized <- normalize_export_path(path)
  grepl(
    "/05_celltype_enrichment_EWCE/EWCE_E9_comparison/",
    normalized,
    fixed = TRUE
  )
}

is_exportable_result_path <- function(path) {
  normalized <- normalize_export_path(path)
  # Regex, not fixed: the trailing (/|$) alternation must be interpreted.
  # With fixed = TRUE this pattern never matched, so the manuscript payload
  # was silently exportable and a rebuild could re-ingest its stale contents.
  !grepl("/results/manuscript(/|$)", normalized) &
    !is_noncanonical_ewce_export_path(normalized)
}

# Legacy dataset-suffixed figure aliases.
#
# Current WGCNA figure code writes an unsuffixed canonical name. Older runs
# also wrote "<stem>_<dataset>.<ext>" aliases. Where both survive in the same
# directory the suffixed file is a stale leftover of a superseded run and must
# not reach the manuscript payload. Removing it from the *export selection*
# leaves the scientific figure tree untouched.
#
# The canonical sibling must exist for a file to count as an alias, so
# descriptive names that merely end in a dataset word (for example
# "..._fingerprint_with_soma.svg", which has no "..._fingerprint_with.svg"
# sibling) are deliberately kept.
legacy_dataset_figure_suffixes <- function() {
  c("neuron_neuropil", "neuron_soma", "neuropil", "soma", "microglia")
}

is_legacy_dataset_suffixed_alias <- function(paths) {
  vapply(paths, function(path) {
    ext <- tools::file_ext(path)
    if (!nzchar(ext)) return(FALSE)
    stem <- tools::file_path_sans_ext(basename(path))
    directory <- dirname(path)
    for (suffix in legacy_dataset_figure_suffixes()) {
      tail <- paste0("_", suffix)
      if (!endsWith(stem, tail)) next
      canonical_stem <- substr(stem, 1L, nchar(stem) - nchar(tail))
      if (!nzchar(canonical_stem)) next
      canonical <- file.path(directory, paste0(canonical_stem, ".", ext))
      if (file.exists(canonical)) return(TRUE)
    }
    FALSE
  }, logical(1), USE.NAMES = FALSE)
}

drop_legacy_dataset_suffixed_aliases <- function(paths) {
  if (!length(paths)) return(paths)
  paths[!is_legacy_dataset_suffixed_alias(paths)]
}

# Figure families with no current producer.
#
# Entries are exact stems that were each individually confirmed to be written
# by no script, config or contract in this repository -- i.e. leftovers of
# removed code. They are enumerated deliberately rather than inferred:
# "basename does not appear in the source tree" is NOT a safe rule, because a
# figure name can legitimately be assembled at runtime. Add a stem here only
# after confirming for that exact family that no producer exists.
#
# all_supermodule_eigengene_spatial_group_plot: the stem appears nowhere in
# R/, the numbered stage directories, tools/, pipeline.yml, config/ or docs/.
# Four files survive on disk (an unsuffixed microglia file plus _microglia,
# _neuropil and _soma variants) from a superseded naming scheme. The whole
# family is excluded, including the unsuffixed member: none of it is a
# current canonical output, so keeping any member would ship an orphan.
orphan_figure_family_stems <- function() {
  c("all_supermodule_eigengene_spatial_group_plot")
}

# Manuscript figure target naming, budgeted against MAX_PATH.
#
# Windows MAX_PATH is 260 characters and file.copy() merely returns FALSE beyond
# it, so flattening a deep source-relative path into a single filename could
# silently fail to copy (measured: every target at 260+ chars failed, every
# target at 259 or fewer succeeded). Targets are therefore built against a
# conservative full-path budget well under the limit.
#
# When a name does not fit, the readable source-relative prefix is truncated and
# a short deterministic hash of the COMPLETE source-relative path is appended.
# Hashing the complete path (never the truncated form) is what keeps two sources
# that share a truncated prefix distinct.
#
# This is deliberately export-local: it does not use \\?\ long-path prefixes and
# does not alter the global safe_filename(), which is called here with an
# effectively unlimited max_chars so that its own truncation cannot collapse two
# names before the hash is computed.
manuscript_figure_path_budget <- function() 240L

manuscript_figure_target_hash <- function(rel_path) {
  substr(
    digest::digest(enc2utf8(as.character(rel_path)), algo = "sha1", serialize = FALSE),
    1L, 10L
  )
}

# `identity` is the string the disambiguating hash is taken over. It defaults to
# rel_paths, which is what the figure export needs and preserves its behaviour
# exactly. The table export passes the full root-qualified relative path so that
# two files sharing a within-root relative path (one under results/source_data,
# one under results/tables) can never derive the same hash.
manuscript_figure_target_paths <- function(rel_paths, target_dir,
                                           budget = manuscript_figure_path_budget(),
                                           identity = rel_paths) {
  if (!length(rel_paths)) return(character(0))
  if (length(identity) != length(rel_paths)) {
    stop("manuscript_figure_target_paths: identity and rel_paths differ in length.",
         call. = FALSE)
  }
  target_dir <- normalize_export_path(target_dir)
  prefix_chars <- nchar(target_dir) + 1L
  rel_paths <- as.character(rel_paths)
  identity <- as.character(identity)

  vapply(seq_along(rel_paths), function(i) {
    rel <- rel_paths[[i]]
    ext <- tolower(tools::file_ext(rel))
    stem <- safe_filename(tools::file_path_sans_ext(rel), max_chars = .Machine$integer.max)
    room <- budget - prefix_chars - (nchar(ext) + 1L)
    hash <- manuscript_figure_target_hash(identity[[i]])
    min_room <- nchar(hash) + 9L
    if (room < min_room) {
      stop(
        "Manuscript figure target directory is too deep for the ", budget,
        "-character path budget: ", target_dir,
        " leaves only ", room, " characters for a filename (need at least ",
        min_room, ").", call. = FALSE
      )
    }
    name <- if (nchar(stem) <= room) {
      paste0(stem, ".", ext)
    } else {
      paste0(substr(stem, 1L, room - nchar(hash) - 1L), "_", hash, ".", ext)
    }
    file.path(target_dir, name)
  }, character(1), USE.NAMES = FALSE)
}

# Fail closed on duplicate targets: two sources copied to one target would
# silently overwrite, leaving a payload that looks complete but is not.
assert_unique_export_targets <- function(targets, sources = NULL) {
  dup_targets <- unique(targets[duplicated(targets)])
  if (length(dup_targets)) {
    detail <- vapply(head(dup_targets, 5L), function(t) {
      src <- if (is.null(sources)) character(0) else sources[targets == t]
      paste0(basename(t), " <- ", paste(basename(src), collapse = ", "))
    }, character(1))
    stop(
      "Export target collision: ", length(dup_targets),
      " target path(s) claimed by more than one source. ",
      paste(detail, collapse = " | "), call. = FALSE
    )
  }
  invisible(TRUE)
}

# Copy and verify. Any FALSE return or any absent target aborts before the
# caller can write a manifest, so a partial payload is never reported as good.
copy_export_targets <- function(sources, targets) {
  if (!length(sources)) return(invisible(TRUE))
  if (length(sources) != length(targets)) {
    stop("copy_export_targets: sources and targets differ in length.", call. = FALSE)
  }
  copied <- file.copy(sources, targets, overwrite = TRUE)
  failed <- which(!copied)
  if (length(failed)) {
    stop(
      "Export copy failed for ", length(failed), " of ", length(sources),
      " file(s); refusing to write a manifest for a partial payload. First: ",
      paste(sprintf("%s (target %d chars)", basename(targets[failed]),
                    nchar(targets[failed]))[seq_len(min(3L, length(failed)))],
            collapse = " | "), call. = FALSE
    )
  }
  absent <- which(!file.exists(targets))
  if (length(absent)) {
    stop(
      "Export copy reported success but ", length(absent),
      " target(s) are absent: ",
      paste(basename(targets[absent])[seq_len(min(3L, length(absent)))],
            collapse = " | "), call. = FALSE
    )
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Manuscript / journal source-data scope.
#
# SCOPE OF THESE RULES: they filter the *manuscript source-data export* only
# (09_export_source_data.R). They never delete or alter an analysis output, and
# they are deliberately NOT applied to PRIDE/deposition packaging, whose
# selectors live in R/pride_helpers.R and remain untouched.
#
# Every rule names an explicit path or family and records why it is diagnostic,
# superseded, proposed or intermediate. There is intentionally NO size-based,
# timestamp-pattern, "proposed appears in the name" or "unreferenced in code"
# heuristic: each of those would fail open on legitimate manuscript source data.
# Where an exclusion is justified by the existence of a canonical sibling, that
# sibling is checked and its absence is a hard error rather than a silent drop.
# ---------------------------------------------------------------------------

# Repo-relative, forward-slashed path used for matching ("results/...").
source_data_scope_relpath <- function(paths) {
  normalized <- gsub("\\\\", "/", as.character(paths))
  sub("^.*/(results/)", "\\1", normalized)
}

# (1) DIAGNOSTIC: internal per-term GSEA direction audits. Manuscript tables
#     cite enrichment results, never this per-term direction bookkeeping.
source_data_excluded_diagnostic <- function(rel) {
  grepl(paste0(
    "^results/tables/04_differential_expression_enrichment/",
    "01b_gsea_protein_direction_audit/.*/BP/gsea_term_direction_audit\\.csv$"
  ), rel)
}

# (2) SUPERSEDED REMOVED PIPELINE: outputs of
#     04_neuropil_contamination_annotation.r, which docs/NAMING_MIGRATION.md
#     records as removed and replaced by 04_neuropil_reference_annotation.r.
#     Scoped to this one tree; not generalised.
source_data_excluded_removed_pipeline <- function(rel) {
  startsWith(rel, paste0(
    "results/tables/04_differential_expression_enrichment/",
    "neuropil_contamination_annotation/"
  ))
}

# (3) SUPERSEDED TIMESTAMPED VARIANTS: only the two named
#     neuropil_reference_annotation stems, and only when the explicit canonical
#     "<stem>_latest.csv" pointer exists in the same directory. No generic
#     "timestamped filename is superseded" rule.
source_data_timestamped_variant_stems <- function() {
  c("microglia_neuropil_annotation", "microglia_neuropil_annotation_summary")
}

source_data_excluded_timestamped_variant <- function(rel, results_root = path_results()) {
  family_dir <- paste0(
    "results/tables/04_differential_expression_enrichment/",
    "neuropil_reference_annotation/"
  )
  vapply(rel, function(one) {
    if (!startsWith(one, family_dir)) return(FALSE)
    base <- basename(one)
    for (stem in source_data_timestamped_variant_stems()) {
      pattern <- paste0("^", stem, "_[0-9]{8}_[0-9]{6}\\.csv$")
      if (!grepl(pattern, base)) next
      canonical_rel <- file.path(dirname(one), paste0(stem, "_latest.csv"))
      canonical <- file.path(dirname(results_root), canonical_rel)
      # Fail closed by KEEPING. The exclusion is conditional on an explicit
      # "<stem>_latest.csv" pointer, so without one these timestamped files may
      # be the only copy of that output and are retained. This is a normal,
      # permanent state for some stems (microglia_neuropil_annotation has a
      # _latest pointer; microglia_neuropil_annotation_summary does not), so it
      # must not abort the export.
      if (!file.exists(canonical)) return(FALSE)
      return(TRUE)
    }
    FALSE
  }, logical(1), USE.NAMES = FALSE)
}

# (4) PROPOSED / NON-CANONICAL VARIANTS: two explicitly named trees, excluded
#     only because a canonical sibling tree exists. Listed as
#     (proposed prefix, canonical prefix) pairs relative to the repository.
source_data_proposed_tree_pairs <- function() {
  list(
    c("results/tables/04_differential_expression_enrichment/compareGO_spatial_atlas_validation_proposed/",
      "results/tables/04_differential_expression_enrichment/compareGO_spatial_atlas"),
    c("results/source_data/04_differential_expression_enrichment/compareGO_spatial_atlas_validation_proposed/",
      "results/source_data/04_differential_expression_enrichment/compareGO_spatial_atlas"),
    c("results/tables/04_differential_expression_enrichment/microglia_targeted_signature_enrichment/microglia_validation_proposed/",
      "results/tables/04_differential_expression_enrichment/microglia_targeted_signature_enrichment/microglia"),
    c("results/source_data/04_differential_expression_enrichment/microglia_targeted_signature_enrichment/microglia_validation_proposed/",
      "results/source_data/04_differential_expression_enrichment/microglia_targeted_signature_enrichment/microglia")
  )
}

source_data_excluded_proposed_tree <- function(rel, results_root = path_results()) {
  repo <- dirname(results_root)
  out <- logical(length(rel))
  for (pair in source_data_proposed_tree_pairs()) {
    hit <- startsWith(rel, pair[[1]])
    if (!any(hit)) next
    canonical <- file.path(repo, pair[[2]])
    # Fail closed: the exclusion is only justified by the canonical tree.
    if (!dir.exists(canonical)) {
      stop(
        "Refusing to exclude proposed source-data tree without its canonical ",
        "sibling: ", pair[[1]], " (expected ", pair[[2]], ")", call. = FALSE
      )
    }
    out <- out | hit
  }
  out
}

# (5) INTERMEDIATE TERM-GENE EXPANSIONS: per-term gene provenance/evidence
#     expansions. They stay in analysis results for reproducibility but are not
#     per-figure journal source data.
source_data_excluded_term_gene_expansion <- function(rel) {
  provenance <- grepl(paste0(
    "^results/tables/04_differential_expression_enrichment/compareGO/[^/]+/BP/",
    "phenotype_within_unit/all_route_units/compareGO_term_gene_provenance\\.csv$"
  ), rel)
  evidence <- grepl(paste0(
    "^results/source_data/04_differential_expression_enrichment/",
    "biological_program_summary/(neuron_neuropil|neuron_soma|microglia)/",
    "program_term_gene_evidence\\.csv$"
  ), rel)
  provenance | evidence
}

# (6) STAGE-11 INTERMEDIATE: the joint control-pair GO FDR expansion. The
#     derived Stage-11 manuscript summaries, geometry tables, trajectories and
#     provenance are kept.
source_data_excluded_stage11_intermediate <- function(rel) {
  basename(rel) == "control_pair_joint_GO_FDR.csv" &
    grepl("/stress_response_biological_audit/", rel, fixed = TRUE)
}

# (7) CONCORDANCE INTERMEDIATE: the all-contrasts theme-assignment expansion,
#     in both the tables and source_data copies. Compact Figure-3/concordance
#     manuscript tables are kept.
source_data_excluded_concordance_intermediate <- function(rel) {
  basename(rel) == "ontology_aware_gsea_theme_assignments_all_contrasts.csv" &
    grepl("/gsea_wgcna_concordance/", rel, fixed = TRUE)
}

# (8) INTERMEDIATE LEGACY REPLAY: one exact prefix. This tree is a legacy
#     replay / validation output set, not a manuscript source-data family, and
#     no required manuscript-facing keep artifact lives in it. It also holds the
#     nine sources whose paths are exactly 260 characters and are therefore
#     unreadable through normal R file operations on Windows.
#
#     Scoped to this single prefix on purpose. It is NOT a rule about "legacy"
#     appearing in a name, nor about comparison trees, nor about sensitivity
#     analyses: the sibling animal_level/, repro1/, repro2/, rnga/ and rngb/
#     trees under the same comparison root remain selected, and a test pins
#     that.
source_data_legacy_replay_prefix <- function() {
  paste0(
    "results/source_data/04_differential_expression_enrichment_comparison/",
    "legacy_replay/"
  )
}

source_data_excluded_legacy_replay <- function(rel) {
  startsWith(rel, source_data_legacy_replay_prefix())
}

source_data_scope_exclusion_reasons <- function(paths, results_root = path_results()) {
  rel <- source_data_scope_relpath(paths)
  reason <- rep(NA_character_, length(rel))
  assign_reason <- function(reason, hit, label) {
    take <- hit & is.na(reason)
    reason[take] <- label
    reason
  }
  reason <- assign_reason(reason, source_data_excluded_diagnostic(rel), "diagnostic")
  reason <- assign_reason(reason, source_data_excluded_legacy_replay(rel), "intermediate_legacy_replay")
  reason <- assign_reason(reason, source_data_excluded_removed_pipeline(rel), "superseded_removed_pipeline")
  reason <- assign_reason(reason, source_data_excluded_timestamped_variant(rel, results_root), "superseded_timestamped_variant")
  reason <- assign_reason(reason, source_data_excluded_proposed_tree(rel, results_root), "proposed_non_canonical")
  reason <- assign_reason(reason, source_data_excluded_term_gene_expansion(rel), "intermediate_term_gene_expansion")
  reason <- assign_reason(reason, source_data_excluded_stage11_intermediate(rel), "intermediate_stage11")
  reason <- assign_reason(reason, source_data_excluded_concordance_intermediate(rel), "intermediate_concordance")
  reason
}

# Manuscript source-data / supplementary-table target naming.
#
# Same contract as the figure targets: a <= 240 character complete target path,
# no \\?\ prefixes, the global safe_filename() untouched and called with an
# unlimited max_chars so its own truncation cannot collapse two names before the
# hash decision, a readable source-relative prefix, and a deterministic
# 10-character hash of the COMPLETE root-qualified source identity when
# shortening is required. The real extension is preserved.
#
# Candidates arrive from two roots (results/source_data and results/tables) and
# are routed to two different output directories. Routing and naming are derived
# from the SAME repo-relative prefix so they cannot disagree, and the hash
# identity includes the root, so two files sharing a within-root relative path
# stay distinct.
manuscript_table_export_root <- function(rel_full) {
  ifelse(
    startsWith(rel_full, "results/source_data/"), "source_data",
    ifelse(startsWith(rel_full, "results/tables/"), "tables", NA_character_)
  )
}

manuscript_table_target_paths <- function(paths, source_data_dir, supplementary_dir,
                                          budget = manuscript_figure_path_budget()) {
  if (!length(paths)) return(character(0))
  rel_full <- source_data_scope_relpath(paths)
  root <- manuscript_table_export_root(rel_full)
  unknown <- which(is.na(root))
  if (length(unknown)) {
    stop(
      "Cannot route ", length(unknown), " source-data candidate(s): expected a ",
      "results/source_data/ or results/tables/ prefix. First: ",
      rel_full[[unknown[[1]]]], call. = FALSE
    )
  }
  within <- ifelse(root == "source_data",
                   sub("^results/source_data/", "", rel_full),
                   sub("^results/tables/", "", rel_full))
  out <- character(length(paths))
  for (r in c("source_data", "tables")) {
    idx <- which(root == r)
    if (!length(idx)) next
    out[idx] <- manuscript_figure_target_paths(
      within[idx],
      if (r == "source_data") source_data_dir else supplementary_dir,
      budget = budget,
      identity = rel_full[idx]
    )
  }
  out
}

# Fail closed on a target claimed by more than one source, and name the roots
# involved so a cross-root clash is obvious rather than silently resolved by a
# new naming policy.
assert_unique_table_export_targets <- function(targets, sources) {
  dup <- unique(targets[duplicated(targets)])
  if (!length(dup)) return(invisible(TRUE))
  detail <- vapply(head(dup, 5L), function(t) {
    src <- source_data_scope_relpath(sources[targets == t])
    paste0(basename(t), " <- ", paste(src, collapse = " + "))
  }, character(1))
  stop(
    "Source-data export target collision: ", length(dup),
    " target path(s) claimed by more than one source. ",
    paste(detail, collapse = " | "), call. = FALSE
  )
}

# Every selected source must be readable through ordinary R file operations.
# Runs after journal-scope filtering so excluded trees cannot block the export.
assert_export_sources_accessible <- function(paths) {
  if (!length(paths)) return(invisible(TRUE))
  info <- file.info(paths)
  bad <- which(!file.exists(paths) | is.na(info$size))
  if (length(bad)) {
    stop(
      length(bad), " selected source(s) are inaccessible (missing, or ",
      "unreadable file.info -- typically a path at or beyond the ",
      "260-character Windows limit). First: ",
      paste(sprintf("%s (%d chars)", source_data_scope_relpath(paths[bad]),
                    nchar(paths[bad]))[seq_len(min(3L, length(bad)))],
            collapse = " | "), call. = FALSE
    )
  }
  invisible(TRUE)
}

is_manuscript_source_data_path <- function(paths, results_root = path_results()) {
  is.na(source_data_scope_exclusion_reasons(paths, results_root))
}

apply_manuscript_source_data_scope <- function(paths, results_root = path_results()) {
  if (!length(paths)) return(paths)
  paths[is_manuscript_source_data_path(paths, results_root)]
}

is_orphan_figure_family_path <- function(paths) {
  families <- orphan_figure_family_stems()
  vapply(paths, function(path) {
    stem <- tools::file_path_sans_ext(basename(path))
    any(stem == families | startsWith(stem, paste0(families, "_")))
  }, logical(1), USE.NAMES = FALSE)
}

drop_orphan_figure_families <- function(paths) {
  if (!length(paths)) return(paths)
  paths[!is_orphan_figure_family_path(paths)]
}

canonical_ewce_figure_root <- function() {
  path_results("figures", "05_celltype_enrichment_EWCE", "EWCE_E9")
}

pg_matrix_input_paths <- function(config) {
  glob_paths(config$canonical_inputs$pg_matrix$paths)
}

metadata_input_paths <- function(config) {
  glob_paths(config$canonical_inputs$sample_metadata$paths)
}

dataset_matrix_patterns <- function(dataset) {
  dataset <- validate_dataset(dataset)
  c(
    sprintf("^\\d{8}_pgmatrix_imputed_%s_[0-9]+samples_missing70pct(\\.xlsx|_with_metadata\\.xlsx)$", dataset),
    sprintf("^\\d{8}_pgmatrix_imputed_%s_[0-9]+samples_missing70pct_with_metadata\\.gct$", dataset)
  )
}

processed_files_for_dataset <- function(dataset, config, include_derived = FALSE) {
  dataset <- validate_dataset(dataset)
  files <- character()

  prep_root <- repo_path(config$canonical_inputs$processed_preprocessing$root)
  subdirs <- config$canonical_inputs$processed_preprocessing$per_dataset_subdirs
  for (sub in subdirs) {
    base <- file.path(prep_root, sub)
    if (sub %in% c("gct_extractR", "protigy_output")) {
      ds_dir <- file.path(base, dataset)
      if (dir.exists(ds_dir)) {
        files <- c(files, list_files_shallow(ds_dir, max_depth = 3L))
      }
      next
    }
    if (!dir.exists(base)) next
    if (sub %in% c("impute", "excel_convert", "morpheus")) {
      hits <- list.files(
        base,
        pattern = paste0("pgmatrix_imputed_", dataset),
        full.names = TRUE,
        recursive = FALSE,
        ignore.case = TRUE
      )
      files <- c(files, hits)
    }
  }

  map_root <- file.path(repo_path(config$canonical_inputs$processed_id_mapping$root), dataset)
  if (dir.exists(map_root)) {
    files <- c(files, list.files(map_root, pattern = "\\.(csv|tsv|xlsx|yml|yaml)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
  }

  de_root <- repo_path(config$canonical_inputs$processed_de_enrichment$root)
  de_ds <- file.path(de_root, c("clusterProfiler", "compareGO"), dataset)
  for (d in de_ds) {
    if (dir.exists(d)) {
      files <- c(files, list.files(d, pattern = "\\.(csv|tsv|xlsx|yml|yaml)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
    }
  }

  qc_root <- file.path(repo_path(config$canonical_inputs$qc_reports$root), dataset)
  if (dir.exists(qc_root)) {
    files <- c(files, list.files(qc_root, pattern = "\\.(csv|tsv|md|txt|yml|yaml)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE))
  }

  if (isTRUE(include_derived)) {
    supp_globs <- config$supplementary_table_globs
    supp_globs <- supp_globs[vapply(supp_globs, function(g) grepl(dataset, g, fixed = TRUE) || grepl("\\*", g, fixed = TRUE), logical(1))]
    files <- c(files, glob_paths(supp_globs))
  }

  files <- unique(normalizePath(files[file.exists(files)], winslash = "/", mustWork = FALSE))
  files[!file.info(files)$isdir]
}

classify_export_category <- function(path, config) {
  p <- tolower(gsub("\\\\", "/", path))
  ext <- tolower(tools::file_ext(p))
  if (any(grepl(paste0("/", gsub("\\.", "\\\\.", basename(pg_matrix_input_paths(config)[1])), "$"), p, fixed = FALSE)) ||
      grepl("/pg_matrix/", p)) return("pg_matrix_input")
  if (grepl("/metadata/", p) || grepl("sample_metadata", p)) return("sample_metadata")
  if (grepl("/01_preprocessing/", p)) return("processed_preprocessing")
  if (grepl("/02_id_mapping/", p)) return("protein_id_mapping")
  if (grepl("/04_differential_expression_enrichment/", p)) return("differential_abundance_enrichment")
  if (grepl("/03_qc_exploration/", p) || grepl("/reports/03_qc", p)) return("qc_report")
  if (grepl("/results/tables/", p) || grepl("/results/source_data/", p)) return("derived_results")
  if (grepl("/pride_submission/", p)) return("pride_staging")
  if (grepl("/methods/", p)) return("provenance")
  if (ext %in% c("raw", "wiff", "mzml", "mzxml", "d", "tdf", "baf")) return("external_raw_ms")
  if (grepl("fasta|\\.fa$", p)) return("external_fasta")
  if (grepl("mqpar|search|parameter", p)) return("external_search_parameters")
  "other"
}

export_scope_label <- function(export_level = "pg_matrix_onward") {
  switch(
    export_level,
    pg_matrix_onward = "PRIDE/journal processed-data export from pg_matrix onward",
    .default = paste0("export_level=", export_level)
  )
}

dataset_for_manifest_file <- function(f, datasets) {
  f <- as.character(f)[[1]]
  datasets <- as.character(datasets)
  datasets <- datasets[nzchar(datasets)]
  for (d in datasets) {
    if (grepl(d, f, fixed = TRUE)) return(d)
  }
  "all"
}

manifest_rows_for_files <- function(files, datasets, config, hash_algo = "sha256") {
  if (!length(files)) {
    return(data.frame(
      relative_path = character(),
      dataset = character(),
      export_category = character(),
      size_bytes = numeric(),
      modified_time = character(),
      sha256 = character(),
      intended_for_PRIDE = logical(),
      intended_for_supplement = logical(),
      stringsAsFactors = FALSE
    ))
  }
  info <- file.info(files)
  rel <- vapply(files, relative_to, character(1), root = repo_root())
  ds <- vapply(files, dataset_for_manifest_file, character(1), datasets = datasets)
  cats <- vapply(files, classify_export_category, character(1), config = config)
  hashes <- if (identical(hash_algo, "sha256")) vapply(files, sha256_file, character(1)) else unname(tools::md5sum(files))
  data.frame(
    relative_path = rel,
    dataset = ds,
    export_category = cats,
    size_bytes = as.numeric(info$size),
    modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    sha256 = hashes,
    intended_for_PRIDE = cats %in% c("pg_matrix_input", "sample_metadata", "processed_preprocessing", "protein_id_mapping", "pride_staging", "provenance"),
    intended_for_supplement = cats %in% c("differential_abundance_enrichment", "derived_results", "qc_report"),
    stringsAsFactors = FALSE
  )
}

copy_export_file <- function(src, dest, dry_run = FALSE) {
  dest_dir <- dirname(dest)
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(dry_run)) {
    message("[DRY-RUN] copy ", src, " -> ", dest)
    return(invisible(FALSE))
  }
  src_norm <- normalizePath(src, winslash = "/", mustWork = FALSE)
  dest_norm <- normalizePath(dest, winslash = "/", mustWork = FALSE)
  if (identical(src_norm, dest_norm)) {
    return(invisible(TRUE))
  }
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!ok) warning("Failed to copy: ", src, " -> ", dest, call. = FALSE)
  invisible(ok)
}

supplementary_candidate_files <- function(config, datasets, include_derived = TRUE) {
  globs <- config$supplementary_table_globs
  hits <- glob_paths(globs)
  hits <- hits[is_exportable_result_path(hits)]
  if (!isTRUE(include_derived)) return(hits)
  hits
}

write_validation_summary_md <- function(checks, out_path, export_level, config) {
  n_fail <- sum(checks$status == "FAIL")
  n_warn <- sum(checks$status == "WARN")
  n_info <- sum(checks$status == "INFO")
  lines <- c(
    "# PRIDE export validation summary",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    paste0("Export level: ", export_level),
    paste0("Scope: ", export_scope_label(export_level)),
    "",
    config$scope_note,
    "",
    "## Summary",
    "",
    paste0("- FAIL: ", n_fail),
    paste0("- WARN: ", n_warn),
    paste0("- INFO: ", n_info),
    paste0("- PASS: ", sum(checks$status == "PASS")),
    "",
    "## Checks",
    ""
  )
  for (i in seq_len(nrow(checks))) {
    lines <- c(lines, paste0("- **", checks$check[[i]], "**: ", checks$status[[i]], " — ", checks$detail[[i]]))
  }
  lines <- c(
    lines,
    "",
    "## Reproducibility boundary",
    "",
    "The committed workflow is reproducible from the pg_matrix stage onward.",
    "Missing raw/vendor MS or search-engine files are expected in partial/processed PRIDE depositions.",
    ""
  )
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, out_path)
  invisible(out_path)
}
