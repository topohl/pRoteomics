# Canonical and historical locations for preprocessing outputs.
#
# Path semantics only. No preprocessing or scientific logic belongs in this
# file: it answers "where does this object live" and nothing else, so the
# fifty-odd analyses that read preprocessing output can resolve a path without
# loading the mapping, dataset-input or joint-QC libraries.
#
# Preprocessing is the root of the dependency graph - seventy runtime consumer
# files across nine analysis domains - and its historical layout is the least
# uniform in the repository. It spans two stage namespaces (01_preprocessing
# and 02_id_mapping) and two lifecycles (data/processed for the derived
# scientific objects, results/<kind> for their summaries and manifests), which
# is why the candidate helpers below take the historical family explicitly
# rather than deriving it.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

PREPROCESSING_DOMAIN <- "preprocessing"

# Historical families a preprocessing object could have lived in. "processed"
# is data/processed/<stage>/...; the rest are results/<kind>/<stage>/...
PREPROCESSING_LEGACY_FAMILIES <- c("processed", "tables", "logs", "reports", "figures")

# Normalized children, as config/output_layout.yml declares them.
PREPROCESSING_CHILDREN <- c("tables", "plots", "models", "manifests", "reports")

# Drop zero-length arguments before file.path() sees them. file.path() returns
# character(0) if any argument has length zero, which silently turns a real
# destination into an empty path.
pp_join <- function(...) {
  parts <- Filter(function(x) length(x) && nzchar(x), list(...))
  if (!length(parts)) return(character(0))
  do.call(file.path, parts)
}

# Does this directory actually hold a contract file?
#
# Regular files only, and defined exactly once. dir.exists() is not evidence
# that data exist: a dry run, or dirs(create = TRUE), leaves a skeleton behind,
# and one of those skeletons contains a tables/source_data subdirectory - so
# counting directory entries made an empty normalized directory look populated
# and shadow the real historical data. That bug survived in a second copy of
# this predicate for two phases, which is why there is only one copy here.
pp_has_files <- function(d) {
  if (!length(d) || is.na(d) || !nzchar(d) || !dir.exists(d)) return(FALSE)
  e <- list.files(d, full.names = TRUE, no.. = TRUE)
  if (!length(e)) return(FALSE)
  any(!file.info(e)$isdir)
}

# Where preprocessing results live after Phase 6G.7:
#   results/preprocessing/<analysis_id>/<scope>/<child>/
preprocessing_dirs <- function(analysis_id, scope = "global", suffix = NULL,
                               create = FALSE) {
  if (is.null(scope) || !length(scope) || !nzchar(scope)) scope <- "global"
  canonical_module_dirs(PREPROCESSING_DOMAIN, analysis_id, scope = scope,
                        suffix = suffix, create = create)
}

# Normalized-first candidates for one named preprocessing artifact.
#
# `owner` is the preprocessing analysis_id that declares the file;
# `legacy_stage` is 01_preprocessing or 02_id_mapping and `legacy_substep` the
# directory inside it. Both are checked against the registry by
# tests/testthat/test-preprocessing-writer-namespace.R, so a wrong owner is a
# test failure rather than a silent mis-resolution.
#
# Trailing `...` segments are preserved in both candidates, because
# config/output_layout.yml keeps meaningful post-scope segments: a contrast
# direction (forward/reverse) or a mapping branch is part of the object's
# identity, not chronology.
#
# Existence is tested on the file, never on its directory, so the empty
# skeleton a dry run leaves behind cannot shadow real historical data.
preprocessing_artifact_candidates <- function(filename, owner, legacy_stage,
                                              legacy_substep, scope = "global",
                                              child = "tables",
                                              legacy_family = "processed",
                                              legacy_scoped = TRUE, ...) {
  child <- match.arg(child, PREPROCESSING_CHILDREN)
  legacy_family <- match.arg(legacy_family, PREPROCESSING_LEGACY_FAMILIES)
  if (is.null(scope) || !length(scope) || !nzchar(scope)) scope <- "global"

  norm <- canonical_result_path(PREPROCESSING_DOMAIN, owner, scope, child,
                                ...) |> pp_join(filename)

  ## The historical layout put the scope in different places, and two substeps
  ## took no dataset segment at all, so both shapes are offered rather than
  ## assumed.
  legacy_root <- if (identical(legacy_family, "processed")) {
    pp_join(path_processed(legacy_stage), legacy_substep)
  } else {
    pp_join(path_results(legacy_family, legacy_stage), legacy_substep)
  }
  legacy_a <- if (isTRUE(legacy_scoped)) {
    pp_join(legacy_root, scope, ..., filename)
  } else {
    pp_join(legacy_root, ..., filename)
  }
  legacy_b <- if (isTRUE(legacy_scoped)) {
    pp_join(legacy_root, ..., filename)
  } else {
    pp_join(legacy_root, scope, ..., filename)
  }

  c(normalized = unname(norm), legacy = unname(legacy_a),
    legacy_alt = unname(legacy_b))
}

# Resolve one named preprocessing artifact, normalized first.
#
#   normalized present + historical present -> normalized
#   normalized absent  + historical present -> historical
#   normalized present + historical absent  -> normalized
#   neither                                 -> the normalized path, so the
#                                              caller's own missing-input
#                                              message names the future home
preprocessing_find <- function(filename, owner, legacy_stage, legacy_substep,
                               scope = "global", child = "tables",
                               legacy_family = "processed",
                               legacy_scoped = TRUE, ...) {
  cand <- preprocessing_artifact_candidates(filename, owner, legacy_stage,
                                            legacy_substep, scope, child,
                                            legacy_family, legacy_scoped, ...)
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(unname(hit[1]))
  unname(cand[["normalized"]])
}

# A directory, normalized first, for a reader that wants the whole folder
# rather than one named file.
#
# The normalized directory is preferred only when it actually holds regular
# files. Counting directory entries is not enough: preprocessing_dirs(create =
# TRUE) leaves a tables/source_data subdirectory behind, so an entry count made
# an empty normalized directory look populated and shadow the real historical
# data. Anchoring on a single probe filename is not enough either - that is how
# the joint compartment renderer once reported all seven of its inputs missing.
preprocessing_dir_any <- function(owner, legacy_stage, legacy_substep,
                                  scope = "global", child = "tables",
                                  legacy_family = "processed",
                                  legacy_scoped = TRUE, ...) {
  child <- match.arg(child, PREPROCESSING_CHILDREN)
  legacy_family <- match.arg(legacy_family, PREPROCESSING_LEGACY_FAMILIES)
  if (is.null(scope) || !length(scope) || !nzchar(scope)) scope <- "global"


  norm <- canonical_result_path(PREPROCESSING_DOMAIN, owner, scope, child, ...)
  if (pp_has_files(norm)) return(unname(norm))

  legacy_root <- if (identical(legacy_family, "processed")) {
    pp_join(path_processed(legacy_stage), legacy_substep)
  } else {
    pp_join(path_results(legacy_family, legacy_stage), legacy_substep)
  }
  ordered <- if (isTRUE(legacy_scoped)) {
    list(pp_join(legacy_root, scope, ...), pp_join(legacy_root, ...))
  } else {
    list(pp_join(legacy_root, ...), pp_join(legacy_root, scope, ...))
  }
  for (d in ordered) if (pp_has_files(d)) return(unname(d))
  unname(norm)
}

# ------------------------------------------------- named objects readers want
#
# Three preprocessing objects carry almost the whole downstream dependency
# graph between them, and their readers build the path themselves rather than
# going through one helper. These give them a chokepoint, so a reader is one
# line rather than six resolver arguments repeated across nine domains.

# The merged sample-metadata workbook: the canonical AnimalID / StressGroup /
# spatial-unit assignment. Nineteen analyses read it.
preprocessing_module_score_metadata <- function(dataset) {
  preprocessing_find("sample_metadata_merged_clean_for_module_scores.xlsx",
                     owner = "build_module_score_metadata",
                     legacy_stage = "01_preprocessing",
                     legacy_substep = "06_merged_metadata_module_score",
                     scope = dataset, child = "tables",
                     legacy_family = "processed")
}

# The serialised joint-compartment QC bundle, a versioned contract object read
# by four QC analyses and one WGCNA audit.
preprocessing_joint_qc_bundle <- function() {
  preprocessing_find("joint_compartment_qc_matrices.rds",
                     owner = "build_joint_protigy_input",
                     legacy_stage = "01_preprocessing",
                     legacy_substep = "joint_compartment_qc",
                     scope = "global", child = "models",
                     legacy_family = "processed")
}

# The audit-table directory that sits beside that bundle.
preprocessing_joint_qc_tables <- function() {
  preprocessing_dir_any("build_joint_protigy_input", "01_preprocessing",
                        "joint_compartment_qc", "global", "tables", "processed")
}

# The mapped contrast directory: docs/file_contracts.tsv object
# mapped_contrast_csv, read by six analysis domains.
preprocessing_mapped_contrast_dir <- function(dataset, direction = "forward") {
  preprocessing_mapping_dir_find("mapped", dataset, direction)
}

# ---------------------------------------------------------------- output roots
#
# Two preprocessing writers take a configurable output root from the
# environment, so their destination cannot be a fixed path. Both keep their
# override semantics exactly: an explicit root is a deliberate branch replay
# and still gets the historical root-relative layout. Only the default moves
# into the normalized namespace, and only these two functions know that, so no
# literal normalized path is scattered across the writers (Phase 6G.7 section
# 23).

pp_paths_equal <- function(left, right) {
  left <- normalizePath(left, winslash = "/", mustWork = FALSE)
  right <- normalizePath(right, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows") {
    left <- tolower(left)
    right <- tolower(right)
  }
  identical(left, right)
}

# TRUE when a resolved root is still the historical default, i.e. no deliberate
# override is in force and the normalized namespace should be used instead.
preprocessing_root_is_default <- function(root, default_root) {
  if (!length(root) || !nzchar(root)) return(TRUE)
  pp_paths_equal(root, default_root)
}

# Destination for one split ProTigy contrast direction.
#
# Historical: data/processed/01_preprocessing/gct_extractR/<dataset>/<direction>
# Normalized: results/preprocessing/extract_protigy_contrasts/<dataset>/tables/<direction>
#
# `root` is the resolved PROTEOMICS_GCT_OUTPUT_ROOT. When it is the historical
# default the normalized destination is returned; when it is an explicit
# override - including the legacy comparison replay branch - the caller's root
# is honoured unchanged.
preprocessing_gct_extract_dir <- function(dataset, direction = NULL,
                                          root = NULL,
                                          default_root = path_processed("01_preprocessing", "gct_extractR")) {
  if (!preprocessing_root_is_default(root, default_root)) {
    return(unname(pp_join(root, dataset, direction)))
  }
  unname(canonical_result_path(PREPROCESSING_DOMAIN, "extract_protigy_contrasts",
                               dataset, "tables") |> pp_join(direction))
}

# Destination for the extraction index, which is a manifest of the split rather
# than a contrast table.
preprocessing_gct_extract_manifest_dir <- function(dataset, root = NULL,
                                                   default_root = path_processed("01_preprocessing", "gct_extractR")) {
  if (!preprocessing_root_is_default(root, default_root)) {
    return(unname(pp_join(root, dataset)))
  }
  unname(canonical_result_path(PREPROCESSING_DOMAIN, "extract_protigy_contrasts",
                               dataset, "manifests"))
}

# The extraction provenance manifest, as a write destination.
#
# Two twin functions computed this path independently before Phase 6G.7 -
# gct_extract_contract_manifest_path() on the writer side and
# canonical_gct_extract_manifest() on the reader side - so both now delegate
# here and cannot drift apart.
preprocessing_gct_extract_manifest_path <- function(dataset, root = NULL,
                                                    default_root = path_processed("01_preprocessing", "gct_extractR")) {
  unname(pp_join(preprocessing_gct_extract_manifest_dir(dataset, root, default_root),
                 "canonical_gct_extract_manifest.csv"))
}

# The same manifest, for a reader: normalized first, historical fallback.
#
# map_protein_identifiers() hard-gates on this file, so a reader that only
# looked at the normalized path would refuse to run against the historical
# extraction that exists on disk today.
preprocessing_gct_extract_manifest_find <- function(dataset, root = NULL,
                                                    default_root = path_processed("01_preprocessing", "gct_extractR")) {
  norm <- preprocessing_gct_extract_manifest_path(dataset, root, default_root)
  if (file.exists(norm)) return(norm)
  base <- if (length(root) && nzchar(root)) root else default_root
  hist <- unname(pp_join(base, dataset, "canonical_gct_extract_manifest.csv"))
  if (file.exists(hist)) return(hist)
  norm
}

# One split-contrast direction, for a reader: normalized first, historical
# fallback. Existence is decided on regular files, not on the directory.
preprocessing_gct_extract_dir_find <- function(dataset, direction, root = NULL,
                                               default_root = path_processed("01_preprocessing", "gct_extractR")) {
  norm <- preprocessing_gct_extract_dir(dataset, direction, root, default_root)
  if (pp_has_files(norm)) return(norm)
  base <- if (length(root) && nzchar(root)) root else default_root
  hist <- unname(pp_join(base, dataset, direction))
  if (pp_has_files(hist)) return(hist)
  norm
}

# Destination for one identifier-mapping branch.
#
# Historical: data/processed/02_id_mapping/<branch>/<dataset>/<direction>/per_file
# Normalized: results/preprocessing/map_protein_identifiers/<dataset>/tables/<branch>/<direction>/per_file
#
# The dataset moves from the third segment to the scope position, so this is
# not a prefix swap and the writer cannot build it by concatenating a root.
preprocessing_mapping_dir <- function(branch, dataset, direction,
                                      leaf = "per_file", root = NULL,
                                      default_root = path_processed("02_id_mapping")) {
  if (!preprocessing_root_is_default(root, default_root)) {
    return(unname(pp_join(root, branch, dataset, direction, leaf)))
  }
  unname(canonical_result_path(PREPROCESSING_DOMAIN, "map_protein_identifiers",
                               dataset, "tables") |>
           pp_join(branch, direction, leaf))
}

# The summary/manifest/report side of identifier mapping.
#
# Historical: results/<kind>/<namespace>/MapThatProt_batch/<dataset>/...
# Normalized: results/preprocessing/map_protein_identifiers/<dataset>/<child>/...
#
# Returned already scoped to the dataset, so callers append only the segments
# that carry meaning - branch, direction, and the leaf kind.
preprocessing_mapping_result_dir <- function(kind, dataset, root = NULL,
                                             default_root = path_processed("02_id_mapping"),
                                             namespace = "02_id_mapping") {
  kind <- match.arg(kind, c("tables", "logs", "reports", "figures"))
  if (!preprocessing_root_is_default(root, default_root)) {
    return(unname(pp_join(path_results(kind, namespace, "MapThatProt_batch"), dataset)))
  }
  child <- switch(kind, tables = "tables", logs = "manifests",
                  reports = "reports", figures = "plots")
  unname(canonical_result_path(PREPROCESSING_DOMAIN, "map_protein_identifiers",
                               dataset, child))
}

# One identifier-mapping branch, for a reader: normalized first, historical
# fallback.
#
# This is the read-side chokepoint for the most heavily depended-on object in
# the repository - the mapped contrast CSVs are a versioned file contract
# (docs/file_contracts.tsv, mapped_contrast_csv) with consumers in six analysis
# domains - so the readers convert to this rather than each growing its own
# two-candidate fallback.
preprocessing_mapping_dir_find <- function(branch, dataset, direction,
                                           leaf = "per_file", root = NULL,
                                           default_root = path_processed("02_id_mapping")) {
  norm <- preprocessing_mapping_dir(branch, dataset, direction, leaf, root, default_root)
  if (pp_has_files(norm)) return(norm)
  base <- if (length(root) && nzchar(root)) root else default_root
  hist <- unname(pp_join(base, branch, dataset, direction, leaf))
  if (pp_has_files(hist)) return(hist)
  norm
}
