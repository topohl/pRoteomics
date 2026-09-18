# Canonical and historical locations for QC results.
#
# Kept in its own file, separate from qc_exploration_utils.R, because forty
# analyses outside analysis/qc read QC output and they should not have to
# load the whole QC utility library - protein mapping, dataset inputs and
# the rest - just to resolve a path.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

# Where QC results live, before and after Phase 6G.5.
#
# Historical layout:
#   results/<kind>/03_qc_exploration/<substep>/<dataset>/
# Normalized layout, keyed on analysis identity:
#   results/qc/<analysis_id>/<scope>/<child>/
#
# The historical substep names were almost pure chronology: 00_, 00b_, 00c_,
# 01_ through 07_, with 05_ and 06_ each used by two unrelated analyses. None
# of that survives, because the analysis_id replaces it, and nothing is lost:
# 05_pca_confounding_qc and 05_empirical_roi_marker_discovery are
# assess_pca_confounding and discover_empirical_roi_markers.
#
# qc_paths() was the single factory that eleven of the fifteen writers used,
# so this is the chokepoint: the writers change only in what they pass.
qc_dirs <- function(analysis_id, scope = "global", create = FALSE) {
  if (is.null(scope) || !length(scope) || !nzchar(scope)) scope <- "global"
  canonical_module_dirs("qc", analysis_id, scope = scope, create = create)
}

# Normalized-first resolution for a QC artifact that something else reads.
#
# QC is upstream of every other domain: forty analyses outside analysis/qc read
# its outputs, and they build those paths themselves rather than through one
# helper, so there is no read-side chokepoint to convert. This gives them one.
#
# `owner` is the QC analysis_id that declares the file and `legacy_substep` is
# the directory it used to live in. Both are checked against the registry by
# tests/testthat/test-qc-writer-namespace.R, so a wrong owner is a test
# failure rather than a silent mis-resolution.
#
# Existence is tested on the file, never on its directory, so the empty
# skeleton a dry run leaves behind cannot shadow real data.
qc_artifact_candidates <- function(filename, owner, legacy_substep,
                                   scope = "global",
                                   kind = c("tables", "figures", "reports",
                                            "logs", "source_data")) {
  kind <- match.arg(kind)
  child <- switch(kind, tables = "tables", figures = "plots",
                  reports = "reports", logs = "manifests",
                  source_data = "tables")
  norm <- if (identical(kind, "source_data")) {
    canonical_result_path("qc", owner, scope, child, "source_data", filename)
  } else {
    canonical_result_path("qc", owner, scope, child, filename)
  }
  ## The historical layout is not uniform: qc_paths() appended the dataset, so
  ## most substeps have a scope segment, but the two that took no dataset wrote
  ## straight into the substep directory. Both shapes are offered rather than
  ## assumed, which is what the empirical ROI marker sets need - eight readers
  ## across four domains address them with no scope segment at all.
  legacy_scoped <- path_results(kind, "03_qc_exploration", legacy_substep,
                                scope, filename)
  legacy_flat <- path_results(kind, "03_qc_exploration", legacy_substep, filename)
  c(normalized = norm, legacy_scoped = legacy_scoped, legacy_flat = legacy_flat)
}

# A directory, normalized first, for a reader that wants the whole folder
# rather than one named file.
#
# The normalized directory is preferred only when it is non-empty. Anchoring a
# folder on one probe filename does not work: a dry run creates the normalized
# skeleton, and if the probe happens to be absent every path built from the
# returned root points into an empty directory. That is exactly how the joint
# compartment renderer first reported all seven of its inputs missing.
qc_dir_any <- function(owner, legacy_substep, scope = "global",
                       kind = c("tables", "figures", "reports", "logs",
                                "source_data")) {
  kind <- match.arg(kind)
  child <- switch(kind, tables = "tables", figures = "plots",
                  reports = "reports", logs = "manifests",
                  source_data = "tables")
  norm <- if (identical(kind, "source_data")) {
    canonical_result_path("qc", owner, scope, child, "source_data")
  } else {
    canonical_result_path("qc", owner, scope, child)
  }
  ## Regular files only. list.files() also reports subdirectories, and
  ## qc_dirs(create = TRUE) leaves a tables/source_data subdirectory behind, so
  ## counting entries made an empty normalized directory look populated and
  ## shadow the real historical data.
  has_files <- function(d) {
    if (!dir.exists(d)) return(FALSE)
    e <- list.files(d, full.names = TRUE, no.. = TRUE)
    if (!length(e)) return(FALSE)
    any(!file.info(e)$isdir)
  }
  if (has_files(norm)) return(norm)
  scoped <- path_results(kind, "03_qc_exploration", legacy_substep, scope)
  if (has_files(scoped)) return(scoped)
  flat <- path_results(kind, "03_qc_exploration", legacy_substep)
  if (has_files(flat)) return(flat)
  norm
}

qc_find <- function(filename, owner, legacy_substep, scope = "global",
                    kind = c("tables", "figures", "reports", "logs",
                             "source_data")) {
  kind <- match.arg(kind)
  cand <- qc_artifact_candidates(filename, owner, legacy_substep, scope, kind)
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(unname(hit[1]))
  unname(cand[["normalized"]])
}
