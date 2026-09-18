# Where spatial_validation's results live, before and after Phase 6G.3.
#
# The historical layout grouped these outputs by result family under one stage
# directory: results/tables/11_spatial_systems/<family>/<file>. Several
# analyses wrote into one family, which is what made those families
# "coordinating" in config/results_ownership.csv.
#
# The normalized layout is keyed on analysis identity instead:
# results/spatial_validation/<analysis_id>/<scope>/<child>/<file>. The family
# segment is not reproduced, because it is now redundant: the owner is the
# analysis_id, and the coordination relationship it used to express is
# recorded in config/results_ownership.csv, which is the registry of record
# for ownership. Nothing is lost, it moves from the path to the registry.
#
# One consequence worth stating plainly. Because the normalized path contains
# the producing analysis, two different analyses cannot produce the same
# normalized destination. The collision the ownership registry exists to
# prevent is therefore unreachable by construction after migration, and the
# only place a genuine duplicate writer could exist is the historical layout.
# It does not exist there either: all 71 declared files have distinct names
# and no filename is written by two scripts.
#
# --- why reads need a resolver ---------------------------------------------
#
# These analyses read each other. summarize_spatial_atlas reads five families;
# build_module_spatial_atlas reads three. Those reads used to be plain
# path_results("tables", "11_spatial_systems", family, file) calls, which name
# the historical tree directly.
#
# After the writers migrate, the historical tree still holds the last run's
# output and the normalized tree holds nothing until something is rerun. A
# reader pinned to the historical path would therefore keep reading the old
# object forever, and the migration would look complete while changing
# nothing. A reader pinned to the normalized path would break immediately,
# because nothing has been written there yet and this phase does not rerun
# any analysis.
#
# So reads resolve normalized-first with a legacy fallback: today the legacy
# copy answers, and the first normalized run silently takes over. That
# ordering is the whole point and is asserted in
# tests/testthat/test-spatial-validation-writer-namespace.R across all four
# existence combinations.

# The canonical output directories for one spatial_validation analysis.
# scope is global for every analysis in this domain: none of them takes
# --dataset, they use integration_cli(default_dataset = "all") and loop
# internally, and none of the 75 declared outputs carries a dataset segment.
# create = FALSE by default and the caller creates what it writes. Section 7
# asks for only the lifecycle directories an analysis actually uses, and
# canonical_module_dirs() would otherwise create all six for every writer,
# including plots/ and models/ for an analysis that only writes tables. The
# existing accessors already create lazily, so this also keeps their behaviour
# unchanged.
spatial_systems_dirs <- function(analysis_id, create = FALSE) {
  canonical_module_dirs("spatial_validation", analysis_id, scope = "global",
                        create = create)
}

# Candidates for one artifact, normalized first.
#
# `owner` is the analysis_id that declares the file in pipeline.yml, and
# `legacy_family` is the directory it used to live in under
# 11_spatial_systems. Both are verified against the registry by
# test-spatial-validation-writer-namespace.R, so a wrong owner here is a test
# failure rather than a silent mis-resolution.
#
# `kind` is the historical top-level results directory: "tables" or "figures".
# figures maps onto the normalized child plots/, because the layout contract
# reserves "figures" for the manuscript repository's assembled journal figures.
spatial_systems_artifact_candidates <- function(filename, owner,
                                                legacy_family = NULL,
                                                kind = c("tables", "figures")) {
  kind <- match.arg(kind)
  child <- if (identical(kind, "figures")) "plots" else "tables"
  legacy <- if (is.null(legacy_family) || !length(legacy_family) ||
                !nzchar(legacy_family)) {
    path_results(kind, "11_spatial_systems", filename)
  } else {
    path_results(kind, "11_spatial_systems", legacy_family, filename)
  }
  c(
    normalized = canonical_result_path("spatial_validation", owner, "global",
                                       child, filename),
    legacy = legacy
  )
}

# Resolve one artifact, preferring the normalized copy.
#
# When neither exists the normalized path is returned, so a failure message
# names where the object is supposed to be rather than where it used to be.
resolve_spatial_systems_artifact <- function(filename, owner,
                                             legacy_family = NULL,
                                             kind = c("tables", "figures"),
                                             env = "PROTEOMICS_SPATIAL_SYSTEMS_ROOT") {
  kind <- match.arg(kind)
  override <- Sys.getenv(env, unset = "")
  if (nzchar(override)) {
    return(normalizePath(file.path(override, filename), winslash = "/",
                         mustWork = FALSE))
  }
  cand <- spatial_systems_artifact_candidates(filename, owner, legacy_family, kind)
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(unname(hit[1]))
  unname(cand[["normalized"]])
}

# Resolve an artifact by filename alone, without naming its owner.
#
# The historical families mix owners: the networks family holds files from
# build_animal_spatial_networks, test_network_group_organization and
# validate_network_workbook. A reader that pulls thirteen files out of one
# family would otherwise have to name a different owner for each, and would go
# stale the moment ownership moved.
#
# Keying on the filename is unambiguous here because all 71 files declared
# under 11_spatial_systems have distinct names and no filename is written by
# two scripts. That invariant is what makes the lookup sound, so the function
# enforces it rather than assuming it: two normalized producers of one
# filename is the collision the ownership registry exists to prevent, and it
# stops here instead of silently picking one.
spatial_systems_find <- function(filename, legacy_family = NULL,
                                 kind = c("tables", "figures", "models", "reports")) {
  kind <- match.arg(kind)
  child <- switch(kind, figures = "plots", tables = "tables",
                  models = "models", reports = "reports")
  norm <- Sys.glob(repo_path("results", "spatial_validation", "*", "global",
                             child, filename))
  norm <- norm[file.exists(norm)]
  if (length(norm) > 1L) {
    stop("two spatial_validation analyses produced '", filename, "': ",
         paste(basename(dirname(dirname(dirname(norm)))), collapse = ", "),
         ". One canonical writer per artifact; see config/results_ownership.csv.",
         call. = FALSE)
  }
  if (length(norm) == 1L) return(unname(norm))
  ## Nothing normalized: fall back to the historical family.
  ##
  ## When neither copy exists this returns the historical path rather than the
  ## canonical destination, which is the opposite of what
  ## resolve_spatial_systems_artifact() does. The difference is deliberate:
  ## that function is given an owner and can name the owner's directory, this
  ## one is keyed on a filename alone and cannot. Every call site is a read
  ## guarded by file.exists(), so the path returned in the nothing-exists case
  ## only shapes a diagnostic. The case where authority actually matters, both
  ## copies present, resolves to the normalized copy.
  hist_kind <- if (identical(kind, "figures")) "figures" else "tables"
  if (is.null(legacy_family) || !length(legacy_family) || !nzchar(legacy_family)) {
    path_results(hist_kind, "11_spatial_systems", filename)
  } else {
    path_results(hist_kind, "11_spatial_systems", legacy_family, filename)
  }
}

# A directory, normalized first.
#
# For a reader that lists a directory rather than naming one file, such as a
# test that checks every figure has its source data. The normalized directory
# is preferred only when it is non-empty, so an empty skeleton left by a dry
# run does not hide a populated historical directory.
spatial_systems_dir_any <- function(owner, legacy_family = NULL,
                                    kind = c("tables", "figures", "models", "reports")) {
  kind <- match.arg(kind)
  child <- switch(kind, figures = "plots", tables = "tables",
                  models = "models", reports = "reports")
  norm <- canonical_result_path("spatial_validation", owner, "global", child)
  ## Regular files only: list.files() also reports subdirectories, so a
  ## normalized directory holding nothing but a source_data subdirectory
  ## would otherwise look populated and shadow real historical data. Found
  ## in Phase 6G.5 on the equivalent QC resolver.
  norm_files <- if (dir.exists(norm)) list.files(norm, full.names = TRUE, no.. = TRUE) else character(0)
  if (length(norm_files) && any(!file.info(norm_files)$isdir)) return(norm)
  hist_kind <- if (identical(kind, "figures")) "figures" else "tables"
  if (is.null(legacy_family) || !length(legacy_family) || !nzchar(legacy_family)) {
    path_results(hist_kind, "11_spatial_systems")
  } else {
    path_results(hist_kind, "11_spatial_systems", legacy_family)
  }
}

# Convenience for a reader that wants the directory rather than one file:
# returns the normalized directory when it holds the named file, else the
# historical family directory. A directory is only preferred when it actually
# contains the artifact, so an empty skeleton created by a dry run cannot
# shadow real data.
resolve_spatial_systems_dir <- function(required_file, owner,
                                        legacy_family = NULL,
                                        kind = c("tables", "figures")) {
  kind <- match.arg(kind)
  p <- resolve_spatial_systems_artifact(required_file, owner, legacy_family, kind)
  dirname(p)
}
