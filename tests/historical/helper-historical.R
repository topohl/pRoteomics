# HISTORICAL SNAPSHOTS - not part of current acceptance.
#
# The files in tests/historical/ pin literal facts about the network built
# before 2026-10-19 (original centring, old output layout, microglia
# supermodule cut height 0.40) and about one-off audits made on it. They are
# kept as evidence, not as acceptance criteria: tests/testthat.R runs
# tests/testthat only, and current consistency (model files against the
# contracts that recorded them) is tested there.
#
# They read ONLY the historical root: PROTEOMICS_HISTORICAL_ROOT, or else the
# preserved all_observed_median build declared in config/normalisation_versions.yml
# (the hash-verified archive Analysis/_archive_proteomics_pre_renormalisation_2026-10-19).
# Without an available root every test skips with the reason. Run with
#
#   Rscript -e 'testthat::test_dir("tests/historical")'
#   PROTEOMICS_HISTORICAL_ROOT=<root> Rscript -e 'testthat::test_dir("tests/historical")'

if (!exists("repo_path", mode = "function")) {
  source(testthat::test_path("..", "..", "R", "paths.R"))
}

historical_root <- function() {
  root <- Sys.getenv("PROTEOMICS_HISTORICAL_ROOT", unset = "")
  if (nzchar(root)) return(sub("/+$", "", gsub("\\\\", "/", root)))
  versions <- tryCatch(
    yaml::read_yaml(repo_path("config", "normalisation_versions.yml"))$versions,
    error = function(e) list()
  )
  hit <- Filter(function(v) identical(v$sample_centring, "all_observed_median") &&
                  identical(v$status, "preserved"), versions)
  if (length(hit) == 1L && length(hit[[1]]$root)) sub("/+$", "", hit[[1]]$root) else ""
}

# Use inside test_that(): the root, or a skip with the reason.
local_historical_root <- function() {
  root <- historical_root()
  testthat::skip_if(
    !nzchar(root) || !dir.exists(root),
    paste0("Historical archive not available (",
           if (nzchar(root)) root else "no preserved all_observed_median build declared",
           "); set PROTEOMICS_HISTORICAL_ROOT to run the snapshot checks.")
  )
  root
}

historical_path <- function(...) file.path(local_historical_root(), ...)

skip_unless_historical <- function(path) {
  testthat::skip_if_not(file.exists(path) || dir.exists(path),
                        paste("Historical artefact not present:", path))
}
