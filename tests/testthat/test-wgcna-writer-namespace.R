source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "wgcna_paths.R"))

# Phase 6G.8: the WGCNA migration.
#
# The largest domain in the repository: 24 registered writers, 123 declared
# outputs, 24,015 lines of analysis code, and 5,926 historical files across 21
# roots. It is also the one whose migration needed the legacy-write detector
# generalized first, because results/reviewer_audit/ carries no stage number
# and five writers wrote there while their declarations already said
# results/wgcna/... with the split-brain gate reporting zero.

WRITERS <- local({
  reg <- read_pipeline_registry(repo_path("pipeline.yml"))
  s <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                      include_unsupported = TRUE)
  s <- s[startsWith(s$script, "analysis/wgcna/"), , drop = FALSE]
  sort(unique(sub("[.]R$", "", basename(s$script))))
})
# build_module_spatial_networks.R sits in pipeline.yml's legacy: section, so it
# is not a registered writer and is deliberately not migrated.
LEGACY_ONLY <- "build_module_spatial_networks"
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")

sp <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}
registry_steps <- function() {
  reg <- read_pipeline_registry(repo_path("pipeline.yml"))
  s <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                      include_unsupported = TRUE)
  s[!duplicated(s$script), , drop = FALSE]
}
writer_file <- function(aid) {
  repo_path(file.path("analysis", "wgcna", paste0(aid, ".R")))
}

is_skippable <- function(x) {
  tryCatch(is.null(x) || (is.symbol(x) && !nzchar(as.character(x))),
           error = function(...) TRUE)
}
walk <- function(e, fn) {
  if (!tryCatch({ e; TRUE }, error = function(...) FALSE)) return(invisible(NULL))
  fn(e)
  if (is.call(e) || is.expression(e) || is.list(e)) {
    for (i in seq_along(e)) {
      el <- tryCatch(e[[i]], error = function(...) NULL)
      if (is_skippable(el)) next
      walk(el, fn)
    }
  }
  invisible(NULL)
}
cname <- function(e) {
  if (!is.call(e)) return(NA_character_)
  fn <- e[[1]]
  if (is.name(fn)) return(as.character(fn))
  if (is.call(fn) && length(fn) == 3L && is.name(fn[[1]]) &&
      as.character(fn[[1]]) %in% c("::", ":::")) return(as.character(fn[[3]]))
  NA_character_
}

# Derived from the resolver file, never hard-coded.
#
# A hard-coded list went stale three times in this phase. The last time, the
# eight family helpers added for the previously unserved families were absent
# from it, so both the source-adder and the ordering check reported "all clean"
# while four freshly converted files had no source line at all. Deriving the
# names means a new helper is covered the moment it is written.
RESOLVERS <- local({
  src <- readLines(repo_path("R", "wgcna_paths.R"), warn = FALSE)
  decl <- grep("^wgcna_[A-Za-z0-9_.]* <- function", src, value = TRUE)
  sort(unique(sub("^(wgcna_[A-Za-z0-9_.]*) <- function.*$", "\\1", decl)))
})

# --- sections 17 and 34: writer truth, all 24, no sampling ---------------

testthat::test_that("every registered WGCNA writer declares only normalized destinations", {
  s <- registry_steps()
  testthat::expect_identical(length(WRITERS), 24L)
  for (aid in WRITERS) {
    i <- which(s$script == paste0("analysis/wgcna/", aid, ".R"))
    testthat::expect_length(i, 1L)
    declared <- sp(s$produces[i[1]])
    testthat::expect_gt(length(declared), 0L)
    root <- paste0("results/wgcna/", aid, "/")
    testthat::expect_true(all(startsWith(declared, root)),
      info = paste(aid, "declares outside its namespace:",
                   paste(declared[!startsWith(declared, root)], collapse = ", ")))
    ## none of the historical roots may survive in a declaration
    for (ns in c("06_modules_WGCNA", "reviewer_audit",
                 "04_differential_expression_enrichment")) {
      testthat::expect_false(any(grepl(ns, declared, fixed = TRUE)),
        info = paste(aid, "still declares", ns))
    }
    testthat::expect_false(any(startsWith(declared, "data/processed/")),
      info = paste(aid, "still declares a data/processed destination"))
  }
})

testthat::test_that("every declared child is one the output layout allows", {
  s <- registry_steps()
  allowed <- output_layout_children()
  for (aid in WRITERS) {
    declared <- sp(s$produces[which(s$script == paste0("analysis/wgcna/", aid, ".R"))[1]])
    for (p in declared) {
      seg <- strsplit(p, "/", fixed = TRUE)[[1]]
      testthat::expect_true(length(seg) >= 5L, info = p)
      testthat::expect_true(identical(seg[[2]], "wgcna"), info = p)
      testthat::expect_true(identical(seg[[3]], aid), info = p)
      testthat::expect_true(seg[[5]] %in% allowed,
                            info = paste(p, "has child", seg[[5]]))
    }
  }
})

testthat::test_that("no normalized WGCNA destination collides with another", {
  s <- registry_steps()
  s <- s[startsWith(s$script, "analysis/wgcna/"), , drop = FALSE]
  all_out <- unlist(lapply(seq_len(nrow(s)), function(i) sp(s$produces[i])))
  concrete <- all_out[grepl("[.][A-Za-z0-9]{2,5}$", basename(all_out))]
  testthat::expect_gt(length(concrete), 0L)
  testthat::expect_identical(sum(duplicated(concrete)), 0L,
    info = paste("colliding:", paste(unique(concrete[duplicated(concrete)]),
                                     collapse = ", ")))
})

testthat::test_that("no registered writer calls a legacy directory factory", {
  for (aid in WRITERS) {
    ex <- parse(writer_file(aid))
    hit <- character(0)
    for (e in ex) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% c("create_module_dirs", "module_paths",
                                  "qc_paths", "wgcna_downstream_paths")) {
        hit <<- c(hit, nm)
      }
    })
    testthat::expect_identical(unique(hit), character(0),
      info = paste(aid, "still calls", paste(unique(hit), collapse = ", ")))
  }
})

testthat::test_that("wgcna_downstream_paths is gone, so it cannot come back", {
  ## It wrapped create_module_dirs() and had six call sites. Removing it keeps
  ## the legacy factory from returning through a helper, the same way qc_paths()
  ## was removed in Phase 6G.5.
  lib <- repo_path("R", "wgcna_downstream_utils.R")
  testthat::expect_true(file.exists(lib))
  txt <- paste(readLines(lib, warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("wgcna_downstream_paths <- function", txt, fixed = TRUE))
})

# --- sections 10 and 29: PARSE, LOAD, RESOLVE, REACHABLE -----------------

testthat::test_that("every registered WGCNA writer parses", {
  for (aid in WRITERS) testthat::expect_silent(invisible(parse(writer_file(aid))))
})

testthat::test_that("the resolver source precedes EVERY resolver call, not just the first", {
  ## build_wgcna_modules.R has an early dry-run block with its own nested
  ## bootstrap that runs entirely above the main one. A source appended after
  ## the last top-level source() landed at line 197 while the first call sat at
  ## line 38: present, but too late. Reachability alone would not catch that,
  ## so ordering is asserted here.
  ## Repo-wide, not just analysis/wgcna. The 174 live reads span nine consumer
  ## areas, and a directory-scoped version of this check missed
  ## analysis/qc/summarize_marker_detectability.R calling wgcna_dirs() with no
  ## source at all.
  tracked <- system2("git", c("-C", repo_root(), "ls-files"), stdout = TRUE)
  rel <- tracked[grepl("^(analysis|R|tools|audits)/.*[.][Rr]$", tracked)]
  rel <- rel[basename(rel) != "wgcna_paths.R"]
  files <- repo_path(rel)
  files <- files[file.exists(files)]
  offenders <- character(0)
  for (f in files) {
    ln <- readLines(f, warn = FALSE)
    code <- ln
    code[grepl("^[[:space:]]*#", code)] <- ""
    call_lines <- sort(unique(unlist(lapply(RESOLVERS, function(r)
      grep(paste0(r, "("), code, fixed = TRUE)))))
    if (!length(call_lines)) next
    src_lines <- grep("wgcna_paths.R", code, fixed = TRUE)
    for (L in call_lines) {
      if (!length(src_lines[src_lines < L])) {
        offenders <- c(offenders, paste0(basename(f), ":", L))
      }
    }
  }
  testthat::expect_identical(offenders, character(0),
    info = paste("resolver called before it is sourced:",
                 paste(offenders, collapse = ", ")))
})

testthat::test_that("every file calling the WGCNA resolver also loads it", {
  tracked <- system2("git", c("-C", repo_root(), "ls-files"), stdout = TRUE)
  code <- tracked[grepl("^(analysis|R)/.*[.][Rr]$", tracked)]
  code <- code[basename(code) != "wgcna_paths.R"]
  offenders <- character(0)
  for (rel in code) {
    f <- repo_path(rel)
    if (!file.exists(f)) next
    ln <- readLines(f, warn = FALSE)
    ln[grepl("^[[:space:]]*#", ln)] <- ""
    if (!any(vapply(RESOLVERS, function(r)
      any(grepl(paste0(r, "("), ln, fixed = TRUE)), logical(1)))) next
    if (!any(grepl("wgcna_paths.R", ln, fixed = TRUE))) offenders <- c(offenders, rel)
  }
  testthat::expect_identical(offenders, character(0),
    info = paste("call the resolver without sourcing it:",
                 paste(offenders, collapse = ", ")))
})

# --- sections 30 and 31: dry-run capability ------------------------------

testthat::test_that("the dry-run capability audit covers every entrypoint and is honest", {
  f <- repo_path("audits", "phase6g_wgcna_dry_run_capability.csv")
  testthat::expect_true(file.exists(f))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  reg <- d[as.logical(d$registered), , drop = FALSE]
  testthat::expect_identical(nrow(reg), 24L)
  testthat::expect_true(all(reg$classification %in%
    c("PATH_TESTABLE_NO_COMPUTE", "PATH_VERIFIED_STRUCTURALLY_ONLY")))
  ## the legacy-only helper is recorded but not counted as a registered writer
  testthat::expect_true(any(!as.logical(d$registered)))
  testthat::expect_true(all(nzchar(reg$reason)))

  ## Every structurally-only writer must have a stated reason, because the rule
  ## is that it is NOT executed. An empty reason would mean the classification
  ## was asserted rather than derived.
  so <- reg[reg$classification == "PATH_VERIFIED_STRUCTURALLY_ONLY", ]
  testthat::expect_true(all(nzchar(so$reason)))
})

# --- section 2: the state carriers, including the fourth -----------------

testthat::test_that("the state-carrier audit records all four carriers with roles", {
  f <- repo_path("audits", "phase6g_wgcna_state_carriers.csv")
  testthat::expect_true(file.exists(f))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(d), 4L)
  testthat::expect_identical(sum(d$frozen_state_role == "CANONICAL_FROZEN_STATE"), 3L)
  testthat::expect_identical(sum(d$frozen_state_role == "FAILED_RUN_PROVENANCE"), 1L)
  ## the three canonical hashes the phase brief pins
  h <- setNames(substr(d$sha256, 1, 8), d$carrier)
  testthat::expect_identical(unname(h[["neuron_neuropil"]]), "0fca92a6")
  testthat::expect_identical(unname(h[["neuron_soma"]]), "e218b20c")
  testthat::expect_identical(unname(h[["microglia"]]), "62fba847")
  ## the failed carrier is not exposed through any resolver
  fr <- d[d$frozen_state_role == "FAILED_RUN_PROVENANCE", ]
  testthat::expect_identical(fr$resolver_exposure, "NONE - unreachable by construction")
  testthat::expect_false(as.logical(fr$active_producer))
})

testthat::test_that("the failed microglia state cannot be reached as a dataset scope", {
  ## Section 8: a resolver must never substitute the failed state for the
  ## canonical microglia one. Two independent barriers are asserted.
  testthat::expect_error(validate_dataset("microglia_failed_20260720_133211"),
                         "Unsupported dataset")
  testthat::expect_false("microglia_failed_20260720_133211" %in% valid_datasets())

  ## and the resolver, given the real dataset, names only the canonical path
  p <- wgcna_final_state("microglia")
  testthat::expect_true(grepl("/01_WGCNA/microglia/wgcna_final_model_state[.]rds$", p) ||
                          grepl("/wgcna/build_wgcna_modules/microglia/models/", p))
  testthat::expect_false(grepl("microglia_failed", p, fixed = TRUE))
})

# --- sections 9 and 26: reader precedence, all four presence states ------

testthat::test_that("resolution prefers normalized, falls back to historical", {
  root <- withr::local_tempdir()
  withr::local_envvar(c(PROTEOMICS_PROJECT_ROOT = root))
  FN <- "wgcna_final_model_state.rds"
  norm <- file.path(root, "results", "wgcna", "build_wgcna_modules",
                    "neuron_neuropil", "models")
  hist <- file.path(root, "data", "processed", "06_modules_WGCNA", "01_WGCNA",
                    "neuron_neuropil")
  dir.create(norm, recursive = TRUE); dir.create(hist, recursive = TRUE)
  np <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

  ## neither -> the normalized path, so the caller's own missing-input message
  ## names the future home
  testthat::expect_identical(np(wgcna_final_state("neuron_neuropil")),
                             np(file.path(norm, FN)))
  ## historical only -> historical
  writeLines("h", file.path(hist, FN))
  testthat::expect_identical(np(wgcna_final_state("neuron_neuropil")),
                             np(file.path(hist, FN)))
  ## both -> normalized
  writeLines("n", file.path(norm, FN))
  testthat::expect_identical(np(wgcna_final_state("neuron_neuropil")),
                             np(file.path(norm, FN)))
  ## normalized only -> normalized
  unlink(file.path(hist, FN))
  testthat::expect_identical(np(wgcna_final_state("neuron_neuropil")),
                             np(file.path(norm, FN)))
})

testthat::test_that("an empty normalized directory cannot shadow historical data", {
  ## Availability is decided on the required artifact, never on directory
  ## presence: wgcna_dirs(create = TRUE) and a dry run both leave skeletons,
  ## and one of those skeletons contains a tables/source_data subdirectory, so
  ## counting entries is not enough either.
  root <- withr::local_tempdir()
  withr::local_envvar(c(PROTEOMICS_PROJECT_ROOT = root))
  norm <- file.path(root, "results", "wgcna", "build_wgcna_modules",
                    "neuron_neuropil", "tables", "modules")
  hist <- file.path(root, "results", "tables", "06_modules_WGCNA", "01_WGCNA",
                    "neuron_neuropil", "modules")
  dir.create(norm, recursive = TRUE); dir.create(hist, recursive = TRUE)
  writeLines("x", file.path(hist, "WGCNA_modules_long.csv"))
  np <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

  ## empty normalized directory -> historical
  testthat::expect_identical(np(wgcna_modules_dir("neuron_neuropil")), np(hist))
  ## a subdirectory only is still empty of contract files
  dir.create(file.path(norm, "source_data"))
  testthat::expect_identical(np(wgcna_modules_dir("neuron_neuropil")), np(hist))
  ## one regular file flips it
  writeLines("y", file.path(norm, "WGCNA_modules_long.csv"))
  testthat::expect_identical(np(wgcna_modules_dir("neuron_neuropil")), np(norm))
})

# --- section 8: object identity ------------------------------------------

testthat::test_that("each dataset resolves to its own state and no other", {
  for (ds in DATASETS) {
    p <- wgcna_final_state(ds)
    testthat::expect_true(grepl(paste0("/", ds, "/"), p, fixed = TRUE),
                          info = paste(ds, "resolved to", p))
    for (other in setdiff(DATASETS, ds)) {
      testthat::expect_false(grepl(paste0("/", other, "/"), p, fixed = TRUE),
                             info = paste(ds, "leaked into", other))
    }
  }
})

testthat::test_that("the adjudicated final label lookup is not confused with raw labels", {
  ## Phase 6F adjudicated "final label" as scientifically meaningful: the
  ## adjudicated Stage-07 label as opposed to the raw Stage-01 label. The
  ## filename is therefore preserved verbatim and must not be renamed.
  p <- wgcna_final_label_lookup("neuron_neuropil")
  testthat::expect_true(grepl("WGCNA_final_label_lookup[.]csv$", p))
  testthat::expect_true(grepl("interpretable_summary|summarize_module_interpretation", p))
  ## and it is a different object from the module membership table
  testthat::expect_false(identical(p, wgcna_modules_long("neuron_neuropil")))
})

testthat::test_that("the scope argument is never defaulted away", {
  ## A resolver that silently defaulted the scope to global would return one
  ## dataset's modules for another.
  a <- wgcna_dirs("build_wgcna_modules", "microglia")$tables
  b <- wgcna_dirs("build_wgcna_modules", "neuron_soma")$tables
  testthat::expect_false(identical(a, b))
  testthat::expect_true(grepl("/microglia/", a, fixed = TRUE))
  ## an empty scope becomes global explicitly, not silently some dataset
  g <- wgcna_dirs("build_wgcna_modules", "")$tables
  testthat::expect_true(grepl("/global/", g, fixed = TRUE))
})

# --- section 15 and 32: the resolver's own contract ----------------------

testthat::test_that("wgcna_dirs produces the contract shape and creates nothing by default", {
  d <- wgcna_dirs("build_wgcna_modules", "neuron_neuropil")
  testthat::expect_true(grepl("results/wgcna/build_wgcna_modules/neuron_neuropil/tables$", d$tables))
  testthat::expect_true(grepl("results/wgcna/build_wgcna_modules/neuron_neuropil/models$", d$models))
  ## logs is an alias of manifests, and figures of plots, which is why 266 uses
  ## of $figures/$logs/$source_data needed no change during the migration
  testthat::expect_identical(d$logs, d$manifests)
  testthat::expect_identical(d$figures, d$plots)
  testthat::expect_true(grepl("/tables/source_data$", d$source_data))
})

testthat::test_that("an unknown child or legacy family is refused, not silently accepted", {
  testthat::expect_error(
    wgcna_artifact_candidates("x.csv", "build_wgcna_modules", "01_WGCNA",
                              "neuron_neuropil", child = "processed"))
  testthat::expect_error(
    wgcna_artifact_candidates("x.csv", "build_wgcna_modules", "01_WGCNA",
                              "neuron_neuropil", child = "tables",
                              legacy_family = "not_a_family"))
})

testthat::test_that("candidates are normalized first and the reviewer family is flat", {
  cand <- wgcna_artifact_candidates("wgcna_final_model_state.rds",
                                    "build_wgcna_modules", "01_WGCNA",
                                    "neuron_neuropil", "models", "processed")
  testthat::expect_identical(names(cand)[1], "normalized")
  testthat::expect_true(grepl("results/wgcna/build_wgcna_modules/neuron_neuropil/models/",
                              cand[["normalized"]], fixed = TRUE))
  testthat::expect_true(grepl("data/processed/06_modules_WGCNA/01_WGCNA/neuron_neuropil/",
                              cand[["legacy"]], fixed = TRUE))

  ## the reviewer_audit family has no 06_modules_WGCNA segment
  rv <- wgcna_artifact_candidates("WGCNA_module_adjudication.csv",
                                  "adjudicate_module_labels",
                                  "wgcna_label_adjudication", "global",
                                  "tables", "reviewer_audit",
                                  legacy_scoped = FALSE)
  testthat::expect_true(grepl("results/reviewer_audit/wgcna_label_adjudication/",
                              rv[["legacy"]], fixed = TRUE))
  testthat::expect_false(grepl("06_modules_WGCNA", rv[["legacy"]], fixed = TRUE))
})

testthat::test_that("a zero-length argument cannot collapse a path to nothing", {
  p <- wgcna_dir_any("build_wgcna_modules", "01_WGCNA", "neuron_neuropil",
                     "tables", "tables", TRUE, NULL)
  testthat::expect_length(p, 1L)
  testthat::expect_true(nzchar(p))
})

# --- section 15: the historical fallback must remain addressable ---------

testthat::test_that("the resolver still names every historical family", {
  txt <- paste(readLines(repo_path("R", "wgcna_paths.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl("06_modules_WGCNA", txt, fixed = TRUE))
  testthat::expect_true(grepl("reviewer_audit", txt, fixed = TRUE))
  for (fam in c("tables", "figures", "source_data", "logs", "reports",
                "processed", "reviewer_audit")) {
    testthat::expect_true(fam %in% WGCNA_LEGACY_FAMILIES, info = fam)
  }
})

testthat::test_that("the empty-directory predicate exists exactly once", {
  ## Three inline copies of this guard is how spatial_systems_dir_any stayed
  ## broken for two phases after its sibling was fixed.
  src <- readLines(repo_path("R", "wgcna_paths.R"), warn = FALSE)
  defs <- grep("has_files <- function", src, fixed = TRUE)
  testthat::expect_identical(length(defs), 1L)
})

# --- section 12: publication_source_data stays unmigrated ---------------

testthat::test_that("publication_source_data keeps its own write sites", {
  s <- registry_steps()
  psd <- s[startsWith(s$script, "analysis/publication_source_data/"), , drop = FALSE]
  testthat::expect_gt(nrow(psd), 0L)
  for (i in seq_len(nrow(psd))) {
    prod <- sp(psd$produces[i])
    ## no publication_source_data output may have been moved into results/wgcna
    testthat::expect_false(any(startsWith(prod, "results/wgcna/")),
      info = paste(psd$script[i], "gained a WGCNA output contract"))
  }
  ## and its two Stage13 path-level artifacts are untouched
  i <- which(psd$script == "analysis/publication_source_data/build_biological_claims_table.R")
  testthat::expect_length(i, 1L)
  declared <- sp(psd$produces[i[1]])
  testthat::expect_true(any(grepl("wgcna_stage13_claim_cardinality_audit.csv",
                                  declared, fixed = TRUE)))
  testthat::expect_true(any(grepl("microglia_wgcna_overlap_stage13_identity_audit.csv",
                                  declared, fixed = TRUE)))
})

testthat::test_that("WGCNA owns no stage13 path-level artifact", {
  s <- registry_steps()
  s <- s[startsWith(s$script, "analysis/wgcna/"), , drop = FALSE]
  for (i in seq_len(nrow(s))) {
    testthat::expect_false(any(grepl("stage13", sp(s$produces[i]), fixed = TRUE)),
      info = paste(s$script[i], "declares a stage13 output"))
  }
})

# --- section 13: contracts agree with code ------------------------------

testthat::test_that("results ownership records the normalized WGCNA families", {
  d <- utils::read.csv(repo_path("config", "results_ownership.csv"),
                       stringsAsFactors = FALSE)
  missing <- character(0)
  for (aid in WRITERS) {
    fam <- d[grepl(paste0("results/wgcna/", aid), d$result_family, fixed = TRUE), , drop = FALSE]
    if (!nrow(fam)) { missing <- c(missing, aid); next }
    testthat::expect_true(all(fam$canonical_owner == paste0("analysis/wgcna/", aid, ".R")),
      info = paste(aid, "ownership names another writer"))
  }
  testthat::expect_identical(missing, character(0),
    info = paste("no ownership row for:", paste(missing, collapse = ", ")))
})

testthat::test_that("the preprocessing domain is declared alongside wgcna in the layout", {
  layout <- yaml::read_yaml(repo_path("config", "output_layout.yml"))
  testthat::expect_true("wgcna" %in% unlist(layout$domains))
})
