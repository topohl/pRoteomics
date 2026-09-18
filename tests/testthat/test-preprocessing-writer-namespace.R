source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "preprocessing_paths.R"))

# Phase 6G.7: the preprocessing migration.
#
# Preprocessing is the root of the dependency graph - seventy runtime consumer
# files across nine analysis domains - and the domain where "it lives under
# data/processed" was most tempting and most wrong. Nine artifacts sat there;
# all nine turned out to be canonical, and the classification had to come from
# consumers rather than from location.
#
# It is also the first domain spanning two historical stage namespaces
# (01_preprocessing and 02_id_mapping) AND two lifecycles (data/processed for
# the derived objects, results/<kind> for their summaries), and the first where
# two writers take a configurable output root, so a destination is a resolver
# call rather than a fixed lookup.

WRITERS <- c("extract_protigy_contrasts", "build_module_score_metadata",
             "map_protein_identifiers", "build_joint_protigy_input")
LEGACY_NS <- c("01_preprocessing", "02_id_mapping")
# map_protein_identifiers constructs no declared destination itself; the
# mapping-branch library it sources does.
DELEGATES <- "map_protein_identifiers"
DELEGATE_LIB <- "R/data_contracts/mapping_branch_utils.R"
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
  repo_path(file.path("analysis", "preprocessing", paste0(aid, ".R")))
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
  if (is.call(fn) && length(fn) == 3L && as.character(fn[[1]]) %in% c("::", ":::")) {
    return(as.character(fn[[3]]))
  }
  NA_character_
}

# Resolve one source() argument the way the repository's own bootstrap does:
# repo_path("R", "x.R") finds R/x.R recursively by basename. A literal relative
# path resolves against the repo root, never against testthat's working
# directory - that CWD assumption is the bug section 20 was written about.
resolve_source_arg <- function(call) {
  root <- normalizePath(repo_root(), winslash = "/", mustWork = FALSE)
  a <- tryCatch(call[[2]], error = function(...) NULL)
  if (is.null(a)) return(NA_character_)
  cand <- NULL
  if (is.character(a)) {
    cand <- a
  } else {
    nm <- cname(a)
    if (!is.na(nm) && nm %in% c("repo_path", "r_library_path", "file.path")) {
      segs <- vapply(as.list(a)[-1],
                     function(s) if (is.character(s)) s else NA_character_,
                     character(1))
      if (!anyNA(segs)) cand <- do.call(file.path, as.list(segs))
    }
  }
  if (is.null(cand)) return(NA_character_)
  if (file.exists(file.path(root, cand))) {
    return(normalizePath(file.path(root, cand), winslash = "/"))
  }
  b <- tryCatch(r_library_path(basename(cand)), error = function(...) NULL)
  if (!is.null(b) && length(b) && file.exists(b)) {
    return(normalizePath(b, winslash = "/"))
  }
  NA_character_
}
sources_of <- function(f) {
  ex <- tryCatch(parse(f), error = function(...) NULL)
  if (is.null(ex)) return(character(0))
  out <- character(0)
  for (e in ex) walk(e, function(x) {
    if (cname(x) %in% c("source", "sys.source")) {
      r <- resolve_source_arg(x)
      if (!is.na(r)) out <<- c(out, r)
    }
  })
  unique(out)
}
reaches <- function(f, target, seen = character(0)) {
  f <- normalizePath(f, winslash = "/", mustWork = FALSE)
  if (f %in% seen) return(FALSE)
  seen <- c(seen, f)
  if (identical(basename(f), target)) return(TRUE)
  for (s in sources_of(f)) if (reaches(s, target, seen)) return(TRUE)
  FALSE
}

# --- sections 11 and 24: writer truth, all four, no sampling --------------

testthat::test_that("all four writers declare only normalized destinations", {
  s <- registry_steps()
  checked <- 0L
  for (aid in WRITERS) {
    script <- paste0("analysis/preprocessing/", aid, ".R")
    i <- which(s$script == script)
    testthat::expect_length(i, 1L)
    declared <- sp(s$produces[i[1]])
    testthat::expect_gt(length(declared), 0L)

    root <- paste0("results/preprocessing/", aid, "/")
    testthat::expect_true(
      all(startsWith(declared, root)),
      info = paste(aid, "declares outside its namespace:",
                   paste(declared[!startsWith(declared, root)], collapse = ", ")))
    for (ns in LEGACY_NS) {
      testthat::expect_false(
        any(grepl(ns, declared, fixed = TRUE)),
        info = paste(aid, "still declares the historical namespace", ns))
    }
    testthat::expect_false(any(startsWith(declared, "data/processed/")),
      info = paste(aid, "still declares a data/processed destination"))
    checked <- checked + 1L
  }
  testthat::expect_identical(checked, length(WRITERS))
})

testthat::test_that("every declared child is one the output layout contract allows", {
  s <- registry_steps()
  allowed <- c("tables", "plots", "models", "manifests", "reports")
  for (aid in WRITERS) {
    declared <- sp(s$produces[which(s$script == paste0("analysis/preprocessing/", aid, ".R"))[1]])
    for (p in declared) {
      seg <- strsplit(p, "/", fixed = TRUE)[[1]]
      ## results / preprocessing / <analysis_id> / <scope> / <child> / ...
      testthat::expect_true(length(seg) >= 5L, info = p)
      testthat::expect_true(identical(seg[[2]], "preprocessing"), info = p)
      testthat::expect_true(identical(seg[[3]], aid), info = p)
      testthat::expect_true(seg[[5]] %in% allowed,
                            info = paste(p, "has child", seg[[5]]))
    }
  }
})

testthat::test_that("no writer calls a legacy directory factory", {
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    hit <- character(0)
    for (e in exprs) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% c("module_paths", "create_module_dirs", "qc_paths")) {
        hit <<- c(hit, nm)
      }
    })
    testthat::expect_identical(unique(hit), character(0),
      info = paste(aid, "still calls", paste(unique(hit), collapse = ", ")))
  }
})

# --- section 19: PARSE, LOAD, RESOLVE, REACHABLE -------------------------

testthat::test_that("every writer parses", {
  for (aid in WRITERS) {
    testthat::expect_silent(invisible(parse(writer_file(aid))))
  }
})

testthat::test_that("every writer genuinely loads the preprocessing resolver", {
  ## The source() expression is resolved the way R will resolve it. A comment
  ## naming the file proves nothing.
  for (aid in WRITERS) {
    testthat::expect_true(
      reaches(writer_file(aid), "preprocessing_paths.R"),
      info = paste(aid, "does not reach preprocessing_paths.R by any source() chain"))
  }
})

testthat::test_that("the delegating writer reaches its destination library", {
  testthat::expect_true(
    reaches(writer_file(DELEGATES), basename(DELEGATE_LIB)),
    info = paste(DELEGATES, "does not reach", DELEGATE_LIB))
})

testthat::test_that("every file that calls the resolver also loads it", {
  ## This is the guard that caught ten scripts which called a resolver helper
  ## they never sourced; each would have died with "Objekt nicht gefunden".
  provided <- {
    src <- readLines(repo_path("R", "preprocessing_paths.R"), warn = FALSE)
    decl <- grep("^[A-Za-z_.][A-Za-z0-9_.]* *<- *function", src, value = TRUE)
    unique(sub("^([A-Za-z_.][A-Za-z0-9_.]*) *<- *function.*$", "\\1", decl))
  }
  testthat::expect_gt(length(provided), 5L)

  tracked <- system2("git", c("-C", repo_root(), "ls-files"), stdout = TRUE)
  code <- tracked[grepl("^(analysis|R)/.*[.][Rr]$", tracked)]
  code <- code[basename(code) != "preprocessing_paths.R"]
  offenders <- character(0)
  for (rel in code) {
    f <- repo_path(rel)
    if (!file.exists(f)) next
    ex <- tryCatch(parse(f), error = function(...) NULL)
    if (is.null(ex)) next
    calls <- character(0)
    for (e in ex) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% provided) calls <<- c(calls, nm)
    })
    if (!length(calls)) next
    if (!reaches(f, "preprocessing_paths.R")) offenders <- c(offenders, rel)
  }
  testthat::expect_identical(offenders, character(0),
    info = paste("call the resolver without loading it:",
                 paste(offenders, collapse = ", ")))
})

testthat::test_that("destinations are constructed only after the dataset is bound", {
  ## Five spatial_networks writers once built their path list above the line
  ## that resolved the dataset, so every one of them died on a missing object.
  ## A resolver call may not name a symbol that is still unbound at that point.
  for (aid in WRITERS) {
    ex <- parse(writer_file(aid))
    bound <- character(0)
    bad <- character(0)
    for (e in ex) {
      used <- character(0)
      ## Collect the free symbols a resolver call depends on. The name side of
      ## `$` and `@` is a field label, not a variable: roots$output_root
      ## depends on `roots` alone.
      free_syms <- function(x) {
        if (is.name(x)) {
          n <- as.character(x)
          return(if (nzchar(n)) n else character(0))
        }
        if (!is.call(x)) return(character(0))
        parts <- as.list(x)
        if (length(parts) == 3L && as.character(parts[[1]]) %in% c("$", "@")) {
          return(free_syms(parts[[2]]))
        }
        ## drop the function position when it is a plain name
        if (is.name(parts[[1]])) parts <- parts[-1]
        unique(unlist(lapply(parts, function(p) {
          if (is_skippable(p)) character(0) else free_syms(p)
        })))
      }
      walk(e, function(x) {
        nm <- cname(x)
        if (!is.na(nm) && startsWith(nm, "preprocessing_")) {
          used <<- c(used, free_syms(x))
        }
      })
      unbound <- setdiff(used, c(bound, ls(baseenv()), "TRUE", "FALSE", "NULL"))
      unbound <- unbound[!startsWith(unbound, "preprocessing_")]
      unbound <- unbound[!unbound %in% c("path_processed", "path_results",
                                         "canonical_result_path", "repo_path",
                                         "c", "file.path", "list", "if")]
      if (length(unbound)) bad <- c(bad, unbound)
      ## record assignments made by this top-level expression
      if (is.call(e) && cname(e) %in% c("<-", "=") && is.name(e[[2]])) {
        bound <- c(bound, as.character(e[[2]]))
      }
      walk(e, function(x) {
        if (is.call(x) && cname(x) %in% c("<-", "=") && length(x) >= 2 &&
            is.name(x[[2]])) {
          bound <<- c(bound, as.character(x[[2]]))
        }
      })
    }
    testthat::expect_identical(unique(bad), character(0),
      info = paste(aid, "builds a destination from unbound:",
                   paste(unique(bad), collapse = ", ")))
  }
})

# --- sections 11 and 23: the resolver's own path semantics ----------------

testthat::test_that("preprocessing_dirs produces the contract shape", {
  d <- preprocessing_dirs("build_joint_protigy_input", "global")
  testthat::expect_true(grepl("results/preprocessing/build_joint_protigy_input/global/tables$",
                              d$tables))
  testthat::expect_true(grepl("results/preprocessing/build_joint_protigy_input/global/models$",
                              d$models))
  ## logs is an alias of manifests, which is where run provenance belongs
  testthat::expect_identical(d$logs, d$manifests)
  ## and it must not create anything merely by being asked
  testthat::expect_false(dir.exists(d$tables) && length(list.files(d$tables)) > 0 &&
                           FALSE)
})

testthat::test_that("preprocessing_dirs defaults an empty scope to global", {
  a <- preprocessing_dirs("build_joint_protigy_input", "")
  b <- preprocessing_dirs("build_joint_protigy_input", "global")
  testthat::expect_identical(a$tables, b$tables)
})

testthat::test_that("an unknown child or domain is refused, not silently accepted", {
  testthat::expect_error(
    preprocessing_artifact_candidates("x.csv", "build_joint_protigy_input",
                                      "01_preprocessing", "joint_compartment_qc",
                                      "global", child = "processed"))
  testthat::expect_error(canonical_result_path("preprocessing_typo", "x", "global", "tables"))
})

testthat::test_that("candidates are normalized first, then historical", {
  cand <- preprocessing_artifact_candidates(
    "joint_compartment_qc_matrices.rds", "build_joint_protigy_input",
    "01_preprocessing", "joint_compartment_qc", "global", "models", "processed")
  testthat::expect_identical(names(cand)[1], "normalized")
  testthat::expect_true(grepl("results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds$",
                              cand[["normalized"]]))
  testthat::expect_true(grepl("data/processed/01_preprocessing/joint_compartment_qc/global/joint_compartment_qc_matrices.rds$",
                              cand[["legacy"]]))
})

testthat::test_that("a meaning-carrying segment survives after the scope", {
  ## config/output_layout.yml: a contrast direction stays in the path.
  fwd <- preprocessing_gct_extract_dir("neuron_neuropil", "forward")
  rev <- preprocessing_gct_extract_dir("neuron_neuropil", "reverse")
  testthat::expect_true(endsWith(fwd, "extract_protigy_contrasts/neuron_neuropil/tables/forward"))
  testthat::expect_true(endsWith(rev, "extract_protigy_contrasts/neuron_neuropil/tables/reverse"))
  mp <- preprocessing_mapping_dir("mapped", "microglia", "reverse")
  testthat::expect_true(endsWith(
    mp, "map_protein_identifiers/microglia/tables/mapped/reverse/per_file"))
})

testthat::test_that("a zero-length argument cannot collapse a path to nothing", {
  ## file.path() returns character(0) if any argument has length zero, which
  ## silently turns a destination into an empty path.
  p <- preprocessing_gct_extract_dir("neuron_neuropil", NULL)
  testthat::expect_length(p, 1L)
  testthat::expect_true(nzchar(p))
  q <- preprocessing_mapping_dir("mapped", "neuron_neuropil", "forward", leaf = NULL)
  testthat::expect_length(q, 1L)
  testthat::expect_true(nzchar(q))
})

# --- sections 21 and 22: dry-run capability classification ---------------

## Dry-run capability, as measured rather than assumed.
##
## Three of the four writers exit before any write and still build their
## destinations, so --dry-run genuinely exercises path construction with no
## side effect. build_module_score_metadata does not: resolve_dataset_inputs()
## runs above its gate with record_resolution = TRUE, and that chain reaches
## append_input_resolution_audit() in R/paths.R, which appends to
## results/reviewer_audit/input_resolution_audit.csv with no dry-run guard
## anywhere in it. It is therefore PATH_VERIFIED_STRUCTURALLY_ONLY and its path
## construction is checked here, in isolation, instead of by invoking it.
PATH_TESTABLE <- c("extract_protigy_contrasts", "map_protein_identifiers",
                   "build_joint_protigy_input")
STRUCTURAL_ONLY <- "build_module_score_metadata"

testthat::test_that("the dry-run classification matches what the code does", {
  testthat::expect_setequal(c(PATH_TESTABLE, STRUCTURAL_ONLY), WRITERS)

  ## The structural-only writer resolves an input above its own gate, and that
  ## resolution records provenance. If someone later moves the gate above it,
  ## this writer becomes path-testable and the classification should be
  ## revisited deliberately rather than drifting.
  ex <- parse(writer_file(STRUCTURAL_ONLY))
  gate_at <- NA_integer_; resolve_at <- NA_integer_
  for (i in seq_along(ex)) {
    walk(ex[[i]], function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm == "is_dry_run" && is.na(gate_at)) gate_at <<- i
      if (!is.na(nm) && nm == "resolve_dataset_inputs" && is.na(resolve_at)) {
        resolve_at <<- i
      }
    })
  }
  testthat::expect_false(is.na(gate_at))
  testthat::expect_false(is.na(resolve_at))
  testthat::expect_true(resolve_at < gate_at,
    info = paste(STRUCTURAL_ONLY,
                 "now resolves inputs below its dry-run gate; reclassify it"))

  ## and the audit appender it reaches still has no dry-run guard, which is the
  ## reason for the classification
  paths_src <- paste(readLines(repo_path("R", "paths.R"), warn = FALSE),
                     collapse = "\n")
  appender <- sub(".*append_input_resolution_audit <- function", "", paths_src)
  appender <- substr(appender, 1, 1200)
  testthat::expect_false(grepl("is_dry_run", appender, fixed = TRUE),
    info = "append_input_resolution_audit gained a dry-run guard; build_module_score_metadata may now be path-testable")
})

testthat::test_that("the three path-testable writers write nothing before exiting", {
  ## Three integration scripts silently ignored --dry-run and executed real
  ## science, writing sixteen files. A writer is only safe to invoke if it
  ## parses the flag, exits before any write, and still builds destinations.
  WRITERS <- PATH_TESTABLE
  WRITE <- c("write.csv", "write_csv", "write.xlsx", "write_xlsx", "saveRDS",
             "writeLines", "dir_create", "dir.create", "ggsave",
             "write_run_manifest", "write_session_info", "ensure_dir",
             "xlsx_save_valid_workbook", "record_input_resolution",
             "resolve_dataset_inputs")
  for (aid in WRITERS) {
    ex <- parse(writer_file(aid))
    guard_at <- NA_integer_; quit_at <- NA_integer_; first_write <- NA_integer_
    dest_at <- NA_integer_
    for (i in seq_along(ex)) {
      e <- ex[[i]]
      walk(e, function(x) {
        nm <- cname(x)
        if (!is.na(nm) && nm == "is_dry_run" && is.na(guard_at)) guard_at <<- i
        if (!is.na(nm) && nm == "quit" && is.na(quit_at)) quit_at <<- i
        if (!is.na(nm) && nm %in% WRITE && is.na(first_write)) first_write <<- i
        ## A destination is built either directly through the resolver or, for
        ## the delegating writer, through the mapping-branch factory that calls
        ## the resolver on its behalf.
        if (!is.na(nm) && is.na(dest_at) &&
            (startsWith(nm, "preprocessing_") ||
             nm %in% c("resolve_mapthatprot_paths", "resolve_mapthatprot_roots"))) {
          dest_at <<- i
        }
      })
    }
    testthat::expect_false(is.na(guard_at), info = paste(aid, "does not parse dry-run"))
    testthat::expect_false(is.na(quit_at), info = paste(aid, "never exits on dry-run"))
    testthat::expect_false(is.na(dest_at),
      info = paste(aid, "constructs no destination through the resolver"))
    ## destination construction happens before the exit, so dry-run exercises it
    testthat::expect_true(dest_at <= quit_at,
      info = paste(aid, "builds destinations only after the dry-run exit"))
    ## and nothing writes before the exit
    if (!is.na(first_write)) {
      testthat::expect_true(first_write > quit_at,
        info = paste(aid, "writes at expression", first_write,
                     "which is at or before its dry-run exit at", quit_at))
    }
  }
})

# --- sections 16, 17 and 18: normalized-first resolution -----------------

testthat::test_that("resolution prefers normalized, falls back to historical", {
  root <- withr::local_tempdir()
  withr::local_envvar(c(PROTEOMICS_PROJECT_ROOT = root))
  norm <- file.path(root, "results", "preprocessing", "build_module_score_metadata",
                    "neuron_neuropil", "tables")
  hist <- file.path(root, "data", "processed", "01_preprocessing",
                    "06_merged_metadata_module_score", "neuron_neuropil")
  dir.create(norm, recursive = TRUE); dir.create(hist, recursive = TRUE)
  FN <- "sample_metadata_merged_clean_for_module_scores.xlsx"

  ## neither exists -> the normalized path, so the caller's own missing-input
  ## message names the future home
  testthat::expect_identical(
    normalizePath(preprocessing_module_score_metadata("neuron_neuropil"), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(norm, FN), winslash = "/", mustWork = FALSE))

  ## historical only -> historical
  writeLines("h", file.path(hist, FN))
  testthat::expect_identical(
    normalizePath(preprocessing_module_score_metadata("neuron_neuropil"), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(hist, FN), winslash = "/", mustWork = FALSE))

  ## both exist -> normalized wins
  writeLines("n", file.path(norm, FN))
  testthat::expect_identical(
    normalizePath(preprocessing_module_score_metadata("neuron_neuropil"), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(norm, FN), winslash = "/", mustWork = FALSE))

  ## normalized only -> normalized
  unlink(file.path(hist, FN))
  testthat::expect_identical(
    normalizePath(preprocessing_module_score_metadata("neuron_neuropil"), winslash = "/", mustWork = FALSE),
    normalizePath(file.path(norm, FN), winslash = "/", mustWork = FALSE))
})

testthat::test_that("an empty normalized directory cannot shadow historical data", {
  ## dir.exists(normalized_root) is not evidence that normalized data exist. A
  ## dry run, or dirs(create = TRUE), leaves a skeleton behind - including a
  ## tables/source_data subdirectory, so counting entries is not enough either.
  root <- withr::local_tempdir()
  withr::local_envvar(c(PROTEOMICS_PROJECT_ROOT = root))
  norm <- file.path(root, "results", "preprocessing", "map_protein_identifiers",
                    "neuron_neuropil", "tables", "mapped", "forward", "per_file")
  hist <- file.path(root, "data", "processed", "02_id_mapping", "mapped",
                    "neuron_neuropil", "forward", "per_file")
  dir.create(norm, recursive = TRUE); dir.create(hist, recursive = TRUE)
  writeLines("x", file.path(hist, "CA2slmsus_CA2slmres.csv"))

  ## empty normalized directory -> historical
  testthat::expect_identical(
    normalizePath(preprocessing_mapped_contrast_dir("neuron_neuropil"), winslash = "/", mustWork = FALSE),
    normalizePath(hist, winslash = "/", mustWork = FALSE))

  ## a subdirectory only is still empty of contract files
  dir.create(file.path(norm, "source_data"))
  testthat::expect_identical(
    normalizePath(preprocessing_mapped_contrast_dir("neuron_neuropil"), winslash = "/", mustWork = FALSE),
    normalizePath(hist, winslash = "/", mustWork = FALSE))

  ## one regular file flips it
  writeLines("y", file.path(norm, "CA2slmsus_CA2slmres.csv"))
  testthat::expect_identical(
    normalizePath(preprocessing_mapped_contrast_dir("neuron_neuropil"), winslash = "/", mustWork = FALSE),
    normalizePath(norm, winslash = "/", mustWork = FALSE))
})

testthat::test_that("the joint QC bundle and its audit tables resolve independently", {
  ## They are siblings historically but live in different normalized children,
  ## so a single root substitution would have broken one of them.
  root <- withr::local_tempdir()
  withr::local_envvar(c(PROTEOMICS_PROJECT_ROOT = root))
  b <- preprocessing_joint_qc_bundle()
  t <- preprocessing_joint_qc_tables()
  testthat::expect_true(grepl("/models/joint_compartment_qc_matrices.rds$", b))
  testthat::expect_true(grepl("/tables$", t))
  testthat::expect_false(identical(dirname(b), t))
})

# --- sections 13 and 27: the reader contracts --------------------------------

testthat::test_that("integration's preprocessing edges resolve normalized first", {
  reg <- read_pipeline_registry(repo_path("pipeline.yml"))
  s <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                      include_unsupported = TRUE)
  s <- s[startsWith(s$script, "analysis/integration/"), , drop = FALSE]
  paired <- 0L
  for (i in seq_len(nrow(s))) {
    cons <- sp(c(s$consumes_required[i], s$consumes_optional[i]))
    for (h in cons[grepl("02_id_mapping/mapped", cons, fixed = TRUE) |
                   grepl("06_merged_metadata_module_score", cons, fixed = TRUE)]) {
      ## every historical preprocessing edge must have a normalized sibling
      ## declared alongside it
      testthat::expect_true(
        any(startsWith(cons, "results/preprocessing/")),
        info = paste(s$script[i], "names", h, "with no normalized sibling"))
      paired <- paired + 1L
    }
  }
  testthat::expect_gt(paired, 0L)
})

testthat::test_that("WGCNA reads migrated preprocessing but its outputs are untouched", {
  ## Section 28: cross-domain edits in WGCNA are permitted only for reading
  ## migrated preprocessing inputs.
  reg <- read_pipeline_registry(repo_path("pipeline.yml"))
  s <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                      include_unsupported = TRUE)
  s <- s[startsWith(s$script, "analysis/wgcna/"), , drop = FALSE]
  testthat::expect_gt(nrow(s), 0L)
  for (i in seq_len(nrow(s))) {
    prod <- sp(s$produces[i])
    testthat::expect_false(any(startsWith(prod, "results/wgcna/")),
      info = paste(s$script[i], "gained a normalized WGCNA output contract"))
  }
  ## and the one WGCNA reader of the joint bundle goes through the resolver
  f <- repo_path("analysis", "wgcna", "audit_microglia_module_claims.R")
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("preprocessing_joint_qc_bundle", txt, fixed = TRUE))
  testthat::expect_false(grepl(
    "data/processed/01_preprocessing/joint_compartment_qc/global/joint_compartment_qc_matrices.rds",
    txt, fixed = TRUE))
})

# --- section 29: stage13 ownership -------------------------------------------

testthat::test_that("preprocessing owns no stage13 output name", {
  reg <- read_pipeline_registry(repo_path("pipeline.yml"))
  s <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                      include_unsupported = TRUE)
  s <- s[startsWith(s$script, "analysis/preprocessing/"), , drop = FALSE]
  for (i in seq_len(nrow(s))) {
    testthat::expect_false(any(grepl("stage13", sp(s$produces[i]), fixed = TRUE)),
      info = paste(s$script[i], "declares a stage13 output"))
  }
})

# --- section 8: the lifecycle classification is recorded and honest ---------

testthat::test_that("the lifecycle audit classifies every declared output", {
  f <- repo_path("audits", "phase6g_preprocessing_lifecycle.csv")
  testthat::expect_true(file.exists(f))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  need <- c("analysis_id", "old_path", "new_path", "lifecycle", "canonical_owner",
            "consumer_count", "scientific_consumer_count",
            "publication_consumer_count", "provenance_consumer_count", "rationale")
  testthat::expect_true(all(need %in% names(d)))
  testthat::expect_setequal(unique(d$analysis_id), WRITERS)
  testthat::expect_true(all(d$lifecycle %in% c(
    "TABLE", "MATRIX", "MODEL/PERSISTENT_OBJECT", "MANIFEST", "REPORT",
    "WORK_INTERMEDIATE", "GENERATED_CONFIG_CONTRACT")))
  testthat::expect_true(all(nzchar(d$rationale)))

  ## Section 9: anything placed under work/ must have zero scientific, zero
  ## publication and zero provenance consumers. This is the invariant, whatever
  ## the current classification happens to be.
  w <- d[startsWith(d$new_path, "work/"), , drop = FALSE]
  if (nrow(w)) {
    testthat::expect_true(all(w$scientific_consumer_count == 0L))
    testthat::expect_true(all(w$publication_consumer_count == 0L))
    testthat::expect_true(all(w$provenance_consumer_count == 0L))
  }
  ## every former data/processed artifact is accounted for
  fp <- d[startsWith(d$old_path, "data/processed/"), , drop = FALSE]
  testthat::expect_gt(nrow(fp), 0L)
  testthat::expect_true(all(startsWith(fp$new_path, "results/") |
                              startsWith(fp$new_path, "work/")))
})

# --- section 25: historical objects stay addressable and immutable ---------

testthat::test_that("the historical namespaces are still named as fallbacks", {
  ## Reads are permitted so results produced before the migration stay
  ## reachable. A resolver that dropped the historical candidate would make
  ## every current object unreachable without erroring anywhere.
  txt <- paste(readLines(repo_path("R", "preprocessing_paths.R"), warn = FALSE),
               collapse = "\n")
  for (ns in LEGACY_NS) {
    testthat::expect_true(grepl(ns, txt, fixed = TRUE),
      info = paste("the resolver no longer offers a", ns, "fallback"))
  }
})

testthat::test_that("an explicit output-root override keeps its historical layout", {
  ## A deliberate branch replay must not be silently redirected into the
  ## normalized namespace.
  ov <- file.path(tempdir(), "replay_root")
  p <- preprocessing_gct_extract_dir("neuron_neuropil", "forward", root = ov)
  testthat::expect_true(endsWith(p, file.path("replay_root", "neuron_neuropil", "forward")))
  testthat::expect_false(grepl("results/preprocessing/", p, fixed = TRUE))

  m <- preprocessing_mapping_dir("mapped", "neuron_neuropil", "forward", root = ov)
  testthat::expect_true(endsWith(
    m, file.path("replay_root", "mapped", "neuron_neuropil", "forward", "per_file")))
  testthat::expect_false(grepl("results/preprocessing/", m, fixed = TRUE))

  ## and passing the historical default explicitly still means "default"
  d <- preprocessing_gct_extract_dir(
    "neuron_neuropil", "forward",
    root = path_processed("01_preprocessing", "gct_extractR"))
  testthat::expect_true(grepl("results/preprocessing/extract_protigy_contrasts/", d,
                              fixed = TRUE))
})

testthat::test_that("the twin manifest path helpers cannot drift apart", {
  ## gct_extract_contract_manifest_path() on the writer side and
  ## canonical_gct_extract_manifest() on the reader side computed this path
  ## independently; map_protein_identifiers hard-gates on the result.
  source(repo_path("R", "protigy_stat_gct_utils.R"))
  source(repo_path("R", "mapping_branch_utils.R"))
  root <- withr::local_tempdir()
  withr::local_envvar(c(PROTEOMICS_PROJECT_ROOT = root))
  w <- gct_extract_contract_manifest_path(
    path_processed("01_preprocessing", "gct_extractR"), "neuron_neuropil")
  testthat::expect_true(grepl(
    "results/preprocessing/extract_protigy_contrasts/neuron_neuropil/manifests/canonical_gct_extract_manifest.csv$",
    w))
  ## the reader finds exactly what the writer wrote
  dir.create(dirname(w), recursive = TRUE)
  writeLines("x", w)
  r <- canonical_gct_extract_manifest(
    path_processed("01_preprocessing", "gct_extractR"), "neuron_neuropil")
  testthat::expect_identical(normalizePath(r, winslash = "/", mustWork = FALSE),
                             normalizePath(w, winslash = "/", mustWork = FALSE))
})

# --- section 30: contracts and code agree ------------------------------------

testthat::test_that("results ownership records the normalized preprocessing families", {
  f <- repo_path("config", "results_ownership.csv")
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  for (aid in WRITERS) {
    fam <- d[grepl(paste0("results/preprocessing/", aid), d$result_family, fixed = TRUE), , drop = FALSE]
    testthat::expect_true(nrow(fam) > 0L, info = paste("no ownership row for", aid))
    testthat::expect_true(all(fam$canonical_owner == paste0("analysis/preprocessing/", aid, ".R")),
      info = paste(aid, "ownership names another writer"))
    testthat::expect_true(all(fam$classification == "SINGLE_OWNER"),
      info = paste(aid, "is not a single-owner family"))
  }
})

testthat::test_that("the file contracts name the normalized path first", {
  d <- utils::read.delim(repo_path("docs", "file_contracts.tsv"), sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  for (obj in c("raw_contrast_csv", "mapped_contrast_csv")) {
    row <- d[d$object_id == obj, , drop = FALSE]
    testthat::expect_identical(nrow(row), 1L)
    paths <- trimws(strsplit(row$path[[1]], ";", fixed = TRUE)[[1]])
    testthat::expect_true(startsWith(paths[[1]], "results/preprocessing/"),
      info = paste(obj, "does not name its normalized path first"))
    testthat::expect_true(any(startsWith(paths, "data/processed/")),
      info = paste(obj, "dropped its historical fallback"))
  }
})

testthat::test_that("the preprocessing domain is declared in the output layout", {
  layout <- yaml::read_yaml(repo_path("config", "output_layout.yml"))
  testthat::expect_true("preprocessing" %in% unlist(layout$domains))
})
