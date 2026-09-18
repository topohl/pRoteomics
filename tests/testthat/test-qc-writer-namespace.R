source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "qc_result_paths.R"))
source(repo_path("R", "migration_gate_utils.R"))

# Phase 6G.5: the QC migration.
#
# QC is upstream of everything: forty analyses outside analysis/qc read its
# outputs, and each built those paths itself, so this is the first domain where
# the read side had to be given a resolver rather than converted at a
# chokepoint.
#
# It is also the domain that held the last coordinating family in the
# repository, and the only one that legitimately writes a file outside
# results/: a generated configuration contract.

WRITERS <- c("assess_dataset_quality", "assess_joint_compartment_quality",
             "assess_marker_rank_abundance", "assess_pca_confounding",
             "assess_replicate_consistency", "assess_sample_quality",
             "build_reference_marker_registry", "discover_empirical_roi_markers",
             "export_marker_traits", "partition_variance",
             "render_compartment_abundance_figures",
             "render_joint_compartment_qc_figures",
             "summarize_marker_detectability", "summarize_missingness",
             "summarize_qc_confounding")
LEGACY_NS <- "03_qc_exploration"
GENERATED_CONFIG <- "config/marker_panels/wgcna_reference_marker_sets.csv"
# render_compartment_abundance_figures constructs no destination itself; the
# workflow library it sources does.
DELEGATES <- "render_compartment_abundance_figures"
DELEGATE_LIB <- "R/qc/control_compartment_abundance_workflow.R"

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
  repo_path(file.path("analysis", "qc", paste0(aid, ".R")))
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
syms <- function(e) {
  out <- character(0)
  walk(e, function(x) if (is.name(x)) {
    n <- as.character(x); if (nzchar(n)) out <<- c(out, n)
  })
  unique(out)
}

# --- section 20: writer truth, all fifteen, no sampling -------------------

testthat::test_that("all fifteen writers declare only normalized destinations", {
  s <- registry_steps()
  checked <- 0L
  for (aid in WRITERS) {
    script <- paste0("analysis/qc/", aid, ".R")
    i <- which(s$script == script)
    testthat::expect_length(i, 1L)
    declared <- sp(s$produces[i[1]])
    testthat::expect_gt(length(declared), 0L)

    ## the adjudicated generated configuration contract is exempt
    results_declared <- setdiff(declared, GENERATED_CONFIG)
    root <- paste0("results/qc/", aid, "/")
    testthat::expect_true(
      all(startsWith(results_declared, root)),
      info = paste(aid, "declares outside its namespace:",
                   paste(results_declared[!startsWith(results_declared, root)],
                         collapse = ", ")))
    testthat::expect_false(any(grepl(LEGACY_NS, results_declared, fixed = TRUE)),
                           info = paste(aid, "still declares the historical namespace"))
    checked <- checked + 1L
  }
  testthat::expect_identical(checked, length(WRITERS))
})

testthat::test_that("no writer calls a legacy directory factory", {
  ## qc_paths() was the shared factory and no longer exists; module_paths()
  ## and create_module_dirs() are the other two.
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    hit <- character(0)
    for (e in exprs) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% c("qc_paths", "module_paths", "create_module_dirs")) {
        hit <<- c(hit, nm)
      }
    })
    testthat::expect_identical(unique(hit), character(0),
      info = paste(aid, "still calls", paste(unique(hit), collapse = ", ")))
  }
})

testthat::test_that("qc_paths is gone, so nothing can quietly keep using it", {
  txt <- paste(readLines(repo_path("R", "qc_result_paths.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_false(grepl("\nqc_paths <- function", txt, fixed = TRUE))
  utils <- paste(readLines(repo_path("R", "qc_exploration_utils.R"), warn = FALSE),
                 collapse = "\n")
  testthat::expect_false(grepl("\nqc_paths <- function", utils, fixed = TRUE))
})

# --- sections 16 and 17: LOAD and REACHABLE ------------------------------

testthat::test_that("every writer that resolves paths actually loads the resolver", {
  ## The source() expression is evaluated the way R will evaluate it. A comment
  ## naming the file proves nothing.
  users <- 0L
  for (aid in WRITERS) {
    f <- writer_file(aid)
    exprs <- parse(f)
    uses <- FALSE
    for (e in exprs) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% c("qc_dirs", "qc_find", "qc_dir_any",
                                  "qc_artifact_candidates")) uses <<- TRUE
    })
    if (!uses) next
    users <- users + 1L

    loaded <- FALSE
    for (e in exprs) walk(e, function(x) {
      if (loaded || !is.call(x)) return(invisible(NULL))
      nm <- cname(x)
      if (is.na(nm) || !nm %in% c("source", "sys.source")) return(invisible(NULL))
      p <- tryCatch(eval(x[[2]], envir = globalenv()), error = function(...) NA_character_)
      if (length(p) == 1L && !is.na(p) && file.exists(p) &&
          basename(p) %in% c("qc_result_paths.R", "qc_exploration_utils.R")) {
        loaded <<- TRUE
      }
      invisible(NULL)
    })
    testthat::expect_true(loaded,
      info = paste(aid, "calls a QC path resolver without loading it"))
  }
  ## fourteen of fifteen; the delegating writer constructs nothing itself
  testthat::expect_identical(users, length(WRITERS) - 1L)
})

testthat::test_that("the delegating writer's library resolves normalized", {
  ## render_compartment_abundance_figures has no write call and no path
  ## construction. Its destinations come from the workflow library, so that is
  ## where the assertion belongs.
  f <- repo_path(DELEGATE_LIB)
  testthat::expect_true(file.exists(f))
  exprs <- parse(f)
  uses_api <- FALSE
  names_legacy <- FALSE
  for (e in exprs) walk(e, function(x) {
    nm <- cname(x)
    if (!is.na(nm) && nm %in% c("qc_dirs", "canonical_module_dirs",
                                "canonical_result_path")) uses_api <<- TRUE
    if (!is.na(nm) && nm %in% c("path_results", "path_processed")) {
      lits <- unlist(lapply(as.list(x)[-1], function(a) if (is.character(a)) a else NULL))
      if (any(lits == LEGACY_NS)) names_legacy <<- TRUE
    }
  })
  testthat::expect_true(uses_api)
  testthat::expect_false(names_legacy)
  ## and the entrypoint really does source it
  src <- paste(readLines(writer_file(DELEGATES), warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("control_compartment_abundance_workflow.R", src, fixed = TRUE))
})

testthat::test_that("each destination call runs after the variables it needs", {
  DEST <- c("qc_dirs", "canonical_module_dirs", "canonical_result_path")
  checked <- 0L
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    assigned_at <- list()
    for (k in seq_along(exprs)) {
      e <- exprs[[k]]
      if (is.call(e) && length(e) == 3L &&
          as.character(e[[1]])[1] %in% c("<-", "=", "<<-") && is.name(e[[2]])) {
        nm <- as.character(e[[2]])
        if (is.null(assigned_at[[nm]])) assigned_at[[nm]] <- k
      }
    }
    idx <- which(vapply(seq_along(exprs),
                        function(k) any(DEST %in% syms(exprs[[k]])), logical(1)))
    for (j in idx) {
      deps <- setdiff(syms(exprs[[j]]),
                      c(DEST, "PATHS", "CANONICAL_PATHS", "out_root", "early_paths"))
      for (d in deps) {
        at <- assigned_at[[d]]
        if (is.null(at)) next
        ## at == j is the enclosing assignment or a function definition
        testthat::expect_false(
          at > j,
          info = paste0(aid, ": '", d, "' assigned at expression ", at,
                        " but the destination call runs at ", j))
        checked <- checked + 1L
      }
    }
  }
  testthat::expect_gt(checked, 0L)
})

# --- section 23: normalized-first, all four existence states -------------

testthat::test_that("qc_find prefers normalized and falls back to historical", {
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("qc_", Sys.getpid()))
  owner <- "assess_pca_confounding"
  sub <- "05_pca_confounding_qc"
  ds <- "neuron_neuropil"
  fn <- "PCA_confounding_summary.csv"
  nd <- file.path(tmp, "results", "qc", owner, ds, "tables")
  ld <- file.path(tmp, "results", "tables", LEGACY_NS, sub, ds)
  dir.create(nd, recursive = TRUE, showWarnings = FALSE)
  dir.create(ld, recursive = TRUE, showWarnings = FALSE)
  nf <- file.path(nd, fn); lf <- file.path(ld, fn)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  combos <- list(
    list(n = FALSE, l = TRUE,  want = "legacy"),
    list(n = TRUE,  l = TRUE,  want = "normalized"),
    list(n = TRUE,  l = FALSE, want = "normalized"),
    list(n = FALSE, l = FALSE, want = "normalized"))

  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    for (cb in combos) {
      unlink(c(nf, lf))
      if (cb$n) writeLines("x", nf)
      if (cb$l) writeLines("x", lf)
      got <- qc_find(fn, owner = owner, legacy_substep = sub, scope = ds)
      is_norm <- grepl("/results/qc/", got, fixed = TRUE)
      if (identical(cb$want, "normalized")) {
        testthat::expect_true(is_norm,
          info = paste("normalized =", cb$n, "legacy =", cb$l, "->", got))
      } else {
        testthat::expect_false(is_norm,
          info = paste("normalized =", cb$n, "legacy =", cb$l, "->", got))
      }
    }
  })
})

testthat::test_that("the flat historical shape is found too", {
  ## qc_paths() appended the dataset, but the two substeps that took no dataset
  ## wrote straight into the substep directory. Eight readers across four
  ## domains address the empirical ROI marker sets with no scope segment.
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("qcflat_", Sys.getpid()))
  sub <- "05_empirical_roi_marker_discovery"
  fn <- "empirical_roi_marker_sets.csv"
  flat <- file.path(tmp, "results", "tables", LEGACY_NS, sub)
  dir.create(flat, recursive = TRUE, showWarnings = FALSE)
  writeLines("x", file.path(flat, fn))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    got <- qc_find(fn, owner = "discover_empirical_roi_markers", legacy_substep = sub)
    testthat::expect_true(grepl(sub, got, fixed = TRUE))
    testthat::expect_false(grepl("/results/qc/", got, fixed = TRUE))
  })
})

# --- section 33: an empty skeleton must not shadow real data -------------

testthat::test_that("a normalized directory holding only subdirectories does not shadow", {
  ## This is the concrete defect: qc_dirs(create = TRUE) leaves a
  ## tables/source_data subdirectory, and counting directory entries made the
  ## empty normalized directory look populated. Seven real inputs were then
  ## reported missing.
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("qcshadow_", Sys.getpid()))
  owner <- "assess_joint_compartment_quality"
  sub <- "00b_joint_compartment_qc"
  nd <- file.path(tmp, "results", "qc", owner, "global", "tables")
  dir.create(file.path(nd, "source_data"), recursive = TRUE, showWarnings = FALSE)
  ld <- file.path(tmp, "results", "tables", LEGACY_NS, sub, "global")
  dir.create(ld, recursive = TRUE, showWarnings = FALSE)
  writeLines("x", file.path(ld, "joint_primary_pca_scores.csv"))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    got <- qc_dir_any(owner = owner, legacy_substep = sub)
    testthat::expect_true(grepl(sub, got, fixed = TRUE),
      info = paste("resolved to", got, "but the only real file is historical"))
  })

  ## and once a real file lands in the normalized directory it wins
  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    writeLines("x", file.path(nd, "joint_primary_pca_scores.csv"))
    got2 <- qc_dir_any(owner = owner, legacy_substep = sub)
    testthat::expect_true(grepl("/results/qc/", got2, fixed = TRUE))
  })
})

# --- section 11: the generated configuration contract --------------------

testthat::test_that("the generated marker registry stays configuration", {
  s <- registry_steps()
  i <- which(s$script == "analysis/qc/build_reference_marker_registry.R")
  declared <- sp(s$produces[i[1]])
  testthat::expect_true(GENERATED_CONFIG %in% declared)

  ## it is committed, which is what makes it a build-then-promote contract
  tracked <- system2("git", c("-C", repo_root(), "ls-files", GENERATED_CONFIG),
                     stdout = TRUE)
  testthat::expect_length(tracked, 1L)

  ## it sits beside hand-maintained panels rather than among results
  siblings <- system2("git", c("-C", repo_root(), "ls-files", "config/marker_panels"),
                      stdout = TRUE)
  testthat::expect_gt(length(siblings), 1L)

  ## and it is consumed as configuration, not as a result
  reg <- read_pipeline_registry(repo_path("pipeline.yml"))
  st <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                       include_unsupported = TRUE)
  consumers <- sum(vapply(seq_len(nrow(st)), function(k) {
    any(c(sp(st$consumes_required[k]), sp(st$consumes_optional[k])) == GENERATED_CONFIG)
  }, logical(1)))
  testthat::expect_gt(consumers, 1L)
})

# --- sections 26 and 27: coordination, replayed from preserved evidence --

testthat::test_that("QC has no coordinating family left, and that is the mechanism", {
  own <- utils::read.csv(repo_path("config", "results_ownership.csv"),
                         stringsAsFactors = FALSE)
  co <- own[own$classification == "LEGITIMATE_MULTI_STAGE_COORDINATION", , drop = FALSE]
  qc <- co[grepl("analysis/qc/", co$canonical_owner), , drop = FALSE]
  testthat::expect_identical(nrow(qc), 0L)
  mine <- own[grepl("^analysis/qc/", own$canonical_owner), , drop = FALSE]
  testthat::expect_gt(nrow(mine), 0L)
  for (i in seq_len(nrow(mine))) {
    testthat::expect_identical(mine$classification[i], "SINGLE_OWNER")
  }
})

testthat::test_that("the pre-migration QC family still classifies the same way", {
  ## Section 27: the circumstances that exposed the classifier are gone, so the
  ## evidence is kept and replayed. The QC family was one owner writing
  ## publication-style figures and one contributor writing diagnostics into a
  ## sibling directory: distinct files, distinct parents, no collision.
  f <- repo_path("audits", "phase6g_coordinating_families_qc_premigration.csv")
  testthat::skip_if_not(file.exists(f))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  testthat::expect_gt(nrow(d), 0L)
  testthat::expect_identical(unique(d$sharing_class), "SHARED_DIRECTORY_ONLY")

  ## several writers really did share the family, with exactly one owner
  testthat::expect_gt(length(unique(d$script)), 1L)
  testthat::expect_identical(
    length(unique(d$script[d$role == "CANONICAL_OWNER"])), 1L)
  ## and no artifact had two writers
  testthat::expect_false(any(d$shared_with_another_writer %in% TRUE))

  own <- utils::read.csv(repo_path("config", "results_ownership.csv"),
                         stringsAsFactors = FALSE)
  proposals <- split(d$current_output, d$script)
  proposals <- lapply(proposals, function(p) unique(p[!is.na(p) & nzchar(p)]))
  collisions <- migration_destination_collisions(proposals, own)
  testthat::expect_length(collisions, 0L)
  testthat::expect_length(migration_gate_blockers(collisions), 0L)
})

testthat::test_that("the synthetic shared-artifact fixture still hard-stops", {
  ## Section 27: keep the collision detector testable now that no real
  ## coordinating family remains anywhere in the repository.
  shared <- "results/qc/assess_pca_confounding/global/tables/x.csv"
  collisions <- migration_destination_collisions(
    list(assess_pca_confounding = shared, partition_variance = shared), NULL)
  testthat::expect_length(collisions, 1L)
  testthat::expect_false(collisions[[1]]$resolvable)
  testthat::expect_gt(length(migration_gate_blockers(collisions)), 0L)
})

testthat::test_that("no QC artifact has two writers", {
  s <- registry_steps()
  q <- s[grepl("^analysis/qc/", s$script), , drop = FALSE]
  out <- do.call(rbind, lapply(seq_len(nrow(q)), function(i) {
    p <- sp(q$produces[i])
    p <- p[grepl("[.][A-Za-z0-9]{2,5}$", basename(p))]
    if (!length(p)) return(NULL)
    data.frame(script = q$script[i], out = p, stringsAsFactors = FALSE)
  }))
  testthat::expect_gt(nrow(out), 0L)
  testthat::expect_identical(sum(duplicated(out$out)), 0L)
})

# --- sections 10 and 28: lifecycle and naming ----------------------------

testthat::test_that("no QC output is declared under work/ and none carries chronology", {
  s <- registry_steps()
  banned <- c(LEGACY_NS, "Stage13", "stage13", "01_Tables")
  for (i in seq_len(nrow(s))) {
    if (!grepl("^analysis/qc/", s$script[i])) next
    for (p in sp(s$produces[i])) {
      testthat::expect_false(startsWith(p, "work/"), info = p)
      if (identical(p, GENERATED_CONFIG)) next
      for (b in banned) {
        testthat::expect_false(grepl(b, p, fixed = TRUE), info = paste(p, "contains", b))
      }
      seg <- strsplit(p, "/", fixed = TRUE)[[1]]
      testthat::expect_false(any(grepl("^[0-9]{2}[a-z]?_", seg)), info = p)
    }
  }
})

testthat::test_that("the lifecycle audit covers every declared output", {
  f <- repo_path("audits", "phase6g_qc_lifecycle.csv")
  testthat::skip_if_not(file.exists(f))
  lc <- utils::read.csv(f, stringsAsFactors = FALSE)
  s <- registry_steps()
  q <- s[grepl("^analysis/qc/", s$script), , drop = FALSE]
  declared <- unlist(lapply(q$produces, sp))
  testthat::expect_identical(nrow(lc), length(declared))
  testthat::expect_identical(sum(lc$lifecycle == "WORK_INTERMEDIATE"), 0L)
  testthat::expect_identical(sum(lc$lifecycle == "GENERATED_CONFIG_CONTRACT"), 1L)
})

# --- section 24: the cross-domain readers -------------------------------

testthat::test_that("every cross-domain QC reader loads the resolver it calls", {
  tracked <- system2("git", c("-C", repo_root(), "ls-files"), stdout = TRUE)
  code <- tracked[grepl("^(analysis|R)/.*[.][Rr]$", tracked)]
  code <- code[!grepl("^analysis/qc/", code)]
  callers <- 0L
  for (rel in code) {
    f <- repo_path(rel)
    if (!file.exists(f)) next
    exprs <- tryCatch(parse(f), error = function(e) NULL)
    if (is.null(exprs)) next
    uses <- FALSE
    for (e in exprs) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% c("qc_find", "qc_dir_any")) uses <<- TRUE
    })
    if (!uses) next
    callers <- callers + 1L
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    testthat::expect_true(
      grepl("qc_result_paths.R", txt, fixed = TRUE) ||
        grepl("qc_exploration_utils.R", txt, fixed = TRUE),
      info = paste(rel, "calls a QC resolver without loading it"))
  }
  ## twelve analysis scripts plus two libraries were repointed
  testthat::expect_gte(callers, 12L)
})

testthat::test_that("the archived fidelity artifact is deliberately not repointed", {
  ## Its producer is archive/02_qc/04d_compartment_marker_fidelity.r, so no
  ## normalized copy will ever exist. Its reader guards with
  ## read_csv_if_exists, and giving it a normalized candidate would invent a
  ## destination nothing writes.
  txt <- paste(readLines(repo_path("analysis", "publication_source_data",
                                   "build_biological_claims_table.R"),
                         warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("04d_compartment_marker_fidelity", txt, fixed = TRUE))
  reg <- paste(readLines(repo_path("pipeline.yml"), warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("archive/02_qc/04d_compartment_marker_fidelity.r",
                              reg, fixed = TRUE))
})
