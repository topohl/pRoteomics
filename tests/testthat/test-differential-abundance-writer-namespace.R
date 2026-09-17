source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "differential_abundance_paths.R"))

# Phase 6G.4: the differential_abundance migration.
#
# This domain is the first where the writers use several different idioms for
# the same job: a legacy directory factory, explicit path_results() blocks, a
# do.call over assembled parts, a conditional substep, and one writer that
# supports redirecting its whole output root. So the writer-truth check has to
# be on every writer, not a sample.
#
# It also holds two artifacts that other domains read as state, the
# clusterProfiler and compareGO manifests. Those resolve through two central
# functions, so the whole cross-domain surface turns on those two resolving
# normalized-first and still finding the historical copy today.

WRITERS <- c("annotate_neuropil_reference", "audit_gsea_protein_direction",
             "audit_stress_response_biology", "build_go_program_atlas",
             "build_sus_res_dap_atlas", "compare_external_stress_signatures",
             "compare_go_enrichment", "run_clusterprofiler_enrichment",
             "summarize_biological_programs", "test_microglia_targeted_signatures",
             "validate_control_spatial_identity")
LEGACY_NS <- "04_differential_expression_enrichment"

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
  repo_path(file.path("analysis", "differential_abundance", paste0(aid, ".R")))
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

# --- section 23: writer truth, every writer, no sampling ------------------

testthat::test_that("all eleven writers agree: registry, code and normalized path", {
  s <- registry_steps()
  checked <- 0L
  for (aid in WRITERS) {
    script <- paste0("analysis/differential_abundance/", aid, ".R")
    i <- which(s$script == script)
    testthat::expect_length(i, 1L)

    declared <- sp(s$produces[i[1]])
    testthat::expect_gt(length(declared), 0L)

    root <- paste0("results/differential_abundance/", aid, "/")
    testthat::expect_true(
      all(startsWith(declared, root)),
      info = paste(aid, "declares outside its namespace:",
                   paste(declared[!startsWith(declared, root)], collapse = ", ")))
    testthat::expect_false(any(grepl(LEGACY_NS, declared, fixed = TRUE)),
                           info = paste(aid, "still declares the historical namespace"))
    checked <- checked + 1L
  }
  testthat::expect_identical(checked, length(WRITERS))
})

testthat::test_that("no writer constructs a destination in the historical namespace", {
  ## the parse tree, so a comment or a recorded input path cannot count
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    hit <- FALSE
    for (e in exprs) walk(e, function(x) {
      if (hit || !is.call(x)) return(invisible(NULL))
      nm <- cname(x)
      if (is.na(nm) || !nm %in% c("create_module_dirs", "module_paths")) return(invisible(NULL))
      hit <<- TRUE
    })
    testthat::expect_false(hit, info = paste(aid, "still calls the legacy directory factory"))
  }
})

testthat::test_that("every writer resolves through the central path API", {
  API <- c("differential_abundance_dirs", "differential_abundance_relative_path",
           "canonical_module_dirs", "canonical_result_path")
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    used <- FALSE
    for (e in exprs) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% API) used <<- TRUE
    })
    testthat::expect_true(used, info = paste(aid, "builds no destination through the path API"))
  }
})

# --- section 17: LOAD and REACHABLE --------------------------------------

testthat::test_that("every writer actually sources the domain path helper", {
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    loaded <- FALSE
    for (e in exprs) walk(e, function(x) {
      if (loaded || !is.call(x)) return(invisible(NULL))
      nm <- cname(x)
      if (is.na(nm) || !nm %in% c("source", "sys.source")) return(invisible(NULL))
      p <- tryCatch(eval(x[[2]], envir = globalenv()), error = function(...) NA_character_)
      if (length(p) == 1L && !is.na(p) &&
          identical(basename(p), "differential_abundance_paths.R") && file.exists(p)) {
        loaded <<- TRUE
      }
      invisible(NULL)
    })
    testthat::expect_true(loaded, info = paste(aid, "never sources the domain path helper"))
  }
})

testthat::test_that("each destination call runs after the variables it needs", {
  DEST <- c("differential_abundance_dirs", "canonical_module_dirs",
            "differential_abundance_relative_path")
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
    ## only top-level expressions constrain order; a use inside a function body
    ## is evaluated when that function is called
    idx <- which(vapply(seq_along(exprs),
                        function(k) any(DEST %in% syms(exprs[[k]])), logical(1)))
    for (j in idx) {
      deps <- setdiff(syms(exprs[[j]]), c(DEST, "CANONICAL_PATHS", "PATHS", "paths"))
      for (d in deps) {
        at <- assigned_at[[d]]
        if (is.null(at)) next
        ## at == j means the symbol is assigned by this very expression, which
        ## is the case when the destination call sits inside a function
        ## definition: validate_control_spatial_identity builds its paths in
        ## out(), nested in control_spatial_identity_main. That call runs when
        ## the function is invoked, so it constrains nothing. Only an
        ## assignment strictly *after* the destination expression is a fault.
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

# --- section 18 fallback: isolated path resolution ------------------------

testthat::test_that("the two writers without a usable dry-run still resolve", {
  ## run_clusterprofiler_enrichment stops at a pre-existing MAX_PATH guard and
  ## compare_go_enrichment at a pre-existing bootstrap ordering bug, both of
  ## which reproduce identically before this migration. Their destinations are
  ## therefore checked directly against the registry instead.
  s <- registry_steps()
  for (aid in c("run_clusterprofiler_enrichment", "compare_go_enrichment")) {
    d <- differential_abundance_dirs(aid, scope = "neuron_neuropil", create = FALSE)
    testthat::expect_true(is.list(d))
    testthat::expect_true(all(c("tables", "models", "manifests") %in% names(d)))
    ## the manifest each one declares must sit under the models directory the
    ## helper returns
    declared <- sp(s$produces[s$script == paste0("analysis/differential_abundance/", aid, ".R")])
    man <- declared[grepl("_manifest[.]csv$", declared)]
    testthat::expect_gte(length(man), 1L)
    for (m in man) {
      expected <- file.path("results", "differential_abundance", aid,
                            "<dataset>", "models", basename(m))
      testthat::expect_identical(m, expected)
    }
    ## and the helper's models directory matches once the dataset is concrete
    testthat::expect_true(grepl(
      file.path("results", "differential_abundance", aid, "neuron_neuropil", "models"),
      d$models, fixed = TRUE))
  }
})

testthat::test_that("the dropped processed alias is gone and unused", {
  ## canonical_module_dirs() deliberately offers no `processed` member, so a
  ## writer that used to write into data/processed has to choose models/ or
  ## work/. Leaving a $processed reference behind yields NULL and fails at
  ## runtime, which is how it was caught.
  d <- differential_abundance_dirs("run_clusterprofiler_enrichment",
                                   scope = "neuron_neuropil", create = FALSE)
  testthat::expect_null(d$processed)
  for (aid in WRITERS) {
    txt <- paste(readLines(writer_file(aid), warn = FALSE), collapse = "\n")
    testthat::expect_false(grepl("CANONICAL_PATHS$processed", txt, fixed = TRUE),
                           info = aid)
  }
})

# --- section 22: normalized-first for the cross-domain state -------------

testthat::test_that("the shared manifests resolve normalized-first", {
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("da_", Sys.getpid()))
  ds <- "neuron_neuropil"
  norm_dir <- file.path(tmp, "results", "differential_abundance",
                        "run_clusterprofiler_enrichment", ds, "models")
  leg_dir <- file.path(tmp, "data", "processed", LEGACY_NS, "clusterProfiler", ds)
  dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(leg_dir, recursive = TRUE, showWarnings = FALSE)
  nf <- file.path(norm_dir, "clusterProfiler_manifest.csv")
  lf <- file.path(leg_dir, "clusterProfiler_manifest.csv")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  combos <- list(
    list(norm = FALSE, leg = TRUE,  want = "legacy"),
    list(norm = TRUE,  leg = TRUE,  want = "normalized"),
    list(norm = TRUE,  leg = FALSE, want = "normalized"),
    list(norm = FALSE, leg = FALSE, want = "normalized"))

  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    for (cb in combos) {
      unlink(c(nf, lf))
      if (cb$norm) writeLines("x", nf)
      if (cb$leg) writeLines("x", lf)
      got <- resolve_differential_abundance_state(
        "clusterProfiler_manifest.csv", "run_clusterprofiler_enrichment", ds,
        "clusterProfiler")
      is_norm <- grepl("results/differential_abundance/", got, fixed = TRUE)
      if (identical(cb$want, "normalized")) {
        testthat::expect_true(is_norm,
          info = paste("normalized =", cb$norm, "legacy =", cb$leg, "->", got))
      } else {
        testthat::expect_false(is_norm,
          info = paste("normalized =", cb$norm, "legacy =", cb$leg, "->", got))
      }
    }
  })
})

testthat::test_that("an empty normalized directory does not shadow the historical manifest", {
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("dashadow_", Sys.getpid()))
  ds <- "neuron_neuropil"
  norm_dir <- file.path(tmp, "results", "differential_abundance",
                        "compare_go_enrichment", ds, "models")
  leg_dir <- file.path(tmp, "data", "processed", LEGACY_NS, "compareGO", ds)
  dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(leg_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("x", file.path(leg_dir, "compareGO_input_manifest.csv"))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    got <- resolve_differential_abundance_state(
      "compareGO_input_manifest.csv", "compare_go_enrichment", ds, "compareGO")
    testthat::expect_true(grepl("data/processed", got, fixed = TRUE))
  })
})

testthat::test_that("today both shared manifests still resolve to the historical copy", {
  source(repo_path("R", "enrichment_io.R"))
  for (ds in c("neuron_neuropil", "neuron_soma", "microglia")) {
    p <- canonical_clusterprofiler_manifest_path(ds)
    q <- canonical_comparego_manifest_path(ds)
    ## nothing has been rerun, so the historical copy is what exists
    if (file.exists(p)) {
      testthat::expect_true(grepl("data/processed", p, fixed = TRUE), info = ds)
    }
    if (file.exists(q)) {
      testthat::expect_true(grepl("data/processed", q, fixed = TRUE), info = ds)
    }
  }
})

# --- sections 10, 12 and 29: naming ---------------------------------------

testthat::test_that("normalized paths carry no chronology and no misleading identifier", {
  banned <- c(LEGACY_NS, "01b_", "Stage13", "stage13", "01_Tables")
  misleading <- c("hotspot", "ca2_slm_hotspot", "hotspot_DAP")
  s <- registry_steps()
  for (i in seq_len(nrow(s))) {
    if (!grepl("^analysis/differential_abundance/", s$script[i])) next
    for (p in sp(s$produces[i])) {
      for (b in banned) {
        testthat::expect_false(grepl(b, p, fixed = TRUE), info = paste(p, "contains", b))
      }
      for (m in misleading) {
        testthat::expect_false(grepl(m, tolower(p), fixed = TRUE),
                               info = paste(p, "reintroduces", m))
      }
      seg <- strsplit(p, "/", fixed = TRUE)[[1]]
      testthat::expect_false(any(grepl("^[0-9]{2}[a-z]?_", seg)), info = p)
    }
  }
})

testthat::test_that("the DAP hierarchy and contrast vocabulary are untouched", {
  ## section 11 and 12: this phase changes addresses, not schemas. The
  ## normalized paths must still name the same objects.
  s <- registry_steps()
  da <- s[grepl("^analysis/differential_abundance/", s$script), , drop = FALSE]
  all_out <- unlist(lapply(da$produces, sp))
  ## the SUS-RES atlas and theme summary keep their names
  testthat::expect_true(any(grepl("sus_res_manuscript_theme_summary.csv", all_out, fixed = TRUE)))
  testthat::expect_true(any(grepl("build_sus_res_dap_atlas/", all_out, fixed = TRUE)))
  ## no contrast table was split or merged: the count of declared outputs is
  ## the same as the historical declaration count recorded in the lifecycle audit
  f <- repo_path("audits", "phase6g_differential_abundance_lifecycle.csv")
  testthat::skip_if_not(file.exists(f))
  lc <- utils::read.csv(f, stringsAsFactors = FALSE)
  testthat::expect_identical(length(all_out), nrow(lc))
})

# --- sections 6 and 35: coordination taxonomy ----------------------------

testthat::test_that("no differential_abundance artifact has two writers", {
  s <- registry_steps()
  da <- s[grepl("^analysis/differential_abundance/", s$script), , drop = FALSE]
  out <- do.call(rbind, lapply(seq_len(nrow(da)), function(i) {
    p <- sp(da$produces[i])
    p <- p[grepl("[.][A-Za-z0-9]{2,5}$", basename(p))]
    if (!length(p)) return(NULL)
    data.frame(script = da$script[i], out = p, stringsAsFactors = FALSE)
  }))
  testthat::expect_gt(nrow(out), 0L)
  testthat::expect_identical(sum(duplicated(out$out)), 0L)
})

testthat::test_that("the enrichment-branch sandbox keeps its historical namespace", {
  ## It is not canonical output, nothing declares or consumes it, and the
  ## publication journal-scope rules exclude it by that name. Repointing it
  ## would mean editing those rules, which belongs to another domain.
  txt <- paste(readLines(writer_file("run_clusterprofiler_enrichment"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl("04_differential_expression_enrichment_comparison",
                              txt, fixed = TRUE))
  helpers <- paste(readLines(repo_path("R", "export_helpers.R"), warn = FALSE),
                   collapse = "\n")
  testthat::expect_true(grepl("04_differential_expression_enrichment_comparison",
                              helpers, fixed = TRUE))
})
