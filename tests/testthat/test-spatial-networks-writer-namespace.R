source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "spatial_network_utils.R"))
source(repo_path("R", "migration_gate_utils.R"))

# Phase 6G.2: the spatial_networks migration.
#
# The enrichment pilot proved a writer's destination by static extraction plus
# a dry run. This domain needs two things that phase did not:
#
#   * runtime-load correctness. Three scripts in this migration called the
#     shared resolver without sourcing it. Every one of them parsed cleanly,
#     because R resolves functions at call time, so a parse check said "ok"
#     while the script would have died on invocation. A parse check is not a
#     load check.
#   * resolver precedence under all four existence combinations, because
#     today's correct answer is still the historical object and will be until a
#     normalized run exists. A migration that silently kept preferring the
#     historical path would look identical to a working one.

DOMAIN <- "spatial_networks"
WRITERS <- c("build_spatial_networks", "build_differential_networks",
             "test_network_stability", "test_differential_network_stability",
             "render_differential_network_figures", "render_network_chord_diagram")
LEGACY_NS <- "07_spatial_networks"

registry_steps <- function() {
  reg <- read_pipeline_registry(repo_path("pipeline.yml"))
  s <- pipeline_steps(reg, pipeline_stage_names(reg), dataset = "all",
                      include_unsupported = TRUE)
  s[!duplicated(s$script), , drop = FALSE]
}
sp <- function(x) {
  p <- trimws(unlist(strsplit(paste(x, collapse = "|"), "[|;]")))
  unique(p[nzchar(p)])
}
active_code <- function(f) {
  l <- readLines(repo_path(f), warn = FALSE)
  paste(l[!grepl("^\\s*#", l)], collapse = "\n")
}

# --- section 3: declared == code-resolved == normalized -------------------

testthat::test_that("all six writers agree: registry, code and normalized path", {
  s <- registry_steps()
  ds <- "neuron_neuropil"
  checked <- 0L

  for (aid in WRITERS) {
    f <- file.path("analysis", DOMAIN, paste0(aid, ".R"))
    row <- s[s$script == f, , drop = FALSE]
    testthat::expect_identical(nrow(row), 1L, info = aid)

    declared <- sp(row$produces)
    testthat::expect_gt(length(declared), 0L)

    # every declared path must be under this analysis's normalized root
    root <- paste0("results/", DOMAIN, "/", aid, "/")
    testthat::expect_true(all(startsWith(declared, root)),
                          info = paste(aid, ":", paste(declared, collapse = ", ")))

    # and the code must build that same root through the resolver
    code <- active_code(f)
    testthat::expect_match(code, "canonical_module_dirs\\(\"spatial_networks\", ANALYSIS_ID",
                           info = aid)
    testthat::expect_match(code, paste0('ANALYSIS_ID <- "', aid, '"'), fixed = TRUE)

    resolved <- relative_to(canonical_result_path(DOMAIN, aid, ds, "tables"))
    testthat::expect_identical(resolved,
                               paste0("results/", DOMAIN, "/", aid, "/", ds, "/tables"))
    checked <- checked + 1L
  }
  testthat::expect_identical(checked, 6L)
})

testthat::test_that("no spatial_networks writer targets the historical namespace", {
  for (aid in WRITERS) {
    code <- active_code(file.path("analysis", DOMAIN, paste0(aid, ".R")))
    testthat::expect_false(grepl(LEGACY_NS, code, fixed = TRUE), info = aid)
    testthat::expect_false(grepl("create_module_dirs", code, fixed = TRUE), info = aid)
    # the helper deliberately offers no processed alias, so a writer cannot
    # avoid choosing models/ or work/ for what used to be data/processed
    testthat::expect_false(grepl("CANONICAL_PATHS$processed", code, fixed = TRUE), info = aid)
  }
})

# --- section 4: parse, then load, then resolve ---------------------------

testthat::test_that("every resolver caller can actually load the resolver", {
  # A textual mention of the library is not proof: the delegation comment
  # inserted during this migration names spatial_network_utils.R, which
  # defeated a substring guard and left three scripts unable to load it.
  # The source() call is found in the parse tree instead.
  callers <- character(0)
  # git ls-files is relative to the working directory, and testthat runs from
  # tests/testthat, so an unanchored call lists only the test files and finds
  # no callers at all - a silently empty search.
  tracked <- system2("git", c("-C", repo_root(), "ls-files"), stdout = TRUE)
  rfiles <- tracked[grepl("^analysis/.*[.][Rr]$", tracked)]
  for (f in rfiles) {
    p <- repo_path(f)
    if (!file.exists(p)) next
    if (!grepl("resolve_spatial_network_object|resolve_bootstrap_differential_tables",
               active_code(f))) next
    callers <- c(callers, f)
  }
  testthat::expect_gt(length(callers), 0L)

  sourced_libs <- function(f) {
    ex <- tryCatch(parse(repo_path(f)), error = function(e) NULL)
    if (is.null(ex)) return(character(0))
    libs <- character(0)
    walk <- function(e) {
      if (!tryCatch({ e; TRUE }, error = function(...) FALSE)) return(invisible(NULL))
      if (is.call(e)) {
        fn <- e[[1]]
        if (is.name(fn) && identical(as.character(fn), "source")) {
          lits <- character(0)
          collect <- function(x) {
            if (!tryCatch({ x; TRUE }, error = function(...) FALSE)) return(invisible(NULL))
            if (is.character(x) && length(x) == 1L) lits <<- c(lits, x)
            if (is.call(x)) for (k in seq_along(x)) collect(tryCatch(x[[k]], error = function(...) NULL))
            invisible(NULL)
          }
          collect(e)
          libs <<- c(libs, lits)
        }
        for (k in seq_along(e)) walk(tryCatch(e[[k]], error = function(...) NULL))
      }
      invisible(NULL)
    }
    for (e in ex) walk(e)
    unique(libs)
  }

  missing <- character(0)
  for (f in callers) {
    testthat::expect_error(parse(repo_path(f)), NA, info = paste("parses:", f))
    libs <- sourced_libs(f)
    if (!any(grepl("spatial_network_utils", libs, fixed = TRUE))) {
      missing <- c(missing, f)
    }
  }
  if (length(missing)) {
    testthat::fail(paste0("resolver called without a source() for its library:\n  ",
                          paste(missing, collapse = "\n  ")))
  }
  testthat::expect_length(missing, 0L)
})

testthat::test_that("the resolver is reachable after loading only what a caller loads", {
  # emulate the caller's own bootstrap, then require the function to exist
  env <- new.env(parent = globalenv())
  sys.source(repo_path("R", "paths.R"), envir = env)
  sys.source(repo_path("R", "spatial_network_utils.R"), envir = env)
  for (fn in c("resolve_spatial_network_object", "spatial_network_object_candidates",
               "resolve_bootstrap_differential_tables",
               "bootstrap_differential_tables_candidates")) {
    testthat::expect_true(exists(fn, envir = env, inherits = FALSE), info = fn)
    testthat::expect_true(is.function(get(fn, envir = env)), info = fn)
  }
})

# --- section 5: resolver precedence, dynamically ------------------------

testthat::test_that("resolver precedence holds in all four existence combinations", {
  root <- withr::local_tempdir("sn_resolve_")
  old <- Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = NA_character_)
  Sys.setenv(PROTEOMICS_PROJECT_ROOT = root)
  on.exit({
    if (is.na(old)) Sys.unsetenv("PROTEOMICS_PROJECT_ROOT") else Sys.setenv(PROTEOMICS_PROJECT_ROOT = old)
  }, add = TRUE)
  testthat::skip_if(!identical(normalizePath(repo_root(), winslash = "/"),
                               normalizePath(root, winslash = "/")),
                    "project root is not overridable in this checkout")

  ds <- "neuron_neuropil"
  unit <- "CA1_SLM"
  cand <- spatial_network_object_candidates(ds, unit)
  mk <- function(p) {
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    writeLines("x", p)
  }

  # neither present -> the canonical destination is named, so a failure points
  # at where the object belongs rather than where it used to be
  testthat::expect_identical(resolve_spatial_network_object(ds, unit),
                             unname(cand[["normalized"]]))

  # legacy only -> legacy selected, which is today's correct behaviour
  mk(cand[["legacy_scoped"]])
  testthat::expect_identical(resolve_spatial_network_object(ds, unit),
                             unname(cand[["legacy_scoped"]]))

  # both present -> normalized wins, or the migration never takes effect
  mk(cand[["normalized"]])
  testthat::expect_identical(resolve_spatial_network_object(ds, unit),
                             unname(cand[["normalized"]]))

  # normalized only
  unlink(cand[["legacy_scoped"]])
  testthat::expect_identical(resolve_spatial_network_object(ds, unit),
                             unname(cand[["normalized"]]))
})

testthat::test_that("the candidate order puts normalized first and legacy last", {
  cand <- spatial_network_object_candidates("microglia", "region")
  testthat::expect_identical(names(cand),
                             c("normalized", "legacy_scoped", "legacy_unscoped"))
  testthat::expect_true(grepl("results/spatial_networks/build_spatial_networks/",
                              cand[["normalized"]], fixed = TRUE))
  testthat::expect_true(grepl("data/processed/07_spatial_networks/",
                              cand[["legacy_scoped"]], fixed = TRUE))

  boot <- bootstrap_differential_tables_candidates("microglia")
  testthat::expect_identical(names(boot),
                             c("normalized", "legacy_flat", "legacy_01_tables"))
  testthat::expect_true(grepl("results/spatial_networks/test_differential_network_stability/",
                              boot[["normalized"]], fixed = TRUE))
})

# --- section 6: the models lifecycle ------------------------------------

testthat::test_that("the network object is a model and never a work intermediate", {
  ds <- "neuron_neuropil"
  cand <- spatial_network_object_candidates(ds, "CA1_SLM")
  norm <- relative_to(cand[["normalized"]])

  testthat::expect_true(grepl("^results/spatial_networks/build_spatial_networks/", norm))
  testthat::expect_true(grepl("/models/", norm, fixed = TRUE))
  testthat::expect_false(grepl("^work/", norm))
  testthat::expect_false(grepl("/work/", norm, fixed = TRUE))

  # the writer must place it in models/, not work/
  code <- active_code("analysis/spatial_networks/build_spatial_networks.R")
  testthat::expect_match(code, "processed = CANONICAL_PATHS$models", fixed = TRUE)
  # and its scratch network files must go to work/
  testthat::expect_match(code, 'file.path(CANONICAL_PATHS$work, "network_files")', fixed = TRUE)
})

testthat::test_that("every consumer of the network object resolves the model path first", {
  # the four in-repo readers all delegate, so one precedence rule covers them
  readers <- c("analysis/spatial_networks/build_differential_networks.R",
               "analysis/spatial_networks/test_network_stability.R",
               "analysis/spatial_networks/test_differential_network_stability.R",
               "analysis/integration/test_network_behaviour_coupling.R",
               "analysis/wgcna/build_module_spatial_networks.R")
  for (f in readers) {
    code <- active_code(f)
    testthat::expect_match(code, "resolve_spatial_network_object(", fixed = TRUE, info = f)
    # no reader may still build the historical path itself
    testthat::expect_false(grepl('path_processed("07_spatial_networks"', code, fixed = TRUE),
                           info = f)
  }
})

# --- section 8: 01_Tables is gone from active destinations --------------

testthat::test_that("no active spatial_networks destination contains 01_Tables", {
  s <- registry_steps()
  sn <- s[startsWith(s$script, paste0("analysis/", DOMAIN, "/")), , drop = FALSE]
  testthat::expect_identical(sum(grepl("01_Tables", sn$produces)), 0L)
  testthat::expect_identical(sum(grepl("01_Tables", sn$consumes_required)), 0L)

  # and no writer constructs it in active code
  for (aid in WRITERS) {
    code <- active_code(file.path("analysis", DOMAIN, paste0(aid, ".R")))
    testthat::expect_false(grepl("01_Tables", code, fixed = TRUE), info = aid)
  }

  # it survives only as the last compatibility candidate, which is deliberate
  boot <- bootstrap_differential_tables_candidates("microglia")
  testthat::expect_true(grepl("01_Tables", boot[["legacy_01_tables"]], fixed = TRUE))
})

# --- section 9: the collision gate, on synthetic cases ------------------

testthat::test_that("the collision gate fails two writers sharing a destination", {
  dest <- "results/spatial_validation/atlas/global/tables/shared.csv"
  proposals <- list(writer_A = dest, writer_B = dest)

  col <- migration_destination_collisions(proposals, ownership = NULL)
  testthat::expect_length(col, 1L)
  testthat::expect_setequal(col[[dest]]$writers, c("writer_A", "writer_B"))
  testthat::expect_false(col[[dest]]$resolvable)

  blockers <- migration_gate_blockers(col)
  testthat::expect_gt(length(blockers), 0L)
  testthat::expect_true(any(grepl("share a destination", blockers)))
})

testthat::test_that("a shared destination is not excused by a declared owner", {
  # the registry names the owner of a family; it never licenses two writers to
  # race for one file
  dest <- "results/spatial_validation/atlas/global/tables/shared.csv"
  own <- data.frame(
    result_family = "results/spatial_validation/atlas/global/tables",
    canonical_owner = "analysis/spatial_validation/writer_A.R",
    stringsAsFactors = FALSE)
  col <- migration_destination_collisions(list(writer_A = dest, writer_B = dest), own)
  testthat::expect_length(col, 1L)
  testthat::expect_true(col[[dest]]$owner_is_one_of_the_writers)
  testthat::expect_false(col[[dest]]$resolvable)
  testthat::expect_gt(length(migration_gate_blockers(col)), 0L)
})

testthat::test_that("the gate passes an owner plus a contributor writing elsewhere", {
  own <- data.frame(
    result_family = "results/spatial_validation/atlas/global/tables",
    canonical_owner = "analysis/spatial_validation/writer_A.R",
    stringsAsFactors = FALSE)
  proposals <- list(
    writer_A = "results/spatial_validation/atlas/global/tables/canonical.csv",
    writer_B = "results/spatial_validation/atlas/global/tables/contributor_intermediate.csv")
  col <- migration_destination_collisions(proposals, own)
  testthat::expect_length(col, 0L)
  testthat::expect_length(migration_gate_blockers(col), 0L)
})

testthat::test_that("the gate blocks on an unclassified dependency or a failed inventory", {
  testthat::expect_true(any(grepl("unknown dependency",
                                  migration_gate_blockers(list(), unknown_dependency_kinds = 3L))))
  testthat::expect_true(any(grepl("could not be produced",
                                  migration_gate_blockers(list(), inventories_failed = 1L))))
  testthat::expect_length(migration_gate_blockers(list(), 0L, 0L), 0L)
})

# --- resolve: the destination call must be reachable, not merely present ---

# The third leg of parse -> load -> resolve. canonical_module_dirs() is scoped
# by dataset, so it cannot run before the dataset is resolved;
# create_module_dirs(MODULE_ID, SUBSTEP_ID) could, and the migration first
# left the call where that one sat. Five of the six writers referenced
# NETWORK_DATASET one line before it was assigned. Each parsed cleanly and
# each would have failed immediately on invocation.
#
# The rule is general: any symbol the call depends on that this file assigns
# at top level must be assigned in an earlier top-level expression. A symbol
# the file never assigns is a function or a global and is not this test's
# business.
testthat::test_that("every writer defines its scope before building destinations", {
  checked <- 0L

  for (aid in WRITERS) {
    f <- repo_path(file.path("analysis", "spatial_networks", paste0(aid, ".R")))
    testthat::expect_true(file.exists(f))
    exprs <- parse(f)

    syms <- function(e) {
      out <- character(0)
      rec <- function(x) {
        if (!tryCatch({ x; TRUE }, error = function(...) FALSE)) return(invisible(NULL))
        if (is.name(x)) {
          nm <- as.character(x)
          if (nzchar(nm)) out <<- c(out, nm)
          return(invisible(NULL))
        }
        if (is.call(x) || is.expression(x) || is.list(x)) {
          for (i in seq_along(x)) {
            el <- tryCatch(x[[i]], error = function(...) NULL)
            skip <- tryCatch(is.null(el) || (is.symbol(el) && !nzchar(as.character(el))),
                             error = function(...) TRUE)
            if (skip) next
            rec(el)
          }
        }
        invisible(NULL)
      }
      rec(e)
      unique(out)
    }

    ## top-level assignment targets, in order of appearance
    assigned_at <- list()
    for (k in seq_along(exprs)) {
      e <- exprs[[k]]
      if (is.call(e) && length(e) == 3L &&
          as.character(e[[1]])[1] %in% c("<-", "=", "<<-") && is.name(e[[2]])) {
        nm <- as.character(e[[2]])
        if (is.null(assigned_at[[nm]])) assigned_at[[nm]] <- k
      }
    }

    ## the expression that builds the destinations
    idx <- which(vapply(seq_along(exprs), function(k) {
      "canonical_module_dirs" %in% syms(exprs[[k]])
    }, logical(1)))
    testthat::expect_length(idx, 1L)

    deps <- setdiff(syms(exprs[[idx]]), c("canonical_module_dirs", "CANONICAL_PATHS"))
    for (d in deps) {
      at <- assigned_at[[d]]
      if (is.null(at)) next          # a function or a global, not assigned here
      testthat::expect_true(at < idx,
        info = paste0(aid, ": '", d, "' is assigned at top-level expression ", at,
                      " but canonical_module_dirs() runs at ", idx,
                      ", so the writer cannot build its destinations"))
      checked <- checked + 1L
    }
  }

  testthat::expect_gt(checked, 0L)
})
