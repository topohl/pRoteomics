source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "spatial_systems_paths.R"))
source(repo_path("R", "migration_gate_utils.R"))

# Phase 6G.3: the spatial_validation migration.
#
# This domain is the first with real coordinating result families, and the
# first where the writers read each other: summarize_spatial_atlas reads five
# families, build_module_spatial_atlas three. So beyond the checks the earlier
# domains needed, this file has to prove two things they did not:
#
#   * that a family shared by several analyses is one owner plus contributors
#     writing distinct artifacts, not two writers racing for one file, and
#     that the gate can tell those apart on real registry data;
#   * that an intra-domain read prefers the normalized copy but still finds
#     the historical one today, because nothing has been rerun and the
#     historical tree is the only place the objects exist.
#
# The second is the one that would fail silently. A reader pinned to the
# historical path keeps working, so the migration would look complete while
# changing nothing at all.

WRITERS <- c("annotate_module_celltypes", "audit_ca2_slm_robustness",
             "audit_stress_identity_robustness", "build_animal_spatial_networks",
             "build_module_spatial_atlas", "build_protein_spatial_atlas",
             "build_spatial_data_contract", "decompose_bilateral_variance",
             "quantify_bilateral_spatial_identity", "quantify_empirical_compartments",
             "quantify_module_bilateral_identity", "quantify_neuropil_detection_context",
             "quantify_neuropil_precision", "summarize_ca2_slm_robustness",
             "summarize_spatial_atlas", "test_network_group_organization",
             "validate_network_workbook", "validate_spatial_foundations")
LEGACY_NS <- "11_spatial_systems"
REPORTS <- c("spatial_systems_atlas.xlsx", "spatial_systems_networks.xlsx",
             "CA2_SLM_robustness_audit.xlsx")

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
  repo_path(file.path("analysis", "spatial_validation", paste0(aid, ".R")))
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

# --- section 17: registry == code == normalized, for every writer ----------

testthat::test_that("all eighteen writers agree: registry, code and normalized path", {
  s <- registry_steps()
  checked <- 0L
  for (aid in WRITERS) {
    script <- paste0("analysis/spatial_validation/", aid, ".R")
    i <- which(s$script == script)
    testthat::expect_length(i, 1L)

    declared <- sp(s$produces[i[1]])
    testthat::expect_gt(length(declared), 0L)

    ## every declared output is in this analysis's own normalized namespace
    expected_root <- paste0("results/spatial_validation/", aid, "/")
    testthat::expect_true(
      all(startsWith(declared, expected_root)),
      info = paste(aid, "declares outputs outside its namespace:",
                   paste(setdiff(declared[!startsWith(declared, expected_root)],
                                 character(0)), collapse = ", ")))

    ## and none of them names the historical stage namespace
    testthat::expect_false(any(grepl(LEGACY_NS, declared, fixed = TRUE)),
                           info = paste(aid, "still declares the historical namespace"))
    checked <- checked + 1L
  }
  testthat::expect_identical(checked, length(WRITERS))
})

testthat::test_that("no spatial_validation writer constructs a legacy destination", {
  ## the parse tree, not the text: a comment naming the old namespace is fine
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    legacy_write <- FALSE
    for (e in exprs) walk(e, function(x) {
      if (legacy_write || !is.call(x)) return(invisible(NULL))
      nm <- cname(x)
      if (is.na(nm) || !nm %in% c("path_results", "path_processed")) return(invisible(NULL))
      lits <- unlist(lapply(as.list(x)[-1], function(a) if (is.character(a)) a else NULL))
      ## a read of the historical tree is permitted; what must not exist is a
      ## destination built there. dir_create/ file.path around it is the tell,
      ## and those are covered by the writer audit; here we assert the domain's
      ## own stage namespace appears in no path construction at all.
      if (any(lits == LEGACY_NS)) legacy_write <<- TRUE
      invisible(NULL)
    })
    testthat::expect_false(
      legacy_write,
      info = paste(aid, "constructs a path in the historical", LEGACY_NS, "namespace"))
  }
})

# --- section 11: the resolver must actually be loaded ----------------------

testthat::test_that("every caller of the shared resolver actually sources it", {
  ## A textual mention is not proof: R resolves functions at call time, so a
  ## script can call a helper it never loaded and still parse cleanly. The
  ## check is on the source() expression, evaluated the way R will evaluate it.
  callers <- 0L
  for (aid in WRITERS) {
    f <- writer_file(aid)
    exprs <- parse(f)

    uses <- FALSE
    for (e in exprs) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% c("spatial_systems_dirs", "spatial_systems_find",
                                  "resolve_spatial_systems_artifact",
                                  "spatial_systems_artifact_candidates",
                                  "resolve_spatial_systems_dir")) uses <<- TRUE
    })
    if (!uses) next
    callers <- callers + 1L

    loaded <- FALSE
    for (e in exprs) walk(e, function(x) {
      if (loaded || !is.call(x)) return(invisible(NULL))
      nm <- cname(x)
      if (is.na(nm) || !nm %in% c("source", "sys.source")) return(invisible(NULL))
      p <- tryCatch(eval(x[[2]], envir = globalenv()), error = function(...) NA_character_)
      if (length(p) == 1L && !is.na(p) &&
          identical(basename(p), "spatial_systems_paths.R") && file.exists(p)) {
        loaded <<- TRUE
      }
      invisible(NULL)
    })
    testthat::expect_true(
      loaded,
      info = paste(aid, "calls the spatial_systems resolver without sourcing",
                   "R/spatial/spatial_systems_paths.R"))
  }
  testthat::expect_identical(callers, length(WRITERS))
})

testthat::test_that("the resolver is reachable after loading only what a caller loads", {
  env <- new.env(parent = globalenv())
  sys.source(repo_path("R", "paths.R"), envir = env)
  sys.source(repo_path("R", "spatial_systems_paths.R"), envir = env)
  for (fn in c("spatial_systems_dirs", "spatial_systems_find",
               "resolve_spatial_systems_artifact",
               "spatial_systems_artifact_candidates")) {
    testthat::expect_true(exists(fn, envir = env, inherits = FALSE), info = fn)
    testthat::expect_true(is.function(get(fn, envir = env)))
  }
})

# --- section 12: the destination call must be reachable --------------------

testthat::test_that("every writer binds what its destination call needs first", {
  ## Phase 6G.2 shipped five writers that built their destinations one line
  ## before the variable they were scoped by existed. Each parsed; each died on
  ## invocation. Parsing proves nothing about order.
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
    idx <- which(vapply(seq_along(exprs), function(k) {
      "spatial_systems_dirs" %in% syms(exprs[[k]])
    }, logical(1)))
    testthat::expect_gte(length(idx), 1L)

    for (j in idx) {
      deps <- setdiff(syms(exprs[[j]]), c("spatial_systems_dirs", "CANONICAL_PATHS"))
      for (d in deps) {
        at <- assigned_at[[d]]
        if (is.null(at)) next     # a function or a global, not assigned here
        testthat::expect_true(
          at < j,
          info = paste0(aid, ": '", d, "' is assigned at expression ", at,
                        " but the destination call runs at ", j))
        checked <- checked + 1L
      }
    }
  }
  testthat::expect_gt(checked, 0L)
})

# --- section 16: normalized-first, all four existence combinations ---------

testthat::test_that("resolution prefers normalized and falls back to legacy", {
  ## repo_root() reads PROTEOMICS_PROJECT_ROOT on every call, so the real
  ## resolver can be pointed at a temporary tree and exercised end to end
  ## rather than re-implemented in the test.
  testthat::skip_if_not_installed("withr")

  tmp <- file.path(tempdir(), paste0("sv_", as.integer(Sys.time()), "_", Sys.getpid()))
  owner <- "build_protein_spatial_atlas"
  fname <- "protein_spatial_cell_affinity.csv"
  norm_dir <- file.path(tmp, "results", "spatial_validation", owner, "global", "tables")
  leg_dir <- file.path(tmp, "results", "tables", "11_spatial_systems", "atlas")
  dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(leg_dir, recursive = TRUE, showWarnings = FALSE)
  nf <- file.path(norm_dir, fname)
  lf <- file.path(leg_dir, fname)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  ## The fourth case is the documented contract of the filename-keyed
  ## resolver: with no owner it cannot name a normalized destination, so it
  ## names the historical one. Every call site guards with file.exists(), so
  ## this only shapes a diagnostic. The owner-keyed resolver is asserted
  ## separately below and does point forward.
  combos <- list(
    list(norm = FALSE, leg = TRUE,  want = "legacy"),
    list(norm = TRUE,  leg = TRUE,  want = "normalized"),
    list(norm = TRUE,  leg = FALSE, want = "normalized"),
    list(norm = FALSE, leg = FALSE, want = "legacy"))

  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    for (cb in combos) {
      unlink(c(nf, lf))
      if (cb$norm) writeLines("x", nf)
      if (cb$leg) writeLines("x", lf)

      got <- spatial_systems_find(fname, "atlas")
      is_norm <- grepl("results/spatial_validation/", got, fixed = TRUE)
      is_leg <- grepl("11_spatial_systems", got, fixed = TRUE)

      if (identical(cb$want, "normalized")) {
        testthat::expect_true(is_norm,
          info = paste("normalized =", cb$norm, "legacy =", cb$leg,
                       "resolved to", got))
        testthat::expect_false(is_leg)
      } else {
        testthat::expect_true(is_leg,
          info = paste("normalized =", cb$norm, "legacy =", cb$leg,
                       "resolved to", got))
      }
    }

    ## the owner-keyed resolver does point at the canonical destination when
    ## nothing exists, so a failure there names where the object belongs
    unlink(c(nf, lf))
    fwd <- resolve_spatial_systems_artifact(fname, owner, "atlas")
    testthat::expect_true(grepl("results/spatial_validation/", fwd, fixed = TRUE))
    testthat::expect_false(grepl("11_spatial_systems", fwd, fixed = TRUE))

    ## and it prefers the normalized copy when both exist
    writeLines("x", nf); writeLines("x", lf)
    both <- resolve_spatial_systems_artifact(fname, owner, "atlas")
    testthat::expect_true(grepl("results/spatial_validation/", both, fixed = TRUE))
  })
})

testthat::test_that("an empty normalized directory does not shadow real data", {
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("svshadow_", Sys.getpid()))
  owner <- "build_protein_spatial_atlas"
  fname <- "protein_spatial_cell_affinity.csv"
  norm_dir <- file.path(tmp, "results", "spatial_validation", owner, "global", "tables")
  leg_dir <- file.path(tmp, "results", "tables", "11_spatial_systems", "atlas")
  dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(leg_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines("x", file.path(leg_dir, fname))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  ## the normalized directory exists but is empty, which is exactly what a
  ## dry run leaves behind
  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    got <- spatial_systems_find(fname, "atlas")
    testthat::expect_true(grepl("11_spatial_systems", got, fixed = TRUE))
  })
})

testthat::test_that("two normalized producers of one filename is an error", {
  ## the one-canonical-writer invariant, enforced where it is relied upon
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("svdup_", Sys.getpid()))
  fname <- "protein_spatial_cell_affinity.csv"
  for (o in c("build_protein_spatial_atlas", "build_module_spatial_atlas")) {
    d <- file.path(tmp, "results", "spatial_validation", o, "global", "tables")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    writeLines("x", file.path(d, fname))
  }
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    testthat::expect_error(spatial_systems_find(fname, "atlas"),
                           "One canonical writer per artifact")
  })
})

# --- sections 8 and 25: lifecycle from evidence, not extension -------------

testthat::test_that("the persistent network object is declared under models", {
  s <- registry_steps()
  i <- which(s$script == "analysis/spatial_validation/build_animal_spatial_networks.R")
  declared <- sp(s$produces[i[1]])
  rds <- declared[grepl("animal_network_objects[.]rds$", declared)]
  testthat::expect_length(rds, 1L)
  ## it is read by test_network_group_organization, so it is a model and not a
  ## work intermediate: work/ may never be named by a downstream contract
  testthat::expect_true(grepl("/global/models/", rds, fixed = TRUE), info = rds)
  testthat::expect_false(startsWith(rds, "work/"))
})

testthat::test_that("the summary workbooks are declared under reports", {
  s <- registry_steps()
  found <- 0L
  for (i in seq_len(nrow(s))) {
    if (!grepl("^analysis/spatial_validation/", s$script[i])) next
    for (p in sp(s$produces[i])) {
      if (!basename(p) %in% REPORTS) next
      testthat::expect_true(grepl("/global/reports/", p, fixed = TRUE), info = p)
      found <- found + 1L
    }
  }
  testthat::expect_identical(found, length(REPORTS))
})

testthat::test_that("no declared spatial_validation output lands in work/", {
  s <- registry_steps()
  for (i in seq_len(nrow(s))) {
    if (!grepl("^analysis/spatial_validation/", s$script[i])) next
    for (p in sp(s$produces[i])) {
      testthat::expect_false(startsWith(p, "work/"), info = p)
    }
  }
})

# --- section 26: no historical layout chronology in the new paths ----------

testthat::test_that("normalized paths carry no chronology tokens", {
  banned <- c("01_Tables", "02_Figures", "Stage13", "11_spatial_systems")
  s <- registry_steps()
  for (i in seq_len(nrow(s))) {
    if (!grepl("^analysis/spatial_validation/", s$script[i])) next
    for (p in sp(s$produces[i])) {
      for (b in banned) {
        testthat::expect_false(grepl(b, p, fixed = TRUE),
                               info = paste(p, "contains", b))
      }
      ## and no bare numeric stage prefix segment
      seg <- strsplit(p, "/", fixed = TRUE)[[1]]
      testthat::expect_false(any(grepl("^[0-9]{2}[a-z]?_", seg)), info = p)
    }
  }
})

# --- sections 4, 6 and 18: the real coordinating families -----------------

testthat::test_that("migration left this domain with no coordinating family", {
  ## Before the migration four families here were written by several analyses
  ## each: atlas, bilateral, ca2_slm_robustness and networks. None survives,
  ## and that is the mechanism working rather than a gap in the registry. A
  ## normalized destination contains its producing analysis, so two analyses
  ## cannot write into one, and each former contributor now owns its own tree.
  ## The coordination relationship did not vanish, it moved into
  ## config/results_ownership.csv, which is the registry of record for it.
  own <- utils::read.csv(repo_path("config", "results_ownership.csv"),
                         stringsAsFactors = FALSE)
  co <- own[own$classification == "LEGITIMATE_MULTI_STAGE_COORDINATION", , drop = FALSE]
  sv <- co[grepl("spatial_validation", co$canonical_owner), , drop = FALSE]
  testthat::expect_identical(nrow(sv), 0L)

  ## every spatial_validation family that remains is single-owner
  mine <- own[grepl("^analysis/spatial_validation/", own$canonical_owner), , drop = FALSE]
  testthat::expect_gt(nrow(mine), 0L)
  for (i in seq_len(nrow(mine))) {
    testthat::expect_length(sp(mine$canonical_owner[i]), 1L)
    testthat::expect_identical(mine$classification[i], "SINGLE_OWNER")
  }
})

testthat::test_that("the gate finds no collision in the real pre-migration families", {
  ## Section 18 asks for the gate to be exercised on actual repository data,
  ## not only synthetically as in Phase 6G.2. The four coordinating families
  ## no longer exist after migration, so the evidence is preserved: the
  ## pre-migration owner/contributor sets and their declared destinations are
  ## recorded in the _premigration audit, and the gate is run against them
  ## here. This is the case the gate has to get right - several writers in one
  ## family, none of them sharing a destination - and it is real data.
  f <- repo_path("audits",
                 "phase6g_coordinating_families_spatial_validation_premigration.csv")
  testthat::skip_if_not(file.exists(f))
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  own <- utils::read.csv(repo_path("config", "results_ownership.csv"),
                         stringsAsFactors = FALSE)

  fams <- unique(d$result_family)
  testthat::expect_gte(length(fams), 4L)

  for (fam in fams) {
    rows <- d[d$result_family == fam, , drop = FALSE]
    ## several writers really did share this family
    testthat::expect_gt(length(unique(rows$script)), 1L)
    ## exactly one of them was the canonical owner
    testthat::expect_identical(
      length(unique(rows$script[rows$role == "CANONICAL_OWNER"])), 1L)

    proposals <- split(rows$current_output, rows$script)
    proposals <- lapply(proposals, function(p) unique(p[!is.na(p) & nzchar(p)]))
    collisions <- migration_destination_collisions(proposals, own)
    testthat::expect_length(collisions, 0L)
    testthat::expect_length(migration_gate_blockers(collisions), 0L)

    ## and no contributor wrote an artifact another one also wrote
    testthat::expect_false(any(rows$shared_with_another_writer %in% TRUE))
  }
})

testthat::test_that("no declared artifact in this domain has two writers", {
  s <- registry_steps()
  sv <- s[grepl("^analysis/spatial_validation/", s$script), , drop = FALSE]
  all_out <- do.call(rbind, lapply(seq_len(nrow(sv)), function(i) {
    p <- sp(sv$produces[i])
    p <- p[grepl("[.][A-Za-z0-9]{2,5}$", basename(p))]
    if (!length(p)) return(NULL)
    data.frame(script = sv$script[i], out = p, stringsAsFactors = FALSE)
  }))
  testthat::expect_gt(nrow(all_out), 0L)
  testthat::expect_identical(sum(duplicated(all_out$out)), 0L)
  ## and filenames are unique too, which is what makes the filename-keyed
  ## resolver sound
  testthat::expect_identical(sum(duplicated(basename(all_out$out))), 0L)
})

testthat::test_that("the gate still fails a genuine shared destination", {
  ## a gate that only ever passes is not evidence, so the same rule is shown
  ## rejecting the shape it exists to reject
  shared <- "results/spatial_validation/build_protein_spatial_atlas/global/tables/x.csv"
  collisions <- migration_destination_collisions(
    list(build_protein_spatial_atlas = shared, build_module_spatial_atlas = shared),
    NULL)
  testthat::expect_length(collisions, 1L)
  testthat::expect_false(collisions[[1]]$resolvable)
  testthat::expect_gt(length(migration_gate_blockers(collisions)), 0L)
})

# --- the legacy comparison oracle must stay pinned to the legacy tree -----

testthat::test_that("the legacy network comparison still addresses the old tree", {
  ## validate_network_workbook compares the animal-level pipeline against the
  ## historical 07_spatial_networks one and asserts that tree is unmodified.
  ## Resolving those reads normalized-first would destroy the comparison, so
  ## this read is deliberately not migrated.
  f <- writer_file("validate_network_workbook")
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl('LEG <- function(...) path_results("tables", "07_spatial_networks", ...)',
                              txt, fixed = TRUE))
})

# --- every legacy fallback family named at a call site is a real family ----

testthat::test_that("each resolver call site names a real filename and family", {
  s <- registry_steps()
  sv <- s[grepl("^analysis/spatial_validation/", s$script), , drop = FALSE]
  declared <- do.call(rbind, lapply(seq_len(nrow(sv)), function(i) {
    p <- sp(sv$produces[i])
    p <- p[grepl("[.][A-Za-z0-9]{2,5}$", basename(p))]
    if (!length(p)) return(NULL)
    data.frame(file = basename(p), stringsAsFactors = FALSE)
  }))
  known_files <- unique(declared$file)
  known_families <- c("atlas", "bilateral", "ca2_slm_robustness", "networks",
                      "precision", "celltype_annotation", "data_contract")

  sites <- 0L
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    for (e in exprs) walk(e, function(x) {
      if (!is.call(x) || !identical(cname(x), "spatial_systems_find")) return(invisible(NULL))
      a <- as.list(x)[-1]
      if (length(a) < 1L || !is.character(a[[1]])) return(invisible(NULL))
      fn <- a[[1]]
      testthat::expect_true(fn %in% known_files,
        info = paste(aid, "resolves", fn, "which no spatial_validation analysis declares"))
      if (length(a) >= 2L && is.character(a[[2]])) {
        testthat::expect_true(a[[2]] %in% known_families,
          info = paste(aid, "names legacy family", a[[2]]))
      }
      sites <<- sites + 1L
      invisible(NULL)
    })
  }
  testthat::expect_gt(sites, 0L)
})
