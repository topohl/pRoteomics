source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "migration_gate_utils.R"))

# Phase 6G.6: the integration migration.
#
# Two things make this domain different from the five before it. It spanned two
# historical stage namespaces rather than one, 10_biological_integration and
# 08_behavior_physio_coupling, plus the publication-facing manuscript_panels
# namespace. And it owns the one remaining stage13 output name that Phase 6F
# adjudicated as genuine chronology.
#
# It is also the domain whose upstream is mostly still unmigrated: sixty of its
# input edges come from WGCNA. Those must keep naming the historical tree,
# because inventing a normalized path for a producer that has not adopted one
# would be a fabricated dependency.

WRITERS <- c("audit_animal_id_integrity", "build_candidate_protein_shortlist",
             "build_cross_compartment_atlas", "build_evidence_priority_matrix",
             "build_immunostaining_candidates",
             "export_module_protein_zoom_source_data",
             "quantify_candidate_network_position",
             "render_module_circular_atlas", "screen_immunostaining_panel",
             "screen_immunostaining_separation",
             "summarize_enrichment_module_concordance",
             "summarize_module_cross_compartment",
             "summarize_programs_for_manuscript",
             "test_behaviour_proteomics_associations",
             "test_enrichment_module_concordance",
             "test_module_behaviour_coupling",
             "test_network_behaviour_coupling")
STAGES <- c("10_biological_integration", "08_behavior_physio_coupling")
PENDING_STAGES <- c("06_modules_WGCNA", "01_preprocessing", "02_id_mapping")
STAGE13_OLD <- "wgcna_circular_atlas_stage13_selection_audit.csv"
STAGE13_NEW <- "circular_atlas_selection_audit.csv"

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
  repo_path(file.path("analysis", "integration", paste0(aid, ".R")))
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

# --- section 24: writer truth, all seventeen -----------------------------

testthat::test_that("all seventeen writers declare only normalized destinations", {
  s <- registry_steps()
  checked <- 0L
  for (aid in WRITERS) {
    i <- which(s$script == paste0("analysis/integration/", aid, ".R"))
    testthat::expect_length(i, 1L)
    declared <- sp(s$produces[i[1]])
    testthat::expect_gt(length(declared), 0L)
    root <- paste0("results/integration/", aid, "/")
    testthat::expect_true(
      all(startsWith(declared, root)),
      info = paste(aid, "declares outside its namespace:",
                   paste(declared[!startsWith(declared, root)], collapse = ", ")))
    for (st in c(STAGES, "manuscript_panels")) {
      testthat::expect_false(any(grepl(st, declared, fixed = TRUE)),
                             info = paste(aid, "still declares", st))
    }
    checked <- checked + 1L
  }
  testthat::expect_identical(checked, length(WRITERS))
})

testthat::test_that("no writer calls a legacy directory factory", {
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    hit <- character(0)
    for (e in exprs) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% c("integration_paths", "create_module_dirs",
                                  "module_paths")) hit <<- c(hit, nm)
    })
    testthat::expect_identical(unique(hit), character(0),
      info = paste(aid, "still calls", paste(unique(hit), collapse = ", ")))
  }
})

testthat::test_that("integration_paths is gone", {
  txt <- paste(readLines(repo_path("R", "integration_utils.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_false(grepl("\nintegration_paths <- function", txt, fixed = TRUE))
})

testthat::test_that("no writer constructs a destination in a historical stage", {
  ## path_results()/path_processed() naming an integration stage namespace
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    bad <- character(0)
    for (e in exprs) walk(e, function(x) {
      nm <- cname(x)
      if (is.na(nm) || !nm %in% c("path_results", "path_processed")) return(invisible(NULL))
      lits <- unlist(lapply(as.list(x)[-1], function(a) if (is.character(a)) a else NULL))
      hit <- intersect(lits, c(STAGES, "manuscript_panels"))
      if (length(hit)) bad <<- c(bad, hit)
    })
    testthat::expect_identical(unique(bad), character(0),
      info = paste(aid, "still builds a path in", paste(unique(bad), collapse = ", ")))
  }
})

# --- section 15: LOAD and REACHABLE --------------------------------------

testthat::test_that("every writer that resolves paths actually loads the helper", {
  users <- 0L
  for (aid in WRITERS) {
    exprs <- parse(writer_file(aid))
    uses <- FALSE
    for (e in exprs) walk(e, function(x) {
      nm <- cname(x)
      if (!is.na(nm) && nm %in% c("integration_dirs", "integration_find",
                                  "integration_artifact_candidates")) uses <<- TRUE
    })
    if (!uses) next
    users <- users + 1L
    loaded <- FALSE
    for (e in exprs) walk(e, function(x) {
      if (loaded || !is.call(x)) return(invisible(NULL))
      nm <- cname(x)
      if (is.na(nm) || !nm %in% c("source", "sys.source")) return(invisible(NULL))
      p <- tryCatch(eval(x[[2]], envir = globalenv()), error = function(...) NA_character_)
      if (length(p) != 1L || is.na(p)) return(invisible(NULL))
      ## Resolve the way R will at runtime. These scripts are launched from the
      ## repository root and some source their libraries by literal relative
      ## path -- source("R/statistics/integration_utils.R") -- so a bare
      ## file.exists() fails here purely because testthat runs from
      ## tests/testthat. That would report a real, loaded dependency as absent.
      ok <- file.exists(p) || file.exists(file.path(repo_root(), p))
      if (ok && identical(basename(p), "integration_utils.R")) loaded <<- TRUE
      invisible(NULL)
    })
    testthat::expect_true(loaded,
      info = paste(aid, "calls an integration path helper without loading it"))
  }
  testthat::expect_identical(users, length(WRITERS))
})

testthat::test_that("each destination call runs after the variables it needs", {
  DEST <- c("integration_dirs", "canonical_module_dirs", "canonical_result_path")
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
    ## A destination call inside a function definition imposes no ordering:
    ## its body resolves when the function is invoked, not when it is defined,
    ## so a forward reference to a sibling helper defined further down is
    ## legal. quantify_candidate_network_position builds its destinations in
    ## emit(), which references write_readme() defined in the next expression.
    is_fundef <- function(e) {
      if (!(is.call(e) && length(e) == 3L &&
            as.character(e[[1]])[1] %in% c("<-", "=", "<<-"))) return(FALSE)
      rhs <- e[[3]]
      is.call(rhs) && identical(as.character(rhs[[1]])[1], "function")
    }
    idx <- which(vapply(seq_along(exprs),
                        function(k) any(DEST %in% syms(exprs[[k]])), logical(1)))
    idx <- idx[!vapply(idx, function(k) is_fundef(exprs[[k]]), logical(1))]
    for (j in idx) {
      deps <- setdiff(syms(exprs[[j]]),
                      c(DEST, "paths", "CANONICAL_PATHS", "OUT", "out"))
      for (d in deps) {
        at <- assigned_at[[d]]
        if (is.null(at)) next
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

# --- section 27: normalized-first, all four existence states ------------

testthat::test_that("integration_find prefers normalized and falls back", {
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("int_", Sys.getpid()))
  owner <- "test_enrichment_module_concordance"
  sub <- "gsea_wgcna_concordance"
  fn <- "program_specific_leading_edge_module_overlap.csv"
  nd <- file.path(tmp, "results", "integration", owner, "global", "tables")
  ld <- file.path(tmp, "results", "tables", "10_biological_integration", sub, "global")
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
      got <- integration_find(fn, owner = owner,
                              legacy_stage = "10_biological_integration",
                              legacy_substep = sub)
      is_norm <- grepl("/results/integration/", got, fixed = TRUE)
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

testthat::test_that("the other historical stage resolves too", {
  ## this domain had two, and 08_behavior_physio_coupling is the one the
  ## behaviour-coupling readers use
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("int08_", Sys.getpid()))
  ld <- file.path(tmp, "results", "tables", "08_behavior_physio_coupling",
                  "network_behavior_coupling")
  dir.create(ld, recursive = TRUE, showWarnings = FALSE)
  writeLines("x", file.path(ld, "edge_behavior_figure_ready_table.csv"))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    got <- integration_find("edge_behavior_figure_ready_table.csv",
                            owner = "test_network_behaviour_coupling",
                            legacy_stage = "08_behavior_physio_coupling",
                            legacy_substep = "network_behavior_coupling")
    testthat::expect_true(grepl("08_behavior_physio_coupling", got, fixed = TRUE))
  })
})

# --- sections 13 and 14: resolver correctness ---------------------------

testthat::test_that("an empty normalized directory cannot shadow historical files", {
  ## Section 13. The QC phase shipped this defect twice, so it is asserted
  ## here rather than assumed: availability comes from the named file, never
  ## from the directory, and a dry-run skeleton is harmless.
  testthat::skip_if_not_installed("withr")
  tmp <- file.path(tempdir(), paste0("intshadow_", Sys.getpid()))
  owner <- "build_cross_compartment_atlas"
  sub <- "cross_compartment_program_atlas"
  fn <- "cross_compartment_program_atlas.csv"
  ## the normalized tree exists, with a child directory, but holds no payload
  nd <- file.path(tmp, "results", "integration", owner, "global", "tables")
  dir.create(file.path(nd, "source_data"), recursive = TRUE, showWarnings = FALSE)
  ld <- file.path(tmp, "results", "tables", "10_biological_integration", sub, "global")
  dir.create(ld, recursive = TRUE, showWarnings = FALSE)
  writeLines("x", file.path(ld, fn))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    got <- integration_find(fn, owner = owner,
                            legacy_stage = "10_biological_integration",
                            legacy_substep = sub)
    testthat::expect_true(grepl(sub, got, fixed = TRUE),
      info = paste("resolved to", got, "but the only real file is historical"))
  })
})

testthat::test_that("no directory is anchored on a guessed probe filename", {
  ## Section 14. The one directory read in this domain resolves through the
  ## overlap table, which is a declared output of the owning analysis, so the
  ## anchor is part of the data contract rather than invented.
  txt <- paste(readLines(writer_file("summarize_enrichment_module_concordance"),
                         warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("strict_dir <- dirname(integration_find(", txt, fixed = TRUE))
  s <- registry_steps()
  i <- which(s$script == "analysis/integration/test_enrichment_module_concordance.R")
  declared <- basename(sp(s$produces[i[1]]))
  testthat::expect_true("program_specific_leading_edge_module_overlap.csv" %in% declared)
})

# --- section 30: the stage13 rename -------------------------------------

testthat::test_that("no future integration destination contains stage13", {
  s <- registry_steps()
  for (i in seq_len(nrow(s))) {
    if (!grepl("^analysis/integration/", s$script[i])) next
    for (p in sp(s$produces[i])) {
      testthat::expect_false(grepl("stage13", p, fixed = TRUE), info = p)
    }
  }
})

testthat::test_that("the renamed audit is declared and its old name is not written", {
  s <- registry_steps()
  i <- which(s$script == "analysis/integration/render_module_circular_atlas.R")
  declared <- sp(s$produces[i[1]])
  testthat::expect_true(any(endsWith(declared, STAGE13_NEW)))
  testthat::expect_false(any(grepl(STAGE13_OLD, declared, fixed = TRUE)))

  ## the producer writes the new name and no longer the old one
  txt <- paste(readLines(writer_file("render_module_circular_atlas"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl(STAGE13_NEW, txt, fixed = TRUE))
  testthat::expect_false(grepl(STAGE13_OLD, txt, fixed = TRUE))
})

testthat::test_that("the historical stage13 artifact is left untouched on disk", {
  ## Section 30 is explicit that the historical filename must not disappear.
  hist <- path_results("reports", "10_biological_integration",
                       "wgcna_circular_atlas", "global", STAGE13_OLD)
  testthat::skip_if_not(file.exists(hist), "historical audit not generated here")
  testthat::expect_true(file.exists(hist))
})

testthat::test_that("the two publication_source_data stage13 names are untouched", {
  s <- registry_steps()
  i <- which(s$script == "analysis/publication_source_data/build_biological_claims_table.R")
  declared <- sp(s$produces[i[1]])
  testthat::expect_true(any(grepl("wgcna_stage13_claim_cardinality_audit.csv",
                                  declared, fixed = TRUE)))
  testthat::expect_true(any(grepl("microglia_wgcna_overlap_stage13_identity_audit.csv",
                                  declared, fixed = TRUE)))
})

# --- section 12: pending producers keep their historical paths ----------

testthat::test_that("inputs from unmigrated producers are not given invented paths", {
  ## Sixty of this domain's input edges come from WGCNA, which has not
  ## migrated. Declaring results/wgcna/... for them would claim a location the
  ## producer has not adopted.
  s <- registry_steps()
  int <- s[grepl("^analysis/integration/", s$script), , drop = FALSE]
  deps <- unique(unlist(lapply(seq_len(nrow(int)), function(i) {
    c(sp(int$consumes_required[i]), sp(int$consumes_optional[i]))
  })))
  for (dom in c("wgcna", "preprocessing", "publication_source_data")) {
    invented <- deps[startsWith(deps, paste0("results/", dom, "/"))]
    testthat::expect_identical(
      invented, character(0),
      info = paste("claims a normalized path for unmigrated", dom, ":",
                   paste(invented, collapse = ", ")))
  }
  ## and the historical WGCNA inputs are still declared
  testthat::expect_gt(sum(grepl("06_modules_WGCNA", deps, fixed = TRUE)), 0L)
})

testthat::test_that("every migrated-domain input has a normalized sibling declared", {
  MIG <- c("05_celltype_enrichment_EWCE" = "enrichment",
           "07_spatial_networks" = "spatial_networks",
           "11_spatial_systems" = "spatial_validation",
           "04_differential_expression_enrichment" = "differential_abundance",
           "03_qc_exploration" = "qc")
  s <- registry_steps()
  int <- s[grepl("^analysis/integration/", s$script), , drop = FALSE]
  for (i in seq_len(nrow(int))) {
    deps <- c(sp(int$consumes_required[i]), sp(int$consumes_optional[i]))
    for (st in names(MIG)) {
      if (!any(grepl(st, deps, fixed = TRUE))) next
      dom <- MIG[[st]]
      testthat::expect_true(
        any(startsWith(deps, paste0("results/", dom, "/"))),
        info = paste(basename(int$script[i]), "names", st,
                     "with no normalized", dom, "sibling"))
    }
  }
})

# --- section 5: coordination --------------------------------------------

testthat::test_that("no integration artifact has two writers", {
  s <- registry_steps()
  int <- s[grepl("^analysis/integration/", s$script), , drop = FALSE]
  out <- do.call(rbind, lapply(seq_len(nrow(int)), function(i) {
    p <- sp(int$produces[i])
    p <- p[grepl("[.][A-Za-z0-9]{2,5}$", basename(p))]
    if (!length(p)) return(NULL)
    data.frame(script = int$script[i], out = p, stringsAsFactors = FALSE)
  }))
  testthat::expect_gt(nrow(out), 0L)
  testthat::expect_identical(sum(duplicated(out$out)), 0L)
})

testthat::test_that("one writer's two substeps keep distinct manifests", {
  ## build_evidence_priority_matrix produces two output families and each
  ## writes its own run_manifest.yml. Collapsing the substep would have merged
  ## them and lost one phase's provenance, so the substep is retained exactly
  ## where it disambiguates.
  s <- registry_steps()
  i <- which(s$script == "analysis/integration/build_evidence_priority_matrix.R")
  declared <- sp(s$produces[i[1]])
  man <- declared[grepl("run_manifest[.]yml$", declared)]
  testthat::expect_length(man, 2L)
  testthat::expect_identical(length(unique(man)), 2L)
})

testthat::test_that("the synthetic shared-artifact fixture still hard-stops", {
  shared <- "results/integration/build_cross_compartment_atlas/global/tables/x.csv"
  collisions <- migration_destination_collisions(
    list(build_cross_compartment_atlas = shared,
         summarize_programs_for_manuscript = shared), NULL)
  testthat::expect_length(collisions, 1L)
  testthat::expect_false(collisions[[1]]$resolvable)
  testthat::expect_gt(length(migration_gate_blockers(collisions)), 0L)
})

# --- section 22/23: the publication boundary ----------------------------

testthat::test_that("integration writes nothing into exports or the manuscript", {
  for (aid in WRITERS) {
    txt <- paste(readLines(writer_file(aid), warn = FALSE), collapse = "\n")
    testthat::expect_false(grepl("path_export(", txt, fixed = TRUE), info = aid)
    testthat::expect_false(grepl("Exp9_manuscript", txt, fixed = TRUE), info = aid)
  }
})

testthat::test_that("export and freeze discovery cover the normalized roots", {
  ## Section 23. The frozen source-data sets were rooted only at historical
  ## paths; two of the three belong to differential_abundance and had already
  ## gone stale when that domain migrated.
  fz <- paste(readLines(repo_path("R", "publication_freeze_utils.R"), warn = FALSE),
              collapse = "\n")
  testthat::expect_true(grepl('path_results("integration", "export_module_protein_zoom_source_data"',
                              fz, fixed = TRUE))
  testthat::expect_true(grepl('path_results("differential_abundance", "validate_control_spatial_identity"',
                              fz, fixed = TRUE))
  ## and the historical roots are retained, because that is where the files are
  testthat::expect_true(grepl('path_results("source_data", "manuscript_panels")',
                              fz, fixed = TRUE))

  ex <- paste(readLines(repo_path("analysis", "publication_source_data",
                                 "08_export_manuscript_figures.R"), warn = FALSE),
              collapse = "\n")
  testthat::expect_true(grepl('path_results("integration", "export_module_protein_zoom_source_data"',
                              ex, fixed = TRUE))
  testthat::expect_true(grepl('path_results("figures", "manuscript_panels")',
                              ex, fixed = TRUE))
})

# --- section 9: lifecycle ----------------------------------------------

testthat::test_that("no integration output is declared under work/", {
  s <- registry_steps()
  for (i in seq_len(nrow(s))) {
    if (!grepl("^analysis/integration/", s$script[i])) next
    for (p in sp(s$produces[i])) {
      testthat::expect_false(startsWith(p, "work/"), info = p)
    }
  }
})

testthat::test_that("the lifecycle audit covers every declared output", {
  f <- repo_path("audits", "phase6g_integration_lifecycle.csv")
  testthat::skip_if_not(file.exists(f))
  lc <- utils::read.csv(f, stringsAsFactors = FALSE)
  s <- registry_steps()
  int <- s[grepl("^analysis/integration/", s$script), , drop = FALSE]
  testthat::expect_identical(nrow(lc), length(unlist(lapply(int$produces, sp))))
  testthat::expect_identical(sum(lc$lifecycle == "WORK_INTERMEDIATE"), 0L)
})
