source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "wgcna_label_adjudication_utils.R"))
source(repo_path("R", "wgcna_label_activation_utils.R"))

testthat::local_edition(3)

# In provisional label mode (a rebuilt network) the reviewed registries in config/wgcna_labels/ name the
# entities of the previous network, so no consumer may attach them by entity ID.

plm_fx_registry_dir <- function() {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  utils::write.csv(data.frame(dataset = "d", level = c("module", "supermodule"),
                              entity_id = c("WGCNA_m01", "SM01"),
                              reviewed_biological_label = c("old network one", "old network sm"),
                              confidence = "high", adjudication_status = "reviewed",
                              reviewer = "A Human", stringsAsFactors = FALSE),
                   file.path(dir, "d.csv"), row.names = FALSE)
  dir
}
plm_fx_stage07 <- function(ds) {
  data.frame(level = c("module", "supermodule"), entity_id = c("WGCNA_m01", "SM01"),
             canonical_biological_label = c("this run one", "this run sm"), stringsAsFactors = FALSE)
}
plm_fx_claims <- function() {
  data.frame(dataset = "d", wgcna_level = c("module", "supermodule", NA),
             wgcna_entity_id = c("WGCNA_m01", "SM01", NA), claim = c("a", "b", "c"),
             stringsAsFactors = FALSE)
}

testthat::test_that("the label mode is reviewed by default and only reviewed or provisional", {
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = NA)
  testthat::expect_identical(wgcna_label_mode(), "reviewed")
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = "Provisional")
  testthat::expect_identical(wgcna_label_mode(), "provisional")
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = "draft")
  testthat::expect_error(wgcna_label_mode(), "reviewed or provisional")
})

testthat::test_that("the reviewed registries are read in reviewed mode and never in provisional mode", {
  dir <- plm_fx_registry_dir()
  rev <- wal_reviewed_registries("d", dir = dir, label_mode = "reviewed")
  testthat::expect_identical(nrow(rev), 2L)
  testthat::expect_identical(unique(rev$dataset), "d")
  testthat::expect_null(wal_reviewed_registries("d", dir = dir, label_mode = "provisional"))
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = "provisional")
  testthat::expect_null(wal_reviewed_registries("d", dir = dir))
})

testthat::test_that("the claims table shows the reviewed label in reviewed mode and this run's label in provisional mode", {
  dir <- plm_fx_registry_dir()
  reviewed <- attach_canonical_wgcna_display_label(plm_fx_claims(), reviewed_dir = dir,
                                                   stage07_reader = plm_fx_stage07, label_mode = "reviewed")
  testthat::expect_identical(reviewed$canonical_display_label, c("old network one", "old network sm", NA))
  testthat::expect_identical(reviewed$canonical_label_source[1:2], rep("active_reviewed_registry", 2))

  provisional <- attach_canonical_wgcna_display_label(plm_fx_claims(), reviewed_dir = dir,
                                                      stage07_reader = plm_fx_stage07, label_mode = "provisional")
  testthat::expect_identical(provisional$canonical_display_label, c("this run one", "this run sm", NA))
  testthat::expect_identical(provisional$canonical_label_source[1:2], rep("stage07_canonical", 2))
  testthat::expect_identical(provisional$claim, plm_fx_claims()$claim)

  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = "provisional")
  from_env <- attach_canonical_wgcna_display_label(plm_fx_claims(), reviewed_dir = dir, stage07_reader = plm_fx_stage07)
  testthat::expect_identical(from_env$canonical_display_label, provisional$canonical_display_label)
})

testthat::test_that("every script that attaches registry labels reads the mode through wgcna_label_mode()", {
  scripts <- c("analysis/wgcna/build_module_label_registry.R", "analysis/wgcna/summarize_module_interpretation.R",
               "analysis/wgcna/audit_microglia_module_claims.R")
  for (s in scripts) {
    src <- paste(readLines(repo_path(s), warn = FALSE), collapse = "\n")
    testthat::expect_match(src, "wgcna_label_mode()", fixed = TRUE, info = s)
    testthat::expect_false(grepl("Sys.getenv(\"PROTEOMICS_WGCNA_LABEL_MODE\"", src, fixed = TRUE), info = s)
  }
  registry <- paste(readLines(repo_path("analysis/wgcna", "build_module_label_registry.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(registry, "wal_reviewed_registries(", fixed = TRUE)
})

testthat::test_that("the microglia claim audit never falls back to the reviewed registry in provisional mode", {
  src <- readLines(repo_path("analysis/wgcna", "audit_microglia_module_claims.R"), warn = FALSE)
  first <- grep("^registry_path <- if \\(identical\\(wgcna_label_mode\\(\\), \"provisional\"\\)\\)", src)
  testthat::expect_length(first, 1L)
  last <- first + grep("^\\} else file.path\\(audit_root, \"config/wgcna_labels/microglia.csv\"\\)", src[first:length(src)])[1] - 1L
  env <- new.env()
  env$audit_root <- "root"
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = "provisional", PROTEOMICS_WGCNA_PROVISIONAL_REGISTRY = NA)
  testthat::expect_error(eval(parse(text = src[first:last]), envir = env), "PROTEOMICS_WGCNA_PROVISIONAL_REGISTRY")
  withr::local_envvar(PROTEOMICS_WGCNA_PROVISIONAL_REGISTRY = "provisional.csv")
  eval(parse(text = src[first:last]), envir = env)
  testthat::expect_identical(env$registry_path, "provisional.csv")
  withr::local_envvar(PROTEOMICS_WGCNA_LABEL_MODE = NA)
  eval(parse(text = src[first:last]), envir = env)
  testthat::expect_identical(env$registry_path, file.path("root", "config/wgcna_labels/microglia.csv"))
})
