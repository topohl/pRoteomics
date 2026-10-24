source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "empirical_roi_marker_utils.R"))

# A regeneration tree without legacy results showed that
# quantify_empirical_compartments read a crosswalk only the proposed-only
# branch of discover_empirical_roi_markers writes. It now rebuilds it from the
# canonical feature tables. Only the script's helper definitions are evaluated.

QEC <- repo_path("analysis", "spatial_validation", "quantify_empirical_compartments.R")
VSF <- repo_path("analysis", "spatial_validation", "validate_spatial_foundations.R")

is_assign_to <- function(e, nm) {
  is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name(nm))
}
literals <- function(e) {
  if (is.character(e)) return(e)
  if (!is.call(e)) return(character())
  unlist(lapply(as.list(e)[-1], function(x) tryCatch(literals(x), error = function(...) NULL)))
}

testthat::test_that("the compartment crosswalk is rebuilt from the canonical feature tables", {
  testthat::skip_if_not_installed("withr")
  exprs <- parse(QEC, keep.source = FALSE)
  testthat::expect_false(any(grepl("crosswalk_proposed", unlist(lapply(exprs, literals)),
                                   fixed = TRUE)))
  env <- new.env(parent = globalenv())
  for (nm in c("feature_table_path", "read_feature_table", "FEATURE_DATASETS")) {
    hit <- Filter(function(e) is_assign_to(e, nm), exprs)
    testthat::expect_length(hit, 1L)
    eval(hit[[1]], env)
  }
  cw_expr <- Filter(function(e) is_assign_to(e, "cw") &&
                      "build_empirical_roi_protein_index" %in% all.names(e), exprs)
  testthat::expect_length(cw_expr, 1L)

  ft <- list(
    neuron_neuropil = data.frame(ProteinGroupID = c("PG:n:1", "PG:n:2"),
                                 member_accessions = c("P2;P1", NA),
                                 mapping_status = c("mapped", "unmapped")),
    neuron_soma = data.frame(ProteinGroupID = "PG:s:1", member_accessions = "P1;P2",
                             mapping_status = "mapped"),
    microglia = data.frame(ProteinGroupID = "PG:m:1", member_accessions = "P1; P2",
                           mapping_status = "mapped"))
  tmp <- tempfile("qec_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  withr::with_envvar(c(PROTEOMICS_PROJECT_ROOT = tmp), {
    for (ds in names(ft)) {
      p <- env$feature_table_path(ds)
      testthat::expect_true(grepl(paste0("results/qc/discover_empirical_roi_markers/global/tables/",
                                         ds, "/canonical_feature_table.csv"), p, fixed = TRUE))
      dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(cbind(ft[[ds]], extra = 1), p, row.names = FALSE)  # as qc_write_csv
    }
    eval(cw_expr[[1]], env)
    want <- as.data.frame(build_empirical_roi_protein_index(ft))
    testthat::expect_identical(env$cw, want)
    testthat::expect_identical(sum(env$cw$cross_dataset_model_eligible), 3L)

    # a missing canonical input stops the script; there is no legacy fallback
    unlink(env$feature_table_path("microglia"))
    testthat::expect_error(eval(cw_expr[[1]], env),
                           "missing_required_input: .*microglia/canonical_feature_table[.]csv")
  })
})

testthat::test_that("foundation check 7 names the missing compartment summary and its producer", {
  hit <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("add")) &&
                  identical(e[[2]], "compartment_validation_phenotype_independent"),
                parse(VSF, keep.source = FALSE))
  testthat::expect_length(hit, 1L)
  run <- function(comp, leak) {
    env <- new.env(parent = baseenv())
    env$comp <- comp; env$leak <- leak
    env$add <- function(check_id, critical, status, detail) {
      list(critical = critical, status = status, detail = detail)
    }
    eval(hit[[1]], env)
  }
  miss <- run(NULL, FALSE)
  testthat::expect_true(miss$critical)
  testthat::expect_identical(miss$status, "FAIL")
  testthat::expect_identical(miss$detail, paste("bilateral_empirical_compartment_summary.csv",
                                                "missing (quantify_empirical_compartments)"))
  ok <- run(data.frame(x = 1), FALSE)
  testthat::expect_identical(ok$status, "PASS")
  testthat::expect_match(ok$detail, "^design is ~ AnimalID \\+ dataset")
  testthat::expect_identical(run(data.frame(x = 1), TRUE)$status, "FAIL")
})
