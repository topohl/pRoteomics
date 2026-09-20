source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "output_namespace_classification.R"))

# Phase 6G.8: the generalized legacy-write classification.
#
# The rule this replaces classified a destination as historical only when one
# of its literal segments matched a numbered stage namespace. That is a proxy,
# and it had a real cost: results/reviewer_audit/ carries no stage number, so
# five WGCNA writers wrote there while their pipeline.yml declarations already
# said results/wgcna/..., and the split-brain gate reported zero. The fixture
# in the first test below is exactly that case.
#
# These tests exist so the rule can be checked on fixtures without running the
# whole audit, which is the only way a detector correction can be trusted.

q <- function(txt) parse(text = txt)[[1]]

# --- the fixture the correction was made for -----------------------------

testthat::test_that("a historical root with no stage number is still historical", {
  ## results/reviewer_audit/example.csv
  testthat::expect_true(
    is_legacy_path_construction(q('path_results("reviewer_audit", "example.csv")')))
  testthat::expect_true(is_legacy_output_path("results/reviewer_audit/example.csv"))

  ## and the nested WGCNA-owned families under it
  testthat::expect_true(
    is_legacy_path_construction(q('path_results("reviewer_audit", "wgcna_label_approval")')))
  testthat::expect_true(
    is_legacy_path_construction(q('path_results("reviewer_audit", "wgcna_label_review", dataset)')))
})

testthat::test_that("the numbered historical roots are still caught", {
  for (e in c('path_results("tables", "06_modules_WGCNA", "01_WGCNA", ds, "modules")',
              'path_results("figures", "03_qc_exploration", "00b_joint_compartment_qc")',
              'path_results("logs", "11_spatial_systems", "networks")',
              'path_results("source_data", "04_differential_expression_enrichment")')) {
    testthat::expect_true(is_legacy_path_construction(q(e)), info = e)
  }
})

testthat::test_that("data/processed is not a canonical output lifecycle", {
  ## config/output_layout.yml declares exactly three: work, results, exports.
  testthat::expect_true(
    is_legacy_path_construction(q('path_processed("06_modules_WGCNA", "01_WGCNA", ds)')))
  testthat::expect_true(
    is_legacy_path_construction(q('path_processed("01_preprocessing", "impute")')))
  testthat::expect_true(is_legacy_output_path("data/processed/06_modules_WGCNA/01_WGCNA"))
})

# --- what must NOT be flagged -------------------------------------------

testthat::test_that("a canonical results path is not historical", {
  domains <- onc_canonical_domains()
  testthat::expect_gt(length(domains), 0L)
  for (d in domains) {
    e <- sprintf('path_results("%s", "some_analysis", "global", "tables")', d)
    testthat::expect_false(is_legacy_path_construction(q(e)), info = e)
    testthat::expect_false(
      is_legacy_output_path(sprintf("results/%s/some_analysis/global/tables/x.csv", d)),
      info = d)
  }
})

testthat::test_that("the adjudicated non-results destinations are exempt", {
  ## config/ is a generated configuration contract (adjudicated in 6G.5);
  ## exports/ is the frozen bundle; pride_submission/ is export staging.
  for (r in ONC_ALLOWED_NONRESULT_ROOTS) {
    testthat::expect_false(is_legacy_output_path(paste0(r, "/whatever/x.csv")), info = r)
  }
  testthat::expect_false(is_legacy_output_path("work/wgcna/build_wgcna_modules/global/x.rds"))
})

testthat::test_that("a read is not a destination", {
  ## The predicate answers "is this destination historical", not "is this a
  ## write". A migrated writer reading a historical input is not a legacy
  ## write, which is why the audit consults this only for destination
  ## arguments of write calls.
  testthat::expect_false(is_legacy_path_construction(q('readRDS("x.rds")')))
  testthat::expect_false(is_legacy_path_construction(q('list.files(d)')))
})

testthat::test_that("the legacy directory factories are historical by definition", {
  for (f in ONC_LEGACY_FACTORIES) {
    e <- sprintf('%s("06_modules_WGCNA", "sub")', f)
    testthat::expect_true(is_legacy_path_construction(q(e)), info = f)
  }
})

# --- the write-call inventory the audit uses -----------------------------

testthat::test_that("the write-call list knows the spellings this codebase uses", {
  ## write.csv and write.xlsx were known; write_csv (readr) and write_xlsx
  ## (writexl) were not, and between them they have hundreds of call sites in
  ## analysis/, every one invisible to the audit before Phase 6G.8.
  src <- readLines(repo_path("tools", "audit_writer_namespaces.R"), warn = FALSE)
  txt <- paste(src, collapse = "\n")
  for (nm in c("write_csv", "write_xlsx", "write_json", "write_delim",
               "write_lines", "write_rds", "svglite")) {
    testthat::expect_true(grepl(paste0('"', nm, '"'), txt, fixed = TRUE),
                          info = paste(nm, "missing from the write-call list"))
  }
  ## and each has a destination spec, so its payload arguments are not tested
  for (nm in c("write_csv", "write_xlsx", "write_json", "write_delim",
               "write_lines", "write_rds", "svglite")) {
    testthat::expect_true(grepl(paste0(nm, " = list("), txt, fixed = TRUE),
                          info = paste(nm, "has no DEST_ARG entry"))
  }
})

testthat::test_that("cat is a write only when it is given a destination", {
  ## cat() has hundreds of call sites that print to the console. Counting it
  ## unconditionally would make every legacy path named in a message a legacy
  ## write, which is the false-positive class Phase 6G.4 removed.
  txt <- paste(readLines(repo_path("tools", "audit_writer_namespaces.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl("CONDITIONAL_WRITE_CALLS", txt, fixed = TRUE))
  testthat::expect_true(grepl('cat = "file"', txt, fixed = TRUE))
})

# --- the audit and the library agree ------------------------------------

testthat::test_that("the audit delegates to this library rather than copying it", {
  txt <- paste(readLines(repo_path("tools", "audit_writer_namespaces.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl("output_namespace_classification.R", txt, fixed = TRUE))
  testthat::expect_true(grepl("is_legacy_path_construction(", txt, fixed = TRUE))
  ## the old proxy must not still be the deciding rule on its own
  testthat::expect_false(
    grepl("return(any(grepl(STAGE_NS, lits)))", txt, fixed = TRUE))
})

testthat::test_that("every LEGACY_READ_ONLY root is classified historical", {
  reg <- repo_path("config", "legacy_output_registry.csv")
  testthat::skip_if_not(file.exists(reg))
  d <- utils::read.csv(reg, stringsAsFactors = FALSE)
  testthat::expect_gt(nrow(d), 0L)
  ro <- d[d$policy == "LEGACY_READ_ONLY", , drop = FALSE]
  testthat::expect_gt(nrow(ro), 0L)
  ## results/publication_source_data is the one LEGACY_READ_ONLY root whose
  ## segment is also a declared domain name - the pre-6F bundle location - so
  ## the shape check cannot distinguish it from that domain's canonical root.
  ## It is excluded deliberately, and that domain is still unmigrated anyway.
  ro <- ro[ro$legacy_path != "results/publication_source_data", , drop = FALSE]
  bad <- ro$legacy_path[!vapply(ro$legacy_path, is_legacy_output_path, logical(1))]
  testthat::expect_identical(bad, character(0),
    info = paste("registered legacy roots not classified historical:",
                 paste(utils::head(bad, 5), collapse = ", ")))
})

testthat::test_that("an ACTIVE_NOT_LEGACY domain root is not treated as historical", {
  reg <- repo_path("config", "legacy_output_registry.csv")
  testthat::skip_if_not(file.exists(reg))
  d <- utils::read.csv(reg, stringsAsFactors = FALSE)
  act <- d[d$policy == "ACTIVE_NOT_LEGACY", , drop = FALSE]
  testthat::skip_if(nrow(act) == 0L)
  dom <- onc_canonical_domains()
  seg <- vapply(strsplit(act$legacy_path, "/", fixed = TRUE),
                function(p) if (length(p) >= 2L) p[[2]] else NA_character_,
                character(1))
  ## the normalized domain roots appear here with live writers; none of them
  ## may leak into the historical-segment fallback
  leaked <- intersect(seg[!is.na(seg) & seg %in% dom],
                      onc_legacy_registered_segments())
  testthat::expect_identical(leaked, character(0),
    info = paste("canonical domain roots leaked into the legacy segments:",
                 paste(leaked, collapse = ", ")))
})
