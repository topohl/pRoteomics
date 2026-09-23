# Phase 6H.9 curation adjudication, pinned.
#
# Two decisions were taken and deliberately NOT implemented:
#   code     C2 - archive the unreachable compare_go_enrichment tail
#   artifact P2 - retain the orphaned bootstrap workbook as internal provenance
#                 only, i.e. out of the outward PRIDE deposition selection
#
# These tests pin the facts the decisions rest on, so that if any of them stops
# being true the adjudication is revisited rather than silently inherited. They
# assert the CURRENT (un-implemented) state; implementing C2/P2 is expected to
# change several of them, and those changes should be deliberate.

source(testthat::test_path("..", "..", "R", "paths.R"))

SCRIPT <- repo_path("analysis", "differential_abundance", "compare_go_enrichment.R")
WORKBOOK <- repo_path(
  "pride_submission", "supplementary_tables",
  paste0("results_tables_04_differential_expression_enrichment_compareGO_",
         "neuron_neuropil_BP_phenotype_within_unit_08_Bootstrap_Stability_Summary.xlsx"))
WORKBOOK_SHA <- "024d3671f2cd6026e8bfd7feb6d9839c7ff23eeb41c98d8f66235d854185baba"
INVENTORY <- repo_path("audits", "phase6h_orphan_bootstrap_artifact_inventory.csv")
TAIL_CSV <- repo_path("audits", "phase6h_unreachable_comparego_tail.csv")

# ---- the unreachability proof ---------------------------------------------

testthat::test_that("the exit is a bare unconditional top-level quit()", {
  testthat::skip_if_not(file.exists(SCRIPT), "compare_go_enrichment.R absent")
  ex <- parse(SCRIPT, keep.source = TRUE)
  is_quit <- vapply(seq_along(ex), function(i) {
    e <- ex[[i]]
    is.call(e) && identical(as.character(e[[1]])[1], "quit")
  }, logical(1))
  testthat::expect_true(any(is_quit))
  idx <- min(which(is_quit))
  e <- ex[[idx]]
  # constant arguments only - a computed status could in principle be a no-op
  args <- as.list(e)[-1]
  testthat::expect_true(all(vapply(args, function(a) !is.call(a) && !is.name(a), logical(1))))
  # Until Phase 6H.10 this also asserted that the bootstrap's slice_sample sat
  # AFTER the exit. Decision C2 archived that region, so the exit is now the
  # final expression and there is no code behind it to order against.
  testthat::expect_identical(idx, length(ex))
  pd <- utils::getParseData(ex)
  testthat::expect_identical(
    sum(pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "slice_sample"), 0L)
})

testthat::test_that("nothing can bypass the exit by masking quit()", {
  # The only structural bypass in R: `quit` rebound to something that does not
  # exit. Checked in the script and in every library it sources before the exit.
  testthat::skip_if_not(file.exists(SCRIPT), "compare_go_enrichment.R absent")
  lines <- readLines(SCRIPT, warn = FALSE)
  patterns <- c("quit <-", "quit<-", "quit =", "assign(\"quit\"", "assign('quit'")
  testthat::expect_identical(
    length(unlist(lapply(patterns, function(p) grep(p, lines, fixed = TRUE)))), 0L)
  testthat::expect_identical(length(grep("\\btrace\\(", lines)), 0L)

  ex <- parse(SCRIPT, keep.source = TRUE)
  sr <- attr(ex, "srcref")
  is_quit <- vapply(seq_along(ex), function(i) {
    e <- ex[[i]]; is.call(e) && identical(as.character(e[[1]])[1], "quit")
  }, logical(1))
  exit_line <- sr[[min(which(is_quit))]][1]

  head_lines <- lines[seq_len(exit_line)]
  srcs <- unlist(regmatches(head_lines, gregexpr('repo_path\\("R"[^)]*\\)', head_lines)))
  libs <- unique(gsub('.*"([^"]+\\.[Rr])".*', "\\1", grep('\\.[Rr]"', srcs, value = TRUE)))
  checked <- 0L
  for (lb in libs) {
    cand <- c(repo_path("R", lb), Sys.glob(file.path(repo_path("R"), "*", lb)))
    cand <- cand[file.exists(cand)]
    if (!length(cand)) next
    checked <- checked + 1L
    l <- readLines(cand[1], warn = FALSE)
    testthat::expect_identical(
      length(unlist(lapply(patterns, function(p) grep(p, l, fixed = TRUE)))), 0L,
      info = paste("quit() masked in", basename(cand[1])))
  }
  testthat::expect_gt(checked, 0L)
})

testthat::test_that("nothing sources, parses or evals the script", {
  # If any file sourced it, the tail's definitions would be importable and both
  # the unreachability argument and the archive decision would change.
  roots <- c("analysis", "R", "tools", "tests")
  files <- unlist(lapply(roots, function(d) {
    p <- repo_path(d)
    if (dir.exists(p)) list.files(p, pattern = "[.][Rr]$", recursive = TRUE,
                                  full.names = TRUE) else character(0)
  }))
  offenders <- character(0)
  for (g in files) {
    if (normalizePath(g, winslash = "/") == normalizePath(SCRIPT, winslash = "/")) next
    l <- readLines(g, warn = FALSE)
    hit <- grep("compare_go_enrichment", l, fixed = TRUE)
    for (h in hit) {
      if (grepl("source\\s*\\(", l[h]) || grepl("\\beval\\s*\\(", l[h]))
        offenders <- c(offenders, paste0(basename(g), ":", h))
    }
  }
  testthat::expect_identical(offenders, character(0),
    info = paste("script is executed elsewhere:", paste(offenders, collapse = ", ")))
})

testthat::test_that("no tail-defined function is depended on from outside", {
  # jaccard, mode_value and optional_read_csv share names with definitions in
  # other files. Those are collisions, not dependencies - each other file must
  # define its own, or archiving the tail would break it.
  for (nm in c("jaccard", "mode_value", "optional_read_csv")) {
    roots <- c(repo_path("analysis"), repo_path("R"), repo_path("tools"))
    files <- unlist(lapply(roots[dir.exists(roots)], list.files,
                           pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE))
    users <- character(0); definers <- character(0)
    for (g in files) {
      if (normalizePath(g, winslash = "/") == normalizePath(SCRIPT, winslash = "/")) next
      l <- readLines(g, warn = FALSE)
      if (any(grepl(paste0("^\\s*", nm, "\\s*(<-|=)\\s*function"), l))) definers <- c(definers, g)
      pd <- tryCatch(utils::getParseData(parse(g, keep.source = TRUE)), error = function(e) NULL)
      if (!is.null(pd) && any(pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == nm))
        users <- c(users, g)
    }
    # every external user also defines it locally
    testthat::expect_identical(setdiff(users, definers), character(0),
      info = paste(nm, "is called outside without a local definition"))
  }
})

# ---- the tail's shape ------------------------------------------------------

testthat::test_that("the unreachable tail metrics are as adjudicated", {
  testthat::skip_if_not(file.exists(TAIL_CSV), "tail audit table absent")
  a <- utils::read.csv(TAIL_CSV, stringsAsFactors = FALSE)
  get <- function(k) a$value[a$item == k][1]
  testthat::expect_identical(get("exit_line"), "569")
  testthat::expect_identical(get("unreachable_region_start_line"), "572")
  testthat::expect_identical(get("unreachable_region_end_line"), "4321")
  testthat::expect_identical(get("quit_masked_in_file"), "FALSE")
  testthat::expect_identical(get("quit_masked_in_sourced_libs"), "FALSE")
  testthat::expect_identical(get("tail_functions_used_outside"), "none")
})

testthat::test_that("the tests that drew assurance from dead code were rewritten", {
  # Phase 6H.9 recorded four tests coupled to the tail's literal text and
  # predicted they would break under C2. Phase 6H.10 rewrote each against a
  # live invariant. This now asserts the repair rather than the coupling, so a
  # regression to text-matching dead code is caught.
  t1 <- repo_path("tests", "testthat", "test-comparego-canonical-contract.R")
  testthat::skip_if_not(file.exists(t1), "canonical contract test absent")
  # it no longer splits the script at the obsolete marker
  testthat::expect_false(any(grepl("marker <- grep", readLines(t1, warn = FALSE), fixed = TRUE)))

  t2 <- repo_path("tests", "testthat", "test-protein-group-enrichment-utils.R")
  testthat::skip_if_not(file.exists(t2), "protein group test absent")
  s2 <- readLines(t2, warn = FALSE)
  # it no longer proves runtime behaviour by grepping the script for tokens
  # that lived only below the exit; it executes the manifest contract instead
  testthat::expect_true(any(grepl("clusterprofiler_runtime_required_fields", s2, fixed = TRUE)))
  testthat::expect_true(any(grepl("validate_clusterprofiler_manifest_contract", s2, fixed = TRUE)))

  # and the tokens themselves are gone from the active script, as expected
  testthat::skip_if_not(file.exists(SCRIPT), "compare_go_enrichment.R absent")
  lines <- readLines(SCRIPT, warn = FALSE)
  for (tok in c("comparison_input_file", "GeneSymbol"))
    testthat::expect_identical(length(grep(tok, lines, fixed = TRUE)), 0L, info = tok)
})

# ---- the artifact ---------------------------------------------------------

testthat::test_that("the orphaned workbook is unchanged and still a 1x5 sheet", {
  # Follows the artifact to where Phase 6H.10 put it. Pointing this at the old
  # outward path would have turned the whole block into a silent skip the moment
  # P2 was implemented - the same defect this phase criticised in the 6H.8 guard
  # and fixed there.
  wb <- path_results("manuscript", "_curated",
    paste0("results_tables_04_differential_expression_enrichment_compareGO_",
           "neuron_neuropil_BP_phenotype_within_unit_",
           "08_Bootstrap_Stability_Summary.xlsx"))
  testthat::expect_true(file.exists(wb),
    info = "the curated provenance copy is missing")
  testthat::skip_if_not(file.exists(wb))
  testthat::expect_identical(file.size(wb), 5141)
  testthat::expect_identical(unname(tools::sha256sum(wb)), WORKBOOK_SHA)
  testthat::skip_if_not_installed("readxl")
  sh <- readxl::excel_sheets(wb)
  testthat::expect_identical(length(sh), 1L)
  d <- readxl::read_excel(wb, sheet = sh[1])
  testthat::expect_identical(dim(d), c(1L, 5L))
  testthat::expect_identical(names(d), c("Mean_Recovery_Rate", "SD_Recovery_Rate",
                                         "Min_Recovery", "Max_Recovery", "Total_TopTerms"))
  # no metadata or provenance sheet - central to the misleading-in-deposited-form finding
  testthat::expect_false(any(grepl("meta|provenance|readme|session|param", sh,
                                   ignore.case = TRUE)))
})

testthat::test_that("the artifact has no manuscript or source-data dependency", {
  testthat::skip_if_not(file.exists(INVENTORY), "artifact inventory absent")
  a <- utils::read.csv(INVENTORY, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(a), 5L)
  # uniquely identified
  testthat::expect_identical(anyDuplicated(a$sha256), 0L)
  # zero manuscript references is the fact P2 rests on
  testthat::expect_identical(unique(a$manuscript_references), 0L)
  testthat::expect_identical(unique(a$doc_references), 0L)
  testthat::expect_identical(unique(a$statistic_name_refs_in_manuscript), 0L)
  # exactly one manifest reference, and it is the rule-generated PRIDE row
  testthat::expect_identical(unique(a$manifest_references), 1L)
})

testthat::test_that("deposition inclusion is rule-based, not curated", {
  # DEPOSIT_UNDECLARED_EXTRA rests on this: the flag comes from a path match.
  eh <- repo_path("R", "utilities", "export_helpers.R")
  testthat::skip_if_not(file.exists(eh), "export_helpers.R absent")
  src <- readLines(eh, warn = FALSE)
  testthat::expect_true(any(grepl('grepl("/pride_submission/", p)', src, fixed = TRUE)))
  testthat::expect_true(any(grepl("pride_staging", src, fixed = TRUE)))
  # and the manifest builder rescans the directory, re-endorsing whatever is there
  bp <- repo_path("analysis", "publication_source_data", "build_pride_manifest.R")
  if (file.exists(bp)) {
    b <- readLines(bp, warn = FALSE)
    testthat::expect_true(any(grepl("pride_submission_dir\\(\\)|list_files_shallow", b)))
  }
})

# ---- the closed-form structure -------------------------------------------

testthat::test_that("the statistic ignores effect size, p-value magnitude and direction", {
  # The block consults only the boolean p.adjust < 0.05, so identical
  # significance patterns with wildly different biology must score identically.
  # This is the basis for classifying the artifact MISLEADING_OR_UNINTERPRETABLE.
  testthat::skip_if_not_installed("dplyr")
  set.seed(11)
  comps <- paste0("C", 1:4)
  terms <- paste0("GO_", sprintf("%03d", 1:80))
  base <- do.call(rbind, lapply(comps, function(cc) {
    k <- sample(seq_along(terms), 50L)
    data.frame(Comparison = cc, Description = terms[k],
               is_sig = stats::runif(50L) < 0.12, stringsAsFactors = FALSE)
  }))
  top <- utils::head(unique(base$Description[base$is_sig]), 12L)

  stat <- function(df, seed = 1L, n_boot = 30L) {
    set.seed(seed)
    res <- integer(n_boot)
    for (b in seq_len(n_boot)) {
      s <- df |>
        dplyr::group_by(.data$Comparison) |>
        dplyr::slice_sample(prop = 1, replace = TRUE) |>
        dplyr::ungroup()
      res[b] <- sum(top %in% unique(s$Description[s$p.adjust < 0.05]))
    }
    mean(res) / length(top)
  }

  weak <- base; weak$p.adjust <- ifelse(weak$is_sig, 0.0499, 0.51)
  huge <- base; huge$p.adjust <- ifelse(huge$is_sig, 1e-30, 0.99)
  testthat::expect_identical(stat(weak), stat(huge))
})

testthat::test_that("the statistic tracks 1 - exp(-M) and nothing else", {
  testthat::skip_if_not_installed("dplyr")
  comps <- paste0("C", 1:5)
  terms <- paste0("GO_", sprintf("%03d", 1:60))
  top <- terms[1:12]
  mk <- function(mult) {
    rows <- do.call(rbind, lapply(comps, function(cc)
      data.frame(Comparison = cc, Description = terms, stringsAsFactors = FALSE)))
    rows$p.adjust <- 0.5
    for (t in top) {
      idx <- which(rows$Description == t)
      rows$p.adjust[utils::head(idx, mult)] <- 0.001
    }
    rows
  }
  stat <- function(df, n_boot = 30L) {
    set.seed(7)
    res <- integer(n_boot)
    for (b in seq_len(n_boot)) {
      s <- df |>
        dplyr::group_by(.data$Comparison) |>
        dplyr::slice_sample(prop = 1, replace = TRUE) |>
        dplyr::ungroup()
      res[b] <- sum(top %in% unique(s$Description[s$p.adjust < 0.05]))
    }
    mean(res) / length(top)
  }
  for (m in 1:3) {
    observed <- stat(mk(m))
    asymptote <- 1 - exp(-m)
    # finite stratum size puts the observed value slightly ABOVE the asymptote
    testthat::expect_gt(observed, asymptote - 0.12)
    testthat::expect_lt(observed, asymptote + 0.12)
  }
})

testthat::test_that("the recorded precision exceeds the Monte Carlo resolution", {
  # Every surviving copy reports more decimal places than its own
  # se = SD/sqrt(100) can support. Supports classification D/E.
  vals <- data.frame(
    mean_rate = c(0.8738182, 0.8846032, 0.7208, 0.650625, 0.6535714),
    sd_rate = c(0.04100906, 0.0322909, 0.09006933, 0.124714, 0.1348304))
  se <- vals$sd_rate / sqrt(100)
  # none of the five is resolved beyond the 3rd decimal place
  testthat::expect_true(all(se > 0.003))
  testthat::expect_true(all(se < 0.02))
})

# ---- integrity ------------------------------------------------------------

testthat::test_that("this adjudication changed no payload", {
  mp <- path_results("manuscript", "figure_export_manifest.csv")
  testthat::skip_if_not(file.exists(mp), "figure manifest absent")
  testthat::expect_identical(
    unname(tools::sha256sum(mp)),
    "0fd0c9ed9febbc05ed6928b7fc6bfffdab845daba6255d8b47c607dee5d3c02c")
  fz <- repo_path("docs", "publication_freeze_manifest.yml")
  testthat::skip_if_not(file.exists(fz), "freeze absent")
  testthat::expect_identical(
    unname(tools::sha256sum(fz)),
    "b4d37250360e2e07136ccdfbca64bd0a629e22946ced1f5ffd1fb730b445c49a")
  # Both decisions are now implemented (Phase 6H.10), so these two assertions
  # were inverted. Their detailed gates live in
  # tests/testthat/test-comparego-tail-archival.R; what is checked here is only
  # that the adjudication did not drag the publication payload with it.
  testthat::expect_false(file.exists(WORKBOOK),
    info = "the workbook is back in pride_submission - P2 has been reverted")
  # The line count was a proxy for "the tail is still gone". It is checked
  # directly instead, because the script is live code and Phase 6I.2 had to
  # edit it: it read the clusterProfiler manifest raw and hard-stopped on the
  # unmounted P:// root the manifest records. The byte gates for the split now
  # live on the archival commit, in test-comparego-tail-archival.R.
  ex <- parse(SCRIPT, keep.source = TRUE)
  called <- utils::getParseData(ex)
  called <- called$text[called$token == "SYMBOL_FUNCTION_CALL"]
  testthat::expect_false("slice_sample" %in% called)
  # the tail is gone if no code survives past the final exit
  lines <- readLines(SCRIPT, warn = FALSE)
  testthat::expect_identical(
    as.integer(attr(ex, "srcref")[[length(ex)]])[3],
    as.integer(max(which(nzchar(trimws(lines))))))
})
