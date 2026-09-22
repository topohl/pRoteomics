# Phase 6H.10 — implementation of decisions C2 and P2.
#
#   C2  archive the unreachable compare_go_enrichment tail, byte-for-byte
#   P2  move the orphan bootstrap workbook out of the outward PRIDE payload
#       while retaining it internally as provenance
#
# These tests protect the two properties that make the change safe: the
# archived code is exactly the bytes that were removed, and the curated
# workbook is exactly the bytes that left the payload. Everything else here
# guards against the change quietly undoing itself.

source(testthat::test_path("..", "..", "R", "paths.R"))

ACTIVE  <- repo_path("analysis", "differential_abundance", "compare_go_enrichment.R")
ARCHIVE <- repo_path("archive", "04_differential_expression_enrichment", "legacy",
                     "02_compareGO_superseded_tail.r")
SIDECAR <- repo_path("archive", "04_differential_expression_enrichment", "legacy",
                     "02_compareGO_superseded_tail.PROVENANCE.md")

# Recorded before the split, in Phase 6H.10's snapshot.
ORIGINAL_SHA <- "dc9e41d98a309ba7b582576984a7a7fc0e83712deeb8bc3d079e3cb91aef9a91"
PREFIX_SHA   <- "e5bce156c78365e2f393a78cdb2aac9ea7a48ec0476a8f0ab9ec1f0d7450f5d7"
TAIL_SHA     <- "cde1bccf1e4f6381b7a40c9cf66389dfb8d22a58b46bf8f548cc33dcd04403ba"

WORKBOOK_SHA <- "024d3671f2cd6026e8bfd7feb6d9839c7ff23eeb41c98d8f66235d854185baba"
WORKBOOK_NAME <- paste0("results_tables_04_differential_expression_enrichment_compareGO_",
                        "neuron_neuropil_BP_phenotype_within_unit_",
                        "08_Bootstrap_Stability_Summary.xlsx")
# results/manuscript/_curated/ rather than a deeper dated directory: the
# workbook's basename is 139 characters, and a longer parent pushed the
# absolute path to 271 - past the 260-character wall this repository hit in
# Phases 6H.1-6H.3, which would have left the retained provenance copy
# unreadable by R. At 237 it is inside both the 260 wall and the conservative
# 240 write guard. The curation date lives in the sidecar instead.
CURATED <- path_results("manuscript", "_curated", WORKBOOK_NAME)
OUTWARD <- repo_path("pride_submission", "supplementary_tables", WORKBOOK_NAME)

# ---- C2: byte preservation -------------------------------------------------

testthat::test_that("the archived tail is byte-identical to what was removed", {
  testthat::skip_if_not(file.exists(ARCHIVE), "archived tail absent")
  testthat::expect_identical(unname(tools::sha256sum(ARCHIVE)), TAIL_SHA)
})

testthat::test_that("the retained active prefix is byte-identical to the original head", {
  # Section 6 of the brief: this is dead-code removal, not a rewrite. The live
  # code must be the SAME BYTES it was before, not merely equivalent.
  testthat::skip_if_not(file.exists(ACTIVE), "active script absent")
  testthat::expect_identical(unname(tools::sha256sum(ACTIVE)), PREFIX_SHA)
})

testthat::test_that("prefix and archive reconstruct the original file exactly", {
  # The strongest gate available: concatenating what was kept with what was
  # archived must reproduce the pre-change file bit for bit. If either side
  # drifts - a reformat, a lint, a line-ending change - this fails.
  testthat::skip_if_not(file.exists(ACTIVE) && file.exists(ARCHIVE),
                        "split artefacts absent")
  tmp <- withr::local_tempfile()
  con <- file(tmp, open = "wb")
  for (f in c(ACTIVE, ARCHIVE)) writeBin(readBin(f, "raw", file.size(f)), con)
  close(con)
  testthat::expect_identical(unname(tools::sha256sum(tmp)), ORIGINAL_SHA)
})

testthat::test_that("the archived payload carries no added commentary", {
  # Section 5: explanatory metadata belongs in the sidecar, so the payload can
  # stay byte-exact. A stray banner added later would break the hash above, but
  # this states the intent explicitly.
  testthat::skip_if_not(file.exists(ARCHIVE), "archived tail absent")
  src <- readLines(ARCHIVE, warn = FALSE)
  testthat::expect_false(any(grepl("Phase 6H", src, fixed = TRUE)))
  testthat::expect_false(any(grepl("ARCHIVED|PROVENANCE ONLY|DO NOT RUN", src)))
  testthat::expect_identical(length(src), 3752L)

  testthat::expect_true(file.exists(SIDECAR))
  side <- readLines(SIDECAR, warn = FALSE)
  for (needle in c(ORIGINAL_SHA, TAIL_SHA, PREFIX_SHA, "NON_RUNNABLE", "C2"))
    testthat::expect_true(any(grepl(needle, side, fixed = TRUE)),
                          info = paste("sidecar missing", needle))
})

# ---- C2: the active surface ------------------------------------------------

testthat::test_that("the active script ends at its real execution boundary", {
  testthat::skip_if_not(file.exists(ACTIVE), "active script absent")
  ex <- parse(ACTIVE, keep.source = TRUE)
  testthat::expect_identical(length(ex), 103L)
  last <- ex[[length(ex)]]
  testthat::expect_true(is.call(last) && identical(as.character(last[[1]])[1], "quit"))
  testthat::expect_identical(length(readLines(ACTIVE, warn = FALSE)), 569L)
})

testthat::test_that("no former-tail material survives in the active script", {
  testthat::skip_if_not(file.exists(ACTIVE), "active script absent")
  pd <- utils::getParseData(parse(ACTIVE, keep.source = TRUE))
  calls <- function(x) sum(pd$token == "SYMBOL_FUNCTION_CALL" & pd$text %in% x)

  # the statistical draw that occupied two audit phases
  testthat::expect_identical(calls("slice_sample"), 0L)
  # the repel layer that could never render
  testthat::expect_identical(calls(c("geom_text_repel", "geom_label_repel")), 0L)
  # the dead writers: 56 of them lived in the tail
  testthat::expect_identical(
    calls(c("write_xlsx", "write_raw_xlsx", "ggsave", "saveRDS", "svg",
            "sink", "saveWorkbook", "write.table", "png", "pdf")), 0L)

  src <- readLines(ACTIVE, warn = FALSE)
  # the obsolete marker, and the two tokens a test used to match in dead code
  for (tok in c("LEGACY_COMPAREGO_TAIL_DISABLED_BY_CANONICAL_EXIT",
                "comparison_input_file", "GeneSymbol"))
    testthat::expect_identical(length(grep(tok, src, fixed = TRUE)), 0L, info = tok)
})

testthat::test_that("the archive is excluded from the active universe by policy", {
  # Not by luck: the repository declares archive/ non-runnable, and the active
  # scanners take their roots from that declaration.
  source(repo_path("R", "utilities", "pipeline_registry.R"))
  testthat::expect_true("archive" %in% pipeline_analysis_script_exclusions()$roots)
  # and the archived file sits under that root
  rel <- sub(paste0("^", repo_path(), "/"), "", gsub("\\\\", "/", ARCHIVE))
  testthat::expect_true(startsWith(rel, "archive/"))
})

# ---- P2: the workbook ------------------------------------------------------

testthat::test_that("the orphan workbook left the outward payload intact", {
  testthat::expect_false(file.exists(OUTWARD),
    info = "the workbook is back in pride_submission - P2 has been undone")
  testthat::expect_true(file.exists(CURATED),
    info = "the curated copy is missing - the bytes must be retained internally")
  testthat::expect_identical(unname(tools::sha256sum(CURATED)), WORKBOOK_SHA)
  testthat::expect_identical(file.size(CURATED), 5141)
})

testthat::test_that("exactly one artifact left the outward payload", {
  # Guards against a curation that quietly took something else with it.
  supp <- repo_path("pride_submission", "supplementary_tables")
  testthat::skip_if_not(dir.exists(supp), "pride supplementary tables absent")
  files <- list.files(supp, recursive = TRUE)
  testthat::expect_identical(length(files), 495L)
  testthat::expect_identical(
    length(grep("Bootstrap_Stability", files, fixed = TRUE)), 0L)
})

testthat::test_that("the outward manifest lost exactly the curated row", {
  man <- repo_path("pride_submission", "manifests", "pride_file_manifest.tsv")
  testthat::skip_if_not(file.exists(man), "pride manifest absent")
  lines <- readLines(man, warn = FALSE)
  testthat::expect_identical(length(lines) - 1L, 1720L)
  testthat::expect_identical(length(grep(WORKBOOK_SHA, lines, fixed = TRUE)), 0L)
  testthat::expect_identical(length(grep("Bootstrap_Stability", lines, fixed = TRUE)), 0L)
})

testthat::test_that("the validator count agrees with the directory", {
  rep <- repo_path("pride_submission", "validation", "validation_report.tsv")
  testthat::skip_if_not(file.exists(rep), "validation report absent")
  line <- grep("supplementary_tables_present", readLines(rep, warn = FALSE), value = TRUE)[1]
  counted <- as.integer(sub(".*\\t([0-9]+) supplementary.*", "\\1", line))
  supp <- repo_path("pride_submission", "supplementary_tables")
  testthat::expect_identical(counted, length(list.files(supp, recursive = TRUE)))
  testthat::expect_identical(counted, 495L)
  # and the curation introduced no validation failure
  rows <- utils::read.delim(rep, stringsAsFactors = FALSE)
  testthat::expect_identical(sum(rows$status == "FAIL"), 0L)
})

testthat::test_that("future packaging cannot re-select the curated workbook", {
  # Section 11: the mechanism must be semantic, not a filename hack. Deposition
  # eligibility already excludes everything under results/manuscript/, which is
  # where the workbook now lives - so the existing contract does the work and no
  # new rule was needed.
  source(repo_path("R", "utilities", "export_helpers.R"))
  testthat::expect_false(is_exportable_result_path(CURATED))
  # the rule is general, not aimed at this file
  testthat::expect_false(is_exportable_result_path(
    repo_path("results", "manuscript", "anything", "else.xlsx")))
  # and it does not over-reach into the live scientific tree
  testthat::expect_true(is_exportable_result_path(
    repo_path("results", "tables", "04_differential_expression_enrichment", "x.csv")))
  # no "exclude anything called bootstrap" rule was introduced
  eh <- readLines(repo_path("R", "utilities", "export_helpers.R"), warn = FALSE)
  testthat::expect_identical(length(grep("Bootstrap_Stability", eh, fixed = TRUE)), 0L)
})

# ---- immutability ----------------------------------------------------------

testthat::test_that("the frozen superseded snapshot was not touched", {
  # Section 12 hard gate. Its stale reference to the workbook is correct
  # provenance for the state it records, and must not be "cleaned".
  f <- path_results("manuscript", "_superseded_20260622", "source_data_export_manifest.csv")
  testthat::skip_if_not(file.exists(f), "superseded manifest absent")
  testthat::expect_identical(
    unname(tools::sha256sum(f)),
    "97b08c92149a6efbcfae32601222d42c5559c55424317fecdb708d9023fe693a")
  # the stale reference is still there, on purpose
  testthat::expect_gt(
    length(grep("Bootstrap_Stability", readLines(f, warn = FALSE), fixed = TRUE)), 0L)
})

testthat::test_that("the other four historical copies are untouched", {
  base <- path_results("manuscript", "_superseded_20260622", "supplementary_tables")
  testthat::skip_if_not(dir.exists(base), "superseded bundle absent")
  expected <- c(
    "08_Bootstrap_Stability_Summary.xlsx" = 15431,
    "04_differential_expression_enrichment_compareGO_microglia_BP_phenotype_within_unit_08_Bootstrap_Stability_Summary.xlsx" = 5135,
    "04_differential_expression_enrichment_compareGO_neuron_soma_BP_phenotype_within_unit_08_Bootstrap_Stability_Summary.xlsx" = 5136,
    "04_differential_expression_enrichment_compareGO_neuron_neuropil_BP_phenotype_within_unit_08_Bootstrap_Stability_Summary.xlsx" = 5143)
  for (nm in names(expected)) {
    f <- file.path(base, nm)
    testthat::skip_if_not(file.exists(f), paste(nm, "absent"))
    testthat::expect_identical(file.size(f), expected[[nm]], info = nm)
  }
})

testthat::test_that("the figure package and freeze are unchanged", {
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
})
