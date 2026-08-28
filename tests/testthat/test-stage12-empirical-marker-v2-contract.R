testthat::local_edition(3)

# Stage 12 consumes the empirical ROI marker panels DESCRIPTIVELY
# (interpretation_scope = "ROI/cellular context only; not evidence of
# cell-intrinsic process"), so the canonical identity is `mapped_gene_symbol`.
# The producer sets GeneSymbol = if_else(claim_allowed, official_gene_symbol,
# NA), making GeneSymbol claim-gated; keying descriptive panels on it silently
# deletes descriptive-only rows. This mirrors the ec958ec regression contract
# already enforced in R/wgcna_downstream_utils.R.
#
# The reader lives inside a self-executing audit script, so it is extracted
# from the parsed source rather than sourced.

AUDIT_SCRIPT <- testthat::test_path(
  "..", "..", "06_modules_WGCNA", "12_microglia_wgcna_nature_readiness_audit.R"
)
CONTRACT <- "empirical_roi_marker_v2_animal_paired_limma"

stage12_marker_env <- function() {
  exprs <- parse(AUDIT_SCRIPT, keep.source = FALSE)
  env <- new.env(parent = globalenv())
  wanted <- c(
    "EMPIRICAL_MARKER_CONTRACT", "EMPIRICAL_MARKER_IDENTITY_FIELD",
    "read_empirical_roi_markers_v2"
  )
  for (e in exprs) {
    parts <- as.list(e)
    if (length(parts) == 3L && identical(as.character(parts[[1]]), "<-") &&
        is.name(parts[[2]]) && as.character(parts[[2]]) %in% wanted) {
      assign(as.character(parts[[2]]), eval(parts[[3]], envir = env), envir = env)
    }
  }
  env
}

write_fixture <- function(df) {
  # Plain tempfile(): withr::local_tempfile() would delete the file when this
  # helper returns, before the reader under test can open it.
  path <- tempfile(fileext = ".csv")
  utils::write.csv(df, path, row.names = FALSE, na = "")
  path
}

# Parallels the ec958ec fixture: a descriptive-only shared panel whose
# GeneSymbol is claim-gated to NA, an enriched claim-allowed panel, and a
# mapping-blocked row.
marker_fixture <- function(contract = CONTRACT) {
  rbind(
    data.frame(
      marker_set = "empirical_microglia_neuropil_shared",
      GeneSymbol = NA_character_,
      mapped_gene_symbol = c("Shared1", "Shared2"),
      marker_evidence_type = "descriptive_only",
      claim_allowed = FALSE,
      marker_contract_version = contract,
      stringsAsFactors = FALSE
    ),
    data.frame(
      marker_set = "empirical_microglia_roi_enriched",
      GeneSymbol = "Enr1",
      mapped_gene_symbol = "Enr1",
      marker_evidence_type = "inferential",
      claim_allowed = TRUE,
      marker_contract_version = contract,
      stringsAsFactors = FALSE
    ),
    data.frame(
      marker_set = "empirical_soma_enriched",
      GeneSymbol = NA_character_,
      mapped_gene_symbol = NA_character_,
      marker_evidence_type = "descriptive_only",
      claim_allowed = FALSE,
      marker_contract_version = contract,
      stringsAsFactors = FALSE
    )
  )
}

testthat::test_that("Stage 12 declares the v2 contract and mapped identity", {
  env <- stage12_marker_env()
  testthat::expect_identical(env$EMPIRICAL_MARKER_CONTRACT, CONTRACT)
  testthat::expect_identical(
    env$EMPIRICAL_MARKER_IDENTITY_FIELD, "mapped_gene_symbol"
  )
  testthat::expect_true(is.function(env$read_empirical_roi_markers_v2))
})

testthat::test_that("1/2. the current v2 contract is accepted and keyed on mapped_gene_symbol", {
  env <- stage12_marker_env()
  out <- env$read_empirical_roi_markers_v2(write_fixture(marker_fixture()))
  testthat::expect_true("descriptive_gene_symbol" %in% names(out))
  testthat::expect_identical(out$descriptive_gene_symbol, out$mapped_gene_symbol)
})

testthat::test_that("3. a descriptive-only row survives when GeneSymbol is NA", {
  env <- stage12_marker_env()
  out <- env$read_empirical_roi_markers_v2(write_fixture(marker_fixture()))
  shared <- out[out$marker_set == "empirical_microglia_neuropil_shared", ]
  testthat::expect_equal(nrow(shared), 2L)
  # Claim-gated: no usable GeneSymbol (NA on read, "" after a CSV round-trip).
  testthat::expect_true(all(is.na(shared$GeneSymbol) | !nzchar(shared$GeneSymbol)))
  testthat::expect_false(any(shared$claim_allowed %in% TRUE))
  testthat::expect_setequal(shared$descriptive_gene_symbol, c("Shared1", "Shared2"))
})

testthat::test_that("4. a mapping-blocked row is excluded and accounted for", {
  env <- stage12_marker_env()
  out <- env$read_empirical_roi_markers_v2(write_fixture(marker_fixture()))
  testthat::expect_false("empirical_soma_enriched" %in% out$marker_set)
  testthat::expect_identical(attr(out, "rows_mapping_blocked"), 1L)
  testthat::expect_equal(nrow(out), 3L)
})

testthat::test_that("5. a missing mapped_gene_symbol column fails closed", {
  env <- stage12_marker_env()
  fx <- marker_fixture(); fx$mapped_gene_symbol <- NULL
  testthat::expect_error(
    env$read_empirical_roi_markers_v2(write_fixture(fx)),
    "missing required column"
  )
  for (col in c("marker_set", "marker_evidence_type", "claim_allowed",
                "marker_contract_version")) {
    fx2 <- marker_fixture(); fx2[[col]] <- NULL
    testthat::expect_error(
      env$read_empirical_roi_markers_v2(write_fixture(fx2)),
      "missing required column", info = col
    )
  }
})

testthat::test_that("6. a non-v2 marker_contract_version fails closed", {
  env <- stage12_marker_env()
  testthat::expect_error(
    env$read_empirical_roi_markers_v2(
      write_fixture(marker_fixture("empirical_roi_marker_v1_legacy"))
    ),
    "contract mismatch"
  )
  mixed <- marker_fixture()
  mixed$marker_contract_version[1] <- "empirical_roi_marker_v1_legacy"
  testthat::expect_error(
    env$read_empirical_roi_markers_v2(write_fixture(mixed)), "contract mismatch"
  )
})

testthat::test_that("7. there is no silent GeneSymbol fallback or row-wise coalesce", {
  env <- stage12_marker_env()
  fx <- marker_fixture()
  # A row whose GeneSymbol is populated but whose mapped identity is blocked
  # must be dropped, never promoted via GeneSymbol.
  fx$GeneSymbol[3] <- "BlockedButClaimed"
  out <- env$read_empirical_roi_markers_v2(write_fixture(fx))
  testthat::expect_false("BLOCKEDBUTCLAIMED" %in% toupper(out$descriptive_gene_symbol))
  testthat::expect_false("BlockedButClaimed" %in% out$descriptive_gene_symbol)
  testthat::expect_identical(attr(out, "rows_mapping_blocked"), 1L)
  # And the descriptive identity never borrows from the claim-gated field.
  testthat::expect_identical(out$descriptive_gene_symbol, out$mapped_gene_symbol)
})

testthat::test_that("8. Stage-12 panel grepl selection semantics are unchanged", {
  src <- readLines(AUDIT_SCRIPT, warn = FALSE)
  for (pat in c(
    'empirical_microglia_ROI = unique(toupper(empirical_markers$descriptive_gene_symbol[grepl("microglia", empirical_markers$marker_set, ignore.case = TRUE)]))',
    'empirical_neuropil = unique(toupper(empirical_markers$descriptive_gene_symbol[grepl("neuropil", empirical_markers$marker_set, ignore.case = TRUE)]))',
    'empirical_soma = unique(toupper(empirical_markers$descriptive_gene_symbol[grepl("soma", empirical_markers$marker_set, ignore.case = TRUE)]))'
  )) {
    testthat::expect_true(any(grepl(pat, src, fixed = TRUE)), info = pat)
  }
  # The claim-gated field must no longer drive any empirical panel.
  panel_lines <- grep("empirical_(microglia_ROI|neuropil|soma) =", src, value = TRUE)
  testthat::expect_gt(length(panel_lines), 0L)
  testthat::expect_false(any(grepl("empirical_markers$GeneSymbol", panel_lines, fixed = TRUE)))
})

testthat::test_that("9. inferential and statistical thresholds are unchanged", {
  src <- paste(readLines(AUDIT_SCRIPT, warn = FALSE), collapse = "\n")
  # marker significance rule used by the context classifier
  testthat::expect_true(grepl("me$FDR <= 0.10 & me$overlap_count >= 2", src, fixed = TRUE))
  # hypergeometric test and BH family retained
  testthat::expect_true(grepl("stats::phyper(", src, fixed = TRUE))
  testthat::expect_true(grepl("FUN = p.adjust", src, fixed = TRUE))
  testthat::expect_true(grepl('method = "BH"', src, fixed = TRUE))
  # context classification thresholds retained
  testthat::expect_true(grepl("micro_vs_neuropil > 0.25", src, fixed = TRUE))
  testthat::expect_true(grepl("micro_vs_neuropil < -0.25", src, fixed = TRUE))
})
