# The executed GSEA ranked order must not become deletable by inattention.
#
# Phase 6I.2 established that rank_statistic_sensitivity_audit.csv stores the
# ordered vector passed to gseGO() - BYTE_EXACT_STORED_ORDER - and that the
# Figure 3 g/h/i seven-protein selection reproduces exactly from it. Phase 6I.3
# fixed the retention semantics, because three separate things made the file
# look disposable and none of them was a decision about the file:
#
#   * the clusterProfiler manifest has no column for it, so no registry knew
#     it existed;
#   * nothing reads it, so a reader census returns its own write site;
#   * the Phase 6H over-wall inventory recorded its tree as REGENERABLE, which
#     was derived from the root directory, not from what the file is.
#
# These tests assert the SEMANTIC REGISTRY STATUS rather than blacklisting a
# basename: the artifact is declared, the declared role is one the repository
# already uses for hashed-but-unconsumed artifacts, and the declaration
# actually resolves onto the instances that exist. A rename that kept the
# declaration honest would still pass; quietly dropping the declaration, or
# flipping it to cleanup-eligible, would not.
#
# What this file does NOT do is pin the artifact's bytes. Byte and semantic
# stability are tested in test-gsea-rank-provenance.R, against the manifest's
# own n_genes; duplicating a 54-file hash list here would fail the first time
# the enrichment is legitimately rerun, which is not what "protected" means.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "enrichment", "enrichment_io.R"))

CP <- repo_path("data", "processed", "04_differential_expression_enrichment",
                "clusterProfiler")
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")

gsea_bp_manifest <- function(dataset) {
  mf <- file.path(CP, dataset, "clusterProfiler_manifest.csv")
  if (!file.exists(mf)) return(NULL)
  m <- utils::read.csv(mf, stringsAsFactors = FALSE)
  m[m$result_type == "GSEA_GO" & m$ontology == "BP", , drop = FALSE]
}

testthat::test_that("the executed ranked order is a declared artifact, not an unowned one", {
  d <- clusterprofiler_protected_reference_artifacts()
  testthat::expect_s3_class(d, "data.frame")
  testthat::expect_gte(nrow(d), 1L)
  for (col in c("artifact", "role", "records", "producer", "cleanup_eligible", "note")) {
    testthat::expect_true(col %in% names(d), info = col)
  }
  row <- d[d$artifact == "rank_statistic_sensitivity_audit.csv", , drop = FALSE]
  testthat::expect_identical(nrow(row), 1L,
    info = "the executed ranked-order record is no longer declared")
  testthat::expect_identical(row$records[[1]], "BYTE_EXACT_STORED_ORDER")
  testthat::expect_false(row$cleanup_eligible[[1]],
    info = "the executed ranked-order record has been marked cleanup-eligible")
  testthat::expect_true(nzchar(row$producer[[1]]))
  testthat::expect_true(grepl("gseGO", row$note[[1]], fixed = TRUE),
    info = "the declaration no longer says why the file is not reconstructible-and-disposable")
})

testthat::test_that("the declared role is one the repository already uses", {
  # Not a category invented for this artifact. The same role token classifies
  # the Stage-11 SUS-RES workbook, which is likewise hashed and never read.
  testthat::expect_identical(CLUSTERPROFILER_PROTECTED_REFERENCE_ROLE,
                             "protected_reference_not_consumed")
  src <- paste(readLines(testthat::test_path("..", "..", "R", "statistics",
                                             "stress_response_biological_audit_utils.R"),
                         warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("protected_reference_not_consumed", src, fixed = TRUE),
    info = "the shared role token vanished from its original definition site")

  contracts <- paste(readLines(repo_path("docs", "OUTPUT_CONTRACTS.md"), warn = FALSE),
                     collapse = "\n")
  testthat::expect_true(grepl("clusterProfiler executed ranked order", contracts, fixed = TRUE),
    info = "the output contract for the executed ranked order was removed")
  testthat::expect_true(grepl("protected_reference_not_consumed", contracts, fixed = TRUE))
})

testthat::test_that("the predicate answers for the declared artifact and defaults open", {
  # A declaration only protects what it covers. Anything undeclared stays
  # ordinary, so this cannot quietly become a blanket exemption.
  testthat::expect_false(
    clusterprofiler_artifact_is_cleanup_eligible("rank_statistic_sensitivity_audit.csv"))
  testthat::expect_false(
    clusterprofiler_artifact_is_cleanup_eligible(
      file.path("a", "b", "rank_statistic_sensitivity_audit.csv")))
  testthat::expect_true(
    clusterprofiler_artifact_is_cleanup_eligible("per_contrast_aggregate_audit.csv"))
  testthat::expect_true(
    clusterprofiler_artifact_is_cleanup_eligible("ora_universe_audit.csv"))
})

testthat::test_that("the declaration resolves onto the instances that exist", {
  # A registry entry naming a path nothing occupies protects nothing. This is
  # the failure the Phase 6I.3 audit found in docs/file_contracts.tsv, whose
  # clusterProfiler_protein_group_audits row points at a results/ path holding
  # zero files while the 54 real instances sit under data/processed/.
  declared <- character(0)
  for (ds in DATASETS) {
    m <- gsea_bp_manifest(ds)
    testthat::skip_if(is.null(m), paste("manifest absent:", ds))
    declared <- c(declared, clusterprofiler_protected_reference_paths(m))
  }
  testthat::skip_if(!length(declared), "no clusterProfiler manifests in this checkout")
  testthat::expect_identical(length(declared), 54L)

  status <- input_addressability(declared)
  # Not one of them may be absent. Over the wall is fine - that is an
  # addressability state, and the Phase 6H.3 staging contract resolves it.
  testthat::expect_identical(sum(status == INPUT_STATUS_ABSENT), 0L,
    info = paste("declared but absent:",
                 paste(utils::head(declared[status == INPUT_STATUS_ABSENT], 3),
                       collapse = ", ")))
  testthat::expect_identical(sum(status == INPUT_STATUS_PRESENT), 36L)
  testthat::expect_identical(sum(status == INPUT_STATUS_OVER_LIMIT), 18L)
})

testthat::test_that("the declaration is derived from the manifest, not hardcoded", {
  # If the path list were a literal, it would rot silently the moment a
  # comparison was added or renamed. It is built from the manifest's own
  # collapsed_gene_input_file column, so it tracks the run that happened.
  m <- gsea_bp_manifest("microglia")
  testthat::skip_if(is.null(m), "microglia manifest absent")
  testthat::expect_identical(length(clusterprofiler_protected_reference_paths(m)), 12L)

  # a manifest with one row yields one path, in the same directory
  one <- m[1, , drop = FALSE]
  p <- clusterprofiler_protected_reference_paths(one)
  testthat::expect_identical(length(p), 1L)
  testthat::expect_identical(basename(p), "rank_statistic_sensitivity_audit.csv")
  testthat::expect_identical(basename(dirname(p)), "protein_group_audits")

  # and an empty or column-less manifest yields nothing rather than guessing
  testthat::expect_identical(clusterprofiler_protected_reference_paths(m[0, , drop = FALSE]),
                             character(0))
  testthat::expect_identical(
    clusterprofiler_protected_reference_paths(data.frame(x = 1)), character(0))
})

testthat::test_that("the protected set stops at the executed run", {
  # The basename occurs 60 times under data/processed/, not 54. Six of them
  # live in the comparison trees - animal_level and legacy_replay, both
  # neuron_soma/DG_sg, three comparisons each. Those record the ranked order of
  # *a* GSEA execution, not of the one that produced Figure 3, and calling them
  # publication provenance would overstate what they are.
  #
  # So this asserts the boundary in both directions: the declaration reaches
  # every canonical instance, and reaches none of the comparison ones. A future
  # change that widened it to "every file with this basename" would fail here,
  # which is the point - same name, different tree.
  canonical_root <- repo_path("data", "processed",
                              "04_differential_expression_enrichment")
  comparison_root <- repo_path("data", "processed",
                               "04_differential_expression_enrichment_comparison")
  testthat::skip_if_not(dir.exists(canonical_root), "enrichment tree absent")

  declared <- character(0)
  for (ds in DATASETS) {
    m <- gsea_bp_manifest(ds)
    testthat::skip_if(is.null(m), paste("manifest absent:", ds))
    declared <- c(declared, clusterprofiler_protected_reference_paths(m))
  }
  declared <- normalizePath(declared, winslash = "/", mustWork = FALSE)
  testthat::expect_identical(length(declared), 54L)

  # every declared instance is inside the canonical tree
  in_canonical <- startsWith(declared,
                             normalizePath(canonical_root, winslash = "/", mustWork = FALSE))
  testthat::expect_identical(sum(!in_canonical), 0L)

  # and none is inside a comparison tree
  testthat::skip_if_not(dir.exists(comparison_root), "comparison tree absent")
  cmp_norm <- normalizePath(comparison_root, winslash = "/", mustWork = FALSE)
  testthat::expect_identical(sum(startsWith(declared, cmp_norm)), 0L)

  # the comparison instances do exist, so this is an exclusion rather than a
  # vacuous assertion about an empty tree
  cmp_instances <- list.files(comparison_root,
                              pattern = "^rank_statistic_sensitivity_audit[.]csv$",
                              recursive = TRUE)
  testthat::expect_gte(length(cmp_instances), 6L)
  testthat::expect_true(all(grepl("^(animal_level|legacy_replay)/", cmp_instances)),
    info = paste("an unexpected comparison branch appeared:",
                 paste(utils::head(unique(dirname(cmp_instances)), 3), collapse = ", ")))
})

testthat::test_that("frozen Phase 6H snapshots are left as the history they record", {
  # The over-wall inventory says REGENERABLE for this tree. That was true of
  # the classification at the time, and a frozen snapshot that is edited to
  # agree with the present stops being evidence of anything. The current
  # declaration supersedes it prospectively instead.
  inv <- repo_path("audits", "phase6h_live_over_maxpath_inventory.csv")
  testthat::skip_if_not(file.exists(inv), "frozen inventory absent")
  d <- utils::read.csv(inv, stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(d), 724L)
  testthat::expect_setequal(unique(d$lifecycle), c("REGENERABLE", "RESULTS"))
  testthat::expect_gt(sum(grepl("rank_statistic_sensitivity_audit", d$rel, fixed = TRUE)), 0L)
})
