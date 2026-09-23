# The active contracts must describe what the repository actually produces.
#
# Two contracts had drifted, in opposite directions, and both were invisible to
# the tests that were supposed to cover them.
#
#   1. inst/schemas/biological_claims_table.yml declared 70 columns for a table
#      that has 75. The five it omitted are the WGCNA display label and its
#      provenance, attached by attach_canonical_wgcna_display_label(). The
#      producer keeps the superseded stage labels ON PURPOSE - "No historical
#      label is discarded" (R/statistics/wgcna_label_activation_utils.R:406) -
#      and test-wgcna-label-activation.R already asserts they are never
#      dropped. So the output was right and the contract was stale: the fix was
#      to declare them, never to stop writing them.
#
#   2. docs/file_contracts.tsv declared two clusterProfiler objects under a
#      results/ root that exists but holds zero files, while the 54 executed
#      per-comparison directories sit under data/processed/. Nothing noticed,
#      because no validator checks that a declared path resolves to anything.
#
# These tests close that gap: they check agreement between contract and
# reality, not the shape of today's data.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "utilities", "schema_validation.R"))

CLAIMS <- path_results("tables", "biological_claims_table.csv")
CONTRACTS <- repo_path("docs", "file_contracts.tsv")

read_schema_columns <- function(name) {
  y <- yaml::read_yaml(repo_path("inst", "schemas", paste0(name, ".yml")))
  list(required = unlist(y$required_columns %||% character(), use.names = FALSE),
       declared = names(y$columns %||% list()))
}

testthat::test_that("the claims table validates against its own schema", {
  testthat::skip_if_not(file.exists(CLAIMS), "claims table absent")
  d <- utils::read.csv(CLAIMS, stringsAsFactors = FALSE, check.names = FALSE, nrows = 5)
  testthat::expect_silent(validate_table_schema(d, "biological_claims_table", strict = TRUE))
})

testthat::test_that("no active field of the claims table is undeclared", {
  # The failure mode that went unnoticed for months: the producer grows a
  # column and the schema does not hear about it.
  testthat::skip_if_not(file.exists(CLAIMS), "claims table absent")
  obs <- names(utils::read.csv(CLAIMS, nrows = 1, stringsAsFactors = FALSE,
                               check.names = FALSE))
  s <- read_schema_columns("biological_claims_table")
  testthat::expect_identical(setdiff(obs, s$declared), character(0),
    info = paste("columns present but undeclared:",
                 paste(setdiff(obs, s$declared), collapse = ", ")))
  testthat::expect_identical(setdiff(s$required, obs), character(0),
    info = paste("declared required but absent:",
                 paste(setdiff(s$required, obs), collapse = ", ")))
})

testthat::test_that("WGCNA label provenance stays declared, and stays optional", {
  # Two distinct guards. Declared, so a future producer change cannot silently
  # drop the superseded stage labels the submission artifact is supposed to
  # record. Optional, because attach_canonical_wgcna_display_label() returns
  # early when a claims table has no WGCNA entity rows, so requiring them would
  # be a false statement about a legitimate table.
  s <- read_schema_columns("biological_claims_table")
  provenance <- c("canonical_display_label", "canonical_label_source",
                  "Stage01_ModuleLabel_Final", "Stage06_label", "Stage07_label")
  for (p in provenance) {
    testthat::expect_true(p %in% s$declared, info = paste(p, "is no longer declared"))
    testthat::expect_false(p %in% s$required,
      info = paste(p, "was made mandatory; a non-WGCNA claims table legitimately lacks it"))
  }

  # and the producer still says why it keeps them
  src <- paste(readLines(repo_path("R", "statistics", "wgcna_label_activation_utils.R"),
                         warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("No historical label is discarded", src, fixed = TRUE),
    info = "the retention intent for the stage labels was removed from the producer")
  for (p in provenance) {
    testthat::expect_true(grepl(p, src, fixed = TRUE), info = paste(p, "no longer written"))
  }
})

testthat::test_that("the schema change moved no claim", {
  # A contract edit must not touch data. These are properties of the artifact
  # as Phase 6I.5 found it; the schema was the only thing that changed.
  testthat::skip_if_not(file.exists(CLAIMS), "claims table absent")
  d <- utils::read.csv(CLAIMS, stringsAsFactors = FALSE, check.names = FALSE)
  testthat::expect_identical(ncol(d), 75L)
  testthat::expect_gt(nrow(d), 0L)
  # claim identity is unique, which is the invariant a merge accident would break
  if ("claim_id" %in% names(d)) {
    testthat::expect_identical(anyDuplicated(d$claim_id), 0L)
  }
})

# ---- file contracts -------------------------------------------------------

testthat::test_that("every declared file contract resolves to something that exists", {
  # Nothing checked this before, which is how two rows came to name a root
  # holding zero files. Rows whose producer has not run in a given checkout are
  # exempt by necessity - the test asserts the ROOT is real and, where the
  # producer has run, that the declaration finds its artifacts.
  testthat::skip_if_not(file.exists(CONTRACTS), "file contracts absent")
  d <- utils::read.delim(CONTRACTS, sep = "\t", stringsAsFactors = FALSE,
                         check.names = FALSE)
  testthat::expect_gt(nrow(d), 0L)

  # A path spec is one or more ;-separated alternatives. An alternative either
  # carries a <placeholder> or glob, in which case the fixed prefix before the
  # first such segment must be a real directory, or it is literal and must
  # itself exist. A row resolves if any alternative does. Paths are
  # repo-relative and testthat runs from tests/testthat/, so they are anchored.
  resolves <- function(spec) {
    parts <- trimws(strsplit(spec, ";", fixed = TRUE)[[1]])
    for (p in parts[nzchar(parts)]) {
      segs <- strsplit(p, "/", fixed = TRUE)[[1]]
      wild <- grepl("[<*{]", segs)
      if (!any(wild)) {
        if (file.exists(repo_path(p)) || dir.exists(repo_path(p))) return(TRUE)
        next
      }
      prefix <- paste(segs[seq_len(which(wild)[1] - 1L)], collapse = "/")
      if (nzchar(prefix) && dir.exists(repo_path(prefix))) return(TRUE)
    }
    FALSE
  }
  bad <- d$object_id[!vapply(d$path, resolves, logical(1), USE.NAMES = FALSE)]
  testthat::expect_identical(bad, character(0),
    info = paste("contract rows whose declared location does not exist:",
                 paste(bad, collapse = ", ")))
})

testthat::test_that("the clusterProfiler audit contracts point at the executed artifacts", {
  # Phase 6I.5 adjudication F2 WRONG_PATH: the family is active, the producer
  # generates it, and it is consumed - it simply was not where the contract
  # said. Asserting "more than zero" rather than a file count, because the
  # number of comparisons is legitimately variable; what must hold is that the
  # declaration finds the executed run.
  testthat::skip_if_not(file.exists(CONTRACTS), "file contracts absent")
  d <- utils::read.delim(CONTRACTS, sep = "\t", stringsAsFactors = FALSE,
                         check.names = FALSE)
  for (id in c("clusterProfiler_protein_group_audits", "enrichment_term_gene_provenance")) {
    r <- d[d$object_id == id, , drop = FALSE]
    testthat::expect_identical(nrow(r), 1L, info = id)
    testthat::expect_false(grepl("results/tables/04_differential_expression_enrichment/clusterProfiler",
                                 r$path[1], fixed = TRUE),
      info = paste(id, "points back at the empty results/ root"))
    root <- repo_path(sub("/<dataset>.*$", "", trimws(strsplit(r$path[1], ";", fixed = TRUE)[[1]])[1]))
    testthat::skip_if_not(dir.exists(root), paste("root absent in this checkout:", root))
    # The stronger half. A root-exists check would NOT have caught the original
    # defect: results/tables/04_.../clusterProfiler did exist, with three empty
    # dataset subdirectories and zero files. So for the rows this phase
    # adjudicated, require that the declared tree actually holds the artifact
    # family - not a file count, which varies with the number of comparisons.
    dirs <- list.dirs(root, recursive = TRUE, full.names = FALSE)
    found <- sum(basename(dirs) == "protein_group_audits")
    testthat::expect_true(found > 0L,
      info = paste(id, "declared root contains no protein_group_audits directory"))
  }
})

testthat::test_that("the canonical imputed matrices are declared, with an honest producer", {
  # Phase 6I.7. These are a REQUIRED active input - dataset_inputs.R resolves
  # them with required = TRUE for both the WGCNA and Protigy paths - whose only
  # producer is archived and is not a registered pipeline step. The contract
  # has to say that rather than imply an active producer, so the producer
  # string is asserted, not just the row's existence.
  testthat::skip_if_not(file.exists(CONTRACTS), "file contracts absent")
  d <- utils::read.delim(CONTRACTS, sep = "\t", stringsAsFactors = FALSE,
                         check.names = FALSE)
  for (id in c("imputed_protein_matrix", "imputation_seed_provenance")) {
    r <- d[d$object_id == id, , drop = FALSE]
    testthat::expect_identical(nrow(r), 1L, info = id)
    testthat::expect_true(grepl("archive/01_preprocessing/01_impute.r", r$created_by[1], fixed = TRUE),
      info = paste(id, "no longer names the archived producer"))
    testthat::expect_true(grepl("ARCHIVED", r$created_by[1], fixed = TRUE),
      info = paste(id, "stopped flagging that no active step regenerates these"))
  }
})

testthat::test_that("the declared matrix pattern is the one the resolver actually uses", {
  # The cross-check that matters. A contract naming a pattern the resolver does
  # not use would be decorative: the resolver is what decides which file is the
  # active input, so the two must describe the same object.
  testthat::skip_if_not(file.exists(CONTRACTS), "file contracts absent")
  d <- utils::read.delim(CONTRACTS, sep = "\t", stringsAsFactors = FALSE,
                         check.names = FALSE)
  declared <- d$path[d$object_id == "imputed_protein_matrix"][1]
  testthat::expect_true(grepl("data/processed/01_preprocessing/impute/", declared, fixed = TRUE))
  testthat::expect_true(grepl("pgmatrix_imputed", declared, fixed = TRUE))
  testthat::expect_true(grepl("missing70pct.xlsx", declared, fixed = TRUE))

  src <- paste(readLines(repo_path("R", "data_contracts", "dataset_inputs.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl("pgmatrix_imputed", src, fixed = TRUE),
    info = "the resolver no longer addresses these matrices")
  testthat::expect_true(grepl("missing70pct", src, fixed = TRUE))
  # and it is still a newest-wins resolution, which is why the contract states
  # cardinality as a pattern rather than a file list
  testthat::expect_true(grepl("latest_pattern", src, fixed = TRUE) ||
                          grepl("latest_matching_file", src, fixed = TRUE),
    info = "generation selection stopped being newest-wins; the contract's cardinality wording is now wrong")
})

testthat::test_that("every dataset has at least one imputed matrix, without pinning a count", {
  # Cardinality is one per dataset per dated generation, newest wins. Two
  # generations coexist today and that is legitimate, so the assertion is
  # per-dataset presence rather than a file total, which would break the next
  # time the imputation is rerun.
  root <- repo_path("data", "processed", "01_preprocessing", "impute")
  testthat::skip_if_not(dir.exists(root), "impute directory absent")
  for (ds in c("microglia", "neuron_neuropil", "neuron_soma")) {
    pat <- paste0("^[0-9]{8}_pgmatrix_imputed_", ds, "_[0-9]+samples_missing70pct[.]xlsx$")
    testthat::expect_gt(length(list.files(root, pattern = pat)), 0L)
  }
  # the seed record that makes the recovery path verifiable
  qc <- file.path(root, "imputation_qc.csv")
  testthat::skip_if_not(file.exists(qc), "imputation QC absent")
  q <- utils::read.csv(qc, stringsAsFactors = FALSE)
  testthat::expect_setequal(q$celltype_layer, c("microglia", "neuron_neuropil", "neuron_soma"))
  testthat::expect_true(all(q$base_seed == 42L))
})

testthat::test_that("the audit contract and the retention declaration agree on a location", {
  # The retention declaration added in Phase 6I.3 protects the executed rank
  # record inside protein_group_audits/. If the file contract named a different
  # tree, one of the two would be wrong. This is the cross-check.
  testthat::skip_if_not(file.exists(CONTRACTS), "file contracts absent")
  source(testthat::test_path("..", "..", "R", "enrichment", "enrichment_io.R"))
  d <- utils::read.delim(CONTRACTS, sep = "\t", stringsAsFactors = FALSE,
                         check.names = FALSE)
  contract_root <- normalizePath(repo_path(sub("/<dataset>.*$", "",
                       d$path[d$object_id == "clusterProfiler_protein_group_audits"][1])),
                       winslash = "/", mustWork = FALSE)

  mf <- file.path(repo_path("data", "processed",
                            "04_differential_expression_enrichment", "clusterProfiler"),
                  "microglia", "clusterProfiler_manifest.csv")
  testthat::skip_if_not(file.exists(mf), "microglia manifest absent")
  m <- utils::read.csv(mf, stringsAsFactors = FALSE)
  m <- m[m$result_type == "GSEA_GO" & m$ontology == "BP", , drop = FALSE]
  protected <- clusterprofiler_protected_reference_paths(m)
  testthat::expect_gt(length(protected), 0L)
  testthat::expect_true(all(grepl(contract_root, normalizePath(protected, winslash = "/",
                                                               mustWork = FALSE), fixed = TRUE)),
    info = "the retention declaration and the file contract name different trees")
})
