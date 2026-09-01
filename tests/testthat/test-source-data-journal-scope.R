source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "export_helpers.R"))

repo <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)
exporter <- file.path(repo, "09_export_pride_journal", "09_export_source_data.R")

# Build an isolated results tree so the scope predicates are exercised against
# fixtures rather than the live analysis outputs.
make_scope_fixture <- function() {
  root <- withr::local_tempdir("sd_scope_", .local_envir = parent.frame())
  results <- file.path(root, "results")
  touch <- function(rel) {
    p <- file.path(root, rel)
    dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    writeLines("x", p)
    p
  }
  list(root = root, results = results, touch = touch)
}

reason_of <- function(paths, results_root) {
  source_data_scope_exclusion_reasons(paths, results_root)
}

# ---------------------------------------------------------------------------
# (1) DIAGNOSTIC
# ---------------------------------------------------------------------------
testthat::test_that("per-term GSEA direction audits are excluded as diagnostic", {
  fx <- make_scope_fixture()
  excluded <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/",
    "01b_gsea_protein_direction_audit/neuron_neuropil/BP/gsea_term_direction_audit.csv"
  ))
  # Same family, different dataset -> also excluded.
  excluded2 <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/",
    "01b_gsea_protein_direction_audit/microglia/BP/gsea_term_direction_audit.csv"
  ))
  # A neighbouring summary in the same directory is NOT excluded: the rule names
  # one file, not the directory.
  kept <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/",
    "01b_gsea_protein_direction_audit/neuron_neuropil/BP/direction_summary.csv"
  ))

  testthat::expect_identical(reason_of(excluded, fx$results), "diagnostic")
  testthat::expect_identical(reason_of(excluded2, fx$results), "diagnostic")
  testthat::expect_true(is.na(reason_of(kept, fx$results)))
})

# ---------------------------------------------------------------------------
# (2) SUPERSEDED REMOVED PIPELINE
# ---------------------------------------------------------------------------
testthat::test_that("the removed neuropil_contamination_annotation tree is excluded", {
  fx <- make_scope_fixture()
  excluded <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/",
    "neuropil_contamination_annotation/microglia/microglia_neuropil_annotation_20260601_104038.csv"
  ))
  # The replacement tree must NOT be caught by this rule.
  kept <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/",
    "neuropil_reference_annotation/microglia/microglia_neuropil_annotation_latest.csv"
  ))

  testthat::expect_identical(reason_of(excluded, fx$results), "superseded_removed_pipeline")
  testthat::expect_true(is.na(reason_of(kept, fx$results)))
})

# ---------------------------------------------------------------------------
# (3) SUPERSEDED TIMESTAMPED VARIANTS (canonical pointer required)
# ---------------------------------------------------------------------------
testthat::test_that("timestamped neuropil_reference variants are excluded only beside a _latest pointer", {
  fx <- make_scope_fixture()
  dir_rel <- paste0(
    "results/tables/04_differential_expression_enrichment/",
    "neuropil_reference_annotation/microglia/"
  )
  latest <- fx$touch(paste0(dir_rel, "microglia_neuropil_annotation_latest.csv"))
  latest_sum <- fx$touch(paste0(dir_rel, "microglia_neuropil_annotation_summary_latest.csv"))
  ts <- fx$touch(paste0(dir_rel, "microglia_neuropil_annotation_20260717_221250.csv"))
  ts_sum <- fx$touch(paste0(dir_rel, "microglia_neuropil_annotation_summary_20260824_154023.csv"))
  # An unrelated timestamped file elsewhere must NOT be swept up: there is no
  # generic "timestamped means superseded" rule.
  other_ts <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/",
    "some_other_family/microglia/some_other_output_20260601_104038.csv"
  ))

  testthat::expect_identical(reason_of(ts, fx$results), "superseded_timestamped_variant")
  testthat::expect_identical(reason_of(ts_sum, fx$results), "superseded_timestamped_variant")
  testthat::expect_true(is.na(reason_of(latest, fx$results)))
  testthat::expect_true(is.na(reason_of(latest_sum, fx$results)))
  testthat::expect_true(is.na(reason_of(other_ts, fx$results)))
})

testthat::test_that("a timestamped variant without its _latest pointer is kept, not dropped", {
  fx <- make_scope_fixture()
  orphan_ts <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/",
    "neuropil_reference_annotation/microglia/microglia_neuropil_annotation_20260717_221250.csv"
  ))
  # No *_latest.csv beside it -> fail closed by keeping what may be the only copy.
  testthat::expect_true(is.na(reason_of(orphan_ts, fx$results)))
})

testthat::test_that("the two stems are judged independently by their own _latest pointer", {
  # Mirrors the real tree: microglia_neuropil_annotation has a _latest pointer,
  # microglia_neuropil_annotation_summary does not. The summary series must
  # therefore be retained while the annotation series is excluded.
  fx <- make_scope_fixture()
  dir_rel <- paste0(
    "results/tables/04_differential_expression_enrichment/",
    "neuropil_reference_annotation/microglia/"
  )
  fx$touch(paste0(dir_rel, "microglia_neuropil_annotation_latest.csv"))
  ann_ts <- fx$touch(paste0(dir_rel, "microglia_neuropil_annotation_20260824_171702.csv"))
  sum_ts <- fx$touch(paste0(dir_rel, "microglia_neuropil_annotation_summary_20260824_171702.csv"))

  testthat::expect_identical(reason_of(ann_ts, fx$results), "superseded_timestamped_variant")
  testthat::expect_true(is.na(reason_of(sum_ts, fx$results)))
})

# ---------------------------------------------------------------------------
# (4) PROPOSED / NON-CANONICAL VARIANTS (canonical sibling tree required)
# ---------------------------------------------------------------------------
testthat::test_that("proposed variant trees are excluded when the canonical tree exists", {
  fx <- make_scope_fixture()
  base <- "results/tables/04_differential_expression_enrichment/"
  fx$touch(paste0(base, "compareGO_spatial_atlas/spatial_atlas_enrichment_long.csv"))
  proposed <- fx$touch(paste0(
    base, "compareGO_spatial_atlas_validation_proposed/spatial_atlas_enrichment_long.csv"
  ))
  fx$touch(paste0(base, "microglia_targeted_signature_enrichment/microglia/microglia_signature_enrichment.csv"))
  proposed2 <- fx$touch(paste0(
    base, "microglia_targeted_signature_enrichment/microglia_validation_proposed/microglia_signature_enrichment.csv"
  ))

  testthat::expect_identical(reason_of(proposed, fx$results), "proposed_non_canonical")
  testthat::expect_identical(reason_of(proposed2, fx$results), "proposed_non_canonical")
})

testthat::test_that("a proposed tree without its canonical sibling fails closed", {
  fx <- make_scope_fixture()
  proposed <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/",
    "compareGO_spatial_atlas_validation_proposed/spatial_atlas_enrichment_long.csv"
  ))
  testthat::expect_error(
    reason_of(proposed, fx$results),
    "canonical"
  )
})

testthat::test_that("the canonical sibling tree itself is always kept", {
  fx <- make_scope_fixture()
  base <- "results/tables/04_differential_expression_enrichment/"
  canonical <- fx$touch(paste0(base, "compareGO_spatial_atlas/spatial_atlas_enrichment_long.csv"))
  fx$touch(paste0(base, "compareGO_spatial_atlas_validation_proposed/x.csv"))
  testthat::expect_true(is.na(reason_of(canonical, fx$results)))
})

# ---------------------------------------------------------------------------
# (5) INTERMEDIATE TERM-GENE EXPANSIONS
# ---------------------------------------------------------------------------
testthat::test_that("term-gene provenance and evidence expansions are excluded", {
  fx <- make_scope_fixture()
  prov <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/compareGO/neuron_neuropil/BP/",
    "phenotype_within_unit/all_route_units/compareGO_term_gene_provenance.csv"
  ))
  evid <- fx$touch(paste0(
    "results/source_data/04_differential_expression_enrichment/",
    "biological_program_summary/neuron_soma/program_term_gene_evidence.csv"
  ))
  # A compareGO provenance file outside all_route_units is NOT excluded.
  kept <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment/compareGO/neuron_neuropil/BP/",
    "phenotype_within_unit/CA1_neuropil/compareGO_term_gene_provenance.csv"
  ))
  # Other biological_program_summary outputs are kept.
  kept2 <- fx$touch(paste0(
    "results/source_data/04_differential_expression_enrichment/",
    "biological_program_summary/neuron_soma/program_summary.csv"
  ))

  testthat::expect_identical(reason_of(prov, fx$results), "intermediate_term_gene_expansion")
  testthat::expect_identical(reason_of(evid, fx$results), "intermediate_term_gene_expansion")
  testthat::expect_true(is.na(reason_of(kept, fx$results)))
  testthat::expect_true(is.na(reason_of(kept2, fx$results)))
})

# ---------------------------------------------------------------------------
# (6) STAGE-11 INTERMEDIATE
# ---------------------------------------------------------------------------
testthat::test_that("the Stage-11 joint control-pair GO FDR expansion is excluded, summaries kept", {
  fx <- make_scope_fixture()
  base <- paste0(
    "results/source_data/04_differential_expression_enrichment/",
    "stress_response_biological_audit/global/"
  )
  excluded <- fx$touch(paste0(base, "control_pair_joint_GO_FDR.csv"))
  keeps <- c(
    fx$touch(paste0(base, "direct_SUS_RES_DAP_control_geometry.csv")),
    fx$touch(paste0(base, "direct_SUS_RES_DAP_control_geometry_summary_global.csv")),
    fx$touch(paste0(base, "direct_SUS_RES_DAP_control_geometry_summary_by_context.csv")),
    fx$touch(paste0(base, "theme_trajectories.csv")),
    fx$touch(paste0(base, "protein_contrast_algebra_audit.csv"))
  )
  testthat::expect_identical(reason_of(excluded, fx$results), "intermediate_stage11")
  testthat::expect_true(all(is.na(reason_of(keeps, fx$results))))
})

# ---------------------------------------------------------------------------
# (7) CONCORDANCE INTERMEDIATE
# ---------------------------------------------------------------------------
testthat::test_that("the all-contrasts theme-assignment expansion is excluded in both copies", {
  fx <- make_scope_fixture()
  a <- fx$touch(paste0(
    "results/tables/10_biological_integration/gsea_wgcna_concordance/global/",
    "ontology_aware_gsea_theme_assignments_all_contrasts.csv"
  ))
  b <- fx$touch(paste0(
    "results/source_data/10_biological_integration/gsea_wgcna_concordance/global/",
    "ontology_aware_gsea_theme_assignments_all_contrasts.csv"
  ))
  kept <- fx$touch(paste0(
    "results/source_data/10_biological_integration/gsea_wgcna_concordance/global/",
    "integration__gsea_wgcna_concordance_long.csv"
  ))
  testthat::expect_identical(reason_of(a, fx$results), "intermediate_concordance")
  testthat::expect_identical(reason_of(b, fx$results), "intermediate_concordance")
  testthat::expect_true(is.na(reason_of(kept, fx$results)))
})

# ---------------------------------------------------------------------------
# Canonical keep coverage + absence of broad heuristics
# ---------------------------------------------------------------------------
testthat::test_that("required manuscript-facing source families are never excluded", {
  fx <- make_scope_fixture()
  required <- c(
    fx$touch(paste0("results/source_data/04_differential_expression_enrichment/",
                    "control_spatial_identity_validation/global/figure2f_source_data.csv")),
    fx$touch(paste0("results/source_data/04_differential_expression_enrichment/",
                    "control_spatial_identity_validation/global/figure2f_regions_CA1layers_source_data.csv")),
    fx$touch("results/source_data/manuscript_panels/figure_3/figure3b_stage07_effect_source.csv"),
    fx$touch("results/tables/manuscript_panels/figure_3/figure3_validation.csv"),
    fx$touch(paste0("results/source_data/04_differential_expression_enrichment/",
                    "sus_res_spatial_dap_atlas/global/sus_res_dap_membership.csv")),
    fx$touch("results/tables/06_modules_WGCNA/group_effects/microglia/module_group_effects.csv"),
    fx$touch(paste0("results/tables/06_modules_WGCNA/01_WGCNA/microglia/supermodules/",
                    "wgcna_module_supermodule_annotation.csv")),
    fx$touch(paste0("results/source_data/06_modules_WGCNA/identity_contract/microglia/",
                    "WGCNA_identity_contract_status.csv")),
    fx$touch("results/tables/biological_claims_table.csv")
  )
  testthat::expect_true(all(is.na(reason_of(required, fx$results))))
  testthat::expect_identical(
    apply_manuscript_source_data_scope(required, fx$results), required
  )
})

testthat::test_that("the scope uses no size, timestamp-pattern or name heuristics", {
  src <- readLines(file.path(repo, "R", "export_helpers.R"), warn = FALSE)
  block_start <- grep("Manuscript / journal source-data scope", src, fixed = TRUE)
  testthat::expect_length(block_start, 1L)
  block_end <- grep("^is_orphan_figure_family_path <- function", src)
  block <- src[block_start[[1]]:block_end[[1]]]
  code <- block[!grepl("^\\s*#", block)]

  # No size-based rule.
  testthat::expect_false(any(grepl("file.info", code, fixed = TRUE)))
  testthat::expect_false(any(grepl("$size", code, fixed = TRUE)))
  # No "proposed appears anywhere in the name" rule: the proposed exclusion is a
  # table of explicit (proposed, canonical) prefix pairs.
  testthat::expect_true(any(grepl("source_data_proposed_tree_pairs", code, fixed = TRUE)))
  # Both fail-closed behaviours documented: keep-on-missing-pointer for the
  # timestamped stems, hard error for the proposed trees.
  testthat::expect_true(any(grepl("Fail closed by KEEPING", block, fixed = TRUE)))
  testthat::expect_true(any(grepl("without its canonical", block, fixed = TRUE)))
})

testthat::test_that("empty and no-op inputs are handled", {
  fx <- make_scope_fixture()
  testthat::expect_identical(
    apply_manuscript_source_data_scope(character(0), fx$results), character(0)
  )
  testthat::expect_length(source_data_scope_exclusion_reasons(character(0), fx$results), 0L)
})

# ---------------------------------------------------------------------------
# Wiring, dry-run purity and PRIDE isolation
# ---------------------------------------------------------------------------
testthat::test_that("the source-data exporter applies the scope before the dry-run guard", {
  src <- readLines(exporter, warn = FALSE)
  code <- src[!grepl("^\\s*#", src)]
  scope_at <- grep("source_data_scope_exclusion_reasons(", code, fixed = TRUE)
  guard_at <- grep("isTRUE(dry_run)", code, fixed = TRUE)
  testthat::expect_length(scope_at, 1L)
  testthat::expect_length(guard_at, 1L)
  testthat::expect_lt(scope_at[[1]], guard_at[[1]])

  # Dry-run purity preserved: nothing mutating before the guard.
  mutating <- grep(
    paste("dir_create", "dir\\.create", "file\\.copy", "write\\.csv",
          "write\\.table", "writeLines", "write_run_manifest", sep = "|"),
    code
  )
  testthat::expect_identical(mutating[mutating < guard_at[[1]]], integer(0))

  # The manuscript/results exclusion guard is still applied.
  testthat::expect_true(any(grepl("is_exportable_result_path(candidates)", code, fixed = TRUE)))
})

testthat::test_that("the journal scope is not applied to PRIDE selectors", {
  pride <- c("R/pride_helpers.R",
             "09_export_pride_journal/05_make_pride_manifest.R",
             "09_export_pride_journal/10_validate_pride_submission.R",
             "09_export_pride_journal/03_export_processed_pg_matrix_package.R",
             "09_export_pride_journal/04_make_supplementary_tables.R")
  for (p in pride) {
    joined <- paste(readLines(file.path(repo, p), warn = FALSE), collapse = "\n")
    testthat::expect_false(
      grepl("source_data_scope_exclusion_reasons", joined, fixed = TRUE),
      info = p
    )
    testthat::expect_false(
      grepl("apply_manuscript_source_data_scope", joined, fixed = TRUE),
      info = p
    )
  }
})

# ---------------------------------------------------------------------------
# (8) INTERMEDIATE LEGACY REPLAY: one exact prefix, siblings retained.
# ---------------------------------------------------------------------------
testthat::test_that("the exact legacy_replay prefix is excluded", {
  fx <- make_scope_fixture()
  base <- "results/source_data/04_differential_expression_enrichment_comparison/"
  excluded <- c(
    fx$touch(paste0(base, "legacy_replay/clusterProfiler/microglia/KEGG/a.csv")),
    fx$touch(paste0(base, "legacy_replay/summary.csv"))
  )
  testthat::expect_identical(
    reason_of(excluded, fx$results),
    rep("intermediate_legacy_replay", length(excluded))
  )
})

testthat::test_that("sibling comparison trees are retained, not swept up", {
  fx <- make_scope_fixture()
  base <- "results/source_data/04_differential_expression_enrichment_comparison/"
  kept <- c(
    fx$touch(paste0(base, "animal_level/clusterProfiler/microglia/KEGG/a.csv")),
    fx$touch(paste0(base, "repro1/clusterProfiler/neuron_soma/BP/b.csv")),
    fx$touch(paste0(base, "repro2/clusterProfiler/neuron_soma/BP/c.csv")),
    fx$touch(paste0(base, "rnga/clusterProfiler/neuron_soma/BP/d.csv")),
    fx$touch(paste0(base, "rngb/clusterProfiler/neuron_soma/BP/e.csv")),
    fx$touch(paste0(base, "results/overview.csv"))
  )
  testthat::expect_true(all(is.na(reason_of(kept, fx$results))))

  # No "legacy in the name" rule: a legacy-named file outside the exact prefix
  # stays selected.
  legacy_named <- fx$touch(paste0(
    "results/tables/01_preprocessing/03c_legacy_vs_animal_level_da_audit/",
    "neuron_soma/protein_level_comparison.csv"
  ))
  testthat::expect_true(is.na(reason_of(legacy_named, fx$results)))

  # The tables-root counterpart of the prefix is not excluded either: the rule
  # names the source_data prefix exactly.
  tables_side <- fx$touch(paste0(
    "results/tables/04_differential_expression_enrichment_comparison/",
    "legacy_replay/clusterProfiler/microglia/KEGG/a.csv"
  ))
  testthat::expect_true(is.na(reason_of(tables_side, fx$results)))
})

testthat::test_that("the nine 260-character legacy_replay sources are unselected", {
  fx <- make_scope_fixture()
  base <- paste0(
    "results/source_data/04_differential_expression_enrichment_comparison/",
    "legacy_replay/clusterProfiler/microglia/KEGG/phenotype_within_unit/"
  )
  paths <- character(0)
  for (region in c("CA1", "CA2", "CA3")) {
    for (contrast in c("res_%scon", "sus_%scon", "sus_%sres")) {
      nm <- sprintf(
        "%smicroglia%s_KEGG.csv", region,
        sub("%s", paste0(region, "microglia"), sprintf(contrast, ""), fixed = TRUE)
      )
      paths <- c(paths, fx$touch(paste0(base, region, "_microglia/", nm)))
    }
  }
  testthat::expect_length(paths, 9L)
  testthat::expect_identical(
    reason_of(paths, fx$results),
    rep("intermediate_legacy_replay", 9L)
  )
  testthat::expect_length(apply_manuscript_source_data_scope(paths, fx$results), 0L)
})

testthat::test_that("the legacy_replay rule is a single literal prefix", {
  src <- readLines(file.path(repo, "R", "export_helpers.R"), warn = FALSE)
  fn <- grep("^source_data_excluded_legacy_replay <- function", src)
  testthat::expect_length(fn, 1L)
  body <- src[fn[[1]]:(fn[[1]] + 2L)]
  testthat::expect_true(any(grepl("startsWith(rel, source_data_legacy_replay_prefix())",
                                  body, fixed = TRUE)))
  # No regex, no grepl, no "legacy" pattern matching in the predicate.
  testthat::expect_false(any(grepl("grepl", body, fixed = TRUE)))
})
