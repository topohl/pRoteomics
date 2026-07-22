testthat::local_edition(3)

root <- normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE)

testthat::test_that("Stage 13 is a parseable non-circular microglia handoff", {
  stage13 <- file.path(root, "06_modules_WGCNA", "13_wgcna_claim_readiness.R")
  testthat::expect_silent(parse(file = stage13))
  text <- paste(readLines(stage13, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, "WGCNA_entity_claim_readiness.csv", fixed = TRUE)
  testthat::expect_match(text, "conventional_preservation_claim_gate_eligible", fixed = TRUE)
  testthat::expect_match(text, '.data$contrast == "SUS - RES"', fixed = TRUE)
  testthat::expect_match(text, '.data$effect_scope == "spatial_adjusted_global"', fixed = TRUE)
  testthat::expect_match(text, '.data$spatial_unit == "global_spatial_adjusted"', fixed = TRUE)
  lines <- readLines(stage13, warn = FALSE)
  effect_start <- grep("^effect_from <- function", lines)
  effect_end <- grep("^effects <-", lines)
  effect_text <- paste(lines[effect_start:(effect_end - 1L)], collapse = "\n")
  testthat::expect_false(grepl("dplyr::distinct", effect_text, fixed = TRUE))
  testthat::expect_match(effect_text, "any(id_counts != 1L)", fixed = TRUE)
  for (script in c("05_module_supermodule_group_effects.r", "06_annotate_module_microenvironment.r", "07_wgcna_interpretable_summary.r", "12_microglia_wgcna_nature_readiness_audit.R")) {
    source_text <- paste(readLines(file.path(root, "06_modules_WGCNA", script), warn = FALSE), collapse = "\n")
    testthat::expect_false(grepl("13_wgcna_claim_readiness|claim_readiness", source_text))
  }
})

testthat::test_that("completed Stage 13 output has exact current identities without join multiplication", {
  path <- file.path(root, "results", "tables", "06_modules_WGCNA", "claim_readiness", "microglia", "WGCNA_entity_claim_readiness.csv")
  testthat::skip_if_not(file.exists(path), "Stage 13 has not been run")
  x <- read.csv(path, check.names = FALSE)
  testthat::expect_equal(nrow(x), 22L)
  testthat::expect_equal(anyDuplicated(x[c("dataset", "level", "entity_id")]), 0L)
  testthat::expect_setequal(x$entity_id[x$level == "module"], sprintf("WGCNA_m%02d", 1:13))
  testthat::expect_setequal(x$entity_id[x$level == "supermodule"], sprintf("SM%02d", 1:9))
  testthat::expect_false(any(x$conventional_preservation_claim_gate_eligible %in% TRUE))
  source(file.path(root, "R", "paths.R"), local = TRUE)
  source(file.path(root, "R", "schema_validation.R"), local = TRUE)
  testthat::expect_silent(validate_table_schema(x, "wgcna_entity_claim_readiness", strict = TRUE))
  testthat::expect_setequal(unique(x$claim_entity_role), c("canonical_module", "higher_order_block", "compatibility_alias"))
  aliases <- x[x$claim_entity_role == "compatibility_alias", ]
  testthat::expect_equal(nrow(aliases), 6L)
  testthat::expect_true(all(!aliases$separate_manuscript_claim_allowed))
  testthat::expect_true(all(aliases$compatibility_target_level == "module"))
  testthat::expect_identical(aliases$canonical_claim_entity_id, aliases$compatibility_target_entity_id)
  testthat::expect_true(all(aliases$compatibility_target_entity_id %in% x$entity_id[x$level == "module"]))
  testthat::expect_true(all(aliases$manuscript_placement == "compatibility_only"))
  testthat::expect_true(all(grepl("not a separate biological claim", aliases$allowed_wording, fixed = TRUE)))

  claimable_blocks <- x$entity_id[x$claim_entity_role == "higher_order_block" & x$separate_manuscript_claim_allowed]
  testthat::expect_identical(claimable_blocks, "SM09")
  testthat::expect_false(x$separate_manuscript_claim_allowed[x$entity_id == "SM01"])
  testthat::expect_false(x$separate_manuscript_claim_allowed[x$entity_id == "SM03"])
  testthat::expect_equal(anyDuplicated(x$canonical_claim_entity_id[x$separate_manuscript_claim_allowed]), 0L)
})

testthat::test_that("Stage 13 carries the exact selected Stage 05 endpoint and provenance", {
  output_path <- file.path(root, "results", "tables", "06_modules_WGCNA", "claim_readiness", "microglia", "WGCNA_entity_claim_readiness.csv")
  testthat::skip_if_not(file.exists(output_path), "Stage 13 has not been run")
  out <- read.csv(output_path, check.names = FALSE)
  testthat::expect_true(all(out$selected_contrast == "SUS - RES"))
  testthat::expect_true(all(out$selected_effect_scope == "spatial_adjusted_global"))
  testthat::expect_true(all(out$selected_spatial_unit == "global_spatial_adjusted"))
  testthat::expect_equal(out$group_effect_FDR, out$group_effect_FDR_global, tolerance = 0)

  required_provenance <- c(
    "group_effect_estimate", "group_effect_SE", "group_effect_CI_low", "group_effect_CI_high",
    "group_effect_p_value", "group_effect_FDR_within_dataset_level", "group_effect_FDR_global",
    "group_effect_adjustment_scope", "group_effect_evidence_status", "group_effect_model_formula",
    "group_effect_model_type", "group_effect_n_animals", "group_effect_n_samples",
    "group_effect_biological_replicate_unit", "group_effect_source_file", "group_effect_source_key"
  )
  testthat::expect_true(all(required_provenance %in% names(out)))

  compare_level <- function(level, file_name, id_col) {
    source <- read.csv(file.path(root, "results", "tables", "06_modules_WGCNA", "group_effects", "microglia", file_name), check.names = FALSE)
    selected <- source[
      source$dataset == "microglia" & source$contrast == "SUS - RES" &
        source$effect_scope == "spatial_adjusted_global" & source$spatial_unit == "global_spatial_adjusted",
      , drop = FALSE
    ]
    ids <- as.character(selected[[id_col]])
    expected_ids <- out$entity_id[out$level == level]
    testthat::expect_equal(nrow(selected), length(expected_ids))
    testthat::expect_setequal(ids, expected_ids)
    testthat::expect_true(all(table(ids) == 1L))
    observed <- out[out$level == level, ]
    source_match <- selected[match(observed$entity_id, ids), , drop = FALSE]
    testthat::expect_equal(observed$group_effect_FDR_global, source_match$FDR_global, tolerance = 0)
    testthat::expect_equal(observed$group_effect_FDR_within_dataset_level, source_match$FDR_within_dataset_level, tolerance = 0)
    testthat::expect_equal(observed$group_effect_estimate, source_match$estimate, tolerance = 1e-14)
    testthat::expect_equal(observed$group_effect_SE, source_match$SE, tolerance = 1e-14)
    testthat::expect_equal(observed$group_effect_p_value, source_match$p_value, tolerance = 1e-14)
  }
  compare_level("module", "module_group_effects.csv", "module_id")
  compare_level("supermodule", "supermodule_group_effects.csv", "supermodule_id")
})

testthat::test_that("future cut-height defaults remain distinct from frozen microglia provenance", {
  stage01 <- file.path(root, "06_modules_WGCNA", "01_WGCNA.r")
  expressions <- parse(file = stage01)
  is_default_function <- vapply(expressions, function(expr) {
    is.call(expr) && identical(as.character(expr[[1]]), "<-") && identical(as.character(expr[[2]]), "supermodule_merge_cut_height")
  }, logical(1))
  testthat::expect_equal(sum(is_default_function), 1L)
  env <- new.env(parent = baseenv())
  env$validate_dataset <- function(dataset, source = NULL) match.arg(dataset, c("neuron_neuropil", "neuron_soma", "microglia"))
  eval(expressions[[which(is_default_function)]], envir = env)
  testthat::expect_equal(env$supermodule_merge_cut_height("neuron_neuropil"), 0.55)
  testthat::expect_equal(env$supermodule_merge_cut_height("neuron_soma"), 0.35)
  testthat::expect_equal(env$supermodule_merge_cut_height("microglia"), 0.45)

  state_path <- file.path(root, "data", "processed", "06_modules_WGCNA", "01_WGCNA", "microglia", "wgcna_final_model_state.rds")
  state <- readRDS(state_path)
  testthat::expect_equal(as.numeric(state$parameters$supermodule_merge_cut_height), 0.40)
  stage12_text <- paste(readLines(file.path(root, "06_modules_WGCNA", "12_microglia_wgcna_nature_readiness_audit.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(stage12_text, "configured_default_cut_height", fixed = TRUE)
  testthat::expect_match(stage12_text, "selected_network_cut_height", fixed = TRUE)
  testthat::expect_match(stage12_text, "historical_explicit_override", fixed = TRUE)
  testthat::expect_match(stage12_text, "membership_generation_cut_height", fixed = TRUE)
  testthat::expect_match(stage12_text, "0.25, 0.35, 0.40, 0.45, 0.50, 0.55, 0.65", fixed = TRUE)
})

testthat::test_that("pipeline declares Stage 13 and corrected-renderer contracts", {
  pipeline_text <- paste(readLines(file.path(root, "pipeline.yml"), warn = FALSE), collapse = "\n")
  stage13_required <- c(
    "group_effects/microglia/module_group_effects.csv", "group_effects/microglia/supermodule_group_effects.csv",
    "module_annotation/microglia/WGCNA_module_biological_annotation.csv", "module_annotation/microglia/WGCNA_supermodule_biological_annotation.csv",
    "interpretable_summary/microglia/WGCNA_final_label_lookup.csv", "module_robustness_consensus.csv",
    "higher_order_block_readiness_summary.csv", "wgcna_module_supermodule_annotation.csv",
    "claim_readiness/microglia/WGCNA_entity_claim_readiness.csv",
    "results/source_data/06_modules_WGCNA/claim_readiness/microglia/WGCNA_entity_claim_readiness_source.csv"
  )
  for (value in stage13_required) testthat::expect_match(pipeline_text, value, fixed = TRUE)
  testthat::expect_match(pipeline_text, "corrected_multi_supermodule_member_loadings.pdf", fixed = TRUE)
  testthat::expect_match(pipeline_text, "wgcna_readiness_summary_source.csv", fixed = TRUE)
  testthat::expect_match(pipeline_text, "supermodule_group_effects_standardized.csv", fixed = TRUE)
})

testthat::test_that("protected current state and memberships retain audited hashes", {
  protected_path <- file.path(root, "results", "reviewer_audit", "microglia_wgcna_nature_readiness", "protected_output_hash_audit.csv")
  testthat::skip_if_not(file.exists(protected_path), "Protected-hash audit is unavailable")
  protected <- read.csv(protected_path, check.names = FALSE)
  targets <- c(
    "data/processed/06_modules_WGCNA/01_WGCNA/microglia/wgcna_final_model_state.rds",
    "results/tables/06_modules_WGCNA/01_WGCNA/microglia/supermodules/wgcna_module_supermodule_annotation.csv",
    "results/tables/06_modules_WGCNA/group_effects/microglia/module_group_effects.csv",
    "results/tables/06_modules_WGCNA/group_effects/microglia/supermodule_group_effects.csv"
  )
  rows <- protected[match(targets, protected$protected_path), , drop = FALSE]
  testthat::expect_false(anyNA(rows$protected_path))
  current <- unname(tools::md5sum(file.path(root, targets)))
  testthat::expect_identical(current, rows$md5_after)
  testthat::expect_true(all(rows$unchanged))
})
