source(testthat::test_path("..", "..", "R", "paths.R"))

testthat::local_edition(3)

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("readr")
testthat::skip_if_not_installed("digest")

source(repo_path("R", "wgcna_reviewed_label_registry.R"))

# Current canonical Stage 01 outputs (analysis/wgcna/build_wgcna_modules.R).
member_map_path <- current_member_map_path("microglia")
summary_path <- current_wgcna("build_wgcna_modules", "microglia", "supermodules",
                              "wgcna_supermodule_summary.csv")

testthat::test_that("the microglia registry in force is exact, current, and coherent", {
  reg <- wgcna_registry_in_force("microglia")
  registry <- readr::read_csv(reg$path, show_col_types = FALSE)
  member_map <- readr::read_csv(member_map_path, show_col_types = FALSE)
  shape <- current_network_shape("microglia")
  testthat::expect_silent(wgcna_validate_reviewed_registry(
    registry, "microglia", member_map, allow_provisional = reg$allow_provisional
  ))
  testthat::expect_equal(nrow(registry), shape$n_identities)
  testthat::expect_setequal(registry$entity_id[registry$level == "module"], shape$modules)
  testthat::expect_setequal(registry$entity_id[registry$level == "supermodule"], shape$supermodules)
  testthat::expect_true(all(registry$label_confidence %in% c("high", "moderate", "low")))
  testthat::expect_true(all(nzchar(registry$proposal_prepared_by)))
  testthat::expect_false(anyDuplicated(registry[c("dataset", "level", "entity_id")]) > 0L)
  if (identical(reg$mode, "reviewed")) {
    testthat::expect_true(all(registry$adjudication_status == "reviewed"))
    # adjudication_status = reviewed must be backed by a completed human review
    complete <- reviewed_rows_are_complete(registry)
    testthat::expect_true(all(complete),
                          info = paste("incomplete review:", paste(names(complete)[!complete], collapse = ", ")))
    testthat::expect_false(any(registry$manual_review_required))
  } else {
    testthat::expect_true(all(registry$adjudication_status == "provisional_unreviewed"))
    testthat::expect_true(all(registry$manual_review_required))
  }
})

testthat::test_that("singleton rows inherit labels without confidence penalty", {
  reg <- wgcna_registry_in_force("microglia")
  registry <- readr::read_csv(reg$path, show_col_types = FALSE)
  shape <- current_network_shape("microglia")
  member_map <- wgcna_normalize_current_member_map(readr::read_csv(member_map_path, show_col_types = FALSE), "microglia")
  singleton <- member_map |>
    dplyr::group_by(.data$SupermoduleID) |>
    dplyr::summarise(n = dplyr::n(), member = dplyr::first(.data$ModuleID), .groups = "drop") |>
    dplyr::filter(.data$n == 1L)
  inherited <- singleton |>
    dplyr::left_join(
      registry |> dplyr::filter(.data$level == "module") |>
        dplyr::select(member = "entity_id", module_bio = "reviewed_biological_label", module_short = "reviewed_short_label", module_conf = "label_confidence"),
      by = "member", relationship = "many-to-one"
    ) |>
    dplyr::left_join(
      registry |> dplyr::filter(.data$level == "supermodule") |>
        dplyr::select(SupermoduleID = "entity_id", super_bio = "reviewed_biological_label", super_short = "reviewed_short_label", super_conf = "label_confidence", structural_status, aggregation_evidence_class),
      by = "SupermoduleID", relationship = "many-to-one"
    )
  testthat::expect_equal(nrow(inherited), length(shape$singletons))
  testthat::expect_setequal(inherited$SupermoduleID, names(shape$singletons))
  testthat::expect_identical(inherited$module_bio, inherited$super_bio)
  testthat::expect_identical(inherited$module_short, inherited$super_short)
  testthat::expect_identical(inherited$module_conf, inherited$super_conf)
  testthat::expect_true(all(inherited$structural_status == "singleton"))
  testthat::expect_true(all(inherited$aggregation_evidence_class == "not_applicable_singleton"))
})

testthat::test_that("canonical lookup has exact fields and current fingerprints", {
  reg <- wgcna_registry_in_force("microglia")
  registry <- readr::read_csv(reg$path, show_col_types = FALSE)
  member_map <- readr::read_csv(member_map_path, show_col_types = FALSE)
  summary <- readr::read_csv(summary_path, show_col_types = FALSE)
  shape <- current_network_shape("microglia")
  module_rows <- member_map |>
    dplyr::transmute(module_id = .data$ModuleID, raw_GO_BP_terms = .data$best_GO_BP, raw_GO_CC_terms = .data$best_GO_CC, raw_GO_MF_terms = .data$best_GO_MF)
  super_rows <- summary |>
    dplyr::mutate(n_member_modules = suppressWarnings(as.integer(.data$n_modules)))
  candidates <- tibble::tibble(
    level = character(), entity_id = character(), candidate_source = character(),
    candidate_label = character(), final_label_score = numeric()
  )
  lookup <- wgcna_build_reviewed_canonical_lookup(
    registry, member_map, candidates, module_rows, super_rows, "microglia",
    allow_provisional = reg$allow_provisional
  )
  testthat::expect_silent(wgcna_validate_canonical_lookup(lookup, "microglia", member_map))
  # The Stage 07 lookup on disk carries exactly the labels of the registry in
  # force (a registry edit that was not propagated fails here).
  written <- read_current_csv(current_wgcna("summarize_module_interpretation", "microglia",
                                            "WGCNA_final_label_lookup.csv"))
  key <- c("level", "entity_id")
  cols <- c("canonical_biological_label", "canonical_short_label",
            "biological_label_confidence", "member_modules_fingerprint")
  joined <- merge(as.data.frame(lookup)[c(key, cols)], written[c(key, cols)],
                  by = key, suffixes = c("_registry", "_stage07"))
  testthat::expect_equal(nrow(joined), nrow(lookup))
  for (col in cols) {
    testthat::expect_identical(as.character(joined[[paste0(col, "_stage07")]]),
                               as.character(joined[[paste0(col, "_registry")]]), info = col)
  }
  testthat::expect_identical(names(lookup), wgcna_canonical_lookup_columns())
  testthat::expect_equal(nrow(lookup), shape$n_identities)
  testthat::expect_identical(lookup$canonical_plot_label, lookup$final_plot_label)
  testthat::expect_true(all(nchar(lookup$member_modules_fingerprint) == 64L))
  testthat::expect_true(all(nchar(lookup$annotation_evidence_fingerprint) == 64L))
  testthat::expect_false(any(grepl("\\b(SUS|RES|CON|higher|lower)\\b", lookup$canonical_biological_label, ignore.case = TRUE)))
  singleton <- lookup$level == "supermodule" & lookup$structural_status == "singleton"
  testthat::expect_equal(sum(singleton), length(shape$singletons))
  testthat::expect_true(all(grepl("^SM[0-9]{2}\\n", lookup$canonical_plot_label[singleton])))
  testthat::expect_true(all(lookup$structural_coherence_class[singleton] == "not_applicable_singleton"))
})

testthat::test_that("Stage 06 and publication code enforce stable non-circular contracts", {
  stage06 <- paste(readLines(repo_path("analysis/wgcna", "annotate_module_microenvironment.R"), warn = FALSE), collapse = "\n")
  stage08 <- paste(readLines(repo_path("analysis/wgcna", "render_module_figures.R"), warn = FALSE), collapse = "\n")
  score_publication <- paste(readLines(repo_path("analysis/wgcna", "summarize_module_scores.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(stage06, "duplicate SM09 generation is forbidden", fixed = TRUE)
  testthat::expect_match(stage06, "member_modules_fingerprint", fixed = TRUE)
  testthat::expect_match(stage06, "automatic_biological_label", fixed = TRUE)
  testthat::expect_false(grepl("n_member_modules <= 1L, \"low\"", stage06, fixed = TRUE))
  testthat::expect_match(stage08, "multi_ids <- sort(", fixed = TRUE)
  testthat::expect_match(stage08, "Singleton loading panels are forbidden", fixed = TRUE)
  testthat::expect_match(stage08, "singular_models_have_no_claim_symbols", fixed = TRUE)
  testthat::expect_false(grepl("selected_significant|selected significant supermodules", stage08, ignore.case = TRUE))
  testthat::expect_false(grepl("coalesce\\(", stage08))
  testthat::expect_false(grepl("legacy_multifile|Supermodule_CleanPlotLabel|Supermodule_CompositionDisplayLabel", score_publication))
})

# Current canonical figure sources (analysis/wgcna/render_module_figures.R).
testthat::test_that("current publication sources preserve exact all-supermodule cardinality", {
  source_dir <- current_wgcna("render_module_figures", "microglia", "source_data")
  required <- file.path(source_dir, c(
    "all_supermodule_architecture_source.csv",
    "all_supermodule_global_eigengenes_source.csv",
    "all_supermodule_sus_res_forest_source.csv",
    "all_supermodule_sus_res_spatial_heatmap_source.csv",
    "multi_supermodule_member_loadings_source.csv",
    "supplementary_all_contrast_effect_matrix_source.csv"
  ))
  testthat::skip_if_not(all(file.exists(required)), "Current publication sources are not generated.")
  shape <- current_network_shape("microglia")
  n_super <- length(shape$supermodules)

  sources <- lapply(required, readr::read_csv, show_col_types = FALSE)
  names(sources) <- c("architecture", "global", "forest", "spatial", "loadings", "matrix")
  testthat::expect_true(all(vapply(sources, function(x) identical(unique(x$dataset), "microglia"), logical(1))))
  testthat::expect_true(all(vapply(sources, function(x) !any(c("dataset.x", "dataset.y") %in% names(x)), logical(1))))
  testthat::expect_true(all(vapply(sources, function(x) !any(grepl("^(roi_context|canonical_short_label)\\.[xy]$", names(x))), logical(1))))
  testthat::expect_setequal(unique(sources$architecture$supermodule_id), shape$supermodules)
  testthat::expect_setequal(unique(sources$global$supermodule_id), shape$supermodules)
  testthat::expect_equal(nrow(sources$forest), n_super)
  # the design is fixed (four ROIs, three contrasts); only the supermodule count follows the network
  testthat::expect_setequal(unique(sources$spatial$roi), c("CA1", "CA2", "CA3", "DG"))
  testthat::expect_setequal(unique(sources$matrix$contrast), c("RES - CON", "SUS - CON", "SUS - RES"))
  testthat::expect_equal(nrow(sources$spatial), n_super * 4L)
  testthat::expect_equal(nrow(sources$matrix), n_super * 3L)
  testthat::expect_false(anyDuplicated(sources$forest$supermodule_id) > 0L)
  testthat::expect_false(anyDuplicated(sources$spatial[c("supermodule_id", "roi")]) > 0L)
  testthat::expect_false(anyDuplicated(sources$matrix[c("supermodule_id", "contrast")]) > 0L)
  testthat::expect_setequal(unique(sources$loadings$supermodule_id), shape$blocks)
  testthat::expect_true(all(sources$loadings$n_member_modules > 1L))
  diagnostic <- !sources$spatial$model_stable
  testthat::expect_true(all(sources$spatial$support_symbol[diagnostic] == "cross_diagnostic_only"))
})

# Current canonical corrected sources (analysis/wgcna/render_microglia_module_figures.R).
testthat::test_that("corrected publication effect sources use unsuffixed canonical labels", {
  source_dir <- current_wgcna("render_microglia_module_figures", "microglia", "source_data")
  required <- file.path(source_dir, c(
    "corrected_all_supermodule_sus_res_forest_source.csv",
    "corrected_all_supermodule_sus_res_spatial_heatmap_source.csv"
  ))
  testthat::skip_if_not(all(file.exists(required)), "Current corrected publication sources are not generated.")

  sources <- lapply(required, readr::read_csv, show_col_types = FALSE)
  testthat::expect_true(all(vapply(sources, function(x) "canonical_short_label" %in% names(x), logical(1))))
  testthat::expect_true("roi_context" %in% names(sources[[2]]))
  testthat::expect_true(all(vapply(sources, function(x) !any(grepl("^(roi_context|canonical_short_label)\\.[xy]$", names(x))), logical(1))))
})
