testthat::local_edition(3)

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("readr")
testthat::skip_if_not_installed("digest")

source(testthat::test_path("..", "..", "R", "wgcna_reviewed_label_registry.R"))

registry_path <- testthat::test_path("..", "..", "config", "wgcna_labels", "microglia.csv")
member_map_path <- testthat::test_path(
  "..", "..", "results", "tables", "06_modules_WGCNA", "01_WGCNA", "microglia",
  "supermodules", "wgcna_module_supermodule_annotation.csv"
)
summary_path <- testthat::test_path(
  "..", "..", "results", "tables", "06_modules_WGCNA", "01_WGCNA", "microglia",
  "supermodules", "wgcna_supermodule_summary.csv"
)

testthat::test_that("reviewed microglia registry is exact, current, and fully adjudicated", {
  registry <- readr::read_csv(registry_path, show_col_types = FALSE)
  member_map <- readr::read_csv(member_map_path, show_col_types = FALSE)
  testthat::expect_silent(wgcna_validate_reviewed_registry(registry, "microglia", member_map))
  testthat::expect_equal(nrow(registry), 22L)
  testthat::expect_equal(sum(registry$level == "module"), 13L)
  testthat::expect_equal(sum(registry$level == "supermodule"), 9L)
  testthat::expect_true(all(registry$adjudication_status == "reviewed"))
  testthat::expect_false(any(registry$manual_review_required))
  testthat::expect_setequal(unique(registry$label_confidence), c("high", "moderate", "low"))
  testthat::expect_true(all(registry$reviewer == "Tobias Pohl"))
  testthat::expect_true(all(nzchar(registry$proposal_prepared_by)))
  testthat::expect_false(anyDuplicated(registry[c("dataset", "level", "entity_id")]) > 0L)
})

testthat::test_that("six singleton rows inherit labels without confidence penalty", {
  registry <- readr::read_csv(registry_path, show_col_types = FALSE)
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
  testthat::expect_equal(nrow(inherited), 6L)
  testthat::expect_identical(inherited$module_bio, inherited$super_bio)
  testthat::expect_identical(inherited$module_short, inherited$super_short)
  testthat::expect_identical(inherited$module_conf, inherited$super_conf)
  testthat::expect_true(all(inherited$structural_status == "singleton"))
  testthat::expect_true(all(inherited$aggregation_evidence_class == "not_applicable_singleton"))
})

testthat::test_that("reviewed canonical lookup has exact fields and current fingerprints", {
  registry <- readr::read_csv(registry_path, show_col_types = FALSE)
  member_map <- readr::read_csv(member_map_path, show_col_types = FALSE)
  summary <- readr::read_csv(summary_path, show_col_types = FALSE)
  module_rows <- member_map |>
    dplyr::transmute(module_id = .data$ModuleID, raw_GO_BP_terms = .data$best_GO_BP, raw_GO_CC_terms = .data$best_GO_CC, raw_GO_MF_terms = .data$best_GO_MF)
  super_rows <- summary |>
    dplyr::mutate(n_member_modules = suppressWarnings(as.integer(.data$n_modules)))
  candidates <- tibble::tibble(
    level = character(), entity_id = character(), candidate_source = character(),
    candidate_label = character(), final_label_score = numeric()
  )
  lookup <- wgcna_build_reviewed_canonical_lookup(
    registry, member_map, candidates, module_rows, super_rows, "microglia"
  )
  testthat::expect_silent(wgcna_validate_canonical_lookup(lookup, "microglia", member_map))
  testthat::expect_identical(names(lookup), wgcna_canonical_lookup_columns())
  testthat::expect_equal(nrow(lookup), 22L)
  testthat::expect_identical(lookup$canonical_plot_label, lookup$final_plot_label)
  testthat::expect_true(all(nchar(lookup$member_modules_fingerprint) == 64L))
  testthat::expect_true(all(nchar(lookup$annotation_evidence_fingerprint) == 64L))
  testthat::expect_false(any(grepl("\\b(SUS|RES|CON|higher|lower)\\b", lookup$canonical_biological_label, ignore.case = TRUE)))
  singleton <- lookup$level == "supermodule" & lookup$structural_status == "singleton"
  testthat::expect_true(all(grepl("^SM[0-9]{2}\\n", lookup$canonical_plot_label[singleton])))
  testthat::expect_true(all(lookup$structural_coherence_class[singleton] == "not_applicable_singleton"))
})

testthat::test_that("Stage 06 and publication code enforce stable non-circular contracts", {
  stage06 <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "06_annotate_module_microenvironment.r"), warn = FALSE), collapse = "\n")
  stage08 <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "08_wgcna_publication_figures.R"), warn = FALSE), collapse = "\n")
  score_publication <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "08_wgcna_score_publication_summary.R"), warn = FALSE), collapse = "\n")
  testthat::expect_match(stage06, "duplicate SM09 generation is forbidden", fixed = TRUE)
  testthat::expect_match(stage06, "member_modules_fingerprint", fixed = TRUE)
  testthat::expect_match(stage06, "automatic_biological_label", fixed = TRUE)
  testthat::expect_false(grepl("n_member_modules <= 1L, \"low\"", stage06, fixed = TRUE))
  testthat::expect_match(stage08, "multi_ids <- c(\"SM01\", \"SM03\", \"SM09\")", fixed = TRUE)
  testthat::expect_match(stage08, "Singleton loading panels are forbidden", fixed = TRUE)
  testthat::expect_match(stage08, "singular_models_have_no_claim_symbols", fixed = TRUE)
  testthat::expect_false(grepl("selected_significant|selected significant supermodules", stage08, ignore.case = TRUE))
  testthat::expect_false(grepl("coalesce\\(", stage08))
  testthat::expect_false(grepl("legacy_multifile|Supermodule_CleanPlotLabel|Supermodule_CompositionDisplayLabel", score_publication))
})

testthat::test_that("current publication sources preserve exact all-supermodule cardinality", {
  source_dir <- testthat::test_path(
    "..", "..", "results", "source_data", "06_modules_WGCNA",
    "wgcna_publication_figures", "microglia"
  )
  required <- file.path(source_dir, c(
    "all_supermodule_architecture_source.csv",
    "all_supermodule_global_eigengenes_source.csv",
    "all_supermodule_sus_res_forest_source.csv",
    "all_supermodule_sus_res_spatial_heatmap_source.csv",
    "multi_supermodule_member_loadings_source.csv",
    "supplementary_all_contrast_effect_matrix_source.csv"
  ))
  testthat::skip_if_not(all(file.exists(required)), "Generated publication sources are not present.")

  sources <- lapply(required, readr::read_csv, show_col_types = FALSE)
  names(sources) <- c("architecture", "global", "forest", "spatial", "loadings", "matrix")
  testthat::expect_true(all(vapply(sources, function(x) identical(unique(x$dataset), "microglia"), logical(1))))
  testthat::expect_true(all(vapply(sources, function(x) !any(c("dataset.x", "dataset.y") %in% names(x)), logical(1))))
  testthat::expect_equal(dplyr::n_distinct(sources$architecture$supermodule_id), 9L)
  testthat::expect_equal(dplyr::n_distinct(sources$global$supermodule_id), 9L)
  testthat::expect_equal(nrow(sources$forest), 9L)
  testthat::expect_equal(nrow(sources$spatial), 36L)
  testthat::expect_equal(nrow(sources$matrix), 27L)
  testthat::expect_false(anyDuplicated(sources$forest$supermodule_id) > 0L)
  testthat::expect_false(anyDuplicated(sources$spatial[c("supermodule_id", "roi")]) > 0L)
  testthat::expect_false(anyDuplicated(sources$matrix[c("supermodule_id", "contrast")]) > 0L)
  testthat::expect_setequal(unique(sources$loadings$supermodule_id), c("SM01", "SM03", "SM09"))
  testthat::expect_true(all(sources$loadings$n_member_modules > 1L))
  diagnostic <- !sources$spatial$model_stable
  testthat::expect_true(all(sources$spatial$support_symbol[diagnostic] == "cross_diagnostic_only"))
})
