testthat::test_that("WGCNA final labeling is unique, singleton-aware, and context-separated", {
  source(testthat::test_path("..", "..", "R", "wgcna_labeling_utils.R"))
  testthat::skip_if_not_installed("dplyr")
  testthat::skip_if_not_installed("tibble")

  modules <- tibble::tibble(
    dataset = "microglia", module_id = c("M01", "M02"), supermodule_id = c("SM01", "SM02"),
    cleaned_biological_label = c("shared ROI: synaptic scaffold", "RNA regulation"),
    raw_top_GO_label = c("synaptic membrane", NA_character_), raw_hub_proteins = c("SYN1", NA_character_),
    module_program_primary = c("synaptic", "RNA"), microenvironment_caution_label = "shared ROI",
    GO_label_relevance_flag = c("relevant", "ontology_mismatch"), label_warning = c("", "ontology mismatch")
  )
  supers <- tibble::tibble(
    dataset = "microglia", supermodule_id = c("SM01", "SM02"), n_member_modules = c(1L, 2L),
    DominantMemberTheme = c("synaptic scaffold", "RNA regulation"), DominantMemberThemeFraction = c(1, 0.5),
    MemberThemeFractions = c("synaptic scaffold=1.00", "RNA regulation=0.50; translation=0.50"),
    raw_top_GO_label = c("synaptic membrane", NA_character_), Macroprogram_Display = c("synaptic", "RNA"),
    dominant_microenvironment_class = "shared ROI", GO_label_relevance_flag = c("relevant", "ontology_mismatch"),
    label_warning = c("", "ontology mismatch")
  )

  candidates <- wgcna_make_label_candidates(modules, supers, "microglia") |>
    wgcna_score_label_candidates() |>
    wgcna_select_final_labels()
  lookup <- wgcna_build_final_label_lookup(candidates, modules, supers, "microglia")
  testthat::expect_silent(wgcna_validate_label_lookup(lookup))
  testthat::expect_false(any(duplicated(lookup[c("dataset", "level", "entity_id")])))
  testthat::expect_true(all(c(
    "dataset", "level", "entity_id", "candidate_label", "candidate_source",
    "evidence_strength", "hub_support", "marker_context_support", "genericity_penalty",
    "ontology_mismatch_flag", "conflict_penalty", "final_label_score", "selected_label", "rejection_reason"
  ) %in% names(candidates)))
  singleton <- lookup$level == "supermodule" & lookup$n_member_modules == 1L
  testthat::expect_true(all(grepl("singleton:", lookup$final_plot_label[singleton], fixed = TRUE)))
  testthat::expect_false(any(grepl("mostly", lookup$final_plot_label[singleton], ignore.case = TRUE)))
  testthat::expect_false(any(grepl("shared ROI", lookup$best_data_driven_label, ignore.case = TRUE)))
  mismatch_without_evidence <- lookup$ontology_mismatch_flag & !lookup$independent_evidence_recorded
  testthat::expect_true(all(lookup$manual_review_required[mismatch_without_evidence]))
})

testthat::test_that("joining canonical labels leaves numeric source values unchanged", {
  source(testthat::test_path("..", "..", "R", "wgcna_labeling_utils.R"))
  source_numeric <- tibble::tibble(Module = c("SM01", "SM02"), estimate = c(0.25, -0.5), p_value = c(0.01, 0.2))
  lookup <- tibble::tibble(entity_id = c("SM01", "SM02"), final_plot_label = c("SM01 · singleton: RNA", "SM02 · mixed / unresolved"))
  labelled <- dplyr::left_join(source_numeric, lookup, by = c("Module" = "entity_id"))
  testthat::expect_identical(source_numeric$estimate, labelled$estimate)
  testthat::expect_identical(source_numeric$p_value, labelled$p_value)
})

testthat::test_that("score publication script consumes only the canonical final lookup for active label choice", {
  txt <- readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "08_wgcna_score_publication_summary.R"), warn = FALSE)
  active <- paste(txt[!grepl("^legacy_multifile_", trimws(txt))], collapse = "\n")
  testthat::expect_match(active, "WGCNA_final_label_lookup.csv", fixed = TRUE)
  testthat::expect_match(active, "validate_render_ids_in_lookup", fixed = TRUE)
  testthat::expect_match(active, "canonical_final_label_lookup", fixed = TRUE)
})
