testthat::local_edition(3)

testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tibble")
source(testthat::test_path("..", "..", "R", "module_contracts.R"))

supermodule_test_map <- function() {
  data.frame(
    dataset = "neuron_soma",
    module_eigengene = c("MEred", "MEblue", "MEgreen"),
    ModuleColor = c("red", "blue", "green"),
    SupermoduleID = c("SM01", "SM01", "SM02"),
    BiologicalLabel = c("energy", "energy", "energy"),
    Supermodule_DisplayLabel = c(
      "SM01 · dominant: energy",
      "SM01 · dominant: energy",
      "SM02 · dominant: energy"
    ),
    stringsAsFactors = FALSE
  )
}

testthat::test_that("one-of-two significant module GO support is not recurring", {
  map <- supermodule_test_map()[1:2, ]
  go <- data.frame(
    ModuleProteinSetType = c("all", "all", "hub"),
    ModuleColor = c("red", "blue", "blue"),
    Ontology = "BP",
    Description = "mitochondrial organization",
    p.adjust = c(0.01, 0.20, 0.001),
    stringsAsFactors = FALSE
  )
  out <- wgcna_supermodule_go_support(go, map)
  testthat::expect_equal(nrow(out), 1L)
  testthat::expect_equal(out$ModuleProteinSetType, "all")
  testthat::expect_equal(out$n_modules_supporting_term, 1L)
  testthat::expect_equal(out$fraction_member_modules_supporting_term, 0.5)
  testthat::expect_equal(out$go_support_confidence, "not_recurring")
  testthat::expect_false(out$recurring_significant_term)
})

testthat::test_that("high and medium recurring GO support obey member-module rules", {
  map <- data.frame(
    dataset = "neuron_soma",
    module_eigengene = paste0("ME", 1:4),
    ModuleColor = paste0("c", 1:4),
    SupermoduleID = "SM01",
    stringsAsFactors = FALSE
  )
  go <- dplyr::bind_rows(
    data.frame(ModuleProteinSetType = "all", ModuleColor = paste0("c", 1:4), Ontology = "BP", Description = "translation", p.adjust = c(.01, .02, .03, .04)),
    data.frame(ModuleProteinSetType = "all", ModuleColor = paste0("c", 1:4), Ontology = "BP", Description = "RNA processing", p.adjust = c(.01, .02, .20, .30))
  )
  out <- wgcna_supermodule_go_support(go, map)
  high <- out[out$Description == "translation", ]
  medium <- out[out$Description == "RNA processing", ]
  testthat::expect_equal(high$go_support_confidence, "high")
  testthat::expect_equal(high$worst_member_module_fdr, .04)
  testthat::expect_equal(medium$go_support_confidence, "medium")
  testthat::expect_equal(medium$fraction_member_modules_supporting_term, .5)
})

testthat::test_that("duplicate biological labels remain separate while display collisions fail", {
  map <- supermodule_test_map()
  testthat::expect_silent(validate_supermodule_member_map(
    map,
    expected_modules = map$module_eigengene,
    display_col = "Supermodule_DisplayLabel"
  ))
  testthat::expect_equal(length(unique(map$BiologicalLabel)), 1L)
  testthat::expect_equal(length(unique(map$SupermoduleID)), 2L)

  collision <- map
  collision$Supermodule_DisplayLabel[collision$SupermoduleID == "SM02"] <- "SM01 · dominant: energy"
  testthat::expect_error(
    validate_supermodule_member_map(collision, display_col = "Supermodule_DisplayLabel"),
    "display-label collision"
  )
})

testthat::test_that("singleton GO evidence and naming confidence remain low", {
  map <- supermodule_test_map()[1, ]
  go <- data.frame(
    ModuleProteinSetType = "all", ModuleColor = "red", Ontology = "BP",
    Description = "mitochondrial organization", p.adjust = .001
  )
  out <- wgcna_supermodule_go_support(go, map)
  testthat::expect_equal(out$go_support_confidence, "singleton")
  old_root <- Sys.getenv("PROTEOMICS_PROJECT_ROOT", unset = NA_character_)
  on.exit(if (is.na(old_root)) Sys.unsetenv("PROTEOMICS_PROJECT_ROOT") else Sys.setenv(PROTEOMICS_PROJECT_ROOT = old_root), add = TRUE)
  Sys.setenv(PROTEOMICS_PROJECT_ROOT = normalizePath(testthat::test_path("..", ".."), winslash = "/", mustWork = TRUE))
  source(testthat::test_path("..", "..", "R", "wgcna_downstream_utils.R"), chdir = TRUE)
  testthat::expect_equal(classify_supermodule_label_confidence(1L, "high"), "low")
})

testthat::test_that("hub evidence is audit-only for supermodule labels and confidence", {
  source(testthat::test_path("..", "..", "R", "wgcna_labeling_utils.R"))
  super_rows <- data.frame(
    dataset = "neuron_soma",
    supermodule_id = "SM01",
    n_member_modules = 2L,
    Supermodule_DisplayLabel = "SM01 \u00b7 dominant: mitochondrial / energy metabolism",
    Supermodule_LabelConfidence = "medium",
    cleaned_biological_label = "mitochondrial / energy metabolism",
    MemberEvidenceHubs = "HubA;HubB",
    stringsAsFactors = FALSE
  )
  candidates <- wgcna_make_label_candidates(
    module_rows = data.frame(dataset = character(), module_id = character()),
    supermodule_rows = super_rows,
    dataset = "neuron_soma"
  ) |>
    wgcna_score_label_candidates() |>
    wgcna_select_final_labels()
  testthat::expect_true(all(candidates$hub_support == 0))
  testthat::expect_true(all(candidates$hub_score_contribution == 0))
  lookup <- wgcna_build_final_label_lookup(
    candidates,
    module_rows = data.frame(dataset = character(), module_id = character()),
    supermodule_rows = super_rows,
    dataset = "neuron_soma"
  )
  testthat::expect_equal(lookup$final_plot_label, super_rows$Supermodule_DisplayLabel)
  testthat::expect_equal(lookup$label_confidence, "medium")
})

testthat::test_that("coherence validation retains signed correlations", {
  x <- seq(-2.5, 2.5, length.out = 8)
  mes <- data.frame(MEred = x, MEblue = -x, MEgreen = x^2)
  map <- supermodule_test_map()
  out <- wgcna_supermodule_coherence(mes, map)
  sm1 <- out[out$SupermoduleID == "SM01", ]
  testthat::expect_equal(sm1$signed_min_pairwise_eigengene_correlation, -1, tolerance = 1e-10)
  testthat::expect_equal(sm1$signed_mean_pairwise_eigengene_correlation, -1, tolerance = 1e-10)
  testthat::expect_equal(sm1$correlation_metric, "signed_pearson")
  testthat::expect_match(sm1$construction_basis, "1_minus_signed")
  testthat::expect_match(sm1$pc1_interpretation_note, "not proof")
})

testthat::test_that("authoritative cluster IDs exactly match downstream summaries", {
  map <- supermodule_test_map()
  exact <- unique(map[c("dataset", "SupermoduleID")])
  testthat::expect_silent(validate_supermodule_summary_ids(exact, map))
  testthat::expect_error(validate_supermodule_summary_ids(exact[1, ], map), "missing")
  unexpected <- rbind(exact, data.frame(dataset = "neuron_soma", SupermoduleID = "SM03"))
  testthat::expect_error(validate_supermodule_summary_ids(unexpected, map), "unexpected")
})

testthat::test_that("member-map validation rejects duplicates, missing members, and stale IDs", {
  map <- supermodule_test_map()
  testthat::expect_error(validate_supermodule_member_map(rbind(map, map[1, ])), "duplicate")
  testthat::expect_error(
    validate_supermodule_member_map(map[-1, ], expected_modules = map$module_eigengene),
    "missing"
  )
  stale <- map
  stale$SupermoduleID[stale$SupermoduleID == "SM02"] <- "neuropil"
  testthat::expect_error(validate_supermodule_member_map(stale), "unexpected SupermoduleID")
})

testthat::test_that("stale or missing composition columns fail clearly", {
  stale <- data.frame(dataset = "neuron_soma", SupermoduleID = "SM01")
  testthat::expect_error(require_supermodule_composition_columns(stale), "missing required column")
  empty <- data.frame(
    dataset = "neuron_soma", SupermoduleID = "SM01",
    Supermodule_CompositionLabel = "", Supermodule_CompositionDisplayLabel = "",
    Supermodule_CompositionConfidence = "low", Supermodule_CompositionRationale = "missing"
  )
  testthat::expect_error(require_supermodule_composition_columns(empty), "stale or empty")
})

testthat::test_that("supermodule sensitivity grid is complete and stable", {
  testthat::expect_identical(wgcna_supermodule_sensitivity_grid(), c(.25, .35, .45, .50, .55, .65))
})

testthat::test_that("construction source keys summaries by stable IDs and treats hubs as audit only", {
  script <- paste(readLines(testthat::test_path("..", "..", "06_modules_WGCNA", "01_WGCNA.r"), warn = FALSE), collapse = "\n")
  testthat::expect_match(script, "average_linkage_on_1_minus_signed_module_eigengene_correlation", fixed = TRUE)
  testthat::expect_match(script, "partition_integrity_diagnostic_not_supermodule_support", fixed = TRUE)
  testthat::expect_match(script, "expected_module_pairs", fixed = TRUE)
  testthat::expect_match(script, 'display = paste0(supermodule_id, " \\u00b7 singleton")', fixed = TRUE)
  testthat::expect_false(grepl("group_by\\(\\.data\\$Supermodule\\)", script))
})

testthat::test_that("active downstream consumers do not use labels as supermodule keys", {
  paths <- testthat::test_path("..", "..", "06_modules_WGCNA", c(
    "01_WGCNA.r", "03_score_module_activity.R",
    "05_module_supermodule_group_effects.r",
    "06_annotate_module_microenvironment.r",
    "07_wgcna_interpretable_summary.r"
  ))
  script_lines <- unlist(lapply(paths, readLines, warn = FALSE), use.names = FALSE)
  scripts <- paste(script_lines, collapse = "\n")
  testthat::expect_false(any(grepl("group_by\\(.*(Supermodule_DisplayLabel|supermodule_plot_label|Macroprogram_Display)", script_lines)))
  testthat::expect_false(any(grepl("distinct\\(.*(Supermodule_DisplayLabel|supermodule_plot_label|Macroprogram_Display)", script_lines)))
  testthat::expect_match(scripts, "validate_supermodule_member_map", fixed = TRUE)
  testthat::expect_match(scripts, "Module-to-supermodule inputs conflict for dataset + module key", fixed = TRUE)
  testthat::expect_match(scripts, 'by = c("dataset", "supermodule_id")', fixed = TRUE)
})
