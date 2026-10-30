testthat::local_edition(3)

# Reviewed labels for SOME entities of a dataset (config/wgcna_labels/neuron_soma.csv), applied by
# Stage 07 over the automatic labels: R/statistics/wgcna_reviewed_label_registry.R. Fixtures only.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "wgcna_reviewed_label_registry.R"))

member_map <- data.frame(dataset = "neuron_soma",
                         ModuleID = c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03"),
                         SupermoduleID = c("SM04", "SM03", "SM04"), stringsAsFactors = FALSE)

partial_fixture <- function(...) {
  r <- data.frame(
    dataset = "neuron_soma", level = c("module", "supermodule"), entity_id = c("WGCNA_m01", "SM04"),
    reviewed_biological_label = c("ubiquitin-dependent protein catabolism", "Compound: a (m01) + b (m03)"),
    reviewed_short_label = c("ubiquitin-dependent protein catabolism", "a (m01) + b (m03)"),
    subcellular_context = "context", roi_context = "roi", label_confidence = "moderate",
    manual_review_required = FALSE, aggregation_evidence_class = c("not_applicable_module", "mixed_supported_themes"),
    structural_status = c("module", "multi_module"), rationale = "fixture rationale",
    reviewer = "Reviewer Name", proposal_prepared_by = "fixture", review_date = "2026-10-27",
    adjudication_status = "reviewed", member_ModuleIDs = c("", "WGCNA_m01;WGCNA_m03"),
    stringsAsFactors = FALSE)
  mods <- list(...)
  for (nm in names(mods)) r[[nm]] <- mods[[nm]]
  r
}

lookup_fixture <- data.frame(
  dataset = "neuron_soma", level = c("module", "module", "module", "supermodule", "supermodule"),
  entity_id = c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03", "SM03", "SM04"), parent_entity_id = NA_character_,
  n_member_modules = c(1L, 1L, 1L, 1L, 2L), raw_top_GO_label = "go", best_data_driven_label = "automatic",
  parent_program = "p", microenvironment_context = "c",
  final_plot_label = paste(c("WGCNA_m01", "WGCNA_m02", "WGCNA_m03", "SM03", "SM04"), "· automatic"),
  label_confidence = "high", label_basis = "auto", label_score = 6.1, manual_review_required = FALSE,
  ontology_mismatch_flag = FALSE, independent_evidence_recorded = TRUE, label_rationale = "auto",
  unsafe_interpretation = "u", stringsAsFactors = FALSE)

testthat::test_that("reviewed overrides change exactly the listed entities", {
  out <- wgcna_apply_reviewed_overrides(lookup_fixture, partial_fixture(), member_map, "neuron_soma", "fixture.csv")
  changed <- vapply(seq_len(nrow(out)), function(i) !identical(as.list(out[i, ]), as.list(lookup_fixture[i, ])), logical(1))
  testthat::expect_identical(out$entity_id[changed], c("WGCNA_m01", "SM04"))
  testthat::expect_identical(out$final_plot_label[out$entity_id == "SM04"], "SM04 · a (m01) + b (m03)")
  testthat::expect_true(all(grepl("^reviewed_registry=fixture.csv", out$label_basis[changed])))
  testthat::expect_true(all(is.na(out$label_score[changed])))
})

testthat::test_that("the same members listed in another order are accepted", {
  testthat::expect_true(wgcna_validate_partial_reviewed_registry(
    partial_fixture(member_ModuleIDs = c("", "WGCNA_m03; WGCNA_m01")), "neuron_soma", member_map))
})

testthat::test_that("a supermodule row reviewed for other members is refused", {
  moved <- member_map
  moved$SupermoduleID[moved$ModuleID == "WGCNA_m03"] <- "SM03"
  testthat::expect_error(wgcna_validate_partial_reviewed_registry(partial_fixture(), "neuron_soma", moved),
                         "made for other members")
  testthat::expect_error(wgcna_validate_partial_reviewed_registry(
    partial_fixture(member_ModuleIDs = c("", "")), "neuron_soma", member_map), "made for other members")
})

testthat::test_that("an unfinished review, an unknown entity or a group word is refused", {
  testthat::expect_error(wgcna_validate_partial_reviewed_registry(
    partial_fixture(reviewer = "pending review"), "neuron_soma", member_map), "must name their reviewer")
  testthat::expect_error(wgcna_validate_partial_reviewed_registry(
    partial_fixture(manual_review_required = TRUE), "neuron_soma", member_map), "must not still require")
  testthat::expect_error(wgcna_validate_partial_reviewed_registry(
    partial_fixture(entity_id = c("WGCNA_m09", "SM04")), "neuron_soma", member_map), "absent from the current member map")
  testthat::expect_error(wgcna_validate_partial_reviewed_registry(
    partial_fixture(reviewed_biological_label = c("higher in SUS", "x")), "neuron_soma", member_map), "group-effect direction")
})

testthat::test_that("a full reviewed registry rejects a row marked reviewed but still pending", {
  r <- partial_fixture()
  testthat::expect_identical(wgcna_reviewed_rows_unfinished(r), c(FALSE, FALSE))
  r$reviewer[2] <- "pending review (mechanical carry-over, 2026-10-26)"
  testthat::expect_identical(wgcna_reviewed_rows_unfinished(r), c(FALSE, TRUE))
})

testthat::test_that("Stage 07 applies the overrides in reviewed mode only, for the declared datasets", {
  code <- paste(readLines(repo_path("analysis", "wgcna", "summarize_module_interpretation.R"), warn = FALSE),
                collapse = "\n")
  testthat::expect_match(code, 'identical(label_mode, "reviewed")', fixed = TRUE)
  testthat::expect_match(code, "ds %in% wgcna_partial_registry_datasets()", fixed = TRUE)
  testthat::expect_identical(wgcna_partial_registry_datasets(), c("neuron_soma", "neuron_neuropil"))
})

testthat::test_that("Stage 07 displays the reviewed labels of every partial-registry dataset", {
  n <- 0L
  for (ds in wgcna_partial_registry_datasets()) {
    cfg <- repo_path("config", "wgcna_labels", paste0(ds, ".csv"))
    lk_path <- current_wgcna("summarize_module_interpretation", ds, "WGCNA_final_label_lookup.csv")
    if (!file.exists(cfg) || !file.exists(lk_path)) next
    if (identical(wgcna_declared_label_mode(ds)$mode, "provisional")) next   # Stage 07 does not read it there
    reg <- utils::read.csv(cfg, stringsAsFactors = FALSE, encoding = "UTF-8")
    lk <- utils::read.csv(lk_path, stringsAsFactors = FALSE, encoding = "UTF-8")
    m <- lk[match(paste(reg$level, reg$entity_id), paste(lk$level, lk$entity_id)), , drop = FALSE]
    testthat::expect_false(anyNA(m$entity_id), info = ds)
    testthat::expect_identical(m$best_data_driven_label, reg$reviewed_biological_label, info = ds)
    testthat::expect_identical(m$label_confidence, reg$label_confidence, info = ds)
    testthat::expect_identical(m$final_plot_label, paste0(reg$entity_id, " \u00b7 ", reg$reviewed_short_label), info = ds)
    testthat::expect_true(all(startsWith(m$label_basis, "reviewed_registry=")), info = ds)
    status <- utils::read.csv(wgcna_label_status_manifest(ds), stringsAsFactors = FALSE)
    testthat::expect_identical(unique(status$label_review_status), "automatic_with_reviewed_overrides", info = ds)
    n <- n + 1L
  }
  if (!n) testthat::skip("no partial registry with a Stage 07 lookup in reviewed mode")
})
