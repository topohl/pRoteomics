source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "candidate_figure_utils.R"))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
have <- function(p) file.exists(p)
CAND <- cf_contract()
CANDF <- function(fig, ...) cf_output_paths(fig)
SHARED <- cf_shared_paths()

# =====================================================================
# 1-2. the candidate contract is separate, and its ids cannot collide
# =====================================================================

testthat::test_that("the candidate contract is a separate file from the canonical one", {
  testthat::expect_true(have(cf_contract_path()))
  testthat::expect_true(have(cf_canonical_contract_path()))
  testthat::expect_false(identical(normalizePath(cf_contract_path(), mustWork = FALSE),
                                   normalizePath(cf_canonical_contract_path(), mustWork = FALSE)))
  testthat::expect_identical(CAND$contract_version, cf_contract_version())
  testthat::expect_identical(CAND$contract_version, "manuscript_candidate_figures_v1")
  # and it declares itself as candidate-only
  testthat::expect_identical(CAND$status, "candidate_only_not_promoted")
  testthat::expect_identical(CAND$compares_against, "manuscript_figures_v2")

  # The canonical contract is untouched BY THIS LAYER. It has moved on for its
  # own reasons - Phase 5B promoted final_truth_v9, so it is now v3 with 2a-2h
  # and 3a-3i - but nothing here may be the cause of that, which is what the
  # reference check at the end of this file actually enforces.
  y <- yaml::read_yaml(cf_canonical_contract_path())
  testthat::expect_identical(y$contract_version,
                             "manuscript_figures_v3_final_truth_v9_promoted")
  ids2 <- vapply(y$figures[["02"]]$panels, function(p) as.character(p$id), character(1))
  ids3 <- vapply(y$figures[["03"]]$panels, function(p) as.character(p$id), character(1))
  testthat::expect_identical(ids2, paste0("2", letters[1:8]))
  testthat::expect_identical(ids3, paste0("3", letters[1:9]))
})

testthat::test_that("candidate panel ids cannot collide with canonical panel ids", {
  ids <- vapply(CAND$panels, function(p) as.character(p$id), character(1))
  canonical <- c(paste0("2", letters[1:6]), paste0("3", letters[1:5]))
  testthat::expect_length(intersect(ids, canonical), 0L)
  # no candidate id even LOOKS like a canonical one
  testthat::expect_false(any(grepl("^[23][a-f]$", ids)))
  testthat::expect_true(all(grepl("^[23](x|ref)_", ids)))
  testthat::expect_identical(anyDuplicated(ids), 0L)

  # the loader itself rejects a colliding id rather than trusting the file
  tmp <- tempfile(fileext = ".yml")
  yaml::write_yaml(list(contract_version = "x",
                        panels = list(list(id = "2c", candidate_figure = "02"))), tmp)
  testthat::expect_error(cf_contract(tmp), "collide with the canonical panel namespace")
  unlink(tmp)
})

# =====================================================================
# 3. canonical Figure 2/3 assets unchanged
# =====================================================================

testthat::test_that("the canonical figure layer is not modified by the candidate layer", {
  code_of <- function(f) {
    l <- readLines(f, warn = FALSE)
    paste(sub("#.*$", "", l), collapse = "\n")
  }
  # The candidate engine touches the canonical namespace in exactly one place,
  # cf_canonical_panel_source(), and only to READ an existing asset. Prove that
  # every canonical-namespace reference in the code is a read.
  eng <- code_of(repo_path("R", "candidate_figure_utils.R"))
  testthat::expect_false(grepl("manuscript_panels", eng, fixed = TRUE))
  canon_refs <- grep('figures", "manuscript"',
                     sub("#.*$", "", readLines(repo_path("R", "candidate_figure_utils.R"),
                                               warn = FALSE)),
                     fixed = TRUE, value = TRUE)
  testthat::expect_length(canon_refs, 1L)
  testthat::expect_match(canon_refs[1], "path_results", fixed = TRUE)
  # writes in the engine always target a cf_output_paths()/cf_shared_paths() dir
  testthat::expect_false(grepl('write_csv_safe\\(.*"manuscript"', eng))

  # the canonical engine is never sourced or called by the candidate layer
  for (f in c(repo_path("R", "candidate_figure_utils.R"),
              repo_path("R", "candidate_figure_panels.R"),
              repo_path("figures", "candidate_figure_02.R"),
              repo_path("figures", "candidate_figure_03.R"))) {
    s <- code_of(f)
    testthat::expect_false(grepl("source(repo_path(\"R\", \"manuscript_figure_utils.R\"))",
                                 s, fixed = TRUE), info = f)
    testthat::expect_false(grepl("manuscript_figure_main(", s, fixed = TRUE), info = f)
  }
  # canonical entry points are byte-unchanged relative to their committed form
  st <- suppressWarnings(system2("git", c("status", "--porcelain", "--",
    "figures/figure_contract.yml", "figures/figure_02.R", "figures/figure_03.R",
    "R/manuscript_figure_utils.R"), stdout = TRUE, stderr = FALSE))
  testthat::expect_length(st[nzchar(st)], 0L)
})

# =====================================================================
# 4. renderers cannot create new inference
# =====================================================================

testthat::test_that("candidate renderers cannot invoke model-fitting functions", {
  testthat::expect_true(cf_assert_no_model_fitting())
  toks <- cf_forbidden_tokens()
  testthat::expect_true(all(c("lmFit", "eBayes", "p.adjust", "gseGO",
                              "bootstrap_enrichment_test", "blockwiseModules") %in% toks))
  # the guard really fires: plant a violation in a temp file
  tmp <- tempfile(fileext = ".R")
  writeLines(c("f <- function(x) {", "  limma::lmFit(x)", "}"), tmp)
  testthat::expect_error(cf_assert_no_model_fitting(tmp),
                         "must not create new inference")
  unlink(tmp)
  # a file that only MENTIONS a forbidden call in a comment is fine
  tmp2 <- tempfile(fileext = ".R")
  writeLines(c("# we never call lmFit or p.adjust here", "g <- function(x) x"), tmp2)
  testthat::expect_true(cf_assert_no_model_fitting(tmp2))
  unlink(tmp2)
})

# =====================================================================
# 5-8. panels consume the canonical evidence and carry honest status
# =====================================================================

testthat::test_that("the GSEA panels use canonical ranked-GSEA output", {
  p <- Filter(function(x) identical(as.character(x$id), "3x_gsea_atlas_susres"),
              CAND$panels)[[1]]
  testthat::expect_match(as.character(p$primary_source),
                         "ontology_aware_gsea_theme_assignments_all_contrasts.csv",
                         fixed = TRUE)
  testthat::expect_true(any(grepl("manuscript_go_theme_registry.tsv",
                                  as.character(unlist(p$input_dependencies)), fixed = TRUE)))
  testthat::expect_true(have(repo_path(as.character(p$primary_source))))

  sd <- file.path(cf_output_paths("03")$source_data, "3x_gsea_atlas_susres_source_data.csv")
  testthat::skip_if_not(have(sd), "candidate figure 3 not built")
  z <- rd(sd)
  # a theme cell summarises constituent terms and says so
  testthat::expect_true(all(grepl("NOT itself a test", z$summary_basis)))
  testthat::expect_true(all(c("n_terms", "n_terms_FDR_supported", "median_NES",
                              "has_FDR_support", "theme_role") %in% names(z)))
  testthat::expect_setequal(unique(z$contrast), c("RES - CON", "SUS - CON", "SUS - RES"))
  # FDR support is counted from real constituent terms, never asserted
  testthat::expect_true(all(z$n_terms_FDR_supported <= z$n_terms))
  testthat::expect_identical(z$has_FDR_support, z$n_terms_FDR_supported > 0L)
  # qc_review themes are present but flagged, not silently mixed in as claims
  testthat::expect_true(any(z$theme_role == "qc_review"))
})

testthat::test_that("the WGCNA panel reuses the canonical Stage-07 effect source", {
  p <- Filter(function(x) identical(as.character(x$id), "3x_wgcna_annotated"),
              CAND$panels)[[1]]
  testthat::expect_match(as.character(p$primary_source),
                         "figure3b_stage07_effect_source.csv", fixed = TRUE)
  testthat::expect_identical(as.integer(p$expected_rows), 45L)

  sd <- file.path(cf_output_paths("03")$source_data, "3x_wgcna_annotated_source_data.csv")
  testthat::skip_if_not(have(sd), "candidate figure 3 not built")
  z <- rd(sd)
  testthat::expect_identical(nrow(z), 45L)
  # the effects are REUSED verbatim from the canonical source
  canon <- rd(repo_path(as.character(p$primary_source)))
  m <- match(paste(z$module_id, z$contrast), paste(canon$module_id, canon$contrast))
  testthat::expect_false(any(is.na(m)))
  testthat::expect_equal(z$estimate, canon$estimate[m])
  testthat::expect_equal(z$tier_specific_fdr, canon$tier_specific_fdr[m])
  # FDR support is an explicit separate symbol, not implied by colour
  testthat::expect_true("fdr_symbol" %in% names(z))
  testthat::expect_true(all(grepl("never drawn as statistical support",
                                  z$descriptive_geometry_note)))
})

testthat::test_that("the stress-identity panel uses the robustness-audited subsets", {
  p <- Filter(function(x) identical(as.character(x$id), "3x_stress_identity"),
              CAND$panels)[[1]]
  testthat::expect_match(as.character(p$primary_source),
                         "stress_identity_robustness_comparison.csv", fixed = TRUE)
  sd <- file.path(cf_output_paths("03")$source_data, "3x_stress_identity_source_data.csv")
  testthat::skip_if_not(have(sd), "candidate figure 3 not built")
  z <- rd(sd)
  testthat::expect_true(all(c("all_canonical_FDR_supported_hits", "fully_observed_hits",
                              "CA2_SLM_robustness_qualified") %in% z$subset))
  testthat::expect_true(all(nzchar(z$subset_definition)))
  # the headline is still reproduced faithfully
  h <- z[z$subset == "all_canonical_FDR_supported_hits", ]
  testthat::expect_identical(h$n_hits, 37L)
  testthat::expect_identical(h$effect_outside_baseline_affinity, 35L)
  # the language guard: not redistribution, and the rank-10 tail is de-emphasised
  testthat::expect_true(all(grepl("NOT protein redistribution", z$interpretation_note)))
  testthat::expect_true(all(grepl("did not survive", z$interpretation_note)))
  panel_src <- paste(readLines(repo_path("R", "candidate_figure_panels.R"), warn = FALSE),
                     collapse = "\n")
  testthat::expect_false(grepl("redistribution of proteins", panel_src, fixed = TRUE))
})

testthat::test_that("no QC-sensitive CA2-SLM count is presented as fully claimable", {
  sd <- file.path(cf_output_paths("03")$source_data, "3x_dap_status_source_data.csv")
  testthat::skip_if_not(have(sd), "candidate figure 3 not built")
  z <- rd(sd)
  testthat::expect_true("status" %in% names(z))
  testthat::expect_true(all(z$canonical_counts_preserved %in% TRUE))
  # CA2-SLM must be split by status, never shown as one claimable block
  ca2 <- z[z$spatial_unit == "CA2_slm", , drop = FALSE]
  testthat::expect_gt(nrow(ca2), 1L)
  testthat::expect_true(all(c("claimable", "not_claimable", "not_evaluable") %in%
                              ca2$status))
  testthat::expect_identical(sum(ca2$n_proteins), 28L)
  testthat::expect_identical(ca2$n_proteins[ca2$status == "claimable"], 6L)
  # unaudited units are labelled as such, never as claimable
  other <- z[z$spatial_unit != "CA2_slm", , drop = FALSE]
  testthat::expect_true(all(other$status == "not_audited"))
  # and the canonical total is preserved
  testthat::expect_identical(sum(z$n_proteins), 37L)
})

# =====================================================================
# 9-12. outputs, determinism, contact sheet, namespace
# =====================================================================

testthat::test_that("every candidate panel has source data and an SVG", {
  for (fig in c("02", "03")) {
    paths <- cf_output_paths(fig)
    testthat::skip_if_not(dir.exists(paths$panels), "candidate layer not built")
    ids <- vapply(cf_panels_for(CAND, fig), function(p) as.character(p$id), character(1))
    for (id in ids) {
      testthat::expect_true(have(file.path(paths$panels, paste0(id, ".svg"))), info = id)
      testthat::expect_true(
        have(file.path(paths$source_data, paste0(id, "_source_data.csv"))), info = id)
    }
  }
})

testthat::test_that("all declared whole-figure variants were built", {
  for (fig in c("02", "03")) {
    paths <- cf_output_paths(fig)
    testthat::skip_if_not(dir.exists(paths$assembled), "candidate layer not built")
    for (a in cf_assemblies_for(CAND, fig)) {
      testthat::expect_true(
        have(file.path(paths$assembled, paste0(as.character(a$name), ".svg"))),
        info = as.character(a$name))
      # every panel a variant references must be declared
      for (it in a$layout) {
        testthat::expect_silent(cf_panel_by_id(CAND, it$panel))
      }
    }
  }
})

testthat::test_that("candidate panel source data is deterministic", {
  # the panel source-data CSVs are pure functions of validated inputs, so the
  # same inputs must give the same bytes; this pins the renderers against
  # accidental randomness (sampling, jitter, unordered joins)
  for (fig in c("02", "03")) {
    sd <- cf_output_paths(fig)$source_data
    testthat::skip_if_not(dir.exists(sd), "candidate layer not built")
    files <- list.files(sd, pattern = "[.]csv$", full.names = TRUE)
    testthat::expect_gt(length(files), 0L)
    for (f in files) {
      a <- readLines(f, warn = FALSE)
      b <- readLines(f, warn = FALSE)
      testthat::expect_identical(a, b)
    }
  }
  # no renderer draws a random number
  src <- paste(sub("#.*$", "", readLines(repo_path("R", "candidate_figure_panels.R"),
                                         warn = FALSE)), collapse = "\n")
  for (tok in c("runif", "rnorm", "sample(", "jitter", "Sys.time", "Sys.Date")) {
    testthat::expect_false(grepl(tok, src, fixed = TRUE), info = tok)
  }
})

testthat::test_that("the contact sheet covers every declared variant", {
  inv <- file.path(SHARED$tables, "candidate_variant_inventory.csv")
  testthat::skip_if_not(have(inv), "contact sheet not built")
  v <- rd(inv)
  declared <- vapply(CAND$assemblies, function(a) as.character(a$name), character(1))
  testthat::expect_setequal(v$assembly, declared)
  testthat::expect_true(all(v$exists %in% TRUE))
  testthat::expect_true(have(file.path(SHARED$figures, "candidate_contact_sheet.svg")))

  pin <- rd(file.path(SHARED$tables, "candidate_panel_inventory.csv"))
  ids <- vapply(CAND$panels, function(p) as.character(p$id), character(1))
  testthat::expect_setequal(pin$candidate_panel_id, ids)
  for (col in c("candidate_panel_id", "candidate_figure", "scientific_question",
                "evidence_type", "inferential_status", "source_artifact",
                "overlaps_existing_panel", "likely_role", "caveat")) {
    testthat::expect_true(col %in% names(pin), info = col)
  }
  testthat::expect_true(all(pin$likely_role %in%
    c("likely_main", "candidate_replacement", "possible_extension",
      "extended_data", "supplementary_only", "undecided")))
  # no editorial decision is hard-coded as canonical
  testthat::expect_true(all(pin$promotion_status == "not_promoted"))
})

testthat::test_that("the candidate output namespace is isolated", {
  for (fig in c("02", "03")) {
    p <- cf_output_paths(fig)
    for (d in p) {
      testthat::expect_match(d, "manuscript_candidates", fixed = TRUE)
      testthat::expect_false(grepl("manuscript_panels", d, fixed = TRUE))
      testthat::expect_false(grepl(file.path("figures", "manuscript", "figure"), d,
                                   fixed = TRUE))
    }
  }
  # nothing was written into the canonical namespaces
  for (fig in c("figure_02", "figure_03")) {
    d <- path_results("figures", "manuscript", fig, "panels")
    if (!dir.exists(d)) next
    f <- list.files(d)
    testthat::expect_false(any(grepl("^[23](x|ref)_", f)))
  }
  # The whole candidate layer is removable: no canonical file references any
  # candidate artefact. (The canonical engine's own --allow-incomplete message
  # uses the English word "candidate"; what must be absent is a reference to
  # THIS layer's files or namespace.)
  for (f in c(repo_path("figures", "figure_02.R"), repo_path("figures", "figure_03.R"),
              repo_path("R", "manuscript_figure_utils.R"),
              repo_path("figures", "figure_contract.yml"))) {
    s <- paste(readLines(f, warn = FALSE), collapse = "\n")
    for (tok in c("figure_candidate_contract",
                  "candidate_figure_utils", "candidate_figure_panels",
                  "candidate_figure_02", "candidate_figure_03", "cf_build_figure")) {
      testthat::expect_false(grepl(tok, s, fixed = TRUE),
                             info = paste(basename(f), tok))
    }
    # The results root is shared by every non-canonical generation, so the bare
    # directory name cannot be banned outright any more: Phase 5B promoted
    # final_truth_v9, whose renderers still write there. What must stay true is
    # that the ONLY thing the canonical layer reaches into that root for is the
    # promoted generation. A reference to any other generation - or to this
    # comparison layer - would mean the canonical figures depend on something
    # removable.
    # A declared input_dependency may name an older generation, because that is
    # truthful provenance about what produced the artefact. What the canonical
    # layer READS - figure_source and primary_source - may only ever come from
    # the promoted generation.
    read_lines_ <- grep("^\\s*(figure_source|primary_source):", strsplit(s, "\n")[[1]],
                        value = TRUE)
    hits <- unlist(regmatches(read_lines_,
                              gregexpr("manuscript_candidates/[A-Za-z0-9_]+", read_lines_)))
    testthat::expect_true(all(hits == "manuscript_candidates/final_truth_v9"),
      info = paste(basename(f), "reads from a non-promoted candidate generation:",
                   paste(setdiff(unique(hits), "manuscript_candidates/final_truth_v9"),
                         collapse = ", ")))
  }
})
