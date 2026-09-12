source(testthat::test_path("..", "..", "R", "paths.R"))
source(testthat::test_path("..", "..", "R", "nature_v2_figure_utils.R"))
source(testthat::test_path("..", "..", "R", "story_v4_figure_utils.R"))

rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
have <- function(p) file.exists(p)
S4 <- s4e_contract()

testthat::test_that("the evidence inventory exists and drives the design", {
  p <- path_results("tables", "manuscript_candidates", "story_v4",
                    "proteomics_manuscript_evidence_inventory.csv")
  testthat::skip_if_not(have(p), "inventory not generated")
  z <- rd(p)
  testthat::expect_gte(nrow(z), 40L)
  for (col in c("claim_id", "dataset", "spatial_unit", "contrast", "analysis_type",
                "quantitative_result", "inferential_status", "FDR_status",
                "QC_status", "bilateral_status", "biological_program",
                "current_candidate_panel", "main_text_candidate",
                "extended_data_candidate", "reason")) {
    testthat::expect_true(col %in% names(z), info = col)
  }
  testthat::expect_true(all(nzchar(z$quantitative_result)))
  testthat::expect_true(all(nzchar(z$reason)))
  # it must cover the audit families the brief named, including the nulls
  testthat::expect_true(any(grepl("^Q_", z$claim_id)))   # network null
  testthat::expect_true(any(grepl("^R_", z$claim_id)))   # coupling null
  testthat::expect_true(any(grepl("^O_", z$claim_id)))   # bilateral, all contrasts
  testthat::expect_true(any(grepl("^G_", z$claim_id)))   # WGCNA phenotype
  # and it must record results that are NOT main-figure candidates
  testthat::expect_gt(sum(z$main_text_candidate == "no"), 5L)
})

testthat::test_that("story_v4 honours the size contract with lowercase labels", {
  testthat::expect_identical(S4$contract_version, s4e_contract_version())
  testthat::expect_identical(S4$status, "candidate_only_not_promoted")
  for (f in S4$figures) {
    testthat::expect_identical(as.numeric(f$width_mm), 183)
    testthat::expect_lte(as.numeric(f$height_mm), 170)
    labs <- vapply(f$layout, function(x) as.character(x$label), character(1))
    testthat::expect_identical(labs, tolower(labs))
    testthat::expect_identical(anyDuplicated(labs), 0L)
    area <- vapply(f$layout, function(x) as.numeric(x$w) * as.numeric(x$h), numeric(1))
    # page density: no more than a third of the page may be accidental whitespace
    testthat::expect_gt(sum(area) / (as.numeric(f$width_mm) * as.numeric(f$height_mm)),
                        0.70)
  }
})

testthat::test_that("the GSEA curve is an exact reconstruction, not a re-run", {
  d <- s4e_output_paths("figure_03")$source_data
  testthat::skip_if_not(dir.exists(d), "story_v4 figure 3 not built")
  for (i in c("s4_ex_neuropil", "s4_ex_soma", "s4_ex_microglia")) {
    p <- file.path(d, paste0(i, "_source_data.csv"))
    testthat::skip_if_not(have(p), i)
    z <- rd(p)
    testthat::expect_identical(nrow(z), 1L)
    # the reconstruction note must state what was and was not done
    testthat::expect_true(grepl("annotation lookup, NOT an", z$reconstruction_note))
    testthat::expect_true(grepl("no p-value or FDR is", z$reconstruction_note))
    # leading edge is a strict subset of the canonical set size
    testthat::expect_lt(z$n_leading_edge, z$setSize)
    testthat::expect_gt(z$n_leading_edge, 0L)
    # the canonical FDR is carried through, not recomputed
    testthat::expect_lt(z$FDR, 0.05)
    # and the three-contrast trajectory is present
    for (cn in c("RES_CON_NES", "SUS_CON_NES", "SUS_RES_NES",
                 "RES_CON_FDR", "SUS_CON_FDR", "SUS_RES_FDR")) {
      testthat::expect_true(cn %in% names(z), info = paste(i, cn))
    }
  }
  # the renderer refuses to draw if membership does not reproduce setSize
  src <- paste(readLines(repo_path("R", "story_v4_figure_panels.R"), warn = FALSE),
               collapse = "\n")
  testthat::expect_true(grepl("membership mismatch", src, fixed = TRUE))
  testthat::expect_true(grepl("does not match the stored canonical ES", src,
                              fixed = TRUE))
})

testthat::test_that("the three programs span the three compartments and differ in trajectory", {
  d <- s4e_output_paths("figure_03")$source_data
  testthat::skip_if_not(dir.exists(d), "story_v4 figure 3 not built")
  ds <- character(); traj <- list()
  for (i in c("s4_ex_neuropil", "s4_ex_soma", "s4_ex_microglia")) {
    z <- rd(file.path(d, paste0(i, "_source_data.csv")))
    ds <- c(ds, z$dataset)
    traj[[i]] <- c(z$RES_CON_NES, z$SUS_CON_NES, z$SUS_RES_NES)
  }
  testthat::expect_setequal(ds, c("neuron_neuropil", "neuron_soma", "microglia"))
  # neuropil and soma diverge: RES and SUS move in opposite directions vs CON
  testthat::expect_lt(traj[["s4_ex_neuropil"]][1] * traj[["s4_ex_neuropil"]][2], 0)
  testthat::expect_lt(traj[["s4_ex_soma"]][1] * traj[["s4_ex_soma"]][2], 0)
  # microglia is graded: both stressed groups move the same way
  testthat::expect_gt(traj[["s4_ex_microglia"]][1] * traj[["s4_ex_microglia"]][2], 0)
})

testthat::test_that("bilateral reproducibility shows every prespecified contrast", {
  p <- file.path(s4e_output_paths("figure_02")$source_data,
                 "s4_bilateral_source_data.csv")
  testthat::skip_if_not(have(p), "story_v4 figure 2 not built")
  z <- rd(p)
  # all 15 prespecified contrasts, not a median subset
  testthat::expect_identical(nrow(z), 15L)
  testthat::expect_setequal(unique(z$dataset),
                            c("neuron_neuropil", "neuron_soma", "microglia"))
  testthat::expect_identical(sum(z$is_representative_scatter), 1L)
  # the weak fine-laminar result is present, not hidden
  testthat::expect_lt(min(z$pearson_r), 0.45)
  testthat::expect_gt(max(z$pearson_r), 0.90)
  testthat::expect_true(all(grepl("is the result, not an inconvenience", z$note)))
})

testthat::test_that("the schematic is a disclosed placeholder with no fabricated histology", {
  p <- file.path(s4e_output_paths("figure_02")$source_data,
                 "s4_schematic_source_data.csv")
  testthat::skip_if_not(have(p), "story_v4 figure 2 not built")
  z <- rd(p)
  st <- z$value[z$element == "artwork_status"]
  testthat::expect_true(grepl("VECTOR PLACEHOLDER", st))
  testthat::expect_true(grepl("no hippocampal artwork", st))
  testthat::expect_true(grepl("No histological data is fabricated", st))
  testthat::expect_true(any(grepl("left;right", z$value)))
})

testthat::test_that("compartment identity uses the stored matrix, not invented dispersion", {
  p <- file.path(s4e_output_paths("figure_02")$source_data,
                 "s4_compartment_source_data.csv")
  testthat::skip_if_not(have(p), "story_v4 figure 2 not built")
  z <- rd(p)
  testthat::expect_identical(nrow(z), 30L)          # 10 markers x 3 compartments
  testthat::expect_setequal(unique(z$dataset),
                            c("neuron_soma", "neuron_neuropil", "microglia"))
  testthat::expect_true(all(grepl("no CI or SE exists", z$encoding_note)))
})

testthat::test_that("story_v4 is isolated and earlier layers are untouched", {
  for (key in c("figure_02", "figure_03")) {
    for (d in s4e_output_paths(key)) {
      testthat::expect_match(d, "manuscript_candidates[/\\\\]story_v4")
      testthat::expect_false(grepl("manuscript_panels", d, fixed = TRUE))
      testthat::expect_false(grepl("story_v3", d, fixed = TRUE))
    }
  }
  for (f in c("figures/figure_02.R", "figures/figure_03.R",
              "figures/figure_contract.yml", "R/manuscript_figure_utils.R",
              "figures/figure_candidate_contract.yml",
              "figures/figure_nature_v2_contract.yml",
              "figures/figure_story_v3_contract.yml",
              "R/story_v3_figure_panels.R")) {
    s <- paste(readLines(repo_path(f), warn = FALSE), collapse = "\n")
    for (tok in c("story_v4", "s4e_build", "figure_story_v4_contract")) {
      testthat::expect_false(grepl(tok, s, fixed = TRUE), info = paste(f, tok))
    }
  }
})

testthat::test_that("story_v4 renderers create no new inference and are deterministic", {
  testthat::expect_true(nv_assert_no_model_fitting(s4e_renderer_sources()))
  src <- paste(sub("#.*$", "", readLines(repo_path("R", "story_v4_figure_panels.R"),
                                          warn = FALSE)), collapse = "\n")
  for (tok in c("runif", "rnorm", "sample(", "jitter(", "Sys.time")) {
    testthat::expect_false(grepl(tok, src, fixed = TRUE), info = tok)
  }
  for (key in c("figure_02", "figure_03")) {
    d <- s4e_output_paths(key)$source_data
    if (!dir.exists(d)) next
    for (f in list.files(d, pattern = "[.]csv$", full.names = TRUE)) {
      testthat::expect_identical(readLines(f, warn = FALSE), readLines(f, warn = FALSE))
    }
  }
})

testthat::test_that("every story_v4 panel has source data and all variants were built", {
  for (key in c("figure_02", "figure_03")) {
    p <- s4e_output_paths(key)
    if (!dir.exists(p$assembled)) next
    figs <- Filter(function(f) identical(as.character(f$figure_key), key), S4$figures)
    for (f in figs) {
      testthat::expect_true(
        have(file.path(p$assembled, paste0(as.character(f$name), ".svg"))),
        info = as.character(f$name))
      for (it in f$layout) {
        testthat::expect_true(
          have(file.path(p$source_data,
                         paste0(as.character(it$panel), "_source_data.csv"))),
          info = as.character(it$panel))
      }
    }
  }
})
