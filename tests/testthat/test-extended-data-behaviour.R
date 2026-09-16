# Guards for the behavioural Extended Data figure.
#
# The figure renders two contracts that are computed in topohl/MMMSociability and
# frozen there. Three things therefore have to stay true, and none of them is
# checked by simply re-running the renderer: the imported bytes still match the
# bundle they claim to come from, this repository still computes nothing, and the
# panels still say only what the upstream contract permits.

ed_bridge_path <- function(...) {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  repo_path("manuscript", "figure1_bridge_mmmsociability", ...)
}

testthat::test_that("imported behavioural contracts match the recorded import hashes", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  testthat::skip_if_not_installed("digest")

  manifest <- utils::read.csv(
    repo_path("manuscript", "figure1_bridge_import_manifest.csv"),
    stringsAsFactors = FALSE
  )
  rows <- manifest[basename(manifest$imported_file) %in% c(
    "behavior_prediction_model_ladder.csv",
    "behavior_sex_effect_contract.csv"
  ), , drop = FALSE]

  testthat::expect_identical(nrow(rows), 2L)
  for (i in seq_len(nrow(rows))) {
    p <- repo_path(rows$imported_file[[i]])
    testthat::expect_true(file.exists(p))
    testthat::expect_identical(as.numeric(file.size(p)), as.numeric(rows$bytes[[i]]))
    testthat::expect_identical(
      digest::digest(file = p, algo = "sha256"),
      rows$sha256[[i]]
    )
    # An import that is not byte-identical is not an import.
    testthat::expect_true(as.logical(rows$byte_identical_to_source[[i]]))
  }
})

testthat::test_that("the model ladder still names one headline model and one baseline", {
  ladder <- utils::read.csv(ed_bridge_path("behavior_prediction_model_ladder.csv"),
                            stringsAsFactors = FALSE)

  testthat::expect_identical(nrow(ladder), 5L)
  testthat::expect_identical(
    ladder$model_id[ladder$current_status == "HEADLINE PRIMARY"], "movement_mean"
  )
  testthat::expect_identical(
    ladder$model_id[ladder$primary_or_exploratory == "reference_baseline"], "mean_only"
  )
  # All five are fitted on the same animals, or the panel is comparing
  # performance across different denominators without saying so.
  testthat::expect_identical(length(unique(ladder$n)), 1L)
  testthat::expect_identical(unique(ladder$n), 111L)

  # The permutation was run for the behaviour-only models only. If that ever
  # changes, "not run" stops being the right thing to print.
  has_p <- !is.na(ladder$permutation_p)
  testthat::expect_identical(
    sort(ladder$model_id[has_p]),
    sort(c("movement_mean", "primary_behavior_family"))
  )

  # Every interval must actually contain its own mean, or the bar and the tick
  # in panel a are inconsistent with each other.
  testthat::expect_true(all(
    ladder$repeated_cv_mean_r2 >= ladder$cv_r2_q025 &
      ladder$repeated_cv_mean_r2 <= ladder$cv_r2_q975
  ))
})

testthat::test_that("no feature-by-sex interaction is supported, which is what licenses panel c", {
  sexc <- utils::read.csv(ed_bridge_path("behavior_sex_effect_contract.csv"),
                          stringsAsFactors = FALSE)

  testthat::expect_identical(nrow(sexc), 3L)
  testthat::expect_true(all(sexc$classification == "FORMAL_INTERACTION_NOT_SUPPORTED"))
  testthat::expect_true(all(sexc$interaction_q_bh > 0.05))
  # Zero inside every interval is the visual claim panel b makes.
  testthat::expect_true(all(sexc$interaction_ci_low < 0 & sexc$interaction_ci_high > 0))
})

testthat::test_that("the Extended Data contract declares the panels the renderer writes", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "manuscript_figure_utils.R"))
  testthat::skip_if_not_installed("yaml")

  # Phase 6A split one behavioural Extended Data figure into two, because it was
  # carrying two different arguments: how the early window was measured and how
  # well it predicts, versus whether the association differs by sex. Panels 9b
  # and 9c keep the meaning they already had.
  cov <- manuscript_figure_contract("ED_behaviour_coverage")
  sec <- manuscript_figure_contract("ED_behaviour_secondary")

  testthat::expect_identical(
    vapply(cov$panels, function(x) as.character(x$id), character(1)), c("05a", "05b"))
  testthat::expect_identical(
    vapply(sec$panels, function(x) as.character(x$id), character(1)), c("09a", "09b", "09c"))
  testthat::expect_identical(as.integer(cov$extended_data_number), 5L)
  testthat::expect_identical(as.integer(sec$extended_data_number), 9L)
  testthat::expect_identical(as.character(cov$canonical_publication_id), "extended_data_05")
  testthat::expect_identical(as.character(sec$canonical_publication_id), "extended_data_09")

  for (fig in list(cov, sec)) {
    testthat::expect_identical(fig$contract_version,
                               "manuscript_figures_v3_final_truth_v9_promoted")
    testthat::expect_false(isTRUE(fig$rendering_repository_computes_statistics))
    testthat::expect_identical(as.numeric(fig$width_mm), 183)
    testthat::expect_identical(as.numeric(fig$height_mm), 118)
    for (panel in fig$panels) {
      testthat::expect_identical(as.character(panel$producer_script),
                                 "figures/extended_data_behaviour_panels.R")
      testthat::expect_true(file.exists(repo_path(as.character(panel$primary_source))))
    }
  }
})

testthat::test_that("the secondary-feature panel is honest about what it can draw", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  # The bundle exports per-animal values for Movement_mean only, so RMSSD and
  # entropy are shown as effect sizes. If a per-animal export ever appears this
  # test should be revisited deliberately rather than the panel silently changed.
  assoc <- utils::read.csv(ed_bridge_path("early_behavior_later_outcome_association.csv"),
                           stringsAsFactors = FALSE)
  testthat::expect_identical(nrow(assoc), 3L)
  testthat::expect_true(all(c("Movement_mean", "Movement_rmssd", "Entropy_acf1") %in%
                              assoc$predictor))
  testthat::expect_identical(unique(assoc$n), 111L)

  # The two secondary features and their adjudication, which the panel prints.
  r <- assoc[assoc$predictor == "Movement_rmssd", ]
  testthat::expect_lt(r$q_bh, 0.05)
  testthat::expect_lt(r$ci_high, 0)                      # interval excludes zero
  e <- assoc[assoc$predictor == "Entropy_acf1", ]
  testthat::expect_gt(e$q_bh, 0.05)
  testthat::expect_gt(e$ci_high, 0)                      # interval includes zero

  # No per-animal secondary-feature source exists anywhere in the bridge.
  f1c <- utils::read.csv(ed_bridge_path("source_data", "figure1c_movement_combz_source.csv"),
                         stringsAsFactors = FALSE)
  testthat::expect_false(any(grepl("rmssd|entropy", names(f1c), ignore.case = TRUE)))
})

testthat::test_that("early-window coverage is stated as the aggregate the bundle exports", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  tl <- utils::read.csv(ed_bridge_path("source_data", "figure1a_timeline_source.csv"),
                        stringsAsFactors = FALSE)
  cov <- tl[!is.na(tl$n_animals), , drop = FALSE][1, ]
  testthat::expect_identical(as.integer(cov$n_animals), 111L)
  testthat::expect_identical(as.integer(cov$n_animals_complete_slots), 50L)
  # 111 - 50 = 61 animals missing leading slots only.
  testthat::expect_identical(111L - as.integer(cov$n_animals_complete_slots), 61L)
  testthat::expect_equal(round(100 * as.numeric(cov$mean_coverage_fraction), 1), 98.6)
  testthat::expect_equal(round(100 * as.numeric(cov$min_coverage_fraction), 1), 94.4)
})

testthat::test_that("the Extended Data namespace resolves and still fails closed", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  source(repo_path("R", "output_namespace_utils.R"))

  paths <- output_namespace_manuscript_figure_paths("/tmp/root", "ED_behaviour_secondary")
  testthat::expect_true(grepl("extended_data_09$", paths$figures))
  testthat::expect_true(grepl("extended_data_09$", paths$source_data))

  # The numbered figures are untouched.
  testthat::expect_true(grepl(
    "figure_01$", output_namespace_manuscript_figure_paths("/tmp/root", "01")$figures
  ))
  # And an unrecognised ID is still rejected rather than silently creating a
  # namespace, which is the whole point of the allow-list.
  testthat::expect_error(
    output_namespace_manuscript_figure_paths("/tmp/root", "07"), "must be one of"
  )
  testthat::expect_error(
    output_namespace_manuscript_figure_paths("/tmp/root", "ED_nonsense"), "must be one of"
  )
})

testthat::test_that("the renderer computes no statistic", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  src <- readLines(repo_path("figures", "extended_data_behaviour_panels.R"), warn = FALSE)

  # Match calls, not prose: the file's own commentary names several of these
  # procedures in order to say that it does not perform them.
  code <- sub("#.*$", "", src)
  code <- paste(code, collapse = "\n")

  forbidden <- c("cor", "cor.test", "lm", "glm", "lmer", "t.test", "wilcox.test",
                 "p.adjust", "sample", "boot", "quantile", "aov", "anova",
                 "fisher.test", "chisq.test", "predict", "residuals")
  for (fn in forbidden) {
    pattern <- paste0("(^|[^[:alnum:]._$])", gsub(".", "\\.", fn, fixed = TRUE), "[[:space:]]*\\(")
    testthat::expect_false(
      grepl(pattern, code, perl = TRUE),
      info = paste0("renderer appears to call ", fn, "()")
    )
  }

  # It must read the frozen contracts and nothing else.
  testthat::expect_true(grepl("behavior_prediction_model_ladder.csv", code, fixed = TRUE))
  testthat::expect_true(grepl("behavior_sex_effect_contract.csv", code, fixed = TRUE))
  # The classification guard is the reason panel c is allowed to exist.
  testthat::expect_true(grepl("FORMAL_INTERACTION_NOT_SUPPORTED", code, fixed = TRUE))
})

testthat::test_that("prohibited sex wording appears on no rendered panel", {
  source(testthat::test_path("..", "..", "R", "paths.R"))

  panel_dir <- path_results("figures", "manuscript", "extended_data_behaviour_panels")
  testthat::skip_if_not(dir.exists(panel_dir), "panels not rendered in this environment")
  svgs <- list.files(panel_dir, pattern = "[.]svg$", full.names = TRUE)
  testthat::expect_gt(length(svgs), 0L)

  prohibited <- c("female-specific", "sex-specific", "stronger in females",
                  "driven by females")
  for (f in svgs) {
    text <- paste(readLines(f, warn = FALSE), collapse = " ")
    for (phrase in prohibited) {
      testthat::expect_false(
        grepl(phrase, text, ignore.case = TRUE, fixed = FALSE),
        info = paste0(basename(f), " renders prohibited wording: ", phrase)
      )
    }
  }
})

testthat::test_that("Results cites the Extended Data panels and no retired Figure 1 panel", {
  source(testthat::test_path("..", "..", "R", "paths.R"))
  draft <- readLines(repo_path("manuscript", "manuscript_draft.md"), warn = FALSE)
  text <- paste(draft, collapse = "\n")

  testthat::expect_true(grepl("Extended Data Fig. 9a", text, fixed = TRUE))
  testthat::expect_true(grepl("Extended Data Fig. 9b", text, fixed = TRUE))
  testthat::expect_true(grepl("Extended Data Fig. 9c", text, fixed = TRUE))

  # Figure 1 has four panels. Panels e and f were merged into d, so a reference
  # to either is a dangling pointer rather than a typo.
  testthat::expect_false(grepl("Fig. 1e", text, fixed = TRUE))
  testthat::expect_false(grepl("Fig. 1f", text, fixed = TRUE))
})
