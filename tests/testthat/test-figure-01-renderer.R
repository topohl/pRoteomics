# Guards for the manuscript Figure 1 renderer.
#
# Figure 1 is the one main figure whose science was computed in a different
# repository. The risk it carries is therefore not a wrong model but a wrong
# quote, or a panel that implies something the frozen analysis does not support.
# These guards pin the separation: the renderer may read only the frozen source
# data, it may print only numbers that the frozen statistics index resolves, and
# the panels may not carry the columns that would let them overreach.

repo <- function(...) file.path(testthat::test_path("..", ".."), ...)
SD <- repo("manuscript", "figure1_bridge_mmmsociability", "source_data")
PANELS <- repo("results", "figures", "manuscript", "figure_01_panels")
FIG <- repo("results", "figures", "manuscript", "figure_01")
RENDERER <- repo("figures", "figure_01_panels.R")
ENTRY <- repo("figures", "figure_01.R")
LEGEND <- repo("manuscript", "figure1_legend.md")
rd <- function(p) utils::read.csv(p, stringsAsFactors = FALSE)
have_src <- dir.exists(SD) && file.exists(file.path(SD, "figure1_panel_statistics.csv"))

test_that("the renderer reads only frozen Figure 1 source data", {
  src <- readLines(RENDERER, warn = FALSE)
  code <- src[!grepl("^\\s*#", src)]
  # It must not reach into the analysis layer, the raw data, or the upstream
  # repository. The behavioural analysis is owned elsewhere and is frozen.
  # Paths and repositories only. "RFID" is a legitimate word in a panel label,
  # so the list matches locations the renderer must not reach, not vocabulary.
  FORBIDDEN <- c("MMMSociability", "analysis_ready", "data/raw", "data/processed",
                 "[.]gct", "pg_matrix", "canonical/later_outcome_combz")
  for (f in FORBIDDEN)
    expect_equal(grep(f, code, value = TRUE), character(0),
                 info = paste("renderer reaches outside the frozen bridge:", f))
  # and it must not compute a statistic
  STATS <- c("\\blm\\(", "\\bglm\\(", "cor\\.test\\(", "\\bcor\\(", "p\\.adjust\\(",
             "\\bboot\\(", "replicate\\(", "geom_smooth\\(", "stat_smooth\\(",
             "\\bt\\.test\\(", "wilcox\\.test\\(")
  for (f in STATS)
    expect_equal(grep(f, code, value = TRUE), character(0),
                 info = paste("renderer computes a statistic:", f))
  # the only source directory it names is the frozen bridge
  expect_true(any(grepl("figure1_bridge_mmmsociability", code, fixed = TRUE)))
})

test_that("the entry point is a thin contract-driven assembler", {
  src <- readLines(ENTRY, warn = FALSE)
  code <- src[!grepl("^\\s*#", src)]
  code <- code[nzchar(trimws(code))]
  # Same shape as figure_02.R and figure_03.R: source paths, source the utils,
  # call the shared main. Anything longer means Figure 1 grew its own pipeline.
  expect_lte(length(code), 5L)
  expect_true(any(grepl('manuscript_figure_main("01")', code, fixed = TRUE)))
})

test_that("Figure 1 has no dependency on behavioural analysis code", {
  # No .R file from the upstream repository may exist here, and the bridge must
  # remain a CSV-only evidence interface.
  expect_equal(list.files(SD, pattern = "[.][Rr]$"), character(0))
  upstream <- c("09_early_prediction_model_ladder", "build_later_outcome_combz",
                "27_build_behavior_main_figure", "build_figure1_panel_source_data")
  here <- list.files(repo("."), pattern = "[.][Rr]$", recursive = TRUE)
  here <- here[!grepl("^(\\.git|renv)/", here)]
  for (u in upstream) expect_false(any(grepl(u, basename(here), fixed = TRUE)))
})

test_that("every required panel source file is present and hash-matched", {
  skip_if_not(have_src, "figure 1 source data not imported")
  man <- rd(file.path(SD, "figure1_panel_source_manifest.csv"))
  expect_gt(nrow(man), 0)
  skip_if_not(requireNamespace("digest", quietly = TRUE), "digest unavailable")
  for (i in seq_len(nrow(man))) {
    p <- file.path(SD, man$file[i])
    expect_true(file.exists(p), info = man$file[i])
    expect_equal(file.size(p), man$bytes[i])
    expect_equal(digest::digest(normalizePath(p), algo = "sha256", file = TRUE),
                 man$sha256[i])
  }
  # the freeze must not describe itself as authoritative if it was not
  expect_true(all(man$creation_state == "authoritative"))
  expect_true(all(man$source_worktree_state == "clean"))
})

test_that("the animal-level panels carry one row per animal", {
  skip_if_not(have_src, "figure 1 source data not imported")
  a <- rd(file.path(SD, "figure1c_movement_combz_source.csv"))
  d <- rd(file.path(SD, "figure1d_loao_predictions_source.csv"))
  b <- rd(file.path(SD, "figure1b_combz_classification_source.csv"))
  expect_equal(nrow(a), 111L)
  expect_equal(length(unique(a$AnimalID)), 111L)
  expect_equal(nrow(d), 111L)
  expect_equal(anyDuplicated(d$AnimalID), 0L)
  expect_equal(anyDuplicated(b$AnimalID), 0L)
  expect_setequal(a$AnimalID, d$AnimalID)
})

test_that("canonical values reproduce from the frozen source tables", {
  skip_if_not(have_src, "figure 1 source data not imported")
  st <- rd(file.path(SD, "figure1_panel_statistics.csv"))
  sv <- function(panel, stat) {
    v <- st$value[st$figure_panel == panel & st$statistic == stat]
    expect_length(v, 1L)
    as.numeric(v)
  }
  expect_equal(sv("1c", "Spearman rho"), -0.3902639142724438, tolerance = 1e-12)
  expect_equal(sv("1c", "95% CI lower"), -0.547434593760078, tolerance = 1e-12)
  expect_equal(sv("1c", "95% CI upper"), -0.20917806610402842, tolerance = 1e-12)
  expect_equal(sv("1c", "BH q"), 6.8776554749896e-05, tolerance = 1e-12)
  expect_equal(sv("1d", "LOAO R2"), 0.15939455855319962, tolerance = 1e-12)
  expect_equal(sv("1d", "intercept-only baseline R2"), -0.01826446280991756,
               tolerance = 1e-12)
  expect_equal(sv("1f", "repeated grouped CV mean R2"), 0.15582272536971034,
               tolerance = 1e-12)
  expect_equal(sv("1f", "repeated CV seed"), 521)
  expect_equal(sv("1b", "male susceptibility threshold"), -0.4366416981697792,
               tolerance = 1e-12)
  expect_equal(sv("1b", "female susceptibility threshold"), -0.22239084402294293,
               tolerance = 1e-12)

  # p is recomputed from the persisted draws, never taken on trust
  e <- rd(file.path(SD, "figure1e_permutation_source.csv"))
  obs <- e$performance_value[e$row_role == "observed_statistic"]
  nulls <- e$performance_value[e$row_role == "permutation_draw"]
  expect_length(nulls, 1000L)
  expect_equal(obs, 0.15939455855319962, tolerance = 1e-12)
  expect_equal(sum(nulls >= obs), 0L)
  expect_equal((sum(nulls >= obs) + 1) / (length(nulls) + 1), 1 / 1001)
})

test_that("no panel can imply more than the analysis supports", {
  skip_if_not(have_src, "figure 1 source data not imported")
  d <- rd(file.path(SD, "figure1d_loao_predictions_source.csv"))
  b <- rd(file.path(SD, "figure1b_combz_classification_source.csv"))
  # the prediction panel carries the headline behaviour-only model and no other
  expect_equal(unique(d$model_id), "movement_mean")
  # the components that construct CombZ are not plottable in the definition panel
  COMPONENTS <- c("NOR", "sucrose_pref", "weight_dev", "delta_cort",
                  "adrenal_weight", "spleen_weight")
  expect_false(any(COMPONENTS %in% names(b)))
  # the unsupported secondary features are absent from the main figure entirely
  for (f in c("figure1c_movement_combz_source.csv",
              "figure1d_loao_predictions_source.csv")) {
    nm <- names(rd(file.path(SD, f)))
    expect_false(any(grepl("rmssd|entropy|gamm|hmm", nm, ignore.case = TRUE)))
  }
})

test_that("the rendered panels and assembly are vector, not rasterised", {
  skip_if_not(dir.exists(PANELS), "figure 1 panels not rendered")
  for (id in c("a", "b", "c", "d")) {
    p <- file.path(PANELS, paste0("figure_01", id, ".svg"))
    expect_true(file.exists(p))
    expect_gt(file.size(p), 0)
    txt <- readLines(p, warn = FALSE)
    expect_true(any(grepl("</svg>", txt, fixed = TRUE)))
    expect_true(any(grepl("<text", txt, fixed = TRUE)))
    expect_false(any(grepl("<image", txt, fixed = TRUE)))
  }
  skip_if_not(dir.exists(file.path(FIG, "assembled")), "figure 1 not assembled")
  asm <- file.path(FIG, "assembled", "figure_01.svg")
  expect_true(file.exists(asm))
  txt <- paste(readLines(asm, warn = FALSE), collapse = "")
  # The assembler nests each panel as a data:image/svg+xml payload. That is
  # still vector. What must never appear is a raster payload.
  expect_false(grepl("data:image/png", txt, fixed = TRUE))
  expect_false(grepl("data:image/jpeg", txt, fixed = TRUE))
  expect_true(grepl("data:image/svg+xml", txt, fixed = TRUE))
  for (f in c("figure_01.pdf", "figure_01.png")) {
    p <- file.path(FIG, "assembled", f)
    expect_true(file.exists(p))
    expect_gt(file.size(p), 0)
  }
})

test_that("the figure contract registers Figure 1 with its scientific rules", {
  y <- yaml::read_yaml(repo("figures", "figure_contract.yml"))
  e <- y$figures[["01"]]
  expect_false(is.null(e))
  expect_true(isTRUE(e$is_numbered_manuscript_figure))
  expect_false(isTRUE(e$rendering_repository_computes_statistics))
  expect_equal(length(e$panels), 4L)
  for (p in e$panels) {
    expect_equal(p$biological_unit, "AnimalID")
    expect_equal(as.character(p$producer_script), "figures/figure_01_panels.R")
    expect_true(grepl("^manuscript/figure1_bridge_mmmsociability/source_data/",
                      as.character(p$primary_source)))
  }
  # the row counts that make the figure checkable are pinned in the contract
  rows <- vapply(e$panels, function(p) as.integer(p$expected_rows %||% NA_integer_),
                 integer(1))
  expect_equal(rows[[3]], 111L)
  expect_equal(rows[[4]], 111L)
})

test_that("the legend states what it must and avoids what it may not", {
  skip_if_not(file.exists(LEGEND), "legend absent")
  ln <- readLines(LEGEND, warn = FALSE)
  d <- paste(ln, collapse = " ")
  # a line that prohibits a phrase necessarily contains it; exempt those, as the
  # repository's other wording guards do
  DENIAL <- paste0("never|must not|do not |does not|cannot|prohibited|absent|",
                   "instead of|rather than|avoid|not licensed|\\bNOT\\b|",
                   "is not |are not |unsupported|no ")
  BANNED <- c("predicts susceptibility", "predicts resilience", "movement-only",
              "female-specific", "sex-specific", "external validation",
              "externally validated", "independent validation", "biomarker",
              "strong prediction", "highly accurate")
  for (b in BANNED) {
    hits <- grep(b, ln, ignore.case = TRUE, value = TRUE)
    hits <- hits[!grepl(DENIAL, hits, perl = TRUE)]
    expect_equal(hits, character(0), info = paste("legend uses:", b))
  }
  # and it must define the things a reader cannot infer
  for (must in c("P25", "movement-mean model", "internal", "72",
                 "not a confidence interval", "guaranteed")) {
    expect_true(grepl(must, d, fixed = TRUE), info = paste("legend omits:", must))
  }
  # the sign contract, which is the one error that would invert the result
  expect_true(grepl("[Hh]igher `CombZ` indicates a more resilient-like", d))
})

test_that("the superseded fifth panel leaves no stale artefact", {
  # Panel e's content is now the right half of panel d. A leftover figure_01e.svg
  # would still look canonical to anyone browsing the output tree, and the
  # contract no longer declares it.
  skip_if_not(dir.exists(PANELS), "figure 1 panels not rendered")
  expect_false(file.exists(file.path(PANELS, "figure_01e.svg")))
  expect_false(file.exists(file.path(FIG, "panels", "figure_01e.svg")))
  expect_false(file.exists(repo("results", "source_data", "manuscript",
                               "figure_01", "figure_01e_source_data.csv")))
  y <- yaml::read_yaml(repo("figures", "figure_contract.yml"))
  ids <- vapply(y$figures[["01"]]$panels, function(p) as.character(p$id),
                character(1))
  expect_false("1e" %in% ids)
  expect_setequal(ids, c("1a", "1b", "1c", "1d"))
})

test_that("Figure 1 uses the same group palette as Figures 2 and 3", {
  # The repository carries two group palettes. R/utilities/plotting_nature.R is the one
  # Figures 2 and 3 actually render with, and it is byte-identical to
  # MMM_GROUP_COLOURS upstream. config/manuscript_palette.yml declares a
  # different set that no numbered figure uses; Figure 1 previously obeyed it.
  source(repo("R", "plotting_nature.R"))
  expect_equal(unname(NATURE_SEMANTIC_PALETTES$group[c("CON", "RES", "SUS")]),
               c("#3E3C6F", "#C6C3BB", "#E63A48"))
  src <- readLines(RENDERER, warn = FALSE)
  code <- src[!grepl("^[[:space:]]*#", src)]
  expect_true(any(grepl("NATURE_SEMANTIC_PALETTES$group", code, fixed = TRUE)))
  expect_false(any(grepl("manuscript_palette.yml", code, fixed = TRUE)))

  # and the rendered panels must actually carry those inks
  skip_if_not(dir.exists(PANELS), "figure 1 panels not rendered")
  ink <- paste(readLines(file.path(PANELS, "figure_01d.svg"), warn = FALSE),
               collapse = "")
  for (h in c("#3E3C6F", "#C6C3BB", "#E63A48"))
    expect_true(grepl(h, ink, fixed = TRUE), info = paste("missing group ink", h))
})

test_that("the association panel is not coloured by outcome group", {
  # Susceptible animals sit low on the CombZ axis and resilient animals high by
  # construction, so colouring panel c by group would let the correlation read
  # as group separation. It is deliberately single-colour.
  skip_if_not(dir.exists(PANELS), "figure 1 panels not rendered")
  ink <- paste(readLines(file.path(PANELS, "figure_01c.svg"), warn = FALSE),
               collapse = "")
  for (h in c("#3E3C6F", "#E63A48"))
    expect_false(grepl(h, ink, fixed = TRUE),
                 info = paste("panel c carries group ink", h))
})
