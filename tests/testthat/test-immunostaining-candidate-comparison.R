# Guards for the manuscript-supporting immunostaining candidate comparison.
#
# The figure exists to compare three already-nominated proteins; the risk it
# carries is not a wrong model but a wrong copy. These tests therefore check
# identity and provenance: that the three requested proteins are the three
# present, that the mapping is unambiguous, that the animal-level table is free
# of pseudoreplication, and that every displayed effect equals its canonical
# source exactly.

repo <- function(...) file.path(testthat::test_path("..", ".."), ...)
SRC <- repo("results", "source_data", "10_biological_integration",
            "immunostaining_candidate_comparison")
FIG <- repo("results", "figures", "manuscript", "immunostaining_candidates")
have_src <- dir.exists(SRC) &&
  file.exists(file.path(SRC, "immunostaining_candidates_mapping.csv"))
rd <- function(f) utils::read.csv(file.path(SRC, f), stringsAsFactors = FALSE)

test_that("exactly the three requested proteins are present and unambiguous", {
  skip_if_not(have_src, "candidate source data not built")
  m <- rd("immunostaining_candidates_mapping.csv")
  expect_equal(nrow(m), 3L)
  expect_setequal(toupper(m$GeneSymbol), c("OGA", "SLC22A23", "ANXA2"))
  # requested accessions must match what the canonical mapping resolved
  expect_equal(m$UniProt[toupper(m$GeneSymbol) == "OGA"], "Q9EQQ9")
  expect_equal(m$UniProt[toupper(m$GeneSymbol) == "ANXA2"], "P07356")
  # one protein group each, and each a single-accession single-gene group
  expect_equal(anyDuplicated(m$ProteinGroupID), 0L)
  expect_true(all(m$ambiguity_class == "single_accession_single_gene"))
  expect_true(all(as.logical(m$protein_level_claim_allowed)))
})

test_that("the animal-level table carries no pseudoreplication", {
  skip_if_not(have_src, "candidate source data not built")
  a <- rd("immunostaining_candidates_animal_abundance.csv")
  # exactly one row per protein x animal x spatial unit: no hemisphere rows, no
  # technical replicates, no duplicated acquisitions
  key <- paste(a$ProteinGroupID, a$AnimalID, a$spatial_unit, sep = "|")
  expect_equal(anyDuplicated(key), 0L)
  expect_equal(nrow(a), 3L * 9L * 10L)
  expect_equal(length(unique(a$AnimalID)), 9L)
  expect_setequal(unique(a$ExpGroup), c("CON", "RES", "SUS"))
  # three animals per group, and the biological unit is declared
  per_group <- tapply(a$AnimalID, a$ExpGroup, function(x) length(unique(x)))
  expect_true(all(per_group == 3L))
  expect_true(all(a$biological_unit == "animal"))
  expect_false(any(grepl("hemisphere|_L_|_R_", a$sample_id, ignore.case = TRUE)))
})

test_that("displayed SUS-RES effects equal the canonical DA output exactly", {
  skip_if_not(have_src, "candidate source data not built")
  e <- rd("immunostaining_candidates_sus_res_effects.csv")
  m <- rd("immunostaining_candidates_mapping.csv")
  da_dir <- repo("data", "processed", "02_id_mapping", "mapped",
                 "neuron_neuropil", "forward", "per_file")
  skip_if_not(dir.exists(da_dir), "canonical DA inputs not present")
  for (u in unique(e$spatial_unit)) {
    key <- tolower(gsub("_", "", u))
    f <- file.path(da_dir, sprintf("%ssus_%sres.csv", key, key))
    skip_if_not(file.exists(f), paste("missing canonical file for", u))
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    d <- d[match(m$ProteinGroupID, d$ProteinGroupID), , drop = FALSE]
    ee <- e[e$spatial_unit == u, , drop = FALSE]
    ee <- ee[match(m$ProteinGroupID, ee$ProteinGroupID), , drop = FALSE]
    expect_equal(ee$log2FC, d$log2fc, tolerance = 0)
    expect_equal(ee$p_value, d$pval, tolerance = 0)
    expect_equal(ee$fdr_bh, d$padj, tolerance = 0)
  }
  expect_true(all(e$contrast == "SUS - RES"))
  expect_equal(e$fdr_supported, e$fdr_bh <= 0.05)
})

test_that("all three proteins share one identical spatial ordering", {
  skip_if_not(have_src, "candidate source data not built")
  a <- rd("immunostaining_candidates_animal_abundance.csv")
  e <- rd("immunostaining_candidates_sus_res_effects.csv")
  ord <- function(d) {
    o <- unique(d[order(d$spatial_order), c("spatial_unit", "spatial_order")])
    o$spatial_unit
  }
  per_protein <- lapply(split(a, a$GeneSymbol), ord)
  expect_equal(length(unique(per_protein)), 1L)
  # the two panels must use the same order as each other
  expect_equal(ord(a), ord(e))
  # and it must be the canonical region-major order from the shared contract
  expect_equal(ord(a),
               c("CA1_so", "CA1_sr", "CA1_slm", "CA2_so", "CA2_sr", "CA2_slm",
                 "CA3_so", "CA3_sr", "DG_mo", "DG_po"))
})

test_that("the renderer depends only on the prepared source data", {
  p <- repo("figures", "manuscript_supporting_immunostaining_candidates.R")
  skip_if_not(file.exists(p), "renderer absent")
  src <- readLines(p, warn = FALSE)
  # it must not reach back past the frozen source-data layer
  forbidden <- "data/processed|data/raw|\\.gct|protigy|per_file|lmFit|eBayes|p\\.adjust"
  offending <- grep(forbidden, src, value = TRUE)
  # comments may name those paths when explaining what is NOT done
  offending <- offending[!grepl("^\\s*#", offending)]
  expect_equal(offending, character(0))
})

test_that("the figure is registered as manuscript-supporting, not a numbered figure", {
  y <- yaml::read_yaml(repo("figures", "figure_contract.yml"))
  entry <- y$figures$S_immunostaining_candidates
  expect_false(is.null(entry))
  expect_false(isTRUE(entry$is_numbered_manuscript_figure))
  expect_equal(length(entry$panels), 2L)
  for (p in entry$panels) {
    expect_equal(p$biological_unit, "animal")
    expect_true(grepl("animal_level", p$hemisphere_handling))
    expect_true(all(grepl("^results/source_data/", unlist(p$input_dependencies))))
  }
})

test_that("the vector outputs exist, are non-empty and are not rasterised", {
  skip_if_not(dir.exists(FIG), "figure not rendered")
  svgs <- c("immunostaining_candidates_spatial_comparison.svg",
            "immunostaining_candidates_abundance.svg",
            "immunostaining_candidates_sus_res_log2fc.svg")
  for (s in svgs) {
    p <- file.path(FIG, s)
    expect_true(file.exists(p))
    expect_gt(file.size(p), 0)
    txt <- readLines(p, warn = FALSE)
    expect_true(any(grepl("</svg>", txt, fixed = TRUE)))
    # live text, no embedded bitmap
    expect_true(any(grepl("<text", txt, fixed = TRUE)))
    expect_false(any(grepl("<image", txt, fixed = TRUE)))
  }
})

# ---------------------------------------------- the ten-candidate screen

PSRC <- repo("results", "source_data", "10_biological_integration",
             "immunostaining_candidate_panel")
PFIG <- repo("results", "figures", "manuscript", "immunostaining_candidate_panel")
have_panel <- dir.exists(PSRC) &&
  file.exists(file.path(PSRC, "candidate_panel_selection.csv"))
prd <- function(f) utils::read.csv(file.path(PSRC, f), stringsAsFactors = FALSE)

test_that("the screen returns ten candidates in declared evidence classes", {
  skip_if_not(have_panel, "candidate panel not built")
  s <- prd("candidate_panel_selection.csv")
  expect_equal(nrow(s), 10L)
  expect_equal(anyDuplicated(s$ProteinGroupID), 0L)
  # the three already nominated must not reappear
  expect_false(any(toupper(s$GeneSymbol) %in% c("OGA", "SLC22A23", "ANXA2")))
  # every candidate carries one of the three declared classes
  expect_true(all(s$evidence_class %in% c(
    "1_FDR_supported_outside_CA2SLM",
    "2_FDR_supported_CA2SLM_robustness_qualified",
    "3_no_FDR_support_effect_and_consistency_only")))
  # and every one is unambiguously targetable
  expect_true(all(s$ambiguity_class == "single_accession_single_gene"))
  expect_true(all(as.logical(s$protein_level_claim_allowed)))
})

test_that("no CA2-SLM-only candidate enters without robustness qualification", {
  skip_if_not(have_panel, "candidate panel not built")
  s <- prd("candidate_panel_selection.csv")
  # class 2 means the only FDR support is CA2-SLM, so it must have qualified
  c2 <- s[s$evidence_class == "2_FDR_supported_CA2SLM_robustness_qualified", ]
  expect_true(all(c2$CA2_SLM_robustness_class == "robust_to_missingness_and_QC"))
  # class 1 must genuinely have support outside CA2-SLM
  c1 <- s[s$evidence_class == "1_FDR_supported_outside_CA2SLM", ]
  expect_true(all(c1$n_fdr05_outside_ca2slm > 0))
  # class 3 must have no FDR support at all
  c3 <- s[s$evidence_class == "3_no_FDR_support_effect_and_consistency_only", ]
  expect_true(all(c3$n_fdr05 == 0))
  # the excluded list must be non-empty and must not overlap the panel
  ex <- prd("candidate_panel_excluded_ca2slm.csv")
  expect_gt(nrow(ex), 0)
  expect_equal(length(intersect(ex$ProteinGroupID, s$ProteinGroupID)), 0L)
})

test_that("panel effects equal the canonical DA output exactly", {
  skip_if_not(have_panel, "candidate panel not built")
  e <- prd("candidate_panel_sus_res_effects.csv")
  da_dir <- repo("data", "processed", "02_id_mapping", "mapped",
                 "neuron_neuropil", "forward", "per_file")
  skip_if_not(dir.exists(da_dir), "canonical DA inputs not present")
  for (u in unique(e$spatial_unit)) {
    key <- tolower(gsub("_", "", u))
    f <- file.path(da_dir, sprintf("%ssus_%sres.csv", key, key))
    skip_if_not(file.exists(f), paste("missing canonical file for", u))
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    ee <- e[e$spatial_unit == u, , drop = FALSE]
    i <- match(ee$ProteinGroupID, d$ProteinGroupID)
    expect_equal(ee$log2FC, d$log2fc[i], tolerance = 0)
    expect_equal(ee$fdr_bh, d$padj[i], tolerance = 0)
  }
})

test_that("the panel abundance table carries no pseudoreplication", {
  skip_if_not(have_panel, "candidate panel not built")
  a <- prd("candidate_panel_animal_abundance.csv")
  key <- paste(a$ProteinGroupID, a$AnimalID, a$spatial_unit, sep = "|")
  expect_equal(anyDuplicated(key), 0L)
  expect_equal(nrow(a), 10L * 9L * 10L)
  expect_true(all(a$biological_unit == "animal"))
})

test_that("the panel vector outputs exist and are not rasterised", {
  skip_if_not(dir.exists(PFIG), "panel figure not rendered")
  for (s in c("immunostaining_panel_spatial_comparison.svg",
              "immunostaining_panel_abundance.svg",
              "immunostaining_panel_sus_res_log2fc.svg")) {
    p <- file.path(PFIG, s)
    expect_true(file.exists(p))
    expect_gt(file.size(p), 0)
    txt <- readLines(p, warn = FALSE)
    expect_true(any(grepl("</svg>", txt, fixed = TRUE)))
    expect_true(any(grepl("<text", txt, fixed = TRUE)))
    expect_false(any(grepl("<image", txt, fixed = TRUE)))
  }
})

# ------------------------------------------- the separation screen
#
# This screen exists because the FDR screen answers the wrong question for
# staining: with three animals per group a BH-FDR value is variance dominated,
# so it ranks small tight effects above large visible ones. The failure mode
# here is therefore not a wrong p-value but a flattering one - a separation
# score computed over imputed cells, which are artificially tight. These tests
# pin the gates that prevent that.

SSRC <- repo("results", "source_data", "10_biological_integration",
             "immunostaining_separation_screen")
SFIG <- repo("results", "figures", "manuscript",
             "immunostaining_separation_panel")
have_sep <- dir.exists(SSRC) &&
  file.exists(file.path(SSRC, "separation_candidate_panel.csv"))
sepd <- function(f) utils::read.csv(file.path(SSRC, f), stringsAsFactors = FALSE)

test_that("the separation panel is ranked and unambiguously targetable", {
  skip_if_not(have_sep, "separation screen not built")
  s <- sepd("separation_candidate_panel.csv")
  expect_equal(nrow(s), 10L)
  expect_equal(anyDuplicated(s$ProteinGroupID), 0L)
  expect_equal(sort(s$rank), 1:10)
  # Ranked on breadth then separation strength - never on fold change and never
  # on any p-value - with housekeeping candidates pushed last rather than
  # dropped. Reproduce that exact key from the recorded columns.
  key <- order(as.logical(s$housekeeping_flag), -s$n_units_qualifying,
               -abs(s$best_ssmd))
  expect_equal(s$rank[key], 1:10)
  expect_true(all(s$selection_contract == "immunostaining_separation_screen_v1"))
  # the three already nominated must not reappear
  expect_false(any(toupper(s$GeneSymbol) %in% c("OGA", "SLC22A23", "ANXA2")))
})

test_that("every declared gate actually held", {
  skip_if_not(have_sep, "separation screen not built")
  s <- sepd("separation_candidate_panel.csv")
  p <- sepd("separation_screen_provenance.csv")
  expect_equal(p$min_abs_log2FC, 1)
  expect_equal(p$min_abs_ssmd, 1)
  expect_true(as.logical(p$requires_complete_separation))
  expect_true(as.logical(p$requires_fully_observed))
  # The Z'-factor is reported, never gated on: nothing in this dataset reaches
  # Z' > 0, so gating on it would return an empty panel rather than an honest
  # one. The flag must stay FALSE or the reported ranking means something else.
  expect_false(as.logical(p$z_factor_used_as_gate))
  expect_true(all(abs(s$best_log2FC) >= p$min_abs_log2FC))
  expect_true(all(abs(s$best_ssmd) >= p$min_abs_ssmd))
  # and the qualifying cells behind the panel obey the same gates
  q <- sepd("separation_qualifying_cells.csv")
  expect_gt(nrow(q), 0)
  expect_true(all(abs(q$log2FC) >= p$min_abs_log2FC))
  expect_true(all(abs(q$ssmd) >= p$min_abs_ssmd))
  expect_true(all(s$ProteinGroupID %in% q$ProteinGroupID))
})

test_that("no candidate rests on imputed data", {
  skip_if_not(have_sep, "separation screen not built")
  # Imputation draws from a downshifted normal, which shrinks within-group SD
  # and inflates every variance-based separation score. A candidate whose
  # winning cell is not fully observed would be an artefact of that.
  s <- sepd("separation_candidate_panel.csv")
  e <- sepd("separation_panel_sus_res_effects.csv")
  win <- merge(s[, c("ProteinGroupID", "best_unit")], e,
               by.x = c("ProteinGroupID", "best_unit"),
               by.y = c("ProteinGroupID", "spatial_unit"))
  expect_equal(nrow(win), 10L)
  expect_true(all(as.logical(win$fully_observed)))
  expect_true(all(as.logical(win$complete_separation)))
})

test_that("the separation measures reproduce from the animal-level values", {
  skip_if_not(have_sep, "separation screen not built")
  a <- sepd("separation_panel_animal_abundance.csv")
  e <- sepd("separation_panel_sus_res_effects.csv")
  # recompute SSMD and log2FC independently from the plotted points
  g <- a[a$ExpGroup %in% c("SUS", "RES"), ]
  agg <- stats::aggregate(abundance ~ ProteinGroupID + spatial_unit + ExpGroup,
                          data = g,
                          FUN = function(z) c(m = mean(z), s = stats::sd(z)))
  st <- data.frame(agg[, 1:3], m = agg$abundance[, "m"],
                   s = agg$abundance[, "s"])
  w <- merge(st[st$ExpGroup == "SUS", ], st[st$ExpGroup == "RES", ],
             by = c("ProteinGroupID", "spatial_unit"), suffixes = c("_S", "_R"))
  w$log2FC <- w$m_S - w$m_R
  w$ssmd <- w$log2FC / sqrt(w$s_S^2 + w$s_R^2)
  chk <- merge(w, e, by = c("ProteinGroupID", "spatial_unit"))
  expect_equal(nrow(chk), 100L)
  expect_equal(chk$log2FC.x, chk$log2FC.y, tolerance = 1e-9)
  expect_equal(chk$ssmd.x, chk$ssmd.y, tolerance = 1e-9)
  # and the group-mean difference equals the canonical DA log2FC, so this panel
  # is a different ranking of the same effects, not a different effect estimate
  da_dir <- repo("data", "processed", "02_id_mapping", "mapped",
                 "neuron_neuropil", "forward", "per_file")
  skip_if_not(dir.exists(da_dir), "canonical DA inputs not present")
  for (u in unique(e$spatial_unit)) {
    key <- tolower(gsub("_", "", u))
    f <- file.path(da_dir, sprintf("%ssus_%sres.csv", key, key))
    skip_if_not(file.exists(f), paste("missing canonical file for", u))
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    ee <- e[e$spatial_unit == u, , drop = FALSE]
    expect_equal(ee$log2FC,
                 d$log2fc[match(ee$ProteinGroupID, d$ProteinGroupID)],
                 tolerance = 1e-9)
  }
})

test_that("the separation panel claims no FDR support and flags its cautions", {
  skip_if_not(have_sep, "separation screen not built")
  e <- sepd("separation_panel_sus_res_effects.csv")
  s <- sepd("separation_candidate_panel.csv")
  # by construction this panel is not FDR-ranked; nothing in it may be marked
  # FDR-supported, or the figure would read as a statistical result
  expect_false(any(as.logical(e$fdr_supported)))
  expect_false("fdr_bh" %in% names(e))
  # a housekeeping or low-abundance protein may rank, but must be flagged
  # rather than silently dropped
  expect_true("housekeeping_flag" %in% names(s))
  expect_true("low_abundance_flag" %in% names(s))
  expect_true(any(as.logical(s$housekeeping_flag) |
                    as.logical(s$low_abundance_flag)))
})

test_that("the separation abundance table carries no pseudoreplication", {
  skip_if_not(have_sep, "separation screen not built")
  a <- sepd("separation_panel_animal_abundance.csv")
  key <- paste(a$ProteinGroupID, a$AnimalID, a$spatial_unit, sep = "|")
  expect_equal(anyDuplicated(key), 0L)
  expect_equal(nrow(a), 10L * 9L * 10L)
  expect_true(all(a$biological_unit == "animal"))
})

test_that("the separation panel is registered as manuscript-supporting", {
  y <- yaml::read_yaml(repo("figures", "figure_contract.yml"))
  entry <- y$figures$S_immunostaining_separation
  expect_false(is.null(entry))
  expect_false(isTRUE(entry$is_numbered_manuscript_figure))
  expect_equal(length(entry$panels), 2L)
  for (p in entry$panels) {
    expect_equal(p$biological_unit, "animal")
    expect_true(grepl("animal_level", p$hemisphere_handling))
    expect_true(all(grepl("^results/source_data/", unlist(p$input_dependencies))))
  }
  # the contract must record that Z' did not gate, matching the provenance row
  eff <- entry$panels[[2]]
  expect_false(isTRUE(eff$z_factor_used_as_gate))
  expect_true(all(c("fully_observed_only", "complete_separation") %in%
                    unlist(eff$selection_gates)))
})

test_that("the separation vector outputs exist and are not rasterised", {
  skip_if_not(dir.exists(SFIG), "separation figure not rendered")
  for (s in c("immunostaining_separation_spatial_comparison.svg",
              "immunostaining_separation_abundance.svg",
              "immunostaining_separation_sus_res_log2fc.svg")) {
    p <- file.path(SFIG, s)
    expect_true(file.exists(p))
    expect_gt(file.size(p), 0)
    txt <- readLines(p, warn = FALSE)
    expect_true(any(grepl("</svg>", txt, fixed = TRUE)))
    expect_true(any(grepl("<text", txt, fixed = TRUE)))
    expect_false(any(grepl("<image", txt, fixed = TRUE)))
  }
})
