source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "protigy_input_utils.R"))
source(repo_path("R", "spatial_systems_data_utils.R"))
source(repo_path("R", "spatial_systems_bilateral_utils.R"))
source(repo_path("R", "spatial_systems_evidence_registry.R"))
source(repo_path("R", "ewce_gene_set_engine.R"))
source(repo_path("R", "control_spatial_identity_utils.R"))

# ---------------------------------------------------------------- fixtures

# Two animals x two spatial units x two sides, with a DELIBERATE technical
# duplicate on one side so within-side collapsing can be tested separately from
# across-side averaging. The canonical regression case is A1/U1: L = 2, R = 6.
sps_fx <- function(drop = NULL) {
  meta <- expand.grid(AnimalID = c("A1", "A2"), Region = c("U1", "U2"),
                      ReplicateGroup = c("Left", "Right"),
                      stringsAsFactors = FALSE)
  meta$Sample <- paste(meta$AnimalID, meta$Region, meta$ReplicateGroup, sep = "_")
  meta$StressGroup <- ifelse(meta$AnimalID == "A1", "CON", "SUS")
  meta$Layer <- NA_character_
  # technical duplicate: a second row for A1/U1/Left
  dup <- meta[meta$AnimalID == "A1" & meta$Region == "U1" &
                meta$ReplicateGroup == "Left", , drop = FALSE]
  dup$Sample <- paste0(dup$Sample, "_tech2")
  meta <- rbind(meta, dup)
  if (!is.null(drop)) meta <- meta[meta$Sample != drop, , drop = FALSE]

  v <- rep(0, nrow(meta))
  is_a1u1 <- meta$AnimalID == "A1" & meta$Region == "U1"
  # A1/U1 Left technical rows are 1 and 3 -> within-side mean 2
  v[is_a1u1 & meta$ReplicateGroup == "Left" & !grepl("tech2", meta$Sample)] <- 1
  v[is_a1u1 & grepl("tech2", meta$Sample)] <- 3
  v[is_a1u1 & meta$ReplicateGroup == "Right"] <- 6
  mat <- rbind(PG1 = v, PG2 = v + 10)
  colnames(mat) <- meta$Sample
  list(mat = mat, meta = meta)
}
sps_build_fx <- function(...) {
  fx <- sps_fx(...)
  sps_build_spatial_levels("neuron_soma", mat = fx$mat, metadata = fx$meta)
}

k <- function(...) paste(..., sep = "\037")

# =====================================================================
# Q1-Q4  the hemisphere contract
# =====================================================================

testthat::test_that("left and right remain distinct, and the L=2 / R=6 case holds", {
  lv <- sps_build_fx()
  # THE REGRESSION CASE: a hemisphere-resolved object must return 2 and 6,
  # never 4 and 4. An object returning 4/4 is hemisphere-AVERAGED and would
  # make every bilateral difference exactly zero.
  testthat::expect_equal(unname(lv$level3$left["PG1", k("A1", "U1")]), 2)
  testthat::expect_equal(unname(lv$level3$right["PG1", k("A1", "U1")]), 6)
  testthat::expect_false(isTRUE(all.equal(
    unname(lv$level3$left["PG1", k("A1", "U1")]),
    unname(lv$level3$right["PG1", k("A1", "U1")]))))
  # and the bilateral value is their midpoint
  testthat::expect_equal(unname(lv$level2$mat["PG1", k("A1", "U1")]), 4)
})

testthat::test_that("technical rows average only WITHIN a side first", {
  lv <- sps_build_fx()
  # A1/U1/Left has two technical rows (1 and 3). Their within-side mean is 2.
  # If technical rows were pooled across sides instead, the left value would be
  # contaminated by the right value of 6.
  testthat::expect_equal(unname(lv$level1$mat["PG1", k("A1", "L", "U1")]), 2)
  testthat::expect_equal(
    lv$level1$meta$n_source_rows[lv$level1$meta$key == k("A1", "L", "U1")], 2L)
  testthat::expect_equal(unname(lv$level1$mat["PG1", k("A1", "R", "U1")]), 6)
})

testthat::test_that("the bilateral mean uses equal side weighting", {
  lv <- sps_build_fx()
  # Left contributed 2 source rows and right 1, but each SIDE still gets half
  # the weight. A sample-weighted mean would give (1+3+6)/3 = 3.33, not 4.
  testthat::expect_equal(unname(lv$level2$mat["PG1", k("A1", "U1")]), 4)
  testthat::expect_false(isTRUE(all.equal(
    unname(lv$level2$mat["PG1", k("A1", "U1")]), 10 / 3)))
})

testthat::test_that("a one-sided cell follows the canonical no-imputation policy", {
  lv <- sps_build_fx(drop = "A1_U1_Right")
  key <- k("A1", "U1")
  testthat::expect_equal(
    lv$level2$meta$hemisphere_status[lv$level2$meta$key == key], "left_only_observed")
  testthat::expect_equal(
    lv$level2$meta$aggregation_method[lv$level2$meta$key == key],
    "single_observed_hemisphere_no_imputation")
  # the observed side passes through unchanged and nothing is invented
  testthat::expect_equal(unname(lv$level2$mat["PG1", key]), 2)
  testthat::expect_true(is.na(lv$level3$right["PG1", key]))
})

# =====================================================================
# Q5-Q6  the replicate is the animal
# =====================================================================

testthat::test_that("no side is treated as an independent biological replicate", {
  lv <- sps_build_fx()
  # the raw ReplicateGroup name never survives into the analysis levels
  testthat::expect_false("ReplicateGroup" %in% names(lv$level1$meta))
  testthat::expect_false("ReplicateGroup" %in% names(lv$level2$meta))
  testthat::expect_true("Hemisphere" %in% names(lv$level1$meta))
  testthat::expect_true(all(lv$level1$meta$Hemisphere %in% c("L", "R")))
  # provenance back to the source column is retained
  testthat::expect_true(all(lv$level1$meta$source_hemisphere_field == "ReplicateGroup"))
  # the guard rejects a frame that reintroduces the raw column
  bad <- lv$level2$meta; bad$ReplicateGroup <- "Left"
  testthat::expect_error(sps_assert_animal_is_replicate(bad), "raw ReplicateGroup")
  # level 2 has ONE column per animal x unit, not one per sample
  testthat::expect_equal(ncol(lv$level2$mat), 4L)
})

testthat::test_that("AnimalID is preserved at every level", {
  lv <- sps_build_fx()
  for (nm in c("level0", "level1", "level2")) {
    d <- if (nm == "level0") lv[[nm]] else lv[[nm]]$meta
    testthat::expect_true("AnimalID" %in% names(d), info = nm)
    testthat::expect_setequal(unique(d$AnimalID), c("A1", "A2"))
  }
  testthat::expect_error(sps_assert_animal_is_replicate(data.frame(x = 1)),
                         "no AnimalID column")
})

testthat::test_that("an unknown or missing side is rejected outright", {
  testthat::expect_error(sps_normalize_hemisphere(c("Left", "Dorsal")),
                         "Unrecognised")
  testthat::expect_error(sps_normalize_hemisphere(c("Left", NA)), "missing")
  testthat::expect_identical(as.character(sps_normalize_hemisphere(c("Left", "R"))),
                             c("L", "R"))
})

# =====================================================================
# Q7-Q8  contrast registry and side purity
# =====================================================================

testthat::test_that("the contrast registry exists exactly once", {
  # R/ is organised into domain directories, so the search must be recursive.
  # A non-recursive glob would find zero definitions and pass vacuously.
  files <- list.files(repo_path("R"), pattern = "[.]R$",
                      recursive = TRUE, full.names = TRUE)
  defs <- sum(vapply(files, function(f) {
    any(grepl("^control_spatial_contrast_registry <- function",
              readLines(f, warn = FALSE)))
  }, logical(1)))
  testthat::expect_identical(defs, 1L)
  # and the manuscript-locked count is pinned
  soma <- control_spatial_contrast_registry("neuron_soma", c("CA1", "CA2", "CA3", "DG"))
  np <- control_spatial_contrast_registry(
    "neuron_neuropil", c("CA1_SLM", "CA1_SO", "CA1_SR", "CA2_SLM", "CA2_SO",
                         "CA2_SR", "CA3_SO", "CA3_SR", "DG_MO", "DG_PO"))
  testthat::expect_identical(
    sum(control_spatial_contrast_is_manuscript_locked(soma)) +
      sum(control_spatial_contrast_is_manuscript_locked(np)),
    control_spatial_manuscript_contrast_count())
})

testthat::test_that("microglia contrasts are region-level and not manuscript-locked", {
  mg <- control_spatial_contrast_registry("microglia", c("CA1", "CA2", "CA3", "DG"))
  testthat::expect_gt(length(mg), 0L)
  testthat::expect_false(any(control_spatial_contrast_is_manuscript_locked(mg)))
  # no layer token may appear in a microglia contrast name
  testthat::expect_false(any(grepl("SLM|_SO|_SR|_MO|_PO|strata|layers",
                                   names(mg))))
})

testthat::test_that("a side-specific model receives only the requested hemisphere", {
  lv <- sps_build_fx()
  l <- sps_side_samples(lv, "L")
  r <- sps_side_samples(lv, "R")
  testthat::expect_true(all(l$Hemisphere == "L"))
  testthat::expect_true(all(r$Hemisphere == "R"))
  testthat::expect_length(intersect(l$Sample, r$Sample), 0L)
  # and a contrast may not be silently redefined when a unit is absent
  spec <- control_spatial_contrast_registry("neuron_soma", c("CA1", "CA2", "CA3", "DG"))[[1]]
  testthat::expect_error(
    control_spatial_contrast_vector(spec, c("anatomical_unit_CA1", "anatomical_unit_CA2")),
    "Refusing to silently redefine")
})

# =====================================================================
# Q10-Q11  marker transfer direction purity
# =====================================================================

testthat::test_that("marker transfer never reuses the evaluation side for discovery", {
  set.seed(4)
  ids <- paste0("P", 1:200)
  eff_L <- stats::rnorm(200); eff_R <- eff_L + stats::rnorm(200, sd = 0.1)
  # discovery on L, evaluation on R
  disc_L <- ids[utils::head(order(-eff_L), 25)]
  tr <- sps_rank_transfer(disc_L, eff_R, ids)
  testthat::expect_identical(tr$n_transferred, 25L)
  testthat::expect_true(tr$auc_like > 0.5)
  # the discovered set must be a function of L ONLY: permuting R cannot change it
  disc_L2 <- ids[utils::head(order(-eff_L), 25)]
  testthat::expect_identical(disc_L, disc_L2)
  # and the reverse direction picks a genuinely different discovery set
  disc_R <- ids[utils::head(order(-eff_R), 25)]
  testthat::expect_false(identical(disc_L, disc_R))
  # transferring an unrelated set gives no advantage
  rnd <- sps_rank_transfer(sample(ids, 25), eff_R, ids)
  testthat::expect_lt(abs(rnd$auc_like - 0.5), 0.35)
})

# =====================================================================
# Q13-Q14  pairing and variance-model identity
# =====================================================================

testthat::test_that("paired sides match on the exact same animal and unit", {
  df <- data.frame(
    AnimalID = c("A1", "A1", "A2", "A2"), SpatialUnit = c("U1", "U1", "U2", "U2"),
    Hemisphere = c("L", "R", "L", "R"), value = c(1, 2, 10, 20),
    stringsAsFactors = FALSE)
  w <- sps_pair_sides(df, c("AnimalID", "SpatialUnit"))
  testthat::expect_equal(nrow(w), 2L)
  testthat::expect_equal(w$left, c(1, 10))
  testthat::expect_equal(w$right, c(2, 20))
  # a left value may never be paired with a different animal's right value
  testthat::expect_true(all(w$AnimalID %in% c("A1", "A2")))
  # duplicate id x side is rejected rather than silently averaged
  dup <- rbind(df, df[1, ])
  testthat::expect_error(sps_pair_sides(dup, c("AnimalID", "SpatialUnit")),
                         "More than one value per id")
  # un-normalised sides are rejected
  raw <- df; raw$Hemisphere <- c("Left", "Right", "Left", "Right")
  testthat::expect_error(sps_pair_sides(raw, c("AnimalID", "SpatialUnit")),
                         "normalized to L/R")
})

testthat::test_that("the variance model formula separates hemisphere from animal", {
  source(repo_path("R", "spatial_systems_endpoint_utils.R"))
  # sized like the real design (9 animals x 4 units x 2 sides) so the variance
  # components are estimable rather than singular
  e <- expand.grid(AnimalID = paste0("A", 1:9),
                   SpatialUnit = c("U1", "U2", "U3", "U4"),
                   Hemisphere = c("L", "R"), stringsAsFactors = FALSE)
  e$dataset <- "d"; e$endpoint_class <- "test"; e$endpoint_id <- "e1"
  e$StressGroup <- "CON"
  set.seed(7)
  # all four components genuinely non-zero, or the fit is singular by
  # construction and the reliability guard (correctly) suppresses the estimate
  animal_effect <- stats::setNames(stats::rnorm(9, sd = 1.5), paste0("A", 1:9))
  unit_effect <- stats::setNames(stats::rnorm(4, sd = 0.8), c("U1", "U2", "U3", "U4"))
  cell_effect <- stats::setNames(stats::rnorm(36, sd = 0.4),
                                 unique(paste(e$AnimalID, e$SpatialUnit)))
  e$value <- animal_effect[e$AnimalID] + unit_effect[e$SpatialUnit] +
    cell_effect[paste(e$AnimalID, e$SpatialUnit)] +
    stats::rnorm(nrow(e), sd = 0.3)
  vc <- sps_variance_components(e)
  testthat::expect_true(grepl("Hemisphere", vc$formula_used))
  testthat::expect_true(grepl("[(]1[ ]*[|][ ]*AnimalID[)]", vc$formula_used))
  testthat::expect_true(grepl("[(]1[ ]*[|][ ]*SpatialUnit[)]", vc$formula_used))
  # Hemisphere is FIXED and AnimalID is RANDOM - they must not be interchanged
  testthat::expect_false(grepl("[(]1[ ]*[|][ ]*Hemisphere[)]", vc$formula_used))
  testthat::expect_true(all(c("between_animal_variance",
                              "within_animal_hemispheric_variance",
                              "is_singular", "model_convergence",
                              "assumption_status") %in% names(vc)))
  # variance driven by animal identity must land in the animal component
  testthat::skip_if_not(vc$assumption_status == "assumptions_met")
  testthat::expect_gt(vc$between_animal_variance, vc$within_animal_hemispheric_variance)
})

testthat::test_that("reliability is reported only when the guards pass", {
  source(repo_path("R", "spatial_systems_endpoint_utils.R"))
  vc <- data.frame(
    dataset = "d", endpoint_class = "c", endpoint_id = c("ok", "singular"),
    n_animals = 9L, n_spatial_units = 4L,
    between_animal_variance = c(1, 1),
    within_animal_hemispheric_variance = c(1, 1),
    hemisphere_fixed_effect = 0, hemisphere_fixed_effect_p = 1,
    model_convergence = "converged", is_singular = c(FALSE, TRUE),
    assumption_status = c("assumptions_met", "singular_fit_reliability_not_reported"),
    formula_used = "f", stringsAsFactors = FALSE)
  pg <- sps_precision_gain(vc)
  testthat::expect_equal(pg$ICC_single_side[1], 0.5)
  testthat::expect_equal(pg$ICC_bilateral_mean[1], 2 / 3)
  testthat::expect_true(is.na(pg$ICC_single_side[2]))
  testthat::expect_true(is.na(pg$reliability_formula[2]))
  testthat::expect_equal(pg$relative_measurement_variance_reduction[1], 0.5)
})

# =====================================================================
# Q15-Q18  the EWCE annotation layer
# =====================================================================

testthat::test_that("the EWCE annotation API has no phenotype dependency", {
  a <- names(formals(run_ewce_gene_set_annotation))
  for (bad in c("StressGroup", "group", "contrast", "condition", "direction",
                "DAP", "phenotype")) {
    testthat::expect_false(tolower(bad) %in% tolower(a), info = bad)
  }
  # and a gene-set frame carrying a phenotype column is rejected
  gs <- data.frame(gene_set_id = "m1", gene_symbol = "Gria1",
                   StressGroup = "SUS", stringsAsFactors = FALSE)
  testthat::expect_error(
    run_ewce_gene_set_annotation(gs, background = "Gria1", reference = list(),
                                 dataset = "d"),
    "phenotype-blind by contract")
})

testthat::test_that("module-annotation FDR is unaffected by Differential-arm rows", {
  set.seed(99)
  fam <- ewce_fdr_family_module_annotation("microglia", "all", 1L)
  mod <- data.frame(gene_set_id = rep(paste0("m", 1:5), each = 4),
                    cell_type = rep(paste0("CT", 1:4), 5),
                    p_value = stats::runif(20), fdr_family = fam,
                    stringsAsFactors = FALSE)
  a <- ewce_apply_family_fdr(mod)
  extra <- data.frame(gene_set_id = paste0("d", 1:40),
                      cell_type = rep(paste0("CT", 1:4), 10),
                      p_value = stats::runif(40) / 1000,
                      fdr_family = ewce_fdr_family_differential("microglia", 1L),
                      stringsAsFactors = FALSE)
  b <- ewce_apply_family_fdr(rbind(mod, extra))
  bmod <- b[b$fdr_family == fam, ]
  # BIT-IDENTICAL, not merely close
  testthat::expect_identical(a$p_value, bmod$p_value)
  testthat::expect_identical(a$FDR, bmod$FDR)
  # the families are genuinely different strings
  testthat::expect_false(identical(fam, ewce_fdr_family_differential("microglia", 1L)))
  # a row without a family is refused rather than pooled globally
  nofam <- mod; nofam$fdr_family <- NA_character_
  testthat::expect_error(ewce_apply_family_fdr(nofam), "explicit FDR family")
})

testthat::test_that("module scopes are the canonical three", {
  testthat::expect_identical(ewce_module_scopes(), c("all", "core_kME06", "top25"))
  testthat::expect_error(
    ewce_fdr_family_module_annotation("d", "not_a_scope", 1L))
  # each scope yields its own family
  fams <- vapply(ewce_module_scopes(),
                 function(s) ewce_fdr_family_module_annotation("d", s, 1L),
                 character(1))
  testthat::expect_identical(anyDuplicated(fams), 0L)
})

testthat::test_that("the EWCE engine refuses an empty background", {
  testthat::expect_error(
    ewce_bootstrap_once(hits = letters, background = character(),
                        reference = list(), annot_level = 1L),
    "background is empty")
  testthat::expect_error(
    run_ewce_gene_set_annotation(list(m1 = letters), background = character(),
                                 reference = list(), dataset = "d"),
    "measured background is required")
})

# =====================================================================
# Q  evidence registry
# =====================================================================

testthat::test_that("the evidence registry classifies every stream", {
  reg <- sps_evidence_dependence_registry()
  testthat::expect_gte(nrow(reg), 11L)
  testthat::expect_identical(anyDuplicated(reg$evidence_id), 0L)
  testthat::expect_true(all(reg$independence_class %in% sps_independence_classes()))
  testthat::expect_true(all(nzchar(reg$allowed_interpretation)))
  testthat::expect_true(all(nzchar(reg$prohibited_interpretation)))
  # bilateral streams must be classed as internal reproducibility, never as
  # independent replication
  bil <- reg[grepl("bilateral|cross_hemisphere", reg$evidence_id), , drop = FALSE]
  testthat::expect_gt(nrow(bil), 0L)
  testthat::expect_true(all(bil$independence_class == "internal_reproducibility"))
  testthat::expect_true(all(grepl("not independent", ignore.case = TRUE,
                                  bil$prohibited_interpretation)))
  # the phenotype-blind EWCE stream must not be marked primary phenotype evidence
  mod <- reg[reg$evidence_id == "ewce_module_annotation", ]
  testthat::expect_identical(mod$phenotype_used, "none")
  testthat::expect_identical(mod$independence_class, "external_annotation")
})

# =====================================================================
# generated outputs
# =====================================================================

sps_out <- function(...) path_results("tables", "11_spatial_systems", ...)

testthat::test_that("the foundation validation has no critical failure", {
  p <- sps_out("spatial_systems_foundation_validation.csv")
  testthat::skip_if_not(file.exists(p), "foundation validation not generated")
  v <- utils::read.csv(p, stringsAsFactors = FALSE)
  testthat::expect_gt(nrow(v), 10L)
  testthat::expect_equal(sum(v$critical %in% TRUE & v$status == "FAIL"), 0L)
})

testthat::test_that("generated bilateral tables keep the sides apart", {
  p <- sps_out("bilateral", "bilateral_spatial_identity_protein_level.csv")
  testthat::skip_if_not(file.exists(p), "bilateral spatial identity not generated")
  x <- utils::read.csv(p, stringsAsFactors = FALSE)
  testthat::expect_true(all(c("estimate_L", "estimate_R", "estimate_bilateral",
                              "sign_agreement", "abs_L_minus_R") %in% names(x)))
  # if L and R were the same object every difference would be exactly zero
  testthat::expect_gt(stats::median(x$abs_L_minus_R, na.rm = TRUE), 0)
})

testthat::test_that("module bilateral output separates level from pattern", {
  a <- sps_out("bilateral", "WGCNA_module_bilateral_reproducibility.csv")
  b <- sps_out("bilateral", "WGCNA_module_bilateral_spatial_profile_reproducibility.csv")
  testthat::skip_if_not(file.exists(a) && file.exists(b), "module bilateral not generated")
  x <- utils::read.csv(a, stringsAsFactors = FALSE)
  y <- utils::read.csv(b, stringsAsFactors = FALSE)
  testthat::expect_true(all(c("pearson_r", "MAE", "mean_signed_L_minus_R",
                              "bilateral_reproducibility_class") %in% names(x)))
  testthat::expect_true(all(c("spatial_profile_Pearson", "AnimalID",
                              "n_spatial_units") %in% names(y)))
  # the per-animal profile table must be keyed by animal, not pooled
  testthat::expect_gt(length(unique(y$AnimalID)), 1L)
})

testthat::test_that("the registry reproduces the original inline contrast weights", {
  # The manuscript-locked contrasts were previously built inline inside
  # control_spatial_identity_main(). Moving them into a registry must be a pure
  # refactor: same names, same weights, same count. This reconstructs the
  # original construction and compares it element-wise.
  orig <- function(dataset, unit_levels, design_cols) {
    mk <- function(name, weights) {
      v <- stats::setNames(rep(0, length(design_cols)), design_cols)
      v[paste0("anatomical_unit_", names(weights))] <- weights
      list(name = name, weights = v)
    }
    out <- list()
    if (dataset == "neuron_soma") {
      for (target in sort(unique(unit_levels))) {
        out[[length(out) + 1L]] <- mk(paste0(target, "_vs_mean_other_soma_regions"),
          control_spatial_target_rest_weights(unit_levels, target))
      }
    }
    if (dataset == "neuron_neuropil") {
      regions <- sub("_.*$", "", unit_levels)
      out[[1]] <- mk("DG_neuropil_vs_mean_non_DG_regions",
        control_spatial_region_mean_weights(unit_levels, regions, "DG"))
      if (all(c("CA1_SO", "CA3_SO") %in% unit_levels)) {
        out[[length(out) + 1L]] <- mk("CA1_SO_vs_CA3_SO",
          stats::setNames(c(1, -1), c("CA1_SO", "CA3_SO")))
      }
      ca1 <- unit_levels[grepl("^CA1_", unit_levels)]
      if (length(ca1) >= 3L) for (target in ca1) {
        out[[length(out) + 1L]] <- mk(paste0(target, "_vs_mean_other_CA1_strata"),
          control_spatial_target_rest_weights(ca1, target))
      }
      dg <- unit_levels[grepl("^DG_", unit_levels)]
      if (length(dg) >= 2L) for (target in dg) {
        out[[length(out) + 1L]] <- mk(paste0(target, "_vs_mean_other_DG_layers"),
          control_spatial_target_rest_weights(dg, target))
      }
    }
    stats::setNames(out, vapply(out, function(z) z$name, character(1)))
  }

  cases <- list(
    neuron_soma = c("CA1", "CA2", "CA3", "DG"),
    neuron_neuropil = c("CA1_SLM", "CA1_SO", "CA1_SR", "CA2_SLM", "CA2_SO",
                        "CA2_SR", "CA3_SO", "CA3_SR", "DG_MO", "DG_PO"))
  for (ds in names(cases)) {
    ul <- cases[[ds]]
    dc <- c(paste0("anatomical_unit_", ul), "hemisphere_R")
    o <- orig(ds, ul, dc)
    n <- lapply(control_spatial_contrast_registry(ds, ul),
                function(s) control_spatial_contrast_vector(s, dc))
    testthat::expect_setequal(names(o), names(n))
    for (k in names(o)) {
      testthat::expect_equal(n[[k]], o[[k]]$weights, info = paste(ds, k))
    }
  }
})

testthat::test_that("the canonical Stage-09 script consumes the registry", {
  src <- paste(readLines(repo_path("analysis/04_differential_abundance",
                                   "09_control_spatial_identity_validation.r"),
                         warn = FALSE), collapse = "\n")
  testthat::expect_true(grepl("control_spatial_contrast_registry", src, fixed = TRUE))
  testthat::expect_true(grepl("control_spatial_contrast_vector", src, fixed = TRUE))
  # the inline builder must be gone, or there would be two copies of a
  # manuscript-locked list that can drift apart
  testthat::expect_false(grepl("make_contrast <- function", src, fixed = TRUE))
})
