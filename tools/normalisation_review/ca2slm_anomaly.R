# What makes CA2-SLM SUS vs RES jump from 28 to 464 FDR-supported proteins when A111 is removed? Read-only follow-up
# of ca2slm_original_sensitivity.R (normalisation review 2026-10-22). Both centrings throughout; the original-centred
# animal-level input is rebuilt as current + per-sample shift (asserted equal to the archived input).
#
#   S  sample quality of the 18 CA2-SLM hemispheres against all 180 neuropil samples
#   P  animal profiles: how each SUS/RES animal's CA2-SLM profile relates to the others; hemisphere concordance
#   V  residual variance per unit from the joint fit: is CA2-SLM noisier than the pooled variance the model assumes?
#   M  model variants (pipeline, ordinary t, CA2-SLM-only, unit variance weights, plate covariate, hemisphere level
#      with animal-unit blocking, observed values only) with and without A111
#   R  every relabelling of the six SUS/RES animals in CA2-SLM (3 vs 3, and 2 vs 3 without A111)
#   L  leave-one-animal-out (one unit at a time) for all ten neuropil units
#   O  offsets and 70% filter recomputed without A111; imputation redrawn with 20 seeds
#   T  representative proteins: effect size, raw and adjusted P across variants, per-animal values
# Writes results/reviewer_audit/normalisation_choice_2026-10-21/ca2slm/. Changes nothing else.
# USAGE  Rscript tools/normalisation_review/ca2slm_anomaly.R      (working directory: repo root)

suppressPackageStartupMessages({ library(readr); library(readxl); library(dplyr); library(tidyr); library(limma); library(matrixStats) })
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
REV <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21")
OUT <- file.path(REV, "ca2slm"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(x, name) write_csv(x, file.path(OUT, name), na = "")
ds <- "neuron_neuropil"; U <- "CA2_slm"; CC <- "CA2_slm_3-CA2_slm_2"
SUSA <- c("A111", "A755", "A764"); RESA <- c("A0003", "A135", "A139")
o <- readRDS(file.path(REV, "validation", "validation_state.rds"))$offsets |> filter(dataset == ds) |>
  mutate(missing_fraction = unname(missing_fraction))
read_gct <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dims <- as.integer(strsplit(lines[[2]], "\t", fixed = TRUE)[[1]][1:4])
  hdr <- strsplit(lines[[3]], "\t", fixed = TRUE)[[1]]
  desc <- do.call(rbind, strsplit(lines[4:(3 + dims[4])], "\t", fixed = TRUE))
  body <- do.call(rbind, strsplit(lines[(4 + dims[4]):(3 + dims[4] + dims[1])], "\t", fixed = TRUE))
  sc <- (dims[3] + 2L):length(hdr)
  list(m = matrix(as.numeric(body[, sc]), nrow = nrow(body), dimnames = list(body[, 1], hdr[sc])),
       coldesc = setNames(lapply(seq_len(nrow(desc)), function(i) desc[i, sc]), desc[, 1]))
}
imp <- read_excel(Sys.glob(sprintf("data/processed/01_preprocessing/impute/*_pgmatrix_imputed_%s_*_missing70pct.xlsx", ds)), progress = FALSE)
mn <- as.matrix(imp[, o$sample_id]); rownames(mn) <- imp$`T: Protein.Names`
rawm <- read_tsv("data/raw/pg_matrix/quicksearch.pg_matrix.tsv", show_col_types = FALSE, progress = FALSE)
X_all <- log2(as.matrix(rawm[, -(1:4)])); rownames(X_all) <- rawm$Protein.Group
obs <- !is.na(X_all[imp$Protein.Group, o$sample_id]); rownames(obs) <- rownames(mn)
g <- read_gct(sprintf("data/processed/01_preprocessing/protigy_input_animal_level/%s/%s_animal_level.gct", ds, ds))
g_old <- read_gct(file.path(ARCHIVE, "data/processed/01_preprocessing/protigy_input_animal_level", ds, paste0(ds, "_animal_level.gct")))
col_animal <- g$coldesc$AnimalID; col_unit <- paste(g$coldesc$region, g$coldesc$layer, sep = "_")
col_key <- paste(col_animal, g$coldesc$region, g$coldesc$layer, sep = "|")
smp_key <- paste(o$animal, sub("_.*$", "", o$unit), sub("^[^_]+_", "", o$unit), sep = "|")
cell0 <- g$coldesc$region_layer_ExpGroup
plate_of <- setNames(o$plate, o$animal)[!duplicated(o$animal)]
shift_orig <- o$raw_new_complete_case - o$raw_original_all_observed

# animal-level matrix from a sample-level matrix sm (rows = proteins), keeping samples flagged in keep
to_animal <- function(sm, keep = rep(TRUE, ncol(sm))) {
  m0 <- matrix(NA_real_, nrow(g$m), ncol(g$m), dimnames = dimnames(g$m))
  for (j in seq_along(col_key)) {
    s <- o$sample_id[smp_key == col_key[j] & keep]
    if (!length(s)) next
    v <- rowMeans(sm[, s, drop = FALSE], na.rm = TRUE); v[is.nan(v)] <- NA; v <- v[!duplicated(names(v))]
    m0[, j] <- v[rownames(m0)]
  }
  m0
}
sample_matrix <- function(version, observed_only = FALSE, base = mn) {
  sm <- base; if (observed_only) sm[!obs] <- NA
  if (version == "original") sm <- sweep(sm, 2, shift_orig, `+`)
  sm
}
SM <- list(original = sample_matrix("original"), new = sample_matrix("new"))
AM <- lapply(SM, to_animal)
rid <- intersect(rownames(g$m), rownames(g_old$m))
stopifnot(max(abs(AM$original[rid, ] - g_old$m[rid, colnames(g$m)]), na.rm = TRUE) < 1e-8,
          max(abs(AM$new[rid, ] - g$m[rid, ]), na.rm = TRUE) < 1e-8)
drop_lab <- function(animals, unit = U) o$label %in% paste0(rep(animals, each = 2), c("_L_", "_R_"), unit)
is_hit <- function(t) sum(t$adj.P.Val < 0.05, na.rm = TRUE)

# ---------------------------------------------------------------- fitting
# m: animal-level matrix; cell: cell label per column; weights optional; model variants below
fit_contrast <- function(m, cell = cell0, contrast = CC, model = "pipeline", extra = NULL, weights = NULL) {
  keepc <- colSums(!is.na(m)) > 0
  m <- m[, keepc, drop = FALSE]; cell <- cell[keepc]; if (!is.null(weights)) weights <- weights[, keepc, drop = FALSE]
  unit_of <- sub("_[0-9]+$", "", cell)
  if (model %in% c("ca2_only", "ca2_only_ordinary_t")) { sel <- unit_of == U; m <- m[, sel, drop = FALSE]; cell <- cell[sel]; unit_of <- unit_of[sel] }
  f <- factor(cell, levels = unique(cell)); design <- model.matrix(~0 + f); colnames(design) <- levels(f)
  if (model == "plate") { pl <- plate_of[col_animal[keepc]]; design <- cbind(design, plateC = as.numeric(pl == "C")) }
  cm <- makeContrasts(contrasts = contrast, levels = design)
  if (model == "unit_weights") weights <- arrayWeights(m, design, var.design = model.matrix(~ factor(unit_of)))
  fit <- lmFit(m, design, weights = weights)
  fc <- contrasts.fit(fit, cm)
  if (model %in% c("ordinary_t", "ca2_only_ordinary_t")) {
    t <- fc$coefficients[, 1] / (fc$stdev.unscaled[, 1] * fc$sigma); p <- 2 * pt(-abs(t), fc$df.residual)
    return(tibble(protein = rownames(m), logFC = fc$coefficients[, 1], SE = fc$stdev.unscaled[, 1] * fc$sigma, t = t, P.Value = p,
                  adj.P.Val = p.adjust(p, "BH"), df_total = fc$df.residual, df_prior = 0, s2_prior = NA_real_, sigma = fc$sigma))
  }
  eb <- eBayes(fc, robust = TRUE)
  tibble(protein = rownames(m), logFC = eb$coefficients[, 1], SE = eb$stdev.unscaled[, 1] * sqrt(eb$s2.post), t = eb$t[, 1], P.Value = eb$p.value[, 1],
         adj.P.Val = p.adjust(eb$p.value[, 1], "BH"), df_total = eb$df.total, df_prior = eb$df.prior, s2_prior = eb$s2.prior, sigma = fc$sigma)
}
summ_fit <- function(r, ...) tibble(..., n_tested = sum(!is.na(r$adj.P.Val)), n_fdr05 = is_hit(r), n_fdr05_down = sum(r$adj.P.Val < 0.05 & r$logFC < 0, na.rm = TRUE),
  n_p001 = sum(r$P.Value < 0.001, na.rm = TRUE), pi0 = tryCatch(propTrueNull(r$P.Value[!is.na(r$P.Value)]), error = function(e) NA_real_),
  median_logFC = median(r$logFC, na.rm = TRUE), median_abs_logFC_hits = median(abs(r$logFC[r$adj.P.Val < 0.05]), na.rm = TRUE),
  median_SE = median(r$SE, na.rm = TRUE), median_df_total = median(r$df_total, na.rm = TRUE), df_prior = median(r$df_prior), s2_prior = r$s2_prior[1],
  block_correlation = if ("block_cor" %in% names(r)) r$block_cor[1] else NA_real_)

# ---------------------------------------------------------------- S. sample quality
qv <- c("missing_fraction", "Proteins.Identified", "Precursors.Identified", "MS1.Signal", "MS2.Signal", "Total.Quantity", "Normalisation.Instability",
        "FWHM.RT", "Median.RT.Prediction.Acc", "raw_original_all_observed", "raw_new_complete_case")
sq <- o |> select(label, animal, group, hemisphere, unit, plate, run_order, acquisition_date, all_of(qv)) |>
  mutate(across(c(MS1.Signal, MS2.Signal, Total.Quantity), log10), offset_difference = raw_original_all_observed - raw_new_complete_case)
z_all <- sq |> mutate(across(c(all_of(qv), offset_difference), ~ (.x - mean(.x)) / sd(.x), .names = "z_{.col}"),
                      rank_missing_in_neuropil = rank(-missing_fraction))
ca2_q <- z_all |> filter(unit == U) |> mutate(rank_missing_in_CA2slm = rank(-missing_fraction)) |> arrange(group, animal, hemisphere)
w(ca2_q, "sample_quality_ca2slm.csv")

# ---------------------------------------------------------------- P. animal profiles (animal-level CA2-SLM vectors)
prof <- list()
for (v in names(AM)) {
  A <- AM[[v]][, col_unit == U]; colnames(A) <- col_animal[col_unit == U]
  conm <- rowMeans(A[, c("A127", "A129", "A765")])
  for (a in c(SUSA, RESA)) {
    others_sus <- setdiff(SUSA, a); others_res <- setdiff(RESA, a)
    loo_prof <- rowMeans(A[, others_sus, drop = FALSE]) - rowMeans(A[, others_res, drop = FALSE])   # SUS-RES without a
    dev <- A[, a] - conm
    prof[[length(prof) + 1L]] <- tibble(centring = v, animal = a, group = if (a %in% SUSA) "SUS" else "RES",
      median_deviation_from_CON = median(dev), cor_deviation_with_loo_SUS_minus_RES = cor(dev, loo_prof),
      cor_with_mean_other_SUS = cor(A[, a] - rowMeans(A), rowMeans(A[, setdiff(SUSA, a), drop = FALSE]) - rowMeans(A)),
      cor_with_mean_RES = cor(A[, a] - rowMeans(A), rowMeans(A[, setdiff(RESA, a), drop = FALSE]) - rowMeans(A)))
  }
}
prof <- bind_rows(prof)
hem <- o |> filter(unit == U) |> group_by(animal, group) |>
  summarise(LR_cor_observed = { l <- sample_id[hemisphere == "L"]; r <- sample_id[hemisphere == "R"]; ok <- obs[, l] & obs[, r]; cor(mn[ok, l], mn[ok, r]) },
            LR_share_both_observed = { l <- sample_id[hemisphere == "L"]; r <- sample_id[hemisphere == "R"]; mean(obs[, l] & obs[, r]) }, .groups = "drop")
w(prof |> left_join(hem, by = c("animal", "group")), "animal_profiles_ca2slm.csv")

# ---------------------------------------------------------------- V. residual variance per unit (joint fit)
var_tab <- list()
for (v in names(AM)) for (variant in c("all animals", "without A111")) {
  m <- AM[[v]]; if (variant == "without A111") m[, col_unit == U & col_animal == "A111"] <- NA
  keepc <- colSums(!is.na(m)) > 0; mm <- m[, keepc]; cl <- cell0[keepc]; un <- sub("_[0-9]+$", "", cl)
  f <- factor(cl, levels = unique(cl)); design <- model.matrix(~0 + f)
  fit <- lmFit(mm, design); res <- mm - fit$coefficients %*% t(design)
  for (u in unique(un)) {
    dfu <- sum(un == u) - length(unique(cl[un == u]))
    s2u <- rowSums(res[, un == u, drop = FALSE]^2) / dfu
    var_tab[[length(var_tab) + 1L]] <- tibble(centring = v, variant = variant, unit = u, df_unit = dfu, df_pooled = fit$df.residual[1],
      median_unit_sd = median(sqrt(s2u)), median_pooled_sd = median(fit$sigma), median_ratio_unit_to_pooled_var = median(s2u / fit$sigma^2),
      share_proteins_unit_var_gt_2x_pooled = mean(s2u > 2 * fit$sigma^2))
  }
}
var_tab <- bind_rows(var_tab); w(var_tab, "residual_variance_by_unit.csv")
aw <- lapply(AM, function(m) { un <- sub("_[0-9]+$", "", cell0); f <- factor(cell0, levels = unique(cell0)); d <- model.matrix(~0 + f)
  a <- arrayWeights(m, d, var.design = model.matrix(~ factor(un))); tibble(unit = un, weight = a) |> distinct() })
w(bind_rows(aw, .id = "centring"), "unit_variance_weights.csv")

# ---------------------------------------------------------------- M. model variants
blocks_hemi <- function(sm, keep) {   # hemisphere-level joint model with the two hemispheres of an animal-unit as a block
  sm <- sm[, keep]; oo <- o[keep, ]
  cellh <- paste0(oo$unit, "_", c(CON = 1, RES = 2, SUS = 3)[oo$group])
  f <- factor(cellh, levels = unique(cellh)); design <- model.matrix(~0 + f); colnames(design) <- levels(f)
  blk <- paste(oo$animal, oo$unit); dc <- duplicateCorrelation(sm, design, block = blk)
  fit <- lmFit(sm, design, block = blk, correlation = dc$consensus.correlation)
  eb <- eBayes(contrasts.fit(fit, makeContrasts(contrasts = CC, levels = design)), robust = TRUE)
  tibble(protein = rownames(sm), logFC = eb$coefficients[, 1], SE = eb$stdev.unscaled[, 1] * sqrt(eb$s2.post), t = eb$t[, 1], P.Value = eb$p.value[, 1],
         adj.P.Val = p.adjust(eb$p.value[, 1], "BH"), df_total = eb$df.total, df_prior = eb$df.prior, s2_prior = eb$s2.prior, sigma = fit$sigma,
         block_cor = dc$consensus.correlation)
}
runs <- list()
for (v in names(AM)) for (ex in c("all animals", "without A111", "without A755 L", "without A764 R", "without A755 L and A764 R")) {
  keep <- switch(ex, "all animals" = rep(TRUE, nrow(o)), "without A111" = !drop_lab("A111"), "without A755 L" = o$label != "A755_L_CA2_slm",
                 "without A764 R" = o$label != "A764_R_CA2_slm", "without A755 L and A764 R" = !(o$label %in% c("A755_L_CA2_slm", "A764_R_CA2_slm")))
  m <- to_animal(SM[[v]], keep)
  for (md in c("pipeline", "ordinary_t", "ca2_only", "ca2_only_ordinary_t", "unit_weights", "plate"))
    runs[[paste(v, ex, md, sep = " | ")]] <- fit_contrast(m, model = md)
  runs[[paste(v, ex, "observed_only", sep = " | ")]] <- fit_contrast(to_animal(sample_matrix(v, TRUE), keep))
  runs[[paste(v, ex, "hemisphere_blocked", sep = " | ")]] <- blocks_hemi(SM[[v]], keep)
}
msum <- bind_rows(lapply(names(runs), function(k) { p <- strsplit(k, " | ", fixed = TRUE)[[1]]; summ_fit(runs[[k]], centring = p[1], exclusion = p[2], model = p[3]) }))
w(msum, "model_variants.csv")

# ---------------------------------------------------------------- R. relabelling the six SUS/RES animals in CA2-SLM
relab <- list()
ca2col <- which(col_unit == U)
for (v in names(AM)) for (ex in c("all animals", "without A111")) {
  pool <- if (ex == "all animals") c(SUSA, RESA) else setdiff(c(SUSA, RESA), "A111")
  k <- if (ex == "all animals") 3 else 2
  m <- AM[[v]]; if (ex == "without A111") m[, ca2col[col_animal[ca2col] == "A111"]] <- NA
  for (s in combn(pool, k, simplify = FALSE)) {
    cell <- cell0; cell[ca2col] <- ifelse(col_animal[ca2col] %in% s, "CA2_slm_3", ifelse(col_animal[ca2col] %in% pool, "CA2_slm_2", cell0[ca2col]))
    r <- fit_contrast(m, cell = cell)
    relab[[length(relab) + 1L]] <- tibble(centring = v, exclusion = ex, sus_like = paste(s, collapse = "+"), is_true_labelling = setequal(s, intersect(SUSA, pool)),
      contains_A755_and_A764 = all(c("A755", "A764") %in% s), contains_A755_or_A764_only_one = sum(c("A755", "A764") %in% s) == 1,
      n_fdr05 = is_hit(r), n_fdr05_down = sum(r$adj.P.Val < 0.05 & r$logFC < 0), median_logFC = median(r$logFC))
  }
}
relab <- bind_rows(relab); w(relab, "relabelling.csv")

# ---------------------------------------------------------------- L. leave-one-animal-out, one unit at a time, all neuropil units
loao <- list()
for (v in names(AM)) for (u in unique(col_unit)) {
  cc <- paste0(u, "_3-", u, "_2")
  base <- fit_contrast(AM[[v]], contrast = cc)
  loao[[length(loao) + 1L]] <- tibble(centring = v, unit = u, removed = "none", n_fdr05 = is_hit(base))
  for (a in c(SUSA, RESA)) {
    m <- AM[[v]]; m[, col_unit == u & col_animal == a] <- NA
    loao[[length(loao) + 1L]] <- tibble(centring = v, unit = u, removed = a, n_fdr05 = is_hit(fit_contrast(m, contrast = cc)))
  }
}
loao <- bind_rows(loao) |> group_by(centring, unit) |> mutate(base = n_fdr05[removed == "none"]) |> ungroup()
w(loao, "unit_leave_one_animal_out.csv")

# ---------------------------------------------------------------- O. offsets, filter and imputation
anc_all <- rownames(X_all)[rowSums(is.na(X_all[, o$sample_id])) == 0]
chk_offset <- max(abs(colMedians(X_all[anc_all, o$sample_id]) - o$raw_new_complete_case))
off_rows <- list(); filt_rows <- list()
for (ex in c("without A111 CA2-SLM", "without all A111 samples")) {
  keep <- if (ex == "without A111 CA2-SLM") !drop_lab("A111") else o$animal != "A111"
  anc <- rownames(X_all)[rowSums(is.na(X_all[, o$sample_id[keep]])) == 0]
  off_new <- colMedians(X_all[anc, o$sample_id[keep]])
  delta <- off_new - o$raw_new_complete_case[keep]                 # change of each kept sample's offset
  sm <- mn[, keep]; sm <- sweep(sm, 2, delta, `-`)
  full <- matrix(NA_real_, nrow(mn), ncol(mn), dimnames = dimnames(mn)); full[, keep] <- sm
  r <- fit_contrast(to_animal(full, keep))
  off_rows[[ex]] <- tibble(exclusion = ex, n_anchors_before = length(anc_all), n_anchors_after = length(anc), max_abs_offset_change = max(abs(delta)),
                           sd_offset_change = sd(delta), new_centring_fdr05_recomputed = is_hit(r), check_offset_reproduction = chk_offset)
  miss <- rowMeans(is.na(X_all[, o$sample_id[keep]])); miss0 <- rowMeans(is.na(X_all[, o$sample_id]))
  present <- rowSums(!is.na(X_all[, o$sample_id])) > 0
  filt_rows[[ex]] <- tibble(exclusion = ex, n_pass_before = sum(miss0 <= 0.7 & present), n_pass_after = sum(miss <= 0.7 & present),
                            entering = sum(miss <= 0.7 & miss0 > 0.7 & present), leaving = sum(miss > 0.7 & miss0 <= 0.7 & present))
}
w(bind_rows(off_rows), "offsets_recomputed.csv"); w(bind_rows(filt_rows), "filter_recomputed.csv")

impute_normal <- function(x, seed, width = 0.3, downshift = 1.8) {   # as archive/01_preprocessing/01_impute.r, per column
  set.seed(seed)
  for (j in seq_len(ncol(x))) { mi <- which(is.na(x[, j])); if (!length(mi)) next
    ob <- x[-mi, j]; x[mi, j] <- rnorm(length(mi), mean = mean(ob) - downshift * sd(ob), sd = sd(ob) * width) }
  x
}
wb_cols <- intersect(names(imp), o$sample_id)                        # workbook column order, as the impute step drew them
unimp <- mn[, wb_cols]; unimp[!obs[, wb_cols]] <- NA
re43 <- impute_normal(unimp, 43L)                                      # neuron_neuropil is the 2nd sorted celltype_layer: seed 42 + 1
w(tibble(check = "re-imputation with the pipeline seed reproduces the workbook", max_abs_difference = max(abs(re43 - mn[, wb_cols]))), "imputation_reproduction.csv")
seed_rows <- list()
for (sd_ in c(43L, 1001:1019)) {
  re <- impute_normal(unimp, sd_)[, o$sample_id]
  for (v in names(AM)) for (ex in c("all animals", "without A111")) {
    keep <- if (ex == "all animals") rep(TRUE, nrow(o)) else !drop_lab("A111")
    sm <- if (v == "original") sweep(re, 2, shift_orig, `+`) else re
    r <- fit_contrast(to_animal(sm, keep))
    seed_rows[[length(seed_rows) + 1L]] <- tibble(seed = sd_, centring = v, exclusion = ex, n_fdr05 = is_hit(r), n_fdr05_down = sum(r$adj.P.Val < 0.05 & r$logFC < 0))
  }
}
seed_tab <- bind_rows(seed_rows); w(seed_tab, "imputation_seeds.csv")

# ---------------------------------------------------------------- T. representative proteins
base_o <- runs[["original | all animals | pipeline"]]; noA_o <- runs[["original | without A111 | pipeline"]]
only_noA <- noA_o |> filter(adj.P.Val < 0.05, !protein %in% base_o$protein[base_o$adj.P.Val < 0.05]) |> arrange(P.Value)
reps <- unique(c("MAP7_MOUSE", "DPP8_MOUSE", "LUZP1_MOUSE", "SYT2_MOUSE", "GCC2_MOUSE", "ASAH1_MOUSE", "GGCT_MOUSE", "ABHD6_MOUSE",
                 head(only_noA$protein, 5), only_noA$protein[round(seq(0.25, 0.95, length.out = 4) * nrow(only_noA))]))
rep_tab <- bind_rows(lapply(names(runs), function(k) { p <- strsplit(k, " | ", fixed = TRUE)[[1]]
  runs[[k]] |> filter(protein %in% reps) |> transmute(protein, centring = p[1], exclusion = p[2], model = p[3], logFC, SE, t, P.Value, adj.P.Val) }))
ca2s <- o$unit == U
vals <- bind_rows(lapply(names(AM), function(v) { A <- AM[[v]][reps, col_unit == U, drop = FALSE]; colnames(A) <- col_animal[col_unit == U]
  as_tibble(A, rownames = "protein") |> pivot_longer(-protein, names_to = "animal", values_to = "animal_level_value") |> mutate(centring = v) })) |>
  left_join(tibble(animal = c(SUSA, RESA, "A127", "A129", "A765"), group = rep(c("SUS", "RES", "CON"), each = 3)), by = "animal") |>
  left_join(bind_rows(lapply(reps, function(p) tibble(protein = p, animal = o$animal[ca2s], hemisphere = o$hemisphere[ca2s], observed = obs[p, ca2s]))) |>
              group_by(protein, animal) |> summarise(n_hemispheres_imputed = sum(!observed), .groups = "drop"), by = c("protein", "animal"))
w(rep_tab, "representative_proteins.csv"); w(vals, "representative_protein_values.csv")
hit_detail <- noA_o |> filter(adj.P.Val < 0.05) |> mutate(in_base_hits = protein %in% base_o$protein[base_o$adj.P.Val < 0.05],
  n_imputed_A755_A764 = rowSums(!obs[protein, ca2s & o$animal %in% c("A755", "A764"), drop = FALSE]),
  n_imputed_A755L_A764R = rowSums(!obs[protein, o$label %in% c("A755_L_CA2_slm", "A764_R_CA2_slm"), drop = FALSE]),
  ca2_only_P = runs[["original | without A111 | ca2_only"]]$P.Value[match(protein, runs[["original | without A111 | ca2_only"]]$protein)],
  unit_weight_adjP = runs[["original | without A111 | unit_weights"]]$adj.P.Val[match(protein, runs[["original | without A111 | unit_weights"]]$protein)])
w(hit_detail, "hits_without_A111_original.csv")

options(width = 220)
print(as.data.frame(ca2_q |> select(label, group, plate, missing_fraction, Proteins.Identified, rank_missing_in_neuropil, z_missing_fraction, z_Total.Quantity, offset_difference)))
print(as.data.frame(prof))
print(as.data.frame(var_tab |> filter(unit == U | variant == "all animals") |> select(centring, variant, unit, median_ratio_unit_to_pooled_var, share_proteins_unit_var_gt_2x_pooled)))
print(as.data.frame(msum))
print(as.data.frame(relab |> arrange(centring, exclusion, -n_fdr05)))
print(as.data.frame(loao |> filter(removed != "none") |> group_by(centring, unit) |> summarise(base = first(base), min = min(n_fdr05), max = max(n_fdr05), max_by = removed[which.max(n_fdr05)], .groups = "drop")))
print(bind_rows(off_rows)); print(bind_rows(filt_rows))
print(seed_tab |> group_by(centring, exclusion) |> summarise(seed43 = n_fdr05[seed == 43], min = min(n_fdr05), median = median(n_fdr05), max = max(n_fdr05), .groups = "drop"))
print(hit_detail |> summarise(n = n(), in_base = sum(in_base_hits), any_imp_A755_A764 = sum(n_imputed_A755_A764 > 0), any_imp_lowdepth = sum(n_imputed_A755L_A764R > 0),
                              ca2only_p05 = sum(ca2_only_P < 0.05), unitw_fdr05 = sum(unit_weight_adjP < 0.05)))
cat("done:", OUT, "\n")
