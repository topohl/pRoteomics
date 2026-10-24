# Read-only validation of the per-sample centring applied before imputation (normalisation review 2026-10-21).
#
# Compares, on the raw DIA-NN protein-group matrix, the original centring (median of all observed values of a
# sample), the 2026-10-19 centring (median of the proteins quantified in every sample of the compartment) and
# independent estimators: a median-ratio estimator on the same complete-case proteins, a median-ratio estimator on
# proteins observed in at least half of the samples, the cross-compartment joint shared-core median already used by
# the QC stage, and no additional centring (DIA-NN output as delivered).
#
# Because the 70% filter uses missingness only and the down-shift imputation is relative to each sample's own
# observed mean and SD, every downstream matrix under estimator k equals the current matrix plus one constant per
# sample (o_new - o_k). This is verified against the archived original matrices before any effect is computed, so
# the comparison needs no pipeline rerun.
#
# Writes tables to results/reviewer_audit/normalisation_choice_2026-10-21/validation/. Changes nothing else.
# USAGE  Rscript tools/normalisation_review/validate_sample_centring.R          (working directory: repo root)

suppressPackageStartupMessages({
  library(readr); library(readxl); library(dplyr); library(tidyr); library(limma); library(matrixStats)
})
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
OUT <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "validation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(x, name) write_csv(x, file.path(OUT, name), na = "")
DATASETS <- c("microglia", "neuron_neuropil", "neuron_soma")
N_SPLITS <- 50L
FOCUS <- c("A755_L_CA2_slm", "A764_R_CA2_slm")
EST <- c("original_all_observed", "new_complete_case", "median_ratio_complete_case", "median_ratio_shared50",
         "joint_shared_core", "none_diann_only")

# ------------------------------------------------------------------------------------------------ 1. inputs
raw <- read_tsv("data/raw/pg_matrix/quicksearch.pg_matrix.tsv", show_col_types = FALSE,
                locale = locale(decimal_mark = "."), progress = FALSE)
X_all <- log2(as.matrix(raw[, -(1:4)])); rownames(X_all) <- raw$Protein.Group
stopifnot(!anyDuplicated(rownames(X_all)), !any(is.infinite(X_all)))
meta <- read_excel("data/metadata/TPE9_sample_metadata_males.xlsx") |> filter(is.na(exclude) | exclude != TRUE)
dia <- read_excel("data/raw/pg_matrix/quicksearch.stats.annotated.xlsx") |>
  select(sample_id, Precursors.Identified, Proteins.Identified, Total.Quantity, MS1.Signal, MS2.Signal, FWHM.RT,
         Median.Mass.Acc.MS1.Corrected, Median.Mass.Acc.MS2.Corrected, Normalisation.Instability, Median.RT.Prediction.Acc)
smp_all <- meta |>
  transmute(sample_id, dataset = celltype_layer, animal = AnimalID,
            group = unname(c("1" = "CON", "2" = "RES", "3" = "SUS")[as.character(ExpGroup)]),
            hemisphere = substr(ReplicateGroup, 1, 1), unit = paste(region, layer, sep = "_"),
            plate, well = sample_location, run_order = as.numeric(run_order),
            acquisition_date = sub(".*Bluto_([0-9]{8})_.*", "\\1", sample_id)) |>
  left_join(dia, by = "sample_id") |>
  group_by(plate) |> mutate(plate_position = rank(run_order)) |> ungroup() |>
  mutate(label = paste(animal, hemisphere, unit, sep = "_"))
stopifnot(all(smp_all$sample_id %in% colnames(X_all)))
w(dplyr::count(smp_all, plate, acquisition_date, group, animal) |> pivot_wider(names_from = plate, values_from = n, values_fill = 0),
  "design_plate_by_animal.csv")

# joint shared-core offsets (QC stage): unnormalised minus normalised, constant per sample on observed cells
joint <- readRDS("results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds")
jd <- joint$primary$unnormalized_log2 - joint$primary$matrix
jd[is.na(joint$primary$unnormalized_log2)] <- NA
joint_offset <- setNames(colMedians(jd, na.rm = TRUE), joint$metadata$sample_id)
joint_spread <- max(colMaxs(jd, na.rm = TRUE) - colMins(jd, na.rm = TRUE))

# ------------------------------------------------------------------------------------------------ 2. estimators
estimate_offsets <- function(X, exclude = character(0), joint = NULL) {
  Xk <- X[!(rownames(X) %in% exclude), , drop = FALSE]
  Xk <- Xk[rowSums(!is.na(Xk)) > 0, , drop = FALSE]
  cc <- rowSums(is.na(Xk)) == 0
  sh <- rowMeans(!is.na(Xk)) >= 0.5
  Xc <- Xk[cc, , drop = FALSE]; Xs <- Xk[sh, , drop = FALSE]
  out <- list(
    original_all_observed = colMedians(Xk, na.rm = TRUE),
    new_complete_case = colMedians(Xc),
    median_ratio_complete_case = colMedians(Xc - rowMeans(Xc)),
    median_ratio_shared50 = colMedians(Xs - rowMedians(Xs, na.rm = TRUE), na.rm = TRUE),
    joint_shared_core = if (is.null(joint)) rep(NA_real_, ncol(X)) else unname(joint[colnames(X)]),
    none_diann_only = rep(0, ncol(X)))
  lapply(out, function(v) setNames(v, colnames(X)))
}
rel <- function(v) v - mean(v)

# ------------------------------------------------------------------------------------------------ 3. per compartment
read_imputed <- function(path) {
  x <- read_excel(path, progress = FALSE)
  m <- as.matrix(x[, grepl("[.]d$", names(x))]); rownames(m) <- x$Protein.Group; m
}
read_gct <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dims <- as.integer(strsplit(lines[[2]], "\t", fixed = TRUE)[[1]][1:4])
  hdr <- strsplit(lines[[3]], "\t", fixed = TRUE)[[1]]
  desc <- do.call(rbind, strsplit(lines[4:(3 + dims[4])], "\t", fixed = TRUE))
  body <- do.call(rbind, strsplit(lines[(4 + dims[4]):(3 + dims[4] + dims[1])], "\t", fixed = TRUE))
  sc <- (dims[3] + 2L):length(hdr)
  m <- matrix(as.numeric(body[, sc]), nrow = nrow(body), dimnames = list(body[, 1], hdr[sc]))
  list(m = m, coldesc = setNames(lapply(seq_len(nrow(desc)), function(i) desc[i, sc]), desc[, 1]))
}
fit_da <- function(m, cell) {
  f <- factor(cell, levels = unique(cell)); design <- model.matrix(~0 + f); colnames(design) <- levels(f)
  units <- unique(sub("_[0-9]+$", "", cell)); pairs <- list(c("2", "1"), c("3", "2"), c("3", "1"))
  nm <- unlist(lapply(pairs, function(p) paste0(units, "_", p[1], "_over_", units, "_", p[2])))
  ex <- unlist(lapply(pairs, function(p) paste0(units, "_", p[1], "-", units, "_", p[2])))
  cm <- makeContrasts(contrasts = ex, levels = design); colnames(cm) <- nm
  fit <- lmFit(m, design)
  eb <- eBayes(contrasts.fit(fit, cm), robust = TRUE, trend = FALSE)
  list(eb = eb, sigma = fit$sigma, contrasts = nm,
       tab = lapply(setNames(nm, nm), function(cc) topTable(eb, coef = cc, number = Inf, sort.by = "none", adjust.method = "BH")))
}
contrast_label <- function(cc) {
  u <- sub("_[0-9]_over_.*$", "", cc); p <- regmatches(cc, regexpr("[0-9]_over_.*_[0-9]$", cc))
  a <- substr(p, 1, 1); b <- substr(p, nchar(p), nchar(p)); g <- c("1" = "CON", "2" = "RES", "3" = "SUS")
  data.frame(unit = u, comparison = paste(g[a], "vs", g[b]))
}

offsets_tab <- list(); verify <- list(); heldout <- list(); partial <- list(); qc_r2 <- list(); grp <- list()
da_counts <- list(); da_shift <- list(); var_tab <- list(); hemi <- list(); dist_tab <- list(); anchor_tab <- list()
intensity <- list(); ca2 <- list(); anchor_sets <- list(); anchor_group <- list(); da_ca2_sens <- list()

for (ds in DATASETS) {
  cat("==", ds, "\n")
  smp <- filter(smp_all, dataset == ds)
  X <- X_all[, smp$sample_id, drop = FALSE]
  X <- X[rowSums(!is.na(X)) > 0, , drop = FALSE]
  det <- rowMeans(!is.na(X)); cc <- det == 1
  pair_key <- paste(smp$animal, smp$unit)
  anchor_sets[[ds]] <- rownames(X)[cc]
  smp$n_detected <- colSums(!is.na(X)); smp$missing_fraction <- colMeans(is.na(X))
  off <- estimate_offsets(X, joint = joint_offset)

  # ---- 3a. offsets and their concordance
  ot <- bind_cols(smp, as_tibble(lapply(off, rel)) |> rename_with(~ paste0("rel_", .x)),
                  as_tibble(off) |> rename_with(~ paste0("raw_", .x))) |>
    mutate(original_minus_new = raw_original_all_observed - raw_new_complete_case)
  offsets_tab[[ds]] <- ot

  # ---- 3b. verification: the current and archived matrices differ by exactly these offsets
  f_new <- Sys.glob(sprintf("data/processed/01_preprocessing/impute/*_pgmatrix_imputed_%s_*_missing70pct.xlsx", ds))
  f_old <- Sys.glob(sprintf("%s/data/processed/01_preprocessing/impute/*_pgmatrix_imputed_%s_*_missing70pct.xlsx", ARCHIVE, ds))
  stopifnot(length(f_new) == 1L)
  mn <- read_imputed(f_new)[, smp$sample_id]
  obs <- !is.na(X[rownames(mn), smp$sample_id])
  ref_new <- sweep(X[rownames(mn), smp$sample_id], 2, off$new_complete_case)
  v <- tibble(dataset = ds, check = "current imputed matrix = raw log2 minus complete-case median (observed cells)",
              file = f_new, max_abs_difference = max(abs(mn[obs] - ref_new[obs])))
  for (fo in f_old) {
    mo <- read_imputed(fo)[rownames(mn), smp$sample_id]
    expected <- sweep(mn, 2, off$new_complete_case - off$original_all_observed, `+`)
    d <- mo - mn
    v <- bind_rows(v,
      tibble(dataset = ds, check = "archived original matrix = current matrix + (new - original offset), all cells incl. imputed",
             file = fo, max_abs_difference = max(abs(mo - expected))),
      tibble(dataset = ds, check = "archived minus current is one constant per sample (max within-sample spread)",
             file = fo, max_abs_difference = max(colMaxs(d) - colMins(d))))
  }
  # animal-level GCT: archived = current + mean offset change of the animal's hemispheres
  g_new <- read_gct(sprintf("data/processed/01_preprocessing/protigy_input_animal_level/%s/%s_animal_level.gct", ds, ds))
  g_old_path <- sprintf("%s/data/processed/01_preprocessing/protigy_input_animal_level/%s/%s_animal_level.gct", ARCHIVE, ds, ds)
  col_key <- paste(g_new$coldesc$AnimalID, g_new$coldesc$region, g_new$coldesc$layer, sep = "|")
  smp_key <- paste(smp$animal, sub("_.*$", "", smp$unit), sub("^[^_]+_", "", smp$unit), sep = "|")
  agg <- function(o) vapply(col_key, function(k) mean(o[smp_key == k]), numeric(1))
  stopifnot(all(col_key %in% smp_key))
  if (file.exists(g_old_path)) {
    g_old <- read_gct(g_old_path)
    stopifnot(identical(colnames(g_old$m), colnames(g_new$m)))
    rid <- intersect(rownames(g_new$m), rownames(g_old$m))
    exp_old <- sweep(g_new$m[rid, ], 2, agg(off$new_complete_case - off$original_all_observed), `+`)
    v <- bind_rows(v, tibble(dataset = ds, check = "archived animal-level input = current + mean hemisphere offset change",
                             file = g_old_path, max_abs_difference = max(abs(g_old$m[rid, ] - exp_old))))
  }
  verify[[ds]] <- v

  # ---- 3c. held-out benchmark: see section 3c below the loop

  # ---- 3d. pre-centring QC: what explains the raw sample levels and the original-new difference
  qc <- ot |> mutate(log10_MS2 = log10(MS2.Signal), log10_MS1 = log10(MS1.Signal), unit = factor(unit),
                     animal = factor(animal), plate = factor(plate))
  responses <- c(raw_new_complete_case = "complete-case median before centring (DIA-NN level)",
                 raw_original_all_observed = "all-observed median before centring",
                 original_minus_new = "original minus new offset",
                 missing_fraction = "missing fraction")
  predictors <- c("missing_fraction", "log10_MS2", "log10_MS1", "Proteins.Identified", "Normalisation.Instability",
                  "FWHM.RT", "Median.RT.Prediction.Acc", "plate", "plate_position", "hemisphere", "unit", "animal", "group")
  for (y in names(responses)) {
    for (p in setdiff(predictors, y)) {
      fml <- stats::as.formula(paste(y, "~", p))
      fit1 <- lm(fml, data = qc)
      # partial R2 of the predictor given unit (the region/layer composition), except for unit itself
      if (p != "unit") {
        full <- lm(stats::as.formula(paste(y, "~ unit +", p)), data = qc); red <- lm(stats::as.formula(paste(y, "~ unit")), data = qc)
        pr2 <- (sum(resid(red)^2) - sum(resid(full)^2)) / sum(resid(red)^2)
      } else pr2 <- NA_real_
      qc_r2[[length(qc_r2) + 1L]] <- tibble(dataset = ds, response = y, response_description = responses[[y]],
        predictor = p, r2_alone = summary(fit1)$r.squared, partial_r2_given_unit = pr2,
        slope_if_numeric = if (is.numeric(qc[[p]])) unname(coef(fit1)[2]) else NA_real_)
    }
  }

  # ---- 3e. condition dependence of the offsets (animal level, unit-adjusted)
  for (k in EST) {
    if (all(is.na(off[[k]]))) next
    a <- ot |> mutate(o = rel(off[[k]])) |> group_by(unit) |> mutate(o = o - mean(o)) |> ungroup() |>
      group_by(animal, group) |> summarise(o = mean(o), .groups = "drop")
    mm <- tapply(a$o, a$group, mean)
    grp[[length(grp) + 1L]] <- tibble(dataset = ds, estimator = k,
      CON = mm[["CON"]], RES = mm[["RES"]], SUS = mm[["SUS"]],
      SUS_minus_RES = mm[["SUS"]] - mm[["RES"]], SUS_minus_CON = mm[["SUS"]] - mm[["CON"]], RES_minus_CON = mm[["RES"]] - mm[["CON"]],
      anova_p_animal_level = anova(lm(o ~ group, data = a))[["Pr(>F)"]][1],
      cor_offset_vs_missing = stats::cor(off[[k]], smp$missing_fraction),
      cor_offset_vs_log10_MS2 = stats::cor(off[[k]], log10(smp$MS2.Signal)))
  }
  # alternative anchor sets: complete within one group only
  for (gname in c("CON", "RES", "SUS")) {
    sel <- smp$group == gname
    cc_g <- rowSums(is.na(X[, sel, drop = FALSE])) == 0
    o_g <- colMedians(X[cc_g, , drop = FALSE], na.rm = TRUE)
    anchor_group[[length(anchor_group) + 1L]] <- tibble(dataset = ds, anchor_definition = paste("complete in", gname, "samples"),
      n_anchor = sum(cc_g), n_shared_with_all_sample_anchor = sum(cc_g & cc),
      cor_with_new_offset = stats::cor(o_g, off$new_complete_case),
      max_abs_difference_relative_offset = max(abs(rel(o_g) - rel(off$new_complete_case))))
  }

  # ---- 3f. anchor characteristics
  mean_int <- rowMeans(X, na.rm = TRUE)
  anchor_tab[[ds]] <- tibble(dataset = ds, n_samples = ncol(X), n_proteins_any = nrow(X),
    n_proteins_after_70pct_filter = sum(det >= 0.3), n_anchor = sum(cc),
    anchor_fraction_of_filtered = sum(cc) / sum(det >= 0.3),
    anchor_mean_log2_q10 = quantile(mean_int[cc], 0.1), anchor_mean_log2_median = median(mean_int[cc]),
    anchor_mean_log2_q90 = quantile(mean_int[cc], 0.9),
    nonanchor_mean_log2_q10 = quantile(mean_int[!cc & det >= 0.3], 0.1), nonanchor_mean_log2_median = median(mean_int[!cc & det >= 0.3]),
    nonanchor_mean_log2_q90 = quantile(mean_int[!cc & det >= 0.3], 0.9),
    anchor_intensity_percentile_median = mean(mean_int[!cc] < median(mean_int[cc])),
    median_detection_nonanchor_filtered = median(det[!cc & det >= 0.3]))

  # ---- 3g. intensity dependence of the residual sample effect (is one constant per sample adequate?)
  Yc <- sweep(X[cc, ], 2, off$new_complete_case); Yc <- Yc - rowMeans(Yc)
  terc <- cut(mean_int[cc], quantile(mean_int[cc], c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low", "mid", "high"))
  tm <- sapply(levels(terc), function(t) colMedians(Yc[terc == t, , drop = FALSE]))
  intensity[[ds]] <- tibble(dataset = ds, label = smp$label, sample_id = smp$sample_id, missing_fraction = smp$missing_fraction,
    low_tertile = tm[, "low"], mid_tertile = tm[, "mid"], high_tertile = tm[, "high"],
    high_minus_low = tm[, "high"] - tm[, "low"], original_minus_new_relative = rel(off$original_all_observed - off$new_complete_case))

  # ---- 3h. sample distributions after centring, hemisphere concordance
  for (k in EST) {
    if (all(is.na(off[[k]]))) next
    Y <- sweep(X, 2, off[[k]])
    med <- colMedians(Y, na.rm = TRUE); iqr <- colIQRs(Y, na.rm = TRUE)
    Ycc <- Y[cc, ]
    dist_tab[[length(dist_tab) + 1L]] <- tibble(dataset = ds, estimator = k,
      sd_sample_median_all_observed = sd(med), cor_sample_median_vs_missing = stats::cor(med, smp$missing_fraction),
      sd_sample_median_anchor = sd(colMedians(Ycc)), median_iqr = median(iqr))
    lr <- lapply(split(seq_len(nrow(smp)), pair_key), function(ix) if (length(ix) == 2) ix[order(smp$hemisphere[ix])] else NULL)
    lr <- Filter(Negate(is.null), lr)
    gap_cc <- vapply(lr, function(ix) median(Ycc[, ix[1]] - Ycc[, ix[2]]), numeric(1))
    both <- lapply(lr, function(ix) { b <- !is.na(Y[, ix[1]]) & !is.na(Y[, ix[2]]) & !cc; median(Y[b, ix[1]] - Y[b, ix[2]]) })
    rr <- vapply(lr, function(ix) { b <- !is.na(Y[, ix[1]]) & !is.na(Y[, ix[2]]); stats::cor(Y[b, ix[1]], Y[b, ix[2]]) }, numeric(1))
    hemi[[length(hemi) + 1L]] <- tibble(dataset = ds, estimator = k, n_pairs = length(lr),
      median_abs_LR_gap_anchor = median(abs(gap_cc)), median_abs_LR_gap_nonanchor_observed_in_both = median(abs(unlist(both))),
      p90_abs_LR_gap_nonanchor = unname(quantile(abs(unlist(both)), 0.9)),
      median_LR_pearson = median(rr), note = "Pearson r between hemispheres is invariant to any per-sample constant")
  }

  # ---- 3i. downstream DA under each estimator: animal-level matrix + per-column shift, pipeline limma
  cell <- g_new$coldesc$region_layer_ExpGroup
  base <- fit_da(g_new$m, cell)
  if (ds == "neuron_neuropil") {
    # sample-level rebuild of the animal-level matrix; proves the reconstruction, then drops the focus hemispheres
    keep_focus <- !(smp$label %in% FOCUS)
    pn <- raw$Protein.Names[match(rownames(mn), raw$Protein.Group)]
    rebuild <- function(keep) {
      m <- sapply(col_key, function(kk) rowMeans(mn[, smp$sample_id[smp_key == kk & keep], drop = FALSE]))
      rownames(m) <- pn; m <- m[!duplicated(pn), , drop = FALSE]; m[rownames(g_new$m), , drop = FALSE]
    }
    m_full <- rebuild(rep(TRUE, nrow(smp)))
    verify[[ds]] <- bind_rows(verify[[ds]], tibble(dataset = ds,
      check = "current animal-level input = mean of the animal's hemispheres in the current imputed matrix",
      file = "data/processed/01_preprocessing/protigy_input_animal_level/neuron_neuropil/neuron_neuropil_animal_level.gct",
      max_abs_difference = max(abs(m_full - g_new$m))))
    m_drop <- rebuild(keep_focus)
  }
  for (k in EST) {
    if (anyNA(off[[k]])) next
    shift <- agg(off$new_complete_case - off[[k]])
    m_k <- sweep(g_new$m, 2, shift, `+`)
    r <- if (identical(k, "new_complete_case")) base else fit_da(m_k, cell)
    var_tab[[length(var_tab) + 1L]] <- tibble(dataset = ds, estimator = k,
      median_residual_sd = median(r$sigma), mean_residual_var = mean(r$sigma^2),
      s2_prior = r$eb$s2.prior[1], df_prior_median = median(r$eb$df.prior))
    for (ct in r$contrasts) {
      t1 <- r$tab[[ct]]; t0 <- base$tab[[ct]]
      lab <- contrast_label(ct)
      da_counts[[length(da_counts) + 1L]] <- tibble(dataset = ds, estimator = k, contrast = ct, unit = lab$unit,
        comparison = lab$comparison, n_fdr05 = sum(t1$adj.P.Val < 0.05), n_up = sum(t1$adj.P.Val < 0.05 & t1$logFC > 0),
        n_down = sum(t1$adj.P.Val < 0.05 & t1$logFC < 0), median_logFC = median(t1$logFC))
      dl <- t1$logFC - t0$logFC
      da_shift[[length(da_shift) + 1L]] <- tibble(dataset = ds, estimator = k, contrast = ct, unit = lab$unit, comparison = lab$comparison,
        logFC_shift_vs_new = median(dl), logFC_shift_spread = max(dl) - min(dl),
        cor_t_vs_new = stats::cor(t1$t, t0$t), spearman_t_vs_new = stats::cor(t1$t, t0$t, method = "spearman"),
        n_fdr05_shared_with_new = sum(t1$adj.P.Val < 0.05 & t0$adj.P.Val < 0.05))
    }
    # CA2-SLM sensitivity: A755 and A764 represented by their other hemisphere only
    if (ds == "neuron_neuropil") {
      agg_keep <- vapply(col_key, function(kk) mean((off$new_complete_case - off[[k]])[smp_key == kk & keep_focus]), numeric(1))
      m_s <- sweep(m_drop, 2, agg_keep, `+`)
      rs <- fit_da(m_s, cell)
      for (ct in grep("^CA2_slm_", rs$contrasts, value = TRUE)) {
        t1 <- rs$tab[[ct]]; lab <- contrast_label(ct)
        da_ca2_sens[[length(da_ca2_sens) + 1L]] <- tibble(estimator = k, contrast = ct, comparison = lab$comparison,
          variant = "A755 and A764 represented by their other hemisphere", n_fdr05 = sum(t1$adj.P.Val < 0.05),
          n_fdr05_full = sum(r$tab[[ct]]$adj.P.Val < 0.05), cor_logFC_vs_full = stats::cor(t1$logFC, r$tab[[ct]]$logFC))
      }
    }
  }

  # ---- 3j. CA2-SLM sample detail (neuropil)
  if (ds == "neuron_neuropil") {
    rz <- function(x) (x - median(x)) / mad(x)
    cons <- rowMedians(sweep(X[cc, smp$unit == "CA2_slm"], 2, off$new_complete_case[smp$unit == "CA2_slm"]))
    ca2[[ds]] <- ot |> mutate(
        z_proteins_identified_neuropil = rz(Proteins.Identified), z_missing_neuropil = rz(missing_fraction),
        z_log10_MS2_neuropil = rz(log10(MS2.Signal)), z_norm_instability_neuropil = rz(Normalisation.Instability),
        z_original_offset_neuropil = rz(rel_original_all_observed), z_new_offset_neuropil = rz(rel_new_complete_case),
        z_median_ratio_offset_neuropil = rz(rel_median_ratio_complete_case)) |>
      filter(unit == "CA2_slm") |>
      mutate(r_with_CA2slm_consensus_anchor = vapply(sample_id, function(s) stats::cor(X[cc, s], cons), numeric(1))) |>
      select(label, sample_id, group, plate, run_order, plate_position, Proteins.Identified, Precursors.Identified, MS2.Signal,
             Normalisation.Instability, missing_fraction, starts_with("z_"), starts_with("rel_"), original_minus_new,
             r_with_CA2slm_consensus_anchor)
  }
}

# The held-out loop above appended into a copy; rebuild it explicitly here so the tables are complete.
heldout <- list(); partial <- list(); focus_rows <- list()
for (ds in DATASETS) {
  smp <- filter(smp_all, dataset == ds)
  X <- X_all[, smp$sample_id, drop = FALSE]; X <- X[rowSums(!is.na(X)) > 0, , drop = FALSE]
  det <- rowMeans(!is.na(X)); mf <- colMeans(is.na(X))
  anchors <- anchor_sets[[ds]]; part_rows <- rownames(X)[det >= 0.5 & det < 0.95]
  pair_key <- paste(smp$animal, smp$unit)
  set.seed(20261021 + match(ds, DATASETS))
  for (s in seq_len(N_SPLITS)) {
    B <- sample(anchors, floor(length(anchors) / 2)); PB <- sample(part_rows, floor(length(part_rows) / 2))
    o <- estimate_offsets(X, exclude = c(B, PB), joint = joint_offset)
    for (k in EST) {
      if (all(is.na(o[[k]]))) next
      for (set in c("held_out_anchor", "held_out_partially_observed")) {
        rows <- if (set == "held_out_anchor") B else PB
        Y <- sweep(X[rows, , drop = FALSE], 2, o[[k]]); Y <- Y - rowMeans(Y, na.rm = TRUE)
        e <- colMedians(Y, na.rm = TRUE)
        pc <- if (set == "held_out_anchor") median(apply(Y, 1, function(y) stats::cor(y, mf))) else NA_real_
        hl <- tapply(e, pair_key, function(z) if (length(z) == 2) abs(diff(z)) else NA_real_)
        row <- tibble(dataset = ds, split = s, estimator = k, set = set, residual_sample_sd = sd(e),
                      cor_residual_vs_missing = stats::cor(e, mf), median_protein_cor_vs_missing = pc,
                      mean_hemisphere_pair_gap = mean(hl, na.rm = TRUE))
        if (set == "held_out_anchor") heldout[[length(heldout) + 1L]] <- row else partial[[length(partial) + 1L]] <- row
        fl <- intersect(FOCUS, smp$label)
        if (length(fl)) focus_rows[[length(focus_rows) + 1L]] <- tibble(dataset = ds, split = s, estimator = k, set = set,
          label = fl, residual = unname(e[match(fl, smp$label)]), residual_z = unname(((e - median(e)) / mad(e))[match(fl, smp$label)]))
      }
    }
  }
}

# ------------------------------------------------------------------------------------------------ 4. DIA-NN level checks
gcc <- rowSums(is.na(X_all[, smp_all$sample_id])) == 0
diann <- smp_all |> mutate(global_anchor_median = colMedians(X_all[gcc, sample_id]),
                           all_observed_median = colMedians(X_all[, sample_id], na.rm = TRUE),
                           n_detected_global = colSums(!is.na(X_all[, sample_id])))
diann_summary <- diann |> group_by(dataset) |>
  summarise(n = n(), n_global_anchor = sum(gcc), sd_global_anchor_median = sd(global_anchor_median),
            sd_all_observed_median = sd(all_observed_median),
            cv_total_quantity = sd(Total.Quantity) / mean(Total.Quantity), cv_MS2_signal = sd(MS2.Signal) / mean(MS2.Signal),
            cv_MS1_signal = sd(MS1.Signal) / mean(MS1.Signal),
            cor_anchor_median_vs_log10_MS2 = stats::cor(global_anchor_median, log10(MS2.Signal)),
            cor_total_quantity_vs_MS2 = stats::cor(Total.Quantity, MS2.Signal),
            median_norm_instability = median(Normalisation.Instability), max_norm_instability = max(Normalisation.Instability),
            .groups = "drop") |>
  bind_rows(diann |> summarise(dataset = "all compartments", n = n(), n_global_anchor = sum(gcc),
            sd_global_anchor_median = sd(global_anchor_median), sd_all_observed_median = sd(all_observed_median),
            cv_total_quantity = sd(Total.Quantity) / mean(Total.Quantity), cv_MS2_signal = sd(MS2.Signal) / mean(MS2.Signal),
            cv_MS1_signal = sd(MS1.Signal) / mean(MS1.Signal),
            cor_anchor_median_vs_log10_MS2 = stats::cor(global_anchor_median, log10(MS2.Signal)),
            cor_total_quantity_vs_MS2 = stats::cor(Total.Quantity, MS2.Signal),
            median_norm_instability = median(Normalisation.Instability), max_norm_instability = max(Normalisation.Instability)))

# ------------------------------------------------------------------------------------------------ 5. write
w(bind_rows(offsets_tab), "sample_offsets_all_estimators.csv")
off_long <- bind_rows(offsets_tab)
conc <- bind_rows(lapply(DATASETS, function(ds) {
  o <- filter(off_long, dataset == ds)
  cm <- expand.grid(a = EST, b = EST, stringsAsFactors = FALSE) |> filter(a < b)
  bind_rows(lapply(seq_len(nrow(cm)), function(i) {
    x <- o[[paste0("rel_", cm$a[i])]]; y <- o[[paste0("rel_", cm$b[i])]]
    if (all(is.na(x)) || all(is.na(y))) return(NULL)
    tibble(dataset = ds, estimator_a = cm$a[i], estimator_b = cm$b[i], pearson_r = stats::cor(x, y),
           max_abs_difference = max(abs(x - y)), sd_difference = sd(x - y))
  }))
}))
w(conc, "offset_concordance_between_estimators.csv")
w(bind_rows(verify) |> bind_rows(tibble(dataset = "global", check = "joint shared-core offset is one constant per sample",
                                        file = "results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds",
                                        max_abs_difference = joint_spread)), "verification_additive_offsets.csv")
ho <- bind_rows(heldout, partial)
w(ho, "heldout_benchmark_splits.csv")
w(ho |> group_by(dataset, set, estimator) |>
    summarise(splits = n(), residual_sample_sd = mean(residual_sample_sd), cor_residual_vs_missing = mean(cor_residual_vs_missing),
              median_protein_cor_vs_missing = mean(median_protein_cor_vs_missing), mean_hemisphere_pair_gap = mean(mean_hemisphere_pair_gap),
              .groups = "drop") |> arrange(dataset, set, residual_sample_sd), "heldout_benchmark_summary.csv")
w(bind_rows(focus_rows) |> group_by(dataset, set, estimator, label) |>
    summarise(mean_residual = mean(residual), mean_residual_z = mean(residual_z), .groups = "drop"), "heldout_focus_samples.csv")
w(bind_rows(qc_r2), "precentring_qc_variance_explained.csv")
w(bind_rows(grp), "offset_condition_dependence.csv")
w(bind_rows(anchor_group), "anchor_definition_by_group.csv")
w(bind_rows(anchor_tab), "anchor_set_characteristics.csv")
w(bind_rows(intensity), "intensity_dependence_of_sample_effect.csv")
w(bind_rows(intensity) |> group_by(dataset) |>
    summarise(sd_high_minus_low = sd(high_minus_low), cor_high_minus_low_vs_missing = stats::cor(high_minus_low, missing_fraction),
              cor_high_minus_low_vs_original_minus_new = stats::cor(high_minus_low, original_minus_new_relative),
              sd_original_minus_new = sd(original_minus_new_relative), .groups = "drop"), "intensity_dependence_summary.csv")
w(bind_rows(dist_tab), "sample_distribution_after_centring.csv")
w(bind_rows(hemi), "hemisphere_concordance.csv")
w(bind_rows(var_tab), "residual_variance_by_estimator.csv")
w(bind_rows(da_counts), "da_counts_by_estimator.csv")
w(bind_rows(da_shift), "da_shift_vs_new_by_estimator.csv")
w(bind_rows(da_ca2_sens), "ca2slm_sensitivity_focus_hemispheres.csv")
w(bind_rows(ca2), "ca2slm_sample_detail.csv")
w(diann_summary, "diann_level_summary.csv")
w(diann |> select(sample_id, label, dataset, plate, run_order, global_anchor_median, all_observed_median, n_detected_global,
                  Proteins.Identified, MS1.Signal, MS2.Signal, Total.Quantity, Normalisation.Instability), "diann_level_samples.csv")
saveRDS(list(anchor_sets = anchor_sets, offsets = off_long), file.path(OUT, "validation_state.rds"))
cat("done:", OUT, "\n")
