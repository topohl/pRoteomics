# Addendum to validate_sample_centring.R after the adversarial check of 2026-10-21 (read-only for the pipeline).
#
#   1  Held-out benchmark under three evaluation functionals (the first run used only the row-centred median, which is
#      the median-ratio statistic itself): row-centred median, raw median and row-centred mean of the held-out proteins.
#      Every estimator is scored re-estimated without the held-out proteins AND as deployed on all proteins (in-sample),
#      because re-estimating the all-observed median on a depleted protein set inflates its error. The joint shared-core
#      offset is fixed (estimated by the QC stage on all proteins), so its held-out score is in-sample: flagged.
#   2  Offset vs missing fraction within region-layer units (partial r), the way the group test is done.
#   3  Intensity dependence of the residual sample effect within units, and from hemisphere pairs.
#   4  Group differences of the offsets at animal level with plate as covariate (plate is unbalanced across groups).
#   5  Microglia: do the anchors carry neuronal/myelin carry-over? Offsets vs a neuron/myelin marker score, and offsets
#      from anchors restricted to microglia-enriched proteins (microglia vs neuropil in the CON animals).
#
# Writes validation/addendum_*.csv. USAGE  Rscript tools/normalisation_review/validate_sample_centring_addendum.R

suppressPackageStartupMessages({ library(readr); library(readxl); library(dplyr); library(tidyr); library(matrixStats) })
OUT <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "validation")
w <- function(x, name) write_csv(x, file.path(OUT, name), na = "")
DATASETS <- c("microglia", "neuron_neuropil", "neuron_soma")
N_SPLITS <- 50L
state <- readRDS(file.path(OUT, "validation_state.rds")); OFF <- state$offsets
raw <- read_tsv("data/raw/pg_matrix/quicksearch.pg_matrix.tsv", show_col_types = FALSE, progress = FALSE)
X_all <- log2(as.matrix(raw[, -(1:4)])); rownames(X_all) <- raw$Protein.Group
gene <- setNames(toupper(sub(";.*$", "", raw$Genes)), raw$Protein.Group)

estimate <- function(X, exclude = character(0)) {
  Xk <- X[!(rownames(X) %in% exclude), , drop = FALSE]; Xk <- Xk[rowSums(!is.na(Xk)) > 0, , drop = FALSE]
  cc <- rowSums(is.na(Xk)) == 0; sh <- rowMeans(!is.na(Xk)) >= 0.5
  Xc <- Xk[cc, , drop = FALSE]; Xs <- Xk[sh, , drop = FALSE]
  list(original_all_observed = colMedians(Xk, na.rm = TRUE), new_complete_case = colMedians(Xc),
       median_ratio_complete_case = colMedians(Xc - rowMeans(Xc)),
       median_ratio_shared50 = colMedians(Xs - rowMedians(Xs, na.rm = TRUE), na.rm = TRUE),
       none_diann_only = rep(0, ncol(X)))
}
functionals <- list(
  row_centred_median = function(Y) { Y <- Y - rowMeans(Y, na.rm = TRUE); colMedians(Y, na.rm = TRUE) },
  raw_median = function(Y) { e <- colMedians(Y, na.rm = TRUE); e - mean(e) },
  row_centred_mean = function(Y) { Y <- Y - rowMeans(Y, na.rm = TRUE); colMeans(Y, na.rm = TRUE) })

ho <- list(); part2 <- list(); part3 <- list(); part4 <- list()
for (ds in DATASETS) {
  o <- OFF[OFF$dataset == ds, ]
  X <- X_all[, o$sample_id, drop = FALSE]; X <- X[rowSums(!is.na(X)) > 0, , drop = FALSE]
  det <- rowMeans(!is.na(X)); mf <- colMeans(is.na(X)); anchors <- rownames(X)[det == 1]
  part_rows <- rownames(X)[det >= 0.5 & det < 0.95]
  pair_key <- paste(o$animal, o$unit)
  deployed <- estimate(X); deployed$joint_shared_core <- o$raw_joint_shared_core
  set.seed(20261021 + match(ds, DATASETS))
  for (s in seq_len(N_SPLITS)) {
    B <- sample(anchors, floor(length(anchors) / 2)); PB <- sample(part_rows, floor(length(part_rows) / 2))
    re <- estimate(X, exclude = c(B, PB)); re$joint_shared_core <- o$raw_joint_shared_core
    for (mode in c("re-estimated without held-out proteins", "as deployed (in-sample)")) {
      est <- if (startsWith(mode, "re")) re else deployed
      for (k in names(est)) for (set in c("anchor", "partially_observed")) {
        rows <- if (set == "anchor") B else PB
        Y <- sweep(X[rows, , drop = FALSE], 2, est[[k]])
        for (fn in names(functionals)) {
          e <- functionals[[fn]](Y)
          gap <- tapply(e, pair_key, function(z) if (length(z) == 2) abs(diff(z)) else NA_real_)
          ho[[length(ho) + 1L]] <- tibble(dataset = ds, split = s, mode = mode, estimator = k, set = set, functional = fn,
            residual_sample_sd = sd(e), cor_residual_vs_missing = stats::cor(e, mf), hemisphere_gap = mean(gap, na.rm = TRUE),
            in_sample = k == "joint_shared_core" || startsWith(mode, "as deployed"))
        }
      }
    }
  }
  # 2. partial r given unit
  for (k in c("original_all_observed", "new_complete_case", "median_ratio_complete_case", "median_ratio_shared50", "joint_shared_core")) {
    v <- o[[paste0("raw_", k)]]
    ru <- resid(lm(v ~ factor(o$unit))); mu <- resid(lm(mf ~ factor(o$unit)))
    part2[[length(part2) + 1L]] <- tibble(dataset = ds, estimator = k, r_pooled = stats::cor(v, mf), r_within_unit = stats::cor(ru, mu))
  }
  # 3. intensity dependence within unit and from hemisphere pairs
  Yc <- sweep(X[anchors, ], 2, o$raw_new_complete_case)
  mean_int <- rowMeans(X[anchors, ])
  Yw <- Yc
  for (u in unique(o$unit)) { j <- o$unit == u; Yw[, j] <- Yc[, j] - rowMeans(Yc[, j, drop = FALSE]) }
  terc <- cut(mean_int, quantile(mean_int, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low", "mid", "high"))
  hl <- colMedians(Yw[terc == "high", ]) - colMedians(Yw[terc == "low", ])
  pairs <- split(seq_len(nrow(o)), pair_key); pairs <- pairs[lengths(pairs) == 2]
  pair_sd <- sd(vapply(pairs, function(ix) diff(hl[ix]), numeric(1))) / sqrt(2)
  dn <- o$raw_original_all_observed - o$raw_new_complete_case
  part3[[length(part3) + 1L]] <- tibble(dataset = ds, sd_high_minus_low_within_unit = sd(hl), sd_from_hemisphere_pairs = pair_sd,
    cor_high_minus_low_vs_missing_within_unit = stats::cor(hl, resid(lm(mf ~ factor(o$unit)))),
    sd_original_minus_new_within_unit = sd(resid(lm(dn ~ factor(o$unit)))))
  # 4. group differences of offsets with plate as covariate (animal level, unit-centred)
  for (k in c("original_all_observed", "new_complete_case", "median_ratio_complete_case")) {
    a <- o |> mutate(v = .data[[paste0("raw_", k)]]) |> group_by(unit) |> mutate(v = v - mean(v)) |> ungroup() |>
      group_by(animal, group, plate) |> summarise(v = mean(v), .groups = "drop")
    f1 <- lm(v ~ group, a); f2 <- lm(v ~ plate + group, a)
    cf2 <- coef(f2)
    part4[[length(part4) + 1L]] <- tibble(dataset = ds, estimator = k,
      SUS_minus_RES = unname(coef(f1)["groupSUS"] - coef(f1)["groupRES"]),
      SUS_minus_RES_plate_adjusted = unname(cf2["groupSUS"] - cf2["groupRES"]),
      RES_minus_CON_plate_adjusted = unname(cf2["groupRES"]), SUS_minus_CON_plate_adjusted = unname(cf2["groupSUS"]),
      plate_effect = unname(cf2["plateC"]), p_group_given_plate = anova(lm(v ~ plate, a), f2)[["Pr(>F)"]][2])
  }
}
hos <- bind_rows(ho) |> group_by(dataset, set, functional, mode, estimator, in_sample) |>
  summarise(residual_sample_sd = mean(residual_sample_sd), cor_residual_vs_missing = mean(cor_residual_vs_missing),
            hemisphere_gap = mean(hemisphere_gap), .groups = "drop") |> arrange(dataset, set, functional, mode, residual_sample_sd)
w(hos, "addendum_heldout_by_functional.csv")
w(bind_rows(part2), "addendum_offset_vs_missing_within_unit.csv")
w(bind_rows(part3), "addendum_intensity_dependence_within_unit.csv")
w(bind_rows(part4), "addendum_group_offsets_plate_adjusted.csv")

# 5. microglia anchors and neuronal/myelin carry-over
o <- OFF[OFF$dataset == "microglia", ]
X <- X_all[, o$sample_id]; X <- X[rowSums(!is.na(X)) > 0, ]
anchors <- rownames(X)[rowSums(is.na(X)) == 0]
markers <- c("SNAP25", "SYT1", "STX1A", "GRIN1", "GRIA2", "SYN1", "NEFL", "INA", "MAP2", "CAMK2A", "MBP", "PLP1", "MOG", "MAG", "CNP")
mk <- intersect(rownames(X)[gene[rownames(X)] %in% markers], anchors)
Y <- sweep(X, 2, o$raw_new_complete_case)
score <- colMedians(Y[mk, ] - rowMeans(Y[mk, ]))
# uncentred carry-over: marker median minus the all-protein median on raw log2 (no dependence on any offset)
score_raw <- colMedians(X[mk, ]) - colMedians(X, na.rm = TRUE)
# microglia-enriched anchors: higher in microglia than in neuropil of the same CON animals (joint shared-core normalised)
joint <- readRDS("results/preprocessing/build_joint_protigy_input/global/models/joint_compartment_qc_matrices.rds")
jm <- joint$primary$unnormalized_log2 - matrix(colMedians(joint$primary$unnormalized_log2 - joint$primary$matrix, na.rm = TRUE),
                                               nrow(joint$primary$matrix), ncol(joint$primary$matrix), byrow = TRUE)
md <- joint$metadata
con <- md$experimental_group == "1"
fid <- joint$primary$feature_ids
jrows <- joint$feature_table$source_feature_id[match(fid, joint$feature_table$ProteinGroupID)]
cat("joint features mapped to raw protein groups:", sum(jrows %in% rownames(X_all)), "of", length(fid), "
")
enr <- rowMeans(jm[, con & md$dataset == "microglia"], na.rm = TRUE) - rowMeans(jm[, con & md$dataset == "neuron_neuropil"], na.rm = TRUE)
mg_enriched <- if (!is.null(jrows)) intersect(jrows[!is.na(enr) & enr > 1], anchors) else character()
o_mg <- if (length(mg_enriched) >= 50) colMedians(X[mg_enriched, ]) else rep(NA_real_, ncol(X))
part5 <- tibble(n_markers_in_anchors = length(mk), markers = paste(sort(unique(gene[mk])), collapse = " "),
  n_microglia_enriched_anchors = length(mg_enriched),
  cor_marker_score_vs_new_offset = stats::cor(score, o$raw_new_complete_case),
  cor_marker_score_vs_original_minus_new = stats::cor(score, o$raw_original_all_observed - o$raw_new_complete_case),
  cor_marker_score_vs_missing = stats::cor(score, o$missing_fraction),
  cor_uncentred_marker_score_vs_new_offset = stats::cor(score_raw, o$raw_new_complete_case),
  cor_uncentred_marker_score_vs_original_offset = stats::cor(score_raw, o$raw_original_all_observed),
  cor_uncentred_marker_score_vs_missing = stats::cor(score_raw, o$missing_fraction),
  cor_microglia_enriched_anchor_offset_vs_new = if (all(is.na(o_mg))) NA_real_ else stats::cor(o_mg, o$raw_new_complete_case),
  max_abs_diff_microglia_enriched_vs_new = if (all(is.na(o_mg))) NA_real_ else max(abs((o_mg - mean(o_mg)) - (o$raw_new_complete_case - mean(o$raw_new_complete_case)))))
w(part5, "addendum_microglia_anchor_carryover.csv")
print(as.data.frame(part5))
print(as.data.frame(bind_rows(part2) |> mutate(across(where(is.numeric), ~ round(.x, 2)))))
print(as.data.frame(bind_rows(part3) |> mutate(across(where(is.numeric), ~ round(.x, 3)))))
print(as.data.frame(bind_rows(part4) |> mutate(across(where(is.numeric), ~ signif(.x, 2)))))
cat("done\n")
