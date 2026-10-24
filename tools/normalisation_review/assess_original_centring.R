# Is the original centring defensible? Read-only checks added on 2026-10-22 to the normalisation review.
#
#   Q  acquisition QC and raw (uncentred, DIA-NN-normalised) levels by group, at animal level: plate-adjusted, within
#      plate, and within the two balanced staining/LCM processing triplets (A764/A135/A127, A111/A139/A129; from the
#      proteome-state audit of 2026-10-17)
#   D  detection by abundance decile and raw intensity quantiles by group: is the depth difference a loss of
#      low-abundance values only, with high-abundance proteins unchanged?
#   U  per region-layer unit: SUS-RES difference in missing fraction against the SUS-RES offset difference each
#      centring imposes, and the median SUS-RES logFC over all proteins (archived original and current animal-level
#      inputs)
#   A  does the shared-protein (anchor) reference contain group-regulated proteins? Direction balance of nominal anchor
#      effects under both centrings, and offsets recomputed on anchors without any nominal effect under either centring
#
# Writes results/reviewer_audit/normalisation_choice_2026-10-21/original_assessment/. Changes nothing else.
# USAGE  Rscript tools/normalisation_review/assess_original_centring.R      (working directory: repo root)

suppressPackageStartupMessages({ library(readr); library(readxl); library(dplyr); library(tidyr); library(limma); library(matrixStats) })
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
REV <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21")
OUT <- file.path(REV, "original_assessment"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(x, name) write_csv(x, file.path(OUT, name), na = "")
state <- readRDS(file.path(REV, "validation", "validation_state.rds"))
OFF <- state$offsets |> mutate(missing_fraction = unname(missing_fraction), n_detected = unname(n_detected))
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
TRIPLETS <- tibble(animal = c("A764", "A135", "A127", "A111", "A139", "A129"), triplet = rep(c("T3", "T4"), each = 3))
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
gct_path <- function(root, ds) file.path(root, "data", "processed", "01_preprocessing", "protigy_input_animal_level", ds, paste0(ds, "_animal_level.gct"))

# ------------------------------------------------------------------------------------------------ Q. QC by group
qc_vars <- c(proteins_identified = "Proteins.Identified", precursors_identified = "Precursors.Identified", missing_fraction = "missing_fraction",
             log10_MS1_signal = "MS1.Signal", log10_MS2_signal = "MS2.Signal", log10_total_quantity = "Total.Quantity",
             normalisation_instability = "Normalisation.Instability", raw_all_observed_median = "raw_original_all_observed",
             raw_anchor_median = "raw_new_complete_case")
qc <- list()
for (ds in DATASETS) {
  o <- OFF |> filter(dataset == ds)
  for (nm in names(qc_vars)) {
    v <- o[[qc_vars[[nm]]]]; if (grepl("^log10", nm)) v <- log10(v)
    # remove the region-layer unit, then average per animal (hemispheres and units)
    a <- o |> mutate(v = v - ave(v, unit)) |> group_by(animal, group, plate) |> summarise(v = mean(v), .groups = "drop") |>
      left_join(TRIPLETS, by = "animal")
    gm <- tapply(a$v, a$group, mean)
    fit <- lm(v ~ group + plate, data = a |> mutate(group = relevel(factor(group), "RES")))
    cf <- summary(fit)$coefficients
    wp <- a |> group_by(plate) |> summarise(d = mean(v[group == "SUS"]) - mean(v[group == "RES"]), .groups = "drop")
    tr <- a |> filter(!is.na(triplet)) |> group_by(triplet) |> summarise(d = v[group == "SUS"] - v[group == "RES"], .groups = "drop")
    qc[[length(qc) + 1L]] <- tibble(dataset = ds, metric = nm, CON = gm[["CON"]], RES = gm[["RES"]], SUS = gm[["SUS"]],
      SUS_minus_RES_plate_adjusted = cf["groupSUS", "Estimate"], p_SUS_vs_RES_plate_adjusted = cf["groupSUS", "Pr(>|t|)"],
      SUS_minus_RES_within_plate_B = wp$d[wp$plate == "B"], SUS_minus_RES_within_plate_C = wp$d[wp$plate == "C"],
      SUS_minus_RES_within_triplet_T3 = tr$d[tr$triplet == "T3"], SUS_minus_RES_within_triplet_T4 = tr$d[tr$triplet == "T4"],
      between_hemisphere_sd = sd(o |> mutate(v = v) |> group_by(animal, unit) |> filter(n() == 2) |> summarise(d = diff(v[order(hemisphere)]) / sqrt(2), .groups = "drop") |> pull(d)))
  }
}
qc <- bind_rows(qc) |> mutate(note = "animal-level values after removing the region-layer mean; raw medians are DIA-NN-normalised log2 before any centring")
w(qc, "acquisition_qc_by_group.csv")

# ------------------------------------------------------------------------------------------------ D. detection by abundance
raw <- read_tsv("data/raw/pg_matrix/quicksearch.pg_matrix.tsv", show_col_types = FALSE, progress = FALSE)
X_all <- log2(as.matrix(raw[, -(1:4)])); rownames(X_all) <- raw$Protein.Group
dec <- list(); qnt <- list()
for (ds in DATASETS) {
  o <- OFF |> filter(dataset == ds)
  X <- X_all[, o$sample_id]; X <- X[rowSums(!is.na(X)) > 0, , drop = FALSE]
  mu <- rowMeans(X, na.rm = TRUE); d10 <- cut(rank(mu, ties.method = "first"), 10, labels = FALSE)
  miss <- sapply(1:10, function(k) colMeans(is.na(X[d10 == k, , drop = FALSE])))   # samples x deciles
  a <- as_tibble(miss, .name_repair = ~ paste0("d", 1:10)) |> mutate(animal = o$animal, group = o$group, unit = o$unit) |>
    pivot_longer(starts_with("d"), names_to = "decile", values_to = "missing") |> group_by(decile, unit) |> mutate(missing_c = missing - mean(missing)) |>
    group_by(decile, animal, group) |> summarise(missing = mean(missing), missing_c = mean(missing_c), .groups = "drop")
  dec[[ds]] <- a |> group_by(decile) |> summarise(dataset = ds, missing_CON = mean(missing[group == "CON"]), missing_RES = mean(missing[group == "RES"]),
                                                   missing_SUS = mean(missing[group == "SUS"]),
                                                   SUS_minus_RES_unit_adjusted = mean(missing_c[group == "SUS"]) - mean(missing_c[group == "RES"]), .groups = "drop") |>
    mutate(decile = as.integer(sub("d", "", decile))) |> arrange(decile) |>
    mutate(share_of_SUS_excess_missing = pmax(SUS_minus_RES_unit_adjusted, 0) / sum(pmax(SUS_minus_RES_unit_adjusted, 0)))
  # raw quantiles of the observed values and the median of the top half of proteins (by overall abundance)
  top <- d10 >= 6
  s <- tibble(animal = o$animal, group = o$group, unit = o$unit, plate = o$plate,
              q25 = colQuantiles(X, probs = 0.25, na.rm = TRUE), q50 = colMedians(X, na.rm = TRUE), q75 = colQuantiles(X, probs = 0.75, na.rm = TRUE),
              q90 = colQuantiles(X, probs = 0.90, na.rm = TRUE), top_half_median = colMedians(X[top, ], na.rm = TRUE),
              anchor_median = colMedians(X[rowSums(is.na(X)) == 0, ]))
  qnt[[ds]] <- s |> pivot_longer(c(q25, q50, q75, q90, top_half_median, anchor_median), names_to = "statistic") |>
    group_by(statistic, unit) |> mutate(value = value - mean(value)) |> group_by(statistic, animal, group, plate) |> summarise(value = mean(value), .groups = "drop") |>
    group_by(statistic) |> summarise(dataset = ds, SUS_minus_RES = mean(value[group == "SUS"]) - mean(value[group == "RES"]),
                                     SUS_minus_RES_within_plate_B = mean(value[group == "SUS" & plate == "B"]) - mean(value[group == "RES" & plate == "B"]),
                                     SUS_minus_RES_within_plate_C = mean(value[group == "SUS" & plate == "C"]) - mean(value[group == "RES" & plate == "C"]), .groups = "drop")
}
w(bind_rows(dec), "missing_by_abundance_decile.csv")
w(bind_rows(qnt) |> mutate(note = "raw DIA-NN log2 (before centring), unit-adjusted, animal level; q = quantiles of the observed values; top_half_median = median of the proteins in the upper five abundance deciles"),
  "raw_intensity_quantiles_by_group.csv")

# ------------------------------------------------------------------------------------------------ U. per-unit depth artefact
gs <- list(); unit_tab <- list(); anchors_res <- list(); inv_off <- list()
fit_da <- function(m, cell) {
  f <- factor(cell, levels = unique(cell)); design <- model.matrix(~0 + f); colnames(design) <- levels(f)
  units <- unique(sub("_[0-9]+$", "", cell)); pairs <- list(c("2", "1"), c("3", "2"), c("3", "1"))
  nm <- unlist(lapply(pairs, function(p) paste0(units, "_", p[1], "_over_", units, "_", p[2])))
  ex <- unlist(lapply(pairs, function(p) paste0(units, "_", p[1], "-", units, "_", p[2])))
  cm <- makeContrasts(contrasts = ex, levels = design); colnames(cm) <- nm
  eb <- eBayes(contrasts.fit(lmFit(m, design), cm), robust = TRUE, trend = FALSE)
  lapply(setNames(nm, nm), function(cc) topTable(eb, coef = cc, number = Inf, sort.by = "none", adjust.method = "BH"))
}
GRP <- c("1" = "CON", "2" = "RES", "3" = "SUS")
for (ds in DATASETS) {
  o <- OFF |> filter(dataset == ds)
  g_new <- read_gct(gct_path(".", ds)); g_old <- read_gct(gct_path(ARCHIVE, ds))
  rid <- intersect(rownames(g_new$m), rownames(g_old$m))
  cell <- g_new$coldesc$region_layer_ExpGroup; units <- unique(sub("_[0-9]+$", "", cell))
  tabs <- list(original = fit_da(g_old$m[rid, ], cell), new = fit_da(g_new$m[rid, ], cell))
  for (u in units) {
    ou <- o |> filter(sub("_microglia$", "", unit) == u | unit == u) |> group_by(animal, group) |>
      summarise(missing = mean(missing_fraction), off_o = mean(rel_original_all_observed), off_n = mean(rel_new_complete_case),
                off_m = mean(rel_median_ratio_complete_case), .groups = "drop")
    d <- function(v) mean(v[ou$group == "SUS"]) - mean(v[ou$group == "RES"])
    cc <- paste0(u, "_3_over_", u, "_2")
    unit_tab[[length(unit_tab) + 1L]] <- tibble(dataset = ds, unit = u, SUS_minus_RES_missing_fraction = d(ou$missing),
      SUS_minus_RES_offset_original = d(ou$off_o), SUS_minus_RES_offset_new = d(ou$off_n), SUS_minus_RES_offset_median_ratio = d(ou$off_m),
      median_SUS_RES_logFC_original = median(tabs$original[[cc]]$logFC), median_SUS_RES_logFC_new = median(tabs$new[[cc]]$logFC),
      share_negative_logFC_original = mean(tabs$original[[cc]]$logFC < 0), share_negative_logFC_new = mean(tabs$new[[cc]]$logFC < 0),
      n_fdr05_original = sum(tabs$original[[cc]]$adj.P.Val < 0.05), n_fdr05_new = sum(tabs$new[[cc]]$adj.P.Val < 0.05),
      share_down_among_fdr05_original = { t <- tabs$original[[cc]]; if (sum(t$adj.P.Val < 0.05)) mean(t$logFC[t$adj.P.Val < 0.05] < 0) else NA_real_ },
      share_down_among_fdr05_new = { t <- tabs$new[[cc]]; if (sum(t$adj.P.Val < 0.05)) mean(t$logFC[t$adj.P.Val < 0.05] < 0) else NA_real_ })
  }
  # ---- A. anchor reference: direction balance of nominal effects (all units and contrasts), under both centrings
  imp <- read_excel(Sys.glob(sprintf("data/processed/01_preprocessing/impute/*_pgmatrix_imputed_%s_*_missing70pct.xlsx", ds)), progress = FALSE,
                    col_types = "text", .name_repair = "minimal")
  pg2name <- setNames(imp$`T: Protein.Names`, imp$Protein.Group)
  anc_names <- intersect(unname(pg2name[state$anchor_sets[[ds]]]), rid)
  for (v in names(tabs)) for (cc in names(tabs[[v]])) {
    t <- tabs[[v]][[cc]][anc_names, ]; lab <- strsplit(sub("^(.*)_([0-9])_over_(.*)_([0-9])$", "\\1|\\2|\\4", cc), "|", fixed = TRUE)[[1]]
    anchors_res[[length(anchors_res) + 1L]] <- tibble(dataset = ds, version = v, unit = lab[1], comparison = paste(GRP[lab[2]], "vs", GRP[lab[3]]),
      n_anchors = nrow(t), nominal_up = mean(t$P.Value < 0.05 & t$logFC > 0), nominal_down = mean(t$P.Value < 0.05 & t$logFC < 0),
      median_logFC = median(t$logFC), n_fdr05 = sum(t$adj.P.Val < 0.05))
  }
  # anchors without any nominal effect (p < 0.05) in any unit or contrast under either centring
  nominal_any <- Reduce(`|`, lapply(unlist(tabs, recursive = FALSE), function(t) t[anc_names, "P.Value"] < 0.05))
  inv <- anc_names[!nominal_any]
  inv_pg <- names(pg2name)[match(inv, pg2name)]
  X <- X_all[, o$sample_id]
  oi <- colMedians(X[inv_pg, ], na.rm = TRUE); oi <- oi - mean(oi)
  gs[[ds]] <- tibble(dataset = ds, n_anchors = length(anc_names), n_invariant = length(inv),
                     cor_invariant_vs_new = cor(oi, o$rel_new_complete_case), cor_invariant_vs_original = cor(oi, o$rel_original_all_observed),
                     max_abs_invariant_minus_new = max(abs(oi - o$rel_new_complete_case)))
  inv_off[[ds]] <- o |> mutate(rel_invariant = oi) |> group_by(animal, group, unit) |>
    summarise(across(c(rel_invariant, rel_new_complete_case, rel_original_all_observed), mean), .groups = "drop") |>
    group_by(unit) |> summarise(dataset = ds,
      SUS_minus_RES_invariant = mean(rel_invariant[group == "SUS"]) - mean(rel_invariant[group == "RES"]),
      SUS_minus_RES_new = mean(rel_new_complete_case[group == "SUS"]) - mean(rel_new_complete_case[group == "RES"]),
      SUS_minus_RES_original = mean(rel_original_all_observed[group == "SUS"]) - mean(rel_original_all_observed[group == "RES"]), .groups = "drop")
}
ut <- bind_rows(unit_tab)
w(ut, "depth_artefact_by_unit.csv")
w(tibble(statistic = c("cor(SUS-RES missing, SUS-RES original offset) across units", "cor(SUS-RES missing, SUS-RES new offset) across units",
                       "cor(SUS-RES original offset, median SUS-RES logFC original)", "cor(SUS-RES missing, median SUS-RES logFC original)",
                       "cor(SUS-RES missing, median SUS-RES logFC new)"),
         r = c(cor(ut$SUS_minus_RES_missing_fraction, ut$SUS_minus_RES_offset_original), cor(ut$SUS_minus_RES_missing_fraction, ut$SUS_minus_RES_offset_new),
               cor(ut$SUS_minus_RES_offset_original, ut$median_SUS_RES_logFC_original), cor(ut$SUS_minus_RES_missing_fraction, ut$median_SUS_RES_logFC_original),
               cor(ut$SUS_minus_RES_missing_fraction, ut$median_SUS_RES_logFC_new)), n_units = nrow(ut)),
  "depth_artefact_correlations.csv")
w(bind_rows(anchors_res), "anchor_effect_direction_balance.csv")
w(bind_rows(gs), "invariant_anchor_offsets.csv")
w(bind_rows(inv_off), "invariant_anchor_group_offsets_by_unit.csv")
cat("done:", OUT, "\n")
