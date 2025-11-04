# Layer enrichment (controls) + explicit within-layer contrasts + microglia analyses
# FDR alpha = 0.05
# - Species filter: keep only _MOUSE proteins
# - Duplicates: keep first per (protein_id, sample_id)
# - Gene symbols: keep first symbol per protein_id (deterministic)
# - Volcano dots: nonsig grey, sig blue + red overlay
# - Top lists: save top up/down; microglia volcano labels top 5 up and top 5 down

suppressPackageStartupMessages({
  if (!require("pacman", quietly = TRUE)) install.packages("pacman")
  pacman::p_load(readxl, dplyr, tidyr, tibble, stringr, purrr,
                 limma, ggplot2, ggrepel, readr, forcats)
})
set.seed(17)

# ------------------------------- Parameters -----------------------------
region_of_interest <- "DG"               # Set here to "CA1", "CA2", "CA3", etc.
layer_levels <- c("mo", "sg", "po") # Set layers for the chosen region here
alpha_layer  <- 0.05
alpha_within <- 0.05
extra_layer  <- "microglia"

# Optional protein_id -> gene_symbol map
idmap_file <- NULL  # CSV with columns: protein_id, gene_symbol (optional)

# ------------------------------- Paths ---------------------------------
base <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data"
male_xlsx  <- file.path(base, "male.data.xlsx")
sinfo_xlsx <- file.path(base, "sample_info.xlsx")

out_base <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results",
                      tolower(region_of_interest), "layer_and_within_group_fdr005")
dirs <- list(
  figs_layer    = file.path(out_base, "figures", "layer_enrichment"),
  figs_within   = file.path(out_base, "figures", "within_layer"),
  figs_micro    = file.path(out_base, "figures", "microglia_layer"),
  tabs_layer    = file.path(out_base, "tables",  "layer_enrichment"),
  tabs_within   = file.path(out_base, "tables",  "within_layer"),
  tabs_micro    = file.path(out_base, "tables",  "microglia_layer"),
  logs          = file.path(out_base, "logs"),
  rds           = file.path(out_base, "rds"),
  qc            = file.path(out_base, "qc")
)
within_sub <- list(
  con_vs_sus = list(figs = file.path(dirs$figs_within, "con_vs_sus"),
                    tabs = file.path(dirs$tabs_within, "con_vs_sus")),
  con_vs_res = list(figs = file.path(dirs$figs_within, "con_vs_res"),
                    tabs = file.path(dirs$tabs_within, "con_vs_res")),
  sus_vs_res = list(figs = file.path(dirs$figs_within, "sus_vs_res"),
                    tabs = file.path(dirs$tabs_within, "sus_vs_res")),
  res_vs_con = list(figs = file.path(dirs$figs_within, "res_vs_con"),
                    tabs = file.path(dirs$tabs_within, "res_vs_con"))
)
invisible(lapply(c(dirs, within_sub$con_vs_sus, within_sub$con_vs_res, within_sub$sus_vs_res, within_sub$res_vs_con),
                 dir.create, showWarnings = FALSE, recursive = TRUE))
ts <- format(Sys.time(), "%Y%m%d_%H%M")

# ------------------------------ Helpers --------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

to_numeric_matrix_hard <- function(x) {
  stopifnot(is.matrix(x))
  if (is.numeric(x)) return(x)
  nr <- nrow(x); nc <- ncol(x)
  out <- matrix(NA_real_, nrow = nr, ncol = nc, dimnames = dimnames(x))
  for (j in seq_len(nc)) {
    colj <- x[, j]
    if (is.list(colj)) {
      colj <- vapply(colj, function(v) if (length(v)) v[[1]] else NA_character_, character(1))
    }
    colj[colj %in% c("NA","", NA)] <- NA
    out[, j] <- suppressWarnings(as.numeric(colj))
  }
  out
}
row_median_base <- function(m) {
  apply(m, 1, function(v) { v <- v[is.finite(v)]; if (!length(v)) NA_real_ else stats::median(v) })
}
row_mad_base <- function(m, constant = 1) {
  apply(m, 1, function(v) { v <- v[is.finite(v)]; if (!length(v)) return(NA_real_); stats::mad(v, constant = constant, na.rm = TRUE) })
}
impute_layer_robust_base <- function(m, groups,
                                    center_fun = row_median_base,
                                    scale_fun  = function(x) row_mad_base(x, constant = 1)) {
  stopifnot(is.matrix(m))
  out <- if (is.numeric(m)) m else to_numeric_matrix_hard(m)
  layers <- droplevels(groups)
  for (lvl in levels(layers)) {
    idx <- which(layers == lvl); if (!length(idx)) next
    sub <- out[, idx, drop = FALSE]
    sub <- if (is.numeric(sub)) sub else to_numeric_matrix_hard(sub)
    c_row <- center_fun(sub)
    s_row <- scale_fun(sub); s_row[!is.finite(s_row) | s_row == 0] <- 1
    z <- sweep(sub, 1, c_row, FUN = "-"); z <- sweep(z, 1, s_row, FUN = "/")
    z_meds <- row_median_base(z)
    miss <- is.na(z)
    if (any(miss)) {
      coords <- which(miss, arr.ind = TRUE)
      z[cbind(coords[,1], coords[,2])] <- z_meds[coords[,1]]
    }
    sub_imputed <- sweep(z, 1, s_row, FUN = "*"); sub_imputed <- sweep(sub_imputed, 1, c_row, FUN = "+")
    out[, idx] <- sub_imputed
  }
  out
}
impute_layer_mnar_finalize <- function(m, groups, downshift = 1.8) {
  stopifnot(is.matrix(m))
  out <- m
  layers <- droplevels(groups)
  for (lvl in levels(layers)) {
    idx <- which(layers == lvl); if (!length(idx)) next
    sub <- out[, idx, drop = FALSE]
    obs_vals <- as.numeric(sub[is.finite(sub)])
    if (length(obs_vals) < 2) {
      all_obs <- as.numeric(out[is.finite(out)])
      small <- stats::median(all_obs, na.rm = TRUE) - downshift * stats::sd(all_obs, na.rm = TRUE)
    } else {
      small <- mean(obs_vals, na.rm = TRUE) - downshift * stats::sd(obs_vals, na.rm = TRUE)
    }
    miss <- is.na(sub)
    if (any(miss)) sub[miss] <- small
    out[, idx] <- sub
  }
  out
}
keep_idx <- function(m, groups, min_prop = 0.5) {
  by_grp <- split(seq_len(ncol(m)), groups)
  apply(m, 1, function(x) all(vapply(by_grp, function(idx) mean(!is.na(x[idx])) >= min_prop, logical(1))))
}
run_within_layer_contrast <- function(expr_mat_imp, samp_meta, group_a, group_b) {
  layers_present <- levels(droplevels(samp_meta$layer))
  lapply(layers_present, function(lyr) {
    idx <- which(samp_meta$layer == lyr)
    cols <- samp_meta$sample_id[idx]
    grp  <- droplevels(samp_meta$ExpGroup[idx])
    if (!all(c(group_a, group_b) %in% levels(grp))) return(NULL)
    if (any(table(grp) < 2)) return(NULL)
    sel <- grp %in% c(group_a, group_b)
    cols <- cols[sel]; grp <- droplevels(grp[sel])
    X <- expr_mat_imp[, cols, drop = FALSE]
    grp <- relevel(grp, ref = group_b)  # A - B
    design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
    if (!(group_a %in% colnames(design))) return(NULL)
    contr <- limma::makeContrasts(contrasts = paste0(group_a, " - ", group_b), levels = design)
    fit  <- limma::lmFit(X, design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, contr), trend = TRUE, robust = TRUE)
    tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
    tt$protein_id <- rownames(tt); tt$layer <- lyr; tt$group_contrast <- paste0(group_a, " vs ", group_b)
    tt
  }) %>% Filter(Negate(is.null), .)
}
standardize_cols <- function(df) {
  out <- df
  if ("P.Value" %in% names(out)) out <- out %>% rename(pval = P.Value)
  if ("adj.P.Val" %in% names(out)) out <- out %>% rename(padj = adj.P.Val)
  out
}
normalize_idmap_first <- function(idmap_raw) {
  idmap_raw %>%
    distinct(protein_id, gene_symbol) %>%
    filter(!is.na(protein_id), protein_id != "", !is.na(gene_symbol), gene_symbol != "") %>%
    arrange(protein_id, gene_symbol) %>%  # adjust ordering rule if needed
    group_by(protein_id) %>%
    slice_head(n = 1) %>%
    ungroup()
}

save_result_tables <- function(df, out_dir, stem, alpha = 0.05, idmap_norm = NULL) {
  if (!nrow(df)) return(invisible(NULL))
  df2 <- standardize_cols(df)
  if (!is.null(idmap_norm)) df2 <- df2 %>% left_join(idmap_norm, by = "protein_id")
  readr::write_csv(df2, file.path(out_dir, paste0(stem, "_full_", ts, ".csv")))
  sig <- df2 %>% filter(padj <= alpha) %>% arrange(if_any(all_of(intersect("layer", names(df2))), ~.x), padj, desc(abs(logFC)))
  readr::write_csv(sig, file.path(out_dir, paste0(stem, "_significant_", ts, ".csv")))
  if ("layer" %in% names(df2)) {
    top5 <- df2 %>% group_by(layer) %>% arrange(padj, desc(abs(logFC)), .by_group = TRUE) %>% slice_head(n = 5) %>% ungroup()
  } else {
    top5 <- df2 %>% arrange(padj, desc(abs(logFC))) %>% slice_head(n = 5)
  }
  readr::write_csv(top5, file.path(out_dir, paste0(stem, "_top5_", ts, ".csv")))
  invisible(list(full = df2, sig = sig, top5 = top5))
}

save_top_updown_tables <- function(df, out_dir, stem, alpha = 0.05, idmap_norm = NULL, by_layer = TRUE, n = 50) {
  if (!nrow(df)) return(invisible(NULL))
  df2 <- standardize_cols(df)
  if (!is.null(idmap_norm)) df2 <- df2 %>% left_join(idmap_norm, by = "protein_id")
  sig <- df2 %>% filter(padj <= alpha)
  if (!nrow(sig)) return(invisible(NULL))
  if (by_layer && "layer" %in% names(sig)) {
    top_up <- sig %>% filter(logFC > 0) %>% group_by(layer) %>% arrange(desc(logFC), padj, .by_group = TRUE) %>% slice_head(n = n) %>% ungroup()
    top_dn <- sig %>% filter(logFC < 0) %>% group_by(layer) %>% arrange(logFC, padj, .by_group = TRUE) %>% slice_head(n = n) %>% ungroup()
  } else {
    top_up <- sig %>% filter(logFC > 0) %>% arrange(desc(logFC), padj) %>% slice_head(n = n)
    top_dn <- sig %>% filter(logFC < 0) %>% arrange(logFC, padj) %>% slice_head(n = n)
  }
  readr::write_csv(top_up, file.path(out_dir, paste0(stem, "_top_up_", n, "_", ts, ".csv")))
  readr::write_csv(top_dn, file.path(out_dir, paste0(stem, "_top_down_", n, "_", ts, ".csv")))
  invisible(list(top_up = top_up, top_down = top_dn))
}

build_volcano <- function(df, xlab, out_svg, out_png, alpha = 0.05, labels_df = NULL) {
  if (!nrow(df)) return(invisible(NULL))
  plot_df <- standardize_cols(df) %>%
    mutate(minus_log10_padj = -log10(padj),
           sig = padj <= alpha,
           sig_flag = ifelse(sig, "sig", "nonsig"))
  ymax <- max(plot_df$minus_log10_padj[is.finite(plot_df$minus_log10_padj)], na.rm = TRUE) * 1.05

  p <- ggplot(plot_df, aes(x = logFC, y = minus_log10_padj)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_point(aes(color = sig_flag), alpha = 0.55, size = 1.8, na.rm = TRUE) +
    geom_point(data = ~ subset(., sig), color = "red", fill = "red", alpha = 0.95, size = 2.1, na.rm = TRUE) +
    facet_wrap(~ layer, ncol = 1, scales = "fixed") +
    scale_y_continuous(limits = c(0, ymax), name = "-log10 adj. p (within-layer)") +
    scale_x_continuous(name = xlab) +
    scale_color_manual(values = c(nonsig = "grey80", sig = "#2C6BED"), guide = "none") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 13),
          strip.text = element_text(face = "bold"))

  if (!is.null(labels_df) && nrow(labels_df)) {
    p <- p + ggrepel::geom_text_repel(
      data = labels_df,
      aes(x = logFC, y = -log10(padj), label = label),
      color = "black", size = 3, max.overlaps = 200,
      box.padding = 0.25, point.padding = 0.25, segment.color = "grey50", segment.size = 0.3
    )
  }

  ggsave(out_svg, p, width = 4, height = 12, device = "svg")
  ggsave(out_png, p, width = 4, height = 12, dpi = 300)
  invisible(p)
}

# --------------------------- Load and prepare ---------------------------
male <- readxl::read_excel(male_xlsx, col_names = TRUE)
if (is.na(names(male)[1]) || names(male)[1] == "") names(male)[1] <- "protein_id"
male <- male %>% rename(protein_id = 1) %>% mutate(protein_id = as.character(protein_id))

# NEW: species filter — keep only entries with "_MOUSE" in protein_id
male <- male %>%
  filter(!is.na(protein_id), protein_id != "") %>%
  filter(grepl("_MOUSE\\b", protein_id, perl = TRUE))

male_sample_cols <- grep("^s\\d+$", names(male), value = TRUE); stopifnot(length(male_sample_cols) > 0)

sinfo <- readxl::read_excel(sinfo_xlsx, col_names = TRUE)
names(sinfo)[1] <- "sample_id"
sinfo <- sinfo %>% mutate(sample_id = as.character(sample_id))

present_samples <- intersect(male_sample_cols, sinfo$sample_id); stopifnot(length(present_samples) > 0)

# Keep-first deduplication before pivoting
expr_long <- male %>%
  pivot_longer(all_of(present_samples), names_to = "sample_id", values_to = "log2_intensity") %>%
  mutate(log2_intensity = suppressWarnings(as.numeric(log2_intensity))) %>%
  arrange(protein_id, sample_id) %>%
  distinct(protein_id, sample_id, .keep_all = TRUE)

# Optional: load normalized ID map (first symbol only)
idmap_norm <- NULL
if (!is.null(idmap_file) && file.exists(idmap_file)) {
  idmap_raw <- readr::read_csv(idmap_file, show_col_types = FALSE)
  stopifnot(all(c("protein_id","gene_symbol") %in% names(idmap_raw)))
  idmap_norm <- normalize_idmap_first(idmap_raw)
}

samp <- sinfo %>%
  filter(sample_id %in% present_samples) %>%
  mutate(layer = tolower(trimws(layer)),
         region = trimws(region),
         ExpGroup = as.character(ExpGroup))

# -------------------- 1) Layer enrichment (controls) -------------------
samp_con <- samp %>% filter(ExpGroup == "con")
expr_long_con <- expr_long %>% semi_join(tibble(sample_id = samp_con$sample_id), by = "sample_id")

samp_region_con <- samp_con %>%
  filter(region == region_of_interest, layer %in% layer_levels) %>%
  mutate(layer = factor(layer, levels = layer_levels))
expr_region_con <- expr_long_con %>% semi_join(samp_region_con, by = "sample_id")

expr_wide <- expr_region_con %>% select(protein_id, sample_id, log2_intensity) %>%
  pivot_wider(names_from = sample_id, values_from = log2_intensity) %>%
  arrange(protein_id) %>% filter(!is.na(protein_id), protein_id != "")

expr_mat <- expr_wide %>% select(-protein_id) %>% as.matrix()
rownames(expr_mat) <- expr_wide$protein_id
expr_mat <- to_numeric_matrix_hard(expr_mat)

groups <- samp_region_con$layer[match(colnames(expr_mat), samp_region_con$sample_id)]
keep <- keep_idx(expr_mat, groups, 0.5)
expr_mat <- expr_mat[keep, , drop = FALSE]

qc_counts <- samp_region_con %>% count(layer, name = "n_samples")
qc_missing <- data.frame(protein_id = rownames(expr_mat), frac_missing = rowMeans(is.na(expr_mat)))
readr::write_csv(qc_counts,  file.path(dirs$qc, paste0(region_of_interest, "_sample_counts_controls_", ts, ".csv")))
readr::write_csv(qc_missing, file.path(dirs$qc, paste0(region_of_interest, "_protein_missingness_controls_", ts, ".csv")))

expr_stage1 <- impute_layer_robust_base(expr_mat, groups)
expr_mat_imp <- impute_layer_mnar_finalize(expr_stage1, groups, downshift = 1.8)

design <- model.matrix(~ 0 + factor(groups, levels = layer_levels))
colnames(design) <- toupper(layer_levels)
other_avg <- function(level, levels) {
  others <- setdiff(levels, level)
  sprintf("%s - (%s)/%d", toupper(level), paste(toupper(others), collapse = " + "), length(others))
}
# Create contrast names for all layers
contrast_strings <- setNames(
  lapply(layer_levels, function(level) {
    others <- setdiff(layer_levels, level)
    sprintf("%s - (%s)/%d", toupper(level), paste(toupper(others), collapse = " + "), length(others))
  }),
  toupper(layer_levels)
)

# Build the contrast matrix dynamically based on existing layers
contrast_matrix <- do.call(
  limma::makeContrasts,
  c(contrast_strings, list(levels = design))
)

fit  <- limma::lmFit(expr_mat_imp, design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix), trend = TRUE, robust = TRUE)

# Get uppercase layer names dynamically
layer_coefs <- toupper(layer_levels)

tt_list <- lapply(layer_coefs, function(coef_nm) {
  tt <- limma::topTable(fit2, coef = coef_nm, number = Inf, sort.by = "none")
  tt$protein_id <- rownames(tt); tt$layer <- tolower(coef_nm); tt
})

stats_df <- dplyr::bind_rows(tt_list) %>%
  left_join(
    data.frame(protein_id = rownames(expr_mat_imp), x_mean = rowMeans(expr_mat_imp, na.rm = TRUE), check.names = FALSE),
    by = "protein_id"
  ) %>% mutate(layer = factor(layer, levels = layer_levels)) %>%
  standardize_cols()

# Labels for layer enrichment: top 5 positive per layer
labels_df <- stats_df %>%
  { if (!is.null(idmap_norm)) left_join(., idmap_norm, by = "protein_id") else mutate(., gene_symbol = NA_character_) } %>%
  group_by(layer) %>%
  group_modify(~{
    .x %>%
      filter(logFC > 0, padj <= alpha_layer) %>%
      arrange(padj, desc(abs(logFC))) %>%
      mutate(label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, protein_id)) %>%
      slice_head(n = 5)
  }) %>% ungroup()

# Save layer enrichment and top up/down
save_result_tables(stats_df, dirs$tabs_layer, stem = "layer_enrichment_controls", alpha = alpha_layer, idmap_norm = idmap_norm)
save_top_updown_tables(stats_df, dirs$tabs_layer, stem = "layer_enrichment_controls", alpha = alpha_layer, idmap_norm = idmap_norm, by_layer = TRUE, n = 50)

# Plot layer enrichment
min_p <- min(stats_df$padj[is.finite(stats_df$padj)], na.rm = TRUE)
padj_breaks <- c(min_p, alpha_layer)
padj_labels <- c(sprintf("0.0001"), sprintf("%.3f", alpha_layer))
stats_plot <- stats_df %>% mutate(logFC_pos = pmax(logFC, 0), padj_color = pmin(padj, alpha_layer))
ymax <- max(stats_plot$logFC_pos[is.finite(stats_plot$logFC_pos)], na.rm = TRUE) * 1.05

p_layer <- ggplot2::ggplot(stats_plot, aes(x = x_mean, y = logFC_pos, color = padj_color)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_point(alpha = 0.8, size = 1.4, na.rm = TRUE) +
  facet_wrap(~ layer, ncol = 1, scales = "fixed") +
  scale_y_continuous(limits = c(0, ymax), name = "log2 fold enrichment in layer") +
  scale_x_continuous(name = "log2 mean intensity across layers") +
  scale_color_gradient(low = "#5858ac", high = "#d3d3d4", limits = c(min_p, alpha_layer), breaks = padj_breaks, labels = padj_labels) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 13), strip.text = element_text(face = "bold")) +
  ggrepel::geom_text_repel(
    data = labels_df %>% mutate(logFC_pos = pmax(logFC, 0)),
    aes(x = x_mean, y = logFC_pos, label = label),
    color = "black", size = 3, max.overlaps = 200,
    box.padding = 0.25, point.padding = 0.25, segment.color = "grey50", segment.size = 0.3
  ) + labs(color = "adj. p")

ggsave(file.path(dirs$figs_layer, paste0(region_of_interest, "_layers_enrichment_", ts, ".svg")), p_layer, width = 4, height = 12, device = "svg")
ggsave(file.path(dirs$figs_layer, paste0(region_of_interest, "_layers_enrichment_", ts, ".png")), p_layer, width = 4, height = 12, dpi = 300)

# -------------------- 1b) Microglia summary (controls only) -------------------
samp_micro_con <- samp_con %>%
  filter(tolower(layer) == tolower(extra_layer) | grepl("microglia", tolower(samp_con$celltype %||% ""), fixed = TRUE))
if (nrow(samp_micro_con)) {
  expr_micro_con <- expr_long_con %>% semi_join(samp_micro_con, by = "sample_id")
  expr_wide_micro <- expr_micro_con %>% select(protein_id, sample_id, log2_intensity) %>%
    pivot_wider(names_from = sample_id, values_from = log2_intensity) %>%
    arrange(protein_id) %>% filter(!is.na(protein_id), protein_id != "")
  expr_micro_mat <- expr_wide_micro %>% select(-protein_id) %>% as.matrix()
  rownames(expr_micro_mat) <- expr_wide_micro$protein_id
  expr_micro_mat <- to_numeric_matrix_hard(expr_micro_mat)

  if (ncol(expr_micro_mat) >= 2) {
    micro_mean <- rowMeans(expr_micro_mat, na.rm = TRUE)
    micro_missing <- rowMeans(is.na(expr_micro_mat))
    micro_tab <- tibble(
      protein_id = rownames(expr_micro_mat),
      log2_mean = as.numeric(micro_mean),
      frac_missing = as.numeric(micro_missing)
    ) %>% arrange(desc(log2_mean)) %>%
      { if (!is.null(idmap_norm)) left_join(., idmap_norm, by = "protein_id") else . }
    readr::write_csv(micro_tab, file.path(dirs$tabs_micro, paste0("microglia_controls_summary_", ts, ".csv")))
    p_micro <- ggplot(micro_tab, aes(x = rank(-log2_mean), y = log2_mean)) +
      geom_point(size = 0.8, alpha = 0.8) +
      labs(x = "rank (descending intensity)", y = "log2 mean intensity", title = "Microglia (controls)") +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank(), axis.text = element_text(size = 12))
    ggsave(file.path(dirs$figs_micro, paste0("microglia_controls_rank_", ts, ".svg")), p_micro, width = 6, height = 4, device = "svg")
    ggsave(file.path(dirs$figs_micro, paste0("microglia_controls_rank_", ts, ".png")), p_micro, width = 6, height = 4, dpi = 300)
  }
}


# 1c) Microglia vs region average EXCLUDING SP or SG (controls only)
samp_micro_con <- samp_con %>%
  filter(tolower(layer) == tolower(extra_layer) | grepl("microglia", tolower(celltype %||% ""), fixed = TRUE))

# Set the layer to exclude for noSP group dynamically:
exclude_layer <- if (region_of_interest == "DG") "sg" else "sp"

# Define layers for noSP as all layers except the excluded one:
layer_levels_no_sp <- setdiff(layer_levels, exclude_layer)

samp_region_noSP_con <- samp_con %>%
  filter(region == region_of_interest, layer %in% layer_levels_no_sp) %>%
  mutate(layer = factor(layer, levels = layer_levels_no_sp))

if (nrow(samp_micro_con) && nrow(samp_region_noSP_con)) {
  samp_region_or_micro_con <- bind_rows(
    samp_region_noSP_con %>% mutate(layer2 = paste0(tolower(region_of_interest), "_noSP")),
    samp_micro_con      %>% mutate(layer2 = "microglia")
  ) %>% select(sample_id, layer2) %>% mutate(layer2 = factor(layer2, levels = c(paste0(tolower(region_of_interest), "_noSP"), "microglia")))

  expr_con_all <- expr_long_con %>% semi_join(samp_region_or_micro_con, by = "sample_id")
  expr_wide_mg <- expr_con_all %>% select(protein_id, sample_id, log2_intensity) %>%
    pivot_wider(names_from = sample_id, values_from = log2_intensity) %>%
    arrange(protein_id) %>% filter(!is.na(protein_id), protein_id != "")
  expr_mg_mat <- expr_wide_mg %>% select(-protein_id) %>% as.matrix()
  rownames(expr_mg_mat) <- expr_wide_mg$protein_id
  expr_mg_mat <- to_numeric_matrix_hard(expr_mg_mat)

  groups_mg <- samp_region_or_micro_con$layer2[match(colnames(expr_mg_mat), samp_region_or_micro_con$sample_id)]
  keep_mg <- keep_idx(expr_mg_mat, groups_mg, 0.5)
  expr_mg_mat <- expr_mg_mat[keep_mg, , drop = FALSE]

  expr_mg_stage1 <- impute_layer_robust_base(expr_mg_mat, groups_mg)
  expr_mg_imp    <- impute_layer_mnar_finalize(expr_mg_stage1, groups_mg, downshift = 1.8)

  design_mg <- model.matrix(~ 0 + factor(groups_mg, levels = c(paste0(tolower(region_of_interest), "_noSP"), "microglia")))
  colnames(design_mg) <- c(paste0(toupper(region_of_interest), "noSP"), "MICROGLIA")
  contrast_mg <- limma::makeContrasts(MG_vs_regionnoSP = MICROGLIA - get(paste0(toupper(region_of_interest), "noSP")), levels = design_mg)

  fit_mg  <- limma::lmFit(expr_mg_imp, design_mg)
  fit_mg2 <- limma::eBayes(limma::contrasts.fit(fit_mg, contrast_mg), trend = TRUE, robust = TRUE)

  tt_mg <- limma::topTable(fit_mg2, coef = "MG_vs_regionnoSP", number = Inf, sort.by = "none")
  tt_mg$protein_id <- rownames(tt_mg)
  mg_df <- tt_mg %>% standardize_cols() %>% transmute(protein_id, logFC, pval = pval, padj = padj)

  # Top 5 up and top 5 down significant (labels)
  mg_top_sig <- mg_df %>% filter(padj <= alpha_layer)
  mg_top5_up <- mg_top_sig %>% filter(logFC > 0) %>% arrange(padj, desc(logFC)) %>% slice_head(n = 5)
  mg_top5_dn <- mg_top_sig %>% filter(logFC < 0) %>% arrange(padj, logFC)       %>% slice_head(n = 5)
  mg_top10_labels <- bind_rows(
    mg_top5_up %>% mutate(direction = "up"),
    mg_top5_dn %>% mutate(direction = "down")
  ) %>%
    { if (!is.null(idmap_norm)) left_join(., idmap_norm, by = "protein_id") else mutate(., gene_symbol = NA_character_) } %>%
    mutate(label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, protein_id))

  # Save labeled lists and top up/down tables (global)
  readr::write_csv(mg_top5_up, file.path(dirs$tabs_micro, paste0("microglia_vs_", region_of_interest, "avg_noSP_top5_up_", ts, ".csv")))
  readr::write_csv(mg_top5_dn, file.path(dirs$tabs_micro, paste0("microglia_vs_", region_of_interest, "avg_noSP_top5_down_", ts, ".csv")))
  save_top_updown_tables(mg_df, dirs$tabs_micro, stem = paste0("microglia_vs_", region_of_interest, "avg_noSP"), alpha = alpha_layer, idmap_norm = idmap_norm, by_layer = FALSE, n = 50)

  # Volcano with top-5 up and down labels
  plot_mg <- mg_df %>% mutate(minus_log10_padj = -log10(padj), sig = padj <= alpha_layer, sig_flag = ifelse(sig, "sig", "nonsig"))
  ymax_mg <- max(plot_mg$minus_log10_padj[is.finite(plot_mg$minus_log10_padj)], na.rm = TRUE) * 1.05

  p_mg <- ggplot(plot_mg, aes(x = logFC, y = minus_log10_padj)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_point(aes(color = sig_flag), alpha = 0.55, size = 1.8, na.rm = TRUE) +
    geom_point(data = ~ subset(., sig), color = "red", alpha = 0.95, size = 2.1, na.rm = TRUE) +
    ggrepel::geom_text_repel(
      data = mg_top10_labels %>% mutate(minus_log10_padj = -log10(padj)),
      aes(x = logFC, y = minus_log10_padj, label = label),
      color = "black", size = 3, max.overlaps = 200,
      box.padding = 0.25, point.padding = 0.25, segment.color = "grey50", segment.size = 0.3
    ) +
    scale_y_continuous(limits = c(0, ymax_mg), name = "-log10 adj. p") +
    scale_x_continuous(name = paste0("log2 fold change (microglia − ", region_of_interest, "avg(no SP))")) +
    scale_color_manual(values = c(nonsig = "grey80", sig = "#2C6BED"), guide = "none") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 12))

  ggsave(file.path(dirs$figs_micro, paste0("microglia_vs_", region_of_interest, "avg_noSP_", ts, ".svg")), p_mg, width = 5, height = 4, device = "svg")
  ggsave(file.path(dirs$figs_micro, paste0("microglia_vs_", region_of_interest, "avg_noSP_", ts, ".png")), p_mg, width = 5, height = 4, dpi = 300)
}

# ---------------- Microglia within-celltype explicit contrasts ----------------
samp_micro_all <- samp %>%
  filter(tolower(layer) == "microglia" | grepl("microglia", tolower(celltype %||% ""), fixed = TRUE))
if (nrow(samp_micro_all)) {
  expr_micro_all <- expr_long %>% semi_join(samp_micro_all, by = "sample_id")
  expr_wide_micro_all <- expr_micro_all %>%
    select(protein_id, sample_id, log2_intensity) %>%
    pivot_wider(names_from = sample_id, values_from = log2_intensity) %>%
    arrange(protein_id) %>% filter(!is.na(protein_id), protein_id != "")
  expr_micro_all_mat <- expr_wide_micro_all %>% select(-protein_id) %>% as.matrix()
  rownames(expr_micro_all_mat) <- expr_wide_micro_all$protein_id
  expr_micro_all_mat <- to_numeric_matrix_hard(expr_micro_all_mat)

  groups_micro <- factor(rep("microglia", ncol(expr_micro_all_mat)))
  keep_micro <- keep_idx(expr_micro_all_mat, groups_micro, 0.5)
  expr_micro_all_mat <- expr_micro_all_mat[keep_micro, , drop = FALSE]

  expr_micro_stage1 <- impute_layer_robust_base(expr_micro_all_mat, groups_micro)
  expr_micro_all_imp <- impute_layer_mnar_finalize(expr_micro_stage1, groups_micro, downshift = 1.8)

  meta_micro <- samp_micro_all %>%
    select(sample_id, ExpGroup) %>%
    mutate(ExpGroup = factor(ExpGroup)) %>%
    filter(sample_id %in% colnames(expr_micro_all_imp))

  run_micro_contrast <- function(group_a, group_b) {
    grp <- droplevels(meta_micro$ExpGroup)
    if (!all(c(group_a, group_b) %in% levels(grp))) return(NULL)
    if (any(table(grp) < 2)) return(NULL)
    sel <- meta_micro$ExpGroup %in% c(group_a, group_b)
    cols <- meta_micro$sample_id[sel]
    grp2 <- droplevels(meta_micro$ExpGroup[sel])
    X <- expr_micro_all_imp[, cols, drop = FALSE]
    grp2 <- relevel(grp2, ref = group_b)  # A - B
    design <- model.matrix(~ 0 + grp2); colnames(design) <- levels(grp2)
    if (!(group_a %in% colnames(design))) return(NULL)
    contr <- limma::makeContrasts(contrasts = paste0(group_a, " - ", group_b), levels = design)
    fit  <- limma::lmFit(X, design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, contr), trend = TRUE, robust = TRUE)
    tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
    tt$protein_id <- rownames(tt)
    standardize_cols(tt)
  }

  micro_within <- list(
    con_vs_sus = list(figs = file.path(dirs$figs_micro, "within_group", "con_vs_sus"),
                      tabs = file.path(dirs$tabs_micro, "within_group", "con_vs_sus")),
    con_vs_res = list(figs = file.path(dirs$figs_micro, "within_group", "con_vs_res"),
                      tabs = file.path(dirs$tabs_micro, "within_group", "con_vs_res")),
    sus_vs_res = list(figs = file.path(dirs$figs_micro, "within_group", "sus_vs_res"),
                      tabs = file.path(dirs$tabs_micro, "within_group", "sus_vs_res")),
    res_vs_con = list(figs = file.path(dirs$figs_micro, "within_group", "res_vs_con"),
                      tabs = file.path(dirs$tabs_micro, "within_group", "res_vs_con"))
  )
  invisible(lapply(unlist(micro_within, recursive = FALSE), dir.create, showWarnings = FALSE, recursive = TRUE))

  write_micro_tables <- function(df, out_dir, stem, alpha = alpha_within) {
    if (!nrow(df)) return(invisible(NULL))
    df2 <- df %>% { if (!is.null(idmap_norm)) left_join(., idmap_norm, by = "protein_id") else . }
    df2_ranked <- df2 %>% arrange(padj, desc(abs(logFC)))
    df2_first  <- df2_ranked %>% distinct(protein_id, .keep_all = TRUE)  # keep-first per protein_id
    readr::write_csv(df2_first, file.path(out_dir, paste0(stem, "_full_firstOnly_", ts, ".csv")))
    sig <- df2_first %>% filter(padj <= alpha)
    readr::write_csv(sig, file.path(out_dir, paste0(stem, "_significant_firstOnly_", ts, ".csv")))
    top_up  <- sig %>% filter(logFC > 0) %>% arrange(padj, desc(logFC)) %>% slice_head(n = 50)
    top_dn  <- sig %>% filter(logFC < 0) %>% arrange(padj, logFC)       %>% slice_head(n = 50)
    readr::write_csv(top_up, file.path(out_dir, paste0(stem, "_top_up_50_firstOnly_", ts, ".csv")))
    readr::write_csv(top_dn, file.path(out_dir, paste0(stem, "_top_down_50_firstOnly_", ts, ".csv")))
    invisible(list(full = df2_first, sig = sig, top_up = top_up, top_down = top_dn))
  }

  micro_volcano <- function(df, figs_dir, stem, xlab, alpha = alpha_within) {
    if (!nrow(df)) return(invisible(NULL))
    df2 <- df %>% mutate(minus_log10_padj = -log10(padj), sig = padj <= alpha, sig_flag = ifelse(sig, "sig", "nonsig"))
    ymax <- max(df2$minus_log10_padj[is.finite(df2$minus_log10_padj)], na.rm = TRUE) * 1.05
    top_up <- df %>% filter(padj <= alpha, logFC > 0) %>% arrange(padj, desc(logFC)) %>% slice_head(n = 5)
    top_dn <- df %>% filter(padj <= alpha, logFC < 0) %>% arrange(padj, logFC)       %>% slice_head(n = 5)
    lab_df <- bind_rows(top_up, top_dn) %>%
      { if (!is.null(idmap_norm)) left_join(., idmap_norm, by = "protein_id") else mutate(., gene_symbol = NA_character_) } %>%
      mutate(label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, protein_id),
             minus_log10_padj = -log10(padj))
    p <- ggplot(df2, aes(x = logFC, y = minus_log10_padj)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
      geom_point(aes(color = sig_flag), alpha = 0.55, size = 1.8, na.rm = TRUE) +
      geom_point(data = ~ subset(., sig), color = "red", alpha = 0.95, size = 2.1, na.rm = TRUE) +
      ggrepel::geom_text_repel(
        data = lab_df, aes(label = label),
        color = "black", size = 3, max.overlaps = 200,
        box.padding = 0.25, point.padding = 0.25, segment.color = "grey50", segment.size = 0.3
      ) +
      scale_y_continuous(limits = c(0, ymax), name = "-log10 adj. p") +
      scale_x_continuous(name = xlab) +
      scale_color_manual(values = c(nonsig = "grey80", sig = "#2C6BED"), guide = "none") +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank(), axis.text = element_text(size = 12))
    ggsave(file.path(figs_dir, paste0(stem, "_", ts, ".svg")), p, width = 5, height = 4, device = "svg")
    ggsave(file.path(figs_dir, paste0(stem, "_", ts, ".png")), p, width = 5, height = 4, dpi = 300)
    invisible(p)
  }

  micro_contrasts <- list(
    list(a = "con", b = "sus", slot = "con_vs_sus", xlab = "log2 fold change (con − sus)"),
    list(a = "con", b = "res", slot = "con_vs_res", xlab = "log2 fold change (con − res)"),
    list(a = "sus", b = "res", slot = "sus_vs_res", xlab = "log2 fold change (sus − res)"),
    list(a = "res", b = "con", slot = "res_vs_con", xlab = "log2 fold change (res − con)")
  )

  for (cc in micro_contrasts) {
    df <- run_micro_contrast(cc$a, cc$b)
    if (is.null(df) || !nrow(df)) next
    write_micro_tables(df, micro_within[[cc$slot]]$tabs, stem = cc$slot, alpha = alpha_within)
    micro_volcano(df, micro_within[[cc$slot]]$figs, stem = cc$slot, xlab = cc$xlab, alpha = alpha_within)
  }
}

# ---------------- 2) Within-layer explicit contrasts (all samples in region) ----------------
samp_all_region <- samp %>% filter(region == region_of_interest, layer %in% layer_levels) %>% mutate(layer = factor(layer, levels = layer_levels))
message("Counts per layer x ExpGroup:"); print(table(samp_all_region$layer, samp_all_region$ExpGroup))

expr_all_region <- expr_long %>% semi_join(samp_all_region, by = "sample_id")
expr_wide_all <- expr_all_region %>% select(protein_id, sample_id, log2_intensity) %>%
  pivot_wider(names_from = sample_id, values_from = log2_intensity) %>%
  arrange(protein_id) %>% filter(!is.na(protein_id), protein_id != "")
expr_mat_all <- expr_wide_all %>% select(-protein_id) %>% as.matrix()
rownames(expr_mat_all) <- expr_wide_all$protein_id
expr_mat_all <- to_numeric_matrix_hard(expr_mat_all)

groups_all <- samp_all_region$layer[match(colnames(expr_mat_all), samp_all_region$sample_id)]
keep_all <- keep_idx(expr_mat_all, groups_all, 0.5)
expr_mat_all <- expr_mat_all[keep_all, , drop = FALSE]

expr_stage1_all <- impute_layer_robust_base(expr_mat_all, groups_all)
expr_mat_all_imp <- impute_layer_mnar_finalize(expr_stage1_all, groups_all, downshift = 1.8)

samp_sub <- samp_all_region %>% select(sample_id, layer, ExpGroup) %>% mutate(ExpGroup = factor(ExpGroup)) %>% filter(sample_id %in% colnames(expr_mat_all_imp))

do_within <- function(group_a, group_b, slot, xlab) {
  tables <- run_within_layer_contrast(expr_mat_all_imp, samp_sub, group_a, group_b)
  df <- if (length(tables)) bind_rows(tables) %>% standardize_cols() else tibble()
  if (!nrow(df)) return(invisible(NULL))
  save_result_tables(df, within_sub[[slot]]$tabs, stem = slot, alpha = alpha_within, idmap_norm = idmap_norm)
  save_top_updown_tables(df, within_sub[[slot]]$tabs, stem = slot, alpha = alpha_within, idmap_norm = idmap_norm, by_layer = TRUE, n = 50)
  lab_df <- df %>%
    { if (!is.null(idmap_norm)) left_join(., idmap_norm, by = "protein_id") else mutate(., gene_symbol = NA_character_) } %>%
    group_by(layer) %>% arrange(padj, desc(abs(logFC)), .by_group = TRUE) %>% slice_head(n = 5) %>%
    mutate(label = ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, protein_id)) %>% ungroup()
  build_volcano(df, xlab = xlab,
                out_svg = file.path(within_sub[[slot]]$figs, paste0(slot, "_", ts, ".svg")),
                out_png = file.path(within_sub[[slot]]$figs, paste0(slot, "_", ts, ".png")),
                alpha = alpha_within,
                labels_df = lab_df)
  invisible(df)
}

cs_df <- do_within("con","sus","con_vs_sus","log2 fold change (con − sus)")
cr_df <- do_within("con","res","con_vs_res","log2 fold change (con − res)")
sr_df <- do_within("sus","res","sus_vs_res","log2 fold change (sus − res)")
rc_df <- do_within("res","con","res_vs_con","log2 fold change (res − con)")

# ------------------------------- Save RDS -------------------------------
saveRDS(list(
  expr_controls = expr_region_con, samp_controls = samp_region_con,
  expr_mat_controls = expr_mat, expr_stage1_controls = expr_stage1, expr_mat_imp_controls = expr_mat_imp,
  layer_enrichment_stats = stats_df, labels_df = labels_df,
  expr_all = expr_all_region, samp_all = samp_all_region,
  expr_mat_all = expr_mat_all, expr_stage1_all = expr_stage1_all, expr_mat_all_imp = expr_mat_all_imp,
  con_vs_sus = cs_df, con_vs_res = cr_df, sus_vs_res = sr_df, res_vs_con = rc_df,
  layer_levels = layer_levels, alpha_layer = alpha_layer, alpha_within = alpha_within,
  region_of_interest = region_of_interest
), file = file.path(dirs$rds, paste0(region_of_interest, "_layer_within_micro_fdr005_", ts, ".rds")))

sink(file.path(dirs$logs, paste0("sessionInfo_", ts, ".txt"))); print(sessionInfo()); sink()
