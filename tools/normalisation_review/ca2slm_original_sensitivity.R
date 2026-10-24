# CA2-SLM SUS vs RES under the original centring: sample exclusion, missingness/imputation and variance moderation
# (normalisation review 2026-10-22). Read-only. The original-centred animal-level input is rebuilt as the current input
# plus the mean offset change of the kept hemispheres (verified equal to the archived input, validation/), so every
# exclusion is applied to exactly the pipeline's aggregation. Writes original_assessment/ca2slm_*.csv.
#
#   exclusions  none; single hemispheres A755 L, A764 R; both; whole SUS animals (A755, A764, A111) and RES animals
#               (A0003, A135, A139) removed from CA2-SLM
#   models      pipeline (limma, all units, ~0 + unit x group, eBayes robust, no trend); eBayes without robust; with
#               intensity trend; CA2-SLM-only model; no moderation (ordinary t from the same fit)
#   imputation  observed-only: imputed cells set to missing before the hemisphere mean, limma on the observed values
# USAGE  Rscript tools/normalisation_review/ca2slm_original_sensitivity.R      (working directory: repo root)

suppressPackageStartupMessages({ library(readr); library(readxl); library(dplyr); library(tidyr); library(limma) })
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
REV <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21")
OUT <- file.path(REV, "original_assessment"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(x, name) write_csv(x, file.path(OUT, name), na = "")
state <- readRDS(file.path(REV, "validation", "validation_state.rds"))
ds <- "neuron_neuropil"
o <- state$offsets |> filter(dataset == ds)
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
imp_path <- Sys.glob(sprintf("data/processed/01_preprocessing/impute/*_pgmatrix_imputed_%s_*_missing70pct.xlsx", ds))
imp <- read_excel(imp_path, progress = FALSE)
mn <- as.matrix(imp[, o$sample_id]); rownames(mn) <- imp$`T: Protein.Names`
rawm <- read_tsv("data/raw/pg_matrix/quicksearch.pg_matrix.tsv", show_col_types = FALSE, progress = FALSE)
obs <- !is.na(as.matrix(rawm[match(imp$Protein.Group, rawm$Protein.Group), o$sample_id])); rownames(obs) <- rownames(mn)
g <- read_gct(sprintf("data/processed/01_preprocessing/protigy_input_animal_level/%s/%s_animal_level.gct", ds, ds))
g_old <- read_gct(file.path(ARCHIVE, "data/processed/01_preprocessing/protigy_input_animal_level", ds, paste0(ds, "_animal_level.gct")))
col_key <- paste(g$coldesc$AnimalID, g$coldesc$region, g$coldesc$layer, sep = "|")
smp_key <- paste(o$animal, sub("_.*$", "", o$unit), sub("^[^_]+_", "", o$unit), sep = "|")
cell <- g$coldesc$region_layer_ExpGroup
shift_orig <- o$raw_new_complete_case - o$raw_original_all_observed   # current + shift = original centring, per sample

# animal-level matrix under a centring ("original" or "new") with some CA2-SLM hemispheres removed; optional observed-only
animal_matrix <- function(drop_labels = character(), version = "original", observed_only = FALSE) {
  keep <- !(o$label %in% drop_labels)
  sm <- mn; if (observed_only) sm[!obs] <- NA
  if (version == "original") sm <- sweep(sm, 2, shift_orig, `+`)
  m0 <- matrix(NA_real_, nrow(g$m), ncol(g$m), dimnames = dimnames(g$m))
  for (j in seq_along(col_key)) {
    s <- o$sample_id[smp_key == col_key[j] & keep]
    if (!length(s)) next
    v <- rowMeans(sm[, s, drop = FALSE], na.rm = TRUE); v[is.nan(v)] <- NA; v <- v[!duplicated(names(v))]
    m0[, j] <- v[rownames(m0)]
  }
  m0
}
# sanity: the rebuilt original-centred input equals the archived input
rid <- intersect(rownames(g$m), rownames(g_old$m))
stopifnot(max(abs(animal_matrix()[rid, ] - g_old$m[rid, colnames(g$m)]), na.rm = TRUE) < 1e-8)

fit_ca2 <- function(m, model = "pipeline") {
  keepc <- colSums(!is.na(m)) > 0; m <- m[, keepc, drop = FALSE]; cl <- cell[keepc]
  if (model == "ca2_only") { sel <- grepl("^CA2_slm_", cl); m <- m[, sel, drop = FALSE]; cl <- cl[sel] }
  f <- factor(cl, levels = unique(cl)); design <- model.matrix(~0 + f); colnames(design) <- levels(f)
  cm <- makeContrasts(contrasts = "CA2_slm_3-CA2_slm_2", levels = design)
  fit <- contrasts.fit(lmFit(m, design), cm)
  if (model == "no_moderation") {
    t <- fit$coefficients[, 1] / (fit$stdev.unscaled[, 1] * fit$sigma); p <- 2 * pt(-abs(t), fit$df.residual)
    return(tibble(protein = rownames(m), logFC = fit$coefficients[, 1], t = t, P.Value = p, adj.P.Val = p.adjust(p, "BH")))
  }
  eb <- eBayes(fit, robust = model != "not_robust", trend = model == "trend")
  tt <- topTable(eb, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
  tibble(protein = rownames(tt), logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, df_prior = eb$df.prior[1])
}
excl <- list("none" = character(), "without A755 L" = "A755_L_CA2_slm", "without A764 R" = "A764_R_CA2_slm",
             "without A755 L and A764 R" = c("A755_L_CA2_slm", "A764_R_CA2_slm"),
             "without animal A755 (SUS)" = c("A755_L_CA2_slm", "A755_R_CA2_slm"), "without animal A764 (SUS)" = c("A764_L_CA2_slm", "A764_R_CA2_slm"),
             "without animal A111 (SUS)" = c("A111_L_CA2_slm", "A111_R_CA2_slm"), "without animal A0003 (RES)" = c("A0003_L_CA2_slm", "A0003_R_CA2_slm"),
             "without animal A135 (RES)" = c("A135_L_CA2_slm", "A135_R_CA2_slm"), "without animal A139 (RES)" = c("A139_L_CA2_slm", "A139_R_CA2_slm"))
runs <- list()
for (v in c("original", "new")) {
  for (en in names(excl)) runs[[paste(v, en, "pipeline", "imputed", sep = " | ")]] <- fit_ca2(animal_matrix(excl[[en]], v), "pipeline")
  for (md in c("not_robust", "trend", "ca2_only", "no_moderation")) runs[[paste(v, "none", md, "imputed", sep = " | ")]] <- fit_ca2(animal_matrix(character(), v), md)
  for (en in c("none", "without A755 L and A764 R")) runs[[paste(v, en, "pipeline", "observed only", sep = " | ")]] <- fit_ca2(animal_matrix(excl[[en]], v, TRUE), "pipeline")
}
summ <- bind_rows(lapply(names(runs), function(k) { p <- strsplit(k, " | ", fixed = TRUE)[[1]]; r <- runs[[k]]
  tibble(centring = p[1], exclusion = p[2], model = p[3], values = p[4], n_tested = sum(!is.na(r$adj.P.Val)), n_fdr05 = sum(r$adj.P.Val < 0.05, na.rm = TRUE),
         n_fdr05_down = sum(r$adj.P.Val < 0.05 & r$logFC < 0, na.rm = TRUE), n_fdr10 = sum(r$adj.P.Val < 0.10, na.rm = TRUE),
         median_logFC = median(r$logFC, na.rm = TRUE), df_prior = r$df_prior[1] %||% NA_real_) }))
w(summ, "ca2slm_original_sensitivity_summary.csv")

# the 28 original hits across every variant
base <- runs[["original | none | pipeline | imputed"]]; hits <- base$protein[base$adj.P.Val < 0.05]
ca2 <- o$unit == "CA2_slm"
per_hit <- bind_rows(lapply(names(runs), function(k) runs[[k]] |> filter(protein %in% hits) |> mutate(run = k))) |>
  select(protein, run, logFC, adj.P.Val) |>
  pivot_wider(names_from = run, values_from = c(logFC, adj.P.Val))
hit_tab <- tibble(protein = hits) |>
  mutate(logFC_original = base$logFC[match(protein, base$protein)], padj_original = base$adj.P.Val[match(protein, base$protein)],
         n_imputed_SUS_CA2slm = rowSums(!obs[protein, ca2 & o$group == "SUS", drop = FALSE]), n_imputed_RES_CA2slm = rowSums(!obs[protein, ca2 & o$group == "RES", drop = FALSE]),
         imputed_A755L = !obs[protein, o$label == "A755_L_CA2_slm"], imputed_A764R = !obs[protein, o$label == "A764_R_CA2_slm"]) |>
  left_join(per_hit, by = "protein")
lf <- select(hit_tab, starts_with("logFC_original |")); pa <- select(hit_tab, starts_with("adj.P.Val_original |"))
hit_tab <- hit_tab |> mutate(n_variants_fdr05 = rowSums(pa < 0.05, na.rm = TRUE), n_variants = ncol(pa),
                             sign_stable_all_variants = apply(sign(lf) == sign(logFC_original), 1, all, na.rm = TRUE),
                             same_sign_new_centring = sign(`logFC_new | none | pipeline | imputed`) == sign(logFC_original),
                             fdr05_new_centring = `adj.P.Val_new | none | pipeline | imputed` < 0.05,
                             exceeds_offset_shift = abs(logFC_original) > abs(median(base$logFC)))
w(hit_tab, "ca2slm_original_hits_across_variants.csv")
print(as.data.frame(summ))
print(hit_tab |> summarise(n = n(), down = sum(logFC_original < 0), sign_stable = sum(sign_stable_all_variants), same_sign_new = sum(same_sign_new_centring),
                           fdr_new = sum(fdr05_new_centring), exceeds_shift = sum(exceeds_offset_shift), any_imputed = sum(n_imputed_SUS_CA2slm + n_imputed_RES_CA2slm > 0)))
cat("done\n")
