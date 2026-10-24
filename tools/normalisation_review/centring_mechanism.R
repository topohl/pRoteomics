# Mechanism versus biology in the centring choice (normalisation review 2026-10-22). Read-only.
#
#   C  in-silico censoring on the real data: RES samples get SUS's excess missing values (lowest observed values
#      removed, or removed at random within abundance deciles), then both centrings are recomputed as the pipeline
#      does (anchor set over all samples of the compartment). How much of the observed SUS-RES offset difference does
#      censoring alone reproduce? Offset changes are relative to the uncensored samples (a smaller anchor set moves
#      every sample's anchor median)
#   R  relabelling: across the 18 units, does the original's median logFC track the missing-fraction difference for
#      every split of the six SUS/RES animals (a property of the method) or only for the true labelling (biology)?
#   W  within-animal association of the offset difference (original - new) with the missing fraction
#   S  simulation: censoring only, and a true global decrease that DIA-NN normalisation removes; what each centring
#      returns
# Writes results/reviewer_audit/normalisation_choice_2026-10-21/centring_mechanism/. Changes nothing else.
# USAGE  Rscript tools/normalisation_review/centring_mechanism.R      (working directory: repo root)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr); library(matrixStats) })
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
REV <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21")
OUT <- file.path(REV, "centring_mechanism"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
w <- function(x, name) write_csv(x, file.path(OUT, name), na = "")
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
SUSA <- c("A111", "A755", "A764"); RESA <- c("A0003", "A135", "A139")
OFF <- readRDS(file.path(REV, "validation", "validation_state.rds"))$offsets |> mutate(missing_fraction = unname(missing_fraction))
raw <- read_tsv("data/raw/pg_matrix/quicksearch.pg_matrix.tsv", show_col_types = FALSE, progress = FALSE)
X_all <- log2(as.matrix(raw[, -(1:4)])); rownames(X_all) <- raw$Protein.Group
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
grp_diff <- function(v, grp, a = "SUS", b = "RES") mean(v[grp == a]) - mean(v[grp == b])
# unit-adjusted animal-level SUS-RES difference of a per-sample quantity
animal_diff <- function(o, v) { a <- o |> mutate(v = v - ave(v, unit)) |> group_by(animal, group) |> summarise(v = mean(v), .groups = "drop"); grp_diff(a$v, a$group) }

# ---------------------------------------------------------------- C. in-silico censoring
cens <- list()
for (ds in DATASETS) {
  o <- OFF |> filter(dataset == ds); X <- X_all[, o$sample_id]; X <- X[rowSums(!is.na(X)) > 0, , drop = FALSE]
  mu <- rowMeans(X, na.rm = TRUE); d10 <- cut(rank(mu, ties.method = "first"), 10, labels = FALSE)
  # SUS excess per unit (total) and per decile (shape), from the data
  ex_unit <- o |> group_by(unit, animal, group) |> summarise(m = mean(missing_fraction), .groups = "drop") |> group_by(unit) |>
    summarise(excess = mean(m[group == "SUS"]) - mean(m[group == "RES"]), .groups = "drop")
  dec_shape <- sapply(1:10, function(k) { mk <- colMeans(is.na(X[d10 == k, , drop = FALSE])); a <- tapply(mk - ave(mk, o$unit), o$animal, mean)
    g <- o$group[match(names(a), o$animal)]; max(grp_diff(a, g), 0) })
  dec_shape <- dec_shape / sum(dec_shape)
  offs <- function(Xc) { anc <- rowSums(is.na(Xc)) == 0; list(orig = colMedians(Xc, na.rm = TRUE), new = colMedians(Xc[anc, , drop = FALSE]), n_anc = sum(anc)) }
  base <- offs(X)
  obs_diff <- c(original = animal_diff(o, base$orig), new = animal_diff(o, base$new))
  for (mech in c("lowest observed values", "random within abundance deciles")) for (rep in 1:(if (mech == "lowest observed values") 1 else 10)) {
    set.seed(100 + rep); Xc <- X
    for (j in which(o$group == "RES")) {
      e <- max(ex_unit$excess[ex_unit$unit == o$unit[j]], 0); n_rm <- round(e * nrow(X)); if (n_rm < 1) next
      ob <- which(!is.na(Xc[, j]))
      rm_idx <- if (mech == "lowest observed values") ob[order(Xc[ob, j])][seq_len(n_rm)] else
        unlist(lapply(1:10, function(k) { cand <- ob[d10[ob] == k]; nk <- min(length(cand), round(n_rm * dec_shape[k])); if (nk > 0) cand[sample.int(length(cand), nk)] }))
      Xc[rm_idx, j] <- NA
    }
    cc <- offs(Xc); oc <- o |> mutate(missing_fraction = colMeans(is.na(Xc)))
    cens[[length(cens) + 1L]] <- tibble(dataset = ds, mechanism = mech, replicate = rep,
      added_missing_RES_mean = mean(colMeans(is.na(Xc))[o$group == "RES"] - colMeans(is.na(X))[o$group == "RES"]),
      observed_SUS_minus_RES_missing = animal_diff(o, o$missing_fraction), censored_SUS_minus_RES_missing = animal_diff(oc, oc$missing_fraction),
      RES_offset_change_original = mean((cc$orig - base$orig)[o$group == "RES"]) - mean((cc$orig - base$orig)[o$group != "RES"]),
      RES_offset_change_new = mean((cc$new - base$new)[o$group == "RES"]) - mean((cc$new - base$new)[o$group != "RES"]),
      common_offset_change_new = mean((cc$new - base$new)[o$group != "RES"]),
      n_anchors_before = base$n_anc, n_anchors_after = cc$n_anc,
      observed_SUS_minus_RES_offset_original = obs_diff[["original"]], observed_SUS_minus_RES_offset_new = obs_diff[["new"]],
      share_of_observed_original_difference_reproduced = (mean((cc$orig - base$orig)[o$group == "RES"]) - mean((cc$orig - base$orig)[o$group != "RES"])) / obs_diff[["original"]])
  }
}
cens <- bind_rows(cens)
w(cens, "in_silico_censoring.csv")
cens_s <- cens |> group_by(dataset, mechanism) |> summarise(across(c(added_missing_RES_mean, RES_offset_change_original, RES_offset_change_new,
  observed_SUS_minus_RES_offset_original, observed_SUS_minus_RES_offset_new, share_of_observed_original_difference_reproduced, n_anchors_after, common_offset_change_new), mean), .groups = "drop")
w(cens_s, "in_silico_censoring_summary.csv")

# ---------------------------------------------------------------- R. relabelling the depth-logFC relation
rel <- list(); per_unit <- list()
splits <- combn(c(SUSA, RESA), 3, simplify = FALSE); splits <- splits[sapply(splits, function(s) "A111" %in% s)]   # 10 unordered splits
GS <- lapply(setNames(DATASETS, DATASETS), function(ds) list(new = read_gct(gct_path(".", ds)), original = read_gct(gct_path(ARCHIVE, ds))))
for (s in splits) {
  rows <- list()
  for (ds in DATASETS) {
    o <- OFF |> filter(dataset == ds)
    for (v in names(GS[[ds]])) {
      g <- GS[[ds]][[v]]; an <- g$coldesc$AnimalID; un <- sub("_[0-9]+$", "", g$coldesc$region_layer_ExpGroup)
      for (u in unique(un)) {
        A <- g$m[, un == u & an %in% s, drop = FALSE]; B <- g$m[, un == u & an %in% setdiff(c(SUSA, RESA), s), drop = FALSE]
        om <- o |> filter(unit == u | sub("_microglia$", "", unit) == u) |> group_by(animal) |> summarise(m = mean(missing_fraction), .groups = "drop")
        stopifnot(nrow(om) == 9)
        rows[[length(rows) + 1L]] <- tibble(dataset = ds, version = v, unit = u,
          missing_diff = mean(om$m[om$animal %in% s]) - mean(om$m[om$animal %in% setdiff(c(SUSA, RESA), s)]),
          median_logFC = median(rowMeans(A) - rowMeans(B), na.rm = TRUE))
      }
    }
  }
  rows <- bind_rows(rows) |> mutate(split = paste(s, collapse = "+"), is_true_labelling = setequal(s, SUSA))
  per_unit[[length(per_unit) + 1L]] <- rows
  rel[[length(rel) + 1L]] <- rows |> group_by(split, is_true_labelling, version) |>
    summarise(n_units = n(), r_missing_vs_median_logFC = cor(missing_diff, median_logFC), slope = coef(lm(median_logFC ~ missing_diff))[2],
              sd_missing_diff = sd(missing_diff), .groups = "drop")
}
rel <- bind_rows(rel); w(rel, "relabelling_depth_logFC.csv"); w(bind_rows(per_unit), "relabelling_depth_logFC_by_unit.csv")

# ---------------------------------------------------------------- W. within-animal association
wa <- OFF |> group_by(dataset) |> mutate(d = raw_original_all_observed - raw_new_complete_case) |>
  group_by(dataset, unit) |> mutate(d_u = d - mean(d), m_u = missing_fraction - mean(missing_fraction)) |>
  group_by(dataset, animal) |> mutate(d_w = d_u - mean(d_u), m_w = m_u - mean(m_u)) |> group_by(dataset) |>
  summarise(r_unit_adjusted = cor(d_u, m_u), r_within_animal = cor(d_w, m_w), slope_within_animal = coef(lm(d_w ~ m_w))[2],
            n_samples = n(), n_animals_positive_within = sum(tapply(seq_along(d_w), animal, function(i) cor(d_w[i], m_w[i]) > 0)), .groups = "drop")
w(wa, "within_animal_offset_vs_missing.csv")

# ---------------------------------------------------------------- S. simulation (normal abundance, detection threshold)
set.seed(7); n <- 5000; reps <- 200
sim <- list()
for (scen in c("censoring only (detection limit raised)", "true global decrease, removed by DIA-NN normalisation")) for (delta in c(0.05, 0.1, 0.2)) {
  r <- replicate(reps, {
    x <- rnorm(n, 0, 1.5); L <- quantile(x, 0.2)                         # 20% missing in the reference group
    a <- ifelse(x + rnorm(n, 0, 0.2) > L, x, NA)
    if (startsWith(scen, "censoring")) { b <- ifelse(x + rnorm(n, 0, 0.2) > L + delta, x, NA) }
    else { b0 <- x - delta; b <- ifelse(b0 + rnorm(n, 0, 0.2) > L, b0, NA)
           ok <- !is.na(a) & !is.na(b); b <- b + median(a[ok] - b[ok]) }   # DIA-NN-like rescaling on well-detected values
    anc <- !is.na(a) & !is.na(b)
    c(missing_excess = mean(is.na(b)) - mean(is.na(a)), original_logFC_shift = -(median(b, na.rm = TRUE) - median(a, na.rm = TRUE)),
      new_logFC_shift = -(median(b[anc]) - median(a[anc])))
  })
  sim[[length(sim) + 1L]] <- tibble(scenario = scen, delta = delta, true_global_logFC = if (startsWith(scen, "censoring")) 0 else -delta,
    missing_excess = mean(r["missing_excess", ]), original_recovers = mean(r["original_logFC_shift", ]), new_recovers = mean(r["new_logFC_shift", ]))
}
sim <- bind_rows(sim); w(sim, "simulation.csv")

options(width = 220)
print(as.data.frame(cens_s)); print(as.data.frame(rel |> arrange(version, r_missing_vs_median_logFC))); print(as.data.frame(wa)); print(as.data.frame(sim))
cat("done:", OUT, "\n")
