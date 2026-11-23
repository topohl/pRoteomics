#!/usr/bin/env Rscript
#
#
# Where each WGCNA module sits in the control hippocampus (added 2026-11-23; Figure 4 candidate).
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/describe_module_control_geography.R
#   Rscript analysis/spatial_validation/describe_module_control_geography.R --dry-run
# Script: analysis/spatial_validation/describe_module_control_geography.R
# Stage: networks
# Scope: global
# Consumes: required results/wgcna/test_module_phenotypes/<dataset>/tables/WGCNA_group_effect_hemisphere_values.csv; results/spatial_validation/map_input_territories/global/tables/input_territory_carryover_qc.csv
# Produces: results/spatial_validation/describe_module_control_geography/global/tables/module_control_spatial_class.csv; module_control_profiles.csv; neuropil_module_territory_plane.csv; and the same three tables with all acquisitions (check)
# Dataset behavior: global; loops over neuron_neuropil, neuron_soma and microglia.
# Notes: Descriptive and phenotype-blind: control animals only, no group effect, no WGCNA refit.
#
# VALUES: the stored module eigengene of every acquisition (Stage 05 hemisphere values; one acquisition per animal,
# hemisphere and unit). The DG po acquisitions with tissue carry-over are left out as in Figure 2j-l
# (input_territory_carryover_qc.csv, column excluded); the tables named *_all_acquisitions keep them (check).
#
# SPATIAL CLASS (every module): over the control animals, the unit with the highest mean of the animal means is the
# peak. The unit effect is the partial R2 of unit in lm(eigengene ~ animal + unit) over the control acquisitions; its
# P comes from shuffling the unit labels within each animal and hemisphere (4,000 shuffles; 50,000 when the P lies
# between 0.04 and 0.06, where the Monte Carlo error could change the class). Classes: no structure (P >= 0.05);
# dominant (the peak is every control animal's peak and at least five of the six hemispheres peak there; only
# hemispheres with every unit vote, so a hemisphere that lost an acquisition to the leave-out counts as not peaking);
# otherwise shared top (units within two standard errors of the peak). tau (specificity of the control profile,
# 0 = flat, 1 = one unit) is reported with its shuffle null (tau_null95, p_tau); the class does not use it.
#
# PROFILES: per unit, the control mean of the animal means, its standard error and that mean minus the module's mean
# over units in eigengene SD. The SD is that of all acquisitions of the compartment (every eigengene is a unit vector
# over the acquisitions the network was built on, so the SD is the same for every module of a compartment).
#
# TERRITORIES (neuropil only): perforant path (PP: CA1 SLM, CA2 SLM, DG MO) and mossy fibre (MF: DG PO, CA3 SR) as in
# Figure 2j, each minus the mean of the other five laminae, per control animal, with a t interval (df = animals - 1),
# in eigengene SD. The territory and region R2 of the control profile are tested by exact enumeration of every
# distinct assignment of the unit labels (2,520 and 25,200); Benjamini-Hochberg over the modules.
#
#

source("R/paths.R")
source("R/statistics/integration_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))
source("R/spatial/neuropil_lamina_utils.R")

ANALYSIS_ID <- "describe_module_control_geography"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)
SCRIPT_ID <- "analysis/spatial_validation/describe_module_control_geography.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
HEMI <- stats::setNames(vapply(DATASETS, function(d) repo_path("results", "wgcna", "test_module_phenotypes", d, "tables",
                                                              "WGCNA_group_effect_hemisphere_values.csv"), character(1)), DATASETS)
CARRYOVER <- repo_path("results", "spatial_validation", "map_input_territories", "global", "tables", "input_territory_carryover_qc.csv")
tab <- function(f) file.path(CANONICAL_PATHS$tables, f)
OUT <- list(
  class = tab("module_control_spatial_class.csv"),
  profiles = tab("module_control_profiles.csv"),
  territory = tab("neuropil_module_territory_plane.csv"),
  class_all = tab("module_control_spatial_class_all_acquisitions.csv"),
  profiles_all = tab("module_control_profiles_all_acquisitions.csv"),
  territory_all = tab("neuropil_module_territory_plane_all_acquisitions.csv")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")

N_SHUFFLE <- 4000L
N_SHUFFLE_BORDER <- 50000L
BORDER <- c(0.04, 0.06)
SEED <- 77031L
PP <- c("ca1_slm", "ca2_slm", "dg_mo")
MF <- c("dg_po", "ca3_sr")
OTHER <- c("ca1_so", "ca1_sr", "ca2_so", "ca2_sr", "ca3_so")
LAMINAE <- c("ca1_so", "ca1_sr", "ca1_slm", "ca2_so", "ca2_sr", "ca2_slm", "ca3_so", "ca3_sr", "dg_mo", "dg_po")

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Control spatial class, control profiles and neuropil territory plane of every WGCNA module.\n")
  dry_run_inputs(SCRIPT_ID, c(as.list(HEMI), list(carryover = CARRYOVER)))
  quit(save = "no", status = 0L)
}
for (f in c(HEMI, CARRYOVER)) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))

carry <- utils::read.csv(CARRYOVER, stringsAsFactors = FALSE)
left_out <- carry$Sample[carry$excluded %in% TRUE]
if (!length(left_out)) stop("The carry-over table leaves no acquisition out; check its 'excluded' column.", call. = FALSE)

read_hemi <- function(d) {
  h <- utils::read.csv(HEMI[[d]], stringsAsFactors = FALSE)
  h <- h[h$level == "module", c("endpoint_id", "AnimalID", "StressGroup", "SpatialUnit", "Hemisphere",
                                "contributing_sample_ids", "source_sample_count", "hemisphere_value")]
  if (any(h$source_sample_count != 1L)) stop(d, ": a hemisphere value pools several acquisitions.", call. = FALSE)
  h$SpatialUnit <- tolower(h$SpatialUnit)
  h
}
hemi <- lapply(stats::setNames(DATASETS, DATASETS), read_hemi)
control_animals <- sort(unique(unlist(lapply(hemi, function(h) h$AnimalID[h$StressGroup == "CON"]))))
for (d in DATASETS) if (!setequal(unique(hemi[[d]]$AnimalID[hemi[[d]]$StressGroup == "CON"]), control_animals))
  stop(d, ": the control animals differ between compartments.", call. = FALSE)
if (length(control_animals) < 3L) stop("Fewer than three control animals.", call. = FALSE)
n_left_out_control <- sum(hemi$neuron_neuropil$contributing_sample_ids %in% left_out &
                          hemi$neuron_neuropil$StressGroup == "CON") / length(unique(hemi$neuron_neuropil$endpoint_id))

tau_of <- function(p) { w <- p - min(p); if (max(w) <= 0) return(0); sum(1 - w / max(w)) / (length(p) - 1) }
# control profile: mean over the animals of each animal's mean over its hemispheres
animal_means <- function(v, a, u, units) {
  m <- tapply(v, list(a, factor(u, levels = units)), mean)
  if (anyNA(m)) stop("A control animal lacks a unit in both hemispheres.", call. = FALSE)
  m
}

spatial_class <- function(d, h, sd_unit) {
  x0 <- h[h$StressGroup == "CON", , drop = FALSE]
  units <- sort(unique(x0$SpatialUnit))
  rows <- list(); prof <- list()
  mods <- sort(unique(x0$endpoint_id))
  for (i in seq_along(mods)) {
    m <- mods[[i]]
    x <- x0[x0$endpoint_id == m, , drop = FALSE]
    x <- x[order(x$AnimalID, x$Hemisphere, x$SpatialUnit), , drop = FALSE]
    B <- animal_means(x$hemisphere_value, x$AnimalID, x$SpatialUnit, units)
    p <- colMeans(B)
    o <- order(-p); peak <- units[o[1]]; second <- units[o[2]]
    animal_peak <- units[apply(B, 1, which.max)]
    hemis <- split(x, paste(x$AnimalID, x$Hemisphere))
    hemis <- hemis[vapply(hemis, function(z) setequal(z$SpatialUnit, units), logical(1))]   # complete hemispheres vote
    hp <- vapply(hemis, function(z) z$SpatialUnit[which.max(z$hemisphere_value)], character(1))
    margin <- B[, peak] - apply(B[, setdiff(units, peak), drop = FALSE], 1, max)
    se_peak <- stats::sd(B[, peak]) / sqrt(nrow(B))
    peak_set <- units[p >= p[peak] - 2 * se_peak]
    Xa <- stats::model.matrix(~ AnimalID, x)
    rss_a <- sum(.lm.fit(Xa, x$hemisphere_value)$residuals^2)
    lev <- units[-1]
    pr2 <- function(uu) {
      Xu <- cbind(Xa, vapply(lev, function(l) as.numeric(uu == l), numeric(length(uu))))
      (rss_a - sum(.lm.fit(Xu, x$hemisphere_value)$residuals^2)) / rss_a
    }
    obs <- pr2(x$SpatialUnit)
    blocks <- split(seq_len(nrow(x)), paste(x$AnimalID, x$Hemisphere))
    shuffle <- function(n) {
      set.seed(SEED + 1000L * match(d, DATASETS) + i)   # one stream per compartment and module
      r2 <- numeric(n); tau <- numeric(n)
      for (k in seq_len(n)) {
        uu <- x$SpatialUnit
        for (b in blocks) uu[b] <- uu[b][sample.int(length(b))]
        r2[k] <- pr2(uu)
        tau[k] <- tau_of(colMeans(animal_means(x$hemisphere_value, x$AnimalID, uu, units)))
      }
      list(r2 = r2, tau = tau, n = n)
    }
    null <- shuffle(N_SHUFFLE)
    p_unit <- (1 + sum(null$r2 >= obs - 1e-12)) / (null$n + 1)
    if (p_unit > BORDER[1] && p_unit < BORDER[2]) {
      null <- shuffle(N_SHUFFLE_BORDER)
      p_unit <- (1 + sum(null$r2 >= obs - 1e-12)) / (null$n + 1)
    }
    tau <- tau_of(p)
    dom <- all(animal_peak == peak) && sum(hp == peak) >= 2L * nrow(B) - 1L && all(margin > 0)
    cls <- if (p_unit >= 0.05) "no_structure" else if (dom) "dominant" else "shared_top"
    rows[[m]] <- data.frame(dataset = d, module_id = m, peak_unit = peak, second_unit = second,
      animals_peak = sum(animal_peak == peak), n_animals = nrow(B), animal_peaks = paste(animal_peak, collapse = ";"),
      hemis_peak = sum(hp == peak), n_hemispheres = length(hp), margin_min = min(margin),
      shared_top_units = paste(peak_set, collapse = ";"), unit_partial_R2 = obs, p_unit = p_unit, n_shuffles = null$n,
      class_rule = dom, spatial_class = cls, tau = tau, tau_null95 = unname(stats::quantile(null$tau, 0.95)),
      p_tau = (1 + sum(null$tau >= tau - 1e-12)) / (null$n + 1), stringsAsFactors = FALSE)
    con_se <- apply(B, 2, stats::sd) / sqrt(nrow(B))
    prof[[m]] <- data.frame(dataset = d, module_id = m, spatial_unit = units, con_mean = unname(p), con_se = unname(con_se),
      n_con_animals = nrow(B), con_z = unname((p - mean(p)) / stats::sd(p)), sd_unit = unname(sd_unit[[m]]),
      con_sd = unname((p - mean(p)) / sd_unit[[m]]), stringsAsFactors = FALSE)
  }
  list(class = do.call(rbind, rows), profiles = do.call(rbind, prof))
}

# every distinct assignment of a label vector to positions (multinomial enumeration)
enumerate_assignments <- function(lab) {
  sizes <- table(lab); grp <- names(sizes); n <- length(lab); out <- list()
  rec <- function(avail, k, cur) {
    if (k == length(grp)) { cur[avail] <- grp[k]; out[[length(out) + 1L]] <<- cur; return(invisible()) }
    for (cc in utils::combn(avail, sizes[[grp[k]]], simplify = FALSE)) {
      c2 <- cur; c2[cc] <- grp[k]; rec(setdiff(avail, cc), k + 1L, c2)
    }
  }
  rec(seq_len(n), 1L, rep(NA_character_, n))
  do.call(rbind, out)
}
r2_groups <- function(p, g) { tot <- sum((p - mean(p))^2); sum(tapply(p, g, function(z) length(z) * (mean(z) - mean(p))^2)) / tot }

territory_plane <- function(h, sd_unit) {
  x0 <- h[h$StressGroup == "CON", , drop = FALSE]
  if (!setequal(unique(x0$SpatialUnit), LAMINAE)) stop("The neuropil laminae differ from the ten of Figure 2.", call. = FALSE)
  ter <- stats::setNames(ifelse(LAMINAE %in% PP, "PP", ifelse(LAMINAE %in% MF, "MF", "other")), LAMINAE)
  reg <- stats::setNames(sub("_.*", "", LAMINAE), LAMINAE)
  A_ter <- enumerate_assignments(unname(ter)); A_reg <- enumerate_assignments(unname(reg))
  ci <- function(z) { n <- length(z); hw <- stats::qt(0.975, n - 1) * stats::sd(z) / sqrt(n); c(mean(z), mean(z) - hw, mean(z) + hw) }
  rows <- list()
  for (m in sort(unique(x0$endpoint_id))) {
    x <- x0[x0$endpoint_id == m, , drop = FALSE]
    B <- animal_means(x$hemisphere_value, x$AnimalID, x$SpatialUnit, LAMINAE)
    p <- colMeans(B); s <- sd_unit[[m]]
    pp <- rowMeans(B[, PP, drop = FALSE]) - rowMeans(B[, OTHER, drop = FALSE])
    mf <- rowMeans(B[, MF, drop = FALSE]) - rowMeans(B[, OTHER, drop = FALSE])
    cp <- ci(pp) / s; cm <- ci(mf) / s
    o_t <- r2_groups(p, ter); o_r <- r2_groups(p, reg)
    nt <- apply(A_ter, 1, function(g) r2_groups(p, g)); nr <- apply(A_reg, 1, function(g) r2_groups(p, g))
    d0 <- data.frame(a = rep(rownames(B), times = ncol(B)), u = rep(colnames(B), each = nrow(B)), v = as.vector(B))
    d0$t <- ter[d0$u]
    f0 <- stats::lm(v ~ a, d0); f1 <- stats::lm(v ~ a + t, d0); f2 <- stats::lm(v ~ a + u, d0)
    rows[[m]] <- data.frame(dataset = "neuron_neuropil", module_id = m,
      pp_minus_other = cp[1], pp_ci_low = cp[2], pp_ci_high = cp[3],
      mf_minus_other = cm[1], mf_ci_low = cm[2], mf_ci_high = cm[3],
      pp_n_animals_positive = sum(pp > 0), mf_n_animals_positive = sum(mf > 0), n_animals = nrow(B),
      reference_laminae = paste(c("CA1_so", "CA1_sr", "CA2_so", "CA2_sr", "CA3_so"), collapse = ";"),
      territory_R2 = o_t, p_territory_exact = mean(nt >= o_t - 1e-12), n_territory_assignments = nrow(A_ter),
      region_R2 = o_r, p_region_exact = mean(nr >= o_r - 1e-12), n_region_assignments = nrow(A_reg),
      territory_share_of_unit = (stats::deviance(f0) - stats::deviance(f1)) / (stats::deviance(f0) - stats::deviance(f2)),
      sd_unit = s, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$q_territory <- stats::p.adjust(out$p_territory_exact, "BH")
  out$q_region <- stats::p.adjust(out$p_region_exact, "BH")
  out
}

# the eigengene SD is that of every acquisition the network was built on, also where some are left out of the means
sd_all <- lapply(hemi, function(h) tapply(h$hemisphere_value, h$endpoint_id, stats::sd))

run_variant <- function(drop) {
  hh <- lapply(hemi, function(h) if (drop) h[!h$contributing_sample_ids %in% left_out, , drop = FALSE] else h)
  sc <- lapply(DATASETS, function(d) spatial_class(d, hh[[d]], sd_all[[d]]))
  cls <- do.call(rbind, lapply(sc, `[[`, "class"))
  prof <- do.call(rbind, lapply(sc, `[[`, "profiles"))
  prof <- merge(prof, cls[c("dataset", "module_id", "spatial_class", "peak_unit")], by = c("dataset", "module_id"))
  prof$is_dominant_unit <- prof$spatial_class == "dominant" & prof$spatial_unit == prof$peak_unit
  terr <- territory_plane(hh$neuron_neuropil, sd_all$neuron_neuropil)
  list(class = cls[order(cls$dataset, cls$module_id), ], profiles = prof[order(prof$dataset, prof$module_id, prof$spatial_unit), ],
       territory = terr[order(terr$module_id), ])
}

primary <- run_variant(TRUE)
check <- run_variant(FALSE)
note <- "control animals only; DG po acquisitions with tissue carry-over left out as in Figure 2j-l"
note_all <- "check: control animals only, every acquisition (DG po carry-over acquisitions kept)"
for (k in names(primary)) { primary[[k]]$values <- note; check[[k]]$values <- note_all }
utils::write.csv(primary$class, OUT$class, row.names = FALSE)
utils::write.csv(primary$profiles, OUT$profiles, row.names = FALSE)
utils::write.csv(primary$territory, OUT$territory, row.names = FALSE)
utils::write.csv(check$class, OUT$class_all, row.names = FALSE)
utils::write.csv(check$profiles, OUT$profiles_all, row.names = FALSE)
utils::write.csv(check$territory, OUT$territory_all, row.names = FALSE)

write_run_manifest(MANIFEST, inputs = c(as.list(HEMI), list(carryover = CARRYOVER)), outputs = OUT,
                   parameters = list(control_animals = control_animals, left_out_acquisitions = left_out,
                                     control_acquisitions_left_out = n_left_out_control,
                                     shuffles = N_SHUFFLE, shuffles_near_0.05 = N_SHUFFLE_BORDER, border = BORDER, seed = SEED,
                                     perforant_path = PP, mossy_fibre = MF, reference_laminae = OTHER,
                                     provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = "Descriptive control geography of the WGCNA modules for the Figure 4 candidate; phenotype-blind.")
message("Control geography: ", nrow(primary$class), " modules; classes ",
        paste(names(table(primary$class$spatial_class)), table(primary$class$spatial_class), collapse = ", "),
        "; territory q < 0.05: ", sum(primary$territory$q_territory < 0.05))
