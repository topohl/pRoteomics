#!/usr/bin/env Rscript
#
# Cross-hemisphere validation of empirical compartment identity.
#
# Rebuilds the empirical compartment contrast from the TRUE hemisphere-resolved
# matrices (R/data_contracts/spatial_systems_data_utils.R), then runs it twice:
#   LEFT discovery  -> RIGHT evaluation
#   RIGHT discovery -> LEFT evaluation
#
# WHY THIS IS NOT A PHENOTYPE ANALYSIS
#   The canonical design is reused verbatim: ~ AnimalID + dataset. Each animal
#   belongs to exactly one StressGroup, so AnimalID as a fixed effect FULLY
#   ABSORBS group membership. The dataset (compartment) coefficients are
#   therefore estimated within animal and are independent of phenotype by
#   construction - no group term is needed or wanted.
#
# WHY "CROSS-HEMISPHERE VALIDATION" AND NOT "REPLICATION"
#   The same nine animals contribute both sides. This measures within-animal
#   reproducibility of a compartment effect, not independent replication.
#
# MARKER-SET SIZES ARE PRE-SPECIFIED (25 / 50 / 100) and are NOT chosen by
# whichever value produces the best agreement.
#
# USAGE
#   Rscript analysis/03_spatial_validation/03_bilateral_empirical_compartment.R
#   Rscript analysis/03_spatial_validation/03_bilateral_empirical_compartment.R --dry-run

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/qc/empirical_roi_marker_utils.R")
source("R/data_contracts/spatial_systems_data_utils.R")
source("R/spatial/spatial_systems_bilateral_utils.R")
source("R/data_contracts/spatial_systems_endpoint_utils.R")

suppressPackageStartupMessages({ library(readr); library(dplyr); library(tidyr) })

SCRIPT_ID <- "analysis/03_spatial_validation/03_bilateral_empirical_compartment.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

TRANSFER_SIZES <- c(25L, 50L, 100L)   # pre-specified, not tuned

OUT <- function() {
  d <- path_results("tables", "11_spatial_systems", "bilateral"); dir_create(d); d
}
crosswalk_path <- function() {
  path_results("tables", "03_qc_exploration", "05_empirical_roi_marker_discovery",
               "empirical_roi_protein_identity_crosswalk_proposed.csv")
}
marker_sets_path <- function() {
  path_results("tables", "03_qc_exploration", "05_empirical_roi_marker_discovery",
               "empirical_roi_marker_sets.csv")
}
DATASETS <- valid_datasets()

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Cross-hemisphere empirical compartment validation.\n")
  dry_run_inputs(SCRIPT_ID, list(protein_identity_crosswalk = crosswalk_path(),
                                 empirical_marker_sets = marker_sets_path()))
  cat("[DRY-RUN] Design reused verbatim: ~ AnimalID + dataset (phenotype absorbed).\n")
  cat("[DRY-RUN] Pre-specified transfer sizes: ",
      paste(TRANSFER_SIZES, collapse = ", "), "\n", sep = "")
  quit(save = "no", status = 0L)
}
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("missing_required_input: limma is required.", call. = FALSE)
}

# ----------------------------------------------- cross-dataset alignment

cw <- as.data.frame(readr::read_csv(crosswalk_path(), show_col_types = FALSE,
                                    progress = FALSE, guess_max = Inf))
cw <- cw[cw$cross_dataset_model_eligible %in% TRUE, , drop = FALSE]
key_by_ds <- split(cw, cw$dataset)
common_keys <- Reduce(intersect, lapply(key_by_ds, function(z) unique(z$empirical_protein_group_key)))
message(sprintf("Cross-dataset aligned protein groups: %d", length(common_keys)))
if (length(common_keys) < 100L) {
  stop("Too few cross-dataset aligned protein groups to model compartment identity.",
       call. = FALSE)
}

# Build, per side, one value per aligned protein x AnimalID x Region x dataset.
build_side_observations <- function(side) {
  mats <- list(); metas <- list()
  for (ds in DATASETS) {
    lv <- sps_levels_for_dataset(ds)
    m1 <- lv$level1$mat
    meta <- lv$level1$meta
    sel <- meta$Hemisphere == side
    # SIDE PURITY: only this hemisphere's columns may enter.
    stopifnot(all(meta$Hemisphere[sel] == side))
    sub <- m1[, sel, drop = FALSE]
    sm <- meta[sel, , drop = FALSE]
    # collapse neuropil layers to Region so all three datasets share the unit
    cellkey <- paste(sm$AnimalID, sm$Region, sep = "\037")
    agg <- sps_aggregate_columns(sub, cellkey)
    parts <- do.call(rbind, strsplit(colnames(agg$mat), "\037", fixed = TRUE))

    # map this dataset's ProteinGroupIDs onto the shared alignment key
    k <- key_by_ds[[ds]]
    lut <- stats::setNames(k$empirical_protein_group_key, k$canonical_ProteinGroupID)
    keys <- unname(lut[rownames(agg$mat)])
    ok <- !is.na(keys) & keys %in% common_keys & !duplicated(keys)
    x <- agg$mat[ok, , drop = FALSE]
    rownames(x) <- keys[ok]
    x <- x[common_keys[common_keys %in% rownames(x)], , drop = FALSE]

    colnames(x) <- paste(ds, parts[, 1], parts[, 2], sep = "::")
    mats[[ds]] <- x
    metas[[ds]] <- data.frame(observation_id = colnames(x), dataset = ds,
                              AnimalID = parts[, 1], Region = parts[, 2],
                              stringsAsFactors = FALSE)
  }
  shared <- Reduce(intersect, lapply(mats, rownames))
  expr <- do.call(cbind, lapply(mats, function(z) z[shared, , drop = FALSE]))
  meta <- dplyr::bind_rows(metas)
  meta <- meta[match(colnames(expr), meta$observation_id), , drop = FALSE]
  list(expr = expr, meta = meta)
}

# Fit the CANONICAL compartment design on one side only.
fit_compartment_side <- function(obs, side) {
  meta <- obs$meta
  meta$AnimalID <- factor(meta$AnimalID,
                          levels = sort(unique(meta$AnimalID), method = "radix"))
  meta$dataset <- factor(meta$dataset, levels = DATASETS)
  design <- stats::model.matrix(~ AnimalID + dataset, data = meta)
  if (qr(design)$rank != ncol(design)) {
    stop("AnimalID + dataset design is not full rank for side ", side, ".",
         call. = FALSE)
  }
  fit <- limma::eBayes(limma::lmFit(obs$expr, design), robust = TRUE, trend = TRUE)
  coefs <- grep("^dataset", colnames(design), value = TRUE)
  out <- list()
  for (cf in coefs) {
    tt <- limma::topTable(fit, coef = cf, number = Inf, sort.by = "none")
    out[[cf]] <- data.frame(
      compartment_contrast = sub("^dataset", "", cf),
      side = side,
      empirical_protein_group_key = rownames(tt),
      effect = tt$logFC, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
      n_observations = ncol(obs$expr),
      n_animals = nlevels(meta$AnimalID),
      model_description = "limma eBayes, ~ AnimalID + dataset, one hemisphere only; AnimalID absorbs StressGroup",
      stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(out)
}

sides <- list(L = build_side_observations("L"), R = build_side_observations("R"))
for (s in names(sides)) {
  message(sprintf("side %s: %d aligned proteins x %d observations (%d animals)",
                  s, nrow(sides[[s]]$expr), ncol(sides[[s]]$expr),
                  length(unique(sides[[s]]$meta$AnimalID))))
}
fits <- dplyr::bind_rows(lapply(names(sides),
                                function(s) fit_compartment_side(sides[[s]], s)))

# ---------------------------------------------------------- protein level

protein_level <- fits %>%
  dplyr::select("compartment_contrast", "empirical_protein_group_key", "side",
                "effect", "adj.P.Val") %>%
  tidyr::pivot_wider(names_from = "side", values_from = c("effect", "adj.P.Val")) %>%
  dplyr::rename(effect_L = "effect_L", effect_R = "effect_R") %>%
  dplyr::mutate(sign_agreement = sign(.data$effect_L) == sign(.data$effect_R),
                abs_L_minus_R = abs(.data$effect_L - .data$effect_R))

summary_tbl <- protein_level %>%
  dplyr::group_by(.data$compartment_contrast) %>%
  dplyr::group_modify(~ {
    ag <- sps_paired_agreement(.x$effect_L, .x$effect_R)
    ag$rank_agreement_spearman <- sps_safe_cor(rank(.x$effect_L), rank(.x$effect_R),
                                               "pearson")
    ag$n_evaluable_proteins <- ag$n_pairs
    ag
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    independence_note = paste0(
      "Cross-hemisphere validation, NOT independent replication: the same ",
      "animals contribute both sides. AnimalID is in the design, so the ",
      "compartment effect is estimated within animal and is phenotype-independent."))

# ------------------------------------------------------- marker transfer

transfer_rows <- list()
for (cc in unique(protein_level$compartment_contrast)) {
  z <- protein_level[protein_level$compartment_contrast == cc, , drop = FALSE]
  for (n in TRANSFER_SIZES) {
    if (nrow(z) < n) next
    for (dir in c("L_to_R", "R_to_L")) {
      disc_col <- if (dir == "L_to_R") "effect_L" else "effect_R"
      eval_col <- if (dir == "L_to_R") "effect_R" else "effect_L"
      # DISCOVERY RANKING USES THE DISCOVERY SIDE ONLY.
      ord <- order(-z[[disc_col]])
      discovered <- z$empirical_protein_group_key[utils::head(ord, n)]
      tr <- sps_rank_transfer(discovered, z[[eval_col]],
                              z$empirical_protein_group_key,
                              higher_is_better = TRUE)
      transfer_rows[[length(transfer_rows) + 1L]] <- cbind(
        data.frame(compartment_contrast = cc, direction = dir, top_n = n,
                   discovery_side = if (dir == "L_to_R") "L" else "R",
                   evaluation_side = if (dir == "L_to_R") "R" else "L",
                   stringsAsFactors = FALSE), tr)
    }
  }
}

# also evaluate the EXISTING canonical bilateral marker sets on each side
if (file.exists(marker_sets_path())) {
  ms <- as.data.frame(readr::read_csv(marker_sets_path(), show_col_types = FALSE,
                                      progress = FALSE, guess_max = Inf))
  if ("empirical_protein_group_key" %in% names(ms)) {
    for (set in sort(unique(ms$marker_set))) {
      keys <- unique(ms$empirical_protein_group_key[ms$marker_set == set])
      for (cc in unique(protein_level$compartment_contrast)) {
        z <- protein_level[protein_level$compartment_contrast == cc, , drop = FALSE]
        for (side in c("L", "R")) {
          col <- if (side == "L") "effect_L" else "effect_R"
          tr <- sps_rank_transfer(keys, z[[col]], z$empirical_protein_group_key,
                                  higher_is_better = TRUE)
          transfer_rows[[length(transfer_rows) + 1L]] <- cbind(
            data.frame(compartment_contrast = cc,
                       direction = paste0("canonical_bilateral_set:", set),
                       top_n = NA_integer_, discovery_side = "bilateral_canonical",
                       evaluation_side = side, stringsAsFactors = FALSE), tr)
        }
      }
    }
  }
}
transfer <- dplyr::bind_rows(transfer_rows)

root <- OUT()
write_csv_safe(protein_level, file.path(root, "bilateral_empirical_compartment_protein_level.csv"))
write_csv_safe(summary_tbl, file.path(root, "bilateral_empirical_compartment_summary.csv"))
write_csv_safe(transfer, file.path(root, "bilateral_empirical_marker_transfer.csv"))

cat("\n===== Cross-hemisphere empirical compartment validation =====\n")
cat(sprintf("aligned proteins: %d   design: ~ AnimalID + dataset (phenotype absorbed)\n",
            nrow(sides$L$expr)))
cat("\n--- effect concordance ---\n")
cat(sprintf("  %-34s %6s %7s %7s %8s\n", "compartment contrast", "n", "r", "rho", "signAgr"))
for (i in seq_len(nrow(summary_tbl))) {
  cat(sprintf("  %-34s %6d %7.3f %7.3f %8.3f\n",
              substr(summary_tbl$compartment_contrast[i], 1, 34),
              summary_tbl$n_evaluable_proteins[i], summary_tbl$pearson_r[i],
              summary_tbl$spearman_rho[i], summary_tbl$sign_agreement_fraction[i]))
}
cat("\n--- marker transfer (pre-specified sizes) ---\n")
tt <- transfer[transfer$discovery_side %in% c("L", "R"), , drop = FALSE]
for (i in seq_len(nrow(tt))) {
  cat(sprintf("  %-34s %-8s top%-4d  medianRankFrac=%.3f  AUC=%.3f\n",
              substr(tt$compartment_contrast[i], 1, 34), tt$direction[i],
              tt$top_n[i], tt$median_eval_rank_fraction[i], tt$auc_like[i]))
}
cat("\nOutputs:", relative_to(root), "\n")
cat("Cross-hemisphere validation; the same animals contribute both sides.\n")
