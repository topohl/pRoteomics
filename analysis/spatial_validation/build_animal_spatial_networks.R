#!/usr/bin/env Rscript
#
#
# One AnimalID = one independent network replicate. For every animal this builds
# a Left network, a Right network and an equal-weight bilateral network from the
# validated hemisphere-resolved constructor.
#
# THE REPRESENTATION IS CHOSEN PHENOTYPE-BLIND
#   Four transformations are declared prospectively. The winner is selected on
#   CON animals using bilateral reproducibility, inter-animal stability,
#   anatomical plausibility and robustness. RES/SUS separation is never a
#   selection criterion and is not computed during selection.
#
# NODE = an anatomical sampling unit. EDGE = similarity of the molecular spatial
# profiles of two anatomical units within one animal. This is NOT connectivity
# of any kind and NOT a coexpression network.
#
# USAGE
#   Rscript analysis/spatial_validation/build_animal_spatial_networks.R
#   Rscript analysis/spatial_validation/build_animal_spatial_networks.R --dry-run
# Script: analysis/spatial_validation/build_animal_spatial_networks.R
# Stage: networks
# Scope: global
# Consumes: required data/processed/01_preprocessing/06_merged_metadata_module_score/<dataset>/sample_metadata_merged_clean_for_module_scores.xlsx; optional results/spatial_validation/build_spatial_data_contract/global/tables/spatial_systems_hemisphere_inventory.csv; results/tables/11_spatial_systems/data_contract/spatial_systems_hemisphere_inventory.csv
# Produces: results/spatial_validation/build_animal_spatial_networks/global/tables/network_representation_selection.csv; results/spatial_validation/build_animal_spatial_networks/global/tables/animal_network_edges.csv; results/spatial_validation/build_animal_spatial_networks/global/tables/animal_network_nodes.csv; +4 more
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Animal-level spatial molecular-similarity networks.
#  
#  

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/data_contracts/spatial_systems_data_utils.R")
source("R/spatial/spatial_systems_bilateral_utils.R")
source("R/spatial/spatial_atlas_utils.R")
source("R/networks/animal_spatial_network_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))

# Phase 6G.3: destinations resolve through the normalized output contract,
# addressed by this analysis's own identity rather than by the historical
# 11_spatial_systems stage directory. Outputs already written there stay
# exactly where they are and are read, never rewritten.
ANALYSIS_ID <- "build_animal_spatial_networks"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "analysis/spatial_validation/build_animal_spatial_networks.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

EDGE_METHOD <- "spearman"   # prospective choice; audited in the selection table
OUT <- function() {
  d <- CANONICAL_PATHS$tables; dir_create(d); d
}
DATASETS <- valid_datasets()

if (isTRUE(cli$dry_run)) {
  inputs <- list()
  for (ds in DATASETS) {
    inputs[[paste0("canonical_metadata__", ds)]] <-
      preprocessing_module_score_metadata(ds)
  }
  cat("[DRY-RUN] Animal-level spatial molecular-similarity networks.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Representation selected phenotype-blind on CON animals.\n")
  cat("[DRY-RUN] One network per AnimalID; full weighted matrices, no thresholding.\n")
  quit(save = "no", status = 0L)
}

# ------------------------------------------------ per-animal side matrices

message("Building per-animal side-resolved matrices")
animal_matrices <- list()
node_sets <- list()
animal_meta <- list()
for (ds in DATASETS) {
  inputs <- resolve_dataset_inputs(ds, purpose = "wgcna", script = SCRIPT_ID,
                                   stage = "networks")
  md <- preprocessing_module_score_metadata(ds)
  canonical <- qc_load_canonical_expression(inputs$expression_file, md,
                                            dataset = ds, strict = TRUE)
  lv <- sps_build_spatial_levels(ds, canonical = canonical)

  # canonical node vocabulary, fixed and shared by every comparison
  nodes <- sort(unique(sat_canonical_spatial_unit(lv$level2$meta$SpatialUnit, ds)))
  node_sets[[ds]] <- nodes

  l1 <- lv$level1$mat; m1 <- lv$level1$meta
  m1$unit <- sat_canonical_spatial_unit(m1$SpatialUnit, ds)
  l2 <- lv$level2$mat; m2 <- lv$level2$meta
  m2$unit <- sat_canonical_spatial_unit(m2$SpatialUnit, ds)

  for (an in sort(unique(m2$AnimalID))) {
    grp <- m2$StressGroup[m2$AnimalID == an][1]
    pick_side <- function(side) {
      sel <- m1$AnimalID == an & m1$Hemisphere == side
      if (!any(sel)) return(NULL)
      # SIDE PURITY: only this hemisphere's columns may enter
      stopifnot(all(m1$Hemisphere[sel] == side))
      x <- l1[, sel, drop = FALSE]
      colnames(x) <- m1$unit[sel]
      x[, intersect(nodes, colnames(x)), drop = FALSE]
    }
    sel2 <- m2$AnimalID == an
    bil <- l2[, sel2, drop = FALSE]
    colnames(bil) <- m2$unit[sel2]
    bil <- bil[, intersect(nodes, colnames(bil)), drop = FALSE]

    key <- paste(ds, an, sep = "\037")
    animal_matrices[[key]] <- list(L = pick_side("L"), R = pick_side("R"),
                                   bilateral = bil)
    animal_meta[[key]] <- data.frame(dataset = ds, AnimalID = an,
                                     StressGroup = grp, stringsAsFactors = FALSE)
  }
  message(sprintf("  %-16s %d animals, %d nodes", ds,
                  length(unique(m2$AnimalID)), length(nodes)))
}
animal_meta <- dplyr::bind_rows(animal_meta)

# ============================ PART 5: phenotype-blind representation selection

message("Selecting the network representation on CON animals only")
sel_rows <- list()
for (ds in DATASETS) {
  nodes <- node_sets[[ds]]
  con_keys <- names(animal_matrices)[
    grepl(paste0("^", ds, "\037"), names(animal_matrices))]
  con_keys <- con_keys[vapply(con_keys, function(k)
    animal_meta$StressGroup[animal_meta$dataset == ds &
      animal_meta$AnimalID == sub(".*\037", "", k)][1] == "CON", logical(1))]

  for (rep in asn_representations()) {
    bil_r <- c(); stab <- c(); within_gt_across <- c(); retained <- c()
    edge_mats <- list()
    for (k in con_keys) {
      am <- animal_matrices[[k]]
      tb <- asn_transform(am$bilateral, rep)
      if (nrow(tb$mat) < 10L) next
      sb <- asn_similarity_matrix(tb$mat, nodes, EDGE_METHOD)
      edge_mats[[k]] <- asn_edge_vector(sb)
      retained <- c(retained, tb$n_proteins_retained / tb$n_proteins_input)
      # 1. bilateral reproducibility of the network itself
      if (!is.null(am$L) && !is.null(am$R)) {
        tl <- asn_transform(am$L, rep); tr <- asn_transform(am$R, rep)
        common <- intersect(rownames(tl$mat), rownames(tr$mat))
        if (length(common) >= 10L) {
          el <- asn_edge_vector(asn_similarity_matrix(tl$mat[common, , drop = FALSE], nodes, EDGE_METHOD))
          er <- asn_edge_vector(asn_similarity_matrix(tr$mat[common, , drop = FALSE], nodes, EDGE_METHOD))
          bil_r <- c(bil_r, sps_safe_cor(el, er, "pearson"))
        }
      }
      # 3. anatomical plausibility: within-region layer pairs should be more
      #    similar than arbitrary across-region pairs (neuropil only)
      if (identical(ds, "neuron_neuropil")) {
        v <- edge_mats[[k]]
        parts <- do.call(rbind, strsplit(names(v), "__", fixed = TRUE))
        same_region <- sub("_.*$", "", parts[, 1]) == sub("_.*$", "", parts[, 2])
        within_gt_across <- c(within_gt_across,
                              mean(v[same_region], na.rm = TRUE) -
                                mean(v[!same_region], na.rm = TRUE))
      }
    }
    # 2. inter-animal stability among CON
    if (length(edge_mats) >= 2L) {
      em <- do.call(rbind, edge_mats)
      cm <- suppressWarnings(stats::cor(t(em), use = "pairwise.complete.obs"))
      stab <- cm[upper.tri(cm)]
    }
    sel_rows[[length(sel_rows) + 1L]] <- data.frame(
      dataset = ds, representation = rep,
      description = asn_representation_description(rep),
      n_con_animals = length(edge_mats),
      median_bilateral_edge_r = if (length(bil_r)) stats::median(bil_r, na.rm = TRUE) else NA_real_,
      median_inter_animal_stability = if (length(stab)) stats::median(stab, na.rm = TRUE) else NA_real_,
      anatomical_within_minus_across = if (length(within_gt_across))
        stats::median(within_gt_across, na.rm = TRUE) else NA_real_,
      median_fraction_proteins_retained = if (length(retained)) stats::median(retained) else NA_real_,
      selection_criteria = "CON only; bilateral reproducibility, inter-animal stability, anatomical plausibility, protein retention. No phenotype contrast is computed here.",
      stringsAsFactors = FALSE)
  }
}
selection <- dplyr::bind_rows(sel_rows)

# SELECTION RULE. Ordered as validity, then reliability, then stability.
#
#   1. ANATOMICAL VALIDITY FLOOR. A representation must actually separate
#      within-region layer pairs from across-region pairs. This is a validity
#      requirement, not a tiebreaker: a transformation that cannot distinguish
#      anatomy is not measuring spatial organisation at all, whatever else it
#      scores. Raw abundance fails here (0.009 versus 0.12-0.13 for the others)
#      because similarity is dominated by the global abundance rank that every
#      unit shares - which is exactly why it also posts the HIGHEST inter-animal
#      stability. Ranking on stability first would reward that degeneracy.
#   2. RELIABILITY. Among valid representations, maximise the median bilateral
#      edge correlation - the one criterion that measures whether the network
#      itself reproduces within an animal.
#   3. Ties by inter-animal stability, then by interpretability.
#
# Every criterion is computed on CON animals only. No phenotype contrast is
# calculated anywhere in this block, so the choice cannot be steered by group
# separation.
ANATOMICAL_FLOOR <- 0.05
selection$meets_bilateral_floor <- selection$median_bilateral_edge_r >= 0.5
selection$meets_retention_floor <- selection$median_fraction_proteins_retained >= 0.5
selection$meets_anatomical_floor <- is.na(selection$anatomical_within_minus_across) |
  selection$anatomical_within_minus_across >= ANATOMICAL_FLOOR
pick_for <- function(ds) {
  z <- selection[selection$dataset == ds, , drop = FALSE]
  ok <- z[z$meets_bilateral_floor %in% TRUE & z$meets_retention_floor %in% TRUE &
            z$meets_anatomical_floor %in% TRUE, , drop = FALSE]
  if (!nrow(ok)) ok <- z
  ok <- ok[order(-ok$median_bilateral_edge_r,
                 -ok$median_inter_animal_stability), , drop = FALSE]
  ok$representation[1]
}
chosen <- vapply(DATASETS, pick_for, character(1))
# One representation is frozen across datasets for interpretability; the
# neuropil choice governs because it is the only dataset with a layer axis.
PRIMARY <- unname(chosen[["neuron_neuropil"]])
selection$is_primary <- selection$representation == PRIMARY
selection$selection_rule <- paste0(
  "CON only. VALIDITY FLOORS (all three required): bilateral edge r >= 0.5, ",
  "protein retention >= 0.5, and anatomical discrimination ",
  "(within-region minus across-region edge similarity) >= ",
  format(ANATOMICAL_FLOOR), ". Among representations that clear all three, ",
  "maximise the median bilateral edge correlation; ties by median ",
  "inter-animal stability. Anatomical discrimination is a FLOOR, not a ",
  "tiebreaker: raw abundance posts the highest inter-animal stability ",
  "precisely because similarity is dominated by the global abundance rank ",
  "every unit shares, so ranking on stability first would reward that ",
  "degeneracy. No phenotype contrast enters this rule.")
selection$edge_metric <- EDGE_METHOD

message("Primary representation: ", PRIMARY)

# ================================= PART 6: one network per animal

message("Building one network per AnimalID")
edge_rows <- list(); node_rows <- list(); matrices <- list()
for (k in names(animal_matrices)) {
  ds <- sub("\037.*$", "", k); an <- sub(".*\037", "", k)
  grp <- animal_meta$StressGroup[animal_meta$dataset == ds & animal_meta$AnimalID == an][1]
  nodes <- node_sets[[ds]]
  am <- animal_matrices[[k]]
  for (mode in c("L", "R", "bilateral")) {
    x <- am[[mode]]
    if (is.null(x) || !ncol(x)) next
    tt <- asn_transform(x, PRIMARY)
    if (nrow(tt$mat) < 10L) next
    if (!all(nodes %in% colnames(tt$mat))) next
    s <- asn_similarity_matrix(tt$mat, nodes, EDGE_METHOD)
    matrices[[paste(k, mode, sep = "\037")]] <- s
    et <- asn_edge_table(s)
    edge_rows[[length(edge_rows) + 1L]] <- cbind(
      data.frame(dataset = ds, AnimalID = an, StressGroup = grp,
                 hemisphere_mode = mode, representation = PRIMARY,
                 edge_metric = EDGE_METHOD,
                 n_proteins_contributing = tt$n_proteins_retained,
                 stringsAsFactors = FALSE), et)
    ns <- asn_node_strength(s)
    node_rows[[length(node_rows) + 1L]] <- data.frame(
      dataset = ds, AnimalID = an, StressGroup = grp, hemisphere_mode = mode,
      node = names(ns), node_strength = unname(ns),
      n_proteins_contributing = tt$n_proteins_retained, stringsAsFactors = FALSE)
  }
}
edges <- dplyr::bind_rows(edge_rows)
nodes_tbl <- dplyr::bind_rows(node_rows)

# ======================= PART 8 + 9: bilateral network reproducibility

message("Assessing bilateral network reproducibility")
bil_rows <- list(); edge_bil_rows <- list()
for (ds in DATASETS) {
  nodes <- node_sets[[ds]]
  for (an in sort(unique(animal_meta$AnimalID[animal_meta$dataset == ds]))) {
    kl <- paste(ds, an, "L", sep = "\037"); kr <- paste(ds, an, "R", sep = "\037")
    if (is.null(matrices[[kl]]) || is.null(matrices[[kr]])) next
    el <- asn_edge_vector(matrices[[kl]]); er <- asn_edge_vector(matrices[[kr]])
    sl <- asn_node_strength(matrices[[kl]]); sr <- asn_node_strength(matrices[[kr]])
    grp <- animal_meta$StressGroup[animal_meta$dataset == ds & animal_meta$AnimalID == an][1]
    ag <- sps_paired_agreement(el, er)
    bil_rows[[length(bil_rows) + 1L]] <- cbind(
      data.frame(dataset = ds, AnimalID = an, StressGroup = grp,
                 n_edges = length(el), stringsAsFactors = FALSE), ag,
      data.frame(node_strength_pearson = sps_safe_cor(sl, sr, "pearson"),
                 node_strength_spearman = sps_safe_cor(sl, sr, "spearman"),
                 stringsAsFactors = FALSE))
    edge_bil_rows[[length(edge_bil_rows) + 1L]] <- data.frame(
      dataset = ds, AnimalID = an, StressGroup = grp, edge_id = names(el),
      left = unname(el), right = unname(er), abs_difference = abs(el - er),
      stringsAsFactors = FALSE)
  }
}
bilateral <- dplyr::bind_rows(bil_rows)

# PRESPECIFIED quality bands, fixed before inspection and never phenotype-tuned.
bilateral$bilateral_network_quality <- dplyr::case_when(
  !is.finite(bilateral$pearson_r) ~ "insufficient_data",
  bilateral$pearson_r >= 0.70 ~ "high_bilateral_reproducibility",
  bilateral$pearson_r >= 0.40 ~ "moderate_bilateral_reproducibility",
  TRUE ~ "low_bilateral_reproducibility")
bilateral$retention_policy <- paste0(
  "An animal is RETAINED regardless of its bilateral network quality. Exclusion ",
  "would require an independent QC criterion; a weak group result is not one.")

edge_bilateral <- dplyr::bind_rows(edge_bil_rows) %>%
  dplyr::group_by(.data$dataset, .data$edge_id) %>%
  dplyr::summarise(
    n_animals = dplyr::n(),
    left_right_pearson = sps_safe_cor(.data$left, .data$right, "pearson"),
    left_right_spearman = sps_safe_cor(.data$left, .data$right, "spearman"),
    median_abs_LR_difference = stats::median(.data$abs_difference, na.rm = TRUE),
    mean_left = mean(.data$left, na.rm = TRUE),
    mean_right = mean(.data$right, na.rm = TRUE),
    .groups = "drop")

# ------------------------------------------------------------------ write

root <- OUT()
write_csv_safe(selection, file.path(root, "network_representation_selection.csv"))
write_csv_safe(edges, file.path(root, "animal_network_edges.csv"))
write_csv_safe(nodes_tbl, file.path(root, "animal_network_nodes.csv"))
write_csv_safe(bilateral, file.path(root, "animal_network_bilateral_reproducibility.csv"))
write_csv_safe(edge_bilateral, file.path(root, "edge_bilateral_reproducibility.csv"))
# full weighted matrices, retained independently of any visualisation threshold
mat_dir <- file.path(root, "matrices"); dir_create(mat_dir)
for (k in names(matrices)) {
  parts <- strsplit(k, "\037", fixed = TRUE)[[1]]
  f <- file.path(mat_dir, sprintf("%s__%s__%s.csv", parts[1], parts[2], parts[3]))
  m <- matrices[[k]]
  write_csv_safe(cbind(data.frame(node = rownames(m), stringsAsFactors = FALSE),
                       as.data.frame(m, check.names = FALSE)), f)
}
saveRDS(list(matrices = matrices, nodes = node_sets, representation = PRIMARY,
             edge_metric = EDGE_METHOD, contract_version = asn_contract_version()),
        file.path(root, "animal_network_objects.rds"))

cat("\n===== Animal-level spatial molecular-similarity networks =====\n")
cat("primary representation:", PRIMARY, " edge metric:", EDGE_METHOD, "\n\n")
cat("--- representation selection (CON only, phenotype-blind) ---\n")
cat(sprintf("  %-16s %-22s %8s %9s %9s %8s\n", "dataset", "representation",
            "bilat_r", "stability", "anat", "retain"))
for (i in seq_len(nrow(selection))) {
  cat(sprintf("  %-16s %-22s %8.3f %9.3f %9.3f %8.2f%s\n",
              selection$dataset[i], selection$representation[i],
              selection$median_bilateral_edge_r[i],
              selection$median_inter_animal_stability[i],
              ifelse(is.na(selection$anatomical_within_minus_across[i]), NA_real_,
                     selection$anatomical_within_minus_across[i]),
              selection$median_fraction_proteins_retained[i],
              ifelse(selection$is_primary[i], "  <= PRIMARY", "")))
}
cat("\n--- bilateral network reproducibility per animal ---\n")
for (ds in DATASETS) {
  z <- bilateral[bilateral$dataset == ds, , drop = FALSE]
  if (!nrow(z)) next
  cat(sprintf("\n  %s\n", ds))
  for (i in seq_len(nrow(z))) {
    cat(sprintf("    %-6s %-4s edge_r=%6.3f rho=%6.3f RMSE=%6.3f nodeStr_r=%6.3f  %s\n",
                z$AnimalID[i], z$StressGroup[i], z$pearson_r[i], z$spearman_rho[i],
                z$RMSE[i], z$node_strength_pearson[i],
                sub("_bilateral_reproducibility", "", z$bilateral_network_quality[i])))
  }
}
cat("\nquality classes:\n"); print(table(bilateral$bilateral_network_quality))
cat("\nnetworks built:", length(matrices), " edges:", nrow(edges), "\n")
cat("Outputs:", relative_to(root), "\n")
cat("One AnimalID = one network replicate. Full weighted matrices retained.\n")
