#!/usr/bin/env Rscript
#
#
# Unit specificity of the stress effects: one model over all units of a
# dataset with the animal as a blocking factor, consensus effects across
# units and direct unit-by-group contrasts (added 2026-10-29).
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/test_unit_specificity.R
#   Rscript analysis/spatial_validation/test_unit_specificity.R --dry-run
# Script: analysis/spatial_validation/test_unit_specificity.R
# Stage: networks
# Scope: global
# Consumes: required data/processed/01_preprocessing/protigy_input_animal_level/<dataset>/<dataset>_animal_level.gct; results/preprocessing/map_protein_identifiers/<dataset>/tables/mapped/forward/per_file/<unit>sus_<unit>res.csv
# Produces: results/spatial_validation/test_unit_specificity/global/tables/unit_specificity_contrasts.csv; results/spatial_validation/test_unit_specificity/global/tables/unit_specificity_summary.csv; +1 more
# Dataset behavior: global; reads the animal-level matrices of neuron_neuropil, neuron_soma and microglia.
# Notes: Direct unit-by-group contrasts (is CA2 SLM's SUS-RES effect larger than in other laminae?) and consensus effects, limma with the animal blocked (duplicateCorrelation); protein-level.
#
# WHY: Figure 3a finds most SUS-RES proteins in CA2 SLM. Whether an effect is
# specific to a unit can only be shown by a direct contrast between units, not
# by the unit where it reaches significance.
#
# MODEL: per dataset the animal-level matrices of the Figure 3 model
# (fit_animal_level_limma.R; hemispheres averaged, imputed), ~0 + unit x group
# cell means, the animal as a block (limma duplicateCorrelation, consensus
# correlation), eBayes(robust = TRUE, trend = FALSE), BH within each contrast.
# Protein-level inference, labelled as such.
#
# CONTRASTS (declared):
#   per unit       SUS - RES and stressed (RES + SUS)/2 - CON, in the blocked model
#   consensus      the mean over a dataset's units of SUS - RES, and of stressed - CON
#   neuropil       CA2 SLM SUS - RES minus the same in CA1 SLM (same stratum, other region),
#                  minus CA2 SR (same region, other stratum), and minus the mean of the other
#                  nine laminae; the same three for stressed - CON
#   prior          CA3 SR minus CA3 SO, SUS - RES and stressed - CON (apical dendrite
#                  atrophy of CA3 under chronic stress)
#   microglia ROI  CA2 minus the mean of the other three ROIs, SUS - RES and stressed - CON
# The CA2 SLM proteins of Figure 3a (official SUS - RES BH < 0.05) get their own
# table with the three CA2 SLM specificity contrasts, BH over those proteins.
#
#

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))
source("R/spatial/neuropil_lamina_utils.R")
source("R/spatial/stress_effect_utils.R")
suppressPackageStartupMessages(library(limma))

ANALYSIS_ID <- "test_unit_specificity"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

SCRIPT_ID <- "analysis/spatial_validation/test_unit_specificity.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- list(
  contrasts = file.path(CANONICAL_PATHS$tables, "unit_specificity_contrasts.csv"),
  summary   = file.path(CANONICAL_PATHS$tables, "unit_specificity_summary.csv"),
  ca2_slm   = file.path(CANONICAL_PATHS$tables, "unit_specificity_ca2_slm_proteins.csv")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")
LAMINAE <- c("CA1so", "CA1sr", "CA1slm", "CA2so", "CA2sr", "CA2slm", "CA3so", "CA3sr", "DGmo", "DGpo")
ROIS <- c("CA1microglia", "CA2microglia", "CA3microglia", "DGmicroglia")

gct_paths <- vapply(STRESS_DATASETS, animal_level_gct_path, "")
if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Unit specificity of the stress effects (one blocked model per dataset).\n")
  dry_run_inputs(SCRIPT_ID, as.list(gct_paths))
  quit(save = "no", status = 0L)
}
for (f in gct_paths) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))

# a contrast vector over the cell means: weights on unit x group cells
cw <- function(cols, weights) {
  w <- setNames(numeric(length(cols)), cols)
  for (k in names(weights)) {
    if (!k %in% cols) stop("Cell ", k, " not in the design.", call. = FALSE)
    w[[k]] <- w[[k]] + weights[[k]]
  }
  w
}
sus_res <- function(u, f = 1) setNames(c(f, -f), paste0(u, c("_SUS", "_RES")))
stressed_con <- function(u, f = 1) setNames(c(f / 2, f / 2, -f), paste0(u, c("_RES", "_SUS", "_CON")))
combine <- function(...) { l <- list(...); n <- unique(unlist(lapply(l, names))); setNames(vapply(n, function(k) sum(unlist(lapply(l, function(x) if (k %in% names(x)) x[[k]] else 0))), 0), n) }

rows <- list(); summ <- list(); mapped_inputs <- character()
for (ds in STRESS_DATASETS) {
  g <- read_animal_level_gct(ds)
  units <- unique(g$unit)
  cell <- factor(g$cell, levels = unique(g$cell))
  design <- model.matrix(~ 0 + cell); colnames(design) <- levels(cell)
  corfit <- duplicateCorrelation(g$m, design, block = g$animal)
  fit <- lmFit(g$m, design, block = g$animal, correlation = corfit$consensus)
  spec <- list()
  for (u in units) {
    spec[[paste0("per_unit|", u, "|SUS_minus_RES")]] <- sus_res(u)
    spec[[paste0("per_unit|", u, "|stressed_minus_CON")]] <- stressed_con(u)
  }
  spec[["consensus|all units|SUS_minus_RES"]] <- do.call(combine, lapply(units, sus_res, f = 1 / length(units)))
  spec[["consensus|all units|stressed_minus_CON"]] <- do.call(combine, lapply(units, stressed_con, f = 1 / length(units)))
  if (ds == "neuron_neuropil") {
    others <- setdiff(LAMINAE, "CA2slm")
    for (eff in c("SUS_minus_RES", "stressed_minus_CON")) {
      fx <- if (eff == "SUS_minus_RES") sus_res else stressed_con
      spec[[paste0("specificity|CA2slm vs CA1slm|", eff)]] <- combine(fx("CA2slm"), fx("CA1slm", -1))
      spec[[paste0("specificity|CA2slm vs CA2sr|", eff)]] <- combine(fx("CA2slm"), fx("CA2sr", -1))
      spec[[paste0("specificity|CA2slm vs other laminae|", eff)]] <- do.call(combine, c(list(fx("CA2slm")), lapply(others, fx, f = -1 / length(others))))
      spec[[paste0("prior|CA3sr vs CA3so|", eff)]] <- combine(fx("CA3sr"), fx("CA3so", -1))
    }
  }
  if (ds == "microglia") {
    for (eff in c("SUS_minus_RES", "stressed_minus_CON")) {
      fx <- if (eff == "SUS_minus_RES") sus_res else stressed_con
      spec[[paste0("specificity|CA2microglia vs other ROIs|", eff)]] <- do.call(combine, c(list(fx("CA2microglia")), lapply(setdiff(ROIS, "CA2microglia"), fx, f = -1 / 3)))
    }
  }
  cm <- sapply(spec, function(w) cw(colnames(design), w))
  eb <- eBayes(contrasts.fit(fit, cm), robust = TRUE, trend = FALSE)
  ident <- canonical_protein_identity(ds, units[1], nrow(g$m))
  mapped_inputs <- c(mapped_inputs, mapped_unit_table_path(ds, units[1]))
  for (k in names(spec)) {
    tt <- topTable(eb, coef = k, number = Inf, sort.by = "none", adjust.method = "BH")
    parts <- strsplit(k, "|", fixed = TRUE)[[1]]
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = ds, kind = parts[1], comparison = parts[2], effect = parts[3], ProteinGroupID = ident$ProteinGroupID,
      official_gene_symbol = ident$official_gene_symbol, gene_level_claim_allowed = ident$gene_level_claim_allowed,
      estimate = tt$logFC, t = tt$t, p_value = tt$P.Value, p_adjust = tt$adj.P.Val, stringsAsFactors = FALSE)
    summ[[length(summ) + 1L]] <- data.frame(
      dataset = ds, kind = parts[1], comparison = parts[2], effect = parts[3], n_proteins = nrow(tt),
      n_fdr_005 = sum(tt$adj.P.Val < 0.05), n_fdr_010 = sum(tt$adj.P.Val < 0.1), share_p_below_005 = mean(tt$P.Value < 0.05),
      prop_true_null = limma::propTrueNull(tt$P.Value), animal_block_correlation = corfit$consensus, stringsAsFactors = FALSE)
  }
}
contrasts <- do.call(rbind, rows)
summary_tab <- do.call(rbind, summ)
rownames(contrasts) <- rownames(summary_tab) <- NULL

# the Figure 3a CA2 SLM proteins: their per-unit effect in the blocked model and the three specificity contrasts
official <- utils::read.csv(mapped_unit_table_path("neuron_neuropil", "CA2slm"), stringsAsFactors = FALSE)
daps <- official$ProteinGroupID[official$padj < 0.05]
pick <- function(comparison, effect) {
  x <- contrasts[contrasts$dataset == "neuron_neuropil" & contrasts$comparison == comparison & contrasts$effect == effect &
                   contrasts$ProteinGroupID %in% daps, , drop = FALSE]
  x[match(daps, x$ProteinGroupID), , drop = FALSE]
}
ca2 <- data.frame(ProteinGroupID = daps, official_gene_symbol = official$official_gene_symbol[match(daps, official$ProteinGroupID)],
                  official_log2fc = official$log2fc[match(daps, official$ProteinGroupID)],
                  official_p_adjust = official$padj[match(daps, official$ProteinGroupID)], stringsAsFactors = FALSE)
b <- pick("CA2slm", "SUS_minus_RES")
ca2$blocked_log2fc <- b$estimate; ca2$blocked_p_value <- b$p_value
for (cmpn in c("CA2slm vs CA1slm", "CA2slm vs CA2sr", "CA2slm vs other laminae")) {
  x <- pick(cmpn, "SUS_minus_RES")
  key <- gsub("[^a-z0-9]+", "_", tolower(sub("^CA2slm vs ", "", cmpn)))
  ca2[[paste0("vs_", key, "_estimate")]] <- x$estimate
  ca2[[paste0("vs_", key, "_p_value")]] <- x$p_value
  ca2[[paste0("vs_", key, "_p_adjust_within_these_proteins")]] <- stats::p.adjust(x$p_value, method = "BH")
  ca2[[paste0("vs_", key, "_same_direction")]] <- sign(x$estimate) == sign(ca2$official_log2fc)
}

utils::write.csv(contrasts, OUT$contrasts, row.names = FALSE)
utils::write.csv(summary_tab, OUT$summary, row.names = FALSE)
utils::write.csv(ca2, OUT$ca2_slm, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = c(as.list(gct_paths), as.list(stats::setNames(unique(mapped_inputs), basename(unique(mapped_inputs))))),
                   outputs = OUT,
                   parameters = list(model = "per dataset: ~0 + unit x group cell means, animal blocked (duplicateCorrelation consensus), eBayes(robust = TRUE, trend = FALSE), BH within contrast",
                                     contrasts = unique(paste(summary_tab$kind, summary_tab$comparison, summary_tab$effect, sep = " | ")),
                                     ca2_slm_proteins = length(daps), provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = "Protein-level inference (residual variance pooled over units, animal as block).")
message("Unit specificity complete: ", nrow(summary_tab), " contrasts")
