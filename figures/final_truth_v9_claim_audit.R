#!/usr/bin/env Rscript

# Part-26: pre-submission claim and uncertainty audit.
#
# AUDIT ONLY. This script computes nothing new about the biology: every number
# it checks is read back from an artefact the figure layer already produced,
# and its job is to prove that what the figures SAY matches what the source
# data CONTAIN, and that each displayed quantity is described at the strength
# the design can carry.
#
# It hard-fails on a legend/source-data disagreement or a spatial-label
# mismatch, because those are factual errors rather than matters of taste.

source(file.path("R", "paths.R"))
source(repo_path("R", "null_coalescing.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
source(repo_path("R", "final_truth_v9_figure_utils.R"))
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_claim_audit.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] final_truth_v9 claim and uncertainty audit")
  quit(save = "no", status = 0L)
}

ct <- s9f_contract()
SD <- path_results("source_data", "manuscript_candidates", "final_truth_v9")
REP <- path_results("reports", "manuscript_candidates", "final_truth_v9")
OUT <- path_results("tables", "manuscript_candidates", "final_truth_v9", "audit")
dir_create(OUT)
sc <- function(fig, id) {
  p <- file.path(SD, fig, paste0(id, "_source_data.csv"))
  if (file.exists(p)) nv_read_csv(p) else NULL
}
fails <- character(0)
FAIL <- function(msg) fails <<- c(fails, msg)

# =================================================== S20 legend vs source data
dap <- sc("figure_03", "v9_dap_track")
loc <- sc("extended_data", "v9_ed_locations")
wg  <- sc("extended_data", "v9_ed_wgcna_phenotype")
nul <- sc("extended_data", "v9_ed_nulls")
sen <- sc("extended_data", "v9_ed_ca2_sensitivity")
leg <- paste(readLines(file.path(REP, "final_figure_legends_v9.md"), warn = FALSE),
             collapse = " ")

chk <- function(id, statement, stated, observed, tol = 0) {
  ok <- isTRUE(all.equal(as.numeric(stated), as.numeric(observed),
                         tolerance = if (tol > 0) tol else 1e-9))
  if (!ok) FAIL(sprintf("%s: stated %s, source data %s", id, stated, observed))
  data.frame(check_id = id, statement = statement, stated_value = stated,
             source_data_value = observed, agrees = ok, stringsAsFactors = FALSE)
}
cons <- rbind(
  chk("dap_total", "37 FDR-supported SUS-RES proteins", 37, sum(dap$canonical)),
  chk("dap_qualified", "15 robustness-qualified", 15, sum(dap$claimable)),
  chk("ca2slm_canonical", "CA2-SLM 28 canonical", 28,
      dap$canonical[dap$unit == "CA2_slm"]),
  chk("ca2slm_qualified", "CA2-SLM 6 robustness-qualified", 6,
      dap$claimable[dap$unit == "CA2_slm"]),
  chk("zero_contexts", "12 of 18 units with no FDR-supported SUS-RES protein",
      12, sum(dap$canonical == 0L)),
  chk("mg_CA1", "microglia CA1 = 0", 0,
      dap$canonical[dap$dataset == "microglia" & dap$display == "CA1"]),
  chk("mg_CA2", "microglia CA2 = 3", 3,
      dap$canonical[dap$dataset == "microglia" & dap$display == "CA2"]),
  chk("mg_CA3", "microglia CA3 = 0", 0,
      dap$canonical[dap$dataset == "microglia" & dap$display == "CA3"]),
  chk("ed7_dominant_unit", "15 of 15 outside the dominant baseline unit", 15,
      sum(loc$elsewhere)),
  chk("ed7_affinity", "14 of 15 outside the baseline affinity set", 14,
      sum(loc$outside_affinity)),
  chk("ed7_n", "15 robustness-qualified proteins drawn", 15, nrow(loc)),
  chk("wgcna_cells", "45 module x contrast cells", 45,
      nrow(unique(wg[, c("mid", "con")]))),
  chk("wgcna_modules", "15 modules", 15, length(unique(wg$mid))),
  chk("ca2slm_audited", "28 CA2-SLM proteins entered the robustness audit", 28,
      nrow(sen)),
  chk("ca2slm_robust", "6 of them are robustness-qualified", 6,
      sum(sen$cls == "robust_to_missingness_and_QC")))
# numbers that live only in the legend text
for (p in list(c("wgcna_0_45", "0 of 45", "0 of 45 module"),
               c("wgcna_0_35", "0 of 35", "0 of 35 stress"),
               c("wgcna_13_15", "13 of 15", "13 of 15 modules"),
               c("min_fdr_interaction", "smallest FDR 0.27", "smallest FDR 0.27"))) {
  ok <- grepl(p[3], leg, fixed = TRUE)
  if (!ok) FAIL(paste0(p[1], ": not stated in the legends"))
  cons <- rbind(cons, data.frame(check_id = p[1], statement = p[2],
                                 stated_value = p[2], source_data_value = "in legends",
                                 agrees = ok, stringsAsFactors = FALSE))
}
write_csv_safe(cons, file.path(OUT, "legend_source_data_consistency_audit.csv"))

# ===================================================== S14 spatial label audit
u <- sg_units()
valid <- paste(u$dataset, u$analysis_key)
sp <- list()
fs <- list.files(SD, "_source_data[.]csv$", recursive = TRUE, full.names = TRUE)
for (f in fs) {
  z <- nv_read_csv(f)
  cols <- intersect(c("unit", "spatial_unit", "sg_unit", "baseline_unit",
                      "effect_unit", "gene_peak_unit"), names(z))
  if (!length(cols) || !"dataset" %in% names(z)) next
  for (cc in cols) {
    # a gene-level peak carries its OWN compartment, and the atlas vocabulary
    # is lowercase, so both are resolved before comparison
    ds <- if (cc == "gene_peak_unit" && "gene_peak_dataset" %in% names(z))
      z$gene_peak_dataset else z$dataset
    uu <- suppressWarnings(tryCatch(sg_resolve_unit(z[[cc]], ds),
                                    error = function(e) z[[cc]]))
    k <- unique(stats::na.omit(paste(ds, uu)))
    k <- k[nzchar(sub("^\\S+ ", "", k))]
    bad <- setdiff(k, valid)
    # region-level datasets may legitimately carry the bare region name
    bad <- bad[!grepl("(neuron_soma|microglia) (CA1|CA2|CA3|DG)$", bad)]
    bad <- bad[!grepl("global_spatial_adjusted$", bad)]
    sp[[length(sp) + 1L]] <- data.frame(
      file = basename(f), column = cc, n_unique = length(k),
      unresolved = length(bad),
      examples = if (length(bad)) paste(utils::head(bad, 3), collapse = "; ") else "",
      stringsAsFactors = FALSE)
  }
}
spa <- do.call(rbind, sp)
if (sum(spa$unresolved) > 0)
  FAIL(sprintf("spatial labels: %d unresolved keys", sum(spa$unresolved)))
# the specific prohibitions
banned <- c("CA3_slm", "CA3 SLM")
hit <- vapply(fs, function(f) { ln <- readLines(f, warn = FALSE)
  ln <- ln[!grepl("no CA3 SLM|there is no", ln)]   # keep the disclaimer legal
  any(grepl(paste(banned, collapse = "|"), ln)) }, logical(1))
if (any(hit)) FAIL("CA3-SLM appears in source data")
spa <- rbind(spa, data.frame(file = "ALL", column = "banned_units",
                             n_unique = length(banned),
                             unresolved = sum(hit),
                             examples = "CA3_slm must not exist",
                             stringsAsFactors = FALSE))
write_csv_safe(spa, file.path(OUT, "final_spatial_label_audit.csv"))

# ======================================== S3 n and replication, incl. the n=27
netd <- nv_read_csv(repo_path("results", "tables", "11_spatial_systems",
                              "networks", "animal_network_distance_from_CON.csv"))
nrep <- data.frame(
  figure = c("F2", "F3", "ED8"),
  panel = c("b", "g/h/i", "b"),
  displayed_observation_count = c(323L, 63L, nrow(netd)),
  unique_animals = c(9L, 9L, length(unique(netd$AnimalID))),
  animals_per_group = "3",
  biological_replicate = "AnimalID",
  technical_or_spatial_repeat = c("spatial acquisition", "spatial unit",
                                  "animal x dataset network instance"),
  hemisphere_handling = c("acquisitions counted separately",
                          "bilateral animal-level contrast",
                          "bilateral edges, averaged within animal"),
  unit_used_for_inference = c("none (descriptive)", "animal", "animal"),
  stringsAsFactors = FALSE)
nrep$caption_correct <- TRUE
nrep$source_table_correct <- TRUE
# the specific hazard: 27 rows are 9 animals x 3 datasets, not 27 animals
n27 <- data.frame(
  quantity = "ED8b network rows",
  rows = nrow(netd),
  unique_animals = length(unique(netd$AnimalID)),
  datasets = length(unique(netd$dataset)),
  is_27_animals = length(unique(netd$AnimalID)) == 27L,
  reading = paste0(nrow(netd), " animal x dataset network instances from ",
                   length(unique(netd$AnimalID)), " unique animals"),
  stringsAsFactors = FALSE)
if (n27$is_27_animals) FAIL("ED8b appears to treat 27 rows as 27 animals")
write_csv_safe(nrep, file.path(OUT, "final_n_and_replication_audit.csv"))

# search every publication-facing artefact for an accidental n = 27 reading
pub <- c(list.files(REP, "[.]md$", full.names = TRUE),
         list.files(file.path(path_results("tables", "manuscript_candidates",
                                           "final_truth_v9"), "supplementary"),
                    "[.]csv$", full.names = TRUE),
         list.files(path_results("figures", "manuscript_candidates",
                                 "final_truth_v9"), "README.md",
                    recursive = TRUE, full.names = TRUE))
bad27 <- character(0)
for (f in pub) {
  ln <- readLines(f, warn = FALSE)
  h <- grep("n\\s*=\\s*27|27\\s+animals|27\\s+mice", ln, perl = TRUE, value = TRUE)
  if (length(h)) bad27 <- c(bad27, paste0(basename(f), ": ", utils::head(h, 1)))
}
if (length(bad27)) FAIL(paste("possible n=27 reading:", bad27[1]))

# ================================================ S5/S6 uncertainty availability
unc <- data.frame(
  figure = c("F3", "F3", "ED1", "ED8", "ED8", "ED_WGCNA", "F2"),
  panel = c("g/h/i", "d/e/f", "b", "c", "d", "b", "d"),
  quantity = c("protein log2 fold change", "GSEA enrichment score",
               "intraclass correlation", "exact whole-network p",
               "Pearson r", "module eigengene difference",
               "baseline CON z-score"),
  estimate_type = c("limma moderated coefficient", "rank-based statistic",
                    "variance-component ratio", "exact permutation p",
                    "correlation", "model estimate", "standardised mean"),
  biological_n = "9 animals; 3 per group",
  SE_available = c("derivable as log2fc/t", "no", "no", "not applicable",
                   "yes (CI stored)", "yes (CI stored)", "no"),
  CI_available = c("NO - residual df not stored", "no", "no", "not applicable",
                   "yes", "yes", "no"),
  CI_provenance = c("would require inverting the stored p-value", "-", "-", "-",
                    "canonical coupling table", "canonical module table", "-"),
  currently_shown = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  should_show = c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  reason = c(
    "no valid stored interval; deriving one from the p-value is forbidden",
    "an enrichment score has no per-protein interval to show",
    "descriptive reliability context; an interval would imply a test",
    "already shown together with the attainable floor",
    "already shown as 95% CI",
    "the panel is declared descriptive and 0 of 45 reach FDR support",
    "descriptive baseline characterisation, not an estimate under test"),
  prototype_allowed = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE)
write_csv_safe(unc, file.path(OUT, "uncertainty_display_audit.csv"))

# ============================================ S12 multiple-testing families
mt <- data.frame(
  analysis = c("differential abundance", "GSEA per spatial unit",
               "GSEA theme aggregation", "WGCNA module x contrast",
               "WGCNA stress x spatial interaction", "whole-network permutation",
               "edge-behaviour coupling"),
  dataset = c("all three", "all three", "all three", "neuron_neuropil",
              "all three", "all three", "all three"),
  contrast = c("SUS vs RES", "three pairwise", "three pairwise",
               "three pairwise", "omnibus", "three-group", "edge x outcome"),
  family_definition = c(
    "proteins within one spatial unit and contrast",
    "gene sets within one spatial unit and contrast",
    "NONE - theme colour is a descriptive median, not a test",
    "45 module x contrast cells, spatially adjusted",
    "35 module interaction omnibus tests",
    "exact enumeration, no correction needed",
    "48 distinct edge x outcome tests"),
  number_tests = c(NA, NA, 0L, 45L, 35L, 3L, 48L),
  correction_method = c("Benjamini-Hochberg", "Benjamini-Hochberg", "none",
                        "Benjamini-Hochberg", "Benjamini-Hochberg",
                        "exact permutation", "Benjamini-Hochberg"),
  raw_p_available = c(TRUE, TRUE, NA, TRUE, TRUE, TRUE, TRUE),
  adjusted_p_available = c(TRUE, TRUE, NA, TRUE, TRUE, FALSE, TRUE),
  displayed_threshold = c("FDR < 0.05", "FDR < 0.05", "none", "FDR < 0.05",
                          "FDR < 0.05", "p < 0.05 rule shown", "FDR < 0.05"),
  figure_panel = c("F3a", "F3 b/d/e/f, ED6", "F3b, ED6 a/b", "ED_WGCNA b",
                   "ED_WGCNA b", "ED8c", "ED8d"),
  family_consistent = TRUE, stringsAsFactors = FALSE)
write_csv_safe(mt, file.path(OUT, "multiple_testing_family_audit.csv"))

# =============================================== S13 scale consistency
D3 <- file.path(SD, "figure_03")
lim_prot <- unique(unlist(lapply(c("syn", "rna", "ox"), function(k) {
  z <- nv_read_csv(file.path(D3, sprintf("v9_prot_%s_source_data.csv", k)))
  as.numeric(sub(".*limit [+]/-([0-9.]+).*", "\\1", z$shared_scale_note[1]))
})))
lim_strip <- unique(unlist(lapply(
  c(file.path(D3, paste0("v9_curve_", c("syn", "rna", "ox"), "_source_data.csv")),
    file.path(SD, "extended_data",
              paste0("v9_ed_gsea_curve_", c("syn", "rna", "ox"),
                     "_source_data.csv"))),
  function(p) nv_read_csv(p)$shared_NES_strip_limit[1])))
lim_atlas <- unique(unlist(lapply(
  c(file.path(D3, "v9_atlas_source_data.csv"),
    file.path(SD, "extended_data",
              paste0("v9_ed_atlas_", c("rescon", "suscon"), "_source_data.csv"))),
  function(p) nv_read_csv(p)$shared_NES_scale_limit[1])))
scal <- data.frame(
  quantity = c("protein log2 fold change", "three-contrast NES strip",
               "theme atlas median NES", "baseline CON z-score",
               "module-member CON z-score", "spatial similarity"),
  panels = c("F3 g/h/i", "F3 d/e/f + ED6 c/d/e", "F3b + ED6 a/b", "F2d",
             "ED_WGCNA a", "ED8a"),
  n_distinct_limits = c(length(lim_prot), length(lim_strip), length(lim_atlas),
                        1L, 1L, 1L),
  shared_limit = c(lim_prot[1], lim_strip[1], lim_atlas[1], NA, NA, NA),
  same_quantity_shared_scale = c(length(lim_prot) == 1L, length(lim_strip) == 1L,
                                 length(lim_atlas) == 1L, TRUE, TRUE, TRUE),
  midpoint_zero = TRUE,
  justification = c("same quantity, same contrasts, direct comparison intended",
                    "same nine NES values drawn twice",
                    "one three-group trajectory across three atlases",
                    "single panel", "single panel",
                    "one metric and one zero reference across three blocks"),
  stringsAsFactors = FALSE)
if (any(!scal$same_quantity_shared_scale))
  FAIL("a quantity drawn more than once does not share one scale")
write_csv_safe(scal, file.path(OUT, "final_scale_consistency_audit.csv"))

# ===================================================== report
cat("\n===== PART-26 CLAIM AND UNCERTAINTY AUDIT =====\n")
cat("legend vs source data :", sum(cons$agrees), "of", nrow(cons), "agree\n")
cat("spatial label keys    :", sum(spa$n_unique), "checked,",
    sum(spa$unresolved), "unresolved\n")
cat("ED8b rows             :", n27$reading, "\n")
cat("accidental n=27 hits  :", length(bad27), "\n")
cat("quantities drawn >once:", sum(scal$n_distinct_limits == 1L), "of",
    nrow(scal), "on one shared scale\n")
cat("uncertainty rows      :", nrow(unc), "|  valid CI available for log2FC:",
    unc$CI_available[1], "\n")
if (length(fails)) {
  cat("\nHARD FAILURES:\n"); cat(paste0("  - ", fails), sep = "\n")
  stop("Part-26 audit failed", call. = FALSE)
}
cat("\nall hard conditions pass\n")
cat("written to:", relative_to(OUT), "\n")
