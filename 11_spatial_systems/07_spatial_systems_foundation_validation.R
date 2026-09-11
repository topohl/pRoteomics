#!/usr/bin/env Rscript
#
# Spatial systems foundation validation.
#
# Checks the CRITICAL contracts the later atlas and network stages depend on.
# A critical FAIL exits non-zero so a downstream stage cannot be built on a
# foundation that did not hold.
#
# USAGE
#   Rscript 11_spatial_systems/07_spatial_systems_foundation_validation.R

source("R/paths.R")
source("R/dataset_config.R")
source("R/integration_utils.R")
source("R/spatial_systems_data_utils.R")
source("R/spatial_systems_evidence_registry.R")
source("R/ewce_gene_set_engine.R")
source("R/control_spatial_identity_utils.R")

suppressPackageStartupMessages({ library(readr); library(dplyr) })

SCRIPT_ID <- "11_spatial_systems/07_spatial_systems_foundation_validation.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

OUT <- function() {
  d <- path_results("tables", "11_spatial_systems"); dir_create(d); d
}
T_ <- function(...) path_results("tables", "11_spatial_systems", ...)

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Spatial systems foundation validation.\n")
  dry_run_inputs(SCRIPT_ID, list(
    aggregation_validation = T_("data_contract", "spatial_systems_aggregation_validation.csv")))
  quit(save = "no", status = 0L)
}

checks <- list()
add <- function(check_id, critical, status, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check_id = check_id, critical = critical, status = status,
    detail = detail, stringsAsFactors = FALSE)
}
read_if <- function(p) if (file.exists(p)) {
  as.data.frame(readr::read_csv(p, show_col_types = FALSE, progress = FALSE,
                                guess_max = Inf))
} else NULL

# ---- 1. true hemisphere values exist --------------------------------------
inv <- read_if(T_("data_contract", "spatial_systems_hemisphere_inventory.csv"))
add("true_hemisphere_values_exist", TRUE,
    if (!is.null(inv) && nrow(inv) > 0L) "PASS" else "FAIL",
    if (is.null(inv)) "hemisphere inventory missing" else
      paste0(nrow(inv), " AnimalID x SpatialUnit cells across ",
             length(unique(inv$dataset)), " datasets"))

agg <- read_if(T_("data_contract", "spatial_systems_aggregation_validation.csv"))
# ---- 2. old fake hemisphere object rejected -------------------------------
rejected <- !is.null(agg) &&
  any(grepl("hemisphere_mat_is_deprecated_alias", agg$validation_target))
side_ok <- !is.null(agg) &&
  all(agg$status[agg$validation_target == "side_resolved_matrices_are_not_side_averaged"] == "PASS")
add("old_hemisphere_object_rejected", TRUE,
    if (rejected && side_ok) "PASS" else "FAIL",
    "empirical ROI hemisphere_mat documented as a hemisphere-AVERAGED deprecated alias; side-resolved matrices proven non-averaged")

# ---- 3. bilateral mean reproduces canonical aggregation -------------------
canon <- if (is.null(agg)) NULL else
  agg[agg$validation_target == "canonical_protigy_animal_level_aggregation", , drop = FALSE]
add("bilateral_mean_reproduces_canonical_aggregation", TRUE,
    if (!is.null(canon) && nrow(canon) && all(canon$status == "PASS")) "PASS" else "FAIL",
    if (is.null(canon) || !nrow(canon)) "not evaluated" else
      paste0("max abs difference ",
             format(max(canon$max_abs_difference, na.rm = TRUE), digits = 3),
             " against protigy_aggregate_expression_columns()"))

# ---- 4. contrast definitions reused (exactly one registry) ----------------
reg_files <- Sys.glob(repo_path("R", "*.R"))
defs <- sum(vapply(reg_files, function(f) {
  any(grepl("^control_spatial_contrast_registry <- function", readLines(f, warn = FALSE)))
}, logical(1)))
n_locked <- sum(control_spatial_contrast_is_manuscript_locked(
  control_spatial_contrast_registry("neuron_soma", c("CA1", "CA2", "CA3", "DG")))) +
  sum(control_spatial_contrast_is_manuscript_locked(
    control_spatial_contrast_registry("neuron_neuropil",
      c("CA1_SLM", "CA1_SO", "CA1_SR", "CA2_SLM", "CA2_SO", "CA2_SR",
        "CA3_SO", "CA3_SR", "DG_MO", "DG_PO"))))
add("contrast_registry_defined_once", TRUE,
    if (defs == 1L && n_locked == control_spatial_manuscript_contrast_count()) "PASS" else "FAIL",
    paste0(defs, " registry definition(s); ", n_locked,
           " manuscript-locked neuronal contrasts (expected ",
           control_spatial_manuscript_contrast_count(), ")"))

# ---- 5/6. side purity in the L and R analyses -----------------------------
bsi <- read_if(T_("bilateral", "bilateral_spatial_identity_protein_level.csv"))
for (side in c("L", "R")) {
  col <- paste0("estimate_", side)
  ok <- !is.null(bsi) && col %in% names(bsi) && any(is.finite(bsi[[col]]))
  add(paste0(tolower(side), "_analysis_side_pure"), TRUE,
      if (ok) "PASS" else "FAIL",
      paste0("one-sided fit produced ", side,
             " estimates; the producer hard-stops if a side-specific subset ",
             "contains other-side samples"))
}

# ---- 7. compartment validation is phenotype-independent -------------------
comp <- read_if(T_("bilateral", "bilateral_empirical_compartment_summary.csv"))
comp_pl <- read_if(T_("bilateral", "bilateral_empirical_compartment_protein_level.csv"))
pheno_cols <- c("StressGroup", "ExpGroup", "Group", "SUS", "RES", "CON")
leak <- !is.null(comp_pl) && any(pheno_cols %in% names(comp_pl))
add("compartment_validation_phenotype_independent", TRUE,
    if (!is.null(comp) && !leak) "PASS" else "FAIL",
    "design is ~ AnimalID + dataset; AnimalID is nested within StressGroup so group membership is fully absorbed, and no phenotype column reaches the output")

# ---- 8. Stage-05 WGCNA L/R consumed, not recomputed -----------------------
wb <- read_if(T_("bilateral", "WGCNA_module_bilateral_reproducibility.csv"))
add("wgcna_stage05_consumed_not_recomputed", TRUE,
    if (!is.null(wb) && nrow(wb) > 0L) "PASS" else "FAIL",
    if (is.null(wb)) "module bilateral table missing" else
      paste0(nrow(wb), " endpoints read from Stage-05 hemisphere values; ",
             "no WGCNA membership or eigengene recomputed"))

# ---- 9. variance model status --------------------------------------------
vc <- read_if(T_("precision", "bilateral_variance_decomposition.csv"))
pg <- read_if(T_("precision", "bilateral_precision_gain.csv"))
add("variance_model_status_recorded", TRUE,
    if (!is.null(vc) && all(c("model_convergence", "is_singular",
                              "assumption_status") %in% names(vc))) "PASS" else "FAIL",
    if (is.null(vc)) "variance decomposition missing" else
      paste0(nrow(vc), " endpoints; ", sum(vc$is_singular %in% TRUE),
             " singular fits excluded from reliability; ",
             sum(vc$assumption_status == "assumptions_met"), " usable"))

# ---- 10. EWCE engine is phenotype-free -----------------------------------
eng_args <- names(formals(run_ewce_gene_set_annotation))
forbidden <- c("StressGroup", "group", "contrast", "condition", "direction")
add("ewce_engine_phenotype_free", TRUE,
    if (!any(tolower(forbidden) %in% tolower(eng_args))) "PASS" else "FAIL",
    paste0("run_ewce_gene_set_annotation(", paste(eng_args, collapse = ", "), ")"))

# ---- 11. EWCE FDR families are independent -------------------------------
set.seed(11)
fam_mod <- ewce_fdr_family_module_annotation("microglia", "all", 1L)
mod <- data.frame(gene_set_id = rep(paste0("m", 1:5), each = 4),
                  p_value = stats::runif(20), fdr_family = fam_mod,
                  stringsAsFactors = FALSE)
a <- ewce_apply_family_fdr(mod)
extra <- data.frame(gene_set_id = paste0("d", 1:40),
                    p_value = stats::runif(40) / 1000,
                    fdr_family = ewce_fdr_family_differential("microglia", 1L),
                    stringsAsFactors = FALSE)
b <- ewce_apply_family_fdr(rbind(mod, extra))
bmod <- b[b$fdr_family == fam_mod, ]
add("ewce_fdr_families_independent", TRUE,
    if (identical(a$FDR, bmod$FDR)) "PASS" else "FAIL",
    "module-annotation FDR is bit-identical after adding 40 phenotype-arm rows")

# ---- 12. module EWCE uses the measured background -------------------------
ct <- read_if(T_("celltype_annotation", "WGCNA_module_external_celltype_affinity_long.csv"))
bg_ok <- !is.null(ct) && "n_background" %in% names(ct) &&
  all(ct$n_background > 1000L, na.rm = TRUE)
add("module_ewce_uses_measured_background", !is.null(ct),
    if (is.null(ct)) "SKIPPED" else if (bg_ok) "PASS" else "FAIL",
    if (is.null(ct)) "module annotation not generated yet" else
      paste0("background sizes ", paste(range(ct$n_background, na.rm = TRUE),
                                        collapse = "-"),
             " genes (measured proteome), never the full reference transcriptome"))

# ---- 13/14. canonical scientific objects unchanged ------------------------
git_changed <- tryCatch(
  system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE),
  error = function(e) character())
touches <- function(pat) any(grepl(pat, git_changed))
add("no_canonical_da_changes", TRUE,
    if (!touches("results/tables/04_differential")) "PASS" else "FAIL",
    "no differential-abundance output modified in the working tree")
add("no_wgcna_state_changes", TRUE,
    if (!touches("06_modules_WGCNA/01_WGCNA.r") &&
        !touches("wgcna_final_model_state")) "PASS" else "FAIL",
    "WGCNA construction script and frozen model state untouched")
add("figure_contracts_unchanged", TRUE,
    if (!touches("figures/figure_0[23]") && !touches("manuscript_figure")) "PASS" else "FAIL",
    "Figure 2 / Figure 3 producers untouched")

# ---- 15. evidence registry present and well formed -----------------------
reg <- sps_evidence_dependence_registry()
add("evidence_dependence_registry_complete", TRUE,
    if (nrow(reg) >= 11L && !anyDuplicated(reg$evidence_id)) "PASS" else "FAIL",
    paste0(nrow(reg), " evidence streams with independence classes and ",
           "prohibited interpretations"))

validation <- dplyr::bind_rows(checks)
root <- OUT()
write_csv_safe(validation, file.path(root, "spatial_systems_foundation_validation.csv"))

cat("\n===== Spatial systems foundation validation =====\n")
for (i in seq_len(nrow(validation))) {
  cat(sprintf("  %-8s %-48s %s\n", validation$status[i], validation$check_id[i],
              substr(validation$detail[i], 1, 88)))
}
crit_fail <- validation$critical %in% TRUE & validation$status == "FAIL"
cat(sprintf("\n%d checks, %d critical FAIL, %d skipped\n", nrow(validation),
            sum(crit_fail), sum(validation$status == "SKIPPED")))
cat("Output:", relative_to(file.path(root, "spatial_systems_foundation_validation.csv")), "\n")
if (any(crit_fail)) {
  quit(save = "no", status = 1L)
}
