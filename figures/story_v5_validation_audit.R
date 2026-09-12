#!/usr/bin/env Rscript

# Part-20 step 1: structural audit of what canonical Figure 2f and 2g actually
# contain, BEFORE anything is re-plotted.
#
# The finding that drives the redesign: the sparsity in both panels is not
# missing data, it is anatomical hierarchy. A regional contrast cannot be
# validated against a laminar signature, and the internal GO panel is a paired
# list of two terms per contrast rather than a term-by-contrast matrix.
# Structurally inapplicable combinations are recorded as such and are never
# drawn as empty cells.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
suppressPackageStartupMessages({ library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/story_v5_validation_audit.R")

B <- path_results("source_data", "04_differential_expression_enrichment",
                  "control_spatial_identity_validation", "global")
OUT <- function(...) {
  d <- path_results("tables", "manuscript_candidates", "story_v5"); dir_create(d)
  file.path(d, ...)
}

k <- utils::read.csv(file.path(B, "figure2e_source_data.csv"), stringsAsFactors = FALSE)
g <- utils::read.csv(file.path(B, "figure2f_regions_CA1layers_source_data.csv"),
                     stringsAsFactors = FALSE)

# validation_domain is the authoritative hierarchy field: parsing contrast
# names misclassifies CA1_SO_vs_CA3_SO, which is a REGIONAL comparison holding
# the layer constant despite carrying a layer token in its name.
lvl_from_domain <- function(x) {
  ifelse(grepl("strata", x), "laminar",
         ifelse(grepl("region", x), "regional", "other"))
}

kr <- data.frame(
  source_panel = "2f_external_kaulich",
  dataset = k$dataset,
  contrast = k$internal_contrast,
  anatomical_level = lvl_from_domain(k$validation_domain),
  validation_domain = k$validation_domain,
  target = k$external_signature,
  evidence_type = "external Kaulich signature",
  term_or_signature = k$external_signature,
  effect_NES = k$NES,
  canonical_FDR = k$p_adjust,
  FDR_supported = k$p_adjust < 0.05,
  set_size = k$mapped_unique_genes,
  is_expected_pairing = k$expected_match,
  match_type = k$match_type,
  structurally_applicable = TRUE,
  stringsAsFactors = FALSE)

gr <- data.frame(
  source_panel = "2g_internal_anatomical_GO",
  dataset = g$dataset,
  contrast = g$contrast,
  anatomical_level = ifelse(grepl("CA1_strata", g$contrast), "laminar", "regional"),
  validation_domain = ifelse(grepl("CA1_strata", g$contrast),
                             "CA1 strata - internal GO", "region - internal GO"),
  target = g$contrast,
  evidence_type = "internal anatomical GO/GSEA",
  term_or_signature = g$Description,
  effect_NES = g$NES,
  canonical_FDR = g$p_adjust,
  FDR_supported = g$p_adjust < 0.05,
  set_size = g$setSize,
  is_expected_pairing = NA,
  match_type = "internal_program_identity",
  structurally_applicable = TRUE,
  stringsAsFactors = FALSE)

# the inapplicable combinations, recorded explicitly so nobody later reads the
# gaps as failed or missing measurements
sig_level <- c(CA1 = "regional", "CA2/3" = "regional", DG = "regional",
               SLM = "laminar", SO = "laminar", SP = "laminar", SR = "laminar")
grid <- expand.grid(contrast = unique(k$internal_contrast),
                    target = unique(k$external_signature),
                    stringsAsFactors = FALSE)
have <- paste(k$internal_contrast, k$external_signature)
grid <- grid[!(paste(grid$contrast, grid$target) %in% have), , drop = FALSE]
dom <- k$validation_domain[match(grid$contrast, k$internal_contrast)]
inapp <- data.frame(
  source_panel = "2f_external_kaulich", dataset = NA_character_,
  contrast = grid$contrast, anatomical_level = lvl_from_domain(dom),
  validation_domain = dom, target = grid$target,
  evidence_type = "external Kaulich signature",
  term_or_signature = grid$target, effect_NES = NA_real_,
  canonical_FDR = NA_real_, FDR_supported = NA, set_size = NA_integer_,
  is_expected_pairing = NA, match_type = "not_applicable",
  structurally_applicable = FALSE, stringsAsFactors = FALSE)
inapp$reason_inapplicable <- sprintf(
  "a %s contrast cannot be validated against a %s reference signature",
  inapp$anatomical_level, unname(sig_level[inapp$target]))

kr$reason_inapplicable <- NA_character_
gr$reason_inapplicable <- NA_character_
audit <- rbind(kr, gr, inapp)
write_csv_safe(audit, OUT("figure2_spatial_validation_hierarchy_audit.csv"))

cat("\n===== Figure 2f/2g structural audit =====\n")
cat("rows:", nrow(audit), "  applicable:", sum(audit$structurally_applicable),
    "  structurally inapplicable:", sum(!audit$structurally_applicable), "\n")

cat("\n--- 2f: occupancy is block-diagonal by anatomical level ---\n")
kk <- audit[audit$source_panel == "2f_external_kaulich" &
              audit$structurally_applicable, ]
kk$tlevel <- unname(sig_level[kk$target])
print(table(contrast_level = kk$anatomical_level, signature_level = kk$tlevel))
cat("present cells:", nrow(kk), "of", length(unique(k$internal_contrast)) *
      length(unique(k$external_signature)),
    sprintf("(%.0f%% of a full matrix would be structurally meaningless)",
            100 * sum(!audit$structurally_applicable) /
              (length(unique(k$internal_contrast)) * length(unique(k$external_signature)))), "\n")
cat("expected pairings:", sum(kk$is_expected_pairing %in% TRUE),
    " specificity comparisons:", sum(kk$is_expected_pairing %in% FALSE), "\n")
cat("FDR-supported:", sum(kk$FDR_supported, na.rm = TRUE), "of", nrow(kk), "\n")

cat("\n--- 2g: a paired list, not a matrix ---\n")
gg <- audit[audit$source_panel == "2g_internal_anatomical_GO", ]
cat("rows:", nrow(gg), "  contrasts:", length(unique(gg$contrast)),
    "  distinct terms:", length(unique(gg$term_or_signature)), "\n")
cat("terms per contrast:", paste(unique(table(gg$contrast)), collapse = ","), "\n")
cat("terms unique to a single contrast:",
    sum(table(gg$term_or_signature) == 1), "of", length(unique(gg$term_or_signature)), "\n")
cat("a term x contrast matrix would be",
    sprintf("%.0f%% blank", 100 * (1 - nrow(gg) /
      (length(unique(gg$term_or_signature)) * length(unique(gg$contrast))))), "\n")
cat("FDR-supported:", sum(gg$FDR_supported), "of", nrow(gg), "\n")
cat("set_size range:", paste(range(gg$set_size), collapse = " - "),
    "(stored, so it can drive point size)\n")
cat("\nwritten:", relative_to(OUT("figure2_spatial_validation_hierarchy_audit.csv")), "\n")
