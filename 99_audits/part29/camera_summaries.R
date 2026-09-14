#!/usr/bin/env Rscript

# Part-29 sections 23-25: the CAMERA/GSEA concordance summaries.
# Reads the concordance table written by camera_sensitivity.R so the expensive
# cameraPR pass is not repeated. No inference is added here.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
EXEMPLARS <- c("GO:0099536", "GO:0006397", "GO:0006119")

j <- utils::read.csv(file.path(AUD, "gsea_camera_concordance.csv"),
                     stringsAsFactors = FALSE)
cat("concordance rows:", nrow(j), "| comparisons:",
    length(unique(j$source_comparison)), "\n")

j$gsea_direction <- ifelse(j$NES > 0, "Up", "Down")
j$direction_concordant <- j$gsea_direction == j$camera_Direction
j$gsea_supported <- is.finite(j$GSEA_FDR) & j$GSEA_FDR < 0.05
j$camera_supported <- is.finite(j$camera_FDR) & j$camera_FDR < 0.05

signed_camera <- ifelse(j$camera_Direction == "Up", 1, -1) *
  (-log10(pmax(j$camera_PValue, .Machine$double.xmin)))

conc <- do.call(rbind, lapply(split(seq_len(nrow(j)), j$source_comparison),
  function(ix) {
    z <- j[ix, , drop = FALSE]; sc <- signed_camera[ix]
    ok <- is.finite(z$NES) & is.finite(sc)
    data.frame(
      source_comparison = z$source_comparison[1], dataset = z$dataset_gsea[1],
      spatial_unit = z$spatial_unit[1], contrast = z$contrast[1],
      n_terms = nrow(z),
      spearman_NES_vs_signed_camera = if (sum(ok) > 2)
        suppressWarnings(stats::cor(z$NES[ok], sc[ok], method = "spearman"))
        else NA_real_,
      directional_concordance = mean(z$direction_concordant, na.rm = TRUE),
      n_gsea_supported = sum(z$gsea_supported),
      frac_gsea_supported_dir_concordant = if (any(z$gsea_supported))
        mean(z$direction_concordant[z$gsea_supported], na.rm = TRUE) else NA_real_,
      frac_gsea_supported_also_camera_supported = if (any(z$gsea_supported))
        mean(z$camera_supported[z$gsea_supported], na.rm = TRUE) else NA_real_,
      n_camera_supported = sum(z$camera_supported),
      stringsAsFactors = FALSE)
  }))
utils::write.csv(conc, file.path(AUD,
  "gsea_camera_concordance_by_comparison.csv"), row.names = FALSE)

# ------------------------------------------------------ S24 the three exemplars
ex <- j[j$GO_ID %in% EXEMPLARS, , drop = FALSE]
ex$classification <- with(ex, ifelse(
  !is.finite(camera_PValue), "NOT_EVALUABLE",
  ifelse(direction_concordant & camera_supported, "CONCORDANT_STRONG",
  ifelse(direction_concordant, "CONCORDANT_DIRECTION_ONLY", "DISCORDANT"))))
ex$inter_gene_cor_assumption <- 0.01
ex <- ex[order(ex$GO_ID, ex$contrast, ex$dataset_gsea, ex$spatial_unit), ,
         drop = FALSE]
utils::write.csv(ex, file.path(AUD, "camera_three_exemplar_audit.csv"),
                 row.names = FALSE)

# ----------------------------------------------- S25 per-theme, no theme p-value
TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)
elig <- unique(TH[TH$theme_claim_eligible %in% TRUE, c("theme_id", "GO_ID")])
desc <- stats::setNames(TH$GO_description[!duplicated(TH$GO_ID)],
                        TH$GO_ID[!duplicated(TH$GO_ID)])
thm <- do.call(rbind, lapply(sort(unique(elig$theme_id)), function(t) {
  ids <- elig$GO_ID[elig$theme_id == t]
  z <- j[j$GO_ID %in% ids, , drop = FALSE]
  data.frame(theme_id = t, n_constituent_GO_terms = length(ids),
             n_camera_tested = length(unique(z$GO_ID)),
             n_term_comparisons = nrow(z),
             n_direction_concordant = sum(z$direction_concordant, na.rm = TRUE),
             fraction_direction_concordant = mean(z$direction_concordant, na.rm = TRUE),
             n_camera_FDR_supported = sum(z$camera_supported),
             n_gsea_FDR_supported = sum(z$gsea_supported),
             fraction_gsea_supported_also_camera = if (any(z$gsea_supported))
               mean(z$camera_supported[z$gsea_supported]) else NA_real_,
             representative_exact_GO_terms = paste(utils::head(
               unname(desc[sort(unique(z$GO_ID))]), 4), collapse = "; "),
             theme_level_p_value = "NOT COMPUTED - themes are descriptive umbrellas",
             stringsAsFactors = FALSE)
}))
utils::write.csv(thm, file.path(AUD, "camera_theme_concordance.csv"),
                 row.names = FALSE)

cat("\n===== CAMERA / GSEA CONCORDANCE =====\n")
cat("term-comparisons matched:", nrow(j), "\n")
cat("directional concordance overall:",
    sprintf("%.4f", mean(j$direction_concordant, na.rm = TRUE)), "\n")
cat("GSEA FDR-supported:", sum(j$gsea_supported),
    "| of those direction-concordant:",
    sprintf("%.4f", mean(j$direction_concordant[j$gsea_supported], na.rm = TRUE)),
    "| also CAMERA FDR-supported:",
    sprintf("%.4f", mean(j$camera_supported[j$gsea_supported], na.rm = TRUE)), "\n")
cat("CAMERA FDR-supported overall:", sum(j$camera_supported), "\n")
cat("median per-comparison Spearman:",
    sprintf("%.4f", stats::median(conc$spearman_NES_vs_signed_camera, na.rm = TRUE)),
    "| range", sprintf("%.3f", min(conc$spearman_NES_vs_signed_camera, na.rm = TRUE)),
    "to", sprintf("%.3f", max(conc$spearman_NES_vs_signed_camera, na.rm = TRUE)), "\n")
cat("\nexemplar classification:\n"); print(table(ex$classification))
cat("\nthree exemplars in their own SUS-RES cells:\n")
k <- ex[ex$contrast == "SUS - RES" &
          ((ex$GO_ID == "GO:0099536" & ex$spatial_unit == "CA3_sr") |
           (ex$GO_ID == "GO:0006397" & ex$spatial_unit == "CA2_sp") |
           (ex$GO_ID == "GO:0006119" & ex$spatial_unit == "CA1" &
              ex$dataset_gsea == "microglia")), ]
print(k[, c("GO_ID", "dataset_gsea", "spatial_unit", "NES", "GSEA_FDR",
            "camera_Direction", "camera_PValue", "camera_FDR",
            "classification")], row.names = FALSE)
cat("\nper-theme:\n")
print(thm[, c("theme_id", "n_camera_tested", "fraction_direction_concordant",
              "fraction_gsea_supported_also_camera")], row.names = FALSE)
