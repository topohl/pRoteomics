#!/usr/bin/env Rscript

# Part-26: user-facing, publication-ready statistics tables.
#
# The per-panel *_source_data.csv files are machine-oriented: abbreviated
# headers, and a prose provenance column repeated on every row. That is right
# for reproducing a panel and wrong for a reader. This script derives a
# Supplementary Table set from those same files - nothing is recomputed, no
# test is run - with human column names, explicit units, the n each statistic
# rests on, the test that produced it and the multiple-testing family it was
# corrected in. Provenance moves out of the data rows into a data dictionary.
#
# Every value here is a READ of an existing stored result. If a quantity is not
# already in the figure layer, it does not appear.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
source(repo_path("R", "final_truth_v9_figure_utils.R"))
source(repo_path("R", "final_truth_v9_fidelity_panels.R"))
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_supplementary_tables.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] final_truth_v9 supplementary statistics tables")
  quit(save = "no", status = 0L)
}

SD <- path_results("source_data", "manuscript_candidates", "final_truth_v9")
OUT <- path_results("tables", "manuscript_candidates", "final_truth_v9",
                    "supplementary")
dir_create(OUT)

sidecar <- function(fig, id) {
  p <- file.path(SD, fig, paste0(id, "_source_data.csv"))
  if (!file.exists(p)) return(NULL)
  nv_read_csv(p)
}
# keep only columns that exist, then rename; never invent a column
pick <- function(z, map) {
  have <- names(map)[names(map) %in% names(z)]
  if (!length(have)) return(NULL)
  out <- z[, have, drop = FALSE]
  names(out) <- unname(map[have])
  out
}
dict <- list()
emit <- function(tbl, file, title, note, defs) {
  if (is.null(tbl) || !nrow(tbl)) {
    message("[supp] skipped (no source): ", file); return(invisible(NULL))
  }
  write_csv_safe(tbl, file.path(OUT, file))
  dict[[length(dict) + 1L]] <<- data.frame(
    table_file = file, table_title = title,
    column = names(tbl),
    definition = unname(defs[names(tbl)]),
    stringsAsFactors = FALSE)
  attr(dict, "notes") <<- c(attr(dict, "notes"), stats::setNames(note, file))
  message(sprintf("[supp] %-46s %3d rows x %2d cols", file, nrow(tbl),
                  ncol(tbl)))
  invisible(NULL)
}

# ---------------------------------------------------------------- ST1
b <- sidecar("extended_data", "v9_ed_bilateral_full")
st1 <- pick(b, c(contrast = "Anatomical contrast", lab = "Contrast (display)",
                 level = "Identity level", n_pairs = "Paired proteins (n)",
                 pearson_r = "Pearson r (left vs right)",
                 spearman_rho = "Spearman rho (left vs right)",
                 sign_agreement_fraction = "Sign agreement (fraction)"))
if (!is.null(st1)) {
  st1$`Biological replicates (n animals)` <- 9L
  st1$`Statistical test` <- "none; stored descriptive concordance metrics"
  st1$Sidedness <- "not applicable"
  st1$`Degrees of freedom` <- NA_integer_
  st1$`Multiple-testing family` <- "not applicable"
}
emit(st1, "ST1_bilateral_reproducibility.csv",
     "Left-right reproducibility of every prespecified anatomical contrast",
     paste0("Descriptive concordance between hemispheres. No hypothesis test ",
            "is performed on these metrics, so no p-value or FDR is reported."),
     c("Anatomical contrast" = "prespecified CON-only anatomical contrast identifier",
       "Contrast (display)" = "display form of the same contrast",
       "Identity level" = "regional identity or fine/laminar identity",
       "Paired proteins (n)" = "proteins with a measured value in both hemispheres",
       "Pearson r (left vs right)" = "Pearson correlation of the per-protein effect between hemispheres",
       "Spearman rho (left vs right)" = "rank correlation of the same quantity",
       "Sign agreement (fraction)" = "fraction of proteins with the same sign in both hemispheres",
       "Biological replicates (n animals)" = "animals contributing; the animal is the biological replicate",
       "Statistical test" = "test that produced the statistic, or none",
       "Multiple-testing family" = "family the p-value was corrected within",
       "Sidedness" = "one- or two-sided, or not applicable",
       "Degrees of freedom" = "degrees of freedom of the test, where defined"))

# ---------------------------------------------------------------- ST2
ex <- sidecar("extended_data", "v9_ed_external_full")
# MT-04. This table used to map the column named p_adjust onto a column called
# "BH-adjusted p" and declare the FDR method to be Benjamini-Hochberg. In this
# analysis p_adjust is the RAW single-set value - Benjamini-Hochberg over a
# family of one is a no-op - so the table positively relabelled an uncorrected
# p as adjusted, in a manuscript-facing artefact with a matching data
# dictionary. The adjusted statistic is signature_FDR, corrected within the
# signature family, and that is what the Results claim rests on. Both are now
# reported, each under its true name.
st2a <- pick(ex, c(internal_contrast = "Internal anatomical contrast",
                   external_signature = "External hippocampal signature",
                   level = "Identity level", kind = "Pairing type",
                   NES = "Normalised enrichment score",
                   signature_FDR = "Signature-family FDR",
                   single_set_p_unadjusted = "Single-set p (uncorrected)"))
if (!is.null(st2a)) {
  st2a$`Statistical test` <- "gene set enrichment (stored)"
  st2a$`Multiple-testing family` <- "signature family"
  st2a$`FDR method` <- "Benjamini-Hochberg within the signature family"
}
emit(st2a, "ST2_external_signature_validation.csv",
     "All 30 tested internal contrast by external signature pairings",
     paste0("Expected pairings and specificity comparisons are both reported, ",
            "so specificity can be judged rather than assumed. The 9 contrasts ",
            "and 7 signatures form a 63-cell grid; only the 30 structurally ",
            "applicable pairs were tested."),
     c("Internal anatomical contrast" = "CON-only contrast defined in this study",
       "External hippocampal signature" = "published reference signature",
       "Identity level" = "regional or CA1 laminar identity",
       "Pairing type" = "expected pairing, or specificity comparison",
       "Normalised enrichment score" = "GSEA NES; positive = enriched in the first side of the contrast",
       "Signature-family FDR" = "Benjamini-Hochberg adjusted within the signature family; this is the statistic the Results claim rests on",
       "Single-set p (uncorrected)" = "the raw single-set value; Benjamini-Hochberg over a family of one is a no-op, so this is NOT an adjusted p and must not be read as one",
       "Statistical test" = "test that produced the statistic",
       "Multiple-testing family" = "family the p-value was corrected within",
       "FDR method" = "multiple-testing correction applied"))

ig <- sidecar("extended_data", "v9_ed_internal_full")
st2b <- pick(ig, c(contrast = "Anatomical contrast", Description = "GO term",
                   level = "Identity level", NES = "Normalised enrichment score",
                   p_adjust = "BH-adjusted p", setSize = "Gene set size (n proteins)"))
if (!is.null(st2b)) {
  st2b$`Statistical test` <- "gene set enrichment (stored)"
  st2b$`Multiple-testing family` <- "within the internal anatomical inventory"
  st2b$`FDR method` <- "Benjamini-Hochberg"
}
# This table is generated from the same 14-row sidecar as the Extended Data
# panel that was withheld in Phase 6A for claiming completeness it does not
# have. Withholding the picture and leaving the table under the same sentence
# would have moved the misrepresentation rather than removed it, so the title
# now states the selection. The genuinely complete inventory is a different,
# already released artefact and is named here so a reader can find it.
emit(st2b, "ST3_internal_anatomical_programs.csv",
     "Selected canonical GO terms: the strongest per anatomical contrast, for 7 of the 11 contrasts",
     paste0("CON-only. No stress information enters term selection. This is a ",
            "SELECTED SUBSET, not an inventory: 14 rows across 7 contrasts. The ",
            "complete canonical result is the released ",
            "control_anatomical_go_bp_gsea supplementary table, 40,680 rows ",
            "across all 11 contrasts, of which 2,826 are FDR-supported in the ",
            "positive direction."),
     c("Anatomical contrast" = "CON-only contrast identifier",
       "GO term" = "Gene Ontology biological process term",
       "Identity level" = "regional or CA1 laminar identity",
       "Normalised enrichment score" = "GSEA NES",
       "BH-adjusted p" = "Benjamini-Hochberg adjusted p-value",
       "Gene set size (n proteins)" = "measured proteins in the term",
       "Statistical test" = "test that produced the statistic",
       "Multiple-testing family" = "family the p-value was corrected within",
       "FDR method" = "multiple-testing correction applied"))

# ---------------------------------------------------------------- ST4
dp <- sidecar("figure_03", "v9_dap_track")
st4 <- pick(dp, c(display = "Spatial unit", dataset = "Compartment",
                  canonical = "FDR-supported SUS-RES proteins (n)",
                  claimable = "Robustness-qualified proteins (n)"))
if (!is.null(st4)) {
  st4$Compartment <- sg_compartment_label(st4$Compartment)
  st4$`Statistical test` <- "differential abundance, SUS vs RES (stored)"
  st4$`FDR method` <- "Benjamini-Hochberg"
  st4$`Biological replicates (n animals)` <- 9L
}
emit(st4, "ST4_differential_abundance_by_spatial_unit.csv",
     "FDR-supported and robustness-qualified SUS-RES proteins per spatial unit",
     paste0("Robustness qualification restricts CA2-SLM hits to those that ",
            "survive the missingness and QC audit; proteins never at risk in ",
            "CA2-SLM are retained unchanged."),
     c("Spatial unit" = "display name of the spatial unit",
       "Compartment" = "neuropil, neuronal soma, or microglia-enriched ROI",
       "FDR-supported SUS-RES proteins (n)" = "proteins reaching BH FDR < 0.05 for SUS vs RES in that unit",
       "Robustness-qualified proteins (n)" = "subset surviving the CA2-SLM robustness audit",
       "Statistical test" = "test that produced the counts",
       "FDR method" = "multiple-testing correction applied",
       "Biological replicates (n animals)" = "animals contributing",
       "Sidedness" = "one- or two-sided, or not applicable",
       "Degrees of freedom" = "degrees of freedom of the test, where defined"))

# ---------------------------------------------------------------- ST5
se <- sidecar("extended_data", "v9_ed_ca2_sensitivity")
st5 <- pick(se, c(gene = "Gene symbol",
                  canonical = "log2 fold change, SUS - RES (all data)",
                  hemi = "log2 fold change, SUS - RES (QC-failed hemispheres dropped)",
                  shrink = "Magnitude lost (log2 units)",
                  cls = "CA2-SLM robustness class"))
if (!is.null(st5)) {
if (!is.null(st5)) st5$`CA2-SLM robustness class` <-
  f9_qc_class_label(st5$`CA2-SLM robustness class`)
  st5$`Statistical test` <- "differential abundance, SUS vs RES (stored)"
  st5$`Biological replicates (n animals)` <- 9L
}
emit(st5, "ST5_ca2slm_robustness_sensitivity.csv",
     "CA2-SLM effect sizes with and without the QC-failed hemispheres",
     paste0("Both QC failures are a single hemisphere of an animal whose other ",
            "hemisphere passes, so the comparison is a within-animal ",
            "sensitivity analysis rather than an exclusion of animals."),
     c("Gene symbol" = "protein-coding gene symbol",
       "log2 fold change, SUS - RES (all data)" = "canonical stored effect size",
       "log2 fold change, SUS - RES (QC-failed hemispheres dropped)" = "same effect recomputed upstream without the two QC-failed acquisitions",
       "Magnitude lost (log2 units)" = "absolute canonical effect minus absolute recomputed effect",
       "CA2-SLM robustness class" = "stored robustness classification",
       "Statistical test" = "test that produced the effect sizes",
       "Biological replicates (n animals)" = "animals contributing",
       "Sidedness" = "one- or two-sided, or not applicable",
       "Degrees of freedom" = "degrees of freedom of the test, where defined"))

# ---------------------------------------------------------------- ST6
wg <- sidecar("extended_data", "v9_ed_wgcna_phenotype")
st6 <- pick(wg, c(mid = "Module", peak_label = "Peak spatial unit",
                  cell = "External cell-type affinity", con = "Contrast",
                  val = "Module eigengene difference",
                  tier_specific_fdr = "Tier-specific FDR"))
if (!is.null(st6)) {
  st6$`Statistical test` <- "module eigengene group difference, spatially adjusted (stored)"
  st6$`FDR method` <- "Benjamini-Hochberg, tier-specific family"
  st6$`FDR-supported at 0.05` <- is.finite(st6$`Tier-specific FDR`) &
    st6$`Tier-specific FDR` < 0.05
  st6$`Biological replicates (n animals)` <- 9L
}
emit(st6, "ST6_wgcna_module_phenotype.csv",
     "WGCNA module eigengene differences by contrast, with the inferential null",
     paste0("DESCRIPTIVE. No module x contrast cell reaches FDR support, and ",
            "no stress-by-spatial-unit interaction omnibus test does either. ",
            "The effect column must not be read as evidence of a group ",
            "difference."),
     c("Module" = "WGCNA module identifier",
       "Peak spatial unit" = "spatial unit where the module's mean CON abundance is highest",
       "External cell-type affinity" = "closest external reference cell type",
       "Contrast" = "group comparison",
       "Module eigengene difference" = "stored effect estimate",
       "Tier-specific FDR" = "Benjamini-Hochberg FDR within the tier-specific family",
       "Statistical test" = "test that produced the estimate",
       "FDR method" = "multiple-testing correction applied",
       "FDR-supported at 0.05" = "whether the cell reaches FDR < 0.05",
       "Biological replicates (n animals)" = "animals contributing",
       "Sidedness" = "one- or two-sided, or not applicable",
       "Degrees of freedom" = "degrees of freedom of the test, where defined"))

# ---------------------------------------------------------------- ST7
nl <- sidecar("extended_data", "v9_ed_nulls")
st7a <- pick(nl, c(lab = "Compartment", n_edges = "Network edges (n)",
                   p = "Exact whole-network p",
                   floor = "Smallest attainable p"))
cp <- sidecar("extended_data", "v9_ed_coupling")
st7b <- pick(cp, c(edge_lab = "Spatial unit pair", out_lab = "Behavioural outcome",
                   estimate = "Pearson r", conf.low = "95% CI lower",
                   conf.high = "95% CI upper", p.value = "p",
                   p.adj_BH_all_edge_phenotype_tests = "BH-adjusted p"))
if (!is.null(st7a)) {
  st7a$`Statistical test` <- "exact permutation over group labels (stored)"
  st7a$`Biological replicates (n animals)` <- 9L
}
if (!is.null(st7b)) {
  st7b$`Statistical test` <- "Pearson correlation across animals (stored)"
  st7b$`FDR method` <- "Benjamini-Hochberg over all edge-phenotype tests"
  st7b$Sidedness <- "two-sided"
  st7b$`Degrees of freedom` <- 7L
  st7b$`Biological replicates (n animals)` <- 9L
}
emit(st7a, "ST7_network_whole_graph_nulls.csv",
     "Whole-network group comparison with its attainable resolution",
     paste0("The attainable floor is the smallest p the exact enumeration ",
            "could return at this sample size; it is reported so the null is ",
            "informative rather than merely non-significant."),
     c("Compartment" = "neuropil, neuronal soma, or microglia-enriched ROI",
       "Network edges (n)" = "edges in the spatial molecular network",
       "Exact whole-network p" = "exact permutation p-value for a whole-network group difference",
       "Smallest attainable p" = "resolution floor of the enumeration at this n",
       "Statistical test" = "test that produced the p-value",
       "Biological replicates (n animals)" = "animals contributing",
       "Sidedness" = "one- or two-sided, or not applicable",
       "Degrees of freedom" = "degrees of freedom of the test, where defined"))
emit(st7b, "ST8_edge_behaviour_coupling.csv",
     "Every tested neuropil spatial-unit-pair by behavioural-outcome correlation",
     paste0("With 9 animals a single correlation has very little resolution; ",
            "the complete inventory is reported so no subset can be mined."),
     c("Spatial unit pair" = "the two spatial units forming the network edge",
       "Behavioural outcome" = "per-animal behavioural or physiological outcome",
       "Pearson r" = "correlation across animals",
       "95% CI lower" = "lower bound of the 95% confidence interval",
       "95% CI upper" = "upper bound of the 95% confidence interval",
       "p" = "uncorrected p-value",
       "BH-adjusted p" = "Benjamini-Hochberg adjusted p-value",
       "Statistical test" = "test that produced the statistic",
       "FDR method" = "multiple-testing correction applied",
       "Biological replicates (n animals)" = "animals contributing",
       "Sidedness" = "one- or two-sided, or not applicable",
       "Degrees of freedom" = "degrees of freedom of the test, where defined"))

# ------------------------------------------------------------ dictionary
dd <- do.call(rbind, dict)
notes <- attr(dict, "notes")
dd$table_note <- unname(notes[dd$table_file])
write_csv_safe(dd, file.path(OUT, "ST0_data_dictionary.csv"))
message(sprintf("[supp] %-46s %3d rows", "ST0_data_dictionary.csv", nrow(dd)))

cat("\n===== SUPPLEMENTARY TABLES =====\n")
inv <- unique(dd[, c("table_file", "table_title")])
print(inv, row.names = FALSE)
cat("\nwritten to:", relative_to(OUT), "\n")
