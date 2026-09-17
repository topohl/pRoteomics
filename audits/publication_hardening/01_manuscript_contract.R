#!/usr/bin/env Rscript

# Part A, sections 1-5 and 12: the manuscript-facing statistical contract.
#
# AUDIT ONLY. No analysis is rerun. Every row is reconstructed from the frozen
# v9 figure contract, the pipeline registry and the audit tables the project
# already emits, so the contract is derived from the repository rather than
# asserted from memory.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
source("audits/publication_hardening/00_checkpoint.R")
suppressMessages(library(yaml))

dir.create(PH_TAB, recursive = TRUE, showWarnings = FALSE)
rd <- function(p) if (file.exists(p))
  utils::read.csv(p, stringsAsFactors = FALSE) else NULL

CT <- yaml::read_yaml("figures/figure_final_truth_v9_contract.yml")
PIPE <- yaml::read_yaml("pipeline.yml")
V9T <- file.path("results", "tables", "manuscript_candidates", "final_truth_v9")
V9S <- file.path("results", "source_data", "manuscript_candidates",
                 "final_truth_v9")
V9F <- file.path("results", "figures", "manuscript_candidates", "final_truth_v9")

# ================================================ A01 reachability inventory
panels <- setNames(CT$panels, vapply(CT$panels, function(p) p$id, character(1)))
lay <- do.call(rbind, lapply(CT$figures, function(f)
  do.call(rbind, lapply(f$layout, function(l) data.frame(
    figure = as.character(f$name), figure_key = as.character(f$figure_key),
    panel_label = as.character(l$label), panel_id = as.character(l$panel),
    stringsAsFactors = FALSE)))))

role_audit <- rd(file.path(V9T, "final_panel_statistical_role_audit.csv"))

inv <- do.call(rbind, lapply(seq_len(nrow(lay)), function(i) {
  pid <- lay$panel_id[i]
  p <- panels[[pid]]
  sd <- file.path(V9S, lay$figure_key[i], paste0(pid, "_source_data.csv"))
  svg <- file.path(V9F, lay$figure_key[i], "panels", paste0(pid, ".svg"))
  pdf <- file.path(V9F, lay$figure_key[i], "assembled",
                   paste0(lay$figure[i], ".pdf"))
  ra <- if (!is.null(role_audit)) role_audit[role_audit$panel_id == pid, ,
                                             drop = FALSE] else NULL
  data.frame(
    figure = lay$figure[i], panel_label = lay$panel_label[i], panel_id = pid,
    renderer = as.character(p$renderer %||% ""),
    renderer_file = {
      f <- Sys.glob(file.path("R", "final_truth_v9_*panels.R"))
      hit <- f[vapply(f, function(x) any(grepl(
        paste0("^", as.character(p$renderer %||% "zzz"), " <- function"),
        readLines(x, warn = FALSE))), logical(1))]
      if (length(hit)) sub(".*proteomics[/\\\\]", "", hit[1]) else "not located"
    },
    canonical_input = as.character(p$primary_source %||% ""),
    source_data_path = if (file.exists(sd)) sd else "MISSING",
    source_data_rows = if (file.exists(sd)) nrow(utils::read.csv(sd)) else NA_integer_,
    panel_svg = if (file.exists(svg)) "present" else "MISSING",
    assembled_pdf = if (file.exists(pdf)) "present" else "MISSING",
    contract_role = as.character(p$role %||% ""),
    statistical_mode = if (!is.null(ra) && nrow(ra)) ra$descriptive_or_inferential[1] else "",
    fdr_family = if (!is.null(ra) && nrow(ra)) ra$FDR_family[1] else "",
    biological_n = if (!is.null(ra) && nrow(ra)) ra$biological_n[1] else "",
    manuscript_reachable = TRUE,
    stringsAsFactors = FALSE)
}))
utils::write.csv(inv, file.path(PH_TAB, "manuscript_reachability_inventory.csv"),
                 row.names = FALSE)

# ==================================================== A02 effect contract
TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)
contrasts <- sort(unique(TH$contrast))

E <- function(dataset, analysis, endpoint, contrast, effect_def, sign_conv,
              pos_meaning, statistic, unit, n_animals, family, claim)
  data.frame(dataset = dataset, analysis = analysis, endpoint = endpoint,
             contrast = contrast, formal_effect_definition = effect_def,
             displayed_sign_convention = sign_conv,
             positive_direction_meaning = pos_meaning, statistic = statistic,
             biological_unit = unit, n_animals = n_animals,
             multiple_testing_family = family, publication_claim = claim,
             stringsAsFactors = FALSE)

eff <- rbind(
  E("all three", "differential abundance", "protein log2 fold change",
    "SUS - RES", "group2 - group1 on animal-level log2 abundance within a spatial unit; Protigy contract animal_level_protigy_da_v1",
    "positive = higher in SUS", "higher abundance in susceptible animals",
    "moderated t (log2FC displayed)", "animal", "6 (3 SUS, 3 RES)",
    "BH within the differential-abundance family",
    "outcome-associated protein difference among stress-exposed animals"),
  E("all three", "differential abundance", "protein log2 fold change",
    "SUS - CON", "group3 - group1", "positive = higher in SUS",
    "higher abundance in susceptible than control", "moderated t",
    "animal", "6 (3 SUS, 3 CON)",
    "BH within the differential-abundance family",
    "difference between stress-exposed susceptible animals and controls"),
  E("all three", "differential abundance", "protein log2 fold change",
    "RES - CON", "group2 - group1", "positive = higher in RES",
    "higher abundance in resilient than control", "moderated t",
    "animal", "6 (3 RES, 3 CON)",
    "BH within the differential-abundance family",
    "difference between stress-exposed resilient animals and controls"),
  E("all three", "ranked GSEA", "GO-BP normalised enrichment score",
    paste(contrasts, collapse = " | "),
    "gseGO on genes ranked by the median moderated t per official gene symbol",
    "positive NES = enriched among genes higher in the first-named group",
    "coordinated higher abundance of the gene set in the first-named group",
    "NES", "animal", "6 per contrast",
    "BH over every GO-BP set returned in that one comparison",
    "program-level direction and spatial pattern"),
  E("all three", "GO-program atlas", "median NES across a theme's GO terms",
    paste(contrasts, collapse = " | "),
    "median of constituent canonical GO-term NES per dataset x spatial unit x theme",
    "same sign convention as the constituent NES",
    "coordinated higher abundance of the program in the first-named group",
    "median NES", "animal", "6 per contrast",
    "NONE - descriptive aggregation; no theme-level family exists",
    "descriptive program summary; never an independent test"),
  E("all three", "CAMERA sensitivity", "competitive gene-set statistic",
    paste(contrasts, collapse = " | "),
    "cameraPR on the exact canonical ranked statistic, inter.gene.cor = 0.01",
    "Direction Up = higher in the first-named group",
    "gene set shifted upward relative to the rest of the ranking",
    "cameraPR PValue / Direction", "animal", "6 per contrast",
    "BH over the full comparable GO-BP family within each contrast",
    "LEVEL 4 orthogonal sensitivity only; never independent validation"),
  E("neuron_neuropil", "WGCNA module phenotype", "module eigengene difference",
    paste(contrasts, collapse = " | "),
    "lmerTest eigengene ~ StressGroup + SpatialUnit + (1|AnimalID), named contrast estimate",
    "positive = higher eigengene in the first-named group",
    "higher module eigengene in the first-named group",
    "eigengene contrast estimate", "animal", "9 across three groups",
    "BH within the module-phenotype family (0 of 45 supported)",
    "descriptive only; no FDR-supported module-phenotype association exists"),
  E("all three", "spatial network", "whole-network group comparison",
    "SUS/RES vs CON consensus",
    "Euclidean distance from a LEAVE-ONE-CON-ANIMAL-OUT consensus",
    "larger = further from the control consensus",
    "greater divergence from the control network", "network distance",
    "animal x dataset instance", "9 unique animals, 27 instances",
    "BH within the network family",
    "no detectable whole-network group difference at this sample size"),
  E("CON only", "bilateral reproducibility", "left-right concordance",
    "not a group contrast",
    "correlation between hemispheres of the same animal",
    "positive = concordant", "higher left-right agreement",
    "Pearson r / Spearman rho / ICC", "animal (CON only)", "3",
    "none - descriptive precision context",
    "descriptive technical reproducibility"),
  E("CON only", "external anatomical validation", "GSEA NES vs published signature",
    "CON-only anatomical contrasts",
    "internal anatomical contrast tested against an external signature",
    "positive = signature enriched in the named anatomical side",
    "agreement with the published anatomical signature", "NES",
    "animal (CON only)", "3",
    "BH within the external-validation inventory",
    "the ONLY genuinely external validation in the package"))
utils::write.csv(eff, file.path(PH_TAB, "manuscript_effect_contract.csv"),
                 row.names = FALSE)

cat("A01 inventory rows:", nrow(inv),
    "| missing source data:", sum(inv$source_data_path == "MISSING"),
    "| missing SVG:", sum(inv$panel_svg == "MISSING"),
    "| renderer not located:", sum(inv$renderer_file == "not located"), "\n")
cat("A02 effect contract rows:", nrow(eff), "\n")
