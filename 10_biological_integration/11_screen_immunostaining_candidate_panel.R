#!/usr/bin/env Rscript
# ================================================================
# Script: 10_biological_integration/11_screen_immunostaining_candidate_panel.R
# Stage: integration
# Scope: global (neuron_neuropil only)
# Consumes: required data/processed/02_id_mapping/mapped/neuron_neuropil/forward/per_file/<unit>sus_<unit>res.csv;
#   results/tables/11_spatial_systems/ca2_slm_robustness/CA2_SLM_DAP_robustness.csv;
#   results/tables/10_biological_integration/wgcna_candidate_protein_shortlist/neuron_neuropil/wgcna_candidate_proteins_shortlist.csv;
#   config/manuscript_spatial_order.yml.
# Produces: results/source_data/10_biological_integration/immunostaining_candidate_panel/.
# Notes: Candidate SCREEN. Selection only; no new statistics.
# ================================================================
#
# WHAT THIS SCRIPT DOES
#   Screens the canonical SUS-RES protein results for additional immunostaining
#   candidates beyond the three already nominated, and writes the selected panel
#   as frozen source data.
#
# WHY IT IS NOT A SIMPLE "TOP N BY FDR"
#   Of the 31 FDR-supported SUS-RES protein x unit cells in the whole neuropil,
#   28 are in CA2-SLM, and the CA2-SLM robustness audit qualified only 6 of
#   those 28. Ranking by FDR alone therefore returns mostly proteins the
#   repository already classifies as insufficient_observed_data or
#   not_claimable_due_to_QC. This screen consumes that existing classification
#   instead of re-deriving it, and records an explicit evidence class per
#   candidate.
#
# WHAT THIS SCRIPT DOES NOT DO
#   No model is refitted, no threshold is tuned, no statistic is recomputed. The
#   CA2-SLM robustness classes and the WGCNA candidate tiers are read as given.
#   Effect sizes are copied from the canonical files and verified against them.

suppressPackageStartupMessages({
  library(dplyr)
})

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

MODULE_ID <- "10_biological_integration"
SUBSTEP_ID <- "immunostaining_candidate_panel"
DATASET <- "neuron_neuropil"
CONTRACT_VERSION <- "immunostaining_candidate_panel_v1"
N_TARGET <- 10L

OUT <- path_results("source_data", MODULE_ID, SUBSTEP_ID)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ALREADY_NOMINATED <- c("OGA", "SLC22A23", "ANXA2")

DA_DIR <- repo_path("data", "processed", "02_id_mapping", "mapped", DATASET,
                    "forward", "per_file")
ROBUST <- repo_path("results", "tables", "11_spatial_systems",
                    "ca2_slm_robustness", "CA2_SLM_DAP_robustness.csv")
SHORT <- repo_path("results", "tables", "10_biological_integration",
                   "wgcna_candidate_protein_shortlist", DATASET,
                   "wgcna_candidate_proteins_shortlist.csv")
for (p in c(DA_DIR, ROBUST, SHORT))
  if (!file.exists(p) && !dir.exists(p))
    stop("missing_required_input: ", p, call. = FALSE)

# canonical anatomical order
units_yml <- yaml::read_yaml(repo_path("config", "manuscript_spatial_order.yml"))$units
npx <- Filter(function(z) identical(z$dataset, DATASET), units_yml)
SPATIAL <- data.frame(
  spatial_unit = vapply(npx, function(z) z$analysis_key, character(1)),
  spatial_unit_display = vapply(npx, function(z) z$display, character(1)),
  spatial_order = vapply(npx, function(z) as.integer(z$order), integer(1)),
  stringsAsFactors = FALSE)
SPATIAL <- SPATIAL[order(SPATIAL$spatial_order), ]

# --------------------------------------------- canonical SUS-RES, all units
da_file_for <- function(u) {
  k <- tolower(gsub("_", "", u))
  file.path(DA_DIR, sprintf("%ssus_%sres.csv", k, k))
}
da <- do.call(rbind, lapply(SPATIAL$spatial_unit, function(u) {
  d <- utils::read.csv(da_file_for(u), stringsAsFactors = FALSE)
  data.frame(
    spatial_unit = u, ProteinGroupID = d$ProteinGroupID,
    GeneSymbol = d$official_gene_symbol,
    UniProt = d$representative_accession,
    original_identifier = d$original_identifier,
    ambiguity_class = d$protein_group_ambiguity_class,
    protein_level_claim_allowed = d$protein_level_claim_allowed,
    log2FC = d$log2fc, p_value = d$pval, fdr_bh = d$padj,
    canonical_source_key = basename(da_file_for(u)),
    stringsAsFactors = FALSE)
}))
da <- da %>% inner_join(SPATIAL, by = "spatial_unit")

# per-protein summary across the ten contexts
summ <- da %>%
  group_by(ProteinGroupID, GeneSymbol, UniProt, original_identifier,
           ambiguity_class, protein_level_claim_allowed) %>%
  summarise(
    n_contexts = n(),
    n_fdr05 = sum(fdr_bh <= 0.05, na.rm = TRUE),
    n_fdr05_outside_ca2slm = sum(fdr_bh <= 0.05 & spatial_unit != "CA2_slm",
                                 na.rm = TRUE),
    min_fdr = min(fdr_bh, na.rm = TRUE),
    strongest_unit = spatial_unit[which.max(abs(log2FC))],
    strongest_log2FC = log2FC[which.max(abs(log2FC))],
    strongest_fdr = fdr_bh[which.max(abs(log2FC))],
    max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    median_abs_log2FC = stats::median(abs(log2FC), na.rm = TRUE),
    n_negative = sum(log2FC < 0, na.rm = TRUE),
    n_positive = sum(log2FC > 0, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(
    majority_direction = ifelse(n_negative >= n_positive, "negative", "positive"),
    n_matching_majority = pmax(n_negative, n_positive),
    directionally_consistent = n_matching_majority >= 8L)

# ------------------------------------- existing CA2-SLM robustness verdicts
rb <- utils::read.csv(ROBUST, stringsAsFactors = FALSE) %>%
  select(ProteinGroupID, CA2_SLM_robustness_class, classification_reason)
summ <- summ %>% left_join(rb, by = "ProteinGroupID")

# ------------------------------------------- existing WGCNA candidate tiers
sl <- utils::read.csv(SHORT, stringsAsFactors = FALSE) %>%
  transmute(ProteinGroupID, candidate_tier, ModuleID,
            module_display_label, abs_kME, clean_mapping) %>%
  distinct(ProteinGroupID, .keep_all = TRUE)
summ <- summ %>% left_join(sl, by = "ProteinGroupID")

# ------------------------------------------------------- evidence classing
# Three classes, in descending strength. A CA2-SLM hit only counts when the
# robustness audit qualified it; otherwise the protein is excluded from the
# panel and the reason is recorded.
summ <- summ %>%
  mutate(
    ca2_qualified = !is.na(CA2_SLM_robustness_class) &
      CA2_SLM_robustness_class == "robust_to_missingness_and_QC",
    evidence_class = case_when(
      n_fdr05_outside_ca2slm > 0 ~ "1_FDR_supported_outside_CA2SLM",
      n_fdr05 > 0 & ca2_qualified ~ "2_FDR_supported_CA2SLM_robustness_qualified",
      n_fdr05 > 0 & !ca2_qualified ~ "X_FDR_supported_CA2SLM_not_qualified",
      TRUE ~ "3_no_FDR_support_effect_and_consistency_only"),
    excluded_reason = ifelse(
      evidence_class == "X_FDR_supported_CA2SLM_not_qualified",
      paste0("only FDR-supported context is CA2-SLM and the robustness audit classified it ",
             CA2_SLM_robustness_class), NA_character_))

# A protein an antibody cannot unambiguously target is not a candidate.
eligible <- summ %>%
  filter(!toupper(GeneSymbol) %in% ALREADY_NOMINATED,
         !is.na(GeneSymbol), nzchar(GeneSymbol),
         as.logical(protein_level_claim_allowed),
         ambiguity_class == "single_accession_single_gene",
         evidence_class != "X_FDR_supported_CA2SLM_not_qualified")

# Epidermal / cornified-envelope proteins are excluded from the effect-size-only
# class. Ranking by raw effect size alone surfaces DSP, TGM1 and S100A14, and
# this neuropil matrix also carries FLG2, KRT26, KRT222, TGM2 and TGM3 - the
# classic skin-contamination signature for laser-captured brain tissue. The
# repository already treats this branch as QC rather than biology: the GO theme
# registry carries epithelial_epidermal_qc (epidermis development, skin
# development, keratinocyte differentiation) with theme_role = qc_review. A
# large, directionally consistent effect in one of these is far more likely to
# be variable dissection carry-over than a hippocampal phenotype, and it would
# be a poor immunostaining target. This is a screen-level QC exclusion, applied
# only to the effect-size-only class; it changes no statistic.
# The filter is ontology-based rather than a keyword blacklist. A keyword list
# is whack-a-mole: the first pass caught DSP/TGM1/S100A14 and then surfaced PKP1
# and NCCRP1, which are equally epithelial. The GO anchors below are the ones the
# repository's own registry already uses for the epithelial_epidermal_qc theme,
# plus desmosome and cornified envelope, so the exclusion is defined by
# annotation rather than by spelling.
EPIDERMAL_QC_GO <- c("GO:0008544",  # epidermis development
                     "GO:0030855",  # epithelial cell differentiation
                     "GO:0043588",  # skin development
                     "GO:0030057",  # desmosome
                     "GO:0001533")  # cornified envelope
epidermal_qc_symbols <- local({
  ok <- requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)
  if (!ok) {
    warning("org.Mm.eg.db unavailable; epidermal QC exclusion cannot be applied",
            call. = FALSE)
    return(character(0))
  }
  ids <- unique(unlist(suppressMessages(AnnotationDbi::mapIds(
    org.Mm.eg.db::org.Mm.eg.db, keys = EPIDERMAL_QC_GO, keytype = "GOALL",
    column = "SYMBOL", multiVals = "list")), use.names = FALSE))
  sort(unique(ids[!is.na(ids)]))
})
# GO alone is not sufficient: S100A14 and NCCRP1 carry essentially no mouse
# epidermal GO annotation yet are keratinocyte proteins, and both re-entered the
# ranking after the ontology filter. Rather than keep widening a pattern until
# the answer looks right - which would be tuning the rule against the result -
# the residue is handled by one short, stated list of epithelial/cornified
# proteins with no expected hippocampal expression. Both criteria are recorded
# per dropped protein so the decision is auditable.
EPIDERMAL_QC_LITERATURE <- c(
  "S100a14", "Nccrp1", "Pkp1", "Pkp3", "Dsp", "Dsg1a", "Dsc2", "Tgm1", "Tgm3",
  "Flg", "Flg2", "Lor", "Ivl", "Cdsn", "Evpl", "Ppl", "Sbsn", "Krt1", "Krt5",
  "Krt10", "Krt14", "Krt26", "Krt222", "Krtap", "Serpinb3a", "Casp14")
is_epidermal_qc <- function(sym)
  sym %in% epidermal_qc_symbols | sym %in% EPIDERMAL_QC_LITERATURE
epidermal_qc_basis <- function(sym)
  ifelse(sym %in% epidermal_qc_symbols & sym %in% EPIDERMAL_QC_LITERATURE,
         "GO_and_literature",
         ifelse(sym %in% epidermal_qc_symbols, "GO_annotation",
                ifelse(sym %in% EPIDERMAL_QC_LITERATURE, "literature_list", NA)))

# Class 3 mirrors how ANXA2 was nominated: no FDR support, but a large effect
# that points the same way across contexts, in a protein with high module
# membership. Requiring all three keeps it a stated, reproducible rule.
class3_pool <- eligible %>%
  filter(evidence_class == "3_no_FDR_support_effect_and_consistency_only",
         directionally_consistent,
         !is.na(abs_kME), abs_kME >= 0.60,
         !is_epidermal_qc(GeneSymbol)) %>%
  arrange(desc(max_abs_log2FC))

epidermal_dropped <- eligible %>%
  filter(evidence_class == "3_no_FDR_support_effect_and_consistency_only",
         directionally_consistent, !is.na(abs_kME), abs_kME >= 0.60,
         is_epidermal_qc(GeneSymbol)) %>%
  arrange(desc(max_abs_log2FC)) %>%
  mutate(exclusion_basis = epidermal_qc_basis(GeneSymbol)) %>%
  select(GeneSymbol, UniProt, ProteinGroupID, max_abs_log2FC,
         strongest_unit, n_matching_majority, abs_kME, exclusion_basis)

strong <- eligible %>%
  filter(evidence_class %in% c("1_FDR_supported_outside_CA2SLM",
                               "2_FDR_supported_CA2SLM_robustness_qualified")) %>%
  arrange(evidence_class, min_fdr)

n_fill <- max(0L, N_TARGET - nrow(strong))
panel <- bind_rows(strong, utils::head(class3_pool, n_fill)) %>%
  mutate(panel_rank = row_number(),
         selection_contract = CONTRACT_VERSION)

# ------------------------------------------------------ per-unit long table
effects <- da %>%
  semi_join(panel, by = "ProteinGroupID") %>%
  left_join(panel %>% select(ProteinGroupID, evidence_class, panel_rank),
            by = "ProteinGroupID") %>%
  mutate(contrast = "SUS - RES", fdr_supported = fdr_bh <= 0.05) %>%
  arrange(panel_rank, spatial_order)

# verification: copied effects must equal their canonical source exactly
for (u in unique(effects$spatial_unit)) {
  d <- utils::read.csv(da_file_for(u), stringsAsFactors = FALSE)
  e <- effects[effects$spatial_unit == u, , drop = FALSE]
  i <- match(e$ProteinGroupID, d$ProteinGroupID)
  if (!isTRUE(all.equal(e$log2FC, d$log2fc[i], tolerance = 0)) ||
      !isTRUE(all.equal(e$fdr_bh, d$padj[i], tolerance = 0)))
    stop("da_reproduction_failed at ", u, call. = FALSE)
}

# ------------------------------------------------------- animal abundance
# Same canonical bilateral animal-level matrix as the three-candidate figure.
# The GCT keys rows by original identifier, and ExpGroup is numeric, so the
# coding is asserted against the known rosters rather than assumed.
GCT <- repo_path("data", "processed", "01_preprocessing",
                 "protigy_input_animal_level", DATASET,
                 paste0(DATASET, "_animal_level.gct"))
gl <- readLines(GCT, n = 12L)
dims <- as.integer(strsplit(gl[2], "\t")[[1]])
n_col_meta <- dims[4]
hdr <- strsplit(gl[3], "\t")[[1]]
mrows <- lapply(seq_len(n_col_meta), function(i) strsplit(gl[3L + i], "\t")[[1]])
mnames <- vapply(mrows, function(x) x[1], character(1))
mval <- function(k) mrows[[which(mnames == k)]][-(1:2)]
smeta <- data.frame(sample_id = hdr[-(1:2)], AnimalID = mval("AnimalID"),
                    ExpGroupCode = mval("ExpGroup"), region = mval("region"),
                    layer = mval("layer"), stringsAsFactors = FALSE)
smeta$ExpGroup <- unname(c(`1` = "CON", `2` = "RES", `3` = "SUS")[smeta$ExpGroupCode])
smeta$spatial_unit <- paste(smeta$region, smeta$layer, sep = "_")
ROSTER <- list(CON = c("A127", "A129", "A765"), RES = c("A0003", "A135", "A139"),
               SUS = c("A111", "A755", "A764"))
for (g in names(ROSTER))
  if (!identical(sort(unique(smeta$AnimalID[smeta$ExpGroup == g])),
                 sort(ROSTER[[g]])))
    stop("ExpGroup coding mismatch for ", g, call. = FALSE)

ex <- utils::read.delim(GCT, skip = 3L + n_col_meta, header = FALSE,
                        stringsAsFactors = FALSE, check.names = FALSE)
names(ex) <- hdr
ridx <- match(panel$original_identifier, ex[[1]])
if (anyNA(ridx))
  stop("panel protein absent from the animal-level matrix: ",
       paste(panel$GeneSymbol[is.na(ridx)], collapse = ", "), call. = FALSE)

abundance <- do.call(rbind, lapply(seq_len(nrow(panel)), function(i)
  data.frame(
    ProteinGroupID = panel$ProteinGroupID[i], GeneSymbol = panel$GeneSymbol[i],
    UniProt = panel$UniProt[i], panel_rank = panel$panel_rank[i],
    evidence_class = panel$evidence_class[i],
    sample_id = smeta$sample_id, AnimalID = smeta$AnimalID,
    ExpGroup = smeta$ExpGroup, spatial_unit = smeta$spatial_unit,
    abundance = as.numeric(unlist(ex[ridx[i], smeta$sample_id])),
    stringsAsFactors = FALSE))) %>%
  inner_join(SPATIAL, by = "spatial_unit") %>%
  mutate(dataset = DATASET, biological_unit = "animal",
         hemisphere_handling = "bilateral animal-level; equal-weight L/R mean where both present, single hemisphere retained otherwise, no imputation") %>%
  arrange(panel_rank, spatial_order, ExpGroup, AnimalID)

if (nrow(abundance %>% count(ProteinGroupID, AnimalID, spatial_unit) %>%
         filter(n > 1L)))
  stop("pseudoreplication in the panel abundance table", call. = FALSE)
utils::write.csv(abundance,
  file.path(OUT, "candidate_panel_animal_abundance.csv"), row.names = FALSE)

git_sha <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
  error = function(e) NA_character_)

utils::write.csv(panel, file.path(OUT, "candidate_panel_selection.csv"),
                 row.names = FALSE)
utils::write.csv(effects, file.path(OUT, "candidate_panel_sus_res_effects.csv"),
                 row.names = FALSE)
utils::write.csv(
  summ %>% filter(evidence_class == "X_FDR_supported_CA2SLM_not_qualified") %>%
    select(GeneSymbol, UniProt, ProteinGroupID, candidate_tier, min_fdr,
           strongest_unit, strongest_log2FC, CA2_SLM_robustness_class,
           excluded_reason) %>% arrange(min_fdr),
  file.path(OUT, "candidate_panel_excluded_ca2slm.csv"), row.names = FALSE)
utils::write.csv(epidermal_dropped,
  file.path(OUT, "candidate_panel_excluded_epidermal_qc.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(repository_commit = git_sha, contract_version = CONTRACT_VERSION,
             dataset = DATASET, n_target = N_TARGET,
             canonical_da_input = "data/processed/02_id_mapping/mapped/neuron_neuropil/forward/per_file/",
             ca2_robustness_input = "results/tables/11_spatial_systems/ca2_slm_robustness/CA2_SLM_DAP_robustness.csv",
             wgcna_shortlist_input = "results/tables/10_biological_integration/wgcna_candidate_protein_shortlist/neuron_neuropil/wgcna_candidate_proteins_shortlist.csv",
             spatial_order_contract = "config/manuscript_spatial_order.yml",
             biological_unit = "animal", stringsAsFactors = FALSE),
  file.path(OUT, "candidate_panel_provenance.csv"), row.names = FALSE)

cat("candidate panel ->", OUT, "\n")
cat("FDR-supported protein x unit cells in the whole neuropil:",
    sum(da$fdr_bh <= 0.05, na.rm = TRUE), "| in CA2-SLM:",
    sum(da$fdr_bh <= 0.05 & da$spatial_unit == "CA2_slm", na.rm = TRUE), "\n")
cat("excluded as CA2-SLM not-qualified:",
    sum(summ$evidence_class == "X_FDR_supported_CA2SLM_not_qualified"), "\n\n")
cat("selected panel:", nrow(panel), "\n")
print(panel %>%
        transmute(rank = panel_rank, GeneSymbol, UniProt,
                  class = substr(evidence_class, 1, 1),
                  n_fdr05, strongest_unit,
                  log2FC = round(strongest_log2FC, 3),
                  fdr = signif(min_fdr, 3),
                  dir = paste0(n_matching_majority, "/10 ", majority_direction),
                  kME = round(abs_kME, 2), tier = candidate_tier) %>%
        as.data.frame(), row.names = FALSE)
print(table(panel$evidence_class))
