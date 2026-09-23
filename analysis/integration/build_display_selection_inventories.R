#!/usr/bin/env Rscript
# ================================================================
# Script: analysis/integration/build_display_selection_inventories.R
# Stage: integration
# Scope: global
# Consumes: required results/tables/10_biological_integration/gsea_wgcna_concordance/global/ontology_aware_gsea_theme_assignments_all_contrasts.csv
# Produces: results/integration/build_display_selection_inventories/global/tables/pathway_enrichment_inventory.csv; results/integration/build_display_selection_inventories/global/tables/pathway_enrichment_inventory_fdr_supported.csv; results/integration/build_display_selection_inventories/global/tables/leading_edge_protein_inventory.csv; results/integration/build_display_selection_inventories/global/tables/leading_edge_protein_recurrence.csv; results/integration/build_display_selection_inventories/global/tables/display_selection_disclosure.csv; results/integration/build_display_selection_inventories/global/tables/inventory_data_dictionary.csv
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Releases the denominators behind the selected pathways and proteins shown in Figures 2 and 3.
# ================================================================
#
# WHY THIS EXISTS
#
# Figure 3 shows seven atlas rows, three exemplar GSEA curves and 21 leading-
# edge proteins. Each of those is a SELECTION, and until now the set it was
# selected from was not released. A reader could see the numerator and not the
# denominator, which makes an honest selection indistinguishable from a mined
# one.
#
# The package already applies the opposite standard elsewhere: ST2 releases all
# 30 external pairings "so specificity can be judged rather than assumed", and
# ST8 releases every tested edge-behaviour correlation "so no subset can be
# mined". This script extends that standard to the pathway and protein
# evidence behind Figure 3.
#
# WHAT THIS SCRIPT DOES NOT DO
#   * No GSEA is rerun. NES, raw p and BH FDR are copied unchanged from the
#     frozen theme-assignment table.
#   * No enrichment test, model or multiple-testing correction is computed.
#   * No figure, no panel, no manuscript claim is modified.
#   * No new statistic is invented. Every column is either copied, or a COUNT
#     of copied rows, or a flag derived from a written rule.
#
# The one derived quantity is the recurrence classification, which applies the
# leading-edge protein rule already written in
# docs/ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md. Counting is not inference:
# no p-value, FDR or effect estimate is produced for a protein here, and the
# standing caveat that no leading-edge protein is individually FDR-supported is
# carried once in inventory_data_dictionary.csv.
# ================================================================

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))

ANALYSIS_ID <- "build_display_selection_inventories"
SCRIPT_ID <- "analysis/integration/build_display_selection_inventories.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

run <- integration_cli(default_dataset = "all", allow_all = TRUE)
paths <- integration_dirs(ANALYSIS_ID, "global", create = !isTRUE(run$dry_run))

inputs <- list(
  themes = integration_find(
    "ontology_aware_gsea_theme_assignments_all_contrasts.csv",
    owner = "test_enrichment_module_concordance",
    legacy_stage = "10_biological_integration",
    legacy_substep = "gsea_wgcna_concordance")
)

if (isTRUE(run$dry_run)) {
  dry_run_inputs(SCRIPT_ID, inputs)
  quit(status = 0, save = "no")
}

loaded <- read_csv_optional(inputs$themes, "global", "integration",
                            "theme_assignments", required = TRUE)
th <- loaded$data
if (is.null(th) || !nrow(th)) {
  stop("missing_required_input: theme assignment table is empty: ",
       inputs$themes, call. = FALSE)
}

FDR_ALPHA <- 0.05

# ---------------------------------------------------------------- constants
#
# The three Figure 3 exemplars. This is a LITERAL copy of the renderer's
# s4_programs() table in Exp9_manuscript, reproduced here only so the
# disclosure table can state where each exemplar sits in the inventory. It
# selects nothing; changing it would change a label, not a result.
EXEMPLARS <- data.frame(
  program = c("synaptic signalling", "mRNA processing",
              "oxidative phosphorylation"),
  dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
  spatial_unit = c("CA3_sr", "CA2_sp", "CA1"),
  GO_ID = c("GO:0099536", "GO:0006397", "GO:0006119"),
  stringsAsFactors = FALSE)

# The seven claim-eligible atlas rows, in registry display order.
ATLAS_ROWS <- c("rna_processing_splicing_rnp", "ribosome_translation",
                "chromatin_organization", "mitochondrial_respiration_oxphos",
                "synaptic_signaling_vesicle", "neuron_projection_development",
                "autophagy_lysosome_endosome")

compartment_label <- function(dataset) {
  unname(c(neuron_neuropil = "Neuropil", neuron_soma = "Neuronal soma",
           microglia = "Microglia-enriched ROI")[dataset])
}

# collapse a character vector to a stable, deduplicated, semicolon list
collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) return(NA_character_)
  paste(sort(x), collapse = ";")
}

# ============================================================== INVENTORY 1
# Every GO term tested in every spatial unit and contrast.
#
# The stored table carries one row per THEME ASSIGNMENT, so the 180 terms that
# legitimately sit in two themes appear twice. The inventory is keyed on the
# term occurrence instead, and the theme columns are collapsed, so a row here
# is one GSEA result and the row count is the true denominator.
key <- paste(th$dataset, th$spatial_unit, th$contrast, th$GO_ID, sep = "\r")
ord <- order(key)
th <- th[ord, , drop = FALSE]
key <- key[ord]
grp <- !duplicated(key)

inv <- data.frame(
  dataset = th$dataset[grp],
  compartment = compartment_label(th$dataset[grp]),
  spatial_unit = th$spatial_unit[grp],
  contrast = th$contrast[grp],
  GO_ID = th$GO_ID[grp],
  GO_description = th$GO_description[grp],
  NES = th$NES[grp],
  raw_p = th$raw_p[grp],
  BH_FDR = th$GSEA_FDR[grp],
  stringsAsFactors = FALSE)
inv$fdr_supported <- is.finite(inv$BH_FDR) & inv$BH_FDR < FDR_ALPHA

split_key <- split(seq_len(nrow(th)), key)
split_key <- split_key[match(key[grp], names(split_key))]
inv$theme_id <- vapply(split_key, function(ix) collapse_unique(th$theme_id[ix]),
                       character(1))
inv$manuscript_theme <- vapply(split_key,
                               function(ix) collapse_unique(th$manuscript_theme[ix]),
                               character(1))
inv$theme_role <- vapply(split_key, function(ix) collapse_unique(th$theme_role[ix]),
                         character(1))
inv$theme_claim_eligible <- vapply(
  split_key, function(ix) any(th$theme_claim_eligible[ix] %in% TRUE), logical(1))
inv$assignment_status <- vapply(split_key,
                                function(ix) collapse_unique(th$assignment_status[ix]),
                                character(1))
rownames(inv) <- NULL

# what the figure actually shows
#
# Restricted to the contrast the panel draws. This was theme_claim_eligible
# alone, which is contrast-blind: it marked all 12,598 claim-eligible
# occurrences (4,200 RES - CON, 4,199 SUS - CON, 4,199 SUS - RES) as
# contributing to a drawn cell, so 8,399 of them - two thirds - claimed to
# contribute to a panel they are not in. docs/FIGURE_SELECTION_RULES.md states
# the panel's own figure: "4,199 of them for the SUS - RES contrast the panel
# draws". theme_claim_eligible remains available, unchanged, for the
# denominator; this column now means what its name says.
FIG3B_CONTRAST <- "SUS - RES"
inv$contributes_to_figure_3b_cell <- inv$theme_claim_eligible &
  inv$contrast == FIG3B_CONTRAST
inv$is_figure_3_exemplar_term <- paste(inv$dataset, inv$spatial_unit, inv$GO_ID) %in%
  paste(EXEMPLARS$dataset, EXEMPLARS$spatial_unit, EXEMPLARS$GO_ID)

inv <- inv[order(inv$contrast, inv$dataset, inv$spatial_unit, -abs(inv$NES)), ,
           drop = FALSE]
rownames(inv) <- NULL

invisible(write_integration_table(inv, paths, "pathway_enrichment_inventory.csv"))
message(sprintf("[inventory] %-46s %7d rows", "pathway_enrichment_inventory.csv",
                nrow(inv)))

sup <- inv[inv$fdr_supported, , drop = FALSE]
rownames(sup) <- NULL
invisible(write_integration_table(
  sup, paths, "pathway_enrichment_inventory_fdr_supported.csv"))
message(sprintf("[inventory] %-46s %7d rows",
                "pathway_enrichment_inventory_fdr_supported.csv", nrow(sup)))

# ============================================================== INVENTORY 2
# Every leading-edge protein of every FDR-supported term.
#
# Figure 3 g/h/i show 7 proteins per program. This is the set those 7 were
# drawn from, for every supported term in the dataset rather than only the
# three exemplars.
le_src <- th[!duplicated(key), , drop = FALSE]
le_src <- le_src[is.finite(le_src$GSEA_FDR) & le_src$GSEA_FDR < FDR_ALPHA, ,
                 drop = FALSE]
genes <- strsplit(as.character(le_src$leading_edge_genes), ";", fixed = TRUE)
n_le <- lengths(genes)
keep <- n_le > 0

le <- data.frame(
  dataset = rep(le_src$dataset[keep], n_le[keep]),
  compartment = rep(compartment_label(le_src$dataset[keep]), n_le[keep]),
  spatial_unit = rep(le_src$spatial_unit[keep], n_le[keep]),
  contrast = rep(le_src$contrast[keep], n_le[keep]),
  GO_ID = rep(le_src$GO_ID[keep], n_le[keep]),
  GO_description = rep(le_src$GO_description[keep], n_le[keep]),
  NES = rep(le_src$NES[keep], n_le[keep]),
  BH_FDR = rep(le_src$GSEA_FDR[keep], n_le[keep]),
  theme_id = rep(le_src$theme_id[keep], n_le[keep]),
  manuscript_theme = rep(le_src$manuscript_theme[keep], n_le[keep]),
  theme_claim_eligible = rep(le_src$theme_claim_eligible[keep] %in% TRUE, n_le[keep]),
  gene = trimws(unlist(genes[keep], use.names = FALSE)),
  leading_edge_size_of_term = rep(n_le[keep], n_le[keep]),
  stringsAsFactors = FALSE)
le <- le[nzchar(le$gene), , drop = FALSE]
le$is_figure_3_exemplar_term <- paste(le$dataset, le$spatial_unit, le$GO_ID) %in%
  paste(EXEMPLARS$dataset, EXEMPLARS$spatial_unit, EXEMPLARS$GO_ID)
# The caveat that no protein here is individually FDR-supported belongs in the
# data dictionary, not repeated on 282,296 rows. Carrying prose in a data
# column is the machine-oriented habit the supplementary-table layer exists to
# undo, and here it would have tripled the file.
le <- le[order(le$contrast, le$dataset, le$spatial_unit, le$GO_ID, le$gene), ,
         drop = FALSE]
rownames(le) <- NULL

invisible(write_integration_table(le, paths, "leading_edge_protein_inventory.csv"))
message(sprintf("[inventory] %-46s %7d rows", "leading_edge_protein_inventory.csv",
                nrow(le)))

# ============================================================== INVENTORY 3
# The recurrence classification from the atlas naming rules.
#
# docs/ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md, leading-edge protein rule:
#
#   "A protein counts as recurrent core only if it appears in >=3 supported GO
#    terms AND >=3 spatial contexts. Everything else is recorded separately as
#    INTERMEDIATE or SINGLE_APPEARANCE."
#
# The rule is written but was never computed into an artefact; it has been used
# only as a naming sanity check. It is applied here verbatim.
#
# The rule names three classes and defines the boundary of one. SINGLE_APPEARANCE
# is taken to mean exactly one supported term in exactly one spatial context,
# and INTERMEDIATE everything between. That reading is stated on the table
# rather than left implicit, because the source document does not fix it.
#
# Counted within a contrast, over CLAIM-ELIGIBLE supported terms, which is the
# atlas scope the rule was written for. The counts over all supported terms
# regardless of theme are carried alongside so the narrower basis is visible.
RECUR_MIN_TERMS <- 3L
RECUR_MIN_CONTEXTS <- 3L

recurrence_for <- function(z, basis) {
  if (!nrow(z)) return(NULL)
  g <- paste(z$contrast, z$gene, sep = "\r")
  do.call(rbind, lapply(split(seq_len(nrow(z)), g), function(ix) {
    w <- z[ix, , drop = FALSE]
    data.frame(
      contrast = w$contrast[1], gene = w$gene[1],
      n_supported_go_terms = length(unique(w$GO_ID)),
      n_spatial_contexts = length(unique(paste(w$dataset, w$spatial_unit))),
      n_compartments = length(unique(w$dataset)),
      spatial_contexts = collapse_unique(paste0(w$dataset, ":", w$spatial_unit)),
      themes = collapse_unique(w$manuscript_theme),
      basis = basis, stringsAsFactors = FALSE)
  }))
}

rec_atlas <- recurrence_for(le[le$theme_claim_eligible, , drop = FALSE],
                            "claim_eligible_supported_terms")
rec_all <- recurrence_for(le, "all_supported_terms")

rec <- rec_atlas
rec$recurrence_class <- ifelse(
  rec$n_supported_go_terms >= RECUR_MIN_TERMS &
    rec$n_spatial_contexts >= RECUR_MIN_CONTEXTS, "RECURRENT_CORE",
  ifelse(rec$n_supported_go_terms == 1L & rec$n_spatial_contexts == 1L,
         "SINGLE_APPEARANCE", "INTERMEDIATE"))

m <- match(paste(rec$contrast, rec$gene), paste(rec_all$contrast, rec_all$gene))
rec$n_supported_go_terms_any_theme <- rec_all$n_supported_go_terms[m]
rec$n_spatial_contexts_any_theme <- rec_all$n_spatial_contexts[m]

# genes that carry a supported term outside the seven atlas rows only
only_all <- rec_all[!paste(rec_all$contrast, rec_all$gene) %in%
                      paste(rec$contrast, rec$gene), , drop = FALSE]
if (nrow(only_all)) {
  only_all$recurrence_class <- "NOT_IN_A_CLAIM_ELIGIBLE_THEME"
  only_all$n_supported_go_terms_any_theme <- only_all$n_supported_go_terms
  only_all$n_spatial_contexts_any_theme <- only_all$n_spatial_contexts
  only_all$n_supported_go_terms <- 0L
  only_all$n_spatial_contexts <- 0L
  rec <- rbind(rec, only_all[, names(rec), drop = FALSE])
}

rec <- rec[order(rec$contrast, -rec$n_supported_go_terms,
                 -rec$n_spatial_contexts, rec$gene), , drop = FALSE]
rownames(rec) <- NULL

invisible(write_integration_table(rec, paths, "leading_edge_protein_recurrence.csv"))
message(sprintf("[inventory] %-46s %7d rows", "leading_edge_protein_recurrence.csv",
                nrow(rec)))
for (ct in sort(unique(rec$contrast))) {
  z <- rec[rec$contrast == ct, , drop = FALSE]
  message(sprintf("             %-11s core %4d  intermediate %5d  single %5d",
                  ct, sum(z$recurrence_class == "RECURRENT_CORE"),
                  sum(z$recurrence_class == "INTERMEDIATE"),
                  sum(z$recurrence_class == "SINGLE_APPEARANCE")))
}

# ============================================================== DISCLOSURE
# One row per displayed selection in Figures 2 and 3: the rule, the set it
# selected from, and whether the rule looked at the result it displays.
#
# The Figure 2 counts are properties of that figure's own source tables and are
# stated as verified literals. The Figure 3 exemplar ranks are computed here
# from the inventory above.
exemplar_rank <- function(i) {
  z <- inv[inv$contrast == "SUS - RES" & inv$dataset == EXEMPLARS$dataset[i] &
             inv$fdr_supported & inv$theme_claim_eligible, , drop = FALSE]
  zu <- z[z$spatial_unit == EXEMPLARS$spatial_unit[i], , drop = FALSE]
  zu <- zu[order(-abs(zu$NES)), , drop = FALSE]
  k <- which(zu$GO_ID == EXEMPLARS$GO_ID[i])[1]
  sprintf("%d of %d FDR-supported claim-eligible terms in %s",
          k, nrow(zu), EXEMPLARS$spatial_unit[i])
}
ex_ranks <- vapply(seq_len(nrow(EXEMPLARS)), exemplar_rank, character(1))

D <- function(panel, displayed, n_displayed, rule, universe, n_universe,
              outcome_dependent, caveat) {
  data.frame(panel = panel, displayed = displayed, n_displayed = n_displayed,
             selection_rule = rule, selected_from = universe,
             universe_size = n_universe,
             rule_used_the_displayed_result = outcome_dependent,
             caveat = caveat, stringsAsFactors = FALSE)
}

disc <- rbind(
  D("Figure 2d", "spatial fingerprint proteins", "19 unique genes",
    paste0("The top 2 genes by BH-adjusted p on the enriched side (logFC > 0) ",
           "of each prespecified CON-only anatomical contrast. Ties break on ",
           "-|logFC| then gene symbol, so the pick is deterministic. Rows are ",
           "then ordered by baseline CON peak unit."),
    "11 prespecified CON-only anatomical contrasts, fitted on CON animals only",
    "11 contrasts x 2 = 22 picks, 19 unique genes",
    "no - the contrast set was fixed before any stress contrast was computed, and no stress group enters selection or row order",
    "Descriptive characterisation of baseline spatial structure, not validation."),
  D("Figure 2e", "compartment marker proteins", "10 markers",
    paste0("A written list of canonical compartment markers, 3-4 per ",
           "compartment, each of which must then pass primary AND strict ",
           "detection eligibility in its intended compartment or the build ",
           "fails."),
    "curated canonical markers, gated by the detection-eligibility audit",
    "10 requested, 10 eligible",
    "no - selection is recorded as not using observed cross-compartment direction or effect magnitude; expected direction is post-selection validation only",
    "Chosen for recognisability, not for performance."),
  D("Figure 2g", "external signature pairings", "10 pairings",
    paste0("Every pairing declared an expected anatomical correspondence ",
           "(CA1 to CA1, CA2/3 to CA2/3, target stratum to its own signature). ",
           "Expectation is a property of anatomy, written before results."),
    "all 30 tested internal-contrast by external-signature pairings",
    "10 of 30 expected; all 30 released in ST2",
    "no - expectation is declared from anatomy, not from the result",
    "The one externally anchored panel in the package."),
  D("Figure 2h", "internal anatomical GO terms", "7 terms",
    "One GO term per anatomical contrast: the maximum |NES|, ties broken on term description.",
    "the canonical GO BP GSEA of each anatomical contrast",
    "7 contrasts, 1 term each",
    "yes - the term shown is the maximum of the quantity shown, within a prespecified contrast",
    "Characterisation of the same proteomic data, explicitly NOT independent validation."),
  D("Figure 3b", "atlas program rows", "7 themes",
    paste0("Themes flagged claim-eligible in a version-controlled registry, ",
           "drawn in registry display order. Admission requires all five of ",
           "biological coherence, recurrence across spatial contexts, ",
           "distinctness from ontology ancestor ladders, non-QC status and ",
           "material representation. Selection must never use SUS/RES ",
           "direction, NES magnitude or the leading-edge proteins."),
    "all GO BP terms assigned to a registry theme",
    sprintf("%d claim-eligible term-occurrences of %d tested",
            sum(inv$theme_claim_eligible), nrow(inv)),
    "no - ordering and admission are phenotype-independent; the rows are never sorted by NES, FDR, direction or support count",
    "Descriptive aggregation. No theme-level p-value or FDR exists."),
  D("Figure 3 d/e/f", "exemplar GSEA programs", "3 programs",
    paste0("EDITORIAL. One FDR-supported program per compartment, fixed as a ",
           "literal table, chosen to represent distinct biology in each of the ",
           "three compartments. Not an algorithmic maximum. Chosen after the ",
           "phenotype-contrast results were known."),
    "claim-eligible themes FDR-supported for SUS - RES within each compartment",
    paste(ex_ranks, collapse = " | "),
    "yes - chosen after seeing the results, by judgement rather than by a rule",
    paste0("Illustrative, never 'strongest'. The complete evidence they ",
           "illustrate is the atlas and the inventories released here.")),
  D("Figure 3 g/h/i", "leading-edge proteins", "21 proteins (7 per program)",
    paste0("The leading edge of the exemplar's own enrichment, ranked by ",
           "|stored rank statistic|, top 7. The cut at 7 is a display ",
           "constraint chosen for legibility, not a statistical threshold."),
    "the complete leading edge of the three exemplar terms, for SUS - RES",
    ## Restricted to the panel's own contrast. is_figure_3_exemplar_term is
    ## keyed on dataset + spatial_unit + GO_ID and deliberately ignores
    ## contrast, so summing it counts the same three terms once per contrast:
    ## 36 (RES - CON) + 416 (SUS - CON) + 415 (SUS - RES) = 867. Figure 3 g/h/i
    ## draws SUS - RES, so 867 overstated the panel's universe by 2.09x and
    ## contradicted docs/FIGURE_SELECTION_RULES.md, which reports the correct
    ## 255 + 112 + 48 = 415. Both files ship, so the disagreement shipped too.
    sprintf("%d leading-edge memberships across the 3 exemplar terms (SUS - RES)",
            sum(le$is_figure_3_exemplar_term & le$contrast == "SUS - RES")),
    "yes - both the parent program and the within-program ranking are outcome-dependent",
    paste0("No protein is individually FDR-supported (smallest BH FDR 0.53). ",
           "Not validated, not significant, not drivers.")))

## What is ACTUALLY released, checked against the export bundle rather than
## asserted. Two entries were overstated and are corrected here.
disc$complete_inventory_released_as <- c(
  "ST1 / v9_ed_fingerprint_full", "figure2d marker provenance table",
  "ST2_external_signature_validation.csv",
  ## Figure 2h: this table is NOT in the publication bundle - 0 rows in
  ## exports/publication_source_data/manifest.csv and 0 files under exports/
  ## or pride_submission/. It exists only as an internal results table, so
  ## claiming a released complete inventory for this panel was wrong.
  paste0("NOT RELEASED - internal only: results/tables/",
         "04_differential_expression_enrichment/control_spatial_identity_validation/",
         "global/control_anatomical_go_bp_gsea.csv"),
  "pathway_enrichment_inventory.csv", "pathway_enrichment_inventory.csv",
  ## Figure 3 g/h/i: the released inventories establish MEMBERSHIP and the
  ## universe size, but not the ordering. The panel's rule ranks by the stored
  ## per-protein rank statistic, and that statistic is not carried here: the
  ## upstream theme table supplies leading_edge_genes as an ALPHABETICALLY
  ## sorted string, so neither the value nor its order survives. A reader can
  ## therefore confirm which proteins were eligible and how many there were,
  ## but cannot reproduce which 7 were drawn.
  paste0("leading_edge_protein_inventory.csv + leading_edge_protein_recurrence.csv ",
         "(membership and universe size only; the per-protein rank statistic that ",
         "orders the top-7 cut is not carried - upstream supplies leading-edge ",
         "genes alphabetically, so the displayed ordering cannot be reproduced ",
         "from these tables)"))

invisible(write_integration_table(disc, paths, "display_selection_disclosure.csv"))
message(sprintf("[inventory] %-46s %7d rows", "display_selection_disclosure.csv",
                nrow(disc)))

# ============================================================== DICTIONARY
# Column definitions and the standing caveats, carried once instead of on
# every row of a 282,296-row table.
DD <- function(table_file, column, definition) {
  data.frame(table_file = table_file, column = column, definition = definition,
             stringsAsFactors = FALSE)
}
SHARED <- c(
  dataset = "Acquisition compartment: neuron_neuropil, neuron_soma or microglia.",
  compartment = "Reader-facing label for dataset.",
  spatial_unit = "Prespecified hippocampal spatial unit. Neuropil is region x layer; soma and the microglia-enriched ROI are region-level only.",
  contrast = "Stress-group contrast. The biological replicate is the animal, n = 3 per group.",
  GO_ID = "Gene Ontology biological-process identifier.",
  GO_description = "GO term name as stored in the enrichment result.",
  NES = "Normalised enrichment score, copied unchanged from the stored GSEA result.",
  raw_p = "Permutation p-value, floored at eps = 1e-10 by the upstream method.",
  BH_FDR = "Benjamini-Hochberg FDR within the GSEA family for that dataset, spatial unit and contrast.",
  fdr_supported = "TRUE when BH_FDR < 0.05. Absence of support is not absence of effect at n = 3 per group.",
  theme_id = "Registry identifier of the atlas theme this term is assigned to. Semicolon-separated when a term legitimately sits in two themes.",
  manuscript_theme = "Reader-facing theme label.",
  theme_role = "primary for the seven atlas rows, qc_review for the two excluded technical themes.",
  theme_claim_eligible = "TRUE when the term belongs to one of the seven claim-eligible atlas rows drawn in Figure 3b.",
  assignment_status = "single_theme, multi_theme, qc_review or unclassified.",
  contributes_to_figure_3b_cell = "TRUE when this term contributes to a drawn Figure 3b cell: claim-eligible AND in the SUS - RES contrast the panel draws. The cell value is the median NES of its contributing terms. Use theme_claim_eligible for the contrast-independent denominator (12,598 occurrences); this column is the 4,199 that reach the panel.",
  is_figure_3_exemplar_term = "TRUE for the three GO terms shown as exemplar curves in Figure 3 d/e/f.",
  gene = "Official mouse gene symbol of a leading-edge member of that enrichment.",
  leading_edge_size_of_term = "Number of leading-edge genes in that term, spatial unit and contrast. Figure 3 g/h/i display 7 of these.")

dd <- rbind(
  do.call(rbind, lapply(names(inv), function(c_)
    DD("pathway_enrichment_inventory.csv", c_, SHARED[[c_]]))),
  do.call(rbind, lapply(names(sup), function(c_)
    DD("pathway_enrichment_inventory_fdr_supported.csv", c_, SHARED[[c_]]))),
  do.call(rbind, lapply(names(le), function(c_)
    DD("leading_edge_protein_inventory.csv", c_, SHARED[[c_]]))),
  DD("leading_edge_protein_recurrence.csv",
     c("contrast", "gene", "n_supported_go_terms", "n_spatial_contexts",
       "n_compartments", "spatial_contexts", "themes", "basis",
       "recurrence_class", "n_supported_go_terms_any_theme",
       "n_spatial_contexts_any_theme"),
     c(SHARED[["contrast"]], SHARED[["gene"]],
       "Number of distinct claim-eligible FDR-supported GO terms whose leading edge contains this gene, within the contrast.",
       "Number of distinct dataset x spatial_unit combinations in which that holds.",
       "Number of distinct compartments in which that holds.",
       "The spatial contexts themselves, as dataset:spatial_unit.",
       "The atlas themes involved.",
       "The term set the counts were taken over.",
       paste0("RECURRENT_CORE requires >= ", RECUR_MIN_TERMS,
              " supported GO terms AND >= ", RECUR_MIN_CONTEXTS,
              " spatial contexts, the rule written in docs/ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md.",
              " SINGLE_APPEARANCE is one term in one context; INTERMEDIATE is everything between;",
              " NOT_IN_A_CLAIM_ELIGIBLE_THEME means the gene's supported terms all fall outside the seven atlas rows.",
              " The source document names the three classes and fixes only the RECURRENT_CORE boundary;",
              " the other two boundaries are this table's stated reading."),
       "The same count taken over all supported terms regardless of theme.",
       "The same context count taken over all supported terms regardless of theme.")),
  DD("display_selection_disclosure.csv",
     c("panel", "displayed", "n_displayed", "selection_rule", "selected_from",
       "universe_size", "rule_used_the_displayed_result", "caveat",
       "complete_inventory_released_as"),
     c("Manuscript figure panel.", "What that panel shows.",
       "How many items are drawn.",
       "The rule that chose them, stated so it could be re-executed.",
       "The set the rule chose from.", "How large that set is.",
       "Whether the rule looked at the quantity the panel displays. This is the question a reader is entitled to ask and the reason this table exists.",
       "The standing interpretive limit on that panel.",
       "Where the complete set is released.")))
dd <- dd[!is.na(dd$definition), , drop = FALSE]

NOTES <- data.frame(
  table_file = "ALL",
  column = c("standing_caveat_leading_edge", "standing_caveat_theme_level",
             "standing_caveat_absence", "standing_caveat_eps_floor"),
  definition = c(
    "Leading-edge membership is not individual protein significance. No protein in these tables is individually FDR-supported by the differential-abundance analysis; BH_FDR is the FDR of the ENRICHMENT, not of the protein. Never describe these as validated, significant, key or driver proteins.",
    "No theme-level p-value or FDR exists. Atlas themes are descriptive aggregations of canonical GO terms and carry no multiple-testing family of their own.",
    "Absence of FDR support is not absence of effect. At three animals per group write 'did not survive correction', never 'no difference' or 'unchanged'.",
    "GSEA p-values are floored at eps = 1e-10. A term at the floor has a true p the method does not resolve, so its FDR bounds the evidence rather than measuring it."),
  stringsAsFactors = FALSE)

invisible(write_integration_table(rbind(dd, NOTES), paths,
                                  "inventory_data_dictionary.csv"))
message(sprintf("[inventory] %-46s %7d rows", "inventory_data_dictionary.csv",
                nrow(dd) + nrow(NOTES)))

write_csv_safe(loaded$status, file.path(paths$reports, "input_status.csv"))
write_integration_manifest(
  paths, inputs,
  list(tables = paths$tables, source_data = paths$source_data),
  list(dataset = run$dataset, fdr_alpha = FDR_ALPHA,
       recurrence_min_terms = RECUR_MIN_TERMS,
       recurrence_min_contexts = RECUR_MIN_CONTEXTS),
  paste0("Release of the pathway and leading-edge protein denominators behind ",
         "the selected results shown in Figures 2 and 3. Nothing is recomputed."))

message("Display selection inventories complete: ", paths$tables)
