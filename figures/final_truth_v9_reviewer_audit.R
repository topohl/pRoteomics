#!/usr/bin/env Rscript

# Part-26: claim crosswalk, selection dependency, provenance and the
# reviewer-vulnerability report.
#
# AUDIT ONLY. Nothing here changes a figure. Several of these tables record
# EDITORIAL facts about how the package was built - which examples were chosen
# and why, what a panel can and cannot support - because those facts are not
# derivable from the data and are exactly what a reviewer will ask about.

source(file.path("R", "paths.R"))
source(repo_path("R", "null_coalescing.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
source(repo_path("R", "story_v4_figure_panels.R"))
source(repo_path("R", "final_truth_v9_figure_utils.R"))
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_reviewer_audit.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] final_truth_v9 reviewer audit")
  quit(save = "no", status = 0L)
}

OUT <- path_results("tables", "manuscript_candidates", "final_truth_v9", "audit")
REP <- path_results("reports", "manuscript_candidates", "final_truth_v9")
dir_create(OUT); dir_create(REP)
N <- "9 animals; 3 per group"

# ============================================== S2 claim-evidence crosswalk
C <- function(id, claim, panel, etype, ds, lvl, mode, sel, indep, test, fam,
              eff, unc, strength, overread, wording)
  data.frame(claim_id = id, claim_text = claim, figure_panel = panel,
             evidence_type = etype, dataset = ds, spatial_level = lvl,
             biological_n = N, replicate_unit = "AnimalID",
             descriptive_or_inferential = mode, selection_dependency = sel,
             independent_validation = indep, test = test,
             multiple_testing_family = fam, effect_size_shown = eff,
             uncertainty_shown = unc, claim_strength_supported = strength,
             potential_overread = overread, recommended_wording = wording,
             status = "supported as worded", stringsAsFactors = FALSE)

cw <- rbind(
 C("F2.1", "compartments have distinct molecular identities", "F2e",
   "marker abundance", "all three", "compartment", "descriptive",
   "prespecified marker panel", "no", "none", "none", "yes", "no",
   "descriptive characterisation",
   "could be read as a test of compartment purity",
   "prespecified compartment markers are most abundant in their own compartment"),
 C("F2.2", "regional molecular architecture is reproducible", "F2f, ED1a",
   "left-right concordance", "all three", "region", "descriptive",
   "prespecified anatomical contrasts", "no", "none", "none", "yes", "no",
   "descriptive reproducibility",
   "reproducible could be read as validated against truth",
   "the regional anatomical effect is concordant between hemispheres"),
 C("F2.3", "fine laminar architecture is structured but less reproducible",
   "F2f, ED1a", "left-right concordance", "neuron_neuropil", "region x layer",
   "descriptive", "prespecified anatomical contrasts", "no", "none", "none",
   "yes", "no", "descriptive, and the weaker laminar result is left visible",
   "could be read as laminar structure being absent",
   "laminar concordance is lower than regional concordance"),
 C("F2.4", "phenotype-blind spatial fingerprints exist", "F2d",
   "baseline abundance", "all three", "region x layer + region", "descriptive",
   "top-RANKED proteins of CON-only contrasts; same contrasts are visualised",
   "no - selection and display share the CON structure", "none", "none",
   "yes", "no", "descriptive characterisation only",
   "could be read as independent validation of spatial structure",
   "phenotype-blind descriptive characterisation of baseline spatial architecture"),
 C("F2.5", "external signatures support anatomical validity", "F2g, ED2b",
   "GSEA against published signatures", "all three", "region + CA1 layers",
   "inferential", "expected pairings shown; full inventory in ED2",
   "YES - external reference", "gene set enrichment",
   "BH within the external-validation inventory", "yes (NES)", "no",
   "the one genuinely external support in the package",
   "the only panel a reviewer should accept as validation",
   "internal anatomical contrasts recover independently published signatures"),
 C("F2.6", "anatomical contrasts map onto coherent functional programs",
   "F2h, ED2c", "GSEA against canonical GO", "all three",
   "region + CA1 layers", "inferential", "one canonical term per contrast",
   "no - same proteomic data", "gene set enrichment",
   "BH within the internal anatomical inventory", "yes (NES)", "no",
   "functional characterisation of the same data",
   "could be mistaken for a second, independent validation",
   "functional characterisation, not independent validation"),
 C("F3.1", "individual-protein differential abundance is sparse", "F3a",
   "counts of FDR-supported proteins", "all three", "18 spatial units",
   "inferential (counts)", "all FDR-supported SUS-RES proteins", "no",
   "differential abundance", "BH within the DA family", "no (counts)", "no",
   "12 of 18 units carry no FDR-supported SUS-RES protein",
   "absence of DAPs could be read as absence of biology",
   "sparse individual-protein differences at this sample size"),
 C("F3.2", "CA2-SLM concentration attenuates after QC", "F3a, ED3",
   "robustness qualification", "neuron_neuropil", "CA2_slm",
   "descriptive QC", "the 28 CA2-SLM FDR-supported proteins", "no",
   "missingness and QC audit", "not a new test family", "yes", "no",
   "28 to 6 after qualification",
   "could be read as a second statistical test",
   "not excluded by the CA2-SLM missingness audit"),
 C("F3.3", "coordinated program differences are richer than DAP burden",
   "F3a vs F3b", "GSEA theme atlas", "all three", "18 spatial units",
   "descriptive aggregation over inferential inputs",
   "claim-eligible themes only", "no", "gene set enrichment",
   "BH within each constituent GO family; no theme family", "yes (median NES)",
   "no", "program-level structure where protein-level structure is sparse",
   "theme colour could be read as an independently tested quantity",
   "descriptive theme aggregation; support marks are constituent-term FDR"),
 C("F3.4", "CA3-SR synaptic program is susceptibility-associated", "F3c, F3d",
   "GSEA", "neuron_neuropil", "CA3_sr", "inferential",
   "editorially chosen exemplar, one per compartment", "no",
   "gene set enrichment", "BH within the GSEA family", "yes (NES)", "no",
   "susceptibility-associated; RES-CON is not FDR-supported",
   "could be read as an established RES-versus-SUS divergence",
   "susceptibility-associated, not divergent"),
 C("F3.5", "CA2 soma mRNA-processing program is susceptibility-associated",
   "F3c, F3e", "GSEA", "neuron_soma", "CA2_sp", "inferential",
   "editorially chosen exemplar", "no", "gene set enrichment",
   "BH within the GSEA family", "yes (NES)", "no",
   "susceptibility-associated", "same as F3.4",
   "susceptibility-associated, not divergent"),
 C("F3.6", "CA1 microglia-enriched OXPHOS shows a graded stress-associated pattern",
   "F3c, F3f", "GSEA", "microglia", "CA1", "inferential",
   "editorially chosen exemplar", "no", "gene set enrichment",
   "BH within the GSEA family", "yes (NES)", "no",
   "the only exemplar whose three-group ordering is fully FDR-supported",
   "microglia-enriched ROI could be read as cell-intrinsic microglia",
   "microglia-enriched local microenvironment, program level"),
 C("F3.7", "displayed proteins contribute to those enrichment results",
   "F3 g/h/i", "leading-edge decomposition", "all three",
   "3 spatial units", "descriptive",
   "leading edge of the SAME enrichment, by stored rank statistic",
   "no - decomposition, not corroboration", "differential abundance",
   "BH within the DA family", "yes (log2FC)",
   "no - no valid stored interval",
   "descriptive decomposition; none is individually FDR-supported",
   "could be read as individual-protein confirmation of the program",
   "selected leading-edge proteins contributing to the enrichment"),
 C("ED.1", "bilateral averaging improves reliability", "ED1b",
   "intraclass correlation", "all three", "endpoint classes", "descriptive",
   "all stored endpoints", "no", "none", "none", "yes (ICC)", "no",
   "precision/reliability of the animal-level bilateral estimate",
   "could be read as improved accuracy against a truth",
   "improves reliability of the animal-level estimate; hemispheres are repeated tissue"),
 C("ED.2", "WGCNA modules are spatially organised", "ED_WGCNA a",
   "module-member abundance", "all three", "18 spatial units", "descriptive",
   "all modules", "no", "none", "none", "yes (CON z)", "no",
   "descriptive spatial organisation",
   "could be read as a tested spatial effect",
   "modules show descriptive spatial organisation"),
 C("ED.3", "WGCNA phenotype effects are not FDR-supported", "ED_WGCNA b",
   "module eigengene difference", "neuron_neuropil", "spatially adjusted",
   "inferential (negative)", "all 15 modules", "no",
   "module eigengene group difference", "BH over 45 module x contrast cells",
   "yes", "no", "0 of 45 reach FDR support",
   "could be read as evidence that module effects are absent",
   "no module x contrast cell survived multiple-testing correction"),
 C("ED.4", "stress x space WGCNA interactions are not FDR-supported",
   "ED_WGCNA b", "interaction omnibus", "all three", "module x unit",
   "inferential (negative)", "all modules", "no", "interaction omnibus",
   "BH over 35 tests", "no", "no",
   "0 of 35, smallest FDR 0.27",
   "could be read as proving spatial homogeneity",
   "no interaction survived correction; absence is not established"),
 C("ED.5", "robustness-qualified effects commonly occur away from the baseline dominant location",
   "ED7", "spatial cross-tabulation", "all three", "18 spatial units",
   "descriptive within a selected set",
   "proteins already FDR-supported AND robustness-qualified",
   "no - selected set", "none", "none", "yes (counts)", "no",
   "15 of 15 outside the dominant unit; 14 of 15 outside the affinity set",
   "could be read as an inferential spatial-null result",
   "among the robustness-qualified proteins, effects were observed outside"),
 C("ED.6", "no whole-network group difference is detected", "ED8c",
   "exact permutation", "all three", "network", "inferential (negative)",
   "all edges", "no", "exact permutation over animal labels",
   "exact enumeration", "no", "yes (attainable floor)",
   "informative null: resolution to 0.0036",
   "could be read as proving networks are identical",
   "no whole-network difference detected at this sample size"),
 C("ED.7", "no edge-behaviour association survives correction", "ED8d",
   "Pearson correlation", "all three", "edge x outcome",
   "inferential (negative)", "all 48 distinct tests", "no",
   "Pearson correlation across animals", "BH over 48 tests", "yes (r)",
   "yes (95% CI)", "limited-power negative result",
   "could be read as evidence of no coupling",
   "no association survived correction; power is very limited at n = 9"))
write_csv_safe(cw, file.path(OUT, "final_claim_evidence_crosswalk.csv"))

# ==================================== S4 selection dependency
sd <- data.frame(
  figure_panel = c("F2d", "F2g", "F2h", "F3 d/e/f", "F3 g/h/i", "ED7",
                   "ED_WGCNA d (m11)"),
  displayed_result = c("baseline spatial fingerprint",
                       "external signature enrichment",
                       "canonical GO program enrichment",
                       "three representative GSEA programs",
                       "leading-edge protein effects",
                       "baseline versus effect location",
                       "m11 worked example"),
  how_selected = c(
    "top-RANKED proteins of prespecified CON-only anatomical contrasts",
    "prespecified expected pairings against external signatures",
    "one canonical GO term per anatomical contrast",
    "EDITORIAL: one exemplar per compartment from FDR-supported GSEA results",
    "leading edge of the same enrichment, ranked by stored rank statistic",
    "proteins already FDR-supported and robustness-qualified",
    "chosen as a worked example; external EWCE is context only"),
  selected_from_same_analysis = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
  constitutes_independent_validation = FALSE,
  consequence = c(
    "descriptive characterisation, not validation of spatial structure",
    "the one externally anchored panel in the package",
    "characterisation of the same data",
    "exemplars illustrate; the atlas carries the complete evidence",
    "decomposition of the enrichment, not corroboration of it",
    "selected-set descriptive pattern, not a spatial-null test",
    "label is not changed by the external evidence"),
  caption_states_this = TRUE, stringsAsFactors = FALSE)
sd$constitutes_independent_validation[sd$figure_panel == "F2g"] <- TRUE
write_csv_safe(sd, file.path(OUT, "selection_dependency_audit.csv"))

# ==================================== S10 F3 program example selection
pr <- s4_programs()
th <- nv_read_csv(repo_path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"))
# The relevant denominator is TWO-fold: how the exemplar ranks among every
# FDR-supported claim-eligible theme in its compartment, and how it ranks
# within its OWN spatial unit. The second is the one a reviewer will ask about.
rank_in_compartment <- function(i) {
  z <- th[th$dataset == pr$dataset[i] & th$contrast == "SUS - RES" &
            th$theme_claim_eligible %in% TRUE, , drop = FALSE]
  z <- z[is.finite(z$GSEA_FDR) & z$GSEA_FDR < 0.05, , drop = FALSE]
  if (!nrow(z)) return(c(NA_integer_, NA_integer_, NA_real_))
  z <- z[order(-abs(z$NES)), , drop = FALSE]
  k <- which(z$GO_ID == pr$term[i] & z$spatial_unit == pr$unit[i])[1]
  zu <- z[z$spatial_unit == pr$unit[i], , drop = FALSE]
  zu <- zu[order(-abs(zu$NES)), , drop = FALSE]
  ku <- which(zu$GO_ID == pr$term[i])[1]
  c(if (length(k) && !is.na(k)) k else NA_integer_, nrow(z),
    if (length(k) && !is.na(k)) z$NES[k] else NA_real_,
    if (!is.na(ku)) ku else NA_integer_, nrow(zu))
}
rk <- t(vapply(seq_len(nrow(pr)), rank_in_compartment, numeric(5)))
sel <- data.frame(
  program = pr$label, dataset = pr$dataset, spatial_unit = pr$unit,
  go_term = pr$term,
  candidate_universe = "claim-eligible themes, FDR-supported for SUS - RES, within the compartment",
  n_candidates = rk[, 2],
  rank_by_abs_NES_in_compartment = rk[, 1], NES = rk[, 3],
  rank_by_abs_NES_within_own_unit = rk[, 4],
  n_fdr_supported_terms_in_own_unit = rk[, 5],
  criterion = "EDITORIAL: one exemplar per compartment; not an algorithmic maximum",
  one_per_compartment = TRUE,
  strongest_means = "NOT the maximum by any measure; the table is a hardcoded literal in s4_programs()",
  selected_after_seeing_phenotype_direction = TRUE,
  honest_statement = paste0(
    "The three exemplars are an editorial choice of one FDR-supported program ",
    "per compartment, fixed as a literal table. They were chosen after the ",
    "phenotype-contrast results were known. They are illustrative, and the ",
    "complete evidence they illustrate is the theme atlas and ED6, which are ",
    "not selected. They must NOT be described as the strongest programs: the ",
    "synaptic exemplar ranks 51st of 60 FDR-supported terms within CA3-SR ",
    "itself, where RNA-processing terms carry larger NES. It is a valid ",
    "FDR-supported result chosen to represent distinct biology per ",
    "compartment, not a maximum."),
  stringsAsFactors = FALSE)
write_csv_safe(sel, file.path(OUT, "f3_program_example_selection.csv"))

# ==================================== S11 GSEA provenance
gp <- data.frame(
  analysis = c(rep("direct GSEA example", 3), "theme atlas"),
  program = c(pr$label, "all claim-eligible themes"),
  dataset = c(pr$dataset, "all three"),
  spatial_unit = c(pr$unit, "18 spatial units"),
  biological_n = N,
  ranking_statistic = "stored per-gene contrast statistic, collapsed from protein groups by the prespecified median rule",
  universe = "genes measured in that spatial unit and contrast",
  gene_set_database = "Gene Ontology, biological process",
  enrichment_method = "clusterProfiler GSEA (stored canonical result)",
  testing_method = "permutation-based, as implemented upstream",
  fdr_family = "gene sets within one spatial unit and contrast",
  software = "clusterProfiler; membership reconstructed via org.Mm.eg.db",
  annotation_pinned = TRUE,
  exact_reconstruction = c(rep("membership reproduces the canonical setSize exactly", 3),
                           "not reconstructed; read from the stored theme table"),
  fdr_interpretation = paste0(
    "gene-set enrichment evidence conditional on the ranked protein-level ",
    "contrast; NOT a count of independent biological observations"),
  stringsAsFactors = FALSE)
write_csv_safe(gp, file.path(OUT, "gsea_statistical_provenance_audit.csv"))

# ==================================== S21/S22 main-vs-ED role and prominence
ctr <- s9f_contract()
rows <- list()
CLASS <- c(
  v9_schematic = "technical_QC", v9_depth = "technical_QC",
  v9_pca = "central_descriptive", v9_fingerprint = "central_descriptive",
  v9_compartment = "supporting_validation", v9_bilateral_main = "technical_QC",
  v9_external_main = "supporting_validation",
  v9_internal_main = "central_descriptive",
  v9_dap_track = "central_inferential", v9_atlas = "central_inferential",
  v9_bridge = "central_descriptive", v9_curve_syn = "central_inferential",
  v9_curve_rna = "central_inferential", v9_curve_ox = "central_inferential",
  v9_prot_syn = "central_descriptive", v9_prot_rna = "central_descriptive",
  v9_prot_ox = "central_descriptive")
for (f in ctr$figures) for (it in f$layout) {
  id <- as.character(it$panel)
  p <- Filter(function(q) identical(as.character(q$id), id), ctr$panels)[[1]]
  main <- grepl("figure_0", as.character(f$figure_key))
  rows[[length(rows) + 1L]] <- data.frame(
    figure = as.character(f$name), panel = as.character(it$label), panel_id = id,
    main_or_ed = if (main) "main" else "extended_data",
    declared_role = as.character(p$main_or_ed %||% ""),
    area_mm2 = as.numeric(it$w) * as.numeric(it$h),
    evidence_class = unname(CLASS[id] %||% NA) %||%
      if (main) "central_descriptive" else "negative_constraint",
    stringsAsFactors = FALSE)
}
rl <- do.call(rbind, rows)
rl$evidence_class[is.na(rl$evidence_class)] <-
  ifelse(rl$main_or_ed[is.na(rl$evidence_class)] == "main",
         "central_descriptive", "negative_constraint")
rl$share_of_page <- round(100 * rl$area_mm2 /
  ave(rl$area_mm2, rl$figure, FUN = sum), 1)
rl$prominence_appropriate <- TRUE
write_csv_safe(rl, file.path(OUT, "final_main_ed_role_audit.csv"))

cat("\n===== PART-26 REVIEWER AUDIT =====\n")
cat("claims audited              :", nrow(cw), "\n")
cat("  selection-dependent       :", sum(sd$selected_from_same_analysis), "of",
    nrow(sd), "\n")
cat("  independently validated   :", sum(sd$constitutes_independent_validation), "\n")
cat("F3 exemplar rank by |NES| within compartment:",
    paste(sprintf("%s=%s/%s", pr$key, rk[, 1], rk[, 2]), collapse = "  "), "\n")
cat("panels classified           :", nrow(rl), "\n")
print(table(rl$main_or_ed, rl$evidence_class))
cat("\nwritten to:", relative_to(OUT), "\n")

# ================================================= S23 reviewer vulnerabilities
V <- function(n, obj, answered, where, residual, wording_fixes, needs_experiment)
  data.frame(id = n, objection = obj, figures_answer_it = answered,
             where = where, residual_weakness = residual,
             solvable_by_wording = wording_fixes,
             new_experiment_required = needs_experiment,
             stringsAsFactors = FALSE)
vl <- rbind(
 V(1, "n = 3 per group is too small for proteomics", "partially",
   "F2b states 9 animals; every legend states biological n",
   "3 per group genuinely limits power; several nulls are power-limited rather than informative",
   "yes - state power honestly and avoid equivalence language", "no"),
 V(2, "spatial samples are pseudoreplicated within animal", "yes",
   "F2b separates 323 acquisitions from 9 animals; DA is a 3v3 animal-level design",
   "the acquisition-level PCA and depth panels are descriptive only, which is stated",
   "already stated", "no"),
 V(3, "hemispheres are treated as independent replicates", "yes",
   "ED1 and ED8b; bilateral averaging is framed as precision, not independence",
   "none identified", "already stated", "no"),
 V(4, "CA2-SLM results are a missingness artefact", "yes",
   "F3a, ED3 a-e; 28 canonical fall to 6 robustness-qualified",
   "the 6 survivors rest on a small, QC-screened set",
   "yes - never call CA2-SLM a hotspot", "no"),
 V(5, "program examples were chosen after seeing the result", "partially",
   "f3_program_example_selection.csv records the rule",
   "the three exemplars ARE editorial and were chosen after the contrasts were known; the synaptic exemplar ranks 51st of 60 FDR-supported terms within CA3-SR",
   "yes - call them illustrative exemplars, never the strongest programs", "no"),
 V(6, "GSEA FDR values look implausibly small for n = 3", "yes",
   "every GSEA legend states that the FDR is conditional on the ranked protein statistic",
   "small FDRs still invite over-reading",
   "yes - keep the conditional-on-ranking wording", "no"),
 V(7, "leading-edge proteins are not individually significant", "yes",
   "F3 g/h/i legends; smallest BH FDR among the 63 displayed values is 0.53",
   "none - this is stated plainly", "already stated", "no"),
 V(8, "WGCNA phenotype nulls may simply be underpowered", "yes",
   "ED_WGCNA b states 0/45 and 0/35 with the smallest FDR",
   "a negative interaction result cannot establish spatial homogeneity",
   "yes - no detectable X at this sample size", "no"),
 V(9, "microglia-enriched ROI is treated as cell-intrinsic microglia", "yes",
   "the artwork headers and the F2a, F3f, F3i and ED6e legends state that the microglia-enriched ROI is an enriched context, not sorted cells",
   "the compartment is an enriched measurement context, not sorted cells",
   "done - the caveat now travels with every microglia panel", "would need sorted or single-cell microglia"),
 V(10, "internal characterisation is presented as validation", "yes",
   "F2h is functional characterisation; only F2g is external",
   "the package rests on ONE external anchor",
   "yes", "an independent external dataset would strengthen this"),
 V(11, "there is no independent proteomics replication", "no",
   "not addressed by any panel",
   "this is the single largest structural limitation of the package",
   "no - wording cannot fix it", "YES - an independent cohort"),
 V(12, "ED7 15/15 is presented as a spatial principle", "yes",
   "ED7 legends say among the robustness-qualified proteins and descriptive",
   "the observation is within a set already selected for being phenotype-associated",
   "yes - keep selected-set wording",
   "a prespecified spatial-null test would be needed"),
 V(13, "network nulls with n = 9 prove nothing", "partially",
   "ED8c shows the attainable floor (0.0036) beside the observed p",
   "the floor makes the whole-network null informative; the coupling null remains power-limited",
   "yes - limited-power negative result", "no"))
write_csv_safe(vl, file.path(OUT, "nature_reviewer_vulnerabilities.csv"))

md <- c("# Nature-reviewer vulnerability report", "",
  "Candidate layer; audit only. Ordered roughly by how hard each objection is",
  "to answer with the present data.", "")
for (i in seq_len(nrow(vl))) md <- c(md,
  sprintf("## %d. %s", vl$id[i], vl$objection[i]), "",
  sprintf("- **Do the figures answer it?** %s", vl$figures_answer_it[i]),
  sprintf("- **Where:** %s", vl$where[i]),
  sprintf("- **Residual weakness:** %s", vl$residual_weakness[i]),
  sprintf("- **Fixable by wording:** %s", vl$solvable_by_wording[i]),
  sprintf("- **New experiment required:** %s", vl$new_experiment_required[i]), "")

# ================================================= S24 highest-value validation
md <- c(md, "# Highest-value missing validation", "",
  paste0("Audit only; this changes no figure. Ranked by how much each would ",
         "move a reviewer on the CENTRAL claim, which is that stress outcome ",
         "is associated with spatially restricted coordinated program ",
         "differences rather than with broad protein-level change."), "",
  "## 1. An independent animal cohort, same spatial workflow", "",
  paste0("**Impact: highest. Feasibility: lowest.** Vulnerability 11 is the ",
         "only objection in this audit that no wording can address and that no ",
         "existing panel touches. Every current result rests on one cohort of ",
         "9 animals. A second cohort, even a small one, processed through the ",
         "same spatial workflow and asked only whether the three exemplar ",
         "programs move in the same direction, would convert the central claim ",
         "from internally consistent to replicated. Nothing else in this list ",
         "does that."), "",
  "## 2. Spatially matched orthogonal validation of one exemplar program", "",
  paste0("**Impact: high. Feasibility: moderate.** Immunostaining or targeted ",
         "proteomics for a small number of leading-edge proteins in the ",
         "matching spatial unit would test whether a program-level enrichment ",
         "corresponds to a measurable local protein difference. The CA1 ",
         "microglia-enriched OXPHOS program is the best candidate: it is the ",
         "only exemplar whose three-group ordering is fully FDR-supported, and ",
         "it is where the compartment interpretation (enriched ROI versus ",
         "cell-intrinsic) is most exposed, so a sorted or immunolabelled ",
         "measurement would answer vulnerabilities 9 and 10 at once. It would ",
         "NOT address replication."), "",
  "## 3. A second external spatial reference", "",
  paste0("**Impact: moderate. Feasibility: highest.** The package rests on a ",
         "single external anchor (F2g). Adding a second independent spatial ",
         "proteomic or transcriptomic reference would strengthen the ",
         "anatomical validity argument at essentially no experimental cost, ",
         "since it is a reanalysis rather than new data collection. It ",
         "strengthens Figure 2 rather than the central Figure 3 claim, which ",
         "is why it ranks third despite being the easiest."), "",
  paste0("_Not recommended:_ deriving any additional statistic from the ",
         "present data. This audit found no unexploited valid uncertainty, and ",
         "manufacturing one would add apparent rigour without adding ",
         "information."), "",
  "_Generated by `figures/final_truth_v9_reviewer_audit.R`. Do not edit._")
writeLines(md, file.path(REP, "nature_reviewer_vulnerabilities.md"))
cat("vulnerabilities:", nrow(vl), "| needing new data:",
    sum(!vl$new_experiment_required %in% "no"), "\n")
