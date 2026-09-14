#!/usr/bin/env Rscript

# Part-27: terminology registry, claim-verb registry and the semantic rules.
#
# SEMANTIC AUDIT ONLY. No analysis is rerun and no figure is redesigned here.
# These tables fix, once, what each thing is CALLED at each level of evidence -
# exact canonical GO term, reviewed umbrella theme, figure label, manuscript
# phrase - so that a familiar shorthand cannot quietly claim more than the
# analysis establishes.
#
# The exact terms and umbrella memberships below are READ from the canonical
# ontology mapping, not asserted, so the registry cannot drift from it.

source(file.path("R", "paths.R"))
source(repo_path("R", "null_coalescing.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
source(repo_path("R", "story_v4_figure_panels.R"))
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_semantics.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] final_truth_v9 semantic registry")
  quit(save = "no", status = 0L)
}

OUT <- path_results("tables", "manuscript_candidates", "final_truth_v9", "audit")
REP <- path_results("reports", "manuscript_candidates", "final_truth_v9")
dir_create(OUT); dir_create(REP)

TH <- nv_read_csv(repo_path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"))
elig <- TH[TH$theme_claim_eligible %in% TRUE, , drop = FALSE]
go_name <- function(id) unique(TH$GO_description[TH$GO_ID == id])[1]
theme_name <- function(tid) unique(elig$manuscript_theme[elig$theme_id == tid])[1]
theme_n <- function(tid) length(unique(elig$GO_ID[elig$theme_id == tid]))
theme_members <- function(tid)
  sort(unique(elig$GO_description[elig$theme_id == tid]))

# ==================================================== S1/S30 terminology registry
T <- function(canon, src, figlab, manu, short, spat, cell, infer, may, maynot,
              notes)
  data.frame(canonical_analysis_term = canon, canonical_source = src,
             exact_figure_label = figlab, preferred_manuscript_label = manu,
             allowed_short_label = short, spatial_resolution = spat,
             cell_type_or_compartment_status = cell, inferential_status = infer,
             permitted_interpretation = may, prohibited_interpretation = maynot,
             notes = notes, stringsAsFactors = FALSE)

GO_EXACT <- "exact canonical GO biological-process term"
UMB <- "reviewed umbrella theme over multiple canonical GO terms"

reg <- rbind(
 T(go_name("GO:0006119"), paste0(GO_EXACT, " GO:0006119"),
   "oxidative phosphorylation", "oxidative-phosphorylation program",
   "OXPHOS", "microglia-enriched ROI, CA1 (region level)",
   "enriched ROI, not sorted cells", "FDR-supported (NES -2.11)",
   "reduced oxidative-phosphorylation program in the CA1 microglia-enriched ROI",
   "mitochondrial dysfunction; impaired mitochondrial function; defective OXPHOS; microglial metabolic dysfunction",
   "OXPHOS is acceptable shorthand ONLY for this exact term, never for the umbrella"),
 T(theme_name("mitochondrial_respiration_oxphos"), paste0(UMB, "; ",
     theme_n("mitochondrial_respiration_oxphos"), " GO terms"),
   "Mitochondrial respiration", "mitochondrial respiration / oxidative phosphorylation",
   "mitochondrial respiration", "18 spatial units",
   "compartment-level", "descriptive median NES; no theme-level FDR",
   "a coordinated mitochondrial respiration / OXPHOS program summary",
   "OXPHOS as the name of the whole theme; mitochondrial dysfunction",
   paste0("the umbrella also contains 4 GLYCOLYTIC terms, which are cytosolic; ",
          "within the CA1 microglia cell those 4 are the non-supported ones, so ",
          "the FDR-supported signal in this theme is mitochondrial respiration")),
 T(go_name("GO:0006397"), paste0(GO_EXACT, " GO:0006397"),
   "mRNA processing", "mRNA-processing program", "mRNA processing",
   "neuronal soma, CA2 (region level)", "soma-enriched capture",
   "FDR-supported (NES +2.13)",
   "increased mRNA-processing program in the CA2 neuronal-soma compartment",
   "RNA biology; transcriptomic change; splicing dysfunction",
   "do not silently broaden the exact term to RNA processing"),
 T(theme_name("rna_processing_splicing_rnp"), paste0(UMB, "; ",
     theme_n("rna_processing_splicing_rnp"), " GO terms"),
   "RNA processing", "RNA processing / splicing / RNP organization",
   "RNA processing", "18 spatial units", "compartment-level",
   "descriptive median NES; no theme-level FDR",
   "a coordinated RNA-processing program summary",
   "mRNA processing as the name of the whole theme",
   paste0("the umbrella genuinely includes rRNA maturation, miRNA processing ",
          "and mitochondrial RNA processing, so it is broader than mRNA processing")),
 T(go_name("GO:0099536"), paste0(GO_EXACT, " GO:0099536"),
   "synaptic signalling", "synaptic-signalling program", "synaptic signalling",
   "neuropil, CA3-SR (region x layer)", "neuropil compartment",
   "FDR-supported (NES -1.72)",
   "reduced synaptic-signalling program in CA3-SR neuropil",
   "synaptic dysfunction; impaired synapses; reduced synaptic transmission",
   paste0("GO spells this 'synaptic signaling'; the figure uses UK spelling. ",
          "A protein-abundance program is not a physiological measurement")),
 T(theme_name("synaptic_signaling_vesicle"), paste0(UMB, "; ",
     theme_n("synaptic_signaling_vesicle"), " GO terms"),
   "Synaptic signalling", "synaptic signalling / vesicle-mediated transport",
   "synaptic / vesicle", "18 spatial units", "compartment-level",
   "descriptive median NES; no theme-level FDR",
   "a coordinated synaptic and vesicle-transport program summary",
   "equating the theme with the exact GO:0099536 term",
   paste0("the atlas label and the exact term differ only in spelling, so the ",
          "legend must state that the atlas row is the broader umbrella")),
 T("WGCNA module", "co-expression module from the WGCNA stage",
   "m01 ... m15", "co-abundance module", "module", "compartment-specific",
   "defined on protein co-abundance, not on cell identity",
   "descriptive; 0/45 phenotype cells FDR-supported",
   "a group of co-abundant proteins with a reviewed label",
   "a cell type; a pathway; a causal unit",
   "module identity, cell-type affinity and spatial peak are three separate things"),
 T("module eigengene", "first principal component of a module",
   "Module eigengene difference", "module eigengene difference", "eigengene",
   "spatially adjusted global model", "compartment-level",
   "inferential but FDR-negative (0/45)", "a per-animal module summary score",
   "module abundance; module activity",
   "distinct from module-member abundance, which is the ED_WGCNA panel a quantity"),
 T("module-member abundance", "mean CON z of a module's member proteins",
   "Mean module-member abundance (CON z-score)",
   "mean abundance of module member proteins", "member abundance",
   "18 spatial units", "compartment-level", "descriptive",
   "where a module's proteins are most abundant in control animals",
   "the module eigengene; module activity",
   "must never share a name with the eigengene quantity"),
 T("external cell-type affinity", "EWCE against external reference data",
   "external cell type", "external cell-type affinity", "EWCE affinity",
   "module level", "EXTERNAL annotation, not a measured cell state",
   "descriptive context", "the closest external reference cell type",
   "cell-intrinsic identity; automatic relabelling of the module",
   "context only; it never becomes the module label"),
 T("microglia-enriched ROI", "sampling design",
   "Microglia ROI", "microglia-enriched ROI", "microglia-enriched ROI",
   "region level only", "ENRICHED measurement compartment, not sorted cells",
   "not applicable", "a local microenvironment measurement",
   "microglial proteome; microglia-specific; cell-intrinsic; cell-autonomous",
   "this is the single most attackable compartment claim in the package"),
 T("neuronal soma", "sampling design", "Soma", "neuronal soma compartment",
   "soma", "region level only", "soma-enriched capture", "not applicable",
   "a region-level soma-enriched compartment",
   "a laminar soma experiment equivalent to a neuropil layer contrast",
   "CA2_sp names the anatomical source, not an independently contrasted layer"),
 T("neuropil", "sampling design", "Neuropil", "neuropil compartment",
   "neuropil", "region x layer", "neuropil compartment", "not applicable",
   "the only compartment with genuine laminar resolution", "",
   "10 units across CA1/CA2/CA3/DG; there is no CA3-SLM"),
 T("dominant control spatial unit", "CON-only spatial atlas",
   "most abundant HERE at baseline", "dominant control spatial unit",
   "dominant control unit", "18 spatial units", "compartment-level",
   "descriptive", "where a protein is most abundant in control animals",
   "a pre-stress longitudinal baseline",
   "nothing was measured pre-stress; prefer control over baseline"),
 T("control affinity set", "canonical effect_identity_relationship",
   "canonical baseline affinity set", "control-reference affinity set",
   "control affinity set", "18 spatial units", "compartment-level",
   "descriptive classification",
   "the top-2 control-reference units for a protein",
   "a tested spatial null", "authoritative column for the 14/15 statement"),
 T("robustness-qualified protein", "CA2-SLM missingness and QC audit",
   "Robustness-qualified (15 total)", "robustness-qualified protein",
   "robustness-qualified", "18 spatial units", "compartment-level",
   "not a uniform second-stage test",
   "not excluded by the CA2-SLM missingness audit",
   "passed an additional significance test",
   "28 CA2-SLM proteins were audited and 6 qualified; the 9 outside were never at risk"),
 T("FDR-supported DAP", "differential abundance, SUS vs RES",
   "FDR-supported DAPs (37 total)", "FDR-supported differentially abundant protein",
   "FDR-supported DAP", "18 spatial units", "compartment-level",
   "inferential, BH within the DA family",
   "a protein reaching BH FDR < 0.05 for SUS versus RES in that unit",
   "a stress effect; a marker; a driver",
   "SUS-RES compares outcomes AMONG stress-exposed animals"),
 T("GSEA theme", "reviewed ontology-aware theme registry",
   "atlas row labels", "molecular program theme", "theme", "18 spatial units",
   "compartment-level", "DESCRIPTIVE aggregation; no theme-level FDR family",
   "a median NES over mapped canonical GO terms",
   "FDR-significant theme; theme p-value; theme-level significance",
   "a support dot means >= 1 constituent GO term passed its own FDR"),
 T("canonical GO term", "Gene Ontology biological process",
   "panel titles in F3 d/e/f and ED6 c/d/e", "canonical GO term", "GO term",
   "one spatial unit", "compartment-level",
   "inferential, BH within its own family",
   "the exact tested gene set", "a pathway, unless genuinely pathway-defined",
   "GO BP terms are processes, not pathways"),
 T("leading-edge protein", "stored GSEA leading edge",
   "gene labels in F3 g/h/i", "leading-edge protein", "leading edge",
   "one spatial unit", "compartment-level",
   "DESCRIPTIVE decomposition of the same enrichment",
   "a protein contributing to the enrichment signal",
   "validated protein; significant protein; driver protein; key protein",
   "none of the 63 displayed values is individually FDR-supported (min BH 0.53)"),
 T("susceptibility-associated", "SUS-RES and SUS-CON contrasts",
   "susceptibility-associated", "susceptibility-associated",
   "susceptibility-associated", "varies", "compartment-level",
   "inferential where FDR-supported",
   "a difference associated with the susceptible outcome",
   "susceptibility-specific; divergent; a resilience program",
   "RES-CON is not FDR-supported for the synaptic and RNA exemplars"),
 T("graded stress-associated pattern", "three pairwise contrasts",
   "graded stress-associated", "graded, shared stress-associated pattern",
   "graded pattern", "microglia-enriched ROI CA1", "compartment-level",
   "inferential; the only exemplar with all three contrasts FDR-supported",
   "a directional pattern consistent across the three pairwise contrasts",
   "a trajectory; three independent replications",
   "the three pairwise contrasts are algebraically related"),
 T("stress x spatial-unit interaction", "WGCNA stage-05 omnibus",
   "not drawn", "stress-by-spatial-unit interaction", "interaction",
   "module x unit", "compartment-level", "inferential, FDR-negative (0/35)",
   "no interaction survived multiple-testing correction",
   "spatial homogeneity is established; module effects are absent",
   "smallest FDR 0.27"),
 T("spatial molecular network", "per-animal edge correlation network",
   "Molecular similarity between spatial proteomic profiles",
   "spatial molecular similarity network", "molecular network",
   "unit-by-unit", "compartment-level", "descriptive network construction",
   "correlation between two units' protein profiles",
   "anatomical connectivity; a projection; a circuit",
   "the panel title already carries the NOT-connectivity disclaimer"),
 T("network distance from the control consensus",
   "euclidean distance in Fisher-z edge space",
   "network distance from CON consensus", "distance from the control consensus",
   "network distance", "whole network", "compartment-level", "descriptive",
   "how far one animal's network sits from the control consensus",
   "network rewiring; reorganization",
   "CON animals use a leave-one-CON-out consensus"),
 T("whole-network group difference", "exact permutation over animal labels",
   "exact whole-network p", "whole-network group difference",
   "network difference", "whole network", "compartment-level",
   "inferential, negative; attainable floor 0.0036",
   "no whole-network difference was detected at this sample size",
   "networks are unchanged; no reorganization occurred",
   "the attainable floor makes this an informative rather than empty null"))
write_csv_safe(reg, file.path(OUT, "final_terminology_registry.csv"))

# =================================================== S25 claim-verb registry
V <- function(verb, ev, ok, notok)
  data.frame(verb = verb, allowed_evidence_type = ev, examples_allowed = ok,
             examples_not_allowed = notok, stringsAsFactors = FALSE)
vb <- rbind(
 V("shows", "any displayed quantity",
   "Figure 3b shows the theme-level enrichment atlas",
   "Figure 3b shows that stress reprograms the hippocampus"),
 V("indicates", "FDR-supported inferential result",
   "the FDR-supported enrichment indicates a coordinated program difference",
   "indicates a mechanism"),
 V("supports", "external or independent evidence",
   "the external hippocampal signatures support the anatomical assignment",
   "panel h supports the spatial claim (same data)"),
 V("is associated with", "any phenotype contrast, FDR-supported or descriptive",
   "the CA3-SR synaptic program is associated with susceptibility",
   "-"),
 V("predicts", "out-of-sample prediction only",
   "-", "the module predicts susceptibility (no out-of-sample test exists)"),
 V("distinguishes", "a classifier or a tested separation",
   "-", "the program distinguishes SUS from RES (no classifier was fit)"),
 V("enriches / is enriched for", "gene-set enrichment",
   "the ranked contrast is enriched for oxidative-phosphorylation proteins",
   "the compartment is enriched for microglia (that is the sampling design)"),
 V("corresponds to", "descriptive correspondence",
   "the module corresponds to an external oligodendrocyte affinity", "-"),
 V("validates", "EXTERNAL independent evidence only",
   "the external signature comparison validates the anatomical assignment",
   "the leading-edge proteins validate the enrichment"),
 V("drives", "causal or perturbation evidence", "-",
   "the program drives susceptibility"),
 V("mediates", "mediation analysis", "-", "microglia mediate the stress effect"),
 V("causes", "intervention", "-", "stress causes reduced OXPHOS"),
 V("reorganizes", "the tested network metric, in its null form", "-",
   "stress reorganizes the molecular network"),
 V("rewires", "not supported by any analysis here", "-",
   "stress rewires hippocampal networks"),
 V("redistributes / relocates", "not supported; ED7 is a selected-set description",
   "-", "proteins redistribute to other layers"),
 V("impairs", "functional or physiological measurement", "-",
   "impairs synaptic transmission"),
 V("restores", "intervention with a rescue arm", "-", "resilience restores OXPHOS"))
write_csv_safe(vb, file.path(OUT, "claim_verb_registry.csv"))

# =================================================== S12 spatial resolution
u <- sg_units()
sr <- data.frame(
  dataset = u$dataset, display_label = u$display,
  region = sub("_.*$", "", u$analysis_key),
  layer = ifelse(u$dataset == "neuron_neuropil",
                 toupper(sub("^[^_]*_", "", u$analysis_key)), ""),
  independent_layer_resolution = u$dataset == "neuron_neuropil",
  stringsAsFactors = FALSE)
sr$recommended_manuscript_phrase <- ifelse(
  sr$dataset == "neuron_neuropil",
  sprintf("%s %s neuropil", sr$region, sr$layer),
  ifelse(sr$dataset == "neuron_soma",
         sprintf("%s neuronal soma", sr$region),
         sprintf("%s microglia-enriched ROI", sr$region)))
sr$potential_overread <- ifelse(
  sr$independent_layer_resolution, "none",
  paste0("the label names the anatomical source of a REGION-LEVEL compartment; ",
         "it must not be read as an independently contrasted layer"))
write_csv_safe(sr, file.path(OUT, "spatial_resolution_semantics.csv"))

# =================================================== S13 validation language
vl <- data.frame(
  panel = c("F2d", "F2g", "F2h", "F3 g/h/i", "ED2b", "ED2c",
            "ED_WGCNA c (EWCE)"),
  what_it_is = c("phenotype-blind descriptive spatial characterisation",
                 "external validation against published signatures",
                 "internal functional characterisation of the same data",
                 "decomposition of the same GSEA signal",
                 "complete external validation inventory",
                 "complete functional characterisation inventory",
                 "external cell-type context"),
  may_be_called_validation = c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
  required_wording = c("descriptive characterisation",
                       "external validation / external support",
                       "functional characterisation / annotation",
                       "leading-edge decomposition",
                       "external validation inventory",
                       "functional characterisation inventory",
                       "external cell-type context or support"),
  stringsAsFactors = FALSE)
write_csv_safe(vl, file.path(OUT, "validation_language_audit.csv"))

# =================================================== S21 WGCNA module naming
appr <- nv_read_csv(repo_path(
  "results", "reviewer_audit", "wgcna_label_approval",
  "WGCNA_final_label_approval_table.csv"))
appr <- appr[appr$level == "module", , drop = FALSE]
mn <- data.frame(
  dataset = appr$dataset, module_id = sub("^WGCNA_", "", appr$entity_id),
  active_reviewed_label = trimws(sub("^WGCNA_m[0-9]+\\s*·\\s*", "",
                                     appr$current_active_label)),
  proposed_label = appr$proposed_final_label,
  proposal_status = appr$adjudication_action,
  confidence = appr$confidence, stringsAsFactors = FALSE)
mn$allowed_manuscript_label <- ifelse(
  mn$proposal_status %in% c("MIXED", "UNRESOLVED"),
  "mixed / unresolved - do not label", mn$active_reviewed_label)
mn$prohibited_label <- ifelse(
  mn$dataset == "neuron_neuropil" & mn$module_id == "m11",
  "oligodendrocyte / myelin - PROPOSED ONLY, not activated",
  "any proposed_final_label that has not been activated")
write_csv_safe(mn, file.path(OUT, "wgcna_module_naming_audit.csv"))

cat("\n===== PART-27 SEMANTIC REGISTRIES =====\n")
cat("terminology entries  :", nrow(reg), "\n")
cat("claim verbs          :", nrow(vb), " of which never allowed:",
    sum(vb$examples_allowed == "-"), "\n")
cat("spatial labels       :", nrow(sr), " with genuine layer resolution:",
    sum(sr$independent_layer_resolution), "\n")
cat("panels that may be called validation:", sum(vl$may_be_called_validation),
    "of", nrow(vl), "\n")
cat("modules              :", nrow(mn), " un-labelable (mixed/unresolved):",
    sum(grepl("^mixed", mn$allowed_manuscript_label)), "\n")
cat("\nexact vs umbrella, the three exemplars:\n")
for (g in c("GO:0099536", "GO:0006397", "GO:0006119"))
  cat(sprintf("  %-12s exact = %-28s\n", g, go_name(g)))
for (t in c("synaptic_signaling_vesicle", "rna_processing_splicing_rnp",
            "mitochondrial_respiration_oxphos"))
  cat(sprintf("  umbrella %-34s %2d GO terms\n", theme_name(t), theme_n(t)))
cat("\nwritten to:", relative_to(OUT), "\n")

# ============================================ S26 claim-chain, S27 story test
CC <- function(sent, panel, obj, mode, sel, verb_ok, fix)
  data.frame(sentence = sent, figure_panel = panel, statistical_object = obj,
             biological_n = "9 animals; 3 per group",
             descriptive_or_inferential = mode,
             selected_from_same_analysis = sel,
             verb_strength_appropriate = verb_ok, corrected_sentence = fix,
             stringsAsFactors = FALSE)
cc <- rbind(
 CC("Compartments carry distinct molecular identities.", "F2e",
    "median CON z of prespecified markers", "descriptive", FALSE, TRUE,
    "NONE - accurate as written"),
 CC("The regional molecular architecture is reproducible.", "F2f, ED1a",
    "left-right Pearson/Spearman/sign agreement", "descriptive", FALSE, FALSE,
    paste0("The regional molecular architecture is reproducible between ",
           "hemispheres. (Hemispheres are repeated tissue within an animal, ",
           "so this is internal consistency, not independent replication.)")),
 CC("Laminar architecture is structured but less reproducible than regional.",
    "F2f, ED1a", "left-right concordance by level", "descriptive", FALSE, TRUE,
    "NONE - the weaker laminar result is deliberately left visible"),
 CC("A phenotype-blind spatial fingerprint exists.", "F2d",
    "CON-only z across 18 units", "descriptive", TRUE, FALSE,
    paste0("Control animals alone define a spatial molecular fingerprint. ",
           "(Descriptive characterisation: the proteins shown were ranked from ",
           "the same CON contrasts being displayed, so this is not independent ",
           "validation of spatial structure.)")),
 CC("External signatures validate the anatomical assignment.", "F2g, ED2b",
    "GSEA NES against published signatures", "inferential", FALSE, TRUE,
    "NONE - this is the one externally anchored panel"),
 CC("Anatomical contrasts are validated against functional programs.",
    "F2h, ED2c", "GSEA NES against canonical GO", "inferential", TRUE, FALSE,
    paste0("Anatomical contrasts correspond to coherent functional programs. ",
           "(Functional characterisation of the same proteomic data, not ",
           "independent validation.)")),
 CC("Individual-protein differential abundance is sparse.", "F3a",
    "counts of BH FDR < 0.05 SUS-RES proteins", "inferential", FALSE, TRUE,
    paste0("Individual-protein differences are sparse: 37 proteins reach FDR ",
           "support across 18 spatial units, and 12 of 18 units contain none.")),
 CC("CA2-SLM is a hotspot of protein-level change.", "F3a, ED3",
    "28 canonical falling to 6 after the missingness audit", "descriptive",
    FALSE, FALSE,
    paste0("The apparent concentration of FDR-supported proteins in CA2-SLM ",
           "weakens after the missingness and QC audit, from 28 to 6.")),
 CC("Program-level differences are richer than the protein-level burden.",
    "F3a versus F3b", "median NES per theme x unit", "descriptive", FALSE,
    FALSE,
    paste0("Coordinated program-level differences are detectable in spatial ",
           "units where no individual protein reaches FDR support.")),
 CC("The CA3-SR synaptic program diverges between resilient and susceptible animals.",
    "F3c, F3d", "GSEA NES, GO:0099536", "inferential", TRUE, FALSE,
    paste0("A synaptic-signalling program in CA3-SR neuropil is ",
           "susceptibility-associated. (RES-CON is not FDR-supported, so the ",
           "RES and SUS arms must not be described as divergent.)")),
 CC("CA2 soma shows increased RNA processing.", "F3c, F3e",
    "GSEA NES, GO:0006397", "inferential", TRUE, FALSE,
    paste0("An mRNA-processing program in the CA2 neuronal-soma compartment is ",
           "susceptibility-associated. (The exact term is mRNA processing; ",
           "RNA processing is the broader atlas umbrella.)")),
 CC("Microglia show reduced OXPHOS.", "F3c, F3f", "GSEA NES, GO:0006119",
    "inferential", TRUE, FALSE,
    paste0("An oxidative-phosphorylation program in the CA1 ",
           "microglia-enriched ROI is reduced, with a graded pattern across ",
           "the three pairwise contrasts. (An enriched ROI is not a ",
           "cell-intrinsic microglial measurement.)")),
 CC("Leading-edge proteins validate the enrichment results.", "F3 g/h/i",
    "log2FC of selected leading-edge proteins", "descriptive", TRUE, FALSE,
    paste0("Selected leading-edge proteins contributing to each enrichment are ",
           "shown. None is individually FDR-supported (smallest BH FDR 0.53), ",
           "so they decompose the signal rather than confirm it.")),
 CC("Bilateral averaging improves accuracy.", "ED1b", "ICC single vs bilateral",
    "descriptive", FALSE, FALSE,
    paste0("Bilateral averaging improves the reliability of the animal-level ",
           "estimate. (No external truth is available, so accuracy is not ",
           "established.)")),
 CC("WGCNA modules show no phenotype effect.", "ED_WGCNA b",
    "module eigengene difference, 45 cells", "inferential (negative)", FALSE,
    FALSE,
    paste0("No module-by-contrast cell survived multiple-testing correction ",
           "(0 of 45), and no stress-by-spatial-unit interaction test did ",
           "either (0 of 35, smallest FDR 0.27).")),
 CC("Stress does not alter the molecular network.", "ED8c",
    "exact permutation, whole network", "inferential (negative)", FALSE, FALSE,
    paste0("We did not detect a whole-network group difference; the exact ",
           "enumeration had resolution to p = 0.0036.")),
 CC("Proteins redistribute away from their baseline locations.", "ED7",
    "cross-tabulation of control peak versus effect unit", "descriptive",
    TRUE, FALSE,
    paste0("Among the robustness-qualified proteins, the strongest ",
           "phenotype-associated difference occurred outside the dominant ",
           "control spatial unit in 15 of 15 cases, and outside the broader ",
           "control affinity set in 14 of 15. Nothing moves: these are two ",
           "independent measurements of the same protein.")))
write_csv_safe(cc, file.path(OUT, "claim_chain_audit.csv"))

md <- c("# Claim-chain audit", "",
  paste0("Sentence-level logic check. Each row is an intended Results ",
         "sentence, the panel it rests on, and a corrected version where the ",
         "original says more than the analysis establishes. Biological n is 9 ",
         "animals, 3 per group, throughout."), "")
for (i in seq_len(nrow(cc))) md <- c(md,
  sprintf("### %d. %s", i, cc$sentence[i]), "",
  sprintf("- **Panel:** %s", cc$figure_panel[i]),
  sprintf("- **Statistical object:** %s", cc$statistical_object[i]),
  sprintf("- **Mode:** %s | **Selected from the same analysis:** %s",
          cc$descriptive_or_inferential[i], cc$selected_from_same_analysis[i]),
  sprintf("- **Verb strength appropriate:** %s",
          cc$verb_strength_appropriate[i]),
  sprintf("- **Corrected:** %s", cc$corrected_sentence[i]), "")

# ------------------------------------------------------------- S27 story test
md <- c(md, "# The core story, tested sentence by sentence", "",
 "## 1. \"The hippocampus exhibits a reproducible multi-resolution molecular geography.\"", "",
 paste0("**reproducible** - supported, but the evidence is left-right ",
        "concordance WITHIN animals (Pearson 0.50-0.92). Hemispheres are ",
        "repeated tissue, so this is internal consistency, not independent ",
        "replication. Qualify it."),
 paste0("**multi-resolution** - supported, but not uniformly: only neuropil ",
        "carries region x layer. Soma and microglia-enriched ROI are ",
        "region-level."),
 paste0("**geography** - a metaphor. The rest of the package says ",
        "architecture; use one word consistently."), "",
 "## 2. \"Later stress outcome is associated with sparse individual-protein changes but richer spatially restricted coordinated molecular-program differences.\"", "",
 paste0("**stress outcome** - correct. SUS-RES compares outcomes among ",
        "stress-exposed animals."),
 paste0("**sparse** - supported: 37 FDR-supported proteins over 18 units, 12 ",
        "of 18 units with none."),
 paste0("**richer** - a loose comparison between two different statistical ",
        "objects, protein counts against enrichment statistics. State the ",
        "actual asymmetry instead."),
 "**spatially restricted** - supported.",
 "**coordinated** - fair; gene-set enrichment is about coordinated sets.",
 paste0("**changes** - implies a before-and-after. SUS-RES is a between-group ",
        "difference. Use differences."), "",
 "## 3. \"Representative effects span CA3-SR neuropil synaptic signalling, CA2 neuronal-soma RNA processing and CA1 microglia-enriched mitochondrial respiration.\"", "",
 "**representative** - correct and important; keep it.",
 paste0("**RNA processing** - WRONG LEVEL. The exact term for the CA2 soma ",
        "panel is mRNA processing (GO:0006397). RNA processing is the broader ",
        "44-term atlas umbrella."),
 paste0("**mitochondrial respiration** - WRONG LEVEL for the direct panel. The ",
        "exact term is oxidative phosphorylation (GO:0006119). Mitochondrial ",
        "respiration is the 20-term umbrella, which also contains glycolytic ",
        "terms."),
 "**spatial levels** - CA3-SR is region x layer; CA2 soma and CA1 microglia-enriched ROI are region-level.", "",
 "## 4. \"These effects occur without FDR-supported WGCNA phenotype effects or detectable global molecular-network reorganization.\"", "",
 paste0("**without** - overstates. Absence of FDR support is not absence of ",
        "effect, and with 3 animals per group power is limited."),
 "**detectable** - helps, and should be applied to both halves.",
 paste0("**global** - the test is a whole-network comparison per compartment. ",
        "Whole-network is the exact word."),
 paste0("**reorganization** - acceptable ONLY as the name of the tested ",
        "network metric, in its null form."), "",
 "# Corrected preferred story", "",
 paste0("> The hippocampus exhibits a bilaterally reproducible, ",
        "multi-resolution spatial molecular architecture, resolved to region ",
        "and layer in neuropil and to region in the neuronal-soma and ",
        "microglia-enriched compartments. Later stress outcome is associated ",
        "with sparse individual-protein differences - 37 proteins reach FDR ",
        "support across 18 spatial units, and 12 of those units contain none - ",
        "but with coordinated molecular-program differences that remain ",
        "detectable in units where no individual protein does. Representative ",
        "FDR-supported examples, one per measurement compartment, are a ",
        "synaptic-signalling program in CA3-SR neuropil, an mRNA-processing ",
        "program in the CA2 neuronal soma, and an oxidative-phosphorylation ",
        "program in the CA1 microglia-enriched ROI. These differences occur ",
        "alongside module-phenotype tests that do not survive multiple-testing ",
        "correction (0 of 45 module-by-contrast cells; 0 of 35 ",
        "stress-by-spatial-unit interaction tests, smallest FDR 0.27) and no ",
        "detectable whole-network group difference."), "",
 paste0("_Every number in that paragraph is checked against the source data by ",
        "`figures/final_truth_v9_claim_audit.R`._"))
writeLines(md, file.path(REP, "claim_chain_audit.md"))

# ------------------------------------------------------------- S28 rules
rules <- c("# Manuscript semantic rules", "",
  paste0("Derived from the canonical ontology mapping and the frozen figure ",
         "layer. Each entry is USE / AVOID / ONLY USE WHEN."), "")
RULE <- function(topic, use, avoid, only)
  c(sprintf("## %s", topic), "", sprintf("**USE:** %s", use),
    sprintf("**AVOID:** %s", avoid), sprintf("**ONLY USE WHEN:** %s", only), "")
rules <- c(rules,
 RULE("baseline", "control spatial profile; dominant control spatial unit; abundance in control animals",
      "baseline, where it could imply a pre-stress longitudinal measurement",
      "baseline is explicitly defined in the same sentence as the control-reference state. Nothing was measured pre-stress."),
 RULE("stress effect", "outcome-associated; susceptibility-associated; associated with later stress outcome",
      "the stress effect, for SUS-RES",
      "the contrast really is RES-CON or SUS-CON, which compare stress-exposed animals with controls."),
 RULE("trajectory", "three-contrast pattern; directional pattern across the three pairwise contrasts",
      "trajectory, molecular trajectory, three independent trajectories",
      "describing the ordering informally, and the algebraic relation between the contrasts is stated nearby."),
 RULE("divergent", "susceptibility-associated; directionally distinct RES and SUS contrasts",
      "divergent, RES and SUS diverge, divergent resilience program",
      "both arms are FDR-supported in opposite directions, which is not the case for the synaptic or mRNA exemplars."),
 RULE("susceptibility-specific", "susceptibility-associated; resilience-associated",
      "susceptibility-specific, resilience-specific",
      "a specificity analysis was actually performed."),
 RULE("pathway", "molecular program; biological program; GO-defined process",
      "pathway, for GO biological-process enrichment",
      "the gene set really is a curated pathway, for example KEGG or Reactome."),
 RULE("validation", "external validation (F2g and ED2b only); functional characterisation; decomposition; external context",
      "validation used generically for F2d, F2h, F3 g/h/i or EWCE",
      "the comparison is against genuinely external, independent evidence."),
 RULE("microglial", "microglia-enriched ROI; microglia-enriched local microenvironment",
      "microglial proteome, microglia-specific, cell-intrinsic, cell-autonomous",
      "referring to external annotation or to the literature, never to a measured intrinsic cell state."),
 RULE("mitochondrial dysfunction", "reduced oxidative-phosphorylation program; reduced mitochondrial respiration program",
      "mitochondrial dysfunction, impaired mitochondrial function, defective OXPHOS, metabolic dysfunction",
      "a functional or respirometric measurement exists. Protein abundance is not function."),
 RULE("synaptic dysfunction", "reduced synaptic-signalling program",
      "synaptic dysfunction, impaired synapses, reduced synaptic transmission",
      "a physiological measurement exists."),
 RULE("reprogramming / rewiring", "program-level differences; spatially restricted molecular differences",
      "reprogramming, rewiring, redistribution, relocation, remodelling",
      "never, on the present evidence. Reorganization may name the tested network metric, in its null form."),
 RULE("stable", "largely preserved spatial molecular architecture; spatial organisation remained evident across groups",
      "unchanged architecture; stable, presented as a tested result",
      "pointing at the bilateral reproducibility and the FDR-negative network comparison, both named."),
 RULE("significant theme", "theme-level summary; theme containing one or more FDR-supported canonical GO terms",
      "FDR-significant theme, theme p-value, theme-level significance",
      "never. Theme aggregation is descriptive and has no FDR family."),
 RULE("leading-edge proteins", "selected leading-edge proteins; leading-edge contributors; proteins contributing to the enrichment",
      "validated, significant, key or driver proteins",
      "never for F3 g/h/i. None of the 63 displayed values is individually FDR-supported."),
 RULE("strongest program", "representative; illustrative; selected FDR-supported example; compartment-specific exemplar",
      "strongest, top, dominant, most altered, major response",
      "never for the three exemplars. The synaptic exemplar ranks 51st of 60 FDR-supported terms within CA3-SR itself."),
 RULE("OXPHOS", "oxidative phosphorylation, for the exact GO:0006119 term",
      "OXPHOS as the name of the 20-term atlas theme",
      "referring to GO:0006119 specifically. The theme is mitochondrial respiration / OXPHOS and also contains glycolytic terms."),
 RULE("mRNA versus RNA processing", "mRNA processing for the exact GO:0006397 term; RNA processing / splicing / RNP organization for the umbrella",
      "silently broadening mRNA processing to RNA biology",
      "the umbrella is meant, and the broadening is stated. The umbrella genuinely includes rRNA and miRNA processing."),
 RULE("null results", "no ... survived multiple-testing correction; we did not detect ...; no detectable ... at the present sample size",
      "no effect, no difference, unchanged, absent, equivalent",
      "never claim absence or equivalence. With 3 animals per group, power is limited."))
writeLines(rules, file.path(REP, "manuscript_semantic_rules.md"))

cat("claim-chain sentences:", nrow(cc), " corrected:",
    sum(cc$corrected_sentence != "NONE - accurate as written"), "\n")
cat("semantic rules        :", sum(grepl("^## ", rules)), "\n")

# ======================================== S5/S6/S9/S29 reader-facing term scan
#
# The registries above declare what each thing MAY be called. This section is
# the enforcement half: it reads every reader-facing artefact this layer emits
# and reports each place where an audited term appears WITHOUT the qualifier
# that licenses it. It is a gate, not a report - a surviving P0 stops the build.
#
# Prohibition text is exempt by construction. The rules file, the claim-chain
# audit and any line that tells the reader NOT to say something necessarily
# contain the banned phrase; matching those is the same false positive the
# Part-26 banned-unit check made against the schematic's own disclaimer.

corpus_files <- function() {
  roots <- c(path_results("reports", "manuscript_candidates", "final_truth_v9"),
             path_results("tables", "manuscript_candidates", "final_truth_v9"),
             path_results("figures", "manuscript_candidates", "final_truth_v9"))
  f <- unlist(lapply(roots[dir.exists(roots)], list.files, pattern = "[.]md$",
                     recursive = TRUE, full.names = TRUE))
  sort(c(f, repo_path("figures", "figure_final_truth_v9_contract.yml")))
}

# files whose whole purpose is to quote the prohibited wording
EXEMPT_FILE <- c("manuscript_semantic_rules.md", "claim_chain_audit.md",
                 "semantic_search_hits.csv")
# a line that instructs against the phrase rather than using it
PROHIBITION <- paste0("never|must not|do not |does not|cannot|prohibited|",
                      "instead of|not be read|not be interpreted|",
                      "rather than|avoid|no hypothesis|banned|->|forbidden|",
                      "over-read|overread|misread|incorrect|wrong|",
                      "deliberately|is treated as|no claim|not a claim")
# the repository capitalises NOT and NEVER exactly where it disclaims
EMPHATIC <- "\\bNOT\\b|\\bNEVER\\b"
disclaimed <- function(ctx)
  grepl(PROHIBITION, ctx, ignore.case = TRUE, perl = TRUE) ||
  grepl(EMPHATIC, ctx, perl = TRUE)

scan_corpus <- function(spec) {
  hits <- list()
  for (f in corpus_files()) {
    if (basename(f) %in% EXEMPT_FILE) next
    txt <- readLines(f, warn = FALSE, encoding = "UTF-8")
    rel <- sub(".*proteomics[/\\]", "", f)
    for (k in seq_len(nrow(spec))) {
      i <- grep(spec$pattern[k], txt, ignore.case = TRUE, perl = TRUE)
      heads <- grep("^#", txt)
      for (ln in i) {
        line <- txt[ln]
        h <- heads[heads <= ln]
        ctx <- paste(c(if (length(h)) txt[max(h)],
                       txt[max(1L, ln - 2L):min(length(txt), ln + 1L)]),
                     collapse = " ")
        if (disclaimed(ctx)) next
        # a term directly preceded by a negation is being disclaimed, not used
        pos <- regexpr(spec$pattern[k], line, ignore.case = TRUE, perl = TRUE)
        pre <- substr(line, max(1L, pos - 40L), max(1L, pos - 1L))
        if (grepl("(^|[^[:alpha:]])(not|never|no|neither|nor|without)([^[:alpha:]]|$)",
                  pre, ignore.case = TRUE, perl = TRUE)) next
        licensed <- nzchar(spec$licence[k]) &&
          grepl(spec$licence[k], line, ignore.case = TRUE, perl = TRUE)
        hits[[length(hits) + 1L]] <- data.frame(
          file = rel, line_or_field = ln, term = spec$term[k],
          context = substr(trimws(sub("^[#>*[:space:]-]+", "", line)), 1, 160),
          licensed_by = if (licensed) spec$licence_name[k] else "NONE",
          severity = if (licensed) "OK" else spec$severity[k],
          recommended_fix = if (licensed) "none - qualifier present in the same sentence"
                            else spec$fix[k],
          stringsAsFactors = FALSE)
      }
    }
  }
  empty <- data.frame(file = character(), line_or_field = integer(),
                      term = character(), context = character(),
                      licensed_by = character(), severity = character(),
                      recommended_fix = character(), stringsAsFactors = FALSE)
  if (!length(hits)) return(empty)
  out <- do.call(rbind, hits)
  out[order(out$severity, out$file, out$line_or_field, out$term), , drop = FALSE]
}

# a clean scan legitimately returns zero rows; add columns at the frame length
addcol <- function(df, name, value) {
  df[[name]] <- if (nrow(df)) value else vector(mode = mode(value), 0L)
  df
}
S <- function(term, pattern, licence_name, licence, severity, fix)
  data.frame(term = term, pattern = pattern, licence_name = licence_name,
             licence = licence, severity = severity, fix = fix,
             stringsAsFactors = FALSE)

# ------------------------------------------------------------ S5 stress language
#
# Both RES and SUS animals were stressed. A SUS-RES contrast therefore separates
# OUTCOMES within stress exposure, not stressed from unstressed. Only RES-CON
# and SUS-CON compare a stress-exposed group with an unexposed one.
stress_spec <- rbind(
  S("stress effect", "stress[ -]?(effect|response|signature|change)",
    "contrast named as CON-referenced", "RES[^a-z]?(-|\\u2212|vs)[^a-z]?CON|SUS[^a-z]?(-|\\u2212|vs)[^a-z]?CON|versus control|vs control",
    "P0 factual",
    "SUS-RES separates outcomes within stress exposure; say outcome-associated or susceptibility-associated"),
  S("stress-induced", "stress[ -]?induced", "none", "",
    "P0 factual",
    "induction was not measured longitudinally; say associated with later stress outcome"),
  S("stressed/unstressed dichotomy", "un[ -]?stressed|non[ -]?stressed",
    "explicit CON definition", "CON|control",
    "P1 overclaim",
    "name the group (CON) rather than implying an unexposed comparison within the stressed arms"),
  S("effect of stress", "effect of (chronic )?stress|impact of stress",
    "contrast named as CON-referenced", "RES|SUS.{0,12}CON|control",
    "P1 overclaim",
    "attribute the difference to the named contrast, not to stress as an agent"))
sl <- scan_corpus(stress_spec)
# the S5 column set, emitted whether or not anything is flagged
sl <- addcol(sl, "phrase", sl$term)
sl <- addcol(sl, "location", paste0(sl$file, ":", sl$line_or_field))
sl$underlying_contrast <- ifelse(
  grepl("SUS.{0,3}(\u2212|-).{0,3}RES", sl$context),
  "SUS-RES: both arms stress-exposed, so this is an OUTCOME contrast",
  ifelse(grepl("RES.{0,3}(\u2212|-).{0,3}CON|SUS.{0,3}(\u2212|-).{0,3}CON",
                sl$context),
         "CON-referenced: genuinely compares stress-exposed with unexposed",
         "no contrast named in the sentence"))
sl <- addcol(sl, "current_wording", sl$context)
sl <- addcol(sl, "accurate", sl$severity == "OK")
sl <- addcol(sl, "replacement", ifelse(sl$accurate, "", sl$recommended_fix))
write_csv_safe(sl, file.path(OUT, "stress_language_audit.csv"))

# ---------------------------------------------------------- S6 baseline language
#
# No pre-stress measurement exists. Every animal was sampled once, after the
# paradigm. "Baseline" is only legal as a synonym for the control-group state
# and only when that identity is stated in the same sentence.
baseline_spec <- rbind(
  S("baseline", "baseline",
    "control identity stated in the same sentence", "CON\\b|control",
    "P1 overclaim",
    "define baseline as the control-group state in the same sentence, or say in control animals"),
  S("pre-stress", "pre[ -]?stress|before stress|prior to stress", "none", "",
    "P0 factual",
    "no pre-stress sample exists; remove the temporal reference entirely"),
  S("at baseline (temporal)", "at baseline",
    "control identity stated in the same sentence", "CON\\b|control",
    "P0 factual",
    "at baseline reads as a timepoint; say in control animals"))
bl <- scan_corpus(baseline_spec)
bl <- addcol(bl, "design_fact",
             "single terminal sampling; no longitudinal within-animal reference")
bl <- addcol(bl, "measured_pre_stress", FALSE)
write_csv_safe(bl, file.path(OUT, "baseline_language_audit.csv"))

# --------------------------------------------- S9 phenotype specificity language
#
# No specificity test was performed anywhere in this layer. Establishing that an
# effect is specific to SUS requires testing the SUS contrast against the RES
# contrast and showing the difference is itself supported; what exists is a set
# of separately FDR-assessed contrasts.
spec_spec <- rbind(
  S("susceptibility/resilience-specific",
    "susceptibility[ -]specific|resilience[ -]specific|SUS[ -]specific|RES[ -]specific",
    "none", "", "P0 factual",
    "no specificity test was performed; say susceptibility-associated or resilience-associated"),
  S("only in / unique to", "only in (SUS|RES|susceptible|resilient)|unique to (SUS|RES)",
    "count framing", "FDR-supported|of 18|of 15|survived",
    "P1 overclaim",
    "state the count of units in which the contrast was FDR-supported rather than exclusivity"),
  S("selectively / exclusively", "selectively|exclusively",
    "none", "", "P1 overclaim",
    "absence of support is not evidence of absence at n = 3 per group"),
  S("divergent", "divergen",
    "both arms supported in opposite directions", "opposite direction|both.{0,20}FDR",
    "P2 consistency",
    "reserve divergent for cases where both arms are FDR-supported in opposite directions"),
  S("dysfunction", "dysfunction|impair|deficit|defective",
    "none", "", "P0 factual",
    "protein abundance is not function; say reduced <program> program"))
ps <- scan_corpus(spec_spec)
ps <- addcol(ps, "test_that_would_license_it", ifelse(
  grepl("specific", ps$term),
  "an explicit SUS-versus-RES contrast-of-contrasts with its own FDR control",
  "a functional or physiological measurement, which this study does not contain"))
write_csv_safe(ps, file.path(OUT, "phenotype_specificity_language_audit.csv"))

# -------------------------------------------------------- S29 consolidated hits
KEEP <- c("file", "line_or_field", "term", "context", "severity",
          "recommended_fix")
hits <- rbind(sl[, KEEP, drop = FALSE], bl[, KEEP, drop = FALSE],
              ps[, KEEP, drop = FALSE])
hits <- hits[hits$severity != "OK", , drop = FALSE]
hits <- hits[order(hits$severity, hits$file, hits$line_or_field), , drop = FALSE]
write_csv_safe(hits, file.path(OUT, "semantic_search_hits.csv"))

p0 <- hits[grepl("^P0", hits$severity), , drop = FALSE]
if (nrow(p0)) {
  print(p0)
  stop("semantic scan: ", nrow(p0), " P0 factual hit(s) survive in ",
       length(unique(p0$file)), " reader-facing artefact(s)")
}

# ------------------------------------------- S29b what is PRINTED on the panel
#
# The scan above reads the prose a reader is given alongside the figure. This
# one reads the text a reader actually sees ON it, extracted from the panel SVG
# text nodes, so an axis label cannot quietly contradict its own legend.
#
# Severity is rescaled for printed text. A short label has no room to carry its
# own qualifier and legitimately borrows it from the legend, so an unqualified
# term is a CONSISTENCY issue there; but a label that asserts a fact the design
# cannot support - a pre-stress reference, a specificity claim, a functional
# claim - is printed on the page and stays P0.
svg_text <- function(f) {
  x <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  n <- unlist(regmatches(x, gregexpr("<text[^>]*>[^<]*</text>", x)))
  t <- sub(".*>([^<]*)</text>", "\\1", n)
  t <- gsub("&amp;", "&", gsub("&lt;", "<", gsub("&gt;", ">", t)))
  trimws(t[nzchar(trimws(t))])
}

scan_panels <- function(spec) {
  root <- path_results("figures", "manuscript_candidates", "final_truth_v9")
  svgs <- sort(list.files(root, pattern = "[.]svg$", recursive = TRUE,
                          full.names = TRUE))
  svgs <- svgs[grepl("[/\\]panels[/\\]", svgs)]
  hits <- list()
  for (f in svgs) {
    labs <- svg_text(f)
    rel <- sub(".*proteomics[/\\]", "", f)
    for (k in seq_len(nrow(spec))) {
      j <- grep(spec$pattern[k], labs, ignore.case = TRUE, perl = TRUE)
      for (m in j) {
        lab <- labs[m]
        if (disclaimed(lab)) next
        # a multi-line label is emitted as consecutive text nodes, so the
        # qualifier of "Baseline abundance / (CON z-score)" sits in the NEXT
        # node; licence against that window, not the single node.
        win <- paste(labs[max(1L, m - 1L):min(length(labs), m + 1L)],
                     collapse = " ")
        if (nzchar(spec$licence[k]) &&
            grepl(spec$licence[k], win, ignore.case = TRUE, perl = TRUE)) next
        sev <- if (grepl("^P0", spec$severity[k])) spec$severity[k]
               else "P2 consistency"
        hits[[length(hits) + 1L]] <- data.frame(
          file = rel, line_or_field = paste0("text node ", m),
          term = spec$term[k], context = substr(lab, 1, 160), severity = sev,
          recommended_fix = if (identical(sev, "P2 consistency"))
              paste0("printed label; qualifier lives in the legend - ",
                     spec$fix[k])
            else spec$fix[k],
          stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(hits))
    return(data.frame(file = character(), line_or_field = character(),
                      term = character(), context = character(),
                      severity = character(), recommended_fix = character(),
                      stringsAsFactors = FALSE))
  out <- do.call(rbind, hits)
  out[order(out$severity, out$file, out$term), , drop = FALSE]
}

pn <- rbind(scan_panels(stress_spec), scan_panels(baseline_spec),
            scan_panels(spec_spec))
pn <- pn[!duplicated(pn[, c("file", "line_or_field", "term")]), , drop = FALSE]
write_csv_safe(pn, file.path(OUT, "printed_panel_language_audit.csv"))

hits <- rbind(hits, pn[, KEEP, drop = FALSE])
hits <- hits[order(hits$severity, hits$file, hits$line_or_field), , drop = FALSE]
write_csv_safe(hits, file.path(OUT, "semantic_search_hits.csv"))

p0 <- hits[grepl("^P0", hits$severity), , drop = FALSE]
if (nrow(p0)) {
  print(p0)
  stop("semantic scan: ", nrow(p0), " P0 factual hit(s) survive")
}
cat("printed panel labels :", nrow(pn), "unresolved on-figure hits\n")
cat("semantic scan (total) :", nrow(hits), "unresolved (",
    sum(grepl("^P1", hits$severity)), "P1,", sum(grepl("^P2", hits$severity)),
    "P2 ); 0 P0\n")

# ================================================= S10 program versus pathway
#
# "Pathway" is only honest for a gene set that a curated pathway database
# defines as one. The enrichment evidence in this project is read here rather
# than assumed, so the verdict cannot drift if a collection is ever added.
ev_families <- sort(unique(TH$evidence_source_family))
go_frac <- mean(grepl("^GO:", TH$GO_ID))
pw <- data.frame(
  gene_set_collection = c(ev_families, "KEGG", "Reactome", "WikiPathways"),
  present_in_this_project = c(rep(TRUE, length(ev_families)), FALSE, FALSE,
                              FALSE),
  definition = c(
    rep("GO biological process: a term in a structured ontology of processes, not a curated reaction network",
        length(ev_families)),
    rep("curated reaction network with defined membership and topology", 3L)),
  pathway_label_justified = c(rep(FALSE, length(ev_families)),
                              rep(NA, 3L)),
  required_label = c(
    rep("molecular program / biological program / GO-defined process / enriched biological process",
        length(ev_families)),
    rep("pathway would be justified here, but no such collection is used",
        3L)),
  stringsAsFactors = FALSE)
pw$evidence <- sprintf(
  "%.1f%% of the %d theme-assignment rows carry a GO: identifier; no KEGG, Reactome or WikiPathways collection appears anywhere in pipeline.yml",
  100 * go_frac, nrow(TH))
write_csv_safe(pw, file.path(OUT, "program_vs_pathway_audit.csv"))

path_hits <- scan_corpus(rbind(
  S("pathway", "pathway", "none", "", "P1 overclaim",
    "the enrichment evidence is GO biological process, not a curated pathway database; say molecular program")))
path_hits <- path_hits[path_hits$severity != "OK", , drop = FALSE]
if (nrow(path_hits)) {
  hits <- rbind(hits, path_hits[, KEEP, drop = FALSE])
  hits <- hits[order(hits$severity, hits$file, hits$line_or_field), ,
               drop = FALSE]
  write_csv_safe(hits, file.path(OUT, "semantic_search_hits.csv"))
  if (any(grepl("^P0", path_hits$severity)))
    stop("semantic scan: P0 pathway hit survives", call. = FALSE)
}
cat("pathway/program      : GO rows", sprintf("%.1f%%", 100 * go_frac),
    "| curated pathway collections used: 0 | unjustified 'pathway' in prose:",
    nrow(path_hits), "\n")

# ====================================================== S27 the core story
#
# Each clause of the proposed four-sentence story is checked against the one
# artefact that would have to support it. A clause with no such artefact is not
# softened here - it is removed, because softening an unsupported clause keeps
# the claim and only hides the evidence gap.
CS <- function(sentence, clause, supported, evidence, verdict, corrected)
  data.frame(sentence = sentence, clause = clause, supported = supported,
             supporting_artefact = evidence, verdict = verdict,
             corrected_clause = corrected, stringsAsFactors = FALSE)

cs <- rbind(
 CS(1, "reproducible", TRUE,
    "F2f / ED1 bilateral audit: stored left-right concordance for every prespecified anatomical contrast",
    "keep", "reproducible"),
 CS(1, "multi-resolution", TRUE,
    "F2a schematic and spatial_resolution_semantics.csv: 10 of 18 units carry region x layer resolution, 8 are region-level",
    "keep, but the asymmetry must be stated once",
    "multi-resolution, laminar in the neuropil and region-level in the neuronal soma and microglia-enriched ROI"),
 CS(1, "geography", TRUE,
    "F2d fingerprint: CON-only, phenotype-blind, ordered by peak spatial unit in control animals",
    "keep as a descriptive metaphor", "molecular geography"),
 CS(2, "later stress outcome", TRUE,
    "the SUS/RES split is a behavioural outcome measured after the paradigm",
    "keep; it is the outcome, not an exposure contrast",
    "later stress outcome"),
 CS(2, "sparse", TRUE,
    "F3a DAP track: 37 FDR-supported SUS-RES proteins, 12 of 18 spatial units with none",
    "keep, with the count", "sparse individual-protein differences"),
 CS(2, "richer", FALSE,
    "no test compares protein-level and program-level evidence on a common scale; the two have different multiple-testing families",
    "REMOVE - this is a cross-family comparison that was never made",
    "and with differences at the level of molecular programs, with both inventories stated and neither claimed to be the stronger"),
 CS(2, "spatially restricted", TRUE,
    "the theme atlas and ED6: FDR-supported terms occur in a subset of units",
    "keep", "spatially restricted"),
 CS(2, "coordinated", TRUE,
    "GSEA operates on a ranked gene set, so coordination is what the statistic measures",
    "keep", "coordinated"),
 CS(3, "CA3-SR neuropil synaptic signalling", TRUE,
    "GO:0099536 synaptic signaling, FDR-supported in the CA3-SR neuropil unit",
    "keep the exact term; the unit is genuinely laminar",
    "CA3 stratum radiatum neuropil synaptic signalling"),
 CS(3, "CA2 neuronal-soma RNA processing", TRUE,
    "GO:0006397 mRNA processing, FDR-supported in the CA2 neuronal-soma compartment",
    "narrow to the exact term; RNA processing is the umbrella",
    "CA2 neuronal-soma mRNA processing"),
 CS(3, "CA1 microglia-enriched mitochondrial respiration", TRUE,
    "GO:0006119 oxidative phosphorylation, FDR-supported in the CA1 microglia-enriched ROI",
    "keep the umbrella at story level; the exact term is narrower",
    "CA1 microglia-enriched ROI mitochondrial respiration / oxidative phosphorylation"),
 CS(3, "representative", TRUE,
    "f3_program_example_selection.csv: a hardcoded one-per-compartment choice made after the contrasts were known; the synaptic exemplar is not the top-ranked term in its own unit",
    "keep ONLY as representative; strongest is prohibited",
    "representative FDR-supported examples"),
 CS(4, "without", FALSE,
    "0 of 45 and 0 of 35 WGCNA tests survived correction; a non-detection at n = 3 per group is not an absence",
    "REMOVE - without asserts absence",
    "and were not accompanied by"),
 CS(4, "detectable", TRUE,
    "ST6 / ST7 report the attainable resolution alongside the null",
    "keep; it is what converts the null into a limited-power statement",
    "any detectable"),
 CS(4, "global", FALSE,
    "ST7 tests the whole graph; ST8 tested 8 neuropil edges only, not every unit pair",
    "QUALIFY - global is true of the whole-graph test but not of the edge-coupling test",
    "whole-network"))
cs$sentence_text <- c(
  "The hippocampus exhibits a reproducible multi-resolution molecular geography.",
  "Later stress outcome is associated with sparse individual-protein changes but richer spatially restricted coordinated molecular-program differences.",
  "Representative effects span CA3-SR neuropil synaptic signalling, CA2 neuronal-soma RNA processing and CA1 microglia-enriched mitochondrial respiration.",
  "These effects occur without FDR-supported WGCNA phenotype effects or detectable global molecular-network reorganization."
)[cs$sentence]
write_csv_safe(cs, file.path(OUT, "core_story_audit.csv"))

STORY <- c(
"# Corrected core story",
"",
"Audited clause by clause against the artefact that would have to support it;",
"see core_story_audit.csv. Three clauses were not supportable and are not",
"softened but removed: 'richer' compares two multiple-testing families that",
"were never placed on a common scale, 'without' asserts an absence that a",
"non-detection at three animals per group cannot establish, and 'global'",
"is true of the whole-graph test but not of the edge-coupling test, which",
"covered eight neuropil edges rather than every spatial-unit pair.",
"",
"## Preferred version",
"",
"The hippocampal proteome is organised as a reproducible molecular geography",
"that is resolved at the laminar level in the neuropil and at the region level",
"in the neuronal soma and microglia-enriched ROI. Later stress outcome is",
"associated with sparse individual-protein differences - 37 FDR-supported",
"SUS-RES proteins, with none at all in 12 of the 18 spatial units - and with",
"coordinated, spatially restricted differences at the level of molecular",
"programs. The two are assessed in separate multiple-testing families and are",
"not placed on a common scale, so neither is claimed to be the stronger.",
"Representative FDR-supported examples span synaptic signalling in",
"the CA3 stratum radiatum neuropil, mRNA processing in the CA2 neuronal soma,",
"and mitochondrial respiration / oxidative phosphorylation in the CA1",
"microglia-enriched ROI; these are illustrative, one per measurement",
"compartment, and are not the strongest program in their own unit. No WGCNA",
"module-phenotype association survived multiple-testing correction (0 of 45",
"and 0 of 35), and we did not detect a whole-network group difference or an",
"FDR-supported edge-behaviour association at the present sample size.",
"",
"## What the preferred version deliberately does not say",
"",
"- that program-level evidence is stronger than protein-level evidence;",
"- that WGCNA or network architecture is unchanged;",
"- that every spatial-unit pair was tested for behavioural coupling;",
"- that any exemplar is the dominant or strongest program;",
"- that any measured quantity is a function rather than an abundance.")
writeLines(STORY, file.path(REP, "core_story_corrected.md"))
cat("core story           :", nrow(cs), "clauses audited,",
    sum(!cs$supported), "removed or qualified\n")

# ============================== S11/S14/S15/S16/S17/S18/S19 the remaining rules
#
# The rules document declares these; this scan enforces them, so a rule cannot
# be stated in one artefact and broken in another.
other_spec <- rbind(
  # S11 cell type: every compartment is an enriched ROI or a soma-enriched
  # capture, never a sorted or single-cell population
  S("cell-intrinsic", "cell[ -]intrinsic|cell[ -]autonomous|microglia[ -]specific",
    "none", "", "P0 factual",
    "an enriched ROI cannot establish a cell-intrinsic mechanism; say microglia-enriched ROI"),
  S("microglial proteome", "microglial proteome|neuronal proteome",
    "none", "", "P1 overclaim",
    "say microglia-enriched ROI or neuronal-soma compartment"),
  # S14 exemplar logic: the three Figure-3 programs were chosen one per
  # compartment after the contrasts were known, and the synaptic exemplar is
  # not the top-ranked FDR-supported term in its own unit
  S("strongest program", "strongest (program|pathway|theme)|top (program|pathway|theme)|dominant program|most altered|major response",
    "none", "", "P0 factual",
    "these are illustrative exemplars, not the strongest; say representative FDR-supported example"),
  # S15 theme aggregation is descriptive and carries no FDR family of its own
  S("theme-level significance",
    "FDR[ -]significant theme|theme[ -]level significance|theme p[ -]?value|significant theme",
    "none", "", "P0 factual",
    "no theme-level test exists; say theme containing one or more FDR-supported canonical GO terms"),
  # S16 none of the displayed leading-edge proteins is individually supported
  S("leading-edge overclaim",
    "validated protein|significant protein|key protein|driver protein",
    "none", "", "P0 factual",
    "say selected leading-edge proteins; none is individually FDR-supported"),
  # S17/S18 a non-detection at three animals per group is not an absence
  S("absence claim",
    "no effect\\b|no difference\\b|unchanged|\\bequivalent\\b|\\babsent\\b",
    "a resolution or correction is named, or the subject is an artefact rather than a result",
    paste0("survived|correction|did not detect|detectable|attainable|power|",
           "sample size|label|table|canonical|file|artefact|unactivated|",
           "historical|untouched|promoted"),
    "P1 overclaim",
    "say no ... survived multiple-testing correction, or we did not detect ... at the present sample size"),
  S("stable architecture", "stabl[ey][a-z ]{0,40}(architecture|organisation|organization)|unchanged architecture",
    "named as reproducibility or as a null", "reproducib|bilateral|did not detect|survived",
    "P1 overclaim",
    "say largely preserved spatial molecular architecture, and name the evidence"),
  # S19 movement verbs describe a redistribution no design here can observe
  S("rewiring", "rewir|reprogramm|redistribut|relocat|remodell",
    "none", "", "P0 factual",
    "nothing was observed to move; say program-level or spatially restricted differences"),
  S("reorganization", "reorgani",
    "the tested network metric, in its null form",
    "did not detect|no detectable|survived|null",
    "P1 overclaim",
    "reorganization is admissible only for the tested network metric in its null form"))

oth <- scan_corpus(other_spec)
oth_bad <- oth[oth$severity != "OK", , drop = FALSE]
oth_panel <- scan_panels(other_spec)
write_csv_safe(oth, file.path(OUT, "claim_strength_language_audit.csv"))

if (nrow(oth_bad) || nrow(oth_panel)) {
  hits <- rbind(hits, oth_bad[, KEEP, drop = FALSE],
                oth_panel[, KEEP, drop = FALSE])
  hits <- hits[order(hits$severity, hits$file, hits$line_or_field), ,
               drop = FALSE]
  write_csv_safe(hits, file.path(OUT, "semantic_search_hits.csv"))
}
p0 <- hits[grepl("^P0", hits$severity), , drop = FALSE]
if (nrow(p0)) {
  print(p0[, c("file", "line_or_field", "term", "context")])
  stop("semantic scan: ", nrow(p0), " P0 factual hit(s) survive", call. = FALSE)
}
cat("claim-strength scan  :", nrow(oth_bad), "prose +", nrow(oth_panel),
    "printed unresolved\n")
cat("SEMANTIC SCAN FINAL  :", nrow(hits), "unresolved (",
    sum(grepl("^P1", hits$severity)), "P1,", sum(grepl("^P2", hits$severity)),
    "P2 ); 0 P0\n")
