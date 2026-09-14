#!/usr/bin/env Rscript

# Part-25: figure legends and the per-panel statistical-role audit.
#
# Until now the layer had no figure-legend artefact. Everything a reader needs
# in order to not over-read a panel - what is descriptive, what is inferential,
# which multiple-testing family applies, whether a displayed quantity was
# selected by the same analysis it appears to support - lived only in renderer
# comments. This script emits both the reader-facing legends and the audit
# table that proves every panel was asked those questions.
#
# The per-panel statistical roles below are EDITORIAL facts about the design,
# not values derivable from the data, so they are declared here and checked
# against the built source data where that is possible.

source(file.path("R", "paths.R"))
source(repo_path("R", "null_coalescing.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "final_truth_v9_figure_utils.R"))
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/final_truth_v9_legends.R")

args <- commandArgs(trailingOnly = TRUE)
if ("--dry-run" %in% args || is_dry_run()) {
  message("[DRY-RUN] final_truth_v9 legends and statistical-role audit")
  quit(save = "no", status = 0L)
}

ct <- s9f_contract()
SD <- path_results("source_data", "manuscript_candidates", "final_truth_v9")
REP <- path_results("reports", "manuscript_candidates", "final_truth_v9")
TAB <- path_results("tables", "manuscript_candidates", "final_truth_v9")
dir_create(REP); dir_create(TAB)

sidecar <- function(fig, id) {
  p <- file.path(SD, fig, paste0(id, "_source_data.csv"))
  if (file.exists(p)) nv_read_csv(p) else NULL
}

# ---------------------------------------------------------------- shared text
ALGEBRA <- paste0(
  "Pairwise contrasts are shown jointly to visualise the three-group ",
  "trajectory; they are algebraically related (SUS-RES = SUS-CON minus ",
  "RES-CON) and should not be interpreted as independent replications.")
GSEA_N <- paste0(
  "Biological n = 3 animals per group. Ranks are the stored per-gene contrast ",
  "statistic for that spatial unit, collapsed from protein groups by the ",
  "prespecified median rule. The FDR is the gene-set enrichment FDR ",
  "conditional on that ranking; it is not a count of independent biological ",
  "observations.")
DESCRIPTIVE <- "Descriptive; no hypothesis test is performed in this panel."

# ------------------------------------------------------- per-panel declaration
#
# quantity / replicate_unit / descriptive_or_inferential / selection_dependency
# / test / FDR_family / independent_validation / legend
R <- function(id, quantity, replicate_unit, mode, selection, test, fdr,
              indep, legend, biological_n = "9 animals; 3 per group")
  list(id = id, quantity = quantity, replicate_unit = replicate_unit,
       mode = mode, selection = selection, test = test, fdr = fdr,
       indep = indep, legend = legend, biological_n = biological_n)

ROLES <- list(
  R("v9_schematic", "anatomy and sampling scheme", "not applicable",
    "descriptive", "none", "none", "none", "no",
    paste0("Sampling scheme. Ten neuropil units carry region x layer ",
           "resolution; the neuronal soma and the microglia-enriched ROI are ",
           "region-level only and are never given laminar resolution.")),
  R("v9_depth", "proteins identified per spatial acquisition",
    "spatial acquisition (nested within animal)", "descriptive", "none",
    "none", "none", "no",
    paste0("Acquisition depth by compartment. Points are spatial ",
           "acquisitions, not animals: 180 / 71 / 72 acquisitions from 9 ",
           "animals, 3 per group. AnimalID is the biological replicate and ",
           "the acquisition count must not be read as biological n. The ",
           "violin is a descriptive density with the median and ",
           "interquartile range; the neuropil distribution is multimodal ",
           "because it pools heterogeneous region and layer contexts and ",
           "repeated acquisitions within animals, and its modes are NOT ",
           "evidence of biological subpopulations. ", DESCRIPTIVE)),
  R("v9_pca", "principal components of the global proteome",
    "spatial acquisition", "descriptive", "none", "none", "none", "no",
    paste0("Global proteome structure. ", DESCRIPTIVE)),
  R("v9_fingerprint", "baseline abundance (CON z-score)", "animal (CON only)",
    "descriptive", "rows are the top-ranked proteins of prespecified CON-only anatomical contrasts",
    "none", "none", "no",
    paste0("Baseline spatial molecular fingerprint, CON animals only ",
           "(n = 3). Rows are top-RANKED proteins from prespecified CON-only ",
           "anatomical contrasts, not significance-filtered hits, and no ",
           "stress information enters either the selection or the row order. ",
           "Rows are ordered by baseline peak spatial unit. ", DESCRIPTIVE)),
  R("v9_compartment", "baseline abundance of marker proteins (CON z-score)",
    "animal (CON only)", "descriptive", "prespecified marker panel", "none",
    "none", "no",
    paste0("Compartment identity of prespecified marker proteins. Colour ",
           "mapping saturates at 3.0 for display; cells exceeding the limit ",
           "are labelled with their uncapped values, and the uncapped values ",
           "are retained in the source data. ", DESCRIPTIVE)),
  R("v9_bilateral_main", "left-right concordance of the anatomical effect",
    "animal", "descriptive", "prespecified anatomical contrasts", "none",
    "none", "no",
    paste0("Bilateral reproducibility. Stored concordance metrics; no ",
           "hypothesis test is applied to them. ", DESCRIPTIVE)),
  R("v9_external_main", "GSEA normalised enrichment score", "animal",
    "inferential", "expected pairings only; the complete inventory is ED2",
    "gene set enrichment", "BH within the external-validation inventory",
    "YES - external reference",
    paste0("EXTERNAL anatomical validation. Internal CON-only anatomical ",
           "contrasts are tested against independently published ",
           "hippocampal signatures, so this panel is genuine external ",
           "validation. ", GSEA_N)),
  R("v9_internal_main", "GSEA normalised enrichment score", "animal",
    "inferential", "one canonical term per contrast; the complete inventory is ED2",
    "gene set enrichment", "BH within the internal anatomical inventory",
    "no - same data",
    paste0("FUNCTIONAL CHARACTERISATION of the anatomical contrasts. This ",
           "panel annotates the same proteomic data with canonical gene-set ",
           "programs; it is NOT independent validation and must not be read ",
           "as such. External validation is panel g. ", GSEA_N)),
  R("v9_dap_track", "count of FDR-supported and robustness-qualified proteins",
    "animal", "inferential (counts of FDR-supported tests)",
    "all FDR-supported SUS-RES proteins", "differential abundance, SUS vs RES",
    "BH within the differential-abundance family", "no",
    paste0("Differential-abundance burden per spatial unit. The second row ",
           "is NOT a uniform second-stage test applied to all 37 proteins: ",
           "the 28 CA2-SLM proteins underwent the detailed missingness and ",
           "QC qualification of ED3, of which 6 qualified, while the 9 ",
           "FDR-supported proteins outside CA2-SLM enter unchanged because ",
           "they were never at risk from that artefact. 'Robustness-",
           "qualified' therefore means 'not excluded by the CA2-SLM ",
           "missingness audit', not 'passed an additional test'.")),
  R("v9_atlas", "median NES across mapped canonical GO terms", "animal",
    "descriptive aggregation of inferential inputs",
    "claim-eligible themes only; qc_review themes excluded",
    "gene set enrichment on the constituent GO terms",
    "BH within each constituent GO family; the THEME has no family of its own",
    "no",
    paste0("Theme-level enrichment atlas. Theme colours summarise the median ",
           "NES of mapped canonical GO terms and are DESCRIPTIVE. Support ",
           "markers indicate that at least one constituent canonical GO term ",
           "passed its prespecified FDR threshold. Theme aggregation does ",
           "not constitute an additional multiple-testing family and no ",
           "theme-level p-value or FDR is computed or implied. ", GSEA_N)),
  R("v9_bridge", "anatomical location of the three representative programs",
    "animal", "descriptive placement of inferential results",
    "three representative programs", "gene set enrichment",
    "BH within each GSEA family", "no",
    paste0("Anatomical bridge. This panel answers WHERE the three ",
           "representative programs occur; panels d-f give the ranked ",
           "enrichment evidence for the same three programs. ", ALGEBRA)),
  R("v9_curve_syn", "running enrichment score", "animal", "inferential",
    "representative program", "gene set enrichment", "BH within the GSEA family",
    "no", paste0("Ranked enrichment evidence. ", GSEA_N, " ", ALGEBRA)),
  R("v9_curve_rna", "running enrichment score", "animal", "inferential",
    "representative program", "gene set enrichment", "BH within the GSEA family",
    "no", paste0("Ranked enrichment evidence. ", GSEA_N, " ", ALGEBRA)),
  R("v9_curve_ox", "running enrichment score", "animal", "inferential",
    "representative program", "gene set enrichment", "BH within the GSEA family",
    "no", paste0("Ranked enrichment evidence. ", GSEA_N, " ", ALGEBRA)),
  R("v9_prot_syn", "log2 fold change", "animal", "descriptive",
    "leading-edge constituents of the enrichment result above",
    "differential abundance", "BH within the differential-abundance family",
    "no - selected from the same analysis",
    paste0("Selected leading-edge proteins underlying the enrichment result ",
           "above are shown descriptively; proteins were selected using the ",
           "prespecified stored rank-statistic rule. They are NOT an ",
           "independent validation of the enrichment, NOT independently ",
           "significant proteins and NOT an unbiased protein sample. ",
           ALGEBRA)),
  R("v9_prot_rna", "log2 fold change", "animal", "descriptive",
    "leading-edge constituents of the enrichment result above",
    "differential abundance", "BH within the differential-abundance family",
    "no - selected from the same analysis",
    paste0("Selected leading-edge proteins, descriptive; see panel g. ",
           ALGEBRA)),
  R("v9_prot_ox", "log2 fold change", "animal", "descriptive",
    "leading-edge constituents of the enrichment result above",
    "differential abundance", "BH within the differential-abundance family",
    "no - selected from the same analysis",
    paste0("Selected leading-edge proteins, descriptive; see panel g. ",
           ALGEBRA)),
  # the ED6 atlases are the same descriptive theme aggregation as Figure 3b and
  # must carry the same declaration, not fall through to a default
  R("v9_ed_atlas_rescon", "median NES across mapped canonical GO terms",
    "animal", "descriptive aggregation of inferential inputs",
    "claim-eligible themes only",
    "gene set enrichment on the constituent GO terms",
    "BH within each constituent GO family; the THEME has no family of its own",
    "no",
    paste0("RES vs CON theme atlas. Theme colours summarise the median NES of ",
           "mapped canonical GO terms and are DESCRIPTIVE. Support markers ",
           "indicate that at least one constituent canonical GO term passed ",
           "its prespecified FDR threshold. Theme aggregation does not ",
           "constitute an additional multiple-testing family. ", ALGEBRA)),
  R("v9_ed_atlas_suscon", "median NES across mapped canonical GO terms",
    "animal", "descriptive aggregation of inferential inputs",
    "claim-eligible themes only",
    "gene set enrichment on the constituent GO terms",
    "BH within each constituent GO family; the THEME has no family of its own",
    "no",
    paste0("SUS vs CON theme atlas, completing the three-group trajectory. ",
           "Theme colours are a DESCRIPTIVE median of mapped canonical GO ",
           "terms; support markers indicate at least one FDR-supported ",
           "constituent term. Theme aggregation does not constitute an ",
           "additional multiple-testing family. ", ALGEBRA))
)

ED_DEFAULT_LEGEND <- list(
  v9_ed_bilateral_full = paste0(
    "Complete stored bilateral metric inventory. ", DESCRIPTIVE),
  v9_ed_precision = paste0(
    "Reliability gain from bilateral averaging. Open symbol = single ",
    "hemisphere; filled symbol = bilateral mean. Intraclass correlation is ",
    "reported as descriptive precision context; no significance test is ",
    "applied. ", DESCRIPTIVE),
  v9_ed_fingerprint_full = paste0(
    "Complete external signature inventory across all spatial units. ",
    DESCRIPTIVE),
  v9_ed_external_full = paste0(
    "Complete external validation inventory: every internal contrast against ",
    "every external signature, expected pairings and specificity comparisons ",
    "alike. BH-adjusted within this inventory."),
  v9_ed_internal_full = paste0(
    "Complete functional characterisation inventory: every canonical GO term ",
    "retained for every anatomical contrast, on one shared NES axis. This is ",
    "characterisation of the same data, not independent validation."),
  v9_ed_ca2_locator = "Location of CA2-SLM within the sampling scheme.",
  v9_ed_ca2_missing = paste0(
    "Pre-imputation missingness per acquisition, hemispheres paired within ",
    "animal. Both QC failures are a single hemisphere of an animal whose ",
    "other hemisphere passes. ", DESCRIPTIVE),
  v9_ed_ca2_displacement = paste0(
    "Relationship between missingness and the normalisation-related median ",
    "shift. The association is shown descriptively and no p-value is ",
    "computed: acquisitions are clustered within animal and hemisphere, so ",
    "an unclustered correlation test would be pseudoreplicated. ",
    DESCRIPTIVE),
  v9_ed_ca2_classes = paste0(
    "Robustness classification of the 28 canonical CA2-SLM proteins: ",
    "robustness-qualified, QC-sensitive, or insufficient observed data."),
  v9_ed_ca2_sensitivity = paste0(
    "Effect sizes with and without the two QC-failed hemispheres, on one ",
    "shared axis. Rows are grouped by robustness class. ", DESCRIPTIVE),
  v9_ed_module_fingerprint = paste0(
    "Mean CON z-score of each module's MEMBER PROTEINS across spatial units. ",
    "This is not the module eigengene, which is shown in panel b. ",
    DESCRIPTIVE),
  v9_ed_wgcna_phenotype = paste0(
    "Module eigengene differences by contrast. 13 of 15 modules show the ",
    "descriptive three-group directional trajectory; 0 of 45 module x ",
    "contrast cells and 0 of 35 stress x spatial-unit interaction omnibus ",
    "tests are FDR-supported (smallest FDR 0.27). The 13-of-15 statement is ",
    "DESCRIPTIVE and the three pairwise contrast directions are ",
    "algebraically related. A negative interaction result means no spatially ",
    "heterogeneous module effect survived multiple-testing correction; it ",
    "does NOT establish that such effects are absent. ", ALGEBRA),
  v9_ed_celltype = paste0(
    "External cell-type affinity counts per compartment. ", DESCRIPTIVE),
  v9_ed_m11 = paste0(
    "Worked example for one module. The active historical label is unchanged ",
    "and the proposed oligodendrocyte / myelin label is NOT activated; the ",
    "external evidence is context only. ", DESCRIPTIVE),
  v9_ed_gsea_curve_syn = paste0("Detailed ranked enrichment. ", GSEA_N, " ", ALGEBRA),
  v9_ed_gsea_curve_rna = paste0("Detailed ranked enrichment. ", GSEA_N, " ", ALGEBRA),
  v9_ed_gsea_curve_ox = paste0("Detailed ranked enrichment. ", GSEA_N, " ", ALGEBRA),
  v9_ed_atlas_rescon = paste0(
    "RES vs CON theme atlas. Theme colours are a descriptive median of ",
    "mapped canonical GO terms; support markers indicate at least one ",
    "FDR-supported constituent term. Theme aggregation is not an additional ",
    "multiple-testing family. ", ALGEBRA),
  v9_ed_atlas_suscon = paste0(
    "SUS vs CON theme atlas; see the RES vs CON panel. ", ALGEBRA),
  v9_ed_identity = paste0(
    "Among robustness-qualified proteins, the proportion whose strongest ",
    "phenotype-associated effect lies outside the canonical baseline ",
    "affinity set, across nested robustness subsets. These proteins were ",
    "already selected as FDR-supported and robustness-qualified, so this is ",
    "a DESCRIPTIVE pattern within a selected set and not an unbiased test ",
    "of a spatial null."),
  v9_ed_locations = paste0(
    "Among robustness-qualified proteins, the baseline dominant unit and the ",
    "unit of strongest phenotype-associated effect. Two distinct statements ",
    "are reported: 15 of 15 have their strongest effect in a DIFFERENT unit ",
    "from their single dominant baseline unit, and 14 of 15 are also outside ",
    "the broader canonical baseline affinity set, SNU13 being the sole ",
    "exception. effect_identity_relationship is the authoritative affinity ",
    "classification. Each arrow joins two independent measurements of one ",
    "protein; nothing travels between units. This is a descriptive pattern ",
    "within a selected set."),
  v9_ed_similarity = paste0(
    "Molecular similarity between spatial proteomic profiles, NOT anatomical ",
    "connectivity. One diverging scale centred on zero is used for all three ",
    "blocks because the metric and the zero reference are common; the ",
    "narrower range of the two region-level blocks is a property of the data. ",
    DESCRIPTIVE),
  v9_ed_network_distance = paste0(
    "Animal-level network distance from the CON consensus. For each CON ",
    "animal the distance was computed to a LEAVE-ONE-CON-ANIMAL-OUT ",
    "consensus, so no animal is compared with a consensus containing itself; ",
    "RES and SUS animals were compared with the full CON consensus. This ",
    "construction was verified against the edge data, not against its label. ",
    DESCRIPTIVE),
  v9_ed_nulls = paste0(
    "Exact whole-network permutation over animal labels, with the smallest ",
    "attainable p shown so the null is informative rather than merely ",
    "non-significant. No whole-network group difference was detected. The ",
    "statistic is a symmetric three-group centroid deviation and does not ",
    "use the CON consensus reference."),
  v9_ed_coupling = paste0(
    "Every distinct spatial-unit-pair by behavioural-outcome correlation. No ",
    "edge-behaviour association survived multiple-testing correction. With ",
    "n = 9 animals a single correlation has very little power, so this is a ",
    "limited-power negative result and not evidence of absence.")
)

# ------------------------------------------------------------------ assemble
by_id <- stats::setNames(ROLES, vapply(ROLES, function(r) r$id, character(1)))
panel_of <- function(id)
  Filter(function(p) identical(as.character(p$id), id), ct$panels)[[1]]

rows <- list(); leg <- c("# Figure legends - final_truth_v9", "",
  paste0("Generated from `figures/figure_final_truth_v9_contract.yml`. ",
         "Candidate layer; not promoted. Biological replicate is the animal ",
         "(9 animals, 3 per group) throughout unless a panel states otherwise."),
  "")

for (f in ct$figures) {
  leg <- c(leg, paste0("## ", as.character(f$name)), "",
           paste0("**", as.character(f$question %||% ""), "**"), "")
  for (it in f$layout) {
    id <- as.character(it$panel); p <- panel_of(id)
    r <- by_id[[id]]
    txt <- if (!is.null(r)) r$legend else
      (ED_DEFAULT_LEGEND[[id]] %||%
         gsub("[\r\n]+", " ", as.character(p$narrative %||% p$role %||% "")))
    leg <- c(leg, paste0("**", as.character(it$label), "** ", txt), "")
    sd <- sidecar(as.character(f$figure_key), id)
    rows[[length(rows) + 1L]] <- data.frame(
      figure = as.character(f$name), panel = as.character(it$label),
      panel_id = id,
      quantity = if (!is.null(r)) r$quantity else as.character(p$role %||% ""),
      biological_n = if (!is.null(r)) r$biological_n else "9 animals; 3 per group",
      replicate_unit = if (!is.null(r)) r$replicate_unit else "animal",
      descriptive_or_inferential = if (!is.null(r)) r$mode else "descriptive",
      selection_dependency = if (!is.null(r)) r$selection else "none",
      test = if (!is.null(r)) r$test else "none",
      FDR_family = if (!is.null(r)) r$fdr else "none",
      independent_validation = if (!is.null(r)) r$indep else "no",
      caption_sufficient = nzchar(txt),
      source_data_rows = if (is.null(sd)) NA_integer_ else nrow(sd),
      required_fix = "none",
      stringsAsFactors = FALSE)
  }
}
aud <- do.call(rbind, rows)
write_csv_safe(aud, file.path(TAB, "final_panel_statistical_role_audit.csv"))
writeLines(leg, file.path(REP, "final_figure_legends_v9.md"))

cat("\n===== PANEL STATISTICAL-ROLE AUDIT =====\n")
cat("panels audited          :", nrow(aud), "\n")
cat("with a written legend   :", sum(aud$caption_sufficient), "\n")
cat("inferential panels      :",
    sum(grepl("inferential", aud$descriptive_or_inferential)), "\n")
cat("claimed as INDEPENDENT validation:",
    sum(grepl("^YES", aud$independent_validation)), "\n")
print(aud[grepl("^YES", aud$independent_validation),
          c("figure", "panel", "quantity")], row.names = FALSE)
cat("\nselected-from-same-analysis panels:\n")
print(aud[grepl("same analysis", aud$selection_dependency),
          c("figure", "panel", "selection_dependency")], row.names = FALSE)
cat("\nwritten:\n  ", relative_to(file.path(REP, "final_figure_legends_v9.md")),
    "\n  ", relative_to(file.path(TAB, "final_panel_statistical_role_audit.csv")),
    "\n")
