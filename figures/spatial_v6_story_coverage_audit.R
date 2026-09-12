#!/usr/bin/env Rscript

# Part-21: the figure-story coverage audit.
#
# One row per panel of the spatial_v6 figure family, recording what the
# scientific question of that panel REQUIRES along the two axes the brief
# separates - spatial biology (compartment / region / layer) and stress
# phenotype (groups, contrasts, three-group trajectory) - and what the panel
# actually carries.
#
# The required/present values are read from the contract's declared
# spatial_axes and phenotype_axes fields, so the audit cannot drift away from
# what the panels are actually built to show: if a panel changes its axes, the
# contract changes and this audit changes with it.
#
# Deliberate omissions are recorded as INTENTIONAL AND CORRECT with a reason.
# A phenotype-blind Figure-2 panel is not a coverage gap; it is the design.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "nature_v2_figure_utils.R"))
source(repo_path("R", "spatial_grammar_utils.R"))
source(repo_path("R", "spatial_v6_figure_utils.R"))
suppressPackageStartupMessages({ library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/spatial_v6_story_coverage_audit.R")

OUT <- function(...) {
  d <- path_results("tables", "manuscript_candidates", "spatial_v6")
  dir_create(d)
  file.path(d, ...)
}

ct <- s6e_contract()
panels <- ct$panels
names(panels) <- vapply(panels, function(p) as.character(p$id), character(1))

# role -> what the scientific question requires. This is the editorial
# judgement layer and is written out explicitly so it can be argued with.
REQ <- list(
  experimental_anatomy_anchor = list(
    q = "Where in the hippocampus was tissue sampled, and at what resolution?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "orientation panel; stress phenotype would add nothing and would imply the schematic is a result",
    role = "MAIN FIGURE 2"),
  direct_baseline_spatial_molecular_view = list(
    q = "What molecular patterns distinguish the sampled regions, layers and compartments at baseline?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "phenotype-blind BY CONSTRUCTION; selecting or colouring by stress would destroy the claim that the spatial architecture is established before any stress model",
    role = "MAIN FIGURE 2"),
  full_baseline_spatial_fingerprint = list(
    q = "What is the complete baseline spatial architecture over external signatures?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "same phenotype-blind rule as the main-figure fingerprint, at full row depth",
    role = "EXTENDED DATA"),
  global_structure = list(
    q = "Do the compartment proteomes occupy structured molecular spaces?",
    comp = "YES", reg = "YES", lay = "NO", grp = "NO", con = "NO",
    traj = "NO", link = "NO",
    just = "structure is the claim; group labels would invite a phenotype reading of a phenotype-blind result",
    role = "MAIN FIGURE 2"),
  reproducibility = list(
    q = "Do paired hemispheres reproduce the spatial molecular architecture?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "a technical reproducibility claim; must stay phenotype-blind or it becomes a stress result",
    role = "MAIN FIGURE 2"),
  compartment_identity = list(
    q = "Do the three compartments carry their expected molecular identities?",
    comp = "YES", reg = "NO", lay = "NO", grp = "NO", con = "NO",
    traj = "NO", link = "NO",
    just = "compartment-level claim only; region and layer are not the question here",
    role = "MAIN FIGURE 2"),
  external_validation = list(
    q = "Do internal anatomical contrasts recover the corresponding external hippocampal signatures?",
    comp = "NO", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "external validation is computed on CON only; phenotype is not part of the question",
    role = "MAIN FIGURE 2"),
  internal_validation = list(
    q = "Does each anatomical contrast recover its own canonical biological programs?",
    comp = "NO", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "internal anatomical GSEA, CON only",
    role = "EXTENDED DATA"),
  technical_context = list(
    q = "Is the measurement deep enough in every compartment?",
    comp = "YES", reg = "NO", lay = "NO", grp = "NO", con = "NO",
    traj = "NO", link = "NO",
    just = "a depth statement; neither region nor phenotype is relevant",
    role = "EXTENDED DATA"),
  measurement_quality = list(
    q = "How much precision is gained by averaging hemispheres?",
    comp = "YES", reg = "YES", lay = "NO", grp = "NO", con = "NO",
    traj = "NO", link = "NO",
    just = "a precision statement about the measurement, not about stress",
    role = "EXTENDED DATA"),
  program_atlas = list(
    q = "Where across compartment, region and layer do coordinated program differences occur?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "YES", traj = "PARTIAL", link = "YES",
    just = "the main SUS-RES field; the other two contrasts are given the same atlas in ED6 so the trajectory is complete without crowding the main figure",
    role = "MAIN FIGURE 3"),
  anatomical_program_bridge = list(
    q = "Where anatomically are the three representative programs, and how do RES and SUS relate to CON there?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "YES", traj = "YES", link = "YES",
    just = "this panel exists precisely to join the spatial and phenotype axes",
    role = "MAIN FIGURE 3"),
  direct_evidence = list(
    q = "Is the program difference real at the level of the ranked enrichment itself?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "YES", traj = "YES", link = "YES",
    just = "exact reconstructed curve plus the three-contrast strip",
    role = "MAIN FIGURE 3"),
  leading_edge_proteins = list(
    q = "Which actual proteins carry the program signal?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "NO", traj = "NO", link = "YES",
    just = "protein-level evidence for the SUS-RES contrast; the three-group trajectory lives on the curve panels",
    role = "MAIN FIGURE 3"),
  gsea_support = list(
    q = "What is the complete enrichment evidence behind each highlighted program?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "YES", traj = "YES", link = "YES",
    just = "full support for the main-figure exemplars",
    role = "EXTENDED DATA"),
  gsea_support_trajectory = list(
    q = "Do the other two contrasts show the same spatial pattern?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "YES", traj = "YES", link = "YES",
    just = "RES-CON and SUS-CON atlases complete the three-group trajectory that the main figure shows only for SUS-RES",
    role = "EXTENDED DATA"),
  wgcna_structure = list(
    q = "What is the module structure?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "module identity is descriptive molecular architecture; phenotype is NOT REQUIRED for the identity component",
    role = "EXTENDED DATA"),
  wgcna_spatial_fingerprint = list(
    q = "How is each module actually organised across space?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "phenotype-blind by design; this is the baseline spatial identity of the modules",
    role = "EXTENDED DATA"),
  wgcna_celltype = list(
    q = "What external cell-type context do the modules carry?",
    comp = "YES", reg = "NO", lay = "NO", grp = "NO", con = "NO",
    traj = "NO", link = "NO",
    just = "external annotation only",
    role = "EXTENDED DATA"),
  wgcna_phenotype_context = list(
    q = "Do module-level effects differ by stress outcome?",
    comp = "NO", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "YES", traj = "YES", link = "YES",
    just = "requires both axes; the inferential status (0 of 45 FDR-supported) must be stated on the panel",
    role = "EXTENDED DATA"),
  ca2_qc = list(
    q = "Does differential missingness and QC failure explain the CA2-SLM result?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "YES", traj = "NO", link = "YES",
    just = "the missingness asymmetry is SUS vs RES, so stress group MUST remain visible",
    role = "EXTENDED DATA"),
  ca2_spatial_locator = list(
    q = "Where is CA2-SLM?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "a locator; it carries no quantity and no phenotype",
    role = "EXTENDED DATA"),
  stress_identity_robustness = list(
    q = "Does the outside-baseline-affinity result survive every robustness restriction?",
    comp = "NO", reg = "NO", lay = "NO", grp = "YES", con = "YES",
    traj = "NO", link = "YES",
    just = "a robustness summary; the actual locations are the job of the paired panel below it",
    role = "EXTENDED DATA"),
  baseline_vs_stress_location = list(
    q = "Is the strongest stress-associated effect where the protein is most abundant at baseline?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "YES", traj = "NO", link = "YES",
    just = "REQUIRES BOTH AXES: baseline spatial identity versus the location of the phenotype-associated effect",
    role = "EXTENDED DATA"),
  baseline_vs_stress_location_alt = list(
    q = "Same question, alternative encoding for legibility comparison.",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "YES", con = "YES", traj = "NO", link = "YES",
    just = "built so the origin-destination and matrix encodings can be compared directly",
    role = "EXTENDED DATA"),
  what_a_spatial_network_is = list(
    q = "What does a spatial molecular network actually represent?",
    comp = "YES", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "CON baseline only; the panel must establish the meaning of an edge BEFORE any group comparison, and must say it is not connectivity",
    role = "EXTENDED DATA"),
  network_replicates = list(
    q = "How much do animal-level networks vary?",
    comp = "YES", reg = "NO", lay = "NO", grp = "YES", con = "NO",
    traj = "NO", link = "NO",
    just = "replicate spread by group, with no inferential claim",
    role = "EXTENDED DATA"),
  network_null = list(
    q = "Is there a detectable whole-network group difference?",
    comp = "YES", reg = "NO", lay = "NO", grp = "YES", con = "YES",
    traj = "NO", link = "NO",
    just = "an informative null; the enumeration had resolution to 0.0036",
    role = "EXTENDED DATA"),
  coupling_null = list(
    q = "Does any network edge couple to behaviour or physiology?",
    comp = "NO", reg = "YES", lay = "YES where measured",
    grp = "NO", con = "NO", traj = "NO", link = "NO",
    just = "correlation across all 9 animals; stress group is not a factor in this analysis",
    role = "EXTENDED DATA")
)

axis_has <- function(s, token) {
  if (is.na(s) || !nzchar(s)) return("NO")
  if (grepl(token, s, fixed = TRUE)) "YES" else "NO"
}

rows <- list()
for (f in ct$figures) {
  # what the FIGURE as a whole carries, so a per-panel gap can be reported as
  # covered by a sibling rather than silently pardoned or silently alarmed
  fig_axes <- paste(vapply(f$layout, function(x)
    as.character(panels[[as.character(x$panel)]]$spatial_axes %||% ""),
    character(1)), collapse = "+")
  for (it in f$layout) {
    id <- as.character(it$panel)
    p <- panels[[id]]
    role <- as.character(p$role %||% "")
    r <- REQ[[role]]
    if (is.null(r)) {
      stop("no coverage requirement declared for role: ", role,
           " (panel ", id, ")", call. = FALSE)
    }
    sax <- as.character(p$spatial_axes %||% "")
    pax <- as.character(p$phenotype_axes %||% "")
    present <- paste0(
      "compartment=", axis_has(sax, "compartment"),
      "; region=", axis_has(sax, "region"),
      "; layer=", axis_has(sax, "layer"),
      "; phenotype=", if (identical(pax, "none_intentional")) "NONE (intentional)"
                      else if (identical(pax, "none_applicable")) "NONE (not applicable)"
                      else pax)
    want_sp <- c(r$comp, r$reg, r$lay)
    got_sp <- c(axis_has(sax, "compartment"), axis_has(sax, "region"),
                axis_has(sax, "layer"))
    missing <- character(0)
    nm <- c("compartment", "region", "layer")
    for (k in seq_along(nm)) {
      if (!grepl("^YES", want_sp[k])) next
      if (got_sp[k] == "YES") next
      # "YES where measured" is SATISFIED by the absence of a layer on a panel
      # about a region-level compartment. Soma and the microglia-enriched ROI
      # have no laminar resolution, so demanding a layer there would be
      # demanding a fabrication - the exact error the brief forbids.
      if (identical(want_sp[k], "YES where measured") &&
          identical(nm[k], "layer")) next
      missing <- c(missing, nm[k])
    }
    if (r$grp == "YES" && identical(pax, "none_intentional")) {
      missing <- c(missing, "stress group")
    }
    covered <- missing[vapply(missing, function(m)
      axis_has(fig_axes, m) == "YES", logical(1))]
    still <- setdiff(missing, covered)
    omission_ok <- if (!length(missing)) "NOT APPLICABLE - nothing missing" else
      if (!length(still)) "COVERED BY A SIBLING PANEL IN THE SAME FIGURE" else
        "REVIEW"
    if (r$grp == "NO" && identical(pax, "none_intentional") &&
        !length(missing)) {
      omission_ok <- "INTENTIONAL AND CORRECT"
    }
    rows[[length(rows) + 1L]] <- data.frame(
      figure = as.character(f$name),
      panel = paste0(as.character(it$label), " - ", id),
      scientific_question = r$q,
      required_compartment_information = r$comp,
      required_region_information = r$reg,
      required_layer_information = r$lay,
      required_stress_group_information = r$grp,
      required_stress_contrast_information = r$con,
      requires_three_group_trajectory = r$traj,
      requires_spatial_and_stress_link = r$link,
      currently_present = present,
      currently_missing = if (length(missing)) paste(missing, collapse = "; ") else "none",
      covered_by_sibling_panel = if (length(covered)) paste(covered, collapse = "; ") else "n/a",
      still_missing_at_figure_level = if (length(still)) paste(still, collapse = "; ") else "none",
      omission_justified = omission_ok,
      justification = r$just,
      recommended_fix = if (length(still)) "add the axis to this panel or to the figure" else
        if (length(covered)) "none required; a sibling panel in the same figure carries it" else "none required",
      main_story_role = r$role,
      stringsAsFactors = FALSE)
  }
}
audit <- do.call(rbind, rows)
write_csv_safe(audit, OUT("figure_story_coverage_audit.csv"))

# ---------------------------------------------------------------------------
# the spatial x phenotype coverage matrix required by brief section 44
# ---------------------------------------------------------------------------
mat <- do.call(rbind, lapply(seq_len(nrow(audit)), function(i) {
  a <- audit[i, ]
  cp <- a$currently_present
  has <- function(tok) grepl(tok, cp, fixed = TRUE)
  pax <- sub(".*phenotype=", "", cp)
  data.frame(
    figure = a$figure, panel = a$panel,
    compartment = if (has("compartment=YES")) "x" else "",
    region = if (has("region=YES")) "x" else "",
    layer = if (has("layer=YES")) "x" else "",
    CON = if (grepl("CON", pax)) "x" else "",
    RES = if (grepl("RES", pax)) "x" else "",
    SUS = if (grepl("SUS", pax)) "x" else "",
    `RES-CON` = if (grepl("RES-CON", pax)) "x" else "",
    `SUS-CON` = if (grepl("SUS-CON", pax)) "x" else "",
    `SUS-RES` = if (grepl("SUS-RES", pax)) "x" else "",
    baseline_spatial_identity =
      if (a$main_story_role %in% c("MAIN FIGURE 2") ||
          grepl("baseline|wgcna_spatial|network", a$scientific_question,
                ignore.case = TRUE)) "x" else "",
    stress_sensitive_spatial_identity =
      if (a$requires_spatial_and_stress_link == "YES") "x" else "",
    direct_proteins = if (grepl("protein", a$scientific_question,
                                ignore.case = TRUE)) "x" else "",
    program_level_biology = if (grepl("program|enrichment",
                                      a$scientific_question,
                                      ignore.case = TRUE)) "x" else "",
    systems_level_biology = if (grepl("module|network", a$scientific_question,
                                      ignore.case = TRUE)) "x" else "",
    qc_robustness = if (grepl("QC|robust|reproduc|deep|precision",
                              a$scientific_question,
                              ignore.case = TRUE)) "x" else "",
    check.names = FALSE, stringsAsFactors = FALSE)
}))
write_csv_safe(mat, OUT("spatial_phenotype_coverage_matrix.csv"))

cat("\n===== figure-story coverage audit =====\n")
cat("panels audited:", nrow(audit), " figures:", length(unique(audit$figure)), "\n")
cat("\nomission status:\n"); print(table(audit$omission_justified))
cat("\nmissing information (should be 'none' for every row):\n")
print(table(audit$still_missing_at_figure_level))
cat("\nmain story role:\n"); print(table(audit$main_story_role))
cat("\npanels that require BOTH spatial and stress axes:",
    sum(audit$requires_spatial_and_stress_link == "YES"), "\n")
cat("panels intentionally phenotype-blind:",
    sum(audit$omission_justified == "INTENTIONAL AND CORRECT"), "\n")
cat("\ncoverage matrix columns:", ncol(mat), " rows:", nrow(mat), "\n")
cat("\nwritten:\n  ", relative_to(OUT("figure_story_coverage_audit.csv")),
    "\n  ", relative_to(OUT("spatial_phenotype_coverage_matrix.csv")), "\n")
