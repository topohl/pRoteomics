#!/usr/bin/env Rscript
# Script: analysis/integration/test_theme_module_composition.R
# Stage: integration
# Consumes: required results/integration/test_enrichment_module_concordance/global/tables/ontology_aware_gsea_theme_assignments_all_contrasts.csv; results/tables/10_biological_integration/gsea_wgcna_concordance/global/ontology_aware_gsea_theme_assignments_all_contrasts.csv; config/manuscript_go_theme_registry.tsv; results/wgcna/build_wgcna_modules/<dataset>/tables/modules/WGCNA_feature_universe.csv; results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_feature_universe.csv; results/wgcna/build_wgcna_modules/<dataset>/tables/modules/WGCNA_modules_long.csv; results/tables/06_modules_WGCNA/01_WGCNA/<dataset>/modules/WGCNA_modules_long.csv; results/wgcna/build_module_identity_contract/<dataset>/tables/WGCNA_module_supermodule_membership_contract.csv; results/tables/06_modules_WGCNA/identity_contract/<dataset>/WGCNA_module_supermodule_membership_contract.csv; results/differential_abundance/run_clusterprofiler_enrichment/<dataset>/models/clusterProfiler_manifest.csv; data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/clusterProfiler_manifest.csv; optional config/gsea_wgcna_theme_module_mapping.csv; results/spatial_validation/build_module_spatial_atlas/global/tables/WGCNA_module_spatial_cell_affinity.csv; results/tables/11_spatial_systems/atlas/WGCNA_module_spatial_cell_affinity.csv
# Produces: results/integration/test_theme_module_composition/global/tables/theme_entity_composition_long.csv; results/integration/test_theme_module_composition/global/tables/theme_composition_summary.csv; results/integration/test_theme_module_composition/global/tables/theme_partition_agreement.csv; +6 more
# Scope: global
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Phenotype-blind composition check: do WGCNA modules and multi-module supermodules carry the curated GSEA themes? Structural evidence only; no theme-level phenotype test.
#
# Question (Tobias, 2026-10-08): do the WGCNA modules and supermodules fit or support the GSEA theme grouping?
#
# What this does, per compartment (microglia-enriched ROI, neuropil, soma):
#   theme protein sets  - proteins annotated (org.Mm.eg.db GOALL) to the theme's claim-eligible GO terms
#                         tested in that dataset, mapped through the dataset's GSEA-ranked genes to
#                         ProteinGroupIDs. Phenotype-independent, as the atlas membership rule requires.
#   entities            - WGCNA modules (identity contract) and the multi-module supermodules.
#   universe            - proteins in WGCNA AND the GSEA ranking AND annotated to a tested GO term.
#   statistic           - one-sided Fisher per theme x entity, conditional-MLE odds ratio with its 95% CI,
#                         fold enrichment, shares; BH (BY sensitivity) within dataset x level x variant.
#   support             - BH-FDR <= 0.05 and odds ratio >= 2.
#   sensitivities       - single-theme proteins, FDR<0.05 leading edges (any contrast and SUS - RES only),
#                         the full WGCNA universe, and core module members (|kME| >= 0.6).
#   partition agreement - Cramer's V and NMI between the theme and WGCNA partitions of the single-theme
#                         proteins, against 1000 label permutations.
#   curated mapping     - config/gsea_wgcna_theme_module_mapping.csv classified against the result (overlay).
#   wheel table         - one row per module with its supermodule, size, spatial structure and the themes it
#                         carries, for the module wheel figure.
#
# Downstream only: no WGCNA, GSEA or group-effect model is refitted; no source value is changed.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "wgcna_paths.R"))
source(repo_path("R", "enrichment_io.R"))
source(repo_path("R", "spatial_systems_paths.R"))
source(repo_path("R", "theme_module_composition_utils.R"))

ANALYSIS_ID <- "test_theme_module_composition"
SCRIPT_ID <- "analysis/integration/test_theme_module_composition.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
run <- integration_cli(default_dataset = "all", allow_all = TRUE)
datasets <- integration_datasets(run$dataset)
paths <- integration_dirs(ANALYSIS_ID, "global", create = !isTRUE(run$dry_run))

required_packages <- c("data.table", "AnnotationDbi", "org.Mm.eg.db", "ggplot2", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)

PRIMARY_SET <- "go_membership"
PRIMARY_UNIVERSE <- "wgcna_ranked_annotated"
PRIMARY_ENTITIES <- "all_members"

theme_table_file <- integration_find(
  "ontology_aware_gsea_theme_assignments_all_contrasts.csv",
  owner = "test_enrichment_module_concordance",
  legacy_stage = "10_biological_integration", legacy_substep = "gsea_wgcna_concordance")
registry_file <- repo_path("config", "manuscript_go_theme_registry.tsv")
mapping_file <- repo_path("config", "gsea_wgcna_theme_module_mapping.csv")
spatial_file <- spatial_systems_find("WGCNA_module_spatial_cell_affinity.csv", "atlas")
dataset_inputs <- lapply(stats::setNames(datasets, datasets), function(ds) c(
  wgcna_universe = wgcna_modules_artifact("WGCNA_feature_universe.csv", ds, child = "tables", "modules"),
  wgcna_modules_long = wgcna_modules_artifact("WGCNA_modules_long.csv", ds, child = "tables", "modules"),
  membership_contract = wgcna_identity_contract_artifact("WGCNA_module_supermodule_membership_contract.csv", ds),
  clusterprofiler_manifest = canonical_clusterprofiler_manifest_path(ds)))
required_inputs <- c(theme_assignments = theme_table_file, theme_registry = registry_file,
                     unlist(lapply(names(dataset_inputs), function(ds)
                       stats::setNames(dataset_inputs[[ds]], paste0(ds, "_", names(dataset_inputs[[ds]]))))))
optional_inputs <- c(theme_entity_mapping = mapping_file, module_spatial_atlas = spatial_file)

output_files <- c("theme_entity_composition_long.csv", "theme_composition_summary.csv", "theme_partition_agreement.csv",
                  "curated_mapping_composition_audit.csv", "wgcna_module_wheel_structure.csv",
                  "theme_protein_sets.csv", "composition_universe_audit.csv")

if (isTRUE(run$dry_run)) {
  dry_run_inputs(SCRIPT_ID, c(required_inputs, optional_inputs))
  for (f in output_files) dry_run_line("Output table", file.path(paths$tables, f))
  dry_run_line("Diagnostic figure", file.path(paths$plots, "theme_module_composition.svg"))
  dry_run_line("Interpretation README", file.path(paths$reports, "README.md"))
  quit(save = "no", status = 0L)
}

missing_required <- required_inputs[!file.exists(required_inputs)]
if (length(missing_required)) {
  stop("Required input(s) missing: ", paste(names(missing_required), missing_required, sep = " = ", collapse = "; "), call. = FALSE)
}
suppressPackageStartupMessages({ library(data.table) })

# ---------------------------------------------------------------- themes
registry <- data.table::fread(registry_file, showProgress = FALSE)
primary <- unique(registry[theme_role == "primary", .(theme_id, display_label, display_order, registry_version)])
data.table::setorder(primary, display_order)
THEMES <- primary$theme_id
if (length(THEMES) != 7L) stop("Expected 7 primary themes in the registry, found ", length(THEMES), call. = FALSE)

themes_dt <- data.table::fread(theme_table_file, showProgress = FALSE,
  select = c("dataset", "contrast", "GO_ID", "GSEA_FDR", "leading_edge_proteins", "theme_id", "theme_role",
             "theme_claim_eligible", "registry_version"))
if (!identical(unique(themes_dt$registry_version), unique(primary$registry_version))) {
  stop("Theme table registry version (", paste(unique(themes_dt$registry_version), collapse = ","),
       ") does not match config/manuscript_go_theme_registry.tsv (", unique(primary$registry_version), ").", call. = FALSE)
}
eligible <- themes_dt[theme_claim_eligible %in% TRUE & theme_id %in% THEMES]
tested_go <- unique(themes_dt[, .(dataset, GO_ID)])

go_members <- local({
  keys <- unique(tested_go$GO_ID)
  ann <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = keys, keytype = "GOALL", columns = "SYMBOL"))
  unique(data.table::as.data.table(ann)[!is.na(SYMBOL), .(GO_ID = GOALL, SYMBOL)])
})

# ---------------------------------------------------------------- per dataset inputs
ranked_universe <- function(ds) {
  man <- read_canonical_clusterprofiler_manifest(dataset_inputs[[ds]][["clusterprofiler_manifest"]], ds)
  man <- man[man$result_type == "GSEA_GO" & man$ontology == "BP", , drop = FALSE]
  if (!nrow(man)) stop("No GSEA_GO BP rows in the clusterProfiler manifest for ", ds, call. = FALSE)
  files <- unique(man$collapsed_gene_input_file)   # the resolved (rebased or staged) runtime path
  rows <- data.table::rbindlist(lapply(files, function(f)
    data.table::fread(f, select = c("GeneSymbol", "contributing_ProteinGroupIDs"), showProgress = FALSE)))
  pg <- unique(rows[, .(ProteinGroupID = unlist(strsplit(contributing_ProteinGroupIDs, ";", fixed = TRUE))), by = GeneSymbol])
  list(map = pg, n_genes_per_comparison = range(as.integer(man$n_genes), na.rm = TRUE), n_comparisons = length(files))
}

comp_all <- list(); sets_long <- list(); universe_audit <- list(); partition <- list()
membership_all <- list(); sizes_all <- list()
for (ds in datasets) {
  inp <- dataset_inputs[[ds]]
  fu <- data.table::fread(inp[["wgcna_universe"]], select = c("ProteinGroupID", "included_in_wgcna", "ModuleID"), showProgress = FALSE)
  fu <- fu[included_in_wgcna %in% TRUE]
  ml <- data.table::fread(inp[["wgcna_modules_long"]], select = c("ModuleID", "ProteinGroupID", "is_core_kME_0.6"), showProgress = FALSE)
  contract <- data.table::fread(inp[["membership_contract"]], showProgress = FALSE)
  if (!all(contract$membership_status == "current_authoritative_membership")) {
    stop("Membership contract for ", ds, " is not current_authoritative_membership.", call. = FALSE)
  }
  if (!setequal(unique(fu$ModuleID), contract$module_id)) {
    stop("Feature-universe modules and the membership contract disagree for ", ds, ".", call. = FALSE)
  }
  rk <- ranked_universe(ds)
  U_wgcna <- unique(fu$ProteinGroupID)
  U_ranked <- unique(rk$map$ProteinGroupID)
  annotated_genes <- unique(go_members[GO_ID %in% tested_go[dataset == ds]$GO_ID]$SYMBOL)
  U_annot <- intersect(intersect(U_wgcna, U_ranked), rk$map[GeneSymbol %in% annotated_genes]$ProteinGroupID)

  # theme sets
  genes_of <- function(th) unique(go_members[GO_ID %in% unique(eligible[dataset == ds & theme_id == th]$GO_ID)]$SYMBOL)
  le_of <- function(th, contrast = NULL) {
    x <- eligible[dataset == ds & theme_id == th & GSEA_FDR < 0.05]
    if (!is.null(contrast)) x <- x[x$contrast == contrast]
    intersect(tmc_split_ids(x$leading_edge_proteins), U_wgcna)
  }
  go_sets <- lapply(stats::setNames(THEMES, THEMES), function(th) intersect(rk$map[GeneSymbol %in% genes_of(th)]$ProteinGroupID, U_wgcna))
  set_variants <- list(
    go_membership = go_sets,
    go_membership_single_theme = tmc_single_theme_sets(go_sets),
    leading_edge_fdr = lapply(stats::setNames(THEMES, THEMES), le_of),
    leading_edge_fdr_sus_res = lapply(stats::setNames(THEMES, THEMES), le_of, contrast = "SUS - RES"))
  for (v in names(set_variants)) for (th in THEMES) {
    if (!length(set_variants[[v]][[th]])) next          # an empty set is recorded by n_theme = 0 in the composition table
    sets_long[[length(sets_long) + 1L]] <- data.table::data.table(dataset = ds, theme_id = th, theme_set_definition = v,
      n_terms = length(unique(eligible[dataset == ds & theme_id == th]$GO_ID)), ProteinGroupID = set_variants[[v]][[th]])
  }

  # entities
  mods_all <- split(fu$ProteinGroupID, fu$ModuleID)
  core_pg <- unique(ml[`is_core_kME_0.6` %in% TRUE]$ProteinGroupID)
  mods_core <- lapply(mods_all, function(p) intersect(p, core_pg))
  multi <- contract[, .(n = .N, members = paste(sort(module_id), collapse = ";")), by = supermodule_id][n > 1]
  sm_of <- function(mods) stats::setNames(lapply(multi$members, function(mm) unique(unlist(mods[strsplit(mm, ";", fixed = TRUE)[[1]]]))),
                                          multi$supermodule_id)
  members_of_sm <- stats::setNames(multi$members, multi$supermodule_id)

  runs <- list(
    list(set = "go_membership", universe = "wgcna_ranked_annotated", entities = "all_members"),
    list(set = "go_membership_single_theme", universe = "wgcna_ranked_annotated", entities = "all_members"),
    list(set = "leading_edge_fdr", universe = "wgcna_ranked_annotated", entities = "all_members"),
    list(set = "leading_edge_fdr_sus_res", universe = "wgcna_ranked_annotated", entities = "all_members"),
    list(set = "go_membership", universe = "wgcna_all", entities = "all_members"),
    list(set = "go_membership", universe = "wgcna_ranked_annotated", entities = "core_kME_0.6"))
  for (r in runs) {
    U <- if (r$universe == "wgcna_all") U_wgcna else U_annot
    mods <- if (r$entities == "core_kME_0.6") mods_core else mods_all
    comp_all[[length(comp_all) + 1L]] <- tmc_composition(set_variants[[r$set]], mods, U, ds, "module", r$set, r$universe,
                                                         entity_definition = r$entities)
    if (nrow(multi)) comp_all[[length(comp_all) + 1L]] <- tmc_composition(set_variants[[r$set]], sm_of(mods), U, ds, "supermodule",
                                                                          r$set, r$universe, entity_members = members_of_sm,
                                                                          entity_definition = r$entities)
  }

  # partition agreement on single-theme proteins in the primary universe
  st <- data.table::rbindlist(lapply(THEMES, function(th) data.table::data.table(ProteinGroupID = intersect(set_variants$go_membership_single_theme[[th]], U_annot), theme_id = th)))
  st <- merge(st, fu[, .(ProteinGroupID, module_id = ModuleID)], by = "ProteinGroupID")
  st <- merge(st, contract[, .(module_id, supermodule_id)], by = "module_id")
  for (lev in c("module", "supermodule")) {
    lab <- if (lev == "module") st$module_id else st$supermodule_id
    partition[[length(partition) + 1L]] <- cbind(data.frame(dataset = ds, level = lev, protein_set = "single_theme_go_membership",
                                                            universe_definition = PRIMARY_UNIVERSE),
                                                 tmc_partition_agreement(st$theme_id, lab))
  }

  membership_all[[ds]] <- data.frame(dataset = ds, module_id = contract$module_id, supermodule_id = contract$supermodule_id,
                                     stringsAsFactors = FALSE)
  sizes_all[[ds]] <- data.frame(dataset = ds, module_id = names(mods_all), n_proteins = lengths(mods_all),
                                n_core_kME_0.6 = lengths(mods_core)[names(mods_all)],
                                n_in_primary_universe = vapply(mods_all, function(p) length(intersect(p, U_annot)), integer(1)),
                                stringsAsFactors = FALSE)
  universe_audit[[ds]] <- data.frame(dataset = ds, n_wgcna_included = length(U_wgcna), n_ranked_protein_groups = length(U_ranked),
    n_ranked_genes = length(unique(rk$map$GeneSymbol)), n_gsea_go_bp_comparisons = rk$n_comparisons,
    n_genes_per_comparison_min = rk$n_genes_per_comparison[1], n_genes_per_comparison_max = rk$n_genes_per_comparison[2],
    n_wgcna_and_ranked = length(intersect(U_wgcna, U_ranked)), n_primary_universe = length(U_annot),
    n_tested_go_terms = nrow(tested_go[dataset == ds]), n_modules = length(mods_all), n_multi_module_supermodules = nrow(multi),
    theme_registry_version = unique(primary$registry_version),
    org_Mm_eg_db_version = as.character(utils::packageVersion("org.Mm.eg.db")), stringsAsFactors = FALSE)
}

comp <- do.call(rbind, comp_all)
is_primary <- comp$theme_set_definition == PRIMARY_SET & comp$universe_definition == PRIMARY_UNIVERSE &
  comp$entity_definition == PRIMARY_ENTITIES
comp$is_primary_variant <- is_primary
comp_primary <- comp[is_primary, , drop = FALSE]
summary_all <- tmc_coherence_summary(comp)
summary_primary <- summary_all[summary_all$theme_set_definition == PRIMARY_SET & summary_all$universe_definition == PRIMARY_UNIVERSE &
                                 summary_all$entity_definition == PRIMARY_ENTITIES, , drop = FALSE]
summary_primary <- tmc_carrier_robustness(summary_primary, comp, PRIMARY_SET, PRIMARY_UNIVERSE, PRIMARY_ENTITIES)
n_ds_with_carrier <- tapply(summary_primary$n_carriers > 0, summary_primary$theme_id, sum)
summary_primary$n_compartments_with_carrier <- as.integer(n_ds_with_carrier[summary_primary$theme_id])
summary_primary$theme_label <- primary$display_label[match(summary_primary$theme_id, primary$theme_id)]

mapping_audit <- if (file.exists(mapping_file)) {
  tmc_mapping_audit(utils::read.csv(mapping_file, stringsAsFactors = FALSE, colClasses = "character"),
                    comp_primary[comp_primary$level %in% c("module", "supermodule"), , drop = FALSE])
} else NULL

spatial <- if (file.exists(spatial_file)) {
  s <- data.table::fread(spatial_file, showProgress = FALSE)
  as.data.frame(s[, .(dataset, module_id = ModuleID, canonical_display_label, peak_unit, spatial_tau,
                      bilateral_reproducibility_class)])
} else data.frame(dataset = character(), module_id = character())
wheel <- tmc_wheel_structure(do.call(rbind, membership_all), do.call(rbind, sizes_all), spatial, comp_primary, THEMES)

# ---------------------------------------------------------------- write
tables <- list(
  theme_entity_composition_long.csv = comp,
  theme_composition_summary.csv = summary_primary,
  theme_partition_agreement.csv = do.call(rbind, partition),
  curated_mapping_composition_audit.csv = mapping_audit %||% data.frame(),
  wgcna_module_wheel_structure.csv = wheel,
  theme_protein_sets.csv = as.data.frame(data.table::rbindlist(sets_long)),
  composition_universe_audit.csv = do.call(rbind, universe_audit))
for (nm in names(tables)) invisible(write_integration_table(tables[[nm]], paths, nm))
write_csv_safe(summary_all, file.path(paths$tables, "theme_composition_summary_all_variants.csv"))

# diagnostic figure (not a manuscript panel)
dir_create(paths$plots)
pd <- comp_primary
pd$theme <- factor(primary$display_label[match(pd$theme_id, THEMES)], levels = rev(primary$display_label))
pd$entity <- ifelse(pd$level == "module", sub("^WGCNA_", "", pd$entity_id), paste0(pd$entity_id, " (", gsub("WGCNA_", "", pd$member_modules), ")"))
pd$facet <- paste(pd$dataset, pd$level)
pd$log2_fold <- pmin(pmax(log2(pd$fold_enrichment), 0), 3)
g <- ggplot2::ggplot(pd, ggplot2::aes(entity, theme)) +
  ggplot2::geom_point(ggplot2::aes(size = frac_theme_in_entity, fill = log2_fold, colour = composition_supported), shape = 21) +
  ggplot2::scale_size_area(max_size = 5, limits = c(0, 1), name = "share of theme") +
  ggplot2::scale_fill_gradient(low = "white", high = "#1F3D52", limits = c(0, 3), name = "log2 fold") +
  ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey80"), name = "BH<=0.05 & OR>=2") +
  ggplot2::facet_grid(. ~ facet, scales = "free_x", space = "free_x") +
  ggplot2::labs(x = NULL, y = NULL, title = "Theme x WGCNA composition (primary variant; structural, not phenotype evidence)") +
  ggplot2::theme_minimal(base_size = 7) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
svg_path <- file.path(paths$plots, "theme_module_composition.svg")
ggplot2::ggsave(svg_path, g, width = 320, height = 90, units = "mm", device = svglite::svglite)
ggplot2::ggsave(sub("[.]svg$", ".pdf", svg_path), g, width = 320, height = 90, units = "mm", device = grDevices::cairo_pdf)

# README
sp <- summary_primary
readme <- c(
  "# Theme x WGCNA module composition",
  "",
  paste0("Generated: ", format(Sys.time(), tz = "Europe/Berlin", usetz = TRUE)),
  "",
  "## Question",
  "",
  "Do the phenotype-blind WGCNA co-abundance modules and multi-module supermodules group proteins the way the curated GSEA themes (manuscript_go_themes registry) do?",
  "",
  "## Evidence boundary",
  "",
  "Structural composition only. Theme sets come from GO annotation of the theme's claim-eligible terms; modules come from co-abundance across all samples with group labels unused. A supported cell says WGCNA carries the theme's proteins together. It is not a group effect, and no theme-level phenotype p-value or FDR is created. The GO annotation that defines the themes also informed the module labels and the curated theme-module mapping, so the curated mapping is only an overlay here; the primary analysis tests every theme against every module.",
  "",
  "## Method",
  "",
  "Primary variant: theme set = proteins annotated (org.Mm.eg.db GOALL) to the theme's claim-eligible GO terms tested in the dataset, mapped through the GSEA-ranked genes to ProteinGroupIDs; universe = proteins in WGCNA and the GSEA ranking with a GO annotation to any tested term; entities = WGCNA modules (identity contract) and multi-module supermodules. One-sided Fisher test; conditional-MLE odds ratio with exact 95% CI; BH within dataset x level x variant (BY reported). Support = BH-FDR <= 0.05 and odds ratio >= 2. Sensitivities: single-theme proteins, FDR<0.05 leading edges (any contrast; SUS - RES only), the full WGCNA universe, core members (|kME| >= 0.6). Classes by the share of the theme its supporting modules capture: coherent >= 50%, core 20-50%, minor < 20%, scattered = no supporting module. Partition agreement: Cramer's V and NMI between theme and WGCNA labels of single-theme proteins against 1000 permutations of the WGCNA labels.",
  "",
  "## Primary result (module level)",
  "",
  paste0("- ", sp$dataset, " | ", sp$theme_id, ": ", sp$composition_class, "; carriers ",
         ifelse(nzchar(sp$carriers), sp$carriers, "none"), "; share captured ", sprintf("%.0f%%", 100 * sp$share_captured_by_carriers),
         ifelse(nzchar(sp$shares_carrier_with), paste0("; shares a carrier with ", sp$shares_carrier_with), "")),
  "",
  "## Limitations",
  "",
  "Very large modules dilute fold enrichment and inflate the share captured; odds ratio, fold and share are reported together, with a core-member sensitivity. Themes overlap (shared GO terms and proteins) and modules partition the universe, so tests in one family are dependent (BY reported). The three compartments are separately built networks from the same animals, not independent replications. Microglia refers to a microglia-enriched ROI.",
  "")
dir_create(paths$reports)
writeLines(readme, file.path(paths$reports, "README.md"), useBytes = TRUE)

write_integration_manifest(
  paths, c(required_inputs, optional_inputs),
  list(tables = paths$tables, source_data = paths$source_data, plots = paths$plots,
       report = file.path(paths$reports, "README.md")),
  list(datasets = datasets, primary_theme_set = PRIMARY_SET, primary_universe = PRIMARY_UNIVERSE,
       primary_entities = PRIMARY_ENTITIES, support_rule = "BH-FDR <= 0.05 and odds ratio >= 2",
       family = "dataset x level x theme_set_definition x universe_definition x entity_definition",
       class_breaks = "coherent >= 0.5, core 0.2-0.5, minor < 0.2, scattered = no supported module",
       partition_permutations = 1000L, theme_registry_version = unique(primary$registry_version)),
  "Structural composition of curated GSEA themes in WGCNA modules and supermodules. Phenotype-blind; no theme-level phenotype test; no source value recomputed.")

message("Theme x WGCNA composition complete: ", paths$tables)
