#!/usr/bin/env Rscript
# Export the tables the Figure 4 WGCNA candidate draws (added 2026-11-23, eight animals).
# Reads finished results only; fits nothing. Run from the repository root:
#   Rscript tools/export_figure_04_wgcna_candidate.R [new output folder; default exports/publication_source_data/figure_04_wgcna_candidate]
#
# modules.csv            every module: size, reviewed label and short label (config/wgcna_labels), bilateral
#                        reproducibility (the stability column), control spatial class, top carried family, and
#                        the bridge rows of panel a
# module_profiles.csv    control profile per module and unit (describe_module_control_geography)
# territory_plane.csv    neuropil perforant-path and mossy-fibre contrasts (describe_module_control_geography)
# figure2_links.csv      Figure 2 home terms whose best carrier module peaks at home (link_home_terms_to_modules)
# family_composition.csv module x Figure 3 family composition, primary variant (test_theme_module_composition)
# family_near_misses.csv per compartment and family without a bridge row, the modules that came closest
# exemplar_modules.csv   share of each Figure 3 exemplar's SUS-RES leading edge in each module
# input_manifest.csv     every input with its SHA-256
#
# BRIDGE ROWS (rule approved 2026-10-30, kept for eight animals 2026-11-23): every module that holds at least 20% of a
# Figure 3 family's proteins and carries it (composition test: BH FDR <= 0.05 and odds ratio >= 2), grouped by its
# largest such family; then every module that carries a Figure 2 home-term link and is not already a row.
#
# STABILITY: bilateral reproducibility (quantify_module_bilateral_identity): the median over the animals of the
# Pearson r between the left and right spatial profile; rows of class poor_bilateral_reproducibility are drawn faded.

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "pipeline.yml"))) stop("Run from the pRoteomics repository root.", call. = FALSE)
local({  # waits for the normalisation decision (config/normalisation.yml)
  source(file.path(root, "R", "paths.R"), local = TRUE)
  source(repo_path("R", "normalisation_contract.R"), local = TRUE)
  assert_normalisation_decided("tools/export_figure_04_wgcna_candidate.R")
})
rel <- function(...) file.path(root, ...)
args <- commandArgs(trailingOnly = TRUE)   # optional: another, new output folder (a trial export outside exports/)
out <- if (length(args)) args[[1]] else rel("exports", "publication_source_data", "figure_04_wgcna_candidate")
if (file.exists(out)) stop("Figure 4 WGCNA export already exists: ", out, call. = FALSE)

DS <- c("neuron_neuropil", "neuron_soma", "microglia")
GROUP_OF_FAMILY <- c(rna_processing_splicing_rnp = "rna_chromatin", ribosome_translation = "rna_chromatin",
                     chromatin_organization = "rna_chromatin", mitochondrial_respiration_oxphos = "mitochondria",
                     synaptic_signaling_vesicle = "synaptic_projection", neuron_projection_development = "synaptic_projection",
                     autophagy_lysosome_endosome = "autophagy_endolysosomal")
GROUPS <- c(unique(unname(GROUP_OF_FAMILY)), "figure2_links")
FAMILY_SHARE_MIN <- 0.20
NEAR_SHARE_MIN <- 0.10

inputs <- data.frame(role = character(), path = character(), stringsAsFactors = FALSE)
rd <- function(role, ...) {
  p <- rel(...)
  if (!file.exists(p)) stop("missing input: ", p, call. = FALSE)
  inputs <<- rbind(inputs, data.frame(role = role, path = p, stringsAsFactors = FALSE))
  utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
}
module_number <- function(id) as.integer(sub("^WGCNA_m", "", id))

# ---- labels and sizes
reg <- do.call(rbind, lapply(DS, function(d) rd("reviewed_label_registry", "config", "wgcna_labels", paste0(d, ".csv"))[
  , c("dataset", "level", "entity_id", "reviewed_biological_label", "reviewed_short_label", "label_confidence")]))
reg <- reg[reg$level == "module", , drop = FALSE]
mods <- do.call(rbind, lapply(DS, function(d) {
  l <- rd("wgcna_module_membership", "results", "wgcna", "build_wgcna_modules", d, "tables", "modules", "WGCNA_modules_long.csv")
  n <- table(l$ModuleID)
  data.frame(dataset = d, module_id = names(n), n_proteins = as.integer(n), stringsAsFactors = FALSE)
}))
mods <- merge(mods, data.frame(dataset = reg$dataset, module_id = reg$entity_id, display_label = reg$reviewed_biological_label,
                               short_label = reg$reviewed_short_label, label_confidence = reg$label_confidence,
                               stringsAsFactors = FALSE), by = c("dataset", "module_id"), all.x = TRUE)
if (anyNA(mods$display_label)) stop("A module has no reviewed label: ", paste(mods$dataset[is.na(mods$display_label)],
                                    mods$module_id[is.na(mods$display_label)], collapse = ", "), call. = FALSE)

# ---- stability: bilateral reproducibility
rep_tab <- rd("bilateral_reproducibility", "results", "spatial_validation", "quantify_module_bilateral_identity", "global", "tables",
              "WGCNA_module_bilateral_reproducibility.csv")
rep_tab <- rep_tab[rep_tab$level == "module", , drop = FALSE]
mods <- merge(mods, data.frame(dataset = rep_tab$dataset, module_id = rep_tab$endpoint_id,
                               stability_profile_r = rep_tab$median_profile_Pearson, stability_level_r = rep_tab$pearson_r,
                               stability_n_animals = rep_tab$n_animals, stability_class = rep_tab$bilateral_reproducibility_class,
                               stringsAsFactors = FALSE), by = c("dataset", "module_id"), all.x = TRUE)
if (anyNA(mods$stability_class)) stop("A module lacks a bilateral reproducibility row.", call. = FALSE)
mods$fragile <- mods$stability_class == "poor_bilateral_reproducibility"
mods$stability_label <- sprintf("%.2f", mods$stability_profile_r)

# ---- control spatial class
geo <- c("results", "spatial_validation", "describe_module_control_geography", "global", "tables")
cls <- do.call(rd, c(list("control_spatial_class"), as.list(geo), list("module_control_spatial_class.csv")))
mods <- merge(mods, cls[c("dataset", "module_id", "peak_unit", "animals_peak", "hemis_peak", "shared_top_units",
                          "unit_partial_R2", "p_unit", "class_rule", "spatial_class", "tau")], by = c("dataset", "module_id"), all.x = TRUE)
if (anyNA(mods$spatial_class) || nrow(mods) != nrow(cls)) stop("Control spatial class and modules disagree.", call. = FALSE)
prof <- do.call(rd, c(list("control_profiles"), as.list(geo), list("module_control_profiles.csv")))
prof <- prof[c("dataset", "module_id", "spatial_unit", "con_mean", "con_se", "n_con_animals", "con_z", "sd_unit", "con_sd",
               "spatial_class", "peak_unit", "is_dominant_unit")]
terr <- do.call(rd, c(list("territory_plane"), as.list(geo), list("neuropil_module_territory_plane.csv")))
terr <- terr[c("dataset", "module_id", "pp_minus_other", "pp_ci_low", "pp_ci_high", "mf_minus_other", "mf_ci_low", "mf_ci_high",
               "pp_n_animals_positive", "mf_n_animals_positive", "n_animals", "reference_laminae", "territory_R2", "p_territory_exact",
               "q_territory", "region_R2", "p_region_exact", "q_region", "sd_unit")]

# ---- Figure 2 links
lk <- rd("home_term_links", "results", "spatial_validation", "link_home_terms_to_modules", "global", "tables", "home_term_module_links.csv")
links <- lk[lk$linked %in% TRUE, c("dataset", "module_id", "spatial_unit", "figure2_panel", "GO_ID", "GO_term", "n_term",
                                   "null_p_link", "peaks_without_term_proteins")]

# ---- family composition (module level, primary variant)
long <- rd("theme_module_composition", "results", "integration", "test_theme_module_composition", "global", "tables",
           "theme_entity_composition_long.csv")
long <- long[long$level == "module" & long$is_primary_variant %in% TRUE, , drop = FALSE]
if (!all(long$theme_id %in% names(GROUP_OF_FAMILY))) stop("Unknown Figure 3 family in the composition table.", call. = FALSE)
comp <- data.frame(dataset = long$dataset, module_id = long$entity_id, theme_id = long$theme_id, n_theme = long$n_theme,
                   n_overlap = long$n_overlap, share_of_family = long$frac_theme_in_entity,
                   module_share_of_universe = long$n_entity / long$n_universe, odds_ratio = long$odds_ratio, BH_FDR = long$BH_FDR,
                   carried = long$composition_supported %in% TRUE, stringsAsFactors = FALSE)
if (anyDuplicated(comp[c("dataset", "module_id", "theme_id")])) stop("Duplicated composition rows.", call. = FALSE)
car <- comp[comp$carried, , drop = FALSE]
car <- car[order(car$dataset, car$module_id, -car$share_of_family), , drop = FALSE]
top <- car[!duplicated(car[c("dataset", "module_id")]), c("dataset", "module_id", "theme_id")]
names(top)[3] <- "top_family"
mods <- merge(mods, top, by = c("dataset", "module_id"), all.x = TRUE)

# ---- bridge rows
fam_rows <- car[car$share_of_family >= FAMILY_SHARE_MIN, , drop = FALSE]
has_row <- unique(paste(fam_rows$dataset, fam_rows$theme_id))   # every family a row carries, not only its largest
fam_rows <- fam_rows[!duplicated(fam_rows[c("dataset", "module_id")]), , drop = FALSE]   # largest qualifying family first
fam_rows$bridge_group <- unname(GROUP_OF_FAMILY[fam_rows$theme_id])
fam_rows$bridge_reason <- sprintf("carries %s (%.0f%% of the family, odds ratio %.2f, BH FDR %.2g)", fam_rows$theme_id,
                                  100 * fam_rows$share_of_family, fam_rows$odds_ratio, fam_rows$BH_FDR)
link_mods <- unique(links[c("dataset", "module_id")])
link_only <- link_mods[!paste(link_mods$dataset, link_mods$module_id) %in% paste(fam_rows$dataset, fam_rows$module_id), , drop = FALSE]
if (nrow(link_only)) { link_only$bridge_group <- "figure2_links"; link_only$bridge_reason <- "carries a Figure 2 home-term link" }
rows <- rbind(fam_rows[c("dataset", "module_id", "bridge_group", "bridge_reason")], link_only)
rows <- rows[order(match(rows$bridge_group, GROUPS), match(rows$dataset, DS), module_number(rows$module_id)), , drop = FALSE]
rows$bridge_order <- seq_len(nrow(rows))
mods <- merge(mods, rows, by = c("dataset", "module_id"), all.x = TRUE)
mods$in_bridge <- !is.na(mods$bridge_order)
mods <- mods[order(match(mods$dataset, DS), module_number(mods$module_id)), , drop = FALSE]

# ---- near misses: families without a row in a compartment, and the modules closest to the rule
near <- comp[comp$BH_FDR <= 0.05 & comp$share_of_family >= NEAR_SHARE_MIN &
             !(comp$carried & comp$share_of_family >= FAMILY_SHARE_MIN), , drop = FALSE]
near <- near[!paste(near$dataset, near$theme_id) %in% has_row, , drop = FALSE]
near$fails <- ifelse(!near$carried & near$share_of_family >= FAMILY_SHARE_MIN, "odds ratio below 2",
              ifelse(near$carried, "below 20% of the family", "odds ratio below 2 and below 20% of the family"))
near <- near[order(match(near$dataset, DS), near$theme_id, -near$share_of_family), , drop = FALSE]

# ---- Figure 3 exemplars: share of the SUS-RES leading edge in each module
sel <- rd("figure_03_selection", "exports", "publication_source_data", "figure_03_adaptation", "selection.csv")
curves <- rd("figure_03_running_curves", "exports", "publication_source_data", "figure_03_adaptation", "running_enrichment_curves.csv")
ga <- rd("gsea_leading_edges", "results", "integration", "test_enrichment_module_concordance", "global", "tables",
         "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
ex <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
  s <- sel[i, ]
  e <- unique(ga$leading_edge_proteins[ga$dataset == s$dataset & ga$contrast == "SUS - RES" & ga$GO_ID == s$GO_ID &
                                       tolower(ga$spatial_unit) == tolower(s$spatial_unit)])
  if (length(e) != 1L) stop("Exemplar ", s$exemplar, ": no unique leading edge.", call. = FALSE)
  pg <- strsplit(e, ";", fixed = TRUE)[[1]]
  z <- curves[curves$exemplar == s$exemplar, , drop = FALSE]
  pk <- z$rank[z$peak %in% TRUE][1]; es <- z$running_ES[z$peak %in% TRUE][1]
  n_curve <- if (es > 0) sum(z$hit & z$rank <= pk) else sum(z$hit & z$rank >= pk)
  if (n_curve != length(pg)) stop("Exemplar ", s$exemplar, ": leading edge differs from the Figure 3 curve.", call. = FALSE)
  l <- utils::read.csv(rel("results", "wgcna", "build_wgcna_modules", s$dataset, "tables", "modules", "WGCNA_modules_long.csv"),
                       stringsAsFactors = FALSE)[c("ProteinGroupID", "ModuleID")]
  m <- l[l$ProteinGroupID %in% pg, , drop = FALSE]
  n <- table(m$ModuleID)
  data.frame(module_id = names(n), n_in_module = as.integer(n), exemplar = s$exemplar, dataset = s$dataset,
             spatial_unit = s$spatial_unit, GO_ID = s$GO_ID, n_leading_edge = length(pg), n_in_network = nrow(m),
             share_of_leading_edge = as.integer(n) / nrow(m), stringsAsFactors = FALSE)
}))
ex <- ex[order(ex$exemplar, -ex$n_in_module), , drop = FALSE]

# ---- checks and write
if (!all(paste(links$dataset, links$module_id) %in% paste(mods$dataset, mods$module_id))) stop("A link names an unknown module.", call. = FALSE)
n_units <- tapply(prof$spatial_unit, prof$dataset, function(u) length(unique(u)))
n_mods <- table(mods$dataset)
if (anyDuplicated(prof[c("dataset", "module_id", "spatial_unit")]) ||
    !identical(as.integer(table(prof$dataset)[DS]), as.integer(n_mods[DS] * n_units[DS])))
  stop("The profile table does not have one row per module and unit.", call. = FALSE)
win_long <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows" && nchar(path) >= 240L) paste0("\\\\?\\", gsub("/", "\\\\", path)) else path
}
inputs <- unique(inputs)
inputs$sha256 <- unname(tools::sha256sum(vapply(inputs$path, win_long, character(1))))
inputs$path <- substring(gsub("\\\\", "/", inputs$path), nchar(root) + 2L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
w <- function(x, f) utils::write.csv(x, file.path(out, f), row.names = FALSE)
w(mods, "modules.csv")
w(prof[order(match(prof$dataset, DS), module_number(prof$module_id), prof$spatial_unit), ], "module_profiles.csv")
w(terr, "territory_plane.csv")
w(links, "figure2_links.csv")
w(comp[order(match(comp$dataset, DS), module_number(comp$module_id), comp$theme_id), ], "family_composition.csv")
w(near[c("dataset", "theme_id", "module_id", "share_of_family", "odds_ratio", "BH_FDR", "carried", "fails")], "family_near_misses.csv")
w(ex, "exemplar_modules.csv")
w(inputs, "input_manifest.csv")
cat("Exported the Figure 4 WGCNA candidate to ", out, "\n", sep = "")
print(mods[mods$in_bridge, c("dataset", "module_id", "short_label", "bridge_group", "spatial_class", "peak_unit", "stability_label",
                             "fragile", "top_family")][order(mods$bridge_order[mods$in_bridge]), ], row.names = FALSE)
