# Inventory of every analysis against the two sample-centring versions (normalisation review 2026-10-21). Read-only.
#
# Rows: every pipeline.yml step; the legacy preprocessing chain that produced the analysis inputs (impute, metadata
# join, animal-level input); and the steps outside the registry that read or export results (tools/ exporters and
# review tools, the export trees, the Exp9_manuscript side). For each row and dataset:
#   - outputs in the working tree (complete-case centring, rerun 2026-10-19/20) and in the archive (original
#     centring): by path, and by file name anywhere in the archive manifest, because most archived originals sit at
#     legacy paths;
#   - whether the step depends on the centring, with the reason; reproducibility of each version;
#   - historical dependencies: label registries, curated mappings, outcome-derived definitions, literal module IDs
#     and hard-coded legacy result paths in the code;
#   - what the rerun read (input-resolution ledger, rerun window only, last read per input, by dataset) and how much of
#     the step the ledger covers;
#   - an ordering check: a step whose oldest output predates the regeneration of what it reads (network, labels,
#     differential tables, enrichment) ran on the previous inputs;
#   - last run status, resulting status and action.
# Also writes analysis_families.csv (each analysis family requested for review mapped to its steps and status).
#
# USAGE  Rscript tools/normalisation_review/build_normalisation_inventory.R        (working directory: repo root)

suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(purrr) })
ARCHIVE <- if (dir.exists("Q:/results")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
OUT <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "inventory")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RERUN_START <- as.POSIXct("2026-10-19 12:00:00")
LEDGER_RERUN_ROWS <- c(57858L, 59249L)   # the ledger rows written by the rerun of 2026-10-19/20
ORIGINAL_COMMIT <- "b478a13"
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
ftime <- function(x) ifelse(is.na(x), NA_character_, format(x, "%Y-%m-%d %H:%M %Z"))

reg <- yaml::read_yaml("pipeline.yml")
steps <- imap_dfr(reg$stages, function(st, stage_name) map_dfr(st$scripts, function(s) tibble(
  stage = stage_name, script = s$script, required = isTRUE(s$required),
  datasets = list(unlist(s$datasets %||% "global")), produces = list(unlist(s$produces %||% character())),
  consumes = list(c(unlist(s$consumes_required %||% character()), unlist(s$consumes_optional %||% character()))),
  notes = s$notes %||% "")))
steps <- distinct(steps, stage, script, .keep_all = TRUE)
# legacy preprocessing chain (rerun manually on 2026-10-19) and the steps outside the registry
extra <- tribble(
  ~stage, ~script, ~datasets, ~produces, ~notes,
  "preprocessing (legacy, run manually)", "archive/01_preprocessing/01_impute.r", "global", "data/processed/01_preprocessing/impute", "the centring step itself",
  "preprocessing (legacy, run manually)", "archive/01_preprocessing/02_excel_convert.r", "global", "data/processed/01_preprocessing/excel_convert", "metadata-joined workbooks",
  "preprocessing (legacy, run manually)", "analysis/preprocessing/build_animal_level_protigy_input.R", "dataset", "data/processed/01_preprocessing/protigy_input_animal_level/<dataset>", "animal-level input of the differential model",
  "outside registry: exporters", "tools/export_publication_source_data.R", "global", "exports/publication_source_data", "frozen source-data bundle",
  "outside registry: exporters", "tools/export_figure_03_adaptation_source.R", "global", "exports/publication_source_data/figure_03_adaptation", "Figure 3 adaptation classes source",
  "outside registry: exporters", "tools/export_figure_03_go_atlas_appendix.R", "global", "exports/publication_source_data/figure_03_go_atlas_appendix", "Figure 3 GO atlas appendix",
  "outside registry: exporters", "tools/export_supplementary_go_pathways_candidate.R", "global", "exports/publication_source_data/supplementary_go_pathways_candidate", "candidate",
  "outside registry: exporters", "tools/export_supplementary_go_protein_zooms_candidate.R", "global", "exports/publication_source_data/supplementary_go_protein_zooms_candidate", "candidate",
  "outside registry: exporters", "tools/export_wgcna_theme_composition_candidate.R", "global", "exports/publication_source_data/wgcna_theme_composition_candidate", "candidate",
  "outside registry: exporters", "tools/generate_publication_freeze_manifest.R", "global", "docs/publication_freeze_manifest.yml", "freeze manifest",
  "outside registry: exporters", "tools/migrate_publication_export_bundle.R", "global", "exports/publication_source_data", "bundle migration",
  "outside registry: review tools", "tools/map_wgcna_module_lineage.R", "global", "results/wgcna/map_wgcna_module_lineage", "previous-to-new network lineage (reads both versions by design)",
  "outside registry: review tools", "tools/attribute_wgcna_change.R", "global", "results/wgcna/attribute_wgcna_change", "network change attribution (reads both versions by design)",
  "outside registry: review tools", "tools/build_renormalisation_workbooks.R", "global", "results/reviewer_audit/renormalisation_2026-10-19/workbooks", "review workbooks (both versions by design)",
  "outside registry: export trees", "results/manuscript (08 export)", "global", "results/manuscript", "frozen manuscript figure tree",
  "outside registry: export trees", "pride_submission (PRIDE export)", "global", "pride_submission", "deposit tree",
  "outside registry: Exp9_manuscript", "Exp9_manuscript import tools and figure/ED builds", "global", "", "reads P:/exports; not checked by the normalisation guard (other repository)",
  "outside registry: discovery audits", "Proteome_state/code/proteome_state.R", "global", "Proteome_state", "proteome-state discovery audit 2026-10-17/19 (CAMERA, programme tests; untracked)",
  "outside registry: discovery audits", "audits/part29/camera_summaries.R", "global", "audits/part29", "CAMERA summaries (audit)",
  "outside registry: discovery audits", "Behavior_coupling/code/proteomics_scores.R", "global", "Behavior_coupling", "proteomics-behaviour coupling scores (built by Behavior_coupling/code/rebuild_proteomics_coupling.R; untracked)") |>
  mutate(required = NA, datasets = map(datasets, ~ if (.x == "dataset") DATASETS else "global"), produces = map(produces, ~ if (nzchar(.x)) .x else character()),
         consumes = list(character()))
steps <- bind_rows(steps, extra)
steps$id <- seq_len(nrow(steps))

INDEPENDENT <- c(
  "analysis/preprocessing/build_joint_protigy_input.R" = "own joint shared-core normalisation of the raw matrix",
  "analysis/qc/assess_joint_compartment_quality.R" = "joint shared-core bundle only",
  "analysis/qc/render_joint_compartment_qc_figures.R" = "joint shared-core bundle only",
  "analysis/qc/render_compartment_abundance_figures.R" = "joint shared-core offsets (04e), not the per-compartment centring",
  "analysis/qc/build_reference_marker_registry.R" = "marker configuration only",
  "analysis/qc/assess_sample_quality.R" = "DIA-NN acquisition statistics only",
  "analysis/publication_source_data/build_sample_metadata.R" = "sample metadata only",
  "analysis/preprocessing/build_module_score_metadata.R" = "sample metadata only (the protein matrix is resolved but not used)",
  "analysis/integration/audit_animal_id_integrity.R" = "sample and animal identifiers only")
# Stale reads found by the adversarial check of 2026-10-21 that the code-text markers below do not catch.
KNOWN_STALE <- tribble(~script, ~dataset, ~reason,
  "analysis/qc/summarize_marker_detectability.R", "neuron_neuropil", "ran on 2026-10-19 18:19, before the rebuilt network existed; module-marker tables computed on the previous network",
  "analysis/qc/summarize_marker_detectability.R", "neuron_soma", "ran on 2026-10-19 18:22, before the rebuilt network existed; module-marker tables computed on the previous network",
  "analysis/differential_abundance/build_sus_res_dap_atlas.R", "global", "reads the original-centring ranked GO atlas and supported-GO audit (results/source_data/04_differential_expression_enrichment/compareGO_spatial_atlas, 2026-08-26)",
  "analysis/integration/render_module_circular_atlas.R", "global", "claim-status track reads results/tables/biological_claims_table.csv (2026-09-11, original centring)",
  "analysis/differential_abundance/test_microglia_targeted_signatures.R", "microglia", "ran 2026-10-19 20:23 and read the original-centring EWCE results (data/processed/05_celltype_enrichment_EWCE/EWCE_E9/microglia/EWCE_results_full.rds, 2026-08-21) before the new EWCE run (2026-10-20 00:40)",
  "Proteome_state/code/proteome_state.R", "global", "discovery audit run 2026-10-19 09:27, before the rerun (original centring; CAMERA, programme tests)",
  "audits/part29/camera_summaries.R", "global", "CAMERA summaries of the original-centring results (audit, not part of the pipeline)",
  "Behavior_coupling/code/proteomics_scores.R", "global", "built 2026-10-19 12:07, before the rebuilt networks; reads results/tables/06_modules_WGCNA (original centring) eigengenes and theme scores")
# notes on steps that write into the legacy result trees (these files are new-version files inside old-version trees)
WRITES_LEGACY <- c(
  "analysis/spatial_validation/summarize_ca2_slm_robustness.R" = "writes CA2_SLM tables into results/tables/11_spatial_systems/ca2_slm_robustness",
  "analysis/preprocessing/build_animal_level_protigy_input.R" = "writes its aggregation tables into results/tables/01_preprocessing (read by the normalisation guard)")
# a step is not downstream of what it produces itself or feeds
NOT_UPSTREAM <- list("analysis/wgcna/build_wgcna_modules.R" = c("labels"))
ORIGINAL_ONLY <- c(
  "analysis/wgcna/build_curated_overlap_programs.R" = "defined from recurrent compareGO proteins of the original-centring results (legacy table); not regenerated",
  "analysis/wgcna/compare_recurrent_module_proteins.R" = "reads legacy compareGO top-gene workbooks of the original-centring results; not regenerated")
steps <- steps |> mutate(
  depends_on_centring = !script %in% names(INDEPENDENT),
  dependence_path = case_when(script %in% names(INDEPENDENT) ~ paste("independent:", INDEPENDENT[script]),
                              script %in% names(ORIGINAL_ONLY) ~ paste("original centring only:", ORIGINAL_ONLY[script]),
                              script == "archive/01_preprocessing/01_impute.r" ~ "applies the centring",
                              TRUE ~ "reads the imputed matrices or results derived from them"))

# ---- historical dependencies
hist_patterns <- c("config/wgcna_labels" = "reviewed WGCNA label registry",
  "gsea_wgcna_theme_module_mapping" = "theme-module mapping (curated on the previous network)",
  "gsea_wgcna_program_module_mapping" = "program-module mapping (curated on the previous network)",
  "wgcna_supermodules" = "supermodule configuration", "curated_overlap" = "curated overlap programmes (defined from previous DA results)",
  "compareGO/overlapping_proteins" = "recurrent compareGO proteins (previous DA results)",
  "manuscript_go_theme_registry" = "theme registry (predefined, not outcome-selected)")
LEGACY_CODE <- 'path_results\\(\\s*"tables"|"results/tables/|results/source_data/|"data/processed/0[2-9]_|path_processed\\(\\s*"0[2-9]_|"results/tables"|biological_claims_table[.]csv"'
code_of <- function(s) if (file.exists(s)) { x <- readLines(s, warn = FALSE); x[!grepl("^\\s*#", x)] } else character()
changed_files <- system2("git", c("diff", "--name-only", ORIGINAL_COMMIT), stdout = TRUE)
untracked <- system2("git", c("ls-files", "--others", "--exclude-standard", "--", "analysis", "R", "tools", "archive"), stdout = TRUE)
steps <- steps |> mutate(
  code = map(script, code_of),
  historical_inputs = map2_chr(consumes, code, function(cs, cd) paste(unique(unlist(lapply(names(hist_patterns), function(p)
    if (any(grepl(p, c(cs, cd), fixed = TRUE))) hist_patterns[[p]]))), collapse = "; ")),
  literal_module_ids = map_chr(code, function(cd) { ids <- unique(unlist(regmatches(cd, gregexpr("WGCNA_m[0-9]{2}|\\bSM[0-9]{2}\\b", cd))))
    if (length(ids)) paste(head(sort(ids), 12), collapse = " ") else "" }),
  legacy_paths_in_code = map_chr(code, function(cd) { h <- unique(trimws(cd[grepl(LEGACY_CODE, cd, perl = TRUE)]))
    if (length(h)) paste(substr(head(h, 3), 1, 110), collapse = " || ") else "" }),
  code_changed_since_original = script %in% c(changed_files, untracked),
  upstream_markers = map(code, function(cd) {
    j <- paste(cd, collapse = "\n")
    c(network = grepl("wgcna_final_model_state|WGCNA_modules_long|modules_long|wgcna_modules_artifact", j),
      labels = grepl("final_label_lookup|WGCNA_final_label_lookup|summarize_module_interpretation", j),
      differential = grepl("map_protein_identifiers|preprocessing_mapped_contrast_dir|02_id_mapping", j),
      enrichment = grepl("clusterProfiler_manifest|run_clusterprofiler_enrichment", j)) }))
PRODUCER_OF <- c(network = "analysis/wgcna/build_wgcna_modules.R", labels = "analysis/wgcna/summarize_module_interpretation.R",
                 differential = "analysis/preprocessing/map_protein_identifiers.R", enrichment = "analysis/differential_abundance/run_clusterprofiler_enrichment.R")
upstream_time <- function(kind, ds) {
  dss <- if (ds == "global") DATASETS else ds
  f <- switch(kind,
    network = file.path("results/wgcna/build_wgcna_modules", dss, "models/wgcna_final_model_state.rds"),
    labels = file.path("results/wgcna/summarize_module_interpretation", dss, "tables/WGCNA_final_label_lookup.csv"),
    differential = unlist(lapply(dss, function(d) Sys.glob(file.path("results/preprocessing/map_protein_identifiers", d, "tables/mapped/forward/per_file/*.csv"))[1])),
    enrichment = file.path("results/differential_abundance/run_clusterprofiler_enrichment", dss, "models/clusterProfiler_manifest.csv"))
  f <- f[!is.na(f) & file.exists(f)]
  if (length(f)) max(file.info(f)$mtime) else as.POSIXct(NA)
}

# ---- outputs in each version
arch <- gsub("\\\\", "/", read_csv(file.path(ARCHIVE, "SHA256_MANIFEST.csv"), show_col_types = FALSE, col_select = 1)[[1]])
arch_names <- table(basename(arch))
# archived files matching a registry output pattern, from the archive manifest (listing the archive is slow)
arch_count <- function(pattern, ds) {
  p <- gsub("<dataset>", ds, pattern, fixed = TRUE)
  rx <- gsub("([.+^$(){}|])", "\\\\\\1", p); rx <- gsub("*", "[^/]*", rx, fixed = TRUE); rx <- gsub("?", "[^/]", rx, fixed = TRUE)
  sum(grepl(paste0("^", sub("/$", "", rx), "($|/)"), arch))
}
GENERIC <- "^(run_manifest|manifest|sessionInfo|README|session_info|checksums|index)[.]|[.]log$"
expand_one <- function(pattern, ds, root) {
  p <- file.path(root, gsub("<dataset>", ds, pattern, fixed = TRUE))
  if (grepl("[*?\\[]", p)) return(Sys.glob(p))
  if (dir.exists(p)) return(list.files(p, recursive = TRUE, full.names = TRUE))
  if (file.exists(p)) p else character()
}
rows <- list()
for (i in seq_len(nrow(steps))) {
  dss <- steps$datasets[[i]]; if (identical(dss, "global") || !length(dss)) dss <- "global"
  for (ds in dss) {
    pr <- steps$produces[[i]]
    nf <- unique(as.character(unlist(lapply(pr, expand_one, ds = ds, root = "."))))
    n_old <- sum(vapply(pr, arch_count, numeric(1), ds = ds))
    nm <- unique(basename(nf)); nm <- nm[!grepl(GENERIC, nm)]
    info <- if (length(nf)) file.info(nf) else NULL
    up <- steps$upstream_markers[[i]]; up <- names(up)[up & PRODUCER_OF[names(up)] != steps$script[i]]
    up <- setdiff(up, NOT_UPSTREAM[[steps$script[i]]] %||% character())
    up_t <- vapply(up, function(k) as.numeric(upstream_time(k, ds)), numeric(1))
    oldest <- if (length(nf)) min(info$mtime, na.rm = TRUE) else as.POSIXct(NA)
    before <- up[!is.na(up_t) & !is.na(oldest) & as.numeric(oldest) < up_t]
    rows[[length(rows) + 1L]] <- tibble(id = steps$id[i], dataset = ds, n_outputs_working_tree = length(nf),
      oldest_output_working_tree = oldest, newest_output_working_tree = if (length(nf)) max(info$mtime, na.rm = TRUE) else as.POSIXct(NA),
      n_outputs_archive_same_path = n_old,
      archive_has_same_named_outputs = if (length(nm)) sprintf("%d of %d file names", sum(nm %in% names(arch_names)), length(nm)) else "",
      ran_before_upstream_regenerated = paste(before, collapse = ", "))
  }
}
outs <- bind_rows(rows)

# ---- what the rerun read (ledger, rerun window), by script and dataset
ledger <- read_csv("results/reviewer_audit/input_resolution_audit.csv", show_col_types = FALSE, guess_max = 1e5)
ledger <- ledger[LEDGER_RERUN_ROWS[1]:LEDGER_RERUN_ROWS[2], ]
ledger$row <- seq_len(nrow(ledger)); ledger$key <- sub("[.][Rr]$", "", basename(ledger$script))
unattributed <- sum(is.na(ledger$script))
last <- ledger |> filter(!is.na(script)) |> group_by(script, dataset, input_name) |> slice_max(row, n = 1, with_ties = FALSE) |> ungroup()
mt <- as.POSIXct(substr(last$file_mtime, 1, 19), format = "%Y-%m-%d %H:%M:%S")
last$stale_derived <- grepl("/(results|data/processed)/", last$resolved_path) & !is.na(mt) & mt < RERUN_START &
  !grepl("06_merged_metadata_module_score", last$resolved_path, fixed = TRUE)
# every last read of a pre-rerun derived input; merged sample metadata carries no protein values (resolver record only)
pre_any <- grepl("/(results|data/processed)/", last$resolved_path) & !is.na(mt) & mt < RERUN_START
write_csv(last[pre_any, ] |> transmute(script, dataset, input_name, file_mtime, resolved_path = sub("^.*proteomics/", "", resolved_path),
                                       note = ifelse(grepl("06_merged_metadata_module_score", resolved_path, fixed = TRUE),
                                                     "sample metadata (no protein values); resolver record", "original-centring derived input")),
          file.path(OUT, "last_reads_of_pre_rerun_derived_inputs.csv"), na = "")
reads <- last |> mutate(dataset = ifelse(dataset %in% DATASETS, dataset, "global")) |> group_by(key, dataset) |> summarise(
  ledger_inputs = n(), ledger_inputs_legacy = sum(grepl("legacy", resolution_mode)),
  last_read_pre_rerun_derived_inputs = sum(stale_derived),
  pre_rerun_derived_inputs = paste(head(unique(sub("^.*/(results|data/processed)/", "\\1/", resolved_path[stale_derived])), 4), collapse = " | "),
  .groups = "drop")
man <- list.files("results/logs/pipeline", pattern = "^pipeline_manifest_.*[.]csv$", full.names = TRUE)
man <- man[file.info(man)$mtime >= RERUN_START & file.info(man)$mtime < as.POSIXct("2026-10-21 00:00:00")]
runs <- bind_rows(lapply(man, function(f) read_csv(f, show_col_types = FALSE, col_types = cols(.default = "c")) |>
  transmute(script, dataset, status, finished_at))) |> arrange(finished_at) |> group_by(script, dataset) |> slice_tail(n = 1) |> ungroup()

inv <- outs |> left_join(steps |> select(-produces, -consumes, -datasets, -code, -upstream_markers), by = "id") |>
  mutate(key = sub("[.][Rr]$", "", basename(script))) |>
  left_join(reads, by = c("key", "dataset")) |>
  left_join(runs |> rename(last_run_status = status, last_run_finished = finished_at), by = c("script", "dataset")) |>
  left_join(KNOWN_STALE |> rename(known_stale_reason = reason), by = c("script", "dataset")) |>
  mutate(
    regenerated_in_rerun = !is.na(newest_output_working_tree) & newest_output_working_tree >= RERUN_START,
    original_reproducible = case_when(!depends_on_centring ~ "not affected",
      script %in% names(ORIGINAL_ONLY) ~ "not regenerable through the setting (its input is missing; only a superseded copy exists)",
      grepl("export trees|Exp9_manuscript", stage) ~ "not applicable (frozen tree or other repository)",
      grepl("discovery audits", stage) ~ "regenerable only by rerunning the audit",
      code_changed_since_original ~ "regenerable through the setting; not byte-identical (code changed since b478a13)",
      TRUE ~ "regenerable through the setting with unchanged code"),
    new_reproducible = case_when(!depends_on_centring ~ "not affected", grepl("export trees|Exp9_manuscript", stage) ~ "not applicable (rebuilt by the export stage)",
                                 TRUE ~ "regenerable with the current (uncommitted) code"),
    writes_into_legacy_tree = unname(ifelse(script %in% names(WRITES_LEGACY), WRITES_LEGACY[script], "")),
    status = case_when(
      !depends_on_centring ~ "independent of centring",
      script %in% names(ORIGINAL_ONLY) ~ "original centring only (outcome-derived definitions; decision needed)",
      grepl("Exp9_manuscript", stage) ~ "held: not guarded (other repository); add a centring check before the next import",
      grepl("export trees|exporters", stage) | stage == "export" ~ "held: holds or produces the original-version export bundle; blocked until the decision",
      grepl("review tools", stage) ~ "review tool (reads both versions by design)",
      !is.na(known_stale_reason) & grepl("discovery audits", stage) ~ paste("stale (audit):", known_stale_reason),
      grepl("^failed", coalesce(last_run_status, "")) & !(regenerated_in_rerun & n_outputs_working_tree > 0 &
        newest_output_working_tree > as.POSIXct(substr(last_run_finished, 1, 19), format = "%Y-%m-%d %H:%M:%S")) ~ "pending (last rerun failed)",
      n_outputs_working_tree == 0 ~ "pending (no outputs)",
      !regenerated_in_rerun ~ "stale (not regenerated)",
      !is.na(known_stale_reason) ~ paste("stale:", known_stale_reason),
      nzchar(ran_before_upstream_regenerated) ~ "stale: ran before its inputs were regenerated",
      coalesce(last_read_pre_rerun_derived_inputs, 0L) > 0 ~ "rerun, but read an original-centring input",
      nzchar(legacy_paths_in_code) ~ "rerun; check hard-coded legacy paths",
      grepl("label|theme|module mapping|curated|registry", historical_inputs) | nzchar(literal_module_ids) ~ "provisional (depends on labels/curated definitions)",
      grepl("^failed", coalesce(last_run_status, "")) ~ "complete (rerun outside the launcher after a failed launcher attempt)",
      TRUE ~ "complete (rerun)"),
    action = case_when(
      grepl("^stale", status) ~ "rerun under the chosen centring (after its upstream), or retire",
      grepl("^original centring only", status) ~ "decide: retire, or keep as a predefined descriptive analysis labelled as defined on the original results",
      grepl("original-centring input", status) ~ "decide the fate of the input listed, then rerun",
      grepl("legacy paths", status) ~ "check whether the hard-coded path is a live read; repoint and rerun if so",
      grepl("^pending", status) ~ "resolve and rerun after the decision",
      grepl("^held", status) ~ "rework export and import after the decision",
      grepl("^provisional", status) ~ "rerun after label/definition review",
      TRUE ~ ""),
    across(c(oldest_output_working_tree, newest_output_working_tree), ftime)) |>
  select(stage, script, dataset, required, depends_on_centring, dependence_path, n_outputs_archive_same_path, archive_has_same_named_outputs,
         n_outputs_working_tree, oldest_output_working_tree, newest_output_working_tree, regenerated_in_rerun, ran_before_upstream_regenerated,
         original_reproducible, new_reproducible, last_run_status, ledger_inputs, ledger_inputs_legacy, last_read_pre_rerun_derived_inputs,
         pre_rerun_derived_inputs, historical_inputs, literal_module_ids, legacy_paths_in_code, writes_into_legacy_tree, code_changed_since_original, status, action, notes)
write_csv(inv, file.path(OUT, "normalisation_dependence_inventory.csv"), na = "")

# ---- analysis families requested for review
FAMILIES <- list(
  "differential abundance" = "fit_animal_level_limma|extract_protigy_contrasts|map_protein_identifiers|dap_atlas|differential_abundance/(summarize|audit|quantify)",
  "GSEA (GO)" = "run_clusterprofiler_enrichment|audit_gsea", "compareGO" = "compare_go_enrichment|compareGO",
  "theme summaries / programme atlas" = "build_go_program_atlas|summarize_biological_programs|audit_stress_response_biology|theme",
  "adaptation classes (Figure 3c)" = "adaptation", "atlas / PCA" = "atlas|pca|spatial_atlas|empirical_compartment",
  "WGCNA network" = "build_wgcna_modules|module_robustness|module_complex", "module / supermodule enrichment and labels" = "annotate_module|label|interpretation|adjudicate",
  "eigengenes / module scores" = "score_module_activity|summarize_module_scores|eigengene",
  "group and behaviour associations" = "test_module_phenotypes|behaviour|behavior|coupling",
  "CAMERA" = "camera", "leading edges" = "leading|core_enrichment|enrichment_module_concordance|theme_module_composition",
  "EWCE / cell-type enrichment" = "ewce|celltype", "overlap scores" = "curated_overlap|overlap",
  "QC summaries" = "^analysis/qc/", "figures and Extended Data" = "render_|08_export|Exp9_manuscript|figure",
  "spatial networks" = "spatial_network|network", "export and deposit" = "export|pride|freeze|supplementary_tables|methods_summary")
FAMILY_EXTRA <- list(
  "differential abundance" = c("audit_ca2_slm_robustness", "audit_stress_identity_robustness", "validate_control_spatial_identity", "build_animal_level_protigy_input", "summarize_ca2_slm_robustness"),
  "module / supermodule enrichment and labels" = c("build_wgcna_modules"),
  "leading edges" = c("build_go_program_atlas"),
  "overlap scores" = c("score_module_activity", "compare_recurrent_module_proteins", "build_curated_overlap_programs"),
  "eigengenes / module scores" = c("proteomics_scores"),
  "group and behaviour associations" = c("proteomics_scores"),
  "CAMERA" = c("proteome_state", "camera_summaries"))
camera_scripts <- unlist(lapply(unique(inv$script), function(s) if (any(grepl("camera(", code_of(s), fixed = TRUE) | grepl("cameraPR", code_of(s), fixed = TRUE))) s))
fam <- bind_rows(lapply(names(FAMILIES), function(f) {
  extra_hit <- sub("[.][Rr]$", "", basename(inv$script)) %in% (FAMILY_EXTRA[[f]] %||% character())
  hit <- inv |> filter(grepl(FAMILIES[[f]], script, ignore.case = TRUE) | extra_hit | (f == "CAMERA" & script %in% camera_scripts))
  tibble(family = f, n_steps = n_distinct(hit$script), steps = paste(sort(unique(basename(hit$script))), collapse = "; "),
         statuses = paste(sprintf("%s (%d)", names(table(hit$status)), as.integer(table(hit$status))), collapse = "; "))
}))
fam <- bind_rows(fam, tibble(family = "manuscript statements", n_steps = NA_integer_,
  steps = "claims register (Exp9_manuscript docx register: 53 statements, 41 proteomics-related; see Statements sheet)", statuses = "re-evaluated 2026-10-20/21; no wording proposed until the decision"))
write_csv(fam, file.path(OUT, "analysis_families.csv"), na = "")
write_csv(tibble(ledger_rows_in_rerun = nrow(ledger), unattributed_rows = unattributed,
                 dependent_steps = n_distinct(inv$script[inv$depends_on_centring]),
                 dependent_steps_with_ledger_reads = n_distinct(inv$script[inv$depends_on_centring & coalesce(inv$ledger_inputs, 0L) > 0])),
          file.path(OUT, "ledger_coverage.csv"))
print(count(inv, depends_on_centring, status))
cat("done:", file.path(OUT, "normalisation_dependence_inventory.csv"), "\n")
