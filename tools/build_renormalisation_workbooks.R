#!/usr/bin/env Rscript
# Reader-facing Excel workbooks for the re-normalised proteomics results (2026-10-19).
#
# Each workbook opens with a "Read me" sheet and shows the current analysis next to the previous
# median-centring as a sensitivity analysis. Style: Arial 10, dark-blue header, light row rules,
# p values as "<0.001", values below 0.05 highlighted, a blue-white-yellow scale on signed effects,
# frozen panes, no gridlines (the house style of the behaviour and physiology workbooks).
#
#   Proteomics_differential_abundance.xlsx   protein-level results, all compartments
#   (GSEA and WGCNA workbooks are added by later sections once those results exist)
#
# Previous centring = the frozen results before 2026-10-19 (hash-verified archive).
# USAGE  Rscript tools/build_renormalisation_workbooks.R [--sections da,gsea,wgcna] [--out <dir>]

source(if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R"))
source(repo_path("R", "preprocessing_paths.R"))
suppressPackageStartupMessages({ library(openxlsx); library(dplyr); library(readr); library(tidyr) })

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1L]] }
SECTIONS <- strsplit(arg("--sections", "da"), ",", fixed = TRUE)[[1]]
OUT <- arg("--out", path_results("reviewer_audit", "renormalisation_2026-10-19", "workbooks"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
# The archive root is long enough to push some files past R's 260-character path limit on this
# host; `subst Q: <archive>` gives a short alias, used when present.
ARCHIVE <- if (dir.exists("Q:/results")) "Q:" else "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
SUPERSEDED <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_superseded_proteomics_inputs_2026-10-19"
DATASETS <- c(neuron_neuropil = "Neuronal neuropil", neuron_soma = "Neuronal soma", microglia = "Microglia")
FDR_CUT <- 0.05

# ------------------------------------------------------------------ house style
source(repo_path("R", "review_workbook_style.R"))

# ------------------------------------------------------------------ contrast parsing
GROUP <- c(con = "CON", res = "RES", sus = "SUS")
parse_contrast <- function(stem, ds) {
  m <- regmatches(stem, regexec("^(CA1|CA2|CA3|DG)([a-z]+?)(sus|res)_(CA1|CA2|CA3|DG)([a-z]+?)(res|con)$", stem))[[1]]
  if (length(m) != 7) stop("Unrecognised contrast file: ", stem, call. = FALSE)
  layer <- m[3]
  unit <- if (identical(layer, "microglia")) m[2] else paste(m[2], layer)
  tibble(spatial_unit = unit, comparison = paste(GROUP[[m[4]]], "vs", GROUP[[m[7]]]))
}
read_da <- function(dir, ds) {
  files <- sort(list.files(dir, pattern = "\\.csv$", full.names = TRUE))
  bind_rows(lapply(files, function(f) {
    x <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE, guess_max = 1e5)
    bind_cols(parse_contrast(sub("\\.csv$", "", basename(f)), ds)[rep(1, nrow(x)), ],
              transmute(x, ProteinGroupID, gene = coalesce(official_gene_symbol, representative_gene_symbol, gene_symbol),
                        accession = representative_accession, log2fc, aveExpr, t, pval, padj))
  }))
}
UNIT_ORDER <- c("CA1 slm", "CA1 so", "CA1 sr", "CA2 slm", "CA2 so", "CA2 sr", "CA3 so", "CA3 sr", "DG mo", "DG po",
                "CA1 sp", "CA2 sp", "CA3 sp", "DG sg", "CA1", "CA2", "CA3", "DG")
CMP_ORDER <- c("RES vs CON", "SUS vs CON", "SUS vs RES")

# ================================================================== differential abundance
if ("da" %in% SECTIONS) {
  cat("Differential abundance workbook\n")
  da <- bind_rows(lapply(names(DATASETS), function(ds) {
    cur <- read_da(preprocessing_mapped_contrast_dir(ds, "forward"), ds)
    prev <- read_da(file.path(ARCHIVE, "data", "processed", "02_id_mapping", "mapped", ds, "forward", "per_file"), ds) |>
      select(spatial_unit, comparison, ProteinGroupID, prev_log2fc = log2fc, prev_padj = padj)
    full_join(cur, prev, by = c("spatial_unit", "comparison", "ProteinGroupID")) |> mutate(dataset = ds, .before = 1)
  })) |>
    mutate(compartment = factor(DATASETS[dataset], levels = DATASETS),
           spatial_unit = factor(spatial_unit, levels = UNIT_ORDER), comparison = factor(comparison, levels = CMP_ORDER)) |>
    arrange(compartment, spatial_unit, comparison, padj)

  overview <- da |> group_by(compartment, spatial_unit, comparison) |>
    summarise(proteins_tested = sum(!is.na(padj)),
              fdr_sig = sum(padj < FDR_CUT, na.rm = TRUE),
              higher = sum(padj < FDR_CUT & log2fc > 0, na.rm = TRUE),
              lower = sum(padj < FDR_CUT & log2fc < 0, na.rm = TRUE),
              prev_fdr_sig = sum(prev_padj < FDR_CUT, na.rm = TRUE),
              in_both = sum(padj < FDR_CUT & prev_padj < FDR_CUT, na.rm = TRUE), .groups = "drop")
  ov <- overview |> transmute(Compartment = as.character(compartment), `Spatial unit` = as.character(spatial_unit),
                              Comparison = as.character(comparison), `Proteins tested` = proteins_tested,
                              `FDR < 0.05` = fdr_sig, `Higher in first group` = higher, `Lower in first group` = lower,
                              `FDR < 0.05, previous centring` = prev_fdr_sig, `FDR < 0.05 in both` = in_both)
  sig <- da |> filter(padj < FDR_CUT | prev_padj < FDR_CUT) |>
    transmute(Compartment = as.character(compartment), `Spatial unit` = as.character(spatial_unit), Comparison = as.character(comparison),
              Gene = gene, `UniProt accession` = accession, `log2 fold change` = log2fc, t = t, P = pval, FDR = padj,
              `log2 fold change, previous centring` = prev_log2fc, `FDR, previous centring` = prev_padj,
              `Called at FDR < 0.05` = case_when(padj < FDR_CUT & prev_padj < FDR_CUT ~ "both",
                                                 padj < FDR_CUT ~ "current only", TRUE ~ "previous centring only"),
              `Protein group ID` = ProteinGroupID)

  # Per-sample centring offsets and detection depth.
  meta <- readxl::read_excel(path_metadata("TPE9_sample_metadata_males.xlsx"))
  raw <- readr::read_tsv(path_raw("pg_matrix", "quicksearch.pg_matrix.tsv"), show_col_types = FALSE, progress = FALSE)
  detected <- vapply(raw[-(1:4)], function(v) sum(!is.na(v) & v > 0), integer(1))
  offsets <- bind_rows(lapply(names(DATASETS), function(ds) {
    old <- readxl::read_excel(Sys.glob(file.path(SUPERSEDED, "data", "processed", "01_preprocessing", "impute", sprintf("20260601_pgmatrix_imputed_%s_*.xlsx", ds))))
    new <- readxl::read_excel(Sys.glob(path_processed("01_preprocessing", "impute", sprintf("*_pgmatrix_imputed_%s_*.xlsx", ds)))[1])
    num <- names(new)[vapply(new, is.numeric, TRUE)]
    tibble(dataset = ds, sample_id = num, offset = vapply(num, function(s) median(new[[s]] - old[[s]]), numeric(1)))
  })) |>
    left_join(meta |> mutate(sample_id = as.character(sample_id)), by = "sample_id") |>
    mutate(detected = unname(detected[sample_id]))
  qc <- readr::read_csv(path_processed("01_preprocessing", "impute", "imputation_qc.csv"), show_col_types = FALSE)
  smp <- offsets |>
    transmute(Compartment = unname(DATASETS[dataset]), Animal = AnimalID,
              Group = unname(c(`1` = "CON", `2` = "RES", `3` = "SUS")[as.character(ExpGroup)]),
              Hemisphere = unname(c(L = "left", R = "right")[sub("^.*_(L|R)_.*$", "\\1", sample_id)]),
              Region = region, Layer = if_else(layer %in% c("none", "microglia"), "", layer),
              `Proteins detected` = detected,
              `Centring offset (current minus previous, log2)` = offset,
              `Sample file` = basename(gsub("\\\\", "/", sample_id))) |>
    arrange(Compartment, Group, Animal)

  wb <- createWorkbook()
  add_readme(wb, "Hippocampal spatial proteomics: differential protein abundance", list(
    list(heading = "Design", lines = c(
      "Male mice after social stress, classified as susceptible (SUS) or resilient (RES), and unstressed controls (CON); n = 3 mice per group.",
      "Laser-captured hippocampal regions (CA1, CA2, CA3, dentate gyrus) and layers, profiled in three compartments: neuronal neuropil (10 region-layer units), neuronal soma (4 units) and microglia (4 regions).",
      "The two hemispheres of each mouse were averaged, so every test has one value per mouse.")),
    list(heading = "Processing", lines = c(
      "Protein-group intensities were log2 transformed and each sample was median-centred on the proteins quantified in every sample of its compartment (3,125 neuropil, 4,192 soma, 3,652 microglia proteins).",
      "Proteins missing in more than 70% of a compartment's samples were removed; remaining missing values were imputed per sample from a down-shifted normal distribution (width 0.3, down-shift 1.8 SD; fixed seeds).",
      "Model per compartment: limma, one coefficient per spatial unit and group (~0 + unit:group); comparisons within each unit; robust empirical Bayes moderation; Benjamini-Hochberg FDR within each comparison.")),
    list(heading = "Previous centring (sensitivity analysis)", lines = c(
      "Previously each sample was centred on the median of all proteins it detected. Samples that detected fewer low-abundance proteins then had a higher median, which shifted every protein in that sample down.",
      "Centring on the proteins detected in every sample removes this dependence on detection depth. Because the change is one constant per sample, rows, missing values and imputed values are otherwise identical.",
      "The previous results are shown in the columns marked 'previous centring' and in the per-sample offsets on the sheet 'Sample centring'.")),
    list(heading = "Sheets", lines = c(
      "Overview: number of proteins tested and at FDR < 0.05 per compartment, spatial unit and comparison, now and with the previous centring.",
      "Significant proteins: every protein at FDR < 0.05 in at least one version, with both versions side by side.",
      "All results (one sheet per compartment): every protein, unit and comparison.",
      "Sample centring: per-sample detection depth and centring offset.")),
    list(heading = "Columns", lines = c(
      "log2 fold change: first group minus second group, in log2 units (positive = higher in the first group). Colour scale: blue lower, yellow higher, full colour at |log2 FC| = 1.",
      "t: moderated t statistic. P: moderated P value. FDR: Benjamini-Hochberg adjusted P within the comparison; values below 0.05 are highlighted.",
      "Mean log2 abundance: average centred log2 intensity of the protein in the compartment.",
      "Protein group ID: stable identifier of the protein group used across all analyses."))
  ))
  put_table(wb, "Overview", ov, int_cols = 4:9, left_cols = 1:3, freeze_col = 4, widths = c(20, 12, 13, rep(14, 6)))
  put_table(wb, "Significant proteins", sig, d_cols = col_idx(sig, c("log2 fold change", "log2 fold change, previous centring")),
            p_cols = col_idx(sig, c("P", "FDR", "FDR, previous centring")), sig_cols = col_idx(sig, c("FDR", "FDR, previous centring")),
            n2_cols = col_idx(sig, "t"), left_cols = 1:5, freeze_col = 6, widths = c(18, 11, 12, 12, 13, 12, 9, 10, 10, 14, 13, 16, 34))
  for (ds in names(DATASETS)) {
    tab <- da |> filter(dataset == ds) |>
      transmute(Gene = gene, `UniProt accession` = accession, `Spatial unit` = as.character(spatial_unit), Comparison = as.character(comparison),
                `log2 fold change` = log2fc, `Mean log2 abundance` = aveExpr, t = t, P = pval, FDR = padj,
                `log2 fold change, previous centring` = prev_log2fc, `FDR, previous centring` = prev_padj, `Protein group ID` = ProteinGroupID)
    put_table(wb, paste("All results -", c(neuron_neuropil = "neuropil", neuron_soma = "soma", microglia = "microglia")[[ds]]), tab,
              d_cols = col_idx(tab, c("log2 fold change", "log2 fold change, previous centring")),
              p_cols = col_idx(tab, c("P", "FDR", "FDR, previous centring")), sig_cols = col_idx(tab, c("FDR", "FDR, previous centring")),
              n2_cols = col_idx(tab, c("Mean log2 abundance", "t")), left_cols = 1:2, freeze_col = 5,
              widths = c(13, 13, 11, 12, 12, 12, 9, 10, 10, 14, 13, 34))
  }
  put_table(wb, "Sample centring", smp, n2_cols = col_idx(smp, "Centring offset (current minus previous, log2)"),
            int_cols = col_idx(smp, "Proteins detected"), left_cols = c(1, 9), freeze_col = 3, widths = c(18, 10, 9, 11, 9, 9, 12, 16, 70))
  f <- file.path(OUT, "Proteomics_differential_abundance.xlsx")
  saveWorkbook(wb, f, overwrite = TRUE)
  cat("  wrote", f, "\n"); print(as.data.frame(overview |> filter(fdr_sig > 0 | prev_fdr_sig > 0)), row.names = FALSE)
}

# ================================================================== GO gene-set enrichment (GSEA)
# Rows of the clusterProfiler manifests: GSEA on GO biological process, full (unsimplified) term
# tables, one per within-unit comparison. Old manifests record P:// paths of the old checkout; those
# files are read from the archive.
read_gsea <- function(manifest, root_map = NULL) {
  m <- readr::read_csv(manifest, show_col_types = FALSE, progress = FALSE) |>
    filter(.data$result_type == "GSEA_GO", .data$ontology == "BP", !as.logical(.data$simplified))
  bind_rows(lapply(seq_len(nrow(m)), function(i) {
    f <- m$output_table[[i]]
    if (!is.null(root_map)) f <- sub("^P:/+", paste0(root_map, "/"), f)
    if (is.na(f) || !file.exists(f)) stop("GSEA table listed in ", manifest, " not readable: ", f, call. = FALSE)
    x <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
    bind_cols(parse_contrast(m$contrast[[i]], m$dataset[[i]])[rep(1, nrow(x)), ],
              transmute(x, ID, Description, setSize, NES, pvalue, p.adjust, core_enrichment))
  }))
}
if ("gsea" %in% SECTIONS) {
  cat("GSEA workbook\n")
  gs <- bind_rows(lapply(names(DATASETS), function(ds) {
    cur <- read_gsea(path_results("differential_abundance", "run_clusterprofiler_enrichment", ds, "models", "clusterProfiler_manifest.csv"))
    prev <- read_gsea(file.path(ARCHIVE, "data", "processed", "04_differential_expression_enrichment", "clusterProfiler", ds, "clusterProfiler_manifest.csv"),
                      root_map = ARCHIVE) |>
      select(spatial_unit, comparison, ID, prev_NES = NES, prev_padj = p.adjust)
    full_join(cur, prev, by = c("spatial_unit", "comparison", "ID")) |> mutate(dataset = ds, .before = 1)
  })) |>
    mutate(compartment = factor(DATASETS[dataset], levels = DATASETS),
           spatial_unit = factor(spatial_unit, levels = UNIT_ORDER), comparison = factor(comparison, levels = CMP_ORDER)) |>
    arrange(compartment, spatial_unit, comparison, p.adjust)
  gov <- gs |> group_by(compartment, spatial_unit, comparison) |>
    summarise(`Terms tested` = sum(!is.na(p.adjust)), `FDR < 0.05` = sum(p.adjust < FDR_CUT, na.rm = TRUE),
              `Enriched in first group` = sum(p.adjust < FDR_CUT & NES > 0, na.rm = TRUE),
              `Enriched in second group` = sum(p.adjust < FDR_CUT & NES < 0, na.rm = TRUE),
              `FDR < 0.05, previous centring` = sum(prev_padj < FDR_CUT, na.rm = TRUE),
              `FDR < 0.05 in both` = sum(p.adjust < FDR_CUT & prev_padj < FDR_CUT, na.rm = TRUE), .groups = "drop") |>
    mutate(across(c(compartment, spatial_unit, comparison), as.character)) |>
    rename(Compartment = compartment, `Spatial unit` = spatial_unit, Comparison = comparison)
  lead <- function(x, n = 25) vapply(strsplit(coalesce(x, ""), "/", fixed = TRUE), function(v) {
    v <- v[nzchar(v)]; if (length(v) > n) paste0(paste(v[seq_len(n)], collapse = ", "), ", ... (", length(v), " in total)") else paste(v, collapse = ", ")
  }, character(1))
  gsig <- gs |> filter(p.adjust < FDR_CUT | prev_padj < FDR_CUT) |>
    transmute(Compartment = as.character(compartment), `Spatial unit` = as.character(spatial_unit), Comparison = as.character(comparison),
              `GO term` = Description, `GO ID` = ID, `Genes in set` = setSize, NES = NES, P = pvalue, FDR = p.adjust,
              `NES, previous centring` = prev_NES, `FDR, previous centring` = prev_padj,
              `Called at FDR < 0.05` = case_when(p.adjust < FDR_CUT & prev_padj < FDR_CUT ~ "both",
                                                 p.adjust < FDR_CUT ~ "current only", TRUE ~ "previous centring only"),
              `Leading-edge genes` = lead(core_enrichment))
  wb <- createWorkbook()
  add_readme(wb, "Hippocampal spatial proteomics: GO gene-set enrichment", list(
    list(heading = "Method", lines = c(
      "Gene-set enrichment analysis (clusterProfiler gseGO) on GO biological process for every within-unit comparison of the differential abundance analysis (same design, normalisation and model).",
      "Genes are ranked by the moderated t statistic of the comparison (median over the protein groups of a gene; protein groups that map to several genes are left out). Gene sets of 10 to 800 genes; fgsea with 100,000 simple permutations and multilevel estimation of small P values, with a fixed seed per comparison; Benjamini-Hochberg FDR within each comparison.",
      "NES: normalised enrichment score. Positive = the term's proteins rank towards higher abundance in the first group of the comparison. Colour scale: blue negative, yellow positive, full colour at |NES| = 2.")),
    list(heading = "Previous centring (sensitivity analysis)", lines = c(
      "The columns marked 'previous centring' repeat the result with the earlier per-sample median centring on all detected proteins, which shifted samples with lower detection depth downwards.",
      "A term that is significant only with the previous centring depended on that shift.")),
    list(heading = "Sheets", lines = c(
      "Overview: terms tested and at FDR < 0.05 per compartment, spatial unit and comparison.",
      "Significant terms: every term at FDR < 0.05 in at least one version, with both versions side by side and the leading-edge genes of the current analysis.",
      "All terms (one sheet per compartment): every tested term."))
  ))
  put_table(wb, "Overview", gov, int_cols = 4:9, left_cols = 1:3, freeze_col = 4, widths = c(20, 12, 13, rep(14, 6)))
  put_table(wb, "Significant terms", gsig, d_cols = col_idx(gsig, c("NES", "NES, previous centring")),
            p_cols = col_idx(gsig, c("P", "FDR", "FDR, previous centring")), sig_cols = col_idx(gsig, c("FDR", "FDR, previous centring")),
            int_cols = col_idx(gsig, "Genes in set"), left_cols = c(1:5, 13), freeze_col = 5, lim = 2,
            widths = c(18, 11, 12, 46, 12, 10, 9, 10, 10, 13, 13, 16, 90))
  for (ds in names(DATASETS)) {
    tab <- gs |> filter(dataset == ds) |>
      transmute(`GO term` = Description, `GO ID` = ID, `Spatial unit` = as.character(spatial_unit), Comparison = as.character(comparison),
                `Genes in set` = setSize, NES = NES, P = pvalue, FDR = p.adjust,
                `NES, previous centring` = prev_NES, `FDR, previous centring` = prev_padj)
    put_table(wb, paste("All terms -", c(neuron_neuropil = "neuropil", neuron_soma = "soma", microglia = "microglia")[[ds]]), tab,
              d_cols = col_idx(tab, c("NES", "NES, previous centring")), p_cols = col_idx(tab, c("P", "FDR", "FDR, previous centring")),
              sig_cols = col_idx(tab, c("FDR", "FDR, previous centring")), int_cols = col_idx(tab, "Genes in set"),
              left_cols = 1:2, freeze_col = 3, lim = 2, widths = c(46, 12, 11, 12, 10, 9, 10, 10, 13, 13))
  }
  f <- file.path(OUT, "Proteomics_GO_gene_set_enrichment.xlsx")
  saveWorkbook(wb, f, overwrite = TRUE)
  cat("  wrote", f, "\n"); print(as.data.frame(filter(gov, `FDR < 0.05` > 0 | `FDR < 0.05, previous centring` > 0)), row.names = FALSE)
}

# ================================================================== WGCNA network review
# Reads the attribution (tools/attribute_wgcna_change.R) and lineage (tools/map_wgcna_module_lineage.R)
# tables and the rebuilt Stage 01 outputs. The checklist has empty decision columns for review.
if ("wgcna" %in% SECTIONS) {
  cat("WGCNA review workbook\n")
  rdq <- function(...) { p <- file.path(...); if (file.exists(p)) readr::read_csv(p, show_col_types = FALSE, progress = FALSE, guess_max = 1e5) else NULL }
  att <- function(ds, f) rdq(path_results("wgcna", "attribute_wgcna_change", ds, "tables", f))
  lin <- function(ds, f) rdq(path_results("wgcna", "map_wgcna_module_lineage", ds, "tables", f))
  newt <- function(ds, ...) rdq(path_results("wgcna", "build_wgcna_modules", ds, "tables", ...))
  oldt <- function(ds, ...) rdq(ARCHIVE, "results", "tables", "06_modules_WGCNA", ...)
  have <- names(DATASETS)[vapply(names(DATASETS), function(ds) !is.null(lin(ds, "old_module_correspondence.csv")), logical(1))]
  lab <- function(ds) unname(DATASETS[ds])
  top_go <- function(go, n = 5) {
    if (is.null(go)) return(tibble(ModuleID = character(), top_GO_terms = character()))
    go |> filter(ModuleProteinSetType == "all", p.adjust < FDR_CUT) |> arrange(p.adjust) |>
      group_by(ModuleID) |> summarise(top_GO_terms = paste(head(paste0(Description, " (", Ontology, ")"), n), collapse = "; "), .groups = "drop")
  }
  old_labels <- function(ds) {
    rev <- rdq(repo_path("config", "wgcna_labels", "reference_pre_renormalisation_2026-10-19", paste0(ds, ".csv")))
    auto <- oldt(ds, "interpretable_summary", ds, "WGCNA_final_label_lookup.csv")
    a <- if (is.null(auto)) tibble(level = character(), entity_id = character(), old_label = character()) else
      transmute(auto, level, entity_id, old_label = coalesce(if ("canonical_biological_label" %in% names(auto)) canonical_biological_label else NA_character_,
                                                             if ("best_data_driven_label" %in% names(auto)) best_data_driven_label else NA_character_))
    r <- if (is.null(rev)) tibble(level = character(), entity_id = character(), reviewed = character()) else
      transmute(rev, level, entity_id = as.character(entity_id), reviewed = reviewed_biological_label)
    full_join(a, r, by = c("level", "entity_id")) |>
      transmute(level, entity_id, old_label = coalesce(reviewed, old_label), old_label_status = if_else(!is.na(reviewed), "reviewed", "automatic"))
  }
  new_labels <- function(ds) {
    reg <- lin(ds, "provisional_label_registry.csv")
    if (!is.null(reg)) return(transmute(reg, level, entity_id, provisional_label = reviewed_biological_label,
                                        label_status = if_else(grepl("^Carried from", rationale) | grepl("Carried from", rationale), "carried from predecessor (provisional)", "provisionally unannotated")))
    auto <- rdq(path_results("wgcna", "summarize_module_interpretation", ds, "tables", "WGCNA_final_label_lookup.csv"))
    if (is.null(auto)) return(tibble(level = character(), entity_id = character(), provisional_label = character(), label_status = character()))
    transmute(auto, level, entity_id, provisional_label = best_data_driven_label, label_status = "automatic (from the new GO evidence; not reviewed)")
  }

  steps <- bind_rows(lapply(have, function(ds) att(ds, "attribution_steps.csv")))
  checks <- bind_rows(lapply(have, function(ds) left_join(att(ds, "inclusion_and_offset_check.csv"), att(ds, "correlation_change_summary.csv"), by = "dataset")))
  corr <- bind_rows(lapply(have, function(ds) {
    oc <- lin(ds, "old_module_correspondence.csv"); ol <- old_labels(ds) |> filter(level == "module")
    sz <- oldt(ds, "01_WGCNA", ds, "modules", "WGCNA_modules_long.csv") |> count(ModuleID, name = "n")
    oc |> left_join(ol, by = c("old_module" = "entity_id")) |> left_join(sz, by = c("old_module" = "ModuleID")) |>
      transmute(Compartment = lab(ds), `Previous module` = old_module, `Previous label` = old_label, `Label status` = old_label_status,
                Proteins = n, Class = class, `Best new module` = best_new, Jaccard = best_jaccard,
                `Share in best new module` = share_of_old_in_best, `Reciprocal best match` = reciprocal_best,
                `Where its proteins went` = targets, `Share unassigned` = share_unassigned,
                `Top-10 GO terms significant again` = paste0(n_old_top_still_significant, " of ", n_old_top_terms),
                `GO evidence supports previous label` = go_supports_old,
                `Previous top GO terms` = old_top_terms, `Top GO terms of best new module` = new_top_terms)
  }))
  origin <- bind_rows(lapply(have, function(ds) {
    no <- lin(ds, "new_module_origin.csv"); nl <- new_labels(ds) |> filter(level == "module")
    sz <- newt(ds, "modules", "WGCNA_modules_long.csv") |> count(ModuleID, name = "n")
    sm <- newt(ds, "supermodules", "wgcna_module_supermodule_annotation.csv") |> distinct(ModuleID, SupermoduleID)
    tg <- top_go(newt(ds, "modules", "WGCNA_module_GO_enrichment_long.csv"))
    no |> left_join(sz, by = c("new_module" = "ModuleID")) |> left_join(sm, by = c("new_module" = "ModuleID")) |>
      left_join(nl, by = c("new_module" = "entity_id")) |> left_join(tg, by = c("new_module" = "ModuleID")) |>
      transmute(Compartment = lab(ds), `New module` = new_module, Supermodule = SupermoduleID, Proteins = n, Origin = class,
                `Main predecessor` = main_old, `Share from main predecessor` = share_from_main_old,
                `Composition (previous modules)` = sources, `Provisional label` = provisional_label, `Label status` = label_status,
                `Top GO terms (new evidence)` = top_GO_terms)
  }))
  smt <- bind_rows(lapply(have, function(ds) {
    sl <- lin(ds, "supermodule_lineage.csv"); ol <- old_labels(ds) |> filter(level == "supermodule")
    sl |> left_join(ol, by = c("old_supermodule" = "entity_id")) |>
      transmute(Compartment = lab(ds), `Previous supermodule` = old_supermodule, `Previous label` = old_label, Class = class,
                `Best new supermodule` = best_new, Jaccard = best_jaccard, `Share in best new supermodule` = share_of_old_in_best,
                `Where its proteins went` = targets, `GO evidence supports previous label` = go_supports_old,
                `Previous top GO terms` = old_top_terms, `Top GO terms of best new supermodule` = new_top_terms)
  }))
  carry <- bind_rows(lapply(have, function(ds) { x <- lin(ds, "label_carry_decisions.csv"); if (is.null(x)) NULL else mutate(x, Compartment = lab(ds)) }))
  carry_t <- if (nrow(carry)) carry |> transmute(Compartment, Level = level, `Previous entity` = old_entity_id, `Reviewed label` = reviewed_biological_label,
                                                 Class = class, `New entity` = new_entity_id, Decision = decision, Reason = reason) else tibble()
  themes <- bind_rows(lapply(have, function(ds) { x <- lin(ds, "theme_mapping_carry.csv"); if (is.null(x)) NULL else mutate(x, Compartment = lab(ds)) }))
  themes_t <- if (nrow(themes)) themes |> transmute(Compartment, Theme = theme_id, Level = entity_level, `Previous entity` = old_entity_id,
                                                    `Previous status` = mapping_status, Class = class, `New entity` = new_entity_id,
                                                    `Theme supported by new GO evidence` = theme_supported_in_new, Carried = carry, Reason = reason) else tibble()
  checklist <- bind_rows(
    origin |> transmute(Compartment, Level = "module", Entity = `New module`, Proteins, Origin, `Main predecessor`,
                        `Provisional label`, `Label status`, `Top GO terms (new evidence)`),
    bind_rows(lapply(have, function(ds) {
      nl <- new_labels(ds) |> filter(level == "supermodule")
      sm <- newt(ds, "supermodules", "wgcna_module_supermodule_annotation.csv") |> distinct(ModuleID, SupermoduleID)
      sl <- lin(ds, "supermodule_lineage.csv")
      sizes <- newt(ds, "modules", "WGCNA_modules_long.csv") |> count(ModuleID, name = "n") |>
        inner_join(sm, by = "ModuleID") |> group_by(SupermoduleID) |> summarise(Proteins = sum(n), .groups = "drop")
      tibble(Compartment = lab(ds), Level = "supermodule", Entity = sort(unique(sm$SupermoduleID))) |>
        left_join(sizes, by = c("Entity" = "SupermoduleID")) |>
        mutate(Origin = paste("member modules:", vapply(Entity, function(e) paste(sort(sm$ModuleID[sm$SupermoduleID == e]), collapse = ", "), character(1))),
               `Main predecessor` = vapply(Entity, function(e) paste(sl$old_supermodule[sl$best_new %in% e], collapse = ", "), character(1))) |>
        left_join(transmute(nl, Entity = entity_id, `Provisional label` = provisional_label, `Label status` = label_status), by = "Entity") |>
        mutate(`Top GO terms (new evidence)` = NA_character_)
    }))
  ) |> mutate(Decision = "", `Reviewed label` = "", Notes = "")

  wb <- createWorkbook()
  add_readme(wb, "Hippocampal spatial proteomics: co-abundance network (WGCNA) after re-normalisation", list(
    list(heading = "What changed", lines = c(
      "The networks were rebuilt from the re-normalised protein matrices with the same settings (signed network, biweight midcorrelation, topological overlap, dynamic tree cut with deepSplit 2 and minimum module size 30, merge height 0.25; soft-thresholding power chosen by the same rule, scale-free fit R2 >= 0.8).",
      "The re-normalised matrix equals the previous one plus one constant per sample, so the proteins, the missing values and the imputed values are unchanged; the network changes only through that offset (sheet 'Attribution').",
      "The new modules are the reference. Labels of the previous network are not transferred by module number: numbers were re-assigned and need not denote the same proteins.")),
    list(heading = "Attribution", lines = c(
      "Steps A to E rebuild the network to separate the sources of change: A previous data and settings (reproduces the previous network); B new data, previous samples and power (effect of the offset alone); C + samples retained by the sample-clustering rule; D + the rule-selected power (reproduces the new network); E previous data at the new power (power alone).",
      "Agreement is the adjusted Rand index over shared proteins (1 = identical partition, 0 = chance).")),
    list(heading = "Correspondence between previous and new modules", lines = c(
      "Proteins are matched by protein group. Jaccard = shared proteins / union. For each previous module the best new module is the one with the highest Jaccard.",
      "retained: the two are each other's best match and Jaccard >= 0.5. merged: >= 50% of the previous module sits in one new module that also holds >= 50% of another previous module. changed: >= 50% in one new module otherwise. split: no new module holds 50%, two or more hold >= 20% each. dissolved: otherwise.",
      "Supermodules are compared the same way on the proteins of their member modules.")),
    list(heading = "Labels", lines = c(
      "A previous label is carried only when the module is retained AND at least half of its ten most significant GO terms (biological process, cellular component, molecular function; BH < 0.05) are significant again in the new module. Overlap alone never carries a label.",
      "All other modules are provisionally unannotated and flagged for review. Every label in this workbook is provisional until reviewed; carried labels are marked as such.",
      "Theme-to-module mappings used for the gene-set/module concordance are carried under the same rule, with support defined as a significant GO biological-process term of the theme in the new module.")),
    list(heading = "Sheets", lines = c(
      "Attribution, Network checks: why the modules changed. Previous modules, New modules, Supermodules: the correspondence. Label decisions, Theme mappings: what was carried and why.",
      "Review checklist: one row per new module and supermodule with its provisional label and new GO evidence; the columns Decision (approve / edit / reject), Reviewed label and Notes are for the reviewer."))
  ))
  if (nrow(steps)) put_table(wb, "Attribution", steps |> transmute(Compartment = lab(dataset), Step = step, Description = description, Power = power,
                                                                   Samples = n_samples, Proteins = n_proteins, Modules = n_modules,
                                                                   `ARI vs previous network` = ari_vs_old_network, `ARI vs new network` = ari_vs_new_network),
                             n2_cols = 8:9, int_cols = 4:7, left_cols = 1:3, freeze_col = 4, widths = c(18, 7, 60, 8, 9, 9, 9, 13, 13))
  if (nrow(checks)) put_table(wb, "Network checks", checks |> transmute(Compartment = lab(dataset), `Proteins (previous / new / shared)` = paste(proteins_old, proteins_new, proteins_shared, sep = " / "),
                                                                        `Samples (previous / new / shared)` = paste(samples_old, samples_new, samples_shared, sep = " / "),
                                                                        `Samples only in new network` = samples_only_new,
                                                                        `Largest deviation from a per-sample constant` = max_residual_after_removing_offset,
                                                                        `Power (previous / new)` = paste(old_power, new_power, sep = " / "),
                                                                        `Median change in protein-protein correlation` = median_abs_change,
                                                                        `95th percentile of change` = p95_abs_change,
                                                                        `Correlation of previous and new correlations` = r_old_vs_new),
                              n2_cols = 7:9, left_cols = 1:4, freeze_col = 2, widths = c(18, 18, 18, 40, 16, 12, 16, 14, 16))
  put_table(wb, "Previous modules", corr, n2_cols = col_idx(corr, c("Jaccard", "Share in best new module", "Share unassigned")),
            int_cols = col_idx(corr, "Proteins"), left_cols = c(1:3, 11, 15, 16), freeze_col = 3,
            widths = c(18, 12, 34, 11, 9, 11, 12, 9, 11, 11, 46, 10, 12, 12, 60, 60))
  put_table(wb, "New modules", origin, n2_cols = col_idx(origin, "Share from main predecessor"), int_cols = col_idx(origin, "Proteins"),
            left_cols = c(1:2, 8:11), freeze_col = 3, widths = c(18, 12, 11, 9, 18, 13, 12, 46, 34, 26, 70))
  put_table(wb, "Supermodules", smt, n2_cols = col_idx(smt, c("Jaccard", "Share in best new supermodule")), left_cols = c(1:3, 8, 10, 11),
            freeze_col = 3, widths = c(18, 12, 34, 11, 12, 9, 11, 40, 12, 60, 60))
  if (nrow(carry_t)) put_table(wb, "Label decisions", carry_t, left_cols = c(1, 4, 8), freeze_col = 4, widths = c(18, 11, 12, 46, 11, 11, 24, 70))
  if (nrow(themes_t)) put_table(wb, "Theme mappings", themes_t, left_cols = c(1:2, 10), freeze_col = 3, widths = c(18, 32, 11, 12, 20, 11, 12, 14, 9, 60))
  put_table(wb, "Review checklist", checklist, int_cols = col_idx(checklist, "Proteins"), left_cols = c(1, 5:9, 11:12), freeze_col = 4,
            widths = c(18, 11, 11, 9, 30, 16, 34, 30, 70, 12, 34, 40))
  f <- file.path(OUT, "Proteomics_WGCNA_network_review.xlsx")
  saveWorkbook(wb, f, overwrite = TRUE)
  cat("  wrote", f, "for", paste(have, collapse = ", "), "\n")
}

# ================================================================== conclusions review
# One row per proteomics-related statement of the manuscript claims register, re-evaluated against
# the current results (input: --claims <csv>, assembled outside the pipeline from the register).
if ("conclusions" %in% SECTIONS) {
  cat("Conclusions review workbook\n")
  cl <- readr::read_csv(arg("--claims", "C:/Users/topohl/Documents/proteomics_renorm_refreeze_2026-10-19/conclusions_review.csv"),
                        show_col_types = FALSE, progress = FALSE)
  status_label <- c(
    supported_unchanged = "supported, unchanged", supported_new_value = "supported, new value",
    changed_scope_or_direction = "changed (scope or direction)", no_longer_supported = "no longer supported",
    not_affected_by_renormalisation = "not affected", not_evaluable = "not yet evaluable")
  tab <- cl |> transmute(`Claim` = claim_id, `Where in the text` = location, `Statement` = claim,
                         `Previous value` = previous_value, `Current value` = current_value,
                         `Status` = unname(status_label[status]), `Note` = note,
                         `Suggested wording` = suggested_wording, `Evidence` = evidence)
  summ <- tab |> count(Status, name = "Statements") |> arrange(desc(Statements))
  wb <- createWorkbook()
  add_readme(wb, "Hippocampal spatial proteomics: statements re-checked after re-normalisation", list(
    list(heading = "What this is", lines = c(
      "Every proteomics-related statement of the manuscript claims register, re-checked against the results computed with the current normalisation (each sample median-centred on the proteins quantified in all samples of its compartment).",
      "Previous value: as recorded in the register. Current value: recomputed from the current results with the same definition; where a quantity was recomputed for both versions with one script, that script reproduced the previous value exactly.")),
    list(heading = "Status", lines = c(
      "supported, unchanged: the statement and its number hold.",
      "supported, new value: the statement holds; the number in the text changes (see Current value and Suggested wording).",
      "changed (scope or direction): the statement needs to be narrowed or reworded.",
      "no longer supported: the statement should be removed.",
      "not affected: the statement does not depend on the proteomics quantification (behaviour, design, methods, text).",
      "not yet evaluable: the producing step is being rebuilt; the statement is re-checked when it exists.")),
    list(heading = "Sheets", lines = c("Summary: number of statements per status.", "Statements: one row per statement."))
  ))
  put_table(wb, "Summary", summ, int_cols = 2, freeze_col = 2, widths = c(34, 14))
  put_table(wb, "Statements", tab, left_cols = c(2:5, 7:9), freeze_col = 2, widths = c(9, 16, 50, 46, 56, 24, 60, 56, 50))
  bodywrap <- createStyle(wrapText = TRUE, valign = "top")
  addStyle(wb, "Statements", bodywrap, rows = 2:(nrow(tab) + 1), cols = seq_len(ncol(tab)), gridExpand = TRUE, stack = TRUE)
  f <- file.path(OUT, "Proteomics_statements_review.xlsx")
  saveWorkbook(wb, f, overwrite = TRUE)
  cat("  wrote", f, "\n"); print(as.data.frame(summ), row.names = FALSE)
}
