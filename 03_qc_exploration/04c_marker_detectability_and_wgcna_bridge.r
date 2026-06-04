# Marker-panel detectability and optional WGCNA marker-enrichment bridge.
# This script is QC/interpreter support, not a formal purity or deconvolution analysis.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
SUBSTEP_ID <- "04c_marker_detectability_and_wgcna_bridge"
PATHS <- qc_paths(SUBSTEP_ID, DATASET)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_MARKER_DETECTABILITY_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_MARKER_DETECTABILITY_METADATA_FILE")
marker_file <- Sys.getenv(
  "PROTEOMICS_MARKER_PANEL_FILE",
  unset = repo_path("config", "marker_panels", "wgcna_reference_marker_sets.csv")
)

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/04c_marker_detectability_and_wgcna_bridge.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c(
      paste0("Marker panel file: ", marker_file, if (file.exists(marker_file)) " [PASS]" else " [fallback internal panels]"),
      "Writes marker detectability/missingness tables and SVG figures.",
      "Optionally writes WGCNA module x marker-panel enrichment if module assignments exist."
    )
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

if (!file.exists(matrix_file)) {
  stop("Marker-detectability matrix not found for dataset '", DATASET, "': ", matrix_file,
       ". Set PROTEOMICS_MARKER_DETECTABILITY_MATRIX_FILE or PROTEOMICS_QC_MATRIX_FILE.", call. = FALSE)
}

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
mat <- expr$mat
meta <- expr$meta

internal_marker_sets <- function() {
  list(
    neuronal_synaptic_neuropil = c("Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2", "Snap25", "Syn1", "Syp", "Dlg4", "Camk2a"),
    nuclear_soma = c("H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3", "Matr3", "Srsf3", "Ddx39b"),
    microglia_homeostatic = c("Aif1", "Tmem119", "P2ry12", "Cx3cr1", "Csf1r", "Hexb", "Fcrls", "Gpr34", "Mertk"),
    microglia_phagolysosomal = c("Tyrobp", "Trem2", "Apoe", "Lpl", "Cst7", "Ctsb", "Ctsd", "Lgals3", "Itgax", "Axl", "C1qa", "C1qb", "C1qc"),
    astrocyte = c("Gfap", "Aqp4", "Aldh1l1", "Slc1a2", "Slc1a3", "Aldoc"),
    oligodendrocyte_myelin = c("Mbp", "Mog", "Plp1", "Cnp", "Mag", "Mobp"),
    endothelial_pericyte_vascular = c("Pecam1", "Cldn5", "Kdr", "Rgs5", "Pdgfrb", "Vtn"),
    extracellular_matrix_basement_membrane = c("Agrn", "Lamb2", "Lamc1", "Col4a1", "Col4a2", "Nid1", "Nid2", "Lama2", "Lama4", "Lama5", "Itga1", "Itgb1"),
    mitochondrial_oxphos = c("Ndufs1", "Ndufa9", "Sdha", "Uqcrc2", "Cox4i1", "Atp5f1a", "Atp5f1b"),
    ribosomal_translation = c("Rpl3", "Rpl4", "Rpl5", "Rps3", "Rps6", "Eef1a1", "Eef2"),
    rnp_rna_processing = c("Hnrnpa2b1", "Hnrnpc", "Sfpq", "Snrnp70", "Ddx5", "Ddx17", "Pabpc1")
  )
}

gene_norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", sub("_MOUSE$", "", as.character(x), ignore.case = TRUE)))
collapse_unique <- function(x) paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse = ";")

read_marker_sets <- function(path) {
  if (!file.exists(path)) return(internal_marker_sets())
  df <- tryCatch(qc_read_table(path), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(internal_marker_sets())
  panel_col <- qc_first_col(df, c("marker_panel", "panel", "celltype", "cell_type", "signature", "set", "source_panel"))
  gene_col <- qc_first_col(df, c("gene", "gene_symbol", "marker", "symbol", "protein", "protein_id", "Gene"))
  if (is.na(panel_col) || is.na(gene_col)) return(internal_marker_sets())
  df <- df[!is.na(df[[panel_col]]) & nzchar(as.character(df[[panel_col]])) & !is.na(df[[gene_col]]) & nzchar(as.character(df[[gene_col]])), , drop = FALSE]
  split(as.character(df[[gene_col]]), as.character(df[[panel_col]]))
}

marker_sets <- read_marker_sets(marker_file)
protein_ids <- rownames(mat)
protein_key <- gene_norm(protein_ids)
universe <- unique(protein_key[!is.na(protein_key) & nzchar(protein_key)])

long <- as.data.frame(mat, check.names = FALSE) |>
  tibble::rownames_to_column("Protein") |>
  tidyr::pivot_longer(-Protein, names_to = "Sample", values_to = "Log2Intensity") |>
  dplyr::left_join(meta, by = "Sample")

rank_group_cols <- intersect(c("Region", "region", "Layer", "layer", "Group", "group", "ExpGroup", "Sex", "sex", "plate", "batch"), names(meta))
rank_group_cols <- rank_group_cols[!duplicated(tolower(rank_group_cols))]
if (!length(rank_group_cols)) rank_group_cols <- "Sample"
rank_data <- long |>
  dplyr::mutate(RankGroup = do.call(paste, c(dplyr::across(dplyr::all_of(rank_group_cols)), sep = " | "))) |>
  dplyr::filter(!is.na(Log2Intensity), !is.na(RankGroup), nzchar(RankGroup)) |>
  dplyr::group_by(RankGroup, Protein) |>
  dplyr::summarise(MeanLog2 = mean(Log2Intensity, na.rm = TRUE), .groups = "drop") |>
  dplyr::group_by(RankGroup) |>
  dplyr::arrange(dplyr::desc(MeanLog2), .by_group = TRUE) |>
  dplyr::mutate(Rank = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(protein_key = gene_norm(Protein))

marker_lookup <- dplyr::bind_rows(lapply(names(marker_sets), function(panel) {
  data.frame(marker_panel = panel, requested_marker = unique(marker_sets[[panel]]), marker_key = gene_norm(unique(marker_sets[[panel]])), stringsAsFactors = FALSE)
})) |>
  dplyr::filter(!is.na(marker_key), nzchar(marker_key)) |>
  dplyr::distinct(marker_panel, requested_marker, marker_key)

protein_stats <- data.frame(Protein = protein_ids, protein_key = protein_key, stringsAsFactors = FALSE) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_samples = ncol(mat),
    n_nonmissing = sum(is.finite(mat[Protein, ])),
    fraction_nonmissing = n_nonmissing / n_samples,
    mean_log2_abundance = ifelse(n_nonmissing > 0, mean(mat[Protein, ], na.rm = TRUE), NA_real_),
    median_log2_abundance = ifelse(n_nonmissing > 0, median(mat[Protein, ], na.rm = TRUE), NA_real_),
    sd_log2_abundance = ifelse(n_nonmissing > 1, stats::sd(mat[Protein, ], na.rm = TRUE), NA_real_),
    cv_like = ifelse(is.finite(sd_log2_abundance) && is.finite(mean_log2_abundance) && abs(mean_log2_abundance) > .Machine$double.eps,
                     sd_log2_abundance / abs(mean_log2_abundance), NA_real_)
  ) |>
  dplyr::ungroup()

rank_stats <- rank_data |>
  dplyr::group_by(protein_key) |>
  dplyr::summarise(mean_rank = mean(Rank, na.rm = TRUE), median_rank = median(Rank, na.rm = TRUE), .groups = "drop")

marker_detectability_by_protein <- marker_lookup |>
  dplyr::left_join(protein_stats, by = c("marker_key" = "protein_key")) |>
  dplyr::group_by(dataset = DATASET, marker_panel, requested_marker, marker_key) |>
  dplyr::summarise(
    matched_protein_id = collapse_unique(Protein),
    detected = any(!is.na(Protein)),
    n_matched_protein_ids = sum(!is.na(Protein)),
    n_samples = dplyr::first(n_samples[!is.na(n_samples)]) %||% ncol(mat),
    n_nonmissing = ifelse(all(is.na(n_nonmissing)), 0L, sum(n_nonmissing, na.rm = TRUE)),
    fraction_nonmissing = ifelse(n_matched_protein_ids > 0, max(fraction_nonmissing, na.rm = TRUE), 0),
    mean_log2_abundance = ifelse(n_matched_protein_ids > 0, mean(mean_log2_abundance, na.rm = TRUE), NA_real_),
    median_log2_abundance = ifelse(n_matched_protein_ids > 0, median(median_log2_abundance, na.rm = TRUE), NA_real_),
    sd_log2_abundance = ifelse(n_matched_protein_ids > 0, mean(sd_log2_abundance, na.rm = TRUE), NA_real_),
    cv_like = ifelse(n_matched_protein_ids > 0, mean(cv_like, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) |>
  dplyr::left_join(rank_stats, by = c("marker_key" = "protein_key")) |>
  dplyr::mutate(
    matched_protein_id = ifelse(nzchar(matched_protein_id), matched_protein_id, NA_character_),
    fraction_nonmissing = ifelse(is.finite(fraction_nonmissing), fraction_nonmissing, 0),
    mean_log2_abundance = ifelse(is.nan(mean_log2_abundance), NA_real_, mean_log2_abundance),
    median_log2_abundance = ifelse(is.nan(median_log2_abundance), NA_real_, median_log2_abundance),
    sd_log2_abundance = ifelse(is.nan(sd_log2_abundance), NA_real_, sd_log2_abundance),
    cv_like = ifelse(is.nan(cv_like), NA_real_, cv_like)
  )

sample_scores <- lapply(names(marker_sets), function(panel) {
  genes <- gene_norm(marker_sets[[panel]])
  idx <- which(protein_key %in% genes)
  if (!length(idx)) {
    return(data.frame(dataset = DATASET, Sample = colnames(mat), marker_panel = panel, marker_score = NA_real_,
                      n_detected_markers = 0L, fraction_detected_markers = 0, stringsAsFactors = FALSE))
  }
  data.frame(
    dataset = DATASET,
    Sample = colnames(mat),
    marker_panel = panel,
    marker_score = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
    n_detected_markers = length(unique(protein_key[idx])),
    fraction_detected_markers = length(unique(protein_key[idx])) / length(unique(genes)),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::left_join(meta, by = "Sample")

marker_detectability_by_panel <- marker_detectability_by_protein |>
  dplyr::group_by(dataset, marker_panel) |>
  dplyr::summarise(
    n_markers_requested = dplyr::n_distinct(requested_marker),
    n_markers_detected = sum(detected),
    fraction_markers_detected = n_markers_detected / n_markers_requested,
    median_fraction_nonmissing = median(fraction_nonmissing, na.rm = TRUE),
    mean_marker_score = mean(sample_scores$marker_score[sample_scores$marker_panel == dplyr::first(marker_panel)], na.rm = TRUE),
    median_marker_score = median(sample_scores$marker_score[sample_scores$marker_panel == dplyr::first(marker_panel)], na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    mean_marker_score = ifelse(is.nan(mean_marker_score), NA_real_, mean_marker_score),
    median_marker_score = ifelse(is.nan(median_marker_score), NA_real_, median_marker_score)
  )

qc_write_csv(marker_detectability_by_protein, file.path(PATHS$tables, "marker_detectability_by_protein.csv"))
qc_write_csv(marker_detectability_by_panel, file.path(PATHS$tables, "marker_detectability_by_panel.csv"))
qc_write_csv(sample_scores, file.path(PATHS$tables, "marker_detectability_by_sample.csv"))
qc_write_xlsx(
  list(
    marker_detectability_by_sample = sample_scores,
    marker_detectability_by_protein = marker_detectability_by_protein,
    marker_detectability_by_panel = marker_detectability_by_panel
  ),
  file.path(PATHS$tables, "marker_detectability_qc.xlsx")
)

p_detection <- ggplot(marker_detectability_by_panel, aes(x = reorder(marker_panel, fraction_markers_detected), y = fraction_markers_detected)) +
  geom_col(fill = "grey55") +
  geom_text(aes(label = paste0(n_markers_detected, "/", n_markers_requested)), hjust = -0.05, size = 2) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(x = NULL, y = "Detected marker fraction") +
  theme_classic(base_size = 8)
ggsave(file.path(PATHS$figures, "marker_panel_detection_barplot.svg"), p_detection, width = 140, height = 100, units = "mm", device = svglite::svglite)

plot_missing <- marker_detectability_by_protein |>
  dplyr::mutate(marker_label = ifelse(detected, requested_marker, paste0(requested_marker, " (not detected)")))
p_missing <- ggplot(plot_missing, aes(x = marker_panel, y = marker_label, size = fraction_nonmissing, fill = mean_log2_abundance)) +
  geom_point(shape = 21, color = "grey30", alpha = 0.85) +
  scale_size_continuous(range = c(0.4, 3.2), limits = c(0, 1)) +
  labs(x = NULL, y = NULL, size = "Nonmissing fraction", fill = "Mean log2 abundance") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(file.path(PATHS$figures, "marker_missingness_dotplot.svg"), p_missing, width = 180, height = max(90, min(260, 3.2 * nrow(plot_missing))), units = "mm", device = svglite::svglite)

group_col <- intersect(c("Dataset", "dataset", "RegionLayer", "Region", "region", "Layer", "layer", "Group", "group", "ExpGroup", "Sex", "sex", "Batch", "batch", "plate"), names(sample_scores))
group_col <- group_col[!duplicated(tolower(group_col))]
if (length(group_col)) {
  p_box <- ggplot(sample_scores, aes(x = marker_panel, y = marker_score, color = .data[[group_col[[1]]]])) +
    geom_boxplot(outlier.shape = NA, color = "grey35") +
    geom_point(position = position_jitter(width = 0.15), alpha = 0.75, size = 0.9) +
    coord_flip() +
    labs(x = NULL, y = "Mean marker abundance", color = group_col[[1]]) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom")
  ggsave(file.path(PATHS$figures, "marker_score_sample_boxplot.svg"), p_box, width = 150, height = 105, units = "mm", device = svglite::svglite)

  heat <- sample_scores |>
    dplyr::mutate(marker_group = as.character(.data[[group_col[[1]]]])) |>
    dplyr::filter(!is.na(marker_group), nzchar(marker_group)) |>
    dplyr::group_by(marker_panel, marker_group) |>
    dplyr::summarise(median_marker_score = median(marker_score, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(median_marker_score = ifelse(is.nan(median_marker_score), NA_real_, median_marker_score))
  qc_write_csv(heat, file.path(PATHS$tables, "marker_detectability_heatmap_source.csv"))
  p_heat <- ggplot(heat, aes(x = marker_group, y = marker_panel, fill = median_marker_score)) +
    geom_tile(color = "white", linewidth = 0.2) +
    labs(x = group_col[[1]], y = NULL, fill = "Median marker score") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  ggsave(file.path(PATHS$figures, "marker_detectability_heatmap.svg"), p_heat, width = 150, height = 95, units = "mm", device = svglite::svglite)
}

find_wgcna_modules <- function(dataset) {
  override <- Sys.getenv("PROTEOMICS_WGCNA_MODULE_ASSIGNMENT_FILE", unset = "")
  if (nzchar(override)) return(override)
  roots <- c(
    path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules"),
    path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset),
    path_results("tables", "06_modules_WGCNA", dataset),
    path_processed("06_modules_WGCNA", "01_WGCNA", dataset, "modules")
  )
  patterns <- c("module.*assignment.*\\.(csv|tsv)$", "module.*membership.*\\.(csv|tsv)$", "protein.*module.*\\.(csv|tsv)$", ".*modules.*\\.(csv|tsv)$")
  files <- unique(unlist(lapply(roots[dir.exists(roots)], function(root) {
    unique(unlist(lapply(patterns, function(pat) list.files(root, pattern = pat, recursive = TRUE, full.names = TRUE, ignore.case = TRUE))))
  })))
  if (length(files)) files[[1]] else NA_character_
}

read_wgcna_modules <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  df <- tryCatch(qc_read_table(path), error = function(e) NULL)
  if (is.null(df) || !nrow(df)) return(NULL)
  protein_col <- qc_first_col(df, c("Protein", "protein", "ProteinID", "protein_id", "gene", "gene_symbol", "Gene", "feature", "Feature", "id"))
  module_col <- qc_first_col(df, c("module", "Module", "module_label", "module_color", "color", "WGCNA_module"))
  if (is.na(protein_col) || is.na(module_col)) return(NULL)
  df |>
    dplyr::transmute(module = as.character(.data[[module_col]]), Protein = as.character(.data[[protein_col]]), protein_key = gene_norm(.data[[protein_col]])) |>
    dplyr::filter(!is.na(module), nzchar(module), !is.na(protein_key), nzchar(protein_key)) |>
    dplyr::distinct(module, Protein, protein_key)
}

module_file <- find_wgcna_modules(DATASET)
module_df <- read_wgcna_modules(module_file)
if (!is.null(module_df) && nrow(module_df)) {
  marker_panel_keys <- marker_lookup |>
    dplyr::filter(marker_key %in% universe) |>
    dplyr::group_by(marker_panel) |>
    dplyr::summarise(panel_keys = list(unique(marker_key)), .groups = "drop")
  modules <- split(module_df$protein_key, module_df$module)
  enrich <- dplyr::bind_rows(lapply(names(modules), function(mod) {
    mod_keys <- unique(modules[[mod]][modules[[mod]] %in% universe])
    dplyr::bind_rows(lapply(seq_len(nrow(marker_panel_keys)), function(i) {
      panel <- marker_panel_keys$marker_panel[[i]]
      panel_keys <- unique(marker_panel_keys$panel_keys[[i]])
      overlap <- intersect(mod_keys, panel_keys)
      a <- length(overlap)
      b <- length(setdiff(mod_keys, panel_keys))
      c <- length(setdiff(panel_keys, mod_keys))
      d <- length(setdiff(universe, union(mod_keys, panel_keys)))
      ft <- tryCatch(stats::fisher.test(matrix(c(a, b, c, d), nrow = 2)), error = function(e) NULL)
      data.frame(
        dataset = DATASET,
        module = mod,
        marker_panel = panel,
        module_size = length(mod_keys),
        panel_size_detected = length(panel_keys),
        overlap_n = a,
        overlap_genes = paste(sort(overlap), collapse = ";"),
        odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
        p_value = if (is.null(ft)) NA_real_ else ft$p.value,
        stringsAsFactors = FALSE
      )
    }))
  })) |>
    dplyr::group_by(marker_panel) |>
    dplyr::mutate(padj_BH = stats::p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup() |>
    dplyr::mutate(enrichment_status = dplyr::case_when(
      !is.na(padj_BH) & padj_BH < 0.05 ~ "FDR_significant",
      !is.na(p_value) & p_value < 0.05 ~ "nominal",
      TRUE ~ "not_supported"
    ))
  qc_write_csv(enrich, file.path(PATHS$tables, "wgcna_module_marker_panel_enrichment.csv"))

  heat_wgcna <- enrich |>
    dplyr::mutate(score = ifelse(is.finite(padj_BH) & padj_BH > 0, -log10(padj_BH), NA_real_),
                  score = ifelse(is.finite(odds_ratio) & odds_ratio < 1, -score, score))
  p_wgcna <- ggplot(heat_wgcna, aes(x = marker_panel, y = module, fill = score)) +
    geom_tile(color = "white", linewidth = 0.2) +
    labs(x = NULL, y = "WGCNA module", fill = "signed -log10 FDR") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
  ggsave(file.path(PATHS$figures, "wgcna_module_marker_panel_enrichment_heatmap.svg"), p_wgcna,
         width = 170, height = max(90, min(240, 4 * length(unique(heat_wgcna$module)))), units = "mm", device = svglite::svglite)
} else {
  message("No WGCNA module assignment table found for dataset ", DATASET, "; skipping WGCNA marker-panel enrichment bridge.")
}

notes <- c(
  "# Marker detectability and WGCNA marker-panel bridge", "",
  "This QC report evaluates whether canonical and reference marker proteins are detectable and consistently quantified in the active dataset.", "",
  "Interpretation constraints:",
  "- Marker-panel scoring is QC and interpretive support, not formal cell-type purity estimation or deconvolution.",
  "- Low detectability or high missingness of canonical microglia markers can explain why WGCNA modules are driven by better-detected ECM, synaptic, vascular, mitochondrial, ribosomal, or microenvironment-associated proteins.",
  "- Marker abundance and WGCNA module assignment answer different questions: sample identity/abundance versus covariance/module topology.",
  "- Agreement in directionality between WGCNA eigengenes and previous GSEA/leading-edge results is robustness evidence even when module labels differ from marker-panel labels.", "",
  paste0("Dataset: ", DATASET),
  paste0("Expression matrix: ", matrix_file),
  paste0("Marker panel source: ", if (file.exists(marker_file)) marker_file else "internal fallback panels"),
  paste0("WGCNA module assignment source: ", if (!is.na(module_file) && file.exists(module_file)) module_file else "not found")
)
writeLines(notes, file.path(PATHS$reports, "marker_detectability_interpretation_notes.md"))

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(matrix = matrix_file, metadata = metadata_file, marker_panels = marker_file, wgcna_modules = module_file),
  outputs = list(figures = PATHS$figures, tables = PATHS$tables, reports = PATHS$reports),
  parameters = list(dataset = DATASET, marker_sets = names(marker_sets)),
  notes = "Marker detectability/missingness QC plus optional WGCNA module x marker-panel enrichment bridge; not formal purity/deconvolution."
)

message("Marker detectability QC complete for dataset: ", DATASET)
