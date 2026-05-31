# Dataset-aware rank-abundance and marker abundance sanity checks.
# Marker panels are abundance/compartment checks, not definitive purity estimates.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
SUBSTEP_ID <- "04_marker_rank_abundance_qc"
PATHS <- qc_paths(SUBSTEP_ID, DATASET)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_RANK_ABUNDANCE_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_RANK_ABUNDANCE_METADATA_FILE")

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/04_marker_rank_abundance_qc.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c("Writes rank tables, marker score tables, and SVG marker sanity-check plots.")
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "ggrepel", "scales", "svglite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))

if (!file.exists(matrix_file)) {
  stop("Rank-abundance matrix not found for dataset '", DATASET, "': ", matrix_file,
       ". Set PROTEOMICS_RANK_ABUNDANCE_MATRIX_FILE or PROTEOMICS_QC_MATRIX_FILE.", call. = FALSE)
}

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
mat <- expr$mat
meta <- expr$meta

marker_sets <- list(
  neuronal_synaptic_neuropil = c("Stxbp1", "Gpm6a", "Nptn", "Sh3gl2", "Atp6v1g2", "Snap25", "Syn1", "Syp", "Dlg4", "Camk2a"),
  nuclear_soma = c("H2ac1", "H4c1", "H3-3a", "H1-4", "H1-3", "Matr3", "Srsf3", "Ddx39b"),
  microglia = c("Aif1", "Tmem119", "P2ry12", "Cx3cr1", "Csf1r", "C1qa", "C1qb", "Hexb", "Mertk"),
  astrocyte = c("Gfap", "Aqp4", "Aldh1l1", "Slc1a2", "Slc1a3", "Aldoc"),
  oligodendrocyte_myelin = c("Mbp", "Mog", "Plp1", "Cnp", "Mag", "Mobp"),
  endothelial_pericyte_vascular = c("Pecam1", "Cldn5", "Kdr", "Rgs5", "Pdgfrb", "Vtn"),
  mitochondrial_oxphos = c("Ndufs1", "Ndufa9", "Sdha", "Uqcrc2", "Cox4i1", "Atp5f1a", "Atp5f1b"),
  ribosomal_translation = c("Rpl3", "Rpl4", "Rpl5", "Rps3", "Rps6", "Eef1a1", "Eef2"),
  rnp_rna_processing = c("Hnrnpa2b1", "Hnrnpc", "Sfpq", "Snrnp70", "Ddx5", "Ddx17", "Pabpc1")
)

gene_norm <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
protein_ids <- rownames(mat)
protein_key <- gene_norm(sub("_MOUSE$", "", protein_ids, ignore.case = TRUE))

sample_scores <- lapply(names(marker_sets), function(panel) {
  genes <- marker_sets[[panel]]
  idx <- which(protein_key %in% gene_norm(genes))
  if (!length(idx)) {
    return(data.frame(Sample = colnames(mat), marker_panel = panel, n_markers_detected = 0L,
                      marker_score = NA_real_, stringsAsFactors = FALSE))
  }
  data.frame(
    Sample = colnames(mat),
    marker_panel = panel,
    n_markers_detected = length(idx),
    marker_score = colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
sample_scores <- dplyr::bind_rows(sample_scores) |>
  dplyr::left_join(meta, by = "Sample")

rank_group_cols <- intersect(c("Region", "region", "Layer", "layer", "Group", "group", "ExpGroup", "plate"), names(meta))
if (!length(rank_group_cols)) rank_group_cols <- "Sample"
rank_group_cols <- unique(c(intersect(c("Region", "region", "Layer", "layer"), rank_group_cols), rank_group_cols[1]))

long <- as.data.frame(mat, check.names = FALSE) |>
  tibble::rownames_to_column("Protein") |>
  tidyr::pivot_longer(-Protein, names_to = "Sample", values_to = "Log2Intensity") |>
  dplyr::left_join(meta, by = "Sample")

rank_data <- long |>
  dplyr::mutate(RankGroup = do.call(paste, c(dplyr::across(dplyr::all_of(rank_group_cols)), sep = " | "))) |>
  dplyr::filter(!is.na(Log2Intensity), !is.na(RankGroup), nzchar(RankGroup)) |>
  dplyr::group_by(RankGroup, Protein) |>
  dplyr::summarise(MeanLog2 = mean(Log2Intensity, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(LinearValue = 2^MeanLog2) |>
  dplyr::group_by(RankGroup) |>
  dplyr::arrange(dplyr::desc(LinearValue), .by_group = TRUE) |>
  dplyr::mutate(Rank = dplyr::row_number()) |>
  dplyr::ungroup()

marker_lookup <- dplyr::bind_rows(lapply(names(marker_sets), function(panel) {
  data.frame(marker_panel = panel, marker = marker_sets[[panel]], marker_key = gene_norm(marker_sets[[panel]]))
}))
rank_data <- rank_data |>
  dplyr::mutate(marker_key = gene_norm(sub("_MOUSE$", "", Protein, ignore.case = TRUE))) |>
  dplyr::left_join(marker_lookup, by = "marker_key") |>
  dplyr::mutate(marker_panel = ifelse(is.na(marker_panel), "none", marker_panel))

qc_write_csv(rank_data, file.path(PATHS$tables, "rank_abundance_table.csv"))
qc_write_csv(sample_scores, file.path(PATHS$tables, "marker_scores_by_sample.csv"))
qc_write_xlsx(list(rank_abundance = rank_data, marker_scores = sample_scores),
              file.path(PATHS$tables, "rank_abundance_marker_qc.xlsx"))

plot_data <- rank_data |>
  dplyr::mutate(is_marker = marker_panel != "none")
label_data <- plot_data |>
  dplyr::filter(is_marker) |>
  dplyr::group_by(RankGroup, marker_panel) |>
  dplyr::slice_min(Rank, n = 3, with_ties = FALSE) |>
  dplyr::ungroup()

p_rank <- ggplot(plot_data, aes(Rank, LinearValue)) +
  geom_point(color = "grey78", alpha = 0.08, size = 0.2) +
  geom_point(data = label_data, aes(color = marker_panel), size = 1.2) +
  ggrepel::geom_label_repel(
    data = label_data,
    aes(label = Protein, fill = marker_panel),
    color = "white", size = 2, label.size = 0, max.overlaps = 100
  ) +
  scale_y_log10(labels = scales::label_number()) +
  facet_wrap(~RankGroup, scales = "free_x") +
  labs(x = "Protein rank", y = "Intensity (log10)", color = "Marker panel", fill = "Marker panel") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave(file.path(PATHS$figures, "rank_abundance_marker_sanity_check.svg"), p_rank,
       width = 180, height = 140, units = "mm", device = svglite::svglite)

summary_vars <- intersect(c("Group", "group", "ExpGroup", "Region", "region", "Layer", "layer", "plate"), names(sample_scores))
summary_vars <- summary_vars[!duplicated(tolower(summary_vars))]
if (length(summary_vars)) {
  score_summary <- sample_scores |>
    tidyr::pivot_longer(dplyr::all_of(summary_vars), names_to = "metadata_term", values_to = "metadata_value") |>
    dplyr::filter(!is.na(metadata_value), nzchar(as.character(metadata_value))) |>
    dplyr::group_by(marker_panel, metadata_term, metadata_value) |>
    dplyr::summarise(n = dplyr::n(), mean_score = mean(marker_score, na.rm = TRUE),
                     median_score = median(marker_score, na.rm = TRUE), .groups = "drop")
  qc_write_csv(score_summary, file.path(PATHS$tables, "marker_score_summary_by_metadata.csv"))

  p_score <- ggplot(sample_scores, aes(x = marker_panel, y = marker_score, color = .data[[summary_vars[[1]]]])) +
    geom_boxplot(outlier.shape = NA, color = "grey35") +
    geom_point(position = position_jitter(width = 0.16), alpha = 0.75, size = 1) +
    coord_flip() +
    labs(x = NULL, y = "Mean marker abundance", color = summary_vars[[1]]) +
    theme_classic(base_size = 8) +
    theme(legend.position = "bottom")
  ggsave(file.path(PATHS$figures, "marker_abundance_summary.svg"), p_score,
         width = 150, height = 100, units = "mm", device = svglite::svglite)
}

write_run_manifest(
  file.path(PATHS$logs, "run_manifest.yml"),
  inputs = list(matrix = matrix_file, metadata = metadata_file),
  outputs = list(figures = PATHS$figures, tables = PATHS$tables),
  parameters = list(dataset = DATASET, marker_sets = names(marker_sets)),
  notes = "Marker abundance/compartment sanity checks only; not cell-type purity estimates."
)

message("Rank-abundance and marker QC complete for dataset: ", DATASET)
