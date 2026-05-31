# Dataset-aware replicate and animal-level consistency QC.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
PATHS <- qc_paths("03_replicate_consistency", DATASET)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_REPLICATE_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_REPLICATE_METADATA_FILE")

if (run$dry_run) {
  status <- qc_dry_run_contract("03_qc_exploration/03_replicate_consistency.r", DATASET,
                                matrix_file = matrix_file, metadata_file = metadata_file, paths = PATHS,
                                extra = c("Requires AnimalID and/or ReplicateGroup for strongest replicate interpretation."))
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "tibble", "ggplot2", "svglite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))
if (!file.exists(matrix_file)) stop("Replicate matrix not found: ", matrix_file, call. = FALSE)

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
mat <- qc_impute_for_pca(expr$mat)
meta <- expr$meta[colnames(mat), , drop = FALSE]
if (nrow(mat) < 2L || ncol(mat) < 2L) stop("Matrix too small for replicate consistency after filtering.", call. = FALSE)

cor_mat <- stats::cor(mat, use = "pairwise.complete.obs", method = "pearson")
pairs <- utils::combn(colnames(mat), 2)
pair_df <- data.frame(Sample1 = pairs[1, ], Sample2 = pairs[2, ],
                      correlation = cor_mat[cbind(pairs[1, ], pairs[2, ])],
                      stringsAsFactors = FALSE)
meta1 <- meta[pair_df$Sample1, , drop = FALSE]
meta2 <- meta[pair_df$Sample2, , drop = FALSE]
pair_df$same_AnimalID <- if ("AnimalID" %in% names(meta)) meta1$AnimalID == meta2$AnimalID else NA
pair_df$same_ReplicateGroup <- if ("ReplicateGroup" %in% names(meta)) meta1$ReplicateGroup == meta2$ReplicateGroup else NA
pair_df$same_plate <- if ("plate" %in% names(meta)) meta1$plate == meta2$plate else NA
for (term in intersect(c("Group", "group", "ExpGroup", "Region", "region", "Layer", "layer"), names(meta))) {
  pair_df[[paste0("same_", term)]] <- meta1[[term]] == meta2[[term]]
}

summary_terms <- grep("^same_", names(pair_df), value = TRUE)
pair_summary <- dplyr::bind_rows(lapply(summary_terms, function(term) {
  pair_df |>
    dplyr::filter(!is.na(.data[[term]])) |>
    dplyr::group_by(status = .data[[term]]) |>
    dplyr::summarise(term = term, n_pairs = dplyr::n(),
                     median_correlation = median(correlation, na.rm = TRUE),
                     mean_correlation = mean(correlation, na.rm = TRUE), .groups = "drop")
}))

animal_agg <- NULL
if ("AnimalID" %in% names(meta)) {
  animals <- unique(meta$AnimalID[!is.na(meta$AnimalID) & nzchar(as.character(meta$AnimalID))])
  if (length(animals) >= 2L) {
    animal_agg <- sapply(animals, function(a) rowMeans(mat[, meta$AnimalID == a, drop = FALSE], na.rm = TRUE))
    animal_agg <- as.data.frame(animal_agg, check.names = FALSE) |>
      tibble::rownames_to_column("Protein")
    qc_write_csv(animal_agg, file.path(PATHS$tables, "animal_aggregated_matrix.csv"))
  }
}

qc_write_csv(pair_df, file.path(PATHS$tables, "pairwise_sample_correlations.csv"))
qc_write_csv(pair_summary, file.path(PATHS$tables, "replicate_correlation_summary.csv"))
qc_write_xlsx(list(pairwise_correlations = pair_df, summary = pair_summary),
              file.path(PATHS$tables, "replicate_consistency.xlsx"))

p <- ggplot(pair_df, aes(x = factor(same_AnimalID), y = correlation)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey30") +
  geom_point(position = position_jitter(width = 0.12), alpha = 0.25, size = 0.6) +
  labs(x = "Same AnimalID", y = "Pairwise sample correlation") +
  theme_classic(base_size = 8)
ggsave(file.path(PATHS$figures, "within_vs_across_animal_correlations.svg"), p,
       width = 100, height = 80, units = "mm", device = svglite::svglite)

flag <- "PASS"
if ("same_AnimalID" %in% names(pair_df) && any(pair_df$same_AnimalID, na.rm = TRUE)) {
  within <- median(pair_df$correlation[pair_df$same_AnimalID %in% TRUE], na.rm = TRUE)
  across <- median(pair_df$correlation[pair_df$same_AnimalID %in% FALSE], na.rm = TRUE)
  if (is.finite(within) && is.finite(across) && within < across) flag <- "WARN"
} else {
  flag <- "WARN"
}
writeLines(c("# Replicate Consistency", "", paste0("Dataset: ", DATASET),
             paste0(flag, ": within-animal replicate structure checked where metadata allowed."),
             "The biological unit is animal; downstream models should avoid treating tissue punches as independent animals."),
           file.path(PATHS$reports, "replicate_consistency_summary.md"))

write_run_manifest(file.path(PATHS$logs, "run_manifest.yml"),
                   inputs = list(matrix = matrix_file, metadata = metadata_file),
                   outputs = list(tables = PATHS$tables, figures = PATHS$figures, reports = PATHS$reports),
                   parameters = list(dataset = DATASET),
                   notes = "Pairwise sample correlation and animal aggregation QC.")

message("Replicate consistency QC complete for dataset: ", DATASET)
