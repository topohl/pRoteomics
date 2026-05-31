# Dataset-aware missingness diagnostics before imputation when possible.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
PATHS <- qc_paths("02_missingness_diagnostics", DATASET)
matrix_file <- qc_resolve_matrix(DATASET, env = "PROTEOMICS_MISSINGNESS_MATRIX_FILE")
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_MISSINGNESS_METADATA_FILE")

if (run$dry_run) {
  status <- qc_dry_run_contract(
    "03_qc_exploration/02_missingness_diagnostics.r",
    DATASET,
    matrix_file = matrix_file,
    metadata_file = metadata_file,
    paths = PATHS,
    extra = c("Prefer raw/non-imputed matrix via PROTEOMICS_MISSINGNESS_MATRIX_FILE.")
  )
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "svglite", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))
if (!file.exists(matrix_file)) stop("Missingness matrix not found: ", matrix_file, call. = FALSE)

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
mat <- expr$mat
meta <- expr$meta
input_is_imputed <- grepl("imput", basename(matrix_file), ignore.case = TRUE)

sample_missing <- data.frame(Sample = colnames(mat), missing_fraction = colMeans(is.na(mat))) |>
  dplyr::left_join(meta, by = "Sample")
protein_missing <- data.frame(Protein = rownames(mat), missing_fraction = rowMeans(is.na(mat)))
terms <- qc_metadata_terms(meta)

by_term <- dplyr::bind_rows(lapply(terms, function(term) {
  sample_missing |>
    dplyr::filter(!is.na(.data[[term]]), nzchar(as.character(.data[[term]]))) |>
    dplyr::group_by(value = .data[[term]]) |>
    dplyr::summarise(term = term, n = dplyr::n(),
                     mean_missing = mean(missing_fraction, na.rm = TRUE),
                     median_missing = median(missing_fraction, na.rm = TRUE), .groups = "drop")
}))

assoc <- dplyr::bind_rows(lapply(terms, function(term) {
  data.frame(term = term,
             p_value = qc_safe_lm_p(sample_missing$missing_fraction, sample_missing[[term]]),
             eta2 = qc_eta2(sample_missing$missing_fraction, sample_missing[[term]]))
}))
assoc$q_value <- p.adjust(assoc$p_value, method = "BH")
assoc <- assoc[order(-assoc$eta2), ]

group_terms <- intersect(c("Group", "group", "ExpGroup", "group2"), assoc$term)
tech_terms <- intersect(c("plate", "batch", "run", "order", "sample_prep"), assoc$term)
finite_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  max(x)
}
group_eta <- finite_max(assoc$eta2[assoc$term %in% group_terms])
tech_eta <- finite_max(assoc$eta2[assoc$term %in% tech_terms])
flag <- if (is.finite(group_eta) && group_eta >= 0.15) "WARN" else "PASS"
confounded <- if (is.finite(group_eta) && is.finite(tech_eta) && tech_eta >= 0.15) "possible" else "not_detected"

qc_write_csv(sample_missing, file.path(PATHS$tables, "missing_fraction_by_sample.csv"))
qc_write_csv(protein_missing, file.path(PATHS$tables, "missing_fraction_by_protein.csv"))
qc_write_csv(by_term, file.path(PATHS$tables, "missing_fraction_by_metadata.csv"))
qc_write_csv(assoc, file.path(PATHS$tables, "missingness_metadata_associations.csv"))
qc_write_xlsx(list(sample = sample_missing, protein = protein_missing, by_metadata = by_term, associations = assoc),
              file.path(PATHS$tables, "missingness_diagnostics.xlsx"))

p_sample <- ggplot(sample_missing, aes(x = reorder(Sample, missing_fraction), y = missing_fraction)) +
  geom_col(fill = "#4C78A8") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = NULL, y = "Missing fraction") +
  theme_classic(base_size = 8)
ggsave(file.path(PATHS$figures, "missing_fraction_by_sample.svg"), p_sample,
       width = 140, height = 120, units = "mm", device = svglite::svglite)

if (length(terms)) {
  color_term <- terms[[1]]
  p_box <- ggplot(sample_missing, aes(x = .data[[color_term]], y = missing_fraction)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", color = "grey30") +
    geom_point(position = position_jitter(width = 0.12), alpha = 0.75, size = 1) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = color_term, y = "Missing fraction") +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(PATHS$figures, "missingness_by_primary_metadata.svg"), p_box,
         width = 120, height = 90, units = "mm", device = svglite::svglite)
}

summary_lines <- c(
  "# Missingness Diagnostics",
  "",
  paste0("Dataset: ", DATASET),
  paste0("Input: ", relative_to(matrix_file)),
  paste0("Input appears imputed: ", input_is_imputed),
  if (input_is_imputed) "WARN: Missingness was assessed on an imputed-looking file; pre-imputation missingness may be underestimated." else "PASS: Input path does not look imputed.",
  paste0(flag, ": Group-associated missingness eta2=", ifelse(is.finite(group_eta), signif(group_eta, 3), "NA")),
  paste0("Group/technical missingness confounding: ", confounded)
)
writeLines(summary_lines, file.path(PATHS$reports, "missingness_summary.md"))

write_run_manifest(file.path(PATHS$logs, "run_manifest.yml"),
                   inputs = list(matrix = matrix_file, metadata = metadata_file),
                   outputs = list(tables = PATHS$tables, figures = PATHS$figures, reports = PATHS$reports),
                   parameters = list(dataset = DATASET, input_is_imputed = input_is_imputed),
                   notes = "Missingness diagnostics and Group confounding screen.")

message("Missingness diagnostics complete for dataset: ", DATASET)
