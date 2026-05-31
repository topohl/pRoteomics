# Dataset-aware variance partitioning for spatial proteomics.

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))
source(repo_path("R", "qc_exploration_utils.R"))

run <- qc_args()
DATASET <- run$dataset
PATHS <- qc_paths("06_variance_partitioning", DATASET)
matrix_file <- path_or_env("PROTEOMICS_VARPART_MATRIX_FILE", qc_resolve_matrix(DATASET), must_exist = FALSE)
metadata_file <- qc_resolve_metadata(DATASET, env = "PROTEOMICS_VARPART_METADATA_FILE")

if (run$dry_run) {
  status <- qc_dry_run_contract("03_qc_exploration/06_variance_partitioning.r", DATASET,
                                matrix_file = matrix_file, metadata_file = metadata_file, paths = PATHS,
                                extra = c("Formula is adapted to available non-collinear metadata terms.",
                                          "Uses (1|AnimalID) only when repeated samples per animal exist."))
  quit(status = status, save = "no")
}

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "svglite", "variancePartition")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(required_pkgs, library, character.only = TRUE)))
if (!file.exists(matrix_file)) stop("VariancePartition matrix not found: ", matrix_file, call. = FALSE)

expr <- qc_read_expression(matrix_file, metadata_file, DATASET)
mat <- qc_impute_for_pca(expr$mat)
meta <- expr$meta[colnames(mat), , drop = FALSE]
meta <- as.data.frame(meta, stringsAsFactors = FALSE)

candidate_terms <- c("Group", "group", "ExpGroup", "Sex", "sex", "Region", "region", "Layer", "layer",
                     "ReplicateGroup", "plate", "batch", "sample_prep", "run", "order")
candidate_terms <- intersect(candidate_terms, names(meta))
usable <- candidate_terms[vapply(candidate_terms, function(term) {
  x <- meta[[term]]
  sum(!is.na(x) & nzchar(as.character(x))) >= 3L && length(unique(x[!is.na(x) & nzchar(as.character(x))])) >= 2L
}, logical(1))]
usable <- usable[!duplicated(tolower(usable))]

random_terms <- character()
if ("AnimalID" %in% names(meta) && any(table(meta$AnimalID) > 1L)) random_terms <- c(random_terms, "(1|AnimalID)")
for (term in intersect(c("Region", "region", "Layer", "layer", "ReplicateGroup", "plate", "batch"), usable)) {
  if (length(unique(meta[[term]])) >= 3L) random_terms <- c(random_terms, paste0("(1|", term, ")"))
}
fixed_terms <- setdiff(usable, gsub("^\\(1\\|(.+)\\)$", "\\1", random_terms))
formula_terms <- c(fixed_terms, random_terms)
if (!length(formula_terms)) stop("No usable metadata terms for variance partitioning.", call. = FALSE)
form <- as.formula(paste("~", paste(formula_terms, collapse = " + ")))

cca <- try(variancePartition::canCorPairs(form, meta), silent = TRUE)
if (!inherits(cca, "try-error")) {
  cca_df <- as.data.frame(as.table(as.matrix(cca)))
  names(cca_df) <- c("term_1", "term_2", "canonical_correlation")
  qc_write_csv(cca_df, file.path(PATHS$tables, "metadata_term_canonical_correlations.csv"))
}

vp <- variancePartition::fitExtractVarPartModel(mat, form, meta)
vp <- variancePartition::sortCols(vp)
vp_df <- as.data.frame(vp)
vp_df$Protein <- rownames(vp_df)
vp_long <- vp_df |>
  tidyr::pivot_longer(-Protein, names_to = "term", values_to = "variance_fraction")
median_vp <- vp_long |>
  dplyr::group_by(term) |>
  dplyr::summarise(median_variance_fraction = median(variance_fraction, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(median_variance_fraction))
top_driven <- vp_long |>
  dplyr::group_by(term) |>
  dplyr::slice_max(variance_fraction, n = 50, with_ties = FALSE) |>
  dplyr::ungroup()

qc_write_csv(vp_long, file.path(PATHS$tables, "variance_fraction_by_protein.csv"))
qc_write_csv(median_vp, file.path(PATHS$tables, "median_variance_by_term.csv"))
qc_write_csv(top_driven, file.path(PATHS$tables, "top_term_driven_proteins.csv"))
qc_write_xlsx(list(variance_fraction = vp_long, median_by_term = median_vp, top_term_driven = top_driven),
              file.path(PATHS$tables, "variance_partitioning.xlsx"))

p <- ggplot(vp_long, aes(x = reorder(term, variance_fraction, median, na.rm = TRUE), y = variance_fraction, fill = term)) +
  geom_violin(scale = "width", linewidth = 0.2, alpha = 0.8) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2) +
  coord_flip() +
  labs(x = NULL, y = "Fraction of variance explained") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none")
ggsave(file.path(PATHS$figures, "variance_partition_summary.svg"), p,
       width = 130, height = 95, units = "mm", device = svglite::svglite)

group_terms <- intersect(c("Group", "group", "ExpGroup", "group2"), names(meta))
technical_terms <- intersect(c("plate", "batch", "sample_prep", "run", "order"), names(meta))
confound_rows <- dplyr::bind_rows(lapply(group_terms, function(g) {
  dplyr::bind_rows(lapply(technical_terms, function(t) {
    data.frame(group_term = g, technical_term = t,
               eta2_group_by_technical = qc_eta2(seq_len(nrow(meta)), interaction(meta[[g]], meta[[t]], drop = TRUE)),
               group_technical_eta2 = qc_eta2(as.numeric(factor(meta[[g]])), meta[[t]]))
  }))
}))
if (nrow(confound_rows)) qc_write_csv(confound_rows, file.path(PATHS$tables, "group_technical_confounding_screen.csv"))
flag <- if (nrow(confound_rows) && any(confound_rows$group_technical_eta2 >= 0.8, na.rm = TRUE)) "WARN" else "PASS"

writeLines(c("# Variance Partitioning", "", paste0("Dataset: ", DATASET),
             paste0("Formula: ", deparse(form)),
             paste0(flag, ": Group/batch/region/layer confounding screen written where metadata allowed.")),
           file.path(PATHS$reports, "variance_partitioning_summary.md"))

write_run_manifest(file.path(PATHS$logs, "run_manifest.yml"),
                   inputs = list(matrix = matrix_file, metadata = metadata_file),
                   outputs = list(tables = PATHS$tables, figures = PATHS$figures, reports = PATHS$reports),
                   parameters = list(dataset = DATASET, formula = deparse(form)),
                   notes = "Canonical variance partitioning with adaptive formula and confounding checks.")

message("Variance partitioning complete for dataset: ", DATASET)

