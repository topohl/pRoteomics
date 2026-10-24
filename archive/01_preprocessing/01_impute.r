# Consumes:
#   - raw pg_matrix TSV and sample metadata workbook
# Produces:
#   - imputed per-celltype_layer workbooks plus imputation QC and run manifest
# File contract:
#   - sample IDs in metadata$sample_id must match sample columns after the first four annotation columns

library(readxl)
library(readr)
library(dplyr)
library(writexl)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "normalisation_contract.R"))

# Paths/config
MODULE_ID <- "01_preprocessing"
SUBSTEP_ID <- "impute"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)
IMPUTATION_SEED <- 42L
metadata_path <- path_or_env("PROTEOMICS_IMPUTE_METADATA", path_metadata("TPE9_sample_metadata_males.xlsx"), must_exist = TRUE)
input_path <- path_or_env("PROTEOMICS_IMPUTE_INPUT", path_raw("pg_matrix", "quicksearch.pg_matrix.tsv"), must_exist = TRUE)
output_dir <- path_or_env("PROTEOMICS_IMPUTE_OUTPUT_DIR", CANONICAL_PATHS$processed, kind = "dir")
ensure_dir(output_dir)
# config/normalisation.yml selects the centring. Writing into the working tree under a setting other than the one it
# holds is refused unless that version is preserved; an override directory keeps its own manifest and stamp.
NORMALISATION <- normalisation_setting()
same_dir <- function(a, b) identical(tolower(sub("/+$", "", normalizePath(a, winslash = "/", mustWork = FALSE))),
                                     tolower(sub("/+$", "", normalizePath(b, winslash = "/", mustWork = FALSE))))
writes_working_tree <- same_dir(output_dir, CANONICAL_PATHS$processed)
# Any directory that already holds a stamped version (however its path is spelled) is protected the same way.
if (writes_working_tree && !identical(NORMALISATION$source, "config/normalisation.yml")) {
    stop("PROTEOMICS_SAMPLE_CENTRING may only redirect a run into a separate output directory (PROTEOMICS_IMPUTE_OUTPUT_DIR); ",
         "the working tree follows config/normalisation.yml.", call. = FALSE)
}
if (writes_working_tree || dir_holds_version(output_dir)) assert_working_version_preserved(NORMALISATION$sample_centring, output_dir)
manifest_dir <- if (writes_working_tree) CANONICAL_PATHS$logs else file.path(output_dir, "run_manifest")
if (!writes_working_tree) Sys.setenv(PROTEOMICS_NORMALISATION_STAMP_DIR = output_dir)
message("Sample centring: ", NORMALISATION$sample_centring, " (", NORMALISATION$source, ")")

# Read metadata
metadata <- read_excel(metadata_path)
if (!"sample_id" %in% names(metadata)) stop("Metadata must contain a sample_id column.", call. = FALSE)

# Exclude rows where 'exclude' is TRUE
metadata <- metadata %>% filter(is.na(exclude) | exclude != TRUE)

# Read data to impute, treating "." as decimal separator
df <- read_tsv(input_path, locale = locale(decimal_mark = "."))

# === Rename Protein.Names to T: Protein.Names ===
if ("Protein.Names" %in% names(df)) {
    names(df)[names(df) == "Protein.Names"] <- "T: Protein.Names"
}
# ===============================================

# Identify annotation columns (first 4 columns, adjust if needed)
annotation_cols <- 1:4
if (ncol(df) <= length(annotation_cols)) {
    stop("Input matrix has no sample columns after annotation columns.", call. = FALSE)
}
annotation_names <- names(df)[annotation_cols]

# Check sample_id matching
sample_ids_metadata <- metadata$sample_id
sample_ids_data <- names(df)[-annotation_cols]
missing_in_data <- setdiff(sample_ids_metadata, sample_ids_data)
missing_in_metadata <- setdiff(sample_ids_data, sample_ids_metadata)
if (length(missing_in_data) > 0) {
    warning("These sample_ids from metadata are missing in data: ", paste(missing_in_data, collapse = ", "))
}
if (length(missing_in_metadata) > 0) {
    warning("These sample_ids from data are missing in metadata: ", paste(missing_in_metadata, collapse = ", "))
}

# Only keep sample_ids present in both
common_sample_ids <- intersect(sample_ids_metadata, sample_ids_data)
if (length(common_sample_ids) == 0) {
    stop("No common sample IDs between metadata and input matrix.", call. = FALSE)
}

# Function to impute missing values (expects log2 scale)
impute_normal <- function(df, numeric_cols, width = 0.3, downshift = 1.8, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    imputed_df <- df
    for (col in numeric_cols) {
        col_data <- df[[col]]
        missing_idx <- which(is.na(col_data))
        if (length(missing_idx) > 0) {
            observed <- col_data[!is.na(col_data)]
            if (length(observed) == 0) {
                stop("Cannot impute column with no observed values: ", names(df)[col], call. = FALSE)
            }
            mean_obs <- mean(observed)
            sd_obs <- sd(observed)
            if (!is.finite(sd_obs) || sd_obs == 0) sd_obs <- 0
            impute_mean <- mean_obs - downshift * sd_obs
            impute_sd <- sd_obs * width
            imputed_values <- rnorm(length(missing_idx), mean = impute_mean, sd = impute_sd)
            imputed_df[missing_idx, col] <- imputed_values
        }
    }
    return(imputed_df)
}

# Helper function to get current date string
get_date_str <- function() {
    format(Sys.Date(), "%Y%m%d")
}

# Helper function to create scientific filenames
make_filename <- function(celltype_layer, n_samples, n_proteins, method = "normal", missing_thresh = 0.7) {
    # Clean celltype_layer for filename
    celltype_layer_clean <- gsub("[^A-Za-z0-9]+", "_", celltype_layer)
    date_str <- get_date_str()
    paste0(
        date_str, "_",
        "pgmatrix_imputed_",
        celltype_layer_clean,
        "_", n_samples, "samples",
        "_missing", missing_thresh*100, "pct.xlsx"
    )
}

qc_rows <- list()
# Consumers read the newest matching workbook; never leave a second version beside a new one. Checked for every
# celltype_layer before anything is written, so a refusal cannot leave a half-regenerated directory.
today_names <- character()
metadata <- metadata %>%
    mutate(celltype_layer_for_impute = ifelse(is.na(celltype_layer) | !nzchar(as.character(celltype_layer)), "missing_celltype_layer", as.character(celltype_layer)))
celltype_layers <- sort(unique(as.character(metadata$celltype_layer_for_impute)))
existing <- list.files(output_dir, pattern = "^[0-9]{8}_pgmatrix_imputed_.+_[0-9]+samples_missing70pct[.]xlsx$")
stale <- existing[!startsWith(existing, paste0(get_date_str(), "_"))]
if (length(stale)) {
    stop("Imputed workbook(s) from another run in ", output_dir, ": ", paste(stale, collapse = ", "),
         ". Move them to a superseded folder (keeping a hash record) before regenerating.", call. = FALSE)
}

# Split by celltype_layer and process each subset
for (idx in seq_along(celltype_layers)) {
    celltype_layer <- celltype_layers[[idx]]
    # Get sample_ids for this celltype_layer
    subset_sample_ids <- metadata %>%
        filter(celltype_layer_for_impute == !!celltype_layer) %>%
        pull(sample_id)
    subset_sample_ids <- intersect(subset_sample_ids, common_sample_ids)

    # If no samples, skip
    if (length(subset_sample_ids) == 0) next

    # Subset data: annotation columns + sample columns
    subset_cols <- c(annotation_names, subset_sample_ids)
    df_subset <- df[, subset_cols]

    # Identify numeric columns for this subset (sample columns only)
    numeric_cols <- seq.int(length(annotation_names) + 1L, ncol(df_subset))
    if (length(numeric_cols) == 0) {
        stop("No numeric sample columns available for celltype_layer: ", celltype_layer, call. = FALSE)
    }

    # Log2 transform numeric columns
    df_subset[, numeric_cols] <- log2(df_subset[, numeric_cols])

    # Median centre each sample column (config/normalisation.yml): on the proteins quantified in every
    # sample of this celltype_layer (complete_case_median, 2026-10-19) or on all observed values of the
    # sample (all_observed_median, the original analysis). The observed-value median moves with
    # detection depth: a sample that missed more low-abundance proteins has a higher observed median,
    # so every protein in it is shifted down.
    offsets <- sample_centring_offsets(as.data.frame(df_subset[, numeric_cols]), NORMALISATION$sample_centring)
    for (j in seq_along(numeric_cols)) {
        df_subset[[numeric_cols[j]]] <- df_subset[[numeric_cols[j]]] - offsets[[j]]
    }

    # Filter out proteins (rows) with >70% missingness in numeric columns
    n_proteins_before_filter <- nrow(df_subset)
    missing_prop <- apply(df_subset[, numeric_cols], 1, function(x) mean(is.na(x)))
    df_filtered <- df_subset[missing_prop <= 0.7, ]
    n_missing_before_imputation <- sum(is.na(df_filtered[, numeric_cols]))

    # Impute missing values (on log2 scale)
    subset_seed <- IMPUTATION_SEED + idx - 1L
    imputed_df <- impute_normal(df_filtered, numeric_cols, seed = subset_seed)
    n_missing_after_imputation <- sum(is.na(imputed_df[, numeric_cols]))

    # Move annotation columns to the end
    anno_idx <- match(annotation_names, names(imputed_df))
    sample_idx <- setdiff(seq_along(imputed_df), anno_idx)
    imputed_df <- imputed_df[, c(sample_idx, anno_idx)]

    # Create scientific filename
    n_samples <- length(subset_sample_ids)
    n_proteins <- nrow(imputed_df)
    filename <- make_filename(celltype_layer, n_samples, n_proteins)
    output_path <- file.path(output_dir, filename)

    # Write output
    write_xlsx(imputed_df, output_path)

    qc_rows[[length(qc_rows) + 1L]] <- data.frame(
        base_seed = IMPUTATION_SEED,
        subset_seed = subset_seed,
        celltype_layer = celltype_layer,
        n_samples = n_samples,
        sample_centring = NORMALISATION$sample_centring,
        centring = if (identical(NORMALISATION$sample_centring, "complete_case_median")) "median of proteins quantified in all samples of the celltype_layer" else "median of all observed values of the sample",
        n_centring_reference_proteins = attr(offsets, "n_reference_proteins"),
        n_proteins_before_filter = n_proteins_before_filter,
        n_proteins_after_filter = n_proteins,
        n_proteins_filtered_out = n_proteins_before_filter - n_proteins,
        n_missing_before_imputation = n_missing_before_imputation,
        n_missing_after_imputation = n_missing_after_imputation,
        n_values_imputed = n_missing_before_imputation - n_missing_after_imputation,
        output_path = output_path,
        stringsAsFactors = FALSE
    )
}

if (length(qc_rows) > 0) {
    imputation_qc <- do.call(rbind, qc_rows)
    write_csv(imputation_qc, file.path(output_dir, "imputation_qc.csv"))
    writeLines(c(
        paste("IMPUTATION_SEED:", IMPUTATION_SEED),
        paste("metadata_path:", metadata_path),
        paste("input_path:", input_path),
        paste("output_dir:", output_dir),
        "",
        capture.output(sessionInfo())
    ), file.path(output_dir, "sessionInfo.txt"))
    write_normalisation_stamp(output_dir, NORMALISATION, imputation_qc$output_path)
    write_run_manifest(
        file.path(manifest_dir, "run_manifest.yml"),
        inputs = list(metadata = metadata_path, matrix = input_path),
        outputs = list(imputed_outputs = imputation_qc$output_path, qc = file.path(output_dir, "imputation_qc.csv"),
                       normalisation_stamp = normalisation_stamp_path(output_dir)),
        parameters = list(imputation_seed = IMPUTATION_SEED, missing_threshold = 0.7, width = 0.3, downshift = 1.8,
                          sample_centring = NORMALISATION$sample_centring, sample_centring_source = NORMALISATION$source),
        notes = "Deterministic per-celltype_layer seeds are assigned after sorting celltype_layer labels."
    )
}
