library(readxl)
library(readr)
library(dplyr)
library(writexl)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

# Paths/config
IMPUTATION_SEED <- 42L
metadata_path <- Sys.getenv("PROTEOMICS_IMPUTE_METADATA", unset = path_metadata("TPE9_sample_metadata_males.xlsx"))
input_path <- Sys.getenv("PROTEOMICS_IMPUTE_INPUT", unset = path_raw("pg_matrix", "quicksearch.pg_matrix.tsv"))
output_dir <- Sys.getenv("PROTEOMICS_IMPUTE_OUTPUT_DIR", unset = path_processed("pg_matrix", "imputed"))

if (!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path, call. = FALSE)
if (!file.exists(input_path)) stop("Input file not found: ", input_path, call. = FALSE)
ensure_dir(output_dir)

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
celltype_layers <- unique(metadata$celltype_layer)

# Split by celltype_layer and process each subset
for (idx in seq_along(celltype_layers)) {
    celltype_layer <- celltype_layers[[idx]]
    # Get sample_ids for this celltype_layer
    subset_sample_ids <- metadata %>%
        filter(celltype_layer == !!celltype_layer) %>%
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

    # Median center each sample column
    for (col in numeric_cols) {
        df_subset[[col]] <- df_subset[[col]] - median(df_subset[[col]], na.rm = TRUE)
    }

    # Filter out proteins (rows) with >70% missingness in numeric columns
    missing_prop <- apply(df_subset[, numeric_cols], 1, function(x) mean(is.na(x)))
    df_filtered <- df_subset[missing_prop <= 0.7, ]

    # Impute missing values (on log2 scale)
    subset_seed <- IMPUTATION_SEED + idx - 1L
    imputed_df <- impute_normal(df_filtered, numeric_cols, seed = subset_seed)

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
        celltype_layer = celltype_layer,
        seed = subset_seed,
        n_samples = n_samples,
        n_proteins = n_proteins,
        output_file = basename(output_path),
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
}
