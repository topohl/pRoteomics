library(readxl)
library(readr)
library(dplyr)
library(writexl)

# Paths
metadata_path <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/metadata/TPE9_sample_metadata_males.xlsx"
input_path <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/pg_matrix/raw/quicksearch.pg_matrix.tsv"
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/pg_matrix/imputed/"

# Read metadata
metadata <- read_excel(metadata_path)

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

# Function to impute missing values (expects log2 scale)
impute_normal <- function(df, numeric_cols, width = 0.3, downshift = 1.8) {
    imputed_df <- df
    for (col in numeric_cols) {
        col_data <- df[[col]]
        missing_idx <- which(is.na(col_data))
        if (length(missing_idx) > 0) {
            observed <- col_data[!is.na(col_data)]
            mean_obs <- mean(observed)
            sd_obs <- sd(observed)
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

# Split by celltype_layer and process each subset
for (celltype_layer in unique(metadata$celltype_layer)) {
    # Get sample_ids for this celltype_layer
    subset_sample_ids <- metadata %>%
        filter(celltype_layer == !!celltype_layer) %>%
        pull(sample_id)
    subset_sample_ids <- intersect(subset_sample_ids, common_sample_ids)
    
    # If no samples, skip
    if (length(subset_sample_ids) == 0) next
    
    # Subset data
    subset_cols <- c(annotation_names, subset_sample_ids)
    df_subset <- df[, subset_cols]
    
    # Identify numeric columns for this subset
    numeric_cols <- (length(annotation_names) + 1):ncol(df_subset)
    
    # Log2 transform numeric columns (like Perseus)
    df_subset[, numeric_cols] <- log2(df_subset[, numeric_cols])
    
    # Filter out proteins (rows) with >70% missingness in numeric columns
    missing_prop <- apply(df_subset[, numeric_cols], 1, function(x) mean(is.na(x)))
    df_filtered <- df_subset[missing_prop <= 0.7, ]
    
    # Impute missing values (on log2 scale)
    imputed_df <- impute_normal(df_filtered, numeric_cols)
    
    # === Move annotation columns to the end ===
    sample_idx <- setdiff(seq_along(imputed_df), match(annotation_names, names(imputed_df)))
    anno_idx <- match(annotation_names, names(imputed_df))
    imputed_df <- imputed_df[, c(sample_idx, anno_idx)]
    # =========================================
    
    # Create scientific filename
    n_samples <- length(subset_sample_ids)
    n_proteins <- nrow(imputed_df)
    filename <- make_filename(celltype_layer, n_samples, n_proteins)
    output_path <- file.path(output_dir, filename)
    
    # Write output
    write_xlsx(imputed_df, output_path)
}
