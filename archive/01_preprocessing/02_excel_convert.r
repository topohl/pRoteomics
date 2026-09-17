# Consumes:
#   - imputed Excel workbook or folder of imputed Excel workbooks
#   - sample metadata workbook
# Produces:
#   - metadata-augmented Excel workbooks and strict GCT v1.3 files
# File contract:
#   - ProTigy-targeted GCT v1.3 writer emits reserved id plus one explicit Description row descriptor, column metadata rows, and numeric sample columns

library(readxl)
library(dplyr)
library(writexl)
library(tools)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "protigy_input_utils.R"))
MODULE_ID <- "01_preprocessing"
SUBSTEP_ID <- "excel_convert"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

# ---- CONFIGURATION ----
# Set mode: "excel" for sheets in one Excel file, "folder" for all Excel files in a folder
mode <- Sys.getenv("PROTEOMICS_EXCEL_CONVERT_MODE", unset = "folder")

# Define file paths. Environment variables allow local overrides without
# committing machine-specific paths.
file_path <- Sys.getenv(
    "PROTEOMICS_EXCEL_CONVERT_FILE",
    unset = path_processed("01_preprocessing", "impute", "20260526_pgmatrix_imputed_neuron_neuropil_180samples_missing70pct.xlsx")
)
folder_path <- Sys.getenv(
    "PROTEOMICS_EXCEL_CONVERT_FOLDER",
    unset = path_processed("01_preprocessing", "impute")
)
metadata_path <- Sys.getenv(
    "PROTEOMICS_EXCEL_CONVERT_METADATA",
    unset = path_metadata("TPE9_sample_metadata_males.xlsx")
)
output_dir <- path_or_env("PROTEOMICS_EXCEL_CONVERT_OUTPUT_DIR", CANONICAL_PATHS$processed, kind = "dir")

if (!file.exists(metadata_path)) stop("Metadata file not found: ", metadata_path, call. = FALSE)
if (mode == "excel" && !file.exists(file_path)) stop("Excel input file not found: ", file_path, call. = FALSE)
if (mode == "folder" && !dir.exists(folder_path)) stop("Excel input folder not found: ", folder_path, call. = FALSE)
ensure_dir(output_dir)

# Read metadata
metadata <- read_excel(metadata_path)

# ---- READ DATA FILES ----
sheet_dfs <- list()
if (mode == "excel") {
    # Read all sheets from one Excel file
    sheet_names <- excel_sheets(file_path)
    sheet_dfs <- setNames(
        lapply(sheet_names, function(sheet) read_excel(file_path, sheet = sheet)),
        sheet_names
    )
} else if (mode == "folder") {
    # Read all Excel files in folder (ignore sheets, use file name as key)
    excel_files <- list.files(folder_path, pattern = "\\.xlsx$", full.names = TRUE)
    if (length(excel_files) == 0) stop("No .xlsx files found in folder: ", folder_path, call. = FALSE)
    sheet_dfs <- setNames(
        lapply(excel_files, function(f) read_excel(f)),
        basename(file_path_sans_ext(excel_files))
    )
}

# ---- PROCESS EACH DATA FRAME ----
new_dfs <- lapply(sheet_dfs, function(df) {
    # Preserve the historical ProTigy id/truncation semantics in the shared helper.
    df <- protigy_prepare_legacy_expression(df)
    sample_cols <- intersect(colnames(df), metadata$sample_id)
    if (length(sample_cols) == 0) return(df)
    # Build metadata rows
    meta_rows <- lapply(names(metadata), function(meta_col) {
        row <- rep(NA, ncol(df))
        names(row) <- colnames(df)
        row["id"] <- meta_col
        if (meta_col %in% colnames(metadata) && length(sample_cols) > 0) {
            idxs <- match(sample_cols, metadata$sample_id)
            valid <- !is.na(idxs) & idxs >= 1 & idxs <= nrow(metadata)
            if (any(valid)) {
                row[sample_cols[valid]] <- metadata[[meta_col]][idxs[valid]]
            }
        }
        row
    })
    # ---- Add combined row: region_layer_ExpGroup ----
    combined_row <- rep(NA, ncol(df))
    names(combined_row) <- colnames(df)
    combined_row["id"] <- "region_layer_ExpGroup"
    if (all(c("region", "layer", "ExpGroup") %in% colnames(metadata))) {
        idxs <- match(sample_cols, metadata$sample_id)
        valid <- !is.na(idxs) & idxs >= 1 & idxs <= nrow(metadata)
        if (any(valid)) {
            region <- as.character(metadata$region[idxs[valid]])
            layer <- as.character(metadata$layer[idxs[valid]])
            expgroup <- as.character(metadata$ExpGroup[idxs[valid]])
            combined <- paste(region, layer, expgroup, sep = "_")
            combined_row[sample_cols[valid]] <- combined
        }
    }
    meta_df <- as.data.frame(do.call(rbind, meta_rows), stringsAsFactors = FALSE)
    meta_df <- rbind(meta_df, as.data.frame(t(combined_row), stringsAsFactors = FALSE))
    # Combine metadata rows and original data
    final_df <- rbind(meta_df, df)
    # ---- Replace NA in layer and celltype_layer rows with "microglia" ----
    layer_idx <- which(final_df$id == "layer")
    celltype_layer_idx <- which(final_df$id == "celltype_layer")
    if (length(layer_idx) == 1) {
        final_df[layer_idx, sample_cols] <- as.character(final_df[layer_idx, sample_cols])
        na_idx <- which(is.na(final_df[layer_idx, sample_cols]) | final_df[layer_idx, sample_cols] == "NA")
        if (length(na_idx) > 0) {
            final_df[layer_idx, sample_cols[na_idx]] <- "microglia"
        }
    }
    if (length(celltype_layer_idx) == 1) {
        final_df[celltype_layer_idx, sample_cols] <- as.character(final_df[celltype_layer_idx, sample_cols])
        na_idx <- which(is.na(final_df[celltype_layer_idx, sample_cols]) | final_df[celltype_layer_idx, sample_cols] == "NA")
        if (length(na_idx) > 0) {
            final_df[celltype_layer_idx, sample_cols[na_idx]] <- "microglia"
        }
    }
    # ---- Add phenotypeWithinUnit row after NA replacement ----
    phenotype_row <- rep(NA, ncol(final_df))
    names(phenotype_row) <- colnames(final_df)
    phenotype_row["id"] <- "phenotypeWithinUnit"
    # Use combined_row values for phenotypeWithinUnit
    if (all(c("region", "layer", "ExpGroup") %in% colnames(metadata))) {
        idxs <- match(sample_cols, metadata$sample_id)
        valid <- !is.na(idxs) & idxs >= 1 & idxs <= nrow(metadata)
        if (any(valid)) {
            region <- as.character(metadata$region[idxs[valid]])
            layer <- as.character(final_df[layer_idx, sample_cols[valid]])
            expgroup <- as.character(metadata$ExpGroup[idxs[valid]])
            combined <- paste(region, layer, expgroup, sep = "_")
            phenotype_row[sample_cols[valid]] <- combined
        }
    }
    # Insert phenotypeWithinUnit row after layer and celltype_layer rows
    insert_idx <- max(layer_idx, celltype_layer_idx, na.rm = TRUE)
    final_df <- rbind(
        final_df[seq_len(insert_idx), , drop = FALSE],
        as.data.frame(t(phenotype_row), stringsAsFactors = FALSE),
        final_df[(insert_idx + 1):nrow(final_df), , drop = FALSE]
    )
    # ---- Remove columns where "exclude" metadata row is TRUE ----
    exclude_row_idx <- which(final_df$id == "exclude")
    if (length(exclude_row_idx) == 1) {
        exclude_row <- final_df[exclude_row_idx, , drop = FALSE]
        exclude_cols <- which(tolower(as.character(exclude_row)) == "true")
        exclude_cols <- setdiff(exclude_cols, which(colnames(final_df) == "id"))
        if (length(exclude_cols) > 0) {
            final_df <- final_df[, -exclude_cols, drop = FALSE]
        }
    }
    final_df
})

# ---- SAVE EACH DATA FRAME AS EXCEL ----
for (sheet in names(new_dfs)) {
    out_path <- file.path(output_dir, paste0(sheet, "_with_metadata.xlsx"))
    write_xlsx(new_dfs[[sheet]], out_path)
}

# ---- SAVE EACH DATA FRAME AS GCT v1.3 ----
self_test_write_gct_v1.3()

for (sheet in names(new_dfs)) {
    out_path <- file.path(output_dir, paste0(sheet, "_with_metadata.gct"))
    write_gct_v1.3(new_dfs[[sheet]], out_path, metadata)
}

write_run_manifest(
    file.path(CANONICAL_PATHS$logs, "run_manifest.yml"),
    inputs = list(file_path = file_path, folder_path = folder_path, metadata = metadata_path),
    outputs = list(output_dir = output_dir, sheets = names(new_dfs)),
    parameters = list(mode = mode, gct_format = "strict GCT v1.3"),
    notes = paste(
      "Excel export behavior is preserved; GCT writer self-test runs before writing outputs.",
      "ProTigy handoffs use one explicit Description row descriptor; source descriptions are retained",
      "when present and id is the compatibility fallback when absent."
    )
)
