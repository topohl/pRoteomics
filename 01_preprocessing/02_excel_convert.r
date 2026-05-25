# Consumes:
#   - imputed Excel workbook or folder of imputed Excel workbooks
#   - sample metadata workbook
# Produces:
#   - metadata-augmented Excel workbooks and strict GCT v1.3 files
# File contract:
#   - GCT v1.3 writer emits #1.3, four-field dimensions, one row metadata column, column metadata rows, and numeric sample columns

library(readxl)
library(dplyr)
library(writexl)
library(tools)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
MODULE_ID <- "01_preprocessing"
SUBSTEP_ID <- "excel_convert"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

# ---- CONFIGURATION ----
# Set mode: "excel" for sheets in one Excel file, "folder" for all Excel files in a folder
mode <- "folder" # or "folder"

# Define file paths. Environment variables allow local overrides without
# committing machine-specific paths.
file_path <- Sys.getenv(
    "PROTEOMICS_EXCEL_CONVERT_FILE",
    unset = path_processed("pg_matrix", "imputed", "pg_matrix_groupfiltered_70percentvalid_imputed_updated_Feb2026.xlsx")
)
folder_path <- Sys.getenv(
    "PROTEOMICS_EXCEL_CONVERT_FOLDER",
    unset = path_processed("pg_matrix", "imputed", "grouped")
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
    # Rename "T: Protein.Names" to "id"
    # Rename "T: Protein.Names" or "Protein.Names" to "id"
    name_col <- intersect(c("T: Protein.Names", "Protein.Names"), colnames(df))
    if (length(name_col) == 1) {
        colnames(df)[colnames(df) == name_col] <- "id"
        df <- df[, c("id", setdiff(colnames(df), "id"))]
    }
    # Remove columns starting from "T: Protein.Group" or "Protein.Group" and after
    group_col <- intersect(c("T: Protein.Group", "Protein.Group"), colnames(df))
    if (length(group_col) == 1) {
        idx <- which(colnames(df) == group_col)
        if (idx > 1) {
            df <- df[, seq_len(idx - 1)]
        } else {
            # If group_col is the first column, remove all columns (or handle as needed)
            df <- df[, FALSE]
        }
    }
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
validate_gct_v1.3 <- function(file, expected_rows, expected_sample_cols, expected_row_metadata_cols, expected_col_metadata_rows) {
    first_lines <- readLines(file, n = 2, warn = FALSE)
    if (length(first_lines) < 2 || !identical(first_lines[[1]], "#1.3")) {
        stop("Invalid GCT marker in written file: ", file, call. = FALSE)
    }
    dims <- strsplit(first_lines[[2]], "\t", fixed = TRUE)[[1]]
    if (length(dims) != 4 || any(is.na(suppressWarnings(as.integer(dims))))) {
        stop("Invalid GCT v1.3 dimension line in written file: ", file, call. = FALSE)
    }
    dims <- as.integer(dims)
    expected <- c(expected_rows, expected_sample_cols, expected_row_metadata_cols, expected_col_metadata_rows)
    if (!identical(dims, expected)) {
        stop(
            "GCT v1.3 dimension mismatch in written file: ", file,
            ". Expected ", paste(expected, collapse = "\t"),
            " but found ", paste(dims, collapse = "\t"),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

write_gct_v1.3 <- function(df, file, metadata) {
    meta_row_idx <- which(df$id %in% c(names(metadata), "region_layer_ExpGroup", "phenotypeWithinUnit"))
    data_rows <- df[-meta_row_idx, , drop = FALSE]
    data_rows <- data_rows[!is.na(data_rows$id) & data_rows$id != "", , drop = FALSE]
    sample_cols <- intersect(setdiff(colnames(df), "id"), metadata$sample_id)
    sample_cols <- intersect(sample_cols, colnames(data_rows))
    if (length(sample_cols) == 0 || nrow(data_rows) == 0) {
        warning(sprintf("No data to write for file: %s", file))
        return()
    }
    is_numeric_row <- function(row) {
        all(suppressWarnings(!is.na(as.numeric(row[sample_cols]))))
    }
    numeric_rows_idx <- apply(data_rows, 1, is_numeric_row)
    data_rows <- data_rows[numeric_rows_idx, , drop = FALSE]
    for (col in sample_cols) {
        data_rows[[col]] <- as.numeric(data_rows[[col]])
    }
    if (nrow(data_rows) == 0) {
        warning(sprintf("No numeric data to write for file: %s", file))
        return()
    }
    meta_rows <- df[meta_row_idx, c("id", sample_cols), drop = FALSE]
    meta_rows <- meta_rows[!is.na(meta_rows$id) & meta_rows$id != "", , drop = FALSE]
    row_metadata_cols <- "id"
    con <- file(file, "wt")
    on.exit(close(con))
    writeLines("#1.3", con)
    writeLines(sprintf("%d\t%d\t%d\t%d", nrow(data_rows), length(sample_cols), length(row_metadata_cols), nrow(meta_rows)), con)
    writeLines(paste(c(row_metadata_cols, sample_cols), collapse = "\t"), con)
    if (nrow(meta_rows) > 0) {
        for (i in seq_len(nrow(meta_rows))) {
            writeLines(paste(unlist(meta_rows[i, c(row_metadata_cols, sample_cols)]), collapse = "\t"), con)
        }
    }
    for (i in seq_len(nrow(data_rows))) {
        writeLines(paste(unlist(data_rows[i, c(row_metadata_cols, sample_cols)]), collapse = "\t"), con)
    }
    close(con)
    on.exit()
    validate_gct_v1.3(file, nrow(data_rows), length(sample_cols), length(row_metadata_cols), nrow(meta_rows))
}

self_test_write_gct_v1.3 <- function() {
    tmp <- tempfile(fileext = ".gct")
    test_metadata <- data.frame(
        sample_id = c("sample1", "sample2", "sample3"),
        group = c("A", "B", "A"),
        stringsAsFactors = FALSE
    )
    test_df <- data.frame(
        id = c("group", "protein_a", "protein_b"),
        sample1 = c("A", 1, 4),
        sample2 = c("B", 2, 5),
        sample3 = c("A", 3, 6),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    write_gct_v1.3(test_df, tmp, test_metadata)
    lines <- readLines(tmp, warn = FALSE)
    stopifnot(identical(lines[[1]], "#1.3"))
    stopifnot(identical(lines[[2]], "2\t3\t1\t1"))
    invisible(TRUE)
}

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
    notes = "Excel export behavior is preserved; GCT writer self-test runs before writing outputs."
)
