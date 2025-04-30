#' UniProt ID Mapping for ClusterProfiler Analysis
#'
#' This script processes gene data with UniProtKB IDs (e.g., "GENE_MOUSE") and maps them
#' to UniProt Accessions using an official UniProt ID mapping file.
#'
#' @details
#' The workflow includes:
#' \enumerate{
#'   \item Setting up the environment and loading necessary libraries.
#'   \item Defining project directories.
#'   \item Checking for the existence of input data and the mapping file.
#'   \item Preprocessing gene identifiers by extracting the primary ID when multiple are provided.
#'   \item Reading and filtering the official UniProt mapping file (MOUSE_10090_idmapping.dat) to obtain an Entry Name to Accession map.
#'   \item Merging the mapping data with the gene information.
#'   \item Saving the mapped results for downstream analysis.
#' }
#'
#' @note Ensure that the mapping file is downloaded, unzipped, and placed in the specified directory before running this script.
#'
#' @author Tobias Pohl

# -----------------------------------------------------
# Load Libraries
# -----------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE))
    install.packages("pacman")

# Load required libraries with pacman
pacman::p_load(dplyr, stringr, tidyr, purrr, readr)

# -----------------------------------------------------
# Define Paths and Project Directories
# -----------------------------------------------------

working_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler"
mapped_dir <- file.path(working_dir, "Datasets", "mapped")
unmapped_dir <- file.path(working_dir, "Datasets", "unmapped")

# Create directories if they do not exist
if (!dir.exists(mapped_dir)) {
    dir.create(mapped_dir, recursive = TRUE)
}
if (!dir.exists(unmapped_dir)) {
    dir.create(unmapped_dir, recursive = TRUE)
}

setwd(working_dir)

# Specify cell types which determine dataset file name
cell_types <- c("mcherry1", "mcherry3")
file_name <- paste0(paste(cell_types, collapse = "_"), ".csv")
data_path <- file.path(working_dir, "Datasets", "raw", file_name)

# IMPORTANT: Set path for the downloaded UniProt mapping file
uniprot_mapping_file_path <- file.path(working_dir, "Datasets", "MOUSE_10090_idmapping.dat")

# -----------------------------------------------------
# Verify Input Files
# -----------------------------------------------------

if (!file.exists(data_path)) {
    stop("Input data file does not exist: ", data_path)
}
if (!file.exists(uniprot_mapping_file_path)) {
    cat("UniProt mapping file not found at:", uniprot_mapping_file_path, "\n")
    cat("Attempting to download the file...\n")
    
    # Increase timeout to 3600 seconds (1 hour) to prevent download timeout
    options(timeout = 3600)
    
    gz_url <- "https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz"
    gz_file <- paste0(uniprot_mapping_file_path, ".gz")
    tryCatch({
        download.file(gz_url, gz_file, mode = "wb")
        if (!requireNamespace("R.utils", quietly = TRUE)) {
            install.packages("R.utils")
        }
        R.utils::gunzip(gz_file, destname = uniprot_mapping_file_path, remove = TRUE)
        cat("Downloaded and unzipped the UniProt mapping file successfully.\n")
    }, error = function(e) {
        stop("Failed to download or unzip the UniProt mapping file: ", e$message)
    })
}

# -----------------------------------------------------
# Preprocess Gene Data
# -----------------------------------------------------

cat("Reading gene data from:", data_path, "\n")
df <- readr::read_csv2(data_path, col_names = TRUE, show_col_types = FALSE) %>%
    mutate(gene_symbol = as.character(gene_symbol)) %>%
    # In the case of semicolon-separated IDs, retain the primary (first) identifier
    mutate(gene_symbol = map_chr(str_split(gene_symbol, ";"), ~str_trim(.x[1]))) %>%
    filter(gene_symbol != "" & !is.na(gene_symbol)) %>%
    distinct(gene_symbol, .keep_all = TRUE)

cat("Total genes loaded:", nrow(df), "\n")
if (nrow(df) == 0) {
    stop("No valid gene symbols were identified after preprocessing.")
}
cat("Preview of gene data BEFORE mapping:\n")
print(head(df))

# -----------------------------------------------------
# Parse UniProt Mapping File
# -----------------------------------------------------

cat("Loading UniProt mapping file from:", uniprot_mapping_file_path, "\n")
uniprot_mapping <- readr::read_tsv(
    uniprot_mapping_file_path,
    col_names = c("UniProt_Accession", "Type", "Value"),
    col_types = "ccc",  # reading all columns as characters
    quote = ""         # avoid handling quotes in the identifiers
)

cat("Total entries in UniProt mapping file:", nrow(uniprot_mapping), "\n")

# Extract mappings for UniProtKB-ID to Accession
entry_name_to_accession <- uniprot_mapping %>%
    dplyr::filter(Type == "UniProtKB-ID") %>%
    dplyr::select(UniProt_Accession, UniProtKB_ID = Value) %>% 
    dplyr::distinct(UniProtKB_ID, .keep_all = TRUE)

cat("Total unique UniProtKB-ID mappings identified:", nrow(entry_name_to_accession), "\n")
if (nrow(entry_name_to_accession) == 0) {
    stop("Failed to obtain UniProtKB-ID mappings from the file.")
}

# -----------------------------------------------------
# Map Gene IDs to UniProt Accessions
# -----------------------------------------------------

cat("Mapping gene identifiers to UniProt Accessions...\n")
df_mapped <- df %>%
    left_join(entry_name_to_accession, by = c("gene_symbol" = "UniProtKB_ID")) %>%
    # Replace the original gene_symbol with the Accession if mapping exists
    mutate(
        gene_symbol = if_else(
            !is.na(UniProt_Accession) & UniProt_Accession != "",
            UniProt_Accession,
            gene_symbol
        )
    ) %>%
    # Retain and order the columns for downstream analysis (adjust if necessary)
    dplyr::select(gene_symbol, pval, padj, log2fc)

cat("Preview of gene data AFTER mapping:\n")
print(head(df_mapped))

cat("Final gene count:", nrow(df_mapped), "\n")

# report proteins that werent mapped
unmapped_proteins <- df %>%
  left_join(entry_name_to_accession, by = c("gene_symbol" = "UniProtKB_ID")) %>%
  filter(is.na(UniProt_Accession) | UniProt_Accession == "") %>%
  select(gene_symbol)
cat("Total unmapped proteins:", nrow(unmapped_proteins), "\n")
cat("Unmapped proteins:\n")
print(unmapped_proteins)

# -----------------------------------------------------
# Save Mapped Results
# -----------------------------------------------------

mapped_file <- file.path(mapped_dir, paste0(tools::file_path_sans_ext(file_name), ".csv"))
readr::write_csv(df_mapped, mapped_file)
cat("Mapped results saved to:", mapped_file, "\n")

unmapped_file <- file.path(unmapped_dir, paste0(tools::file_path_sans_ext(file_name), ".csv"))
readr::write_csv(unmapped_proteins, unmapped_file)
cat("Unmapped proteins saved to:", unmapped_file, "\n")

cat("Pipeline completed successfully.\n")
