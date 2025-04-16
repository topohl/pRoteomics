#' @title Gene Symbol Mapping for ClusterProfiler Analysis
#'
#' @description
#' This script processes a gene dataset by cleaning gene symbols and mapping them to corresponding UniProt 
#' SwissProt identifiers using the biomaRt package. It handles entries with multiple gene symbols by 
#' retaining only the primary symbol and then integrates UniProt database data to update these symbols when available.
#'
#' @details
#' The script performs the following steps:
#' \itemize{
#'   \item Loads necessary libraries including dplyr, stringr, tidyr, purrr, biomaRt, and readr.
#'   \item Sets the working directory and defines directories for results and datasets.
#'   \item Constructs the file name dynamically based on predefined cell types and verifies that the
#'         dataset exists.
#'   \item Reads the gene data from a CSV file, transforms gene symbol entries with multiple entries by 
#'         splitting on ";" and retaining only the first entry.
#'   \item Filters out empty or duplicate gene symbols to ensure a clean dataset.
#'   \item Queries the UniProt database via biomaRt to get corresponding SwissProt identifiers for gene names.
#'   \item Joins the gene dataset with the UniProt data, replacing gene symbols with their SwissProt
#'         identifiers when available.
#' }
#'
#' @note 
#' Ensure that the CSV data file is present at the designated path. The script will terminate with an error
#' if the file does not exist.
#'
#' @examples
#' \dontrun{
#'   # Set up working directory and run the script to map gene symbols to SwissProt identifiers.
#' }
#'
#' @author
#' Tobias Pohl

library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(biomaRt)
library(readr)

# Set working directory and results directory paths
working_dir <- "/Users/tobiaspohl/Documents/clusterProfiler"
results_dir <- file.path(working_dir, "Results")
setwd(working_dir)

# Define cell types and construct the CSV file name accordingly
cell_types <- c("mcherryG1", "mcherryG2")
file_name <- paste0(paste(cell_types, collapse = "_"), ".csv")
data_path <- file.path(working_dir, "Datasets", file_name)

# Verify that the data file exists before proceeding
if (!file.exists(data_path)) {
    stop("Data file does not exist: ", data_path)
}

# Load the gene dataset, convert gene symbols to character vectors,
# split entries with multiple genes at the semicolon, and use only the first gene.
# Also, remove empty entries and duplicate gene symbols.
df <- read_csv(data_path) %>%
    mutate(gene_symbol = as.character(gene_symbol)) %>%
    mutate(gene_symbol = str_split(gene_symbol, ";")) %>%
    mutate(gene_symbol = map_chr(gene_symbol, ~str_trim(.x[1]))) %>%
    filter(gene_symbol != "") %>%
    distinct(gene_symbol, .keep_all = TRUE)

# Connect to the UniProt database using biomaRt
uniprot <- useMart("uniprot", dataset = "uniprotkb")

# Retrieve UniProt SwissProt identifiers and associated gene names,
# group by gene name to consolidate multiple identifiers,
# and filter out any invalid or missing gene names.
gene_symbols_db <- getBM(attributes = c("uniprot_swissprot", "gene_name"), mart = uniprot) %>%
    group_by(gene_name) %>%
    summarise(uniprot_swissprot = paste(unique(uniprot_swissprot), collapse = ";"),
              .groups = "drop") %>%
    filter(!is.na(gene_name) & gene_name != "")

# Integrate the UniProt data with the gene dataset by performing a left join.
# Replace the gene symbol with the corresponding SwissProt identifier when available.
df <- df %>%
    left_join(gene_symbols_db, by = c("gene_symbol" = "gene_name")) %>%
    mutate(gene_symbol = if_else(uniprot_swissprot != "" & !is.na(uniprot_swissprot),
                                 uniprot_swissprot,
                                 gene_symbol)) %>%
    select(gene_symbol, everything(), -uniprot_swissprot)
