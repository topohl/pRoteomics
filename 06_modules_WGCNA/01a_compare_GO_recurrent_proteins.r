# ================================================================
# Script: 06_modules_WGCNA/01a_compare_GO_recurrent_proteins.r
# Stage: modules_downstream
# Scope: global
# Consumes: required data/processed/04_differential_expression_enrichment/compareGO/neuron_neuropil/compareGO_input_manifest.csv; data/processed/04_differential_expression_enrichment/compareGO/neuron_soma/compareGO_input_manifest.csv; +1 more; optional none.
# Produces: results/tables/06_modules_WGCNA/01a_compare_GO_recurrent_proteins/.
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Global recurrent GO/protein comparison used by curated overlap programs.
# ================================================================

# ================================================================
# Compare overlapping proteins across GSEA top-gene Excel files
# Keeps proteins present in at least N datasets
# ================================================================

library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(openxlsx)
library(tibble)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))

# ------------------------------------------------
# 1) SETTINGS
# ------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[[1]] == length(args)) return(default)
  args[[hit[[1]] + 1L]]
}

DATASET <- validate_dataset(value_after("--dataset", Sys.getenv("PROTEOMICS_DATASET", unset = NA_character_)), source = "--dataset")
Sys.setenv(PROTEOMICS_DATASET = DATASET)
dry_run <- is_dry_run()
target_file   <- "07_Top_Genes_Driving_TopTerms.xlsx"

min_datasets <- 3

comparego_dir <- path_results(
  "tables",
  "04_differential_expression_enrichment",
  "compareGO",
  DATASET,
  "BP",
  "phenotype_within_unit"
)

aggregate_input_file <- file.path(comparego_dir, target_file)

output_dir <- path_results(
  "tables",
  "04_differential_expression_enrichment",
  "compareGO",
  DATASET,
  "recurrent_top_gene_drivers"
)
dir_create(output_dir)

output_file <- file.path(output_dir, paste0("overlapping_proteins_min", min_datasets, "_spatial_units.xlsx"))

expected_spatial_units <- switch(
  DATASET,
  neuron_neuropil = c(
    "CA1_slm", "CA1_so", "CA1_sr",
    "CA2_slm", "CA2_so", "CA2_sr",
    "CA3_so", "CA3_sr",
    "DG_mo", "DG_po"
  ),
  neuron_soma = c("CA1", "CA2", "CA3", "DG"),
  microglia = c("CA1", "CA2", "CA3", "DG")
)

if (dry_run) {
  dry_run_line("Script", "06_modules_WGCNA/01a_compare_GO_recurrent_proteins.r")
  dry_run_line("Dataset", DATASET)
  dry_run_line("compareGO aggregate top-gene file", aggregate_input_file, if (file.exists(aggregate_input_file)) "PASS" else "FAIL")
  dry_run_line("Expected spatial units", paste(expected_spatial_units, collapse = ", "), "INFO")
  dry_run_line("Output workbook", output_file, "INFO")
  quit(status = if (file.exists(aggregate_input_file)) 0L else 1L, save = "no")
}

# ------------------------------------------------
# 2) FIND INPUT FILES
# ------------------------------------------------

input_files <- tibble(
  Input_type = "compareGO_aggregate_top_gene_drivers",
  File = aggregate_input_file,
  File_exists = file.exists(File)
)

if (!file.exists(aggregate_input_file)) {
  stop("No input file found. Expected compareGO output: ", aggregate_input_file, call. = FALSE)
}

message("Found compareGO input file: ", aggregate_input_file)

# ------------------------------------------------
# 3) READ AND STANDARDIZE DATA
# ------------------------------------------------

compact_dataset_key <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("_", "")
}

dataset_lookup <- tibble(
  Dataset = expected_spatial_units,
  Dataset_key = compact_dataset_key(expected_spatial_units)
)

parse_dataset_from_comparison <- function(x) {
  x_key <- compact_dataset_key(x)
  detected <- dataset_lookup$Dataset[stringr::str_detect(x_key, dataset_lookup$Dataset_key)]
  if (length(detected) == 0) NA_character_ else detected[[1]]
}

read_top_genes <- function(file) {
  df <- readxl::read_excel(file)

  required_cols <- c("Description", "Gene", "Freq", "Mean_NES", "Comparisons")
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      "File missing required columns: ",
      basename(file), " | Missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df %>%
    mutate(
      Description = as.character(Description),
      Gene = as.character(Gene),
      Gene = str_trim(Gene),
      Freq = suppressWarnings(as.numeric(Freq)),
      Mean_NES = suppressWarnings(as.numeric(Mean_NES)),
      Comparisons = as.character(Comparisons)
    ) %>%
    filter(!is.na(Gene), Gene != "") %>%
    tidyr::separate_rows(Comparisons, sep = ";\\s*") %>%
    mutate(
      Comparison = stringr::str_trim(Comparisons),
      Dataset = purrr::map_chr(Comparison, parse_dataset_from_comparison)
    ) %>%
    filter(!is.na(Dataset), Dataset %in% expected_spatial_units) %>%
    group_by(Dataset, Description, Gene) %>%
    summarise(
      Freq = n_distinct(Comparison),
      Mean_NES = mean(Mean_NES, na.rm = TRUE),
      Comparisons = paste(sort(unique(Comparison)), collapse = "; "),
      .groups = "drop"
    )
}

all_data <- read_top_genes(aggregate_input_file)

found_datasets <- sort(unique(all_data$Dataset))
missing_datasets <- setdiff(expected_spatial_units, found_datasets)
message("Detected datasets: ", paste(found_datasets, collapse = ", "))
if (length(missing_datasets) > 0) {
  message("Missing expected spatial units in compareGO comparisons: ", paste(missing_datasets, collapse = ", "))
}

if (nrow(all_data) == 0) {
  stop("No dataset-level rows could be parsed from the compareGO Comparisons column.", call. = FALSE)
}

# ------------------------------------------------
# 4) DATASET-LEVEL OVERLAP
# ------------------------------------------------

protein_dataset_summary <- all_data %>%
  distinct(Gene, Dataset) %>%
  group_by(Gene) %>%
  summarise(
    N_datasets = n_distinct(Dataset),
    Datasets = paste(sort(unique(Dataset)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(N_datasets >= min_datasets) %>%
  arrange(desc(N_datasets), Gene)

# ------------------------------------------------
# 5) TERM-LEVEL DETAILS FOR OVERLAPPING PROTEINS
# ------------------------------------------------

overlap_details <- all_data %>%
  semi_join(protein_dataset_summary, by = "Gene") %>%
  group_by(Gene, Dataset) %>%
  summarise(
    N_terms = n_distinct(Description),
    Terms = paste(sort(unique(Description)), collapse = "; "),
    Mean_Freq = mean(Freq, na.rm = TRUE),
    Max_Freq = max(Freq, na.rm = TRUE),
    Mean_NES_mean = mean(Mean_NES, na.rm = TRUE),
    Mean_NES_max_abs = Mean_NES[which.max(abs(Mean_NES))],
    .groups = "drop"
  ) %>%
  left_join(protein_dataset_summary, by = "Gene") %>%
  arrange(desc(N_datasets), Gene, Dataset)

# ------------------------------------------------
# 6) WIDE PRESENCE / ABSENCE MATRIX
# ------------------------------------------------

presence_matrix <- all_data %>%
  distinct(Gene, Dataset) %>%
  mutate(Present = 1L) %>%
  semi_join(protein_dataset_summary, by = "Gene") %>%
  pivot_wider(
    names_from = Dataset,
    values_from = Present,
    values_fill = 0
  ) %>%
  left_join(protein_dataset_summary, by = "Gene") %>%
  relocate(Gene, N_datasets, Datasets) %>%
  arrange(desc(N_datasets), Gene)

# ------------------------------------------------
# 7) DESCRIPTION-LEVEL OVERLAPS
# Optional: proteins shared within similar GO descriptions
# ------------------------------------------------

description_overlap <- all_data %>%
  distinct(Description, Gene, Dataset) %>%
  group_by(Description, Gene) %>%
  summarise(
    N_datasets = n_distinct(Dataset),
    Datasets = paste(sort(unique(Dataset)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(N_datasets >= min_datasets) %>%
  arrange(desc(N_datasets), Description, Gene)

# ------------------------------------------------
# 8) WRITE EXCEL OUTPUT
# ------------------------------------------------

wb <- createWorkbook()

addWorksheet(wb, "Overlap_summary")
writeData(wb, "Overlap_summary", protein_dataset_summary)

addWorksheet(wb, "Overlap_details")
writeData(wb, "Overlap_details", overlap_details)

addWorksheet(wb, "Presence_matrix")
writeData(wb, "Presence_matrix", presence_matrix)

addWorksheet(wb, "Description_overlap")
writeData(wb, "Description_overlap", description_overlap)

addWorksheet(wb, "Input_files")
writeData(wb, "Input_files", input_files)

header_style <- createStyle(
  textDecoration = "bold",
  border = "Bottom",
  fgFill = "#D9EAF7"
)

for (sheet in names(wb)) {
  addStyle(
    wb,
    sheet = sheet,
    style = header_style,
    rows = 1,
    cols = 1:200,
    gridExpand = TRUE
  )
  freezePane(wb, sheet = sheet, firstRow = TRUE)
  setColWidths(wb, sheet = sheet, cols = 1:200, widths = "auto")
}

saveWorkbook(wb, output_file, overwrite = TRUE)

message("Done.")
message("Output written to: ", output_file)
