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

# ------------------------------------------------
# 1) SETTINGS
# ------------------------------------------------

main_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/compareGO/BP/phenotype_within_unit"

target_subdir <- "01_Tables_and_Data"
target_file   <- "07_Top_Genes_Driving_TopTerms.xlsx"

min_datasets <- 3

output_file <- file.path(
  main_dir,
  paste0("overlapping_proteins_min", min_datasets, "_datasets.xlsx")
)

# Optional: restrict to expected datasets
expected_datasets <- c(
  "CA1_slm", "CA1_so", "CA1_sr",
  "CA2_slm", "CA2_so", "CA2_sr",
  "CA3_sr",
  "DG_mo", "DG_po"
)

# ------------------------------------------------
# 2) FIND INPUT FILES
# ------------------------------------------------

dataset_dirs <- list.dirs(main_dir, recursive = FALSE, full.names = TRUE)

input_files <- tibble(
  Dataset = basename(dataset_dirs),
  File = file.path(dataset_dirs, target_subdir, target_file)
) %>%
  filter(Dataset %in% expected_datasets) %>%
  mutate(File_exists = file.exists(File))

missing_files <- input_files %>%
  filter(!File_exists)

input_files_found <- input_files %>%
  filter(File_exists)

if (nrow(input_files_found) == 0) {
  stop("No input files found. Check main_dir, folder names, and target_file.")
}

message("Found ", nrow(input_files_found), " input files.")
if (nrow(missing_files) > 0) {
  message("Missing files in: ", paste(missing_files$Dataset, collapse = ", "))
}

# ------------------------------------------------
# 3) READ AND STANDARDIZE DATA
# ------------------------------------------------

read_top_genes <- function(file, dataset) {
  df <- readxl::read_excel(file)

  required_cols <- c("Description", "Gene", "Freq", "Mean_NES")
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
      Dataset = dataset,
      Description = as.character(Description),
      Gene = as.character(Gene),
      Gene = str_trim(Gene),
      Freq = suppressWarnings(as.numeric(Freq)),
      Mean_NES = suppressWarnings(as.numeric(Mean_NES))
    ) %>%
    filter(!is.na(Gene), Gene != "") %>%
    distinct(Dataset, Description, Gene, .keep_all = TRUE)
}

all_data <- purrr::map2_dfr(
  input_files_found$File,
  input_files_found$Dataset,
  read_top_genes
)

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