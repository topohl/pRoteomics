# Consumes:
#   - raw sample metadata workbook from data/metadata/sample_metadata_R.xlsx
# Produces:
#   - processed sample metadata workbook under data/metadata/
# File contract:
#   - docs/active_script_io_audit.tsv object 01_preprocessing/04_format_metadata.r

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

file_path <- path_metadata("sample_metadata_R.xlsx")
output_file <- path_metadata("sample_metadata_R_processed.xlsx")
qc_csv_file <- path_metadata("sample_metadata_R_processed_QC.csv")
qc_txt_file <- path_metadata("sample_metadata_R_processed_QC.txt")
MODULE_ID <- "01_preprocessing"
SUBSTEP_ID <- "format_metadata"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

if (is_dry_run()) {
  dry_run_line("Script", "01_preprocessing/04_format_metadata.r")
  dry_run_line("Metadata workbook", file_path, if (file.exists(file_path)) "PASS" else "FAIL")
  dry_run_line("Output file", output_file)
  quit(status = if (file.exists(file_path)) 0 else 1, save = "no")
}
if (!file.exists(file_path)) stop("Metadata workbook not found: ", file_path, call. = FALSE)

# Install and load required packages only when explicitly requested.
packages <- c("readxl", "writexl", "dplyr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "),
       ". Install them explicitly before running this script.", call. = FALSE)
}
lapply(packages, library, character.only = TRUE)

# Get available sheet names
sheet_names <- excel_sheets(file_path)
print(sheet_names)

# Load the sheet named "samples" into a data frame
df <- read_excel(file_path, sheet = "samples")
if (!"sample_id" %in% names(df)) stop("Metadata workbook sheet 'samples' must contain sample_id.", call. = FALSE)
head(df)

# extract info from column "sample_id"; example basename:
# Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A0003_L_CA1_Microglia_S113_S2-A3_1_13026.d
# in Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A0003_L_CA1_so_Neuron_S113_S2-A3_1_13026.d 
# A0003 is AnimalID, _L_ or _R_ indicate left or right replicates and go into "ReplicateGroup".
# _13026 is the run number and goes into "run_order"
df <- df %>%
  mutate(
    AnimalID = sub(".*_(A[0-9]+)_.*", "\\1", sample_id),
    ReplicateGroup = ifelse(grepl("_L_", sample_id), "Left", "Right"),
    run_order = sub(".*_(\\d+)\\.d$", "\\1", sample_id),
    group = sub(".*_(?:A[0-9]+)_(?:L|R)_([^_]+(?:_[^_]+)*)_S\\d+_.*", "\\1", sample_id),
    exclude = ifelse(grepl("CA3_slm", sample_id), TRUE, FALSE),
    plate = ifelse(grepl("20250703", sample_id), "B",
             ifelse(grepl("20250707", sample_id), "C", NA)),
    sample_number = sub(".*_(S\\d+)_.*", "\\1", sample_id),
    sample_location = ifelse(grepl("_[S]\\d+-", sample_id),
                             sub(".*_(S\\d+)-([A-Z]\\d+)_.*", "\\2", sample_id),
                             sub(".*_(S\\d+)_.*", "\\1", sample_id)),
    shortname = paste(plate, sample_number, sample_location, run_order, sep="_"),
    group2 = paste(group, ReplicateGroup, sep = "_"),
    region = sub(".*_(L|R)_((?:CA1|CA2|CA3|DG))_.*", "\\2", sample_id),
    # Keep string "NA" for layer/celltype_layer because 02_excel_convert.r
    # historically maps these metadata rows to microglia for Morpheus exports.
    layer = ifelse(grepl(".*_(?:L|R)_[^_]+_((?:so|sp|sr|slm|mo|po|sg))_.*", sample_id),
             sub(".*_(?:L|R)_[^_]+_((?:so|sp|sr|slm|mo|po|sg))_.*", "\\1", sample_id),
             "NA"),
    celltype = ifelse(
      grepl("_(?:L|R)_[^_]+_Microglia", sample_id, perl = TRUE),
      "microglia",
      ifelse(
      grepl("_(?:L|R)_[^_]+_[^_]+_Neuron", sample_id, perl = TRUE),
      "neuron",
      NA
      )
    ),
    celltype_layer = case_when(
      grepl("_(sp|sg)_Neuron", sample_id) ~ "neuron_soma",
      grepl("_Neuron", sample_id) ~ "neuron_neuropil",
      grepl("Microglia", sample_id) ~ "NA",
      TRUE ~ "NA"
    )
  )

if (any(df$AnimalID == df$sample_id, na.rm = TRUE)) {
  failed <- df$sample_id[df$AnimalID == df$sample_id]
  stop("AnimalID regex failed for sample_id(s): ", paste(failed, collapse = ", "), call. = FALSE)
}
if (any(df$run_order == df$sample_id | is.na(df$run_order) | !nzchar(df$run_order), na.rm = TRUE)) {
  failed <- df$sample_id[df$run_order == df$sample_id | is.na(df$run_order) | !nzchar(df$run_order)]
  stop("run_order regex failed for sample_id(s): ", paste(failed, collapse = ", "), call. = FALSE)
}
if (any(df$sample_number == df$sample_id | is.na(df$sample_number) | !nzchar(df$sample_number), na.rm = TRUE)) {
  failed <- df$sample_id[df$sample_number == df$sample_id | is.na(df$sample_number) | !nzchar(df$sample_number)]
  stop("sample_number regex failed for sample_id(s): ", paste(failed, collapse = ", "), call. = FALSE)
}
if (any(df$sample_location == df$sample_id | is.na(df$sample_location) | !nzchar(df$sample_location), na.rm = TRUE)) {
  failed <- df$sample_id[df$sample_location == df$sample_id | is.na(df$sample_location) | !nzchar(df$sample_location)]
  stop("sample_location regex failed for sample_id(s): ", paste(failed, collapse = ", "), call. = FALSE)
}

allowed_regions <- c("CA1", "CA2", "CA3", "DG")
bad_regions <- unique(df$region[!is.na(df$region) & !df$region %in% allowed_regions])
if (length(bad_regions) > 0) {
  stop("Unexpected or failed region value(s): ", paste(bad_regions, collapse = ", "), call. = FALSE)
}

allowed_layers <- c("so", "sp", "sr", "slm", "mo", "po", "sg", "NA")
bad_layers <- unique(df$layer[!is.na(df$layer) & !df$layer %in% allowed_layers])
if (length(bad_layers) > 0) {
  warning("Unexpected layer value(s): ", paste(bad_layers, collapse = ", "), call. = FALSE)
}

ca3_slm_rows <- df %>% filter(region == "CA3", layer == "slm")
ca3_slm_excluded <- nrow(ca3_slm_rows) > 0 && all(ca3_slm_rows$exclude == TRUE, na.rm = FALSE)
if (nrow(ca3_slm_rows) > 0 && !ca3_slm_excluded) {
  stop("CA3_slm rows must be exclude == TRUE for all matching samples.", call. = FALSE)
}

qc_counts <- bind_rows(
  data.frame(qc_section = "n_samples", value = "total", n = nrow(df), stringsAsFactors = FALSE),
  df %>% count(region, name = "n") %>% transmute(qc_section = "region", value = as.character(region), n),
  df %>% count(layer, name = "n") %>% transmute(qc_section = "layer", value = as.character(layer), n),
  df %>% count(celltype, name = "n") %>% transmute(qc_section = "celltype", value = as.character(celltype), n),
  df %>% count(celltype_layer, name = "n") %>% transmute(qc_section = "celltype_layer", value = as.character(celltype_layer), n),
  df %>% count(exclude, name = "n") %>% transmute(qc_section = "exclude", value = as.character(exclude), n)
)

interaction_qc <- df %>%
  count(region, layer, celltype, exclude, name = "n") %>%
  mutate(
    qc_section = "region_x_layer_x_celltype_x_exclude",
    value = paste(region, layer, celltype, exclude, sep = " | ")
  ) %>%
  select(qc_section, value, n)

ca3_slm_qc <- data.frame(
  qc_section = "CA3_slm",
  value = if (nrow(ca3_slm_rows) > 0) "all exclude == TRUE" else "no CA3_slm rows detected",
  n = nrow(ca3_slm_rows),
  stringsAsFactors = FALSE
)

qc_table <- bind_rows(qc_counts, interaction_qc, ca3_slm_qc)
write.csv(qc_table, qc_csv_file, row.names = FALSE)
writeLines(c(
  "Metadata parsing QC",
  paste("Input:", file_path),
  paste("Output:", output_file),
  paste("Samples:", nrow(df)),
  paste("Unexpected regions:", if (length(bad_regions)) paste(bad_regions, collapse = ", ") else "none"),
  paste("Unexpected layers:", if (length(bad_layers)) paste(bad_layers, collapse = ", ") else "none"),
  paste("CA3_slm rows:", nrow(ca3_slm_rows)),
  paste("CA3_slm exclude all TRUE:", if (nrow(ca3_slm_rows) > 0) ca3_slm_excluded else "not applicable"),
  "",
  capture.output(print(qc_table, row.names = FALSE))
), qc_txt_file)

# save as excel file using writexl package
write_xlsx(df, output_file)
write_run_manifest(
  file.path(CANONICAL_PATHS$logs, "run_manifest.yml"),
  inputs = list(metadata_workbook = file_path),
  outputs = list(processed_metadata = output_file, qc_csv = qc_csv_file, qc_txt = qc_txt_file),
  parameters = list(
    allowed_regions = allowed_regions,
    allowed_layers = allowed_layers,
    layer_na_representation = "string NA retained for downstream Morpheus metadata-row compatibility"
  ),
  notes = "Parsing QC fails on AnimalID, region, run_order, sample_number, and sample_location extraction failures."
)
