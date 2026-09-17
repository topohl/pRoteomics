# Consumes:
#   - sample/protein metadata workbooks from data/metadata/
# Produces:
#   - processed protein metadata TSV and long metadata workbook under data/processed/01_preprocessing/metadata_create/
# File contract:
#   - docs/active_script_io_audit.tsv object 01_preprocessing/05_metadata_create.r

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
MODULE_ID <- "01_preprocessing"
SUBSTEP_ID <- "metadata_create"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

work_direction <- path_metadata()
metadata_file <- file.path(work_direction, "TPE9_sample_metadata_males.xlsx")
files <- list.files(path = work_direction, pattern = "TPE9_samples_males\\.xlsx$", full.names = TRUE)
if (is_dry_run()) {
    dry_run_line("Script", "01_preprocessing/05_metadata_create.r")
    dry_run_line("Metadata directory", work_direction, if (dir.exists(work_direction)) "PASS" else "FAIL")
    dry_run_line("Sample metadata workbook", metadata_file, if (file.exists(metadata_file)) "PASS" else "FAIL")
    dry_run_line("Protein workbook count", length(files), if (length(files) > 0) "PASS" else "FAIL")
    dry_run_line("Output folders", paste(unlist(CANONICAL_PATHS), collapse = "; "))
    quit(status = if (file.exists(metadata_file) && length(files) > 0) 0 else 1, save = "no")
}
if (!file.exists(metadata_file)) stop("Sample metadata workbook not found: ", metadata_file, call. = FALSE)
if (length(files) == 0) {
    stop("No matching .xlsx files found in work_direction: ", work_direction)
}

packages <- c("readxl", "dplyr", "tidyr", "stringr", "purrr", "tibble", "writexl")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
    stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "),
         ". Install them explicitly before running this script.", call. = FALSE)
}
invisible(lapply(packages, library, character.only = TRUE))
metadata <- readxl::read_xlsx(metadata_file)
safe_read <- purrr::possibly(readxl::read_xlsx, otherwise = tibble::tibble())
data <- purrr::map_dfr(files, safe_read) %>%
    dplyr::filter(!is.na(`Protein.Group`))

# save data as filename but with _processed.tsv
output_file <- file.path(CANONICAL_PATHS$processed, "TPE9_samples_males_processed.tsv")
# write as tab-separated values, no row names, no quotes, empty strings for NA
write.table(data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# reshape protein table to long and join with metadata, then save
keep_cols <- c("Protein.Group", "Protein.Names", "Genes", "First.Protein.Description")

# identify sample columns (those that are not protein annotation columns)
sample_cols <- setdiff(names(data), keep_cols)

# Prepare a conservative "cleaning" function so column names and metadata IDs can be matched
clean_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    # drop common prefixes/suffixes you may have in column names, adjust if needed
    str_remove_all("^Intensity\\.|^X|^Sample[_-]") %>%
    str_remove_all("\\.raw$|\\.mzML$|_1$|_2$") %>%
    str_replace_all("[^A-Za-z0-9_-]+", "_") %>%
    str_to_lower()
}

# add cleaned IDs to metadata so we can match robustly
if (!"sample_id" %in% names(metadata)) {
  stop("metadata must contain a column named 'sample_id' to match against sample columns")
}
metadata <- metadata %>% mutate(sample_id_clean = clean_id(sample_id))

# create mapping from original column names to cleaned sample_id
sample_map <- tibble(original = sample_cols) %>%
  mutate(sample_id_clean = clean_id(original))

# pivot the sample columns to long form using the original column names, then attach cleaned id and metadata
long_proteins <- data %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "original",
    values_to = "intensity"
  ) %>%
  mutate(intensity = as.numeric(intensity)) %>%
  left_join(sample_map, by = "original") %>%             # add cleaned sample id
  left_join(metadata, by = "sample_id_clean") %>%        # add metadata columns
  filter(!is.na(intensity))

writexl::write_xlsx(long_proteins, file.path(CANONICAL_PATHS$processed, "TPE9_samples_males_long_with_metadata.xlsx"))
write_session_info(file.path(CANONICAL_PATHS$logs, "sessionInfo.txt"))
