# read in xlsx files from directory

if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
}
pacman::p_load(readxl, dplyr, tidyr, stringr, purrr, tibble, writexl)

work_direction <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics"
metadata <- readxl::read_xlsx(file.path(work_direction, "TPE9_sample_metadata_males.xlsx"))

files <- list.files(path = work_direction, pattern = "TPE9_samples_males\\.xlsx$", full.names = TRUE)
if (length(files) == 0) {
    stop("No matching .xlsx files found in work_direction: ", work_direction)
}
safe_read <- purrr::possibly(readxl::read_xlsx, otherwise = tibble::tibble())
data <- purrr::map_dfr(files, safe_read) %>%
    dplyr::filter(!is.na(`Protein.Group`))

# save data as filename but with _processed.tsv
output_file <- file.path(work_direction, "TPE9_samples_males_processed.tsv")
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

writexl::write_xlsx(long_proteins, file.path(work_direction, "TPE9_samples_males_long_with_metadata.xlsx"))