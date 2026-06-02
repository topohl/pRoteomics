# ================================================================
# Build clean merged metadata sheet for module-score analysis
# ================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(openxlsx)
library(janitor)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "dataset_inputs.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = "") {
  hit <- which(args == flag)
  if (!length(hit) || hit[1] == length(args)) return(default)
  args[[hit[1] + 1]]
}
dataset_cli <- arg_value("--dataset", default = "")
if (nzchar(dataset_cli)) Sys.setenv(PROTEOMICS_DATASET = validate_dataset(dataset_cli, source = "--dataset"))
dataset_profile <- current_dataset()
dataset_inputs <- resolve_dataset_inputs(dataset_profile, purpose = "module_score")

# ------------------------------------------------
# 1) PATHS
# ------------------------------------------------

first_existing_path <- function(paths) {
  paths <- unique(normalizePath(paths[nzchar(paths)], winslash = "/", mustWork = FALSE))
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

latest_matching_file <- function(root, pattern) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!dir.exists(root)) return(NA_character_)

  candidates <- list.files(root, pattern = pattern, full.names = TRUE, recursive = TRUE)
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) return(NA_character_)

  info <- file.info(candidates)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

latest_matching_file_anywhere <- function(roots, pattern) {
  candidates <- vapply(roots, latest_matching_file, character(1), pattern = pattern)
  candidates <- candidates[!is.na(candidates) & file.exists(candidates)]
  if (length(candidates) == 0) return(NA_character_)

  info <- file.info(candidates)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}

analysis_root <- function() {
  normalizePath(dirname(repo_root()), winslash = "/", mustWork = FALSE)
}

resolve_module_score_proteomics_file <- function() {
  override <- Sys.getenv("PROTEOMICS_MODULE_SCORE_PROTEOMICS_FILE", unset = "")
  if (nzchar(override)) return(normalizePath(override, winslash = "/", mustWork = FALSE))
  dataset_inputs$expression_file
}

resolve_module_score_auc_all_file <- function() {
  override <- Sys.getenv("PROTEOMICS_MODULE_SCORE_AUC_ALL", unset = "")
  if (nzchar(override)) return(normalizePath(override, winslash = "/", mustWork = FALSE))

  direct <- first_existing_path(c(
    path_external("behavior", "auc_individual_animals_all.csv"),
    file.path(analysis_root(), "Behavior", "RFID", "statistics", "gamm_new", "analyses", "gamm", "tables", "auc_individual_animals_all.csv"),
    file.path(analysis_root(), "Behavior", "RFID", "statistics", "analyses_", "gamm", "tables", "auc_individual_animals_all.csv")
  ))
  if (!is.na(direct)) return(direct)

  latest_matching_file(
    file.path(analysis_root(), "Behavior"),
    "^auc_individual_animals_all\\.csv$"
  )
}

resolve_module_score_auc_first_file <- function() {
  override <- Sys.getenv("PROTEOMICS_MODULE_SCORE_AUC_FIRST", unset = "")
  if (nzchar(override)) return(normalizePath(override, winslash = "/", mustWork = FALSE))

  direct <- first_existing_path(c(
    path_external("behavior", "auc_individual_animals_firstChangeActive.csv"),
    file.path(analysis_root(), "Behavior", "RFID", "statistics", "gamm_new", "analyses", "gamm", "tables", "auc_individual_animals_firstChangeActive.csv"),
    file.path(analysis_root(), "Behavior", "RFID", "statistics", "analyses_", "gamm", "tables", "auc_individual_animals_firstChangeActive.csv")
  ))
  if (!is.na(direct)) return(direct)

  latest_matching_file(
    file.path(analysis_root(), "Behavior"),
    "^auc_individual_animals_firstChangeActive\\.csv$"
  )
}

resolve_module_score_behavior_file <- function() {
  override <- Sys.getenv("PROTEOMICS_MODULE_SCORE_BEHAVIOR_FILE", unset = "")
  if (nzchar(override)) return(normalizePath(override, winslash = "/", mustWork = FALSE))

  direct <- first_existing_path(c(
    path_external("behavior", "E9_Behavior_Data.xlsx"),
    file.path(analysis_root(), "SIS_Analysis", "E9_Behavior_Data.xlsx")
  ))
  if (!is.na(direct)) return(direct)

  latest_matching_file_anywhere(
    c(file.path(analysis_root(), "SIS_Analysis"), file.path(analysis_root(), "Behavior")),
    "^E9_Behavior_Data\\.xlsx$"
  )
}

proteomics_file <- resolve_module_score_proteomics_file()
auc_all_file <- resolve_module_score_auc_all_file()
auc_first_file <- resolve_module_score_auc_first_file()
behavior_file <- resolve_module_score_behavior_file()

input_files <- c(
  "Proteomics workbook" = proteomics_file,
  "AUC all file" = auc_all_file,
  "AUC first active file" = auc_first_file,
  "Behavior z-score workbook" = behavior_file
)
missing_inputs <- names(input_files)[is.na(input_files) | !file.exists(input_files)]
if (length(missing_inputs) > 0) {
  missing_lines <- paste0(missing_inputs, ": ", input_files[missing_inputs])
  stop(
    "Required module-score merge input file(s) not found:\n",
    paste(missing_lines, collapse = "\n"),
    "\n\nIf the proteomics workbook is missing, run source('01_preprocessing/02_excel_convert.r') ",
    "or set PROTEOMICS_MODULE_SCORE_PROTEOMICS_FILE to the metadata-augmented workbook.",
    "\nIf behavior inputs are missing, set PROTEOMICS_MODULE_SCORE_AUC_ALL, ",
    "PROTEOMICS_MODULE_SCORE_AUC_FIRST, or PROTEOMICS_MODULE_SCORE_BEHAVIOR_FILE.",
    call. = FALSE
  )
}

message("Using proteomics workbook: ", proteomics_file)
message("Using AUC all file: ", auc_all_file)
message("Using AUC first active file: ", auc_first_file)
message("Using behavior z-score workbook: ", behavior_file)
message("Resolved dataset: ", dataset_profile)
message("Resolved dataset inputs diagnostics: ", paste(dataset_inputs$diagnostics, collapse = " | "))

out_dir <- Sys.getenv("PROTEOMICS_MODULE_SCORE_OUTPUT_DIR", unset = path_results("module_scores", dataset_profile))
ensure_dir(out_dir)

out_file <- file.path(out_dir, "sample_metadata_merged_clean_for_module_scores.xlsx")
global_compat_out_file <- path_results("module_scores", "sample_metadata_merged_clean_for_module_scores.xlsx")
write_global_compat <- identical(dataset_profile, "neuron_neuropil") || tolower(Sys.getenv("PROTEOMICS_WRITE_GLOBAL_MODULE_SCORE_METADATA_COMPAT", unset = "")) %in% c("1", "true", "yes")

if (is_dry_run()) {
  dry_run_line("Script", "01_preprocessing/06_merged_metadata_module_score.r")
  dry_run_line("Dataset", dataset_profile)
  dry_run_line("Proteomics workbook", proteomics_file, if (file.exists(proteomics_file)) "PASS" else "FAIL")
  dry_run_line("AUC all file", auc_all_file, if (file.exists(auc_all_file)) "PASS" else "FAIL")
  dry_run_line("AUC first active file", auc_first_file, if (file.exists(auc_first_file)) "PASS" else "FAIL")
  dry_run_line("Behavior z-score workbook", behavior_file, if (file.exists(behavior_file)) "PASS" else "FAIL")
  dry_run_line("Dataset-scoped output", out_file)
  dry_run_line("Legacy global compatibility copy", if (write_global_compat) global_compat_out_file else "disabled")
  quit(status = if (all(file.exists(c(proteomics_file, auc_all_file, auc_first_file, behavior_file)))) 0 else 1, save = "no")
}

# ------------------------------------------------
# 2) HELPERS
# ------------------------------------------------

standardize_animal_id <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_to_upper(x)

  numeric_core <- str_extract(x, "\\d+$")

  numeric_core <- ifelse(
    !is.na(numeric_core),
    as.character(as.integer(numeric_core)),
    NA_character_
  )

  ifelse(!is.na(numeric_core), numeric_core, x)
}

standardize_group <- function(x) {
  x <- as.character(x)
  x <- str_to_upper(str_trim(x))

  case_when(
    x %in% c("1", "CON", "CONTROL", "CTRL") ~ "CON",
    x %in% c("2", "RES", "RESILIENT") ~ "RES",
    x %in% c("3", "SUS", "SUSCEPTIBLE") ~ "SUS",
    x %in% c("SIS", "STRESS") ~ "SIS",
    TRUE ~ x
  )
}

# ------------------------------------------------
# 3) READ PROTEOMICS METADATA BLOCK
# ------------------------------------------------

raw <- read_excel(proteomics_file, col_names = FALSE)

first_col <- as.character(raw[[1]])

metadata_keys <- c(
  "id",
  "sample_id",
  "shortmicrogliame",
  "exclude",
  "group",
  "group2",
  "region",
  "layer",
  "celltype",
  "celltype_layer",
  "phenotypeWithinUnit",
  "ExpGroup",
  "plate",
  "sample_location",
  "sample_number",
  "AnimalID",
  "ReplicateGroup",
  "run_order",
  "region_layer_ExpGroup"
)

meta_rows <- which(first_col %in% metadata_keys)

if (length(meta_rows) == 0) {
  stop("No metadata rows detected in proteomics file.")
}

metadata_block <- raw[meta_rows, ]

sample_meta <- metadata_block %>%
  as.data.frame()

colnames(sample_meta) <- paste0("V", seq_len(ncol(sample_meta)))

sample_meta <- sample_meta %>%
  rename(meta_key = V1) %>%
  pivot_longer(
    cols = -meta_key,
    names_to = "sample_col",
    values_to = "value"
  ) %>%
  filter(sample_col != "V1") %>%
  mutate(sample_index = as.integer(str_remove(sample_col, "^V")) - 1L) %>%
  dplyr::select(sample_index, meta_key, value) %>%
  pivot_wider(names_from = meta_key, values_from = value)

sample_meta <- sample_meta %>%
  mutate(
    Sample = as.character(sample_id),
    ShortSample = as.character(shortmicrogliame),
    AnimalID_raw = as.character(AnimalID),
    AnimalID_join = standardize_animal_id(AnimalID),
    Region = as.character(region),
    Layer = as.character(layer),
    RegionLayer = if (dataset_profile %in% c("microglia", "neuron_soma")) Region else paste(Region, Layer, sep = "_"),
    SpatialUnit = if (dataset_profile %in% c("microglia", "neuron_soma")) "region" else "region_layer",
    SpatialLabel = if (dataset_profile %in% c("microglia", "neuron_soma")) Region else RegionLayer,
    CellType = as.character(celltype),
    CellTypeLayer = as.character(celltype_layer),
    ReplicateGroup = as.character(ReplicateGroup),
    Plate = as.character(plate),
    SampleLocation = as.character(sample_location),
    SampleNumber = as.character(sample_number),
    RunOrder = suppressWarnings(as.numeric(run_order)),
    RegionLayerGroup = as.character(region_layer_ExpGroup),
    ExpGroup_raw = as.character(ExpGroup),
    StressGroup_from_proteomics = standardize_group(ExpGroup_raw),
    Exclude = case_when(
      str_to_upper(as.character(exclude)) == "TRUE" ~ TRUE,
      str_to_upper(as.character(exclude)) == "FALSE" ~ FALSE,
      TRUE ~ NA
    )
  )

# ------------------------------------------------
# 4) READ AUC: WHOLE STRESS PARADIGM
# ------------------------------------------------

auc_all <- read_csv(auc_all_file, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    animal_id_join = standardize_animal_id(animal_num),
    group_auc_all = standardize_group(group),
    sex_auc_all = as.character(sex),
    batch_auc_all = as.character(batch)
  ) %>%
  filter(metric == "Movement") %>%
  dplyr::select(
    animal_id_join,
    batch_auc_all,
    group_auc_all,
    sex_auc_all,
    phase,
    change,
    window,
    auc,
    auc_norm
  ) %>%
  rename(
    Phase_all = phase,
    Change_all = change,
    Window_all = window,
    AUC_all = auc,
    AUC_norm_all = auc_norm
  )

auc_all_one <- auc_all %>%
  mutate(
    priority = case_when(
      Phase_all %in% c("Active", "active") & Window_all == "all" ~ 1,
      TRUE ~ 2
    )
  ) %>%
  arrange(animal_id_join, priority) %>%
  group_by(animal_id_join) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-priority)

# ------------------------------------------------
# 5) READ AUC: FIRST ACTIVE PHASE
# ------------------------------------------------

auc_first <- read_csv(auc_first_file, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    animal_id_join = standardize_animal_id(animal_num),
    group_auc_first = standardize_group(group),
    sex_auc_first = as.character(sex),
    batch_auc_first = as.character(batch)
  ) %>%
  filter(metric == "Movement") %>%
  dplyr::select(
    animal_id_join,
    batch_auc_first,
    group_auc_first,
    sex_auc_first,
    phase,
    change,
    window,
    auc,
    auc_norm
  ) %>%
  rename(
    Phase_firstActive = phase,
    Change_firstActive = change,
    Window_firstActive = window,
    AUC_firstActive = auc,
    AUC_norm_firstActive = auc_norm
  )

auc_first_one <- auc_first %>%
  mutate(
    priority = case_when(
      Phase_firstActive %in% c("Active", "active") & Window_firstActive == "all" ~ 1,
      TRUE ~ 2
    )
  ) %>%
  arrange(animal_id_join, priority) %>%
  group_by(animal_id_join) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(-priority)

# ------------------------------------------------
# 6) READ BEHAVIOR / PHYSIOLOGY Z-SCORES
# ------------------------------------------------

behavior <- read_excel(behavior_file, sheet = "zScore") %>%
  clean_names() %>%
  mutate(
    animal_id_join = standardize_animal_id(id),
    group_behavior_raw = as.character(group),
    sex_behavior = as.character(sex),
    batch_behavior = as.character(batch)
  ) %>%
  rename(
    NOR_z = nor,
    SucrosePref_z = sucrose_pref,
    WeightDev_z = weight_dev,
    DeltaCort_z = delta_cort,
    AdrenalWeight_z = adrenal_weight,
    SpleenWeight_z = spleen_weight,
    CombZ = comb_z
  ) %>%
  dplyr::select(
    animal_id_join,
    group_behavior_raw,
    sex_behavior,
    batch_behavior,
    NOR_z,
    SucrosePref_z,
    WeightDev_z,
    DeltaCort_z,
    AdrenalWeight_z,
    SpleenWeight_z,
    CombZ
  )

# ------------------------------------------------
# 7) MERGE
# ------------------------------------------------

merged_meta <- sample_meta %>%
  left_join(auc_all_one, by = c("AnimalID_join" = "animal_id_join")) %>%
  left_join(auc_first_one, by = c("AnimalID_join" = "animal_id_join")) %>%
  left_join(behavior, by = c("AnimalID_join" = "animal_id_join")) %>%
  mutate(
    StressGroup = case_when(
      !is.na(group_auc_all) ~ group_auc_all,
      !is.na(group_auc_first) ~ group_auc_first,
      StressGroup_from_proteomics %in% c("CON", "RES", "SUS") ~ StressGroup_from_proteomics,
      TRUE ~ NA_character_
    ),
    Sex = coalesce(sex_behavior, sex_auc_all, sex_auc_first),
    Batch = coalesce(batch_behavior, batch_auc_all, batch_auc_first, Plate),
    AnimalID = AnimalID_join
  )

# ------------------------------------------------
# 8) CLEAN OUTPUT
# ------------------------------------------------

merged_meta_clean <- merged_meta %>%
  dplyr::select(
    Sample,
    ShortSample,
    AnimalID,
    AnimalID_raw,
    Region,
    Layer,
    RegionLayer,
    SpatialUnit,
    SpatialLabel,
    CellType,
    CellTypeLayer,
    ReplicateGroup,
    StressGroup,
    Sex,
    Batch,
    CombZ,
    NOR_z,
    SucrosePref_z,
    WeightDev_z,
    DeltaCort_z,
    AdrenalWeight_z,
    SpleenWeight_z,
    AUC_all,
    AUC_norm_all,
    AUC_firstActive,
    AUC_norm_firstActive,
    RunOrder,
    Plate,
    SampleLocation,
    SampleNumber,
    Exclude
  ) %>%
  mutate(
    Region = factor(Region),
    Layer = factor(Layer),
    RegionLayer = factor(as.character(RegionLayer)),
    SpatialUnit = factor(SpatialUnit),
    SpatialLabel = factor(SpatialLabel),
    CellType = factor(CellType),
    CellTypeLayer = factor(CellTypeLayer),
    ReplicateGroup = factor(ReplicateGroup),
    StressGroup = factor(StressGroup, levels = c("CON", "RES", "SUS")),
    Sex = factor(Sex),
    Batch = factor(Batch),
    AnimalID = as.character(AnimalID),
    AnimalID_raw = as.character(AnimalID_raw),
    Exclude = as.logical(Exclude)
  ) %>%
  arrange(AnimalID, Region, Layer, ReplicateGroup)

# ------------------------------------------------
# 9) QC TABLES
# ------------------------------------------------

qc_sample_counts <- merged_meta_clean %>%
  count(RegionLayer, StressGroup, Sex, Batch, name = "n_samples") %>%
  arrange(RegionLayer, StressGroup, Sex, Batch)

qc_animal_counts <- merged_meta_clean %>%
  distinct(AnimalID, StressGroup, Sex, Batch, CombZ, AUC_all, AUC_firstActive) %>%
  count(StressGroup, Sex, Batch, name = "n_animals") %>%
  arrange(StressGroup, Sex, Batch)

qc_missingness <- merged_meta_clean %>%
  summarise(
    n_samples = n(),
    n_animals = n_distinct(AnimalID),
    missing_stress_group = sum(is.na(StressGroup)),
    missing_comb_z = sum(is.na(CombZ)),
    missing_auc_all = sum(is.na(AUC_all)),
    missing_auc_first = sum(is.na(AUC_firstActive)),
    missing_sex = sum(is.na(Sex)),
    missing_batch = sum(is.na(Batch))
  )

unmatched_animals <- merged_meta_clean %>%
  filter(
    is.na(StressGroup) |
      is.na(CombZ) |
      is.na(AUC_all) |
      is.na(AUC_firstActive)
  ) %>%
  distinct(
    AnimalID,
    AnimalID_raw,
    StressGroup,
    Sex,
    Batch,
    CombZ,
    AUC_all,
    AUC_firstActive
  ) %>%
  arrange(AnimalID)

id_overlap <- sample_meta %>%
  distinct(AnimalID_raw, AnimalID_join) %>%
  mutate(
    in_auc_all = AnimalID_join %in% auc_all_one$animal_id_join,
    in_auc_first = AnimalID_join %in% auc_first_one$animal_id_join,
    in_behavior = AnimalID_join %in% behavior$animal_id_join
  ) %>%
  arrange(AnimalID_join)

duplicate_id_check_behavior <- behavior %>%
  count(animal_id_join, name = "n") %>%
  filter(n > 1)

duplicate_id_check_proteomics <- sample_meta %>%
  distinct(AnimalID_raw, AnimalID_join) %>%
  count(AnimalID_join, name = "n") %>%
  filter(n > 1)

# ------------------------------------------------
# 10) SAVE
# ------------------------------------------------

write.xlsx(
  list(
    MergedMetadata_Clean = merged_meta_clean,
    QCSampleCounts = qc_sample_counts,
    QCAnimalCounts = qc_animal_counts,
    QCMissingness = qc_missingness,
    UnmatchedAnimals = unmatched_animals,
    IDOverlap = id_overlap,
    DuplicateID_Behavior = duplicate_id_check_behavior,
    DuplicateID_Proteomics = duplicate_id_check_proteomics
  ),
  out_file,
  overwrite = TRUE
)

if (isTRUE(write_global_compat)) {
  ensure_dir(dirname(global_compat_out_file))
  file.copy(out_file, global_compat_out_file, overwrite = TRUE)
  if (!identical(dataset_profile, "neuron_neuropil")) {
    warning(
      "Wrote legacy global compatibility metadata copy from dataset '", dataset_profile,
      "'. Dataset-aware scripts should consume dataset-scoped metadata at: ", out_file,
      call. = FALSE
    )
  }
}

cat("Done.\nResolved dataset:", dataset_profile,
    "\nSaved clean merged metadata to:\n", out_file,
    if (isTRUE(write_global_compat)) paste0("\nLegacy compatibility copy:\n", global_compat_out_file) else "",
    "\n")
