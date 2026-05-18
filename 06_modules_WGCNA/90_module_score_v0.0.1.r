# ================================================================
# General neuropil module score analysis with group-aware replicate QC
# Spatial proteomics: neuron_neuropil, nuclei excluded
# ================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(openxlsx)
library(stringr)
library(purrr)
library(tibble)
library(ggpubr)
library(scales)

select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
first <- dplyr::first

# ------------------------------------------------
# 1) PATHS
# ------------------------------------------------

analysis_primary <- "primary_all_replicates"
analysis_qc_sensitivity <- "sensitivity_flagged_replicates_removed"

analysis_display_label <- function(x) {
  case_when(
    x == analysis_primary ~ "Primary: all technical replicates",
    x == analysis_qc_sensitivity ~ "Sensitivity: flagged technical replicates removed",
    TRUE ~ x
  )
}

protein_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/morpheus/20260218_pgmatrix_imputed_neuron_neuropil_180samples_missing70pct_with_metadata.xlsx"

metadata_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/module_scores/sample_metadata_merged_clean_for_module_scores.xlsx"

mapping_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/MOUSE_10090_idmapping.dat"

saving_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Results/module_scores/general_neuropil_module_score_results/"
dir.create(saving_dir, recursive = TRUE, showWarnings = FALSE)

dir_tables <- file.path(saving_dir, "tables")
dir_qc     <- file.path(saving_dir, "qc")
dir_plots  <- file.path(saving_dir, "plots")
dir_group  <- file.path(dir_plots, "module_group_scores")
dir_cor    <- file.path(dir_plots, "behavior_correlations")
dir_qc_replicate <- file.path(dir_qc, "replicate_qc")
dir_group_primary <- file.path(dir_group, analysis_primary)
dir_group_qc <- file.path(dir_group, analysis_qc_sensitivity)
dir_directional <- file.path(dir_plots, "directional_robustness")

dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_qc, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_group, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_cor, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_qc_replicate, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_group_primary, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_group_qc, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_directional, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# 2) LOAD METADATA
# ------------------------------------------------

metadata <- read_excel(metadata_file, sheet = "MergedMetadata_Clean") %>%
  mutate(
    Sample = as.character(Sample),
    AnimalID = as.character(AnimalID),
    Region = as.character(Region),
    Layer = as.character(Layer),
    RegionLayer = paste(Region, Layer, sep = "_"),
    StressGroup = factor(StressGroup, levels = c("CON", "RES", "SUS")),
    Sex = factor(Sex),
    Batch = factor(Batch),
    ReplicateGroup = factor(ReplicateGroup),
    CellTypeLayer = as.character(CellTypeLayer),
    Exclude = as.logical(Exclude)
  ) %>%
  filter(is.na(Exclude) | Exclude == FALSE) %>%
  filter(CellTypeLayer == "neuron_neuropil")

# ------------------------------------------------
# 3) LOAD PROTEIN MATRIX FROM MORPHEUS FORMAT
# ------------------------------------------------

raw <- read_excel(protein_file, col_names = FALSE)

first_col <- as.character(raw[[1]])

metadata_keys <- c(
  "id", "sample_id", "shortmicrogliame", "exclude",
  "group", "group2", "region", "layer", "celltype",
  "celltype_layer", "phenotypeWithinUnit", "ExpGroup",
  "plate", "sample_location", "sample_number", "AnimalID",
  "ReplicateGroup", "run_order", "region_layer_ExpGroup"
)

metadata_key_rows <- which(first_col %in% metadata_keys)

if (length(metadata_key_rows) == 0) {
  stop("Could not identify the metadata block in the Morpheus-format protein matrix.")
}

protein_start <- max(metadata_key_rows) + 1

protein_df <- raw[protein_start:nrow(raw), ] %>%
  as.data.frame()

colnames(protein_df) <- c("ProteinID", as.character(raw[1, -1]))

protein_df <- protein_df %>%
  mutate(ProteinID = as.character(ProteinID)) %>%
  filter(!is.na(ProteinID), ProteinID != "") %>%
  distinct(ProteinID, .keep_all = TRUE)

sample_cols <- intersect(colnames(protein_df), metadata$Sample)

if (length(sample_cols) == 0) {
  stop("No matching sample names between protein matrix and metadata.")
}

metadata <- metadata %>%
  filter(Sample %in% sample_cols)

protein_df <- protein_df %>%
  select(ProteinID, all_of(metadata$Sample))

mat <- protein_df %>%
  as.data.frame()

rownames(mat) <- toupper(mat$ProteinID)
mat$ProteinID <- NULL
mat <- as.matrix(mat)
storage.mode(mat) <- "numeric"

# ------------------------------------------------
# 4) LOAD UNIPROT MAPPING
# ------------------------------------------------

mapping_raw <- read.delim(
  mapping_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  col.names = c("Accession", "Type", "Value")
)

mapping_table <- mapping_raw %>%
  filter(Type %in% c("UniProtKB-ID", "Gene_Name")) %>%
  mutate(
    Accession = toupper(Accession),
    Value = toupper(Value)
  ) %>%
  pivot_wider(
    names_from = Type,
    values_from = Value,
    values_fn = ~ paste(unique(.x), collapse = ";")
  ) %>%
  rename(
    MatrixID = `UniProtKB-ID`,
    GeneSymbol = Gene_Name
  ) %>%
  mutate(
    MatrixID = str_split(MatrixID, ";") %>% map_chr(1),
    GeneSymbol = str_split(GeneSymbol, ";") %>% map_chr(1)
  )

map_accessions_to_matrix_ids <- function(accessions, mapping_table) {
  mapping_table %>%
    filter(Accession %in% toupper(accessions)) %>%
    pull(MatrixID) %>%
    unique() %>%
    toupper()
}

# ------------------------------------------------
# 5) GENERAL NEUROPIL MODULES
# ------------------------------------------------

rna_ca1_accessions <- c(
  "A0A7N9VR94", "P05132", "Q64213", "Q6PFR5", "Q7TNC4",
  "Q8BMZ5", "Q8C2Q3", "Q8VH51", "Q921F4", "Q99M28",
  "Q8CCS6", "Q9CQQ8", "Q99LI7", "E9Q616", "O88569",
  "P09405", "P17225", "P26369"
)

rna_ca2sr_accessions <- c(
  "O70551", "O89086", "Q505B7", "Q810A7", "Q8BMZ5",
  "Q8BX17", "Q8JZX4", "Q921M3", "Q922X9", "A0A1W2P872",
  "O70133", "P26369", "P60824", "Q3UEB3", "Q6PDM2",
  "Q8BL97", "P62307", "P62311", "Q99K48", "E9Q616",
  "O55135", "Q8C8T8", "Q8CJG1", "Q8K1C0", "Q921F2",
  "Q62376", "Q60823"
)

rna_additions_ca3sr <- c(
  "O35737", "O08583", "O35326", "O70133",
  "P61979", "P62309", "P62996", "P70333",
  "P83870", "P62311", "Q6PDM2", "Q8BL97",
  "P09405", "Q61937", "P05132", "Q64213",
  "P32067", "P35550", "P50580", "P84104",
  "Q00PI9", "Q3UKJ7", "Q60749", "Q8R081",
  "P22366", "P35123", "P63058", "P70372",
  "Q4VGL6", "Q6ZWM4", "Q8BTM8"
)

rna_additions_ca2so <- c(
  "A0A1W2P872", "P62317", "P62320", "P84104",
  "Q64213", "Q7TNC4", "Q8C181", "Q91VC3",
  "Q921M3", "Q9D883", "O35737", "P60229",
  "P62996", "Q3TCH7", "Q3UJB9", "Q62376",
  "Q8BL97", "Q8BP71", "Q8VHM5", "Q91YQ3",
  "Q8BMZ5", "Q9Z1N5", "D3YXK2"
)

mito_accessions <- c(
  "P03899", "P19536", "Q91VD9", "Q91WD5", "Q9DC70",
  "Q9DCT2", "Q3UIU2", "Q06185", "P0DN34", "P03930",
  "P52503", "P97450", "Q7TMF3", "Q8BK30", "Q8K3J1",
  "Q91YT0", "Q99LC3", "Q9CQ54", "P43024", "P56391",
  "P99028", "Q5HZI9", "P17665", "Q9CPX8", "Q9CQ69",
  "Q9CQA3"
)

mito_additions_ca3sr <- c(
  "P03888", "P17665", "P48771", "Q61387",
  "Q9CR68", "Q9D855", "Q9DCJ5", "Q9QX60",
  "P00158", "P97478", "O09111", "O35943",
  "P00397", "P14152", "Q80X68", "Q9CQH3",
  "Q9CQY6", "Q9DCS9", "O55042", "Q06185",
  "Q99LY9", "Q9CQC7", "Q9CR61", "Q9D6J5",
  "Q9WUR9", "P00416", "Q3TXX4", "Q9D6M3",
  "Q9Z1G3", "Q9Z2B2", "P03893", "P03930",
  "P12787"
)

endo_lysosomal_accessions <- c(
  "P97412", "O35114", "O88384", "O89116", "O89053",
  "P62878", "Q5EG47", "Q8R1T1", "Q8R307", "Q9EQH3",
  "Q9R0N7", "Q9Z266", "O88668", "Q9DB50", "Q9DC29",
  "O55102", "Q9Z1T1", "Q9Z1T6", "P39054", "Q2NL51",
  "E9Q137", "O88597", "O88811", "Q5KU39", "Q6DFW0",
  "Q9CR26", "O88983", "P17047", "P62746", "P70280",
  "P11438", "P84078", "P97950", "Q09200", "Q63844",
  "Q64008"
)

chromatin_rnp_accessions <- c(
  "P06537", "P14733", "Q3TEA8", "Q60972", "Q61937",
  "Q64522", "Q6PDQ2", "Q8C050", "Q8VEK3", "Q8VIJ6",
  "Q61315", "Q99JN2", "E9Q7G0", "O54692", "O54879",
  "O88569", "P10922", "P28352", "P43274", "P43276"
)

chromatin_additions_ca3sr <- c(
  "G5E870", "Q3THW5", "P10922", "P14733", "P21619",
  "P25976", "Q61214", "Q61687", "Q64522", "Q6P5D8"
)

# -----------------------------
# Main RNA/RNP consensus module
# recurrent in >=2 RNA-related sources
# -----------------------------

rna_sources <- list(
  CA1 = rna_ca1_accessions,
  CA2_SR = rna_ca2sr_accessions,
  CA3_SR = rna_additions_ca3sr,
  CA2_SO = rna_additions_ca2so
)

rna_recurrence <- table(unlist(rna_sources))

rna_rnp_consensus <- names(rna_recurrence)[rna_recurrence >= 2]

modules_accession <- list(
  Neuropil_RNP_RNA_processing_main = unique(rna_rnp_consensus),

  Neuropil_RNP_RNA_processing_full = unique(c(
    rna_ca1_accessions,
    rna_ca2sr_accessions,
    rna_additions_ca3sr,
    rna_additions_ca2so
  )),

  Neuropil_mito_bioenergetics = unique(c(
    mito_accessions,
    mito_additions_ca3sr
  )),

  Neuropil_endo_lysosomal_proteostasis = unique(
    endo_lysosomal_accessions
  ),

  Neuropil_chromatin_RNP_related_exploratory = unique(c(
    chromatin_rnp_accessions,
    chromatin_additions_ca3sr
  ))
)

modules_mapped <- imap(
  modules_accession,
  ~ map_accessions_to_matrix_ids(.x, mapping_table)
)

modules <- imap(
  modules_mapped,
  ~ intersect(.x, rownames(mat))
)

module_size_check <- tibble(
  Module = names(modules_accession),
  n_accessions = lengths(modules_accession),
  n_mapped = lengths(modules_mapped),
  n_found = lengths(modules)
)

print(module_size_check)

write.xlsx(
  module_size_check,
  file.path(dir_qc, "module_size_check.xlsx"),
  overwrite = TRUE
)

# ------------------------------------------------
# 6) MODULE COVERAGE QC
# ------------------------------------------------

module_coverage <- imap_dfr(modules_accession, function(accessions, module_name) {
  mapped_ids <- modules_mapped[[module_name]]
  found_ids <- modules[[module_name]]

  tibble(
    Module = module_name,
    n_accessions_input = length(accessions),
    n_mapped_to_matrix_id = length(mapped_ids),
    n_found_in_matrix = length(found_ids),
    accessions_input = paste(accessions, collapse = "; "),
    mapped_matrix_ids = paste(mapped_ids, collapse = "; "),
    found_matrix_ids = paste(found_ids, collapse = "; "),
    missing_after_mapping = paste(setdiff(mapped_ids, rownames(mat)), collapse = "; ")
  )
})

write.xlsx(module_coverage, file.path(dir_qc, "module_gene_coverage.xlsx"), overwrite = TRUE)

# ------------------------------------------------
# 7) Z-SCORE PER PROTEIN
# ------------------------------------------------

mat_z <- t(scale(t(mat)))
mat_z <- mat_z[rowSums(is.finite(mat_z)) > 0, , drop = FALSE]

# ------------------------------------------------
# 8) MODULE SCORE FUNCTIONS
# ------------------------------------------------

compute_module_score <- function(mat_z, genes, min_genes = 5) {
  genes_found <- intersect(genes, rownames(mat_z))

  if (length(genes_found) < min_genes) {
    warning("Only ", length(genes_found), " proteins found for module. Returning NA.")
    return(rep(NA_real_, ncol(mat_z)))
  }

  colMeans(mat_z[genes_found, , drop = FALSE], na.rm = TRUE)
}

compute_pca_module_score <- function(mat_z, genes, min_genes = 5) {
  genes_found <- intersect(genes, rownames(mat_z))

  if (length(genes_found) < min_genes) {
    return(rep(NA_real_, ncol(mat_z)))
  }

  submat <- mat_z[genes_found, , drop = FALSE]
  submat <- submat[rowSums(is.finite(submat)) == ncol(submat), , drop = FALSE]

  if (nrow(submat) < min_genes) {
    return(rep(NA_real_, ncol(mat_z)))
  }

  pca <- prcomp(t(submat), center = FALSE, scale. = FALSE)
  score <- pca$x[, 1]
  mean_score <- colMeans(submat, na.rm = TRUE)

  if (cor(score, mean_score, use = "complete.obs") < 0) {
    score <- -score
  }

  score
}

# ------------------------------------------------
# 9) COMPUTE MODULE SCORES
# ------------------------------------------------

scores_df <- imap_dfr(modules, function(genes, module_name) {
  tibble(
    Sample = colnames(mat_z),
    Module = module_name,
    ModuleScore = compute_module_score(mat_z, genes)
  )
}) %>%
  left_join(metadata, by = "Sample")

pca_scores_df <- imap_dfr(modules, function(genes, module_name) {
  tibble(
    Sample = colnames(mat_z),
    Module = module_name,
    PC1_ModuleScore = compute_pca_module_score(mat_z, genes)
  )
}) %>%
  left_join(metadata, by = "Sample")

write.xlsx(scores_df, file.path(dir_tables, "module_scores_per_sample.xlsx"), overwrite = TRUE)
write.xlsx(pca_scores_df, file.path(dir_tables, "module_scores_pca_sensitivity.xlsx"), overwrite = TRUE)

# ------------------------------------------------
# 10) ADVANCED GROUP-AWARE REPLICATE QC
# ------------------------------------------------

library(ggrepel)

qc_palette <- c(
  OK = "grey35",
  Review = "#e63947"
)

group_palette <- c(
  CON = "#3e3c6f",
  RES = "#c6c3bb",
  SUS = "#e63a47"
)

qc_thresholds <- list(
  group_delta_moderate = 1.00,
  group_delta_strong = 1.50,
  replicate_diff_moderate = 0.75,
  replicate_diff_strong = 1.25,
  driver_score_moderate = 1.00,
  driver_score_strict = 1.50
)

qc_df <- scores_df %>%
  filter(
    is.finite(ModuleScore),
    !is.na(StressGroup),
    !is.na(RegionLayer),
    !is.na(ReplicateGroup)
  ) %>%
  mutate(
    StressGroup = factor(StressGroup, levels = c("CON", "RES", "SUS")),
    ReplicateGroup = as.character(ReplicateGroup),
    QC_ID = paste(AnimalID, RegionLayer, ReplicateGroup, Module, sep = " | ")
  )

# ------------------------------------------------
# 10.1 Within-group replicate deviation
# ------------------------------------------------

qc_df <- qc_df %>%
  group_by(RegionLayer, Module, StressGroup) %>%
  mutate(
    n_group_replicates = sum(is.finite(ModuleScore)),
    group_median = median(ModuleScore, na.rm = TRUE),
    group_mean = mean(ModuleScore, na.rm = TRUE),
    group_sd = sd(ModuleScore, na.rm = TRUE),
    group_mad = mad(ModuleScore, na.rm = TRUE),
    
    robust_group_z = case_when(
      is.finite(group_mad) & group_mad > 0 ~
        (ModuleScore - group_median) / group_mad,
      is.finite(group_sd) & group_sd > 0 ~
        (ModuleScore - group_mean) / group_sd,
      TRUE ~ 0
    ),
    
    abs_robust_group_z = abs(robust_group_z),
    group_abs_delta = abs(ModuleScore - group_median),
    group_deviation_flag = case_when(
      n_group_replicates >= 8 ~
        abs_robust_group_z >= 2.5 | group_abs_delta >= qc_thresholds$group_delta_moderate,
      n_group_replicates >= 4 ~
        abs_robust_group_z >= 2.0 | group_abs_delta >= qc_thresholds$group_delta_moderate,
      TRUE ~ FALSE
    ),
    group_strong_deviation_flag = case_when(
      n_group_replicates >= 8 ~
        abs_robust_group_z >= 3.5 | group_abs_delta >= qc_thresholds$group_delta_strong,
      n_group_replicates >= 4 ~
        abs_robust_group_z >= 3.0 | group_abs_delta >= qc_thresholds$group_delta_strong,
      TRUE ~ FALSE
    )
  ) %>%
  ungroup()

# ------------------------------------------------
# 10.2 Left/right replicate pair QC
# ------------------------------------------------

replicate_qc <- qc_df %>%
  select(
    Sample,
    ShortSample,
    AnimalID,
    StressGroup,
    Sex,
    RegionLayer,
    Module,
    ReplicateGroup,
    ModuleScore,
    group_median,
    robust_group_z,
    abs_robust_group_z
  ) %>%
  filter(ReplicateGroup %in% c("Left", "Right")) %>%
  pivot_wider(
    names_from = ReplicateGroup,
    values_from = c(
      Sample,
      ShortSample,
      ModuleScore,
      group_median,
      robust_group_z,
      abs_robust_group_z
    ),
    names_sep = "_"
  ) %>%
  mutate(
    # group median should be identical for left/right, but keep fallback
    group_median_pair = dplyr::coalesce(group_median_Left, group_median_Right),

    replicate_diff = abs(ModuleScore_Left - ModuleScore_Right),
    replicate_mean = rowMeans(cbind(ModuleScore_Left, ModuleScore_Right), na.rm = TRUE),

    left_distance_to_group = abs(ModuleScore_Left - group_median_pair),
    right_distance_to_group = abs(ModuleScore_Right - group_median_pair),

    replicate_direction = case_when(
      ModuleScore_Left > ModuleScore_Right ~ "Left > Right",
      ModuleScore_Right > ModuleScore_Left ~ "Right > Left",
      TRUE ~ "Equal"
    ),

    discordance_driver = case_when(
      abs_robust_group_z_Left > abs_robust_group_z_Right ~ "Left",
      abs_robust_group_z_Right > abs_robust_group_z_Left ~ "Right",
      TRUE ~ "Both"
    ),

    discordance_driver_score = case_when(
      discordance_driver == "Left" ~ abs_robust_group_z_Left - abs_robust_group_z_Right,
      discordance_driver == "Right" ~ abs_robust_group_z_Right - abs_robust_group_z_Left,
      TRUE ~ 0
    ),

    candidate_bad_replicate = case_when(
      discordance_driver == "Left" ~ Sample_Left,
      discordance_driver == "Right" ~ Sample_Right,
      TRUE ~ NA_character_
    ),

    candidate_bad_replicate_short = case_when(
      discordance_driver == "Left" ~ ShortSample_Left,
      discordance_driver == "Right" ~ ShortSample_Right,
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(RegionLayer, Module, StressGroup) %>%
  mutate(
    n_group_pairs = sum(is.finite(replicate_diff)),
    diff_median = median(replicate_diff, na.rm = TRUE),
    diff_mean = mean(replicate_diff, na.rm = TRUE),
    diff_sd = sd(replicate_diff, na.rm = TRUE),
    diff_mad = mad(replicate_diff, na.rm = TRUE),

    replicate_diff_z = case_when(
      is.finite(diff_mad) & diff_mad > 0 ~
        (replicate_diff - diff_median) / diff_mad,
      is.finite(diff_sd) & diff_sd > 0 ~
        (replicate_diff - diff_mean) / diff_sd,
      TRUE ~ 0
    ),

    replicate_disagreement_flag = case_when(
      n_group_pairs >= 4 ~
        replicate_diff_z >= 2.5 | replicate_diff >= qc_thresholds$replicate_diff_moderate,
      n_group_pairs >= 3 ~
        replicate_diff_z >= 1.75 | replicate_diff >= qc_thresholds$replicate_diff_moderate,
      n_group_pairs >= 2 ~
        replicate_diff >= qc_thresholds$replicate_diff_moderate,
      TRUE ~ FALSE
    ),
    replicate_strong_disagreement_flag = case_when(
      n_group_pairs >= 4 ~
        replicate_diff_z >= 3.5 | replicate_diff >= qc_thresholds$replicate_diff_strong,
      n_group_pairs >= 3 ~
        replicate_diff_z >= 2.5 | replicate_diff >= qc_thresholds$replicate_diff_strong,
      n_group_pairs >= 2 ~
        replicate_diff >= qc_thresholds$replicate_diff_strong,
      TRUE ~ FALSE
    )
  ) %>%
  ungroup() %>%
  mutate(
    candidate_single_replicate_exclusion =
      replicate_disagreement_flag &
      discordance_driver %in% c("Left", "Right") &
      discordance_driver_score >= qc_thresholds$driver_score_moderate,

    candidate_single_replicate_exclusion_strict =
      replicate_strong_disagreement_flag &
      discordance_driver %in% c("Left", "Right") &
      discordance_driver_score >= qc_thresholds$driver_score_strict,

    replicate_qc_status = case_when(
      candidate_single_replicate_exclusion_strict ~ "Candidate_exclude_strict",
      candidate_single_replicate_exclusion ~ "Candidate_exclude_moderate",
      replicate_strong_disagreement_flag ~ "Review_strong_disagreement",
      replicate_disagreement_flag ~ "Review_disagreement",
      TRUE ~ "OK"
    )
  )

# ------------------------------------------------
# 10.3 Animal-level leverage QC
# ------------------------------------------------

animal_mean_qc <- qc_df %>%
  group_by(
    AnimalID,
    StressGroup,
    Sex,
    RegionLayer,
    Module
  ) %>%
  summarise(
    Batch = {
      batch_values <- sort(unique(as.character(Batch[!is.na(Batch)])))
      if (length(batch_values) == 0) NA_character_ else paste(batch_values, collapse = ";")
    },
    animal_mean = mean(ModuleScore, na.rm = TRUE),
    animal_sd_between_replicates = sd(ModuleScore, na.rm = TRUE),
    n_replicates_for_animal = n(),
    .groups = "drop"
  ) %>%
  group_by(RegionLayer, Module, StressGroup) %>%
  mutate(
    n_animals_in_group = n_distinct(AnimalID),
    animal_group_mean = mean(animal_mean, na.rm = TRUE),
    animal_group_sd = sd(animal_mean, na.rm = TRUE),
    animal_group_median = median(animal_mean, na.rm = TRUE),
    animal_group_mad = mad(animal_mean, na.rm = TRUE),
    
    animal_leverage_z = case_when(
      is.finite(animal_group_mad) & animal_group_mad > 0 ~
        (animal_mean - animal_group_median) / animal_group_mad,
      is.finite(animal_group_sd) & animal_group_sd > 0 ~
        (animal_mean - animal_group_mean) / animal_group_sd,
      TRUE ~ 0
    ),
    
    abs_animal_leverage_z = abs(animal_leverage_z),
    animal_leverage_flag = n_animals_in_group >= 4 & abs_animal_leverage_z >= 2.0,
    animal_strong_leverage_flag = n_animals_in_group >= 4 & abs_animal_leverage_z >= 2.8
  ) %>%
  ungroup()

# ------------------------------------------------
# 10.4 Merge QC metrics
# ------------------------------------------------

sample_qc_summary <- qc_df %>%
  left_join(
    replicate_qc %>%
      select(
        AnimalID,
        RegionLayer,
        Module,
        replicate_diff,
        replicate_mean,
        replicate_direction,
        discordance_driver,
        discordance_driver_score,
        candidate_bad_replicate,
        candidate_bad_replicate_short,
        candidate_single_replicate_exclusion,
        candidate_single_replicate_exclusion_strict,
        replicate_diff_z,
        replicate_disagreement_flag,
        replicate_strong_disagreement_flag,
        replicate_qc_status,
        n_group_pairs
      ),
    by = c("AnimalID", "RegionLayer", "Module")
  ) %>%
  left_join(
    animal_mean_qc %>%
      select(
        AnimalID,
        RegionLayer,
        Module,
        animal_mean,
        animal_sd_between_replicates,
        animal_leverage_z,
        abs_animal_leverage_z,
        animal_leverage_flag,
        animal_strong_leverage_flag
      ),
    by = c("AnimalID", "RegionLayer", "Module")
  ) %>%
  mutate(
    group_deviation_flag = coalesce(group_deviation_flag, FALSE),
    group_strong_deviation_flag = coalesce(group_strong_deviation_flag, FALSE),
    replicate_disagreement_flag = coalesce(replicate_disagreement_flag, FALSE),
    replicate_strong_disagreement_flag = coalesce(replicate_strong_disagreement_flag, FALSE),
    candidate_single_replicate_exclusion = coalesce(candidate_single_replicate_exclusion, FALSE),
    candidate_single_replicate_exclusion_strict = coalesce(candidate_single_replicate_exclusion_strict, FALSE),
    replicate_qc_status = coalesce(replicate_qc_status, "OK"),
    animal_leverage_flag = coalesce(animal_leverage_flag, FALSE),
    animal_strong_leverage_flag = coalesce(animal_strong_leverage_flag, FALSE),

    flagged_single_replicate =
      Sample == candidate_bad_replicate &
      candidate_single_replicate_exclusion,
    
    flagged_single_replicate_strict =
      Sample == candidate_bad_replicate &
      candidate_single_replicate_exclusion_strict,

    flagged_single_replicate = coalesce(flagged_single_replicate, FALSE),
    flagged_single_replicate_strict = coalesce(flagged_single_replicate_strict, FALSE),
    
    n_flags =
      as.integer(group_deviation_flag) +
      as.integer(replicate_disagreement_flag) +
      as.integer(animal_leverage_flag),
    
    n_strong_flags =
      as.integer(group_strong_deviation_flag) +
      as.integer(replicate_strong_disagreement_flag) +
      as.integer(animal_strong_leverage_flag),
    
    QC_status = case_when(
      flagged_single_replicate_strict ~ "Candidate_exclude_strict",
      flagged_single_replicate ~ "Candidate_exclude_moderate",
      n_strong_flags >= 1 ~ "Review_high_priority",
      n_flags >= 2 ~ "Review_high_priority",
      n_flags == 1 ~ "Review_moderate_priority",
      TRUE ~ "OK"
    ),
    
    QC_plot_status = ifelse(QC_status == "OK", "OK", "Review"),
    
    QC_reason = case_when(
      flagged_single_replicate_strict ~
        "candidate technical replicate drives strong left/right disagreement",
      flagged_single_replicate ~
        "candidate technical replicate drives left/right disagreement",
      group_strong_deviation_flag & replicate_strong_disagreement_flag ~
        "strong group deviation + strong replicate disagreement",
      group_strong_deviation_flag ~
        "strong deviation from same-group technical replicates",
      replicate_strong_disagreement_flag ~
        "strong left/right disagreement",
      animal_strong_leverage_flag ~
        "strong animal-level leverage",
      group_deviation_flag & replicate_disagreement_flag ~
        "group deviation + replicate disagreement",
      group_deviation_flag ~
        "deviation from same-group technical replicates",
      replicate_disagreement_flag ~
        "left/right disagreement",
      animal_leverage_flag ~
        "animal-level leverage",
      TRUE ~ "none"
    ),
    
    QC_priority_score =
      coalesce(abs_robust_group_z, 0) +
      coalesce(pmax(replicate_diff_z, 0, na.rm = TRUE), 0) +
      coalesce(abs_animal_leverage_z, 0) +
      coalesce(pmax(discordance_driver_score, 0, na.rm = TRUE), 0)
  )

qc_review_priority <- sample_qc_summary %>%
  filter(QC_status != "OK") %>%
  arrange(
    desc(flagged_single_replicate_strict),
    desc(flagged_single_replicate),
    desc(n_strong_flags),
    desc(n_flags),
    desc(QC_priority_score)
  )

qc_summary_counts <- sample_qc_summary %>%
  count(RegionLayer, Module, StressGroup, QC_status, QC_reason, name = "n_samples") %>%
  arrange(RegionLayer, Module, StressGroup, QC_status)

replicate_status_counts <- replicate_qc %>%
  count(RegionLayer, Module, StressGroup, replicate_qc_status, name = "n_pairs") %>%
  arrange(RegionLayer, Module, StressGroup, replicate_qc_status)

qc_animal_review <- animal_mean_qc %>%
  mutate(
    QC_status = case_when(
      animal_strong_leverage_flag ~ "Review_high_priority",
      animal_leverage_flag ~ "Review_moderate_priority",
      TRUE ~ "OK"
    )
  ) %>%
  arrange(desc(abs_animal_leverage_z), RegionLayer, Module, StressGroup)

# ------------------------------------------------
# 10.5 Save QC tables
# ------------------------------------------------

write.xlsx(
  list(
    ReviewPriority = qc_review_priority,
    SampleLevelQC = sample_qc_summary,
    ReplicatePairQC = replicate_qc,
    AnimalLevelQC = qc_animal_review,
    SummaryCounts = qc_summary_counts,
    ReplicateStatusCounts = replicate_status_counts,
    Thresholds = tibble(
      threshold = names(qc_thresholds),
      value = unlist(qc_thresholds)
    )
  ),
  file.path(dir_qc_replicate, "00_replicate_QC_tables.xlsx"),
  overwrite = TRUE
)

write.xlsx(qc_review_priority, file.path(dir_qc_replicate, "01_QC_review_priority.xlsx"), overwrite = TRUE)
write.xlsx(sample_qc_summary, file.path(dir_qc_replicate, "02_sample_level_QC_full.xlsx"), overwrite = TRUE)
write.xlsx(replicate_qc, file.path(dir_qc_replicate, "03_replicate_pair_QC.xlsx"), overwrite = TRUE)
write.xlsx(qc_animal_review, file.path(dir_qc_replicate, "04_animal_level_QC.xlsx"), overwrite = TRUE)
write.xlsx(qc_summary_counts, file.path(dir_qc_replicate, "05_QC_summary_counts.xlsx"), overwrite = TRUE)
write.xlsx(replicate_status_counts, file.path(dir_qc_replicate, "06_replicate_status_counts.xlsx"), overwrite = TRUE)

# ------------------------------------------------
# 10.6 QC plots
# ------------------------------------------------

top_flags <- sample_qc_summary %>%
  filter(QC_status != "OK") %>%
  arrange(desc(QC_priority_score)) %>%
  slice_head(n = 40) %>%
  mutate(
    PlotLabel = paste(AnimalID, RegionLayer, Module, ReplicateGroup, sep = " | "),
    QC_status = factor(
      QC_status,
      levels = c(
        "Review_moderate_priority",
        "Review_high_priority",
        "Candidate_exclude_moderate",
        "Candidate_exclude_strict"
      )
    )
  )

if (nrow(top_flags) > 0) {
  p_top_flags <- ggplot(
    top_flags,
    aes(
      x = reorder(PlotLabel, QC_priority_score),
      y = QC_priority_score,
      fill = QC_status
    )
  ) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    coord_flip() +
    labs(
      title = "QC review priority: top flagged samples",
      x = NULL,
      y = "QC priority score"
    ) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      axis.text.y = element_text(size = 5.5),
      plot.title = element_text(size = 9, face = "bold")
    )
  
  ggsave(file.path(dir_qc_replicate, "QC_review_priority_top_flags.svg"), p_top_flags,
         width = 190, height = 145, units = "mm")
}

rep_label_df <- replicate_qc %>%
  filter(replicate_qc_status != "OK") %>%
  mutate(
    Label = paste0(
      AnimalID, "\n",
      discordance_driver, " driver\n",
      replicate_qc_status
    )
  )

p_rep <- ggplot(
  replicate_qc,
  aes(
    x = replicate_diff,
    y = discordance_driver_score,
    color = replicate_qc_status
  )
) +
  geom_point(size = 2.1, alpha = 0.9) +
  geom_vline(xintercept = median(replicate_qc$replicate_diff, na.rm = TRUE), linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = 1.5, linetype = "dashed", linewidth = 0.25) +
  ggrepel::geom_text_repel(
    data = rep_label_df,
    aes(label = Label),
    size = 2,
    max.overlaps = 100,
    segment.linewidth = 0.2,
    show.legend = FALSE
  ) +
  facet_grid(Module ~ RegionLayer, scales = "free") +
  scale_color_manual(
    values = c(
      OK = "grey35",
      Review_disagreement = "#f4a261",
      Review_strong_disagreement = "#e76f51",
      Candidate_exclude_moderate = "#e63947",
      Candidate_exclude_strict = "#8d0801"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Replicate-driver QC",
    x = "|Left - Right| module score",
    y = "driver score\n(distance imbalance)"
  ) +
  theme_classic(base_size = 7) +
  theme(
    strip.text = element_text(size = 5.8),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(size = 9, face = "bold")
  )

ggsave(file.path(dir_qc_replicate, "QC_replicate_driver_labeled.svg"), p_rep,
       width = 270, height = 190, units = "mm")

context_label_df <- sample_qc_summary %>%
  filter(QC_status != "OK") %>%
  arrange(desc(QC_priority_score)) %>%
  group_by(RegionLayer, Module) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    PlotLabel = paste0(
      AnimalID,
      " ",
      ReplicateGroup,
      "\n",
      QC_reason
    )
  )

p_context <- ggplot(
  sample_qc_summary,
  aes(
    x = StressGroup,
    y = ModuleScore
  )
) +
  geom_boxplot(
    aes(fill = StressGroup),
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.32,
    linewidth = 0.25
  ) +
  geom_point(
    aes(
      color = QC_plot_status,
      shape = ReplicateGroup
    ),
    position = position_jitter(width = 0.10, height = 0),
    size = 1.25,
    alpha = 0.82
  ) +
  ggrepel::geom_text_repel(
    data = context_label_df,
    aes(label = PlotLabel),
    size = 1.8,
    max.overlaps = 80,
    min.segment.length = 0,
    segment.linewidth = 0.18,
    box.padding = 0.25,
    point.padding = 0.12,
    show.legend = FALSE
  ) +
  facet_grid(Module ~ RegionLayer, scales = "free_y") +
  scale_fill_manual(values = group_palette, drop = FALSE) +
  scale_color_manual(values = qc_palette, drop = FALSE) +
  labs(
    title = "QC sample group context",
    x = NULL,
    y = "Module score"
  ) +
  theme_classic(base_size = 7) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.text = element_text(size = 5.8),
    plot.title = element_text(size = 9, face = "bold")
  )

ggsave(
  file.path(dir_qc_replicate, "QC_sample_group_context_labeled.svg"),
  p_context,
  width = 285,
  height = 205,
  units = "mm"
)

replicate_ba_df <- replicate_qc %>%
  mutate(
    signed_replicate_diff = ModuleScore_Left - ModuleScore_Right,
    BA_Label = paste0(
      AnimalID,
      "\n",
      discordance_driver,
      " driver\n",
      replicate_qc_status
    )
  )

replicate_ba_label_df <- replicate_ba_df %>%
  filter(replicate_qc_status != "OK")

p_ba <- ggplot(
  replicate_ba_df,
  aes(
    x = replicate_mean,
    y = signed_replicate_diff,
    color = replicate_qc_status
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  geom_point(size = 2.0, alpha = 0.88) +
  ggrepel::geom_text_repel(
    data = replicate_ba_label_df,
    aes(label = BA_Label),
    size = 2,
    max.overlaps = 100,
    segment.linewidth = 0.2,
    show.legend = FALSE
  ) +
  facet_grid(Module ~ RegionLayer, scales = "free") +
  scale_color_manual(
    values = c(
      OK = "grey35",
      Review_disagreement = "#f4a261",
      Review_strong_disagreement = "#e76f51",
      Candidate_exclude_moderate = "#e63947",
      Candidate_exclude_strict = "#8d0801"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Left/right agreement QC",
    x = "Mean of left/right module scores",
    y = "Left - Right module score"
  ) +
  theme_classic(base_size = 7) +
  theme(
    strip.text = element_text(size = 5.8),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(size = 9, face = "bold")
  )

ggsave(
  file.path(dir_qc_replicate, "QC_left_right_bland_altman_labeled.svg"),
  p_ba,
  width = 270,
  height = 190,
  units = "mm"
)

qc_burden <- sample_qc_summary %>%
  group_by(RegionLayer, Module) %>%
  summarise(
    n_samples = n(),
    n_review = sum(QC_status != "OK", na.rm = TRUE),
    n_candidate_exclude = sum(flagged_single_replicate, na.rm = TRUE),
    pct_review = 100 * n_review / n_samples,
    pct_candidate_exclude = 100 * n_candidate_exclude / n_samples,
    .groups = "drop"
  )

p_burden <- ggplot(
  qc_burden,
  aes(
    x = RegionLayer,
    y = Module,
    fill = pct_review
  )
) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = paste0(n_review, "/", n_samples)),
    size = 2.1
  ) +
  scale_fill_gradient(
    low = "grey95",
    high = "#e63947",
    na.value = "grey90"
  ) +
  labs(
    title = "QC review burden",
    x = NULL,
    y = NULL,
    fill = "% review"
  ) +
  theme_classic(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 9, face = "bold")
  )

ggsave(
  file.path(dir_qc_replicate, "QC_review_burden_heatmap.svg"),
  p_burden,
  width = 160,
  height = 90,
  units = "mm"
)

# ------------------------------------------------
# 11) ANALYSIS DATASETS: PRIMARY ALL vs QC SENSITIVITY
# ------------------------------------------------

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

make_scores_animal <- function(df, analysis_label) {
  optional_mean_vars <- c(
    "CombZ", "NOR_z", "SucrosePref_z", "WeightDev_z", "DeltaCort_z",
    "AdrenalWeight_z", "SpleenWeight_z",
    "AUC_all", "AUC_norm_all", "AUC_firstActive", "AUC_norm_firstActive"
  )
  missing_optional_vars <- setdiff(optional_mean_vars, colnames(df))
  if (length(missing_optional_vars) > 0) {
    df[missing_optional_vars] <- NA_real_
  }

  df %>%
    group_by(
      AnimalID,
      Region,
      Layer,
      RegionLayer,
      CellType,
      CellTypeLayer,
      StressGroup,
      Sex,
      Module
    ) %>%
    summarise(
      Batch = {
        batch_values <- sort(unique(as.character(Batch[!is.na(Batch)])))
        if (length(batch_values) == 0) NA_character_ else paste(batch_values, collapse = ";")
      },
      n_replicates = sum(is.finite(ModuleScore)),
      n_samples_averaged = n(),
      ModuleScore = safe_mean(ModuleScore),
      CombZ = safe_mean(CombZ),
      NOR_z = safe_mean(NOR_z),
      SucrosePref_z = safe_mean(SucrosePref_z),
      WeightDev_z = safe_mean(WeightDev_z),
      DeltaCort_z = safe_mean(DeltaCort_z),
      AdrenalWeight_z = safe_mean(AdrenalWeight_z),
      SpleenWeight_z = safe_mean(SpleenWeight_z),
      AUC_all = safe_mean(AUC_all),
      AUC_norm_all = safe_mean(AUC_norm_all),
      AUC_firstActive = safe_mean(AUC_firstActive),
      AUC_norm_firstActive = safe_mean(AUC_norm_firstActive),
      .groups = "drop"
    ) %>%
    mutate(Analysis = analysis_label)
}

scores_animal_all <- make_scores_animal(scores_df, analysis_primary)

scores_df_qc_sensitivity <- sample_qc_summary %>%
  filter(!flagged_single_replicate) %>%
  select(all_of(colnames(scores_df)))

scores_animal_qc_sensitivity <- make_scores_animal(
  scores_df_qc_sensitivity,
  analysis_qc_sensitivity
)

proteomics_for_behavior <- scores_animal_all %>%
  select(
    AnimalID,
    StressGroup,
    Sex,
    RegionLayer,
    Module,
    ModuleScore
  ) %>%
  mutate(
    ProteomicFeature = paste(RegionLayer, Module, sep = "__")
  ) %>%
  select(AnimalID, StressGroup, Sex, ProteomicFeature, ModuleScore) %>%
  pivot_wider(
    names_from = ProteomicFeature,
    values_from = ModuleScore
  ) %>%
  rename(AnimalNum = AnimalID)

dir.create("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/proteomics", recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  proteomics_for_behavior,
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/proteomics/module_scores.csv"
)

write.xlsx(scores_animal_all, file.path(dir_tables, paste0("module_scores_per_animal_regionlayer_", analysis_primary, ".xlsx")), overwrite = TRUE)
write.xlsx(scores_animal_qc_sensitivity, file.path(dir_tables, paste0("module_scores_per_animal_regionlayer_", analysis_qc_sensitivity, ".xlsx")), overwrite = TRUE)
write.xlsx(scores_df_qc_sensitivity, file.path(dir_tables, paste0("module_scores_per_sample_", analysis_qc_sensitivity, ".xlsx")), overwrite = TRUE)

# ------------------------------------------------
# 12) PARAMETRIC AND NONPARAMETRIC GROUP STATISTICS
# ------------------------------------------------

signif_label <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    p <= 0.10  ~ ".",
    TRUE ~ "ns"
  )
}

empty_parametric_result <- function(message) {
  list(
    anova = tibble(
      Effect = character(),
      Df = numeric(),
      `Sum Sq` = numeric(),
      `Mean Sq` = numeric(),
      `F value` = numeric(),
      p.value = numeric(),
      model_type = character(),
      model_note = character(),
      p.signif = character()
    ),
    contrasts = tibble(
      contrast = character(),
      group1 = character(),
      group2 = character(),
      comparison = character(),
      estimate = numeric(),
      SE = numeric(),
      df = numeric(),
      t.ratio = numeric(),
      p.value = numeric(),
      model_type = character(),
      model_note = character(),
      p_adj_within_model_BH = numeric(),
      p.signif = character()
    ),
    emmeans = tibble(
      StressGroup = factor(levels = c("CON", "RES", "SUS")),
      emmean = numeric(),
      SE = numeric(),
      df = numeric(),
      lower.CL = numeric(),
      upper.CL = numeric(),
      model_type = character(),
      model_note = character()
    ),
    diagnostics = tibble(
      n = 0,
      n_groups = 0,
      model_type = NA_character_,
      model_note = message,
      residual_shapiro_p = NA_real_,
      bartlett_p = NA_real_,
      max_cooks_distance = NA_real_
    )
  )
}

fit_parametric_model <- function(df) {
  model <- lm(ModuleScore ~ StressGroup, data = df)

  list(
    model = model,
    model_type = "lm_stress_only",
    model_note = "Batch intentionally not included"
  )
}

run_parametric_stats <- function(df) {
  df <- df %>%
    filter(is.finite(ModuleScore), !is.na(StressGroup)) %>%
    mutate(
      StressGroup = factor(StressGroup, levels = c("CON", "RES", "SUS"))
    )

  group_counts <- df %>%
    count(StressGroup, name = "n") %>%
    filter(n > 0)

  if (nrow(group_counts) < 2 || nrow(df) < 6 || any(group_counts$n < 2)) {
    return(empty_parametric_result("Skipped: fewer than two groups, n < 6, or any group n < 2"))
  }

  fit <- fit_parametric_model(df)
  model <- fit$model

  emm <- emmeans(model, ~ StressGroup)

  emm_df <- confint(emm) %>%
    as.data.frame() %>%
    mutate(
      model_type = fit$model_type,
      model_note = fit$model_note
    )

  contrast_df <- contrast(emm, method = "pairwise", adjust = "BH") %>%
    as.data.frame() %>%
    separate(contrast, into = c("group1", "group2"), sep = " - ", remove = FALSE) %>%
    mutate(
      comparison = paste(group1, group2, sep = "_vs_"),
      model_type = fit$model_type,
      model_note = fit$model_note,
      p_adj_within_model_BH = p.value,
      p.signif = signif_label(p.value)
    )

  anova_df <- anova(model) %>%
    as.data.frame() %>%
    rownames_to_column("Effect") %>%
    rename(p.value = `Pr(>F)`) %>%
    mutate(
      model_type = fit$model_type,
      model_note = fit$model_note,
      p.signif = signif_label(p.value)
    )

  residuals_model <- residuals(model)
  residual_shapiro_p <- if (length(residuals_model) >= 3 && length(residuals_model) <= 5000) {
    tryCatch(shapiro.test(residuals_model)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }

  bartlett_p <- tryCatch(
    bartlett.test(ModuleScore ~ StressGroup, data = df)$p.value,
    error = function(e) NA_real_
  )

  cooks_values <- suppressWarnings(cooks.distance(model))
  max_cooks <- if (any(is.finite(cooks_values))) max(cooks_values, na.rm = TRUE) else NA_real_

  diagnostics <- tibble(
    n = nrow(df),
    n_groups = nrow(group_counts),
    n_CON = group_counts$n[match("CON", as.character(group_counts$StressGroup))],
    n_RES = group_counts$n[match("RES", as.character(group_counts$StressGroup))],
    n_SUS = group_counts$n[match("SUS", as.character(group_counts$StressGroup))],
    model_type = fit$model_type,
    model_note = fit$model_note,
    residual_shapiro_p = residual_shapiro_p,
    bartlett_p = bartlett_p,
    max_cooks_distance = max_cooks
  ) %>%
    mutate(across(starts_with("n_"), ~ replace_na(.x, 0L)))

  list(anova = anova_df, contrasts = contrast_df, emmeans = emm_df, diagnostics = diagnostics)
}

rank_biserial <- function(x, g1, g2) {
  x1 <- x$ModuleScore[x$StressGroup == g1]
  x2 <- x$ModuleScore[x$StressGroup == g2]
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]

  if (length(x1) < 1 || length(x2) < 1) return(NA_real_)

  diffs <- outer(x2, x1, "-")
  mean(sign(diffs))
}

run_nonparametric_stats <- function(df) {
  df <- df %>%
    filter(is.finite(ModuleScore), !is.na(StressGroup)) %>%
    mutate(StressGroup = factor(StressGroup, levels = c("CON", "RES", "SUS")))

  group_counts <- df %>%
    count(StressGroup, name = "n") %>%
    filter(n > 0)

  if (nrow(group_counts) < 2 || nrow(df) < 5) {
    return(list(
      omnibus = tibble(test = "Kruskal-Wallis", statistic = NA_real_, p.value = NA_real_, p.signif = NA_character_),
      contrasts = tibble(
        contrast = character(),
        group1 = character(),
        group2 = character(),
        comparison = character(),
        p.value = numeric(),
        rank_biserial = numeric(),
        p_adj_within_model_BH = numeric(),
        p.signif = character()
      )
    ))
  }

  kw <- tryCatch(kruskal.test(ModuleScore ~ StressGroup, data = df), error = function(e) NULL)

  omnibus <- tibble(
    test = "Kruskal-Wallis",
    statistic = if (is.null(kw)) NA_real_ else unname(kw$statistic),
    p.value = if (is.null(kw)) NA_real_ else kw$p.value,
    p.signif = signif_label(p.value)
  )

  pairwise <- combn(as.character(group_counts$StressGroup), 2, simplify = FALSE) %>%
    map_dfr(function(pair) {
      g1 <- pair[1]
      g2 <- pair[2]
      tmp <- df %>% filter(StressGroup %in% c(g1, g2))

      wt <- tryCatch(
        wilcox.test(ModuleScore ~ StressGroup, data = tmp, exact = FALSE),
        error = function(e) NULL
      )

      tibble(
        contrast = paste(g1, g2, sep = " - "),
        group1 = g1,
        group2 = g2,
        comparison = paste(g1, g2, sep = "_vs_"),
        p.value = if (is.null(wt)) NA_real_ else wt$p.value,
        rank_biserial = rank_biserial(df, g1, g2)
      )
    }) %>%
    mutate(
      p_adj_within_model_BH = p.adjust(p.value, method = "BH"),
      p.signif = signif_label(p_adj_within_model_BH)
    )

  list(omnibus = omnibus, contrasts = pairwise)
}

make_parametric_stats_tables <- function(df, analysis_label) {
  stats_nested <- df %>%
    group_by(RegionLayer, Module) %>%
    nest() %>%
    mutate(stats = map(data, run_parametric_stats))

  anova_out <- stats_nested %>%
    transmute(RegionLayer, Module, anova = map(stats, "anova")) %>%
    unnest(anova) %>%
    mutate(Analysis = analysis_label, .before = 1)

  contrasts_out <- stats_nested %>%
    transmute(RegionLayer, Module, contrasts = map(stats, "contrasts")) %>%
    unnest(contrasts) %>%
    mutate(Analysis = analysis_label, .before = 1)

  emmeans_out <- stats_nested %>%
    transmute(RegionLayer, Module, emmeans = map(stats, "emmeans")) %>%
    unnest(emmeans) %>%
    mutate(Analysis = analysis_label, .before = 1)

  diagnostics_out <- stats_nested %>%
    transmute(RegionLayer, Module, diagnostics = map(stats, "diagnostics")) %>%
    unnest(diagnostics) %>%
    mutate(Analysis = analysis_label, .before = 1)

  list(anova = anova_out, contrasts = contrasts_out, emmeans = emmeans_out, diagnostics = diagnostics_out)
}

make_nonparametric_stats_tables <- function(df, analysis_label) {
  stats_nested <- df %>%
    group_by(RegionLayer, Module) %>%
    nest() %>%
    mutate(stats = map(data, run_nonparametric_stats))

  omnibus_out <- stats_nested %>%
    transmute(RegionLayer, Module, omnibus = map(stats, "omnibus")) %>%
    unnest(omnibus) %>%
    mutate(Analysis = analysis_label, .before = 1)

  contrasts_out <- stats_nested %>%
    transmute(RegionLayer, Module, contrasts = map(stats, "contrasts")) %>%
    unnest(contrasts) %>%
    mutate(Analysis = analysis_label, .before = 1)

  list(omnibus = omnibus_out, contrasts = contrasts_out)
}

stats_parametric_all <- make_parametric_stats_tables(scores_animal_all, analysis_primary)
stats_parametric_qc <- make_parametric_stats_tables(
  scores_animal_qc_sensitivity,
  analysis_qc_sensitivity
)

stats_nonparametric_all <- make_nonparametric_stats_tables(scores_animal_all, analysis_primary)
stats_nonparametric_qc <- make_nonparametric_stats_tables(
  scores_animal_qc_sensitivity,
  analysis_qc_sensitivity
)

anova_out <- bind_rows(stats_parametric_all$anova, stats_parametric_qc$anova) %>%
  group_by(Analysis, Effect) %>%
  mutate(p_adj_global_BH = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(p.signif.global = signif_label(p_adj_global_BH))

contrasts_out <- bind_rows(stats_parametric_all$contrasts, stats_parametric_qc$contrasts) %>%
  group_by(Analysis) %>%
  mutate(p_adj_global_BH = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(p.signif.global = signif_label(p_adj_global_BH))

diagnostics_out <- bind_rows(stats_parametric_all$diagnostics, stats_parametric_qc$diagnostics)

emmeans_out <- bind_rows(stats_parametric_all$emmeans, stats_parametric_qc$emmeans)

nonparametric_omnibus_out <- bind_rows(
  stats_nonparametric_all$omnibus,
  stats_nonparametric_qc$omnibus
) %>%
  group_by(Analysis) %>%
  mutate(p_adj_global_BH = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(p.signif.global = signif_label(p_adj_global_BH))

nonparametric_contrasts_out <- bind_rows(
  stats_nonparametric_all$contrasts,
  stats_nonparametric_qc$contrasts
) %>%
  group_by(Analysis) %>%
  mutate(p_adj_global_BH = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  mutate(p.signif.global = signif_label(p_adj_global_BH))

# ------------------------------------------------
# 13) EFFECT SIZES AND GROUP SUMMARIES
# ------------------------------------------------

cohens_d <- function(x, g1, g2) {
  x1 <- x$ModuleScore[x$StressGroup == g1]
  x2 <- x$ModuleScore[x$StressGroup == g2]
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  n1 <- length(x1)
  n2 <- length(x2)

  if (n1 < 2 || n2 < 2) return(NA_real_)

  pooled_sd <- sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / (n1 + n2 - 2))

  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)

  (mean(x2) - mean(x1)) / pooled_sd
}

make_effect_sizes <- function(df, label) {
  df %>%
    group_by(RegionLayer, Module) %>%
    group_modify(~ tibble(
      Contrast = c("RES-CON", "SUS-CON", "SUS-RES"),
      group1 = c("CON", "CON", "RES"),
      group2 = c("RES", "SUS", "SUS"),
      Cohen_d = c(
        cohens_d(.x, "CON", "RES"),
        cohens_d(.x, "CON", "SUS"),
        cohens_d(.x, "RES", "SUS")
      ),
      rank_biserial = c(
        rank_biserial(.x, "CON", "RES"),
        rank_biserial(.x, "CON", "SUS"),
        rank_biserial(.x, "RES", "SUS")
      )
    )) %>%
    ungroup() %>%
    mutate(Analysis = label, .before = 1)
}

effect_sizes <- bind_rows(
  make_effect_sizes(scores_animal_all, analysis_primary),
  make_effect_sizes(scores_animal_qc_sensitivity, analysis_qc_sensitivity)
)

contrast_effects <- contrasts_out %>%
  left_join(
    effect_sizes,
    by = c("Analysis", "RegionLayer", "Module", "group1", "group2")
  )

summary_group <- bind_rows(scores_animal_all, scores_animal_qc_sensitivity) %>%
  group_by(Analysis, RegionLayer, Module, StressGroup) %>%
  summarise(
    n = sum(is.finite(ModuleScore)),
    mean = mean(ModuleScore, na.rm = TRUE),
    sd = sd(ModuleScore, na.rm = TRUE),
    sem = sd / sqrt(n),
    median = median(ModuleScore, na.rm = TRUE),
    q25 = quantile(ModuleScore, 0.25, na.rm = TRUE),
    q75 = quantile(ModuleScore, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------
# 13.1 DIRECTIONAL ROBUSTNESS AND CONVERGENCE
# ------------------------------------------------

signed_z_from_p <- function(p, direction) {
  ifelse(
    is.finite(p) & is.finite(direction) & direction != 0,
    qnorm(pmin(pmax(1 - p / 2, .Machine$double.eps), 1 - .Machine$double.eps)) * sign(direction),
    NA_real_
  )
}

combine_p_fisher <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (length(p) == 0) return(NA_real_)
  pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}

combine_p_stouffer <- function(z) {
  z <- z[is.finite(z)]
  if (length(z) == 0) return(NA_real_)
  2 * pnorm(-abs(sum(z) / sqrt(length(z))))
}

direction_label <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x > 0 ~ "group2_higher",
    x < 0 ~ "group1_higher",
    TRUE ~ "no_change"
  )
}

directional_effects <- contrast_effects %>%
  mutate(
    contrast_label = paste(group2, "vs", group1),
    estimate_group2_minus_group1 = -estimate,
    estimate_se = SE,
    estimate_df = df,
    estimate_ci_half_width = ifelse(
      is.finite(estimate_se) & is.finite(estimate_df),
      qt(0.975, estimate_df) * estimate_se,
      NA_real_
    ),
    estimate_ci_low = estimate_group2_minus_group1 - estimate_ci_half_width,
    estimate_ci_high = estimate_group2_minus_group1 + estimate_ci_half_width,
    direction = direction_label(Cohen_d),
    signed_z = signed_z_from_p(p.value, Cohen_d),
    p_nominal = p.value,
    p_within_BH = p_adj_within_model_BH,
    p_global_BH = p_adj_global_BH,
    nominal_significant = is.finite(p_nominal) & p_nominal <= 0.05,
    within_BH_significant = is.finite(p_within_BH) & p_within_BH <= 0.05,
    global_BH_significant = is.finite(p_global_BH) & p_global_BH <= 0.05,
    ci_preserves_direction = case_when(
      Cohen_d > 0 ~ estimate_ci_low > 0,
      Cohen_d < 0 ~ estimate_ci_high < 0,
      TRUE ~ FALSE
    ),
    trend_same_direction = is.finite(Cohen_d) & abs(Cohen_d) >= 0.20,
    robustness_note = case_when(
      within_BH_significant & ci_preserves_direction ~ "significant_consistent_CI",
      within_BH_significant ~ "significant_directional",
      trend_same_direction ~ "directional_trend",
      TRUE ~ "weak_or_unclear"
    )
  )

summarise_directional_consistency <- function(df) {
  df %>%
    filter(is.finite(Cohen_d), direction != "no_change") %>%
    summarise(
      n_tests = n(),
      n_positive = sum(Cohen_d > 0, na.rm = TRUE),
      n_negative = sum(Cohen_d < 0, na.rm = TRUE),
      dominant_direction = case_when(
        n_positive > n_negative ~ "group2_higher",
        n_negative > n_positive ~ "group1_higher",
        n_tests > 0 ~ "mixed",
        TRUE ~ NA_character_
      ),
      n_dominant = pmax(n_positive, n_negative),
      direction_consistency = ifelse(n_tests > 0, n_dominant / n_tests, NA_real_),
      mean_Cohen_d = mean(Cohen_d, na.rm = TRUE),
      median_Cohen_d = median(Cohen_d, na.rm = TRUE),
      mean_abs_Cohen_d = mean(abs(Cohen_d), na.rm = TRUE),
      n_nominal_sig = sum(nominal_significant, na.rm = TRUE),
      n_within_BH_sig = sum(within_BH_significant, na.rm = TRUE),
      n_global_BH_sig = sum(global_BH_significant, na.rm = TRUE),
      n_directional_trends = sum(trend_same_direction, na.rm = TRUE),
      n_ci_preserved = sum(ci_preserves_direction, na.rm = TRUE),
      fisher_p = combine_p_fisher(p_nominal),
      stouffer_signed_z = ifelse(
        sum(is.finite(signed_z)) > 0,
        sum(signed_z, na.rm = TRUE) / sqrt(sum(is.finite(signed_z))),
        NA_real_
      ),
      stouffer_p = combine_p_stouffer(signed_z),
      evidence_class = case_when(
        direction_consistency >= 0.80 & n_within_BH_sig > 0 & n_directional_trends == n_tests ~
          "convergent_with_variable_power",
        direction_consistency >= 0.80 & n_directional_trends >= ceiling(0.70 * n_tests) ~
          "directionally_stable_trend",
        direction_consistency <= 0.60 & n_tests >= 3 ~
          "mixed_or_unstable_direction",
        TRUE ~
          "limited_evidence"
      ),
      .groups = "drop"
    )
}

directional_consistency_by_module <- directional_effects %>%
  group_by(Analysis, Module, group1, group2, contrast_label) %>%
  group_modify(~ summarise_directional_consistency(.x)) %>%
  ungroup()

directional_consistency_by_region <- directional_effects %>%
  group_by(Analysis, RegionLayer, group1, group2, contrast_label) %>%
  group_modify(~ summarise_directional_consistency(.x)) %>%
  ungroup()

directional_consistency_overall <- directional_effects %>%
  group_by(Analysis, group1, group2, contrast_label) %>%
  group_modify(~ summarise_directional_consistency(.x)) %>%
  ungroup()

directional_region_module_matrix <- directional_effects %>%
  group_by(Analysis, RegionLayer, Module, group1, group2, contrast_label) %>%
  summarise(
    Cohen_d = first(Cohen_d),
    direction = first(direction),
    p_nominal = first(p_nominal),
    p_within_BH = first(p_within_BH),
    p_global_BH = first(p_global_BH),
    nominal_significant = first(nominal_significant),
    within_BH_significant = first(within_BH_significant),
    global_BH_significant = first(global_BH_significant),
    ci_preserves_direction = first(ci_preserves_direction),
    robustness_note = first(robustness_note),
    .groups = "drop"
  )

run_batch_alignment <- function(df) {
  if (!"Batch" %in% colnames(df)) {
    return(tibble(
      n = sum(is.finite(df$ModuleScore)),
      n_batches = 0,
      batch_p = NA_real_,
      batch_eta2 = NA_real_,
      batch_alignment_flag = FALSE,
      batch_note = "not tested: Batch column not present"
    ))
  }

  tmp <- df %>%
    filter(is.finite(ModuleScore), !is.na(Batch)) %>%
    mutate(Batch = droplevels(factor(Batch)))

  n_batches <- n_distinct(tmp$Batch)

  if (nrow(tmp) < 6 || n_batches < 2) {
    return(tibble(
      n = nrow(tmp),
      n_batches = n_batches,
      batch_p = NA_real_,
      batch_eta2 = NA_real_,
      batch_alignment_flag = FALSE,
      batch_note = "not tested: n < 6 or fewer than two batches"
    ))
  }

  fit <- tryCatch(lm(ModuleScore ~ Batch, data = tmp), error = function(e) NULL)

  if (is.null(fit)) {
    return(tibble(
      n = nrow(tmp),
      n_batches = n_batches,
      batch_p = NA_real_,
      batch_eta2 = NA_real_,
      batch_alignment_flag = FALSE,
      batch_note = "not tested: batch model failed"
    ))
  }

  aov_df <- anova(fit)
  ss_total <- sum(aov_df$`Sum Sq`, na.rm = TRUE)
  ss_batch <- aov_df$`Sum Sq`[1]
  eta2 <- ifelse(is.finite(ss_total) & ss_total > 0, ss_batch / ss_total, NA_real_)
  batch_p <- aov_df$`Pr(>F)`[1]
  batch_flag <- is.finite(batch_p) & is.finite(eta2) & batch_p <= 0.05 & eta2 >= 0.25

  tibble(
    n = nrow(tmp),
    n_batches = n_batches,
    batch_p = batch_p,
    batch_eta2 = eta2,
    batch_alignment_flag = batch_flag,
    batch_note = case_when(
      batch_flag ~ "review: strong batch alignment",
      is.finite(eta2) & eta2 >= 0.25 ~ "review: large batch effect size",
      TRUE ~ "no strong batch alignment"
    )
  )
}

batch_alignment_diagnostics <- bind_rows(scores_animal_all, scores_animal_qc_sensitivity) %>%
  group_by(Analysis, RegionLayer, Module) %>%
  group_modify(~ run_batch_alignment(.x)) %>%
  ungroup() %>%
  group_by(Analysis) %>%
  mutate(batch_p_adj_BH = p.adjust(batch_p, method = "BH")) %>%
  ungroup()

run_module_correlation_structure <- function(df) {
  wide <- df %>%
    filter(is.finite(ModuleScore)) %>%
    select(AnimalID, Module, ModuleScore) %>%
    distinct(AnimalID, Module, .keep_all = TRUE) %>%
    pivot_wider(names_from = Module, values_from = ModuleScore)

  module_cols <- setdiff(colnames(wide), "AnimalID")

  if (length(module_cols) < 2) {
    return(tibble(
      module1 = character(),
      module2 = character(),
      n = integer(),
      rho = numeric(),
      p = numeric()
    ))
  }

  combn(module_cols, 2, simplify = FALSE) %>%
    map_dfr(function(pair) {
      tmp <- wide %>%
        filter(is.finite(.data[[pair[1]]]), is.finite(.data[[pair[2]]]))

      if (nrow(tmp) < 5) {
        return(tibble(
          module1 = pair[1],
          module2 = pair[2],
          n = nrow(tmp),
          rho = NA_real_,
          p = NA_real_
        ))
      }

      ct <- tryCatch(
        cor.test(tmp[[pair[1]]], tmp[[pair[2]]], method = "spearman", exact = FALSE),
        error = function(e) NULL
      )

      tibble(
        module1 = pair[1],
        module2 = pair[2],
        n = nrow(tmp),
        rho = if (is.null(ct)) NA_real_ else unname(ct$estimate),
        p = if (is.null(ct)) NA_real_ else ct$p.value
      )
    })
}

module_correlation_structure <- bind_rows(scores_animal_all, scores_animal_qc_sensitivity) %>%
  group_by(Analysis, RegionLayer) %>%
  group_modify(~ run_module_correlation_structure(.x)) %>%
  ungroup() %>%
  group_by(Analysis) %>%
  mutate(
    p_adj_BH = p.adjust(p, method = "BH"),
    correlation_note = case_when(
      is.finite(rho) & abs(rho) >= 0.70 & p_adj_BH <= 0.05 ~ "strong_module_coherence",
      is.finite(rho) & abs(rho) >= 0.50 ~ "moderate_module_coherence",
      TRUE ~ "limited_or_no_coherence"
    )
  ) %>%
  ungroup()

make_module_member_animal_scores <- function(sample_ids, analysis_label) {
  sample_ids <- intersect(sample_ids, colnames(mat_z))

  imap_dfr(modules, function(genes, module_name) {
    genes_found <- intersect(genes, rownames(mat_z))

    if (length(genes_found) == 0 || length(sample_ids) == 0) {
      return(tibble())
    }

    as.data.frame(t(mat_z[genes_found, sample_ids, drop = FALSE]), check.names = FALSE) %>%
      rownames_to_column("Sample") %>%
      pivot_longer(
        cols = all_of(genes_found),
        names_to = "Protein",
        values_to = "ProteinScore"
      ) %>%
      mutate(Module = module_name)
  }) %>%
    left_join(
      metadata %>%
        select(Sample, AnimalID, RegionLayer, StressGroup),
      by = "Sample"
    ) %>%
    filter(!is.na(AnimalID), !is.na(RegionLayer), !is.na(StressGroup)) %>%
    group_by(AnimalID, RegionLayer, StressGroup, Module, Protein) %>%
    summarise(
      ProteinScore = safe_mean(ProteinScore),
      .groups = "drop"
    ) %>%
    mutate(Analysis = analysis_label, .before = 1)
}

protein_group_difference <- function(df, g1, g2) {
  x1 <- df$ProteinScore[df$StressGroup == g1]
  x2 <- df$ProteinScore[df$StressGroup == g2]
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]

  p <- if (length(x1) >= 2 && length(x2) >= 2) {
    tryCatch(t.test(x2, x1)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }

  tibble(
    group1 = g1,
    group2 = g2,
    contrast_label = paste(g2, "vs", g1),
    n_group1 = length(x1),
    n_group2 = length(x2),
    protein_delta_group2_minus_group1 = mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE),
    protein_p = p
  )
}

make_module_member_effects <- function(df) {
  df %>%
    group_by(Analysis, RegionLayer, Module, Protein) %>%
    group_modify(~ bind_rows(
      protein_group_difference(.x, "CON", "RES"),
      protein_group_difference(.x, "CON", "SUS"),
      protein_group_difference(.x, "RES", "SUS")
    )) %>%
    ungroup() %>%
    group_by(Analysis) %>%
    mutate(protein_p_adj_BH = p.adjust(protein_p, method = "BH")) %>%
    ungroup()
}

module_member_animal_scores <- bind_rows(
  make_module_member_animal_scores(scores_df$Sample, analysis_primary),
  make_module_member_animal_scores(
    scores_df_qc_sensitivity$Sample,
    analysis_qc_sensitivity
  )
)

module_member_effects <- make_module_member_effects(module_member_animal_scores)

module_driver_consistency <- module_member_effects %>%
  left_join(
    directional_effects %>%
      select(
        Analysis,
        RegionLayer,
        Module,
        group1,
        group2,
        contrast_label,
        module_Cohen_d = Cohen_d,
        module_direction = direction,
        module_p_within_BH = p_within_BH
      ),
    by = c("Analysis", "RegionLayer", "Module", "group1", "group2", "contrast_label")
  ) %>%
  mutate(
    protein_same_direction_as_module = case_when(
      is.finite(module_Cohen_d) & module_Cohen_d > 0 ~ protein_delta_group2_minus_group1 > 0,
      is.finite(module_Cohen_d) & module_Cohen_d < 0 ~ protein_delta_group2_minus_group1 < 0,
      TRUE ~ NA
    )
  ) %>%
  group_by(Analysis, RegionLayer, Module, group1, group2, contrast_label) %>%
  summarise(
    module_Cohen_d = first(module_Cohen_d),
    module_direction = first(module_direction),
    module_p_within_BH = first(module_p_within_BH),
    n_module_proteins_tested = sum(is.finite(protein_delta_group2_minus_group1)),
    median_protein_delta = median(protein_delta_group2_minus_group1, na.rm = TRUE),
    mean_abs_protein_delta = mean(abs(protein_delta_group2_minus_group1), na.rm = TRUE),
    fraction_proteins_same_direction_as_module = mean(protein_same_direction_as_module, na.rm = TRUE),
    n_nominal_protein_drivers = sum(protein_p <= 0.05, na.rm = TRUE),
    n_BH_protein_drivers = sum(protein_p_adj_BH <= 0.05, na.rm = TRUE),
    driver_consistency_note = case_when(
      fraction_proteins_same_direction_as_module >= 0.70 ~ "coherent_module_member_direction",
      fraction_proteins_same_direction_as_module <= 0.45 ~ "review_unstable_member_direction",
      TRUE ~ "mixed_module_member_direction"
    ),
    .groups = "drop"
  )

write.xlsx(
  list(
    Scores_sample_primary = scores_df,
    Scores_sample_sensitivity = scores_df_qc_sensitivity,
    Scores_animal_primary = scores_animal_all,
    Scores_animal_sensitivity = scores_animal_qc_sensitivity,
    GroupSummary = summary_group,
    Parametric_ANOVA = anova_out,
    Parametric_Contrasts = contrasts_out,
    Parametric_EMMeans = emmeans_out,
    Parametric_Diagnostics = diagnostics_out,
    Nonparametric_Omnibus = nonparametric_omnibus_out,
    Nonparametric_Contrasts = nonparametric_contrasts_out,
    EffectSizes = effect_sizes,
    ParametricContrastsWithEffects = contrast_effects,
    DirectionalEffects = directional_effects,
    DirectionalConsistencyByModule = directional_consistency_by_module,
    DirectionalConsistencyByRegion = directional_consistency_by_region,
    DirectionalConsistencyOverall = directional_consistency_overall,
    DirectionalRegionModuleMatrix = directional_region_module_matrix,
    BatchAlignmentDiagnostics = batch_alignment_diagnostics,
    ModuleCorrelationStructure = module_correlation_structure,
    ModuleMemberProteinEffects = module_member_effects,
    ModuleDriverConsistency = module_driver_consistency,
    ModuleCoverage = module_coverage,
    QCReview = qc_review_priority,
    SampleQC = sample_qc_summary,
    ReplicateQC = replicate_qc,
    ReplicateStatusCounts = replicate_status_counts,
    AnimalQC = animal_mean_qc
  ),
  file.path(dir_tables, "general_neuropil_module_score_statistics_full.xlsx"),
  overwrite = TRUE
)

# ------------------------------------------------
# 14) GROUP PLOTS: PRIMARY ALL AND QC SENSITIVITY
# ------------------------------------------------

format_p <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "0.000",
    TRUE ~ sprintf("%.3f", p)
  )
}

nature_module_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 15, face = "plain", color = "#4d4d4d", hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.title.y = element_text(size = 12.5, color = "#4d4d4d", margin = margin(r = 7)),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 13, color = "#4d4d4d", margin = margin(t = 2)),
      axis.text.y = element_text(size = 10.5, color = "#4d4d4d"),
      axis.line = element_line(linewidth = 0.75, color = "#4d4d4d"),
      axis.ticks = element_line(linewidth = 0.75, color = "#4d4d4d"),
      axis.ticks.length = unit(2.0, "mm"),
      plot.margin = margin(8, 8, 6, 8)
    )
}

mean_sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(data.frame(y = NA_real_, ymin = NA_real_, ymax = NA_real_))
  }

  sem <- if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
  data.frame(
    y = mean(x),
    ymin = mean(x) - sem,
    ymax = mean(x) + sem
  )
}

make_stat_annotation <- function(
    stats_df,
    plot_df,
    region_layer,
    module_name,
    analysis_label,
    p_column = "p_adj_within_model_BH",
    alpha = 0.05) {
  if (!p_column %in% colnames(stats_df)) {
    stop("Requested plot annotation p-value column is not present in Parametric_Contrasts.")
  }

  ann <- stats_df %>%
    filter(
      Analysis == analysis_label,
      RegionLayer == region_layer,
      Module == module_name,
      is.finite(.data[[p_column]]),
      .data[[p_column]] <= alpha
    ) %>%
    mutate(
      p_for_plot = .data[[p_column]],
      p.label = format_p(p_for_plot)
    ) %>%
    select(group1, group2, p.value, p_adj_within_model_BH, p_adj_global_BH, p_for_plot, p.label)

  if (nrow(ann) == 0) return(NULL)

  y_max <- max(plot_df$ModuleScore, na.rm = TRUE)
  y_min <- min(plot_df$ModuleScore, na.rm = TRUE)
  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range == 0) y_range <- 1

  ann %>% mutate(y.position = y_max + seq_len(n()) * 0.14 * y_range)
}

plot_module_scores <- function(
    df_animal,
    df_sample,
    contrast_stats_df,
    region_layer,
    module_name,
    analysis_label,
    contrast_p_column = "p_adj_within_model_BH") {
  plot_animal <- df_animal %>%
    filter(
      RegionLayer == region_layer,
      Module == module_name,
      is.finite(ModuleScore)
    ) %>%
    mutate(StressGroup = factor(StressGroup, levels = c("CON", "RES", "SUS")))

  plot_sample <- df_sample %>%
    filter(
      RegionLayer == region_layer,
      Module == module_name,
      is.finite(ModuleScore)
    ) %>%
    mutate(StressGroup = factor(StressGroup, levels = c("CON", "RES", "SUS")))

  if (nrow(plot_animal) == 0) return(NULL)

  ann <- make_stat_annotation(
    contrast_stats_df,
    plot_animal,
    region_layer,
    module_name,
    analysis_label,
    p_column = contrast_p_column
  )

  y_max <- max(c(plot_animal$ModuleScore, plot_sample$ModuleScore), na.rm = TRUE)
  y_min <- min(c(plot_animal$ModuleScore, plot_sample$ModuleScore), na.rm = TRUE)
  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  plot_title <- str_wrap(str_replace_all(module_name, "_", " "), width = 26)

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.75, color = "#4d4d4d") +
    geom_point(
      data = plot_sample,
      aes(
        x = StressGroup,
        y = ModuleScore,
        color = StressGroup
      ),
      position = position_jitter(width = 0.10, height = 0),
      shape = 16,
      size = 3.1,
      alpha = 0.48
    ) +
    stat_summary(
      data = plot_animal,
      aes(x = StressGroup, y = ModuleScore),
      fun.data = mean_sem,
      geom = "errorbar",
      width = 0.18,
      linewidth = 1.05,
      color = "#4d4d4d"
    ) +
    stat_summary(
      data = plot_animal,
      aes(x = StressGroup, y = ModuleScore),
      fun = mean,
      geom = "crossbar",
      width = 0.46,
      linewidth = 0.95,
      color = "#4d4d4d"
    ) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
    labs(
      title = plot_title,
      x = NULL,
      y = "Module score (z-score)"
    ) +
    coord_cartesian(
      ylim = c(
        y_min - 0.10 * y_range,
        y_max + ifelse(is.null(ann), 0.16, 0.28 + 0.14 * nrow(ann)) * y_range
      ),
      clip = "off"
    ) +
    nature_module_theme(base_size = 10)

  if (!is.null(ann)) {
    p <- p +
      stat_pvalue_manual(
        ann,
        label = "p.label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.85,
        size = 5.0,
        color = "#4d4d4d",
        hide.ns = TRUE
      )
  }

  p
}

region_layers <- unique(scores_animal_all$RegionLayer)
module_names <- names(modules)

for (rl in region_layers) {
  for (mod in module_names) {
    p_all <- plot_module_scores(
      df_animal = scores_animal_all,
      df_sample = sample_qc_summary,
      contrast_stats_df = contrasts_out,
      region_layer = rl,
      module_name = mod,
      analysis_label = analysis_primary
    )

    if (!is.null(p_all)) {
      ggsave(file.path(dir_group_primary, paste0("module_score_", rl, "_", mod, "_", analysis_primary, ".svg")),
             plot = p_all, width = 74, height = 86, units = "mm")
    }

    p_qc <- plot_module_scores(
      df_animal = scores_animal_qc_sensitivity,
      df_sample = scores_df_qc_sensitivity,
      contrast_stats_df = contrasts_out,
      region_layer = rl,
      module_name = mod,
      analysis_label = analysis_qc_sensitivity
    )

    if (!is.null(p_qc)) {
      ggsave(file.path(dir_group_qc, paste0("module_score_", rl, "_", mod, "_", analysis_qc_sensitivity, ".svg")),
             plot = p_qc, width = 74, height = 86, units = "mm")
    }
  }
}

# ------------------------------------------------
# 14.1 DIRECTIONAL ROBUSTNESS PLOTS - refined manuscript style
# ------------------------------------------------

# -----------------------------
# manuscript-style heatmap setup
# -----------------------------

region_order <- c(
  "CA1_slm", "CA1_so", "CA1_sr",
  "CA2_slm", "CA2_so", "CA2_sr",
  "CA3_sr",
  "DG_mo", "DG_po"
)

module_order <- c(
  "Neuropil_RNP_RNA_processing_main",
  "Neuropil_RNP_RNA_processing_full",
  "Neuropil_chromatin_RNP_related_exploratory",
  "Neuropil_mito_bioenergetics",
  "Neuropil_endo_lysosomal_proteostasis"
)

clean_module_label <- function(x) {
  x %>%
    str_replace("^Neuropil_", "") %>%
    str_replace_all("_", " ") %>%
    str_replace("RNP RNA processing main", "RNP/RNA\nmain") %>%
    str_replace("RNP RNA processing full", "RNP/RNA\nfull") %>%
    str_replace("chromatin RNP related exploratory", "Chromatin/RNP\nexploratory") %>%
    str_replace("mito bioenergetics", "Mitochondrial\nbioenergetics") %>%
    str_replace("endo lysosomal proteostasis", "Endo-lysosomal\nproteostasis")
}

ordered_levels <- function(x, preferred) {
  observed <- unique(as.character(x[!is.na(x)]))
  c(preferred[preferred %in% observed], sort(setdiff(observed, preferred)))
}

contrast_order <- c(
  "RES vs CON",
  "SUS vs CON",
  "SUS vs RES"
)

heat_cols <- c(
  low  = "#3D5A80",   # muted blue
  mid  = "#F7F4EF",   # warm off-white
  high = "#B23A34",   # muted red
  ink  = "#3F3F3F",
  grid = "#FFFFFF"
)

theme_nature_heat <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),

      axis.text.x = element_text(
        size = 6.5,
        color = heat_cols["ink"],
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(
        size = 6.8,
        color = heat_cols["ink"],
        lineheight = 0.9
      ),

      strip.text = element_text(
        size = 7.2,
        color = heat_cols["ink"],
        margin = margin(b = 3)
      ),
      strip.background = element_blank(),

      plot.title = element_text(
        size = 9.2,
        face = "plain",
        color = heat_cols["ink"],
        hjust = 0,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        size = 6.8,
        color = "#6A6A6A",
        hjust = 0,
        margin = margin(b = 4)
      ),

      legend.position = "right",
      legend.title = element_text(size = 6.5, color = heat_cols["ink"]),
      legend.text = element_text(size = 6, color = heat_cols["ink"]),
      legend.key.height = unit(13, "mm"),
      legend.key.width = unit(2.4, "mm"),

      plot.margin = margin(4, 5, 4, 4)
    )
}

scale_effect_fill <- function(name, limits) {
  scale_fill_gradient2(
    low = heat_cols["low"],
    mid = heat_cols["mid"],
    high = heat_cols["high"],
    midpoint = 0,
    limits = limits,
    oob = scales::squish,
    name = name,
    guide = guide_colorbar(
      frame.colour = heat_cols["ink"],
      frame.linewidth = 0.25,
      ticks.colour = heat_cols["ink"],
      ticks.linewidth = 0.25,
      barheight = unit(20, "mm"),
      barwidth = unit(2.5, "mm")
    )
  )
}

# -----------------------------
# 1) Main directional effect heatmap
# -----------------------------

plot_directional_effect_heatmap <- function(df, analysis_label) {

  plot_df <- df %>%
    filter(
      Analysis == analysis_label,
      is.finite(Cohen_d)
    )

  region_levels <- ordered_levels(plot_df$RegionLayer, region_order)
  module_levels <- ordered_levels(plot_df$Module, module_order)
  contrast_levels <- ordered_levels(plot_df$contrast_label, contrast_order)

  plot_df <- plot_df %>%
    mutate(
      RegionLayer = factor(RegionLayer, levels = region_levels),
      Module = factor(Module, levels = rev(module_levels)),
      Module_label = factor(clean_module_label(as.character(Module)), levels = clean_module_label(rev(module_levels))),
      contrast_label = factor(contrast_label, levels = contrast_levels),
      show_label = abs(Cohen_d) >= 0.8,
      d_label = sprintf("%.1f", Cohen_d)
    ) %>%
    filter(!is.na(RegionLayer), !is.na(Module_label), !is.na(contrast_label))

  if (nrow(plot_df) == 0) return(NULL)

  max_abs_d <- quantile(abs(plot_df$Cohen_d), 0.95, na.rm = TRUE)
  max_abs_d <- max(1, max_abs_d)
  max_abs_d <- min(max_abs_d, 2.5)

  ggplot(plot_df, aes(x = RegionLayer, y = Module_label, fill = Cohen_d)) +
    geom_tile(color = heat_cols["grid"], linewidth = 0.55, width = 0.94, height = 0.94) +

    geom_text(
      data = plot_df %>% filter(show_label),
      aes(label = d_label),
      size = 1.85,
      color = heat_cols["ink"]
    ) +

    geom_point(
      data = plot_df %>% filter(within_BH_significant),
      shape = 21,
      size = 1.55,
      color = heat_cols["ink"],
      fill = "white",
      stroke = 0.28
    ) +

    facet_wrap(~ contrast_label, nrow = 1) +
    scale_effect_fill("Cohen's d", limits = c(-max_abs_d, max_abs_d)) +
    coord_equal() +
    labs(
      title = paste0("Neuropil module effect structure"),
      subtitle = paste0(analysis_display_label(analysis_label), " | dots indicate BH-adjusted p <= 0.05; numbers shown for |d| >= 0.8")
    ) +
    theme_nature_heat(base_size = 7)
}

# -----------------------------
# 2) Module correlation heatmap
# -----------------------------

plot_module_correlation_heatmap <- function(df, analysis_label) {

  plot_df <- df %>%
    filter(
      Analysis == analysis_label,
      is.finite(rho)
    )

  if (nrow(plot_df) == 0) return(NULL)

  module_levels <- ordered_levels(c(plot_df$module1, plot_df$module2), module_order)
  region_levels <- ordered_levels(plot_df$RegionLayer, region_order)
  module_label_levels <- clean_module_label(module_levels)

  plot_df <- plot_df %>%
    mutate(
      module1_clean = clean_module_label(module1),
      module2_clean = clean_module_label(module2)
    ) %>%
    transmute(
      RegionLayer,
      module1 = module1_clean,
      module2 = module2_clean,
      rho,
      p_adj_BH
    )

  mirror_df <- plot_df %>%
    transmute(
      RegionLayer,
      module1 = module2,
      module2 = module1,
      rho,
      p_adj_BH
    )

  diag_df <- expand.grid(
    RegionLayer = unique(plot_df$RegionLayer),
    module1 = module_label_levels,
    module2 = module_label_levels,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    filter(module1 == module2) %>%
    mutate(
      rho = 1,
      p_adj_BH = 0
    )

  full_df <- bind_rows(plot_df, mirror_df, diag_df) %>%
    mutate(
      RegionLayer = factor(RegionLayer, levels = region_levels),
      module1 = factor(module1, levels = module_label_levels),
      module2 = factor(module2, levels = rev(module_label_levels))
    ) %>%
    filter(!is.na(RegionLayer), !is.na(module1), !is.na(module2))

  ggplot(full_df, aes(x = module1, y = module2, fill = rho)) +
    geom_tile(color = heat_cols["grid"], linewidth = 0.55, width = 0.94, height = 0.94) +

    geom_point(
      data = full_df %>% filter(module1 != module2, p_adj_BH <= 0.05),
      shape = 21,
      size = 1.45,
      color = heat_cols["ink"],
      fill = "white",
      stroke = 0.28
    ) +

    facet_wrap(~ RegionLayer) +
    scale_effect_fill("Spearman\nrho", limits = c(-1, 1)) +
    coord_equal() +
    labs(
      title = "Module correlation structure",
      subtitle = paste0(analysis_display_label(analysis_label), " | dots indicate BH-adjusted p <= 0.05")
    ) +
    theme_nature_heat(base_size = 6.5) +
    theme(
      axis.text.x = element_text(size = 5.5, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 5.5, lineheight = 0.9)
    )
}

# -----------------------------
# 3) Protein-driver consistency heatmap
# -----------------------------

plot_module_driver_consistency <- function(df, analysis_label) {

  plot_df <- df %>%
    filter(
      Analysis == analysis_label,
      is.finite(fraction_proteins_same_direction_as_module)
    )

  region_levels <- ordered_levels(plot_df$RegionLayer, region_order)
  module_levels <- ordered_levels(plot_df$Module, module_order)
  contrast_levels <- ordered_levels(plot_df$contrast_label, contrast_order)

  plot_df <- plot_df %>%
    mutate(
      RegionLayer = factor(RegionLayer, levels = region_levels),
      Module = factor(Module, levels = rev(module_levels)),
      Module_label = factor(clean_module_label(as.character(Module)), levels = clean_module_label(rev(module_levels))),
      contrast_label = factor(contrast_label, levels = contrast_levels),
      show_label = fraction_proteins_same_direction_as_module >= 0.70 |
        fraction_proteins_same_direction_as_module <= 0.45,
      frac_label = sprintf("%.2f", fraction_proteins_same_direction_as_module)
    ) %>%
    filter(!is.na(RegionLayer), !is.na(Module_label), !is.na(contrast_label))

  if (nrow(plot_df) == 0) return(NULL)

  ggplot(plot_df, aes(x = RegionLayer, y = Module_label, fill = fraction_proteins_same_direction_as_module)) +
    geom_tile(color = heat_cols["grid"], linewidth = 0.55, width = 0.94, height = 0.94) +

    geom_text(
      data = plot_df %>% filter(show_label),
      aes(label = frac_label),
      size = 1.85,
      color = heat_cols["ink"]
    ) +

    facet_wrap(~ contrast_label, nrow = 1) +

    scale_fill_gradient2(
      low = "#7A5A8A",
      mid = heat_cols["mid"],
      high = heat_cols["high"],
      midpoint = 0.50,
      limits = c(0, 1),
      oob = scales::squish,
      name = "Protein fraction\nsame direction",
      guide = guide_colorbar(
        frame.colour = heat_cols["ink"],
        frame.linewidth = 0.25,
        ticks.colour = heat_cols["ink"],
        ticks.linewidth = 0.25,
        barheight = unit(20, "mm"),
        barwidth = unit(2.5, "mm")
      )
    ) +

    coord_equal() +
    labs(
      title = "Module member driver consistency",
      subtitle = paste0(analysis_display_label(analysis_label), " | numbers shown for coherent or unstable modules")
    ) +
    theme_nature_heat(base_size = 7)
}

# -----------------------------
# 4) Directional consistency dot plot
# -----------------------------

plot_directional_consistency <- function(df, analysis_label) {

  plot_df <- df %>%
    filter(
      Analysis == analysis_label,
      is.finite(direction_consistency)
    )

  module_levels <- ordered_levels(plot_df$Module, module_order)
  contrast_levels <- ordered_levels(plot_df$contrast_label, contrast_order)

  plot_df <- plot_df %>%
    mutate(
      Module = factor(Module, levels = rev(module_levels)),
      Module_label = factor(clean_module_label(as.character(Module)), levels = clean_module_label(rev(module_levels))),
      contrast_label = factor(contrast_label, levels = contrast_levels),
      evidence_class = factor(
        evidence_class,
        levels = c(
          "convergent_with_variable_power",
          "directionally_stable_trend",
          "limited_evidence",
          "mixed_or_unstable_direction"
        )
      )
    ) %>%
    filter(!is.na(Module_label), !is.na(contrast_label))

  if (nrow(plot_df) == 0) return(NULL)

  evidence_cols <- c(
    convergent_with_variable_power = heat_cols["high"],
    directionally_stable_trend = heat_cols["low"],
    limited_evidence = "#8A8A84",
    mixed_or_unstable_direction = "#7A5A8A"
  )

  ggplot(
    plot_df,
    aes(
      x = direction_consistency,
      y = Module_label,
      color = evidence_class,
      size = mean_abs_Cohen_d
    )
  ) +
    geom_vline(xintercept = 0.80, linetype = "dashed", linewidth = 0.30, color = heat_cols["ink"]) +
    geom_point(alpha = 0.90, stroke = 0) +
    facet_wrap(~ contrast_label, nrow = 1) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_color_manual(values = evidence_cols, drop = FALSE) +
    scale_size_continuous(range = c(1.4, 3.4), name = "|mean d|") +
    labs(
      title = "Directional consistency across regions",
      subtitle = paste0(analysis_display_label(analysis_label), " | dashed line marks 80% directional agreement"),
      x = "Fraction in dominant direction",
      y = NULL,
      color = NULL
    ) +
    theme_minimal(base_size = 7) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.22, color = "#EAE6DD"),
      panel.spacing = unit(3.5, "mm"),

      plot.title = element_text(size = 9.2, color = heat_cols["ink"], hjust = 0, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 6.8, color = "#6A6A6A", hjust = 0, margin = margin(b = 4)),

      axis.title.x = element_text(size = 7, color = heat_cols["ink"], margin = margin(t = 4)),
      axis.text = element_text(size = 6.5, color = heat_cols["ink"]),

      strip.background = element_blank(),
      strip.text = element_text(size = 7.2, color = heat_cols["ink"], margin = margin(b = 3)),

      legend.position = "right",
      legend.title = element_text(size = 6.5, color = heat_cols["ink"]),
      legend.text = element_text(size = 6, color = heat_cols["ink"])
    )
}

# -----------------------------
# save plots
# -----------------------------

for (analysis_label in unique(directional_effects$Analysis)) {

  p_effect_heatmap <- plot_directional_effect_heatmap(directional_effects, analysis_label)

  if (!is.null(p_effect_heatmap)) {
    ggsave(
      file.path(dir_directional, paste0("NATURE_directional_effect_heatmap_", analysis_label, ".svg")),
      plot = p_effect_heatmap,
      width = 185,
      height = 95,
      units = "mm"
    )
    ggsave(
      file.path(dir_directional, paste0("NATURE_directional_effect_heatmap_", analysis_label, ".pdf")),
      plot = p_effect_heatmap,
      width = 185,
      height = 95,
      units = "mm"
    )
  }

  p_consistency <- plot_directional_consistency(directional_consistency_by_module, analysis_label)

  if (!is.null(p_consistency)) {
    ggsave(
      file.path(dir_directional, paste0("NATURE_directional_consistency_by_module_", analysis_label, ".svg")),
      plot = p_consistency,
      width = 175,
      height = 80,
      units = "mm"
    )
    ggsave(
      file.path(dir_directional, paste0("NATURE_directional_consistency_by_module_", analysis_label, ".pdf")),
      plot = p_consistency,
      width = 175,
      height = 80,
      units = "mm"
    )
  }

  p_module_cor <- plot_module_correlation_heatmap(module_correlation_structure, analysis_label)

  if (!is.null(p_module_cor)) {
    ggsave(
      file.path(dir_directional, paste0("NATURE_module_correlation_structure_", analysis_label, ".svg")),
      plot = p_module_cor,
      width = 195,
      height = 120,
      units = "mm"
    )
    ggsave(
      file.path(dir_directional, paste0("NATURE_module_correlation_structure_", analysis_label, ".pdf")),
      plot = p_module_cor,
      width = 195,
      height = 120,
      units = "mm"
    )
  }

  p_driver_consistency <- plot_module_driver_consistency(module_driver_consistency, analysis_label)

  if (!is.null(p_driver_consistency)) {
    ggsave(
      file.path(dir_directional, paste0("NATURE_module_driver_consistency_", analysis_label, ".svg")),
      plot = p_driver_consistency,
      width = 185,
      height = 95,
      units = "mm"
    )
    ggsave(
      file.path(dir_directional, paste0("NATURE_module_driver_consistency_", analysis_label, ".pdf")),
      plot = p_driver_consistency,
      width = 185,
      height = 95,
      units = "mm"
    )
  }
}

# ------------------------------------------------
# 15) CORRELATIONS
# ------------------------------------------------

correlation_vars <- c(
  "CombZ", "NOR_z", "SucrosePref_z", "WeightDev_z", "DeltaCort_z",
  "AdrenalWeight_z", "SpleenWeight_z",
  "AUC_all", "AUC_norm_all", "AUC_firstActive", "AUC_norm_firstActive"
)

correlation_vars <- intersect(correlation_vars, colnames(scores_animal_all))

run_correlations <- function(df, vars) {
  map_dfr(vars, function(v) {
    tmp <- df %>%
      filter(is.finite(ModuleScore), is.finite(.data[[v]]))

    if (nrow(tmp) < 5) {
      return(tibble(
        Variable = v,
        n = nrow(tmp),
        rho = NA_real_,
        p = NA_real_
      ))
    }

    ct <- cor.test(tmp$ModuleScore, tmp[[v]], method = "spearman", exact = FALSE)

    tibble(
      Variable = v,
      n = nrow(tmp),
      rho = unname(ct$estimate),
      p = ct$p.value
    )
  })
}

correlations_out <- scores_animal_all %>%
  group_by(RegionLayer, Module) %>%
  group_modify(~ run_correlations(.x, correlation_vars)) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p, method = "BH"),
    p.signif = signif_label(p),
    p.signif.adj = signif_label(p_adj_BH)
  )

write.xlsx(correlations_out, file.path(dir_tables, paste0("module_score_correlations_", analysis_primary, ".xlsx")), overwrite = TRUE)

# ------------------------------------------------
# 16) CORRELATION PLOTS
# ------------------------------------------------

plot_correlation <- function(df, cor_df, region_layer, module_name, variable) {
  plot_df <- df %>%
    filter(
      RegionLayer == region_layer,
      Module == module_name,
      is.finite(ModuleScore),
      is.finite(.data[[variable]])
    )

  if (nrow(plot_df) < 5) return(NULL)

  cor_label <- cor_df %>%
    filter(
      RegionLayer == region_layer,
      Module == module_name,
      Variable == variable
    ) %>%
    mutate(
      label = paste0(
        "Spearman rho = ", round(rho, 2),
        "\nBH p = ", signif(p_adj_BH, 2)
      )
    ) %>%
    pull(label)

  if (length(cor_label) == 0) cor_label <- ""

  ggplot(plot_df, aes(x = .data[[variable]], y = ModuleScore, color = StressGroup)) +
    geom_point(size = 2.2, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.4) +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = cor_label,
      hjust = -0.05,
      vjust = 1.15,
      size = 2.6
    ) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    labs(
      title = paste(region_layer, "-", module_name),
      x = variable,
      y = "Neuropil module score"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35)
    )
}

for (rl in region_layers) {
  for (mod in module_names) {
    for (v in correlation_vars) {
      p <- plot_correlation(scores_animal_all, correlations_out, rl, mod, v)

      if (!is.null(p)) {
        ggsave(
          file.path(dir_cor, paste0("cor_", rl, "_", mod, "_", v, "_", analysis_primary, ".svg")),
          plot = p,
          width = 68,
          height = 56,
          units = "mm"
        )
      }
    }
  }
}

cat("Done. Results saved to:\n", saving_dir, "\n")
