# ================================================================
# Script: 08_behavior_physio_coupling/01_correlate_proteomics_with_behavior.r
# Stage: coupling
# Scope: dataset_specific
# Consumes: required data/external/behavior/auc_individual_animals_firstChangeActive.csv; optional data/processed/morpheus/20260218_pgmatrix_imputed_neuron_soma_71samples_missing70pct_with_metadata.xlsx; data/external/MOUSE_10090_idmapping.dat.
# Produces: results/figures/08_behavior_physio_coupling/01_correlate_proteomics_with_behavior/; results/tables/08_behavior_physio_coupling/01_correlate_proteomics_with_behavior/.
# Dataset behavior: runs for neuron_soma according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Legacy direct proteomics-behavior coupling retained as optional soma-only analysis.
# ================================================================

# --- 1. Load Necessary Libraries ---
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(readr)
library(readxl)
library(svglite)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "protein_mapping_utils.R"))
MODULE_ID <- "08_behavior_physio_coupling"
SUBSTEP_ID <- "correlate_proteomics_with_behavior"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)
out_dir <- path_or_env("PROTEOMICS_BEHAVIOR_COR_OUTPUT_DIR", CANONICAL_PATHS$figures, kind = "dir")
ensure_dir(out_dir)

# --- 2. Load Your Data Files ---
proteomics_input_file <- Sys.getenv(
  "PROTEOMICS_BEHAVIOR_COR_PROTEOMICS_FILE",
  unset = path_processed("01_preprocessing/excel_convert", "20260526_pgmatrix_imputed_neuron_soma_71samples_missing70pct_with_metadata.xlsx")
)
behavior_input_file <- Sys.getenv(
  "PROTEOMICS_BEHAVIOR_COR_BEHAVIOR_FILE",
  unset = path_external("behavior", "auc_individual_animals_firstChangeActive.csv")
)
if (!file.exists(proteomics_input_file)) stop("Proteomics behavior-correlation input not found: ", proteomics_input_file, call. = FALSE)
if (!file.exists(behavior_input_file)) stop("Behavior correlation input not found: ", behavior_input_file, call. = FALSE)
prot_data <- readxl::read_xlsx(proteomics_input_file, sheet = 1, col_names = TRUE, col_types = NULL, na = c("", "NA", "N/A"), trim_ws = TRUE, skip = 0) # Individual log2 intensities
beh_data  <- read.csv(behavior_input_file) # Column: 'MouseID', 'Group', 'Simulated_AUC'

# The 13 Leading Edge IDs you identified for "Skin Development"
leading_edge_ids <- c("E9Q557", "Q2VIS4", "Q9CQH5", "P97350", "Q02257", 
                      "P98200", "Q08189", "Q9DC29", "Q9JLF6", "P13516", 
                      "Q62470", "P06797", "Q9D8B6")

# --- 2c. Analysis Options ---
# Set to NULL to keep all. Example for your case: target_layer <- "sr"
target_region <- "ca1"
target_layer  <- "sp"
target_metric <- "Movement"  # Set to NULL to keep all metrics

# ExpGroup coding in proteomics metadata
expgroup_map <- c("1" = "CON", "2" = "RES", "3" = "SUS")

# --- 2b. UniProt Mapping Utilities ---

normalize_mouse_id <- function(x) {
  x <- as.character(x)
  x <- gsub("[^0-9]", "", x)
  x <- sub("^0+", "", x)
  x[x == ""] <- NA_character_
  x
}

normalize_meta_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "", x)
  x
}

get_meta_row <- function(tbl, key, sample_cols) {
  keys <- normalize_meta_key(tbl$raw_id)
  target <- normalize_meta_key(key)

  # Exact normalized match first
  idx <- which(keys == target)

  # Fallback: prefix/contains matching for odd formatting
  if (!length(idx)) idx <- which(startsWith(keys, target))
  if (!length(idx)) idx <- which(grepl(target, keys, fixed = TRUE))

  if (!length(idx)) return(setNames(rep(NA_character_, length(sample_cols)), sample_cols))
  vals <- as.character(tbl[idx[1], sample_cols, drop = TRUE])
  stats::setNames(vals, sample_cols)
}

is_trueish <- function(x) {
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y")
}

normalize_meta_value <- function(x) {
  nms <- names(x)
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("", "na", "n/a", "null")] <- NA_character_
  if (!is.null(nms)) names(x) <- nms
  x
}

diagnostic_id_table <- function(source, ids) {
  tibble::tibble(
    source = rep(source, length(ids)),
    MouseID = as.character(ids)
  )
}


behavior_idmapping_candidates <- c(
  "MOUSE_10090_idmapping.dat",
  file.path("Datasets", "MOUSE_10090_idmapping.dat"),
  path_external("MOUSE_10090_idmapping.dat")
)
behavior_idmapping_path <- behavior_idmapping_candidates[file.exists(behavior_idmapping_candidates)][1]
if (is.na(behavior_idmapping_path)) behavior_idmapping_path <- "MOUSE_10090_idmapping.dat"
uniprot_mapping <- load_mouse_idmapping(behavior_idmapping_path, auto_download = TRUE)
maps <- build_mouse_maps(uniprot_mapping)

leading_edge_accessions <- vapply(
  leading_edge_ids,
  map_token_to_mouse_accession,
  FUN.VALUE = character(1),
  entry_map = maps$entry_map,
  gene_map = maps$gene_map
) %>% unique()

leading_edge_accessions <- leading_edge_accessions[!is.na(leading_edge_accessions) & nzchar(leading_edge_accessions)]

# --- Diagnostic: leading edge ID mapping ---
le_map_tbl <- tibble::tibble(
  input_id  = leading_edge_ids,
  accession = vapply(leading_edge_ids, map_token_to_mouse_accession,
                     FUN.VALUE = character(1),
                     entry_map = maps$entry_map, gene_map = maps$gene_map)
)
cat("=== Leading Edge ID Mapping ===\n")
cat(sprintf("Mapped: %d / %d\n", sum(!is.na(le_map_tbl$accession)), nrow(le_map_tbl)))
print(le_map_tbl)
cat("================================\n")

# --- 3. Extract and Calculate Individual Scores ---

# Step A: Map proteomics row IDs to canonical UniProt accessions and keep the Leading Edge proteins
prot_tbl <- prot_data
if (ncol(prot_tbl) < 1) stop("prot_data has no columns.")
names(prot_tbl)[1] <- "raw_id"

prot_tbl <- prot_tbl %>%
  mutate(raw_id = as.character(raw_id)) %>%
  mutate(
    mapped_accession = vapply(
      raw_id,
      map_token_to_mouse_accession,
      FUN.VALUE = character(1),
      entry_map = maps$entry_map,
      gene_map = maps$gene_map
    )
  )

sample_cols <- setdiff(names(prot_tbl), c("raw_id", "mapped_accession"))
sample_cols <- sample_cols[!is.na(sample_cols) & nzchar(sample_cols)]
sample_cols <- unique(sample_cols)

# Metadata rows that describe each sample column
meta_region <- get_meta_row(prot_tbl, "region", sample_cols)
meta_layer <- get_meta_row(prot_tbl, "layer", sample_cols)
meta_group <- get_meta_row(prot_tbl, "group", sample_cols)
meta_exclude <- get_meta_row(prot_tbl, "exclude", sample_cols)
meta_animal <- get_meta_row(prot_tbl, "AnimalID", sample_cols)
meta_expgroup <- get_meta_row(prot_tbl, "ExpGroup", sample_cols)

meta_region <- normalize_meta_value(meta_region)
meta_layer <- normalize_meta_value(meta_layer)
meta_group <- normalize_meta_value(meta_group)

if (all(is.na(meta_layer))) {
  stop("Layer filtering is configured to use the 'layer' metadata row, but that row is missing or empty.")
}

# Light harmonization for common layer aliases
meta_layer[meta_layer == "stratum radiatum"] <- "sr"
meta_layer[meta_layer == "stratum pyramidale"] <- "sp"

selected_cols <- sample_cols
if (!is.null(target_region)) {
  keep_region <- !is.na(meta_region[selected_cols]) & meta_region[selected_cols] == tolower(target_region)
  selected_cols <- selected_cols[keep_region]
}
if (!is.null(target_layer)) {
  keep_layer <- !is.na(meta_layer[selected_cols]) & meta_layer[selected_cols] == tolower(target_layer)
  selected_cols <- selected_cols[keep_layer]
}

# Respect exclude flag if present
if (any(!is.na(meta_exclude[selected_cols]))) {
  keep_exclude <- !is_trueish(meta_exclude[selected_cols])
  keep_exclude[is.na(keep_exclude)] <- TRUE
  selected_cols <- selected_cols[keep_exclude]
}

selected_cols <- selected_cols[!is.na(selected_cols) & nzchar(selected_cols)]
selected_cols <- intersect(selected_cols, names(prot_tbl))

if (!length(selected_cols)) {
  available_layers <- sort(unique(stats::na.omit(meta_layer[sample_cols])))
  available_regions <- sort(unique(stats::na.omit(meta_region[sample_cols])))
  stop(
    paste0(
      "No proteomics sample columns left after region/layer/exclude filtering. ",
      "Requested region=", ifelse(is.null(target_region), "<all>", target_region),
      ", layer=", ifelse(is.null(target_layer), "<all>", target_layer), ". ",
      "Available layers: ", ifelse(length(available_layers), paste(available_layers, collapse = ", "), "<none>"), ". ",
      "Available regions: ", ifelse(length(available_regions), paste(available_regions, collapse = ", "), "<none>")
    )
  )
}

subset_prot <- prot_tbl %>%
  filter(!is.na(mapped_accession), mapped_accession %in% leading_edge_accessions) %>%
  dplyr::select(all_of(c("mapped_accession", selected_cols))) %>%
  mutate(across(-mapped_accession, ~ suppressWarnings(as.numeric(.)))) %>%
  dplyr::select(where(~ !is.character(.) || any(!is.na(.))))

# --- Diagnostic: leading edge proteins found in data ---
found_accessions <- unique(subset_prot$mapped_accession)
missing_accessions <- setdiff(leading_edge_accessions, found_accessions)
cat("=== Leading Edge Proteins in Proteomics Data ===\n")
cat(sprintf("Found: %d / %d accessions\n", length(found_accessions), length(leading_edge_accessions)))
if (length(found_accessions))  cat("Found:   ", paste(found_accessions,   collapse = ", "), "\n")
if (length(missing_accessions)) cat("Missing: ", paste(missing_accessions, collapse = ", "), "\n")
cat("================================================\n")

subset_prot <- subset_prot %>% dplyr::select(-mapped_accession)

if (nrow(subset_prot) == 0) {
  stop("No overlap found between leading_edge_ids and mapped proteomics IDs. Check mapping input/data source.")
}

if (ncol(subset_prot) == 0) {
  stop("No numeric proteomics sample columns remain after filtering.")
}

# Step B: Calculate the average log2 intensity per mouse (column-wise mean)
mouse_scores <- tibble::tibble(
  sample_col = colnames(subset_prot),
  Structural_Load_Mean = as.numeric(colMeans(subset_prot, na.rm = TRUE)),
  MouseID = normalize_mouse_id(meta_animal[colnames(subset_prot)]),
  ExpGroup = as.character(meta_expgroup[colnames(subset_prot)])
) %>%
  mutate(Group_from_prot = unname(expgroup_map[ExpGroup])) %>%
  filter(!is.na(MouseID) & nzchar(MouseID)) %>%
  group_by(MouseID) %>%
  summarise(
    Structural_Load_Mean = mean(Structural_Load_Mean, na.rm = TRUE),
    Group_from_prot = dplyr::first(na.omit(Group_from_prot)),
    .groups = "drop"
  )

if (nrow(mouse_scores) == 0) {
  stop("No valid MouseID values found in proteomics metadata row 'AnimalID'.")
}

# Standardize behavior columns for joining
beh_tbl <- beh_data

if (!is.null(target_metric)) {
  if (!"metric" %in% names(beh_tbl)) stop("Behavior table does not contain a 'metric' column.")
  beh_tbl <- beh_tbl %>% filter(metric == target_metric)
  if (nrow(beh_tbl) == 0) stop(paste0("No rows in behavior table match metric='", target_metric, "'."))
}

if (!"AUC" %in% names(beh_tbl) && "Simulated_AUC" %in% names(beh_tbl)) {
  beh_tbl <- beh_tbl %>% rename(AUC = Simulated_AUC)
}
if (!"Group" %in% names(beh_tbl)) {
  stop("Behavior table must contain a Group column.")
}

beh_id_col <- dplyr::case_when(
  "AnimalNum" %in% names(beh_tbl) ~ "AnimalNum",
  "MouseID" %in% names(beh_tbl) ~ "MouseID",
  TRUE ~ NA_character_
)
if (is.na(beh_id_col)) stop("Behavior table must contain AnimalNum or MouseID for joining.")

beh_tbl <- beh_tbl %>%
  mutate(MouseID = normalize_mouse_id(.data[[beh_id_col]])) %>%
  filter(!is.na(MouseID) & nzchar(MouseID))

if (!"AUC" %in% names(beh_tbl)) {
  stop("Behavior table must contain AUC (or Simulated_AUC).")
}

# Step C: Merge with Behavioral Data

# --- Diagnostic: show ID overlap before joining ---
cat("=== ID Matching Diagnostics ===\n")
cat("Proteomics MouseIDs (from AnimalID metadata row):\n")
print(sort(unique(mouse_scores$MouseID)))
cat("Behavior MouseIDs (normalized):\n")
print(sort(unique(beh_tbl$MouseID)))
matched <- intersect(mouse_scores$MouseID, beh_tbl$MouseID)
prot_only <- setdiff(mouse_scores$MouseID, beh_tbl$MouseID)
beh_only  <- setdiff(beh_tbl$MouseID, mouse_scores$MouseID)
cat("Matched (", length(matched), "):", paste(matched, collapse = ", "), "\n")
cat("Prot only (", length(prot_only), "):", paste(prot_only, collapse = ", "), "\n")
cat("Beh only (", length(beh_only), "):", paste(beh_only, collapse = ", "), "\n")
cat("================================\n")

final_df <- inner_join(beh_tbl, mouse_scores, by = "MouseID") %>%
  mutate(Group = dplyr::coalesce(Group, Group_from_prot))

if (nrow(final_df) == 0) {
  stop("No overlap between behavior IDs and proteomics AnimalID after normalization.")
}

join_diagnostics <- dplyr::bind_rows(
  data.frame(metric = "proteomics_animals_before_join", value = length(unique(mouse_scores$MouseID))),
  data.frame(metric = "behavior_animals_before_join", value = length(unique(beh_tbl$MouseID))),
  data.frame(metric = "animals_after_join", value = length(unique(final_df$MouseID))),
  data.frame(metric = "animals_lost_from_proteomics", value = length(prot_only)),
  data.frame(metric = "animals_lost_from_behavior", value = length(beh_only))
)
write.csv(join_diagnostics, file.path(CANONICAL_PATHS$tables, "join_diagnostics_summary.csv"), row.names = FALSE)
write.csv(diagnostic_id_table("proteomics_only", prot_only), file.path(CANONICAL_PATHS$tables, "join_diagnostics_proteomics_only.csv"), row.names = FALSE)
write.csv(diagnostic_id_table("behavior_only", beh_only), file.path(CANONICAL_PATHS$tables, "join_diagnostics_behavior_only.csv"), row.names = FALSE)
coverage_cols <- intersect(c("Sex", "sex", "Group", "window"), names(final_df))
if (length(coverage_cols) > 0) {
  coverage <- final_df %>%
    dplyr::count(dplyr::across(dplyr::all_of(coverage_cols)), name = "n")
  write.csv(coverage, file.path(CANONICAL_PATHS$tables, "join_diagnostics_coverage.csv"), row.names = FALSE)
}

# Step D: Z-Score the protein data relative to the CONTROL group
# This makes the Y-axis represent "deviation from normal"
con_mean <- mean(final_df$Structural_Load_Mean[final_df$Group == "CON"], na.rm = TRUE)
con_sd   <- sd(final_df$Structural_Load_Mean[final_df$Group == "CON"], na.rm = TRUE)

final_df <- final_df %>%
  mutate(Structural_ZScore = (Structural_Load_Mean - con_mean) / con_sd)

# --- 4. Statistical Analysis (Correlation) ---

# Run a linear model to see if AUC predicts the protein signature
fit <- lm(Structural_ZScore ~ AUC, data = final_df)
summary_fit <- summary(fit)
r_squared <- summary_fit$r.squared
p_val <- summary_fit$coefficients[2,4]

# --- 5. Publication-style Visualization ---

publication_theme <- theme_classic(base_size = 7, base_family = "sans") +
  theme(
    # Axes
    axis.line        = element_line(colour = "black", linewidth = 0.4),
    axis.ticks       = element_line(colour = "black", linewidth = 0.4),
    axis.ticks.length = unit(2, "pt"),
    axis.text        = element_text(size = 6, colour = "black"),
    axis.title       = element_text(size = 7, colour = "black"),
    # Legend
    legend.position  = "top",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 6),
    legend.key.size  = unit(6, "pt"),
    legend.spacing.x = unit(3, "pt"),
    # Title / subtitle as annotation (not used as a panel title)
    plot.title       = element_text(size = 7, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 6, hjust = 0, colour = "grey40"),
    # Panel
    panel.grid       = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.margin      = margin(4, 6, 4, 4, "pt")
  )

p_label <- paste0("italic(R)^2 == ", round(r_squared, 3),
                  "~'|'~italic(p) == ", format.pval(p_val, digits = 2, eps = 0.001))

plot <- ggplot(final_df, aes(x = AUC, y = Structural_ZScore, color = Group)) +
  geom_smooth(method = "lm", colour = "black", linewidth = 0.5,
              linetype = "dashed", se = TRUE, fill = "grey85", alpha = 0.4) +
  geom_point(size = 1.6, alpha = 0.85, stroke = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = p_label, parse = TRUE,
           hjust = -0.05, vjust = 1.4, size = 2, colour = "grey30") +
  scale_color_manual(values = c("CON" = "#888888", "RES" = "#377eb8", "SUS" = "#e41a1c")) +
  labs(
    title    = "Hippocampal Structural Leading Edge vs. Locomotor Load",
    x        = "12 h locomotor load (AUC)",
    y        = "Structural leading-edge z-score\n(CA1 sp, relative to CON)"
  ) +
  publication_theme

# --- 6. Save Everything ---

# Save as SVG (Publication single-column width: 88 mm)
svglite::svglite(file.path(out_dir, "Figure3_AUC_vs_Proteomics.svg"),
                 width = 1.5, height = 2)
print(plot)
dev.off()

# Save the final data table for your Supplemental Information
write.csv(final_df, file.path(out_dir, "SourceData_Figure3_Correlation.csv"), row.names = FALSE)
write_run_manifest(
  file.path(CANONICAL_PATHS$logs, "run_manifest.yml"),
  inputs = list(proteomics = proteomics_input_file, behavior = behavior_input_file),
  outputs = list(figure_dir = out_dir, join_diagnostics = CANONICAL_PATHS$tables),
  parameters = list(target_region = target_region, target_layer = target_layer, target_metric = target_metric)
)

# Print success message
cat("Analysis complete. R²:", r_squared, "| P-value:", p_val)
