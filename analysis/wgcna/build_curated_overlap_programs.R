# ================================================================
# Script: analysis/wgcna/build_curated_overlap_programs.R
# Stage: modules_downstream
# Scope: global
# Consumes: required results/wgcna/compare_recurrent_module_proteins/neuron_neuropil/tables/recurrent_top_gene_drivers/overlapping_proteins_min3_spatial_units.xlsx; data/external/MOUSE_10090_idmapping.dat; optional none declared in pipeline.yml
# Produces: results/wgcna/build_curated_overlap_programs/global/tables
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
# Notes: Curated overlap programs are global and do not recompute WGCNA.
# ================================================================
# Build curated overlap-derived neuropil biological programs from the genes
# that drive the top GO terms in at least three neuron_neuropil spatial units
# (sheet Overlap_summary of compare_recurrent_module_proteins.R: Gene,
# N_datasets, Datasets, ProteinGroupIDs, member_accessions).
#
# Gene is the official gene symbol and the program rules classify it. One row
# is kept per gene and protein group. UniProt is the protein group's
# representative accession, its first member accession in sorted order (the
# fallback rule of select_representative_member() in
# R/data_contracts/protein_mapping_utils.R); member_accessions keeps all of the
# group's accessions, and the mitochondrially encoded accession list matches
# any of them. The UniProt mapping file supplies protein names and GO text by
# representative accession.
# ================================================================

library(readr)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(openxlsx)

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "module_contracts.R"))
source(repo_path("R", "wgcna_paths.R"))
source(repo_path("R", "top_term_driver_utils.R"))
MODULE_ID <- "06_modules_WGCNA"
SUBSTEP_ID <- file.path("curated_overlap_programs", "global")
CANONICAL_PATHS <- wgcna_dirs("build_curated_overlap_programs", "global")

# ------------------------------------------------
# 1) PATHS
# ------------------------------------------------

OVERLAP_DATASET <- "neuron_neuropil"
overlap_file <- file.path(
  wgcna_dirs("compare_recurrent_module_proteins", OVERLAP_DATASET)$tables,
  "recurrent_top_gene_drivers", "overlapping_proteins_min3_spatial_units.xlsx"
)
overlap_sheet <- "Overlap_summary"

mapping_file <- path_external("MOUSE_10090_idmapping.dat")

saving_dir <- CANONICAL_PATHS$tables

if (is_dry_run()) {
  dry_run_line("Script", "analysis/wgcna/build_curated_overlap_programs.R")
  dry_run_line("Dataset", "global")
  dry_run_line("Module source", "curated_overlap_programs")
  dry_run_line("Overlap file", overlap_file, if (file.exists(overlap_file)) "PASS" else "FAIL")
  dry_run_line("Mapping file", mapping_file, if (file.exists(mapping_file)) "PASS" else "FAIL")
  dry_run_line("Output folders", paste(unlist(CANONICAL_PATHS), collapse = "; "))
  dry_run_line("Expected output workbook", file.path(saving_dir, "curated_overlap_programs.xlsx"))
  dry_run_line("Expected manifest", file.path(CANONICAL_PATHS$logs, "curated_overlap_programs_manifest.yml"))
  quit(status = if (file.exists(overlap_file) && file.exists(mapping_file)) 0 else 1, save = "no")
}
if (!file.exists(overlap_file)) stop("Overlap file not found: ", overlap_file, call. = FALSE)
if (!file.exists(mapping_file)) stop("UniProt mapping file not found: ", mapping_file, call. = FALSE)

dir.create(saving_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(CANONICAL_PATHS$models, recursive = TRUE, showWarnings = FALSE)
dir.create(CANONICAL_PATHS$logs, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# 2) LOAD OVERLAP TABLE
# ------------------------------------------------

overlap_raw <- readxl::read_excel(overlap_file, sheet = overlap_sheet)
overlap_required <- c("Gene", "N_datasets", "Datasets", "ProteinGroupIDs", "member_accessions")
overlap_missing <- setdiff(overlap_required, names(overlap_raw))
if (length(overlap_missing)) {
  stop("Overlap table missing columns: ", paste(overlap_missing, collapse = ", "), call. = FALSE)
}

overlap_genes <- overlap_raw %>%
  mutate(
    Gene = str_trim(as.character(Gene)),
    N_datasets = as.integer(N_datasets),
    Datasets = as.character(Datasets),
    ProteinGroupIDs = as.character(ProteinGroupIDs),
    member_accessions = as.character(member_accessions)
  ) %>%
  filter(!is.na(Gene), Gene != "", N_datasets >= 3)

no_group <- overlap_genes$Gene[is.na(overlap_genes$ProteinGroupIDs) | !nzchar(overlap_genes$ProteinGroupIDs)]
if (length(no_group)) {
  stop("Overlap genes without a protein group: ", paste(utils::head(no_group, 10), collapse = ", "), call. = FALSE)
}

# one row per gene and protein group
overlap_df <- overlap_genes %>%
  mutate(groups = purrr::map2(ProteinGroupIDs, member_accessions, top_term_protein_group_pairs)) %>%
  dplyr::select(Gene, N_datasets, Datasets, groups) %>%
  tidyr::unnest(groups) %>%
  mutate(UniProt = top_term_representative_accession(member_accessions)) %>%
  distinct(Gene, ProteinGroupID, .keep_all = TRUE)

if (any(is.na(overlap_df$UniProt))) {
  stop("Protein groups without a member accession: ",
       paste(utils::head(overlap_df$ProteinGroupID[is.na(overlap_df$UniProt)], 10), collapse = ", "),
       call. = FALSE)
}

# ------------------------------------------------
# 3) LOAD UNIPROT MAPPING
# ------------------------------------------------

mapping <- readr::read_tsv(
  mapping_file,
  col_names = c("UniProt", "Type", "Value"),
  show_col_types = FALSE
)

protein_names <- mapping %>%
  filter(Type %in% c("Protein names", "Protein_Name")) %>%
  group_by(UniProt) %>%
  summarise(ProteinName = dplyr::first(Value), .groups = "drop")

go_terms <- mapping %>%
  filter(str_detect(Type, "^GO")) %>%
  group_by(UniProt) %>%
  summarise(GO_text = paste(unique(Value), collapse = " | "), .groups = "drop")

annot_df <- overlap_df %>%
  left_join(protein_names, by = "UniProt") %>%
  left_join(go_terms, by = "UniProt") %>%
  mutate(
    GeneName = Gene,
    ProteinName = if_else(is.na(ProteinName), "", ProteinName),
    GO_text = if_else(is.na(GO_text), "", GO_text),
    gene_lower = str_to_lower(GeneName),
    annotation_text = str_to_lower(
      paste(GeneName, ProteinName, GO_text, sep = " | ")
    )
  )

# ------------------------------------------------
# 4) GENE-SYMBOL-BASED MODULE RULES
# ------------------------------------------------

is_rnp_rna <- function(x) {
  str_detect(
    x,
    paste(
      c(
        "^hnrnp", "^snrnp", "^srsf", "^sf3", "^u2af", "^prpf",
        "^ddx", "^dhx", "^rbm", "^nono", "^pcbp", "^fus",
        "^matrin", "^pabp", "^pabpc", "^pabpn",
        "^ncl$", "^ncl", "^nol", "^nop", "^fbl$",
        "^rps", "^rpl"
      ),
      collapse = "|"
    )
  )
}

is_ribosome_translation <- function(x) {
  str_detect(
    x,
    paste(
      c(
        "^rps", "^rpl",
        "^eif", "^eef", "^eif3", "^eif4",
        "^mt-", "^mrpl", "^mrps"
      ),
      collapse = "|"
    )
  )
}

mito_encoded_accessions <- c(
  "P03888", "P03893", "P03899", "P03930", "P05064",
  "P0DN34", "P12023", "P12787", "P17665", "P35487",
  "P48771", "P56391", "P62897", "P63330", "P97450",
  "P99028"
)

# member_accessions: the protein group's accessions, ";"-separated; a group
# matches the accession list if any of its members does
is_mito <- function(x, member_accessions) {
  str_detect(
    x,
    paste(
      c(
        "^mt-", "^nduf", "^cox", "^cyc", "^uqcr", "^sdh",
        "^atp5", "^atp6", "^atp8",
        "^mrpl", "^mrps",
        "^idha", "^idhb", "^mdh", "^ogdh", "^cs$",
        "^aco2", "^dlst", "^dlat", "^pdha", "^pdhb",
        "^slc25", "^tomm", "^timm", "^hspa9"
      ),
      collapse = "|"
    )
  ) |
    vapply(member_accessions, function(a) any(top_term_split_accessions(a) %in% mito_encoded_accessions),
           logical(1), USE.NAMES = FALSE)
}

is_proteostasis <- function(x) {
  str_detect(
    x,
    paste(
      c(
        "^psm", "^uba", "^ube", "^ubc", "^ubqln",
        "^hsp", "^hspa", "^hspb", "^hspd", "^hspel",
        "^dna", "^cltc", "^clta",
        "^lamp", "^ctsa", "^ctsb", "^ctsd", "^ctsl",
        "^vps", "^rab", "^atg", "^sqstm"
      ),
      collapse = "|"
    )
  )
}

is_synaptic_cytoskeleton <- function(x) {
  str_detect(
    x,
    paste(
      c(
        "^syn", "^snap", "^stx", "^sv2", "^syt", "^dlg",
        "^camk", "^grin", "^gria", "^shank", "^homer",
        "^map", "^tubb", "^tuba", "^act", "^actb", "^actg",
        "^dyn", "^kif", "^myo", "^sept"
      ),
      collapse = "|"
    )
  )
}

is_chromatin_exploratory <- function(x) {
  str_detect(
    x,
    paste(
      c(
        "^hist", "^h1f", "^h2a", "^h2b", "^h3", "^h4",
        "^hmgb", "^hp1", "^cbx", "^chd", "^smarca",
        "^hdac", "^kat", "^set"
      ),
      collapse = "|"
    )
  )
}

# ------------------------------------------------
# 5) CLASSIFY PROTEINS
# ------------------------------------------------

classified_df <- annot_df %>%
  mutate(
    is_RNP_RNA = is_rnp_rna(gene_lower),
    is_Ribosome_Translation = is_ribosome_translation(gene_lower),
    is_Mito_Bioenergetics = is_mito(gene_lower, member_accessions),
    is_Proteostasis = is_proteostasis(gene_lower),
    is_Synaptic_Cytoskeleton = is_synaptic_cytoskeleton(gene_lower),
    is_Chromatin_Exploratory = is_chromatin_exploratory(gene_lower)
  ) %>%
  dplyr::select(
    UniProt, GeneName, ProteinGroupID, member_accessions, ProteinName, GO_text,
    N_datasets, Datasets, gene_lower, annotation_text,
    starts_with("is_")
  )

# ------------------------------------------------
# 6) BUILD MODULES
# ------------------------------------------------

module_rules <- list(
  Neuropil_RNP_RNA_processing_main = classified_df$N_datasets >= 6 & classified_df$is_RNP_RNA,
  Neuropil_RNP_RNA_processing_full = classified_df$N_datasets >= 3 & classified_df$is_RNP_RNA,
  Neuropil_ribosome_translation = classified_df$N_datasets >= 3 & classified_df$is_Ribosome_Translation,
  Neuropil_mito_bioenergetics = classified_df$N_datasets >= 3 & classified_df$is_Mito_Bioenergetics,
  Neuropil_endo_lysosomal_proteostasis = classified_df$N_datasets >= 3 & classified_df$is_Proteostasis,
  Neuropil_synaptic_cytoskeleton = classified_df$N_datasets >= 3 & classified_df$is_Synaptic_Cytoskeleton,
  Neuropil_chromatin_RNP_related_exploratory = classified_df$N_datasets >= 3 & classified_df$is_Chromatin_Exploratory
)

module_defs <- lapply(module_rules, function(hit) unique(classified_df$UniProt[hit]))

# ------------------------------------------------
# 7) EXPORT LONG MODULE TABLE
# ------------------------------------------------

module_long <- imap_dfr(module_rules, function(hit, module_name) {
  classified_df[hit, , drop = FALSE] %>%
    transmute(
      ModuleSet = "curated_overlap_programs",
      ModuleID = module_name,
      Module = module_name,
      UniProt,
      Source = "neuron_neuropil_recurrent_top_term_drivers",
      GeneName, ProteinName, ProteinGroupID, member_accessions,
      N_datasets, Datasets,
      is_RNP_RNA,
      is_Ribosome_Translation,
      is_Mito_Bioenergetics,
      is_Proteostasis,
      is_Synaptic_Cytoskeleton,
      is_Chromatin_Exploratory
    )
}) %>%
  mutate(GeneSymbol = GeneName) %>%
  arrange(Module, desc(N_datasets), GeneName, UniProt)
validate_curated_overlap_programs(module_long, "curated overlap programs")

module_summary <- module_long %>%
  count(Module, name = "n_proteins") %>%
  arrange(desc(n_proteins))

unassigned_df <- classified_df %>%
  filter(
    !is_RNP_RNA,
    !is_Ribosome_Translation,
    !is_Mito_Bioenergetics,
    !is_Proteostasis,
    !is_Synaptic_Cytoskeleton,
    !is_Chromatin_Exploratory
  ) %>%
  dplyr::select(UniProt, GeneName, ProteinGroupID, member_accessions, ProteinName, N_datasets, Datasets) %>%
  arrange(desc(N_datasets), GeneName)

# ------------------------------------------------
# 8) WRITE EXCEL
# ------------------------------------------------

wb <- createWorkbook()

addWorksheet(wb, "Module_summary")
writeData(wb, "Module_summary", module_summary)

addWorksheet(wb, "Modules_long")
writeData(wb, "Modules_long", module_long)

addWorksheet(wb, "All_classified")
writeData(wb, "All_classified", classified_df)

addWorksheet(wb, "Unassigned")
writeData(wb, "Unassigned", unassigned_df)

for (mod in names(module_defs)) {
  sheet_name <- str_sub(mod, 1, 31)

  mod_df <- module_long %>%
    filter(Module == mod) %>%
    dplyr::select(UniProt, GeneName, ProteinGroupID, member_accessions, ProteinName, N_datasets, Datasets)

  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, mod_df)
}

saveWorkbook(
  wb,
  file.path(saving_dir, "curated_overlap_programs.xlsx"),
  overwrite = TRUE
)
saveWorkbook(
  wb,
  file.path(saving_dir, "Overlap_based_neuropil_modules_classified.xlsx"),
  overwrite = TRUE
)

saveRDS(
  module_defs,
  file.path(CANONICAL_PATHS$models, "curated_overlap_programs.rds")
)
saveRDS(
  module_defs,
  file.path(CANONICAL_PATHS$models, "Overlap_based_neuropil_modules_classified.rds")
)

write_csv(
  module_long,
  file.path(CANONICAL_PATHS$models, "curated_overlap_programs_long.csv")
)
write_csv(
  module_long,
  file.path(CANONICAL_PATHS$models, "Overlap_based_neuropil_modules_classified_long.csv")
)
write_csv(
  module_long,
  file.path(saving_dir, "curated_overlap_programs_long.csv")
)
write_csv(
  module_long,
  file.path(saving_dir, "Overlap_based_neuropil_modules_classified_long.csv")
)
write_run_manifest(
  file.path(CANONICAL_PATHS$logs, "curated_overlap_programs_manifest.yml"),
  inputs = list(
    overlap_file = overlap_file,
    mapping_file = mapping_file
  ),
  outputs = list(
    workbook = file.path(saving_dir, "curated_overlap_programs.xlsx"),
    long_csv = file.path(saving_dir, "curated_overlap_programs_long.csv"),
    processed_rds = file.path(CANONICAL_PATHS$models, "curated_overlap_programs.rds")
  ),
  parameters = list(
    dataset_scope = "global",
    overlap_dataset = OVERLAP_DATASET,
    overlap_sheet = overlap_sheet,
    ModuleSet = "curated_overlap_programs",
    min_datasets = 3,
    uniprot_column = "representative accession of each protein group (first member accession in sorted order); all member accessions in member_accessions"
  ),
  notes = "Recurrent top-term driver genes of neuron_neuropil (>= 3 spatial units), classified by gene symbol into global-scoped curated overlap programs; one row per gene and protein group."
)
write_session_info(file.path(CANONICAL_PATHS$logs, "sessionInfo.txt"))

# ------------------------------------------------
# 9) PRINT RESULTS
# ------------------------------------------------

cat("\nModule counts:\n")
print(module_summary)

cat("\nUnassigned proteins:\n")
print(nrow(unassigned_df))

cat("\nSaved to:\n")
cat(saving_dir, "\n")
