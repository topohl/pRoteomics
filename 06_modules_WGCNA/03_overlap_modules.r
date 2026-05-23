# ================================================================
# Consumes:
#   - compareGO overlap table from canonical results/tables
#   - UniProt mapping file from data/external
# Produces:
#   - overlap-derived module definitions under data/processed and results/tables
# File contract:
#   - docs/active_script_io_audit.tsv object 06_modules_WGCNA/03_overlap_modules.r
# ================================================================
# Build overlap-based neuropil modules from recurrent proteins
# Requires overlap table with columns:
# Gene, N_datasets, Datasets
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
MODULE_ID <- "06_modules_WGCNA"
SUBSTEP_ID <- "overlap_modules"
CANONICAL_PATHS <- create_module_dirs(MODULE_ID, SUBSTEP_ID)

# ------------------------------------------------
# 1) PATHS
# ------------------------------------------------

overlap_file <- path_results("tables", "04_differential_expression_enrichment", "compareGO", "overlapping_proteins_min3_datasets.xlsx")

mapping_file <- path_external("MOUSE_10090_idmapping.dat")

saving_dir <- CANONICAL_PATHS$tables

dir.create(saving_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(CANONICAL_PATHS$processed, recursive = TRUE, showWarnings = FALSE)

if (is_dry_run()) {
  dry_run_line("Script", "06_modules_WGCNA/03_overlap_modules.r")
  dry_run_line("Overlap file", overlap_file, if (file.exists(overlap_file)) "PASS" else "FAIL")
  dry_run_line("Mapping file", mapping_file, if (file.exists(mapping_file)) "PASS" else "FAIL")
  dry_run_line("Output folders", paste(unlist(CANONICAL_PATHS), collapse = "; "))
  quit(status = if (file.exists(overlap_file) && file.exists(mapping_file)) 0 else 1, save = "no")
}
if (!file.exists(overlap_file)) stop("Overlap file not found: ", overlap_file, call. = FALSE)
if (!file.exists(mapping_file)) stop("UniProt mapping file not found: ", mapping_file, call. = FALSE)

# ------------------------------------------------
# 2) LOAD OVERLAP TABLE
# ------------------------------------------------

overlap_df <- if (str_detect(overlap_file, "\\.xlsx$")) {
  readxl::read_excel(overlap_file)
} else {
  readr::read_csv(overlap_file, show_col_types = FALSE)
}

overlap_df <- overlap_df %>%
  { missing <- setdiff(c("Gene", "N_datasets", "Datasets"), names(.)); if (length(missing)) stop("Overlap table missing columns: ", paste(missing, collapse = ", "), call. = FALSE); . } %>%
  rename(UniProt = Gene) %>%
  mutate(
    UniProt = as.character(UniProt),
    N_datasets = as.integer(N_datasets)
  ) %>%
  filter(N_datasets >= 3)

# ------------------------------------------------
# 3) LOAD UNIPROT MAPPING
# ------------------------------------------------

mapping <- readr::read_tsv(
  mapping_file,
  col_names = c("UniProt", "Type", "Value"),
  show_col_types = FALSE
)

gene_names <- mapping %>%
  filter(Type %in% c("Gene_Name", "Gene Names", "Gene_Synonym")) %>%
  group_by(UniProt) %>%
  summarise(GeneName = first(Value), .groups = "drop")

protein_names <- mapping %>%
  filter(Type %in% c("Protein names", "Protein_Name")) %>%
  group_by(UniProt) %>%
  summarise(ProteinName = first(Value), .groups = "drop")

go_terms <- mapping %>%
  filter(str_detect(Type, "^GO")) %>%
  group_by(UniProt) %>%
  summarise(GO_text = paste(unique(Value), collapse = " | "), .groups = "drop")

annot_df <- overlap_df %>%
  left_join(gene_names, by = "UniProt") %>%
  left_join(protein_names, by = "UniProt") %>%
  left_join(go_terms, by = "UniProt") %>%
  mutate(
    GeneName = if_else(is.na(GeneName), "", GeneName),
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

is_mito <- function(x, uniprot) {
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
    uniprot %in% c(
      "P03888", "P03893", "P03899", "P03930", "P05064",
      "P0DN34", "P12023", "P12787", "P17665", "P35487",
      "P48771", "P56391", "P62897", "P63330", "P97450",
      "P99028"
    )
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
    is_Mito_Bioenergetics = is_mito(gene_lower, UniProt),
    is_Proteostasis = is_proteostasis(gene_lower),
    is_Synaptic_Cytoskeleton = is_synaptic_cytoskeleton(gene_lower),
    is_Chromatin_Exploratory = is_chromatin_exploratory(gene_lower)
  )

# ------------------------------------------------
# 6) BUILD MODULES
# ------------------------------------------------

module_defs <- list(
  Neuropil_RNP_RNA_processing_main =
    classified_df %>%
    filter(N_datasets >= 6, is_RNP_RNA) %>%
    pull(UniProt) %>%
    unique(),

  Neuropil_RNP_RNA_processing_full =
    classified_df %>%
    filter(N_datasets >= 3, is_RNP_RNA) %>%
    pull(UniProt) %>%
    unique(),

  Neuropil_ribosome_translation =
    classified_df %>%
    filter(N_datasets >= 3, is_Ribosome_Translation) %>%
    pull(UniProt) %>%
    unique(),

  Neuropil_mito_bioenergetics =
    classified_df %>%
    filter(N_datasets >= 3, is_Mito_Bioenergetics) %>%
    pull(UniProt) %>%
    unique(),

  Neuropil_endo_lysosomal_proteostasis =
    classified_df %>%
    filter(N_datasets >= 3, is_Proteostasis) %>%
    pull(UniProt) %>%
    unique(),

  Neuropil_synaptic_cytoskeleton =
    classified_df %>%
    filter(N_datasets >= 3, is_Synaptic_Cytoskeleton) %>%
    pull(UniProt) %>%
    unique(),

  Neuropil_chromatin_RNP_related_exploratory =
    classified_df %>%
    filter(N_datasets >= 3, is_Chromatin_Exploratory) %>%
    pull(UniProt) %>%
    unique()
)

# ------------------------------------------------
# 7) EXPORT LONG MODULE TABLE
# ------------------------------------------------

module_long <- imap_dfr(module_defs, function(accessions, module_name) {
  tibble(
    Module = module_name,
    UniProt = accessions
  )
}) %>%
  left_join(
    classified_df %>%
      select(
        UniProt, GeneName, ProteinName,
        N_datasets, Datasets,
        is_RNP_RNA,
        is_Ribosome_Translation,
        is_Mito_Bioenergetics,
        is_Proteostasis,
        is_Synaptic_Cytoskeleton,
        is_Chromatin_Exploratory
      ),
    by = "UniProt"
  ) %>%
  arrange(Module, desc(N_datasets), GeneName, UniProt)

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
  select(UniProt, GeneName, ProteinName, N_datasets, Datasets) %>%
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
    select(UniProt, GeneName, ProteinName, N_datasets, Datasets)

  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, mod_df)
}

saveWorkbook(
  wb,
  file.path(saving_dir, "Overlap_based_neuropil_modules_classified.xlsx"),
  overwrite = TRUE
)

saveRDS(
  module_defs,
  file.path(CANONICAL_PATHS$processed, "Overlap_based_neuropil_modules_classified.rds")
)

write_csv(
  module_long,
  file.path(CANONICAL_PATHS$processed, "Overlap_based_neuropil_modules_classified_long.csv")
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
