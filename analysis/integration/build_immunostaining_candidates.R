#!/usr/bin/env Rscript
# ================================================================
# Script: analysis/integration/build_immunostaining_candidates.R
# Stage: integration
# Scope: global (neuron_neuropil only)
# Consumes: required data/processed/01_preprocessing/protigy_input_animal_level/neuron_neuropil/neuron_neuropil_animal_level.gct;
#   data/processed/02_id_mapping/mapped/neuron_neuropil/forward/per_file/<unit>sus_<unit>res.csv;
#   config/manuscript_spatial_order.yml.
# Produces: results/source_data/10_biological_integration/immunostaining_candidate_comparison/.
# Notes: Source-data preparation only. No new statistics.
# ================================================================
#
# WHAT THIS SCRIPT DOES
#   Turns three already-nominated immunostaining candidates into frozen,
#   manuscript-facing source data. It reads the canonical bilateral animal-level
#   neuropil abundance matrix that the differential-abundance analysis was run
#   on, and the canonical per-spatial-unit SUS-RES protein results, and writes
#   tidy tables for exactly those three protein groups.
#
# WHAT THIS SCRIPT DOES NOT DO
#   No aggregation is redone: the GCT is already bilateral animal-level, one
#   value per protein group x animal x spatial unit. No model is refitted, no
#   contrast is recomputed, no candidate is selected or ranked, no threshold is
#   tuned. Every log2FC, p and FDR is copied from its canonical file and then
#   verified against it.
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

MODULE_ID <- "10_biological_integration"
SUBSTEP_ID <- "immunostaining_candidate_comparison"
DATASET <- "neuron_neuropil"
CONTRACT_VERSION <- "immunostaining_candidate_comparison_v1"

OUT <- path_results("source_data", MODULE_ID, SUBSTEP_ID)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# The three candidates, nominated upstream. Requested UniProt accessions are
# recorded so the resolved mapping can be checked against them rather than
# trusted. SLC22A23 was nominated by symbol only.
CANDIDATES <- data.frame(
  requested_symbol = c("OGA", "SLC22A23", "ANXA2"),
  requested_uniprot = c("Q9EQQ9", NA_character_, "P07356"),
  stringsAsFactors = FALSE)

# ---------------------------------------------------------------- inputs
GCT <- repo_path("data", "processed", "01_preprocessing",
                 "protigy_input_animal_level", DATASET,
                 paste0(DATASET, "_animal_level.gct"))
DA_DIR <- repo_path("data", "processed", "02_id_mapping", "mapped", DATASET,
                    "forward", "per_file")
if (!file.exists(GCT)) stop("missing_required_input: ", GCT, call. = FALSE)
if (!dir.exists(DA_DIR)) stop("missing_required_input: ", DA_DIR, call. = FALSE)

# Canonical anatomical order. Region-major, layers deep-to-superficial, taken
# from the shared contract rather than parsed out of the unit string.
ORDER_YML <- repo_path("config", "manuscript_spatial_order.yml")
units_yml <- yaml::read_yaml(ORDER_YML)$units
npx <- Filter(function(z) identical(z$dataset, DATASET), units_yml)
SPATIAL <- data.frame(
  spatial_unit = vapply(npx, function(z) z$analysis_key, character(1)),
  spatial_unit_display = vapply(npx, function(z) z$display, character(1)),
  region = vapply(npx, function(z) z$region, character(1)),
  layer = vapply(npx, function(z) z$layer, character(1)),
  spatial_order = vapply(npx, function(z) as.integer(z$order), integer(1)),
  stringsAsFactors = FALSE)
SPATIAL <- SPATIAL[order(SPATIAL$spatial_order), ]
if (nrow(SPATIAL) != 10L)
  stop("expected 10 neuropil spatial units, got ", nrow(SPATIAL), call. = FALSE)

# ExpGroup is stored numerically in the GCT. The mapping is fixed by the study
# design and is asserted below against the known animal rosters rather than
# assumed.
EXPGROUP_LABEL <- c(`1` = "CON", `2` = "RES", `3` = "SUS")
ROSTER <- list(CON = c("A127", "A129", "A765"),
               RES = c("A0003", "A135", "A139"),
               SUS = c("A111", "A755", "A764"))

# ------------------------------------------------- read the animal-level GCT
gct_lines <- readLines(GCT, n = 12L)
dims <- as.integer(strsplit(gct_lines[2], "\t")[[1]])
n_row_meta <- dims[3]; n_col_meta <- dims[4]
header <- strsplit(gct_lines[3], "\t")[[1]]
sample_ids <- header[-seq_len(2L + n_row_meta - 1L)]
meta_rows <- lapply(seq_len(n_col_meta), function(i)
  strsplit(gct_lines[3L + i], "\t")[[1]])
meta_names <- vapply(meta_rows, function(x) x[1], character(1))
meta_value <- function(key) {
  i <- which(meta_names == key)
  if (!length(i)) stop("GCT column metadata missing: ", key, call. = FALSE)
  meta_rows[[i]][-seq_len(2L + n_row_meta - 1L)]
}
sample_meta <- data.frame(
  sample_id = sample_ids,
  AnimalID = meta_value("AnimalID"),
  ExpGroupCode = meta_value("ExpGroup"),
  region = meta_value("region"),
  layer = meta_value("layer"),
  stringsAsFactors = FALSE)
sample_meta$ExpGroup <- unname(EXPGROUP_LABEL[sample_meta$ExpGroupCode])
sample_meta$spatial_unit <- paste(sample_meta$region, sample_meta$layer, sep = "_")

# assert the numeric ExpGroup coding against the known rosters
for (g in names(ROSTER)) {
  got <- sort(unique(sample_meta$AnimalID[sample_meta$ExpGroup == g]))
  if (!identical(got, sort(ROSTER[[g]])))
    stop("ExpGroup coding mismatch for ", g, ": got ",
         paste(got, collapse = ","), call. = FALSE)
}

# GCT layout: line 1 version, line 2 dims, line 3 header, then n_col_meta
# metadata rows, then the data block.
expr <- utils::read.delim(GCT, skip = 3L + n_col_meta, header = FALSE,
                          stringsAsFactors = FALSE, check.names = FALSE)
names(expr) <- header
# The GCT keys rows by the ORIGINAL identifier (e.g. OGA_MOUSE), not by
# ProteinGroupID, so the mapping table carries that key across.
row_id <- expr[[1]]

# ------------------------------------------------------- resolve the mapping
# Mapping is taken from a canonical DA file, which already carries the frozen
# gene-annotation contract. Ambiguity is a hard failure.
probe <- utils::read.csv(file.path(DA_DIR, "CA2slmsus_CA2slmres.csv"),
                         stringsAsFactors = FALSE)
resolve_one <- function(sym, acc) {
  # which() rather than logical subsetting: official_gene_symbol carries NA for
  # unannotated groups, and `NA == "OGA"` is NA, which would return one NA row
  # per unannotated protein group instead of no rows.
  hit <- probe[which(toupper(probe$official_gene_symbol) == toupper(sym)), ,
               drop = FALSE]
  if (!nrow(hit) && !is.na(acc))
    hit <- probe[which(grepl(acc, probe$member_accessions, fixed = TRUE)), ,
                 drop = FALSE]
  if (nrow(hit) != 1L)
    stop("ambiguous_or_missing_mapping: ", sym, " resolved to ", nrow(hit),
         " protein groups", call. = FALSE)
  if (!is.na(acc) && !grepl(acc, hit$member_accessions[1], fixed = TRUE))
    stop("uniprot_mismatch: ", sym, " requested ", acc, " but canonical is ",
         hit$member_accessions[1], call. = FALSE)
  if (!isTRUE(as.logical(hit$protein_level_claim_allowed[1])))
    stop("protein_level_claim_not_allowed: ", sym, call. = FALSE)
  data.frame(
    requested_symbol = sym,
    requested_uniprot = acc %||% NA_character_,
    ProteinGroupID = hit$ProteinGroupID[1],
    GeneSymbol = hit$official_gene_symbol[1],
    UniProt = hit$representative_accession[1],
    member_accessions = hit$member_accessions[1],
    original_identifier = hit$original_identifier[1],
    ambiguity_class = hit$protein_group_ambiguity_class[1],
    protein_level_claim_allowed = hit$protein_level_claim_allowed[1],
    stringsAsFactors = FALSE)
}
mapping <- do.call(rbind, Map(resolve_one, CANDIDATES$requested_symbol,
                              CANDIDATES$requested_uniprot))
rownames(mapping) <- NULL
if (anyDuplicated(mapping$ProteinGroupID))
  stop("ambiguous_mapping: two candidates resolved to the same protein group",
       call. = FALSE)

# ------------------------------------------------------- animal abundance
keep <- match(mapping$original_identifier, row_id)
if (anyNA(keep))
  stop("candidate absent from the animal-level matrix: ",
       paste(mapping$GeneSymbol[is.na(keep)], collapse = ", "), call. = FALSE)

abundance <- do.call(rbind, lapply(seq_len(nrow(mapping)), function(i) {
  v <- as.numeric(unlist(expr[keep[i], sample_meta$sample_id]))
  data.frame(
    ProteinGroupID = mapping$ProteinGroupID[i],
    GeneSymbol = mapping$GeneSymbol[i],
    UniProt = mapping$UniProt[i],
    sample_id = sample_meta$sample_id,
    AnimalID = sample_meta$AnimalID,
    ExpGroup = sample_meta$ExpGroup,
    spatial_unit = sample_meta$spatial_unit,
    abundance = v,
    stringsAsFactors = FALSE)
}))
abundance <- abundance %>%
  inner_join(SPATIAL, by = "spatial_unit") %>%
  mutate(dataset = DATASET,
         biological_unit = "animal",
         hemisphere_handling = "bilateral animal-level; equal-weight L/R mean where both present, single hemisphere retained otherwise, no imputation",
         abundance_type = "canonical animal-level value as supplied to the differential-abundance analysis") %>%
  arrange(GeneSymbol, spatial_order, ExpGroup, AnimalID) %>%
  select(dataset, ProteinGroupID, GeneSymbol, UniProt, spatial_unit,
         spatial_unit_display, region, layer, spatial_order, AnimalID, ExpGroup,
         abundance, sample_id, biological_unit, hemisphere_handling,
         abundance_type)

# no pseudoreplication: exactly one row per protein x animal x spatial unit
dup <- abundance %>% count(ProteinGroupID, AnimalID, spatial_unit) %>%
  filter(n > 1L)
if (nrow(dup))
  stop("pseudoreplication: duplicate protein x animal x spatial unit rows",
       call. = FALSE)

# ----------------------------------------------------------- SUS-RES effects
da_file_for <- function(unit) {
  key <- tolower(gsub("_", "", unit))
  file.path(DA_DIR, sprintf("%ssus_%sres.csv", key, key))
}
effects <- do.call(rbind, lapply(seq_len(nrow(SPATIAL)), function(j) {
  unit <- SPATIAL$spatial_unit[j]
  f <- da_file_for(unit)
  if (!file.exists(f))
    stop("missing_required_input: canonical SUS-RES file ", f, call. = FALSE)
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d <- d[match(mapping$ProteinGroupID, d$ProteinGroupID), , drop = FALSE]
  if (anyNA(d$ProteinGroupID))
    stop("candidate absent from canonical SUS-RES file for ", unit, call. = FALSE)
  data.frame(
    dataset = DATASET,
    ProteinGroupID = mapping$ProteinGroupID,
    GeneSymbol = mapping$GeneSymbol,
    UniProt = mapping$UniProt,
    spatial_unit = unit,
    spatial_unit_display = SPATIAL$spatial_unit_display[j],
    region = SPATIAL$region[j], layer = SPATIAL$layer[j],
    spatial_order = SPATIAL$spatial_order[j],
    contrast = "SUS - RES",
    log2FC = d$log2fc, aveExpr = d$aveExpr, t_statistic = d$t,
    p_value = d$pval, fdr_bh = d$padj,
    fdr_supported = d$padj <= 0.05,
    canonical_source_key = basename(f),
    stringsAsFactors = FALSE)
}))
effects <- effects %>% arrange(GeneSymbol, spatial_order)

# verification: the copied values must equal their canonical source exactly
for (j in seq_len(nrow(SPATIAL))) {
  unit <- SPATIAL$spatial_unit[j]
  d <- utils::read.csv(da_file_for(unit), stringsAsFactors = FALSE)
  d <- d[match(mapping$ProteinGroupID, d$ProteinGroupID), , drop = FALSE]
  e <- effects[effects$spatial_unit == unit, , drop = FALSE]
  e <- e[match(mapping$ProteinGroupID, e$ProteinGroupID), , drop = FALSE]
  for (col in c(log2FC = "log2fc", p_value = "pval", fdr_bh = "padj")) {
    got <- e[[names(which(c(log2FC = "log2fc", p_value = "pval",
                            fdr_bh = "padj") == col))]]
    if (!isTRUE(all.equal(got, d[[col]], tolerance = 0)))
      stop("da_reproduction_failed: ", unit, " column ", col, call. = FALSE)
  }
}

# -------------------------------------------------------------- provenance
git_sha <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
  error = function(e) NA_character_)
prov <- do.call(rbind, lapply(seq_len(nrow(effects)), function(i) data.frame(
  repository_commit = git_sha,
  contract_version = CONTRACT_VERSION,
  dataset = DATASET,
  ProteinGroupID = effects$ProteinGroupID[i],
  GeneSymbol = effects$GeneSymbol[i],
  UniProt = effects$UniProt[i],
  spatial_unit = effects$spatial_unit[i],
  canonical_da_source_key = effects$canonical_source_key[i],
  canonical_da_input = file.path("data/processed/02_id_mapping/mapped",
                                 DATASET, "forward/per_file",
                                 effects$canonical_source_key[i]),
  canonical_abundance_input = file.path(
    "data/processed/01_preprocessing/protigy_input_animal_level", DATASET,
    paste0(DATASET, "_animal_level.gct")),
  spatial_order_contract = "config/manuscript_spatial_order.yml",
  biological_unit = "animal",
  stringsAsFactors = FALSE)))

utils::write.csv(abundance,
  file.path(OUT, "immunostaining_candidates_animal_abundance.csv"),
  row.names = FALSE)
utils::write.csv(effects,
  file.path(OUT, "immunostaining_candidates_sus_res_effects.csv"),
  row.names = FALSE)
utils::write.csv(mapping,
  file.path(OUT, "immunostaining_candidates_mapping.csv"), row.names = FALSE)
utils::write.csv(prov,
  file.path(OUT, "immunostaining_candidates_provenance.csv"), row.names = FALSE)

cat("immunostaining candidate source data ->", OUT, "\n")
cat("  mapping rows:", nrow(mapping), "\n")
print(mapping[, c("GeneSymbol", "UniProt", "ProteinGroupID", "ambiguity_class")])
cat("  abundance rows:", nrow(abundance), "(",
    length(unique(abundance$ProteinGroupID)), "proteins x",
    length(unique(abundance$AnimalID)), "animals x",
    length(unique(abundance$spatial_unit)), "units )\n")
cat("  effect rows:", nrow(effects), "| FDR-supported:", sum(effects$fdr_supported), "\n")
for (g in unique(effects$GeneSymbol)) {
  e <- effects[effects$GeneSymbol == g, ]
  s <- e[which.max(abs(e$log2FC)), ]
  cat(sprintf("   %-9s strongest %-7s log2FC %+.3f FDR %.4g | negative %d/10 | FDR<=0.05 in %d\n",
              g, s$spatial_unit, s$log2FC, s$fdr_bh, sum(e$log2FC < 0),
              sum(e$fdr_supported)))
}
cat("  DA values verified identical to canonical source: TRUE\n")
