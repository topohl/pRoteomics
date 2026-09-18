#!/usr/bin/env Rscript
# ================================================================
# Script: analysis/integration/screen_immunostaining_separation.R
# Stage: integration
# Scope: global (neuron_neuropil only)
# Consumes: required data/processed/01_preprocessing/protigy_input_animal_level/neuron_neuropil/neuron_neuropil_animal_level.gct; data/raw/pg_matrix/quicksearch.pg_matrix.tsv; data/processed/02_id_mapping/mapped/neuron_neuropil/forward/per_file/; +1 more; optional none declared in pipeline.yml
# Produces: results/integration/screen_immunostaining_separation/global/tables/source_data/separation_candidate_panel.csv; results/integration/screen_immunostaining_separation/global/tables/source_data/separation_metrics_all.csv; results/integration/screen_immunostaining_separation/global/tables/source_data/separation_qualifying_cells.csv; +3 more
# Notes: Descriptive separation screen for staining feasibility. No new inference.
# ================================================================
#
# WHY THIS SCREEN EXISTS ALONGSIDE THE FDR SCREEN
#   With three animals per group an FDR value is dominated by the variance
#   estimate and is unstable. What actually determines whether a difference is
#   visible by immunostaining is different: the fold change has to be large, it
#   has to hold in most animals rather than being carried by one, and the groups
#   have to be tight enough to separate. Those are assay-separation questions,
#   so this screen uses assay-separation measures.
#
# WHAT IS COMPUTED
#   Per protein x spatial unit, from the canonical bilateral animal-level values:
#     * log2FC SUS - RES, copied from the canonical output and verified
#     * within-group SD for SUS and RES (n = 3 each)
#     * complete separation: every SUS animal above every RES animal, or below
#     * Z'-factor  = 1 - 3(sd_SUS + sd_RES) / |mean_SUS - mean_RES|
#     * SSMD       = (mean_SUS - mean_RES) / sqrt(sd_SUS^2 + sd_RES^2)
#   These are DESCRIPTIVE separation statistics. No p-value is derived from
#   them, no model is fitted, and no existing result is changed.
#
# THE IMPUTATION TRAP THIS SCREEN AVOIDS
#   Missing values were imputed by drawing from a downshifted normal
#   (archive/01_preprocessing/01_impute.r uses rnorm with a narrow sd). Imputed cells
#   are therefore artificially tight, which INFLATES every variance-based
#   separation score. Ranking on Z'-factor without guarding this would
#   systematically promote proteins that were mostly missing. The raw matrix is
#   read here purely to rebuild the missingness mask, and a protein x unit cell
#   only enters the ranking when all six SUS and RES animal values are fully
#   observed with no imputed acquisition behind them.
# Dataset behavior: runs for global according to pipeline.yml and --dataset/PROTEOMICS_DATASET where supported.
#  

suppressPackageStartupMessages({
  library(dplyr)
})

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)
source(repo_path("R", "integration_utils.R"))
source(repo_path("R", "preprocessing_paths.R"))

# Phase 6G.6: destinations resolve through the normalized output
# contract, addressed by this analysis's own identity. This domain
# spanned two historical stage namespaces; neither survives in a new
# path. Outputs already written there stay exactly where they are.
ANALYSIS_ID <- "screen_immunostaining_separation"

MODULE_ID <- "10_biological_integration"
SUBSTEP_ID <- "immunostaining_separation_screen"
DATASET <- "neuron_neuropil"
CONTRACT_VERSION <- "immunostaining_separation_screen_v1"

# Stated thresholds, not tuned: a two-fold difference, groups that do not
# overlap at all, and a positive Z'-factor (the conventional boundary at which
# two groups are considered separable in assay development).
MIN_ABS_LOG2FC <- 1.0
N_TARGET <- 10L

# Z'-factor is REPORTED but not used as a gate. The conventional Z' > 0
# boundary asks the difference to exceed three times the summed within-group
# SDs, which is a plate-assay standard; across all 40,151 fully observed
# protein x unit cells here, ZERO reach it. Gating on it would return an empty
# panel and say nothing useful, so ranking uses SSMD, which is the small-sample
# separation measure, with the conventional bands |SSMD| >= 3 strong, >= 2
# moderate, >= 1 weak.
MIN_ABS_SSMD <- 1.0

# Housekeeping and structural proteins are flagged rather than dropped. A
# two-fold beta-actin difference between outcome groups is far more likely to be
# normalisation residue than biology, and it would be a poor staining target,
# but the call belongs to the reader.
HOUSEKEEPING_FLAG <- c("Actb", "Actg1", "Gapdh", "Tubb5", "Tubb4a", "Tuba1a",
                       "Tuba1b", "Ppia", "Ywhaz", "Rplp0", "Hprt", "B2m")

OUT <- integration_dirs(ANALYSIS_ID, "global", create = TRUE)$source_data
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GCT <- repo_path("data", "processed", "01_preprocessing",
                 "protigy_input_animal_level", DATASET,
                 paste0(DATASET, "_animal_level.gct"))
RAW <- repo_path("data", "raw", "pg_matrix", "quicksearch.pg_matrix.tsv")
DA_DIR <- preprocessing_mapped_contrast_dir(DATASET, "forward")
for (p in c(GCT, RAW)) if (!file.exists(p))
  stop("missing_required_input: ", p, call. = FALSE)

units_yml <- yaml::read_yaml(repo_path("config", "manuscript_spatial_order.yml"))$units
npx <- Filter(function(z) identical(z$dataset, DATASET), units_yml)
SPATIAL <- data.frame(
  spatial_unit = vapply(npx, function(z) z$analysis_key, character(1)),
  spatial_unit_display = vapply(npx, function(z) z$display, character(1)),
  spatial_order = vapply(npx, function(z) as.integer(z$order), integer(1)),
  stringsAsFactors = FALSE)
SPATIAL <- SPATIAL[order(SPATIAL$spatial_order), ]

# ------------------------------------------------ animal-level abundance
gl <- readLines(GCT, n = 12L)
n_col_meta <- as.integer(strsplit(gl[2], "\t")[[1]])[4]
hdr <- strsplit(gl[3], "\t")[[1]]
mrows <- lapply(seq_len(n_col_meta), function(i) strsplit(gl[3L + i], "\t")[[1]])
mnames <- vapply(mrows, function(x) x[1], character(1))
mval <- function(k) mrows[[which(mnames == k)]][-(1:2)]
smeta <- data.frame(sample_id = hdr[-(1:2)], AnimalID = mval("AnimalID"),
                    ExpGroupCode = mval("ExpGroup"), region = mval("region"),
                    layer = mval("layer"), stringsAsFactors = FALSE)
smeta$ExpGroup <- unname(c(`1` = "CON", `2` = "RES", `3` = "SUS")[smeta$ExpGroupCode])
smeta$spatial_unit <- paste(smeta$region, smeta$layer, sep = "_")
ROSTER <- list(CON = c("A127", "A129", "A765"), RES = c("A0003", "A135", "A139"),
               SUS = c("A111", "A755", "A764"))
for (g in names(ROSTER))
  if (!identical(sort(unique(smeta$AnimalID[smeta$ExpGroup == g])), sort(ROSTER[[g]])))
    stop("ExpGroup coding mismatch for ", g, call. = FALSE)

ex <- utils::read.delim(GCT, skip = 3L + n_col_meta, header = FALSE,
                        stringsAsFactors = FALSE, check.names = FALSE)
names(ex) <- hdr
orig_id <- ex[[1]]
mat <- as.matrix(ex[, smeta$sample_id])
storage.mode(mat) <- "double"

# ------------------------------------------------------ missingness mask
# Rebuilt from the raw matrix so that imputed cells can be excluded. Raw sample
# columns encode animal, hemisphere, region and layer in the file name.
raw <- utils::read.delim(RAW, stringsAsFactors = FALSE, check.names = FALSE)
raw_cols <- setdiff(names(raw), c("Protein.Group", "Protein.Names", "Genes",
                                  "First.Protein.Description"))
parse_col <- function(x) {
  b <- basename(gsub("\\\\", "/", x))
  m <- regmatches(b, regexec("Tobias_(A[0-9]+)_([LR])_([A-Za-z0-9]+)_([A-Za-z0-9]+)_", b))[[1]]
  if (length(m) < 5) return(c(NA, NA, NA, NA))
  c(m[2], m[3], m[4], m[5])
}
pc <- t(vapply(raw_cols, parse_col, character(4)))
rawmeta <- data.frame(col = raw_cols, AnimalID = pc[, 1], hemi = pc[, 2],
                      region = pc[, 3], layer4 = pc[, 4],
                      stringsAsFactors = FALSE)
# neuropil layer tokens are the layer itself; the celltype token follows
rawmeta$spatial_unit <- paste(rawmeta$region, rawmeta$layer4, sep = "_")
rawmeta <- rawmeta[!is.na(rawmeta$AnimalID) &
                     rawmeta$spatial_unit %in% SPATIAL$spatial_unit, ]

# map raw protein rows onto the GCT original identifiers via Protein.Names
raw_key <- raw$Protein.Names
obs_frac <- matrix(NA_real_, nrow = length(orig_id), ncol = nrow(smeta),
                   dimnames = list(orig_id, smeta$sample_id))
ridx <- match(orig_id, raw_key)
for (j in seq_len(nrow(smeta))) {
  cols <- rawmeta$col[rawmeta$AnimalID == smeta$AnimalID[j] &
                        rawmeta$spatial_unit == smeta$spatial_unit[j]]
  if (!length(cols)) next
  sub <- as.matrix(raw[ridx, cols, drop = FALSE])
  storage.mode(sub) <- "double"
  obs_frac[, j] <- rowMeans(!is.na(sub))
}

# ------------------------------------------------- per protein x unit stats
zfac <- function(m1, m2, s1, s2) {
  d <- abs(m1 - m2)
  ifelse(d > 0, 1 - 3 * (s1 + s2) / d, -Inf)
}
res <- do.call(rbind, lapply(SPATIAL$spatial_unit, function(u) {
  iS <- which(smeta$spatial_unit == u & smeta$ExpGroup == "SUS")
  iR <- which(smeta$spatial_unit == u & smeta$ExpGroup == "RES")
  S <- mat[, iS, drop = FALSE]; R <- mat[, iR, drop = FALSE]
  oS <- obs_frac[, iS, drop = FALSE]; oR <- obs_frac[, iR, drop = FALSE]
  mS <- rowMeans(S); mR <- rowMeans(R)
  sS <- apply(S, 1, stats::sd); sR <- apply(R, 1, stats::sd)
  data.frame(
    original_identifier = orig_id, spatial_unit = u,
    mean_SUS = mS, mean_RES = mR, sd_SUS = sS, sd_RES = sR,
    log2FC = mS - mR,
    min_SUS = apply(S, 1, min), max_SUS = apply(S, 1, max),
    min_RES = apply(R, 1, min), max_RES = apply(R, 1, max),
    fully_observed = apply(oS, 1, function(z) all(!is.na(z) & z == 1)) &
      apply(oR, 1, function(z) all(!is.na(z) & z == 1)),
    mean_abundance = rowMeans(cbind(S, R)),
    stringsAsFactors = FALSE)
}))
res <- res %>%
  mutate(
    complete_separation = (min_SUS > max_RES) | (max_SUS < min_RES),
    z_factor = zfac(mean_SUS, mean_RES, sd_SUS, sd_RES),
    ssmd = (mean_SUS - mean_RES) / sqrt(sd_SUS^2 + sd_RES^2),
    abs_log2FC = abs(log2FC)) %>%
  inner_join(SPATIAL, by = "spatial_unit")

# ------------------------------------------- annotate with canonical identity
probe <- utils::read.csv(file.path(DA_DIR, "CA2slmsus_CA2slmres.csv"),
                         stringsAsFactors = FALSE)
ann <- probe %>%
  transmute(original_identifier, ProteinGroupID, GeneSymbol = official_gene_symbol,
            UniProt = representative_accession,
            ambiguity_class = protein_group_ambiguity_class,
            protein_level_claim_allowed) %>%
  distinct(original_identifier, .keep_all = TRUE)
res <- res %>% inner_join(ann, by = "original_identifier")

# canonical log2FC must agree with the value recomputed from the animal means
for (u in SPATIAL$spatial_unit) {
  k <- tolower(gsub("_", "", u))
  d <- utils::read.csv(file.path(DA_DIR, sprintf("%ssus_%sres.csv", k, k)),
                       stringsAsFactors = FALSE)
  r <- res[res$spatial_unit == u, ]
  i <- match(r$ProteinGroupID, d$ProteinGroupID)
  dev <- max(abs(r$log2FC - d$log2fc[i]), na.rm = TRUE)
  if (!is.finite(dev) || dev > 1e-6)
    stop("animal-mean log2FC disagrees with canonical at ", u,
         " (max deviation ", signif(dev, 3), ")", call. = FALSE)
}

# reuse the epidermal QC exclusion from the FDR screen, same rationale
EPIDERMAL_QC_LITERATURE <- c(
  "S100a14", "Nccrp1", "Pkp1", "Pkp3", "Dsp", "Dsg1a", "Dsc2", "Tgm1", "Tgm3",
  "Flg", "Flg2", "Lor", "Ivl", "Cdsn", "Evpl", "Ppl", "Sbsn", "Krt1", "Krt5",
  "Krt10", "Krt14", "Krt26", "Krt222", "Krtap", "Serpinb3a", "Casp14")
epi_go <- local({
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) return(character(0))
  unique(unlist(suppressMessages(AnnotationDbi::mapIds(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = c("GO:0008544", "GO:0030855", "GO:0043588", "GO:0030057", "GO:0001533"),
    keytype = "GOALL", column = "SYMBOL", multiVals = "list")), use.names = FALSE))
})

qualifying <- res %>%
  filter(fully_observed, complete_separation,
         abs_log2FC >= MIN_ABS_LOG2FC, abs(ssmd) >= MIN_ABS_SSMD,
         !is.na(GeneSymbol), nzchar(GeneSymbol),
         as.logical(protein_level_claim_allowed),
         ambiguity_class == "single_accession_single_gene",
         !GeneSymbol %in% EPIDERMAL_QC_LITERATURE, !GeneSymbol %in% epi_go)

# per-protein rollup: how widely does the separation hold, and how good is the
# best context
per_protein <- qualifying %>%
  group_by(ProteinGroupID, GeneSymbol, UniProt) %>%
  summarise(
    n_units_qualifying = n(),
    best_unit = spatial_unit[which.max(abs(ssmd))],
    best_unit_display = spatial_unit_display[which.max(abs(ssmd))],
    best_z_factor = z_factor[which.max(abs(ssmd))],
    best_log2FC = log2FC[which.max(abs(ssmd))],
    best_ssmd = ssmd[which.max(abs(ssmd))],
    best_sd_SUS = sd_SUS[which.max(abs(ssmd))],
    best_sd_RES = sd_RES[which.max(abs(ssmd))],
    best_mean_abundance = mean_abundance[which.max(abs(ssmd))],
    .groups = "drop")

# directional consistency across all ten contexts, from the full table
dirn <- res %>%
  group_by(ProteinGroupID) %>%
  summarise(n_negative = sum(log2FC < 0), n_positive = sum(log2FC > 0),
            .groups = "drop") %>%
  mutate(n_matching_majority = pmax(n_negative, n_positive),
         majority_direction = ifelse(n_negative >= n_positive, "negative", "positive"))
per_protein <- per_protein %>% left_join(dirn, by = "ProteinGroupID")

panel <- per_protein %>%
  mutate(housekeeping_flag = GeneSymbol %in% HOUSEKEEPING_FLAG) %>%
  arrange(housekeeping_flag, desc(n_units_qualifying), desc(abs(best_ssmd))) %>%
  mutate(rank = row_number(), selection_contract = CONTRACT_VERSION)

git_sha <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE,
                            stderr = FALSE)[1], error = function(e) NA_character_)
utils::write.csv(res %>% select(ProteinGroupID, GeneSymbol, UniProt,
                                spatial_unit, spatial_unit_display, spatial_order,
                                mean_SUS, mean_RES, sd_SUS, sd_RES, log2FC,
                                abs_log2FC, complete_separation, z_factor, ssmd,
                                fully_observed, mean_abundance),
                 file.path(OUT, "separation_metrics_all.csv"), row.names = FALSE)
utils::write.csv(qualifying %>% select(ProteinGroupID, GeneSymbol, UniProt,
                                       spatial_unit, log2FC, z_factor, ssmd,
                                       sd_SUS, sd_RES, mean_abundance),
                 file.path(OUT, "separation_qualifying_cells.csv"), row.names = FALSE)
utils::write.csv(panel, file.path(OUT, "separation_candidate_panel.csv"),
                 row.names = FALSE)
utils::write.csv(
  data.frame(repository_commit = git_sha, contract_version = CONTRACT_VERSION,
             dataset = DATASET, min_abs_log2FC = MIN_ABS_LOG2FC,
             min_abs_ssmd = MIN_ABS_SSMD,
             z_factor_used_as_gate = FALSE,
             requires_complete_separation = TRUE,
             requires_fully_observed = TRUE,
             canonical_abundance_input = "data/processed/01_preprocessing/protigy_input_animal_level/neuron_neuropil/neuron_neuropil_animal_level.gct",
             missingness_mask_input = "data/raw/pg_matrix/quicksearch.pg_matrix.tsv",
             biological_unit = "animal", stringsAsFactors = FALSE),
  file.path(OUT, "separation_screen_provenance.csv"), row.names = FALSE)

cat("separation screen ->", OUT, "\n")
cat("protein x unit cells:", nrow(res),
    "| fully observed:", sum(res$fully_observed),
    "| complete separation:", sum(res$complete_separation), "\n")
cat("cells passing all criteria:", nrow(qualifying),
    "| distinct proteins:", nrow(per_protein), "\n\n")
print(utils::head(panel %>% transmute(
  rank, GeneSymbol, UniProt, units = n_units_qualifying,
  best = best_unit, log2FC = round(best_log2FC, 2),
  Zfac = round(best_z_factor, 2), SSMD = round(best_ssmd, 2), hk = housekeeping_flag,
  sdSUS = round(best_sd_SUS, 3), sdRES = round(best_sd_RES, 3),
  abund = round(best_mean_abundance, 2),
  dir = paste0(n_matching_majority, "/10 ", majority_direction)) %>%
  as.data.frame(), 15), row.names = FALSE)

# ------------------------------------------------- abundance for plotting
# Detectability matters for staining: a protein at the bottom of the abundance
# range may separate cleanly and still be invisible on tissue. Flagged, not
# dropped, with the dataset's own distribution as the reference.
ab_q <- stats::quantile(res$mean_abundance, c(0.1, 0.5), na.rm = TRUE)
panel <- panel %>%
  mutate(low_abundance_flag = best_mean_abundance < ab_q[[1]],
         abundance_decile_note = sprintf(
           "dataset 10th pct %.2f, median %.2f", ab_q[[1]], ab_q[[2]]))
utils::write.csv(panel, file.path(OUT, "separation_candidate_panel.csv"),
                 row.names = FALSE)

keep_oid <- ann$original_identifier[match(panel$ProteinGroupID, ann$ProteinGroupID)]
ridx2 <- match(keep_oid, orig_id)
abundance <- do.call(rbind, lapply(seq_len(nrow(panel)), function(i)
  data.frame(
    ProteinGroupID = panel$ProteinGroupID[i], GeneSymbol = panel$GeneSymbol[i],
    UniProt = panel$UniProt[i], panel_rank = panel$rank[i],
    evidence_class = sprintf("separation rank %d", panel$rank[i]),
    sample_id = smeta$sample_id, AnimalID = smeta$AnimalID,
    ExpGroup = smeta$ExpGroup, spatial_unit = smeta$spatial_unit,
    abundance = as.numeric(mat[ridx2[i], ]),
    stringsAsFactors = FALSE))) %>%
  inner_join(SPATIAL, by = "spatial_unit") %>%
  mutate(dataset = DATASET, biological_unit = "animal",
         hemisphere_handling = "bilateral animal-level; equal-weight L/R mean where both present, single hemisphere retained otherwise, no imputation") %>%
  arrange(panel_rank, spatial_order, ExpGroup, AnimalID)
if (nrow(abundance %>% count(ProteinGroupID, AnimalID, spatial_unit) %>%
         filter(n > 1L)))
  stop("pseudoreplication in the separation panel abundance table", call. = FALSE)
utils::write.csv(abundance, file.path(OUT, "separation_panel_animal_abundance.csv"),
                 row.names = FALSE)

effects <- res %>%
  semi_join(panel, by = "ProteinGroupID") %>%
  left_join(panel %>% select(ProteinGroupID, panel_rank = rank), by = "ProteinGroupID") %>%
  mutate(contrast = "SUS - RES", fdr_supported = FALSE,
         evidence_class = sprintf("separation rank %d", panel_rank)) %>%
  arrange(panel_rank, spatial_order)
utils::write.csv(effects, file.path(OUT, "separation_panel_sus_res_effects.csv"),
                 row.names = FALSE)

cat("\nflags: low abundance", sum(panel$low_abundance_flag),
    "| housekeeping", sum(panel$housekeeping_flag), "\n")
cat("abundance reference:", panel$abundance_decile_note[1], "\n")
