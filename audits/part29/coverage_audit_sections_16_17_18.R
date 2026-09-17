# =====================================================================
# Part-29 audit, sections 16 / 17 / 18
# Effective gene-set coverage of the canonical ranked GSEA across
# datasets and spatial compartments.
#
# READ-ONLY with respect to every canonical analysis output.
# Writes only into results/tables/publication_audits/upstream_enrichment_v10/.
#
# Canonical evidence sources (never regenerated here):
#   - compareGO term tables  (evidence_source_family
#       canonical_compareGO_ranked_GSEA_GO_BP) -> setSize / NES / core_enrichment
#   - per-comparison clusterProfiler GSEA_BP_results_full.csv (cross-check)
#   - per-comparison protein_group_audits/collapsed_gene_input.csv
#       -> the ranked gene universe that was handed to gseGO()
#   - org.Mm.eg.db / GO.db (pinned 3.22.0) -> canonical GO-BP annotation
# =====================================================================

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

suppressMessages({
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(GO.db)
})

OUT_DIR <- "results/tables/publication_audits/upstream_enrichment_v10"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# Windows long-path helper (several canonical audit files sit deeper
# than MAX_PATH). Backslashes built from raw bytes to avoid escaping.
# ---------------------------------------------------------------------
BS <- rawToChar(as.raw(92))
abspath <- function(rel) {
  p <- if (grepl("^[A-Za-z]:", rel)) rel else file.path(getwd(), rel)
  p <- gsub(BS, "/", p, fixed = TRUE)
  if (nchar(p) < 250L) return(p)
  paste0(BS, BS, "?", BS, gsub("/", BS, p, fixed = TRUE))
}
read_csv_lp <- function(rel, ...) utils::read.csv(abspath(rel), stringsAsFactors = FALSE, ...)

log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# ---------------------------------------------------------------------
# 1. Manuscript-relevant GO terms: the 204 BP terms carried by the six
#    claim-eligible primary themes, plus the three exemplars.
# ---------------------------------------------------------------------
THEME_FILE <- file.path("results/tables/10_biological_integration",
                        "gsea_wgcna_concordance/global",
                        "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
theme <- read_csv_lp(THEME_FILE, colClasses = "character")
stopifnot(all(theme$evidence_source_family == "canonical_compareGO_ranked_GSEA_GO_BP"))

elig <- theme[theme$theme_claim_eligible == "TRUE", , drop = FALSE]
PRIMARY_THEMES <- c("synaptic_signaling_vesicle", "rna_processing_splicing_rnp",
                    "ribosome_translation", "autophagy_lysosome_endosome",
                    "mitochondrial_respiration_oxphos", "chromatin_organization")
stopifnot(setequal(unique(elig$theme_id), PRIMARY_THEMES))

EXEMPLARS <- c("GO:0099536", "GO:0006397", "GO:0006119")
go_ids <- sort(unique(c(elig$GO_ID, EXEMPLARS)))
log_line("manuscript-relevant GO terms:", length(go_ids),
         "| exemplars already inside the 204:", all(EXEMPLARS %in% elig$GO_ID))

# theme membership string per GO ID (7 terms sit in two themes)
tm <- unique(elig[, c("GO_ID", "theme_id")])
theme_of <- tapply(tm$theme_id, tm$GO_ID, function(x) paste(sort(unique(x)), collapse = ";"))

# manuscript spatial-unit / contrast vocabulary, keyed by source_comparison
map <- unique(theme[, c("dataset", "spatial_unit", "contrast",
                        "phenotype_contrast", "source_comparison")])
map <- map[order(map$dataset, map$spatial_unit, map$contrast), ]
stopifnot(nrow(map) == 54L, !any(duplicated(map$source_comparison)))

# ---------------------------------------------------------------------
# 2. Canonical GO-BP annotation, reproduced exactly the way
#    clusterProfiler:::get_GO_data() builds it for keyType = "SYMBOL":
#      mapIds(OrgDb, keys = <BP terms>, column = "SYMBOL",
#             keytype = "GOALL", multiVals = "list")
#    i.e. full is_a/part_of ancestor propagation (GOALL), mouse only.
# ---------------------------------------------------------------------
go_ont <- AnnotationDbi::Ontology(GO.db::GOTERM)
stopifnot(all(go_ids %in% names(go_ont[go_ont == "BP"])))
go_term_name <- AnnotationDbi::Term(GO.db::GOTERM)[go_ids]

g2g <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys = go_ids,
                                              column = "SYMBOL", keytype = "GOALL",
                                              multiVals = "list"))
annot <- lapply(g2g, function(x) unique(x[!is.na(x)]))
names(annot) <- go_ids
annot_n <- vapply(annot, length, integer(1))
log_line("GO-BP annotated mouse SYMBOLs per term: min", min(annot_n),
         "median", stats::median(annot_n), "max", max(annot_n))

# ---------------------------------------------------------------------
# 3. Ranked universes + protein-group mapping audits, per comparison.
# ---------------------------------------------------------------------
CP_ROOT <- "data/processed/04_differential_expression_enrichment/clusterProfiler"
DATASETS <- c("microglia", "neuron_neuropil", "neuron_soma")

comp_dir <- function(ds, unit_dir, comp) {
  file.path(CP_ROOT, ds, "phenotype_within_unit", unit_dir, comp)
}

universes  <- list()   # per comparison: measured gene SYMBOLs
multi_pg   <- list()   # per comparison: genes collapsed from > 1 protein group
amb_excl   <- list()   # per comparison: excluded-ambiguous protein group member symbols
comp_index <- list()

for (ds in DATASETS) {
  units <- list.dirs(file.path(CP_ROOT, ds, "phenotype_within_unit"),
                     recursive = FALSE, full.names = FALSE)
  for (u in units) {
    comps <- list.dirs(file.path(CP_ROOT, ds, "phenotype_within_unit", u),
                       recursive = FALSE, full.names = FALSE)
    for (cp in comps) {
      base <- comp_dir(ds, u, cp)
      aud  <- file.path(base, "protein_group_audits")
      gi   <- read_csv_lp(file.path(aud, "collapsed_gene_input.csv"))
      universes[[cp]] <- unique(as.character(gi$official_gene_symbol))
      multi_pg[[cp]]  <- unique(as.character(
        gi$official_gene_symbol[suppressWarnings(as.integer(gi$n_protein_groups_for_gene)) > 1L]))

      tr <- read_csv_lp(file.path(aud, "protein_group_to_gene_transformation_audit.csv"))
      ex <- tr[tr$eligibility_status != "eligible", , drop = FALSE]
      mgs <- ifelse(is.na(ex$member_gene_symbols), "", as.character(ex$member_gene_symbols))
      amb_excl[[cp]] <- list(
        n_excluded_groups = nrow(ex),
        n_excluded_with_symbols = sum(nzchar(mgs)),
        symbols = lapply(strsplit(mgs[nzchar(mgs)], ";", fixed = TRUE), trimws))

      comp_index[[length(comp_index) + 1L]] <- data.frame(
        dataset = ds, unit_dir = u, source_comparison = cp,
        ranked_universe_size = length(universes[[cp]]),
        n_genes_multi_protein_group = length(multi_pg[[cp]]),
        n_excluded_protein_groups = nrow(ex),
        n_excluded_protein_groups_with_symbols = sum(nzchar(mgs)),
        stringsAsFactors = FALSE)
    }
  }
}
comp_index <- do.call(rbind, comp_index)
stopifnot(nrow(comp_index) == 54L)
comp_index <- merge(comp_index, map, by = c("dataset", "source_comparison"), all.x = TRUE)
stopifnot(!any(is.na(comp_index$spatial_unit)))

# The measured universe is constant inside each dataset -> verify, do not assume.
universe_by_dataset <- list()
universe_constant <- logical(0)
for (ds in DATASETS) {
  ks <- comp_index$source_comparison[comp_index$dataset == ds]
  ref <- sort(universes[[ks[1]]])
  ok <- all(vapply(ks, function(k) identical(sort(universes[[k]]), ref), logical(1)))
  universe_constant[ds] <- ok
  universe_by_dataset[[ds]] <- ref
  log_line("dataset", ds, ": n_comparisons", length(ks),
           "| measured universe identical across all comparisons:", ok,
           "| size", length(ref))
}
stopifnot(all(universe_constant))
dataset_universe_size <- vapply(universe_by_dataset, length, integer(1))

# dataset-level mapping-ambiguity gene pools
amb_genes_by_dataset <- list()
multi_pg_by_dataset  <- list()
amb_group_syms_by_dataset <- list()
for (ds in DATASETS) {
  ks <- comp_index$source_comparison[comp_index$dataset == ds]
  syms <- unlist(lapply(ks, function(k) unlist(amb_excl[[k]]$symbols)), use.names = FALSE)
  amb_genes_by_dataset[[ds]] <- unique(syms[nzchar(syms)])
  multi_pg_by_dataset[[ds]]  <- unique(unlist(multi_pg[ks], use.names = FALSE))
  # per-group symbol vectors (deduplicated across identical comparisons)
  gl <- amb_excl[[ks[1]]]$symbols
  amb_group_syms_by_dataset[[ds]] <- gl
}

# ---------------------------------------------------------------------
# 4. Canonical GSEA results (compareGO tables) + cross-check against the
#    per-comparison clusterProfiler tables.
# ---------------------------------------------------------------------
gsea <- list()
for (ds in DATASETS) {
  p <- file.path("results/tables/04_differential_expression_enrichment/compareGO",
                 ds, "BP/phenotype_within_unit/all_route_units",
                 "compareGO_term_comparison.csv")
  d <- read_csv_lp(p)
  d$dataset <- ds
  gsea[[ds]] <- d
}
gsea <- do.call(rbind, gsea)
gsea <- gsea[gsea$result_type == "GSEA_GO" & gsea$ontology == "BP", , drop = FALSE]
log_line("canonical compareGO GSEA rows:", nrow(gsea),
         "| comparisons:", length(unique(gsea$comparison)))
stopifnot(setequal(unique(gsea$comparison), comp_index$source_comparison))

# cross-check one comparison per dataset against GSEA_BP_results_full.csv
xcheck <- list()
for (ds in DATASETS) {
  ci <- comp_index[comp_index$dataset == ds, ][1, ]
  full <- read_csv_lp(file.path(comp_dir(ds, ci$unit_dir, ci$source_comparison),
                                "GO/BP/GSEA_BP_results_full.csv"))
  a <- gsea[gsea$comparison == ci$source_comparison, ]
  m <- merge(a[, c("ID", "setSize", "NES", "core_enrichment")],
             full[, c("ID", "setSize", "NES", "core_enrichment")], by = "ID")
  xcheck[[ds]] <- data.frame(
    dataset = ds, comparison = ci$source_comparison,
    same_terms = setequal(a$ID, full$ID),
    n_terms = nrow(a),
    setSize_identical = all(m$setSize.x == m$setSize.y),
    max_abs_NES_diff = max(abs(m$NES.x - m$NES.y)),
    core_enrichment_identical = all(m$core_enrichment.x == m$core_enrichment.y),
    stringsAsFactors = FALSE)
}
xcheck <- do.call(rbind, xcheck)
print(xcheck)
stopifnot(all(xcheck$same_terms), all(xcheck$setSize_identical),
          all(xcheck$max_abs_NES_diff == 0), all(xcheck$core_enrichment_identical))

le_size <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  vapply(strsplit(x, "/", fixed = TRUE),
         function(p) sum(nzchar(p)), integer(1))
}
gsea$leading_edge_size <- le_size(gsea$core_enrichment)

key <- paste(gsea$comparison, gsea$ID, sep = "||")
stopifnot(!any(duplicated(key)))
gsea_idx <- setNames(seq_len(nrow(gsea)), key)

# ---------------------------------------------------------------------
# 5. SECTION 16 / 18 TABLE: effective set size by dataset x compartment
# ---------------------------------------------------------------------
MIN_GS <- 10L   # analysis_params$min_gs_size
MAX_GS <- 800L  # analysis_params$max_gs_size

obs_by_dataset <- matrix(NA_integer_, nrow = length(go_ids), ncol = length(DATASETS),
                         dimnames = list(go_ids, DATASETS))
amb_by_dataset <- matrix(NA_integer_, nrow = length(go_ids), ncol = length(DATASETS),
                         dimnames = list(go_ids, DATASETS))
amb_grp_by_dataset  <- amb_by_dataset
amb_lost_by_dataset <- amb_by_dataset
amb_mpg_by_dataset  <- amb_by_dataset
for (ds in DATASETS) {
  uni <- universe_by_dataset[[ds]]
  ambg <- amb_genes_by_dataset[[ds]]
  mpg  <- multi_pg_by_dataset[[ds]]
  grps <- amb_group_syms_by_dataset[[ds]]
  for (g in go_ids) {
    a <- annot[[g]]
    obs_by_dataset[g, ds] <- length(intersect(a, uni))
    lost <- setdiff(intersect(a, ambg), uni)
    collapsed <- intersect(a, mpg)
    amb_lost_by_dataset[g, ds] <- length(lost)
    amb_mpg_by_dataset[g, ds]  <- length(collapsed)
    amb_grp_by_dataset[g, ds]  <- sum(vapply(grps, function(s) any(s %in% a), logical(1)))
    amb_by_dataset[g, ds] <- length(unique(c(intersect(a, ambg), collapsed)))
  }
}

n_comp <- nrow(comp_index)
rows <- vector("list", n_comp)
for (i in seq_len(n_comp)) {
  ci <- comp_index[i, ]
  ds <- ci$dataset
  cp <- ci$source_comparison
  k  <- paste(cp, go_ids, sep = "||")
  j  <- gsea_idx[k]
  obs <- obs_by_dataset[, ds]
  eff <- gsea$setSize[j]
  rows[[i]] <- data.frame(
    GO_ID = go_ids,
    GO_description = unname(go_term_name[go_ids]),
    dataset = ds,
    spatial_unit = ci$spatial_unit,
    canonical_annotated_genes = unname(annot_n[go_ids]),
    genes_in_ranked_universe = unname(obs),
    effective_setSize = eff,
    fraction_annotated_observed = unname(obs) / unname(annot_n[go_ids]),
    leading_edge_size = gsea$leading_edge_size[j],
    dataset_measured_universe_size = unname(dataset_universe_size[ds]),
    contrast = ci$contrast,
    phenotype_contrast = ci$phenotype_contrast,
    source_comparison = cp,
    theme_ids = unname(theme_of[go_ids]),
    is_exemplar = go_ids %in% EXEMPLARS,
    tested_in_gsea = !is.na(eff),
    gsea_size_filter_reason = ifelse(!is.na(eff), "tested",
                              ifelse(obs < MIN_GS, "below_min_gs_size_10",
                              ifelse(obs > MAX_GS, "above_max_gs_size_800",
                                     "size_eligible_but_absent_from_stored_result"))),
    NES = gsea$NES[j],
    enrichmentScore = gsea$enrichmentScore[j],
    pvalue = gsea$pvalue[j],
    p_adjust = gsea$p.adjust[j],
    leading_edge_tags = gsea$leading_edge[j],
    fraction_leading_edge_of_setSize = gsea$leading_edge_size[j] / eff,
    ranked_universe_size_this_comparison = ci$ranked_universe_size,
    mapping_ambiguity_count = unname(amb_by_dataset[, ds]),
    mapping_ambiguity_genes_lost_to_excluded_groups = unname(amb_lost_by_dataset[, ds]),
    mapping_ambiguity_genes_collapsed_multi_protein_group = unname(amb_mpg_by_dataset[, ds]),
    mapping_ambiguity_excluded_protein_groups = unname(amb_grp_by_dataset[, ds]),
    stringsAsFactors = FALSE)
}
cov <- do.call(rbind, rows)
rownames(cov) <- NULL

# --- integrity checks -------------------------------------------------
tested <- cov[cov$tested_in_gsea, ]
chk_setsize_identity <- all(tested$effective_setSize == tested$genes_in_ranked_universe)
untested <- cov[!cov$tested_in_gsea, ]
chk_untested_explained <- all(untested$gsea_size_filter_reason %in%
                                c("below_min_gs_size_10", "above_max_gs_size_800"))
log_line("effective_setSize == |GO-annotated INTERSECT ranked universe| for all tested rows:",
         chk_setsize_identity)
log_line("every untested manuscript term explained by the min/max gene-set-size filter:",
         chk_untested_explained)
log_line("untested reasons:", paste(names(table(untested$gsea_size_filter_reason)),
                                    table(untested$gsea_size_filter_reason),
                                    collapse = "; "))
stopifnot(chk_setsize_identity)

# Size-eligible manuscript terms that are nevertheless missing from the stored
# canonical result: record them explicitly instead of silently dropping them.
anom <- cov[cov$gsea_size_filter_reason == "size_eligible_but_absent_from_stored_result",
            c("GO_ID", "GO_description", "theme_ids", "dataset", "spatial_unit", "contrast",
              "source_comparison", "canonical_annotated_genes", "genes_in_ranked_universe")]
if (nrow(anom)) {
  cat("\n--- size-eligible manuscript terms absent from the stored GSEA result ---\n")
  print(anom, row.names = FALSE)
}
# same check over ALL BP terms, not just manuscript terms
all_anom <- do.call(rbind, lapply(DATASETS, function(ds) {
  d <- gsea[gsea$dataset == ds, ]
  nc <- length(unique(d$comparison))
  tt2 <- table(d$ID)
  ids <- names(tt2)[tt2 < nc]
  if (!length(ids)) return(NULL)
  do.call(rbind, lapply(ids, function(m) {
    s <- d[d$ID == m, ]
    data.frame(dataset = ds, GO_ID = m, GO_description = s$Description[1],
               setSize = paste(unique(s$setSize), collapse = ";"),
               present_in = nrow(s), of_comparisons = nc,
               missing_from = paste(setdiff(unique(d$comparison), s$comparison), collapse = ";"),
               in_manuscript_204 = m %in% go_ids,
               stringsAsFactors = FALSE)
  }))
}))
if (!is.null(all_anom)) {
  cat("\n--- ALL GO-BP terms not present in every comparison of their dataset ---\n")
  print(all_anom, row.names = FALSE)
  write.csv(all_anom, file.path(OUT_DIR, "gsea_terms_missing_from_some_comparisons.csv"),
            row.names = FALSE)
}

write.csv(cov, file.path(OUT_DIR, "gsea_effective_setsize_by_dataset.csv"), row.names = FALSE)
log_line("wrote gsea_effective_setsize_by_dataset.csv :", nrow(cov), "rows")

# ---------------------------------------------------------------------
# 6. SECTION 16 TABLE: three-exemplar coverage
# ---------------------------------------------------------------------
ex <- cov[cov$GO_ID %in% EXEMPLARS, ]
ex_ds_stats <- do.call(rbind, lapply(split(ex, list(ex$GO_ID, ex$dataset), drop = TRUE),
  function(d) data.frame(
    GO_ID = d$GO_ID[1], dataset = d$dataset[1],
    leading_edge_count_dataset_min = min(d$leading_edge_size),
    leading_edge_count_dataset_median = stats::median(d$leading_edge_size),
    leading_edge_count_dataset_max = max(d$leading_edge_size),
    stringsAsFactors = FALSE)))
exemplar <- merge(ex, ex_ds_stats, by = c("GO_ID", "dataset"), all.x = TRUE)
exemplar <- data.frame(
  GO_ID = exemplar$GO_ID,
  GO_description = exemplar$GO_description,
  dataset = exemplar$dataset,
  spatial_unit = exemplar$spatial_unit,
  contrast = exemplar$contrast,
  GO_annotated_gene_count = exemplar$canonical_annotated_genes,
  measured_eligible_gene_count = exemplar$genes_in_ranked_universe,
  effective_setSize = exemplar$effective_setSize,
  leading_edge_count = exemplar$leading_edge_size,
  fraction_measured = exemplar$genes_in_ranked_universe / exemplar$canonical_annotated_genes,
  fraction_leading_edge = exemplar$leading_edge_size / exemplar$effective_setSize,
  mapping_ambiguity_count = exemplar$mapping_ambiguity_count,
  mapping_ambiguity_genes_lost_to_excluded_groups =
    exemplar$mapping_ambiguity_genes_lost_to_excluded_groups,
  mapping_ambiguity_genes_collapsed_multi_protein_group =
    exemplar$mapping_ambiguity_genes_collapsed_multi_protein_group,
  mapping_ambiguity_excluded_protein_groups =
    exemplar$mapping_ambiguity_excluded_protein_groups,
  dataset_measured_universe_size = exemplar$dataset_measured_universe_size,
  NES = exemplar$NES, pvalue = exemplar$pvalue, p_adjust = exemplar$p_adjust,
  tested_in_gsea = exemplar$tested_in_gsea,
  leading_edge_count_dataset_min = exemplar$leading_edge_count_dataset_min,
  leading_edge_count_dataset_median = exemplar$leading_edge_count_dataset_median,
  leading_edge_count_dataset_max = exemplar$leading_edge_count_dataset_max,
  stringsAsFactors = FALSE)
exemplar <- exemplar[order(exemplar$GO_ID, exemplar$dataset,
                           exemplar$spatial_unit, exemplar$contrast), ]
write.csv(exemplar, file.path(OUT_DIR, "three_exemplar_gsea_coverage.csv"), row.names = FALSE)
log_line("wrote three_exemplar_gsea_coverage.csv :", nrow(exemplar), "rows")

# ---------------------------------------------------------------------
# 7. SECTION 17: cross-compartment coverage asymmetry vs NES
#    "compartment" = one of the 18 spatial units. Ratios are computed
#    per GO term within a contrast, over the compartments where the term
#    was actually tested. NOTHING IS NORMALISED.
# ---------------------------------------------------------------------
tt <- cov[cov$tested_in_gsea, ]
sp <- split(tt, list(tt$GO_ID, tt$contrast), drop = TRUE)

asym <- do.call(rbind, lapply(sp, function(d) {
  smin <- min(d$effective_setSize); smax <- max(d$effective_setSize)
  imin <- which(d$effective_setSize == smin); imax <- which(d$effective_setSize == smax)
  ds_stats <- function(ds) {
    z <- d[d$dataset == ds, ]
    if (!nrow(z)) return(c(NA_real_, NA_real_))
    c(z$effective_setSize[1], stats::median(z$NES))
  }
  mg <- ds_stats("microglia"); np <- ds_stats("neuron_neuropil"); so <- ds_stats("neuron_soma")
  # within-dataset (same-dataset compartments) set-size spread
  wd <- vapply(split(d$effective_setSize, d$dataset),
               function(v) if (length(v)) max(v) / min(v) else NA_real_, numeric(1))
  data.frame(
    GO_ID = d$GO_ID[1],
    GO_description = d$GO_description[1],
    theme_ids = d$theme_ids[1],
    contrast = d$contrast[1],
    is_exemplar = d$is_exemplar[1],
    n_compartments_tested = nrow(d),
    n_datasets_tested = length(unique(d$dataset)),
    canonical_annotated_genes = d$canonical_annotated_genes[1],
    effective_setSize_min = smin,
    effective_setSize_max = smax,
    setSize_ratio_max_over_min = smax / smin,
    compartments_at_min_setSize = paste(sort(unique(d$spatial_unit[imin])), collapse = ";"),
    compartments_at_max_setSize = paste(sort(unique(d$spatial_unit[imax])), collapse = ";"),
    datasets_at_min_setSize = paste(sort(unique(d$dataset[imin])), collapse = ";"),
    datasets_at_max_setSize = paste(sort(unique(d$dataset[imax])), collapse = ";"),
    NES_at_min_setSize_median = stats::median(d$NES[imin]),
    NES_at_min_setSize_range = paste(sprintf("%.4f", range(d$NES[imin])), collapse = " to "),
    NES_at_max_setSize_median = stats::median(d$NES[imax]),
    NES_at_max_setSize_range = paste(sprintf("%.4f", range(d$NES[imax])), collapse = " to "),
    NES_min = min(d$NES), NES_max = max(d$NES),
    NES_range = max(d$NES) - min(d$NES),
    NES_sd = stats::sd(d$NES),
    NES_sign_flip_across_compartments = any(d$NES > 0) && any(d$NES < 0),
    max_within_dataset_setSize_ratio = max(wd, na.rm = TRUE),
    setSize_microglia = mg[1], NES_median_microglia = mg[2],
    setSize_neuron_neuropil = np[1], NES_median_neuron_neuropil = np[2],
    setSize_neuron_soma = so[1], NES_median_neuron_soma = so[2],
    stringsAsFactors = FALSE)
}))
rownames(asym) <- NULL
asym <- asym[asym$n_compartments_tested > 1L, ]

# --- threshold, fixed from the distribution shape BEFORE inspecting
#     which terms it selects: Tukey upper fence Q3 + 1.5*IQR of the
#     observed max/min set-size ratios.
r <- asym$setSize_ratio_max_over_min
q <- stats::quantile(r, c(0, .10, .25, .50, .75, .90, .95, .99, 1))
IQR_r <- unname(q["75%"] - q["25%"])
THRESH <- unname(q["75%"]) + 1.5 * IQR_r
asym$coverage_asymmetry_threshold_rule <- "Tukey upper fence: Q3 + 1.5*IQR of observed max/min setSize ratios"
asym$coverage_asymmetry_threshold_value <- THRESH
asym$flag_large_coverage_asymmetry <- r >= THRESH
asym$flag_ratio_ge_1p5 <- r >= 1.5
asym$flag_ratio_ge_2p0 <- r >= 2.0
asym <- asym[order(-asym$setSize_ratio_max_over_min, asym$GO_ID, asym$contrast), ]
write.csv(asym, file.path(OUT_DIR, "gsea_cross_compartment_setsize_asymmetry.csv"),
          row.names = FALSE)
log_line("wrote gsea_cross_compartment_setsize_asymmetry.csv :", nrow(asym), "rows")

cat("\n--- set-size ratio distribution across compartments (per GO term x contrast) ---\n")
print(round(q, 4))
cat("Tukey upper fence threshold =", round(THRESH, 4), "\n")
cat("flagged (>= fence):", sum(asym$flag_large_coverage_asymmetry),
    "of", nrow(asym), sprintf("(%.1f%%)", 100 * mean(asym$flag_large_coverage_asymmetry)), "\n")
cat("ratio >= 1.5:", sum(asym$flag_ratio_ge_1p5), " ratio >= 2.0:", sum(asym$flag_ratio_ge_2p0), "\n")
cat("max within-dataset setSize ratio over all terms:",
    max(asym$max_within_dataset_setSize_ratio), "\n")

# pooled coverage-vs-NES association (diagnostic, nothing normalised)
sp_all <- suppressWarnings(stats::cor(tt$effective_setSize, tt$NES, method = "spearman"))
sp_abs <- suppressWarnings(stats::cor(tt$effective_setSize, abs(tt$NES), method = "spearman"))
cat("pooled Spearman(setSize, NES) =", round(sp_all, 4),
    "| Spearman(setSize, |NES|) =", round(sp_abs, 4), "\n")

# ---------------------------------------------------------------------
# 8. Provenance / summary statistics
# ---------------------------------------------------------------------
fr <- cov$fraction_annotated_observed
# unweighted: one value per distinct GO term x dataset pair (612 pairs)
fp <- unique(cov[, c("GO_ID", "dataset", "fraction_annotated_observed")])$fraction_annotated_observed
prov <- data.frame(
  item = c("git_HEAD", "R_version", "org.Mm.eg.db", "GO.db", "clusterProfiler", "DOSE",
           "theme_table", "theme_table_rows", "manuscript_GO_terms",
           "n_datasets", "n_spatial_units", "n_comparisons",
           "universe_microglia", "universe_neuron_neuropil", "universe_neuron_soma",
           "universe_union_all_datasets",
           "min_gs_size", "max_gs_size",
           "rows_gsea_effective_setsize_by_dataset", "rows_three_exemplar_gsea_coverage",
           "rows_cross_compartment_asymmetry",
           "fraction_annotated_observed_min", "fraction_annotated_observed_q25",
           "fraction_annotated_observed_median", "fraction_annotated_observed_mean",
           "fraction_annotated_observed_q75", "fraction_annotated_observed_max",
           "n_distinct_GO_x_dataset_pairs",
           "pairwise_fraction_measured_min", "pairwise_fraction_measured_q25",
           "pairwise_fraction_measured_median", "pairwise_fraction_measured_mean",
           "pairwise_fraction_measured_q75", "pairwise_fraction_measured_max",
           "pct_pairs_fraction_measured_ge_0.50", "pct_pairs_fraction_measured_ge_0.30",
           "pct_pairs_fraction_measured_lt_0.20",
           "setsize_ratio_p50", "setsize_ratio_p75", "setsize_ratio_p90",
           "setsize_ratio_p95", "setsize_ratio_p99", "setsize_ratio_max",
           "coverage_asymmetry_threshold", "n_flagged_large_asymmetry",
           "max_within_dataset_setSize_ratio",
           "pooled_spearman_setSize_NES", "pooled_spearman_setSize_absNES",
           "effective_setSize_equals_annotated_intersect_universe",
           "untested_terms_all_explained_by_size_filter"),
  value = c("f556f08", R.version.string,
            as.character(packageVersion("org.Mm.eg.db")),
            as.character(packageVersion("GO.db")),
            as.character(packageVersion("clusterProfiler")),
            as.character(packageVersion("DOSE")),
            THEME_FILE, nrow(theme), length(go_ids),
            length(DATASETS), length(unique(comp_index$spatial_unit)), nrow(comp_index),
            dataset_universe_size[["microglia"]],
            dataset_universe_size[["neuron_neuropil"]],
            dataset_universe_size[["neuron_soma"]],
            length(unique(unlist(universe_by_dataset))),
            MIN_GS, MAX_GS,
            nrow(cov), nrow(exemplar), nrow(asym),
            sprintf("%.6f", min(fr)), sprintf("%.6f", unname(stats::quantile(fr, .25))),
            sprintf("%.6f", stats::median(fr)), sprintf("%.6f", mean(fr)),
            sprintf("%.6f", unname(stats::quantile(fr, .75))), sprintf("%.6f", max(fr)),
            length(fp),
            sprintf("%.6f", min(fp)), sprintf("%.6f", unname(stats::quantile(fp, .25))),
            sprintf("%.6f", stats::median(fp)), sprintf("%.6f", mean(fp)),
            sprintf("%.6f", unname(stats::quantile(fp, .75))), sprintf("%.6f", max(fp)),
            sprintf("%.4f", 100 * mean(fp >= 0.50)), sprintf("%.4f", 100 * mean(fp >= 0.30)),
            sprintf("%.4f", 100 * mean(fp < 0.20)),
            sprintf("%.6f", unname(q["50%"])), sprintf("%.6f", unname(q["75%"])),
            sprintf("%.6f", unname(q["90%"])), sprintf("%.6f", unname(q["95%"])),
            sprintf("%.6f", unname(q["99%"])), sprintf("%.6f", unname(q["100%"])),
            sprintf("%.6f", THRESH), sum(asym$flag_large_coverage_asymmetry),
            sprintf("%.6f", max(asym$max_within_dataset_setSize_ratio)),
            sprintf("%.6f", sp_all), sprintf("%.6f", sp_abs),
            chk_setsize_identity, chk_untested_explained),
  stringsAsFactors = FALSE)
write.csv(prov, file.path(OUT_DIR, "gsea_coverage_audit_provenance.csv"), row.names = FALSE)

cat("\n--- fraction of annotated genes actually measured (all 204 terms x 54 comparisons) ---\n")
print(round(stats::quantile(fr, c(0, .05, .25, .5, .75, .95, 1)), 4))
cat("mean =", round(mean(fr), 4), "\n")
cat("\n--- per dataset ---\n")
for (ds in DATASETS) {
  f <- cov$fraction_annotated_observed[cov$dataset == ds]
  cat(sprintf("%-16s median %.4f  IQR %.4f-%.4f  min %.4f  max %.4f\n", ds,
              stats::median(f), stats::quantile(f, .25), stats::quantile(f, .75),
              min(f), max(f)))
}

cat("\n--- exemplar coverage (dataset level; setSize is constant within a dataset) ---\n")
exs <- unique(exemplar[, c("GO_ID", "GO_description", "dataset", "GO_annotated_gene_count",
                           "measured_eligible_gene_count", "effective_setSize",
                           "fraction_measured", "mapping_ambiguity_count",
                           "leading_edge_count_dataset_min", "leading_edge_count_dataset_median",
                           "leading_edge_count_dataset_max")])
exs$fraction_measured <- round(exs$fraction_measured, 4)
print(exs, row.names = FALSE)

cat("\n--- top 15 cross-compartment coverage asymmetries ---\n")
print(head(asym[, c("GO_ID", "GO_description", "contrast", "n_compartments_tested",
                    "effective_setSize_min", "effective_setSize_max",
                    "setSize_ratio_max_over_min", "datasets_at_min_setSize",
                    "datasets_at_max_setSize", "NES_at_min_setSize_median",
                    "NES_at_max_setSize_median", "NES_range",
                    "flag_large_coverage_asymmetry")], 15), row.names = FALSE)

cat("\n--- terms NOT tested in every compartment ---\n")
nt <- unique(cov[!cov$tested_in_gsea, c("GO_ID", "GO_description", "dataset",
                                        "genes_in_ranked_universe", "canonical_annotated_genes",
                                        "gsea_size_filter_reason")])
print(nt[order(nt$GO_ID, nt$dataset), ], row.names = FALSE)

log_line("DONE")
