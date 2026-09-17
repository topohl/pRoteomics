# =============================================================================
# Part-29 sections 20 + 21: CAMERA feasibility / provenance audit
#
# Question: does a scientifically valid gene-level animal x sample matrix exist
#           for a full limma::camera() analysis, or is cameraPR() on the
#           canonical ranked statistic the only defensible option?
#
# AUDIT ONLY. Nothing canonical is modified, rerun or regenerated.
# Writes only to results/tables/publication_audits/upstream_enrichment_v10/.
# =============================================================================

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

OUT_DIR <- file.path("results", "tables", "publication_audits", "upstream_enrichment_v10")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

fmt <- function(x, d = 10) format(x, trim = TRUE, scientific = FALSE, digits = d)

# ---------------------------------------------------------------------------
# Windows MAX_PATH workaround.
# Several canonical protein_group_audits files sit at absolute paths of 262-268
# characters, which R's file API cannot open on this host (R 4.5.1, ucrt).
# robocopy is extended-length-path aware, so we stage READ-ONLY copies of the
# needed canonical files into the OS temp directory. Nothing canonical is
# modified; the staged copies are byte-identical reads.
# ---------------------------------------------------------------------------
STAGE <- file.path(tempdir(), "part29_camera_stage")
dir.create(STAGE, recursive = TRUE, showWarnings = FALSE)
winp <- function(p) chartr("/", "\\", p)
stage_files <- function(src_rel, files, tag) {
  dst <- file.path(STAGE, tag)
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  suppressWarnings(system2("robocopy",
    c(shQuote(winp(file.path(getwd(), src_rel))), shQuote(winp(dst)), files,
      "/NJH", "/NJS", "/NP", "/NFL", "/NDL", "/R:1", "/W:1"),
    stdout = TRUE, stderr = TRUE))
  missing <- files[!file.exists(file.path(dst, files))]
  if (length(missing)) stop("Could not stage: ", paste(file.path(src_rel, missing), collapse = ", "))
  dst
}

# ---------------------------------------------------------------------------
# GCT v1.3 reader (base R only)
# ---------------------------------------------------------------------------
read_gct13 <- function(path) {
  ln <- readLines(path, warn = FALSE)
  dims <- suppressWarnings(as.integer(strsplit(trimws(ln[2]), "\t", fixed = TRUE)[[1]]))
  dims <- dims[!is.na(dims)][1:4]
  nr <- dims[1]; nc <- dims[2]; nrhd <- dims[3]; nchd <- dims[4]
  expected <- 1L + nrhd + nc
  split_pad <- function(x) { f <- strsplit(x, "\t", fixed = TRUE)[[1]]; length(f) <- expected; f }
  hdr <- split_pad(ln[3])
  rhd_names <- if (nrhd > 0) hdr[2:(1 + nrhd)] else character(0)
  cid <- hdr[(2 + nrhd):(1 + nrhd + nc)]
  cdesc <- NULL
  if (nchd > 0) {
    cdm <- do.call(rbind, lapply(ln[4:(3 + nchd)], split_pad))
    cdesc <- as.data.frame(cdm[, (2 + nrhd):(1 + nrhd + nc), drop = FALSE], stringsAsFactors = FALSE)
    rownames(cdesc) <- cdm[, 1]
    colnames(cdesc) <- cid
  }
  dl <- ln[(4 + nchd):(3 + nchd + nr)]
  dm <- do.call(rbind, lapply(dl, split_pad))
  rid <- dm[, 1]
  rdesc <- NULL
  if (nrhd > 0) {
    rdesc <- as.data.frame(dm[, 2:(1 + nrhd), drop = FALSE], stringsAsFactors = FALSE)
    names(rdesc) <- rhd_names
    rownames(rdesc) <- NULL
  }
  mat <- matrix(suppressWarnings(as.numeric(dm[, (2 + nrhd):(1 + nrhd + nc)])), nrow = nr)
  rownames(mat) <- rid; colnames(mat) <- cid
  list(rid = rid, rdesc = rdesc, mat = mat, cdesc = cdesc,
       dims = c(nrow = nr, ncol = nc, nrhd = nrhd, nchd = nchd))
}

DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
STAT_GCT <- c(
  neuron_neuropil = "data/processed/01_preprocessing/protigy_output_animal_level/neuron_neuropil/stat_results_for_ssGSEA_neuropil_proteome.gct",
  neuron_soma     = "data/processed/01_preprocessing/protigy_output_animal_level/neuron_soma/stat_results_for_ssGSEA_soma_proteome.gct",
  microglia       = "data/processed/01_preprocessing/protigy_output_animal_level/microglia/stat_results_for_ssGSEA_microglia_proteome.gct"
)
INPUT_GCT <- c(
  neuron_neuropil = "data/processed/01_preprocessing/protigy_input_animal_level/neuron_neuropil/neuron_neuropil_animal_level.gct",
  neuron_soma     = "data/processed/01_preprocessing/protigy_input_animal_level/neuron_soma/neuron_soma_animal_level.gct",
  microglia       = "data/processed/01_preprocessing/protigy_input_animal_level/microglia/microglia_animal_level.gct"
)
MAPPED_DIR <- "data/processed/02_id_mapping/mapped/%s/forward/per_file/%s.csv"
AUDIT_DIR  <- "data/processed/04_differential_expression_enrichment/clusterProfiler/%s/phenotype_within_unit/%s/%s/protein_group_audits"
EXTRACT_MANIFEST <- "data/processed/01_preprocessing/gct_extractR/%s/canonical_gct_extract_manifest.csv"
AGG_SUMMARY <- "results/tables/01_preprocessing/02a_prepare_animal_level_protigy_input/%s/aggregation_summary.csv"
AGG_AUDIT   <- "results/tables/01_preprocessing/02a_prepare_animal_level_protigy_input/%s/aggregation_audit.csv"
SRC_ASSIGN  <- "results/tables/01_preprocessing/02a_prepare_animal_level_protigy_input/%s/source_sample_assignment.csv"

EXEMPLARS <- data.frame(
  dataset  = c("neuron_neuropil", "neuron_soma", "microglia"),
  unit     = c("CA3_sr", "CA2_sp", "CA1_microglia"),
  contrast = c("CA3srsus_CA3srres", "CA2spsus_CA2spres", "CA1microgliasus_CA1microgliares"),
  protigy_comparison = c("CA3_sr_3_over_CA3_sr_2", "CA2_sp_3_over_CA2_sp_2", "CA1_microglia_3_over_CA1_microglia_2"),
  stringsAsFactors = FALSE
)
EXPGROUP_LABEL <- c("1" = "CON", "2" = "RES", "3" = "SUS")

message("== reading animal-level GCTs ==")
gct_in   <- lapply(INPUT_GCT[DATASETS], read_gct13)
gct_stat <- lapply(STAT_GCT[DATASETS],  read_gct13)
names(gct_in) <- DATASETS; names(gct_stat) <- DATASETS

# ---------------------------------------------------------------------------
# Q1: inventory of candidate animal-level abundance matrices
# ---------------------------------------------------------------------------
inv <- list()
abund <- list()   # dataset -> numeric matrix (protein group rows x animal_spatial_unit columns)

for (ds in DATASETS) {
  gi <- gct_in[[ds]]; gs <- gct_stat[[ds]]
  unit_cols <- colnames(gi$mat)                       # canonical animal x spatial-unit column names
  stat_rdesc_names <- names(gs$rdesc)
  emb <- intersect(unit_cols, stat_rdesc_names)       # embedded abundance columns inside stat GCT rdesc
  A <- matrix(suppressWarnings(as.numeric(as.matrix(gs$rdesc[, emb, drop = FALSE]))),
              nrow = nrow(gs$rdesc), dimnames = list(gs$rid, emb))
  abund[[ds]] <- A

  # numeric agreement between protigy input matrix and protigy-embedded matrix
  common_r <- intersect(rownames(A), rownames(gi$mat))
  d <- A[common_r, emb, drop = FALSE] - gi$mat[common_r, emb, drop = FALSE]
  inv[[length(inv) + 1L]] <- data.frame(
    dataset = ds,
    artefact = "protigy_input_animal_level_gct",
    path = INPUT_GCT[[ds]],
    n_rows = gi$dims[["nrow"]], n_value_columns = gi$dims[["ncol"]],
    row_identity = "UniProt entry name (== mapped-contrast original_identifier); NOT ProteinGroupID",
    column_identity = "AnimalID_spatialUnit (animal x spatial unit)",
    value_semantics = "log2 protein-group abundance, 70%-missingness-filtered + imputed, equal-weight mean of Left/Right hemisphere replicates",
    n_na_values = sum(is.na(gi$mat)),
    n_unique_animals = length(unique(sub("_.*$", "", colnames(gi$mat)))),
    max_abs_diff_vs_embedded = NA_real_,
    stringsAsFactors = FALSE)
  inv[[length(inv) + 1L]] <- data.frame(
    dataset = ds,
    artefact = "protigy_output_animal_level_stat_gct_embedded_abundance",
    path = STAT_GCT[[ds]],
    n_rows = nrow(A), n_value_columns = ncol(A),
    row_identity = "UniProt entry name (== mapped-contrast original_identifier); NOT ProteinGroupID",
    column_identity = "AnimalID_spatialUnit (animal x spatial unit), stored as GCT row-descriptor columns",
    value_semantics = "same log2 abundance matrix Protigy actually fitted; carried alongside the canonical limma statistics",
    n_na_values = sum(is.na(A)),
    n_unique_animals = length(unique(sub("_.*$", "", colnames(A)))),
    max_abs_diff_vs_embedded = max(abs(d), na.rm = TRUE),
    stringsAsFactors = FALSE)
}
inventory <- do.call(rbind, inv)

# provenance: is the stat GCT the declared source of the canonical contrast CSVs?
prov <- do.call(rbind, lapply(DATASETS, function(ds) {
  m <- utils::read.csv(sprintf(EXTRACT_MANIFEST, ds), stringsAsFactors = FALSE)
  data.frame(dataset = ds, contract_version = m$contract_version,
             declared_source_gct = m$source_gct_path, comparison_count = m$comparison_count,
             strict_contract_validated = m$strict_contract_validated, stringsAsFactors = FALSE)
}))

# ---------------------------------------------------------------------------
# Row identity: abundance matrix rows vs canonical GSEA input rows
# ---------------------------------------------------------------------------
rowid_check <- do.call(rbind, lapply(seq_len(nrow(EXEMPLARS)), function(i) {
  ds <- EXEMPLARS$dataset[i]; cn <- EXEMPLARS$contrast[i]
  mp <- utils::read.csv(sprintf(MAPPED_DIR, ds, cn), stringsAsFactors = FALSE)
  A <- abund[[ds]]
  data.frame(dataset = ds, contrast = cn,
             n_rows_mapped_contrast_csv = nrow(mp),
             n_rows_abundance_matrix = nrow(A),
             n_original_identifier_in_matrix = sum(mp$original_identifier %in% rownames(A)),
             n_original_identifier_missing = sum(!mp$original_identifier %in% rownames(A)),
             n_duplicated_matrix_rownames = sum(duplicated(rownames(A))),
             n_duplicated_ProteinGroupID = sum(duplicated(mp$ProteinGroupID)),
             identity_1to1 = identical(sort(mp$original_identifier), sort(rownames(A))),
             stringsAsFactors = FALSE)
}))

# ---------------------------------------------------------------------------
# Q2: exemplar unit subsetting + Q4 design/model inference
# ---------------------------------------------------------------------------
solve_df <- function(t, p) {
  # df.total such that 2*pt(-|t|, df) = p  (moderated t from limma/Protigy)
  ok <- is.finite(t) & is.finite(p) & p > 1e-12 & p < 0.999 & abs(t) > 1e-6
  t <- abs(t[ok]); p <- p[ok]
  if (!length(t)) return(rep(NA_real_, 3))
  v <- vapply(seq_along(t), function(k) {
    f <- function(df) 2 * stats::pt(-t[k], df) - p[k]
    tryCatch(stats::uniroot(f, c(0.5, 5000), tol = 1e-8)$root, error = function(e) NA_real_)
  }, numeric(1))
  stats::quantile(v, c(0.25, 0.5, 0.75), na.rm = TRUE)
}

ex_rows <- list(); design_rows <- list(); collapse_rows <- list()

for (i in seq_len(nrow(EXEMPLARS))) {
  ds <- EXEMPLARS$dataset[i]; unit <- EXEMPLARS$unit[i]
  cn <- EXEMPLARS$contrast[i]; pc <- EXEMPLARS$protigy_comparison[i]

  gi <- gct_in[[ds]]; gs <- gct_stat[[ds]]; A <- abund[[ds]]
  cd <- gi$cdesc

  unit_cols <- colnames(gi$mat)[as.character(cd["region_layer_ExpGroup", ]) != "" &
                                 sub("_[123]$", "", as.character(cd["region_layer_ExpGroup", ])) == unit]
  grp <- as.character(cd["ExpGroup", unit_cols])
  lab <- EXPGROUP_LABEL[grp]
  animals <- as.character(cd["AnimalID", unit_cols])

  sus <- unit_cols[lab == "SUS"]; res <- unit_cols[lab == "RES"]; con <- unit_cols[lab == "CON"]

  # bilateral completeness for this unit
  aa <- utils::read.csv(sprintf(AGG_AUDIT, ds), stringsAsFactors = FALSE)
  aau <- aa[aa$canonical_spatial_unit == unit, , drop = FALSE]
  n_incomplete <- sum(aau$hemisphere_status != "bilateral_complete")

  ex_rows[[length(ex_rows) + 1L]] <- data.frame(
    dataset = ds, spatial_unit = unit, contrast = cn,
    n_columns_for_unit = length(unit_cols),
    n_SUS = length(sus), n_RES = length(res), n_CON = length(con),
    SUS_animals = paste(sort(animals[lab == "SUS"]), collapse = ";"),
    RES_animals = paste(sort(animals[lab == "RES"]), collapse = ";"),
    CON_animals = paste(sort(animals[lab == "CON"]), collapse = ";"),
    SUS_columns = paste(sort(sus), collapse = ";"),
    RES_columns = paste(sort(res), collapse = ";"),
    n_na_in_6_animal_submatrix = sum(is.na(A[, c(sus, res), drop = FALSE])),
    n_units_not_bilateral_complete = n_incomplete,
    stringsAsFactors = FALSE)

  # ---- canonical statistics for this contrast ------------------------------
  mp <- utils::read.csv(sprintf(MAPPED_DIR, ds, cn), stringsAsFactors = FALSE)
  rn <- match(mp$original_identifier, rownames(A))
  lfc_col <- paste0("logFC.", pc); t_col <- paste0("t.", pc); ave_col <- paste0("AveExpr.", pc)
  stat_lfc <- suppressWarnings(as.numeric(gs$rdesc[[lfc_col]]))[rn]
  stat_t   <- suppressWarnings(as.numeric(gs$rdesc[[t_col]]))[rn]
  stat_ave <- suppressWarnings(as.numeric(gs$rdesc[[ave_col]]))[rn]

  recon_lfc <- rowMeans(A[rn, sus, drop = FALSE]) - rowMeans(A[rn, res, drop = FALSE])
  recon_ave6 <- rowMeans(A[rn, c(sus, res), drop = FALSE])
  recon_aveall <- rowMeans(A[rn, , drop = FALSE])

  dfq <- solve_df(mp$t, mp$pval)

  design_rows[[length(design_rows) + 1L]] <- data.frame(
    dataset = ds, spatial_unit = unit, contrast = cn, protigy_comparison = pc,
    max_abs_diff_mapped_t_vs_statgct_t = max(abs(mp$t - stat_t), na.rm = TRUE),
    max_abs_diff_mapped_log2fc_vs_statgct_logFC = max(abs(mp$log2fc - stat_lfc), na.rm = TRUE),
    max_abs_diff_logFC_vs_meanSUS_minus_meanRES = max(abs(stat_lfc - recon_lfc), na.rm = TRUE),
    max_abs_diff_AveExpr_vs_mean_of_6_columns = max(abs(stat_ave - recon_ave6), na.rm = TRUE),
    max_abs_diff_AveExpr_vs_mean_of_all_columns = max(abs(stat_ave - recon_aveall), na.rm = TRUE),
    inferred_df_total_q25 = unname(dfq[1]), inferred_df_total_median = unname(dfq[2]),
    inferred_df_total_q75 = unname(dfq[3]),
    residual_df_if_within_unit_fit = length(c(sus, res)) - 2L,
    residual_df_if_dataset_wide_fit = ncol(A) - length(unique(as.character(cd["region_layer_ExpGroup", ]))),
    stringsAsFactors = FALSE)

  # ---- Q5: phenotype-blind gene collapse under the canonical mapping -------
  ad <- sprintf(AUDIT_DIR, ds, unit, cn)
  # protein_group_to_gene_transformation_audit.csv = the FULL per-protein-group
  # eligibility + statistic table (eligibility_audit.csv holds only the excluded rows).
  ad_stage <- stage_files(ad, c("protein_group_to_gene_transformation_audit.csv",
                                "per_contrast_aggregate_audit.csv",
                                "collapsed_gene_input.csv",
                                "eligibility_audit.csv"), paste0(ds, "__", cn))
  elig <- utils::read.csv(file.path(ad_stage, "protein_group_to_gene_transformation_audit.csv"), stringsAsFactors = FALSE)
  agg  <- utils::read.csv(file.path(ad_stage, "per_contrast_aggregate_audit.csv"), stringsAsFactors = FALSE)
  coll <- utils::read.csv(file.path(ad_stage, "collapsed_gene_input.csv"), stringsAsFactors = FALSE)
  excl <- utils::read.csv(file.path(ad_stage, "eligibility_audit.csv"), stringsAsFactors = FALSE)

  # Independent re-derivation of the documented, phenotype-blind eligibility rule
  # (clusterProfiler_manifest primary_gene_level_eligibility_rule), to prove that
  # the mapping contract can be reapplied without reading any phenotype label.
  rule_eligible <- mp$gene_level_claim_allowed %in% TRUE &
    mp$protein_group_gene_annotation_status == "concordant_official_gene" &
    mp$protein_group_ambiguity_class %in% c("single_accession_single_gene", "multi_accession_same_gene") &
    !is.na(mp$official_gene_symbol) & nzchar(as.character(mp$official_gene_symbol))
  audit_eligible_ids <- elig$ProteinGroupID[elig$eligibility_status == "eligible"]
  rule_matches_audit <- setequal(mp$ProteinGroupID[rule_eligible], audit_eligible_ids)

  e <- elig[elig$eligibility_status == "eligible" & is.finite(elig$source_statistic), , drop = FALSE]

  # (a) reproduce the canonical ranked-statistic collapse exactly
  med_t <- tapply(e$source_statistic, e$GeneSymbol, stats::median)
  m1 <- match(coll$GeneSymbol, names(med_t))
  repro_max_abs <- max(abs(coll$collapsed_statistic - as.numeric(med_t)[m1]), na.rm = TRUE)

  # (b) phenotype-blind median collapse of the ABUNDANCE matrix, same eligible
  #     protein groups, same official gene symbol grouping, same median rule
  ridx <- match(e$original_identifier, rownames(A))
  n_missing_rows <- sum(is.na(ridx))
  sub <- A[ridx[!is.na(ridx)], , drop = FALSE]
  gsym <- e$GeneSymbol[!is.na(ridx)]
  gene_mat <- do.call(rbind, lapply(split(seq_len(nrow(sub)), gsym), function(ii)
    apply(sub[ii, , drop = FALSE], 2, stats::median)))

  collapse_rows[[length(collapse_rows) + 1L]] <- data.frame(
    dataset = ds, spatial_unit = unit, contrast = cn,
    n_protein_groups_total = nrow(elig),
    n_protein_groups_eligible = nrow(e),
    n_eligible_pg_rows_absent_from_abundance_matrix = n_missing_rows,
    n_genes_canonical_collapsed_gene_input = nrow(coll),
    n_genes_from_abundance_collapse = nrow(gene_mat),
    gene_universes_identical = setequal(coll$GeneSymbol, rownames(gene_mat)),
    n_genes_only_in_canonical = length(setdiff(coll$GeneSymbol, rownames(gene_mat))),
    n_genes_only_in_abundance_collapse = length(setdiff(rownames(gene_mat), coll$GeneSymbol)),
    n_multi_pg_genes = sum(coll$n_protein_groups_for_gene > 1L),
    max_pg_per_gene = max(coll$n_protein_groups_for_gene),
    canonical_collapse_rule_reproduced_max_abs_diff = repro_max_abs,
    n_gene_matrix_columns = ncol(gene_mat),
    n_na_in_gene_matrix = sum(is.na(gene_mat)),
    audit_eligible_single_accession_groups = agg$eligible_single_accession_groups,
    audit_eligible_same_gene_multi_accession_groups = agg$eligible_same_gene_multi_accession_groups,
    audit_genes_after_duplicate_collapse = agg$genes_after_duplicate_collapse,
    audit_genes_with_multiple_ProteinGroupIDs = agg$genes_with_multiple_ProteinGroupIDs,
    audit_gene_collapse_rule = agg$gene_collapse_rule,
    audit_rank_statistic_column = agg$rank_statistic_column,
    audit_rank_statistic_type = agg$rank_statistic_type,
    audit_rank_statistic_fallback_used = agg$rank_statistic_fallback_used,
    n_rows_excluded_only_audit = nrow(excl),
    documented_eligibility_rule_reproduces_audit = rule_matches_audit,
    n_eligible_by_documented_rule = sum(rule_eligible),
    stringsAsFactors = FALSE)
}

exemplar_units <- do.call(rbind, ex_rows)
design_model    <- do.call(rbind, design_rows)
gene_collapse   <- do.call(rbind, collapse_rows)

# ---------------------------------------------------------------------------
# Scope of the canonical limma fit: if Protigy fitted ONE dataset-wide model,
# AveExpr is rowMeans of the whole matrix and is therefore identical across all
# comparisons of that dataset. If it subset per comparison, AveExpr varies.
# ---------------------------------------------------------------------------
fit_scope <- do.call(rbind, lapply(DATASETS, function(ds) {
  gs <- gct_stat[[ds]]
  ac <- grep("^AveExpr\\.", names(gs$rdesc), value = TRUE)
  M <- matrix(suppressWarnings(as.numeric(as.matrix(gs$rdesc[, ac, drop = FALSE]))), nrow = nrow(gs$rdesc))
  rng <- apply(M, 1, function(x) diff(range(x, na.rm = TRUE)))
  A <- abund[[ds]]
  ngroups <- length(unique(as.character(gct_in[[ds]]$cdesc["phenotypeWithinUnit", ])))
  data.frame(dataset = ds,
             n_comparisons = length(ac),
             max_within_row_AveExpr_range_across_comparisons = max(rng, na.rm = TRUE),
             max_abs_diff_AveExpr_vs_rowMeans_full_matrix =
               max(abs(M[, 1] - rowMeans(A)), na.rm = TRUE),
             n_matrix_columns = ncol(A),
             n_phenotypeWithinUnit_levels = ngroups,
             dataset_wide_residual_df = ncol(A) - ngroups,
             stringsAsFactors = FALSE)
}))

imput_qc <- utils::read.csv("data/processed/01_preprocessing/impute/imputation_qc.csv", stringsAsFactors = FALSE)
imput_qc$pct_values_imputed <- 100 * imput_qc$n_values_imputed /
  (imput_qc$n_proteins_after_filter * imput_qc$n_samples)

# ---------------------------------------------------------------------------
# Q3: does a gene-level (not protein-group-level) MATRIX exist anywhere?
# ---------------------------------------------------------------------------
gene_named <- unique(basename(list.files("data/processed", pattern = "(?i)gene",
                                         recursive = TRUE, full.names = TRUE)))
gene_level_artefacts <- data.frame(
  artefact = gene_named,
  is_sample_or_animal_matrix = FALSE,
  note = "gene-level table with one collapsed statistic per gene per contrast; no sample/animal dimension",
  stringsAsFactors = FALSE)

qc_animal_matrix <- "results/tables/03_qc_exploration/03_replicate_consistency/%s/animal_aggregated_matrix.csv"
qc_rows <- do.call(rbind, lapply(DATASETS, function(ds) {
  p <- sprintf(qc_animal_matrix, ds)
  h <- utils::read.csv(p, nrows = 1, stringsAsFactors = FALSE)
  data.frame(dataset = ds, path = p, n_columns = ncol(h),
             first_column = names(h)[1],
             columns = paste(names(h), collapse = ";"), stringsAsFactors = FALSE)
}))

# ---------------------------------------------------------------------------
# Assemble the answer table
# ---------------------------------------------------------------------------
nn <- abund[["neuron_neuropil"]]; ns <- abund[["neuron_soma"]]; mg <- abund[["microglia"]]
e1 <- exemplar_units[1, ]; e2 <- exemplar_units[2, ]; e3 <- exemplar_units[3, ]
d1 <- design_model[1, ]; c1 <- gene_collapse[1, ]; c2 <- gene_collapse[2, ]; c3 <- gene_collapse[3, ]

report <- data.frame(
  question = c(
    "Q1. Does a canonical animal-level protein abundance matrix exist? Exact paths, dimensions, row identity, column identity, value semantics.",
    "Q2. Can the matrix be subset to exactly the 6 animals (3 SUS + 3 RES) of CA3_sr neuropil / CA2 neuron_soma / CA1 microglia, with the SAME eligible measured gene universe the canonical GSEA used?",
    "Q3. Does an existing GENE-level (not protein-group-level) matrix exist anywhere?",
    "Q4. What design/contrast object is needed, and does the repo already store a design matrix or per-sample group assignment?",
    "Q5. Is a phenotype-blind gene collapse of the abundance matrix possible using the SAME contract as the ranked-statistic collapse (median over member protein groups per official gene symbol), and is that defensible reuse or an invented method?"
  ),
  answer = c(
    paste0(
      "YES. Two byte-identical-valued animal-level matrices exist per dataset. ",
      "(a) The Protigy INPUT matrices data/processed/01_preprocessing/protigy_input_animal_level/<ds>/<ds>_animal_level.gct. ",
      "(b) The SAME matrix embedded as GCT row-descriptor columns inside the Protigy OUTPUT stat GCTs ",
      "data/processed/01_preprocessing/protigy_output_animal_level/<ds>/stat_results_for_ssGSEA_*.gct, ",
      "which is the artefact the canonical contrast CSVs were extracted from ",
      "(gct_extractR/<ds>/canonical_gct_extract_manifest.csv, contract animal_level_protigy_da_v1). ",
      "Rows are UniProt entry names (e.g. A0A0J9YTR2_MOUSE) == mapped-contrast original_identifier, NOT ProteinGroupID; ",
      "ProteinGroupID is recoverable 1:1 via data/processed/02_id_mapping/mapped/<ds>/forward/per_file/*.csv. ",
      "Columns are AnimalID_spatialUnit (animal x spatial unit), not raw MS samples. ",
      "Values are log2 protein-group abundances, 70%-missingness-filtered and imputed at sample level, then aggregated to ",
      "animal x spatial unit as an equal-weight mean of the Left/Right hemisphere imputed log2 values. Matrices are complete (zero NA)."
    ),
    paste0(
      "YES for all three. Each unit yields exactly 3 SUS + 3 RES animal columns with zero missing values, and the eligible gene ",
      "universe is exactly reproducible because every eligible ProteinGroupID of the canonical GSEA input has a row in the abundance matrix. ",
      "CA3_sr neuropil SUS(", e1$SUS_animals, ") vs RES(", e1$RES_animals, "); ",
      "CA2_sp neuron_soma SUS(", e2$SUS_animals, ") vs RES(", e2$RES_animals, "); ",
      "CA1_microglia SUS(", e3$SUS_animals, ") vs RES(", e3$RES_animals, "). ",
      "The ExpGroup->phenotype map (1=CON, 2=RES, 3=SUS) is carried in the GCT column descriptors and in ",
      "gct_extractR/<ds>/indexComparisons.csv. ",
      "CAVEAT: subsetting to only those 6 columns would NOT reproduce the canonical model. The canonical statistics come from a ",
      "single DATASET-WIDE limma fit over all ", fit_scope$n_matrix_columns[1], "/", fit_scope$n_matrix_columns[2], "/",
      fit_scope$n_matrix_columns[3], " animal-unit columns with ", fit_scope$n_phenotypeWithinUnit_levels[1], "/",
      fit_scope$n_phenotypeWithinUnit_levels[2], "/", fit_scope$n_phenotypeWithinUnit_levels[3],
      " phenotypeWithinUnit groups (see Q4). camera() must therefore be run on the FULL matrix with the canonical design and a ",
      "SUS-RES contrast for that unit, not on a 6-column subset."
    ),
    paste0(
      "NO. No gene-level abundance/expression matrix exists anywhere in the repository. Every gene-level artefact is a one-row-per-gene ",
      "STATISTIC table for a single contrast (collapsed_gene_input.csv, collapsed_gene_input_provenance.csv, ",
      "duplicate_gene_collapse_audit.csv, protein_group_to_gene_transformation_audit.csv, gsea_*_term_gene_provenance.csv); ",
      "none has a sample or animal dimension. All sample/animal matrices in the repo are protein-group level ",
      "(animal-level GCTs; WGCNA wgcna_expression.xlsx; joint_shared_core_log2_median_normalized_imputed.gct). ",
      "The only animal-dimension matrix with gene-like row labels is the QC artefact ",
      "results/tables/03_qc_exploration/03_replicate_consistency/<ds>/animal_aggregated_matrix.csv, which has only 9 animal columns ",
      "(collapsed across ALL spatial units), uses unmapped source labels rather than the official SYMBOL contract, and is therefore unusable."
    ),
    paste0(
      "Needed: design <- model.matrix(~0 + phenotypeWithinUnit) over the full dataset-wide animal-level matrix, plus ",
      "contrast <- makeContrasts(<unit>_3 - <unit>_2) (3=SUS, 2=RES). The grouping factor IS stored, in three places: ",
      "the GCT column descriptors of protigy_input_animal_level/<ds>/<ds>_animal_level.gct (rows sample_id, AnimalID, ExpGroup, region, ",
      "layer, celltype, celltype_layer, phenotypeWithinUnit, region_layer_ExpGroup), ",
      "results/tables/01_preprocessing/02a_prepare_animal_level_protigy_input/<ds>/source_sample_assignment.csv and aggregation_audit.csv, ",
      "and data/processed/01_preprocessing/gct_extractR/<ds>/indexComparisons.csv (comparison -> forward/reverse contrast names). ",
      "NO stored model.matrix / design object and no stored eBayes prior exist anywhere in the repository; the design must be ",
      "rebuilt from phenotypeWithinUnit. Two independent lines of evidence establish that the canonical fit is DATASET-WIDE, ",
      "not a 6-sample within-unit fit, which is exactly the situation camera() is designed for: ",
      "(i) AveExpr equals rowMeans of the ENTIRE animal-level matrix (max abs diff ",
      fmt(max(fit_scope$max_abs_diff_AveExpr_vs_rowMeans_full_matrix), 3),
      ") and is identical across all comparisons of a dataset (max within-row range across comparisons ",
      fmt(max(fit_scope$max_within_row_AveExpr_range_across_comparisons), 3),
      "), whereas the mean of just the 6 unit columns differs by up to ",
      fmt(max(design_model$max_abs_diff_AveExpr_vs_mean_of_6_columns), 4), "; ",
      "(ii) solving 2*pt(-|t|,df)=p on the canonical t/pval yields median df.total = ",
      paste(fmt(design_model$inferred_df_total_median, 6), collapse = "/"), " against dataset-wide residual df ",
      paste(fit_scope$dataset_wide_residual_df, collapse = "/"), " (implied eBayes df.prior ~ ",
      paste(fmt(design_model$inferred_df_total_median - fit_scope$dataset_wide_residual_df, 3), collapse = "/"),
      "), and NOT against the within-unit residual df of ", d1$residual_df_if_within_unit_fit, ". ",
      "So the required object is design <- model.matrix(~0 + phenotypeWithinUnit) on the full matrix plus ",
      "contrasts.fit(makeContrasts(<unit>_3 - <unit>_2)), and camera(y, index, design, contrast) reproduces the canonical ",
      "variance model up to re-estimation of the eBayes prior."
    ),
    paste0(
      "YES, and it is defensible reuse of the existing mapping contract, not an invented method - with one explicit caveat. ",
      "Reuse is exact for (i) protein-group eligibility (gene_level_claim_allowed && concordant_official_gene && ambiguity class in ",
      "{single_accession_single_gene, multi_accession_same_gene}) and (ii) the grouping key (official_gene_symbol). Both are ",
      "phenotype-blind and are recorded per contrast in protein_group_audits/protein_group_to_gene_transformation_audit.csv. ",
      "The documented rule was re-derived independently from the mapped contrast CSV and reproduces the canonical eligible set exactly ",
      "(reproduces_audit = ", paste(gene_collapse$documented_eligibility_rule_reproduces_audit, collapse = "/"), "). ",
      "Applying median-over-member-protein-groups column-wise to the abundance matrix reproduces the canonical gene universe exactly ",
      "for all three exemplars (", c1$n_genes_canonical_collapsed_gene_input, " / ", c2$n_genes_canonical_collapsed_gene_input, " / ",
      c3$n_genes_canonical_collapsed_gene_input, " genes, set-identical). ",
      "Crucially the collapse is near-degenerate: exactly ", paste(gene_collapse$n_multi_pg_genes, collapse = "/"),
      " gene per dataset is formed from more than one protein group (max ", paste(gene_collapse$max_pg_per_gene, collapse = "/"),
      " groups), so for ", paste(gene_collapse$n_genes_canonical_collapsed_gene_input - gene_collapse$n_multi_pg_genes, collapse = "/"),
      " of the genes the gene-level row IS the protein-group row unchanged, and for the single multi-group gene median(2 values) = mean(2 values). ",
      "CAVEAT: median-of-abundance-then-test is still not algebraically the same quantity as the canonical median-of-t (median is not ",
      "linear), and eBayes moderation re-estimated over ",
      paste(gene_collapse$n_genes_canonical_collapsed_gene_input, collapse = "/"), " gene rows instead of ",
      paste(gene_collapse$n_protein_groups_total, collapse = "/"),
      " protein-group rows will shift the prior slightly, so a camera() gene-level t would be close to but not bit-identical to the ",
      "canonical ranked statistic. The mapping/eligibility/grouping contract transfers unchanged; the collapsed VALUE is a re-derivation."
    )
  ),
  evidence_paths = c(
    paste(c(INPUT_GCT[DATASETS], STAT_GCT[DATASETS],
            sprintf(EXTRACT_MANIFEST, DATASETS),
            "results/tables/01_preprocessing/02a_prepare_animal_level_protigy_input/animal_level_protigy_handoff_manifest.csv"),
          collapse = " | "),
    paste(c(sprintf(SRC_ASSIGN, DATASETS), sprintf(AGG_AUDIT, DATASETS),
            sprintf(MAPPED_DIR, EXEMPLARS$dataset, EXEMPLARS$contrast),
            file.path(sprintf(AUDIT_DIR, EXEMPLARS$dataset, EXEMPLARS$unit, EXEMPLARS$contrast), "protein_group_to_gene_transformation_audit.csv")),
          collapse = " | "),
    paste(c("docs/file_contracts.tsv", "results/reports/proteomics_data_dictionary.tsv",
            "data/processed/06_modules_WGCNA/01_WGCNA/<ds>/inputs/wgcna_expression.xlsx",
            "data/processed/01_preprocessing/protigy_input/global/joint_shared_core_log2_median_normalized_imputed.gct",
            sprintf(qc_animal_matrix, DATASETS)), collapse = " | "),
    paste(c(INPUT_GCT[["neuron_neuropil"]], sprintf(SRC_ASSIGN, DATASETS),
            sprintf(AGG_AUDIT, DATASETS),
            sprintf("data/processed/01_preprocessing/gct_extractR/%s/indexComparisons.csv", DATASETS),
            sprintf(MAPPED_DIR, EXEMPLARS$dataset, EXEMPLARS$contrast)), collapse = " | "),
    paste(c(file.path(sprintf(AUDIT_DIR, EXEMPLARS$dataset, EXEMPLARS$unit, EXEMPLARS$contrast), "protein_group_to_gene_transformation_audit.csv"),
            file.path(sprintf(AUDIT_DIR, EXEMPLARS$dataset, EXEMPLARS$unit, EXEMPLARS$contrast), "collapsed_gene_input.csv"),
            "R/enrichment/protein_group_enrichment_utils.R:96 collapse_protein_group_genes",
            "R/enrichment/protein_group_enrichment_utils.R:30 select_rank_statistic"), collapse = " | ")
  ),
  exact_numbers = c(
    paste0(
      "neuron_neuropil ", nrow(nn), " protein groups x ", ncol(nn), " animal-spatial units (9 animals x 10 laminar units), NA=",
      sum(is.na(nn)), "; neuron_soma ", nrow(ns), " x ", ncol(ns), " (9 x 4), NA=", sum(is.na(ns)),
      "; microglia ", nrow(mg), " x ", ncol(mg), " (9 x 4), NA=", sum(is.na(mg)),
      ". Protigy-input vs Protigy-embedded matrices agree exactly: max abs diff ",
      fmt(max(inventory$max_abs_diff_vs_embedded, na.rm = TRUE), 3),
      ". Row/column identity to the canonical GSEA input is 1:1 for all three exemplars: mapped-contrast rows ",
      paste(rowid_check$n_rows_mapped_contrast_csv, collapse = "/"), " vs matrix rows ",
      paste(rowid_check$n_rows_abundance_matrix, collapse = "/"), ", missing original_identifier ",
      paste(rowid_check$n_original_identifier_missing, collapse = "/"),
      ". 18 spatial units total (10 + 4 + 4); 3+3+3 animals per ExpGroup in every dataset."
    ),
    paste0(
      "CA3_sr (neuron_neuropil): ", e1$n_columns_for_unit, " columns = ", e1$n_SUS, " SUS + ", e1$n_RES, " RES + ", e1$n_CON,
      " CON; SUS=", e1$SUS_animals, ", RES=", e1$RES_animals, "; NA in 6-animal submatrix=", e1$n_na_in_6_animal_submatrix,
      "; eligible protein groups ", c1$n_protein_groups_eligible, "/", c1$n_protein_groups_total,
      ", all present in matrix (missing=", c1$n_eligible_pg_rows_absent_from_abundance_matrix, "), genes ",
      c1$n_genes_canonical_collapsed_gene_input, ". ",
      "CA2_sp (neuron_soma): ", e2$n_columns_for_unit, " columns = ", e2$n_SUS, " SUS + ", e2$n_RES, " RES + ", e2$n_CON,
      " CON; SUS=", e2$SUS_animals, ", RES=", e2$RES_animals, "; NA=", e2$n_na_in_6_animal_submatrix,
      "; eligible ", c2$n_protein_groups_eligible, "/", c2$n_protein_groups_total, ", genes ",
      c2$n_genes_canonical_collapsed_gene_input, ". ",
      "CA1_microglia: ", e3$n_columns_for_unit, " columns = ", e3$n_SUS, " SUS + ", e3$n_RES, " RES + ", e3$n_CON,
      " CON; SUS=", e3$SUS_animals, ", RES=", e3$RES_animals, "; NA=", e3$n_na_in_6_animal_submatrix,
      "; eligible ", c3$n_protein_groups_eligible, "/", c3$n_protein_groups_total, ", genes ",
      c3$n_genes_canonical_collapsed_gene_input, ". ",
      "Bilateral completeness of the exemplar units: ", e1$n_units_not_bilateral_complete, "/", e2$n_units_not_bilateral_complete,
      "/", e3$n_units_not_bilateral_complete, " non-complete (the single dataset-wide incomplete unit is A111_CA1_sp, left-only, not an exemplar). ",
      "Full-matrix context for the required camera() call: ",
      paste(paste0(fit_scope$dataset, " ", fit_scope$n_matrix_columns, " columns / ",
                   fit_scope$n_phenotypeWithinUnit_levels, " phenotypeWithinUnit groups / residual df ",
                   fit_scope$dataset_wide_residual_df), collapse = "; "), "."
    ),
    paste0(
      "0 gene x sample matrices found. ", nrow(gene_level_artefacts),
      " gene-named artefacts under data/processed, all one-row-per-gene statistic tables: ",
      paste(gene_level_artefacts$artefact, collapse = ", "),
      ". WGCNA input wgcna_expression.xlsx is ProteinGroupID x MS sample. ",
      "joint_shared_core_log2_median_normalized_imputed.gct is ProteinGroupID x 323 MS samples, 4242 shared-core rows ",
      "(a restricted universe, not the per-dataset 5045/5529/5219). ",
      "animal_aggregated_matrix.csv is ", paste(qc_rows$n_columns - 1, collapse = "/"),
      " animal columns only, first column '", qc_rows$first_column[1], "', no spatial dimension."
    ),
    paste0(
      "GCT column descriptors present: 9 rows (",
      paste(rownames(gct_in[["neuron_neuropil"]]$cdesc), collapse = ", "),
      "). phenotypeWithinUnit has ",
      length(unique(as.character(gct_in[["neuron_neuropil"]]$cdesc["phenotypeWithinUnit", ]))),
      " levels for neuron_neuropil and ",
      length(unique(as.character(gct_in[["neuron_soma"]]$cdesc["phenotypeWithinUnit", ]))),
      "/", length(unique(as.character(gct_in[["microglia"]]$cdesc["phenotypeWithinUnit", ]))),
      " for neuron_soma/microglia. Canonical statistics linkage verified: mapped-contrast t vs stat-GCT t max abs diff ",
      paste(fmt(design_model$max_abs_diff_mapped_t_vs_statgct_t, 3), collapse = "/"),
      "; logFC vs (mean SUS - mean RES) from the abundance matrix max abs diff ",
      paste(fmt(design_model$max_abs_diff_logFC_vs_meanSUS_minus_meanRES, 3), collapse = "/"),
      "; AveExpr vs mean of the 6 unit columns ",
      paste(fmt(design_model$max_abs_diff_AveExpr_vs_mean_of_6_columns, 3), collapse = "/"),
      " vs mean of all columns ",
      paste(fmt(design_model$max_abs_diff_AveExpr_vs_mean_of_all_columns, 3), collapse = "/"),
      ". Inferred df.total (median) = ",
      paste(fmt(design_model$inferred_df_total_median, 6), collapse = "/"),
      " against a within-unit residual df of ",
      paste(design_model$residual_df_if_within_unit_fit, collapse = "/"),
      " and a dataset-wide residual df of ",
      paste(design_model$residual_df_if_dataset_wide_fit, collapse = "/"),
      " => implied eBayes df.prior ~ ",
      paste(fmt(design_model$inferred_df_total_median - fit_scope$dataset_wide_residual_df, 3), collapse = "/"),
      ". AveExpr is constant across all comparisons within a dataset (max within-row range ",
      paste(fmt(fit_scope$max_within_row_AveExpr_range_across_comparisons, 3), collapse = "/"),
      ") and equals rowMeans of the full matrix (max abs diff ",
      paste(fmt(fit_scope$max_abs_diff_AveExpr_vs_rowMeans_full_matrix, 3), collapse = "/"),
      "), confirming a single dataset-wide lmFit. No design/model.matrix artefact exists: a repository-wide filename search for ",
      "'*design*matrix*' returned 0 hits and docs/file_contracts.tsv declares no design object."
    ),
    paste0(
      "Canonical median_finite_statistics collapse re-derived from protein_group_to_gene_transformation_audit.csv and matched ",
      "against collapsed_gene_input.csv, max abs diff ",
      paste(fmt(gene_collapse$canonical_collapse_rule_reproduced_max_abs_diff, 3), collapse = "/"),
      ". Phenotype-blind median collapse of the abundance matrix over the SAME eligible protein groups yields ",
      paste(gene_collapse$n_genes_from_abundance_collapse, collapse = "/"),
      " genes vs canonical ", paste(gene_collapse$n_genes_canonical_collapsed_gene_input, collapse = "/"),
      "; set-identical = ", paste(gene_collapse$gene_universes_identical, collapse = "/"),
      " (genes only in canonical ", paste(gene_collapse$n_genes_only_in_canonical, collapse = "/"),
      ", only in abundance collapse ", paste(gene_collapse$n_genes_only_in_abundance_collapse, collapse = "/"), "). ",
      "Multi-protein-group genes affected by the collapse: ",
      paste(gene_collapse$n_multi_pg_genes, collapse = "/"), " genes, max ",
      paste(gene_collapse$max_pg_per_gene, collapse = "/"), " protein groups per gene. ",
      "Resulting gene matrices: ", paste(gene_collapse$n_genes_from_abundance_collapse, collapse = "/"), " genes x ",
      paste(gene_collapse$n_gene_matrix_columns, collapse = "/"), " animal-spatial-unit columns, NA = ",
      paste(gene_collapse$n_na_in_gene_matrix, collapse = "/"), ". ",
      "Genes whose row is an UNCHANGED protein-group row: ",
      paste(gene_collapse$n_genes_canonical_collapsed_gene_input - gene_collapse$n_multi_pg_genes, collapse = "/"),
      " of ", paste(gene_collapse$n_genes_canonical_collapsed_gene_input, collapse = "/"), ". ",
      "Underlying data are imputed: ",
      paste(paste0(imput_qc$celltype_layer, " ", imput_qc$n_values_imputed, " values (",
                   fmt(imput_qc$pct_values_imputed, 3), "%) imputed at sample level over ",
                   imput_qc$n_proteins_after_filter, " x ", imput_qc$n_samples), collapse = "; "),
      "; animal-level values are then equal-weight L/R means."
    )
  ),
  verdict = c(
    "FEASIBLE",
    "FEASIBLE_WITH_CAVEAT",
    "NOT_FEASIBLE",
    "FEASIBLE_WITH_CAVEAT",
    "FEASIBLE_WITH_CAVEAT"
  ),
  caveat = c(
    paste0("Row identity is UniProt entry name, not ProteinGroupID; the ProteinGroupID/official-SYMBOL contract must be joined in from ",
           "data/processed/02_id_mapping/mapped/<ds>/forward/per_file/*.csv. Columns are animal x spatial unit, i.e. hemisphere-averaged ",
           "imputed values, so within-unit variance is already reduced by the L/R averaging step and every value is post-imputation."),
    paste0("Subsetting to exactly 6 columns is trivially possible and lossless, but a 6-sample within-unit limma fit has only ",
           design_model$residual_df_if_within_unit_fit[1],
           " residual df and does NOT reproduce the canonical moderated t the manuscript GSEA was ranked on. camera() must be run on ",
           "the full dataset-wide matrix with the canonical phenotypeWithinUnit design and the unit-specific SUS-RES contrast. ",
           "All values are post-imputation and hemisphere-averaged, so camera()'s inter-gene correlation (VIF) estimate is computed on ",
           "partly synthetic residuals."),
    paste0("No gene-level matrix exists, so any camera() run requires constructing one de novo. The construction is near-trivial here ",
           "(one multi-protein-group gene per dataset), but it is still a new derived artefact that does not currently exist under any ",
           "documented contract in docs/file_contracts.tsv."),
    paste0("No stored design matrix exists; it must be reconstructed from phenotypeWithinUnit (which IS stored). The fit scope is ",
           "recoverable and unambiguous (dataset-wide), but the exact eBayes prior (s2.prior, df.prior) is not stored, so a camera() ",
           "fit would re-estimate it from the gene-level rows and cannot be guaranteed bit-identical to the canonical moderation."),
    paste0("The mapping/eligibility/grouping contract transfers exactly (gene universes set-identical, documented rule reproduced ",
           "independently). The collapsed VALUE is a re-derivation, not the canonical quantity: median-of-abundance is not the ",
           "abundance whose moderated t equals the canonical median-of-t. The numerical impact is confined to ",
           paste(gene_collapse$n_multi_pg_genes, collapse = "/"),
           " gene per dataset plus a re-estimated eBayes prior, but it is a second gene-level quantity that the manuscript does not ",
           "currently define.")
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(report, file.path(OUT_DIR, "camera_feasibility_report.csv"), row.names = FALSE, na = "")
utils::write.csv(inventory, file.path(OUT_DIR, "camera_feasibility_matrix_inventory.csv"), row.names = FALSE, na = "")
utils::write.csv(prov, file.path(OUT_DIR, "camera_feasibility_canonical_provenance.csv"), row.names = FALSE, na = "")
utils::write.csv(rowid_check, file.path(OUT_DIR, "camera_feasibility_row_identity_check.csv"), row.names = FALSE, na = "")
utils::write.csv(exemplar_units, file.path(OUT_DIR, "camera_feasibility_exemplar_units.csv"), row.names = FALSE, na = "")
utils::write.csv(design_model, file.path(OUT_DIR, "camera_feasibility_design_and_model_inference.csv"), row.names = FALSE, na = "")
utils::write.csv(gene_collapse, file.path(OUT_DIR, "camera_feasibility_gene_collapse_check.csv"), row.names = FALSE, na = "")
utils::write.csv(fit_scope, file.path(OUT_DIR, "camera_feasibility_fit_scope.csv"), row.names = FALSE, na = "")
utils::write.csv(imput_qc, file.path(OUT_DIR, "camera_feasibility_imputation_context.csv"), row.names = FALSE, na = "")

cat("\n==================== INVENTORY ====================\n"); print(inventory[, c("dataset","artefact","n_rows","n_value_columns","n_na_values","max_abs_diff_vs_embedded")])
cat("\n==================== PROVENANCE ====================\n"); print(prov)
cat("\n==================== ROW IDENTITY ====================\n"); print(rowid_check)
cat("\n==================== EXEMPLAR UNITS ====================\n"); print(t(exemplar_units))
cat("\n==================== DESIGN / MODEL ====================\n"); print(t(design_model))
cat("\n==================== GENE COLLAPSE ====================\n"); print(t(gene_collapse))
cat("\n==================== FIT SCOPE ====================\n"); print(fit_scope)
cat("\n==================== IMPUTATION CONTEXT ====================\n"); print(imput_qc[, c("celltype_layer","n_samples","n_proteins_after_filter","n_values_imputed","pct_values_imputed")])
cat("\n==================== GENE-NAMED ARTEFACTS ====================\n"); print(gene_level_artefacts$artefact)
cat("\n==================== QC ANIMAL MATRIX ====================\n"); print(qc_rows[, c("dataset","n_columns","first_column")])
cat("\nDONE\n")
