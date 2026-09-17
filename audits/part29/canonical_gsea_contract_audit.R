# ============================================================================
# Part-29 sections 1, 2, 6 -- canonical enrichment contract reconstruction,
# rank-statistic audit, gene-set-size-bound audit.
#
# READ-ONLY AUDIT. No canonical output is modified, rerun or regenerated.
# Writes only under results/tables/publication_audits/upstream_enrichment_v10/.
# ============================================================================

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

OUT_DIR <- "results/tables/publication_audits/upstream_enrichment_v10"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DATASETS <- c("microglia", "neuron_neuropil", "neuron_soma")
THEME_TABLE <- file.path(
  "results/tables/10_biological_integration/gsea_wgcna_concordance/global",
  "ontology_aware_gsea_theme_assignments_all_contrasts.csv"
)

# ---------------------------------------------------------------------------
# Windows MAX_PATH workaround.
# Several microglia per-contrast audit CSVs sit at absolute paths > 259 chars
# (e.g. .../microglia/phenotype_within_unit/CA1_microglia/
#        CA1microgliares_CA1microgliacon/protein_group_audits/
#        per_contrast_aggregate_audit.csv  = 270 chars) and cannot be opened by
# the Win32 file API from this repo root. The pipeline's own manifests record
# this repo root as "P://...", i.e. the canonical run mounted the repo on a
# substituted short drive. Reuse that convention: prefer a short root when one
# is available, create it if it is not, and release it again on exit.
# ---------------------------------------------------------------------------
REPO_ROOT <- getwd()
SHORT_ROOT <- NULL
SHORT_ROOT_CREATED <- FALSE
if (.Platform$OS.type == "windows") {
  if (dir.exists("P:/data/processed/04_differential_expression_enrichment")) {
    SHORT_ROOT <- "P:/"
  } else if (!dir.exists("P:/")) {
    invisible(suppressWarnings(system2("cmd", c("/c", "subst", "P:",
      shQuote(gsub("/", "\\\\", REPO_ROOT))), stdout = NULL, stderr = NULL)))
    if (dir.exists("P:/data/processed/04_differential_expression_enrichment")) {
      SHORT_ROOT <- "P:/"
      SHORT_ROOT_CREATED <- TRUE
    }
  }
  if (isTRUE(SHORT_ROOT_CREATED)) {
    reg.finalizer(environment(), function(e) {
      suppressWarnings(system2("cmd", c("/c", "subst", "P:", "/D"),
                               stdout = NULL, stderr = NULL))
    }, onexit = TRUE)
  }
}
rp <- function(relpath) {
  if (!is.null(SHORT_ROOT)) {
    cand <- paste0(SHORT_ROOT, relpath)
    if (file.exists(cand)) return(cand)
  }
  if (!file.exists(relpath)) {
    stop("Cannot open '", relpath, "' (absolute length ",
         nchar(file.path(REPO_ROOT, relpath)), " chars). Windows MAX_PATH is 259 and ",
         "no short repo root is mounted. Run:  subst P: \"",
         gsub("/", "\\\\", REPO_ROOT), "\"  and re-run this script.", call. = FALSE)
  }
  relpath
}

rd <- function(p) utils::read.csv(rp(p), stringsAsFactors = FALSE, check.names = FALSE)
one <- function(x) {
  u <- unique(as.character(x))
  if (length(u) == 1L) u else paste(u, collapse = "|")
}

# ---------------------------------------------------------------------------
# 0. Software versions (renv.lock pin vs installed)
# ---------------------------------------------------------------------------
lock <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)
pkgs <- c("clusterProfiler", "DOSE", "fgsea", "org.Mm.eg.db", "GO.db",
          "AnnotationDbi", "limma", "withr")
lock_ver <- vapply(pkgs, function(p) {
  v <- lock$Packages[[p]]$Version
  if (is.null(v)) NA_character_ else as.character(v)
}, character(1))
inst_ver <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1))
ver_match <- identical(unname(lock_ver), unname(inst_ver))
SOFTWARE_VERSIONS <- paste0(
  "R ", lock$R$Version, " (installed ", paste(R.version$major, R.version$minor, sep = "."), "); ",
  paste(paste0(pkgs, " ", lock_ver), collapse = "; "),
  "; renv.lock_pin_equals_installed=", ver_match
)

# ---------------------------------------------------------------------------
# 1. Manuscript key set: every (dataset, contrast, spatial_unit) in the theme table
# ---------------------------------------------------------------------------
theme <- rd(THEME_TABLE)
keys <- unique(theme[, c("dataset", "contrast", "phenotype_contrast",
                         "spatial_unit", "source_comparison")])
keys <- keys[order(keys$dataset, keys$contrast, keys$spatial_unit), ]
stopifnot(nrow(keys) == 54L)
EVIDENCE_FAMILY <- one(theme$evidence_source_family)
GO_DB_VERSION <- one(theme$GO_db_package_version)
GO_SOURCE_DATE <- one(theme$GO_source_date)

# ---------------------------------------------------------------------------
# 2. clusterProfiler manifests -- GSEA_GO rows only
# ---------------------------------------------------------------------------
man <- do.call(rbind, lapply(DATASETS, function(d) {
  m <- rd(file.path("data/processed/04_differential_expression_enrichment/clusterProfiler",
                    d, "clusterProfiler_manifest.csv"))
  m[m$result_type == "GSEA_GO", , drop = FALSE]
}))
stopifnot(nrow(man) == 54L)

# ---------------------------------------------------------------------------
# 3. ProTigy comparison index -> effect definition; ExpGroup -> phenotype
# ---------------------------------------------------------------------------
idx <- do.call(rbind, lapply(DATASETS, function(d) {
  i <- rd(file.path("data/processed/01_preprocessing/gct_extractR_animal_level", d,
                    "indexComparisons.csv"))
  i$dataset <- d
  i
}))
grp_tok <- do.call(rbind, lapply(seq_len(nrow(idx)), function(i) {
  lab <- idx$comparison[i]                       # e.g. CA1_sp_2.over.CA1_sp_1
  parts <- strsplit(lab, ".over.", fixed = TRUE)[[1]]
  num <- sub(".*_([0-9]+)$", "\\1", parts)
  fwd <- idx$parsed_forward_comparison[i]        # e.g. CA1spres_CA1spcon
  toks <- strsplit(fwd, "_", fixed = TRUE)[[1]]
  tok <- sub("^.*(con|res|sus)$", "\\1", toks)
  data.frame(ExpGroup = num, phenotype = tok, stringsAsFactors = FALSE)
}))
GROUP_MAP <- unique(grp_tok)
GROUP_MAP <- GROUP_MAP[order(GROUP_MAP$ExpGroup), ]
stopifnot(nrow(GROUP_MAP) == 3L)
GROUP_MAP$phenotype <- toupper(GROUP_MAP$phenotype)
grp2pheno <- stats::setNames(GROUP_MAP$phenotype, GROUP_MAP$ExpGroup)

# ---------------------------------------------------------------------------
# 4. Animal-level design (biological replicates per dataset / unit / group)
# ---------------------------------------------------------------------------
assign_list <- lapply(DATASETS, function(d) {
  a <- rd(file.path("results/tables/01_preprocessing/02a_prepare_animal_level_protigy_input",
                    d, "source_sample_assignment.csv"))
  a[a$inclusion_status == "included", , drop = FALSE]
})
names(assign_list) <- DATASETS

map_route_to_canonical_unit <- function(dataset, route_unit) {
  units <- unique(assign_list[[dataset]]$canonical_spatial_unit)
  if (route_unit %in% units) return(route_unit)
  stripped <- sub("_[^_]+$", "", route_unit)
  if (stripped %in% units) return(stripped)
  stop("Cannot map route_unit '", route_unit, "' for dataset ", dataset)
}

animals_for <- function(dataset, route_unit, groups) {
  a <- assign_list[[dataset]]
  cu <- map_route_to_canonical_unit(dataset, route_unit)
  a <- a[a$canonical_spatial_unit == cu & as.character(a$ExpGroup) %in% groups, , drop = FALSE]
  per <- tapply(a$AnimalID, as.character(a$ExpGroup), function(x) length(unique(x)))
  list(n_unique = length(unique(a$AnimalID)),
       per_group = paste(paste0(grp2pheno[names(per)], "=", as.integer(per)), collapse = ";"),
       canonical_unit = cu)
}

# ---------------------------------------------------------------------------
# 5. Independent reimplementation of the repo's rank/collapse contract
# ---------------------------------------------------------------------------
norm_names <- function(nm) tolower(gsub("[^a-z0-9]", "", tolower(nm)))
find_col <- function(df, candidates) {
  normalized <- norm_names(names(df))
  wanted <- norm_names(candidates)
  hit <- match(wanted, normalized)
  hit <- hit[!is.na(hit)]
  if (length(hit)) names(df)[hit[[1]]] else NA_character_
}
select_rank_statistic_local <- function(df) {
  moderated <- find_col(df, c("t", "statistic", "stat", "moderated_t", "t_statistic"))
  if (!is.na(moderated)) {
    return(list(column = moderated, type = "moderated_or_signed_inferential",
                fallback_used = FALSE))
  }
  fallback <- find_col(df, c("log2fc", "logfc", "log2foldchange", "avg_log2fc", "avg_logfc"))
  if (is.na(fallback)) stop("no supported signed rank statistic")
  list(column = fallback, type = "log_fold_change", fallback_used = TRUE)
}
ALLOWED_CLASSES <- c("single_accession_single_gene", "multi_accession_same_gene")

derive_seed_local <- function(base_seed, comparison, analysis_type) {
  modulus <- 2147483646
  key <- paste0(nchar(comparison, type = "bytes"), ":", comparison, "|",
                nchar(analysis_type, type = "bytes"), ":", analysis_type)
  h <- as.double(base_seed) %% modulus
  for (code in utf8ToInt(enc2utf8(key))) h <- (h * 131 + as.double(code)) %% modulus
  s <- as.integer(h)
  if (s < 1L) s <- 1L
  s
}

# ---------------------------------------------------------------------------
# 6. GO BP gene sets exactly as clusterProfiler::gseGO(keyType="SYMBOL") builds them
# ---------------------------------------------------------------------------
GO_DATA <- clusterProfiler:::get_GO_data(org.Mm.eg.db, "BP", "SYMBOL")
GO_SETS <- DOSE:::getGeneSet(GO_DATA)
N_GO_BP_SETS <- length(GO_SETS)
go_term_name <- function(ids) {
  nm <- suppressWarnings(DOSE:::TERM2NAME(ids, GO_DATA))
  as.character(nm)
}
set_sizes_for <- function(gene_names) {
  vapply(GO_SETS, function(p) length(unique(stats::na.omit(match(p, gene_names)))), integer(1))
}

MIN_GS <- 10L
MAX_GS <- 800L

cfg_lines <- readLines("config/clusterProfiler_config.yml", warn = FALSE)
cfg_get <- function(k) {
  ln <- grep(paste0("^\\s*", k, ":"), cfg_lines, value = TRUE)[1]
  trimws(sub(paste0("^\\s*", k, ":\\s*"), "", ln))
}
CFG_MIN <- as.integer(cfg_get("min_gs_size"))
CFG_MAX <- as.integer(cfg_get("max_gs_size"))
CFG_PCUT <- as.numeric(cfg_get("pvalue_cutoff"))
CFG_QCUT <- as.numeric(cfg_get("qvalue_cutoff"))
CFG_PADJ <- cfg_get("p_adjust_method")
CFG_SEEDBASE <- as.integer(cfg_get("gsea_seed_base"))
CFG_HASH <- unname(tools::md5sum("config/clusterProfiler_config.yml"))
stopifnot(CFG_MIN == MIN_GS, CFG_MAX == MAX_GS)

GSEGO_FORMALS <- formals(clusterProfiler::gseGO)
GSEA_EXPONENT <- eval(GSEGO_FORMALS$exponent)
GSEA_EPS <- eval(GSEGO_FORMALS$eps)
GSEGO_DEFAULT_MIN <- eval(GSEGO_FORMALS$minGSSize)
GSEGO_DEFAULT_MAX <- eval(GSEGO_FORMALS$maxGSSize)

# ---------------------------------------------------------------------------
# 7. Per-result contract reconstruction
# ---------------------------------------------------------------------------
compare_go <- lapply(DATASETS, function(d) {
  rd(file.path("results/tables/04_differential_expression_enrichment/compareGO", d,
               "BP/phenotype_within_unit/all_route_units/compareGO_term_comparison.csv"))
})
names(compare_go) <- DATASETS

universe_cache <- list()
rows <- list()
rank_rows <- list()
dropped_records <- list()

for (i in seq_len(nrow(man))) {
  mr <- man[i, ]
  ds <- mr$dataset
  cmp <- mr$comparison
  ru <- mr$route_unit
  kk <- keys[keys$dataset == ds & keys$source_comparison == cmp, ]
  stopifnot(nrow(kk) == 1L)

  audit_root <- file.path("data/processed/04_differential_expression_enrichment/clusterProfiler",
                          ds, mr$route_category, ru, cmp, "protein_group_audits")
  agg <- rd(file.path(audit_root, "per_contrast_aggregate_audit.csv"))
  diag <- rd(file.path(audit_root, "gsea_input_diagnostics.csv"))
  repro <- rd(file.path(audit_root, "gsea_reproducibility.csv"))
  repro_go <- repro[repro$analysis_type == "gseGO_BP", ]
  coll <- rd(file.path(audit_root, "duplicate_gene_collapse_audit.csv"))
  cgi <- rd(file.path(audit_root, "collapsed_gene_input.csv"))
  qc <- rd(file.path("data/processed/04_differential_expression_enrichment/clusterProfiler",
                     ds, mr$route_category, ru, cmp, "QC_summary.csv"))

  # ---- independent re-derivation from the mapped DA table ----
  src <- file.path("data/processed/02_id_mapping/mapped", ds, "forward/per_file",
                   paste0(cmp, ".csv"))
  dfin <- rd(src)
  stat <- select_rank_statistic_local(dfin)
  has_t <- !is.na(find_col(dfin, "t"))
  has_log2fc <- !is.na(find_col(dfin, "log2fc"))
  vals <- suppressWarnings(as.numeric(dfin[[stat$column]]))
  elig <- dfin$gene_level_claim_allowed %in% c(TRUE, "TRUE") &
    dfin$protein_group_ambiguity_class %in% ALLOWED_CLASSES &
    dfin$protein_group_gene_annotation_status == "concordant_official_gene" &
    !is.na(dfin$official_gene_symbol) & nzchar(as.character(dfin$official_gene_symbol))
  use <- elig & is.finite(vals)
  sp <- split(vals[use], as.character(dfin$official_gene_symbol)[use])
  re_med <- vapply(sp, function(x) stats::median(x[is.finite(x)]), numeric(1))
  re_n <- vapply(sp, length, integer(1))
  re_disc <- vapply(sp, function(x) any(sign(x) > 0) && any(sign(x) < 0), logical(1))

  ref <- stats::setNames(cgi$collapsed_statistic, cgi$GeneSymbol)
  shared <- intersect(names(ref), names(re_med))
  collapse_reproduced <- identical(sort(names(ref)), sort(names(re_med))) &&
    max(abs(ref[shared] - re_med[shared])) < 1e-9
  n_multi_re <- sum(re_n > 1L)
  n_disc_re <- sum(re_disc)

  # ---- gene universe -> GO set sizes (cached per dataset, verified identical) ----
  gene_names <- names(sort(ref, decreasing = TRUE))
  if (is.null(universe_cache[[ds]])) {
    universe_cache[[ds]] <- list(genes = sort(gene_names),
                                 sizes = set_sizes_for(gene_names),
                                 ref_comparison = cmp)
  }
  universe_identical <- identical(universe_cache[[ds]]$genes, sort(gene_names))
  sizes <- if (universe_identical) universe_cache[[ds]]$sizes else set_sizes_for(gene_names)
  n_tested <- sum(sizes >= MIN_GS & sizes <= MAX_GS)

  cg <- compare_go[[ds]]
  cg_i <- cg[cg$comparison == cmp, , drop = FALSE]
  n_reported <- nrow(cg_i)
  stopifnot(as.integer(mr$n_terms) == n_reported)

  # Gene sets that passed the size filter but carry no row in the canonical
  # output: DOSE:::GSEA_fgsea drops rows with NA p-value after BH adjustment.
  tested_ids <- names(sizes)[sizes >= MIN_GS & sizes <= MAX_GS]
  missing_ids <- setdiff(tested_ids, cg_i$ID)
  if (length(missing_ids)) {
    dropped_records[[length(dropped_records) + 1L]] <- data.frame(
      dataset = ds, comparison = kk$contrast, spatial_unit = kk$spatial_unit,
      source_comparison = cmp, GO_ID = missing_ids,
      GO_description = go_term_name(missing_ids),
      setSize_in_universe = as.integer(sizes[missing_ids]),
      stringsAsFactors = FALSE)
  }

  seed_file <- as.integer(repro_go$gsea_seed)
  seed_recomputed <- derive_seed_local(CFG_SEEDBASE, cmp, "gseGO_BP")

  arms <- strsplit(kk$contrast, " - ")[[1]]
  ix <- idx[idx$dataset == ds & idx$parsed_forward_comparison == cmp, ]
  effect_def <- paste0(
    kk$contrast, "; ProTigy contrast ", ix$comparison[1],
    "; log2FC = mean(log2 ", arms[1], ") - mean(log2 ", arms[2],
    "); positive t/log2FC = higher in ", arms[1],
    "; ExpGroup coding 1=", grp2pheno[["1"]], ", 2=", grp2pheno[["2"]], ", 3=", grp2pheno[["3"]])
  an <- animals_for(ds, ru, names(grp2pheno)[grp2pheno %in% arms])

  rows[[length(rows) + 1L]] <- data.frame(
    dataset = ds,
    comparison = kk$contrast,
    spatial_unit = kk$spatial_unit,
    n_unique_animals = an$n_unique,
    animals_per_group = an$per_group,
    input_matrix_or_table = src,
    DA_model = paste0(
      "limma moderated t-test computed by ProTigy on animal-level log2 imputed protein-group ",
      "intensities (Left/Right hemispheres averaged within AnimalID x canonical spatial unit; ",
      "analysis/01_preprocessing/02a_prepare_animal_level_protigy_input.r); two-group comparison fitted ",
      "within one spatial unit; topTable fields logFC/AveExpr/t/P.Value/adj.P.Val/B split out by ",
      "analysis/01_preprocessing/03_gct_extractR.r and carried into 02_id_mapping as log2fc/aveExpr/t/pval/padj/B"),
    effect_definition = effect_def,
    rank_statistic_column = mr$rank_statistic_column,
    rank_statistic_type = mr$rank_statistic_type,
    rank_statistic_fallback_used = mr$rank_statistic_fallback_used,
    n_ranked_protein_groups = sum(coll$n_protein_groups_for_gene),
    n_ranked_unique_genes = nrow(cgi),
    duplicate_gene_collapse_rule = mr$duplicate_gene_collapse_rule,
    n_genes_with_multiple_protein_groups = agg$genes_with_multiple_ProteinGroupIDs,
    n_genes_with_discordant_protein_group_directions = agg$genes_with_discordant_directions,
    GO_ontology = paste0("GO:", mr$ontology),
    GO_annotation_version = paste0(
      "org.Mm.eg.db ", utils::packageVersion("org.Mm.eg.db"),
      "; GO.db ", GO_DB_VERSION, "; GO source date ", GO_SOURCE_DATE,
      "; stage02 orgdb_package_version ", one(dfin$orgdb_package_version)),
    minGSSize = MIN_GS,
    maxGSSize = MAX_GS,
    GSEA_backend = paste0(
      "clusterProfiler::gseGO(ont=BP, keyType=SYMBOL, OrgDb=org.Mm.eg.db, by='",
      repro_go$clusterprofiler_by, "') -> DOSE:::GSEA_internal -> DOSE:::GSEA_fgsea -> ",
      "fgsea::fgsea/fgseaMultilevel; clusterProfiler seed argument = ",
      repro_go$clusterprofiler_seed_argument),
    GSEA_exponent = GSEA_EXPONENT,
    eps = GSEA_EPS,
    nPermSimple = as.integer(repro_go$n_perm_simple),
    seed_base = as.integer(repro_go$gsea_seed_base),
    derived_seed = seed_file,
    RNGkind = repro_go$rng_kind,
    pAdjustMethod = CFG_PADJ,
    pvalueCutoff = CFG_PCUT,
    qvalueCutoff = paste0(
      "not_passed_to_gseGO (config analysis.qvalue_cutoff=", CFG_QCUT,
      " is applied to ORA enrichGO only, 01_clusterProfiler.r:1866/1890)"),
    number_GO_sets_tested = n_tested,
    FDR_family_definition = paste0(
      "BH over all GO:BP gene sets passing minGSSize/maxGSSize inside this single gseGO call ",
      "(n=", n_tested, "); p.adjust applied in DOSE:::GSEA_fgsea to the full size-filtered ",
      "fgsea result before NA-p rows are dropped; family = one (dataset, contrast, spatial unit); ",
      "no pooling across units, contrasts, datasets or themes, and no re-adjustment in ",
      "02_compareGO.r or in the theme table"),
    software_versions = SOFTWARE_VERSIONS,
    # ---------------- provenance / verification extras ----------------
    source_comparison = cmp,
    route_unit = ru,
    analysis_id = mr$analysis_id,
    run_id = mr$run_id,
    analysis_status = mr$analysis_status,
    config_file = mr$config_file,
    config_hash = mr$config_hash,
    config_hash_matches_current_file = identical(as.character(mr$config_hash), as.character(CFG_HASH)),
    output_table = mr$output_table,
    n_terms_reported_in_canonical_output = n_reported,
    n_sets_dropped_na_pvalue = n_tested - n_reported,
    na_pvalue_dropped_GO_IDs = paste(missing_ids, collapse = ";"),
    gene_universe_identical_to_dataset_reference = universe_identical,
    derived_seed_recomputed = seed_recomputed,
    derived_seed_matches_recomputation = identical(seed_file, seed_recomputed),
    derived_seed_matches_QC_summary = identical(seed_file, as.integer(qc$gsea_go_seed)),
    gsea_input_total_ranked_genes = diag$total_ranked_genes,
    gsea_input_orgdb_symbol_match_fraction = diag$orgdb_symbol_match_fraction,
    gsea_input_duplicated_gene_names = diag$duplicated_gene_names,
    gsea_input_min_rank_statistic = diag$min_rank_statistic,
    gsea_input_max_rank_statistic = diag$max_rank_statistic,
    observed_min_setSize = min(cg_i$setSize),
    observed_max_setSize = max(cg_i$setSize),
    total_ProteinGroupIDs_in_input = agg$total_ProteinGroupIDs,
    excluded_multi_gene_groups = agg$excluded_multi_gene_groups,
    excluded_unresolved_groups = agg$excluded_unresolved_groups,
    excluded_mixed_species_or_contaminant_groups = agg$excluded_mixed_species_or_contaminant_groups,
    stringsAsFactors = FALSE
  )

  # ---- rank statistic audit row ----
  expected_col <- "t"
  fb <- as.character(mr$rank_statistic_fallback_used)
  fails <- character(0)
  if (!identical(as.character(mr$rank_statistic_column), expected_col)) {
    fails <- c(fails, "rank_statistic_column_not_t")
  }
  if (!identical(as.character(mr$rank_statistic_type), "moderated_or_signed_inferential")) {
    fails <- c(fails, "rank_statistic_type_not_moderated")
  }
  if (!identical(toupper(fb), "FALSE")) fails <- c(fails, "fallback_used_TRUE")
  if (!identical(as.character(stat$column), as.character(mr$rank_statistic_column))) {
    fails <- c(fails, "manifest_disagrees_with_reselected_column")
  }
  if (isTRUE(stat$fallback_used)) fails <- c(fails, "reselection_fallback_TRUE")
  if (!identical(as.character(mr$duplicate_gene_collapse_rule), "median_finite_statistics")) {
    fails <- c(fails, "collapse_rule_not_median")
  }
  if (!isTRUE(collapse_reproduced)) fails <- c(fails, "median_collapse_not_reproducible")
  if (!identical(as.integer(agg$genes_with_multiple_ProteinGroupIDs), as.integer(n_multi_re))) {
    fails <- c(fails, "n_genes_multi_group_mismatch")
  }
  if (!identical(as.integer(agg$genes_with_discordant_directions), as.integer(n_disc_re))) {
    fails <- c(fails, "n_genes_discordant_mismatch")
  }

  panels <- if (identical(kk$contrast, "SUS - RES")) {
    "Fig3b v9_atlas; Fig3c v9_bridge; Fig3d-f v9_curve_syn/rna/ox; ED6c-e v9_ed_gsea_curve_syn/rna/ox"
  } else if (identical(kk$contrast, "RES - CON")) {
    "ED6a v9_ed_atlas_rescon; Fig3c v9_bridge; Fig3d-f v9_curve_syn/rna/ox; ED6c-e v9_ed_gsea_curve_syn/rna/ox"
  } else {
    "ED6b v9_ed_atlas_suscon; Fig3c v9_bridge; Fig3d-f v9_curve_syn/rna/ox; ED6c-e v9_ed_gsea_curve_syn/rna/ox"
  }

  rank_rows[[length(rank_rows) + 1L]] <- data.frame(
    dataset = ds,
    comparison = kk$contrast,
    spatial_unit = kk$spatial_unit,
    source_comparison = cmp,
    route_unit = ru,
    analysis_id = mr$analysis_id,
    result_type = mr$result_type,
    ontology = mr$ontology,
    evidence_source_family = EVIDENCE_FAMILY,
    feeds_manuscript_theme_table = TRUE,
    figure3_relevant = TRUE,
    ED6_relevant = TRUE,
    figure_panels = panels,
    rank_statistic_column = mr$rank_statistic_column,
    rank_statistic_type = mr$rank_statistic_type,
    rank_statistic_fallback_used = fb,
    rank_statistic_column_reselected_from_input = stat$column,
    rank_statistic_type_reselected = stat$type,
    rank_statistic_fallback_reselected = stat$fallback_used,
    input_has_moderated_t_column = has_t,
    input_has_log2fc_column = has_log2fc,
    expected_rank_statistic_column = expected_col,
    expected_rank_statistic_type = "moderated_or_signed_inferential",
    duplicate_gene_collapse_rule = mr$duplicate_gene_collapse_rule,
    expected_collapse_rule = "median_finite_statistics",
    median_collapse_independently_reproduced = collapse_reproduced,
    n_ranked_unique_genes = nrow(cgi),
    n_genes_with_multiple_protein_groups_manifest = agg$genes_with_multiple_ProteinGroupIDs,
    n_genes_with_multiple_protein_groups_recomputed = n_multi_re,
    n_genes_with_discordant_directions_manifest = agg$genes_with_discordant_directions,
    n_genes_with_discordant_directions_recomputed = n_disc_re,
    documented_in_manifest = TRUE,
    manifest_path = file.path("data/processed/04_differential_expression_enrichment/clusterProfiler",
                              ds, "clusterProfiler_manifest.csv"),
    status = if (length(fails)) "FAIL" else "PASS",
    fail_reason = if (length(fails)) paste(fails, collapse = ";") else "",
    stringsAsFactors = FALSE
  )

  cat(sprintf("[%2d/54] %-16s %-10s %-8s rank=%s fallback=%s collapse_ok=%s tested=%d reported=%d\n",
              i, ds, kk$contrast, kk$spatial_unit, mr$rank_statistic_column, fb,
              collapse_reproduced, n_tested, n_reported))
  utils::flush.console()
}

contract <- do.call(rbind, rows)
rank_audit <- do.call(rbind, rank_rows)
contract <- contract[order(contract$dataset, contract$comparison, contract$spatial_unit), ]
rank_audit <- rank_audit[order(rank_audit$dataset, rank_audit$comparison, rank_audit$spatial_unit), ]

utils::write.csv(contract, file.path(OUT_DIR, "canonical_gsea_contract_audit.csv"), row.names = FALSE)
utils::write.csv(rank_audit, file.path(OUT_DIR, "gsea_rank_statistic_audit.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 8. Gene-set size-bound audit
# ---------------------------------------------------------------------------
prim <- theme[theme$theme_role == "primary" &
                as.character(theme$theme_claim_eligible) %in% c("TRUE", "True", "true"), ]
anchor_raw <- unique(prim$anchor_GO_ID[!is.na(prim$anchor_GO_ID) & nzchar(prim$anchor_GO_ID)])
anchor_ids <- sort(unique(trimws(unlist(strsplit(anchor_raw, ";", fixed = TRUE)))))
anchor_theme <- vapply(anchor_ids, function(g) {
  paste(sort(unique(prim$manuscript_theme[grepl(g, prim$anchor_GO_ID, fixed = TRUE)])), collapse = "|")
}, character(1))
EXEMPLARS <- c("GO:0099536", "GO:0006397", "GO:0006119")
exemplar_note <- c("GO:0099536" = "Fig3d / ED6c exemplar (synaptic signaling)",
                   "GO:0006397" = "Fig3e / ED6d exemplar (mRNA processing)",
                   "GO:0006119" = "Fig3f / ED6e exemplar (oxidative phosphorylation)")
check_ids <- sort(unique(c(anchor_ids, EXEMPLARS)))
check_lbl <- go_term_name(check_ids)

size_rows <- list()
for (ds in DATASETS) {
  u <- universe_cache[[ds]]
  sizes <- u$sizes
  rep_cmp <- u$ref_comparison
  kk <- keys[keys$dataset == ds & keys$source_comparison == rep_cmp, ]
  cg <- compare_go[[ds]]
  cg_rep <- cg[cg$comparison == rep_cmp, , drop = FALSE]
  n_tested <- sum(sizes >= MIN_GS & sizes <= MAX_GS)
  n_cmp_ds <- length(unique(cg$comparison))

  size_rows[[length(size_rows) + 1L]] <- data.frame(
    record_type = "bounds_and_size_distribution",
    dataset = ds,
    representative_comparison = rep_cmp,
    contrast = kk$contrast,
    spatial_unit = kk$spatial_unit,
    minGSSize_config = CFG_MIN,
    maxGSSize_config = CFG_MAX,
    minGSSize_in_force = MIN_GS,
    maxGSSize_in_force = MAX_GS,
    bounds_source = paste0(
      "config/clusterProfiler_config.yml analysis.min_gs_size / analysis.max_gs_size (md5 ",
      CFG_HASH, ") passed to clusterProfiler::gseGO at ",
      "analysis/04_differential_abundance/01_clusterProfiler.r:1502; gseGO package defaults ",
      "minGSSize=", GSEGO_DEFAULT_MIN, " / maxGSSize=", GSEGO_DEFAULT_MAX, " are overridden"),
    n_go_bp_sets_in_orgdb = N_GO_BP_SETS,
    n_ranked_genes = length(u$genes),
    n_sets_with_zero_overlap = sum(sizes == 0L),
    n_sets_excluded_below_min = sum(sizes < MIN_GS),
    n_sets_excluded_below_min_nonzero = sum(sizes > 0L & sizes < MIN_GS),
    n_sets_excluded_above_max = sum(sizes > MAX_GS),
    n_sets_tested_within_bounds = n_tested,
    n_terms_reported_in_canonical_output = nrow(cg_rep),
    n_sets_dropped_na_pvalue = n_tested - nrow(cg_rep),
    observed_min_setSize_in_output = min(cg_rep$setSize),
    observed_max_setSize_in_output = max(cg_rep$setSize),
    n_comparisons_in_dataset = n_cmp_ds,
    term_GO_ID = NA_character_,
    term_label = NA_character_,
    term_role = NA_character_,
    term_theme = NA_character_,
    term_setSize_in_universe = NA_integer_,
    term_within_bounds = NA,
    term_excluded_reason = NA_character_,
    term_present_in_representative_output = NA,
    term_n_comparisons_present_in_dataset = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (j in seq_along(check_ids)) {
    g <- check_ids[j]
    sz <- if (g %in% names(sizes)) as.integer(sizes[[g]]) else NA_integer_
    within <- !is.na(sz) && sz >= MIN_GS && sz <= MAX_GS
    reason <- if (is.na(sz)) {
      "GO_ID_absent_from_org.Mm.eg.db_GO_BP_collection"
    } else if (sz == 0L) {
      "no_ranked_gene_overlap"
    } else if (sz < MIN_GS) {
      paste0("setSize_", sz, "_below_minGSSize_", MIN_GS)
    } else if (sz > MAX_GS) {
      paste0("setSize_", sz, "_above_maxGSSize_", MAX_GS)
    } else {
      "not_excluded"
    }
    role <- paste(c(if (g %in% anchor_ids) "theme_anchor",
                    if (g %in% EXEMPLARS) "figure_exemplar"), collapse = "+")
    th <- paste(c(if (g %in% anchor_ids) unname(anchor_theme[[g]]),
                  if (g %in% EXEMPLARS) unname(exemplar_note[[g]])), collapse = " | ")
    size_rows[[length(size_rows) + 1L]] <- data.frame(
      record_type = "term_check",
      dataset = ds,
      representative_comparison = rep_cmp,
      contrast = kk$contrast,
      spatial_unit = kk$spatial_unit,
      minGSSize_config = CFG_MIN,
      maxGSSize_config = CFG_MAX,
      minGSSize_in_force = MIN_GS,
      maxGSSize_in_force = MAX_GS,
      bounds_source = "config/clusterProfiler_config.yml analysis.min_gs_size / analysis.max_gs_size",
      n_go_bp_sets_in_orgdb = N_GO_BP_SETS,
      n_ranked_genes = length(u$genes),
      n_sets_with_zero_overlap = NA_integer_,
      n_sets_excluded_below_min = NA_integer_,
      n_sets_excluded_below_min_nonzero = NA_integer_,
      n_sets_excluded_above_max = NA_integer_,
      n_sets_tested_within_bounds = n_tested,
      n_terms_reported_in_canonical_output = NA_integer_,
      n_sets_dropped_na_pvalue = NA_integer_,
      observed_min_setSize_in_output = NA_integer_,
      observed_max_setSize_in_output = NA_integer_,
      n_comparisons_in_dataset = n_cmp_ds,
      term_GO_ID = g,
      term_label = check_lbl[j],
      term_role = role,
      term_theme = th,
      term_setSize_in_universe = sz,
      term_within_bounds = within,
      term_excluded_reason = reason,
      term_present_in_representative_output = g %in% cg_rep$ID,
      term_n_comparisons_present_in_dataset = length(unique(cg$comparison[cg$ID == g])),
      stringsAsFactors = FALSE
    )
  }
}
# Gene sets that passed the size bounds but are absent from the canonical output
# (fgsea returned NA p-value; DOSE:::GSEA_fgsea drops those rows). These are
# absences of a result, not non-significant results.
dropped <- if (length(dropped_records)) do.call(rbind, dropped_records) else NULL
if (!is.null(dropped) && nrow(dropped)) {
  tmpl <- size_rows[[1]]
  for (i in seq_len(nrow(dropped))) {
    r <- tmpl
    ds <- dropped$dataset[i]
    u <- universe_cache[[ds]]
    cg <- compare_go[[ds]]
    g <- dropped$GO_ID[i]
    role <- paste(c(if (g %in% anchor_ids) "theme_anchor",
                    if (g %in% EXEMPLARS) "figure_exemplar"), collapse = "+")
    r$record_type <- "size_filtered_but_absent_from_output"
    r$dataset <- ds
    r$representative_comparison <- dropped$source_comparison[i]
    r$contrast <- dropped$comparison[i]
    r$spatial_unit <- dropped$spatial_unit[i]
    r$n_ranked_genes <- length(u$genes)
    r$n_sets_with_zero_overlap <- NA_integer_
    r$n_sets_excluded_below_min <- NA_integer_
    r$n_sets_excluded_below_min_nonzero <- NA_integer_
    r$n_sets_excluded_above_max <- NA_integer_
    r$n_sets_tested_within_bounds <- sum(u$sizes >= MIN_GS & u$sizes <= MAX_GS)
    r$n_terms_reported_in_canonical_output <-
      sum(cg$comparison == dropped$source_comparison[i])
    r$n_sets_dropped_na_pvalue <- r$n_sets_tested_within_bounds -
      r$n_terms_reported_in_canonical_output
    r$observed_min_setSize_in_output <- NA_integer_
    r$observed_max_setSize_in_output <- NA_integer_
    r$term_GO_ID <- g
    r$term_label <- dropped$GO_description[i]
    r$term_role <- if (nzchar(role)) role else "non_anchor_GO_BP_term"
    r$term_theme <- paste(c(if (g %in% anchor_ids) unname(anchor_theme[[g]]),
                            if (g %in% EXEMPLARS) unname(exemplar_note[[g]])), collapse = " | ")
    r$term_setSize_in_universe <- dropped$setSize_in_universe[i]
    r$term_within_bounds <- TRUE
    r$term_excluded_reason <-
      "within_size_bounds_but_no_row_in_canonical_output_fgsea_NA_pvalue"
    r$term_present_in_representative_output <- FALSE
    r$term_n_comparisons_present_in_dataset <- length(unique(cg$comparison[cg$ID == g]))
    size_rows[[length(size_rows) + 1L]] <- r
  }
}

size_audit <- do.call(rbind, size_rows)
utils::write.csv(size_audit, file.path(OUT_DIR, "gsea_geneset_size_audit.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 9. Console summary
# ---------------------------------------------------------------------------
cat("\n================ SUMMARY ================\n")
cat("rows canonical_gsea_contract_audit.csv    :", nrow(contract), "\n")
cat("rows gsea_rank_statistic_audit.csv        :", nrow(rank_audit), "\n")
cat("rows gsea_geneset_size_audit.csv          :", nrow(size_audit), "\n")
cat("rank_statistic_column values              :", paste(unique(contract$rank_statistic_column), collapse = ","), "\n")
cat("rank_statistic_type values                :", paste(unique(contract$rank_statistic_type), collapse = ","), "\n")
cat("fallback_used values                      :", paste(unique(contract$rank_statistic_fallback_used), collapse = ","), "\n")
cat("collapse rule values                      :", paste(unique(contract$duplicate_gene_collapse_rule), collapse = ","), "\n")
cat("median collapse reproduced (all 54)       :", all(rank_audit$median_collapse_independently_reproduced), "\n")
cat("n PASS / n FAIL                           :", sum(rank_audit$status == "PASS"), "/", sum(rank_audit$status == "FAIL"), "\n")
if (any(rank_audit$status == "FAIL")) {
  print(rank_audit[rank_audit$status == "FAIL", c("dataset", "comparison", "spatial_unit", "fail_reason")],
        row.names = FALSE)
}
cat("seeds match recomputation (all 54)        :", all(contract$derived_seed_matches_recomputation), "\n")
cat("seeds match QC_summary (all 54)           :", all(contract$derived_seed_matches_QC_summary), "\n")
cat("config hash matches current file (all 54) :", all(contract$config_hash_matches_current_file), "\n")
cat("gene universe identical within dataset    :", all(contract$gene_universe_identical_to_dataset_reference), "\n")
cat("n_unique_animals values                   :", paste(sort(unique(contract$n_unique_animals)), collapse = ","), "\n")
cat("animals_per_group values                  :", paste(sort(unique(contract$animals_per_group)), collapse = " / "), "\n")
cat("number_GO_sets_tested range               :", paste(range(contract$number_GO_sets_tested), collapse = "-"), "\n")
cat("n_terms reported range                    :", paste(range(contract$n_terms_reported_in_canonical_output), collapse = "-"), "\n")
cat("n_ranked_unique_genes by dataset          :",
    paste(paste0(names(tapply(contract$n_ranked_unique_genes, contract$dataset, one)), "=",
                 tapply(contract$n_ranked_unique_genes, contract$dataset, one)), collapse = " "), "\n")
cat("observed setSize range across all results :",
    min(contract$observed_min_setSize), "-", max(contract$observed_max_setSize), "\n")
cat("GSEA exponent / eps / nPermSimple         :", GSEA_EXPONENT, "/", GSEA_EPS, "/",
    one(contract$nPermSimple), "\n")
cat("n distinct theme anchor GO IDs checked    :", length(anchor_ids), "\n")
excl <- size_audit[size_audit$record_type == "term_check" & !size_audit$term_within_bounds, ]
cat("anchor/exemplar term-checks excluded      :", nrow(excl), "of", sum(size_audit$record_type == "term_check"), "\n")
if (nrow(excl)) {
  print(excl[, c("dataset", "term_GO_ID", "term_label", "term_setSize_in_universe", "term_excluded_reason")],
        row.names = FALSE)
}
dropped_rows <- size_audit[size_audit$record_type == "size_filtered_but_absent_from_output", ]
cat("size-filtered sets absent from output     :", nrow(dropped_rows), "\n")
if (nrow(dropped_rows)) {
  print(dropped_rows[, c("dataset", "contrast", "spatial_unit", "term_GO_ID", "term_label",
                         "term_setSize_in_universe", "term_role")], row.names = FALSE)
}
cat("=========================================\n")
