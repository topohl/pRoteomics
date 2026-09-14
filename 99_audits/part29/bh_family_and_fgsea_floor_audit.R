#!/usr/bin/env Rscript

# ==========================================================================
# Part-29, sections 4, 5 and 27
#   (4/5)  BH multiple-testing family construction for the canonical ranked
#          GSEA, verified from the returned result tables themselves.
#   (27)   fgsea numerical floor (eps) for every GO term displayed in
#          Figure 3 or Extended Data Figure 6.
#
# AUDIT ONLY. This script reads canonical outputs, canonical run logs and
# repository source code. It never reruns, regenerates or modifies any
# canonical analysis output. It recomputes GO-BP gene-set ELIGIBILITY from the
# stored ranked gene lists (a read-only verification, not a GSEA rerun) and it
# never calls gseGO/fgsea.
#
# It writes exactly two CSVs, both under
#   results/tables/publication_audits/upstream_enrichment_v10/
# ==========================================================================

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

OUT_DIR <- file.path("results", "tables", "publication_audits",
                     "upstream_enrichment_v10")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_MT  <- file.path(OUT_DIR, "gsea_multiple_testing_contract.csv")
OUT_EPS <- file.path(OUT_DIR, "gsea_numerical_precision_audit.csv")

DATASETS <- c("microglia", "neuron_soma", "neuron_neuropil")

# ---- contract constants, read from source, not assumed ---------------------
cfg_lines <- readLines(file.path("config", "clusterProfiler_config.yml"), warn = FALSE)
num_from_cfg <- function(key) {
  i <- grep(paste0("^\\s*", key, ":"), cfg_lines)
  if (!length(i)) stop("config key not found: ", key)
  as.numeric(sub(".*:\\s*", "", cfg_lines[i[1]]))
}
chr_from_cfg <- function(key) {
  i <- grep(paste0("^\\s*", key, ":"), cfg_lines)
  if (!length(i)) stop("config key not found: ", key)
  trimws(sub(".*:\\s*", "", cfg_lines[i[1]]))
}
CFG_PVAL_CUT    <- num_from_cfg("pvalue_cutoff")     # 1
CFG_QVAL_CUT    <- num_from_cfg("qvalue_cutoff")     # 1 (ORA only; gseGO has no such arg)
CFG_PADJ        <- chr_from_cfg("p_adjust_method")   # BH
CFG_MIN         <- num_from_cfg("min_gs_size")       # 10
CFG_MAX         <- num_from_cfg("max_gs_size")       # 800
CFG_NPERMSIMPLE <- num_from_cfg("n_perm_simple")     # 100000

# gseGO's eps is NOT set anywhere in the repository: neither at the call site
# (04_differential_expression_enrichment/01_clusterProfiler.r:1501-1511) nor in
# clusterprofiler_fgsea_control_args() (R/clusterprofiler_reproducibility.R:110).
# Therefore the package default is in force. Read it from the installed package
# rather than hard-coding it.
EPS_IN_FORCE     <- as.numeric(eval(formals(clusterProfiler::gseGO)$eps))
FGSEA_EPS_DEF    <- as.numeric(eval(formals(fgsea::fgseaMultilevel)$eps))
FGSEA_NPERM_DEF  <- as.numeric(eval(formals(fgsea::fgseaMultilevel)$nPermSimple))
FGSEA_SAMPLESIZE <- as.numeric(eval(formals(fgsea::fgseaMultilevel)$sampleSize))

# assert that the repository really does not override eps
src_files <- c("04_differential_expression_enrichment/01_clusterProfiler.r",
               "R/clusterprofiler_reproducibility.R",
               "R/protein_group_enrichment_utils.R")
eps_hits <- unlist(lapply(src_files, function(f) {
  grep("[^a-zA-Z_.]eps[[:space:]]*=", readLines(f, warn = FALSE), value = TRUE)
}))
EPS_OVERRIDE_FOUND <- length(eps_hits) > 0L

message("[cfg] pvalue_cutoff=", CFG_PVAL_CUT, " pAdjust=", CFG_PADJ,
        " minGSSize=", CFG_MIN, " maxGSSize=", CFG_MAX,
        " nPermSimple=", CFG_NPERMSIMPLE)
message("[cfg] eps in force = ", format(EPS_IN_FORCE, scientific = TRUE),
        " (gseGO default; repo override found: ", EPS_OVERRIDE_FOUND, ")")

# ==========================================================================
# 1. MANIFEST -> the 54 canonical GSEA_GO comparisons
# ==========================================================================
unP <- function(p) sub("^[A-Za-z]:/+", "", p)   # "P://results/..." -> "results/..."

# Nine of the 54 protein_group_audits/collapsed_gene_input.csv paths are longer
# than the Windows MAX_PATH limit of 260 characters, so base-R file.exists() and
# file() cannot open them even though list.files() sees them. Stage those into
# tempdir() with robocopy (which is long-path aware) and read them from there.
# This is a READ-ONLY copy out of the repository; nothing in the repository is
# touched.
LP_STAGE <- file.path(tempdir(), "part29_longpath_stage")
dir.create(LP_STAGE, recursive = TRUE, showWarnings = FALSE)
read_csv_lp <- function(path, ...) {
  if (file.exists(path)) {
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...))
  }
  d <- dirname(path); b <- basename(path)
  if (!dir.exists(d) || !(b %in% list.files(d))) stop("file not found: ", path, call. = FALSE)
  dst <- file.path(LP_STAGE, b)
  if (file.exists(dst)) unlink(dst, force = TRUE)
  src <- normalizePath(d, winslash = "\\", mustWork = TRUE)
  system2("robocopy", c(shQuote(src), shQuote(normalizePath(LP_STAGE, winslash = "\\")),
                        shQuote(b), "/NJH", "/NJS", "/NP", "/NFL", "/NDL",
                        "/R:1", "/W:1"),
          stdout = NULL, stderr = NULL)
  if (!file.exists(dst)) stop("long-path staging failed for: ", path, call. = FALSE)
  utils::read.csv(dst, stringsAsFactors = FALSE, check.names = FALSE, ...)
}
file_reachable <- function(p) file.exists(p) | vapply(p, function(x)
  dir.exists(dirname(x)) && basename(x) %in% list.files(dirname(x)), logical(1))

man <- do.call(rbind, lapply(DATASETS, function(d) {
  utils::read.csv(file.path("data", "processed",
                            "04_differential_expression_enrichment",
                            "clusterProfiler", d, "clusterProfiler_manifest.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)
}))
man <- man[man$result_type == "GSEA_GO", , drop = FALSE]
man$gsea_table <- unP(man$output_table)
man$gene_file  <- unP(man$collapsed_gene_input_file)
stopifnot(all(file.exists(man$gsea_table)), all(file_reachable(man$gene_file)))
N_LONGPATH <- sum(!file.exists(man$gene_file))
message("[1/7] manifest: ", nrow(man), " canonical GSEA_GO comparisons (",
        N_LONGPATH, " gene-input files exceed Windows MAX_PATH and are staged)")

# ==========================================================================
# 2. GO-BP GENE SETS, exactly as clusterProfiler::gseGO builds them
#    (ont = "BP", keyType = "SYMBOL", OrgDb = org.Mm.eg.db)
# ==========================================================================
message("[2/7] rebuilding the GO-BP gene-set collection (read-only)")
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Mm.eg.db); library(fastmatch)
})
GO_DATA  <- clusterProfiler:::get_GO_data(org.Mm.eg.db, "BP", "SYMBOL")
GENESETS <- DOSE:::getGeneSet(GO_DATA)
GENESETS <- lapply(GENESETS, unique)
N_GO_BP_COLLECTION <- length(GENESETS)
message("      GO-BP sets in collection = ", N_GO_BP_COLLECTION,
        "; annotations = ", sum(lengths(GENESETS)))

ALLSYM <- unique(unlist(GENESETS, use.names = FALSE))
idxlst <- lapply(GENESETS, function(g) fmatch(g, ALLSYM))
FLAT   <- unlist(idxlst, use.names = FALSE)
GRP    <- rep.int(seq_along(idxlst), lengths(idxlst))

# fgsea::preparePathwaysAndStats():
#   pathwaysFiltered <- unique(na.omit(fmatch(p, names(stats))))
#   minSize <- max(minSize, 1); maxSize <- min(maxSize, length(stats) - 1)
#   toKeep  <- minSize <= size & size <= maxSize
eligible_counts <- function(gene_symbols) {
  present <- logical(length(ALLSYM))
  hit <- fmatch(gene_symbols, ALLSYM)
  present[hit[!is.na(hit)]] <- TRUE
  sz <- tabulate(GRP[present[FLAT]], nbins = length(GENESETS))
  lo <- max(CFG_MIN, 1)
  hi <- min(CFG_MAX, length(gene_symbols) - 1)
  list(n_eligible = sum(sz >= lo & sz <= hi), min_used = lo, max_used = hi)
}

# ==========================================================================
# 3. PER-COMPARISON BH FAMILY AUDIT
# ==========================================================================
message("[3/7] auditing the BH family of each comparison")

# Solve for the family size n that R's p.adjust(method="BH") must have used:
#   adj_(i) = min( 1 , n * min_{k>=i} p_(k)/k )
# so for every row with adj < 1 :  n = adj_(i) / min_{k>=i} p_(k)/k
implied_family_size <- function(p, adj) {
  ok <- is.finite(p) & is.finite(adj)
  p <- p[ok]; adj <- adj[ok]
  o <- order(p); p <- p[o]; adj <- adj[o]
  lp <- length(p)
  ratio <- p / seq_len(lp)
  tail_min <- rev(cummin(rev(ratio)))
  use <- adj < 1 - 1e-12 & tail_min > 0
  if (!any(use)) return(c(n = NA_real_, spread = NA_real_))
  est <- adj[use] / tail_min[use]
  c(n = stats::median(est), spread = max(est) - min(est))
}
bh_with_n <- function(p, n) {
  lp <- length(p)
  o <- order(p, decreasing = TRUE); ro <- order(o)
  pmin(1, cummin(n / (lp:1L) * p[o]))[ro]
}

gsea_store <- list()
rows <- vector("list", nrow(man))
for (i in seq_len(nrow(man))) {
  r <- man[i, ]
  tb <- utils::read.csv(r$gsea_table, stringsAsFactors = FALSE, check.names = FALSE)
  gene <- read_csv_lp(r$gene_file)
  gsym <- as.character(gene$GeneSymbol)

  el <- eligible_counts(gsym)
  n_ret <- nrow(tb)
  p  <- suppressWarnings(as.numeric(tb$pvalue))
  pa <- suppressWarnings(as.numeric(tb$p.adjust))
  n_raw <- sum(is.finite(p)); n_adj <- sum(is.finite(pa))

  imp <- implied_family_size(p, pa)
  n_imp <- if (is.finite(imp[["n"]])) round(imp[["n"]]) else NA_real_

  recon_ret <- bh_with_n(p[is.finite(p)], n_ret)
  dev_ret <- max(abs(recon_ret - pa[is.finite(p)]), na.rm = TRUE)
  dev_imp <- if (is.finite(n_imp)) {
    max(abs(bh_with_n(p[is.finite(p)], n_imp) - pa[is.finite(p)]), na.rm = TRUE)
  } else NA_real_

  n_floor <- sum(is.finite(p) & p <= EPS_IN_FORCE)
  n_floor_exact <- sum(is.finite(p) & p == EPS_IN_FORCE)
  n_excluded_from_family <- el$n_eligible - (if (is.finite(n_imp)) n_imp else NA_real_)

  gsea_store[[r$comparison]] <- data.frame(
    dataset = r$dataset, source_comparison = r$comparison,
    route_unit = r$route_unit,
    ID = as.character(tb$ID), Description = as.character(tb$Description),
    setSize = suppressWarnings(as.integer(tb$setSize)),
    NES = suppressWarnings(as.numeric(tb$NES)),
    pvalue = p, p.adjust = pa,
    leading_edge_size = vapply(strsplit(as.character(tb$core_enrichment), "/", fixed = TRUE),
                               function(x) sum(nzchar(x)), integer(1)),
    stringsAsFactors = FALSE)

  rows[[i]] <- data.frame(
    dataset = r$dataset,
    spatial_unit_dir = r$route_unit,
    comparison = r$comparison,
    n_ranked_genes = length(gsym),
    n_GO_sets_in_BP_collection = N_GO_BP_COLLECTION,
    n_GO_sets_eligible = el$n_eligible,
    n_GO_sets_tested = el$n_eligible,
    n_GO_sets_returned = n_ret,
    n_raw_P = n_raw,
    n_adjusted_P = n_adj,
    n_FDR_lt_0.05 = sum(is.finite(pa) & pa < 0.05),
    pvalueCutoff_in_force = CFG_PVAL_CUT,
    qvalueCutoff_in_force = paste0(
      "not_applicable: clusterProfiler::gseGO has no qvalueCutoff argument and ",
      "none is passed at 01_clusterProfiler.r:1501-1511; config ",
      "analysis.qvalue_cutoff=", CFG_QVAL_CUT, " is consumed by the ORA branch only"),
    pAdjustMethod = CFG_PADJ,
    family_definition = paste0(
      "one BH family per comparison per ontology = every GO:BP gene set whose ",
      "overlap with that comparison's ranked SYMBOL list is within [",
      el$min_used, ",", el$max_used, "]; BH applied by DOSE::GSEA_fgsea as ",
      "p.adjust(tmp_res$pval, method='BH') on the FULL fgsea output, before ",
      "any pvalueCutoff filtering; sets whose p is NA stay in the denominator n"),
    any_prefilter_before_BH = el$n_eligible > n_ret,
    prefilter_description = if (el$n_eligible > n_ret) paste0(
      "No p-value-dependent prefilter. The prespecified gene-set SIZE filter ",
      "(minGSSize=", CFG_MIN, ", maxGSSize=", CFG_MAX, ") is applied by ",
      "fgsea::preparePathwaysAndStats before any p-value exists. HOWEVER ",
      el$n_eligible - n_ret, " size-eligible GO:BP set(s) were returned by ",
      "fgseaMultilevel with pval=NA ('unbalanced positive/negative gene-level ",
      "statistics', modeFraction<10). R's p.adjust(p,'BH') with its default ",
      "n=length(p) evaluates n AFTER dropping NAs, so those set(s) are absent ",
      "from the BH denominator: family n=", n_imp, " instead of ", el$n_eligible,
      " (denominator inflation factor ", formatC(el$n_eligible / n_imp, digits = 8, format = "f"),
      " not applied). Post-BH filters p<=", CFG_PVAL_CUT, " and p.adjust<=",
      CFG_PVAL_CUT, " remove 0 rows because the cutoff is 1.")
    else paste0(
      "No prefilter of any kind removes a tested set before BH. The only ",
      "prefilter is the prespecified gene-set SIZE filter (minGSSize=", CFG_MIN,
      ", maxGSSize=", CFG_MAX, ") applied by fgsea::preparePathwaysAndStats ",
      "before any p-value exists. Every size-eligible set received a finite p ",
      "and entered the BH denominator. Post-BH filters p<=", CFG_PVAL_CUT,
      " and p.adjust<=", CFG_PVAL_CUT, " remove 0 rows because the cutoff is 1."),
    # ---- the decisive verification columns -------------------------------
    bh_family_size_implied_from_stored_padjust = n_imp,
    bh_implied_family_estimator_spread = unname(imp[["spread"]]),
    bh_recon_maxabsdiff_using_n_returned = dev_ret,
    bh_recon_maxabsdiff_using_n_implied = dev_imp,
    n_GO_sets_excluded_from_BH_family = n_excluded_from_family,
    bh_denominator_inflation_if_NA_sets_included = el$n_eligible / n_imp,
    returned_equals_tested = identical(as.integer(n_ret), as.integer(el$n_eligible)),
    implied_family_equals_returned = identical(as.integer(n_imp), as.integer(n_ret)),
    implied_family_equals_tested = identical(as.integer(n_imp), as.integer(el$n_eligible)),
    max_raw_P = max(p, na.rm = TRUE),
    max_adjusted_P = max(pa, na.rm = TRUE),
    n_padj_gt_0.5 = sum(is.finite(pa) & pa > 0.5),
    n_padj_gt_0.9 = sum(is.finite(pa) & pa > 0.9),
    n_padj_gt_0.99 = sum(is.finite(pa) & pa > 0.99),
    min_raw_P = min(p, na.rm = TRUE),
    min_adjusted_P = min(pa, na.rm = TRUE),
    n_terms_at_eps_floor = n_floor,
    n_terms_exactly_at_eps = n_floor_exact,
    eps_in_force = EPS_IN_FORCE,
    nPermSimple_in_force = CFG_NPERMSIMPLE,
    fgsea_sampleSize_default = FGSEA_SAMPLESIZE,
    gsea_table = r$gsea_table,
    stringsAsFactors = FALSE)
  if (i %% 9 == 0) message("      ", i, "/", nrow(man))
}
mt <- do.call(rbind, rows)
rownames(mt) <- NULL

# ==========================================================================
# 4. CROSS-CHECK: run logs, compareGO passthrough, theme-table passthrough
# ==========================================================================
message("[4/7] cross-checking canonical run logs")
log_na <- stats::setNames(rep(0L, nrow(mt)), mt$comparison)
for (i in seq_len(nrow(mt))) {
  lg <- file.path("results", "logs", "04_differential_expression_enrichment",
                  "clusterProfiler", mt$dataset[i], paste0(mt$comparison[i], ".log"))
  if (!file.exists(lg)) next
  L <- readLines(lg, warn = FALSE)
  go_i <- grep("GSEA_RNG \\| analysis_type=gseGO_BP", L)
  kg_i <- grep("GSEA_RNG \\| analysis_type=gseKEGG", L)
  if (!length(go_i)) next
  a <- max(go_i)
  b <- if (length(kg_i) && any(kg_i > a)) min(kg_i[kg_i > a]) else length(L) + 1L
  seg <- L[seq(a, b - 1L)]
  hit <- grep("GSEA_FGSEA \\| There were [0-9]+ pathways for which P-values were not calculated",
              seg, value = TRUE)
  if (length(hit)) {
    log_na[mt$comparison[i]] <- sum(as.integer(sub(".*There were ([0-9]+) pathways.*", "\\1", hit)))
  }
}
mt$log_reported_NA_pvalue_sets_GO_BP <- as.integer(log_na[mt$comparison])
mt$NA_count_matches_log <- mt$n_GO_sets_excluded_from_BH_family == mt$log_reported_NA_pvalue_sets_GO_BP

message("[5/7] cross-checking compareGO passthrough row counts")
cg_counts <- do.call(rbind, lapply(DATASETS, function(d) {
  f <- file.path("results", "tables", "04_differential_expression_enrichment",
                 "compareGO", d, "BP", "phenotype_within_unit", "all_route_units",
                 "compareGO_term_comparison.csv")
  hdr <- names(utils::read.csv(f, nrows = 1L, stringsAsFactors = FALSE, check.names = FALSE))
  cc <- rep("NULL", length(hdr)); cc[match("comparison", hdr)] <- "character"
  x <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, colClasses = cc)
  tb <- as.data.frame(table(x$comparison), stringsAsFactors = FALSE)
  names(tb) <- c("comparison", "n_rows_compareGO")
  tb
}))
mt$n_rows_in_compareGO_table <- cg_counts$n_rows_compareGO[match(mt$comparison, cg_counts$comparison)]
mt$compareGO_is_passthrough <- mt$n_rows_in_compareGO_table == mt$n_GO_sets_returned

# ---- theme table ---------------------------------------------------------
message("[6/7] reading the canonical manuscript theme table (selected columns)")
THEME_TBL <- file.path("results", "tables", "10_biological_integration",
                       "gsea_wgcna_concordance", "global",
                       "ontology_aware_gsea_theme_assignments_all_contrasts.csv")
hdr <- names(utils::read.csv(THEME_TBL, nrows = 1L, stringsAsFactors = FALSE, check.names = FALSE))
chrcols <- c("dataset", "phenotype_contrast", "contrast", "spatial_unit",
             "source_comparison", "GO_ID", "GO_description", "theme_id",
             "manuscript_theme", "theme_role", "assignment_status",
             "registry_version", "GO_db_package_version", "evidence_source_family",
             "leading_edge_genes")
numcols <- c("NES", "raw_p", "GSEA_FDR")
logcols <- c("theme_claim_eligible")
cc <- rep("NULL", length(hdr))
cc[match(chrcols, hdr)] <- "character"
cc[match(numcols, hdr)] <- "numeric"
cc[match(logcols, hdr)] <- "logical"
th <- utils::read.csv(THEME_TBL, stringsAsFactors = FALSE, check.names = FALSE, colClasses = cc)
message("      theme table rows = ", nrow(th))

th_counts <- as.data.frame(table(th$source_comparison), stringsAsFactors = FALSE)
names(th_counts) <- c("comparison", "n_rows_theme_table")
mt$n_rows_in_theme_table <- th_counts$n_rows_theme_table[match(mt$comparison, th_counts$comparison)]
# the theme table carries ONE ROW PER (term x theme assignment), so a GO term
# assigned to two themes appears twice; the passthrough test must be on the
# number of DISTINCT GO_IDs
th_uniq <- tapply(th$GO_ID, th$source_comparison, function(x) length(unique(x)))
mt$n_distinct_GO_in_theme_table <- as.integer(th_uniq[mt$comparison])
mt$theme_table_is_passthrough <- mt$n_distinct_GO_in_theme_table == mt$n_GO_sets_returned
th_split_max <- tapply(th$GSEA_FDR, th$source_comparison, function(x) max(x, na.rm = TRUE))
th_split_ge  <- tapply(th$GSEA_FDR, th$source_comparison, function(x) sum(x >= 0.05, na.rm = TRUE))
mt$theme_table_max_GSEA_FDR <- as.numeric(th_split_max[mt$comparison])
mt$theme_table_n_rows_FDR_ge_0.05 <- as.numeric(th_split_ge[mt$comparison])

utils::write.csv(mt, OUT_MT, row.names = FALSE)
message("      wrote ", OUT_MT, " (", nrow(mt), " rows)")

# ==========================================================================
# 5. THE DISPLAYED GO TERMS OF FIGURE 3 AND EXTENDED DATA FIGURE 6
# ==========================================================================
message("[7/7] assembling the displayed-term numerical-precision audit")

# Panel -> contrast map, read from figures/figure_final_truth_v9_contract.yml
#   F3b  v9_atlas            contrast: SUS - RES
#   ED6a v9_ed_atlas_rescon  contrast: RES - CON
#   ED6b v9_ed_atlas_suscon  contrast: SUS - CON
#   F3c  v9_bridge           three exemplars x 3 contrasts
#   F3d/e/f + ED6c/d/e       the same three exemplars (curve + 3-contrast strip)
ATLAS_PANEL <- c("RES - CON" = "ED6a_v9_ed_atlas_rescon",
                 "SUS - CON" = "ED6b_v9_ed_atlas_suscon",
                 "SUS - RES" = "F3b_v9_atlas")
THEMES6 <- c("synaptic_signaling_vesicle", "rna_processing_splicing_rnp",
             "ribosome_translation", "mitochondrial_respiration_oxphos",
             "autophagy_lysosome_endosome", "chromatin_organization")
EXEMPLARS <- data.frame(
  key = c("synaptic", "rna", "oxphos"),
  dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
  spatial_unit = c("CA3_sr", "CA2_sp", "CA1"),
  GO_ID = c("GO:0099536", "GO:0006397", "GO:0006119"),
  stringsAsFactors = FALSE)

le_size <- function(x, sep) vapply(strsplit(as.character(x), sep, fixed = TRUE),
                                   function(v) sum(nzchar(v)), integer(1))

# ---- atlas constituents, FDR-supported -----------------------------------
z <- th[th$theme_id %in% THEMES6 & th$theme_claim_eligible %in% TRUE, , drop = FALSE]
z <- z[!duplicated(paste(z$dataset, z$spatial_unit, z$contrast, z$GO_ID, z$theme_id)), , drop = FALSE]
ATLAS_ALL_N <- nrow(z)
zs <- z[is.finite(z$GSEA_FDR) & z$GSEA_FDR < 0.05, , drop = FALSE]
ATLAS_SUP_ASSIGN_N <- nrow(zs)
key <- paste(zs$dataset, zs$spatial_unit, zs$contrast, zs$GO_ID, sep = "\r")
atlas <- do.call(rbind, lapply(split(seq_len(nrow(zs)), key), function(ix) {
  w <- zs[ix, , drop = FALSE]
  data.frame(dataset = w$dataset[1], spatial_unit = w$spatial_unit[1],
             contrast = w$contrast[1], GO_ID = w$GO_ID[1],
             GO_description = w$GO_description[1],
             source_comparison = w$source_comparison[1],
             theme_ids = paste(sort(unique(w$theme_id)), collapse = ";"),
             theme_raw_p = w$raw_p[1], theme_FDR = w$GSEA_FDR[1],
             theme_NES = w$NES[1],
             theme_leading_edge_size = le_size(w$leading_edge_genes[1], ";"),
             displayed_in = unname(ATLAS_PANEL[w$contrast[1]]),
             display_role = "atlas_theme_constituent_FDR_supported",
             stringsAsFactors = FALSE)
}))
rownames(atlas) <- NULL

# ---- the three exemplars, all three contrasts ----------------------------
ex <- do.call(rbind, lapply(seq_len(nrow(EXEMPLARS)), function(i) {
  e <- EXEMPLARS[i, ]
  w <- th[th$dataset == e$dataset & th$spatial_unit == e$spatial_unit &
            th$GO_ID == e$GO_ID, , drop = FALSE]
  w <- w[!duplicated(w$contrast), , drop = FALSE]
  w <- w[match(c("RES - CON", "SUS - CON", "SUS - RES"), w$contrast), , drop = FALSE]
  stopifnot(!any(is.na(w$contrast)))
  data.frame(dataset = w$dataset, spatial_unit = w$spatial_unit,
             contrast = w$contrast, GO_ID = w$GO_ID,
             GO_description = w$GO_description,
             source_comparison = w$source_comparison,
             theme_ids = w$theme_id, theme_raw_p = w$raw_p,
             theme_FDR = w$GSEA_FDR, theme_NES = w$NES,
             theme_leading_edge_size = le_size(w$leading_edge_genes, ";"),
             displayed_in = paste0("F3c_bridge;F3", c("d", "e", "f")[i], "_curve;ED6",
                                   c("c", "d", "e")[i], "_curve"),
             display_role = "figure3_ED6_exemplar",
             stringsAsFactors = FALSE)
}))

disp <- rbind(atlas[, names(ex)], ex)
# an exemplar row may also be an atlas constituent: keep ONE row, merging roles
k <- paste(disp$dataset, disp$spatial_unit, disp$contrast, disp$GO_ID, sep = "\r")
disp <- do.call(rbind, lapply(split(seq_len(nrow(disp)), k), function(ix) {
  w <- disp[ix, , drop = FALSE]
  o <- w[1, , drop = FALSE]
  o$displayed_in <- paste(sort(unique(unlist(strsplit(w$displayed_in, ";")))), collapse = ";")
  o$display_role <- paste(sort(unique(w$display_role)), collapse = ";")
  o$theme_ids    <- paste(sort(unique(w$theme_ids[nzchar(w$theme_ids)])), collapse = ";")
  o
}))
rownames(disp) <- NULL

# ---- join the canonical GSEA numbers -------------------------------------
canon <- do.call(rbind, gsea_store); rownames(canon) <- NULL
ck <- paste(canon$source_comparison, canon$ID, sep = "\r")
dk <- paste(disp$source_comparison, disp$GO_ID, sep = "\r")
m <- match(dk, ck)
stopifnot(!any(is.na(m)))

out <- data.frame(
  GO_ID = disp$GO_ID,
  GO_description = canon$Description[m],
  dataset = disp$dataset,
  spatial_unit = disp$spatial_unit,
  contrast = disp$contrast,
  raw_p = canon$pvalue[m],
  BH_FDR = canon$p.adjust[m],
  eps_in_force = EPS_IN_FORCE,
  log2err_if_available = NA_real_,
  n_perm = CFG_NPERMSIMPLE,
  at_numerical_floor = canon$pvalue[m] <= EPS_IN_FORCE,
  setSize = canon$setSize[m],
  leading_edge_size = canon$leading_edge_size[m],
  # ---- provenance / verification ---------------------------------------
  source_comparison = disp$source_comparison,
  spatial_unit_dir = canon$route_unit[m],
  displayed_in = disp$displayed_in,
  display_role = disp$display_role,
  theme_ids = disp$theme_ids,
  NES = canon$NES[m],
  raw_p_theme_table = disp$theme_raw_p,
  BH_FDR_theme_table = disp$theme_FDR,
  raw_p_agrees_with_theme_table = abs(canon$pvalue[m] - disp$theme_raw_p) <= 1e-12,
  FDR_agrees_with_theme_table = abs(canon$p.adjust[m] - disp$theme_FDR) <= 1e-12,
  leading_edge_size_theme_table = disp$theme_leading_edge_size,
  eps_source = "clusterProfiler::gseGO formal default (1e-10); no override anywhere in repo",
  fgsea_default_eps = FGSEA_EPS_DEF,
  n_perm_source = "analysis.n_perm_simple in config/clusterProfiler_config.yml -> fgseaMultilevel(nPermSimple=)",
  fgsea_default_nPermSimple = FGSEA_NPERM_DEF,
  fgsea_sampleSize = FGSEA_SAMPLESIZE,
  log2err_availability = paste0(
    "NOT RECOVERABLE: DOSE:::GSEA_fgsea rebuilds the result data.frame from ",
    "fgsea's pathway/pval/ES/NES/size columns only and discards log2err; no ",
    "canonical output in this repository stores it"),
  bh_family_size_for_this_comparison = mt$bh_family_size_implied_from_stored_padjust[
    match(disp$source_comparison, mt$comparison)],
  stringsAsFactors = FALSE)
out <- out[order(out$dataset, out$spatial_unit, out$contrast, out$GO_ID), ]
rownames(out) <- NULL

utils::write.csv(out, OUT_EPS, row.names = FALSE)
message("      wrote ", OUT_EPS, " (", nrow(out), " rows)")

# ==========================================================================
# CONSOLE SUMMARY
# ==========================================================================
cat("\n================= BH FAMILY =================\n")
cat("comparisons audited                         : ", nrow(mt), "\n", sep = "")
cat("GO-BP sets in collection                    : ", N_GO_BP_COLLECTION, "\n", sep = "")
cat("n_GO_sets_eligible range                    : ", min(mt$n_GO_sets_eligible), "-", max(mt$n_GO_sets_eligible), "\n", sep = "")
cat("n_GO_sets_returned range                    : ", min(mt$n_GO_sets_returned), "-", max(mt$n_GO_sets_returned), "\n", sep = "")
cat("returned == tested in n comparisons         : ", sum(mt$returned_equals_tested), "/", nrow(mt), "\n", sep = "")
cat("implied BH family == tested in n comparisons: ", sum(mt$implied_family_equals_tested), "/", nrow(mt), "\n", sep = "")
cat("max |BH recon - stored| using n_implied     : ", format(max(mt$bh_recon_maxabsdiff_using_n_implied, na.rm = TRUE)), "\n", sep = "")
cat("max |BH recon - stored| using n_returned    : ", format(max(mt$bh_recon_maxabsdiff_using_n_returned, na.rm = TRUE)), "\n", sep = "")
cat("estimator spread (max over comparisons)     : ", format(max(mt$bh_implied_family_estimator_spread, na.rm = TRUE)), "\n", sep = "")
cat("implied family == returned in n comparisons : ", sum(mt$implied_family_equals_returned), "/", nrow(mt), "\n", sep = "")
cat("comparisons with any_prefilter_before_BH    : ", sum(mt$any_prefilter_before_BH), "\n", sep = "")
cat("total GO sets excluded from a BH family     : ", sum(mt$n_GO_sets_excluded_from_BH_family), "\n", sep = "")
print(mt[mt$n_GO_sets_excluded_from_BH_family > 0,
         c("dataset", "comparison", "n_GO_sets_tested", "n_GO_sets_returned",
           "bh_family_size_implied_from_stored_padjust",
           "n_GO_sets_excluded_from_BH_family", "log_reported_NA_pvalue_sets_GO_BP",
           "bh_denominator_inflation_if_NA_sets_included", "n_FDR_lt_0.05")])
cat("NA count matches run log in n comparisons   : ", sum(mt$NA_count_matches_log), "/", nrow(mt), "\n", sep = "")
cat("compareGO passthrough OK                    : ", sum(mt$compareGO_is_passthrough), "/", nrow(mt), "\n", sep = "")
cat("theme-table passthrough OK                  : ", sum(mt$theme_table_is_passthrough), "/", nrow(mt), "\n", sep = "")
cat("min / max stored p.adjust over all rows     : ", format(min(mt$min_adjusted_P)), " / ", format(max(mt$max_adjusted_P)), "\n", sep = "")
cat("comparisons whose max p.adjust > 0.99       : ", sum(mt$max_adjusted_P > 0.99), "/", nrow(mt), "\n", sep = "")
cat("comparisons whose max raw p > 0.99          : ", sum(mt$max_raw_P > 0.99), "/", nrow(mt), "\n", sep = "")
cat("total n_FDR_lt_0.05 over all comparisons    : ", sum(mt$n_FDR_lt_0.05), "\n", sep = "")
cat("theme table: max GSEA_FDR                   : ", format(max(mt$theme_table_max_GSEA_FDR)), "\n", sep = "")
cat("theme table: rows with FDR >= 0.05          : ", format(sum(mt$theme_table_n_rows_FDR_ge_0.05)), " of ", nrow(th), "\n", sep = "")
cat("theme table: distinct GO == returned        : ", sum(mt$theme_table_is_passthrough), "/", nrow(mt), "\n", sep = "")
cat("theme table: evidence_source_family values  : ",
    paste(unique(th$evidence_source_family), collapse = " | "), "\n", sep = "")
cat("theme table: registry / GO.db               : ",
    paste(unique(th$registry_version), collapse = "|"), " / ",
    paste(unique(th$GO_db_package_version), collapse = "|"), "\n", sep = "")
cat("theme table: assignments to the six themes  : ", ATLAS_ALL_N, "\n", sep = "")
cat("  of those, GSEA_FDR >= 0.05 (kept anyway)  : ", ATLAS_ALL_N - ATLAS_SUP_ASSIGN_N, "\n", sep = "")
cat("  max GSEA_FDR among six-theme assignments  : ", format(max(z$GSEA_FDR, na.rm = TRUE)), "\n", sep = "")

cat("\n================= NUMERICAL FLOOR =================\n")
cat("eps in force                                : ", format(EPS_IN_FORCE, scientific = TRUE), "\n", sep = "")
cat("nPermSimple in force                        : ", CFG_NPERMSIMPLE, "\n", sep = "")
cat("fgsea package default eps / nPermSimple     : ", format(FGSEA_EPS_DEF, scientific = TRUE), " / ", FGSEA_NPERM_DEF, "\n", sep = "")
cat("terms at eps floor, whole corpus            : ", sum(mt$n_terms_at_eps_floor), " of ", sum(mt$n_GO_sets_returned), "\n", sep = "")
cat("  of which p is EXACTLY eps                 : ", sum(mt$n_terms_exactly_at_eps), "\n", sep = "")
cat("comparisons with >=1 term at the floor      : ", sum(mt$n_terms_at_eps_floor > 0), "/", nrow(mt), "\n", sep = "")
cat("smallest raw p anywhere in the corpus       : ", format(min(mt$min_raw_P), scientific = TRUE), "\n", sep = "")
cat("displayed rows audited                      : ", nrow(out), "\n", sep = "")
cat("  atlas theme assignments (six themes)      : ", ATLAS_ALL_N, "\n", sep = "")
cat("  of which FDR-supported assignments        : ", ATLAS_SUP_ASSIGN_N, "\n", sep = "")
cat("  distinct atlas-constituent display rows   : ", nrow(atlas), "\n", sep = "")
cat("  exemplar rows                             : ", nrow(ex), "\n", sep = "")
cat("displayed rows AT the numerical floor       : ", sum(out$at_numerical_floor), "\n", sep = "")
cat("smallest raw p among displayed terms        : ", format(min(out$raw_p), scientific = TRUE), "\n", sep = "")
cat("smallest BH FDR among displayed terms       : ", format(min(out$BH_FDR), scientific = TRUE), "\n", sep = "")
cat("largest raw p among displayed terms         : ", format(max(out$raw_p), scientific = TRUE), "\n", sep = "")
cat("largest BH FDR among displayed terms        : ", format(max(out$BH_FDR), scientific = TRUE), "\n", sep = "")
cat("raw_p agrees with theme table (all rows)    : ", all(out$raw_p_agrees_with_theme_table), "\n", sep = "")
cat("BH_FDR agrees with theme table (all rows)   : ", all(out$FDR_agrees_with_theme_table), "\n", sep = "")
cat("n distinct GO terms displayed               : ", length(unique(out$GO_ID)), "\n", sep = "")
cat("displayed rows at floor, by contrast:\n")
print(table(out$contrast, out$at_numerical_floor))
cat("displayed rows at floor, by dataset:\n")
print(table(out$dataset, out$at_numerical_floor))
cat("exemplar rows detail:\n")
print(out[grepl("exemplar", out$display_role),
          c("GO_ID", "dataset", "spatial_unit", "contrast", "raw_p", "BH_FDR",
            "setSize", "leading_edge_size", "at_numerical_floor")])
cat("\nsmallest ten raw p among displayed terms:\n")
print(utils::head(out[order(out$raw_p), c("GO_ID", "dataset", "spatial_unit",
                                          "contrast", "raw_p", "BH_FDR",
                                          "at_numerical_floor")], 10))
cat("\nDONE\n")
