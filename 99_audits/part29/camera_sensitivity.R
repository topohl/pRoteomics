#!/usr/bin/env Rscript

# Part-29 sections 19-27: an orthogonal competitive gene-set sensitivity
# analysis, using limma's CAMERA framework, which accounts for within-set
# correlation that GSEA's permutation of gene labels does not.
#
# THIS IS A SENSITIVITY ANALYSIS, NOT A SECOND PRIMARY ENDPOINT (section 39).
# The canonical ranked GSEA is unchanged and remains the program-level evidence.
# Nothing here reselects a Figure-3 exemplar, creates a theme-level p-value, or
# alters a multiple-testing family.
#
# WHY cameraPR AND NOT camera(). A full camera() needs a GENE-level expression
# matrix. The repository stores a PROTEIN-GROUP-level animal matrix and collapses
# to genes only at the level of the STATISTIC (median of moderated t per official
# gene symbol). Collapsing abundances instead and refitting would not reproduce
# the canonical statistic - the median is not linear, and eBayes would be
# re-estimated over a different number of rows - so a gene-level matrix would
# have to be invented for this audit, which section 20 forbids. cameraPR is
# therefore run on the EXACT canonical ranked statistic, which is the ranking
# the manuscript's GSEA actually used, and is documented as the weaker
# preranked variant.
#
# inter.gene.cor is PRESPECIFIED at the limma-recommended 0.01 and is never
# tuned to improve agreement with GSEA. Two further fixed values are reported
# as a sensitivity band, also prespecified.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages({
  library(limma); library(org.Mm.eg.db); library(GO.db); library(AnnotationDbi)
})
source(file.path("R", "protein_group_enrichment_utils.R"))

AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
dir.create(AUD, recursive = TRUE, showWarnings = FALSE)

MAPPED <- file.path("data", "processed", "02_id_mapping_animal_level", "mapped")
MIN_GS <- 10L; MAX_GS <- 800L
ICOR <- c(0.00, 0.01, 0.05)       # prespecified; 0.01 is the primary
ICOR_PRIMARY <- 0.01
EXEMPLARS <- c("GO:0099536", "GO:0006397", "GO:0006119")

TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)

# ------------------------------------------- the canonical ranked statistic
#
# Rebuilt with the repository's own contract functions on the repository's own
# per-comparison mapped inputs, so the ranking cameraPR sees is the ranking
# gseGO saw.
comparisons <- do.call(rbind, lapply(
  list.dirs(MAPPED, recursive = FALSE), function(ds) {
    d <- file.path(ds, "forward", "per_file")
    if (!dir.exists(d)) return(NULL)
    f <- list.files(d, pattern = "[.]csv$", full.names = TRUE)
    data.frame(dataset = basename(ds), file = f,
               comparison = tools::file_path_sans_ext(basename(f)),
               stringsAsFactors = FALSE)
  }))
cat("mapped comparisons found:", nrow(comparisons), "\n")

`%||%` <- function(a, b) if (is.null(a)) b else a
ranked_of <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  st <- select_rank_statistic(df)
  tr <- protein_group_gene_transform(df, statistic = st, strict = TRUE)
  cg <- collapse_protein_group_genes(tr)
  list(stat = stats::setNames(cg$collapsed_statistic %||% cg$statistic,
                              cg$GeneSymbol),
       rank_column = st$column, fallback = st$fallback_used,
       n_groups = nrow(df), n_genes = nrow(cg))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# one ranked list per comparison, and the per-dataset gene universe
ranked <- list(); meta <- list()
for (i in seq_len(nrow(comparisons))) {
  r <- try(ranked_of(comparisons$file[i]), silent = TRUE)
  if (inherits(r, "try-error")) { cat("  SKIP", comparisons$comparison[i], "\n"); next }
  key <- paste(comparisons$dataset[i], comparisons$comparison[i], sep = "|")
  ranked[[key]] <- r$stat
  meta[[key]] <- data.frame(dataset = comparisons$dataset[i],
                            comparison = comparisons$comparison[i],
                            rank_statistic_column = r$rank_column,
                            rank_statistic_fallback_used = r$fallback,
                            n_protein_groups = r$n_groups,
                            n_ranked_genes = r$n_genes,
                            stringsAsFactors = FALSE)
}
meta <- do.call(rbind, meta)
cat("ranked lists built:", length(ranked), "| any fallback used:",
    any(meta$rank_statistic_fallback_used), "\n")

# -------------------------------------------- the comparable GO-BP family
#
# Built the way clusterProfiler builds it for keyType = "SYMBOL": GOALL
# propagation over the pinned org.Mm.eg.db, then restricted to the measured
# universe and to the SAME minGSSize/maxGSSize contract.
go_terms <- sort(unique(TH$GO_ID))
go2sym <- suppressMessages(AnnotationDbi::mapIds(
  org.Mm.eg.db::org.Mm.eg.db, keys = go_terms, keytype = "GOALL",
  column = "SYMBOL", multiVals = "list"))
go2sym <- lapply(go2sym, function(z) unique(z[!is.na(z)]))
cat("GO-BP terms with an annotation:", sum(lengths(go2sym) > 0), "\n")

go_desc <- stats::setNames(TH$GO_description[!duplicated(TH$GO_ID)],
                           TH$GO_ID[!duplicated(TH$GO_ID)])

# --------------------------------------------------------------- run cameraPR
out <- list()
for (key in names(ranked)) {
  s <- ranked[[key]]
  s <- s[is.finite(s)]
  if (!length(s)) next
  idx <- lapply(go2sym, function(g) which(names(s) %in% g))
  n <- lengths(idx)
  keep <- n >= MIN_GS & n <= MAX_GS
  idx <- idx[keep]
  if (!length(idx)) next
  res <- do.call(rbind, lapply(ICOR, function(rho) {
    cp <- limma::cameraPR(statistic = s, index = idx, use.ranks = FALSE,
                          inter.gene.cor = rho, sort = FALSE)
    if (is.null(cp) || !nrow(cp)) return(NULL)
    ids <- rownames(cp)
    if (is.null(ids)) ids <- names(idx)
    data.frame(GO_ID = ids, inter_gene_cor = rho,
               camera_NGenes = as.integer(cp$NGenes),
               camera_Direction = as.character(cp$Direction),
               camera_PValue = as.numeric(cp$PValue),
               camera_FDR = stats::p.adjust(as.numeric(cp$PValue), method = "BH"),
               stringsAsFactors = FALSE)
  }))
  if (is.null(res) || !nrow(res)) { cat("  no sets for", key, "
"); next }
  res$dataset <- sub("[|].*$", "", key)
  res$comparison <- sub("^.*[|]", "", key)
  out[[key]] <- res
}
cam <- do.call(rbind, out)
rownames(cam) <- NULL
cam$GO_description <- unname(go_desc[cam$GO_ID])
utils::write.csv(cam, file.path(AUD, "camera_full_go_bp_results.csv"),
                 row.names = FALSE)

# ------------------------------------------------------------ S22 the contract
contract <- data.frame(
  item = c("method", "why_not_full_camera", "statistic", "statistic_source",
           "gene_collapse", "gene_sets", "annotation", "size_bounds",
           "inter_gene_cor_primary", "inter_gene_cor_sensitivity",
           "multiple_testing", "n_comparisons", "n_gene_sets_tested_median",
           "role_in_claim_hierarchy", "limma_version"),
  value = c(
    "limma::cameraPR (preranked competitive gene-set test)",
    paste0("no canonical GENE-level abundance matrix exists; the contract ",
           "collapses the STATISTIC, not the abundance, so a gene-level ",
           "matrix would have to be invented for this audit"),
    unique(meta$rank_statistic_column)[1],
    "the same per-comparison mapped inputs the canonical gseGO consumed",
    "median of finite moderated t per official gene SYMBOL (canonical rule)",
    "GO biological process, GOALL propagation, restricted to the measured universe",
    paste0("org.Mm.eg.db ", utils::packageVersion("org.Mm.eg.db"),
           ", GO.db ", utils::packageVersion("GO.db")),
    paste0("minGSSize=", MIN_GS, ", maxGSSize=", MAX_GS, " (canonical)"),
    as.character(ICOR_PRIMARY),
    paste(ICOR, collapse = "; "),
    "Benjamini-Hochberg within each comparison, over the full comparable GO-BP family",
    as.character(length(ranked)),
    as.character(stats::median(table(cam$comparison[
      cam$inter_gene_cor == ICOR_PRIMARY]))),
    "LEVEL 4 orthogonal sensitivity; never a second primary endpoint",
    as.character(utils::packageVersion("limma"))),
  stringsAsFactors = FALSE)
utils::write.csv(contract, file.path(AUD, "camera_contract.csv"),
                 row.names = FALSE)

# ------------------------------------------------------- S23 GSEA concordance
g <- TH[, c("dataset", "spatial_unit", "contrast", "source_comparison", "GO_ID",
            "NES", "GSEA_FDR", "raw_p")]
g <- g[!duplicated(paste(g$source_comparison, g$GO_ID)), , drop = FALSE]
cm <- cam[cam$inter_gene_cor == ICOR_PRIMARY, , drop = FALSE]
j <- merge(g, cm, by.x = c("source_comparison", "GO_ID"),
           by.y = c("comparison", "GO_ID"), suffixes = c("_gsea", "_camera"))
j$gsea_direction <- ifelse(j$NES > 0, "Up", "Down")
j$direction_concordant <- j$gsea_direction == j$camera_Direction
j$gsea_supported <- is.finite(j$GSEA_FDR) & j$GSEA_FDR < 0.05
j$camera_supported <- is.finite(j$camera_FDR) & j$camera_FDR < 0.05
utils::write.csv(j, file.path(AUD, "gsea_camera_concordance.csv"),
                 row.names = FALSE)

conc <- do.call(rbind, lapply(split(j, j$source_comparison), function(z) {
  data.frame(source_comparison = z$source_comparison[1],
             dataset = z$dataset[1], spatial_unit = z$spatial_unit[1],
             contrast = z$contrast[1], n_terms = nrow(z),
             spearman_NES_vs_signed_camera = suppressWarnings(stats::cor(
               z$NES, ifelse(z$camera_Direction == "Up", 1, -1) *
                 (-log10(pmax(z$camera_PValue, .Machine$double.xmin))),
               method = "spearman", use = "complete.obs")),
             directional_concordance = mean(z$direction_concordant),
             n_gsea_supported = sum(z$gsea_supported),
             frac_gsea_supported_dir_concordant = if (any(z$gsea_supported))
               mean(z$direction_concordant[z$gsea_supported]) else NA_real_,
             frac_gsea_supported_also_camera_supported = if (any(z$gsea_supported))
               mean(z$camera_supported[z$gsea_supported]) else NA_real_,
             n_camera_supported = sum(z$camera_supported),
             stringsAsFactors = FALSE)
}))
utils::write.csv(conc, file.path(AUD, "gsea_camera_concordance_by_comparison.csv"),
                 row.names = FALSE)

# ------------------------------------------------------- S24 the three exemplars
ex <- j[j$GO_ID %in% EXEMPLARS, , drop = FALSE]
ex$classification <- with(ex, ifelse(
  !is.finite(camera_PValue), "NOT_EVALUABLE",
  ifelse(direction_concordant & camera_supported, "CONCORDANT_STRONG",
  ifelse(direction_concordant, "CONCORDANT_DIRECTION_ONLY", "DISCORDANT"))))
ex$inter_gene_cor_assumption <- ICOR_PRIMARY
ex <- ex[order(ex$GO_ID, ex$contrast, ex$dataset, ex$spatial_unit), ,
         drop = FALSE]
utils::write.csv(ex, file.path(AUD, "camera_three_exemplar_audit.csv"),
                 row.names = FALSE)

# -------------------------------------------------- S25 per-theme concordance
#
# NO theme-level p-value is created. Each theme is described only by counting
# how its CONSTITUENT GO terms behaved.
elig <- TH[TH$theme_claim_eligible %in% TRUE, c("theme_id", "GO_ID")]
elig <- unique(elig)
thm <- do.call(rbind, lapply(sort(unique(elig$theme_id)), function(t) {
  ids <- elig$GO_ID[elig$theme_id == t]
  z <- j[j$GO_ID %in% ids, , drop = FALSE]
  data.frame(theme_id = t,
             n_constituent_GO_terms = length(ids),
             n_camera_tested = length(unique(z$GO_ID)),
             n_term_comparisons = nrow(z),
             n_direction_concordant = sum(z$direction_concordant),
             fraction_direction_concordant = mean(z$direction_concordant),
             n_camera_FDR_supported = sum(z$camera_supported),
             n_gsea_FDR_supported = sum(z$gsea_supported),
             fraction_gsea_supported_also_camera = if (any(z$gsea_supported))
               mean(z$camera_supported[z$gsea_supported]) else NA_real_,
             representative_exact_GO_terms = paste(utils::head(
               unname(go_desc[sort(unique(z$GO_ID))]), 4), collapse = "; "),
             theme_level_p_value = "NOT COMPUTED - themes are descriptive umbrellas",
             stringsAsFactors = FALSE)
}))
utils::write.csv(thm, file.path(AUD, "camera_theme_concordance.csv"),
                 row.names = FALSE)

cat("\n===== PART-29 CAMERA SENSITIVITY =====\n")
cat("comparisons:", length(ranked), "| GO-BP sets tested per comparison (median):",
    stats::median(table(cm$comparison)), "\n")
cat("rank statistic:", unique(meta$rank_statistic_column),
    "| fallback used anywhere:", any(meta$rank_statistic_fallback_used), "\n")
cat("matched GSEA x CAMERA term-comparisons:", nrow(j), "\n")
cat("directional concordance overall:",
    sprintf("%.3f", mean(j$direction_concordant)), "\n")
cat("GSEA FDR-supported:", sum(j$gsea_supported),
    "| of those direction-concordant:",
    sprintf("%.3f", mean(j$direction_concordant[j$gsea_supported])),
    "| also CAMERA FDR-supported:",
    sprintf("%.3f", mean(j$camera_supported[j$gsea_supported])), "\n")
cat("median per-comparison Spearman:",
    sprintf("%.3f", stats::median(conc$spearman_NES_vs_signed_camera,
                                  na.rm = TRUE)), "\n")
cat("\nthree exemplars:\n")
print(ex[, c("GO_ID", "dataset", "spatial_unit", "contrast", "NES", "GSEA_FDR",
             "camera_Direction", "camera_PValue", "camera_FDR",
             "classification")], row.names = FALSE)
cat("\nper-theme concordance:\n")
print(thm[, c("theme_id", "n_camera_tested", "fraction_direction_concordant",
              "fraction_gsea_supported_also_camera")], row.names = FALSE)
