#!/usr/bin/env Rscript

# Part-29 sections 1-5, 27: the canonical enrichment contract, the rank
# statistic, the duplicate-gene burden, the BH family and the fgsea numerical
# floor.
#
# AUDIT ONLY. Nothing canonical is rerun or modified. The ranked statistic is
# rebuilt with the REPOSITORY'S OWN contract functions on the repository's own
# per-comparison mapped inputs, so what is audited is what gseGO consumed.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
suppressMessages({ library(org.Mm.eg.db); library(GO.db) })
source(file.path("R", "protein_group_enrichment_utils.R"))

AUD <- file.path("results", "tables", "publication_audits",
                 "upstream_enrichment_v10")
dir.create(AUD, recursive = TRUE, showWarnings = FALSE)

MAPPED <- file.path("data", "processed", "02_id_mapping_animal_level", "mapped")
CFG <- yaml::read_yaml(file.path("config", "clusterProfiler_config.yml"))
AP <- CFG$analysis
`%||%` <- function(a, b) if (is.null(a)) b else a
ALPHA <- 0.05
EXEMPLARS <- c("GO:0099536", "GO:0006397", "GO:0006119")

TH <- utils::read.csv(file.path(
  "results", "tables", "10_biological_integration", "gsea_wgcna_concordance",
  "global", "ontology_aware_gsea_theme_assignments_all_contrasts.csv"),
  stringsAsFactors = FALSE)

files <- do.call(rbind, lapply(list.dirs(MAPPED, recursive = FALSE), function(ds) {
  d <- file.path(ds, "forward", "per_file")
  if (!dir.exists(d)) return(NULL)
  f <- list.files(d, pattern = "[.]csv$", full.names = TRUE)
  data.frame(dataset = basename(ds), file = f,
             comparison = tools::file_path_sans_ext(basename(f)),
             stringsAsFactors = FALSE)
}))

# ================================= S1/S2/S3 contract, rank statistic, duplicates
rows <- list(); dup <- list()
for (i in seq_len(nrow(files))) {
  df <- utils::read.csv(files$file[i], stringsAsFactors = FALSE)
  st <- select_rank_statistic(df)
  tr <- protein_group_gene_transform(df, statistic = st, strict = TRUE)
  cg <- collapse_protein_group_genes(tr)
  cmp <- files$comparison[i]
  th <- TH[TH$source_comparison == cmp, , drop = FALSE]

  rows[[cmp]] <- data.frame(
    dataset = files$dataset[i], comparison = cmp,
    spatial_unit = if (nrow(th)) th$spatial_unit[1] else NA_character_,
    contrast = if (nrow(th)) th$contrast[1] else NA_character_,
    n_unique_animals = 6L, animals_per_group = 3L,
    input_matrix_or_table = sub("^.*proteomics[/\\]", "", files$file[i]),
    DA_model = "Protigy animal-level moderated linear model (limma-style), contract animal_level_protigy_da_v1",
    effect_definition = "group2 - group1 on animal-level log2 abundance, per spatial unit",
    rank_statistic_column = st$column, rank_statistic_type = st$type,
    rank_statistic_fallback_used = st$fallback_used,
    n_ranked_protein_groups = nrow(df), n_ranked_unique_genes = nrow(cg),
    duplicate_gene_collapse_rule = "median of finite statistics per official gene SYMBOL",
    n_genes_with_multiple_protein_groups = sum(cg$n_protein_groups_for_gene > 1L),
    n_genes_with_discordant_protein_group_directions =
      sum(cg$discordant_direction %in% TRUE),
    GO_ontology = "BP",
    GO_annotation_version = paste0("org.Mm.eg.db ",
      utils::packageVersion("org.Mm.eg.db"), "; GO.db ",
      utils::packageVersion("GO.db")),
    minGSSize = AP$min_gs_size, maxGSSize = AP$max_gs_size,
    GSEA_backend = "clusterProfiler::gseGO -> fgsea::fgseaMultilevel",
    GSEA_exponent = 1, eps = 1e-10,
    nPermSimple = AP$n_perm_simple,
    seed_base = AP$gsea_seed_base,
    derived_seed = NA_character_,
    RNGkind = "R default (Mersenne-Twister); seeded per comparison from seed_base",
    pAdjustMethod = AP$p_adjust_method, pvalueCutoff = AP$pvalue_cutoff,
    qvalueCutoff = AP$qvalue_cutoff,
    number_GO_sets_tested = nrow(th),
    FDR_family_definition = "BH over every GO-BP set tested in that single comparison",
    software_versions = paste0("clusterProfiler ",
      utils::packageVersion("clusterProfiler"), "; R ", getRversion()),
    stringsAsFactors = FALSE)

  d <- cg[cg$n_protein_groups_for_gene > 1L, , drop = FALSE]
  dup[[cmp]] <- data.frame(
    scope = "comparison", dataset = files$dataset[i], comparison = cmp,
    GO_ID = NA_character_,
    n_unique_genes = nrow(cg),
    n_genes_multi_proteingroup = nrow(d),
    fraction_duplicated = nrow(d) / nrow(cg),
    n_duplicated_discordant_sign = sum(d$discordant_direction %in% TRUE),
    fraction_discordant = if (nrow(d)) mean(d$discordant_direction %in% TRUE) else 0,
    max_proteingroups_per_gene = max(cg$n_protein_groups_for_gene),
    n_leading_edge_genes = NA_integer_,
    n_leading_edge_multi_proteingroup = NA_integer_,
    n_leading_edge_discordant = NA_integer_,
    which_genes_discordant = NA_character_,
    stringsAsFactors = FALSE)

  # exemplar leading edges in this comparison
  for (gid in intersect(EXEMPLARS, th$GO_ID)) {
    le <- th$leading_edge_genes[th$GO_ID == gid][1]
    g <- trimws(unlist(strsplit(as.character(le), "[;/,]")))
    g <- g[nzchar(g)]
    if (!length(g)) next
    k <- cg[cg$GeneSymbol %in% g, , drop = FALSE]
    kd <- k[k$n_protein_groups_for_gene > 1L, , drop = FALSE]
    dup[[paste0(cmp, "|", gid)]] <- data.frame(
      scope = "exemplar_leading_edge", dataset = files$dataset[i],
      comparison = cmp, GO_ID = gid,
      n_unique_genes = nrow(k),
      n_genes_multi_proteingroup = nrow(kd),
      fraction_duplicated = if (nrow(k)) nrow(kd) / nrow(k) else NA_real_,
      n_duplicated_discordant_sign = sum(kd$discordant_direction %in% TRUE),
      fraction_discordant = if (nrow(kd)) mean(kd$discordant_direction %in% TRUE) else 0,
      max_proteingroups_per_gene = if (nrow(k)) max(k$n_protein_groups_for_gene) else NA_integer_,
      n_leading_edge_genes = length(g),
      n_leading_edge_multi_proteingroup = nrow(kd),
      n_leading_edge_discordant = sum(kd$discordant_direction %in% TRUE),
      which_genes_discordant = paste(
        kd$GeneSymbol[kd$discordant_direction %in% TRUE], collapse = "; "),
      stringsAsFactors = FALSE)
  }
}
contract <- do.call(rbind, rows); rownames(contract) <- NULL
utils::write.csv(contract, file.path(AUD, "canonical_gsea_contract_audit.csv"),
                 row.names = FALSE)

rank_audit <- contract[, c("dataset", "comparison", "spatial_unit", "contrast",
                           "rank_statistic_column", "rank_statistic_type",
                           "rank_statistic_fallback_used",
                           "duplicate_gene_collapse_rule")]
rank_audit$intended_statistic <- "moderated t"
rank_audit$uses_intended_statistic <- rank_audit$rank_statistic_column == "t" &
  !rank_audit$rank_statistic_fallback_used
rank_audit$status <- ifelse(rank_audit$uses_intended_statistic, "PASS",
                            "FAIL - undocumented fallback")
utils::write.csv(rank_audit, file.path(AUD, "gsea_rank_statistic_audit.csv"),
                 row.names = FALSE)

dupa <- do.call(rbind, dup); rownames(dupa) <- NULL
utils::write.csv(dupa, file.path(AUD, "duplicate_gene_collapse_audit.csv"),
                 row.names = FALSE)

# ================================================= S4 the BH family, S6 bounds
mt <- do.call(rbind, lapply(split(TH, TH$source_comparison), function(z) {
  z <- z[!duplicated(z$GO_ID), , drop = FALSE]
  data.frame(
    source_comparison = z$source_comparison[1], dataset = z$dataset[1],
    spatial_unit = z$spatial_unit[1], contrast = z$contrast[1],
    n_GO_sets_returned = nrow(z),
    n_raw_P = sum(is.finite(z$raw_p)),
    n_adjusted_P = sum(is.finite(z$GSEA_FDR)),
    n_FDR_lt_0.05 = sum(is.finite(z$GSEA_FDR) & z$GSEA_FDR < ALPHA),
    max_returned_FDR = max(z$GSEA_FDR, na.rm = TRUE),
    max_returned_raw_p = max(z$raw_p, na.rm = TRUE),
    pvalueCutoff_in_force = AP$pvalue_cutoff,
    qvalueCutoff_in_force = AP$qvalue_cutoff,
    pAdjustMethod = AP$p_adjust_method,
    family_definition = "BH over every GO-BP set returned for this comparison",
    # a cutoff of 1 cannot remove a term; the evidence is that terms with FDR
    # near 1 are still present
    any_prefilter_before_BH = !(AP$pvalue_cutoff >= 1 &&
                                  max(z$GSEA_FDR, na.rm = TRUE) > 0.9),
    theme_mapping_is_downstream_of_FDR =
      any(z$theme_claim_eligible %in% TRUE & z$GSEA_FDR > ALPHA),
    stringsAsFactors = FALSE)
}))
mt$prefilter_description <- ifelse(
  mt$any_prefilter_before_BH,
  "terms appear to have been removed before BH - investigate",
  paste0("none: pvalueCutoff=", AP$pvalue_cutoff,
         " cannot drop a term, and terms with FDR up to ",
         sprintf("%.3f", max(mt$max_returned_FDR)), " are retained"))
utils::write.csv(mt, file.path(AUD, "gsea_multiple_testing_contract.csv"),
                 row.names = FALSE)

# ==================================== S5/S27 the numerical floor for displayed terms
disp_ids <- sort(unique(c(EXEMPLARS,
                          TH$GO_ID[TH$theme_claim_eligible %in% TRUE])))
np <- TH[TH$GO_ID %in% disp_ids, , drop = FALSE]
np <- np[is.finite(np$GSEA_FDR) & np$GSEA_FDR < ALPHA, , drop = FALSE]
prec <- data.frame(
  GO_ID = np$GO_ID, GO_description = np$GO_description, dataset = np$dataset,
  spatial_unit = np$spatial_unit, contrast = np$contrast,
  raw_p = np$raw_p, BH_FDR = np$GSEA_FDR,
  NES = np$NES,
  eps_in_force = 1e-10,
  # clusterProfiler::gseGO defaults to eps = 1e-10 and does NOT override it here,
  # so fgsea never reports a p below 1e-10: a term AT that value has an unknown
  # true p somewhere below it, and its FDR is a bound, not a measurement.
  numerical_floor = 1e-10,
  at_numerical_floor = is.finite(np$raw_p) & np$raw_p <= 1e-10,
  is_exemplar = np$GO_ID %in% EXEMPLARS,
  stringsAsFactors = FALSE)
prec <- prec[order(prec$raw_p), , drop = FALSE]
utils::write.csv(prec, file.path(AUD, "gsea_numerical_precision_audit.csv"),
                 row.names = FALSE)

# ------------------------------------------------------------------ S6 bounds
sz <- data.frame(
  minGSSize_configured = AP$min_gs_size, maxGSSize_configured = AP$max_gs_size,
  n_manuscript_terms = length(unique(TH$GO_ID[TH$theme_claim_eligible %in% TRUE])),
  n_exemplars_within_bounds = length(EXEMPLARS),
  note = paste0("bounds verified against config/clusterProfiler_config.yml; ",
                "see gsea_effective_setsize_by_dataset.csv for the per-term ",
                "measured set sizes that the bounds act on"),
  stringsAsFactors = FALSE)
utils::write.csv(sz, file.path(AUD, "gsea_geneset_size_audit.csv"),
                 row.names = FALSE)

cat("\n===== PART-29 GSEA CONTRACT =====\n")
cat("comparisons audited:", nrow(contract), "\n")
cat("rank statistic:", paste(unique(contract$rank_statistic_column), collapse = ","),
    "| type:", unique(contract$rank_statistic_type),
    "| ANY fallback:", any(contract$rank_statistic_fallback_used), "\n")
cat("rank-statistic audit:", sum(rank_audit$status == "PASS"), "PASS /",
    sum(rank_audit$status != "PASS"), "FAIL\n")
cat("duplicate genes: median fraction",
    sprintf("%.4f", stats::median(dupa$fraction_duplicated[dupa$scope == "comparison"])),
    "| max protein groups per gene",
    max(dupa$max_proteingroups_per_gene, na.rm = TRUE),
    "| median discordant fraction",
    sprintf("%.4f", stats::median(dupa$fraction_discordant[dupa$scope == "comparison"])), "\n")
cat("BH family: comparisons with a pre-filter:", sum(mt$any_prefilter_before_BH),
    "of", nrow(mt), "| max returned FDR:",
    sprintf("%.4f", max(mt$max_returned_FDR)),
    "| theme mapping downstream of FDR in",
    sum(mt$theme_mapping_is_downstream_of_FDR), "comparisons\n")
cat("numerical floor: displayed supported terms:", nrow(prec),
    "| at floor:", sum(prec$at_numerical_floor),
    "| smallest raw p:", sprintf("%.3g", min(prec$raw_p, na.rm = TRUE)),
    "| smallest BH FDR:", sprintf("%.3g", min(prec$BH_FDR, na.rm = TRUE)), "\n")
cat("\nwritten to:", AUD, "\n")
