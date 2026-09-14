## Part-29 section 3 - duplicate-gene / protein-group collapse burden audit.
## READ-ONLY audit. Does not rerun, regenerate or modify any canonical analysis output.
## Writes exactly one CSV under results/tables/publication_audits/upstream_enrichment_v10/.
##
## Sources (all canonical, all read-only):
##   data/processed/04_differential_expression_enrichment/clusterProfiler/<dataset>/
##     phenotype_within_unit/<unit>/<comparison>/
##       protein_group_audits/collapsed_gene_input.csv                     (collapse_protein_group_genes output)
##       protein_group_audits/protein_group_to_gene_transformation_audit.csv (protein_group_gene_transform output)
##       protein_group_audits/per_contrast_aggregate_audit.csv             (pipeline's own aggregate counts)
##       GO/BP/GSEA_BP_results_full.csv                                    (canonical ranked gseGO BP results)
##   results/tables/10_biological_integration/gsea_wgcna_concordance/global/
##       ontology_aware_gsea_theme_assignments_all_contrasts.csv           (manuscript theme table)

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")

options(stringsAsFactors = FALSE)

OUT_DIR <- "results/tables/publication_audits/upstream_enrichment_v10"
OUT_CSV <- file.path(OUT_DIR, "duplicate_gene_collapse_audit.csv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CP_ROOT <- "data/processed/04_differential_expression_enrichment/clusterProfiler"
THEME_TABLE <- file.path(
  "results/tables/10_biological_integration/gsea_wgcna_concordance/global",
  "ontology_aware_gsea_theme_assignments_all_contrasts.csv"
)

DATASETS <- c("microglia", "neuron_neuropil", "neuron_soma")

REPO <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics"

## Several canonical per-comparison audit paths exceed the Windows MAX_PATH
## limit of 260 characters (e.g. 262 characters for
## .../CA1microgliares_CA1microgliacon/protein_group_audits/collapsed_gene_input.csv),
## which makes plain file.exists()/read.csv() fail silently. Route every read
## through the Windows extended-length prefix. Read-only; nothing is written here.
lp <- function(p) {
  abs <- ifelse(grepl("^([A-Za-z]:|//|\\\\\\\\)", p), p, file.path(REPO, p))
  paste0("\\\\?\\", gsub("/", "\\\\", abs))
}

## ---------------------------------------------------------------- helpers ---

contrast_from_folder <- function(folder) {
  ## folder looks like "<unit><pheno>_<unit><pheno>", e.g. CA3srsus_CA3srres
  parts <- strsplit(folder, "_", fixed = TRUE)[[1]]
  left <- tolower(substr(parts[[1]], nchar(parts[[1]]) - 2, nchar(parts[[1]])))
  right <- tolower(substr(parts[[2]], nchar(parts[[2]]) - 2, nchar(parts[[2]])))
  paste(toupper(left), "-", toupper(right))
}

split_semi <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  v <- trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
  v[nzchar(v)]
}

split_slash <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  v <- trimws(unlist(strsplit(x, "/", fixed = TRUE), use.names = FALSE))
  v[nzchar(v)]
}

fmt <- function(x, digits = 10) {
  if (length(x) == 0) return(NA_character_)
  paste(format(x, trim = TRUE, digits = digits, scientific = FALSE), collapse = ";")
}

## Enumerate every comparison directory that carries the canonical audits.
comparison_index <- function() {
  rows <- list()
  for (ds in DATASETS) {
    unit_root <- file.path(CP_ROOT, ds, "phenotype_within_unit")
    units <- sort(list.dirs(unit_root, recursive = FALSE, full.names = FALSE))
    for (u in units) {
      comps <- sort(list.dirs(file.path(unit_root, u), recursive = FALSE, full.names = FALSE))
      for (cmp in comps) {
        base <- file.path(unit_root, u, cmp)
        rows[[length(rows) + 1L]] <- data.frame(
          dataset = ds,
          spatial_unit = u,
          source_comparison = cmp,
          contrast = contrast_from_folder(cmp),
          collapse_csv = file.path(base, "protein_group_audits", "collapsed_gene_input.csv"),
          transform_csv = file.path(base, "protein_group_audits",
            "protein_group_to_gene_transformation_audit.csv"),
          aggregate_csv = file.path(base, "protein_group_audits", "per_contrast_aggregate_audit.csv"),
          gsea_csv = file.path(base, "GO", "BP", "GSEA_BP_results_full.csv")
        )
      }
    }
  }
  do.call(rbind, rows)
}

idx <- comparison_index()
cat("comparisons found:", nrow(idx), "\n")
stopifnot(all(file.exists(lp(idx$collapse_csv))))
stopifnot(all(file.exists(lp(idx$gsea_csv))))
stopifnot(all(file.exists(lp(idx$transform_csv))))
stopifnot(all(file.exists(lp(idx$aggregate_csv))))

## ------------------------------------------------------------ row template ---

empty_row <- function() {
  data.frame(
    row_type = NA_character_, dataset = NA_character_, spatial_unit = NA_character_,
    contrast = NA_character_, source_comparison = NA_character_,
    go_id = NA_character_, go_description = NA_character_,
    n_unique_genes = NA_integer_, n_genes_multi_proteingroup = NA_integer_,
    fraction_duplicated = NA_real_, n_duplicated_discordant_sign = NA_integer_,
    fraction_discordant = NA_real_, max_proteingroups_per_gene = NA_integer_,
    fraction_discordant_of_all_genes = NA_real_,
    n_eligible_protein_groups = NA_integer_, n_total_protein_groups = NA_integer_,
    n_excluded_protein_groups = NA_integer_,
    n_excluded_gene_level_claim_not_allowed = NA_integer_,
    n_excluded_multi_gene_indistinguishable = NA_integer_,
    n_excluded_unresolved_group = NA_integer_,
    n_excluded_mixed_species_or_contaminant = NA_integer_,
    n_excluded_other = NA_integer_,
    aggregate_audit_enumerated_exclusions = NA_integer_,
    aggregate_audit_exclusion_enumeration_complete = NA_character_,
    multi_proteingroup_genes = NA_character_,
    multi_proteingroup_member_statistics = NA_character_,
    multi_proteingroup_collapsed_statistic = NA_character_,
    multi_proteingroup_direction_pattern = NA_character_,
    multi_proteingroup_contributing_ids = NA_character_,
    rank_statistic_column = NA_character_, rank_statistic_fallback_used = NA_character_,
    gene_collapse_rule = NA_character_,
    n_bp_terms_total = NA_integer_,
    n_bp_terms_with_multi_pg_gene_in_leading_edge = NA_integer_,
    n_leading_edge_genes = NA_integer_, n_leading_edge_multi_proteingroup = NA_integer_,
    n_leading_edge_discordant = NA_integer_, which_genes_discordant = NA_character_,
    term_NES = NA_real_, term_enrichment_score = NA_real_,
    term_pvalue = NA_real_, term_p_adjust = NA_real_,
    term_set_size = NA_integer_, term_leading_edge_tag = NA_character_,
    leading_edge_side = NA_character_,
    leading_edge_boundary_gene = NA_character_,
    leading_edge_boundary_statistic = NA_real_,
    leading_edge_min_abs_margin_to_boundary = NA_real_,
    multi_pg_le_genes = NA_character_,
    multi_pg_le_member_statistics = NA_character_,
    multi_pg_le_collapsed_statistic = NA_character_,
    multi_pg_le_alt_member_worst_case_statistic = NA_character_,
    n_le_genes_flipped_by_alternative_member = NA_integer_,
    genes_flipped_by_alternative_member = NA_character_,
    collapse_materially_determines_leading_edge = NA_character_,
    theme_table_le_gene_count = NA_integer_,
    theme_table_le_matches_canonical_core_enrichment = NA_character_,
    theme_table_manuscript_themes = NA_character_,
    theme_table_claim_eligible = NA_character_,
    theme_table_assignment_status = NA_character_,
    notes = NA_character_,
    stringsAsFactors = FALSE
  )
}

## ----------------------------------------------- per-comparison summary rows ---

summary_rows <- list()
collapse_cache <- list()
gsea_cache <- list()
multi_gene_cache <- list()

for (i in seq_len(nrow(idx))) {
  info <- idx[i, ]
  cg <- utils::read.csv(lp(info$collapse_csv))
  gs <- utils::read.csv(lp(info$gsea_csv))
  key <- paste(info$dataset, info$source_comparison, sep = "|")
  collapse_cache[[key]] <- cg
  gsea_cache[[key]] <- gs

  npg <- as.integer(cg$n_protein_groups_for_gene)
  disc <- as.logical(cg$discordant_direction)
  n_unique <- nrow(cg)
  n_multi <- sum(npg > 1L)
  n_disc <- sum(disc %in% TRUE)
  multi_sel <- which(npg > 1L)
  multi_genes <- cg$GeneSymbol[multi_sel]
  multi_gene_cache[[key]] <- multi_genes

  ## how many GSEA BP terms carry a multi-protein-group gene in their leading edge
  n_terms <- nrow(gs)
  if (length(multi_genes)) {
    le_lists <- lapply(gs$core_enrichment, split_slash)
    hit <- vapply(le_lists, function(v) any(multi_genes %in% v), logical(1))
    n_terms_hit <- sum(hit)
  } else {
    n_terms_hit <- 0L
  }

  agg <- utils::read.csv(lp(info$aggregate_csv))
  tr <- utils::read.csv(lp(info$transform_csv))
  transform_n <- nrow(tr)
  excl <- table(tr$exclusion_reason[tr$eligibility_status == "excluded"])
  getn <- function(nm) if (nm %in% names(excl)) as.integer(excl[[nm]]) else 0L
  n_excl_claim <- getn("gene_level_claim_not_allowed")
  n_excl_multi <- getn("excluded_multi_gene_indistinguishable")
  n_excl_unres <- getn("excluded_unresolved_group")
  n_excl_contam <- getn("excluded_mixed_species_or_contaminant")
  n_excl_total <- sum(tr$eligibility_status == "excluded")
  agg_enum <- as.integer(agg$excluded_multi_gene_groups[[1]]) +
    as.integer(agg$excluded_partially_mapped_groups[[1]]) +
    as.integer(agg$excluded_unresolved_groups[[1]]) +
    as.integer(agg$excluded_mixed_species_or_contaminant_groups[[1]])

  r <- empty_row()
  r$row_type <- "comparison_summary"
  r$dataset <- info$dataset
  r$spatial_unit <- info$spatial_unit
  r$contrast <- info$contrast
  r$source_comparison <- info$source_comparison
  r$n_unique_genes <- n_unique
  r$n_genes_multi_proteingroup <- n_multi
  r$fraction_duplicated <- n_multi / n_unique
  r$n_duplicated_discordant_sign <- n_disc
  r$fraction_discordant <- if (n_multi > 0) n_disc / n_multi else NA_real_
  r$fraction_discordant_of_all_genes <- n_disc / n_unique
  r$max_proteingroups_per_gene <- max(npg)
  r$n_eligible_protein_groups <- sum(npg)
  r$n_total_protein_groups <- as.integer(agg$total_ProteinGroupIDs[[1]])
  r$n_excluded_protein_groups <- n_excl_total
  r$n_excluded_gene_level_claim_not_allowed <- n_excl_claim
  r$n_excluded_multi_gene_indistinguishable <- n_excl_multi
  r$n_excluded_unresolved_group <- n_excl_unres
  r$n_excluded_mixed_species_or_contaminant <- n_excl_contam
  r$n_excluded_other <- n_excl_total - (n_excl_claim + n_excl_multi + n_excl_unres + n_excl_contam)
  r$aggregate_audit_enumerated_exclusions <- agg_enum
  r$aggregate_audit_exclusion_enumeration_complete <-
    if (identical(agg_enum, n_excl_total)) "TRUE" else
      sprintf("FALSE_%d_of_%d_excluded_protein_groups_not_enumerated",
        n_excl_total - agg_enum, n_excl_total)
  r$multi_proteingroup_genes <- paste(multi_genes, collapse = ";")
  r$multi_proteingroup_member_statistics <-
    paste(cg$individual_statistics[multi_sel], collapse = "|")
  r$multi_proteingroup_collapsed_statistic <-
    paste(fmt(cg$collapsed_statistic[multi_sel]), collapse = "|")
  r$multi_proteingroup_direction_pattern <-
    paste(cg$direction_pattern[multi_sel], collapse = "|")
  r$multi_proteingroup_contributing_ids <-
    paste(cg$contributing_ProteinGroupIDs[multi_sel], collapse = "|")
  r$rank_statistic_column <- as.character(agg$rank_statistic_column[[1]])
  r$rank_statistic_fallback_used <- as.character(agg$rank_statistic_fallback_used[[1]])
  r$gene_collapse_rule <- paste(sort(unique(as.character(cg$gene_collapse_rule))), collapse = ";")
  r$n_bp_terms_total <- n_terms
  r$n_bp_terms_with_multi_pg_gene_in_leading_edge <- n_terms_hit

  ## consistency check against the pipeline's own aggregate audit
  chk <- c(
    sprintf("aggregate_genes_after_collapse=%d", as.integer(agg$genes_after_duplicate_collapse[[1]])),
    sprintf("aggregate_genes_multi_PG=%d", as.integer(agg$genes_with_multiple_ProteinGroupIDs[[1]])),
    sprintf("aggregate_genes_discordant=%d", as.integer(agg$genes_with_discordant_directions[[1]])),
    sprintf("transform_rows=%d", transform_n),
    sprintf("recomputed_matches_aggregate=%s",
      identical(n_unique, as.integer(agg$genes_after_duplicate_collapse[[1]])) &&
        identical(n_multi, as.integer(agg$genes_with_multiple_ProteinGroupIDs[[1]])) &&
        identical(n_disc, as.integer(agg$genes_with_discordant_directions[[1]])))
  )
  r$notes <- paste(chk, collapse = "; ")
  summary_rows[[length(summary_rows) + 1L]] <- r
}

summary_df <- do.call(rbind, summary_rows)

## --------------------------------------- leading-edge row builder + sensitivity ---
##
## Descriptive sensitivity (NOT a GSEA rerun): the canonical leading edge is the
## contiguous block of gene-set members at one end of the ranked list. Its inner
## boundary is the set member furthest from that end (largest rank position for
## ES >= 0, smallest for ES < 0). A leading-edge gene that maps to several protein
## groups is held to be "collapse-determined" if substituting any single member
## protein group's own statistic for the median would move it past that boundary
## statistic, i.e. out of the leading-edge block.

build_le_row <- function(info, cg, gs, go_id, row_type) {
  r <- empty_row()
  r$row_type <- row_type
  r$dataset <- info$dataset
  r$spatial_unit <- info$spatial_unit
  r$contrast <- info$contrast
  r$source_comparison <- info$source_comparison
  r$go_id <- go_id

  term <- gs[gs$ID == go_id, , drop = FALSE]
  if (nrow(term) != 1L) {
    r$notes <- "GO term absent from canonical GSEA_BP_results_full.csv for this comparison"
    return(r)
  }

  le_genes <- split_slash(term$core_enrichment[[1]])
  r$go_description <- as.character(term$Description[[1]])
  r$term_NES <- as.numeric(term$NES[[1]])
  r$term_enrichment_score <- as.numeric(term$enrichmentScore[[1]])
  r$term_pvalue <- as.numeric(term$pvalue[[1]])
  r$term_p_adjust <- as.numeric(term$p.adjust[[1]])
  r$term_set_size <- as.integer(term$setSize[[1]])
  r$term_leading_edge_tag <- as.character(term$leading_edge[[1]])
  r$n_leading_edge_genes <- length(le_genes)

  m <- match(le_genes, cg$GeneSymbol)
  stopifnot(!anyNA(m))
  npg_le <- as.integer(cg$n_protein_groups_for_gene[m])
  disc_le <- as.logical(cg$discordant_direction[m])
  r$n_leading_edge_multi_proteingroup <- sum(npg_le > 1L)
  r$n_leading_edge_discordant <- sum(disc_le %in% TRUE)
  disc_names <- le_genes[disc_le %in% TRUE]
  r$which_genes_discordant <- paste(disc_names, collapse = ";")

  ## canonical ranked vector, rebuilt exactly as build_enrichment_gene_inputs() does
  ranked <- stats::setNames(cg$collapsed_statistic, cg$GeneSymbol)
  ranked <- ranked[order(cg$collapsed_statistic, decreasing = TRUE)]
  pos <- match(le_genes, names(ranked))
  stopifnot(!anyNA(pos))
  es <- as.numeric(term$enrichmentScore[[1]])
  r$leading_edge_side <- if (es >= 0) "top_of_ranked_list" else "bottom_of_ranked_list"
  boundary_pos <- if (es >= 0) max(pos) else min(pos)
  boundary_gene <- names(ranked)[boundary_pos]
  boundary_stat <- unname(ranked[boundary_pos])
  r$leading_edge_boundary_gene <- boundary_gene
  r$leading_edge_boundary_statistic <- boundary_stat
  ## closest other leading-edge member to the inner boundary (the boundary gene
  ## itself is excluded, otherwise this is trivially zero)
  inner <- pos[pos != boundary_pos]
  r$leading_edge_min_abs_margin_to_boundary <-
    if (length(inner)) min(abs(unname(ranked[inner]) - boundary_stat)) else NA_real_

  multi_idx <- which(npg_le > 1L)
  if (!length(multi_idx)) {
    r$multi_pg_le_genes <- ""
    r$multi_pg_le_member_statistics <- ""
    r$multi_pg_le_collapsed_statistic <- ""
    r$multi_pg_le_alt_member_worst_case_statistic <- ""
    r$n_le_genes_flipped_by_alternative_member <- 0L
    r$genes_flipped_by_alternative_member <- ""
    r$collapse_materially_determines_leading_edge <-
      "NO_no_leading_edge_gene_maps_to_more_than_one_protein_group"
  } else {
    g_names <- character(0); g_stats <- character(0)
    g_med <- character(0); g_worst <- character(0)
    flipped <- character(0)
    for (k in multi_idx) {
      g <- le_genes[[k]]
      ci <- m[[k]]
      members <- suppressWarnings(as.numeric(split_semi(cg$individual_statistics[[ci]])))
      members <- members[is.finite(members)]
      med <- as.numeric(cg$collapsed_statistic[[ci]])
      ## if the gene itself defines the boundary, use the next leading-edge gene inwards
      b_stat <- boundary_stat
      if (identical(g, boundary_gene)) {
        other <- pos[pos != boundary_pos]
        b_stat <- if (es >= 0) unname(ranked[max(other)]) else unname(ranked[min(other)])
      }
      stays <- if (es >= 0) members >= b_stat else members <= b_stat
      worst <- if (es >= 0) min(members) else max(members)
      g_names <- c(g_names, g)
      g_stats <- c(g_stats, fmt(members))
      g_med <- c(g_med, fmt(med))
      g_worst <- c(g_worst, fmt(worst))
      if (any(!stays)) flipped <- c(flipped, g)
    }
    r$multi_pg_le_genes <- paste(g_names, collapse = ";")
    r$multi_pg_le_member_statistics <- paste(g_stats, collapse = "|")
    r$multi_pg_le_collapsed_statistic <- paste(g_med, collapse = "|")
    r$multi_pg_le_alt_member_worst_case_statistic <- paste(g_worst, collapse = "|")
    r$n_le_genes_flipped_by_alternative_member <- length(flipped)
    r$genes_flipped_by_alternative_member <- paste(flipped, collapse = ";")
    r$collapse_materially_determines_leading_edge <- if (length(flipped)) {
      sprintf("YES_for_%d_of_%d_multi_protein_group_leading_edge_genes",
        length(flipped), length(multi_idx))
    } else {
      sprintf("NO_all_%d_multi_protein_group_leading_edge_genes_stay_under_every_member_statistic",
        length(multi_idx))
    }
  }

  ## carry the comparison-level burden onto the row for context
  srow <- summary_df[summary_df$dataset == info$dataset &
      summary_df$source_comparison == info$source_comparison, , drop = FALSE]
  carry <- c("n_unique_genes", "n_genes_multi_proteingroup", "fraction_duplicated",
    "n_duplicated_discordant_sign", "fraction_discordant", "fraction_discordant_of_all_genes",
    "max_proteingroups_per_gene", "multi_proteingroup_genes",
    "multi_proteingroup_member_statistics", "multi_proteingroup_collapsed_statistic",
    "multi_proteingroup_direction_pattern", "multi_proteingroup_contributing_ids",
    "rank_statistic_column", "rank_statistic_fallback_used", "gene_collapse_rule")
  for (cn in carry) r[[cn]] <- srow[[cn]][[1]]
  r$notes <- sprintf(
    "leading edge from canonical core_enrichment; ranked list n=%d; boundary position=%d; descriptive sensitivity, no GSEA rerun",
    length(ranked), boundary_pos)
  r
}

## ------------------------------------------------------- exemplar GO terms ---

EXEMPLARS <- data.frame(
  dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
  spatial_unit = c("CA3_sr", "CA2_sp", "CA1_microglia"),
  go_id = c("GO:0099536", "GO:0006397", "GO:0006119"),
  stringsAsFactors = FALSE
)

exemplar_rows <- list()
for (e in seq_len(nrow(EXEMPLARS))) {
  ex <- EXEMPLARS[e, ]
  sub <- idx[idx$dataset == ex$dataset & idx$spatial_unit == ex$spatial_unit, , drop = FALSE]
  stopifnot(nrow(sub) == 3L)
  for (j in seq_len(nrow(sub))) {
    info <- sub[j, ]
    key <- paste(info$dataset, info$source_comparison, sep = "|")
    exemplar_rows[[length(exemplar_rows) + 1L]] <-
      build_le_row(info, collapse_cache[[key]], gsea_cache[[key]], ex$go_id,
        "exemplar_leading_edge")
  }
}
exemplar_df <- do.call(rbind, exemplar_rows)

## ------- exhaustive sweep: every BP term whose leading edge contains a multi-PG gene ---

sweep_rows <- list()
for (i in seq_len(nrow(idx))) {
  info <- idx[i, ]
  key <- paste(info$dataset, info$source_comparison, sep = "|")
  multi_genes <- multi_gene_cache[[key]]
  if (!length(multi_genes)) next
  gs <- gsea_cache[[key]]
  cg <- collapse_cache[[key]]
  le_lists <- lapply(gs$core_enrichment, split_slash)
  hit <- which(vapply(le_lists, function(v) any(multi_genes %in% v), logical(1)))
  for (h in hit) {
    sweep_rows[[length(sweep_rows) + 1L]] <-
      build_le_row(info, cg, gs, as.character(gs$ID[[h]]),
        "multi_pg_gene_leading_edge_sensitivity")
  }
}
sweep_df <- if (length(sweep_rows)) do.call(rbind, sweep_rows) else exemplar_df[0, ]

## --------------------------------- cross-check against manuscript theme table ---

all_multi_genes <- sort(unique(unlist(multi_gene_cache, use.names = FALSE)))
theme_scan <- local({
  pat <- paste(c(EXEMPLARS$go_id, all_multi_genes), collapse = "|")
  con <- file(lp(THEME_TABLE), "r")
  on.exit(close(con))
  hdr <- readLines(con, 1L)
  keep <- character(0)
  key_chunks <- list()
  n_lines <- 0L
  repeat {
    x <- readLines(con, 50000L)
    if (!length(x)) break
    n_lines <- n_lines + length(x)
    ## fields 1-6 (dataset, phenotype_contrast, contrast, spatial_unit,
    ## source_comparison, GO_ID) are comma-free, so a positional cut is safe.
    key_chunks[[length(key_chunks) + 1L]] <- sub("^(([^,]*,){5}[^,]*).*$", "\\1", x)
    k <- grepl(pat, x, fixed = FALSE)
    if (any(k)) keep <- c(keep, x[k])
  }
  list(n_lines = n_lines, keys = unlist(key_chunks, use.names = FALSE),
    df = if (length(keep)) utils::read.csv(text = paste(c(hdr, keep), collapse = "\n")) else NULL)
})
cat("theme table data rows scanned:", theme_scan$n_lines, "\n")
theme_hits <- theme_scan$df

## reconcile theme-table row count against the canonical GSEA BP term rows
local({
  keys <- theme_scan$keys
  parts <- do.call(rbind, strsplit(keys, ",", fixed = TRUE))
  tt_ds <- gsub('"', "", parts[, 1]); tt_cmp <- gsub('"', "", parts[, 5])
  tt_go <- gsub('"', "", parts[, 6])
  triple <- paste(tt_ds, tt_cmp, tt_go, sep = "|")
  cat("theme-table rows:", length(keys),
    "| distinct dataset+comparison+GO_ID triples:", length(unique(triple)),
    "| duplicated rows (term in >1 theme):", sum(duplicated(triple)), "\n")
  tt_counts <- table(paste(tt_ds, tt_cmp, sep = "|"))
  can_counts <- stats::setNames(summary_df$n_bp_terms_total,
    paste(summary_df$dataset, summary_df$source_comparison, sep = "|"))
  uniq_counts <- tapply(triple, paste(tt_ds, tt_cmp, sep = "|"),
    function(v) length(unique(v)))
  common <- intersect(names(can_counts), names(uniq_counts))
  cat("comparisons where distinct theme-table GO_IDs == canonical BP term rows:",
    sum(uniq_counts[common] == can_counts[common]), "of", length(common), "\n")
  cat("canonical BP term rows total:", sum(can_counts),
    "| distinct theme-table triples total:", length(unique(triple)), "\n")
})

unit_key_for_theme <- function(dataset, spatial_unit) {
  if (identical(dataset, "microglia")) sub("_microglia$", "", spatial_unit) else spatial_unit
}

annotate_theme <- function(df) {
  if (is.null(theme_hits) || !nrow(df)) return(df)
  for (i in seq_len(nrow(df))) {
    if (is.na(df$go_id[[i]])) next
    uk <- unit_key_for_theme(df$dataset[[i]], df$spatial_unit[[i]])
    th <- theme_hits[theme_hits$dataset == df$dataset[[i]] &
        theme_hits$spatial_unit == uk &
        theme_hits$source_comparison == df$source_comparison[[i]] &
        theme_hits$GO_ID == df$go_id[[i]], , drop = FALSE]
    if (!nrow(th)) {
      df$theme_table_le_matches_canonical_core_enrichment[[i]] <- "absent_from_theme_table"
      next
    }
    le_theme <- sort(unique(split_semi(th$leading_edge_genes[[1]])))
    df$theme_table_le_gene_count[[i]] <- length(le_theme)
    key <- paste(df$dataset[[i]], df$source_comparison[[i]], sep = "|")
    gs <- gsea_cache[[key]]
    le_can <- sort(unique(split_slash(gs$core_enrichment[gs$ID == df$go_id[[i]]][[1]])))
    df$theme_table_le_matches_canonical_core_enrichment[[i]] <-
      if (identical(le_theme, le_can)) "TRUE" else
        sprintf("FALSE_theme_only=%d_canonical_only=%d",
          length(setdiff(le_theme, le_can)), length(setdiff(le_can, le_theme)))
    df$theme_table_manuscript_themes[[i]] <-
      paste(sort(unique(as.character(th$manuscript_theme))), collapse = ";")
    df$theme_table_claim_eligible[[i]] <-
      paste(sort(unique(as.character(th$theme_claim_eligible))), collapse = ";")
    df$theme_table_assignment_status[[i]] <-
      paste(sort(unique(as.character(th$assignment_status))), collapse = ";")
  }
  df
}

exemplar_df <- annotate_theme(exemplar_df)
sweep_df <- annotate_theme(sweep_df)

## ------------------------------------------------------------------ write ---

out <- rbind(summary_df, exemplar_df, sweep_df)
utils::write.csv(out, lp(OUT_CSV), row.names = FALSE, na = "")
cat("wrote:", OUT_CSV, " rows:", nrow(out), "\n")

## ------------------------------------------------------------- console log ---

cat("\n== comparison-level duplicate burden (54 comparisons) ==\n")
print(summary_df[, c("dataset", "spatial_unit", "contrast", "n_unique_genes",
  "n_genes_multi_proteingroup", "fraction_duplicated", "n_duplicated_discordant_sign",
  "fraction_discordant", "max_proteingroups_per_gene", "multi_proteingroup_genes",
  "multi_proteingroup_direction_pattern",
  "n_bp_terms_with_multi_pg_gene_in_leading_edge", "n_bp_terms_total")], row.names = FALSE)

cat("\n-- global duplicate burden --\n")
cat("recomputed counts match pipeline aggregate audit in every comparison:",
  all(grepl("recomputed_matches_aggregate=TRUE", summary_df$notes, fixed = TRUE)), "\n")
cat("n_unique_genes: min", min(summary_df$n_unique_genes), "max", max(summary_df$n_unique_genes), "\n")
cat("n_genes_multi_proteingroup: min", min(summary_df$n_genes_multi_proteingroup),
  "max", max(summary_df$n_genes_multi_proteingroup),
  "sum", sum(summary_df$n_genes_multi_proteingroup), "\n")
cat("fraction_duplicated: min", min(summary_df$fraction_duplicated),
  "max", max(summary_df$fraction_duplicated), "\n")
cat("max_proteingroups_per_gene across all comparisons:",
  max(summary_df$max_proteingroups_per_gene), "\n")
cat("comparisons with a discordant duplicated gene:",
  sum(summary_df$n_duplicated_discordant_sign > 0), "of", nrow(summary_df), "\n")
cat("distinct multi-protein-group genes anywhere:", paste(all_multi_genes, collapse = ";"), "\n")
cat("total BP term rows across comparisons:", sum(summary_df$n_bp_terms_total), "\n")
cat("BP term rows whose leading edge contains a multi-PG gene:",
  sum(summary_df$n_bp_terms_with_multi_pg_gene_in_leading_edge), "\n")

cat("\n-- protein groups per comparison --\n")
print(unique(summary_df[, c("dataset", "n_total_protein_groups", "n_eligible_protein_groups",
  "n_excluded_protein_groups", "n_unique_genes")]), row.names = FALSE)

cat("\n-- protein-group exclusion breakdown (from the canonical transformation audit) --\n")
print(unique(summary_df[, c("dataset", "n_excluded_protein_groups",
  "n_excluded_gene_level_claim_not_allowed", "n_excluded_multi_gene_indistinguishable",
  "n_excluded_unresolved_group", "n_excluded_mixed_species_or_contaminant",
  "n_excluded_other", "aggregate_audit_enumerated_exclusions",
  "aggregate_audit_exclusion_enumeration_complete")]), row.names = FALSE)
cat("comparisons where per_contrast_aggregate_audit.csv enumerates every excluded protein group:",
  sum(summary_df$aggregate_audit_exclusion_enumeration_complete == "TRUE"), "of",
  nrow(summary_df), "\n")

cat("\n== exemplar leading edges ==\n")
print(exemplar_df[, c("dataset", "spatial_unit", "contrast", "go_id", "go_description",
  "term_NES", "term_p_adjust", "term_set_size", "n_leading_edge_genes",
  "n_leading_edge_multi_proteingroup", "n_leading_edge_discordant",
  "which_genes_discordant", "collapse_materially_determines_leading_edge")], row.names = FALSE)
cat("\nexemplar leading-edge cross-check vs manuscript theme table:\n")
print(exemplar_df[, c("go_id", "contrast", "n_leading_edge_genes",
  "theme_table_le_gene_count", "theme_table_le_matches_canonical_core_enrichment")],
  row.names = FALSE)

cat("\n== sweep: every BP term with a multi-PG gene in its leading edge ==\n")
cat("rows:", nrow(sweep_df), "\n")
if (nrow(sweep_df)) {
  cat("terms where an alternative member protein-group statistic would remove the gene",
    "from the leading edge:", sum(sweep_df$n_le_genes_flipped_by_alternative_member > 0), "\n")
  print(table(sweep_df$collapse_materially_determines_leading_edge))
  cat("\naffected comparisons:\n")
  print(table(paste(sweep_df$dataset, sweep_df$spatial_unit, sweep_df$contrast, sep = " / ")))
  cat("\ndistinct GO terms in the sweep:", length(unique(sweep_df$go_id)), "\n")
  cat("GO terms:\n")
  print(unique(sweep_df[, c("go_id", "go_description")]), row.names = FALSE)

  cat("\nmanuscript-theme assignment of the sweep rows:\n")
  print(table(sweep_df$theme_table_manuscript_themes,
    sweep_df$theme_table_claim_eligible, useNA = "ifany"))

  fl <- sweep_df[sweep_df$n_le_genes_flipped_by_alternative_member > 0, , drop = FALSE]
  cat("\nflipped rows:", nrow(fl), "\n")
  if (nrow(fl)) {
    cat("flipped rows assigned to a claim-eligible manuscript theme:",
      sum(grepl("TRUE", fl$theme_table_claim_eligible, fixed = TRUE)), "\n")
    cat("FDR<0.05 among flipped rows:", sum(fl$term_p_adjust < 0.05), "\n")
    cat("\nworked example (first flipped row):\n")
    print(t(fl[1, c("dataset", "spatial_unit", "contrast", "go_id", "go_description",
      "term_NES", "term_enrichment_score", "term_p_adjust", "n_leading_edge_genes",
      "leading_edge_side", "leading_edge_boundary_gene", "leading_edge_boundary_statistic",
      "multi_pg_le_genes", "multi_pg_le_member_statistics", "multi_pg_le_collapsed_statistic",
      "multi_pg_le_alt_member_worst_case_statistic",
      "theme_table_manuscript_themes", "theme_table_claim_eligible")]))
  }

  cat("\nexemplar GO terms present in the sweep (should be none):",
    sum(sweep_df$go_id %in% EXEMPLARS$go_id), "\n")
}

cat("\n", R.version.string, "\n")
