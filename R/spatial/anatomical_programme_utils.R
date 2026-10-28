# Anatomical GO-BP programmes of the control hippocampus: the families of
# sibling units, the contrasts added beyond the 11 manuscript-locked ones, and
# the presentation-only rule that picks the unit-specific programmes shown in
# Figure 2. Used by analysis/spatial_validation/map_anatomical_programmes.R.
#
# The locked contrast definitions stay in control_spatial_contrast_registry()
# (R/spatial/control_spatial_identity_utils.R) and are not changed here. The
# added contrasts are defined below, once, with manuscript_locked = FALSE.

# One row per unit: which contrast tests it against its siblings, and where
# that contrast comes from ("locked" = validate_control_spatial_identity.R,
# "added" = map_anatomical_programmes.R).
anatomical_programme_families <- function() {
  f <- function(family, dataset, source, units, contrast_fmt, label, prefix = "", shown = FALSE) {
    data.frame(family = family, dataset = dataset, source = source,
               unit = units, anatomical_unit = paste0(prefix, units),
               contrast = sprintf(contrast_fmt, units),
               family_label = label, unit_order = seq_along(units),
               shown_in_figure = shown, stringsAsFactors = FALSE)
  }
  rbind(
    f("soma_regions", "neuron_soma", "locked", c("CA1", "CA2", "CA3", "DG"),
      "%s_vs_mean_other_soma_regions", "Neuronal soma, regions", shown = TRUE),
    f("CA1_strata", "neuron_neuropil", "locked", c("SO", "SR", "SLM"),
      "CA1_%s_vs_mean_other_CA1_strata", "CA1 neuropil, strata", "CA1_", TRUE),
    f("CA2_strata", "neuron_neuropil", "added", c("SO", "SR", "SLM"),
      "CA2_%s_vs_mean_other_CA2_strata", "CA2 neuropil, strata", "CA2_", TRUE),
    f("CA3_strata", "neuron_neuropil", "added", c("SO", "SR"),
      "CA3_%s_vs_mean_other_CA3_strata", "CA3 neuropil, strata", "CA3_"),
    f("DG_layers", "neuron_neuropil", "locked", c("MO", "PO"),
      "DG_%s_vs_mean_other_DG_layers", "DG neuropil, layers", "DG_", TRUE),
    f("microglia_regions", "microglia", "added", c("CA1", "CA2", "CA3", "DG"),
      "%s_vs_mean_other_microglia_regions", "Microglia-enriched ROI, regions")
  )
}

# Direct comparisons between sibling units, so that a programme is called
# specific to a unit only when it is enriched there against EACH sibling
# (never from "significant here, not significant there"). One spec per
# unordered pair (first unit +1, second unit -1); a two-unit family's
# one-vs-rest contrast already is that pair and is reused, not refitted.
anatomical_programme_pairwise_name <- function(family, unit_a, unit_b) {
  paste0("pairwise_", family, "_", unit_a, "_vs_", unit_b)
}

anatomical_programme_pairs <- function(families = anatomical_programme_families()) {
  out <- list()
  for (fam in unique(families$family)) {
    x <- families[families$family == fam, , drop = FALSE]
    x <- x[order(x$unit_order), , drop = FALSE]
    if (nrow(x) < 2L) next
    for (i in seq_len(nrow(x) - 1L)) for (j in (i + 1L):nrow(x)) {
      two <- nrow(x) == 2L
      out[[length(out) + 1L]] <- data.frame(
        family = fam, dataset = x$dataset[1], unit_a = x$unit[i], unit_b = x$unit[j],
        anatomical_unit_a = x$anatomical_unit[i], anatomical_unit_b = x$anatomical_unit[j],
        contrast = if (two) x$contrast[i] else anatomical_programme_pairwise_name(fam, x$unit[i], x$unit[j]),
        reuses_one_vs_rest = two, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

anatomical_programme_pairwise_specs <- function(dataset, unit_levels, families = anatomical_programme_families(),
                                                which_families = unique(families$family)) {
  p <- anatomical_programme_pairs(families)
  p <- p[p$dataset == dataset & !p$reuses_one_vs_rest & p$family %in% which_families, , drop = FALSE]
  specs <- list()
  for (i in seq_len(nrow(p))) {
    if (!all(c(p$anatomical_unit_a[i], p$anatomical_unit_b[i]) %in% unit_levels)) {
      stop("Pairwise contrast ", p$contrast[i], " needs units absent from the design.", call. = FALSE)
    }
    w <- stats::setNames(rep(0, length(unit_levels)), unit_levels)
    w[[p$anatomical_unit_a[i]]] <- 1; w[[p$anatomical_unit_b[i]]] <- -1
    specs[[p$contrast[i]]] <- control_spatial_contrast_spec(p$contrast[i], w, paste0(p$family[i], "_pairwise"), FALSE)
  }
  specs
}

# Specs of the added contrasts for one dataset, built from the unit levels
# present in its design. CA2 and CA3 strata are target-vs-mean-of-the-other
# strata of the same region, exactly like the locked CA1 strata; the microglia
# region contrasts are the registry's own (already manuscript_locked = FALSE).
anatomical_programme_added_specs <- function(dataset, unit_levels) {
  unit_levels <- sort(unique(as.character(unit_levels)))
  fam <- anatomical_programme_families()
  fam <- fam[fam$source == "added" & fam$dataset == dataset, , drop = FALSE]
  if (identical(dataset, "microglia")) {
    reg <- control_spatial_contrast_registry("microglia", unit_levels)
    missing <- setdiff(fam$contrast, names(reg))
    if (length(missing)) stop("Microglia registry lacks: ", paste(missing, collapse = ", "), call. = FALSE)
    return(reg[fam$contrast])
  }
  specs <- list()
  for (region in unique(sub("_strata$", "", fam$family))) {
    members <- unit_levels[startsWith(unit_levels, paste0(region, "_"))]
    want <- paste0(region, "_", fam$unit[fam$family == paste0(region, "_strata")])
    if (!setequal(members, want)) {
      stop(region, " strata in the data (", paste(members, collapse = ", "),
           ") differ from the planned set (", paste(want, collapse = ", "), ").", call. = FALSE)
    }
    for (target in want) {
      specs[[length(specs) + 1L]] <- control_spatial_contrast_spec(
        paste0(target, "_vs_mean_other_", region, "_strata"),
        control_spatial_target_rest_weights(members, target),
        paste0("neuropil_", region, "_stratum_vs_rest"), FALSE)
    }
  }
  names(specs) <- vapply(specs, function(s) s$name, character(1))
  specs
}

# GO BP branches whose terms are annotation artefacts when they top a brain
# tissue ranking (e.g. "sperm-egg recognition" is carried by the CCT
# chaperonin, "regulation of heart rate" by neuronal ion channels).
anatomical_programme_artefact_roots <- function() {
  c(reproduction = "GO:0000003", reproductive_process = "GO:0022414",
    circulatory_system_process = "GO:0003013",
    circulatory_system_development = "GO:0072359")
}

anatomical_programme_artefact_flag <- function(ids, roots = anatomical_programme_artefact_roots()) {
  ids <- as.character(ids)
  if (!length(ids)) return(logical())
  if (!requireNamespace("GO.db", quietly = TRUE)) stop("GO.db is required.", call. = FALSE)
  u <- unique(ids[!is.na(ids)])
  anc <- AnnotationDbi::mget(u, GO.db::GOBPANCESTOR, ifnotfound = NA)
  hit <- vapply(u, function(id) id %in% roots || any(anc[[id]] %in% roots), logical(1))
  unname(hit[ids])
}

anatomical_programme_genes <- function(x) {
  g <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  unique(g[nzchar(g)])
}

# Overlap coefficient of two core-enrichment sets: shared / smaller set.
anatomical_programme_overlap <- function(a, b) {
  if (!length(a) || !length(b)) return(0)
  length(intersect(a, b)) / min(length(a), length(b))
}

anatomical_programme_rule <- function(q_max = 0.05, max_size = 300L, overlap_max = 0.3, k = 3L) {
  list(q_max = q_max, max_size = as.integer(max_size), overlap_max = overlap_max, k = as.integer(k),
       text = sprintf(paste0(
         "presentation-only selection from completed control-only GO-BP GSEA: within each family of ",
         "sibling units, a term is a candidate for a unit when q < %s and NES > 0 there and it is not ",
         "q < %s with NES > 0 in any sibling unit; term size <= %d; terms in the reproduction and ",
         "circulatory-system branches of GO BP are left out; candidates ranked by NES (then q, then GO ID); ",
         "a term is skipped when more than %s of the smaller core-enrichment set is shared with a term ",
         "already chosen; at most %d terms per unit"),
         q_max, q_max, as.integer(max_size), overlap_max, as.integer(k)))
}

# go: completed GSEA rows (ID, Description, setSize, NES, p.adjust,
# core_enrichment, contrast); families: anatomical_programme_families().
# Returns every positive term of every unit with its candidate status, rank
# and selection, so the selection can be audited row by row.
anatomical_programme_candidates <- function(go, families = anatomical_programme_families(),
                                            rule = anatomical_programme_rule(),
                                            artefact_ids = NULL) {
  need <- c("ID", "Description", "setSize", "NES", "p.adjust", "core_enrichment", "contrast")
  miss <- setdiff(need, names(go))
  if (length(miss)) stop("GSEA table lacks: ", paste(miss, collapse = ", "), call. = FALSE)
  absent <- setdiff(families$contrast, unique(go$contrast))
  if (length(absent)) stop("GSEA table lacks contrast(s): ", paste(absent, collapse = ", "), call. = FALSE)
  go <- go[go$contrast %in% families$contrast & !is.na(go$NES) & !is.na(go$p.adjust), , drop = FALSE]
  go <- merge(go, families[, c("family", "unit", "contrast", "unit_order")], by = "contrast")
  go$supported_positive <- go$p.adjust < rule$q_max & go$NES > 0
  key <- paste(go$family, go$ID)
  n_pos <- tapply(go$supported_positive, key, sum)
  go$n_family_units_supported <- as.integer(n_pos[key])
  # artefact_ids: GO IDs to treat as artefact-branch terms (tests pass them
  # directly); by default they are looked up in GO.db.
  go$artefact_branch <- if (is.null(artefact_ids)) anatomical_programme_artefact_flag(go$ID) else go$ID %in% artefact_ids
  go$size_ok <- go$setSize <= rule$max_size
  go$unit_specific <- go$supported_positive & go$n_family_units_supported == 1L
  go$candidate <- go$unit_specific & go$size_ok & !go$artefact_branch
  out <- go[go$supported_positive, , drop = FALSE]
  out$candidate_rank <- NA_integer_; out$selected <- FALSE
  out$display_rank <- NA_integer_; out$redundant_with <- NA_character_
  for (ct in unique(out$contrast)) {
    i <- which(out$contrast == ct & out$candidate)
    i <- i[order(-out$NES[i], out$p.adjust[i], out$ID[i])]
    out$candidate_rank[i] <- seq_along(i)
    chosen <- list()
    for (j in i) {
      if (length(chosen) >= rule$k) break
      g <- anatomical_programme_genes(out$core_enrichment[j])
      ov <- vapply(chosen, function(h) anatomical_programme_overlap(g, h), numeric(1))
      if (length(ov) && any(ov > rule$overlap_max)) {
        out$redundant_with[j] <- names(chosen)[which.max(ov)]
        next
      }
      chosen[[out$ID[j]]] <- g
      out$selected[j] <- TRUE
      out$display_rank[j] <- length(chosen)
    }
  }
  out <- out[order(match(out$family, unique(families$family)), out$unit_order,
                   is.na(out$candidate_rank), out$candidate_rank, -out$NES), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# The displayed matrix: every selected term in every unit of its family, with
# the unit it was selected for. Non-supported cells keep their NES and q.
anatomical_programme_display <- function(go, candidates, families = anatomical_programme_families(),
                                         selected_col = "selected", rank_col = "display_rank") {
  sel <- candidates[candidates[[selected_col]] %in% TRUE, c("family", "ID", "unit", rank_col), drop = FALSE]
  names(sel)[names(sel) == "unit"] <- "selected_for_unit"
  names(sel)[names(sel) == rank_col] <- "display_rank"
  go <- merge(go[go$contrast %in% families$contrast,
                 intersect(c("contrast", "ID", "Description", "setSize", "NES", "pvalue", "p.adjust"), names(go)), drop = FALSE],
              families[, c("family", "family_label", "dataset", "source", "unit", "contrast", "unit_order")], by = "contrast")
  out <- merge(go, sel, by = c("family", "ID"))
  out$supported <- out$p.adjust < 0.05
  out$selected_for_order <- families$unit_order[match(paste(out$family, out$selected_for_unit),
                                                      paste(families$family, families$unit))]
  out <- out[order(match(out$family, unique(families$family)), out$selected_for_order,
                   out$display_rank, out$unit_order), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------- robustness of a term in its unit
# Evidence that a programme belongs to a unit, beyond its one-vs-rest q:
#   direct  : enriched in the unit against EACH sibling (pairwise contrasts; NES oriented
#             so that > 0 means higher in the unit)
#   sides   : NES of the unit's one-vs-rest contrast fitted on left-only and on right-only
#             samples (same animals; a consistency check, not replication)
#   proteins: the leading edge is carried by several proteins, not one, and many of them are
#             individually enriched (gene-level collapsed moderated t and FDR of the contrast)
# gene_stats: contrast, GeneSymbol, collapsed_statistic, collapsed_fdr, collapsed_logfc
# pairwise_go: contrast, ID, NES, p.adjust for every pair in anatomical_programme_pairs()
# side_go: contrast, side ("L"/"R"), ID, NES, p.adjust for one-vs-rest contrasts
anatomical_programme_term_metrics <- function(candidates, gene_stats, pairwise_go, side_go,
                                              families = anatomical_programme_families()) {
  x <- candidates
  pairs <- anatomical_programme_pairs(families)
  key_pw <- paste(pairwise_go$contrast, pairwise_go$ID)
  key_side <- paste(side_go$contrast, side_go$side, side_go$ID)
  gs <- split(gene_stats, gene_stats$contrast)
  n <- nrow(x)
  res <- data.frame(n_leading_edge = integer(n), top_leading_protein = NA_character_,
                    top_protein_share = NA_real_, n_leading_edge_supported = NA_integer_,
                    fraction_leading_edge_supported = NA_real_, pairwise_n_siblings = NA_integer_,
                    pairwise_min_NES = NA_real_, pairwise_max_q = NA_real_,
                    pairwise_supported_all = NA, pairwise_detail = NA_character_,
                    NES_left = NA_real_, q_left = NA_real_, NES_right = NA_real_, q_right = NA_real_,
                    both_hemispheres_positive = NA, stringsAsFactors = FALSE)
  for (i in seq_len(n)) {
    le <- anatomical_programme_genes(x$core_enrichment[i])
    res$n_leading_edge[i] <- length(le)
    g <- gs[[x$contrast[i]]]
    if (!is.null(g) && length(le)) {
      g <- g[g$GeneSymbol %in% le, , drop = FALSE]
      s <- pmax(g$collapsed_statistic, 0)
      if (nrow(g) && sum(s) > 0) {
        res$top_leading_protein[i] <- g$GeneSymbol[which.max(s)]
        res$top_protein_share[i] <- max(s) / sum(s)
      }
      sup <- sum(g$collapsed_fdr < 0.05 & g$collapsed_logfc > 0, na.rm = TRUE)
      res$n_leading_edge_supported[i] <- sup
      res$fraction_leading_edge_supported[i] <- sup / length(le)
    }
    p <- pairs[pairs$family == x$family[i] & (pairs$unit_a == x$unit[i] | pairs$unit_b == x$unit[i]), , drop = FALSE]
    if (nrow(p)) {
      j <- match(paste(p$contrast, x$ID[i]), key_pw)
      nes <- pairwise_go$NES[j] * ifelse(p$unit_a == x$unit[i], 1, -1)
      q <- pairwise_go$p.adjust[j]
      other <- ifelse(p$unit_a == x$unit[i], p$unit_b, p$unit_a)
      res$pairwise_n_siblings[i] <- nrow(p)
      res$pairwise_min_NES[i] <- if (all(is.na(nes))) NA_real_ else min(nes, na.rm = TRUE)
      res$pairwise_max_q[i] <- if (all(is.na(q))) NA_real_ else max(q, na.rm = TRUE)
      res$pairwise_supported_all[i] <- !anyNA(nes) && all(nes > 0 & q < 0.05)
      res$pairwise_detail[i] <- paste(sprintf("vs %s: NES %s, q %s", other,
                                              ifelse(is.na(nes), "NA", sprintf("%.2f", nes)),
                                              ifelse(is.na(q), "NA", formatC(q, format = "g", digits = 2))), collapse = "; ")
    }
    for (sd in c("L", "R")) {
      j <- match(paste(x$contrast[i], sd, x$ID[i]), key_side)
      nm <- if (sd == "L") c("NES_left", "q_left") else c("NES_right", "q_right")
      res[[nm[1]]][i] <- side_go$NES[j]; res[[nm[2]]][i] <- side_go$p.adjust[j]
    }
    res$both_hemispheres_positive[i] <- if (anyNA(c(res$NES_left[i], res$NES_right[i]))) NA else
      res$NES_left[i] > 0 && res$NES_right[i] > 0
  }
  cbind(x, res)
}

anatomical_programme_strict_rule <- function(q_max = 0.05, max_size = 300L, min_leading_edge = 5L,
                                             max_top_share = 0.4, overlap_max = 0.3, k = 3L) {
  list(q_max = q_max, max_size = as.integer(max_size), min_leading_edge = as.integer(min_leading_edge),
       max_top_share = max_top_share, overlap_max = overlap_max, k = as.integer(k),
       text = sprintf(paste0(
         "reproducible eligibility: q < %s and NES > 0 in the unit's one-vs-rest contrast; enriched against EACH sibling unit ",
         "in direct pairwise contrasts (NES > 0, q < %s); NES > 0 on left-only and on right-only samples; at least %d ",
         "leading-edge proteins, none carrying %s or more of the leading-edge statistic; term size <= %d; not in the ",
         "reproduction or circulatory-system branches of GO BP. Eligible terms ranked by NES; a term is skipped when ",
         "more than %s of the smaller core-enrichment set is shared with one already chosen; at most %d per unit"),
         q_max, q_max, as.integer(min_leading_edge), max_top_share, as.integer(max_size), overlap_max, as.integer(k)))
}

anatomical_programme_flags <- function(m, rule = anatomical_programme_strict_rule()) {
  m$flag_single_protein <- m$n_leading_edge < 3L | (!is.na(m$top_protein_share) & m$top_protein_share >= 0.5)
  m$flag_dominant_protein <- !is.na(m$top_protein_share) & m$top_protein_share >= 0.25
  m$flag_few_proteins <- m$n_leading_edge < rule$min_leading_edge
  m$flag_not_direct <- !(m$pairwise_supported_all %in% TRUE)
  m$flag_hemisphere <- !(m$both_hemispheres_positive %in% TRUE)
  m$strict_eligible <- m$supported_positive & m$size_ok & !m$artefact_branch & !m$flag_few_proteins &
    !(is.na(m$top_protein_share) | m$top_protein_share >= rule$max_top_share) &
    m$pairwise_supported_all %in% TRUE & m$both_hemispheres_positive %in% TRUE
  m
}

anatomical_programme_strict_select <- function(m, rule = anatomical_programme_strict_rule()) {
  m <- anatomical_programme_flags(m, rule)
  m$strict_rank <- NA_integer_; m$strict_selected <- FALSE; m$strict_display_rank <- NA_integer_
  for (ct in unique(m$contrast)) {
    i <- which(m$contrast == ct & m$strict_eligible)
    i <- i[order(-m$NES[i], m$p.adjust[i], m$ID[i])]
    m$strict_rank[i] <- seq_along(i)
    chosen <- list()
    for (j in i) {
      if (length(chosen) >= rule$k) break
      g <- anatomical_programme_genes(m$core_enrichment[j])
      ov <- vapply(chosen, function(h) anatomical_programme_overlap(g, h), numeric(1))
      if (length(ov) && any(ov > rule$overlap_max)) next
      chosen[[m$ID[j]]] <- g
      m$strict_selected[j] <- TRUE; m$strict_display_rank[j] <- length(chosen)
    }
  }
  m
}

# Hand-curated display picks (config/anatomical_programme_display_picks.csv): one row per
# displayed term with the reason it was chosen. Validated against the metrics; the report
# shows which criteria each pick meets, so a pick that fails one is visible, not hidden.
anatomical_programme_read_picks <- function(path) {
  p <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("family", "unit", "ID", "Description", "reason")
  miss <- setdiff(need, names(p))
  if (length(miss)) stop("Picks file lacks: ", paste(miss, collapse = ", "), call. = FALSE)
  if (anyDuplicated(paste(p$family, p$unit, p$ID))) stop("Duplicate pick in ", path, call. = FALSE)
  if (any(!nzchar(trimws(p$reason)))) stop("Every pick needs a reason.", call. = FALSE)
  p$display_rank <- stats::ave(seq_len(nrow(p)), p$family, p$unit, FUN = seq_along)
  p
}
