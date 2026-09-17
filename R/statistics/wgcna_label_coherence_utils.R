# Pure helpers for the PHENOTYPE-BLIND WGCNA module / supermodule label audit.
#
# QUESTION
#   Is each module's biological NAME actually supported by the proteins at the
#   network centre of that module, and by broader annotation evidence?
#
# STRICT PHENOTYPE BLINDNESS
#   Nothing in this file may read a group contrast, a differential-abundance
#   statistic, a SUS/RES/CON effect, or any candidate tier. Names must never be
#   chosen or scored using outcome information, or the biological annotation
#   becomes circular with the phenotype analysis it is later used to interpret.
#   wcl_assert_phenotype_blind() enforces this on every table this layer emits.
#
# WHAT IT DOES NOT DO
#   No WGCNA, no ORA and no label is recomputed or overwritten. The frozen ORA
#   result stays exactly as it is; this layer only asks whether the term that
#   NAMED the module describes its network centre, and proposes labels for
#   manual adjudication.

# `%||%` comes from the canonical R/null_coalescing.R, loaded via R/paths.R.

# ---------------------------------------------------------------- vocabulary

wcl_contract_version <- function() "wgcna_label_coherence_audit_v1"

# Column-name fragments that must never appear in a phenotype-blind output.
wcl_forbidden_field_patterns <- function() {
  c("sus", "\\bres\\b", "\\bcon\\b", "log2fc", "fold_change", "padj",
    "bh_fdr", "pvalue_da", "contrast", "tier_a", "tier_b", "tier_c", "tier_d",
    "candidate_tier", "group_effect", "phenotype", "stress", "da_", "_da$",
    "dap", "estimate", "effect_size")
}

# Guard: refuse to emit a naming-audit table that carries phenotype columns.
# `allow` lists deliberate exceptions (e.g. a literal blindness-attestation
# column that merely records that no phenotype field was used).
wcl_assert_phenotype_blind <- function(data, label = "Label audit table",
                                       allow = character()) {
  if (!is.data.frame(data)) stop(label, " must be a data frame.", call. = FALSE)
  nms <- setdiff(names(data), allow)
  low <- tolower(nms)
  hits <- character()
  for (pattern in wcl_forbidden_field_patterns()) {
    hits <- c(hits, nms[grepl(pattern, low)])
  }
  hits <- unique(hits)
  if (length(hits)) {
    stop(label, " contains phenotype-derived column(s): ",
         paste(hits, collapse = ", "),
         ". The naming audit must stay phenotype-blind.", call. = FALSE)
  }
  invisible(TRUE)
}

# Diagnostic coherence classes (module level).
wcl_module_coherence_classes <- function() {
  c("strongly_supported", "supported_but_wording_needs_refinement",
    "peripheral_GO_driven", "too_generic", "mixed_biology",
    "contextually_implausible_without_caveat", "unresolved")
}

wcl_module_actions <- function() {
  c("KEEP", "REFINE_WORDING", "RENAME_REVIEW", "CONTEXTUALIZE",
    "MIXED_UNRESOLVED")
}

wcl_supermodule_biological_classes <- function() {
  c("coherent_single_program", "related_program_family",
    "mixed_but_structurally_coherent", "mixed_and_structurally_limited",
    "singleton", "unresolved")
}

# GO Descriptions that describe a very broad process and make a poor module NAME
# even when the enrichment itself is sound.
wcl_generic_label_patterns <- function() {
  c("^protein phosphorylation$", "^phosphorylation$",
    "^protein modification", "^macromolecule (catabolic|metabolic)",
    "modification-dependent", "^cellular (process|metabolic)",
    "^regulation of ", "^positive regulation of ", "^negative regulation of ",
    "^organic substance", "^nitrogen compound", "homeostasis$",
    "^transport$", "^localization$", "^biosynthetic process$")
}

.wcl_stop <- function(...) stop(..., call. = FALSE)

.wcl_require <- function(data, columns, label) {
  if (!is.data.frame(data)) .wcl_stop(label, " must be a data frame.")
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    .wcl_stop(label, " is missing required column(s): ",
              paste(missing, collapse = ", "), ".")
  }
  invisible(TRUE)
}

.wcl_is_true <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
}

.wcl_col <- function(data, column, default) {
  if (!is.null(data[[column]])) return(data[[column]])
  rep(default, nrow(data))
}

# ------------------------------------------------- label-defining GO term

# Recover the GO term that produced the Stage-01 label, exactly as Stage 01
# selected it.
#
# analysis/05_wgcna/build_wgcna_modules.R:2746-2751 does:
#     filter(ModuleProteinSetType == "all") %>%
#     group_by(ModuleID, Ontology) %>%
#     arrange(p.adjust, qvalue, .by_group = TRUE) %>%
#     slice_head(n = 1)
# There is NO significance threshold and NO explicit tie-break, so ties fall
# through to the enrichment table's own row order. This reproduces that,
# including the stable ordering, and additionally reports how many terms were
# tied at the winning p.adjust so an arbitrary tie-break is visible.
wcl_label_defining_terms <- function(go_long, ontology = "BP",
                                     protein_set_type = "all") {
  .wcl_require(
    go_long,
    c("ModuleID", "Ontology", "ModuleProteinSetType", "ID", "Description",
      "p.adjust", "qvalue", "GeneRatio", "Count", "geneID"),
    "GO enrichment table"
  )
  sub <- go_long[
    as.character(go_long$Ontology) == ontology &
      as.character(go_long$ModuleProteinSetType) == protein_set_type, ,
    drop = FALSE
  ]
  if (!nrow(sub)) {
    return(data.frame(ModuleID = character(), stringsAsFactors = FALSE))
  }
  idx <- split(seq_len(nrow(sub)), as.character(sub$ModuleID))
  rows <- lapply(names(idx), function(mod) {
    rows_i <- idx[[mod]]
    go_padj <- suppressWarnings(as.numeric(sub$p.adjust[rows_i]))
    qval <- suppressWarnings(as.numeric(sub$qvalue[rows_i]))
    # stable order on (p.adjust, qvalue), ties keep table order
    ord <- order(go_padj, qval, method = "radix")
    win <- rows_i[ord[[1]]]
    n_tied <- sum(go_padj == go_padj[ord[[1]]], na.rm = TRUE)
    data.frame(
      ModuleID = mod,
      label_go_id = as.character(sub$ID[win]),
      label_go_description = as.character(sub$Description[win]),
      label_go_ontology = ontology,
      label_go_protein_set = protein_set_type,
      label_go_p_adjust = go_padj[ord[[1]]],
      label_go_qvalue = qval[ord[[1]]],
      label_go_gene_ratio = as.character(sub$GeneRatio[win]),
      label_go_count = suppressWarnings(as.integer(sub$Count[win])),
      label_go_is_fdr_significant = is.finite(go_padj[ord[[1]]]) &
        go_padj[ord[[1]]] <= 0.05,
      label_go_n_terms_tied_at_min = as.integer(n_tied),
      label_go_tie_broken_arbitrarily = n_tied > 1L,
      label_go_gene_ids = as.character(sub$geneID[win]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$ModuleID, method = "radix"), , drop = FALSE]
}

# Map the '/'-separated Entrez IDs in a GO geneID field back to canonical
# ProteinGroupIDs using the frozen membership's EntrezID column.
#
# The mapping is reported honestly: an Entrez ID may match several protein
# groups in the module, or none. Nothing is silently dropped.
wcl_map_label_contributors <- function(gene_ids, members) {
  .wcl_require(members, c("ProteinGroupID", "EntrezID"), "Module membership")
  ids <- trimws(strsplit(as.character(gene_ids %||% ""), "/", fixed = TRUE)[[1]])
  ids <- ids[nzchar(ids) & !is.na(ids)]
  member_entrez <- as.character(members$EntrezID)
  matched <- lapply(ids, function(e) which(member_entrez == e))
  n_hits <- vapply(matched, length, integer(1))
  list(
    entrez = ids,
    rows = sort(unique(unlist(matched, use.names = FALSE))),
    n_terms = length(ids),
    n_mapped_unique = sum(n_hits == 1L),
    n_mapped_multiple = sum(n_hits > 1L),
    n_unmapped = sum(n_hits == 0L)
  )
}

# --------------------------------------------- within-module rank enrichment

# Rank-based AUC: probability that a randomly chosen label-contributing protein
# is MORE central than a randomly chosen non-contributor.
#
# 0.5 = contributors sit no differently from the rest of the module.
# > 0.5 = contributors are more central (the label describes the network centre).
# < 0.5 = the term is driven by peripheral proteins, which lowers its
#         suitability as a NAME even though the ORA itself remains correct.
#
# Computed from the Mann-Whitney U statistic; diagnostic only, no p-value.
wcl_centrality_auc <- function(rank_fraction, is_contributor) {
  rf <- suppressWarnings(as.numeric(rank_fraction))
  grp <- .wcl_is_true(is_contributor)
  ok <- is.finite(rf)
  rf <- rf[ok]; grp <- grp[ok]
  n1 <- sum(grp); n0 <- sum(!grp)
  if (!n1 || !n0) return(NA_real_)
  # low rank_fraction = central, so rank ascending and take the contributor side
  r <- rank(rf, ties.method = "average")
  u <- sum(r[grp]) - n1 * (n1 + 1) / 2
  1 - (u / (n1 * n0))
}

# --------------------------------------------------- module coherence audit

# Per-module label-coherence evidence. `members` must already carry
# abs_kME_rank_in_module / n_module_members / rank_fraction.
wcl_module_label_evidence <- function(members, label_terms, top_n_hub = 10L) {
  .wcl_require(
    members,
    c("dataset", "ModuleID", "ProteinGroupID", "GeneSymbol", "EntrezID",
      "abs_kME", "abs_kME_rank_in_module", "n_module_members", "rank_fraction",
      "is_core_kME_0.6", "is_top_hub_25", "is_top10_module_hub",
      "is_top5_module_representative"),
    "Module membership"
  )
  .wcl_require(label_terms, c("ModuleID", "label_go_gene_ids"), "Label terms")

  idx <- split(seq_len(nrow(members)), as.character(members$ModuleID))
  rows <- lapply(names(idx), function(mod) {
    mem <- members[idx[[mod]], , drop = FALSE]
    term <- label_terms[label_terms$ModuleID == mod, , drop = FALSE]
    if (!nrow(term)) term <- label_terms[0, , drop = FALSE]

    mapped <- wcl_map_label_contributors(
      if (nrow(term)) term$label_go_gene_ids[[1]] else "", mem
    )
    contrib <- rep(FALSE, nrow(mem))
    contrib[mapped$rows] <- TRUE

    top10 <- .wcl_is_true(mem$is_top10_module_hub)
    top25 <- .wcl_is_true(mem$is_top_hub_25)
    top5 <- .wcl_is_true(mem$is_top5_module_representative)
    rf <- suppressWarnings(as.numeric(mem$rank_fraction))
    kme <- suppressWarnings(as.numeric(mem$abs_kME))

    hub_symbols <- function(flag) {
      sel <- which(flag)
      sel <- sel[order(suppressWarnings(as.numeric(mem$abs_kME_rank_in_module[sel])))]
      sym <- as.character(mem$GeneSymbol[sel])
      sym[is.na(sym) | !nzchar(sym)] <- as.character(mem$ProteinGroupID[sel])[is.na(sym) | !nzchar(sym)]
      paste(sym, collapse = "; ")
    }

    out <- data.frame(
      dataset = as.character(mem$dataset[[1]]),
      ModuleID = mod,
      module_size = nrow(mem),
      stringsAsFactors = FALSE
    )
    if (nrow(term)) {
      out <- cbind(out, term[, setdiff(names(term), c("ModuleID", "label_go_gene_ids")), drop = FALSE])
    }
    out$n_label_term_genes <- mapped$n_terms
    out$n_label_contributors_mapped <- length(mapped$rows)
    out$n_label_genes_mapped_unique <- mapped$n_mapped_unique
    out$n_label_genes_mapped_multiple <- mapped$n_mapped_multiple
    out$n_label_genes_unmapped <- mapped$n_unmapped
    out$label_mapping_fraction <- if (mapped$n_terms) {
      (mapped$n_mapped_unique + mapped$n_mapped_multiple) / mapped$n_terms
    } else NA_real_

    out$n_contributors_core_kME_0.6 <- sum(contrib & .wcl_is_true(mem[["is_core_kME_0.6"]]))
    out$fraction_contributors_core_kME_0.6 <- if (any(contrib)) {
      mean(.wcl_is_true(mem[["is_core_kME_0.6"]])[contrib])
    } else NA_real_
    out$n_contributors_in_top5 <- sum(contrib & top5)
    out$n_contributors_in_top10 <- sum(contrib & top10)
    out$n_contributors_in_top25 <- sum(contrib & top25)
    out$fraction_top10_hubs_contributing <- if (any(top10)) mean(contrib[top10]) else NA_real_
    out$fraction_top25_hubs_contributing <- if (any(top25)) mean(contrib[top25]) else NA_real_
    out$median_contributor_rank_fraction <- if (any(contrib)) stats::median(rf[contrib], na.rm = TRUE) else NA_real_
    out$median_contributor_abs_kME <- if (any(contrib)) stats::median(kme[contrib], na.rm = TRUE) else NA_real_
    out$median_noncontributor_abs_kME <- if (any(!contrib)) stats::median(kme[!contrib], na.rm = TRUE) else NA_real_
    out$contributor_centrality_auc <- wcl_centrality_auc(rf, contrib)

    out$top10_hub_symbols <- hub_symbols(top10)
    out$top25_hub_symbols <- hub_symbols(top25)
    out$label_contributor_symbols <- {
      sel <- which(contrib)
      sel <- sel[order(suppressWarnings(as.numeric(mem$abs_kME_rank_in_module[sel])))]
      sym <- as.character(mem$GeneSymbol[sel])
      sym[is.na(sym) | !nzchar(sym)] <- as.character(mem$ProteinGroupID[sel])[is.na(sym) | !nzchar(sym)]
      paste(utils::head(sym, 25L), collapse = "; ")
    }
    out
  })
  out <- do.call(rbind, rows)
  out[order(out$dataset, out$ModuleID, method = "radix"), , drop = FALSE]
}

# ------------------------------------------------------- coherence classing

# Assign the diagnostic coherence class and review action from the evidence.
#
# These are REVIEW AIDS, not automatic truth. The thresholds are deliberately
# coarse and are all recorded in the output so a reviewer can disagree.
# Tissue-implausible GO vocabularies: real, coherent enrichments whose wording
# implies a cell type that is not plausibly present in the profiled tissue.
# Such a module should be CONTEXTUALIZED as a signature, not renamed away and
# not asserted as a resident cell population.
wcl_nonneural_label_patterns <- function() {
  c("keratinocyte", "keratinization", "cornified", "epiderm", "skin development",
    "desmosome", "hair follicle", "sperm|fertilization|zona pellucida")
}

wcl_classify_module_labels <- function(evidence,
                                       auc_central = 0.55, auc_peripheral = 0.45,
                                       neural_dataset = TRUE) {
  .wcl_require(
    evidence,
    c("label_go_description", "label_go_is_fdr_significant",
      "contributor_centrality_auc", "fraction_top10_hubs_contributing",
      "median_contributor_rank_fraction"),
    "Module evidence"
  )
  desc <- tolower(as.character(evidence$label_go_description))
  sig <- .wcl_is_true(evidence$label_go_is_fdr_significant)
  auc <- suppressWarnings(as.numeric(evidence$contributor_centrality_auc))
  hub10 <- suppressWarnings(as.numeric(evidence$fraction_top10_hubs_contributing))
  tied <- .wcl_is_true(.wcl_col(evidence, "label_go_tie_broken_arbitrarily", FALSE))

  generic <- Reduce(`|`, lapply(wcl_generic_label_patterns(),
                                function(p) grepl(p, desc)), init = rep(FALSE, length(desc)))
  nonneural <- Reduce(`|`, lapply(wcl_nonneural_label_patterns(),
                                  function(p) grepl(p, desc)), init = rep(FALSE, length(desc)))

  cls <- rep("unresolved", nrow(evidence))
  # peripheral first, then generic, then support level
  cls[sig & is.finite(auc) & auc >= auc_central] <-
    "supported_but_wording_needs_refinement"
  cls[sig & is.finite(auc) & auc >= auc_central &
        is.finite(hub10) & hub10 >= 0.20] <- "strongly_supported"
  cls[sig & is.finite(auc) & auc <= auc_peripheral] <- "peripheral_GO_driven"
  cls[!sig] <- "unresolved"
  cls[generic] <- "too_generic"
  cls[!sig & tied] <- "mixed_biology"
  # A tissue-implausible term that is nonetheless well supported by the network
  # centre is REAL biology that needs a caveat, not a wrong enrichment.
  cls[neural_dataset & nonneural & sig] <- "contextually_implausible_without_caveat"

  action <- rep("MIXED_UNRESOLVED", length(cls))
  action[cls == "strongly_supported"] <- "KEEP"
  action[cls == "supported_but_wording_needs_refinement"] <- "REFINE_WORDING"
  action[cls == "peripheral_GO_driven"] <- "RENAME_REVIEW"
  action[cls == "too_generic"] <- "RENAME_REVIEW"
  action[cls == "mixed_biology"] <- "MIXED_UNRESOLVED"
  action[cls == "unresolved"] <- "RENAME_REVIEW"
  action[cls == "contextually_implausible_without_caveat"] <- "CONTEXTUALIZE"

  evidence$label_is_generic_wording <- generic
  evidence$label_is_tissue_implausible_wording <- nonneural
  evidence$coherence_class <- cls
  evidence$review_action <- action
  evidence$coherence_rationale <- paste0(
    "FDR-significant label term: ", sig,
    "; contributor centrality AUC: ",
    ifelse(is.finite(auc), sprintf("%.2f", auc), "NA"),
    "; fraction of top-10 hubs contributing to the label term: ",
    ifelse(is.finite(hub10), sprintf("%.2f", hub10), "NA"),
    "; generic wording: ", generic,
    "; tissue-implausible wording: ", nonneural,
    "; label term tied at minimum p.adjust: ", tied
  )
  evidence
}

# ---------------------------------------------- supermodule structural class

# Module-balanced supermodule evidence.
#
# Member modules are weighted EQUALLY for theme recurrence: a supermodule
# containing one 1428-member module and one 37-member module must not be named
# after the large one simply because it contributes more proteins. The
# illustrative hub panel is likewise balanced (top N hubs per member module).
wcl_supermodule_evidence <- function(member_map, module_evidence, structural,
                                     hubs_per_module = 3L) {
  .wcl_require(member_map, c("dataset", "SupermoduleID", "ModuleID"),
               "Supermodule member map")
  .wcl_require(module_evidence, c("dataset", "ModuleID"), "Module evidence")

  key <- paste(as.character(member_map$dataset),
               as.character(member_map$SupermoduleID), sep = "\r")
  idx <- split(seq_len(nrow(member_map)), key)
  parts <- do.call(rbind, strsplit(names(idx), "\r", fixed = TRUE))

  rows <- lapply(seq_along(idx), function(i) {
    rows_i <- idx[[i]]
    ds <- parts[i, 1]; sm <- parts[i, 2]
    mods <- sort(as.character(member_map$ModuleID[rows_i]))
    ev <- module_evidence[module_evidence$dataset == ds &
                            module_evidence$ModuleID %in% mods, , drop = FALSE]
    ev <- ev[order(match(ev$ModuleID, mods)), , drop = FALSE]

    themes <- as.character(.wcl_col(ev, "proposed_module_theme", NA_character_))
    themes <- themes[!is.na(themes) & nzchar(themes)]
    tab <- sort(table(themes), decreasing = TRUE)
    dominant <- if (length(tab)) names(tab)[[1]] else NA_character_
    dominant_n <- if (length(tab)) as.integer(tab[[1]]) else 0L
    second <- if (length(tab) > 1L) names(tab)[[2]] else NA_character_
    # equal-weight entropy over member-module themes
    p <- if (length(tab)) as.numeric(tab) / sum(tab) else numeric(0)
    entropy <- if (length(p) > 1L) -sum(p * log(p)) / log(length(p)) else 0

    st <- structural[structural$dataset == ds &
                       structural$SupermoduleID == sm, , drop = FALSE]
    n_mod <- length(mods)
    min_cor <- if (nrow(st)) suppressWarnings(as.numeric(
      st$adjusted_signed_min_pairwise_eigengene_correlation[[1]])) else NA_real_
    pc1 <- if (nrow(st)) suppressWarnings(as.numeric(st$pc1_variance_explained[[1]])) else NA_real_
    stab <- if (nrow(st)) suppressWarnings(as.numeric(st$cut_height_stability_fraction_stable[[1]])) else NA_real_

    structural_class <- if (n_mod <= 1L) "singleton" else if (
      is.finite(min_cor) && min_cor >= 0.5 && is.finite(stab) && stab >= 0.8
    ) "structurally_coherent" else if (
      is.finite(min_cor) && min_cor >= 0.3
    ) "structurally_related" else "structurally_limited"

    # balanced illustrative hub panel: top N per member module
    panel <- unlist(lapply(mods, function(m) {
      e <- module_evidence[module_evidence$dataset == ds &
                             module_evidence$ModuleID == m, , drop = FALSE]
      if (!nrow(e)) return(character())
      syms <- trimws(strsplit(as.character(e$top10_hub_symbols[[1]]), ";", fixed = TRUE)[[1]])
      paste0(m, ":", paste(utils::head(syms[nzchar(syms)], hubs_per_module),
                           collapse = ","))
    }), use.names = FALSE)

    biological_class <- if (n_mod <= 1L) "singleton" else if (
      dominant_n == n_mod && n_mod > 1L
    ) "coherent_single_program" else if (
      dominant_n >= ceiling(n_mod * 2 / 3)
    ) "related_program_family" else if (
      identical(structural_class, "structurally_coherent")
    ) "mixed_but_structurally_coherent" else if (
      is.finite(min_cor)
    ) "mixed_and_structurally_limited" else "unresolved"

    data.frame(
      dataset = ds,
      SupermoduleID = sm,
      n_member_modules = n_mod,
      member_ModuleIDs = paste(mods, collapse = "; "),
      member_module_labels = paste(as.character(
        .wcl_col(ev, "current_stage01_label", NA_character_)), collapse = "; "),
      member_module_proposed_themes = paste(
        as.character(.wcl_col(ev, "proposed_module_theme", NA_character_)),
        collapse = "; "),
      dominant_member_theme = dominant,
      n_modules_supporting_dominant_theme = dominant_n,
      dominant_theme_fraction = if (n_mod) dominant_n / n_mod else NA_real_,
      second_member_theme = second,
      member_theme_entropy = entropy,
      adjusted_signed_min_pairwise_correlation = min_cor,
      pc1_variance_explained = pc1,
      cut_height_stability_fraction_stable = stab,
      structural_coherence_class = structural_class,
      biological_coherence_class = biological_class,
      balanced_hub_panel = paste(panel, collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$dataset, out$SupermoduleID, method = "radix"), , drop = FALSE]
}

# A singleton supermodule must inherit its single member module's reviewed
# identity exactly, rather than being independently renamed.
wcl_apply_singleton_inheritance <- function(supermodule_evidence, module_evidence) {
  .wcl_require(supermodule_evidence,
               c("dataset", "SupermoduleID", "n_member_modules", "member_ModuleIDs"),
               "Supermodule evidence")
  single <- supermodule_evidence$n_member_modules == 1L
  inherited_label <- rep(NA_character_, nrow(supermodule_evidence))
  inherited_action <- rep(NA_character_, nrow(supermodule_evidence))
  for (i in which(single)) {
    mod <- trimws(supermodule_evidence$member_ModuleIDs[[i]])
    e <- module_evidence[module_evidence$dataset == supermodule_evidence$dataset[[i]] &
                           module_evidence$ModuleID == mod, , drop = FALSE]
    if (!nrow(e)) next
    inherited_label[[i]] <- as.character(.wcl_col(e, "proposed_label", NA_character_)[[1]])
    inherited_action[[i]] <- as.character(.wcl_col(e, "review_action", NA_character_)[[1]])
  }
  supermodule_evidence$inherits_from_member_module <- single
  supermodule_evidence$proposed_supermodule_label <- ifelse(
    single, inherited_label, supermodule_evidence$dominant_member_theme
  )
  supermodule_evidence$supermodule_review_action <- ifelse(
    single, inherited_action %||% "MIXED_UNRESOLVED",
    ifelse(supermodule_evidence$biological_coherence_class %in%
             c("coherent_single_program", "related_program_family"),
           "KEEP", "MIXED_UNRESOLVED")
  )
  supermodule_evidence
}
