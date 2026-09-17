# Conservative, PHENOTYPE-BLIND adjudication of WGCNA module / supermodule names.
#
# WHAT THIS IS FOR
#   Establishing, for every module: what the network contains, which annotations
#   are robustly supported, whether the network CENTRE supports them, whether a
#   single coherent program exists, and what the safest human-readable reviewed
#   label would be. It produces PROPOSALS and evidence, never active labels.
#
# FIVE QUESTIONS KEPT SEPARATE (a term can pass A and fail B-E)
#   A. Is the enrichment statistically supported?
#   B. Does the theme describe the CENTRAL portion of the module?
#   C. Does the broader module support the same theme?
#   D. Is the theme specific enough and appropriate as a NAME?
#   E. Is the wording suitable for a manuscript rather than raw ontology syntax?
#
# EVIDENCE INDEPENDENCE (Phase 14 safeguard)
#   A candidate generated FROM a GO term must not be "validated" by that same
#   term. Evidence is tagged by role:
#     generating   - the BP terms that produced the candidate theme
#     supporting   - whether the same annotation describes the network centre
#                    (hub / core overlap with the generating terms)
#     orthogonal   - CC / MF ontology context, and marker/cell-context panels
#   These are NOT four independent experiments and are never counted as such.
#
# STRICT PHENOTYPE BLINDNESS
#   No group contrast, differential-abundance statistic or candidate tier may
#   enter any function here. wla_assert_phenotype_blind() gates every emitted
#   table.

# `%||%` comes from the canonical R/null_coalescing.R, loaded via R/paths.R.

# ---------------------------------------------------------------- vocabulary

wla_contract_version <- function() "wgcna_label_adjudication_v1"

wla_fdr_threshold <- function() 0.05

# Protein sets for which Stage 01 already computed ORA. All three are frozen
# artifacts; nothing is recomputed here.
wla_protein_sets <- function() c("all", "core_kME_0.6", "top_hub_25")

wla_coherence_classes <- function() {
  c("high_confidence_coherent", "coherent_but_wording_poor", "moderate_support",
    "peripheral_enrichment_only", "mixed_biology", "context_sensitive_signature",
    "unresolved")
}

wla_actions <- function() {
  c("KEEP", "REFINE_WORDING", "RENAME", "CONTEXTUALIZE", "MIXED", "UNRESOLVED")
}

wla_confidence_levels <- function() c("high", "moderate", "low")

wla_supermodule_classes <- function() {
  c("singleton", "coherent_single_program", "related_program_family",
    "mixed_but_structurally_coherent", "mixed_and_structurally_limited",
    "unresolved")
}

wla_evidence_roles <- function() c("generating", "supporting", "orthogonal")

# Column-name fragments that must never appear in a naming output.
wla_forbidden_field_patterns <- function() {
  c("\\bsus\\b", "sus_res", "\\bres_\\b", "log2fc", "fold_change", "\\bpadj\\b",
    "bh_fdr", "contrast", "tier_a", "tier_b", "tier_c", "tier_d",
    "candidate_tier", "group_effect", "phenotype_link", "(^|_)dap(s)?($|_)",
    "effect_size", "eigengene_estimate")
}

wla_assert_phenotype_blind <- function(data, label = "Adjudication table",
                                       allow = character()) {
  if (!is.data.frame(data)) stop(label, " must be a data frame.", call. = FALSE)
  nms <- setdiff(names(data), allow)
  low <- tolower(nms)
  hits <- unique(unlist(lapply(wla_forbidden_field_patterns(),
                               function(p) nms[grepl(p, low)])))
  if (length(hits)) {
    stop(label, " contains phenotype-derived column(s): ",
         paste(hits, collapse = ", "),
         ". Naming must remain phenotype-blind.", call. = FALSE)
  }
  invisible(TRUE)
}

.wla_stop <- function(...) stop(..., call. = FALSE)

.wla_require <- function(data, columns, label) {
  if (!is.data.frame(data)) .wla_stop(label, " must be a data frame.")
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    .wla_stop(label, " is missing required column(s): ",
              paste(missing, collapse = ", "), ".")
  }
  invisible(TRUE)
}

.wla_is_true <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
}

.wla_col <- function(data, column, default) {
  if (!is.null(data[[column]])) return(data[[column]])
  rep(default, nrow(data))
}

.wla_split_ids <- function(x) {
  ids <- trimws(strsplit(as.character(x %||% ""), "/", fixed = TRUE)[[1]])
  ids[nzchar(ids) & !is.na(ids)]
}

# ------------------------------------------------------------- theme layer

# Deterministic, transparent consolidation of redundant GO terms into broader
# biological themes.
#
# GO results are massively redundant: one module can carry hundreds of
# significant terms describing one program. Naming from a single lowest-p row is
# exactly the failure mode this layer exists to avoid.
#
# No new ontology dependency is introduced. Consolidation is by an explicit,
# inspectable keyword map over GO Descriptions; every constituent GO ID and
# Description is retained so a reviewer can always see what was merged.
# A term matching no theme keeps its own Description as a singleton theme, so
# nothing is silently discarded.
wla_theme_definitions <- function() {
  list(
    `mitochondrial energy metabolism` = c(
      "mitochondri", "oxidative phosphorylation", "respiratory (chain|electron)",
      "electron transport", "atp synthesis", "atp metabolic", "tricarboxylic",
      "citrate cycle", "acetyl-coa", "nadh", "cellular respiration",
      "generation of precursor metabolites", "proton motive force",
      "oxidoreduction-driven"),
    `cytoplasmic translation / ribosome` = c(
      "cytoplasmic translation", "ribosom", "translational initiation",
      "translational elongation", "peptide biosynthetic", "trna aminoacylation",
      "trna metabolic", "amino acid activation"),
    `RNA processing / RNP` = c(
      "rna processing", "rna splicing", "spliceosom", "mrna (processing|splicing|metabolic)",
      "ribonucleoprotein", "\\brnp\\b", "ncrna", "rrna", "p-body", "stress granule",
      "nuclear speck", "3'-utr", "poly\\(a\\)", "rna localization", "gene silencing by"),
    `protein folding / chaperone` = c(
      "protein folding", "chaperon", "unfolded protein", "heat shock",
      "de novo.*folding", "protein refolding"),
    `proteostasis / ubiquitin-proteasome` = c(
      "proteasom", "ubiquitin", "protein catabolic", "proteolysis",
      "modification-dependent", "erad", "autophag", "lysosom"),
    `synaptic organization / signalling` = c(
      "synap", "postsynap", "presynap", "neurotransmitter", "dendrit", "axon",
      "neuron projection", "trans-synaptic", "glutamatergic", "gabaergic",
      "neuron cellular homeostasis", "learning or memory"),
    `vesicle trafficking / exocytosis` = c(
      "vesicle", "exocyto", "endocyto", "snare", "membrane fusion",
      "golgi", "endosom", "secretion", "transport vesicle"),
    `cytoskeleton / actin-microtubule` = c(
      "cytoskelet", "actin", "microtubule", "tubulin", "myosin", "intermediate filament",
      "cell projection organization", "supramolecular fiber"),
    `myelin / oligodendrocyte ensheathment` = c(
      "ensheathment", "myelin", "oligodendrocyte", "axon ensheathment"),
    `epithelial / keratinization` = c(
      "keratinocyte", "keratinization", "cornified", "epiderm", "skin development",
      "desmosome", "hair follicle", "peptide cross-linking"),
    `immune / inflammatory` = c(
      "immune", "inflammat", "complement", "interferon", "antigen", "phagocyt",
      "cytokine", "leukocyte", "defense response"),
    `vascular / extracellular matrix` = c(
      "extracellular matrix", "basement membrane", "collagen", "laminin",
      "angiogen", "blood vessel", "endotheli", "integrin", "cell-substrate adhesion"),
    `ion transport / membrane potential` = c(
      "ion transport", "cation transport", "channel activity", "membrane potential",
      "ion homeostasis", "atpase-coupled", "transmembrane transport"),
    `signal transduction / phosphorylation` = c(
      "phosphorylation", "kinase", "phosphatase", "signal transduction",
      "gtpase", "second messenger", "signaling pathway"),
    `chromatin / nuclear organization` = c(
      "chromatin", "histone", "nucleosom", "dna repair", "dna replication",
      "transcription", "nuclear matrix")
  )
}

# Assign one theme to each GO Description. First matching theme in definition
# order wins, which makes the assignment deterministic and reproducible.
wla_assign_theme <- function(description) {
  desc <- tolower(as.character(description))
  defs <- wla_theme_definitions()
  out <- rep(NA_character_, length(desc))
  for (theme in names(defs)) {
    pattern <- paste(defs[[theme]], collapse = "|")
    hit <- is.na(out) & grepl(pattern, desc, perl = TRUE)
    out[hit] <- theme
  }
  # unmatched terms become their own theme so nothing is lost
  out[is.na(out)] <- paste0("other: ", desc[is.na(out)])
  out
}

# Per-module theme summary across a given protein set.
#
# `members` supplies the module's frozen kME ranks so contributor centrality can
# be computed; `go` is the frozen enrichment table.
wla_module_themes <- function(go, members, module_id, protein_set = "all",
                              ontology = "BP", fdr = wla_fdr_threshold()) {
  .wla_require(go, c("ModuleID", "Ontology", "ModuleProteinSetType", "ID",
                     "Description", "p.adjust", "Count", "GeneRatio", "geneID"),
               "GO enrichment table")
  sub <- go[as.character(go$ModuleID) == module_id &
              as.character(go$Ontology) == ontology &
              as.character(go$ModuleProteinSetType) == protein_set, , drop = FALSE]
  sub <- sub[is.finite(suppressWarnings(as.numeric(sub$p.adjust))) &
               suppressWarnings(as.numeric(sub$p.adjust)) <= fdr, , drop = FALSE]
  if (!nrow(sub)) {
    return(data.frame(
      ModuleID = character(), protein_set = character(), theme = character(),
      stringsAsFactors = FALSE
    ))
  }
  sub$theme <- wla_assign_theme(sub$Description)

  mem <- members[as.character(members$ModuleID) == module_id, , drop = FALSE]
  entrez <- as.character(mem$EntrezID)
  top10 <- .wla_is_true(mem$is_top10_module_hub)
  top25 <- .wla_is_true(mem$is_top_hub_25)
  core <- .wla_is_true(mem[["is_core_kME_0.6"]])
  kme <- suppressWarnings(as.numeric(mem$abs_kME))
  rf <- suppressWarnings(as.numeric(mem$rank_fraction))

  idx <- split(seq_len(nrow(sub)), sub$theme)
  rows <- lapply(names(idx), function(th) {
    rows_i <- idx[[th]]
    # union of contributing proteins across ALL terms in the theme
    ids <- unique(unlist(lapply(sub$geneID[rows_i], .wla_split_ids),
                         use.names = FALSE))
    contrib <- entrez %in% ids
    padj <- suppressWarnings(as.numeric(sub$p.adjust[rows_i]))
    ord <- order(padj, method = "radix")
    data.frame(
      ModuleID = module_id,
      protein_set = protein_set,
      ontology = ontology,
      theme = th,
      n_supporting_go_terms = length(rows_i),
      best_p_adjust = min(padj, na.rm = TRUE),
      representative_go_id = as.character(sub$ID[rows_i[ord[[1]]]]),
      representative_go_description = as.character(sub$Description[rows_i[ord[[1]]]]),
      supporting_go_ids = paste(as.character(sub$ID[rows_i[ord]]), collapse = "; "),
      supporting_go_descriptions = paste(
        utils::head(as.character(sub$Description[rows_i[ord]]), 12L), collapse = "; "),
      n_contributing_proteins = sum(contrib),
      fraction_module_contributing = if (nrow(mem)) sum(contrib) / nrow(mem) else NA_real_,
      fraction_core_contributing = if (any(core)) mean(contrib[core]) else NA_real_,
      fraction_top10_contributing = if (any(top10)) mean(contrib[top10]) else NA_real_,
      fraction_top25_contributing = if (any(top25)) mean(contrib[top25]) else NA_real_,
      median_contributor_abs_kME = if (any(contrib)) stats::median(kme[contrib], na.rm = TRUE) else NA_real_,
      median_contributor_rank_fraction = if (any(contrib)) stats::median(rf[contrib], na.rm = TRUE) else NA_real_,
      contributor_centrality_auc = wla_centrality_auc(rf, contrib),
      contributing_hub_symbols = {
        sel <- which(contrib & top10)
        sel <- sel[order(suppressWarnings(as.numeric(mem$abs_kME_rank_in_module[sel])))]
        paste(as.character(mem$GeneSymbol[sel]), collapse = "; ")
      },
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$best_p_adjust, method = "radix"), , drop = FALSE]
}

# Rank-based AUC: probability a contributor is more central than a
# non-contributor. 0.5 = indistinguishable; >0.5 = theme describes the centre.
wla_centrality_auc <- function(rank_fraction, is_contributor) {
  rf <- suppressWarnings(as.numeric(rank_fraction))
  grp <- .wla_is_true(is_contributor)
  ok <- is.finite(rf)
  rf <- rf[ok]; grp <- grp[ok]
  n1 <- sum(grp); n0 <- sum(!grp)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(rf, ties.method = "average")
  u <- sum(r[grp]) - n1 * (n1 + 1) / 2
  1 - (u / (n1 * n0))
}

# ------------------------------------------------- centrality support (P4)

# Real centrality support for a theme, replacing the inert
# hub_support = nonempty(hub_text) notion.
#
# Every component stays separately visible; no arbitrary bonus is added and the
# components are deliberately NOT collapsed into one opaque score.
wla_centrality_support <- function(theme_row = NULL) {
  # A module with no significant theme still needs a one-row support record, so
  # every optional field is read length-safely rather than via `$`, which would
  # yield a zero-length vector for an absent column.
  num1 <- function(col) {
    v <- if (is.null(theme_row) || is.null(theme_row[[col]])) NA_real_ else
      suppressWarnings(as.numeric(theme_row[[col]]))
    if (!length(v)) NA_real_ else v[[1]]
  }
  data.frame(
    top10_theme_fraction = num1("fraction_top10_contributing"),
    top25_theme_fraction = num1("fraction_top25_contributing"),
    core_theme_fraction = num1("fraction_core_contributing"),
    contributor_median_abs_kME = num1("median_contributor_abs_kME"),
    contributor_centrality_auc = num1("contributor_centrality_auc"),
    stringsAsFactors = FALSE
  )
}

# Does the network centre support this theme? Deliberately requires BOTH a real
# hub presence and a central rank distribution, so a theme carried by a handful
# of peripheral proteins cannot qualify.
wla_centre_supports_theme <- function(theme_row, min_top10 = 0.20,
                                      min_auc = 0.55) {
  t10 <- suppressWarnings(as.numeric(theme_row$fraction_top10_contributing))
  auc <- suppressWarnings(as.numeric(theme_row$contributor_centrality_auc))
  is.finite(t10) & is.finite(auc) & t10 >= min_top10 & auc >= min_auc
}

# --------------------------------------------------- wording preferences (P6)

# Raw ontology grammar that should not be used verbatim as a manuscript label,
# mapped to readable program wording. Applied ONLY when the theme evidence
# already supports the same biology; it never changes which biology is claimed.
wla_wording_map <- function() {
  c(
    "generation of precursor metabolites and energy" = "mitochondrial energy metabolism",
    "ensheathment of neurons" = "oligodendrocyte / myelin ensheathment",
    "modification-dependent macromolecule catabolic process" = "ubiquitin-proteasome protein catabolism",
    "formation of cytoplasmic translation initiation complex" = "cytoplasmic translation initiation",
    "cytoplasmic translational initiation" = "cytoplasmic translation initiation",
    "trna aminoacylation for protein translation" = "tRNA aminoacylation / translation",
    "keratinocyte differentiation" = "epithelial / keratinization signature",
    "skin development" = "epithelial / keratinization signature",
    "binding of sperm to zona pellucida" = "non-neural adhesion-glycoprotein signature",
    "regulatory ncrna-mediated gene silencing" = "RNA processing / RNP regulation",
    "proton motive force-driven atp synthesis" = "mitochondrial ATP synthesis",
    "proton motive force driven atp synthesis" = "mitochondrial ATP synthesis"
  )
}

# Terms too broad to name a module unless nothing better is supported.
wla_overly_broad_patterns <- function() {
  c("^protein phosphorylation$", "^phosphorylation$", "^protein modification",
    "^macromolecule (catabolic|metabolic)", "^cellular (process|metabolic)",
    "^regulation of ", "^organic substance", "^nitrogen compound",
    "homeostasis$", "^transport$", "^localization$", "^biosynthetic process$",
    "^hydrolase activity", "^catalytic activity", "^binding$")
}

# Themes whose biology is real but does not correspond literally to the sampled
# cell type, and must therefore be described as a SIGNATURE.
wla_context_sensitive_themes <- function() {
  c("epithelial / keratinization", "immune / inflammatory")
}

wla_readable_label <- function(theme, representative_description) {
  desc <- tolower(trimws(as.character(representative_description)))
  map <- wla_wording_map()
  mapped <- unname(map[desc])
  ifelse(!is.na(mapped), mapped,
         ifelse(grepl("^other: ", theme), sub("^other: ", "", theme), theme))
}

# ------------------------------------------------------- module adjudication

# Combine the three protein-set theme tables into one adjudicated record.
#
# Evidence roles are explicit (Phase 14): BP themes GENERATE the candidate,
# hub/core overlap with those same terms SUPPORTS it (not an independent test),
# and CC/MF plus marker context are ORTHOGONAL.
wla_adjudicate_module <- function(module_id, dataset, module_size,
                                  themes_all, themes_core, themes_top25,
                                  cc_terms = character(), mf_terms = character(),
                                  current_label = NA_character_,
                                  marker_context = NA_character_) {
  pick_top <- function(x) if (is.null(x) || !nrow(x)) NULL else x[1, , drop = FALSE]
  top_all <- pick_top(themes_all)
  top_core <- pick_top(themes_core)
  top_hub <- pick_top(themes_top25)

  # The candidate theme is the one with the strongest CONVERGENT support across
  # protein sets, not simply the lowest p-value in the all-protein set.
  candidates <- unique(stats::na.omit(c(
    if (!is.null(top_all)) top_all$theme,
    if (!is.null(top_core)) top_core$theme,
    if (!is.null(top_hub)) top_hub$theme
  )))
  n_sets_supporting <- vapply(candidates, function(th) {
    sum(c(
      !is.null(themes_all) && th %in% themes_all$theme,
      !is.null(themes_core) && th %in% themes_core$theme,
      !is.null(themes_top25) && th %in% themes_top25$theme
    ))
  }, integer(1))
  best_p <- vapply(candidates, function(th) {
    ps <- c(
      if (!is.null(themes_all)) themes_all$best_p_adjust[themes_all$theme == th],
      if (!is.null(themes_core)) themes_core$best_p_adjust[themes_core$theme == th]
    )
    if (!length(ps)) NA_real_ else min(ps, na.rm = TRUE)
  }, numeric(1))
  ord <- order(-n_sets_supporting, best_p, method = "radix")
  dominant <- if (length(candidates)) candidates[ord[[1]]] else NA_character_

  row_for <- function(tbl, th) {
    if (is.null(tbl) || !nrow(tbl) || is.na(th)) return(NULL)
    hit <- tbl[tbl$theme == th, , drop = FALSE]
    if (!nrow(hit)) NULL else hit[1, , drop = FALSE]
  }
  dom_all <- row_for(themes_all, dominant)
  dom_core <- row_for(themes_core, dominant)
  dom_hub <- row_for(themes_top25, dominant)

  # second theme, for mixed detection
  second <- if (length(candidates) > 1L) candidates[ord[[2]]] else NA_character_
  sec_all <- row_for(themes_all, second)

  support <- wla_centrality_support(dom_all)

  centre_ok <- if (!is.null(dom_all)) wla_centre_supports_theme(dom_all) else FALSE
  n_go_all <- if (!is.null(dom_all)) dom_all$n_supporting_go_terms else 0L
  rep_desc <- if (!is.null(dom_all)) dom_all$representative_go_description else
    if (!is.null(dom_core)) dom_core$representative_go_description else NA_character_
  best_padj_all <- if (!is.null(dom_all)) dom_all$best_p_adjust else NA_real_
  best_padj_core <- if (!is.null(dom_core)) dom_core$best_p_adjust else NA_real_
  best_padj_hub <- if (!is.null(dom_hub)) dom_hub$best_p_adjust else NA_real_

  any_sig <- is.finite(best_padj_all) || is.finite(best_padj_core)
  broad <- !is.na(dominant) && any(vapply(
    wla_overly_broad_patterns(),
    function(p) grepl(p, tolower(rep_desc %||% "")), logical(1)))
  context_sensitive <- !is.na(dominant) && dominant %in% wla_context_sensitive_themes()
  sec_share <- if (!is.null(sec_all)) sec_all$fraction_module_contributing else NA_real_
  dom_share <- if (!is.null(dom_all)) dom_all$fraction_module_contributing else NA_real_
  mixed <- is.finite(sec_share) && is.finite(dom_share) &&
    sec_share >= 0.75 * dom_share && !identical(second, dominant)

  # ---- classification (Phase 5)
  cls <- "unresolved"
  if (!any_sig) {
    cls <- "unresolved"
  } else if (context_sensitive && centre_ok) {
    cls <- "context_sensitive_signature"
  } else if (mixed) {
    cls <- "mixed_biology"
  } else if (centre_ok && n_go_all >= 3L && n_sets_supporting[ord[[1]]] >= 2L) {
    cls <- if (broad) "coherent_but_wording_poor" else "high_confidence_coherent"
  } else if (centre_ok) {
    cls <- "coherent_but_wording_poor"
  } else if (is.finite(support$contributor_centrality_auc) &&
             support$contributor_centrality_auc < 0.5) {
    cls <- "peripheral_enrichment_only"
  } else {
    cls <- "moderate_support"
  }

  # ---- confidence (Phase 17): high requires convergence AND centre support
  n_sets_dom <- if (length(candidates)) n_sets_supporting[[ord[[1]]]] else 0L
  confidence <- "low"
  if (cls %in% c("high_confidence_coherent") &&
      n_go_all >= 3L && n_sets_dom >= 2L && centre_ok && !mixed) {
    confidence <- "high"
  } else if (cls %in% c("coherent_but_wording_poor", "context_sensitive_signature",
                        "moderate_support") && any_sig) {
    confidence <- "moderate"
  }
  # A theme carried by a single protein set with no measurable centrality
  # evidence is not moderately confident, however small its p-value. Specificity
  # is not rewarded for its own sake.
  thin_evidence <- n_sets_dom < 2L &&
    !(is.finite(support$contributor_centrality_auc) && centre_ok)
  if (thin_evidence) confidence <- "low"

  # ---- proposed names (Phase 6). Never promoted automatically.
  readable <- if (is.na(dominant)) NA_character_ else
    wla_readable_label(dominant, rep_desc)
  proposed_primary <- switch(
    cls,
    unresolved = "mixed / unresolved",
    mixed_biology = "mixed / unresolved",
    peripheral_enrichment_only = "mixed / unresolved",
    context_sensitive_signature = if (!is.na(readable) && grepl("signature", readable))
      readable else paste0(readable, " signature"),
    readable
  )
  # a broad term is only allowed when nothing more specific is supported
  proposed_alt1 <- if (!is.na(dominant) && !is.null(dom_core)) {
    wla_readable_label(dominant, dom_core$representative_go_description)
  } else NA_character_
  proposed_alt2 <- if (!is.na(second)) wla_readable_label(second, second) else NA_character_
  fallback <- if (!is.na(dominant) && !identical(cls, "unresolved")) dominant else
    "mixed / unresolved"

  action <- switch(
    cls,
    high_confidence_coherent = if (!is.na(current_label) &&
      identical(tolower(trimws(current_label)), tolower(trimws(proposed_primary))))
      "KEEP" else "REFINE_WORDING",
    coherent_but_wording_poor = "REFINE_WORDING",
    context_sensitive_signature = "CONTEXTUALIZE",
    moderate_support = "REFINE_WORDING",
    peripheral_enrichment_only = "RENAME",
    mixed_biology = "MIXED",
    unresolved = "UNRESOLVED"
  )

  data.frame(
    dataset = dataset, ModuleID = module_id, module_size = module_size,
    dominant_theme = dominant %||% NA_character_,
    n_protein_sets_supporting_theme = if (length(candidates)) n_sets_supporting[[ord[[1]]]] else 0L,
    n_supporting_go_terms_all = n_go_all,
    best_p_adjust_all = best_padj_all,
    best_p_adjust_core = best_padj_core,
    best_p_adjust_top25 = best_padj_hub,
    best_go_all = if (!is.null(dom_all)) dom_all$representative_go_description else NA_character_,
    best_go_core = if (!is.null(dom_core)) dom_core$representative_go_description else NA_character_,
    best_go_top25 = if (!is.null(dom_hub)) dom_hub$representative_go_description else NA_character_,
    support,
    centre_supports_theme = centre_ok,
    second_theme = second %||% NA_character_,
    second_theme_module_fraction = sec_share,
    theme_is_overly_broad = broad,
    theme_is_context_sensitive = context_sensitive,
    coherence_class = cls,
    proposed_confidence = confidence,
    recommended_action = action,
    proposed_primary_label = proposed_primary,
    proposed_alternative_1 = proposed_alt1,
    proposed_alternative_2 = proposed_alt2,
    conservative_fallback_label = fallback,
    theme_supporting_top_hubs = if (!is.null(dom_all)) dom_all$contributing_hub_symbols else NA_character_,
    dominant_theme_evidence_is_thin = thin_evidence,
    orthogonal_cc_evidence = paste(utils::head(cc_terms, 5L), collapse = "; "),
    orthogonal_mf_evidence = paste(utils::head(mf_terms, 5L), collapse = "; "),
    orthogonal_marker_context = marker_context %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------- spot check

# Which top hubs support the proposed theme and which do not (Phase 15).
# A label is not required to fit every hub - WGCNA modules contain
# multifunctional proteins - but systematic contradiction downgrades it.
wla_hub_spot_check <- function(members, module_id, themes_all, dominant_theme) {
  mem <- members[as.character(members$ModuleID) == module_id &
                   .wla_is_true(members$is_top10_module_hub), , drop = FALSE]
  mem <- mem[order(suppressWarnings(as.numeric(mem$abs_kME_rank_in_module))), ,
             drop = FALSE]
  if (!nrow(mem)) {
    return(data.frame(supporting_hubs = NA_character_,
                      non_supporting_hubs = NA_character_,
                      n_hubs_supporting = 0L, stringsAsFactors = FALSE))
  }
  # Contributor identity is already resolved per theme during the theme build,
  # so the supporting hubs are read from there rather than re-derived.
  supporting <- if (!is.null(themes_all) && nrow(themes_all) && !is.na(dominant_theme)) {
    hit <- themes_all[themes_all$theme == dominant_theme, , drop = FALSE]
    if (nrow(hit)) trimws(strsplit(hit$contributing_hub_symbols[[1]], ";", fixed = TRUE)[[1]]) else character()
  } else character()
  supporting <- supporting[nzchar(supporting)]
  syms <- as.character(mem$GeneSymbol)
  data.frame(
    supporting_hubs = paste(intersect(syms, supporting), collapse = "; "),
    non_supporting_hubs = paste(setdiff(syms, supporting), collapse = "; "),
    n_hubs_supporting = length(intersect(syms, supporting)),
    stringsAsFactors = FALSE
  )
}

# Stable identifiers for the decisive hubs, so a human can verify against
# UniProt / MGI / GO later without putting external web state in the pipeline.
wla_hub_identifiers <- function(members, module_id, n = 10L) {
  mem <- members[as.character(members$ModuleID) == module_id, , drop = FALSE]
  mem <- mem[order(suppressWarnings(as.numeric(mem$abs_kME_rank_in_module))), ,
             drop = FALSE]
  mem <- utils::head(mem, n)
  paste(sprintf("%s(%s|%s)", as.character(mem$GeneSymbol),
                as.character(.wla_col(mem, "RepresentativeUniProt", NA_character_)),
                as.character(.wla_col(mem, "EntrezID", NA_character_))),
        collapse = "; ")
}

# ------------------------------------------------------ supermodule layer

# Module-balanced supermodule adjudication. Member modules are weighted EQUALLY:
# a 1,400-protein module must not outvote a 40-protein one.
wla_adjudicate_supermodules <- function(member_map, module_adjudication,
                                        structural, hubs_per_module = 3L) {
  .wla_require(member_map, c("dataset", "SupermoduleID", "ModuleID"),
               "Supermodule member map")
  .wla_require(module_adjudication, c("dataset", "ModuleID", "dominant_theme"),
               "Module adjudication")

  key <- paste(member_map$dataset, member_map$SupermoduleID, sep = "\r")
  idx <- split(seq_len(nrow(member_map)), key)
  parts <- do.call(rbind, strsplit(names(idx), "\r", fixed = TRUE))

  rows <- lapply(seq_along(idx), function(i) {
    ds <- parts[i, 1]; sm <- parts[i, 2]
    mods <- sort(as.character(member_map$ModuleID[idx[[i]]]))
    ev <- module_adjudication[module_adjudication$dataset == ds &
                                module_adjudication$ModuleID %in% mods, , drop = FALSE]
    ev <- ev[order(match(ev$ModuleID, mods)), , drop = FALSE]
    n_mod <- length(mods)

    # one vote per MEMBER MODULE, never per protein
    themes <- as.character(ev$dominant_theme)
    themes <- themes[!is.na(themes) & nzchar(themes)]
    tab <- sort(table(themes), decreasing = TRUE)
    dominant <- if (length(tab)) names(tab)[[1]] else NA_character_
    dom_n <- if (length(tab)) as.integer(tab[[1]]) else 0L
    second <- if (length(tab) > 1L) names(tab)[[2]] else NA_character_
    p <- if (length(tab)) as.numeric(tab) / sum(tab) else numeric(0)
    entropy <- if (length(p) > 1L) -sum(p * log(p)) / log(length(p)) else 0

    st <- structural[structural$dataset == ds & structural$SupermoduleID == sm, ,
                     drop = FALSE]
    getn <- function(col) if (nrow(st) && !is.null(st[[col]]))
      suppressWarnings(as.numeric(st[[col]][[1]])) else NA_real_
    min_cor <- getn("adjusted_signed_min_pairwise_eigengene_correlation")
    pc1 <- getn("pc1_variance_explained")
    stab <- getn("cut_height_stability_fraction_stable")

    structural_class <- if (n_mod <= 1L) "singleton" else if (
      is.finite(min_cor) && min_cor >= 0.5 && is.finite(stab) && stab >= 0.8
    ) "structurally_coherent" else if (
      is.finite(min_cor) && min_cor >= 0.3
    ) "structurally_related" else "structurally_limited"

    biological_class <- if (n_mod <= 1L) "singleton" else if (
      dom_n == n_mod
    ) "coherent_single_program" else if (
      dom_n >= ceiling(n_mod * 2 / 3)
    ) "related_program_family" else if (
      identical(structural_class, "structurally_coherent")
    ) "mixed_but_structurally_coherent" else if (
      is.finite(min_cor)
    ) "mixed_and_structurally_limited" else "unresolved"

    panel <- vapply(mods, function(m) {
      e <- ev[ev$ModuleID == m, , drop = FALSE]
      if (!nrow(e)) return(NA_character_)
      hubs <- trimws(strsplit(as.character(.wla_col(e, "top10_hub_symbols", "")[[1]]),
                              ";", fixed = TRUE)[[1]])
      hubs <- hubs[nzchar(hubs)]
      paste0(m, ": ", paste(utils::head(hubs, hubs_per_module), collapse = ", "))
    }, character(1))

    # singletons inherit the member module's proposal exactly
    inherit <- n_mod == 1L
    proposed <- if (inherit) {
      as.character(.wla_col(ev, "proposed_primary_label", NA_character_)[[1]])
    } else if (biological_class %in% c("coherent_single_program",
                                       "related_program_family")) {
      dominant
    } else "mixed / unresolved"
    conf <- if (inherit) {
      as.character(.wla_col(ev, "proposed_confidence", "low")[[1]])
    } else if (identical(biological_class, "coherent_single_program")) "high"
    else if (identical(biological_class, "related_program_family")) "moderate"
    else "low"
    action <- if (inherit) {
      as.character(.wla_col(ev, "recommended_action", "UNRESOLVED")[[1]])
    } else if (biological_class %in% c("coherent_single_program",
                                       "related_program_family")) "REFINE_WORDING"
    else "MIXED"

    data.frame(
      dataset = ds, SupermoduleID = sm, n_member_modules = n_mod,
      member_ModuleIDs = paste(mods, collapse = "; "),
      member_module_proposed_labels = paste(
        as.character(.wla_col(ev, "proposed_primary_label", NA_character_)),
        collapse = "; "),
      member_module_themes = paste(as.character(ev$dominant_theme), collapse = "; "),
      dominant_theme = dominant,
      n_modules_supporting_dominant_theme = dom_n,
      dominant_theme_module_fraction = if (n_mod) dom_n / n_mod else NA_real_,
      second_theme = second,
      theme_entropy = entropy,
      adjusted_signed_min_pairwise_correlation = min_cor,
      pc1_variance_explained = pc1,
      cut_height_stability_fraction_stable = stab,
      structural_coherence_class = structural_class,
      biological_coherence_class = biological_class,
      balanced_hub_panel = paste(stats::na.omit(panel), collapse = " | "),
      inherits_from_member_module = inherit,
      proposed_supermodule_label = proposed,
      proposed_confidence = conf,
      recommended_action = action,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$dataset, out$SupermoduleID, method = "radix"), , drop = FALSE]
}

# ------------------------------------------------- proposal registry (P9)

# A PROPOSED registry, schema-compatible with config/wgcna_labels/<ds>.csv but
# explicitly non-active. No human reviewer is fabricated.
wla_proposal_registry <- function(dataset, module_adj, supermodule_adj) {
  mk <- function(level, ids, label, short, conf, action, rationale) {
    data.frame(
      dataset = dataset, level = level, entity_id = ids,
      reviewed_biological_label = label, reviewed_short_label = short,
      confidence = conf, manual_review_required = TRUE,
      rationale = rationale, proposed_action = action,
      adjudication_status = "proposed",
      proposal_prepared_by = "automated_evidence_audit",
      reviewer = NA_character_, review_date = NA_character_,
      contract_version = wla_contract_version(),
      stringsAsFactors = FALSE
    )
  }
  mod <- mk("module", module_adj$ModuleID, module_adj$proposed_primary_label,
            module_adj$dominant_theme, module_adj$proposed_confidence,
            module_adj$recommended_action,
            paste0("class=", module_adj$coherence_class,
                   "; centre_supports_theme=", module_adj$centre_supports_theme,
                   "; n_GO_terms_all=", module_adj$n_supporting_go_terms_all,
                   "; protein_sets_supporting=", module_adj$n_protein_sets_supporting_theme))
  sm <- mk("supermodule", supermodule_adj$SupermoduleID,
           supermodule_adj$proposed_supermodule_label,
           supermodule_adj$dominant_theme, supermodule_adj$proposed_confidence,
           supermodule_adj$recommended_action,
           paste0("structural=", supermodule_adj$structural_coherence_class,
                  "; biological=", supermodule_adj$biological_coherence_class,
                  "; dominant theme in ",
                  supermodule_adj$n_modules_supporting_dominant_theme, "/",
                  supermodule_adj$n_member_modules, " member modules"))
  rbind(mod, sm)
}

# A proposal must NEVER satisfy the active-registry contract.
wla_validate_proposal <- function(registry, dataset, expected_modules,
                                  expected_supermodules) {
  .wla_require(registry, c("dataset", "level", "entity_id",
                           "reviewed_biological_label", "confidence",
                           "adjudication_status", "reviewer"),
               "Proposed registry")
  if (!all(registry$adjudication_status == "proposed")) {
    .wla_stop("A proposal registry must carry adjudication_status == 'proposed' ",
              "on every row.")
  }
  if (!all(is.na(registry$reviewer))) {
    .wla_stop("A proposal registry must not name a reviewer; human adjudication ",
              "has not happened.")
  }
  mods <- registry$entity_id[registry$level == "module"]
  sms <- registry$entity_id[registry$level == "supermodule"]
  if (!setequal(mods, expected_modules)) {
    .wla_stop("Proposed registry modules do not match the authoritative module set.")
  }
  if (!setequal(sms, expected_supermodules)) {
    .wla_stop("Proposed registry supermodules do not match the authoritative set.")
  }
  if (anyDuplicated(paste(registry$level, registry$entity_id))) {
    .wla_stop("Proposed registry contains duplicate level + entity_id rows.")
  }
  if (!all(registry$confidence %in% wla_confidence_levels())) {
    .wla_stop("Proposed registry confidence must be one of: ",
              paste(wla_confidence_levels(), collapse = ", "), ".")
  }
  invisible(TRUE)
}

# TRUE only for a registry that is genuinely active: adjudicated by a named
# human reviewer. Newly generated proposals can never satisfy this.
wla_is_active_registry <- function(registry) {
  if (!is.data.frame(registry) || !nrow(registry)) return(FALSE)
  status <- if (is.null(registry$adjudication_status)) {
    rep("adjudicated", nrow(registry))
  } else as.character(registry$adjudication_status)
  reviewer <- if (is.null(registry$reviewer)) rep(NA_character_, nrow(registry)) else
    as.character(registry$reviewer)
  all(status %in% c("adjudicated", "reviewed", "active")) &&
    all(!is.na(reviewer) & nzchar(trimws(reviewer)))
}

# ------------------------------------------------- consumer migration plan

# Which downstream scripts display a WGCNA biological label, which naming stage
# that field comes from, and what should happen after human adjudication.
#
# Each row was verified by inspecting the script's label column directly.
# Categories:
#   identity_critical     - keys on ModuleID/SupermoduleID; a label cannot
#                           affect correctness
#   biological_display    - shows a biological name to a human; should
#                           eventually read the reviewed canonical label
#   historical_provenance - deliberately records the Stage-01 label as
#                           provenance and must NOT be migrated
wla_consumer_migration_plan <- function() {
  rows <- list(
    c("analysis/integration/build_candidate_protein_shortlist.R",
      "wgcna_candidate_proteins_all.csv / shortlist workbook", "ModuleLabel_Final",
      "stage01", "biological_display", "canonical_biological_label", "high",
      "Low. Review-facing tables only; no structural contract depends on the label text.",
      "Migrate to the reviewed canonical label once neuropil/soma registries are approved."),
    c("analysis/integration/quantify_candidate_network_position.R",
      "wgcna_sus_res_network_position_protein_level.csv", "ModuleLabel_Final",
      "stage01", "biological_display", "canonical_biological_label", "high",
      "Low. Descriptive audit output; the analysis keys on ModuleID.",
      "Migrate with the candidate shortlist."),
    c("analysis/wgcna/audit_module_label_coherence.R",
      "WGCNA_module_label_coherence_audit.csv", "ModuleLabel_Final",
      "stage01", "historical_provenance", "keep ModuleLabel_Final", "none",
      "None. The audit exists to compare naming stages, so it must retain the raw Stage-01 label.",
      "Do not migrate. Add the reviewed label as an ADDITIONAL column after approval."),
    c("analysis/wgcna/adjudicate_module_labels.R",
      "WGCNA_module_adjudication.csv", "ModuleLabel_Final",
      "stage01", "historical_provenance", "keep ModuleLabel_Final", "none",
      "None. This table is the adjudication record and must show what was there before.",
      "Do not migrate."),
    c("analysis/integration/export_module_protein_zoom_source_data.R",
      "manuscript Figure 3 panels", "ModuleLabel_Final",
      "stage01", "biological_display", "canonical_biological_label", "low",
      "HIGH. Manuscript figure with hard structural contracts (exactly 15 proteins x 15 modules x 45 rows). Changing displayed text requires a deliberate figure re-freeze and caption update.",
      "Do NOT migrate in this pass. Revisit only after labels are approved and with an explicit figure re-freeze."),
    c("analysis/integration/render_module_circular_atlas.R",
      "wgcna_circular_atlas figures", "final_plot_label (plus one ModuleLabel_Final)",
      "stage07_mostly", "biological_display", "canonical_biological_label", "medium",
      "Medium. Already reads the Stage-07 final label in most places; one residual Stage-01 reference should be reconciled.",
      "Reconcile the residual Stage-01 reference, then inherit reviewed labels automatically via Stage 07."),
    c("analysis/wgcna/render_module_figures.R",
      "WGCNA publication figures", "canonical_biological_label",
      "stage07_reviewed", "biological_display", "canonical_biological_label", "none",
      "None. Already consumes the canonical reviewed label; it will pick up approved labels automatically.",
      "No change needed; it currently hard-stops for non-microglia datasets, which the registry generalization now makes unnecessary."),
    c("analysis/publication_source_data/build_biological_claims_table.R",
      "biological_claims_table.csv / .xlsx",
      "mixes ModuleLabel_Final, canonical_biological_label, final_plot_label, safe_display_label",
      "mixed", "biological_display", "canonical_biological_label", "high",
      "Medium. A submission artifact that currently mixes four label fields, so the same module can appear under different names in one export.",
      "Consolidate onto the canonical reviewed label after approval; keep Stage-01 only as an explicit provenance column."),
    c("analysis/wgcna/summarize_module_interpretation.R",
      "WGCNA_final_label_lookup.csv", "reviewed registry or automatic fallback",
      "stage07", "identity_critical", "canonical_biological_label", "medium",
      "Medium. This is where reviewed labels enter. The identical(ds, 'microglia') gates should become a check for an ACTIVE validating registry.",
      "Replace the dataset literal with an active-registry presence check once neuropil/soma registries are approved.")
  )
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- c("script", "output", "current_label_field", "label_stage",
                  "category", "recommended_canonical_field", "migration_priority",
                  "risk", "action_after_adjudication")
  out
}
