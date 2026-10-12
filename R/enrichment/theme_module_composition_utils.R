# Theme x WGCNA module composition: do the phenotype-blind co-abundance structures (WGCNA modules and
# multi-module supermodules) group proteins the way the curated GSEA themes do?
#
# Evidence role: structural composition only. Theme membership comes from GO annotation of the theme's
# claim-eligible terms (manuscript_go_themes registry), restricted to the proteins the GSEA ranked; module
# membership comes from the WGCNA identity contract. Neither side uses a group label, so a supported overlap
# says that WGCNA's co-abundance structure carries the theme's proteins together. It is not evidence of a
# group effect, and no theme-level phenotype p-value or FDR is created here
# (docs/ATLAS_PROGRAM_SELECTION_AND_NAMING_RULES.md).
#
# Statistic per theme x entity: one-sided Fisher test (over-representation) on the 2x2 table within a stated
# protein universe, with the conditional MLE odds ratio and its exact 95% CI, fold enrichment, the share of
# the theme captured by the entity and the share of the entity in the theme. BH (and BY as a sensitivity)
# within dataset x entity level x theme-set definition x universe definition. A cell supports composition
# only with BH-FDR <= 0.05 AND odds ratio >= 2: the odds-ratio floor stops a trivially significant overlap
# with a very large module from counting.

TMC_SUPPORT_FDR <- 0.05
TMC_SUPPORT_OR <- 2
TMC_CLASS_BREAKS <- c(coherent = 0.5, core = 0.2)   # share of the theme captured by its supporting modules
TMC_EVIDENCE_ROLE <- "structural_composition_not_phenotype"

tmc_split_ids <- function(x) {
  x <- as.character(x)
  unique(unlist(strsplit(x[!is.na(x) & nzchar(x)], ";", fixed = TRUE), use.names = FALSE))
}

#' Composition statistics for one theme protein set and one entity protein set within a universe.
tmc_overlap_row <- function(theme_set, entity_set, universe) {
  universe <- unique(universe)
  t <- intersect(unique(theme_set), universe)
  e <- intersect(unique(entity_set), universe)
  a <- length(intersect(t, e)); b <- length(t) - a; c <- length(e) - a; d <- length(universe) - a - b - c
  if (min(a, b, c, d) < 0) stop("tmc_overlap_row: inconsistent 2x2 table", call. = FALSE)
  m <- matrix(c(a, b, c, d), 2L)
  p_greater <- stats::fisher.test(m, alternative = "greater")$p.value
  two_sided <- stats::fisher.test(m)           # conditional MLE odds ratio and its exact 95% CI
  share <- function(num, den) if (den > 0) num / den else NA_real_
  data.frame(
    n_universe = length(universe), n_theme = length(t), n_entity = length(e), n_overlap = a,
    frac_theme_in_entity = share(a, length(t)), frac_entity_in_theme = share(a, length(e)),
    fold_enrichment = if (length(t) > 0 && length(e) > 0) share(a, length(t)) / (length(e) / length(universe)) else NA_real_,
    odds_ratio = unname(two_sided$estimate), odds_ratio_ci_low = two_sided$conf.int[1],
    odds_ratio_ci_high = two_sided$conf.int[2],
    jaccard = share(a, length(union(t, e))), fisher_p_greater = p_greater,
    stringsAsFactors = FALSE)
}

#' All themes x all entities of one dataset and level, adjusted within that family.
#' theme_sets, entities: named lists of protein-group ids; entity_members: named character (entity -> members).
tmc_composition <- function(theme_sets, entities, universe, dataset, level, theme_set_definition,
                            universe_definition, entity_members = NULL, entity_definition = "all_members") {
  if (!length(theme_sets) || !length(entities)) return(NULL)
  rows <- list()
  for (th in names(theme_sets)) for (en in names(entities)) {
    r <- tmc_overlap_row(theme_sets[[th]], entities[[en]], universe)
    rows[[length(rows) + 1L]] <- cbind(data.frame(
      dataset = dataset, level = level, entity_id = en,
      member_modules = if (is.null(entity_members)) en else unname(entity_members[[en]]),
      theme_id = th, theme_set_definition = theme_set_definition, universe_definition = universe_definition,
      entity_definition = entity_definition,
      stringsAsFactors = FALSE), r)
  }
  out <- do.call(rbind, rows)
  out$BH_FDR <- stats::p.adjust(out$fisher_p_greater, "BH")
  out$BY_FDR <- stats::p.adjust(out$fisher_p_greater, "BY")
  out$family_size <- nrow(out)
  out$composition_supported <- out$BH_FDR <= TMC_SUPPORT_FDR & out$odds_ratio >= TMC_SUPPORT_OR
  out$evidence_role <- TMC_EVIDENCE_ROLE
  out
}

#' Proteins carried by exactly one of the themes (drops proteins shared between themes).
tmc_single_theme_sets <- function(theme_sets) {
  ids <- unlist(lapply(theme_sets, unique), use.names = FALSE)
  once <- names(which(table(ids) == 1L))
  lapply(theme_sets, function(s) intersect(s, once))
}

#' Class of a theme from the share its supporting modules capture.
tmc_class <- function(share, n_carriers) {
  ifelse(n_carriers == 0L | is.na(share), "scattered",
         ifelse(share >= TMC_CLASS_BREAKS[["coherent"]], "coherent",
                ifelse(share >= TMC_CLASS_BREAKS[["core"]], "core", "minor")))
}

#' Per dataset x theme summary at module level for one variant: the supporting modules (carriers), the share
#' of the theme they capture (modules partition the universe, so shares add), the top carrier and the class.
tmc_coherence_summary <- function(comp) {
  comp <- comp[comp$level == "module", , drop = FALSE]
  keys <- unique(comp[, c("dataset", "theme_id", "theme_set_definition", "universe_definition", "entity_definition")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    r <- comp[comp$dataset == k$dataset & comp$theme_id == k$theme_id &
                comp$theme_set_definition == k$theme_set_definition &
                comp$universe_definition == k$universe_definition &
                comp$entity_definition == k$entity_definition, , drop = FALSE]
    s <- r[r$composition_supported, , drop = FALSE]
    s <- s[order(-s$frac_theme_in_entity, s$entity_id), , drop = FALSE]
    share <- if (nrow(s)) sum(s$n_overlap) / r$n_theme[1] else 0
    data.frame(k, n_theme = r$n_theme[1], n_carriers = nrow(s),
               carriers = paste(s$entity_id, collapse = ";"),
               share_captured_by_carriers = share,
               top_carrier = if (nrow(s)) s$entity_id[1] else NA_character_,
               top_carrier_share = if (nrow(s)) s$frac_theme_in_entity[1] else NA_real_,
               top_carrier_fold = if (nrow(s)) s$fold_enrichment[1] else NA_real_,
               top_carrier_odds_ratio = if (nrow(s)) s$odds_ratio[1] else NA_real_,
               composition_class = tmc_class(share, nrow(s)), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  # themes whose carriers overlap within a dataset are not separated by WGCNA
  out$shares_carrier_with <- vapply(seq_len(nrow(out)), function(i) {
    if (!nzchar(out$carriers[i])) return("")
    mine <- strsplit(out$carriers[i], ";", fixed = TRUE)[[1]]
    same <- out$dataset == out$dataset[i] & out$theme_set_definition == out$theme_set_definition[i] &
      out$universe_definition == out$universe_definition[i] & out$entity_definition == out$entity_definition[i] &
      out$theme_id != out$theme_id[i]
    hit <- out$theme_id[same][vapply(out$carriers[same], function(cc) any(strsplit(cc, ";", fixed = TRUE)[[1]] %in% mine), logical(1))]
    paste(sort(hit), collapse = ";")
  }, character(1))
  out
}

#' Is the primary top carrier also supported in each sensitivity variant? One column per variant.
tmc_carrier_robustness <- function(summary_primary, comp_all, primary_def, primary_universe,
                                   primary_entities = "all_members") {
  vars <- unique(comp_all[comp_all$level == "module", c("theme_set_definition", "universe_definition", "entity_definition")])
  vars <- vars[!(vars$theme_set_definition == primary_def & vars$universe_definition == primary_universe &
                   vars$entity_definition == primary_entities), , drop = FALSE]
  out <- summary_primary
  for (i in seq_len(nrow(vars))) {
    v <- vars[i, ]
    nm <- paste0("top_carrier_supported_in__", v$theme_set_definition, "__", v$universe_definition, "__", v$entity_definition)
    out[[nm]] <- vapply(seq_len(nrow(out)), function(j) {
      if (is.na(out$top_carrier[j])) return(NA)
      hit <- comp_all$level == "module" & comp_all$dataset == out$dataset[j] & comp_all$theme_id == out$theme_id[j] &
        comp_all$entity_id == out$top_carrier[j] & comp_all$theme_set_definition == v$theme_set_definition &
        comp_all$universe_definition == v$universe_definition & comp_all$entity_definition == v$entity_definition
      if (!any(hit)) return(NA)
      isTRUE(comp_all$composition_supported[hit][1])
    }, logical(1))
  }
  out
}

#' Normalised mutual information (geometric mean of entropies) between two labelings.
tmc_nmi <- function(a, b) {
  tab <- table(a, b); p <- tab / sum(tab); pa <- rowSums(p); pb <- colSums(p)
  nz <- p > 0
  mi <- sum(p[nz] * log(p[nz] / outer(pa, pb)[nz]))
  ha <- -sum(pa[pa > 0] * log(pa[pa > 0])); hb <- -sum(pb[pb > 0] * log(pb[pb > 0]))
  if (ha <= 0 || hb <= 0) return(NA_real_)
  mi / sqrt(ha * hb)
}

#' Cramer's V between two labelings.
tmc_cramers_v <- function(a, b) {
  tab <- table(a, b)
  if (min(dim(tab)) < 2L) return(NA_real_)
  chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
  unname(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1L))))
}

#' Agreement between the theme partition and a WGCNA partition of the same proteins, against a label
#' permutation baseline (entity labels shuffled, theme labels fixed).
tmc_partition_agreement <- function(theme_label, entity_label, n_perm = 1000L, seed = 20261012L) {
  stopifnot(length(theme_label) == length(entity_label))
  obs_v <- tmc_cramers_v(theme_label, entity_label); obs_n <- tmc_nmi(theme_label, entity_label)
  set.seed(seed)
  perm <- vapply(seq_len(n_perm), function(i) {
    s <- sample(entity_label)
    c(tmc_cramers_v(theme_label, s), tmc_nmi(theme_label, s))
  }, numeric(2))
  data.frame(n_proteins = length(theme_label), n_theme_labels = length(unique(theme_label)),
             n_entity_labels = length(unique(entity_label)),
             cramers_v = obs_v, cramers_v_perm_mean = mean(perm[1, ]), cramers_v_perm_max = max(perm[1, ]),
             cramers_v_perm_p = (1 + sum(perm[1, ] >= obs_v)) / (1 + n_perm),
             nmi = obs_n, nmi_perm_mean = mean(perm[2, ]), nmi_perm_max = max(perm[2, ]),
             nmi_perm_p = (1 + sum(perm[2, ] >= obs_n)) / (1 + n_perm),
             n_permutations = n_perm, permutation_seed = seed, stringsAsFactors = FALSE)
}

#' Curated theme -> entity mapping rows classified against the composition (overlay only; never an input).
tmc_mapping_audit <- function(mapping, comp_primary) {
  keep <- mapping[as.character(mapping$approved_for_manuscript_interpretation) %in% c("TRUE", "true", "True"), , drop = FALSE]
  if (!nrow(keep)) return(NULL)
  rows <- lapply(seq_len(nrow(keep)), function(i) {
    k <- keep[i, ]
    hit <- comp_primary[comp_primary$dataset == k$dataset & comp_primary$theme_id == k$theme_id &
                          comp_primary$entity_id == k$entity_id, , drop = FALSE]
    st <- if (!nrow(hit)) "entity_not_tested" else if (isTRUE(hit$composition_supported[1])) "supported" else
      if (hit$BH_FDR[1] <= TMC_SUPPORT_FDR && hit$fold_enrichment[1] > 1) "weak_or_below_2" else "not_enriched"
    data.frame(dataset = k$dataset, theme_id = k$theme_id, entity_level = k$entity_level, entity_id = k$entity_id,
               mapping_role = if ("mapping_role" %in% names(k)) as.character(k$mapping_role) else NA_character_,
               confidence = if ("confidence" %in% names(k)) as.character(k$confidence) else NA_character_,
               n_overlap = if (nrow(hit)) hit$n_overlap[1] else NA_integer_,
               frac_theme_in_entity = if (nrow(hit)) hit$frac_theme_in_entity[1] else NA_real_,
               fold_enrichment = if (nrow(hit)) hit$fold_enrichment[1] else NA_real_,
               odds_ratio = if (nrow(hit)) hit$odds_ratio[1] else NA_real_,
               BH_FDR = if (nrow(hit)) hit$BH_FDR[1] else NA_real_,
               curated_mapping_composition_status = st, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  sup <- comp_primary[comp_primary$composition_supported, , drop = FALSE]
  unmapped <- sup[!paste(sup$dataset, sup$theme_id, sup$entity_id) %in% paste(keep$dataset, keep$theme_id, keep$entity_id), , drop = FALSE]
  if (nrow(unmapped)) {
    out <- rbind(out, data.frame(dataset = unmapped$dataset, theme_id = unmapped$theme_id, entity_level = unmapped$level,
                                 entity_id = unmapped$entity_id, mapping_role = NA_character_, confidence = NA_character_,
                                 n_overlap = unmapped$n_overlap, frac_theme_in_entity = unmapped$frac_theme_in_entity,
                                 fold_enrichment = unmapped$fold_enrichment, odds_ratio = unmapped$odds_ratio,
                                 BH_FDR = unmapped$BH_FDR, curated_mapping_composition_status = "supported_but_unmapped",
                                 stringsAsFactors = FALSE))
  }
  out
}

#' One row per module for the module wheel: identity, supermodule, size, spatial structure, and the themes
#' it carries in the primary variant (supported cells only, largest share first).
tmc_wheel_structure <- function(membership, sizes, spatial, comp_primary, theme_order) {
  m <- merge(membership, sizes, by = c("dataset", "module_id"), all.x = TRUE)
  m <- merge(m, spatial, by = c("dataset", "module_id"), all.x = TRUE)
  n_mod <- stats::aggregate(module_id ~ dataset + supermodule_id, data = m, FUN = length)
  names(n_mod)[3] <- "n_member_modules"
  m <- merge(m, n_mod, by = c("dataset", "supermodule_id"), all.x = TRUE)
  m$is_multi_module_supermodule <- m$n_member_modules > 1L
  sup <- comp_primary[comp_primary$level == "module" & comp_primary$composition_supported, , drop = FALSE]
  sup$theme_rank <- match(sup$theme_id, theme_order)
  carried <- lapply(seq_len(nrow(m)), function(i) {
    s <- sup[sup$dataset == m$dataset[i] & sup$entity_id == m$module_id[i], , drop = FALSE]
    s <- s[order(-s$frac_theme_in_entity, s$theme_rank), , drop = FALSE]
    if (!nrow(s)) return(c(carried_themes = "", top_theme = NA, top_theme_share = NA, top_theme_fold = NA,
                           top_theme_frac_module = NA))
    c(carried_themes = paste(s$theme_id, collapse = ";"), top_theme = s$theme_id[1],
      top_theme_share = s$frac_theme_in_entity[1], top_theme_fold = s$fold_enrichment[1],
      top_theme_frac_module = s$frac_entity_in_theme[1])
  })
  cb <- as.data.frame(do.call(rbind, carried), stringsAsFactors = FALSE)
  for (cl in c("top_theme_share", "top_theme_fold", "top_theme_frac_module")) cb[[cl]] <- as.numeric(cb[[cl]])
  out <- cbind(m, cb)
  out[order(out$dataset, out$supermodule_id, out$module_id), , drop = FALSE]
}
