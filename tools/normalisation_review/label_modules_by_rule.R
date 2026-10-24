# Rule-based classification of WGCNA modules and supermodules for both network versions (normalisation review,
# 2026-10-22). The same code and thresholds are applied to the original-centring network (archive) and to the
# complete-case network (working tree), so the two can be shown side by side. Only the wording of each label is
# written by hand; the threshold test, the signal class, the lineage class and the confidence come from here.
#
#   signals     significant GO BP terms (BH q < 0.05; the module's own over-representation from the network build,
#               against the tested proteome) are grouped by shared proteins: overlap coefficient
#               |A and B| / min(|A|, |B|) between the terms' module proteins, average-linkage clustering, cut at 0.5
#               (the EnrichmentMap default; nested terms such as a complex and its subunit fall into one signal).
#               The 40 most significant terms are grouped; a group whose proteins are mostly (overlap coefficient
#               >= 0.5) those of a larger group is then merged into it, so separate signals are carried by different
#               proteins. Passing signals whose proteins mostly (at least half of each) sit in the same significant
#               CC term of the module are then one signal in that compartment (respiration and ATP synthesis in the
#               mitochondrial inner membrane). Signal size = proteins in its largest term. Signals are
#               processes (BP); CC and MF terms enter through the two-ontology requirement and the wording, so broad
#               location terms (extracellular space, vesicle membrane) do not count as signals. A module without a
#               passing BP signal is classified on the signals over all three ontologies (a location- or complex-defined
#               module, e.g. axon/microtubule/ER membrane), and the basis is recorded; a weak association is also named
#               from all three ontologies.
#   passing     the module needs significant terms in at least two ontologies; a signal then passes the label
#               threshold when its leading (most significant) term is carried by more than 10 proteins. (BP and CC
#               terms of the same biology often share too few proteins to fall into one signal, so the ontology
#               requirement is on the module, not on each signal.)
#   signal class unannotated: no passing signal ('weak' when exactly one signal exists: the label names it as a weak
#               association; 'none' otherwise). Otherwise the largest passing signal is the dominant one, and every
#               passing signal at least half its size is a comparable signal: 'single signal' (only the dominant),
#               'two signals' (both are named in the label) or 'several signals' (three or more).
#   lineage     new network only: carried when one original-network module supplies at least 50% of the proteins and
#               its dominant signal is found again among the comparable signals (its leading term is one of their
#               terms, or the signals share at least 30% of their terms); revised when it supplies 50% but its
#               dominant signal is not found again; new below 50%.
#   confidence  coverage = proteins in the dominant signal's largest term / gene-mapped module proteins; hub support =
#               share of the 10 highest-kME gene-mapped hubs that are in the dominant signal.
#               single signal: high when coverage >= 15% and hub support >= 50%; moderate when coverage >= 8%, or
#               the module holds at least half of a term with at least 50 proteins; low otherwise.
#               two signals: moderate when coverage >= 8% (or the complex exception), low otherwise.
#               several signals and weak associations: low. Unannotated without a signal: none.
#   supermodules singleton (one member: its classes and confidence); unannotated (no member labelled); mixed (one member
#               unannotated, another labelled); concordant (all members labelled, every pair of dominant signals shares
#               at least 30% of its terms); compound (all members labelled, different dominant signals). Confidence of
#               concordant and compound = the lowest member confidence; mixed: low.
#   stability   each module is matched in the 9 leave-one-animal-out rebuilds (validation/wgcna_builds, own power):
#               best Jaccard per rebuild, median and minimum over animals, and the number of rebuilds with Jaccard
#               >= 0.5. Only rebuilds that reproduce the network (ARI >= 0.99) are used. Descriptive, not a filter.
#
# Writes results/reviewer_audit/normalisation_choice_2026-10-21/annotation_review/rule_based/: module_rule_classes.csv,
# module_signals.csv (the grouping that classifies each module, with the reason a signal fails the threshold),
# module_signals_all_ontologies.csv (unannotated modules: the all-ontology grouping that decides weak or unannotated),
# signal_terms.csv (the GO terms of every signal), supermodule_rule_classes.csv, rules.yml.
# USAGE  Rscript tools/normalisation_review/label_modules_by_rule.R      (working directory: repo root)

suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })
source(file.path("R", "paths.R"))
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
REV <- path_results("reviewer_audit", "normalisation_choice_2026-10-21")
OUT <- file.path(REV, "annotation_review", "rule_based"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
BUILDS <- file.path(REV, "validation", "wgcna_builds")
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
RULES <- list(q = 0.05, top_terms = 40L, overlap = 0.5, min_proteins = 10L, min_ontologies = 2L, comparable_ratio = 0.5,
              high_cov = 0.15, mod_cov = 0.08, hub_n = 10L, hub_share = 0.5, complex_share = 0.5, complex_n = 50L,
              term_jaccard = 0.3, carry_share = 0.5, stable_jaccard = 0.5, rebuild_ari = 0.99)
VERSIONS <- list(
  original = list(estimator = "original_all_observed",
                  state = function(ds) file.path(ARCHIVE, "data", "processed", "06_modules_WGCNA", "01_WGCNA", ds, "wgcna_final_model_state.rds"),
                  sm = function(ds) file.path(ARCHIVE, "results", "tables", "06_modules_WGCNA", "module_annotation", ds, "WGCNA_supermodule_biological_annotation.csv")),
  new = list(estimator = "new_complete_case",
             state = function(ds) file.path("results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"),
             sm = function(ds) file.path("results", "wgcna", "annotate_module_microenvironment", ds, "tables", "WGCNA_supermodule_biological_annotation.csv")))
is_module <- function(m) !is.na(m) & !grepl("grey|unassigned|m00", m)
# gene annotation and module GO limits of a WGCNA state (states built before 2026-10-23 used the legacy annotation and 10-500)
state_evidence <- function(s) paste(s$parameters$gene_annotation_contract_version %||% "legacy_uniprot_gene_name",
                                    s$parameters$module_go_min_gs_size %||% 10L, s$parameters$module_go_max_gs_size %||% 500L, sep = " / ")
evidence_seen <- character()
jac <- function(a, b) if (!length(a) && !length(b)) NA_real_ else length(intersect(a, b)) / length(union(a, b))
ids_of <- function(x) if (is.na(x) || !nzchar(x)) character() else strsplit(x, ";", fixed = TRUE)[[1]]

overlap_matrix <- function(sets) {
  k <- length(sets); M <- matrix(1, k, k)
  for (i in seq_len(k)) for (j in seq_len(k)) if (i < j)
    M[i, j] <- M[j, i] <- length(intersect(sets[[i]], sets[[j]])) / min(length(sets[[i]]), length(sets[[j]]))
  M
}

# signals of one module from its significant GO terms; passing signals first, largest first
module_signals <- function(go_m, onts = c("BP", "CC", "MF")) {
  n_onts <- n_distinct(go_m$Ontology[go_m$p.adjust < RULES$q])
  sig <- go_m |> filter(p.adjust < RULES$q, Ontology %in% onts) |> arrange(p.adjust) |> head(RULES$top_terms)
  if (!nrow(sig)) return(tibble())
  sets <- lapply(strsplit(as.character(sig$geneID), "/", fixed = TRUE), unique)
  cl <- if (nrow(sig) == 1) 1L else cutree(hclust(as.dist(1 - overlap_matrix(sets)), method = "average"), h = 1 - RULES$overlap)
  sig$cluster <- cl; sig$genes <- sets
  grp <- bind_rows(lapply(split(seq_len(nrow(sig)), sig$cluster), function(ix) {
    s <- sig[ix, ]; lead <- s[order(s$p.adjust, -s$Count)[1], ]; best <- s[which.max(s$Count), ]   # ties in q: the larger term leads
    tibble(leading_term = lead$Description, leading_id = lead$ID, leading_ontology = lead$Ontology, leading_q = lead$p.adjust,
           leading_count = lead$Count, best_count = best$Count, best_term = best$Description,
           n_terms = nrow(s), ontologies = paste(sort(unique(s$Ontology)), collapse = "+"), n_ontologies = n_distinct(s$Ontology),
           n_proteins = length(unique(unlist(s$genes))), term_ids = list(s$ID), genes = list(unique(unlist(s$genes))),
           top_terms = paste(head(s$Description[order(s$p.adjust)], 4), collapse = " | "),
           complex = any(s$Count >= RULES$complex_n & s$Count / as.integer(sub("/.*", "", s$BgRatio)) >= RULES$complex_share))
  })) |> arrange(desc(best_count), leading_q)
  # second step: a term group whose proteins are mostly those of a larger group (overlap coefficient >= 0.5) is the same
  # biology seen through other terms, and is merged into it; only groups carried by different proteins stay separate
  keep <- list()
  for (i in seq_len(nrow(grp))) {
    hit <- which(vapply(keep, function(k) length(intersect(k$genes[[1]], grp$genes[[i]])) / min(length(k$genes[[1]]), length(grp$genes[[i]])), numeric(1)) >= RULES$overlap)
    if (length(hit)) keep[[hit[1]]] <- merge_signal(keep[[hit[1]]], grp[i, ]) else keep[[length(keep) + 1L]] <- grp[i, ]
  }
  bind_rows(keep) |> mutate(passing = n_onts >= RULES$min_ontologies & leading_count > RULES$min_proteins, n_module_ontologies = n_onts) |>
    arrange(desc(passing), desc(best_count), leading_q) |> mutate(signal = seq_len(n()), .before = 1)
}

merge_signal <- function(k, g) {
  k$term_ids <- list(c(k$term_ids[[1]], g$term_ids[[1]])); k$genes <- list(union(k$genes[[1]], g$genes[[1]]))
  k$n_terms <- k$n_terms + g$n_terms; k$n_proteins <- length(k$genes[[1]]); k$complex <- k$complex || g$complex
  k$ontologies <- paste(sort(unique(c(strsplit(k$ontologies, "+", fixed = TRUE)[[1]], strsplit(g$ontologies, "+", fixed = TRUE)[[1]]))), collapse = "+")
  if (g$leading_q < k$leading_q || (g$leading_q == k$leading_q && g$leading_count > k$leading_count)) { k$leading_term <- g$leading_term; k$leading_id <- g$leading_id; k$leading_q <- g$leading_q
                                   k$leading_count <- g$leading_count; k$leading_ontology <- g$leading_ontology }
  if (g$best_count > k$best_count) { k$best_count <- g$best_count; k$best_term <- g$best_term }
  tt <- unique(c(strsplit(k$top_terms, " | ", fixed = TRUE)[[1]], strsplit(g$top_terms, " | ", fixed = TRUE)[[1]]))
  k$top_terms <- paste(head(tt, 6), collapse = " | ")
  k
}

# third step: passing signals whose proteins mostly sit in the same significant CC term of the module (at least half of
# each signal's proteins) are one signal in that compartment (e.g. respiration and ATP synthesis in mitochondrial inner
# membrane); the compartment is recorded for the wording
group_by_compartment <- function(sg, cc) {
  if (!nrow(sg)) return(sg)
  sg$compartment <- NA_character_
  in_cc <- function(genes) vapply(cc$genes, function(t) length(intersect(genes, t)) / length(genes), numeric(1))
  for (i in seq_len(nrow(sg))) if (sg$passing[i] && length(cc$genes)) {
    f <- in_cc(sg$genes[[i]]); if (any(f >= RULES$overlap)) sg$compartment[i] <- cc$term[which(f >= RULES$overlap)[1]]
  }
  p <- which(sg$passing)
  if (length(p) < 2 || !length(cc$genes)) return(sg)
  share <- lapply(p, function(i) in_cc(sg$genes[[i]]) >= RULES$overlap)
  keep <- list(); kept_share <- list()
  for (j in seq_along(p)) {
    hit <- which(vapply(kept_share, function(ks) any(ks & share[[j]]), logical(1)))
    if (length(hit)) {
      h <- hit[1]; common <- which(kept_share[[h]] & share[[j]])[1]
      keep[[h]] <- merge_signal(keep[[h]], sg[p[j], ]); keep[[h]]$compartment <- cc$term[common]
      kept_share[[h]] <- kept_share[[h]] | share[[j]]
    } else { keep[[length(keep) + 1L]] <- sg[p[j], ]; kept_share[[length(kept_share) + 1L]] <- share[[j]] }
  }
  bind_rows(bind_rows(keep), sg[-p, ]) |> arrange(desc(passing), desc(best_count), leading_q) |> mutate(signal = seq_len(n()))
}

classify_module <- function(sg, n_mapped, hub_ent, n_all_signals) {
  out <- list(n_signals = nrow(sg), n_passing = if (nrow(sg)) sum(sg$passing) else 0L, comparable = integer())
  if (!nrow(sg) || !any(sg$passing)) {
    out$signal_class <- if (n_all_signals == 1) "unannotated (weak)" else "unannotated"
    out$confidence <- if (n_all_signals == 1) "low" else "none"; out$dominant <- NA_integer_
    return(out)
  }
  p <- which(sg$passing); d <- p[1]
  comp <- p[sg$best_count[p] >= RULES$comparable_ratio * sg$best_count[d]]
  out$dominant <- d; out$comparable <- comp
  out$second_ratio <- if (length(p) > 1) sg$best_count[p[2]] / sg$best_count[d] else 0
  out$coverage <- sg$best_count[d] / n_mapped
  out$hub_support <- if (length(hub_ent)) mean(hub_ent %in% sg$genes[[d]]) else NA_real_
  out$signal_class <- c("single signal", "two signals", "several signals")[min(length(comp), 3L)]
  hub_ok <- !is.na(out$hub_support) && out$hub_support >= RULES$hub_share
  mod_ok <- out$coverage >= RULES$mod_cov || sg$complex[d]
  out$confidence <- switch(out$signal_class,
    "single signal" = if (out$coverage >= RULES$high_cov && hub_ok) "high" else if (mod_ok) "moderate" else "low",
    "two signals" = if (mod_ok) "moderate" else "low",
    "several signals" = "low")
  out
}

ari <- function(a, b) {
  ids <- intersect(names(a), names(b)); a <- a[ids]; b <- b[ids]
  tab <- table(a, b); n <- sum(tab); c2 <- function(v) sum(v * (v - 1) / 2)
  sij <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); e <- sa * sb / (n * (n - 1) / 2)
  (sij - e) / ((sa + sb) / 2 - e)
}
# rebuilds that reproduce the pipeline network: the archived-input rebuilds (rebuild_archived_network_loao.R) where they
# exist, otherwise the review's rebuilds of that version
stability <- function(ds, v, est, long) {
  prefixes <- c(if (v == "original") sprintf("%s_archived_original", ds), sprintf("%s_%s", ds, est))
  pipe <- setNames(ifelse(is_module(long$ModuleID), long$ModuleID, "grey"), long$ProteinGroupID)
  full <- NULL; a <- NA_real_; pf <- NA_character_
  for (cand in prefixes) {
    full_f <- file.path(BUILDS, paste0(cand, "_full.rds"))
    if (!file.exists(full_f)) next
    f <- readRDS(full_f); a <- ari(f, pipe)
    if (a >= RULES$rebuild_ari) { full <- f; pf <- cand; break }
    cat(sprintf("  %s: rebuild %s does not reproduce the network (ARI %.3f); not used\n", ds, cand, a))
  }
  if (is.null(full)) return(tibble())
  cat(sprintf("  %s %s: stability from %s (ARI %.3f with the network)\n", v, ds, pf, a))
  map <- tibble(pg = names(full), colour = unname(full)) |> inner_join(transmute(long, pg = ProteinGroupID, ModuleID), by = "pg") |>
    filter(colour != "grey", is_module(ModuleID)) |> count(colour, ModuleID) |> group_by(colour) |> slice_max(n, n = 1, with_ties = FALSE) |> ungroup()
  stopifnot(!anyDuplicated(map$ModuleID))
  loao <- list.files(BUILDS, pattern = sprintf("^%s_loao_A[0-9]+[.]rds$", pf), full.names = TRUE)
  lo_all <- lapply(loao, readRDS)
  bind_rows(lapply(seq_len(nrow(map)), function(i) {
    mod_pg <- names(full)[full == map$colour[i]]
    js <- vapply(lo_all, function(lo) {
      shared <- intersect(mod_pg, names(lo)); if (!length(shared)) return(NA_real_)
      cand <- setdiff(unique(lo[shared]), "grey"); if (!length(cand)) return(0)
      max(vapply(cand, function(cc) jac(shared, intersect(names(lo)[lo == cc], names(full))), numeric(1)))
    }, numeric(1))
    tibble(ModuleID = map$ModuleID[i], stability_median_jaccard = median(js, na.rm = TRUE), stability_min_jaccard = min(js, na.rm = TRUE),
           stability_n_recovered = sum(js >= RULES$stable_jaccard, na.rm = TRUE), stability_n_rebuilds = sum(!is.na(js)), stability_source = pf)
  }))
}

mods <- list(); sigs <- list(); sigs_all <- list(); terms <- list(); sms <- list(); stabs <- list()
for (v in names(VERSIONS)) for (ds in DATASETS) {
  cat("==", v, ds, "\n")
  s <- readRDS(VERSIONS[[v]]$state(ds))
  evidence_seen <- unique(c(evidence_seen, state_evidence(s)))
  if (length(evidence_seen) > 1L) stop("The two networks carry different gene annotation or GO limits (", paste(evidence_seen, collapse = " vs "),
                                       "); regenerate both with the same code before comparing them.", call. = FALSE)
  long <- as.data.frame(s$WGCNA_modules_long)
  go <- as.data.frame(s$GO_enrichment) |> filter(ModuleProteinSetType == "all")
  modules <- sort(unique(long$ModuleID[is_module(long$ModuleID)]))
  stabs[[length(stabs) + 1L]] <- stability(ds, v, VERSIONS[[v]]$estimator, long) |> mutate(version = v, dataset = ds)
  for (m in modules) {
    mm <- long[long$ModuleID == m, , drop = FALSE]
    mapped <- unique(na.omit(as.character(mm$EntrezID)))
    hubs <- mm |> filter(!is.na(EntrezID)) |> arrange(desc(abs_kME)) |> distinct(EntrezID, .keep_all = TRUE) |> head(RULES$hub_n)
    sg_all <- module_signals(go[go$ModuleID == m, , drop = FALSE])
    cc_m <- go |> filter(ModuleID == m, Ontology == "CC", p.adjust < RULES$q) |> arrange(p.adjust) |> head(RULES$top_terms)
    cc <- list(term = cc_m$Description, genes = lapply(strsplit(as.character(cc_m$geneID), "/", fixed = TRUE), unique))
    sg <- group_by_compartment(module_signals(go[go$ModuleID == m, , drop = FALSE], onts = "BP"), cc)
    basis <- "BP"
    if ((!nrow(sg) || !any(sg$passing)) && nrow(sg_all) && any(sg_all$passing)) {   # location- or complex-defined module
      sg <- group_by_compartment(sg_all, cc); basis <- "all ontologies (no passing BP signal)"
    }
    cl <- classify_module(sg, length(mapped), as.character(hubs$EntrezID), nrow(sg_all))
    dom <- if (!is.na(cl$dominant)) sg[cl$dominant, ] else NULL
    weak <- if (cl$signal_class == "unannotated (weak)") sg_all[1, ] else NULL
    comp <- if (length(cl$comparable)) sg[cl$comparable, ] else sg[0, ]
    mods[[length(mods) + 1L]] <- tibble(
      version = v, dataset = ds, ModuleID = m, n_proteins = nrow(mm), n_gene_mapped = length(mapped),
      n_signals = cl$n_signals, n_passing_signals = cl$n_passing, signal_basis = if (grepl("^unannotated", cl$signal_class)) "" else basis,
      signal_class = cl$signal_class, rule_confidence = cl$confidence,
      dominant_signal = if (!is.null(dom)) dom$top_terms else "", dominant_compartment = if (!is.null(dom)) coalesce(dom$compartment, "") else "",
      second_compartment = if (nrow(comp) > 1) coalesce(comp$compartment[2], "") else "", weak_association = if (!is.null(weak)) weak$top_terms else "",
      weak_ontologies = if (!is.null(weak)) weak$ontologies else "", weak_size = if (!is.null(weak)) weak$best_count else NA_integer_,
      dominant_size = if (!is.null(dom)) dom$best_count else NA_integer_, dominant_best_term = if (!is.null(dom)) dom$best_term else "",
      module_ontologies = paste(sort(unique(go$Ontology[go$ModuleID == m & go$p.adjust < RULES$q])), collapse = "+"),
      second_signal = if (nrow(comp) > 1) comp$top_terms[2] else "", third_signal = if (nrow(comp) > 2) comp$top_terms[3] else "",
      other_passing_signals = if (nrow(sg)) paste(sg$top_terms[sg$passing & !(sg$signal %in% c(cl$dominant, cl$comparable))], collapse = " || ") else "",
      second_to_dominant_ratio = cl$second_ratio %||% NA_real_,
      coverage = cl$coverage %||% NA_real_, hub_support = cl$hub_support %||% NA_real_, complex_exception = if (!is.null(dom)) dom$complex else NA,
      hubs = paste(coalesce(hubs$GeneSymbol, ""), collapse = ", "),
      dominant_term_ids = if (!is.null(dom)) paste(dom$term_ids[[1]], collapse = ";") else "",
      dominant_leading_id = if (!is.null(dom)) dom$leading_id else "",
      comparable_term_ids = if (nrow(comp)) paste(unlist(comp$term_ids), collapse = ";") else "")
    sig_rows <- function(x, set, comparable) x |> transmute(version = v, dataset = ds, ModuleID = m, ontology_set = set, signal, passing, comparable = signal %in% comparable,
      leading_term, leading_ontology, leading_q, leading_count, best_term, best_count, n_proteins, compartment = coalesce(if ("compartment" %in% names(x)) compartment else NA_character_, ""),
      n_terms, ontologies, top_terms, complex,
      fails_because = ifelse(passing, "", trimws(paste(ifelse(n_module_ontologies < RULES$min_ontologies, sprintf("module has significant terms in only %d ontology", n_module_ontologies), ""),
                                                       ifelse(leading_count <= RULES$min_proteins, sprintf("leading term carried by %d proteins (needs more than %d)", leading_count, RULES$min_proteins), ""),
                                                       sep = "; "), whitespace = "[; ]")))
    term_rows <- function(x, set) bind_rows(lapply(seq_len(nrow(x)), function(i) tibble(version = v, dataset = ds, ModuleID = m, ontology_set = set, signal = x$signal[i], term_id = x$term_ids[[i]])))
    if (nrow(sg)) { sigs[[length(sigs) + 1L]] <- sig_rows(sg, basis, cl$comparable); terms[[length(terms) + 1L]] <- term_rows(sg, basis) }
    # unannotated modules: the grouping over all three ontologies decides 'weak' (exactly one signal) or plain 'unannotated'
    if (grepl("^unannotated", cl$signal_class) && nrow(sg_all) && basis == "BP") {
      sigs_all[[length(sigs_all) + 1L]] <- sig_rows(sg_all, "all ontologies (decides weak or unannotated)", integer())
      terms[[length(terms) + 1L]] <- term_rows(sg_all, "all ontologies (decides weak or unannotated)") }
  }
  mods_v <- bind_rows(mods) |> filter(version == v, dataset == ds)
  smt <- read_csv(VERSIONS[[v]]$sm(ds), show_col_types = FALSE)
  for (i in seq_len(nrow(smt))) {
    mem <- strsplit(smt$member_modules[i], ";", fixed = TRUE)[[1]]
    mr <- mods_v |> filter(ModuleID %in% mem)
    labelled <- !grepl("^unannotated", mr$signal_class)
    cls <- if (length(mem) == 1) { if (labelled) "singleton" else mr$signal_class } else
      if (!any(labelled)) "unannotated" else if (!all(labelled)) "mixed" else {
        pj <- combn(nrow(mr), 2, function(ix) jac(ids_of(mr$dominant_term_ids[ix[1]]), ids_of(mr$dominant_term_ids[ix[2]])))
        if (all(pj >= RULES$term_jaccard)) "concordant" else "compound" }
    conf <- if (length(mem) == 1) mr$rule_confidence else if (cls == "mixed") "low" else if (cls == "unannotated") {
      if (any(mr$rule_confidence == "low")) "low" else "none" } else c("low", "moderate", "high")[min(match(mr$rule_confidence, c("low", "moderate", "high")))]
    sms[[length(sms) + 1L]] <- tibble(version = v, dataset = ds, SupermoduleID = smt$SupermoduleID[i], member_modules = smt$member_modules[i],
      members = paste(sprintf("%s (%d): %s, %s", sub("WGCNA_", "", mr$ModuleID), mr$n_proteins, mr$signal_class, mr$rule_confidence), collapse = "; "),
      member_dominant_signals = paste(sprintf("%s: %s", sub("WGCNA_", "", mr$ModuleID), ifelse(nzchar(mr$dominant_signal), sub(" [|].*", "", mr$dominant_signal), "none")), collapse = "; "),
      rule_class = cls, rule_confidence = conf, pipeline_label = smt$Supermodule_FinalLabel[i])
  }
}
mod_tab <- bind_rows(mods) |> left_join(bind_rows(stabs), by = c("version", "dataset", "ModuleID"))

# lineage between the versions (map_wgcna_module_lineage): new -> main original source; original -> main new destination
lin <- bind_rows(lapply(DATASETS, function(ds) read_csv(file.path("results", "wgcna", "map_wgcna_module_lineage", ds, "tables", "module_lineage_pairs.csv"), show_col_types = FALSE)))
to_new <- lin |> group_by(dataset, old_module) |> slice_max(share_of_old, n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(version = "original", dataset, ModuleID = old_module, lineage_module = new_module, lineage_share = share_of_old)
from_old <- lin |> group_by(dataset, new_module) |> slice_max(share_of_new, n = 1, with_ties = FALSE) |> ungroup() |>
  transmute(version = "new", dataset, ModuleID = new_module, lineage_module = old_module, lineage_share = share_of_new)
mod_tab <- mod_tab |> left_join(bind_rows(to_new, from_old), by = c("version", "dataset", "ModuleID"))
prev <- mod_tab |> filter(version == "original") |> select(dataset, lineage_module = ModuleID, prev_ids = dominant_term_ids, prev_lead = dominant_leading_id)
new_tab <- mod_tab |> filter(version == "new") |> left_join(prev, by = c("dataset", "lineage_module"), relationship = "many-to-one") |>
  mutate(found = mapply(function(ids, pids, plead) { a <- ids_of(ids); b <- ids_of(pids)
                          (nzchar(coalesce(plead, "")) && plead %in% a) || (length(b) > 0 && length(a) > 0 && jac(a, b) >= RULES$term_jaccard) },
                        comparable_term_ids, prev_ids, prev_lead),
         lineage_class = case_when(grepl("^unannotated", signal_class) ~ NA_character_, lineage_share < RULES$carry_share ~ "new",
                                   found ~ "carried", TRUE ~ "revised")) |>
  select(-prev_ids, -prev_lead, -found)
mod_tab <- bind_rows(mod_tab |> filter(version == "original") |> mutate(lineage_class = NA_character_), new_tab)
stopifnot(!anyDuplicated(mod_tab[c("version", "dataset", "ModuleID")]))
write_csv(mod_tab, file.path(OUT, "module_rule_classes.csv"), na = "")
write_csv(bind_rows(sigs), file.path(OUT, "module_signals.csv"), na = "")
write_csv(bind_rows(sigs_all), file.path(OUT, "module_signals_all_ontologies.csv"), na = "")
write_csv(bind_rows(terms), file.path(OUT, "signal_terms.csv"), na = "")
write_csv(bind_rows(sms), file.path(OUT, "supermodule_rule_classes.csv"), na = "")
yaml::write_yaml(RULES, file.path(OUT, "rules.yml"))
print(mod_tab |> count(version, signal_class) |> pivot_wider(names_from = signal_class, values_from = n, values_fill = 0))
cat("done:", OUT, "\n")
