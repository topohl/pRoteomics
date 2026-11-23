#!/usr/bin/env Rscript
#
#
# Which WGCNA module carries each Figure 2 home term, and does that module peak where the term does (added 2026-11-23;
# Figure 4 candidate).
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/link_home_terms_to_modules.R
#   Rscript analysis/spatial_validation/link_home_terms_to_modules.R --dry-run
# Script: analysis/spatial_validation/link_home_terms_to_modules.R
# Stage: networks
# Scope: global
# Consumes: required results/spatial_validation/map_anatomical_programmes/global/tables/anatomical_programme_display_terms.csv; results/wgcna/build_wgcna_modules/<dataset>/tables/modules/WGCNA_modules_long.csv; results/wgcna/build_wgcna_modules/<dataset>/models/wgcna_final_model_state.rds; results/wgcna/test_module_phenotypes/<dataset>/tables/WGCNA_group_effect_hemisphere_values.csv; results/spatial_validation/map_input_territories/global/tables/input_territory_carryover_qc.csv
# Produces: results/spatial_validation/link_home_terms_to_modules/global/tables/home_term_module_links.csv; home_term_module_overlap.csv; home_term_link_summary.csv; home_term_base_rate.csv
# Dataset behavior: global; neuron_soma (Figure 2f) and neuron_neuropil (Figure 2g-i).
# Notes: Descriptive and phenotype-blind: control animals only, no group effect, no WGCNA refit.
#
# HOME TERMS: the curated GO terms of Figure 2f-i in the unit they were picked for (map_anatomical_programmes,
# selection 'curated', unit = selected_for_unit).
#
# CARRIER: for each home term and module of its compartment, a one-sided Fisher test of the term's proteins in the
# module (GO BP membership with descendants, org.Mm.eg.db; universe = the module proteins with an Entrez ID);
# Benjamini-Hochberg over the modules of the term. Carriers have FDR <= 0.05 and odds ratio >= 2 (conditional maximum
# likelihood); the best carrier holds the largest share of the term's proteins.
#
# LINK: the best carrier's eigengene peaks at the home unit, among the term's sibling units, in every control animal
# (mean of the animal's hemispheres). The DG po acquisitions with tissue carry-over are left out as in Figure 2j-l.
#
# TWO CHECKS on every link:
#   restating: the module's proteins outside the term (unannotated ones included), averaged as per-protein z over the
#     control acquisitions, still peak at home in every control animal (peaks_without_term_proteins);
#   selection-matched null: random protein sets of the term's size with the same number of proteins that peak at home
#     in every control animal as the term itself (a term picked because it is high at home); the carrier rule (with the
#     sample odds ratio, 0.5 added to each cell) and the link are applied to 400 sets per term. null_p_link is the share
#     of sets that link; a link with null_p_link near 1 restates how its term was picked.
# The summary's tail P treats the terms as independent; several home terms are nested or repeated (the same GO term
# in two units), so links are correlated and that P is too small. Read it as descriptive.
#
#

source("R/paths.R")
source("R/statistics/integration_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))
source("R/spatial/neuropil_lamina_utils.R")

ANALYSIS_ID <- "link_home_terms_to_modules"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)
SCRIPT_ID <- "analysis/spatial_validation/link_home_terms_to_modules.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

DATASETS <- c("neuron_soma", "neuron_neuropil")
TERMS <- repo_path("results", "spatial_validation", "map_anatomical_programmes", "global", "tables", "anatomical_programme_display_terms.csv")
CARRYOVER <- repo_path("results", "spatial_validation", "map_input_territories", "global", "tables", "input_territory_carryover_qc.csv")
MODULES <- stats::setNames(vapply(DATASETS, function(d) repo_path("results", "wgcna", "build_wgcna_modules", d, "tables", "modules", "WGCNA_modules_long.csv"), character(1)), DATASETS)
STATE <- stats::setNames(vapply(DATASETS, function(d) repo_path("results", "wgcna", "build_wgcna_modules", d, "models", "wgcna_final_model_state.rds"), character(1)), DATASETS)
HEMI <- stats::setNames(vapply(DATASETS, function(d) repo_path("results", "wgcna", "test_module_phenotypes", d, "tables", "WGCNA_group_effect_hemisphere_values.csv"), character(1)), DATASETS)
tab <- function(f) file.path(CANONICAL_PATHS$tables, f)
OUT <- list(links = tab("home_term_module_links.csv"), overlap = tab("home_term_module_overlap.csv"),
            summary = tab("home_term_link_summary.csv"), base_rate = tab("home_term_base_rate.csv"))
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")

N_SIM <- 400L
SEED <- 20261030L
PANEL <- c(soma_regions = "f", CA1_strata = "g", CA2_strata = "h", DG_layers = "i")
SIBLINGS <- list(soma_regions = c("ca1", "ca2", "ca3", "dg"), CA1_strata = c("ca1_so", "ca1_sr", "ca1_slm"),
                 CA2_strata = c("ca2_so", "ca2_sr", "ca2_slm"), DG_layers = c("dg_mo", "dg_po"))
PREFIX <- c(CA1_strata = "ca1", CA2_strata = "ca2", DG_layers = "dg")

if (isTRUE(cli$dry_run)) {
  cat("[DRY-RUN] Figure 2 home terms against the WGCNA modules: carrier, link, restating check, selection-matched null.\n")
  dry_run_inputs(SCRIPT_ID, c(list(terms = TERMS, carryover = CARRYOVER), as.list(MODULES), as.list(STATE), as.list(HEMI)))
  quit(save = "no", status = 0L)
}
for (f in c(TERMS, CARRYOVER, MODULES, STATE, HEMI)) if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
for (pkg in c("AnnotationDbi", "org.Mm.eg.db"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required annotation package unavailable: ", pkg, call. = FALSE)
invisible(lapply(CANONICAL_PATHS[c("tables", "manifests")], dir_create))

# ------------------------------------------------------------------ home terms
dt <- utils::read.csv(TERMS, stringsAsFactors = FALSE)
home <- dt$selection == "curated" & dt$unit == dt$selected_for_unit
if (any(home & !dt$family %in% names(PANEL)))
  stop("A curated home term lies outside Figure 2f-i: ", paste(unique(dt$family[home & !dt$family %in% names(PANEL)]), collapse = ", "),
       call. = FALSE)
H <- dt[home, , drop = FALSE]
H$panel <- unname(PANEL[H$family])
H$home_unit <- ifelse(H$family == "soma_regions", tolower(H$unit), paste(PREFIX[H$family], tolower(H$unit), sep = "_"))
H$compartment <- ifelse(H$family == "soma_regions", "neuron_soma", "neuron_neuropil")
if (any(H$compartment != H$dataset)) stop("A home term's compartment differs from its dataset column.", call. = FALSE)
H <- H[order(H$panel, H$display_rank, H$ID), , drop = FALSE]
H$row <- seq_len(nrow(H))
if (anyDuplicated(H[c("family", "unit", "ID")])) stop("A home term is listed twice.", call. = FALSE)
go <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = unique(H$ID), keytype = "GOALL", columns = "ENTREZID"))
go <- unique(go[!is.na(go$ENTREZID), c("GOALL", "ENTREZID")])
term_genes <- split(as.character(go$ENTREZID), go$GOALL)
if (!all(H$ID %in% names(term_genes)))
  stop("Home terms without GO annotation: ", paste(setdiff(H$ID, names(term_genes)), collapse = ", "), call. = FALSE)

carry <- utils::read.csv(CARRYOVER, stringsAsFactors = FALSE)
left_out <- carry$Sample[carry$excluded %in% TRUE]
np_samples <- unique(utils::read.csv(HEMI[["neuron_neuropil"]], stringsAsFactors = FALSE)$contributing_sample_ids)
if (!length(left_out) || !all(left_out %in% np_samples))
  stop("The carry-over acquisitions to leave out are not all in the neuropil eigengene table.", call. = FALSE)

# ------------------------------------------------------------------ per compartment
peak_count <- function(prof, home, sibs, animals) {
  sum(vapply(animals, function(a) { p <- prof[prof$animal == a & prof$unit %in% sibs, , drop = FALSE]
    if (nrow(p) != length(sibs)) stop("A control animal lacks a sibling unit.", call. = FALSE)
    p$unit[which.max(p$value)] == home }, logical(1)))
}
ctx <- list()
for (d in DATASETS) {
  L <- utils::read.csv(MODULES[[d]], stringsAsFactors = FALSE, colClasses = c(official_entrez_id = "character"))[c("ProteinGroupID", "ModuleID", "official_entrez_id")]
  L$official_entrez_id[is.na(L$official_entrez_id) | L$official_entrez_id %in% c("", "NA")] <- NA_character_
  h <- utils::read.csv(HEMI[[d]], stringsAsFactors = FALSE)
  h <- h[h$level == "module", , drop = FALSE]
  h$SpatialUnit <- tolower(h$SpatialUnit)
  animals <- sort(unique(h$AnimalID[h$StressGroup == "CON"]))
  hc <- h[h$StressGroup == "CON" & !h$contributing_sample_ids %in% left_out, , drop = FALSE]
  ME <- stats::aggregate(hemisphere_value ~ endpoint_id + AnimalID + SpatialUnit, hc, mean)
  names(ME) <- c("endpoint_id", "animal", "unit", "value")
  # acquisitions of the control animals, as the hemisphere table names them
  acq <- unique(hc[c("contributing_sample_ids", "AnimalID", "SpatialUnit", "Hemisphere")])
  s <- readRDS(STATE[[d]])
  X <- as.matrix(s$expression.data)
  if (!all(acq$contributing_sample_ids %in% rownames(X))) stop(d, ": a control acquisition is missing from the model state.", call. = FALSE)
  Xc <- X[acq$contributing_sample_ids, , drop = FALSE]
  if (!all(L$ProteinGroupID %in% colnames(X))) stop(d, ": a module protein is missing from the model state.", call. = FALSE)
  # the stored eigengenes are those of the model state (one acquisition per hemisphere value)
  me <- as.matrix(s$mergedMEs)
  mm <- unique(as.data.frame(s$WGCNA_modules_long)[c("WGCNAInternalColor", "ModuleID")])
  colnames(me) <- mm$ModuleID[match(sub("^ME", "", colnames(me)), mm$WGCNAInternalColor)]
  key <- match(paste(hc$contributing_sample_ids, hc$endpoint_id), paste(rep(rownames(me), ncol(me)), rep(colnames(me), each = nrow(me))))
  me_diff <- max(abs(as.vector(me)[key] - hc$hemisphere_value))
  if (!is.finite(me_diff) || me_diff > 1e-6) stop(d, ": stored hemisphere values differ from the model-state eigengenes.", call. = FALSE)
  Zc <- scale(Xc)   # per-protein z over the control acquisitions
  ctx[[d]] <- list(L = L, ME = ME, animals = animals, acq = acq, Zc = Zc, me_diff = me_diff, mods = sort(unique(L$ModuleID)))
}
animals <- ctx[[1]]$animals
if (!identical(animals, ctx[[2]]$animals)) stop("The control animals differ between compartments.", call. = FALSE)
n_animals <- length(animals)

# per-protein count of control animals in which the protein peaks at home among the siblings
protein_peak <- function(cx, cols, home, sibs) {
  z <- cx$Zc[, cols, drop = FALSE]; cnt <- integer(length(cols))
  for (a in cx$animals) {
    M <- vapply(sibs, function(u) colMeans(z[cx$acq$AnimalID == a & cx$acq$SpatialUnit == u, , drop = FALSE]), numeric(length(cols)))
    if (is.null(dim(M))) M <- matrix(M, nrow = 1L)
    cnt <- cnt + (max.col(M, ties.method = "first") == match(home, sibs))
  }
  cnt
}
set_profile <- function(cx, cols) {
  v <- rowMeans(cx$Zc[, cols, drop = FALSE])
  a <- stats::aggregate(v ~ AnimalID + SpatialUnit, data.frame(v = v, AnimalID = cx$acq$AnimalID, SpatialUnit = cx$acq$SpatialUnit), mean)
  data.frame(animal = a$AnimalID, unit = a$SpatialUnit, value = a$v)
}

overlap <- list(); links <- list(); base <- list(); nulls <- list()
set.seed(SEED)
for (k in H$row) {
  t <- H[H$row == k, ]; cx <- ctx[[t$compartment]]; sibs <- SIBLINGS[[t$family]]; home <- t$home_unit
  if (!home %in% sibs || !all(sibs %in% cx$acq$SpatialUnit)) stop("Home unit or siblings not found: ", home, call. = FALSE)
  Lu <- cx$L[!is.na(cx$L$official_entrez_id), , drop = FALSE]
  inT <- Lu$official_entrez_id %in% term_genes[[t$ID]]
  modpeak <- vapply(cx$mods, function(m) peak_count(cx$ME[cx$ME$endpoint_id == m, ], home, sibs, cx$animals), integer(1))
  f <- do.call(rbind, lapply(cx$mods, function(m) {
    inM <- Lu$ModuleID == m
    tb <- matrix(c(sum(inT & inM), sum(inT & !inM), sum(!inT & inM), sum(!inT & !inM)), 2)
    ft <- stats::fisher.test(tb, alternative = "greater")
    data.frame(row = k, dataset = t$compartment, spatial_unit = home, GO_ID = t$ID, module_id = m, n_term = sum(inT), n_overlap = tb[1, 1],
               share_of_term = tb[1, 1] / sum(inT), odds_ratio = unname(ft$estimate), p_value = ft$p.value,
               animals_peak_home = modpeak[[m]], stringsAsFactors = FALSE)
  }))
  f$fdr_within_term <- stats::p.adjust(f$p_value, "BH")
  f$carrier <- f$fdr_within_term <= 0.05 & f$odds_ratio >= 2
  overlap[[k]] <- f
  car <- f[f$carrier, , drop = FALSE]; car <- car[order(-car$share_of_term, car$p_value), , drop = FALSE]
  best <- if (nrow(car)) car$module_id[1] else NA_character_
  best_peak <- if (nrow(car)) car$animals_peak_home[1] else NA_integer_
  linked <- isTRUE(best_peak == n_animals)
  # restating check: the carrier's proteins outside the term
  nonterm_peak <- NA_integer_; n_nonterm <- NA_integer_
  if (linked) {
    memb <- cx$L[cx$L$ModuleID == best, , drop = FALSE]
    nonterm <- memb$ProteinGroupID[is.na(memb$official_entrez_id) | !memb$official_entrez_id %in% term_genes[[t$ID]]]
    n_nonterm <- length(nonterm)
    nonterm_peak <- peak_count(set_profile(cx, nonterm), home, sibs, cx$animals)
  }
  termset_peak <- peak_count(set_profile(cx, Lu$ProteinGroupID[inT]), home, sibs, cx$animals)
  # selection-matched null
  hp <- protein_peak(cx, Lu$ProteinGroupID, home, sibs) == n_animals
  modv <- match(Lu$ModuleID, cx$mods); K <- tabulate(modv, length(cx$mods)); N <- nrow(Lu)
  rule <- function(idx) {
    a <- tabulate(modv[idx], length(cx$mods)); n <- length(idx)
    p <- stats::phyper(a - 1, K, N - K, n, lower.tail = FALSE)
    b <- n - a; c <- K - a; dd <- N - K - b
    or <- ((a + 0.5) * (dd + 0.5)) / ((b + 0.5) * (c + 0.5))
    ok <- stats::p.adjust(p, "BH") <= 0.05 & or >= 2
    if (!any(ok)) return(c(carrier = 0, link = 0))
    bi <- which(ok)[which.max((a / n)[ok])]
    c(carrier = 1, link = as.numeric(modpeak[[bi]] == n_animals))
  }
  n <- sum(inT); nh <- sum(inT & hp)
  obs_null_rule <- rule(which(inT))
  sim <- t(replicate(N_SIM, rule(c(sample(which(hp), nh), sample(which(!hp), n - nh)))))
  links[[k]] <- data.frame(row = k, figure2_panel = t$panel, family = t$family, dataset = t$compartment, spatial_unit = home,
    GO_ID = t$ID, GO_term = t$Description, NES = t$NES, n_term = n, n_carriers = nrow(car),
    carriers = paste(sprintf("%s (share %.2f, OR %.1f, FDR %.1e, %d/%d animals)", car$module_id, car$share_of_term, car$odds_ratio,
                             car$fdr_within_term, car$animals_peak_home, n_animals), collapse = "; "),
    module_id = best, animals_peak_home = best_peak, linked = linked, n_nonterm_proteins = n_nonterm,
    nonterm_animals_peak_home = nonterm_peak, peaks_without_term_proteins = linked & isTRUE(nonterm_peak == n_animals),
    termset_animals_peak_home = termset_peak, frac_term_proteins_peak_home = nh / n, frac_universe_peak_home = mean(hp),
    n_modules_peak_home = sum(modpeak == n_animals), linked_under_null_rule = obs_null_rule[["link"]] == 1,
    null_p_carrier = mean(sim[, "carrier"]), null_p_link = mean(sim[, "link"]), n_null_sets = N_SIM, stringsAsFactors = FALSE)
  bkey <- paste(t$compartment, t$family, home)
  if (!bkey %in% names(base))
    base[[bkey]] <- data.frame(dataset = t$compartment, family = t$family, spatial_unit = home, n_modules = length(cx$mods),
                               n_modules_peak_home = sum(modpeak == n_animals),
                               modules_peak_home = paste(cx$mods[modpeak == n_animals], collapse = ";"), stringsAsFactors = FALSE)
}
links <- do.call(rbind, links); overlap <- do.call(rbind, overlap); base <- do.call(rbind, base)
poisson_binomial_tail <- function(p, k) { dd <- 1; for (q in p) dd <- c(dd * (1 - q), 0) + c(0, dd * q); sum(dd[(k + 1):length(dd)]) }
summ <- do.call(rbind, lapply(c(DATASETS, "all"), function(d) {
  x <- if (d == "all") links else links[links$dataset == d, , drop = FALSE]
  data.frame(dataset = d, n_terms = nrow(x), n_with_carrier = sum(x$n_carriers > 0), n_linked = sum(x$linked),
             n_linked_without_term_proteins = sum(x$peaks_without_term_proteins), expected_links_null = sum(x$null_p_link),
             p_links_at_least_observed = poisson_binomial_tail(x$null_p_link, sum(x$linked_under_null_rule)),
             n_linked_under_null_rule = sum(x$linked_under_null_rule), n_control_animals = n_animals,
             note = "tail P treats the terms as independent; nested and repeated home terms make it too small (descriptive)",
             stringsAsFactors = FALSE)
}))
utils::write.csv(links, OUT$links, row.names = FALSE)
utils::write.csv(overlap, OUT$overlap, row.names = FALSE)
utils::write.csv(summ, OUT$summary, row.names = FALSE)
utils::write.csv(base, OUT$base_rate, row.names = FALSE)
write_run_manifest(MANIFEST, inputs = c(list(terms = TERMS, carryover = CARRYOVER), as.list(MODULES), as.list(STATE), as.list(HEMI)),
                   outputs = OUT,
                   parameters = list(control_animals = animals, left_out_acquisitions = left_out, null_sets = N_SIM, seed = SEED,
                                     annotation = paste("org.Mm.eg.db", as.character(utils::packageVersion("org.Mm.eg.db"))),
                                     max_eigengene_difference = lapply(ctx, `[[`, "me_diff"),
                                     provenance = run_output_provenance(OUT, SCRIPT_ID)),
                   notes = "Figure 2 home terms against the WGCNA modules for the Figure 4 candidate; phenotype-blind.")
message("Home-term links: ", nrow(links), " terms, ", sum(links$linked), " linked (", sum(links$peaks_without_term_proteins),
        " without their term proteins); expected under the selection-matched null ", round(sum(links$null_p_link), 1),
        ", P = ", signif(summ$p_links_at_least_observed[summ$dataset == "all"], 2))
