#!/usr/bin/env Rscript
#
#
# Which GO biological-process programmes set each hippocampal region and layer
# apart, in control animals, and how well is each one supported?
#
# WHAT IS ADDED (2026-10-27; NOT among the 11 prespecified contrasts)
#   neuropil:  CA2 strata (SO, SR, SLM, each vs the mean of the other CA2 strata)
#              CA3 strata (SO vs SR, in both directions)
#   microglia-enriched ROI: each region vs the mean of the other regions (the
#              registry's microglia contrasts; regional context only)
#   direct pairwise contrasts between sibling units (soma regions, CA1 strata,
#              CA2 strata; two-unit families already are pairwise), so that no
#              programme is called unit-specific from a difference in significance
#   one-vs-rest contrasts of the families shown in Figure 2 fitted on left-only
#              and on right-only samples (direction consistency; same animals)
# Model, ranking and GO-BP GSEA are those of
# analysis/differential_abundance/validate_control_spatial_identity.R: CON
# animals only; abundance ~ 0 + anatomical_unit + hemisphere (one-sided fits drop
# the hemisphere term); animal as the duplicateCorrelation block; robust trend
# eBayes; the same seed rule, permutations and term sizes. A locked neuropil
# contrast (CA1 SLM) is re-run with its GSEA and must match the locked table; the
# locked soma contrast CA1 must match at protein level; the microglia estimates
# must match quantify_bilateral_spatial_identity.R. A mismatch stops the run.
#
# BILATERAL AGREEMENT: every added contrast is also estimated on left-only and
# right-only samples (Pearson r of the per-protein estimates, as in
# quantify_bilateral_spatial_identity.R). Same animals: not replication.
#
# DISPLAY SELECTION (presentation only, from completed results). Three selections
# are written side by side so that none is hidden:
#   automatic: anatomical_programme_rule() - supported in one unit and not in its
#              siblings, top 3 by NES (kept as a sensitivity comparison)
#   strict:    anatomical_programme_strict_rule() - direct pairwise support against
#              every sibling, both hemispheres positive, several contributing
#              proteins, no single-protein dominance; top 3 by NES
#   curated:   config/anatomical_programme_display_picks.csv - the authors' picks
#              with a reason each, checked against every criterion
# Every positive term of every unit carries the metrics behind these criteria
# (anatomical_programme_candidates.csv).
#
# USAGE (from the repository root)
#   Rscript analysis/spatial_validation/map_anatomical_programmes.R
#   Rscript analysis/spatial_validation/map_anatomical_programmes.R --dry-run
#   PROTEOMICS_ANATOMICAL_PROGRAMMES_FIT=pairwise,hemisphere Rscript ...  (refit only those stages)
#   PROTEOMICS_ANATOMICAL_PROGRAMMES_SELECT_ONLY=true Rscript ...        (re-select and re-plot only)
# Script: analysis/spatial_validation/map_anatomical_programmes.R
# Stage: networks
# Scope: global
# Consumes: required results/preprocessing/build_module_score_metadata/<dataset>/tables/sample_metadata_merged_clean_for_module_scores.xlsx; results/differential_abundance/validate_control_spatial_identity/global/tables/control_anatomical_go_bp_gsea.csv; results/differential_abundance/validate_control_spatial_identity/global/tables/anatomical_protein_contrasts.csv; results/spatial_validation/quantify_bilateral_spatial_identity/global/tables/bilateral_spatial_identity_protein_level.csv; config/clusterProfiler_config.yml; config/anatomical_programme_display_picks.csv
# Produces: results/spatial_validation/map_anatomical_programmes/global/tables/added_anatomical_go_bp_gsea.csv; results/spatial_validation/map_anatomical_programmes/global/tables/anatomical_programme_display_terms.csv; +12 more
# Dataset behavior: fits neuron_soma, neuron_neuropil and microglia; selection reads the locked soma and neuropil GSEA.
# Notes: Control-only anatomical GO programmes beyond the locked contrasts, their robustness, and the Figure 2 display selections.
#
#

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/data_contracts/dataset_inputs.R")
source("R/statistics/integration_utils.R")
source("R/utilities/validation_utils.R")
source("R/qc/qc_exploration_utils.R")
source("R/enrichment/protein_group_enrichment_utils.R")
source("R/enrichment/clusterprofiler_reproducibility.R")
source("R/spatial/control_spatial_identity_utils.R")
source("R/data_contracts/spatial_systems_data_utils.R")
source("R/spatial/anatomical_programme_utils.R")
source(repo_path("R", "spatial_systems_paths.R"))
source(repo_path("R", "differential_abundance_paths.R"))

ANALYSIS_ID <- "map_anatomical_programmes"
CANONICAL_PATHS <- spatial_systems_dirs(ANALYSIS_ID)

SCRIPT_ID <- "analysis/spatial_validation/map_anatomical_programmes.R"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)
cli <- integration_cli(default_dataset = "all")

LOCKED_DIR <- differential_abundance_dirs("validate_control_spatial_identity")$tables
LOCKED <- c(go = file.path(LOCKED_DIR, "control_anatomical_go_bp_gsea.csv"),
            protein = file.path(LOCKED_DIR, "anatomical_protein_contrasts.csv"))
BILATERAL <- spatial_systems_find("bilateral_spatial_identity_protein_level.csv", "bilateral")
GSEA_CONFIG <- repo_path("config", "clusterProfiler_config.yml")
PICKS <- repo_path("config", "anatomical_programme_display_picks.csv")
CHECK_CONTRAST <- c(neuron_neuropil = "CA1_SLM_vs_mean_other_CA1_strata",
                    neuron_soma = "CA1_vs_mean_other_soma_regions")
check_for <- function(ds) if (ds %in% names(CHECK_CONTRAST)) unname(CHECK_CONTRAST[[ds]]) else NA_character_
FIT_DATASETS <- c("neuron_soma", "neuron_neuropil", "microglia")

TABLES <- function(name) file.path(CANONICAL_PATHS$tables, name)
PLOTS <- function(name) file.path(CANONICAL_PATHS$plots, name)
OUT <- list(
  protein     = TABLES("added_anatomical_protein_contrasts.csv"),
  go          = TABLES("added_anatomical_go_bp_gsea.csv"),
  pairwise_go = TABLES("pairwise_anatomical_go_bp_gsea.csv"),
  side_go     = TABLES("hemisphere_anatomical_go_bp_gsea.csv"),
  gene_stats  = TABLES("anatomical_gene_statistics.csv"),
  bilateral   = TABLES("added_contrast_bilateral_agreement.csv"),
  check       = TABLES("method_reproduction_check.csv"),
  status      = TABLES("analysis_status.csv"),
  candidates  = TABLES("anatomical_programme_candidates.csv"),
  display     = TABLES("anatomical_programme_display_terms.csv"),
  picks       = TABLES("anatomical_programme_picks_report.csv"),
  plot        = PLOTS("anatomical_programme_bubbles.svg"),
  plot_strict = PLOTS("anatomical_programme_bubbles_strict.svg"),
  plot_auto   = PLOTS("anatomical_programme_bubbles_automatic.svg")
)
MANIFEST <- file.path(CANONICAL_PATHS$manifests, "run_manifest.yml")

if (isTRUE(cli$dry_run)) {
  inputs <- c(as.list(LOCKED), bilateral = BILATERAL, gsea_config = GSEA_CONFIG, picks = PICKS,
              stats::setNames(lapply(FIT_DATASETS, preprocessing_module_score_metadata),
                              paste0("canonical_metadata__", FIT_DATASETS)))
  cat("[DRY-RUN] Control-only anatomical programmes (added, pairwise and one-sided contrasts; display selections).\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  quit(save = "no", status = 0L)
}

invisible(lapply(CANONICAL_PATHS[c("tables", "plots", "manifests")], dir_create))
suppressPackageStartupMessages({
  library(limma); library(clusterProfiler); library(org.Mm.eg.db); library(ggplot2)
})
for (f in c(LOCKED, BILATERAL, GSEA_CONFIG)) {
  if (!file.exists(f)) stop("missing_required_input: ", f, call. = FALSE)
}
gsea_cfg <- validate_clusterprofiler_gsea_config(yaml::read_yaml(GSEA_CONFIG))
GSEA <- list(
  seed_base = gsea_cfg$analysis$gsea_seed_base,
  n_perm = gsea_cfg$analysis$n_perm_simple,
  min_size = clusterprofiler_integer_scalar(gsea_cfg$analysis$min_gs_size, "min_gs_size"),
  max_size = clusterprofiler_integer_scalar(gsea_cfg$analysis$max_gs_size, "max_gs_size")
)
RULE <- anatomical_programme_rule()
STRICT <- anatomical_programme_strict_rule()
FAMILIES <- anatomical_programme_families()
# Cross-session GSEA tolerance documented in R/enrichment/clusterprofiler_reproducibility.R.
TOL <- list(estimate = 1e-8, nes = 1e-4, p = 1e-4)

select_only <- tolower(Sys.getenv("PROTEOMICS_ANATOMICAL_PROGRAMMES_SELECT_ONLY", "")) %in% c("1", "true", "yes")
fit_stages <- trimws(strsplit(Sys.getenv("PROTEOMICS_ANATOMICAL_PROGRAMMES_FIT", "added,pairwise,hemisphere"), ",")[[1]])
if (!all(fit_stages %in% c("added", "pairwise", "hemisphere"))) stop("Unknown fit stage: ", paste(fit_stages, collapse = ","), call. = FALSE)
if (select_only) fit_stages <- character()

# CON metadata exactly as the locked analysis prepares it. The canonical helper
# reads the hemisphere from the sample name, which does not resolve for the
# microglia ROI; there the ReplicateGroup value is used, as in
# quantify_bilateral_spatial_identity.R.
con_metadata <- function(meta, dataset) {
  x <- tryCatch(control_spatial_prepare_metadata(meta), error = function(e) e)
  if (!inherits(x, "error")) return(x)
  if (!identical(dataset, "microglia") || !grepl("hemisphere", conditionMessage(x))) stop(x)
  x <- meta[meta$StressGroup == "CON" & (is.na(meta$Exclude) | !meta$Exclude), , drop = FALSE]
  if (!nrow(x) || anyNA(x$AnimalID)) stop("No usable CON microglia samples.", call. = FALSE)
  x$hemisphere <- factor(as.character(sps_normalize_hemisphere(x[[sps_hemisphere_source_field()]])),
                         levels = c("L", "R"))
  if (anyNA(x$hemisphere)) stop("Unresolved microglia hemisphere.", call. = FALSE)
  x$anatomical_unit <- ifelse(x$SpatialUnit == "region_layer",
                              paste(toupper(x$Region), toupper(x$Layer), sep = "_"), toupper(x$Region))
  x
}

fit_units <- function(mat, meta) {
  d <- control_spatial_design(meta)
  corfit <- limma::duplicateCorrelation(mat, d$design, block = meta$AnimalID)
  fit <- limma::lmFit(mat, d$design, block = meta$AnimalID, correlation = corfit$consensus)
  list(fit = fit, design = d, correlation = corfit$consensus)
}

contrast_table <- function(model, spec, feature_table, dataset) {
  w <- control_spatial_contrast_vector(spec, colnames(model$design$design))
  cf <- limma::eBayes(limma::contrasts.fit(model$fit, matrix(w, ncol = 1, dimnames = list(names(w), spec$name))),
                      robust = TRUE, trend = TRUE)
  tt <- limma::topTable(cf, number = Inf, sort.by = "none")
  tt$ProteinGroupID <- rownames(tt)
  tt <- cbind(tt, feature_table[match(tt$ProteinGroupID, feature_table$ProteinGroupID),
                                setdiff(names(feature_table), "ProteinGroupID"), drop = FALSE])
  tt$dataset <- dataset; tt$contrast <- spec$name
  tt
}

# Each GSEA result is kept in work/ as soon as it exists, keyed by the exact
# ranked statistics and the GSEA settings, so a crash of this long run loses at
# most the contrast in progress (2026-10-27: a native crash after 36 GSEA runs
# lost all of them). A cached result is used only when the ranking, the seed
# identity and every setting are identical, so it equals a recomputation.
GSEA_CACHE <- file.path(CANONICAL_PATHS$work, "gsea_cache")
run_go <- function(dataset, contrast, ranked, seed_contrast = contrast) {
  repro <- control_spatial_gsea_reproducibility(dataset, seed_contrast, "gseGO_BP", GSEA$seed_base, GSEA$n_perm)
  key <- list(ranked = ranked, repro = repro, settings = GSEA, dataset = dataset, contrast = contrast,
              packages = c(clusterProfiler = as.character(utils::packageVersion("clusterProfiler")),
                           fgsea = as.character(utils::packageVersion("fgsea")),
                           GO.db = as.character(utils::packageVersion("GO.db"))))
  cache_file <- file.path(GSEA_CACHE, paste0(digest::digest(key, algo = "sha256"), ".rds"))
  if (file.exists(cache_file)) {
    hit <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(hit) && identical(hit$key, key)) return(hit$go)
  }
  go <- as.data.frame(run_seeded_clusterprofiler_gsea(
    clusterProfiler::gseGO, gsea_seed = repro$gsea_seed, n_perm_simple = repro$n_perm_simple,
    geneList = ranked, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP",
    minGSSize = GSEA$min_size, maxGSSize = GSEA$max_size, pvalueCutoff = 1, verbose = FALSE))
  go <- if (!nrow(go)) {
    cbind(data.frame(dataset = dataset, contrast = contrast, status = "completed_zero_terms"), repro)
  } else {
    go$dataset <- dataset; go$contrast <- contrast; go$status <- "completed"
    go$p_value <- go$pvalue; go$p_adjust <- go$p.adjust; go$mapped_unique_genes <- go$setSize
    cbind(go, repro[rep(1L, nrow(go)), , drop = FALSE])
  }
  dir_create(GSEA_CACHE)
  tmp <- paste0(cache_file, ".tmp")
  saveRDS(list(key = key, go = go), tmp)
  file.rename(tmp, cache_file)
  go
}

gene_table <- function(tt) {
  g <- build_enrichment_gene_inputs(tt, strict = TRUE)
  list(ranked = g$ranked,
       stats = data.frame(contrast = tt$contrast[1], GeneSymbol = g$collapse$GeneSymbol,
                          collapsed_statistic = g$collapse$collapsed_statistic,
                          collapsed_fdr = g$collapse$collapsed_fdr, collapsed_logfc = g$collapse$collapsed_logfc,
                          stringsAsFactors = FALSE))
}

if (length(fit_stages)) {
  acc <- list(protein = list(), go = list(), pairwise = list(), side = list(), genes = list(),
              agreement = list(), checks = list(), status = list())
  models_info <- list()
  bil <- utils::read.csv(BILATERAL, check.names = FALSE)
  locked_protein <- utils::read.csv(LOCKED[["protein"]], check.names = FALSE)
  locked_go <- utils::read.csv(LOCKED[["go"]], check.names = FALSE)
  for (ds in FIT_DATASETS) {
    message("Loading ", ds)
    inputs <- resolve_dataset_inputs(ds, purpose = "wgcna", script = SCRIPT_ID, stage = "networks")
    canonical <- qc_load_canonical_expression(inputs$expression_file, preprocessing_module_score_metadata(ds),
                                              dataset = ds, strict = TRUE)
    meta <- con_metadata(canonical$meta, ds)
    meta <- meta[match(colnames(canonical$mat)[colnames(canonical$mat) %in% meta$Sample], meta$Sample), , drop = FALSE]
    mat <- canonical$mat[, meta$Sample, drop = FALSE]
    model <- fit_units(mat, meta)
    units <- sub("^anatomical_unit_", "", grep("^anatomical_unit_", colnames(model$design$design), value = TRUE))
    models_info[[ds]] <- list(dataset = ds, n_samples = ncol(mat), n_animals = length(unique(meta$AnimalID)),
                              hemisphere_in_model = model$design$hemisphere_included,
                              duplicate_correlation = model$correlation, units = units)
    fam_ds <- FAMILIES[FAMILIES$dataset == ds, , drop = FALSE]
    locked_specs <- control_spatial_contrast_registry(ds, units)
    added_specs <- if (ds %in% c("neuron_neuropil", "microglia")) anatomical_programme_added_specs(ds, units) else list()
    one_vs_rest <- c(locked_specs[intersect(names(locked_specs), fam_ds$contrast)], added_specs)
    one_vs_rest <- one_vs_rest[!duplicated(names(one_vs_rest))]  # microglia: the added specs are the registry's own

    # gene-level statistics of every family contrast (leading-edge support); locked ones refitted here
    for (nm in names(one_vs_rest)) {
      tt <- contrast_table(model, one_vs_rest[[nm]], canonical$feature_table, ds)
      acc$genes[[nm]] <- gene_table(tt)$stats
      if (identical(nm, check_for(ds))) {
        lp <- locked_protein[locked_protein$contrast == nm, , drop = FALSE]
        acc$checks[[paste(nm, "protein")]] <- data.frame(
          check = "refit of a locked contrast against validate_control_spatial_identity.R (protein level)", contrast = nm,
          proteins_compared = sum(lp$ProteinGroupID %in% tt$ProteinGroupID),
          max_abs_logFC_difference = max(abs(tt$logFC[match(lp$ProteinGroupID, tt$ProteinGroupID)] - lp$logFC)),
          terms_only_one_side = NA_integer_, setSize_mismatches = NA_integer_, core_enrichment_mismatches = NA_integer_,
          max_abs_NES_difference = NA_real_, max_abs_p_difference = NA_real_, q05_crossings = NA_integer_,
          stringsAsFactors = FALSE)
      }
    }

    if ("added" %in% fit_stages && length(added_specs)) {
      run_specs <- added_specs
      if (identical(ds, "neuron_neuropil")) run_specs[[check_for(ds)]] <- locked_specs[[check_for(ds)]]
      sides <- lapply(c(L = "L", R = "R"), function(s) {
        keep <- meta$hemisphere == s
        m <- fit_units(mat[, keep, drop = FALSE], meta[keep, , drop = FALSE])
        vapply(run_specs, function(sp) {
          w <- control_spatial_contrast_vector(sp, colnames(m$design$design))
          limma::contrasts.fit(m$fit, matrix(w, ncol = 1))$coefficients[, 1]
        }, numeric(nrow(mat)))
      })
      for (nm in names(run_specs)) {
        message("  added: ", nm)
        tt <- contrast_table(model, run_specs[[nm]], canonical$feature_table, ds)
        go <- run_go(ds, nm, gene_table(tt)$ranked)
        l <- sides$L[, nm]; r <- sides$R[, nm]; ok <- is.finite(l) & is.finite(r)
        if (identical(nm, check_for(ds))) {
          lp <- locked_protein[locked_protein$contrast == nm, , drop = FALSE]
          lg <- locked_go[locked_go$contrast == nm & locked_go$status == "completed", , drop = FALSE]
          m <- merge(go[, c("ID", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment")],
                     lg[, c("ID", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment")], by = "ID", all = TRUE,
                     suffixes = c("", "_locked"))
          acc$checks[[paste(nm, "gsea")]] <- data.frame(
            check = "re-run of a locked contrast against validate_control_spatial_identity.R", contrast = nm,
            proteins_compared = sum(lp$ProteinGroupID %in% tt$ProteinGroupID),
            max_abs_logFC_difference = max(abs(tt$logFC[match(lp$ProteinGroupID, tt$ProteinGroupID)] - lp$logFC)),
            terms_only_one_side = sum(is.na(m$NES) | is.na(m$NES_locked)),
            setSize_mismatches = sum(m$setSize != m$setSize_locked, na.rm = TRUE),
            core_enrichment_mismatches = sum(m$core_enrichment != m$core_enrichment_locked, na.rm = TRUE),
            max_abs_NES_difference = max(abs(m$NES - m$NES_locked), na.rm = TRUE),
            max_abs_p_difference = max(abs(m$pvalue - m$pvalue_locked), na.rm = TRUE),
            q05_crossings = sum((m$p.adjust < 0.05) != (m$p.adjust_locked < 0.05), na.rm = TRUE),
            stringsAsFactors = FALSE)
        }
        if (identical(nm, check_for(ds)) || identical(ds, "microglia")) {
          b <- bil[bil$contrast == nm, , drop = FALSE]
          if (!nrow(b)) stop("No bilateral reference rows for ", nm, call. = FALSE)
          i <- match(b$ProteinGroupID, tt$ProteinGroupID)
          acc$checks[[paste(nm, "bilateral")]] <- data.frame(
            check = "bilateral and one-sided estimates against quantify_bilateral_spatial_identity.R", contrast = nm,
            proteins_compared = sum(!is.na(i)),
            max_abs_logFC_difference = max(abs(c(tt$logFC[i] - b$estimate_bilateral,
                                                 l[i] - b$estimate_L, r[i] - b$estimate_R)), na.rm = TRUE),
            terms_only_one_side = NA_integer_, setSize_mismatches = NA_integer_, core_enrichment_mismatches = NA_integer_,
            max_abs_NES_difference = NA_real_, max_abs_p_difference = NA_real_, q05_crossings = NA_integer_,
            stringsAsFactors = FALSE)
        }
        if (!identical(nm, check_for(ds))) {
          acc$protein[[nm]] <- tt; acc$go[[nm]] <- go
          acc$agreement[[nm]] <- data.frame(dataset = ds, contrast = nm, n_pairs = sum(ok),
                                            pearson_r = stats::cor(l[ok], r[ok]),
                                            spearman_rho = stats::cor(l[ok], r[ok], method = "spearman"),
                                            sign_agreement_fraction = mean(sign(l[ok]) == sign(r[ok])),
                                            stringsAsFactors = FALSE)
          acc$status[[nm]] <- data.frame(dataset = ds, contrast = nm, family = run_specs[[nm]]$family,
                                         manuscript_locked = FALSE, status = "completed",
                                         message = sprintf("duplicateCorrelation=%.4f; hemisphere=%s", model$correlation,
                                                           if (model$design$hemisphere_included) "included" else model$design$hemisphere_omission_reason),
                                         stringsAsFactors = FALSE)
        }
      }
    }

    if ("pairwise" %in% fit_stages) {
      pw <- anatomical_programme_pairwise_specs(ds, units, FAMILIES,
                                                which_families = unique(FAMILIES$family[FAMILIES$shown_in_figure]))
      for (nm in names(pw)) {
        message("  pairwise: ", nm)
        tt <- contrast_table(model, pw[[nm]], canonical$feature_table, ds)
        acc$pairwise[[nm]] <- run_go(ds, nm, gene_table(tt)$ranked)
      }
    }

    if ("hemisphere" %in% fit_stages) {
      shown <- one_vs_rest[names(one_vs_rest) %in% FAMILIES$contrast[FAMILIES$shown_in_figure]]
      for (s in c("L", "R")) {
        keep <- meta$hemisphere == s
        side_model <- fit_units(mat[, keep, drop = FALSE], meta[keep, , drop = FALSE])
        for (nm in names(shown)) {
          message("  ", s, "-only: ", nm)
          tt <- contrast_table(side_model, shown[[nm]], canonical$feature_table, ds)
          go <- run_go(ds, nm, gene_table(tt)$ranked, seed_contrast = paste0(nm, "__", s, "_only"))
          go$side <- s
          acc$side[[paste(nm, s)]] <- go
        }
      }
    }
  }
  check <- do.call(rbind, acc$checks)
  failed <- check$max_abs_logFC_difference > TOL$estimate |
    (!is.na(check$max_abs_NES_difference) & (check$max_abs_NES_difference > TOL$nes |
                                              check$max_abs_p_difference > TOL$p | check$q05_crossings > 0 |
                                              check$setSize_mismatches > 0 | check$core_enrichment_mismatches > 0 |
                                              check$terms_only_one_side > 0))
  check$passed <- !failed
  utils::write.csv(check, OUT$check, row.names = FALSE)
  if (any(failed)) stop("Method reproduction check failed; see ", OUT$check, call. = FALSE)
  utils::write.csv(do.call(rbind, acc$genes), OUT$gene_stats, row.names = FALSE)
  if ("added" %in% fit_stages) {
    utils::write.csv(control_spatial_bind_rows_fill(acc$protein), OUT$protein, row.names = FALSE)
    utils::write.csv(control_spatial_bind_rows_fill(acc$go), OUT$go, row.names = FALSE)
    utils::write.csv(do.call(rbind, acc$agreement), OUT$bilateral, row.names = FALSE)
    utils::write.csv(do.call(rbind, acc$status), OUT$status, row.names = FALSE)
  }
  if ("pairwise" %in% fit_stages) utils::write.csv(control_spatial_bind_rows_fill(acc$pairwise), OUT$pairwise_go, row.names = FALSE)
  if ("hemisphere" %in% fit_stages) utils::write.csv(control_spatial_bind_rows_fill(acc$side), OUT$side_go, row.names = FALSE)
}

# ---------------------------------------------------------------- selection (presentation only)
for (f in c(OUT$go, OUT$pairwise_go, OUT$side_go, OUT$gene_stats)) {
  if (!file.exists(f)) stop("Fit output missing (run the fit stages first): ", f, call. = FALSE)
}
locked_go <- utils::read.csv(LOCKED[["go"]], check.names = FALSE)
added_go <- utils::read.csv(OUT$go, check.names = FALSE)
cols <- c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust", "core_enrichment", "contrast", "status")
all_go <- rbind(locked_go[locked_go$status == "completed", cols], added_go[added_go$status == "completed", cols])
pairwise_go <- utils::read.csv(OUT$pairwise_go, check.names = FALSE)
pairwise_go <- pairwise_go[pairwise_go$status == "completed", c("contrast", "ID", "NES", "p.adjust")]
two_unit <- anatomical_programme_pairs(FAMILIES)
two_unit <- two_unit$contrast[two_unit$reuses_one_vs_rest]
pairwise_go <- rbind(pairwise_go, all_go[all_go$contrast %in% two_unit, c("contrast", "ID", "NES", "p.adjust")])
side_go <- utils::read.csv(OUT$side_go, check.names = FALSE)
side_go <- side_go[side_go$status == "completed", c("contrast", "side", "ID", "NES", "p.adjust")]
gene_stats <- utils::read.csv(OUT$gene_stats, check.names = FALSE)

candidates <- anatomical_programme_candidates(all_go, FAMILIES, RULE)
candidates <- anatomical_programme_term_metrics(candidates, gene_stats, pairwise_go, side_go, FAMILIES)
candidates <- anatomical_programme_strict_select(candidates, STRICT)

picks <- if (file.exists(PICKS)) anatomical_programme_read_picks(PICKS) else NULL
candidates$curated_selected <- FALSE; candidates$curated_display_rank <- NA_integer_; candidates$curated_reason <- NA_character_
if (!is.null(picks)) {
  k <- match(paste(picks$family, picks$unit, picks$ID), paste(candidates$family, candidates$unit, candidates$ID))
  if (anyNA(k)) stop("Curated pick(s) are not supported positive terms of their unit: ",
                     paste(picks$Description[is.na(k)], collapse = "; "), call. = FALSE)
  candidates$curated_selected[k] <- TRUE
  candidates$curated_display_rank[k] <- picks$display_rank
  candidates$curated_reason[k] <- picks$reason
}
off_rule <- candidates$curated_selected & !candidates$strict_eligible
if (any(off_rule)) message("Curated picks outside the strict support rule (see the picks report flags): ",
                           paste(candidates$unit[off_rule], candidates$Description[off_rule], sep = ": ", collapse = "; "))

shown <- FAMILIES[FAMILIES$shown_in_figure, , drop = FALSE]
display <- list(
  automatic = anatomical_programme_display(all_go, candidates, FAMILIES, "selected", "display_rank"),
  strict = anatomical_programme_display(all_go, candidates, FAMILIES, "strict_selected", "strict_display_rank"))
if (!is.null(picks)) display$curated <- anatomical_programme_display(all_go, candidates, FAMILIES, "curated_selected", "curated_display_rank")
display <- do.call(rbind, lapply(names(display), function(s) cbind(selection = s, display[[s]])))
display$display_selection_rule <- c(automatic = RULE$text, strict = STRICT$text,
                                    curated = "curated by the authors from the supported positive terms of each unit; reason per term in config/anatomical_programme_display_picks.csv; criteria met by each pick in anatomical_programme_picks_report.csv")[display$selection]
display$analytical_source <- ifelse(FAMILIES$source[match(display$contrast, FAMILIES$contrast)] == "locked",
                                    relative_to(LOCKED[["go"]]), relative_to(OUT$go))
report_cols <- c("family", "unit", "ID", "Description", "setSize", "NES", "p.adjust", "n_family_units_supported",
                 "candidate", "candidate_rank", "selected", "display_rank", "strict_eligible", "strict_rank",
                 "strict_selected", "strict_display_rank", "curated_selected", "curated_display_rank", "curated_reason",
                 "pairwise_supported_all", "pairwise_min_NES", "pairwise_max_q", "pairwise_detail",
                 "NES_left", "q_left", "NES_right", "q_right", "both_hemispheres_positive",
                 "n_leading_edge", "top_leading_protein", "top_protein_share", "n_leading_edge_supported",
                 "fraction_leading_edge_supported", "artefact_branch", "flag_single_protein",
                 "flag_dominant_protein", "flag_few_proteins", "flag_not_direct", "flag_hemisphere",
                 "redundant_with", "core_enrichment")
utils::write.csv(candidates[, setdiff(report_cols, "core_enrichment")], OUT$candidates, row.names = FALSE)
utils::write.csv(display, OUT$display, row.names = FALSE)
pr <- candidates[candidates$curated_selected | candidates$strict_selected | candidates$selected, report_cols, drop = FALSE]
pr <- pr[pr$family %in% shown$family, , drop = FALSE]
pr$leading_edge_proteins <- pr$core_enrichment; pr$core_enrichment <- NULL
utils::write.csv(pr[order(match(pr$family, unique(FAMILIES$family)), match(paste(pr$family, pr$unit), paste(FAMILIES$family, FAMILIES$unit)),
                          -pr$curated_selected, pr$curated_display_rank, pr$strict_display_rank, pr$display_rank), ],
                 OUT$picks, row.names = FALSE)

# Review renders (Exp9_manuscript draws the figure panels): one block per shown family.
plot_family <- function(d) {
  d$term <- factor(d$Description, levels = rev(unique(d$Description)))
  d$unit <- factor(d$unit, levels = FAMILIES$unit[FAMILIES$family == d$family[1]])
  ggplot(d, aes(unit, term)) +
    geom_point(aes(size = pmin(-log10(p.adjust), 8), fill = pmax(pmin(NES, 3), -3), alpha = supported),
               shape = 21, stroke = 0) +
    scale_fill_gradient2(low = "#3b6fb6", mid = "grey95", high = "#c0392b", limits = c(-3, 3), name = "NES") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
    scale_size_area(max_size = 4, limits = c(0, 8), name = "-log10 q") +
    labs(title = d$family_label[1], x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(panel.grid = element_line(linewidth = 0.2), plot.title = element_text(face = "bold", size = 7.5))
}
render <- function(sel, path) {
  d <- display[display$selection == sel & display$family %in% shown$family, , drop = FALSE]
  if (!nrow(d) || !requireNamespace("patchwork", quietly = TRUE)) return(invisible(NULL))
  ps <- lapply(split(d, factor(d$family, levels = unique(shown$family))), plot_family)
  ps <- ps[vapply(ps, function(p) nrow(p$data) > 0, logical(1))]
  g <- patchwork::wrap_plots(ps, ncol = 2) + patchwork::plot_layout(guides = "collect")
  ggsave(path, g, width = 260, height = 170, units = "mm", device = "svg")
}
render("curated", OUT$plot); render("strict", OUT$plot_strict); render("automatic", OUT$plot_auto)

write_run_manifest(
  MANIFEST,
  inputs = c(as.list(LOCKED), bilateral = BILATERAL, gsea_config = GSEA_CONFIG, picks = if (file.exists(PICKS)) PICKS),
  outputs = OUT,
  parameters = list(
    select_only = select_only, fit_stages = fit_stages, check_contrasts = as.list(CHECK_CONTRAST), tolerance = TOL,
    go_bp_minGSSize = GSEA$min_size, go_bp_maxGSSize = GSEA$max_size,
    gsea_seed_base = GSEA$seed_base, n_perm_simple = GSEA$n_perm,
    gsea_seed_identity = "control_spatial_gsea_reproducibility(): same comparison string as validate_control_spatial_identity.R (dataset + contrast + enrichment type); one-sided fits seed on '<contrast>__<L|R>_only'",
    added_contrasts = FAMILIES$contrast[FAMILIES$source == "added"],
    pairwise_contrasts = anatomical_programme_pairs(FAMILIES)$contrast,
    automatic_rule = RULE$text, strict_rule = STRICT$text,
    artefact_branches = as.list(anatomical_programme_artefact_roots()),
    models = if (length(fit_stages)) models_info else "not refitted (selection only)"
  ),
  notes = "Control-only (three CON animals, both hemispheres). The added contrasts were defined on 2026-10-27, after the 11 prespecified ones; they and the display selections are functional characterisation of the same data, not validation."
)
message("Anatomical programmes complete")
