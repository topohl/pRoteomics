#!/usr/bin/env Rscript
#
# Do SUS - RES FDR-supported proteins occupy unusual positions within their
# frozen WGCNA modules?
#
# THE QUESTION
#   Given the proteins that were actually TESTED for protein-level SUS - RES
#   differential abundance, are the FDR-supported ones more central or more
#   peripheral within their own modules than an equal number of tested proteins
#   drawn from the same modules?
#
# WHY THE A1/A2/D COUNTS CANNOT ANSWER IT
#   Those tiers are defined BY kME. A module of ~335 members has exactly 10
#   top-10 hubs, so ~3% of any random draw lands there. A small Tier A1 count is
#   therefore the null expectation, not evidence of peripheral targeting. This
#   script uses the CONTINUOUS within-module rank instead and compares against a
#   module-stratified permutation null.
#
# WHAT THIS SCRIPT DOES NOT DO
#   * No WGCNA is recomputed: kME and the frozen hub flags are copied verbatim.
#   * No differential abundance is recomputed: log2FC / p / BH FDR are copied
#     from the manifest-selected Stage 02 mapped contrast files.
#   * Nothing here is causal, and a permutation p-value is not independent
#     replication. Spatial contexts are repeated measurements of the same
#     animals; "FDR-supported in >= 1 context" is treated as a single label.
#
# SCOPE SEPARATION
#   This script is PHENOTYPE-AWARE by design. The companion label-coherence
#   audit (analysis/05_wgcna/14_wgcna_label_coherence_audit.R) is strictly
#   phenotype-blind and must never consume this script's outputs.
#
# USAGE
#   Rscript analysis/08_integration/09_wgcna_sus_res_network_position.R
#   Rscript analysis/08_integration/09_wgcna_sus_res_network_position.R --dataset microglia
#   Rscript analysis/08_integration/09_wgcna_sus_res_network_position.R --permutations 20000
#   Rscript analysis/08_integration/09_wgcna_sus_res_network_position.R --dry-run

source("R/paths.R")
source("R/data_contracts/dataset_config.R")
source("R/statistics/integration_utils.R")
source("R/enrichment/enrichment_io.R")
source("R/statistics/sus_res_spatial_dap_atlas_utils.R")
source("R/statistics/wgcna_candidate_protein_utils.R")
source("R/networks/wgcna_network_position_utils.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

SCRIPT_ID <- "analysis/08_integration/09_wgcna_sus_res_network_position.R"
SUBSTEP <- "wgcna_sus_res_network_position"
Sys.setenv(PROTEOMICS_SCRIPT_ID = SCRIPT_ID)

cli <- integration_cli(default_dataset = "all")
cli_value <- function(flag, default = "") {
  hit <- which(cli$args == flag)
  if (!length(hit) || hit[[1]] == length(cli$args)) return(default)
  cli$args[[hit[[1]] + 1L]]
}
n_permutations <- suppressWarnings(as.integer(
  cli_value("--permutations", as.character(wnp_default_permutations()))
))
if (is.na(n_permutations) || n_permutations < 1L) {
  stop("--permutations must be a positive integer.", call. = FALSE)
}

membership_path <- function(dataset) {
  path_results("tables", "06_modules_WGCNA", "01_WGCNA", dataset, "modules",
               "WGCNA_modules_long.csv")
}
candidate_path <- function(dataset) {
  path_results("tables", "10_biological_integration",
               "wgcna_candidate_protein_shortlist", dataset,
               "wgcna_candidate_proteins_all.csv")
}

datasets <- integration_datasets(cli$dataset)

if (isTRUE(cli$dry_run)) {
  inputs <- list()
  for (dataset in datasets) {
    inputs[[paste0("module_membership__", dataset)]] <- membership_path(dataset)
    inputs[[paste0("clusterprofiler_manifest__", dataset)]] <-
      canonical_clusterprofiler_manifest_path(dataset)
    inputs[[paste0("candidate_summary__", dataset)]] <- candidate_path(dataset)
  }
  cat("[DRY-RUN] WGCNA SUS-RES network-position audit; no outputs written.\n")
  dry_run_inputs(SCRIPT_ID, inputs)
  cat("[DRY-RUN] Would write to results/{tables,figures,source_data,reports,logs}/",
      "10_biological_integration/", SUBSTEP, "/<dataset>/ and /global/.\n", sep = "")
  quit(save = "no", status = 0L)
}

read_required <- function(path, label) {
  if (!file.exists(path)) {
    stop("missing_required_input: ", label, ": ", path, call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE, guess_max = Inf)
}
relpath <- function(path) relative_to(normalizePath(path, winslash = "/", mustWork = FALSE))

# ------------------------------------------------------------ DA eligibility

# The null universe is the set of proteins actually TESTED for SUS - RES DA, not
# every module member. "Tested" reuses the canonical Stage 04 definition from
# R/statistics/sus_res_spatial_dap_atlas_utils.R: a finite BH FDR and a finite effect.
load_sus_res_eligibility <- function(dataset) {
  manifest_path <- canonical_clusterprofiler_manifest_path(dataset)
  if (!file.exists(manifest_path)) {
    stop("missing_required_input: clusterProfiler manifest: ", manifest_path,
         call. = FALSE)
  }
  manifest <- readr::read_csv(
    manifest_path, show_col_types = FALSE, progress = FALSE, guess_max = Inf,
    col_types = readr::cols(.default = readr::col_character())
  ) %>%
    filter(.data$result_type == "GSEA_GO",
           .data$route_category == "phenotype_within_unit") %>%
    mutate(contrast = wcp_contrast_from_comparison(.data$comparison)) %>%
    filter(.data$contrast == wcp_primary_contrast()) %>%
    distinct(.data$comparison, .keep_all = TRUE)
  if (!nrow(manifest)) {
    stop("missing_required_input: no SUS - RES manifest rows for ", dataset,
         call. = FALSE)
  }

  per_context <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
    row <- manifest[i, , drop = FALSE]
    input <- sus_res_resolve_manifest_input(row$input_gene_file[[1]], dataset,
                                            repo_root())
    z <- read_required(input, "manifest-selected protein DA input")
    orientation <- wcp_comparison_orientation(row$comparison[[1]])
    effect <- suppressWarnings(as.numeric(z$log2fc)) *
      orientation$formal_effect_multiplier[[1]]
    fdr <- suppressWarnings(as.numeric(z$padj))
    tested <- is.finite(fdr) & is.finite(effect)
    tibble(
      ProteinGroupID = as.character(z$ProteinGroupID),
      spatial_unit = row$route_unit[[1]],
      is_tested = tested,
      is_hit = tested & fdr <= wcp_fdr_threshold() & effect != 0,
      abs_effect = abs(effect),
      fdr = fdr,
      ave_expr = if ("aveExpr" %in% names(z)) suppressWarnings(as.numeric(z$aveExpr)) else NA_real_
    )
  }))

  per_context %>%
    group_by(.data$ProteinGroupID) %>%
    summarise(
      n_sus_res_contexts_tested = sum(.data$is_tested),
      is_da_tested = any(.data$is_tested),
      is_hit = any(.data$is_hit),
      n_sus_res_contexts_fdr05 = sum(.data$is_hit),
      sus_res_min_BH_FDR = suppressWarnings(min(.data$fdr[.data$is_tested], na.rm = TRUE)),
      sus_res_max_abs_log2FC = suppressWarnings(max(.data$abs_effect[.data$is_tested], na.rm = TRUE)),
      # Mean log2 abundance across tested contexts. aveExpr is the limma
      # per-contrast average expression carried in the frozen mapped DA files;
      # it is used ONLY as a detectability covariate for the sensitivity
      # analysis, never as evidence.
      mean_ave_expr = suppressWarnings(mean(.data$ave_expr[.data$is_tested], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(across(
      c("sus_res_min_BH_FDR", "sus_res_max_abs_log2FC", "mean_ave_expr"),
      ~ ifelse(is.finite(.x), .x, NA_real_)
    ))
}

# ------------------------------------------------------------- build dataset

build_dataset <- function(dataset) {
  membership <- read_required(membership_path(dataset),
                              "frozen WGCNA module membership")
  keep <- c("dataset", "ModuleID", "ModuleLabel_Final", "ProteinGroupID",
            "GeneSymbol", "kME", "abs_kME", "is_core_kME_0.6", "is_top_hub_25")
  membership <- membership[, intersect(keep, names(membership)), drop = FALSE]
  membership <- wcp_rank_module_members(membership)

  eligibility <- load_sus_res_eligibility(dataset)

  # Tier membership is joined for DESCRIPTIVE module composition only. It is
  # never used to define hits, and never enters the permutation.
  tiers <- NULL
  cand_path <- candidate_path(dataset)
  if (file.exists(cand_path)) {
    tiers <- read_required(cand_path, "candidate summary") %>%
      select(any_of(c("ProteinGroupID", "is_tier_A1", "is_tier_A2", "is_tier_D",
                      "candidate_tier", "phenotype_network_class")))
  }

  out <- membership %>%
    inner_join(eligibility, by = "ProteinGroupID") %>%
    filter(.data$is_da_tested)
  if (!is.null(tiers)) out <- left_join(out, tiers, by = "ProteinGroupID")

  # Ranks come from the FULL frozen module, so rank_fraction stays comparable to
  # the frozen membership even though the eligible universe is a subset.
  out$rank_fraction <- wnp_rank_fraction(out$abs_kME_rank_in_module,
                                         out$n_module_members)
  out$is_top10_module_hub <- out$is_top10_module_hub %in% TRUE
  out$membership_source_file <- relpath(membership_path(dataset))
  out$contract_version <- wnp_contract_version()
  out$Source <- SCRIPT_ID
  out
}

# -------------------------------------------------- detectability sensitivity

# Sensitivity permutation stratified by module x abundance tertile.
#
# aveExpr is the only canonical per-protein detectability quantity carried in
# the frozen DA files. If it is unavailable the analysis is skipped and recorded
# as a limitation rather than backed by an invented proxy.
detectability_strata <- function(data, n_strata = 3L) {
  if (!("mean_ave_expr" %in% names(data))) return(NULL)
  v <- suppressWarnings(as.numeric(data$mean_ave_expr))
  if (!any(is.finite(v))) return(NULL)
  out <- rep(NA_character_, length(v))
  for (mod in unique(as.character(data$ModuleID))) {
    sel <- which(as.character(data$ModuleID) == mod & is.finite(v))
    if (length(sel) < n_strata) {
      out[sel] <- paste0(mod, "|all")
      next
    }
    br <- stats::quantile(v[sel], probs = seq(0, 1, length.out = n_strata + 1L),
                          na.rm = TRUE)
    br[1] <- -Inf; br[length(br)] <- Inf
    grp <- cut(v[sel], breaks = unique(br), include.lowest = TRUE, labels = FALSE)
    out[sel] <- paste0(mod, "|", ifelse(is.na(grp), 0L, grp))
  }
  out[is.na(out)] <- paste0(as.character(data$ModuleID)[is.na(out)], "|na")
  out
}

# ----------------------------------------------------------------- figures

save_figure <- function(plot, dir, stem, width_mm, height_mm) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NA_character_)
  dir_create(dir)
  path <- file.path(dir, paste0(stem, ".png"))
  ggplot2::ggsave(path, plot, width = width_mm, height = height_mm,
                  units = "mm", dpi = 300, bg = "white")
  pdf_path <- file.path(dir, paste0(stem, ".pdf"))
  ok <- tryCatch({
    ggplot2::ggsave(pdf_path, plot, width = width_mm, height = height_mm,
                    units = "mm", device = grDevices::cairo_pdf, bg = "white")
    TRUE
  }, error = function(e) FALSE)
  if (!ok && file.exists(pdf_path)) unlink(pdf_path)
  path
}

build_figures <- function(dataset, protein, modules, perm, fig_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(character())
  library(ggplot2)
  base <- theme_classic(base_size = 8) +
    theme(strip.background = element_blank(),
          plot.title = element_text(face = "bold", size = 9))
  written <- character()

  hits <- protein[protein$is_hit %in% TRUE, , drop = FALSE]

  # 1. observed hit rank fractions against the module-stratified null
  null_med <- perm$null[, "median_rank_fraction"]
  obs_med <- perm$observed[["median_rank_fraction"]]
  if (any(is.finite(null_med))) {
    p1 <- ggplot(data.frame(x = null_med[is.finite(null_med)]), aes(x = .data$x)) +
      geom_histogram(bins = 60, fill = "#C6C3BB", colour = NA) +
      geom_vline(xintercept = obs_med, colour = "#E63947", linewidth = 0.7) +
      geom_vline(xintercept = 0.5, colour = "#3D3B6E", linetype = "dashed",
                 linewidth = 0.4) +
      labs(
        title = paste0(dataset, ": median within-module rank fraction of SUS-RES hits"),
        subtitle = sprintf(
          "red = observed (%.3f); grey = %s module-stratified draws; dashed = 0.5. Lower is more central.",
          obs_med, format(perm$n_permutations, big.mark = ",")),
        x = "Median within-module |kME| rank fraction", y = "Null draws"
      ) + base + theme(plot.subtitle = element_text(size = 6.5))
    written <- c(written, save_figure(p1, fig_dir,
      paste0("network_position_null_median_rank_fraction_", dataset), 180, 95))
  }

  # 2. rank-fraction distribution, hits vs eligible non-hits
  d2 <- data.frame(
    rank_fraction = suppressWarnings(as.numeric(protein$rank_fraction)),
    group = ifelse(protein$is_hit %in% TRUE,
                   "SUS-RES FDR-supported", "DA-tested, not FDR-supported"),
    stringsAsFactors = FALSE
  )
  d2 <- d2[is.finite(d2$rank_fraction), , drop = FALSE]
  if (nrow(d2)) {
    p2 <- ggplot(d2, aes(x = .data$rank_fraction, colour = .data$group)) +
      stat_ecdf(linewidth = 0.6) +
      scale_colour_manual(values = c("DA-tested, not FDR-supported" = "#3D3B6E",
                                     "SUS-RES FDR-supported" = "#E63947")) +
      labs(title = paste0(dataset, ": within-module position, hits vs eligible non-hits"),
           subtitle = "Empirical CDF of the within-module |kME| rank fraction. Diagonal = uniform.",
           x = "Within-module |kME| rank fraction (0 = most central)",
           y = "Cumulative fraction", colour = NULL) +
      base + theme(legend.position = "bottom",
                   plot.subtitle = element_text(size = 6.5))
    written <- c(written, save_figure(p2, fig_dir,
      paste0("network_position_rank_fraction_ecdf_", dataset), 150, 110))
  }

  # 3. module-by-module hit topology
  if (nrow(hits)) {
    d3 <- protein
    d3$module <- factor(as.character(d3$ModuleID),
                        levels = sort(unique(as.character(d3$ModuleID))))
    p3 <- ggplot(d3, aes(x = .data$module, y = .data$rank_fraction)) +
      geom_jitter(data = d3[!(d3$is_hit %in% TRUE), , drop = FALSE],
                  width = 0.25, height = 0, colour = "#DDDBD5", size = 0.25) +
      geom_jitter(data = d3[d3$is_hit %in% TRUE, , drop = FALSE],
                  width = 0.15, height = 0, colour = "#E63947", size = 1.1) +
      geom_hline(yintercept = 0.5, linetype = "dashed", colour = "#3D3B6E",
                 linewidth = 0.3) +
      scale_y_reverse() +
      labs(title = paste0(dataset, ": module-by-module hit topology"),
           subtitle = "Red = SUS-RES FDR-supported; grey = DA-tested non-hits. Axis reversed: central at top.",
           x = NULL, y = "Within-module |kME| rank fraction") +
      base + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
                   plot.subtitle = element_text(size = 6.5))
    written <- c(written, save_figure(p3, fig_dir,
      paste0("network_position_module_topology_", dataset), 190, 110))
  }
  written[!is.na(written)]
}

# --------------------------------------------------------------------- main

message("[", SCRIPT_ID, "] permutations = ", n_permutations)
built <- list()
for (dataset in datasets) {
  message("[", SCRIPT_ID, "] building ", dataset)
  built[[dataset]] <- build_dataset(dataset)
}

emit <- function(scope, protein) {
  paths <- create_module_dirs("10_biological_integration",
                              file.path(SUBSTEP, scope))
  perm <- wnp_permute(protein, n_permutations = n_permutations)
  perm_summary <- wnp_permutation_summary(
    perm, dataset = scope, scope = "module_stratified"
  )
  perm_summary$interpretation <- wnp_interpret(perm_summary)

  # detectability sensitivity: module x abundance-tertile strata
  strata <- detectability_strata(protein)
  if (!is.null(strata)) {
    sens_input <- protein
    sens_input$ModuleID <- strata      # stratify on the finer grouping
    sens <- wnp_permute(sens_input, n_permutations = n_permutations)
    sens_summary <- wnp_permutation_summary(
      sens, dataset = scope, scope = "module_x_abundance_tertile_stratified"
    )
    sens_summary$interpretation <- wnp_interpret(sens_summary)
    perm_summary <- rbind(perm_summary, sens_summary)
  }

  modules <- wnp_module_summary(protein)

  written <- list(
    protein = write_integration_table(
      protein, paths, "wgcna_sus_res_network_position_protein_level.csv"),
    modules = write_integration_table(
      modules, paths, "wgcna_sus_res_network_position_module_summary.csv"),
    permutation = write_integration_table(
      perm_summary, paths, "wgcna_sus_res_network_position_permutation_summary.csv")
  )

  figures <- if (!identical(scope, "global")) {
    build_figures(scope, protein, modules, perm, paths$figures)
  } else character()

  readme <- write_readme(paths, scope, protein, modules, perm_summary, strata)

  write_integration_manifest(
    paths,
    inputs = as.list(stats::setNames(
      unique(stats::na.omit(as.character(protein$membership_source_file))),
      "module_membership"
    )),
    outputs = list(
      tables = unlist(lapply(written, `[[`, "table"), use.names = FALSE),
      source_data = unlist(lapply(written, `[[`, "source"), use.names = FALSE),
      figures = figures, report = readme
    ),
    parameters = list(
      contract_version = wnp_contract_version(),
      scope = scope,
      n_permutations = n_permutations,
      permutation_seed = perm$seed,
      null_universe = "proteins with a finite SUS - RES BH FDR and finite effect in >= 1 spatial context",
      stratification = "module; sensitivity adds abundance tertiles within module",
      hit_definition = "SUS - RES BH FDR <= 0.05 with a non-zero effect in >= 1 spatial context",
      rank_fraction = "(within-module |kME| rank - 0.5) / module size"
    ),
    notes = paste(
      "Descriptive positional comparison of SUS - RES FDR-supported proteins",
      "against a module-stratified permutation null. No WGCNA and no",
      "differential abundance is recomputed. Permutation p-values are not",
      "independent replication; spatial contexts are repeated measurements of",
      "the same animals. The A1/A2/D tiers are joined for descriptive module",
      "composition only and never enter the null or the test."
    )
  )
  list(paths = paths, protein = protein, modules = modules,
       permutation = perm_summary, perm = perm)
}

write_readme <- function(paths, scope, protein, modules, perm_summary, strata) {
  primary <- perm_summary[perm_summary$scope == "module_stratified", , drop = FALSE]
  lines <- c(
    paste0("# SUS - RES network position within frozen WGCNA modules (", scope, ")"),
    "",
    paste0("Generated by `", SCRIPT_ID, "` (contract ", wnp_contract_version(), ")."),
    "",
    "## Question",
    "",
    "Given the proteins actually TESTED for protein-level SUS - RES differential",
    "abundance, are the FDR-supported ones positioned unusually centrally or",
    "peripherally within their own WGCNA modules?",
    "",
    "## Why the A1/A2/D counts cannot answer this",
    "",
    "Those tiers are defined by |kME| itself. A module of ~335 members has exactly",
    "10 top-10 hubs, so roughly 3% of ANY random draw falls there. A small Tier A1",
    "count is the null expectation, not evidence of peripheral targeting. This",
    "analysis therefore uses the continuous within-module rank.",
    "",
    "## Design",
    "",
    paste0("- Null universe: DA-TESTED proteins only (finite BH FDR and finite effect",
           " in >= 1 SUS - RES spatial context), not all module members."),
    "- Stratification: within every module, draw exactly as many proteins as that",
    "  module actually contributed hits. This holds module size, module |kME|",
    "  distribution and per-module hit count fixed.",
    paste0("- Permutations: ", format(unique(primary$n_permutations), big.mark = ","),
           "; p-values use the (1 + count) / (B + 1) correction."),
    "- Rank fraction = (within-module |kME| rank - 0.5) / module size; 0 = most",
    "  central. The midpoint form makes modules of different sizes comparable.",
    "",
    "## Result",
    "",
    paste0("Eligible (DA-tested) module members: ", nrow(protein),
           "; SUS - RES FDR-supported: ", sum(protein$is_hit %in% TRUE),
           " across ", sum(modules$n_hits > 0), " modules."),
    "",
    "| Statistic | Observed | Null median | Null 95% | p central | p peripheral | Reading |",
    "| --- | --- | --- | --- | --- | --- | --- |"
  )
  for (i in seq_len(nrow(primary))) {
    lines <- c(lines, sprintf(
      "| %s | %.4f | %.4f | %.4f - %.4f | %.4f | %.4f | %s |",
      primary$statistic_label[i], primary$observed[i], primary$null_median[i],
      primary$null_q025[i], primary$null_q975[i], primary$p_more_central[i],
      primary$p_more_peripheral[i], primary$interpretation[i]
    ))
  }
  lines <- c(
    lines, "",
    "## Caveats",
    "",
    "- This is a DESCRIPTIVE positional comparison. It is not causal, and a",
    "  permutation p-value is not independent replication.",
    "- Spatial contexts are repeated measurements of the SAME animals.",
    "  \"FDR-supported in at least one context\" is treated as a single label; no",
    "  independence across contexts is assumed or claimed.",
    "- The module-level hypergeometric test in the module summary is EXPLORATORY.",
    "  It uses the DA-tested universe and is BH corrected across modules within a",
    "  dataset, but it does not account for differing power across spatial",
    "  contexts or for correlation between them.",
    if (is.null(strata)) {
      "- LIMITATION: no canonical per-protein detectability covariate was available, so the abundance-stratified sensitivity analysis was skipped rather than backed by an invented proxy."
    } else {
      "- Sensitivity: the permutation was repeated stratifying on module x abundance tertile, using `aveExpr` from the frozen mapped DA files as the only canonical detectability covariate. Both scopes are reported in the permutation summary."
    },
    "- Datasets with very few FDR-supported proteins (neuron_soma and microglia",
    "  currently have 3 each) cannot support dataset-specific inference. Their",
    "  permutation results are reported for completeness only.",
    "- A1/A2/D tier columns appear in the module summary for descriptive",
    "  composition. They are never used to define hits or to build the null.",
    ""
  )
  path <- file.path(paths$reports, "README.md")
  dir_create(dirname(path))
  writeLines(lines, path)
  path
}

emitted <- list()
for (dataset in datasets) {
  emitted[[dataset]] <- emit(dataset, built[[dataset]])
}

# Secondary combined analysis: strata are dataset x module, so no module is ever
# compared against a module from another dataset.
combined <- NULL
if (setequal(datasets, valid_datasets())) {
  pooled <- bind_rows(built)
  pooled$ModuleID <- paste(pooled$dataset, pooled$ModuleID, sep = "::")
  combined <- emit("global", pooled)
}

# ------------------------------------------------------------ console report

cat("\n===== SUS-RES network position within WGCNA modules =====\n")
cat("Contract      : ", wnp_contract_version(), "\n", sep = "")
cat("Permutations  : ", format(n_permutations, big.mark = ","),
    " (module-stratified)\n", sep = "")
for (dataset in datasets) {
  e <- emitted[[dataset]]
  p <- e$permutation[e$permutation$scope == "module_stratified", , drop = FALSE]
  cat("\n--- ", dataset, " ---\n", sep = "")
  cat("  DA-tested module members : ", format(nrow(e$protein), big.mark = ","), "\n", sep = "")
  cat("  SUS-RES FDR-supported    : ", sum(e$protein$is_hit %in% TRUE),
      " in ", sum(e$modules$n_hits > 0), " modules\n", sep = "")
  if (sum(e$protein$is_hit %in% TRUE) < 5L) {
    cat("  NOTE: too few hits for dataset-specific inference; reported for completeness.\n")
  }
  for (i in seq_len(nrow(p))) {
    cat(sprintf("    %-42s obs=%7.4f null=%7.4f [%6.4f,%6.4f]  p_cen=%.4f p_per=%.4f  %s\n",
                p$statistic_label[i], p$observed[i], p$null_median[i],
                p$null_q025[i], p$null_q975[i], p$p_more_central[i],
                p$p_more_peripheral[i], p$interpretation[i]))
  }
}
cat("\nOutputs written under:\n")
for (dataset in datasets) cat("  ", relpath(emitted[[dataset]]$paths$tables), "\n", sep = "")
if (!is.null(combined)) cat("  ", relpath(combined$paths$tables), "\n", sep = "")
cat("\nPermutation p-values are descriptive, not independent replication.\n")
cat("Spatial contexts are repeated measurements of the same animals.\n")
