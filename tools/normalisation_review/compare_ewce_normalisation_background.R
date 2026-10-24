# EWCE: separate the effect of the background from the effect of the centring (normalisation review 2026-10-21).
#
#   A  archived run       original centring, reference-transcriptome background (the fallback the old run used)
#   C  comparison branch  original centring, measured-proteome background (PROTEOMICS_EWCE_BRANCH run)
#   B  current run        complete-case centring, measured-proteome background
#
# A vs C isolates the step from run A to run C (background, species space and code since b478a13 together), C vs B the
# centring, A vs B is the total change. Baseline targets (top-N proteins by stratum abundance) are invariant to a
# per-sample constant by construction; their gene lists are not compared here, only their scores and calls (a few calls
# can still differ through the global FDR). Read-only; writes validation/ewce_normalisation_vs_background.csv.
# Archived run: <archive>/results/tables/05_celltype_enrichment_EWCE/EWCE_E9/<dataset>/Supplementary_Table_EWCE.xlsx (2026-08-21).
# USAGE  Rscript tools/normalisation_review/compare_ewce_normalisation_background.R

suppressPackageStartupMessages({ library(dplyr); library(readr); library(openxlsx) })
ARCHIVE <- if (dir.exists("Q:/results")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
OUT <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "validation")
BRANCH <- "original_centring_proteome_background"
DATASETS <- c("neuron_neuropil", "neuron_soma", "microglia")
load_run <- function(path, run) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  read.xlsx(path, sheet = "Full_Results") |>
    transmute(run = run, Target, CellType, annotLevel, AnalysisType, TopN, N_Background, sd_from_mean, q_global, Significant_Global)
}
rows <- list()
for (ds in DATASETS) {
  pB <- file.path("results", "enrichment", "run_ewce_celltype_enrichment", ds, "tables", "Supplementary_Table_EWCE.xlsx")
  pC <- file.path("results", "enrichment", "run_ewce_celltype_enrichment", "comparison", BRANCH, ds, "tables", "Supplementary_Table_EWCE.xlsx")
  pA <- file.path(ARCHIVE, "results", "tables", "05_celltype_enrichment_EWCE", "EWCE_E9", ds, "Supplementary_Table_EWCE.xlsx")
  runs <- list(A = load_run(pA, "A"), B = load_run(pB, "B"), C = load_run(pC, "C"))
  cat(ds, ": A", pA, "| B", file.exists(pB), "| C", file.exists(pC), "\n")
  pairs <- list(c("A", "C", "background (same original centring)"), c("C", "B", "centring (same proteome background)"),
                c("A", "B", "total (as archived vs as now)"))
  for (pr in pairs) {
    x <- runs[[pr[1]]]; y <- runs[[pr[2]]]
    if (is.null(x) || is.null(y)) next
    j <- inner_join(x, y, by = c("Target", "CellType", "annotLevel", "AnalysisType", "TopN"), suffix = c("_1", "_2"))
    for (at in sort(unique(j$AnalysisType))) {
      k <- j[j$AnalysisType == at & j$annotLevel == 2 & j$TopN == 250, ]
      if (!nrow(k)) next
      rows[[length(rows) + 1L]] <- tibble(dataset = ds, analysis_type = at, effect = pr[3], runs = paste(pr[1], "vs", pr[2]),
        n_tests = nrow(k), background_size_1 = median(k$N_Background_1), background_size_2 = median(k$N_Background_2),
        cor_sd_from_mean = suppressWarnings(stats::cor(k$sd_from_mean_1, k$sd_from_mean_2)),
        significant_1 = sum(k$Significant_Global_1 %in% TRUE), significant_2 = sum(k$Significant_Global_2 %in% TRUE),
        significant_in_both = sum(k$Significant_Global_1 %in% TRUE & k$Significant_Global_2 %in% TRUE),
        call_agreement = mean((k$Significant_Global_1 %in% TRUE) == (k$Significant_Global_2 %in% TRUE)))
    }
  }
}
res <- bind_rows(rows)
write_csv(res, file.path(OUT, "ewce_normalisation_vs_background.csv"), na = "")
print(as.data.frame(res |> mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
