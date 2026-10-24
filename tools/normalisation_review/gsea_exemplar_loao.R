# Leave-one-animal-out sensitivity of the three GSEA exemplars, under both centrings (normalisation review 2026-10-22).
# Read-only. Each run drops one animal from every unit of the compartment, refits the joint animal-level limma model
# (~0 + unit x group, eBayes robust), puts the refitted statistics into the pipeline's mapped contrast table, and runs
# the pipeline's gene collapse (build_enrichment_gene_inputs: moderated t, median per gene) and seeded gseGO
# (BP, SYMBOL, min 10 / max 800, BH, per-comparison seed; fgsea nPermSimple 100000 for the full data as in the pipeline,
# 10000 for leave-one-out runs, column n_perm_simple; 4 BiocParallel workers). Rows
# already in gsea_exemplar_loao.csv are kept and skipped. Writes results/reviewer_audit/normalisation_choice_2026-10-21/gsea_exemplars/.
# USAGE  Rscript tools/normalisation_review/gsea_exemplar_loao.R      (working directory: repo root)

suppressPackageStartupMessages({ library(readr); library(dplyr); library(limma); library(clusterProfiler); library(org.Mm.eg.db) })
source("R/enrichment/protein_group_enrichment_utils.R"); source("R/enrichment/clusterprofiler_reproducibility.R")
# fgsea parallelises through the default BiocParallel back end (all cores on this host): cap it at 4 workers
BP <- BiocParallel::SnowParam(workers = 4); BiocParallel::bpstart(BP)        # start the workers once, not per call
BiocParallel::register(BP, default = TRUE)
N_PERM <- function(removed) if (removed == "none") 100000L else 10000L     # pipeline setting for the full data; 10^4 for leave-one-out
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
REV <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21")
OUT <- file.path(REV, "gsea_exemplars"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
read_gct <- function(path) {
  lines <- readLines(path, warn = FALSE)
  dims <- as.integer(strsplit(lines[[2]], "\t", fixed = TRUE)[[1]][1:4])
  hdr <- strsplit(lines[[3]], "\t", fixed = TRUE)[[1]]
  desc <- do.call(rbind, strsplit(lines[4:(3 + dims[4])], "\t", fixed = TRUE))
  body <- do.call(rbind, strsplit(lines[(4 + dims[4]):(3 + dims[4] + dims[1])], "\t", fixed = TRUE))
  sc <- (dims[3] + 2L):length(hdr)
  list(m = matrix(as.numeric(body[, sc]), nrow = nrow(body), dimnames = list(body[, 1], hdr[sc])),
       coldesc = setNames(lapply(seq_len(nrow(desc)), function(i) desc[i, sc]), desc[, 1]))
}
gct_path <- function(root, ds) file.path(root, "data", "processed", "01_preprocessing", "protigy_input_animal_level", ds, paste0(ds, "_animal_level.gct"))
GRP <- c(CON = "1", RES = "2", SUS = "3"); ANIMALS <- list(SUS = c("A111", "A755", "A764"), RES = c("A0003", "A135", "A139"), CON = c("A127", "A129", "A765"))
EX <- tibble::tribble(
  ~dataset,          ~unit,           ~file_unit,      ~a,    ~b,    ~go,
  "neuron_neuropil", "CA3_sr",        "CA3sr",         "SUS", "RES", "GO:0099536",
  "neuron_neuropil", "CA3_sr",        "CA3sr",         "SUS", "CON", "GO:0099536",
  "neuron_soma",     "CA2_sp",        "CA2sp",         "SUS", "RES", "GO:0006397",
  "neuron_soma",     "CA2_sp",        "CA2sp",         "SUS", "CON", "GO:0006397",
  "microglia",       "CA1_microglia", "CA1microglia",  "SUS", "RES", "GO:0006119",
  "microglia",       "CA1_microglia", "CA1microglia",  "SUS", "CON", "GO:0006119",
  "microglia",       "CA1_microglia", "CA1microglia",  "RES", "CON", "GO:0006119")
EX <- EX[order(!(EX$a == "SUS" & EX$b == "RES")), ]                  # SUS vs RES first
fit_t <- function(m, cell, contrast) {
  keepc <- colSums(!is.na(m)) > 0; m <- m[, keepc]; cell <- cell[keepc]
  f <- factor(cell, levels = unique(cell)); design <- model.matrix(~0 + f); colnames(design) <- levels(f)
  eb <- eBayes(contrasts.fit(lmFit(m, design), makeContrasts(contrasts = contrast, levels = design)), robust = TRUE)
  tt <- topTable(eb, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
  tibble(protein = rownames(tt), log2fc = tt$logFC, t = tt$t, pval = tt$P.Value, padj = tt$adj.P.Val)
}
G <- list()
prev_path <- file.path(OUT, "gsea_exemplar_loao.csv")
rows <- if (file.exists(prev_path)) list(read_csv(prev_path, show_col_types = FALSE, col_types = cols(removed = "c"))) else list()
if (length(rows) && !"n_perm_simple" %in% names(rows[[1]])) rows[[1]]$n_perm_simple <- 100000L   # rows written before the switch
done_key <- if (length(rows)) with(rows[[1]], paste(dataset, unit, comparison, centring, removed)) else character()
for (i in seq_len(nrow(EX))) {
  e <- EX[i, ]; ds <- e$dataset
  if (is.null(G[[ds]])) G[[ds]] <- list(new = read_gct(gct_path(".", ds)), original = read_gct(gct_path(ARCHIVE, ds)))
  comp <- paste0(e$file_unit, tolower(e$a), "_", e$file_unit, tolower(e$b))
  mapped <- read_csv(file.path("data/processed/02_id_mapping_animal_level/mapped", ds, "forward/per_file", paste0(comp, ".csv")), show_col_types = FALSE, guess_max = 1e5)
  seed <- derive_clusterprofiler_gsea_seed(20260824L, comp, clusterprofiler_gsea_analysis_types("BP")[["go"]])
  for (v in c("new", "original")) {
    g <- G[[ds]][[v]]; cell <- g$coldesc$region_layer_ExpGroup; an <- g$coldesc$AnimalID
    u <- sub("_[0-9]+$", "", cell[grepl(paste0("^", sub("_microglia$", "", e$unit), "(_microglia)?_[0-9]+$"), cell)][1])
    contrast <- paste0(u, "_", GRP[[e$a]], "-", u, "_", GRP[[e$b]])
    for (rm_a in c("none", ANIMALS[[e$a]], ANIMALS[[e$b]])) {
      if (paste(ds, e$unit, paste(e$a, "vs", e$b), v, rm_a) %in% done_key) next
      m <- g$m; if (rm_a != "none") m[, an == rm_a] <- NA
      st <- fit_t(m, cell, contrast)
      df <- mapped |> dplyr::select(-any_of(c("log2fc", "aveExpr", "t", "pval", "padj", "B", "significant", "logpval", "sign.logP"))) |>
        inner_join(st, by = c("original_identifier" = "protein"))
      gi <- build_enrichment_gene_inputs(as.data.frame(df))
      gse <- run_seeded_clusterprofiler_gsea(clusterProfiler::gseGO, gsea_seed = seed, n_perm_simple = N_PERM(rm_a),
        geneList = gi$ranked, ont = "BP", keyType = "SYMBOL", minGSSize = 10, maxGSSize = 800, pvalueCutoff = 1, verbose = FALSE,
        OrgDb = org.Mm.eg.db, pAdjustMethod = "BH")
      r <- as.data.frame(gse) |> arrange(pvalue); k <- match(e$go, r$ID)
      rows[[length(rows) + 1L]] <- tibble(dataset = ds, unit = e$unit, comparison = paste(e$a, "vs", e$b), go = e$go, centring = v, removed = rm_a,
        n_genes = length(gi$ranked), NES = r$NES[k], pvalue = r$pvalue[k], p_adjust = r$p.adjust[k], rank_by_p = k, n_terms = nrow(r),
        n_terms_fdr05 = sum(r$p.adjust < 0.05), protein_FDR05 = sum(st$padj < 0.05), n_perm_simple = N_PERM(rm_a))
      message(ds, " ", comp, " ", v, " -", rm_a, ": NES ", round(r$NES[k], 2), " FDR ", signif(r$p.adjust[k], 2))
      write_csv(bind_rows(rows), prev_path, na = "")
    }
  }
}
res <- bind_rows(rows)
summ <- res |> group_by(dataset, unit, comparison, go, centring) |>
  summarise(NES_all = NES[removed == "none"], FDR_all = p_adjust[removed == "none"], min_abs_NES = min(abs(NES)), max_FDR = max(p_adjust),
            worst_removed = removed[which.max(p_adjust)], n_loao_fdr05 = sum(p_adjust[removed != "none"] < 0.05), n_loao = sum(removed != "none"),
            sign_stable = all(sign(NES) == sign(NES[removed == "none"])), .groups = "drop")
write_csv(summ, file.path(OUT, "gsea_exemplar_loao_summary.csv"), na = "")
options(width = 220); print(as.data.frame(summ))
cat("done:", OUT, "\n")
