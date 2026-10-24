# Follow-up checks for the centring validation (normalisation review 2026-10-21). Read-only for the pipeline.
#
#   A  CA2-SLM: is the SUS-vs-RES result carried by A755 L and A764 R, or does removing any two hemispheres lose it?
#      Every equivalent removal (one hemisphere of each of two animals of one group) is refitted with the pipeline
#      limma model, for the original, new and median-ratio centrings.
#   B  Functional composition of the anchor proteins (quantified in every sample) against all proteins kept by the
#      70% filter, per compartment; and of the non-anchor proteins (what the anchor set under-represents).
#   C  WGCNA stability: the network rebuilt with the pipeline settings under the original, new and median-ratio
#      centrings on the same samples and proteins, and leave-one-animal-out rebuilds of each. Every build is cached,
#      so the step resumes after an interruption. The pipeline's sample-tree outlier cut is not re-applied (the
#      current network's samples are used for every estimator), so the neuropil original build is not the previous
#      pipeline network (that one dropped A135 L CA1 so).
#   C2 the same at a common soft power (the current network's), because the power rule (lowest power with scale-free
#      R2 >= 0.80) sits on a knife edge and a one-step power change alone moves partitions by ARI 0.6-0.8.
#
# Needs validate_sample_centring.R to have run (reads its validation_state.rds).
# USAGE  Rscript tools/normalisation_review/validate_sample_centring_followup.R [--part A|B|C|all] [--threads 8]

suppressPackageStartupMessages({ library(readr); library(readxl); library(dplyr); library(tidyr); library(limma); library(matrixStats) })
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) { i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[[i + 1L]] }
PART <- arg("--part", "all")
if (identical(PART, "C2")) PART <- c("C", "C2")
ARCHIVE <- if (dir.exists("Q:/data")) "Q:" else
  "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_pre_renormalisation_2026-10-19"
OUT <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "validation")
w <- function(x, name) write_csv(x, file.path(OUT, name), na = "")
state <- readRDS(file.path(OUT, "validation_state.rds"))
OFF <- state$offsets
EST3 <- c("original_all_observed", "new_complete_case", "median_ratio_complete_case")
DATASETS <- c("microglia", "neuron_neuropil", "neuron_soma")

read_imputed <- function(path) {
  x <- read_excel(path, progress = FALSE)
  m <- as.matrix(x[, grepl("[.]d$", names(x))]); rownames(m) <- x$Protein.Group
  list(m = m, names = x$`T: Protein.Names`)
}
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
fit_counts <- function(m, cell, unit = "CA2_slm") {
  f <- factor(cell, levels = unique(cell)); design <- model.matrix(~0 + f); colnames(design) <- levels(f)
  units <- unique(sub("_[0-9]+$", "", cell)); pairs <- list(c("2", "1"), c("3", "2"), c("3", "1"))
  nm <- unlist(lapply(pairs, function(p) paste0(units, "_", p[1], "_over_", units, "_", p[2])))
  ex <- unlist(lapply(pairs, function(p) paste0(units, "_", p[1], "-", units, "_", p[2])))
  cm <- makeContrasts(contrasts = ex, levels = design); colnames(cm) <- nm
  eb <- eBayes(contrasts.fit(lmFit(m, design), cm), robust = TRUE, trend = FALSE)
  sapply(c(SUS_vs_RES = paste0(unit, "_3_over_", unit, "_2"), SUS_vs_CON = paste0(unit, "_3_over_", unit, "_1"),
           RES_vs_CON = paste0(unit, "_2_over_", unit, "_1")),
         function(cc) sum(topTable(eb, coef = cc, number = Inf, sort.by = "none", adjust.method = "BH")$adj.P.Val < 0.05))
}

# ------------------------------------------------------------------------------------------------ A. CA2-SLM control
if (any(PART %in% c("A", "all"))) {
  cat("== A: CA2-SLM removal control\n")
  ds <- "neuron_neuropil"
  o <- OFF[OFF$dataset == ds, ]
  imp <- read_imputed(Sys.glob(sprintf("data/processed/01_preprocessing/impute/*_pgmatrix_imputed_%s_*_missing70pct.xlsx", ds)))
  mn <- imp$m[, o$sample_id]
  g <- read_gct(sprintf("data/processed/01_preprocessing/protigy_input_animal_level/%s/%s_animal_level.gct", ds, ds))
  col_key <- paste(g$coldesc$AnimalID, g$coldesc$region, g$coldesc$layer, sep = "|")
  smp_key <- paste(o$animal, sub("_.*$", "", o$unit), sub("^[^_]+_", "", o$unit), sep = "|")
  pn <- imp$names
  ca2 <- o[o$unit == "CA2_slm", ]
  cell <- g$coldesc$region_layer_ExpGroup
  rows <- list()
  for (grp in c("SUS", "RES", "CON")) {
    an <- sort(unique(ca2$animal[ca2$group == grp]))
    combos <- list()
    for (pr in combn(an, 2, simplify = FALSE)) for (h1 in c("L", "R")) for (h2 in c("L", "R"))
      combos[[length(combos) + 1L]] <- c(paste(pr[1], h1, "CA2_slm", sep = "_"), paste(pr[2], h2, "CA2_slm", sep = "_"))
    for (cb in combos) {
      keep <- !(o$label %in% cb)
      m0 <- g$m
      for (kk in unique(smp_key[o$label %in% cb])) {
        v1 <- rowMeans(mn[, o$sample_id[smp_key == kk & keep], drop = FALSE]); names(v1) <- pn
        v1 <- v1[!duplicated(names(v1))]; m0[, which(col_key == kk)] <- v1[rownames(m0)]
      }
      for (k in EST3) {
        shift <- vapply(col_key, function(kk) mean((o$raw_new_complete_case - o[[paste0("raw_", k)]])[smp_key == kk & keep]), numeric(1))
        n <- fit_counts(sweep(m0, 2, shift, `+`), cell)
        rows[[length(rows) + 1L]] <- tibble(estimator = k, removed_group = grp, removed = paste(cb, collapse = " + "),
          is_A755L_A764R = setequal(cb, c("A755_L_CA2_slm", "A764_R_CA2_slm")),
          removed_missing_fraction = paste(round(o$missing_fraction[match(cb, o$label)], 3), collapse = " + "),
          SUS_vs_RES = n[["SUS_vs_RES"]], SUS_vs_CON = n[["SUS_vs_CON"]], RES_vs_CON = n[["RES_vs_CON"]])
      }
    }
  }
  full <- bind_rows(lapply(EST3, function(k) {
    shift <- vapply(col_key, function(kk) mean((o$raw_new_complete_case - o[[paste0("raw_", k)]])[smp_key == kk]), numeric(1))
    n <- fit_counts(sweep(g$m, 2, shift, `+`), cell)
    tibble(estimator = k, removed_group = "none", removed = "none (all 18 samples)", is_A755L_A764R = FALSE, removed_missing_fraction = "",
           SUS_vs_RES = n[["SUS_vs_RES"]], SUS_vs_CON = n[["SUS_vs_CON"]], RES_vs_CON = n[["RES_vs_CON"]])
  }))
  res <- bind_rows(full, bind_rows(rows))
  w(res, "ca2slm_removal_control.csv")
  print(res |> filter(removed_group != "none") |> group_by(estimator, removed_group) |>
          summarise(n = n(), SUS_vs_RES_median = median(SUS_vs_RES), SUS_vs_RES_min = min(SUS_vs_RES), SUS_vs_RES_max = max(SUS_vs_RES),
                    focus_pair = SUS_vs_RES[is_A755L_A764R][1], .groups = "drop"))
}

# ------------------------------------------------------------------------------------------------ A2. CA2-SLM detail
# Which current hits depend on the two lowest-depth SUS hemispheres; single-hemisphere removals; a CA2-SLM-only model
# (variance from that unit alone); logFC concordance of every removal variant; the padj distribution of the hits.
if (any(PART %in% c("A2", "all"))) {
  cat("== A2: CA2-SLM detail\n")
  ds <- "neuron_neuropil"
  o <- OFF[OFF$dataset == ds, ]
  imp <- read_imputed(Sys.glob(sprintf("data/processed/01_preprocessing/impute/*_pgmatrix_imputed_%s_*_missing70pct.xlsx", ds)))
  mn <- imp$m[, o$sample_id]; rownames(mn) <- imp$names
  rawm <- read_tsv("data/raw/pg_matrix/quicksearch.pg_matrix.tsv", show_col_types = FALSE, progress = FALSE)
  pgid <- read_excel(Sys.glob(sprintf("data/processed/01_preprocessing/impute/*_pgmatrix_imputed_%s_*_missing70pct.xlsx", ds)), progress = FALSE)$Protein.Group
  obs <- !is.na(as.matrix(rawm[match(pgid, rawm$Protein.Group), o$sample_id])); rownames(obs) <- imp$names
  ca2 <- o$unit == "CA2_slm"; focus <- c("A755_L_CA2_slm", "A764_R_CA2_slm")
  da_new <- "results/preprocessing/map_protein_identifiers/neuron_neuropil/tables/mapped/forward/per_file/CA2slmsus_CA2slmres.csv"
  da_old <- file.path(ARCHIVE, "data/processed/02_id_mapping/mapped/neuron_neuropil/forward/per_file/CA2slmsus_CA2slmres.csv")
  hits <- bind_rows(lapply(list(new = da_new, original = da_old), function(f) {
    d <- read_csv(f, show_col_types = FALSE); d <- d[!is.na(d$padj) & d$padj < 0.05, ]
    id <- d$original_identifier; sus <- ca2 & o$group == "SUS"
    others <- sus & !(o$label %in% focus)
    tibble(protein = id, gene_symbol = d$official_gene_symbol %||% NA_character_, log2fc = d$log2fc, padj = d$padj,
           imputed_A755L = !obs[id, o$label == "A755_L_CA2_slm"], imputed_A764R = !obs[id, o$label == "A764_R_CA2_slm"],
           n_imputed_SUS = rowSums(!obs[id, sus, drop = FALSE]), n_imputed_RES = rowSums(!obs[id, ca2 & o$group == "RES", drop = FALSE]),
           detection_fraction_neuropil = rowMeans(obs[id, , drop = FALSE]),
           detection_class = cut(rowMeans(obs[id, , drop = FALSE]), c(-Inf, 0.5, 0.95, 0.9999, Inf), labels = c("under 50%", "50-95% (partially observed)", "95-99.9%", "all samples (anchor)")),
           complete_in_CA2slm = rowSums(!obs[id, ca2, drop = FALSE]) == 0,
           A755L_minus_other_SUS = mn[id, o$label == "A755_L_CA2_slm"] - rowMeans(mn[id, others, drop = FALSE]),
           A764R_minus_other_SUS = mn[id, o$label == "A764_R_CA2_slm"] - rowMeans(mn[id, others, drop = FALSE]))
  }), .id = "version")
  w(hits, "ca2slm_hits_sample_dependence.csv")
  # CA2-SLM-only model, and single/paired removals with logFC concordance, under three centrings
  g <- read_gct(sprintf("data/processed/01_preprocessing/protigy_input_animal_level/%s/%s_animal_level.gct", ds, ds))
  col_key <- paste(g$coldesc$AnimalID, g$coldesc$region, g$coldesc$layer, sep = "|")
  smp_key <- paste(o$animal, sub("_.*$", "", o$unit), sub("^[^_]+_", "", o$unit), sep = "|")
  cell <- g$coldesc$region_layer_ExpGroup
  fit_full <- function(m) {
    f <- factor(cell, levels = unique(cell)); design <- model.matrix(~0 + f); colnames(design) <- levels(f)
    cm <- makeContrasts(contrasts = "CA2_slm_3-CA2_slm_2", levels = design)
    topTable(eBayes(contrasts.fit(lmFit(m, design), cm), robust = TRUE, trend = FALSE), coef = 1, number = Inf, sort.by = "none")
  }
  variant_matrix <- function(keep, k) {
    m0 <- g$m
    for (kk in unique(smp_key[!keep])) {
      v1 <- rowMeans(mn[, o$sample_id[smp_key == kk & keep], drop = FALSE]); v1 <- v1[!duplicated(names(v1))]
      m0[, which(col_key == kk)] <- v1[rownames(m0)]
    }
    shift <- vapply(col_key, function(kk) mean((o$raw_new_complete_case - o[[paste0("raw_", k)]])[smp_key == kk & keep]), numeric(1))
    sweep(m0, 2, shift, `+`)
  }
  variants <- list("all samples" = rep(TRUE, nrow(o)), "without A755 L" = o$label != "A755_L_CA2_slm",
                   "without A764 R" = o$label != "A764_R_CA2_slm", "without A755 R (A755 L alone)" = o$label != "A755_R_CA2_slm",
                   "without A755 L and A764 R" = !(o$label %in% focus))
  rows <- list()
  for (k in EST3) {
    base <- fit_full(variant_matrix(rep(TRUE, nrow(o)), k))
    for (vn in names(variants)) {
      tt <- if (vn == "all samples") base else fit_full(variant_matrix(variants[[vn]], k))
      # CA2-SLM-only model on the same animal-level values
      m <- variant_matrix(variants[[vn]], k)[, grepl("^CA2_slm_", cell)]; gg <- factor(sub("^CA2_slm_", "", cell[grepl("^CA2_slm_", cell)]))
      d <- model.matrix(~0 + gg); colnames(d) <- paste0("g", levels(gg))
      t2 <- topTable(eBayes(contrasts.fit(lmFit(m, d), makeContrasts(g3 - g2, levels = d)), robust = TRUE), coef = 1, number = Inf, sort.by = "none")
      rows[[length(rows) + 1L]] <- tibble(estimator = k, variant = vn, n_fdr05_pipeline_model = sum(tt$adj.P.Val < 0.05),
        n_fdr05_ca2slm_only_model = sum(t2$adj.P.Val < 0.05), median_logFC_all_proteins = median(tt$logFC),
        cor_logFC_vs_all_samples = stats::cor(tt$logFC, base$logFC), spearman_t_vs_all_samples = stats::cor(tt$t, base$t, method = "spearman"),
        n_padj_0.035_0.05 = sum(tt$adj.P.Val >= 0.035 & tt$adj.P.Val < 0.05))
    }
  }
  res <- bind_rows(rows); w(res, "ca2slm_variants_and_unit_only_model.csv"); print(as.data.frame(res))
}

# ------------------------------------------------------------------------------------------------ B. anchor composition
if (any(PART %in% c("B", "all"))) {
  cat("== B: anchor composition\n")
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Mm.eg.db) })
  raw <- read_tsv("data/raw/pg_matrix/quicksearch.pg_matrix.tsv", show_col_types = FALSE, col_select = 1:4, progress = FALSE)
  sym <- setNames(sub(";.*$", "", raw$Genes), raw$Protein.Group)
  res <- list()
  for (ds in DATASETS) {
    imp <- read_imputed(Sys.glob(sprintf("data/processed/01_preprocessing/impute/*_pgmatrix_imputed_%s_*_missing70pct.xlsx", ds)))
    filtered <- rownames(imp$m)
    anchors <- intersect(state$anchor_sets[[ds]], filtered)
    universe <- unique(na.omit(sym[filtered])); universe <- universe[nzchar(universe)]
    for (set in c("anchor", "non_anchor")) {
      genes <- unique(na.omit(sym[if (set == "anchor") anchors else setdiff(filtered, anchors)])); genes <- genes[nzchar(genes)]
      for (ont in c("CC", "BP")) {
        e <- tryCatch(suppressWarnings(enrichGO(genes, universe = universe, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = ont,
                                                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)), error = function(e) NULL)
        if (is.null(e)) next
        d <- as.data.frame(e) |> filter(p.adjust < 0.05) |> arrange(p.adjust) |> head(12)
        if (nrow(d)) res[[length(res) + 1L]] <- d |> transmute(dataset = ds, set, n_set = length(genes), n_universe = length(universe),
          Ontology = ont, ID, Description, GeneRatio, BgRatio, p.adjust,
          fold = (as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))) /
                 (as.numeric(sub("/.*", "", BgRatio)) / as.numeric(sub(".*/", "", BgRatio))))
      }
    }
  }
  w(bind_rows(res), "anchor_functional_composition.csv")
}

# ------------------------------------------------------------------------------------------------ C. WGCNA stability
if (any(PART %in% c("C", "all"))) {
  cat("== C: WGCNA stability\n")
  suppressPackageStartupMessages(library(WGCNA))
  enableWGCNAThreads(nThreads = as.integer(arg("--threads", "8")))
  cache <- file.path(OUT, "wgcna_builds"); dir.create(cache, showWarnings = FALSE)
  build <- function(x, power) {
    adj <- adjacency(x, power = power, type = "signed", corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = 0.05))
    tom <- TOMsimilarity(adj, verbose = 0); rm(adj); invisible(gc())
    diss <- 1 - tom; rm(tom); invisible(gc())
    tree <- hclust(as.dist(diss), method = "average")
    mods <- cutreeDynamic(dendro = tree, distM = diss, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30, verbose = 0)
    rm(diss); invisible(gc())
    merged <- mergeCloseModules(x, labels2colors(mods), cutHeight = 0.25, verbose = 0)$colors
    setNames(as.character(merged), colnames(x))
  }
  ari <- function(a, b) {
    ids <- intersect(names(a), names(b)); a <- a[ids]; b <- b[ids]
    tab <- table(a, b); n <- sum(tab); c2 <- function(v) sum(v * (v - 1) / 2)
    sij <- c2(tab); sa <- c2(rowSums(tab)); sb <- c2(colSums(tab)); e <- sa * sb / (n * (n - 1) / 2)
    (sij - e) / ((sa + sb) / 2 - e)
  }
  n_mod <- function(x) length(setdiff(unique(x), "grey"))
  cached <- function(key, expr) { f <- file.path(cache, paste0(key, ".rds")); if (file.exists(f)) return(readRDS(f)); v <- force(expr); saveRDS(v, f); v }
  summ <- list(); loao <- list()
  for (ds in DATASETS) {
    sn <- readRDS(file.path("results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"))
    so <- readRDS(file.path(ARCHIVE, "data", "processed", "06_modules_WGCNA", "01_WGCNA", ds, "wgcna_final_model_state.rds"))
    xn <- as.matrix(sn$expression.data); storage.mode(xn) <- "double"
    col_new <- setNames(as.character(sn$mergedColors), colnames(xn))
    col_old <- setNames(as.character(so$mergedColors), colnames(so$expression.data))
    o <- OFF[OFF$dataset == ds, ]; o <- o[match(rownames(xn), o$sample_id), ]
    stopifnot(!anyNA(o$sample_id))
    for (k in EST3) {
      xk <- xn + (o$raw_new_complete_case - o[[paste0("raw_", k)]])   # row-wise: one constant per sample
      fi <- cached(sprintf("%s_%s_softpower", ds, k), {
        spt <- pickSoftThreshold(xk, networkType = "signed", corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = 0.05), verbose = 0)
        as.data.frame(spt$fitIndices)
      })
      el <- fi$Power[fi$SFT.R.sq >= 0.80]; pw <- if (length(el)) min(el) else fi$Power[which.max(fi$SFT.R.sq)]
      cat(sprintf("  %s %s power %d\n", ds, k, pw))
      full <- cached(sprintf("%s_%s_full", ds, k), build(xk, pw))
      me <- moduleEigengenes(xk, colors = full)$eigengenes
      me <- me[, colnames(me) != "MEgrey", drop = FALSE]
      rme <- vapply(me, function(v) stats::cor(v, o$missing_fraction), numeric(1))
      animals <- sort(unique(o$animal))
      ar <- vapply(animals, function(a) {
        keep <- o$animal != a
        cl <- cached(sprintf("%s_%s_loao_%s", ds, k, a), build(xk[keep, , drop = FALSE], pw))
        ari(cl, full)
      }, numeric(1))
      loao[[length(loao) + 1L]] <- tibble(dataset = ds, estimator = k, left_out_animal = animals, ari_vs_full_build = ar)
      summ[[length(summ) + 1L]] <- tibble(dataset = ds, estimator = k, power = pw, n_samples = nrow(xk), n_proteins = ncol(xk),
        n_modules = n_mod(full), grey_fraction = mean(full == "grey"),
        ari_vs_current_new_network = ari(full, col_new), ari_vs_previous_network = ari(full, col_old),
        loao_ari_mean = mean(ar), loao_ari_min = min(ar),
        max_abs_cor_ME_vs_missing = max(abs(rme)), median_abs_cor_ME_vs_missing = median(abs(rme)),
        n_modules_ME_cor_missing_gt_0.5 = sum(abs(rme) > 0.5))
    }
  }
  if (!"C2" %in% PART) {
    w(bind_rows(summ), "wgcna_stability_by_estimator.csv")
    w(bind_rows(loao), "wgcna_leave_one_animal_out.csv")
    print(bind_rows(summ))
  }
  if ("C2" %in% PART) {
    # common power = the current (new-centring) network's selected power
    summ2 <- list(); loao2 <- list()
    for (ds in DATASETS) {
      sn <- readRDS(file.path("results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"))
      xn <- as.matrix(sn$expression.data); storage.mode(xn) <- "double"
      col_new <- setNames(as.character(sn$mergedColors), colnames(xn))
      pw <- as.numeric(sn$softPower)
      o <- OFF[OFF$dataset == ds, ]; o <- o[match(rownames(xn), o$sample_id), ]
      for (k in EST3) {
        xk <- xn + (o$raw_new_complete_case - o[[paste0("raw_", k)]])
        own <- readRDS(file.path(cache, sprintf("%s_%s_softpower.rds", ds, k)))
        el <- own$Power[own$SFT.R.sq >= 0.80]; pw_own <- if (length(el)) min(el) else own$Power[which.max(own$SFT.R.sq)]
        tag <- if (pw_own == pw) "" else sprintf("_p%d", pw)
        cat(sprintf("  C2 %s %s common power %d (own %d)
", ds, k, pw, pw_own))
        full <- cached(sprintf("%s_%s_full%s", ds, k, tag), build(xk, pw))
        animals <- sort(unique(o$animal))
        ar <- vapply(animals, function(a) {
          keep <- o$animal != a
          ari(cached(sprintf("%s_%s_loao_%s%s", ds, k, a, tag), build(xk[keep, , drop = FALSE], pw)), full)
        }, numeric(1))
        loao2[[length(loao2) + 1L]] <- tibble(dataset = ds, estimator = k, power = pw, left_out_animal = animals, ari_vs_full_build = ar)
        summ2[[length(summ2) + 1L]] <- tibble(dataset = ds, estimator = k, common_power = pw, own_power = pw_own,
          n_modules = n_mod(full), ari_vs_current_new_network = ari(full, col_new), loao_ari_mean = mean(ar), loao_ari_min = min(ar))
      }
    }
    l2 <- bind_rows(loao2)
    paired <- l2 |> select(dataset, estimator, left_out_animal, ari_vs_full_build) |>
      tidyr::pivot_wider(names_from = estimator, values_from = ari_vs_full_build) |>
      mutate(new_minus_original = new_complete_case - original_all_observed,
             median_ratio_minus_new = median_ratio_complete_case - new_complete_case)
    w(bind_rows(summ2), "wgcna_stability_common_power.csv")
    w(l2, "wgcna_leave_one_animal_out_common_power.csv")
    w(paired, "wgcna_leave_one_animal_out_common_power_paired.csv")
    print(bind_rows(summ2)); print(paired |> group_by(dataset) |> summarise(animals_new_better = sum(new_minus_original > 0), mean_diff = mean(new_minus_original)))
  }
}
cat("done\n")
