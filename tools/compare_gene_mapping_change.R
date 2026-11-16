#!/usr/bin/env Rscript
# Before/after comparison for a change of protein-to-gene annotation (manual gene annotation overrides).
#
#   snapshot: Rscript tools/compare_gene_mapping_change.R snapshot --root <tree> --out <dir> [--datasets a,b,c]
#             copies the outputs compared below (same relative paths) and writes snapshot_manifest.csv
#             (sha256 of every file, git commit of the tree, sha256 of the override file).
#   compare:  Rscript tools/compare_gene_mapping_change.R compare --before <root> --after <root> --out <dir> [--datasets ...]
#             writes one CSV per layer and comparison_summary.csv.
#
# Layers compared, per dataset:
#   1. Stage 02 per-contrast tables: gene-level inclusion per protein group, and the protein-level statistics
#      (log2fc, t, pval, padj), which an annotation change must leave identical.
#   2. clusterProfiler per comparison: gene universe (collapsed_gene_input.csv), GO GSEA (BP), enrichGO and KEGG GSEA:
#      every term's NES/p/adjusted p before and after, and terms gaining or losing FDR significance.
#   3. EWCE: background universe and every test (p, q, fold change), and q-significance changes.
#   4. WGCNA: module membership (must not change) and module GO (universe size, adjusted p, significance changes).
# A term is significant at adjusted p (or q) < 0.05.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else ""
arg <- function(flag, default = NA_character_) { i <- which(args == flag); if (length(i) && i[1] < length(args)) args[[i[1] + 1]] else default }
datasets <- strsplit(arg("--datasets", "neuron_neuropil,neuron_soma,microglia"), ",")[[1]]
alpha <- 0.05

layer_files <- function(root, ds) {
  rel <- function(p) sub(paste0("^", gsub("([.\\\\+*?^$(){}|\\[\\]])", "\\\\\\1", normalizePath(root, winslash = "/")), "/?"), "", normalizePath(p, winslash = "/"))
  f <- c(
    list.files(file.path(root, "results/preprocessing/map_protein_identifiers", ds, "tables/mapped/forward/per_file"), "[.]csv$", full.names = TRUE),
    list.files(file.path(root, "results/differential_abundance/run_clusterprofiler_enrichment", ds, "models"), recursive = TRUE, full.names = TRUE,
               pattern = "^(GSEA_BP_results_full|enrichGO_ALL_results_full|KEGG_GSEA_results|collapsed_gene_input)[.]csv$"),
    list.files(file.path(root, "results/enrichment/run_ewce_celltype_enrichment", ds, "models"), "^EWCE_results_full[.]rds$", full.names = TRUE),
    list.files(file.path(root, "results/wgcna/build_wgcna_modules", ds, "tables/modules"), "^(genes_in_module_.*|WGCNA_module_GO_enrichment_long|WGCNA_module_GO_enrichment_QC)[.]csv$", full.names = TRUE))
  data.table(abs = f, rel = vapply(f, rel, ""))
}

if (identical(mode, "snapshot")) {
  root <- arg("--root"); out <- arg("--out"); stopifnot(!is.na(root), !is.na(out), dir.exists(root))
  if (dir.exists(out) && length(list.files(out))) stop("snapshot folder is not empty: ", out, call. = FALSE)
  files <- rbindlist(lapply(datasets, function(ds) cbind(dataset = ds, layer_files(root, ds))))
  # R on this host cannot open paths of 260+ characters: keep the snapshot root short
  too_long <- nchar(file.path(normalizePath(out, winslash = "/", mustWork = FALSE), files$rel)) >= 255
  if (any(too_long)) stop(sum(too_long), " snapshot paths would exceed 254 characters; use a shorter --out", call. = FALSE)
  for (i in seq_len(nrow(files))) { dst <- file.path(out, files$rel[i]); dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(files$abs[i], dst, copy.date = TRUE)) stop("copy failed: ", files$abs[i], call. = FALSE) }
  sha <- function(f) digest::digest(file = f, algo = "sha256")
  ov <- file.path(root, "data/metadata/manual_gene_annotation_overrides.csv")
  commit <- tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_)
  man <- files[, .(dataset, rel, bytes = file.size(abs), sha256 = vapply(abs, sha, ""))]
  man[, `:=`(source_root = normalizePath(root, winslash = "/"), git_commit = commit[1], override_file_sha256 = if (file.exists(ov)) sha(ov) else "absent",
             snapshot_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))]
  fwrite(man, file.path(out, "snapshot_manifest.csv"))
  cat("snapshot:", nrow(man), "files,", round(sum(man$bytes) / 1e6, 1), "MB ->", out, "\n")
  quit(status = 0, save = "no")
}
if (!identical(mode, "compare")) stop("usage: compare_gene_mapping_change.R snapshot|compare ...", call. = FALSE)
before <- arg("--before"); after <- arg("--after"); out <- arg("--out")
stopifnot(!is.na(before), !is.na(after), !is.na(out)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) {
  if (!file.exists(f)) return(NULL)
  x <- fread(f, showProgress = FALSE)
  if ("official_entrez_id" %in% names(x)) x[, official_entrez_id := as.character(official_entrez_id)]
  x
}
sig_change <- function(b, a, key, pcol) {
  m <- merge(b[, c(key, pcol), with = FALSE], a[, c(key, pcol), with = FALSE], by = key, all = TRUE, suffixes = c("_before", "_after"))
  pb <- m[[paste0(pcol, "_before")]]; pa <- m[[paste0(pcol, "_after")]]
  list(m = m, n_sig_before = sum(pb < alpha, na.rm = TRUE), n_sig_after = sum(pa < alpha, na.rm = TRUE),
       gained = sum((is.na(pb) | pb >= alpha) & pa < alpha, na.rm = TRUE), lost = sum(pb < alpha & (is.na(pa) | pa >= alpha), na.rm = TRUE),
       only_before = sum(!is.na(pb) & is.na(pa)), only_after = sum(is.na(pb) & !is.na(pa)))
}
summ <- list(); s02 <- list(); cp <- list(); ew <- list(); wg <- list()
for (ds in datasets) {
  # 1. Stage 02
  for (f in layer_files(after, ds)[grepl("map_protein_identifiers", rel)]$rel) {
    b <- rd(file.path(before, f)); a <- rd(file.path(after, f)); if (is.null(b) || is.null(a)) next
    m <- merge(b[, .(ProteinGroupID, g_b = gene_level_claim_allowed, s_b = official_gene_symbol, lfc_b = log2fc, t_b = t, p_b = pval, q_b = padj)],
               a[, .(ProteinGroupID, g_a = gene_level_claim_allowed, s_a = official_gene_symbol, lfc_a = log2fc, t_a = t, p_a = pval, q_a = padj)], by = "ProteinGroupID", all = TRUE)
    stat_diff <- suppressWarnings(max(abs(c(m$lfc_a - m$lfc_b, m$t_a - m$t_b, m$p_a - m$p_b, m$q_a - m$q_b)), na.rm = TRUE))
    s02[[length(s02) + 1]] <- data.table(dataset = ds, comparison = sub("[.]csv$", "", basename(f)), pgs_before = nrow(b), pgs_after = nrow(a),
      gene_level_before = sum(b$gene_level_claim_allowed %in% TRUE), gene_level_after = sum(a$gene_level_claim_allowed %in% TRUE),
      newly_gene_level = paste(m[!(g_b %in% TRUE) & g_a %in% TRUE, s_a], collapse = ";"), no_longer_gene_level = paste(m[g_b %in% TRUE & !(g_a %in% TRUE), s_b], collapse = ";"),
      symbol_changed = sum(m$g_b %in% TRUE & m$g_a %in% TRUE & m$s_b != m$s_a, na.rm = TRUE),
      protein_level_stats_max_abs_diff = stat_diff, protein_sets_identical = setequal(b$ProteinGroupID, a$ProteinGroupID))
  }
  # 2. clusterProfiler
  cdirs <- unique(dirname(dirname(dirname(layer_files(after, ds)[grepl("GSEA_BP_results_full", rel)]$rel))))
  for (cd in cdirs) {
    u_b <- rd(file.path(before, cd, "protein_group_audits/collapsed_gene_input.csv")); u_a <- rd(file.path(after, cd, "protein_group_audits/collapsed_gene_input.csv"))
    for (tt in list(c("GO/BP/GSEA_BP_results_full.csv", "GSEA_GO_BP"), c("GO/BP/enrichGO_ALL_results_full.csv", "enrichGO"), c("KEGG/KEGG_GSEA_results.csv", "GSEA_KEGG"))) {
      b <- rd(file.path(before, cd, tt[1])); a <- rd(file.path(after, cd, tt[1])); if (is.null(b) || is.null(a) || !"p.adjust" %in% names(b)) next
      s <- sig_change(b, a, "ID", "p.adjust")
      nes <- if ("NES" %in% names(b)) { mm <- merge(b[, .(ID, NES)], a[, .(ID, NES)], by = "ID"); suppressWarnings(max(abs(mm$NES.x - mm$NES.y), na.rm = TRUE)) } else NA
      cp[[length(cp) + 1]] <- data.table(dataset = ds, comparison = basename(cd), test = tt[2],
        universe_before = if (is.null(u_b)) NA else nrow(u_b), universe_after = if (is.null(u_a)) NA else nrow(u_a),
        terms_before = nrow(b), terms_after = nrow(a), sig_before = s$n_sig_before, sig_after = s$n_sig_after,
        sig_gained = s$gained, sig_lost = s$lost, terms_only_before = s$only_before, terms_only_after = s$only_after, max_abs_delta_NES = nes)
      fwrite(s$m, file.path(out, paste0("terms_", ds, "_", basename(cd), "_", tt[2], ".csv")))
    }
  }
  # 3. EWCE
  fb <- file.path(before, "results/enrichment/run_ewce_celltype_enrichment", ds, "models/EWCE_results_full.rds"); fa <- sub(before, after, fb, fixed = TRUE)
  if (file.exists(fb) && file.exists(file.path(after, "results/enrichment/run_ewce_celltype_enrichment", ds, "models/EWCE_results_full.rds"))) {
    xb <- readRDS(fb); xa <- readRDS(file.path(after, "results/enrichment/run_ewce_celltype_enrichment", ds, "models/EWCE_results_full.rds"))
    key <- intersect(c("CellType", "annotLevel", "Target", "AnalysisType", "Stratum", "Contrast", "Direction", "Metric", "TopN"), names(xb$results_all))
    b <- as.data.table(xb$results_all); a <- as.data.table(xa$results_all)
    s <- sig_change(b, a, key, "q")
    ew[[length(ew) + 1]] <- data.table(dataset = ds, background_before = length(xb$background_universe), background_after = length(xa$background_universe),
      background_added = paste(setdiff(xa$background_universe, xb$background_universe), collapse = ";"),
      background_removed = paste(setdiff(xb$background_universe, xa$background_universe), collapse = ";"),
      tests_before = nrow(b), tests_after = nrow(a), sig_before = s$n_sig_before, sig_after = s$n_sig_after, sig_gained = s$gained, sig_lost = s$lost,
      max_abs_delta_p = suppressWarnings(max(abs(s$m$q_after - s$m$q_before), na.rm = TRUE)))
    fwrite(s$m, file.path(out, paste0("ewce_tests_", ds, ".csv")))
  }
  # 4. WGCNA
  mods <- layer_files(after, ds)[grepl("genes_in_module_", rel)]$rel
  memb <- function(root) rbindlist(lapply(mods, function(f) { x <- rd(file.path(root, f)); if (is.null(x)) NULL else data.table(module = basename(f), pg = x$ProteinGroupID) }))
  mb <- memb(before); ma <- memb(after)
  gob <- rd(file.path(before, "results/wgcna/build_wgcna_modules", ds, "tables/modules/WGCNA_module_GO_enrichment_long.csv"))
  goa <- rd(file.path(after, "results/wgcna/build_wgcna_modules", ds, "tables/modules/WGCNA_module_GO_enrichment_long.csv"))
  if (!is.null(gob) && !is.null(goa)) {
    s <- sig_change(gob, goa, c("ModuleID", "ModuleProteinSetType", "Ontology", "ID"), "p.adjust")
    wg[[length(wg) + 1]] <- data.table(dataset = ds, membership_identical = isTRUE(all.equal(mb[order(module, pg)], ma[order(module, pg)], check.attributes = FALSE)),
      go_universe_before = paste(unique(gob$UniverseSize), collapse = ";"), go_universe_after = paste(unique(goa$UniverseSize), collapse = ";"),
      sig_before = s$n_sig_before, sig_after = s$n_sig_after, sig_gained = s$gained, sig_lost = s$lost)
    fwrite(s$m, file.path(out, paste0("wgcna_module_go_", ds, ".csv")))
  }
}
w <- function(x, name) { d <- rbindlist(x, fill = TRUE); fwrite(d, file.path(out, name)); d }
S02 <- w(s02, "stage02_inclusion_and_protein_stats.csv"); CP <- w(cp, "clusterprofiler_terms.csv"); EW <- w(ew, "ewce.csv"); WG <- w(wg, "wgcna.csv")
summary <- data.table(
  layer = c("Stage 02 protein-level statistics", "Stage 02 gene-level inclusion", "clusterProfiler universes", "clusterProfiler FDR significance",
            "EWCE background", "EWCE q significance", "WGCNA membership", "WGCNA module GO significance"),
  result = c(
    if (nrow(S02)) sprintf("max |difference| %.3g over %d tables; protein sets identical in %d/%d", max(S02$protein_level_stats_max_abs_diff), nrow(S02), sum(S02$protein_sets_identical), nrow(S02)) else "n/a",
    if (nrow(S02)) sprintf("gene-level groups %d -> %d (summed over tables); symbol changes %d", sum(S02$gene_level_before), sum(S02$gene_level_after), sum(S02$symbol_changed)) else "n/a",
    if (nrow(CP)) sprintf("%d comparisons; universe change range %s", uniqueN(CP$comparison), paste(range(CP$universe_after - CP$universe_before, na.rm = TRUE), collapse = " to ")) else "n/a",
    if (nrow(CP)) sprintf("significant terms %d -> %d; gained %d, lost %d", sum(CP$sig_before), sum(CP$sig_after), sum(CP$sig_gained), sum(CP$sig_lost)) else "n/a",
    if (nrow(EW)) paste(sprintf("%s %d -> %d", EW$dataset, EW$background_before, EW$background_after), collapse = "; ") else "n/a",
    if (nrow(EW)) sprintf("significant %d -> %d; gained %d, lost %d", sum(EW$sig_before), sum(EW$sig_after), sum(EW$sig_gained), sum(EW$sig_lost)) else "n/a",
    if (nrow(WG)) paste(sprintf("%s %s", WG$dataset, ifelse(WG$membership_identical, "identical", "CHANGED")), collapse = "; ") else "n/a",
    if (nrow(WG)) sprintf("significant %d -> %d; gained %d, lost %d", sum(WG$sig_before), sum(WG$sig_after), sum(WG$sig_gained), sum(WG$sig_lost)) else "n/a"))
fwrite(summary, file.path(out, "comparison_summary.csv"))
print(summary)
