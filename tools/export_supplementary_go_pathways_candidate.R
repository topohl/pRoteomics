#!/usr/bin/env Rscript
# Downstream-only, one-time source-data export for a supplementary figure candidate.
# Run from the pRoteomics repository root. No enrichment test or DA model is fitted.

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "pipeline.yml")))
  stop("Run from the pRoteomics repository root.", call. = FALSE)
if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
    !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  stop("AnnotationDbi and org.Mm.eg.db are required to verify GO membership.",
       call. = FALSE)
}

rel <- function(...) file.path(root, ...)
inventory_dir <- rel("exports", "publication_source_data",
                     "supplementary_selection_inventories")
supported_path <- file.path(inventory_dir,
                            "pathway_enrichment_inventory_fdr_supported.csv")
all_path <- file.path(inventory_dir, "pathway_enrichment_inventory.csv")
stopifnot(file.exists(supported_path), file.exists(all_path))
out <- rel("exports", "publication_source_data",
           "supplementary_go_pathways_candidate")
if (file.exists(out)) stop("Candidate export already exists: ", out,
                           call. = FALSE)

# Exact terms name the four remaining atlas families directly. GO:0006914 is
# the autophagy child of the registry's GO:0061919 core anchor.
spec <- data.frame(
  display_order = 1:4,
  theme_id = c("ribosome_translation", "chromatin_organization",
               "neuron_projection_development", "autophagy_lysosome_endosome"),
  go_id = c("GO:0006412", "GO:0006325", "GO:0031175", "GO:0006914"),
  stringsAsFactors = FALSE)
supported <- utils::read.csv(supported_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
all <- utils::read.csv(all_path, stringsAsFactors = FALSE,
                       check.names = FALSE)
required <- c("dataset", "spatial_unit", "contrast", "GO_ID",
              "GO_description", "NES", "BH_FDR", "theme_id")
if (!all(required %in% names(supported)) || !all(required %in% names(all)))
  stop("Frozen pathway inventory lacks required columns.", call. = FALSE)

selection <- vector("list", nrow(spec))
curves <- vector("list", nrow(spec))
inputs <- data.frame(role = c("frozen_full_inventory", "frozen_supported_inventory"),
                     path = c(all_path, supported_path), stringsAsFactors = FALSE)
for (i in seq_len(nrow(spec))) {
  s <- spec[i, ]
  candidates <- supported[supported$theme_id == s$theme_id &
                            supported$GO_ID == s$go_id &
                            supported$contrast == "SUS - RES", , drop = FALSE]
  if (!nrow(candidates))
    stop("No FDR-supported SUS - RES occurrence: ", s$go_id, call. = FALSE)
  candidates <- candidates[order(candidates$BH_FDR, candidates$dataset,
                                 candidates$spatial_unit), , drop = FALSE]
  pick <- candidates[1, , drop = FALSE]
  ds <- pick$dataset[[1]]; unit <- pick$spatial_unit[[1]]
  unit_dir <- if (identical(ds, "microglia")) paste0(unit, "_microglia") else unit
  tok <- gsub("_", "", unit)
  if (identical(ds, "microglia")) tok <- paste0(tok, "microglia")
  cmp <- sprintf("%ssus_%sres", tok, tok)
  base <- rel("data", "processed", "04_differential_expression_enrichment",
              "clusterProfiler", ds, "phenotype_within_unit", unit_dir, cmp)
  ranked_path <- file.path(base, "protein_group_audits", "collapsed_gene_input.csv")
  gsea_path <- file.path(base, "GO", "BP", "GSEA_BP_results_full.csv")
  if (!file.exists(ranked_path) || !file.exists(gsea_path))
    stop("Canonical GSEA inputs missing for ", ds, " ", unit, call. = FALSE)
  inputs <- rbind(inputs, data.frame(role = c("ranked_gene_input", "stored_gsea_result"),
                                    path = c(ranked_path, gsea_path)))
  g <- utils::read.csv(ranked_path, stringsAsFactors = FALSE)
  res <- utils::read.csv(gsea_path, stringsAsFactors = FALSE)
  if (!all(c("GeneSymbol", "collapsed_statistic") %in% names(g)) ||
      !all(c("ID", "setSize", "enrichmentScore", "NES", "p.adjust") %in% names(res)))
    stop("Canonical GSEA schema mismatch: ", ds, " ", unit, call. = FALSE)
  if (anyDuplicated(g$GeneSymbol) || anyNA(g$GeneSymbol) ||
      anyNA(g$collapsed_statistic))
    stop("Ranked gene input is not one finite row per symbol.", call. = FALSE)
  r <- res[res$ID == s$go_id, , drop = FALSE]
  if (nrow(r) != 1L) stop("Exact GO result is absent or duplicated: ", s$go_id,
                         call. = FALSE)
  ranked <- stats::setNames(g$collapsed_statistic, g$GeneSymbol)
  ranked <- ranked[order(ranked, decreasing = TRUE)]
  sy <- suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db, keys = s$go_id, keytype = "GOALL",
    columns = "SYMBOL")$SYMBOL)
  members <- intersect(unique(sy[!is.na(sy)]), names(ranked))
  if (length(members) != r$setSize[[1]])
    stop("GO membership/setSize mismatch: ", s$go_id, " ", ds, " ", unit,
         call. = FALSE)
  hit <- names(ranked) %in% members
  n <- length(ranked); nh <- sum(hit)
  phit <- numeric(n); phit[hit] <- abs(ranked[hit])
  running <- cumsum(phit / sum(phit)) -
    cumsum(replace(numeric(n), !hit, 1 / (n - nh)))
  es <- if (abs(max(running)) > abs(min(running))) max(running) else min(running)
  if (!isTRUE(all.equal(es, r$enrichmentScore[[1]], tolerance = 1e-10)))
    stop("Reconstructed enrichment score differs from stored result: ",
         s$go_id, " ", ds, " ", unit, call. = FALSE)
  if (!isTRUE(all.equal(r$NES[[1]], pick$NES[[1]], tolerance = 1e-8)) ||
      !isTRUE(all.equal(r$p.adjust[[1]], pick$BH_FDR[[1]], tolerance = 1e-8)))
    stop("Frozen inventory differs from stored GSEA result: ", s$go_id,
         " ", ds, " ", unit, call. = FALSE)
  peak <- if (es < 0) which.min(running) else which.max(running)
  selection[[i]] <- data.frame(
    display_order = s$display_order, theme_id = s$theme_id, GO_ID = s$go_id,
    GO_description = pick$GO_description, dataset = ds, spatial_unit = unit,
    selection_rule = paste0("fixed exact GO ID; among FDR-supported SUS - RES ",
                            "occurrences choose smallest stored BH_FDR; ",
                            "ties by dataset then spatial_unit"),
    n_eligible_units = nrow(candidates), NES = r$NES[[1]],
    BH_FDR = r$p.adjust[[1]], enrichmentScore = es,
    setSize = r$setSize[[1]], n_ranked = n, peak_rank = peak,
    stringsAsFactors = FALSE)
  curves[[i]] <- data.frame(
    display_order = s$display_order, theme_id = s$theme_id,
    GO_ID = s$go_id, dataset = ds, spatial_unit = unit,
    rank = seq_len(n), running_ES = running, hit = hit,
    peak = seq_len(n) == peak, stringsAsFactors = FALSE)
}
selection <- do.call(rbind, selection)
curves <- do.call(rbind, curves)
regional <- all[all$GO_ID %in% spec$go_id &
                  all$theme_id %in% spec$theme_id &
                  all$contrast %in% c("RES - CON", "SUS - CON", "SUS - RES"),
                required, drop = FALSE]
if (anyDuplicated(regional[c("dataset", "spatial_unit", "contrast", "GO_ID")]))
  stop("Duplicate regional exact-term inventory rows.", call. = FALSE)
inputs$path <- substring(gsub("\\\\", "/", inputs$path), nchar(root) + 2L)
inputs$sha256 <- unname(tools::sha256sum(file.path(root, inputs$path)))
inputs$annotation_package <- as.character(utils::packageVersion("org.Mm.eg.db"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(selection, file.path(out, "selection.csv"), row.names = FALSE)
utils::write.csv(curves, file.path(out, "running_enrichment_curves.csv"),
                 row.names = FALSE)
utils::write.csv(regional, file.path(out, "regional_exact_term_inventory.csv"),
                 row.names = FALSE)
utils::write.csv(inputs, file.path(out, "input_manifest.csv"), row.names = FALSE)
cat("Wrote candidate source data to ", out, "\n", sep = "")
