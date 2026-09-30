#!/usr/bin/env Rscript
# Export the exact three Figure 3 exemplar contexts for downstream rendering.
# Existing GSEA results, ranked inputs and mapped protein DA values are read;
# no enrichment or differential model is fitted.

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "pipeline.yml")))
  stop("Run from the pRoteomics repository root.", call. = FALSE)
for (pkg in c("AnnotationDbi", "org.Mm.eg.db"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Required annotation package unavailable: ", pkg, call. = FALSE)
rel <- function(...) file.path(root, ...)
win_long <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows" && nchar(path) >= 240L)
    paste0("\\\\?\\", gsub("/", "\\\\", path)) else path
}
out <- rel("exports", "publication_source_data", "figure_03_adaptation")
if (file.exists(out)) stop("Adaptation export already exists: ", out,
                           call. = FALSE)

selection_spec <- data.frame(
  exemplar = 1:3,
  programme_id = c("synapse_vesicle", "rna_rnp", "mitochondria_oxphos"),
  dataset = c("neuron_neuropil", "neuron_soma", "microglia"),
  spatial_unit = c("CA3_sr", "CA2_sp", "CA1"),
  GO_ID = c("GO:0099536", "GO:0006397", "GO:0006119"),
  stringsAsFactors = FALSE)

inventory_path <- rel("exports", "publication_source_data",
                      "supplementary_selection_inventories",
                      "pathway_enrichment_inventory.csv")
inventory <- read.csv(inventory_path, stringsAsFactors = FALSE,
                      check.names = FALSE)
inputs <- data.frame(role = "frozen_pathway_inventory", path = inventory_path,
                     stringsAsFactors = FALSE)
selection <- list(); curves <- list(); proteins <- list()

for (i in seq_len(nrow(selection_spec))) {
  sp <- selection_spec[i, , drop = FALSE]
  ds <- sp$dataset[[1]]; unit <- sp$spatial_unit[[1]]; go <- sp$GO_ID[[1]]
  unit_dir <- if (ds == "microglia") paste0(unit, "_microglia") else unit
  token <- gsub("_", "", unit)
  if (ds == "microglia") token <- paste0(token, "microglia")
  contrast_dir <- sprintf("%ssus_%sres", token, token)
  base <- rel("data", "processed", "04_differential_expression_enrichment",
              "clusterProfiler", ds, "phenotype_within_unit", unit_dir,
              contrast_dir)
  ranked_path <- file.path(base, "protein_group_audits",
                           "collapsed_gene_input.csv")
  gsea_path <- file.path(base, "GO", "BP", "GSEA_BP_results_full.csv")
  if (!file.exists(win_long(ranked_path)) || !file.exists(win_long(gsea_path)))
    stop("Canonical GSEA inputs missing for exemplar ", i, call. = FALSE)
  ranked_df <- read.csv(win_long(ranked_path), stringsAsFactors = FALSE)
  gsea <- read.csv(win_long(gsea_path), stringsAsFactors = FALSE)
  ranked <- setNames(ranked_df$collapsed_statistic, ranked_df$GeneSymbol)
  ranked <- ranked[order(ranked, decreasing = TRUE)]
  stored <- gsea[gsea$ID == go, , drop = FALSE]
  inv <- inventory[inventory$dataset == ds &
                     inventory$spatial_unit == unit &
                     inventory$contrast == "SUS - RES" &
                     inventory$GO_ID == go, , drop = FALSE]
  if (nrow(stored) != 1L || nrow(inv) != 1L ||
      !isTRUE(all.equal(stored$NES[[1]], inv$NES[[1]], tolerance = 1e-8)) ||
      !isTRUE(all.equal(stored$p.adjust[[1]], inv$BH_FDR[[1]],
                        tolerance = 1e-8)))
    stop("Stored GSEA/inventory mismatch for exemplar ", i, call. = FALSE)

  symbols <- suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db, keys = go, keytype = "GOALL",
    columns = "SYMBOL")$SYMBOL)
  members <- intersect(unique(symbols[!is.na(symbols)]), names(ranked))
  if (length(members) != stored$setSize[[1]])
    stop("GO membership/setSize mismatch for exemplar ", i, call. = FALSE)
  hit <- names(ranked) %in% members
  n <- length(ranked); nh <- sum(hit)
  phit <- numeric(n); phit[hit] <- abs(ranked[hit])
  running <- cumsum(phit / sum(phit)) -
    cumsum(replace(numeric(n), !hit, 1 / (n - nh)))
  es <- if (abs(max(running)) > abs(min(running))) max(running) else min(running)
  if (!isTRUE(all.equal(es, stored$enrichmentScore[[1]], tolerance = 1e-10)))
    stop("Reconstructed ES mismatch for exemplar ", i, call. = FALSE)
  peak <- if (es < 0) which.min(running) else which.max(running)
  selection[[i]] <- data.frame(
    exemplar = sp$exemplar, programme_id = sp$programme_id,
    dataset = ds, spatial_unit = unit, GO_ID = go,
    GO_description = inv$GO_description, NES = stored$NES,
    BH_FDR = stored$p.adjust, enrichmentScore = es,
    setSize = stored$setSize, n_ranked = n, peak_rank = peak,
    selection_rule = "manuscript-fixed programme/location/constituent-term exemplar",
    stringsAsFactors = FALSE)
  curves[[i]] <- data.frame(
    exemplar = sp$exemplar, GO_ID = go, dataset = ds, spatial_unit = unit,
    rank = seq_len(n), running_ES = running, hit = hit,
    peak = seq_len(n) == peak, stringsAsFactors = FALSE)

  edge_idx <- if (es < 0) peak:n else seq_len(peak)
  leading <- names(ranked)[edge_idx][hit[edge_idx]]
  stats <- ranked[leading]
  stats <- stats[order(-abs(stats))]
  keep <- names(stats)[seq_len(min(7L, length(stats)))]
  p_rows <- list()
  for (pair in list(c("res", "con"), c("sus", "con"), c("sus", "res"))) {
    filename <- sprintf("%s%s_%s%s.csv", token, pair[[1]], token, pair[[2]])
    mapped_path <- rel("data", "processed", "02_id_mapping", "mapped", ds,
                       "forward", "per_file", filename)
    da <- read.csv(mapped_path, stringsAsFactors = FALSE)
    symbol_col <- intersect(c("official_gene_symbol", "gene_symbol"),
                            names(da))[[1]]
    m <- match(keep, da[[symbol_col]])
    z <- data.frame(
      exemplar = sp$exemplar, GO_ID = go, dataset = ds, spatial_unit = unit,
      gene = keep, ProteinGroupID = da$ProteinGroupID[m],
      contrast = paste(toupper(pair), collapse = " - "),
      log2FC = da$log2fc[m], BH_FDR = da$padj[m],
      rank_statistic = unname(stats[keep]), stringsAsFactors = FALSE)
    if (anyNA(z[c("ProteinGroupID", "log2FC", "BH_FDR")]))
      stop("Mapped leading-edge protein missing for exemplar ", i,
           call. = FALSE)
    p_rows[[length(p_rows) + 1L]] <- z
    inputs <- rbind(inputs, data.frame(role = "mapped_protein_da",
                                       path = mapped_path))
  }
  proteins[[i]] <- do.call(rbind, p_rows)
  inputs <- rbind(inputs,
    data.frame(role = c("ranked_gene_input", "stored_gsea_result"),
               path = c(ranked_path, gsea_path)))
}

selection <- do.call(rbind, selection)
curves <- do.call(rbind, curves)
proteins <- do.call(rbind, proteins)
if (nrow(selection) != 3L || anyDuplicated(selection$exemplar) ||
    anyDuplicated(curves[c("exemplar", "rank")]) ||
    anyDuplicated(proteins[c("exemplar", "gene", "contrast")]))
  stop("Adaptation source identities are duplicated or incomplete.", call. = FALSE)

inputs <- unique(inputs)
inputs$path <- substring(gsub("\\\\", "/", inputs$path), nchar(root) + 2L)
hash_paths <- vapply(file.path(root, inputs$path), win_long, character(1))
inputs$sha256 <- unname(tools::sha256sum(hash_paths))
inputs$annotation_package <- as.character(packageVersion("org.Mm.eg.db"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(selection, file.path(out, "selection.csv"), row.names = FALSE)
write.csv(curves, file.path(out, "running_enrichment_curves.csv"),
          row.names = FALSE)
write.csv(proteins, file.path(out, "protein_zoom_values.csv"), row.names = FALSE)
write.csv(inputs, file.path(out, "input_manifest.csv"), row.names = FALSE)
cat("Exported exact Figure 3 adaptation exemplars to ", out, "\n", sep = "")
