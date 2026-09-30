#!/usr/bin/env Rscript
# Figure 3 g-i selection rule applied downstream to the four appendix terms.
# Reads frozen GSEA/rank inputs and stored mapped DA; fits no model.

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "pipeline.yml")))
  stop("Run from the pRoteomics repository root.", call. = FALSE)
rel <- function(...) file.path(root, ...)
source_dir <- rel("exports", "publication_source_data",
                  "supplementary_go_pathways_candidate")
selection <- utils::read.csv(file.path(source_dir, "selection.csv"),
                             stringsAsFactors = FALSE)
curves <- utils::read.csv(file.path(source_dir,
                                   "running_enrichment_curves.csv"),
                          stringsAsFactors = FALSE)
out <- rel("exports", "publication_source_data",
           "supplementary_go_protein_zooms_candidate")
if (file.exists(out)) stop("Candidate export already exists: ", out,
                           call. = FALSE)
if (nrow(selection) != 4L || anyDuplicated(selection$GO_ID))
  stop("Candidate exact-term selection is incomplete.", call. = FALSE)

rows <- list()
inputs <- data.frame(role = c("selection", "running_curves"),
                     path = file.path(source_dir,
                                      c("selection.csv",
                                        "running_enrichment_curves.csv")),
                     stringsAsFactors = FALSE)
for (i in seq_len(nrow(selection))) {
  s <- selection[i, ]
  ds <- s$dataset; unit <- s$spatial_unit
  unit_dir <- if (identical(ds, "microglia")) paste0(unit, "_microglia") else unit
  tok <- gsub("_", "", unit)
  if (identical(ds, "microglia")) tok <- paste0(tok, "microglia")
  cmp <- sprintf("%ssus_%sres", tok, tok)
  ranked_path <- rel("data", "processed", "04_differential_expression_enrichment",
                     "clusterProfiler", ds, "phenotype_within_unit", unit_dir,
                     cmp, "protein_group_audits", "collapsed_gene_input.csv")
  if (!file.exists(ranked_path)) stop("Ranked input missing: ", ranked_path,
                                      call. = FALSE)
  inputs <- rbind(inputs, data.frame(role = "ranked_gene_input",
                                    path = ranked_path))
  g <- utils::read.csv(ranked_path, stringsAsFactors = FALSE)
  if (!all(c("GeneSymbol", "collapsed_statistic") %in% names(g)) ||
      anyDuplicated(g$GeneSymbol) || anyNA(g$collapsed_statistic))
    stop("Ranked input schema/identity mismatch: ", s$GO_ID, call. = FALSE)
  ranked <- stats::setNames(g$collapsed_statistic, g$GeneSymbol)
  ranked <- ranked[order(ranked, decreasing = TRUE)]
  c <- curves[curves$GO_ID == s$GO_ID, , drop = FALSE]
  c <- c[order(c$rank), , drop = FALSE]
  if (nrow(c) != length(ranked) ||
      !identical(c$rank, seq_along(ranked)) ||
      sum(c$peak) != 1L || sum(c$hit) != s$setSize)
    stop("Frozen curve/rank mismatch: ", s$GO_ID, call. = FALSE)
  peak <- which(c$peak)
  edge_idx <- if (s$enrichmentScore < 0) peak:length(ranked) else seq_len(peak)
  leading <- names(ranked)[edge_idx][c$hit[edge_idx]]
  st <- ranked[leading]
  st <- st[order(-abs(st))]
  keep <- names(st)[seq_len(min(7L, length(st)))]
  if (!length(keep)) stop("No leading-edge genes: ", s$GO_ID, call. = FALSE)
  for (cb in list(c("res", "con"), c("sus", "con"), c("sus", "res"))) {
    file <- sprintf("%s%s_%s%s.csv", tok, cb[[1]], tok, cb[[2]])
    mapped_path <- rel("data", "processed", "02_id_mapping", "mapped", ds,
                       "forward", "per_file", file)
    if (!file.exists(mapped_path)) stop("Mapped DA file missing: ", mapped_path,
                                        call. = FALSE)
    inputs <- rbind(inputs, data.frame(role = "mapped_protein_da",
                                      path = mapped_path))
    da <- utils::read.csv(mapped_path, stringsAsFactors = FALSE)
    sym <- intersect(c("official_gene_symbol", "gene_symbol"), names(da))[[1]]
    if (!all(c("ProteinGroupID", "log2fc", "padj") %in% names(da)))
      stop("Mapped DA schema mismatch: ", mapped_path, call. = FALSE)
    m <- match(keep, da[[sym]])
    z <- data.frame(
      display_order = s$display_order, theme_id = s$theme_id,
      GO_ID = s$GO_ID, dataset = ds, spatial_unit = unit,
      gene = keep, ProteinGroupID = da$ProteinGroupID[m],
      contrast = paste(toupper(cb), collapse = " - "),
      log2FC = da$log2fc[m], BH_FDR = da$padj[m],
      rank_statistic = unname(st[keep]),
      stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] <- z[!is.na(z$log2FC), , drop = FALSE]
  }
}
values <- do.call(rbind, rows)
if (anyDuplicated(values[c("GO_ID", "gene", "contrast")]) ||
    anyNA(values$ProteinGroupID) || anyNA(values$BH_FDR))
  stop("Candidate protein values have duplicate or missing provenance.",
       call. = FALSE)
inputs$path <- substring(gsub("\\\\", "/", inputs$path), nchar(root) + 2L)
inputs$sha256 <- unname(tools::sha256sum(file.path(root, inputs$path)))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(values, file.path(out, "protein_zoom_values.csv"),
                 row.names = FALSE)
utils::write.csv(inputs, file.path(out, "protein_input_manifest.csv"),
                 row.names = FALSE)
cat("Wrote candidate protein source data to ", out, "\n", sep = "")
