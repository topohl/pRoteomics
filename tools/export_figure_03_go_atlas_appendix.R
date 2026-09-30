#!/usr/bin/env Rscript
# Frozen, downstream-only source for every FDR-supported SUS-RES GO term in
# the seven Figure 3 atlas themes. No enrichment or protein model is fitted.
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(root, "pipeline.yml")))
  stop("Run from the pRoteomics repository root.", call. = FALSE)
for (pkg in c("AnnotationDbi", "org.Mm.eg.db"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Required annotation package unavailable: ", pkg, call. = FALSE)
rel <- function(...) file.path(root, ...)
src <- rel("exports", "publication_source_data",
           "supplementary_selection_inventories")
full_path <- file.path(src, "pathway_enrichment_inventory.csv")
sup_path <- file.path(src, "pathway_enrichment_inventory_fdr_supported.csv")
out <- rel("exports", "publication_source_data", "figure_03_go_atlas_appendix")
if (file.exists(out)) stop("Appendix export already exists: ", out, call. = FALSE)
all <- read.csv(full_path, stringsAsFactors = FALSE, check.names = FALSE)
sup <- read.csv(sup_path, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("dataset", "spatial_unit", "contrast", "GO_ID",
          "GO_description", "NES", "BH_FDR", "theme_id",
          "theme_claim_eligible")
if (!all(need %in% names(all)) || !all(need %in% names(sup)))
  stop("Frozen inventory schema mismatch.", call. = FALSE)
themes <- c("rna_processing_splicing_rnp", "ribosome_translation",
            "chromatin_organization", "mitochondrial_respiration_oxphos",
            "synaptic_signaling_vesicle", "neuron_projection_development",
            "autophagy_lysosome_endosome")
theme_labels <- c("RNA processing", "Translation / ribosome",
                  "Chromatin / epigenetic regulation",
                  "Mitochondrial respiration", "Synaptic signalling / vesicle",
                  "Neuron projection development",
                  "Autophagy / endolysosomal")
supported <- sup[sup$contrast == "SUS - RES" &
                   sup$theme_claim_eligible %in% TRUE, , drop = FALSE]
membership <- do.call(rbind, lapply(seq_len(nrow(supported)), function(i) {
  ids <- strsplit(supported$theme_id[[i]], ";", fixed = TRUE)[[1]]
  ids <- ids[ids %in% themes]
  data.frame(theme_id = ids, GO_ID = supported$GO_ID[[i]],
             stringsAsFactors = FALSE)
}))
membership <- unique(membership)
if (!nrow(membership) || anyNA(membership) ||
    !setequal(membership$GO_ID, supported$GO_ID))
  stop("Supported term/theme membership is incomplete.", call. = FALSE)
membership <- membership[order(match(membership$theme_id, themes),
                               membership$GO_ID), , drop = FALSE]
term_ids <- unique(membership$GO_ID)
inputs <- data.frame(role = c("frozen_full_inventory", "frozen_supported_inventory"),
                     path = c(full_path, sup_path), stringsAsFactors = FALSE)
selection <- vector("list", length(term_ids))
curves <- vector("list", length(term_ids))
proteins <- vector("list", length(term_ids))
cache <- new.env(parent = emptyenv())
get_context <- function(ds, unit) {
  key <- paste(ds, unit, sep = "|")
  if (exists(key, envir = cache, inherits = FALSE))
    return(get(key, envir = cache))
  unit_dir <- if (ds == "microglia") paste0(unit, "_microglia") else unit
  tok <- gsub("_", "", unit)
  if (ds == "microglia") tok <- paste0(tok, "microglia")
  cmp <- sprintf("%ssus_%sres", tok, tok)
  base <- rel("data", "processed", "04_differential_expression_enrichment",
              "clusterProfiler", ds, "phenotype_within_unit", unit_dir, cmp)
  ranked_path <- file.path(base, "protein_group_audits", "collapsed_gene_input.csv")
  gsea_path <- file.path(base, "GO", "BP", "GSEA_BP_results_full.csv")
  if (!file.exists(ranked_path) || !file.exists(gsea_path))
    stop("Canonical GSEA files missing: ", key, call. = FALSE)
  g <- read.csv(ranked_path, stringsAsFactors = FALSE)
  res <- read.csv(gsea_path, stringsAsFactors = FALSE)
  if (!all(c("GeneSymbol", "collapsed_statistic") %in% names(g)) ||
      !all(c("ID", "setSize", "enrichmentScore", "NES", "p.adjust") %in% names(res)) ||
      anyDuplicated(g$GeneSymbol) || anyNA(g$GeneSymbol) ||
      anyNA(g$collapsed_statistic) || any(!is.finite(g$collapsed_statistic)))
    stop("Canonical GSEA schema/rank mismatch: ", key, call. = FALSE)
  ranked <- setNames(g$collapsed_statistic, g$GeneSymbol)
  ranked <- ranked[order(ranked, decreasing = TRUE)]
  inputs <<- rbind(inputs, data.frame(role = c("ranked_gene_input",
                                             "stored_gsea_result"),
                                      path = c(ranked_path, gsea_path)))
  value <- list(ranked = ranked, gsea = res, token = tok)
  assign(key, value, envir = cache)
  value
}
for (i in seq_along(term_ids)) {
  id <- term_ids[[i]]
  candidates <- supported[supported$GO_ID == id, , drop = FALSE]
  candidates <- candidates[order(candidates$BH_FDR, candidates$dataset,
                                 candidates$spatial_unit), , drop = FALSE]
  pick <- candidates[1, , drop = FALSE]
  ds <- pick$dataset[[1]]; unit <- pick$spatial_unit[[1]]
  ctx <- get_context(ds, unit)
  ranked <- ctx$ranked
  r <- ctx$gsea[ctx$gsea$ID == id, , drop = FALSE]
  if (nrow(r) != 1L ||
      !isTRUE(all.equal(r$NES[[1]], pick$NES[[1]], tolerance = 1e-8)) ||
      !isTRUE(all.equal(r$p.adjust[[1]], pick$BH_FDR[[1]], tolerance = 1e-8)))
    stop("Stored GSEA/frozen inventory mismatch: ", id, call. = FALSE)
  sy <- suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db, keys = id, keytype = "GOALL",
    columns = "SYMBOL")$SYMBOL)
  members <- intersect(unique(sy[!is.na(sy)]), names(ranked))
  if (length(members) != r$setSize[[1]])
    stop("GO membership/setSize mismatch: ", id, " ", ds, " ", unit,
         call. = FALSE)
  hit <- names(ranked) %in% members
  n <- length(ranked); nh <- sum(hit)
  phit <- numeric(n); phit[hit] <- abs(ranked[hit])
  running <- cumsum(phit / sum(phit)) -
    cumsum(replace(numeric(n), !hit, 1 / (n - nh)))
  es <- if (abs(max(running)) > abs(min(running))) max(running) else min(running)
  if (!isTRUE(all.equal(es, r$enrichmentScore[[1]], tolerance = 1e-10)))
    stop("Reconstructed ES differs from stored result: ", id, call. = FALSE)
  peak <- if (es < 0) which.min(running) else which.max(running)
  selection[[i]] <- data.frame(
    GO_ID = id, GO_description = pick$GO_description, dataset = ds,
    spatial_unit = unit, NES = r$NES[[1]], BH_FDR = r$p.adjust[[1]],
    enrichmentScore = es, setSize = r$setSize[[1]], n_ranked = n,
    peak_rank = peak, n_supported_units = nrow(candidates),
    selection_rule = paste0("Smallest stored BH_FDR among FDR-supported ",
                            "SUS - RES occurrences; ties by dataset, spatial_unit"),
    stringsAsFactors = FALSE)
  curves[[i]] <- data.frame(GO_ID = id, rank = seq_len(n),
                            running_ES = running, hit = hit,
                            peak = seq_len(n) == peak)
  edge_idx <- if (es < 0) peak:n else seq_len(peak)
  leading <- names(ranked)[edge_idx][hit[edge_idx]]
  st <- ranked[leading]
  st <- st[order(-abs(st))]
  keep <- names(st)[seq_len(min(7L, length(st)))]
  if (!length(keep)) stop("No leading-edge genes: ", id, call. = FALSE)
  zlist <- list()
  for (cb in list(c("res", "con"), c("sus", "con"), c("sus", "res"))) {
    file <- sprintf("%s%s_%s%s.csv", ctx$token, cb[[1]], ctx$token, cb[[2]])
    mapped_path <- rel("data", "processed", "02_id_mapping", "mapped", ds,
                       "forward", "per_file", file)
    if (!file.exists(mapped_path))
      stop("Mapped DA file missing: ", mapped_path, call. = FALSE)
    da <- read.csv(mapped_path, stringsAsFactors = FALSE)
    sym <- intersect(c("official_gene_symbol", "gene_symbol"), names(da))
    if (!length(sym) ||
        !all(c("ProteinGroupID", "log2fc", "padj") %in% names(da)))
      stop("Mapped DA schema mismatch: ", mapped_path, call. = FALSE)
    m <- match(keep, da[[sym[[1]]]])
    z <- data.frame(GO_ID = id, dataset = ds, spatial_unit = unit,
                    gene = keep, ProteinGroupID = da$ProteinGroupID[m],
                    contrast = paste(toupper(cb), collapse = " - "),
                    log2FC = da$log2fc[m], BH_FDR = da$padj[m],
                    rank_statistic = unname(st[keep]),
                    stringsAsFactors = FALSE)
    if (anyNA(z[c("ProteinGroupID", "log2FC", "BH_FDR")]))
      stop("Mapped leading-edge protein absent: ", id, " ", file,
           call. = FALSE)
    zlist[[length(zlist) + 1L]] <- z
    inputs <- rbind(inputs, data.frame(role = "mapped_protein_da",
                                       path = mapped_path))
  }
  proteins[[i]] <- do.call(rbind, zlist)
  if (i %% 10L == 0L) message("Verified ", i, "/", length(term_ids), " GO terms")
}
selection <- do.call(rbind, selection)
curves <- do.call(rbind, curves)
proteins <- do.call(rbind, proteins)
membership$theme_order <- match(membership$theme_id, themes)
membership$theme_label <- theme_labels[membership$theme_order]
membership$display_order <- seq_len(nrow(membership))
membership <- merge(membership, selection, by = "GO_ID", sort = FALSE)
membership <- membership[order(membership$display_order), , drop = FALSE]
regional <- all[all$GO_ID %in% term_ids &
                  all$contrast %in% c("RES - CON", "SUS - CON", "SUS - RES"),
                c("dataset", "spatial_unit", "contrast", "GO_ID",
                  "GO_description", "NES", "BH_FDR"), drop = FALSE]
if (anyDuplicated(regional[c("dataset", "spatial_unit", "contrast", "GO_ID")]) ||
    anyDuplicated(proteins[c("GO_ID", "gene", "contrast")]) ||
    !setequal(unique(curves$GO_ID), term_ids) ||
    !setequal(unique(proteins$GO_ID), term_ids))
  stop("Appendix source has duplicate or missing term identities.", call. = FALSE)
inputs <- unique(inputs)
inputs$path <- substring(gsub("\\\\", "/", inputs$path), nchar(root) + 2L)
inputs$sha256 <- unname(tools::sha256sum(file.path(root, inputs$path)))
inputs$annotation_package <- as.character(packageVersion("org.Mm.eg.db"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(membership, file.path(out, "theme_term_index.csv"), row.names = FALSE)
write.csv(selection, file.path(out, "selected_contexts.csv"), row.names = FALSE)
write.csv(curves, file.path(out, "running_enrichment_curves.csv"), row.names = FALSE)
write.csv(regional, file.path(out, "regional_exact_term_inventory.csv"),
          row.names = FALSE)
write.csv(proteins, file.path(out, "protein_zoom_values.csv"), row.names = FALSE)
write.csv(inputs, file.path(out, "input_manifest.csv"), row.names = FALSE)
cat("Exported ", length(term_ids), " unique terms in ", nrow(membership),
    " theme memberships to ", out, "\n", sep = "")
