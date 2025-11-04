# ================================ parallel-enabled
# WGCNA with spatial traits + preservation + condition×region×layer×celltype panels
# Outputs organized into subfolders under output_dir
# ================================ parallel-enabled

# Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db", ask = FALSE, update = FALSE)
suppressPackageStartupMessages(
  pacman::p_load(WGCNA, flashClust, curl, readxl, ggplot2, svglite, GO.db,
                 reshape2, gtools, patchwork, cowplot, pheatmap, dplyr, tidyr,
                 httr, jsonlite, purrr, AnnotationDbi, org.Mm.eg.db, readr, stringr, install = TRUE)
)

# Parallel setup
nCores <- tryCatch({
  pc <- parallel::detectCores(logical = FALSE)
  if (is.na(pc) || pc < 2) 2 else pc
}, error = function(e) 2)
enableWGCNAThreads(nThreads = nCores)
WGCNAnThreads()

# --------------------------
# Paths and data load
# --------------------------
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Subfolders
subdirs <- list(
  plots_qc          = file.path(output_dir, "plots_qc"),
  plots_traits      = file.path(output_dir, "plots_traits"),
  network           = file.path(output_dir, "network"),
  tables_modules    = file.path(output_dir, "tables_modules"),
  tables_pres       = file.path(output_dir, "tables_preservation")
)
invisible(lapply(subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Path helpers
fp_qc      <- function(...) file.path(subdirs$plots_qc, ...)
fp_traits  <- function(...) file.path(subdirs$plots_traits, ...)
fp_net     <- function(...) file.path(subdirs$network, ...)
fp_modtab  <- function(...) file.path(subdirs$tables_modules, ...)
fp_prestab <- function(...) file.path(subdirs$tables_pres, ...)

# Optional: safe svg helper
save_svg <- function(path, width, height, expr) {
  svglite::svglite(file = path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

expr_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/male.data.xlsx"
meta_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info.xlsx"

# ================================
# Mouse-only mapping: robust idmapping parser + offline + SYMBOL/ALIAS + Entrez + UniProt gene_primary + QC
# ================================

suppressPackageStartupMessages({
  if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
  pacman::p_load(
    WGCNA, readxl, ggplot2, svglite, dplyr, tidyr, tibble, stringr, readr, purrr,
    AnnotationDbi, org.Mm.eg.db, install = TRUE
  )
})
if (!requireNamespace("UniProt.ws", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("UniProt.ws", ask = FALSE, update = FALSE)
}
library(UniProt.ws)

# --------------------------
# Paths and environment
# --------------------------
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output"
subdirs <- list(
  plots_qc         = file.path(output_dir, "plots_qc"),
  plots_traits     = file.path(output_dir, "plots_traits"),
  network          = file.path(output_dir, "network"),
  tables_modules   = file.path(output_dir, "tables_modules"),
  tables_pres      = file.path(output_dir, "tables_preservation")
)
safe_dir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE); if (file.access(path, 2) != 0) stop(sprintf("Not writable: %s", path)); invisible(normalizePath(path)) }
log_session <- function(out_dir) { safe_dir(out_dir); writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "session_info.txt")) }
invisible(lapply(c(output_dir, unlist(subdirs)), safe_dir)); log_session(output_dir)

# --------------------------
# Inputs
# --------------------------
expr_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/male.data.xlsx"
meta_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info.xlsx"
idmap_dat <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/Datasets/MOUSE_10090_idmapping.dat"

stop_if_missing <- function(path) if (!file.exists(path)) stop(sprintf("Missing file: %s", path))
read_head <- function(path) { df <- readxl::read_excel(path); utils::write.table(utils::head(df, 10), file.path(subdirs$plots_qc, paste0(basename(path), "_head10.tsv")), sep="\t", row.names=FALSE, quote=FALSE); df }

male.data <- { stop_if_missing(expr_xlsx); read_head(expr_xlsx) } %>% dplyr::mutate(.row_id = dplyr::row_number())
meta.data <- { stop_if_missing(meta_xlsx); read_head(meta_xlsx) }

# Keep only tokens with _MOUSE, and drop rows with no such tokens
male.data.filtered <- male.data %>%
  tidyr::separate_rows(gene_symbol, sep = ";") %>%
  dplyr::mutate(gene_symbol = toupper(trimws(gene_symbol))) %>%
  dplyr::filter(grepl("_MOUSE$", gene_symbol)) %>%
  dplyr::group_by(.row_id) %>%
  dplyr::summarise(gene_symbol = paste(unique(gene_symbol), collapse = ";"), .groups = "drop") %>%
  dplyr::right_join(male.data %>% dplyr::select(-gene_symbol, .row_id), by = ".row_id") %>%
  dplyr::filter(!is.na(gene_symbol) & nzchar(gene_symbol))  # drop rows with no _MOUSE tokens

# Replace male.data with the filtered version for all downstream steps
male.data <- male.data.filtered

# --------------------------
# Robust UniProtKB-ID -> Accession map (offline)
# --------------------------
stop_if_missing(idmap_dat)
idmap_tbl <- readr::read_tsv(idmap_dat, col_names = c("ACC","DB","VAL"), col_types = "ccc", progress = FALSE, quote = "", comment = "")
idmap_uid <- idmap_tbl %>%
  dplyr::filter(DB == "UniProtKB-ID" & grepl("_MOUSE\\s*$", VAL) & nzchar(ACC)) %>%
  dplyr::transmute(
    UNIPROT    = toupper(trimws(ACC)),
    entry_full = toupper(trimws(VAL)),
    entry_base = toupper(gsub("_MOUSE$", "", trimws(VAL)))
  )
entry_map <- idmap_uid %>% dplyr::distinct(entry_base, .keep_all = TRUE)
if (!nrow(entry_map)) stop("entry_map is empty after robust parse")

# Sentinel sanity check
sentinels <- c("AIF1","AKAP2","ADCY1","AKAP1","AMPD3","ANXA3","1433S","ACK1","AIP","ADA10")
missing_sentinels <- setdiff(sentinels, entry_map$entry_base)
if (length(missing_sentinels)) {
  warning(sprintf("entry_map missing expected keys: %s", paste(missing_sentinels, collapse=", ")))
  readr::write_tsv(entry_map, file.path(subdirs$tables_modules, "entry_map_debug.tsv"))
}

# --------------------------
# Mouse-only tokenization and classification
# --------------------------
normalize_token <- function(x) { x <- toupper(gsub("\\s+", "", x)); x <- gsub("\\u00A0", "", x); x <- gsub("\\.+", ".", x); x <- gsub("__+", "_", x); x }
to_base_no_iso_mouse <- function(x) { x <- gsub("-\\d+$", "", x); gsub("_MOUSE$", "", x) }

tokenize_mouse_only <- function(male_df) {
  tok <- male_df %>% tidyr::separate_rows(gene_symbol, sep = ";") %>% dplyr::mutate(token_raw = gene_symbol, token_up = normalize_token(gene_symbol))
  dropped_non_mouse <- tok %>% dplyr::filter(!grepl("_MOUSE$", token_up))
  if (nrow(dropped_non_mouse)) readr::write_tsv(dropped_non_mouse, file.path(subdirs$tables_modules, "dropped_non_mouse_tokens.tsv"))
  tok %>%
    dplyr::filter(grepl("_MOUSE$", token_up)) %>%
    dplyr::mutate(
      token_base = to_base_no_iso_mouse(token_up),
      looks_ac   = grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", token_base),
      looks_entry= grepl("^[A-Z0-9][A-Z0-9\\-\\.]+$", token_base),
      id_class = dplyr::case_when(
        looks_ac    ~ "UNIPROT_AC_MOUSE",
        looks_entry ~ "UNIPROT_ENTRY",
        TRUE        ~ "UNKNOWN"
      ),
      Resolved_UNIPROT = NA_character_,
      strategy = NA_character_
    )
}
resolved2 <- tokenize_mouse_only(male.data)

# --------------------------
# Mapping stack
# --------------------------

# 1) Accept accession-like bases
idx_ac <- which(resolved2$id_class == "UNIPROT_AC_MOUSE")
if (length(idx_ac)) {
  resolved2$Resolved_UNIPROT[idx_ac] <- resolved2$token_base[idx_ac]
  resolved2$strategy[idx_ac] <- "accept_accession_base"
}

# 2) Offline entry map (now robust)
idx_en <- which(resolved2$id_class == "UNIPROT_ENTRY" & (is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT)))
if (length(idx_en)) {
  hit <- entry_map$UNIPROT[match(toupper(resolved2$token_base[idx_en]), entry_map$entry_base)]
  ok <- !is.na(hit) & nzchar(hit)
  if (any(ok)) { ii <- idx_en[ok]; resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "entry_local_mouse" }
}

# 3) SYMBOL/ALIAS offline resolver (MGI-first)
need_idx <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
need_ids <- toupper(unique(resolved2$token_base[need_idx]))
is_acc <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", need_ids)
ids_ent <- unique(need_ids[!is_acc])

if (length(ids_ent)) {
  sel_sym <- try(AnnotationDbi::select(org.Mm.eg.db, keys = ids_ent, keytype = "SYMBOL", columns = c("MGIID","ENTREZID","UNIPROT","SYMBOL")), silent = TRUE)
  map_sym <- tibble::tibble()
  if (!inherits(sel_sym, "try-error") && nrow(sel_sym)) {
    map_sym <- tibble::as_tibble(sel_sym) %>%
      dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
      dplyr::group_by(SYMBOL) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup() %>%
      dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT))
  }
  kt <- try(AnnotationDbi::keytypes(org.Mm.eg.db), silent = TRUE)
  map_alias <- tibble::tibble()
  if (!inherits(kt, "try-error") && "ALIAS" %in% kt) {
    sel_alias <- try(AnnotationDbi::select(org.Mm.eg.db, keys = ids_ent, keytype = "ALIAS", columns = c("UNIPROT","ALIAS")), silent = TRUE)
    if (!inherits(sel_alias, "try-error") && nrow(sel_alias)) {
      map_alias <- tibble::as_tibble(sel_alias) %>%
        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
        dplyr::group_by(ALIAS) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup() %>%
        dplyr::transmute(input = toupper(ALIAS), primaryAccession = toupper(UNIPROT))
    }
  }
  map_symall <- dplyr::bind_rows(map_sym, map_alias) %>% dplyr::distinct(input, .keep_all = TRUE)
  if (nrow(map_symall)) {
    need_idx_now <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
    base_need <- toupper(resolved2$token_base[need_idx_now])
    hit <- map_symall$primaryAccession[match(base_need, map_symall$input)]
    ok <- !is.na(hit) & nzchar(hit)
    ii <- need_idx_now[ok]
    if (length(ii)) { resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "orgdb_mgi_symbol_first" }
  }
}

# 4) SYMBOL -> Entrez -> UniProt (offline two-hop)
need_idx <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
need_ids <- toupper(unique(resolved2$token_base[need_idx]))
sym_left <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])

if (length(sym_left)) {
  sym2eg <- try(AnnotationDbi::select(org.Mm.eg.db, keys = sym_left, keytype = "SYMBOL", columns = c("ENTREZID","SYMBOL")), silent = TRUE)
  eg2up  <- tibble::tibble()
  if (!inherits(sym2eg, "try-error") && nrow(sym2eg)) {
    ekeys <- unique(na.omit(sym2eg$ENTREZID))
    if (length(ekeys)) {
      egsel <- try(AnnotationDbi::select(org.Mm.eg.db, keys = ekeys, keytype = "ENTREZID", columns = c("UNIPROT","ENTREZID")), silent = TRUE)
      if (!inherits(egsel, "try-error") && nrow(egsel)) {
        eg2up <- tibble::as_tibble(egsel) %>%
          dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
          dplyr::group_by(ENTREZID) %>% dplyr::arrange(UNIPROT, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup()
      }
    }
    if (nrow(eg2up)) {
      map_sym2up <- tibble::as_tibble(sym2eg) %>%
        dplyr::distinct(SYMBOL, ENTREZID) %>%
        dplyr::left_join(eg2up, by = "ENTREZID") %>%
        dplyr::filter(!is.na(UNIPROT) & nzchar(UNIPROT)) %>%
        dplyr::transmute(input = toupper(SYMBOL), primaryAccession = toupper(UNIPROT)) %>%
        dplyr::distinct(input, .keep_all = TRUE)
      need_idx2 <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
      base_need <- toupper(resolved2$token_base[need_idx2])
      hit <- map_sym2up$primaryAccession[match(base_need, map_sym2up$input)]
      ok <- !is.na(hit) & nzchar(hit)
      ii <- need_idx2[ok]
      if (length(ii)) { resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "orgdb_symbol_entrez_uniprot" }
    }
  }
}

# 5) UniProt gene_primary resolver (Mus musculus), batched with retry, prefer reviewed
need_idx <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
need_ids <- unique(toupper(resolved2$token_base[need_idx]))
sym_left2 <- unique(need_ids[grepl("^[A-Z0-9\\-]{2,}$", need_ids)])

if (length(sym_left2)) {
  batch_vec <- split(sym_left2, ceiling(seq_along(sym_left2)/50))
  picks <- list()
  for (b in batch_vec) {
    q_list <- lapply(b, function(g) list(organism_id = 10090, gene_primary = g))
    query_once <- function(ql) try(UniProt.ws::queryUniProt(query = ql, fields = c("accession","id","gene_primary","reviewed"), collapse = "OR", n = 10, pageSize = 10), silent = TRUE)
    res_list <- lapply(q_list, function(ql) { out <- query_once(ql); if (inherits(out, "try-error") || !is.data.frame(out)) { Sys.sleep(0.8); out <- query_once(ql) }; out })
    ok <- res_list[!vapply(res_list, inherits, logical(1), "try-error")]
    if (length(ok)) {
      tbl <- dplyr::bind_rows(lapply(ok, tibble::as_tibble))
      if (nrow(tbl)) {
        tbl <- tbl %>% dplyr::mutate(gene_primary = toupper(.data$gene_primary), accession = toupper(.data$accession))
        pick <- tbl %>% dplyr::group_by(gene_primary) %>% dplyr::arrange(dplyr::desc(.data$reviewed), accession, .by_group = TRUE) %>% dplyr::slice_head(n=1) %>% dplyr::ungroup() %>% dplyr::transmute(input = gene_primary, primaryAccession = accession)
        picks[[length(picks)+1]] <- pick
      }
    }
  }
  if (length(picks)) {
    map_gene <- dplyr::bind_rows(picks) %>% dplyr::distinct(input, .keep_all = TRUE)
    need_idx3 <- which(is.na(resolved2$Resolved_UNIPROT) | !nzchar(resolved2$Resolved_UNIPROT))
    base_need <- toupper(resolved2$token_base[need_idx3])
    hit <- map_gene$primaryAccession[match(base_need, map_gene$input)]
    ok <- !is.na(hit) & nzchar(hit)
    ii <- need_idx3[ok]
    if (length(ii)) { resolved2$Resolved_UNIPROT[ii] <- hit[ok]; resolved2$strategy[ii] <- "uniprot_gene_primary_retry" }
  }
}

# --------------------------
# Save unmapped lists and audits
# --------------------------
audit_mouse <- resolved2 %>%
  dplyr::mutate(mapped = !is.na(Resolved_UNIPROT) & nzchar(Resolved_UNIPROT)) %>%
  dplyr::count(id_class, strategy, mapped, name = "n") %>%
  dplyr::arrange(desc(n))
readr::write_tsv(audit_mouse, file.path(subdirs$tables_modules, "mapping_audit_mouse_only_robust.tsv"))

unmapped_tokens <- resolved2 %>%
  dplyr::filter(is.na(Resolved_UNIPROT) | !nzchar(Resolved_UNIPROT)) %>%
  dplyr::transmute(
    .row_id, token_raw, token_up, token_base, id_class,
    reason = dplyr::case_when(
      grepl("[^A-Z0-9_\\-\\.]", token_base) ~ "illegal_chars",
      grepl("^[A-Z0-9\\-\\.]+$", token_base) ~ "entry_name_not_in_local_or_query",
      TRUE ~ "unexpected_format_or_na"
    )
  ) %>% dplyr::arrange(id_class, token_base)
readr::write_tsv(unmapped_tokens, file.path(subdirs$tables_modules, "unmapped_mouse_tokens.tsv"))

unmapped_summary <- unmapped_tokens %>% dplyr::count(id_class, reason, token_base, name = "n") %>% dplyr::arrange(dplyr::desc(n))
readr::write_tsv(unmapped_summary, file.path(subdirs$tables_modules, "unmapped_mouse_tokens_summary.tsv"))

# --------------------------
# Collapse to features and build expression matrix
# --------------------------
collapse_ids <- function(x) { x <- unique(x[!is.na(x) & nzchar(x)]); if (!length(x)) return(NA_character_); paste(x, collapse = ";") }
male.norm <- resolved2 %>%
  dplyr::group_by(.row_id) %>%
  dplyr::summarise(gene_symbol = collapse_ids(Resolved_UNIPROT), .groups = "drop") %>%
  dplyr::right_join(male.data %>% dplyr::select(-gene_symbol, .row_id), by = ".row_id") %>%
  dplyr::select(-.row_id) %>%
  dplyr::mutate(gene_symbol = dplyr::na_if(gene_symbol, ""))

to_numeric_matrix <- function(male_norm, qc_dir = subdirs$plots_qc) {
  if (!"gene_symbol" %in% names(male_norm)) stop("male.norm must contain gene_symbol")
  expr <- as.data.frame(lapply(male_norm[, -1, drop = FALSE], function(x) suppressWarnings(as.numeric(x))))
  if (!all(vapply(expr, is.numeric, logical(1)))) stop("Non-numeric columns remain after coercion")
  mat <- as.data.frame(t(expr))
  feat <- male_norm$gene_symbol
  empty <- which(!nzchar(ifelse(is.na(feat), "", feat)))
  if (length(empty)) feat[empty] <- paste0("UNMAPPED_", seq_along(empty))
  feat <- make.unique(feat, sep = "_")
  colnames(mat) <- feat
  if (any(!nzchar(colnames(mat)) | is.na(colnames(mat)))) stop("Empty/NA feature names after repair")
  utils::write.table(utils::head(mat[, 1:min(10, ncol(mat)), drop = FALSE]), file.path(qc_dir, "expression_head10.tsv"), sep = "\t", row.names = TRUE, quote = FALSE)
  mat
}
expression.data <- to_numeric_matrix(male.norm)

# Save core outputs
readr::write_tsv(resolved2, file.path(subdirs$tables_modules, "resolved_tokens_mouse_only_robust.tsv"))
saveRDS(list(expression = expression.data, male.norm = male.norm, mapping = resolved2), file = file.path(subdirs$network, "mouse_only_mapping_outputs_robust.rds"))

# Expression matrix: rows = samples, cols = proteins
#expression.data <- male.norm[, -1]
#expression.data <- as.data.frame(lapply(expression.data, as.numeric))
#expression.data <- as.data.frame(t(expression.data))
#names(expression.data) <- male.norm$gene_symbol

# Force single accession per feature name and ensure uniqueness
fix_feature_ids <- function(nms) {
  # take first token before ';'
  first <- sub(";.*$", "", nms)
  # trim + uppercase
  first <- toupper(trimws(first))
  # keep only accession-like strings; if not accession, keep as-is
  is_acc <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[0-9A-Z]{7}$", first)
  first[!nzchar(first)] <- "UNMAPPED"
  # make unique to avoid downstream errors
  make.unique(first, sep = "_")
}

colnames(expression.data) <- fix_feature_ids(colnames(expression.data))


# --------------------------
# QC and sample clustering
# --------------------------
gsg <- goodSamplesGenes(expression.data)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0) printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", ")))
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE]
}

sampleTree <- hclust(dist(expression.data), method = "average")
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.main = 2); abline(h = 40, col = "red")
svg(file = fp_qc("sample_clustering_outliers.svg"), width = 8, height = 6)
dev.off()

#cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
expression.data <- expression.data[cut.sampleTree == 1, ]

# --------------------------
# Soft-threshold selection
# --------------------------
spt <- pickSoftThreshold(expression.data)

svglite::svglite(file = fp_qc("soft_threshold_scale_independence.svg"), width = 7, height = 5)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
text(spt$fitIndices[,1], spt$fitIndices[,2], labels = spt$fitIndices[,1], col = "red"); abline(h = 0.80, col = "red")
dev.off()

svglite::svglite(file = fp_qc("soft_threshold_mean_connectivity.svg"), width = 7, height = 5)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(spt$fitIndices[,1], spt$fitIndices[,5], labels = spt$fitIndices[,1], col = "red")
dev.off()

softPower <- 6

# --------------------------
# Network construction
# --------------------------
adjacency <- adjacency(expression.data, power = softPower, type = "signed",
                       corFnc = "bicor", corOptions = list(use="p", maxPOutliers=0.05))
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

svg(file = fp_net("gene_dendrogram.svg"), width = 12, height = 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2,
                         pamRespectsDendro = FALSE, minClusterSize = 30)

colorSeq <- c(
  "lemon" = "lemonchiffon", "sage" = "darkseagreen", "bluegray" = "steelblue",
  "mintblue" = "lightsteelblue", "azure" = "deepskyblue", "khaki" = "khaki",
  "skyblue" = "skyblue", "babyblue" = "lightblue", "amber" = "goldenrod",
  "tealgreen" = "darkcyan", "forestgreen" = "forestgreen", "gold" = "gold",
  "violet" = "violet", "seafoam" = "mediumaquamarine", "coral" = "coral",
  "salmonlight" = "lightsalmon", "peach" = "peachpuff", "mint" = "palegreen",
  "lime" = "limegreen", "mauve" = "plum", "freesia" = "lightpink", "cocoa" = "saddlebrown", 
  "lavender" = "lavender", "magenta" = "magenta", "salmon" = "salmon",
  "rose" = "mistyrose", "aquamarine" = "aquamarine", "tomato" = "tomato",
  "plum" = "plum", "hotpink" = "hotpink", "rust" = "sienna"
)

# Map numeric module labels -> colors from colorSeq (recycle palette if needed)
unique_mods <- sort(unique(Modules))
nmods <- length(unique_mods)
palette_vals <- unname(colorSeq)
if (length(palette_vals) < nmods) palette_vals <- rep(palette_vals, length.out = nmods)
mod_colors_map <- setNames(palette_vals[seq_len(nmods)], as.character(unique_mods))
ModuleColors <- as.character(mod_colors_map[as.character(Modules)])
stopifnot(length(ModuleColors) == length(Modules))

svg(file = fp_net("gene_dendrogram_module_colors.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, ModuleColors, "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MElist <- moduleEigengenes(expression.data, colors = ModuleColors)
MEs <- MElist$eigengenes
ME.dissimilarity <- 1 - cor(MEs, use = "p", method = "pearson")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
mergedColors <- merge$colors
mergedMEs <- orderMEs(merge$newMEs)

svg(file = fp_net("gene_dendrogram_modules_merged.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors),
                    c("Original Module","Merged Module"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

# --------------------------
# Eigengene network plots
# --------------------------
MET <- orderMEs(mergedMEs)
svg(file = fp_net("eigengene_dendrogram.svg"), width = 8, height = 6)
plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE, marDendro = c(0, 4, 2, 0))
dev.off()
svg(file = fp_net("eigengene_adjacency_heatmap.svg"), width = 8, height = 6)
par(mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", plotDendrograms = FALSE,
                      marHeatmap = c(5, 5, 2, 2), xLabelsAngle = 90)
dev.off()

# --------------------------
# Spatial + condition trait matrix and heatmaps
# --------------------------
sample_info <- readxl::read_excel(path = meta_xlsx)
stopifnot("row.names" %in% names(sample_info))
rownames(sample_info) <- as.character(sample_info$row.names)

Samples <- rownames(expression.data)
sample_info <- sample_info[Samples, , drop = FALSE]

X_celltype <- model.matrix(~ 0 + celltype, data = sample_info); colnames(X_celltype) <- sub("^celltype", "celltype_", colnames(X_celltype))
X_layer    <- model.matrix(~ 0 + layer,    data = sample_info); colnames(X_layer)    <- sub("^layer",    "layer_",    colnames(X_layer))
X_region   <- model.matrix(~ 0 + region,   data = sample_info); colnames(X_region)   <- sub("^region",   "region_",   colnames(X_region))
X_cond     <- model.matrix(~ 0 + ExpGroup, data = sample_info); colnames(X_cond)     <- sub("^ExpGroup", "cond_",     colnames(X_cond))

datTraits <- as.data.frame(cbind(X_celltype, X_layer, X_region, X_cond), stringsAsFactors = FALSE)
keep_cols <- vapply(datTraits, function(x) sd(as.numeric(x), na.rm = TRUE) > 0, logical(1))
datTraits <- datTraits[, keep_cols, drop = FALSE]
rownames(datTraits) <- Samples

nSamples <- nrow(expression.data)
MEcorr <- cor(mergedMEs, datTraits, use = "p", method = "pearson")
MEp    <- corPvalueStudent(MEcorr, nSamples)

plot_trait_heatmap <- function(matCorr, matP, cols, file) {
  if (length(cols) == 0) return(invisible(NULL))
  textMatrix <- paste(signif(matCorr[, cols, drop=FALSE], 2), "\n
                (", signif(matP[, cols, drop=FALSE], 1), ")", sep = "")
  dim(textMatrix) <- dim(matCorr[, cols, drop=FALSE])
  svg(file = file, width = 6, height = max(4, nrow(matCorr) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = matCorr[, cols, drop=FALSE],
                 xLabels = colnames(matCorr)[cols],
                 yLabels = rownames(matCorr),
                 ySymbols = rownames(matCorr),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.6,
                 zlim = c(-1, 1),
                 main = "Module–trait relationships")
  dev.off()
}

trait_names <- colnames(datTraits)
groups <- list(
  celltype = grep("^celltype_", trait_names),
  layer    = grep("^layer_",   trait_names),
  region   = grep("^region_",  trait_names),
  cond     = grep("^cond_",    trait_names)
)
for (nm in names(groups)) {
  idx <- groups[[nm]]
  if (length(idx) > 0) {
    plot_trait_heatmap(MEcorr, MEp, idx, fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")))
  }
}

# --------------------------
# Pairwise condition contrasts (optional)
# --------------------------
mk_contrast <- function(vec, a, b){v<-rep(NA_real_,length(vec));v[vec==a]<-0;v[vec==b]<-1;v}
grp <- as.character(sample_info$ExpGroup)
contrasts <- list(con_res = mk_contrast(grp, "con", "res"),
                  con_sus = mk_contrast(grp, "con", "sus"),
                  res_sus = mk_contrast(grp, "res", "sus"))
for (nm in names(contrasts)) {
  v <- contrasts[[nm]]; keep <- !is.na(v)
  cmat <- cor(mergedMEs[keep, , drop=FALSE], v[keep], use="p")
  pmat <- corPvalueStudent(cmat, sum(keep))
  txt <- paste(signif(cmat,2), "\n(", signif(pmat,1), ")", sep = "")
  dim(txt) <- dim(cmat)
  svg(file = fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")),
      width = 3, height = max(4, ncol(mergedMEs) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = cmat, xLabels = nm,
                 yLabels = colnames(mergedMEs),
                 ySymbols = colnames(mergedMEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = txt,
                 setStdMargins = FALSE,
                 cex.text = 0.8, zlim = c(-1,1),
                 main = "Module–trait relationships")
  dev.off()
}

# --------------------------
# kME, GS, hubs
# --------------------------
modNames <- substring(colnames(mergedMEs), 3)
geneModuleMembership <- as.data.frame(cor(expression.data,
                                          mergedMEs, use = "p",
                                          method = "pearson"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

cond_numeric <- setNames(c(1, 2, 3), c("con", "res", "sus"))
ExpGroup_num <- as.numeric(cond_numeric[as.character(sample_info$ExpGroup)])
geneTraitSignificance <- as.data.frame(cor(expression.data,
                                           ExpGroup_num,
                                           use = "p",
                                           method = "pearson"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- "GS.ExpGroup"
names(GSPvalue) <- "p.GS.ExpGroup"

modules_of_interest <- unique(mergedColors)
for (module in modules_of_interest) {
  moduleGenes <- mergedColors == module
  if (!any(moduleGenes)) next
  mmcol <- paste0("MM", module); if (!mmcol %in% colnames(geneModuleMembership)) next
  gene_info <- data.frame(
    Gene = colnames(expression.data)[moduleGenes],
    Module = module,
    ModuleMembership = geneModuleMembership[moduleGenes, mmcol],
    GeneSignificance = geneTraitSignificance[moduleGenes, "GS.ExpGroup"]
  )
  write.csv(gene_info, file = fp_modtab(paste0("genes_in_module_", module, ".csv")), row.names = FALSE)
}

# --------------------------
# Preservation across conditions
# --------------------------
idx_con <- which(sample_info$ExpGroup == "con")
idx_res <- which(sample_info$ExpGroup == "res")
idx_sus <- which(sample_info$ExpGroup == "sus")

multiExpr <- list(
  ALL = list(data = expression.data),
  CON = list(data = expression.data[idx_con, , drop = FALSE]),
  RES = list(data = expression.data[idx_res, , drop = FALSE]),
  SUS = list(data = expression.data[idx_sus, , drop = FALSE])
)

good_cols <- lapply(multiExpr, function(e) {
  gsg <- goodSamplesGenes(e$data, verbose = 3)
  which(gsg$goodGenes)
})
common_idx <- Reduce(intersect, good_cols)
common_genes <- colnames(multiExpr[[1]]$data)[common_idx]

has_na <- sapply(common_genes, function(g) any(vapply(multiExpr, function(e) any(is.na(e$data[, g])), logical(1))))
common_genes <- common_genes[!has_na]
stopifnot(length(common_genes) > 0)

multi_expr_clean <- lapply(multiExpr, function(e) {
  dat <- e$data[, common_genes, drop = FALSE]
  dat <- as.data.frame(dat)
  dat[] <- lapply(dat, function(col) as.numeric(col))
  keep_samp <- apply(dat, 1, function(r) sd(r, na.rm = TRUE) > 0)
  dat <- dat[keep_samp, , drop = FALSE]
  list(data = as.matrix(dat))
})
names(multi_expr_clean) <- names(multiExpr)

if (any(vapply(multi_expr_clean, function(e) nrow(e$data) < 2, logical(1)))) {
  stop("One or more sets have < 2 samples after cleaning; reduce filtering or combine sets.")
}

ref_colors <- mergedColors[match(common_genes, colnames(expression.data))]
stopifnot(length(ref_colors) == length(common_genes))

multi_color <- list(
  ALL = ref_colors,
  CON = rep("grey", length(common_genes)),
  RES = rep("grey", length(common_genes)),
  SUS = rep("grey", length(common_genes))
)

set.seed(12345)
mp <- modulePreservation(
  multi_expr_clean,
  multi_color,
  referenceNetworks = 1,
  nPermutations = 200,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = "use = 'p', maxPOutliers = 0.05",
  parallelCalculation = TRUE,  # when supported by your WGCNA build
  verbose = 3
)

Z_CON <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.CON
Z_RES <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.RES
Z_SUS <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.SUS

write.csv(Z_CON, fp_prestab("module_preservation_CON_vs_ALL_Zsummary.csv"), row.names = TRUE)
write.csv(Z_RES, fp_prestab("module_preservation_RES_vs_ALL_Zsummary.csv"), row.names = TRUE)
write.csv(Z_SUS, fp_prestab("module_preservation_SUS_vs_ALL_Zsummary.csv"), row.names = TRUE)

# --------------------------
# Condition × region × layer × celltype panels
# --------------------------
# Build combined label; microglia layer -> "none"
region <- as.character(sample_info$region)
layer  <- as.character(sample_info$layer)
cell   <- as.character(sample_info$celltype)
cond   <- as.character(sample_info$ExpGroup)

combo <- paste(cond, region, layer, cell, sep = "_")
combo <- factor(combo)

# One-hot and correlations
X_combo <- model.matrix(~ 0 + combo)
colnames(X_combo) <- levels(combo)  # no "comb" prefix at all
if (ncol(X_combo) == 0) stop("No combined strata present in X_combo (ncol == 0). Check combo label construction.")

MEcorr_combo <- cor(mergedMEs, X_combo, use = "p")
MEp_combo    <- corPvalueStudent(MEcorr_combo, nrow(expression.data))

df_combo <- reshape2::melt(MEcorr_combo, varnames = c("module","comb"), value.name = "r")
p_combo  <- reshape2::melt(MEp_combo,    varnames = c("module","comb"), value.name = "p")
df_combo$p <- p_combo$p
df_combo$stars <- gtools::stars.pval(df_combo$p)
df_combo$module <- factor(df_combo$module, levels = rownames(MEcorr_combo))

# Robust split of combo into columns
split_combo_labels <- function(labels) {
  parts <- strsplit(labels, "_", fixed = TRUE)
  do.call(rbind, lapply(parts, function(p) {
    p <- as.character(p)
    if (length(p) < 4) p <- c(p, rep("missing", 4 - length(p)))
    p[1:4]
  }))
}
lab_parts <- split_combo_labels(as.character(df_combo$comb))
colnames(lab_parts) <- c("condition","region","layer","cell")
df_combo2 <- cbind(df_combo[, c("module","comb","r","p","stars")], as.data.frame(lab_parts, stringsAsFactors = FALSE))

present_conds <- intersect(c("con","res","sus"), unique(df_combo2$condition))
by_cond <- split(df_combo2, df_combo2$condition)

panel_plot <- function(dfi, panel_title = "") {
  if (!nrow(dfi)) return(NULL)
  dfi$trait <- paste(dfi$region, dfi$layer, dfi$cell, sep = "_")
  ord <- order(dfi$region, dfi$layer, dfi$cell)
  dfi$trait <- factor(dfi$trait, levels = unique(dfi$trait[ord]))
  ggplot(dfi, aes(x = trait, y = module, fill = r)) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes(label = stars), size = 3, color = "black", na.rm = TRUE) +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = "#2166ac", mid = "white", high = "#b2182b") +
    labs(title = panel_title, x = NULL, y = NULL, fill = "Pearson r") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank(),
          plot.margin = margin(2, 2, 2, 2),
          legend.position = "none")
}

plots <- lapply(present_conds, function(cn) panel_plot(by_cond[[cn]], toupper(cn)))
plots <- Filter(function(p) !is.null(p) && inherits(p, "ggplot"), plots)
if (length(plots) == 0) stop("No non-empty condition panels; check combo labels and present conditions.")

# Combine panels; if single, save directly
if (length(plots) == 1) {
  g_combined <- plots[[1]]
  ggsave(fp_traits("panel_ME_vs_condition_region_layer_celltype.svg"),
         g_combined, device = svglite::svglite, width = 10, height = 4)
} else {
  g_combined <- patchwork::wrap_plots(plots, nrow = 1)
  legend_df <- data.frame(x = 1:3, y = 1:3, r = c(-1, 0, 1))
  p_legend <- ggplot(legend_df, aes(x, y, fill = r)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = "#2166ac", mid = "white", high = "#b2182b") +
    theme_void() + theme(legend.position = "right") + labs(fill = "Pearson r")
  legend_only <- cowplot::get_legend(p_legend)
  g <- cowplot::plot_grid(g_combined, legend_only, rel_widths = c(1, 0.08))
  ggsave(fp_traits("panel_ME_vs_condition_region_layer_celltype.svg"),
         g, device = svglite::svglite, width = 14, height = 4.5)
}

# --------------------------
# Ensure required objects are defined for the pheatmap block
# --------------------------
if (!exists("df_all")) {
    if (exists("df_combo2")) {
        df_all <- df_combo2
        if (!"trait" %in% names(df_all)) {
            df_all$trait <- paste(df_all$condition, df_all$region, df_all$layer, df_all$cell, sep = "_")
        }
    } else if (exists("MEcorr_combo") && exists("MEp_combo")) {
        df_all <- reshape2::melt(MEcorr_combo, varnames = c("module", "comb"), value.name = "r")
        p_tmp  <- reshape2::melt(MEp_combo,    varnames = c("module", "comb"), value.name = "p")
        df_all$p <- p_tmp$p
        parts <- strsplit(as.character(df_all$comb), "_", fixed = TRUE)
        parts <- lapply(parts, function(x) { x <- as.character(x); length(x) <- 4; x[is.na(x)] <- "missing"; x })
        parts <- do.call(rbind, parts)
        colnames(parts) <- c("condition", "region", "layer", "cell")
        df_all <- cbind(df_all, as.data.frame(parts, stringsAsFactors = FALSE))
        df_all$trait <- paste(df_all$condition, df_all$region, df_all$layer, df_all$cell, sep = "_")
    } else {
        stop("df_all is missing and cannot be reconstructed: provide df_combo2 or MEcorr_combo/MEp_combo.")
    }
}

# module_levels (ordered modules)
if (!exists("module_levels")) {
    if (exists("MEcorr_combo")) {
        module_levels <- rownames(MEcorr_combo)
    } else {
        module_levels <- unique(df_all$module)
    }
}

# traits_in_order (ordered trait columns)
if (!exists("traits_in_order")) {
    if (exists("traits_in_order") && length(traits_in_order) > 0) {
        # already provided
    } else {
        # reasonable default ordering: condition, region, layer, cell
        df_all$trait <- paste(df_all$condition, df_all$region, df_all$layer, df_all$cell, sep = "_")
        ord <- order(df_all$condition, df_all$region, df_all$layer, df_all$cell)
        traits_in_order <- unique(df_all$trait[ord])
    }
}

# strip_df: module -> color mapping
if (!exists("strip_df")) {
    if (exists("mergedColors")) {
        # mergedColors is gene-level mapping; fall back to simple default
        strip_df <- data.frame(module = module_levels, mod_col = rep("grey80", length(module_levels)), stringsAsFactors = FALSE)
    } else {
        strip_df <- data.frame(module = module_levels, mod_col = rep("grey80", length(module_levels)), stringsAsFactors = FALSE)
    }
}

# Ensure output_dir exists
if (!exists("output_dir")) output_dir <- getwd()
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

library(reshape2)
library(grid)

# 1) Wide matrices and order
r_mat  <- reshape2::acast(df_all, module ~ trait, value.var = "r")
p_mat  <- reshape2::acast(df_all, module ~ trait, value.var = "p")

# FDR matrix
fdr_vec <- p.adjust(as.vector(p_mat), method = "BH")
fdr_mat <- matrix(fdr_vec, nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))

stopifnot(all(module_levels %in% rownames(r_mat)), all(traits_in_order %in% colnames(r_mat)))
r_mat   <- r_mat[module_levels, traits_in_order, drop = FALSE]
p_mat   <- p_mat[module_levels, traits_in_order, drop = FALSE]
fdr_mat <- fdr_mat[module_levels, traits_in_order, drop = FALSE]

# 2) Column clustering and separators
hc_cols <- hclust(as.dist(1 - cor(r_mat, use = "pairwise.complete.obs")), method = "average")
k_clusters <- 6
col_grp <- cutree(hc_cols, k = k_clusters)
grp_ord <- col_grp[hc_cols$order]
gap_pos <- which(grp_ord[-1] != head(grp_ord, -1))
xlines <- gap_pos + 0.5

# 3) Module color mapping for row strip
if (exists("mergedMEs")) {
  me_names <- colnames(mergedMEs)                      # e.g., "MEblue"
  row_me <- rownames(r_mat)
  row_mod_colors <- sub("^ME", "", row_me)
  if (!exists("colorSeq")) {
    colorSeq <- c(
      "turquoise"="turquoise","blue"="blue","brown"="brown","yellow"="yellow",
      "green"="green","red"="red","black"="black","pink"="pink","magenta"="magenta",
      "purple"="purple","greenyellow"="greenyellow","tan"="tan","salmon"="salmon",
      "cyan"="cyan","midnightblue"="midnightblue","lightcyan"="lightcyan",
      "greenyellow2"="#ADFF2F","royalblue"="royalblue","darkred"="darkred",
      "skyblue"="skyblue","orange"="orange","grey"="grey"
    )
  }
  color_lookup <- colorSeq
  missing_cols <- setdiff(unique(row_mod_colors), names(color_lookup))
  if (length(missing_cols)) {
    add <- setNames(missing_cols, missing_cols)
    color_lookup <- c(color_lookup, add)
  }
  row_colors_vec <- unname(color_lookup[row_mod_colors])
  row_colors_vec[is.na(row_colors_vec)] <- "grey80"
  ann_row <- data.frame(ModuleColor = row_colors_vec, row.names = rownames(r_mat))
  ann_colors <- list(ModuleColor = setNames(unique(ann_row$ModuleColor), unique(ann_row$ModuleColor)))
} else {
  ann_row <- data.frame(ModuleColor = rep("grey80", nrow(r_mat)), row.names = rownames(r_mat))
  ann_colors <- list(ModuleColor = c("grey80" = "grey80"))
}

# Column condition annotation
col_condition <- sapply(strsplit(colnames(r_mat), "_", fixed = TRUE), `[`, 1)
ann_col <- data.frame(Condition = col_condition, row.names = colnames(r_mat))
ann_colors$Condition <- c(con = "#eaf0ff", res = "#fff5df", sus = "#ffecef")

# 4) Palette and breaks
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", repos = "https://cloud.r-project.org")
r_lim <- 0.8
bk <- seq(-r_lim, r_lim, length.out = 201)
pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(200)

# 5) Significance dots
display_mat <- matrix("", nrow = nrow(r_mat), ncol = ncol(r_mat), dimnames = dimnames(r_mat))
display_mat[fdr_mat < 0.01] <- "•"

# 6) Render pheatmap once and save
ph <- pheatmap::pheatmap(
  mat = r_mat,
  cluster_rows = FALSE,
  cluster_cols = hc_cols,
  treeheight_row = 0,
  treeheight_col = 50,
  gaps_col = NULL,
  annotation_row = ann_row,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  color = pal, breaks = bk, legend = TRUE,
  border_color = NA,
  cellwidth = 10, cellheight = 10,
  angle_col = 45,
  fontsize_col = 8.5, fontsize_row = 8.8,
  display_numbers = display_mat,
  number_color = "grey20",
  number_cex = 0.7,
  show_rownames = TRUE, show_colnames = TRUE,
  silent = TRUE
)

nr <- nrow(r_mat)

png(fp_traits("ME_trait_pheatmap.png"), width = 3200, height = 1800, res = 300)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
panel_id <- grep("matrix", ph$gtable$layout$name)[1]
seekViewport(ph$gtable$layout$name[panel_id])
for (xl in xlines) {
  grid::grid.lines(x = unit(c(xl, xl), "native"),
                   y = unit(c(0, nr), "native"),
                   gp = gpar(col = "white", lwd = 3))
}
upViewport(0); dev.off()

pdf(fp_traits("ME_trait_pheatmap.pdf"), width = 12, height = 7)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
panel_id <- grep("matrix", ph$gtable$layout$name)[1]
seekViewport(ph$gtable$layout$name[panel_id])
for (xl in xlines) {
  grid::grid.lines(x = unit(c(xl, xl), "native"),
                   y = unit(c(0, nr), "native"),
                   gp = gpar(col = "white", lwd = 2))
}
upViewport(0); dev.off()

svglite::svglite(fp_traits("ME_trait_pheatmap.svg"), width = 12, height = 7)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
panel_id <- grep("matrix", ph$gtable$layout$name)[1]
seekViewport(ph$gtable$layout$name[panel_id])
for (xl in xlines) {
  grid::grid.lines(x = unit(c(xl, xl), "native"),
                   y = unit(c(0, nr), "native"),
                   gp = gpar(col = "white", lwd = 2))
}
upViewport(0); dev.off()

library(patchwork)
library(gtools)
library(cowplot)

# Build blocks: celltype, region, layer, condition (no cellclass)
block_list <- list(
  celltype  = grep("^celltype_", colnames(datTraits)),
  region    = grep("^region_",  colnames(datTraits)),
  layer     = grep("^layer_",   colnames(datTraits)),
  condition = grep("^cond_",    colnames(datTraits))
)

block_cor_df <- function(block_idx, block_name) {
  if (length(block_idx) == 0) return(NULL)
  r <- cor(mergedMEs, datTraits[, block_idx, drop=FALSE], use="p")
  p <- corPvalueStudent(r, nrow(expression.data))
  df <- reshape2::melt(r, varnames = c("module","trait"), value.name = "r")
  p_long <- reshape2::melt(p, varnames = c("module","trait"), value.name = "p")
  df$p <- p_long$p
  df$stars <- gtools::stars.pval(df$p)
  df$block <- block_name
  df$module <- factor(df$module, levels = rownames(r))
  df
}

dfs <- Filter(Negate(is.null),
              mapply(block_cor_df, block_list, names(block_list), SIMPLIFY = FALSE))

panel_plot <- function(dfi, legend = "none") {
  if (nrow(dfi) == 0) return(NULL)
  dfi$trait <- factor(dfi$trait, levels = unique(dfi$trait))

  p_heat <- ggplot(dfi, aes(x = trait, y = module, fill = r)) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes(label = stars), size = 3, color = "black", na.rm = TRUE) +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = "#2166ac", mid = "white", high = "#b2182b") +
    labs(x = NULL, y = NULL, fill = "Pearson r") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = legend,
          panel.grid = element_blank(),
          plot.margin = margin(2, 2, 2, 2))

  if (!identical(legend, "none")) return(p_heat)

  modules <- levels(dfi$module)
  mod_cols <- sub("^ME", "", modules)
  strip_df <- data.frame(module = modules, mod_col = mod_cols, stringsAsFactors = FALSE)
  strip_df$module <- factor(strip_df$module, levels = modules)

  p_strip <- ggplot(strip_df, aes(x = 1, y = module, fill = mod_col)) +
    geom_tile() +
    scale_fill_identity() +
    theme_void() +
    theme(plot.margin = margin(2, 0, 2, 2))

  combined <- p_strip + p_heat + plot_layout(widths = c(0.05, 1))
  combined
}

plots <- lapply(dfs, panel_plot)
p_legend_heat <- panel_plot(dfs[[1]], legend = "right")
legend_only <- cowplot::get_legend(p_legend_heat)

combo <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(title = "Module–trait relationships: cell type, region, layer, condition")

svg(file.path(output_dir, "panel_module_trait_relationships_no_cellclass.svg"), width = 14, height = 4.5)
cowplot::plot_grid(combo, legend_only, rel_widths = c(1, 0.08))
dev.off()

# ==========================================================
# NEW: ME by condition plots with significance and exports
# ==========================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggpubr); library(readr); library(broom)
})

# Long-form eigengenes + metadata
stopifnot(nrow(mergedMEs) == nrow(sample_info))
ME_long <- mergedMEs %>%
  tibble::rownames_to_column("Sample") %>%
  mutate(condition = as.character(sample_info$ExpGroup)) %>%
  pivot_longer(cols = starts_with("ME"),
               names_to = "module", values_to = "ME")

# Ensure order con, res, sus
ME_long$condition <- factor(ME_long$condition, levels = c("con","res","sus"))

# Stats
do_stats <- function(df) {
  aov_fit <- stats::aov(ME ~ condition, data = df)
  aov_p <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
  list(aov_p = aov_p)
}

stat_list <- ME_long %>%
  group_by(module) %>%
  group_map(~{
    st <- do_stats(.x)
    tibble(module = unique(.x$module), aov_p = st$aov_p)
  }) %>% bind_rows() %>%
  mutate(aov_fdr = p.adjust(aov_p, method = "BH")) %>%
  arrange(aov_fdr)

readr::write_csv(stat_list, fp_modtab("ME_by_condition_ANOVA_FDR.csv"))

top_modules <- head(stat_list$module, 12)
comparisons <- list(c("con","res"), c("con","sus"), c("res","sus"))

# Color mapping (full circles)
cond_cols <- c(con = "#457B9D", res = "#C6C3BB", sus = "#E63946")

# Helper: build one dotplot
plot_dot_mod <- function(dfm, mod, aov_fdr = NA_real_, show_subtitle = TRUE, errorbar = c("none","sem","sd")) {
  errorbar <- match.arg(errorbar)
  # Summary for error bars
  summ <- dfm %>%
    group_by(condition) %>%
    summarise(mean = mean(ME, na.rm = TRUE),
              sd = sd(ME, na.rm = TRUE),
              n = dplyr::n(),
              se = sd / sqrt(pmax(n, 1)), .groups = "drop")
  # Base: all sample dots (full circles)
  p <- ggplot(dfm, aes(x = condition, y = ME, color = condition)) +
    geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1),
               size = 2, alpha = 0.85, shape = 16, stroke = 0) +
    scale_color_manual(values = cond_cols, guide = "none")

  # Add mean points (slightly larger) on top
  p <- p + geom_point(data = summ, aes(x = condition, y = mean),
                      inherit.aes = FALSE, size = 3.2, shape = 16, color = "black") +
           geom_point(data = summ, aes(x = condition, y = mean, color = condition),
                      inherit.aes = FALSE, size = 2.8, shape = 16)

  # Optional error bars
  if (errorbar != "none") {
    if (errorbar == "sem") {
      p <- p + geom_errorbar(data = summ,
                             aes(x = condition, ymin = mean - se, ymax = mean + se, color = condition),
                             inherit.aes = FALSE, width = 0.1, alpha = 0.9)
    } else if (errorbar == "sd") {
      p <- p + geom_errorbar(data = summ,
                             aes(x = condition, ymin = mean - sd, ymax = mean + sd, color = condition),
                             inherit.aes = FALSE, width = 0.1, alpha = 0.9)
    }
  }

  subtitle_txt <- if (isTRUE(show_subtitle))
    sprintf("ANOVA FDR=%s", ifelse(is.na(aov_fdr), "NA", signif(aov_fdr, 3))) else NULL

  p <- p +
    labs(title = paste0(mod, " eigengene by condition"),
         subtitle = subtitle_txt,
         x = NULL, y = "Module eigengene") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.x = element_blank())

  # Pairwise significance labels (Wilcoxon BH), shown as p.signif above groups
  p <- p + ggpubr::stat_compare_means(comparisons = comparisons,
                                      method = "wilcox.test",
                                      p.adjust.method = "BH",
                                      label = "p.signif",
                                      hide.ns = TRUE)
  p
}

# Multi-page PDF: all modules as dotplots
pdf(fp_traits("ME_by_condition_all_modules_dotplot.pdf"), width = 7, height = 4.5)
for (mod in unique(ME_long$module)) {
  dfm <- ME_long %>% filter(module == mod)
  aov_fdr <- stat_list$aov_fdr[stat_list$module == mod][1]
  p <- plot_dot_mod(dfm, mod, aov_fdr = aov_fdr, show_subtitle = TRUE, errorbar = "sem")
  print(p)
}
dev.off()

# SVG grid: top N most differential modules (compact)
plot_one_mod_dot <- function(mod) {
  dfm <- ME_long %>% filter(module == mod)
  aov_fdr <- stat_list$aov_fdr[stat_list$module == mod][1]
  summ <- dfm %>% group_by(condition) %>%
    summarise(mean = mean(ME, na.rm = TRUE), se = sd(ME, na.rm = TRUE)/sqrt(n()), .groups = "drop")

  ggplot(dfm, aes(x = condition, y = ME, color = condition)) +
    geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1),
               size = 1.8, alpha = 0.8, shape = 16, stroke = 0) +
    geom_point(data = summ, aes(x = condition, y = mean),
               inherit.aes = FALSE, size = 3, shape = 16, color = "black") +
    geom_point(data = summ, aes(x = condition, y = mean, color = condition),
               inherit.aes = FALSE, size = 2.6, shape = 16) +
    # Optional tiny SEM bars for overview
    geom_errorbar(data = summ, aes(x = condition, ymin = mean - se, ymax = mean + se, color = condition),
                  inherit.aes = FALSE, width = 0.08, alpha = 0.8) +
    scale_color_manual(values = cond_cols, guide = "none") +
    labs(title = paste0(mod, "  (FDR=", signif(aov_fdr, 3), ")"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 9),
          panel.grid.major.x = element_blank())
}

if (length(top_modules) > 0) {
  plots_top <- lapply(top_modules, plot_one_mod_dot)
  ncol_grid <- min(4, ceiling(sqrt(length(plots_top))))
  nrow_grid <- ceiling(length(plots_top)/ncol_grid)
  g <- patchwork::wrap_plots(plots_top, ncol = ncol_grid)
  svglite::svglite(fp_traits("ME_by_condition_top_modules_dotplot.svg"),
                   width = 3.0*ncol_grid, height = 2.5*nrow_grid)
  print(g)
  dev.off()
}

# Optional: export pairwise Wilcoxon BH summary across modules
pw_tables <- ME_long %>%
  group_by(module) %>%
  group_map(~{
    tt <- pairwise.wilcox.test(.x$ME, .x$condition, p.adjust.method = "BH", exact = FALSE)
    broom::tidy(tt) %>% mutate(module = unique(.x$module))
  }) %>% bind_rows()
readr::write_csv(pw_tables, fp_modtab("ME_by_condition_pairwise_Wilcoxon_BH.csv"))


# --------------------------
# End of script
# --------------------------

