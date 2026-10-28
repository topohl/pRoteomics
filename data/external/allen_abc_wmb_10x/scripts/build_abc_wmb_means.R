# Build per-subclass and per-supertype mean expression tables from the ABC Atlas
# whole-mouse-brain 10x taxonomy (Yao et al. 2023; taxonomy CCN20230722, metadata 20231215).
#
# Input (downloaded and verified by download_abc_wmb.R into ../raw/):
#   - MapMyCells precomputed stats (precomputed_stats_ABC_revision_230821.h5). Its `sum` dataset
#     holds the AIT21.0 cluster means of log2(CPM+1) ("src": cl.means.freeze.csv) with
#     n_cells = 1 for every cluster, so sum == cluster mean.
#   - taxonomy tables (cluster -> supertype -> subclass -> class, neurotransmitter),
#     cl.df_CCN202307220.xlsx (anatomical annotation, MERFISH CCF frequencies, markers),
#     gene.csv, cell_metadata.csv (dissection region and 10x chemistry per cell).
#
# Output (this folder):
#   abc_wmb_subclass_mean_log2cpm.csv    genes x 338 subclasses
#   abc_wmb_supertype_mean_log2cpm.csv   genes x 1201 supertypes
#       value = cell-weighted mean of cluster means = sum_k(n_k * mean_k) / sum_k(n_k),
#       n_k = number_of_cells of cluster k in the official taxonomy (10Xv2 + 10Xv3 + multiome)
#   abc_wmb_subclass_types.csv / abc_wmb_supertype_types.csv   one row per type: names,
#       parents, neurotransmitter, n_clusters, n_cells (+ per chemistry), dissection-region
#       fractions, MERFISH CCF regions, curated anatomical annotation, markers
#   abc_wmb_taxonomy.csv                 one row per cluster (5322)
#   abc_wmb_genes.csv                    gene_identifier (Ensembl), gene_symbol, ...
#
# Usage: Rscript build_abc_wmb_means.R   (single-threaded; peak memory ~6 GB)

suppressPackageStartupMessages({
  library(rhdf5)
  library(data.table)
  library(readxl)
  library(Matrix)
  library(jsonlite)
})
setDTthreads(4)

args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- if (length(args_file)) dirname(normalizePath(args_file)) else getwd()
out_dir <- if (basename(here) == "scripts") dirname(here) else here   # the reference folder (data/external/allen_abc_wmb_10x)
raw <- normalizePath(Sys.getenv("ABC_WMB_RAW", file.path(out_dir, "raw")), winslash = "/")   # where download_abc_wmb.R puts them
p <- function(...) file.path(raw, ...)

h5_path <- p("mapmycells/WMB-taxonomy/20240831/precomputed_stats_ABC_revision_230821.h5")
tax_dir <- p("metadata/WMB-taxonomy/20231215")
digits <- 3

# ---- genes -------------------------------------------------------------------------------
genes <- fread(p("metadata/WMB-10X/20241115/gene.csv"))
col_names <- fromJSON(h5read(h5_path, "col_names"))
stopifnot(identical(col_names, genes$gene_identifier))
genes[, symbol_duplicated := gene_symbol %in% gene_symbol[duplicated(gene_symbol)]]

# ---- cluster means -----------------------------------------------------------------------
meta <- fromJSON(h5read(h5_path, "metadata"))
n_cells_h5 <- h5read(h5_path, "n_cells")
stopifnot(all(n_cells_h5 == 1))          # sum == cluster mean
cluster_to_row <- unlist(fromJSON(h5read(h5_path, "cluster_to_row")))
S <- h5read(h5_path, "sum")              # genes x clusters (HDF5 rows = clusters)
stopifnot(nrow(S) == length(col_names), ncol(S) == length(cluster_to_row))
colnames(S) <- names(cluster_to_row)[order(cluster_to_row)]
h5closeAll()

# ---- taxonomy ----------------------------------------------------------------------------
mem <- fread(file.path(tax_dir, "cluster_to_cluster_annotation_membership.csv"))
term <- fread(file.path(tax_dir, "cluster_annotation_term.csv"))
cl_size <- fread(file.path(tax_dir, "cluster.csv"))
level_of <- c(CCN20230722_NEUR = "neurotransmitter", CCN20230722_CLAS = "class",
              CCN20230722_SUBC = "subclass", CCN20230722_SUPT = "supertype",
              CCN20230722_CLUS = "cluster")
mem[, level := level_of[cluster_annotation_term_set_label]]
wide_label <- dcast(mem, cluster_alias ~ level, value.var = "cluster_annotation_term_label")
wide_name <- dcast(mem, cluster_alias ~ level, value.var = "cluster_annotation_term_name")
setnames(wide_label, setdiff(names(wide_label), "cluster_alias"),
         paste0(setdiff(names(wide_label), "cluster_alias"), "_label"))
setnames(wide_name, setdiff(names(wide_name), "cluster_alias"),
         paste0(setdiff(names(wide_name), "cluster_alias"), "_name"))
tax <- merge(wide_label, wide_name, by = "cluster_alias")
tax <- merge(tax, cl_size[, .(cluster_alias, n_cells = number_of_cells)], by = "cluster_alias")
setnames(tax, "neurotransmitter_name", "neurotransmitter")
tax[, neurotransmitter_label := NULL]

cl_df <- as.data.table(read_excel(file.path(tax_dir, "cl.df_CCN202307220.xlsx"),
                                  sheet = "cluster_annotation", .name_repair = "unique_quiet"))
stopifnot("CTX.subclass_id...39" %in% names(cl_df))   # duplicated header: 2nd = Yao 2021 subclass label
cl_df <- cl_df[, .(cluster_label = cell_set_accession.cluster, nt_type_label, nt_type_combo_label,
                   anatomical_annotation, CCF_broad.freq, CCF_acronym.freq,
                   n_cells_10Xv3 = v3.size, n_cells_10Xv2 = v2.size,
                   n_cells_multiome = multiome.size, nt.markers, np.markers,
                   cluster.markers.combo, merfish.markers.combo,
                   yao2021_ctx_hpf_cluster = CTX.cluster_label,
                   yao2021_ctx_hpf_subclass = CTX.subclass_id...39)]
tax <- merge(tax, cl_df, by = "cluster_label", all.x = TRUE)
stopifnot(nrow(tax) == 5322, !anyNA(tax$n_cells),
          all(tax$n_cells == tax$n_cells_10Xv3 + tax$n_cells_10Xv2 + tax$n_cells_multiome),
          setequal(tax$cluster_label, colnames(S)))

# ---- dissection region (region_of_interest) per cluster from cell metadata ---------------
cells <- fread(p("metadata/WMB-10X/20241115/cell_metadata.csv"),
               select = c("cluster_alias", "region_of_interest_acronym", "library_method"))
roi <- cells[, .N, by = .(cluster_alias, roi = region_of_interest_acronym)]
chem <- cells[, .N, by = .(cluster_alias, library_method)]
rm(cells); invisible(gc())
stopifnot(all(chem[, sum(N), by = cluster_alias][order(cluster_alias)]$V1 ==
              tax[order(cluster_alias)]$n_cells))

fmt_freq <- function(names, n, cover = 0.8, min_frac = 0.02) {
  ok <- !is.na(n) & n > 0
  if (!any(ok)) return(NA_character_)
  names <- names[ok]; n <- n[ok]
  o <- order(-n); fr <- n[o] / sum(n); nm <- names[o]
  keep <- seq_len(max(1, min(which(cumsum(fr) >= cover - 1e-9))))
  keep <- keep[fr[keep] >= min_frac | keep == 1]
  paste0(nm[keep], ":", sprintf("%.2f", fr[keep]), collapse = ",")
}
roi_cluster <- roi[, .(dissection_roi_freq = fmt_freq(roi, N)), by = cluster_alias]
tax <- merge(tax, roi_cluster, by = "cluster_alias")

# ---- order clusters as in the h5 and aggregate ----------------------------------------------
setkey(tax, cluster_label)
tax <- tax[colnames(S)]
w <- tax$n_cells

parse_freq <- function(x) {   # "A:0.5,B:0.3" -> data.table(region, frac)
  if (is.na(x) || !nzchar(x)) return(data.table(region = character(), frac = numeric()))
  # split on commas not followed by a space ("CUL4, 5" is one CCF structure name)
  parts <- strsplit(strsplit(x, ",(?! )", perl = TRUE)[[1]], ":", fixed = TRUE)
  data.table(region = vapply(parts, `[`, "", 1),
             frac = as.numeric(vapply(parts, `[`, "", 2)))
}
ccf_long <- tax[, parse_freq(CCF_broad.freq), by = .(cluster_alias, n_cells)]
ccf_fine_long <- tax[, parse_freq(CCF_acronym.freq), by = .(cluster_alias, n_cells)]

collapse_unique <- function(x, sep = "; ") {
  x <- unique(na.omit(x)); paste(x[nzchar(x)], collapse = sep)
}

aggregate_level <- function(level) {
  lab <- paste0(level, "_label"); nm <- paste0(level, "_name")
  groups <- unique(tax[, c(lab, nm), with = FALSE])
  groups <- groups[order(get(nm))]
  idx <- match(tax[[lab]], groups[[lab]])
  W <- sparseMatrix(i = seq_along(idx), j = idx, x = w,
                    dims = c(length(idx), nrow(groups)))
  n_type <- colSums(W)
  M <- as.matrix(S %*% W) / rep(n_type, each = nrow(S))
  colnames(M) <- groups[[nm]]

  # per-type annotation table
  key <- tax[, .(cluster_alias, type = get(lab))]
  roi_t <- merge(roi, key, by = "cluster_alias")[, .(N = sum(N)), by = .(type, roi)]
  roi_t <- roi_t[, .(dissection_roi_freq = fmt_freq(roi, N)), by = type]
  # MERFISH CCF frequencies are only listed until 80% coverage per cluster; weight the
  # listed fractions by cluster size (approximate type-level spatial distribution)
  agg_ccf <- function(long, colname) {
    x <- merge(long, key, by = "cluster_alias")[, .(N = sum(frac * n_cells)), by = .(type, region)]
    setnames(x[, .(v = fmt_freq(region, N, cover = 0.8, min_frac = 0.03)), by = type], "v", colname)
  }
  ccf_t <- agg_ccf(ccf_long, "merfish_ccf_broad_freq")
  ccf_f <- agg_ccf(ccf_fine_long, "merfish_ccf_fine_freq")
  info <- tax[, .(
    n_clusters = .N,
    n_cells = sum(n_cells),
    n_cells_10Xv3 = sum(n_cells_10Xv3),
    n_cells_10Xv2 = sum(n_cells_10Xv2),
    n_cells_multiome = sum(n_cells_multiome),
    neurotransmitter = collapse_unique(neurotransmitter, sep = "|"),
    nt_type_combo = collapse_unique(nt_type_combo_label, sep = "|"),
    anatomical_annotation = collapse_unique(anatomical_annotation),
    yao2021_ctx_hpf_subclass = collapse_unique(yao2021_ctx_hpf_subclass)
  ), by = c(lab, nm, if (level == "supertype") c("subclass_label", "subclass_name"),
            "class_label", "class_name")]
  info[, frac_10Xv3 := round(n_cells_10Xv3 / n_cells, 3)]
  info <- merge(info, roi_t, by.x = lab, by.y = "type")
  info <- merge(info, ccf_t, by.x = lab, by.y = "type", all.x = TRUE)
  info <- merge(info, ccf_f, by.x = lab, by.y = "type", all.x = TRUE)
  stopifnot(all(info[match(groups[[lab]], info[[lab]])]$n_cells == n_type))
  info <- info[match(groups[[lab]], info[[lab]])]
  list(M = M, info = info)
}

markers <- list(
  subclass = as.data.table(read_excel(file.path(tax_dir, "cl.df_CCN202307220.xlsx"),
                                      sheet = "subclass_annotation"))[
    , .(subclass_label = cell_set_accession.subclass, subclass.markers.combo,
        subclass.tf.markers.combo)],
  supertype = as.data.table(read_excel(file.path(tax_dir, "cl.df_CCN202307220.xlsx"),
                                       sheet = "supertype_annotation"))[
    , .(supertype_label = cell_set_accession.supertype, supertype.markers.combo,
        `supertype.markers.combo (within subclass)`)]
)

write_matrix <- function(M, file) {
  dt <- data.table(gene_symbol = genes$gene_symbol, gene_identifier = genes$gene_identifier)
  dt <- cbind(dt, as.data.table(round(M, digits)))
  fwrite(dt, file.path(out_dir, file))
}

for (level in c("subclass", "supertype")) {
  res <- aggregate_level(level)
  lab <- paste0(level, "_label")
  info <- merge(res$info, markers[[level]], by = lab, all.x = TRUE, sort = FALSE)
  info <- info[match(res$info[[lab]], info[[lab]])]
  write_matrix(res$M, sprintf("abc_wmb_%s_mean_log2cpm.csv", level))
  fwrite(info, file.path(out_dir, sprintf("abc_wmb_%s_types.csv", level)))
  message(level, ": ", ncol(res$M), " types, ", nrow(res$M), " genes")
}

# ---- cluster-level taxonomy and gene table -----------------------------------------------
tax_out <- tax[, .(cluster_label, cluster_alias, cluster_name, supertype_label, supertype_name,
                   subclass_label, subclass_name, class_label, class_name, neurotransmitter,
                   nt_type_combo_label, n_cells, n_cells_10Xv3, n_cells_10Xv2, n_cells_multiome,
                   anatomical_annotation, dissection_roi_freq, merfish_ccf_broad_freq = CCF_broad.freq,
                   merfish_ccf_fine_freq = CCF_acronym.freq, nt.markers, np.markers,
                   cluster.markers.combo, merfish.markers.combo,
                   yao2021_ctx_hpf_cluster, yao2021_ctx_hpf_subclass)]
fwrite(tax_out[order(cluster_name)], file.path(out_dir, "abc_wmb_taxonomy.csv"))
fwrite(genes[, .(gene_identifier, gene_symbol, gene_name = name, mapped_ncbi_identifier,
                 symbol_duplicated)], file.path(out_dir, "abc_wmb_genes.csv"))

writeLines(c(
  sprintf("built_at: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  sprintf("h5_source: %s (generated_on %s)", meta$src, meta$generated_on),
  "value: cell-weighted mean over clusters of the cluster mean log2(CPM+1)",
  "weights: cluster number_of_cells (WMB-taxonomy 20231215 cluster.csv; 10Xv2+10Xv3+multiome)",
  sprintf("rounding: %d decimals", digits),
  sprintf("R: %s", R.version.string)
), file.path(out_dir, "build_info.txt"))
message("done")
