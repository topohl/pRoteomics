# ================== PCA analysis and plotting (extended, fixed, organized) ==================
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

pacman::p_load(
  data.table, ggplot2, factoextra, reshape2, stats, ggrepel, tools,
  grid, aricode
)

set.seed(42)

# =============== Config =================
gct_file   <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/E9_pg_matrix_protigy.gct"
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/pca_plots"
#gct_file <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/pg.matrix_filtered_70percent-onegroup_imputed_ANOVA_z-scored.gct"
#output_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/pca_plots"
if (!file.exists(gct_file)) stop(sprintf("Input file not found: %s", gct_file))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============== Helpers =================
trim_ws <- function(x){
  if (is.null(x)) return(x)
  x <- as.character(x)
  x <- gsub("[\u00A0\u2007\u202F]", " ", x, perl = TRUE)
  x <- gsub("^\\s+|\\s+$", "", x, perl = TRUE)
  x
}

# robust writer: prefer fwrite, else write.csv
write_dt <- function(df, path){
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(df, path)
  } else {
    utils::write.csv(df, path, row.names = FALSE)
  }
  invisible(path)
}

# pheatmap saver that forces supported extensions
save_pheatmap <- function(mat, file, width=8, height=10, scale="row"){
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("pheatmap not installed, skipping heatmap: ", file)
    return(invisible(NULL))
  }
  ext <- tools::file_ext(file)
  if (!nzchar(ext) || !(tolower(ext) %in% c("png","pdf","tiff","bmp","jpeg","jpg"))) {
    file <- sub("\\.[A-Za-z0-9]+$", "", file)
    file <- paste0(file, ".png")
  }
  pheatmap::pheatmap(
    mat, show_colnames = FALSE, scale = scale, clustering_method = "complete",
    color = colorRampPalette(c("#2c7fb8","#f7f7f7","#d95f0e"))(101),
    filename = file, width = width, height = height
  )
  invisible(file)
}

# -------- Subfolder helpers ----------
subdir <- function(...){ file.path(output_dir, ...) }
ensure_dir <- function(path){ dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

save_plot <- function(subfolder, filename, plot, width=7.5, height=6.2, dpi=150){
  d <- ensure_dir(subdir(subfolder))
  ggsave(filename = filename, plot = plot, path = d, width = width, height = height, dpi = dpi)
  invisible(file.path(d, filename))
}

save_table <- function(subfolder, filename, df, row.names=FALSE){
  d <- ensure_dir(subdir(subfolder))
  f <- file.path(d, filename)
  if (exists("write_dt")) write_dt(df, f) else utils::write.csv(df, f, row.names = row.names)
  invisible(f)
}

# =============== Reader with audits =================
read_gct_meta_rows_checked <- function(path){
  cat("Reading GCT:", path, "\n")
  dt <- data.table::fread(path, header = FALSE, sep = "\t", data.table = TRUE,
              na.strings = c("NA","NaN","","Inf","-Inf"))
  dropN <- 0L
  while (nrow(dt) > 0 && all(is.na(dt[1,]))) { dt <- dt[-1]; dropN <- dropN + 1L }
  cat("Dropped leading empty rows:", dropN, "\n")

  if (!identical(trim_ws(dt[1, V1][[1]]), "#1.3")) {
    warning("Top marker not '#1.3' after dropping blanks; continuing but verify file format.")
  }

  meta_start <- 3L
  meta_end   <- 19L
  expr_start <- 20L
  if (nrow(dt) < expr_start) stop("Unexpected file length; too short for meta+expr blocks.")

  meta_block <- dt[meta_start:meta_end]
  expr_block <- dt[expr_start:nrow(dt)]
  cat("Meta block rows:", nrow(meta_block), " | Expr block rows:", nrow(expr_block), "\n")

  if (ncol(meta_block) < 2) stop("Metadata block has no value columns (need >=2).")
  keep_cols <- which(colSums(!is.na(meta_block)) > 0)
  if (length(keep_cols) < 2) stop("Metadata has <2 non-empty columns.")
  meta_block <- meta_block[, ..keep_cols]
  ncols_meta <- ncol(meta_block)
  nS <- ncols_meta - 1L
  if (nS < 1L) stop("nS computed from meta width is <1; abort.")
  cat("Detected nS (samples):", nS, "from metadata width\n")

  keys <- trim_ws(meta_block[[1]])
  cat("Meta keys detected:\n")
  print(keys)

  meta_wide <- data.frame(matrix(NA_character_, nrow = nS, ncol = 0), check.names = FALSE)
  for (r in seq_along(keys)) {
    key <- keys[r]
    if (!nzchar(key)) next
    vals <- trim_ws(as.character(meta_block[r, 2:(1+nS), with = FALSE]))
    if (length(vals) != nS) stop(sprintf("Meta row '%s' has %d values; expected %d.", key, length(vals), nS))
    meta_wide[[key]] <- vals
  }

  audit_fields <- intersect(c("AnimalID","ReplicateGroup","celltype","layer","region","group","group2","ExpGroup","plate","sampleNumber","shortname"),
                            names(meta_wide))
  cat("Audit: first 3 samples from meta_wide fields:\n")
  print(head(meta_wide[, audit_fields, drop = FALSE], 3))

  keep_cols_expr <- which(colSums(!is.na(expr_block)) > 0)
  if (length(keep_cols_expr) < (1 + nS)) {
    stop(sprintf("Expr block has only %d non-empty cols; need at least %d for values.", length(keep_cols_expr), 1 + nS))
  }
  expr_block_use <- expr_block[, ..keep_cols_expr]
  if (ncol(expr_block_use) < (1 + nS)) stop("Expr usable columns < 1+nS; abort.")
  proteins <- trim_ws(expr_block_use[[1]])
  expr_vals <- expr_block_use[, 2:(1+nS), with = FALSE]
  expr_mat  <- as.matrix(expr_vals)
  mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- proteins
  cat("Expr matrix dim:", paste(dim(expr_mat), collapse = " x "), "\n")

  pick <- function(nm) if (nm %in% names(meta_wide)) trim_ws(meta_wide[[nm]]) else NULL
  short    <- pick("shortname")
  sampleNo <- pick("sampleNumber")
  animal   <- pick("AnimalID")
  group2   <- pick("group2")
  if (!is.null(short) && all(nzchar(short))) {
    new_names <- make.unique(short); picked <- "shortname"
  } else if (!is.null(sampleNo) && all(nzchar(sampleNo))) {
    new_names <- make.unique(sampleNo); picked <- "sampleNumber"
  } else if (!is.null(animal) && !is.null(group2)) {
    new_names <- make.unique(paste0(animal, "_", group2)); picked <- "AnimalID_group2"
  } else { new_names <- sprintf("S%03d", seq_len(nS)); picked <- "S###" }
  colnames(expr_mat) <- new_names
  cat("Sample name source:", picked, "\n")
  cat("First 3 sample names:", paste(head(new_names, 3), collapse = ", "), "\n")

  field <- function(nm) if (nm %in% names(meta_wide)) trim_ws(meta_wide[[nm]]) else rep(NA_character_, nS)
  meta <- data.frame(
    sample         = new_names,
    AnimalID       = field("AnimalID"),
    group2         = field("group2"),
    ReplicateGroup = field("ReplicateGroup"),
    celltype       = field("celltype"),
    celltype_layer = field("celltype_layer"),
    layer          = field("layer"),
    region         = field("region"),
    group          = field("group"),
    ExpGroup       = field("ExpGroup"),
    plate          = field("plate"),
    sampleNumber   = field("sampleNumber"),
    shortname      = field("shortname"),
    stringsAsFactors = FALSE,
    row.names = new_names
  )

  meta$celltype       <- tolower(trim_ws(meta$celltype))
  meta$layer          <- tolower(trim_ws(meta$layer))
  meta$region         <- toupper(trim_ws(meta$region))
  meta$ReplicateGroup <- tools::toTitleCase(tolower(trim_ws(meta$ReplicateGroup)))

  exp_cols <- c("AnimalID","ReplicateGroup","celltype","layer", "celltype_layer", "region","group","group2","ExpGroup","sampleNumber","shortname")
  present  <- exp_cols[exp_cols %in% names(meta)]
  cat("Meta presence check:\n"); print(present)
  if (nrow(meta) != ncol(expr_mat)) stop("meta rows != expr cols; alignment broken.")

  cat("First sample meta snapshot:\n")
  print(meta[1, present, drop = FALSE])

  list(mat = expr_mat, meta = meta)
}

# ================== Build mat/meta with checks ==================
g <- read_gct_meta_rows_checked(gct_file)
mat  <- g$mat
meta <- g$meta

# ================== PCA pipeline ==================
mat[!is.finite(mat)] <- NA
if (!is.matrix(mat) || nrow(mat) < 2 || ncol(mat) < 2)
  stop(sprintf("Parsed matrix has dim %s; need at least 2x2.", paste(dim(mat), collapse="x")))

keep_cols <- colMeans(is.na(mat)) <= 0.80
if (!all(keep_cols)) {
  message("Dropping high-missing samples: ", paste(colnames(mat)[!keep_cols], collapse = ", "))
  mat  <- mat[, keep_cols, drop = FALSE]
  meta <- meta[colnames(mat), , drop = FALSE]
}

keep_rows <- rowSums(is.na(mat)) < ncol(mat)
mat <- mat[keep_rows, , drop = FALSE]
if (nrow(mat) < 2 || ncol(mat) < 2)
  stop(sprintf("Matrix too small after filtering: %dx%d.", nrow(mat), ncol(mat)))

vz <- apply(mat, 1, stats::var, na.rm = TRUE)
if (any(!is.finite(vz) | vz == 0)) mat <- mat[which(vz > 0 & is.finite(vz)), , drop = FALSE]

if (requireNamespace("impute", quietly = TRUE)) {
  mat <- impute::impute.knn(mat)$data
} else {
  med <- apply(mat, 1, median, na.rm = TRUE)
  idx <- which(is.na(mat), arr.ind = TRUE)
  if (nrow(idx) > 0) mat[idx] <- med[idx[,1]]
}

pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)

# Align names
clean_names <- trim_ws(colnames(mat))
colnames(mat)    <- clean_names
rownames(pca$x)  <- clean_names
rownames(meta)   <- trim_ws(rownames(meta))
meta             <- meta[rownames(pca$x), , drop = FALSE]
stopifnot(identical(rownames(meta), rownames(pca$x)))

# Derived labels
meta$region_layer        <- paste(meta$region, meta$layer, meta$ReplicateGroup, sep = "_")
meta$region_only_layer   <- paste(meta$region, meta$layer, sep = "_")
meta$region_celltype_rep <- paste(meta$region, meta$celltype, meta$ReplicateGroup, sep = "_")

# ================== Minimal modern styling ==================
theme_pca_min <- function() {
  theme_minimal(base_size = 16, base_family = "sans") +
    theme(
      panel.grid.major = element_line(color = "#ECECEC", size = 0.4),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      axis.title = element_text(color = "#444444", size = 18, face = "plain"),
      axis.text  = element_text(color = "#555555", size = 18),
      axis.ticks = element_line(color = "#DDDDDD", size = 0.3),
      axis.line  = element_blank(),
      panel.border = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = "right",
      legend.title = element_text(color = "#444444", size = 12),
      legend.text  = element_text(color = "#666666", size = 13),
      plot.title   = element_text(color = "#222222", face = "bold", size = 16, margin = margin(b = 6)),
      plot.subtitle= element_text(color = "#666666", size = 12, margin = margin(b = 6)),
      plot.caption = element_text(color = "#999999", size = 10),
      strip.text = element_text(color = "#333333", size = 10, face = "bold"),
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(8, 8, 6, 8),
      complete = TRUE
    )
}

base_hex <- c("#F8D247", "#DEC196", "#B2B2B2", "#C7D745", "#C6B18B",
              "#C592C5", "#A89BB0", "#E89369", "#66C1A4", "#FF6F61",
              "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251")
make_modern_palette <- function(n){
  if (n <= length(base_hex)) return(base_hex[seq_len(n)])
  grDevices::colorRampPalette(base_hex, space = "Lab")(n)
}

build_group <- function(key){
  if (!key %in% names(meta)) stop(sprintf("Key '%s' not in meta.", key))
  v <- trim_ws(as.character(meta[[key]]))
  v[v == ""] <- NA
  factor(v)
}

plot_and_save_group <- function(key, title_prefix, out_file, point_size = 8) {
  grp <- build_group(key)
  keep <- !is.na(grp)
  if (sum(keep) < 2) { message(sprintf("Skipping '%s': <2 non-NA samples.", key)); return(invisible(NULL)) }
  grp2 <- droplevels(grp[keep]); nlev <- nlevels(grp2); pal <- make_modern_palette(nlev)

  X <- pca; X$x <- X$x[keep, , drop = FALSE]

  if (!is.null(pca$sdev) && length(pca$sdev) >= 2) {
    varp <- (pca$sdev^2) / sum(pca$sdev^2)
    lab_x <- sprintf("PC1 (%.1f%%)", varp[1] * 100)
    lab_y <- sprintf("PC2 (%.1f%%)", varp[2] * 100)
  } else {
    lab_x <- "PC1"; lab_y <- "PC2"
  }

  p <- fviz_pca_ind(
    X,
    geom = "point",
    habillage = grp2,
    addEllipses = FALSE,
    ellipse.type = "t",
    ellipse.level = 0.5,
    palette = pal,
    repel = TRUE,
    mean.point = FALSE,
    title = sprintf("%s by %s", title_prefix, key)
  )

  if (length(p$layers) && inherits(p$layers[[1]]$geom, "GeomPoint")) {
    p$layers[[1]] <- ggplot2::geom_point(mapping = p$layers[[1]]$mapping, inherit.aes = TRUE, shape = 16, size = point_size, alpha = 0.8)
  } else {
    p <- p + ggplot2::geom_point(size = point_size, shape = 16, alpha = 0.8)
  }

  p <- p +
    theme_pca_min() +
    theme(
      axis.line.x = element_line(color = "#E0E0E0"),
      axis.line.y = element_line(color = "#E0E0E0")
    ) +
    labs(x = lab_x, y = lab_y, subtitle = NULL, caption = NULL)

  # Organized save
  save_plot("plots/base", out_file, p)
  p
}

# ================== Generate and save base plots (SVG) ==================
plot_and_save_group("ReplicateGroup",       "PCA", "pca_by_replicate_group.svg")
plot_and_save_group("celltype",             "PCA", "pca_by_celltype.svg")
plot_and_save_group("layer",                "PCA", "pca_by_layer.svg")
plot_and_save_group("celltype_layer",       "PCA", "pca_by_celltype_layer.svg")
plot_and_save_group("region",               "PCA", "pca_by_region.svg")
plot_and_save_group("region_only_layer",    "PCA", "pca_region_layer.svg")
plot_and_save_group("region_celltype_rep",  "PCA", "pca_region_celltype_replicate.svg")

# Export parsed metadata
save_table("tables/meta", "sample_metadata_parsed.csv", meta, row.names = TRUE)

# ================== Extensions ==================

# 1) Scree and cumulative variance
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
df_scree <- data.frame(PC = seq_along(var_explained),
                       Variance = var_explained,
                       Cumulative = cumsum(var_explained))
save_table("tables/variance", "pca_variance_explained.csv", df_scree)

p_scree <- ggplot(df_scree, aes(PC, Variance)) +
  geom_col(fill="#6B5B95") +
  geom_point() + geom_line(group=1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_pca_min() + labs(title="PCA Scree", x="Principal Component", y="Variance Explained")
save_plot("plots/variance", "pca_scree.svg", p_scree)

p_cum <- ggplot(df_scree, aes(PC, Cumulative)) +
  geom_point(color="#66C1A4") + geom_line(color="#66C1A4") +
  geom_hline(yintercept = 0.8, linetype="dashed", color="#B2B2B2") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,1)) +
  theme_pca_min() + labs(title="Cumulative Variance", x="Principal Component", y="Cumulative Fraction")
save_plot("plots/variance", "pca_cumulative_variance.svg", p_cum)

# 2) PC ~ metadata ANOVA with effect sizes
pc_df <- as.data.frame(pca$x)
pc_df$sample <- rownames(pc_df)
pc_meta <- cbind(pc_df, meta[rownames(pc_df), , drop=FALSE])
fac_vars <- c("region","layer","celltype","ReplicateGroup","plate","group","group2","ExpGroup")
fac_vars <- fac_vars[fac_vars %in% names(pc_meta)]

anova_rows <- list()
for (pc_name in colnames(pca$x)) {
  if (!length(fac_vars)) next
  form <- as.formula(paste(pc_name, "~", paste(fac_vars, collapse = " + ")))
  fit <- try(aov(form, data = pc_meta), silent = TRUE)
  if (inherits(fit, "try-error")) next
  a <- summary(fit)[[1]]
  ss_total <- sum(a[,"Sum Sq"], na.rm=TRUE)
  for (v in fac_vars) {
    if (v %in% rownames(a)) {
      ss <- a[v,"Sum Sq"]; df <- a[v,"Df"]; ms <- a[v,"Mean Sq"]; f <- a[v,"F value"]; p <- a[v,"Pr(>F)"]
      eta2 <- if (!is.na(ss_total) && ss_total>0) ss/ss_total else NA_real_
      anova_rows[[length(anova_rows)+1]] <- data.frame(PC=pc_name, term=v, df=df, ss=ss, ms=ms, F=f, p=p, eta2=eta2, stringsAsFactors = FALSE)
    }
  }
}
anova_tab <- if (length(anova_rows)) do.call(rbind, anova_rows) else data.frame()
if (nrow(anova_tab)) {
  anova_tab$q <- p.adjust(anova_tab$p, method = "BH")
  save_table("tables/associations", "pc_meta_anova.csv", anova_tab)
}

# 3) Top loadings export for PC1/PC2 (+/-) and enrichment-ready lists
rot <- as.data.frame(pca$rotation)
rot$protein <- rownames(pca$rotation)
top_n <- 50
for (k in c("PC1","PC2")) {
  if (!k %in% colnames(rot)) next
  ord_pos <- rot[order(-rot[[k]]), c("protein", k)]
  ord_neg <- rot[order(rot[[k]]),  c("protein", k)]
  save_table("tables/loadings", sprintf("loadings_%s_toppos.csv", k), head(ord_pos, top_n))
  save_table("tables/loadings", sprintf("loadings_%s_topneg.csv", k), head(ord_neg, top_n))
  ranked <- rot[, c("protein", k)]
  ranked <- ranked[order(-ranked[[k]]), ]
  save_table("tables/loadings", sprintf("loadings_%s_ranked_full.csv", k), ranked)
}

# 4) PC–protein correlations (per protein vs PC1/PC2)
pc_scores <- pca$x[, colnames(pca$x)[1:min(2, ncol(pca$x))], drop = FALSE]
cors <- lapply(colnames(pc_scores), function(pc) {
  s <- pc_scores[, pc]
  ct <- apply(mat, 1, function(v) {
    ok <- is.finite(v) & is.finite(s)
    if (sum(ok) < 3) return(c(r=NA_real_, p=NA_real_))
    test <- suppressWarnings(cor.test(v[ok], s[ok], method = "pearson"))
    c(r = unname(test$estimate), p = unname(test$p.value))
  })
  df <- as.data.frame(t(ct))
  df$protein <- rownames(df)
  df$pc <- pc
  df
})
cors_df <- do.call(rbind, cors)
cors_df$q <- ave(cors_df$p, cors_df$pc, FUN = function(x) p.adjust(x, method="BH"))
save_table("tables/correlations", "protein_pc_correlations.csv", cors_df)

# 5) Signed loading heatmap for PC1/PC2 (top 50 +/−)
heat_dir <- ensure_dir(subdir("plots/heatmaps"))
for (k in c("PC1","PC2")) {
  if (!k %in% colnames(rot)) next
  pos <- head(rot[order(-rot[[k]]), "protein"], top_n)
  neg <- head(rot[order(rot[[k]]),  "protein"], top_n)
  sel <- unique(c(pos, neg))
  sel <- sel[sel %in% rownames(mat)]
  if (length(sel) >= 4) {
    ord_samples <- order(pca$x[,k])
    M <- mat[sel, ord_samples, drop = FALSE]
    save_pheatmap(M, file.path(heat_dir, sprintf("heatmap_%s_topposneg.png", k)), width=8, height=10, scale="row")
  }
}

# 6) UMAP on PCs (first 20 PCs)
run_umap <- function(X, n_neighbors = 15, min_dist = 0.2, n_components = 2, metric = "euclidean"){
  if (!requireNamespace("uwot", quietly = TRUE)) {
    message("uwot not installed; skipping UMAP.")
    return(NULL)
  }
  set.seed(42)
  uwot::umap(X, n_neighbors = n_neighbors, min_dist = min_dist, n_components = n_components, metric = metric, verbose = FALSE)
}
npc <- min(20, ncol(pca$x))
um <- run_umap(pca$x[, 1:npc, drop = FALSE], n_neighbors = 15, min_dist = 0.2)
if (!is.null(um)) {
  um_df <- data.frame(UMAP1 = um[,1], UMAP2 = um[,2], meta[rownames(pca$x), , drop = FALSE])
  plot_umap_group <- function(key, out_file){
    if (!key %in% names(um_df)) return(NULL)
    grp <- factor(trim_ws(as.character(um_df[[key]])))
    pal <- make_modern_palette(nlevels(droplevels(grp)))
    p <- ggplot(um_df, aes(UMAP1, UMAP2, color = grp)) +
      geom_point(alpha = 0.8, size = 6) +
      scale_color_manual(values = pal, na.translate = FALSE) +
      theme_pca_min() + labs(title = paste("UMAP by", key), color = key)
    save_plot("plots/umap", out_file, p)
    invisible(p)
  }
  plot_umap_group("ReplicateGroup", "umap_by_replicate_group.svg")
  plot_umap_group("celltype",       "umap_by_celltype.svg")
  plot_umap_group("layer",          "umap_by_layer.svg")
  plot_umap_group("region",         "umap_by_region.svg")
  plot_umap_group("region_only_layer","umap_region_layer.svg")
}

# 7) Clustering on PCs and ARI/NMI vs metadata
pc_for_cluster <- scale(pca$x[, 1:npc, drop = FALSE])
best_k <- 2:8

if (!requireNamespace("cluster", quietly = TRUE)) {
  message("cluster not installed; silhouette may be skipped.")
}
if (!requireNamespace("aricode", quietly = TRUE)) {
  message("aricode not installed; ARI/NMI will be skipped.")
}

sil_tab <- list()
clu_eval <- list()
for (k in best_k) {
  km <- kmeans(pc_for_cluster, centers = k, nstart = 50, iter.max = 100)
  sil_avg <- NA_real_
  if (requireNamespace("cluster", quietly = TRUE)) {
    s <- cluster::silhouette(km$cluster, dist(pc_for_cluster))
    sil_avg <- mean(s[, "sil_width"], na.rm = TRUE)
  }
  sil_tab[[length(sil_tab)+1]] <- data.frame(k=k, silhouette=sil_avg)

  if (requireNamespace("aricode", quietly = TRUE)) {
    for (lab in c("region","layer","celltype","ReplicateGroup")) {
      if (!lab %in% names(meta)) next
      ref <- factor(meta[[lab]])
      pred <- factor(km$cluster, levels = sort(unique(km$cluster)))
      ari <- aricode::ARI(ref, pred)
      nmi <- aricode::NMI(ref, pred)
      clu_eval[[length(clu_eval)+1]] <- data.frame(k=k, label=lab, ARI=ari, NMI=nmi, stringsAsFactors = FALSE)
    }
  }
}
sil_df <- if (length(sil_tab)) do.call(rbind, sil_tab) else data.frame()
if (nrow(sil_df)) save_table("tables/clustering", "cluster_silhouette.csv", sil_df)
clu_df <- if (length(clu_eval)) do.call(rbind, clu_eval) else data.frame()
if (nrow(clu_df)) save_table("tables/clustering", "cluster_ARI_NMI.csv", clu_df)

suggest_k <- if (nrow(sil_df) && is.finite(max(sil_df$silhouette, na.rm=TRUE))) sil_df$k[which.max(sil_df$silhouette)] else 3
km_final <- kmeans(pc_for_cluster, centers = suggest_k, nstart = 100, iter.max = 200)
meta$kmeans_cluster <- factor(km_final$cluster)
pc_grp_plot <- plot_and_save_group("kmeans_cluster", "PCA", "pca_by_kmeans_cluster.svg")
if (inherits(pc_grp_plot, "gg")) save_plot("plots/clustering", "pca_by_kmeans_cluster.svg", pc_grp_plot)

# 8) Batch correction check: pre/post plate (or ReplicateGroup) using limma
run_remove_batch <- function(mat_samples_by_feature, meta_df, batch_key = "plate", design_keys = c("region","layer")) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    message("limma not installed; skipping batch correction.")
    return(NULL)
  }
  if (!batch_key %in% names(meta_df)) {
    message("Batch key not in meta: ", batch_key); return(NULL)
  }
  x <- mat_samples_by_feature
  smp <- colnames(x)
  meta_local <- meta_df[smp, , drop = FALSE]
  if (!identical(rownames(meta_local), smp)) stop("Meta/sample alignment failed for batch removal.")
  batch <- factor(meta_local[[batch_key]])
  if (length(batch) != ncol(x)) stop("Batch length != number of samples in matrix.")
  dk <- intersect(design_keys, names(meta_local))
  if (!length(dk)) {
    design <- matrix(1, ncol(x), 1); colnames(design) <- "Intercept"
  } else {
    design <- stats::model.matrix(as.formula(paste0("~ 0 + ", paste(dk, collapse = "+"))), data = meta_local)
  }
  if (nrow(design) != ncol(x)) stop("Row dimension of design must equal number of samples (columns of x).")
  limma::removeBatchEffect(x, batch = batch, design = design)
}

mat_adj <- try(run_remove_batch(mat_samples_by_feature = mat, meta_df = meta, batch_key = "plate", design_keys = c("region","layer")), silent = TRUE)
if (!inherits(mat_adj, "try-error") && !is.null(mat_adj)) {
  pca_prec <- pca
  pca_post <- prcomp(t(mat_adj), center = TRUE, scale. = TRUE)
  pre_plot <- plot_and_save_group("plate", "PCA (pre-correction)", "pca_by_plate_precorrect.svg")
  if (inherits(pre_plot, "gg")) save_plot("plots/batch", "pca_by_plate_precorrect.svg", pre_plot)
  pca <- pca_post
  post_plot <- plot_and_save_group("plate", "PCA (post-correction)", "pca_by_plate_postcorrect.svg")
  if (inherits(post_plot, "gg")) save_plot("plots/batch", "pca_by_plate_postcorrect.svg", post_plot)
  pca <- pca_prec
} else {
  message("Batch correction skipped or failed; check plate variable presence and design alignment.")
}

# 9) Optional biplot with sparse loadings on PC1/PC2
plot_biplot_sparse <- function(topN = 20, out_file = "pca_biplot_sparse.svg"){
  R <- as.data.frame(pca$rotation)
  R$protein <- rownames(R)
  for (k in c("PC1","PC2")) if (!k %in% names(R)) return(invisible(NULL))
  sel <- unique(c(head(R[order(-R$PC1), "protein"], topN),
                  head(R[order(R$PC1),  "protein"], topN),
                  head(R[order(-R$PC2), "protein"], topN),
                  head(R[order(R$PC2),  "protein"], topN)))
  dd <- data.frame(pca$x[,1:2, drop=FALSE], meta[rownames(pca$x), , drop=FALSE])
  pal <- make_modern_palette(nlevels(factor(dd$region)))
  p <- ggplot(dd, aes(PC1, PC2, color=factor(region))) +
    geom_point(size=3, alpha=0.8) +
    scale_color_manual(values = pal) +
    theme_pca_min() + labs(title="PCA biplot (sparse loadings)", color="region")
  arrows <- R[R$protein %in% sel, c("PC1","PC2","protein")]
  rngx <- diff(range(dd$PC1)); rngy <- diff(range(dd$PC2))
  scale_fac <- 0.5 * min(rngx, rngy)
  arrows$PC1s <- arrows$PC1 * scale_fac
  arrows$PC2s <- arrows$PC2 * scale_fac
  p <- p +
    geom_segment(data=arrows, aes(x=0, y=0, xend=PC1s, yend=PC2s), inherit.aes = FALSE, arrow = arrow(length = unit(0.15,"cm")), color="#444444", alpha=0.7) +
    ggrepel::geom_text_repel(data=arrows, aes(x=PC1s, y=PC2s, label=protein), inherit.aes = FALSE, size=3, max.overlaps = 200)
  save_plot("plots/biplot", out_file, p, width=8, height=6.5)
}
plot_biplot_sparse(20, "pca_biplot_sparse.svg")

# 10) Save session info
ensure_dir(subdir("tables/meta"))
writeLines(c(capture.output(sessionInfo())), con = file.path(subdir("tables/meta"), "sessionInfo.txt"))

# ================== Additional Extensions (organized outputs) ==================

# ----- Subdirectory helpers -----
subdir <- function(...){ file.path(output_dir, ...) }
ensure_dir <- function(path){ dir.create(path, showWarnings = FALSE, recursive = TRUE); path }

# Save plot: create subdir and call ggsave
save_plot <- function(subfolder, filename, plot, width=7.5, height=6.2, dpi=150){
  d <- ensure_dir(subdir(subfolder))
  ggsave(filename = filename, plot = plot, path = d, width = width, height = height, dpi = dpi)
  invisible(file.path(d, filename))
}

# Save table: create subdir and write via write_dt fallback
save_table <- function(subfolder, filename, df){
  d <- ensure_dir(subdir(subfolder))
  f <- file.path(d, filename)
  if (exists("write_dt")) write_dt(df, f) else utils::write.csv(df, f, row.names = FALSE)
  invisible(f)
}

# A) Confidence ellipses and centroids on PCA
add_pca_ellipses <- function(key, fname){
  grp <- build_group(key)
  keep <- !is.na(grp)
  grp2 <- droplevels(grp[keep])
  if (nlevels(grp2) < 2) return(invisible(NULL))
  df <- data.frame(pca$x[keep, 1:2, drop=FALSE], group = grp2, check.names = FALSE)
  varp <- (pca$sdev^2)/sum(pca$sdev^2)
  p <- ggplot(df, aes(PC1, PC2, color = group)) +
    geom_point(alpha=0.9, size=3) +
    stat_ellipse(type="t", level=0.95, linewidth=0.6) +
    stat_summary(aes(group=group, color=group), fun = mean, geom = "point", size=4, shape=4, stroke=1.2) +
    scale_color_manual(values = make_modern_palette(nlevels(grp2))) +
    theme_pca_min() +
    labs(title = paste("PCA with ellipses by", key),
         x = sprintf("PC1 (%.1f%%)", 100*varp[1]),
         y = sprintf("PC2 (%.1f%%)", 100*varp[2]),
         color = key)
  save_plot("plots/ellipses", fname, p)
}
add_pca_ellipses("region",   "pca_ellipses_region.svg")
add_pca_ellipses("layer",    "pca_ellipses_layer.svg")
add_pca_ellipses("celltype", "pca_ellipses_celltype.svg")

# B) Leave-one-batch/plate-out sensitivity (if plate exists)
leave_one_batch_pca <- function(batch_key = "plate"){
  if (!batch_key %in% names(meta)) { message("No batch key ", batch_key); return(invisible(NULL)) }
  res <- list()
  full_rot <- pca$rotation
  if (is.null(full_rot) || nrow(full_rot) == 0) { message("No rotation in base PCA; skipping similarity metrics."); full_rot <- NULL }
  batches <- na.omit(unique(meta[[batch_key]]))
  for (b in batches) {
    keep <- meta[[batch_key]] != b & !is.na(meta[[batch_key]])
    if (sum(keep) < 3) next
    Xsamp <- t(mat[, keep, drop=FALSE])
    pp <- try(prcomp(Xsamp, center=TRUE, scale.=TRUE), silent=TRUE)
    if (inherits(pp, "try-error")) next

    cs <- NA_real_
    if (!is.null(full_rot) && !is.null(pp$rotation)) {
      common <- intersect(rownames(full_rot), rownames(pp$rotation))
      if (length(common) >= 10) {
        F <- full_rot[common, 1:min(5, ncol(full_rot)), drop=FALSE]
        R <- pp$rotation[common, 1:ncol(F), drop=FALSE]
        cs <- sapply(seq_len(ncol(F)), function(i) {
          fi <- F[,i]; ri <- R[,i]
          den <- sqrt(sum(fi^2) * sum(ri^2))
          if (!is.finite(den) || den == 0) return(NA_real_)
          abs(sum(fi*ri) / den)
        })
      } else {
        want <- if (!is.null(full_rot)) min(5, ncol(full_rot)) else 5
        cs <- rep(NA_real_, want)
      }
    }
    var_ex <- (pp$sdev^2)/sum(pp$sdev^2)
    cs_vec <- if (length(cs) < length(var_ex)) c(cs, rep(NA_real_, length(var_ex)-length(cs))) else cs[seq_along(var_ex)]
    res[[length(res)+1]] <- data.frame(batch_left_out = b,
                                       PC = paste0("PC", seq_along(var_ex)),
                                       variance = var_ex,
                                       cos_sim = cs_vec)

    subsamp <- rownames(pp$x)
    subsamp <- intersect(subsamp, rownames(meta))
    df_sc <- data.frame(pp$x[subsamp,1:2, drop=FALSE],
                        region = droplevels(factor(meta[subsamp, "region"])),
                        check.names = FALSE)
    pal <- make_modern_palette(nlevels(droplevels(df_sc$region)))
    p <- ggplot(df_sc, aes(PC1, PC2, color=region)) +
         geom_point(size=3, alpha=0.8) +
         scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) +
         theme_pca_min() + labs(title=sprintf("PCA (-%s)", b))
    save_plot("plots/leave_one_batch", sprintf("pca_region_minus_%s.svg", b), p)
  }
  if (length(res)) {
    tab <- do.call(rbind, res)
    save_table("leave_one_batch", "pca_leave_one_batch_summary.csv", tab)
  }
}
leave_one_batch_pca("plate")

# C) Robust PCA variants
# Assumes run_irlba_pca() is already defined earlier and returns a prcomp-like list

# C1) IRLBA PCA (fast truncated SVD on scaled data), aligned rownames
pca_irlba <- run_irlba_pca()
if (!is.null(pca_irlba)) {
  if (is.null(rownames(pca_irlba$x))) rownames(pca_irlba$x) <- rownames(t(mat))
  samp <- rownames(pca_irlba$x)
  samp <- intersect(samp, rownames(meta))
  if (length(samp) >= 2) {
    df <- data.frame(pca_irlba$x[samp, 1:2, drop=FALSE],
                     region = droplevels(factor(meta[samp, "region"])),
                     check.names = FALSE)
    pal <- make_modern_palette(nlevels(df$region))
    p <- ggplot(df, aes(PC1, PC2, color = region)) +
      geom_point(size=3, alpha=0.8) +
      scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) +
      theme_pca_min() + labs(title="IRLBA PCA by region")
    save_plot("plots/irlba", "irlba_pca_by_region.svg", p)
  } else {
    message("IRLBA: no overlapping samples to plot with metadata.")
  }
}

# C2) ROBPCA (outlier-robust) if available
run_robpca <- function(){
  if (!requireNamespace("rospca", quietly = TRUE)) { message("rospca not installed; skipping ROBPCA."); return(NULL) }
  X <- scale(t(mat), center=TRUE, scale=TRUE)
  rownames_X <- rownames(X)
  rp <- rospca::robpca(X, k = min(10, ncol(X)))
  pcx <- rp$scores
  colnames(pcx) <- paste0("PC", seq_len(ncol(pcx)))
  rownames(pcx) <- rownames_X
  list(x = pcx, sdev = sqrt(rp$eigenvalues), rotation = NULL)
}
pca_ro <- run_robpca()
if (!is.null(pca_ro)) {
  samp <- rownames(pca_ro$x)
  samp <- intersect(samp, rownames(meta))
  df <- data.frame(pca_ro$x[samp,1:2, drop=FALSE],
                   region = droplevels(factor(meta[samp, "region"])),
                   check.names = FALSE)
  pal <- make_modern_palette(nlevels(droplevels(df$region)))
  p <- ggplot(df, aes(PC1, PC2, color = region)) + geom_point(size=3, alpha=0.8) +
       scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) + theme_pca_min() + labs(title="ROBPCA by region")
  save_plot("plots/robpca", "robpca_by_region.svg", p)
}

# D) t-SNE with parameter sweeps (on top PCs), aligned to meta and safe perplexity
run_tsne_and_plot <- function(){
  if (!requireNamespace("Rtsne", quietly = TRUE)) { message("Rtsne not installed; skipping t-SNE."); return(NULL) }
  samp <- rownames(pca$x)
  X <- pca$x[samp, 1:min(20, ncol(pca$x)), drop=FALSE]
  N <- nrow(X)
  if (N < 5) { message("Too few samples for t-SNE."); return(NULL) }
  for (perp in c(10, 30, 50)) {
    p_ok <- max(1, min(perp, floor((N - 1)/3)))
    set.seed(42 + perp)
    ts <- Rtsne::Rtsne(X, perplexity = p_ok, pca = FALSE, theta = 0.5, check_duplicates = FALSE, verbose = FALSE)
    emb <- ts$Y
    rownames(emb) <- samp
    df <- data.frame(tSNE1 = emb[,1], tSNE2 = emb[,2], meta[samp, , drop = FALSE], check.names = FALSE)
    plot_tsne <- function(key, out){
      if (!key %in% names(df)) return(NULL)
      grp <- factor(trim_ws(as.character(df[[key]])))
      pal <- make_modern_palette(nlevels(droplevels(grp)))
      p <- ggplot(df, aes(tSNE1, tSNE2, color = grp)) +
        geom_point(alpha=0.9, size=3) +
        scale_color_manual(values = pal, na.translate = FALSE) +
        theme_pca_min() + labs(title = paste("t-SNE by", key, "(perp", p_ok, ")"), color = key)
      save_plot("plots/tsne", out, p)
    }
    plot_tsne("region",   sprintf("tsne_region_perp%d.svg",  p_ok))
    plot_tsne("layer",    sprintf("tsne_layer_perp%d.svg",   p_ok))
    plot_tsne("celltype", sprintf("tsne_celltype_perp%d.svg", p_ok))
  }
}
run_tsne_and_plot()

# E) Gap statistic for clustering (PC space)
compute_gap <- function(){
  if (!requireNamespace("cluster", quietly = TRUE)) { message("cluster not installed; skipping gap."); return(NULL) }
  X <- scale(pca$x[, 1:min(10, ncol(pca$x)), drop=FALSE])
  set.seed(42)
  gap <- cluster::clusGap(X, FUN = stats::kmeans, K.max = 10, B = 50, nstart = 25, iter.max = 100)
  gap_df <- as.data.frame(gap$Tab)
  gap_df$k <- seq_len(nrow(gap_df))
  save_table("clustering", "cluster_gap_statistic.csv", gap_df)
  p <- ggplot(gap_df, aes(k, gap)) + geom_line() + geom_point() + theme_pca_min() + labs(title="Gap statistic", x="k", y="Gap")
  save_plot("clustering", "cluster_gap_statistic.svg", p)
}
compute_gap()

# F) PC–protein correlation volcano plots (uses existing correlations if available)
if (exists("cors_df") && nrow(cors_df)) {
  for (pc in unique(cors_df$pc)) {
    df <- cors_df[cors_df$pc == pc, , drop=FALSE]
    df$ml10q <- -log10(pmax(df$q, .Machine$double.xmin))
    p <- ggplot(df, aes(r, ml10q)) +
      geom_point(alpha=0.4, size=1.3, color="#6B5B95") +
      geom_vline(xintercept = c(-0.3, 0.3), linetype="dashed", color="#B2B2B2") +
      geom_hline(yintercept = -log10(0.05), linetype="dashed", color="#B2B2B2") +
      theme_pca_min() + labs(title = paste("Protein–", pc, "correlations"), x="Pearson r", y="-log10(q)")
    save_plot("correlations", sprintf("protein_pc_%s_volcano.svg", pc), p)
  }
}

# G) Stratified PCA: by region and by layer (local plotting)
stratified_pca <- function(by = c("region","layer")){
  by <- match.arg(by)
  lv <- levels(droplevels(factor(meta[[by]])))
  for (l in lv) {
    samp <- rownames(meta)[meta[[by]] == l & !is.na(meta[[by]])]
    if (length(samp) < 4) next
    samp <- intersect(samp, colnames(mat))
    if (length(samp) < 4) next
    pp <- prcomp(t(mat[, samp, drop=FALSE]), center=TRUE, scale.=TRUE)
    # by ReplicateGroup
    df1 <- data.frame(pp$x[,1:2, drop=FALSE],
                      ReplicateGroup = droplevels(factor(meta[samp, "ReplicateGroup"])),
                      check.names = FALSE)
    palR <- make_modern_palette(nlevels(droplevels(df1$ReplicateGroup)))
    p1 <- ggplot(df1, aes(PC1, PC2, color = ReplicateGroup)) +
      geom_point(size=3, alpha=0.85) +
      scale_color_manual(values=palR, drop=FALSE, na.translate=FALSE) +
      theme_pca_min() + labs(title = sprintf("PCA (%s=%s) by Replicate", by, l))
    save_plot(file.path("plots/stratified", by), sprintf("pca_by_replicate_%s_%s.svg", by, l), p1)
    # by celltype
    df2 <- data.frame(pp$x[,1:2, drop=FALSE],
                      celltype = droplevels(factor(meta[samp, "celltype"])),
                      check.names = FALSE)
    palC <- make_modern_palette(nlevels(droplevels(df2$celltype)))
    p2 <- ggplot(df2, aes(PC1, PC2, color = celltype)) +
      geom_point(size=3, alpha=0.85) +
      scale_color_manual(values=palC, drop=FALSE, na.translate=FALSE) +
      theme_pca_min() + labs(title = sprintf("PCA (%s=%s) by celltype", by, l))
    save_plot(file.path("plots/stratified", by), sprintf("pca_by_celltype_%s_%s.svg", by, l), p2)
  }
}
stratified_pca("region")
stratified_pca("layer")

# H) Contrast PCA: example contrast between two levels of a factor
contrast_pca <- function(factor_key = "layer", level_a = NULL, level_b = NULL, min_n = 3){
  if (!factor_key %in% names(meta)) return(invisible(NULL))
  fac <- droplevels(factor(meta[[factor_key]]))
  if (is.null(level_a) || is.null(level_b)) {
    if (nlevels(fac) < 2) return(invisible(NULL))
    levs <- levels(fac)[1:2]; level_a <- levs[1]; level_b <- levs[2]
  }
  A <- which(fac == level_a); B <- which(fac == level_b)
  if (length(A) < min_n || length(B) < min_n) return(invisible(NULL))
  M <- mat[, c(A, B), drop=FALSE]
  g <- c(rep(level_a, length(A)), rep(level_b, length(B)))
  Mz <- sweep(M, 1, rowMeans(M, na.rm=TRUE), "-")
  muA <- rowMeans(Mz[, g==level_a, drop=FALSE], na.rm=TRUE)
  muB <- rowMeans(Mz[, g==level_b, drop=FALSE], na.rm=TRUE)
  save_table("contrast", sprintf("contrast_%s_%s_vs_%s.csv", factor_key, level_a, level_b),
             data.frame(protein = rownames(M), muA = muA, muB = muB, diff = muA - muB, check.names = FALSE))
  pp <- prcomp(t(M), center=TRUE, scale.=TRUE)
  df <- data.frame(pp$x[,1:2, drop=FALSE], grp = droplevels(factor(g)), check.names = FALSE)
  pal <- make_modern_palette(nlevels(droplevels(df$grp)))
  p <- ggplot(df, aes(PC1, PC2, color = grp)) + geom_point(size=3, alpha=0.85) +
       scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) + theme_pca_min() +
       labs(title = sprintf("PCA (%s: %s vs %s)", factor_key, level_a, level_b), color = factor_key)
  save_plot("plots/contrast", sprintf("pca_%s_%s_vs_%s.svg", factor_key, level_a, level_b), p)
}
if ("layer" %in% names(meta) && nlevels(factor(meta$layer)) >= 2) {
  levs <- levels(droplevels(factor(meta$layer)))
  contrast_pca("layer", levs[1], levs[2])
}

# I) Optional supervised projections (LDA, PLS-DA) – simple CV metrics
run_lda <- function(label_key = "region"){
  if (!label_key %in% names(meta)) return(NULL)
  if (!requireNamespace("MASS", quietly = TRUE)) { message("MASS not installed; skipping LDA."); return(NULL) }
  y <- factor(meta[[label_key]])
  X <- pca$x[, 1:min(10, ncol(pca$x)), drop=FALSE]
  keep <- !is.na(y) & rowSums(!is.finite(X)) == 0
  y <- droplevels(y[keep]); X <- X[keep, , drop=FALSE]
  if (nlevels(y) < 2) return(NULL)
  fit <- MASS::lda(X, grouping = y, CV = TRUE)
  acc <- mean(fit$class == y)
  save_table("supervised", sprintf("lda_%s_cv_accuracy.csv", label_key), data.frame(metric="CV-accuracy", value=acc))
}
run_lda("region")
run_lda("layer")

run_plsda <- function(label_key = "region"){
  if (!label_key %in% names(meta)) return(NULL)
  if (!requireNamespace("mixOmics", quietly = TRUE)) { message("mixOmics not installed; skipping PLS-DA."); return(NULL) }
  y <- factor(meta[[label_key]])
  X <- t(mat)
  keep <- !is.na(y) & rowSums(!is.finite(X)) == 0
  y <- droplevels(y[keep]); X <- X[keep, , drop=FALSE]
  if (nlevels(y) < 2) return(NULL)
  set.seed(42)
  fit <- mixOmics::plsda(X, y, ncomp = min(3, nlevels(y)-1))
  perf <- mixOmics::perf(fit, validation="Mfold", folds=5, nrepeat=5, progressBar=FALSE)
  acc <- mean(perf$error.rate$overall)
  save_table("supervised", sprintf("plsda_%s_cv_accuracy.csv", label_key), data.frame(metric="Mfold-accuracy", value=acc))
}
run_plsda("region")
run_plsda("layer")

# J) WGCNA module eigengene vs PCs (explicit WGCNA::cor, aligned samples)
wgcna_modules_vs_pcs <- function(){
  if (!requireNamespace("WGCNA", quietly = TRUE)) { message("WGCNA not installed; skipping WGCNA-PC association."); return(NULL) }
  WGCNA::enableWGCNAThreads()
  datExpr <- t(mat) # samples x proteins
  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose=0)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  powers <- 6
  suppressWarnings( WGCNA::pickSoftThreshold(datExpr, powerVector=powers, verbose=0) )
  net <- WGCNA::blockwiseModules(datExpr, power=powers, TOMType="unsigned", minModuleSize=30,
                                 reassignThreshold=0, mergeCutHeight=0.25, numericLabels=TRUE,
                                 pamRespectsDendro=FALSE, saveTOMs=FALSE, verbose=0)
  MEs <- WGCNA::moduleEigengenes(datExpr, colors=net$colors)$eigengenes
  samp <- rownames(datExpr)
  PCs <- pca$x[samp, 1:min(10, ncol(pca$x)), drop=FALSE]
  cc <- WGCNA::cor(MEs, PCs, use="pairwise.complete.obs", method="pearson")
  ps <- WGCNA::corPvalueStudent(cc, nrow(datExpr))
  df <- data.frame(expand.grid(ME=rownames(cc), PC=colnames(cc)),
                   r = as.vector(cc), p = as.vector(ps))
  df$q <- p.adjust(df$p, method="BH")
  save_table("wgcna", "wgcna_module_vs_pc.csv", df)
}
wgcna_modules_vs_pcs()