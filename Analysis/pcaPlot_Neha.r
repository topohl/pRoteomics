# ================== PCA analysis and plotting (extended, fixed, organized) ==================
# Set a reliable CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
}

# Load core packages
pacman::p_load(
    data.table, ggplot2, factoextra, reshape2, stats, ggrepel, tools,
    grid, uwot, RColorBrewer, pheatmap, Rtsne, aricode, rospca, irlba,
    pandoc, treemapify, dplyr
)

# Try to load aricode separately (it's only used for clustering metrics)
if (!requireNamespace("aricode", quietly = TRUE)) {
    message("aricode package not available - clustering ARI/NMI metrics will be skipped")
} else {
    library(aricode)
}

set.seed(42)

# =============== Config =================
#gct_file   <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/E9_pg_matrix_protigy.gct"
#output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/pca_plots"
gct_file <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Datasets/gct/data/pg.matrix_filtered_pcaAdjusted.gct"
output_dir <- "S:/Lab_Member/Tobi/Experiments/Collabs/Neha/clusterProfiler/Results/pca_plots"
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
read_gct <- function(path) {
    cat("Reading GCT file:", path, "\n")
    
    all_lines <- readLines(path)
    
    # Line 1: version marker
    version <- trimws(all_lines[1])
    if (version != "#1.3") {
        warning("Expected '#1.3', found: ", version)
    }
    
    # Line 2: dimensions
    dims_parts <- strsplit(all_lines[2], "\t")[[1]]
    dims <- as.numeric(dims_parts)
    
    # Line 3: sample IDs (first element is "id", skip it)
    line3_parts <- strsplit(all_lines[3], "\t")[[1]]
    sample_ids <- line3_parts[-1]
    
    # Lines 4-12: metadata (9 rows total)
    # sampleNumber, shortname, plate, group2, AnimalID, ReplicateGroup, celltype, ExpGroup, celltype_group
    meta_lines <- all_lines[4:12]
    meta_split <- lapply(meta_lines, function(x) strsplit(x, "\t")[[1]])
    
    meta_keys <- sapply(meta_split, `[`, 1)
    meta_values_list <- lapply(meta_split, function(x) x[-1])
    
    # Build metadata dataframe
    meta_mat <- do.call(rbind, meta_values_list)
    meta_df <- as.data.frame(t(meta_mat), stringsAsFactors = FALSE)
    colnames(meta_df) <- meta_keys
    rownames(meta_df) <- sample_ids
    
    # Convert numeric columns
    for (col in colnames(meta_df)) {
        suppressWarnings({
            x_num <- as.numeric(meta_df[[col]])
            if (!any(is.na(x_num))) meta_df[[col]] <- x_num
        })
    }
    
    # Lines 13+: expression data
    expr_lines <- all_lines[13:length(all_lines)]
    expr_split <- lapply(expr_lines, function(x) strsplit(x, "\t")[[1]])
    
    proteins <- sapply(expr_split, `[`, 1)
    expr_values_list <- lapply(expr_split, function(x) as.numeric(x[-1]))
    
    expr_mat <- do.call(rbind, expr_values_list)
    rownames(expr_mat) <- proteins
    colnames(expr_mat) <- sample_ids
    
    # Filter rows: keep only those containing "_MOUSE"
    mouse_rows <- grepl("_MOUSE", rownames(expr_mat))
    expr_mat <- expr_mat[mouse_rows, , drop = FALSE]
    
    # Keep only first protein name if multiple names separated by ";"
    rownames(expr_mat) <- sapply(strsplit(rownames(expr_mat), ";"), `[`, 1)
    
    # Remove "_MOUSE" suffix from rownames
    rownames(expr_mat) <- sub("_MOUSE$", "", rownames(expr_mat))
    
    cat("Success!\n")
    cat("Expression matrix:", nrow(expr_mat), "proteins x", ncol(expr_mat), "samples\n")
    cat("Metadata:", nrow(meta_df), "samples x", ncol(meta_df), "fields\n")
    cat("Metadata fields:", paste(colnames(meta_df), collapse = ", "), "\n")
    
    list(mat = expr_mat, meta = meta_df)
}


# ================== Build mat/meta with checks ==================
g <- read_gct(gct_file)
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

# Derived labels based on actual metadata fields
# Available fields: sampleNumber, shortname, plate, group2, AnimalID, ReplicateGroup, celltype, ExpGroup, celltype_group
if ("celltype" %in% names(meta) && "group2" %in% names(meta)) {
    meta$celltype_group2 <- paste(meta$celltype, meta$group2, sep = "_")
}
if ("celltype" %in% names(meta) && "ExpGroup" %in% names(meta)) {
    meta$celltype_ExpGroup <- paste(meta$celltype, meta$ExpGroup, sep = "_")
}
if ("ReplicateGroup" %in% names(meta) && "celltype" %in% names(meta)) {
    meta$ReplicateGroup_celltype <- paste(meta$ReplicateGroup, meta$celltype, sep = "_")
}

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

#base_hex <- c("#F8D247", "#E89369", "#B2B2B2", "#C7D745", "#C6B18B",
#              "#C592C5", "#A89BB0", "#DEC196", "#66C1A4", "#FF6F61",
#              "#6B5B95", "#88B04B", "#F7CAC9", "#92A8D1", "#955251")

#base_hex <- c("#FFB6C1", "#EE82EE", "#FF69B4", "#FF1493", "#C71585",
#              "#E6B0E8", "#DDA0DD", "#ff68c0ff", "#DA70D6", "#BA55D3",
#              "#F4C2C2", "#FFD1DC", "#FFDAB9", "#FFE4E1", "#FFF0F5")

base_hex <- c("#1E90FF", "#61d7ffff", "#FFA500", "#FFD700", "#FF4500",
              "#FF6347", "#FFDAB9", "#FFE4B5", "#F0E68C", "#ADFF2F",
              "#98FB98", "#FF69B4", "#FFB6C1", "#FF1493", "#FF4500")

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

# 2) PC ~ metadata ANOVA with effect sizes (celltype and ExpGroup), robust and always writes a file
pc_df <- as.data.frame(pca$x)
pc_df$sample <- rownames(pc_df)
pc_meta <- cbind(pc_df[, "sample", drop = FALSE], 
                 pc_df[, grep("^PC", names(pc_df))], 
                 meta[rownames(pc_df), , drop = FALSE])

# Candidate factors
cand_vars <- intersect(c("celltype","ExpGroup"), names(pc_meta))
if (length(cand_vars)) {
    pc_meta[cand_vars] <- lapply(pc_meta[cand_vars], function(x) {
        x <- trim_ws(as.character(x)); x[x == ""] <- NA; factor(x)
    })
}

anova_rows <- list()
for (pc_name in colnames(pca$x)) {
    if (!length(cand_vars)) break
    dat <- pc_meta[, c(pc_name, cand_vars), drop = FALSE]
    dat <- dat[complete.cases(dat), , drop = FALSE]
    if (nrow(dat) < 3) next
    vars_ok <- cand_vars[sapply(cand_vars, function(v) nlevels(droplevels(dat[[v]])) >= 2)]
    if (!length(vars_ok)) next
    form <- as.formula(sprintf("%s ~ %s", pc_name, paste(vars_ok, collapse = " + ")))
    fit <- try(aov(form, data = dat), silent = TRUE)
    if (inherits(fit, "try-error")) next
    a <- try(summary(fit)[[1]], silent = TRUE)
    if (inherits(a, "try-error")) next
    ss_total <- sum(a[, "Sum Sq"], na.rm = TRUE)
    for (v in vars_ok) {
        if (v %in% rownames(a)) {
            ss <- a[v, "Sum Sq"]; df <- a[v, "Df"]; ms <- a[v, "Mean Sq"]; f <- a[v, "F value"]; p <- a[v, "Pr(>F)"]
            eta2 <- if (is.finite(ss_total) && ss_total > 0) ss / ss_total else NA_real_
            anova_rows[[length(anova_rows)+1L]] <- data.frame(
                PC = pc_name, term = v, df = as.numeric(df), ss = as.numeric(ss),
                ms = as.numeric(ms), F = as.numeric(f), p = as.numeric(p), eta2 = as.numeric(eta2),
                stringsAsFactors = FALSE
            )
        }
    }
}
anova_tab <- if (length(anova_rows)) do.call(rbind, anova_rows) else
    data.frame(PC=character(), term=character(), df=double(), ss=double(), ms=double(), F=double(), p=double(), eta2=double(), stringsAsFactors = FALSE)

if (nrow(anova_tab)) anova_tab$q <- p.adjust(anova_tab$p, method = "BH") else anova_tab$q <- numeric(0)

# Always write an output file (even if empty)
save_table("tables/associations", "pc_meta_anova.csv", anova_tab)

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
      geom_point(shape = 16, alpha = 0.8, size = 6, stroke = 0) +
      scale_color_manual(values = pal, na.translate = FALSE) +
      theme_pca_min() + labs(title = paste("UMAP by", key), color = key)
    save_plot("plots/umap", out_file, p)
    invisible(p)
  }
  plot_umap_group("ReplicateGroup", "umap_by_replicate_group.svg")
  plot_umap_group("celltype",       "umap_by_celltype.svg")
}
# detect switched samples / outliers based on umap
if (!is.null(um)) {
    # Calculate celltype group centers in UMAP space using median
    celltype_centers <- aggregate(um, by = list(celltype = meta[rownames(um), "celltype"]), FUN = median)
    rownames(celltype_centers) <- celltype_centers$celltype
    celltype_centers$celltype <- NULL
    
    # Calculate distance of each sample to its own celltype center
    sample_distances <- sapply(rownames(um), function(sample) {
        sample_celltype <- meta[sample, "celltype"]
        if (is.na(sample_celltype)) return(NA_real_)
        
        center <- as.numeric(celltype_centers[sample_celltype, ])
        sample_coords <- as.numeric(um[sample, ])
        
        # Euclidean distance to center
        sqrt(sum((sample_coords - center)^2))
    })
    
    # Identify outliers per celltype (samples far from their own group center)
    outlier_samples <- character(0)
    outlier_info <- list()
    
    for (ct in unique(meta$celltype)) {
        if (is.na(ct)) next
        
        samples_in_celltype <- rownames(meta)[meta$celltype == ct & !is.na(meta$celltype)]
        samples_in_celltype <- intersect(samples_in_celltype, names(sample_distances))
        
        if (length(samples_in_celltype) < 3) next  # Need at least 3 samples to detect outliers
        
        dists <- sample_distances[samples_in_celltype]
        
        # Use 95th percentile or 2 MADs as threshold
        threshold <- max(quantile(dists, 0.95, na.rm = TRUE), 
                        median(dists, na.rm = TRUE) + 2 * mad(dists, na.rm = TRUE))
        
        outliers_in_group <- samples_in_celltype[dists > threshold]
        
        if (length(outliers_in_group) > 0) {
            outlier_samples <- c(outlier_samples, outliers_in_group)
            for (outlier in outliers_in_group) {
                outlier_info[[outlier]] <- data.frame(
                    sample = outlier,
                    celltype = ct,
                    distance_to_center = dists[outlier],
                    threshold = threshold,
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    
    if (length(outlier_samples) > 0) {
        message("Potential outlier samples detected based on distance to celltype center:")
        print(outlier_samples)
        
        # Find potential switch partners for each outlier
        switch_partners <- list()
        for (outlier in outlier_samples) {
            outlier_celltype <- meta[outlier, "celltype"]
            outlier_coords <- um[outlier, ]
            
            # Find closest celltype center (excluding own celltype)
            other_celltypes <- setdiff(rownames(celltype_centers), outlier_celltype)
            if (length(other_celltypes) > 0) {
                distances_to_centers <- apply(celltype_centers[other_celltypes, , drop = FALSE], 1, function(center) {
                    sqrt(sum((outlier_coords - center)^2))
                })
                closest_celltype <- names(which.min(distances_to_centers))
                
                switch_partners[[outlier]] <- data.frame(
                    outlier = outlier,
                    outlier_celltype = outlier_celltype,
                    distance_to_own_center = sample_distances[outlier],
                    potential_switch_celltype = closest_celltype,
                    distance_to_switch_center = min(distances_to_centers),
                    stringsAsFactors = FALSE
                )
            }
        }
        
        if (length(switch_partners) > 0) {
            switch_df <- do.call(rbind, switch_partners)
            message("\nPotential sample switch partners:")
            print(switch_df)
            save_table("tables/umap", "umap_outlier_samples_with_switches.csv", switch_df)
        }
        
        # Save detailed outlier info
        outlier_detail_df <- do.call(rbind, outlier_info)
        save_table("tables/umap", "umap_outlier_samples_details.csv", outlier_detail_df)
        
        # Label outliers in UMAP plot
        um_df$outlier <- ifelse(rownames(um_df) %in% outlier_samples, "Outlier", "Normal")
        p <- ggplot(um_df, aes(UMAP1, UMAP2, color = outlier)) +
            geom_point(alpha = 0.8, size = 6) +
            scale_color_manual(values = c("Normal" = "grey", "Outlier" = "red"), na.translate = FALSE) +
            theme_pca_min() + labs(title = "UMAP with Outliers Labeled", color = "Sample Status")
        save_plot("plots/umap", "umap_with_outliers.svg", p)
        
        # Additional plot: show distances to centers
        dist_df <- data.frame(
            sample = names(sample_distances),
            celltype = meta[names(sample_distances), "celltype"],
            distance = sample_distances,
            is_outlier = names(sample_distances) %in% outlier_samples,
            stringsAsFactors = FALSE
        )
        p_dist <- ggplot(dist_df, aes(x = celltype, y = distance, color = is_outlier)) +
            geom_jitter(width = 0.2, alpha = 0.7, size = 3) +
            scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), labels = c("Normal", "Outlier")) +
            theme_pca_min() +
            labs(title = "Distance to Celltype Center (Median)", x = "Celltype", y = "UMAP Distance", color = "Status") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        save_plot("plots/umap", "umap_distance_to_centers.svg", p_dist)
        
    } else {
        message("No outlier samples detected based on distance to celltype centers.")
    }
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
    for (lab in c("celltype","ReplicateGroup")) {
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
run_remove_batch <- function(mat_samples_by_feature, meta_df, batch_key = "plate", design_keys = c("celltype")) {
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

mat_adj <- try(run_remove_batch(mat_samples_by_feature = mat, meta_df = meta, batch_key = "plate", design_keys = c("celltype")), silent = TRUE)
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
  pal <- make_modern_palette(nlevels(factor(dd$celltype)))
  p <- ggplot(dd, aes(PC1, PC2, color=factor(celltype))) +
    geom_point(size=3, alpha=0.8) +
    scale_color_manual(values = pal) +
    theme_pca_min() + labs(title="PCA biplot (sparse loadings)", color="celltype")
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
                        celltype = droplevels(factor(meta[subsamp, "celltype"])),
                        check.names = FALSE)
    pal <- make_modern_palette(nlevels(droplevels(df_sc$celltype)))
    p <- ggplot(df_sc, aes(PC1, PC2, color=celltype)) +
         geom_point(size=3, alpha=0.8) +
         scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) +
         theme_pca_min() + labs(title=sprintf("PCA (-%s)", b))
    save_plot("plots/leave_one_batch", sprintf("pca_celltype_minus_%s.svg", b), p)
  }
  if (length(res)) {
    tab <- do.call(rbind, res)
    save_table("leave_one_batch", "pca_leave_one_batch_summary.csv", tab)
  }
}
leave_one_batch_pca("plate")

# C) Robust PCA variants
# Assumes run_irlba_pca() is already defined earlier and returns a prcomp-like list
run_irlba_pca <- function(){
  if (!requireNamespace("irlba", quietly = TRUE)) { message("irlba not installed; skipping IRLBA PCA."); return(NULL) }
  X <- scale(t(mat), center=TRUE, scale=TRUE)
  rownames_X <- rownames(X)
  k <- min(10, ncol(X)-1, nrow(X)-1)
  if (k < 2) { message("Not enough samples/features for IRLBA PCA."); return(NULL) }
  irlba_res <- irlba::irlba(X, nv = k, nu = k)
  pcx <- irlba_res$u %*% diag(irlba_res$d)
  colnames(pcx) <- paste0("PC", seq_len(ncol(pcx)))
  rownames(pcx) <- rownames_X
  list(x = pcx, sdev = irlba_res$d / sqrt(max(1, nrow(X)-1)), rotation = NULL)
}


# C1) IRLBA PCA (fast truncated SVD on scaled data), aligned rownames
pca_irlba <- run_irlba_pca()
if (!is.null(pca_irlba)) {
  if (is.null(rownames(pca_irlba$x))) rownames(pca_irlba$x) <- rownames(t(mat))
  samp <- rownames(pca_irlba$x)
  samp <- intersect(samp, rownames(meta))
  if (length(samp) >= 2) {
    df <- data.frame(pca_irlba$x[samp, 1:2, drop=FALSE],
                     celltype = droplevels(factor(meta[samp, "celltype"])),
                     check.names = FALSE)
    pal <- make_modern_palette(nlevels(droplevels(df$celltype)))
    p <- ggplot(df, aes(PC1, PC2, color = celltype)) +
      geom_point(size=3, alpha=0.8) +
      scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) +
      theme_pca_min() + labs(title="IRLBA PCA by celltype", color="celltype")
    save_plot("plots/irlba", "irlba_pca_by_celltype.svg", p)
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
                   celltype = droplevels(factor(meta[samp, "celltype"])),
                   check.names = FALSE)
  pal <- make_modern_palette(nlevels(droplevels(df$celltype)))
  p <- ggplot(df, aes(PC1, PC2, color = celltype)) +
       geom_point(size=3, alpha=0.8) +
       scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) +
       theme_pca_min() + labs(title="ROBPCA by celltype", color="celltype")
  save_plot("plots/robpca", "robpca_by_celltype.svg", p)
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
#stratified_pca <- function(by = c("region","layer")){
#  by <- match.arg(by)
#  lv <- levels(droplevels(factor(meta[[by]])))
#  for (l in lv) {
#    samp <- rownames(meta)[meta[[by]] == l & !is.na(meta[[by]])]
#    if (length(samp) < 4) next
#    samp <- intersect(samp, colnames(mat))
#    if (length(samp) < 4) next
#    pp <- prcomp(t(mat[, samp, drop=FALSE]), center=TRUE, scale.=TRUE)
#    # by ReplicateGroup
#    df1 <- data.frame(pp$x[,1:2, drop=FALSE],
#                      ReplicateGroup = droplevels(factor(meta[samp, "ReplicateGroup"])),
#                      check.names = FALSE)
#    palR <- make_modern_palette(nlevels(droplevels(df1$ReplicateGroup)))
#    p1 <- ggplot(df1, aes(PC1, PC2, color = ReplicateGroup)) +
#      geom_point(size=3, alpha=0.85) +
#      scale_color_manual(values=palR, drop=FALSE, na.translate=FALSE) +
#      theme_pca_min() + labs(title = sprintf("PCA (%s=%s) by Replicate", by, l))
#    save_plot(file.path("plots/stratified", by), sprintf("pca_by_replicate_%s_%s.svg", by, l), p1)
#    # by celltype
#    df2 <- data.frame(pp$x[,1:2, drop=FALSE],
#                      celltype = droplevels(factor(meta[samp, "celltype"])),
#                      check.names = FALSE)
#    palC <- make_modern_palette(nlevels(droplevels(df2$celltype)))
#    p2 <- ggplot(df2, aes(PC1, PC2, color = celltype)) +
#      geom_point(size=3, alpha=0.85) +
#      scale_color_manual(values=palC, drop=FALSE, na.translate=FALSE) +
#      theme_pca_min() + labs(title = sprintf("PCA (%s=%s) by celltype", by, l))
#    save_plot(file.path("plots/stratified", by), sprintf("pca_by_celltype_%s_%s.svg", by, l), p2)
#  }
#}
#stratified_pca("region")
#stratified_pca("layer")

# H) Contrast PCA: example contrast between two levels of a factor
#contrast_pca <- function(factor_key = "layer", level_a = NULL, level_b = NULL, min_n = 3){
#  if (!factor_key %in% names(meta)) return(invisible(NULL))
#  fac <- droplevels(factor(meta[[factor_key]]))
#  if (is.null(level_a) || is.null(level_b)) {
#    if (nlevels(fac) < 2) return(invisible(NULL))
#    levs <- levels(fac)[1:2]; level_a <- levs[1]; level_b <- levs[2]
#  }
#  A <- which(fac == level_a); B <- which(fac == level_b)
#  if (length(A) < min_n || length(B) < min_n) return(invisible(NULL))
#  M <- mat[, c(A, B), drop=FALSE]
#  g <- c(rep(level_a, length(A)), rep(level_b, length(B)))
#  Mz <- sweep(M, 1, rowMeans(M, na.rm=TRUE), "-")
#  muA <- rowMeans(Mz[, g==level_a, drop=FALSE], na.rm=TRUE)
#  muB <- rowMeans(Mz[, g==level_b, drop=FALSE], na.rm=TRUE)
#  save_table("contrast", sprintf("contrast_%s_%s_vs_%s.csv", factor_key, level_a, level_b),
#             data.frame(protein = rownames(M), muA = muA, muB = muB, diff = muA - muB, check.names = FALSE))
#  pp <- prcomp(t(M), center=TRUE, scale.=TRUE)
#  df <- data.frame(pp$x[,1:2, drop=FALSE], grp = droplevels(factor(g)), check.names = FALSE)
#  pal <- make_modern_palette(nlevels(droplevels(df$grp)))
#  p <- ggplot(df, aes(PC1, PC2, color = grp)) + geom_point(size=3, alpha=0.85) +
#       scale_color_manual(values=pal, drop=FALSE, na.translate=FALSE) + theme_pca_min() +
#       labs(title = sprintf("PCA (%s: %s vs %s)", factor_key, level_a, level_b), color = factor_key)
#  save_plot("plots/contrast", sprintf("pca_%s_%s_vs_%s.svg", factor_key, level_a, level_b), p)
#}
#if ("layer" %in% names(meta) && nlevels(factor(meta$layer)) >= 2) {
#  levs <- levels(droplevels(factor(meta$layer)))
#  contrast_pca("layer", levs[1], levs[2])
#}

# I) Optional supervised projections (LDA, PLS-DA) – simple CV metrics
#run_lda <- function(label_key = "region"){
#  if (!label_key %in% names(meta)) return(NULL)
#  if (!requireNamespace("MASS", quietly = TRUE)) { message("MASS not installed; skipping LDA."); return(NULL) }
#  y <- factor(meta[[label_key]])
#  X <- pca$x[, 1:min(10, ncol(pca$x)), drop=FALSE]
#  keep <- !is.na(y) & rowSums(!is.finite(X)) == 0
#  y <- droplevels(y[keep]); X <- X[keep, , drop=FALSE]
#  if (nlevels(y) < 2) return(NULL)
#  fit <- MASS::lda(X, grouping = y, CV = TRUE)
#  acc <- mean(fit$class == y)
#  save_table("supervised", sprintf("lda_%s_cv_accuracy.csv", label_key), data.frame(metric="CV-accuracy", value=acc))
#}
#run_lda("region")
#run_lda("layer")

#run_plsda <- function(label_key = "region"){
#  if (!label_key %in% names(meta)) return(NULL)
#  if (!requireNamespace("mixOmics", quietly = TRUE)) { message("mixOmics not installed; skipping PLS-DA."); return(NULL) }
#  y <- factor(meta[[label_key]])
#  X <- t(mat)
#  keep <- !is.na(y) & rowSums(!is.finite(X)) == 0
#  y <- droplevels(y[keep]); X <- X[keep, , drop=FALSE]
#  if (nlevels(y) < 2) return(NULL)
#  set.seed(42)
#  fit <- mixOmics::plsda(X, y, ncomp = min(3, nlevels(y)-1))
#  perf <- mixOmics::perf(fit, validation="Mfold", folds=5, nrepeat=5, progressBar=FALSE)
#  acc <- mean(perf$error.rate$overall)
#  save_table("supervised", sprintf("plsda_%s_cv_accuracy.csv", label_key), data.frame(metric="Mfold-accuracy", value=acc))
#}
#run_plsda("region")
#run_plsda("layer")

# J) WGCNA module eigengene vs PCs (explicit WGCNA::cor, aligned samples)
#wgcna_modules_vs_pcs <- function(){
#  if (!requireNamespace("WGCNA", quietly = TRUE)) { message("WGCNA not installed; skipping WGCNA-PC association."); return(NULL) }
#  WGCNA::enableWGCNAThreads()
#  datExpr <- t(mat) # samples x proteins
#  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose=0)
#  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
#  powers <- 6
#  suppressWarnings( WGCNA::pickSoftThreshold(datExpr, powerVector=powers, verbose=0) )
#  net <- WGCNA::blockwiseModules(datExpr, power=powers, TOMType="unsigned", minModuleSize=30,
#                                 reassignThreshold=0, mergeCutHeight=0.25, numericLabels=TRUE,
#                                 pamRespectsDendro=FALSE, saveTOMs=FALSE, verbose=0)
#  MEs <- WGCNA::moduleEigengenes(datExpr, colors=net$colors)$eigengenes
#  samp <- rownames(datExpr)
#  PCs <- pca$x[samp, 1:min(10, ncol(pca$x)), drop=FALSE]
#  cc <- WGCNA::cor(MEs, PCs, use="pairwise.complete.obs", method="pearson")
#  ps <- WGCNA::corPvalueStudent(cc, nrow(datExpr))
#  df <- data.frame(expand.grid(ME=rownames(cc), PC=colnames(cc)),
#                   r = as.vector(cc), p = as.vector(ps))
#  df$q <- p.adjust(df$p, method="BH")
#  save_table("wgcna", "wgcna_module_vs_pc.csv", df)
#}
#wgcna_modules_vs_pcs()

# ================== Additional Supplementary Figure Plots ==================

# K) PC contribution heatmap: show which PCs each sample loads on
pc_contribution_heatmap <- function(n_pcs = 10) {
    pc_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    pc_mat_scaled <- t(scale(t(pc_mat)))  # Scale per sample
    
    # Use sampleNumber as row names if available
    if ("sampleNumber" %in% names(meta)) {
        rownames(pc_mat_scaled) <- meta[rownames(pc_mat_scaled), "sampleNumber"]
    }
    
    # Annotate samples with metadata
    ann_cols <- intersect(c("celltype", "ExpGroup", "ReplicateGroup", "plate"), names(meta))
    if (length(ann_cols) > 0) {
        annotation_row <- meta[rownames(pca$x), ann_cols, drop = FALSE]
        if ("sampleNumber" %in% names(meta)) {
            rownames(annotation_row) <- meta[rownames(pca$x), "sampleNumber"]
        }
        save_pheatmap(pc_mat_scaled, 
                     file.path(ensure_dir(subdir("plots/heatmaps")), "pc_contributions_per_sample.svg"),
                     width = 3, height = 13, scale = "none")
    }
}
pc_contribution_heatmap(10)

# L) Pairwise PC scatterplot matrix for PC1-PC5
pc_pairs_plot <- function(n_pcs = 5) {
    if (!requireNamespace("GGally", quietly = TRUE)) {
        message("GGally not installed; skipping pairs plot")
        return(NULL)
    }
    pc_subset <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    df <- data.frame(pc_subset, celltype = meta[rownames(pc_subset), "celltype"])
    pal <- make_modern_palette(nlevels(droplevels(factor(df$celltype))))
    
    p <- GGally::ggpairs(df, columns = 1:ncol(pc_subset), 
                        aes(color = celltype, alpha = 0.6),
                        upper = list(continuous = "cor"),
                        lower = list(continuous = "points")) +
        scale_color_manual(values = pal) +
        theme_pca_min()
    
    save_plot("plots/pairs", "pca_pairs_matrix.png", p, width = 12, height = 12)
}
pc_pairs_plot(5)

# M) Variance explained bar plot with cumulative line overlay
variance_barplot <- function() {
    var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
    cum_var <- cumsum(var_exp)
    n_show <- min(20, length(var_exp))
    
    df <- data.frame(
        PC = factor(paste0("PC", 1:n_show), levels = paste0("PC", 1:n_show)),
        Variance = var_exp[1:n_show],
        Cumulative = cum_var[1:n_show]
    )
    
    p <- ggplot(df, aes(x = PC)) +
        geom_col(aes(y = Variance), fill = "#6B5B95", alpha = 0.7) +
        geom_line(aes(y = Cumulative, group = 1), color = "#66C1A4", size = 1) +
        geom_point(aes(y = Cumulative), color = "#66C1A4", size = 2) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                          sec.axis = sec_axis(~., name = "Cumulative Variance")) +
        theme_pca_min() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "PCA Variance Explained", x = "Principal Component", y = "Individual Variance")
    
    save_plot("plots/variance", "pca_variance_barplot_with_cumulative.svg", p, width = 10, height = 6)
}
variance_barplot()

# N) Sample-to-sample distance heatmap (Euclidean distance in PC space)
sample_distance_heatmap <- function(n_pcs = 10) {
    pc_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    dist_mat <- as.matrix(dist(pc_mat, method = "euclidean"))
    
    ann_cols <- intersect(c("celltype", "ExpGroup", "plate"), names(meta))
    if (length(ann_cols) > 0) {
        annotation_col <- meta[colnames(dist_mat), ann_cols, drop = FALSE]
        annotation_row <- annotation_col
    }
    
    # Use sampleNumber as row names if available
    if ("sampleNumber" %in% names(meta)) {
        rownames(dist_mat) <- meta[colnames(dist_mat), "sampleNumber"]
    }
    
    save_pheatmap(dist_mat,
                 file.path(ensure_dir(subdir("plots/heatmaps")), "sample_distance_euclidean_pc_space.png"),
                 width = 10, height = 10, scale = "none")
}
sample_distance_heatmap(10)

# O) Loading magnitude plot: show absolute loadings across PCs
loading_magnitude_plot <- function(n_pcs = 5, top_n = 20) {
    rot <- as.data.frame(pca$rotation)
    pc_cols <- paste0("PC", 1:min(n_pcs, ncol(rot)))
    
    # Calculate magnitude across selected PCs
    rot$magnitude <- sqrt(rowSums(rot[, pc_cols]^2))
    rot$protein <- rownames(rot)
    
    top_proteins <- head(rot[order(-rot$magnitude), ], top_n)
    
    # Reshape for plotting
    top_long <- reshape2::melt(top_proteins[, c("protein", pc_cols)], 
                              id.vars = "protein", 
                              variable.name = "PC", 
                              value.name = "loading")
    
    p <- ggplot(top_long, aes(x = reorder(protein, loading), y = loading, fill = PC)) +
        geom_col(position = "dodge") +
        coord_flip() +
        scale_fill_manual(values = make_modern_palette(n_pcs)) +
        theme_pca_min() +
        labs(title = paste("Top", top_n, "Proteins by Loading Magnitude"),
             x = "Protein", y = "Loading Value")
    
    save_plot("plots/loadings", "loading_magnitude_top_proteins.svg", p, width = 8, height = 10)
}
loading_magnitude_plot(5, 20)

# P) PC trajectory plot: connect samples by experimental order if available
pc_trajectory_plot <- function() {
    if (!"sampleNumber" %in% names(meta)) {
        message("sampleNumber not in metadata; skipping trajectory plot")
        return(NULL)
    }
    
    df <- data.frame(pca$x[, 1:2], 
                     sampleNumber = meta[rownames(pca$x), "sampleNumber"],
                     celltype = meta[rownames(pca$x), "celltype"])
    df <- df[order(df$sampleNumber), ]
    
    pal <- make_modern_palette(nlevels(droplevels(factor(df$celltype))))
    
    p <- ggplot(df, aes(PC1, PC2, color = celltype)) +
        geom_path(color = "gray70", alpha = 0.5, size = 0.5) +
        geom_point(size = 4, alpha = 0.8) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "PCA Trajectory by Sample Order", 
             subtitle = "Lines connect consecutive samples")
    
    save_plot("plots/trajectory", "pca_sample_trajectory.svg", p)
}
pc_trajectory_plot()

# Q) Comparison of different normalization methods
compare_normalizations <- function() {
    if (!exists("mat") || is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
        message("Matrix 'mat' is not defined or is empty; skipping normalization comparison.")
        return(NULL)
    }
    
    # Check if matrix has valid data
    if (!is.matrix(mat) || !is.numeric(mat)) {
        message("Matrix 'mat' is not a numeric matrix; skipping normalization comparison.")
        return(NULL)
    }
    
    # Remove rows/columns with all NAs or infinite values
    mat_clean <- mat[rowSums(is.finite(mat)) > 0, , drop = FALSE]
    mat_clean <- mat_clean[, colSums(is.finite(mat_clean)) > 0, drop = FALSE]
    
    if (nrow(mat_clean) < 2 || ncol(mat_clean) < 2) {
        message("Matrix too small after cleaning; skipping normalization comparison.")
        return(NULL)
    }
    
    # Current data (already processed) - use as-is
    pca_current <- prcomp(t(mat_clean), center = TRUE, scale. = TRUE)
    
    # Z-score normalization (per protein)
    mat_zscore <- t(scale(t(mat_clean), center = TRUE, scale = TRUE))
    pca_zscore <- prcomp(t(mat_zscore), center = TRUE, scale. = TRUE)
    
    # Median-centered normalization
    mat_median <- sweep(mat_clean, 1, apply(mat_clean, 1, median, na.rm = TRUE), "-")
    pca_median <- prcomp(t(mat_median), center = TRUE, scale. = TRUE)
    
    # Plot comparison
    plot_list <- list(
        list(pca = pca_current, title = "Current (PCA-adjusted)", file = "norm_current.svg"),
        list(pca = pca_zscore, title = "Z-score normalized", file = "norm_zscore.svg"),
        list(pca = pca_median, title = "Median-centered", file = "norm_median.svg")
    )
    
    for (item in plot_list) {
        # Get common samples between PCA and metadata
        common_samples <- intersect(rownames(item$pca$x), rownames(meta))
        if (length(common_samples) < 2) {
            message(paste("Skipping", item$title, "- no common samples with metadata"))
            next
        }
        
        df <- data.frame(
            item$pca$x[common_samples, 1:2, drop = FALSE],
            celltype = meta[common_samples, "celltype"]
        )
        
        pal <- make_modern_palette(nlevels(droplevels(factor(df$celltype))))
        
        p <- ggplot(df, aes(PC1, PC2, color = celltype)) +
            geom_point(size = 4, alpha = 0.8, stroke = 0) +
            scale_color_manual(values = pal) +
            theme_pca_min() +
            labs(title = paste("PCA:", item$title), color = "celltype")
        
        save_plot("plots/normalization", item$file, p)
    }
    
    message("Normalization comparison completed successfully!")
}
compare_normalizations()

# R) Metadata correlation heatmap with PCs
metadata_pc_correlation <- function(n_pcs = 10) {
    pc_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    
    # Select numeric metadata columns
    numeric_meta <- meta[, sapply(meta, is.numeric), drop = FALSE]
    
    if (ncol(numeric_meta) == 0) {
        message("No numeric metadata found for correlation")
        return(NULL)
    }
    
    # Align samples
    common_samples <- intersect(rownames(pc_mat), rownames(numeric_meta))
    pc_mat <- pc_mat[common_samples, ]
    numeric_meta <- numeric_meta[common_samples, ]
    
    # Calculate correlations
    cor_mat <- cor(pc_mat, numeric_meta, use = "pairwise.complete.obs")
    
    save_pheatmap(t(cor_mat),
                 file.path(ensure_dir(subdir("plots/heatmaps")), "metadata_pc_correlations.png"),
                 width = 8, height = 6, scale = "none")
    
    # Also save as table
    cor_df <- as.data.frame(as.table(cor_mat))
    colnames(cor_df) <- c("PC", "Metadata", "Correlation")
    save_table("tables/correlations", "metadata_pc_correlations.csv", cor_df)
}
metadata_pc_correlation(10)

# S) 3D PCA plot (PC1, PC2, PC3) if plotly available
pca_3d_plot <- function() {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        message("plotly not installed; skipping 3D PCA plot")
        return(NULL)
    }
    
    df <- data.frame(pca$x[, 1:3],
                     celltype = meta[rownames(pca$x), "celltype"],
                     sample = rownames(pca$x))
    
    pal <- make_modern_palette(nlevels(droplevels(factor(df$celltype))))
    
    p <- plotly::plot_ly(df, x = ~PC1, y = ~PC2, z = ~PC3, 
                        color = ~celltype, colors = pal,
                        text = ~sample, type = "scatter3d", mode = "markers",
                        marker = list(size = 5, opacity = 0.8)) %>%
        plotly::layout(title = "3D PCA Plot (PC1-PC3)",
                      scene = list(xaxis = list(title = "PC1"),
                                  yaxis = list(title = "PC2"),
                                  zaxis = list(title = "PC3")))
    
    htmlwidgets::saveWidget(p, 
                           file.path(ensure_dir(subdir("plots/3d")), "pca_3d_interactive.html"),
                           selfcontained = FALSE)
}
pca_3d_plot()

# T) Density plots of PC scores by group
pc_density_plots <- function(n_pcs = 3) {
    for (i in 1:min(n_pcs, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            score = pca$x[, i],
            celltype = meta[rownames(pca$x), "celltype"]
        )
        
        pal <- make_modern_palette(nlevels(droplevels(factor(df$celltype))))
        
        p <- ggplot(df, aes(x = score, fill = celltype)) +
            geom_density(alpha = 0.6) +
            scale_fill_manual(values = pal) +
            theme_pca_min() +
            labs(title = paste(pc_name, "Distribution by Celltype"),
                 x = paste(pc_name, "Score"), y = "Density")
        
        save_plot("plots/density", paste0(pc_name, "_density_by_celltype.svg"), p)
    }
}
pc_density_plots(3)

message("All supplementary plots completed!")

# ================== Additional Main & Supplementary Figure Plots ==================

# U) Boxplots of PC scores by experimental groups
pc_boxplots_by_group <- function(n_pcs = 5, group_key = "ExpGroup") {
    if (!group_key %in% names(meta)) {
        message(paste(group_key, "not in metadata; skipping PC boxplots"))
        return(NULL)
    }
    
    for (i in 1:min(n_pcs, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            score = pca$x[, i],
            group = factor(meta[rownames(pca$x), group_key])
        )
        df <- df[!is.na(df$group), ]
        
        pal <- make_modern_palette(nlevels(df$group))
        
        p <- ggplot(df, aes(x = group, y = score, fill = group)) +
            geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
            scale_fill_manual(values = pal) +
            theme_pca_min() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(title = paste(pc_name, "Scores by", group_key),
                 x = group_key, y = paste(pc_name, "Score"))
        
        save_plot("plots/boxplots", paste0(pc_name, "_boxplot_by_", group_key, ".svg"), p)
    }
}
pc_boxplots_by_group(5, "ExpGroup")
pc_boxplots_by_group(5, "celltype")

# V) Violin plots with overlaid points and statistical tests
pc_violin_plots <- function(n_pcs = 5, group_key = "celltype") {
    if (!group_key %in% names(meta)) return(NULL)
    
    all_normality <- list()
    all_main_tests <- list()
    all_posthocs <- list()
    
    for (i in 1:min(n_pcs, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            score = pca$x[, i],
            group = factor(meta[rownames(pca$x), group_key])
        )
        df <- df[!is.na(df$group), ]
        
        cat("\n=== Processing", pc_name, "for", group_key, "===\n")
        
        pal <- make_modern_palette(nlevels(df$group))
        
        # Create base plot
        p <- ggplot(df, aes(x = group, y = score, fill = group, color = group)) +
            geom_hline(yintercept = 0, linetype = "solid", color = "lightgrey", size = 1) + 
            geom_violin(alpha = 0.4, trim = FALSE, show.legend = FALSE, linetype = "blank") +
            geom_jitter(width = 0.15, alpha = 0.8, size = 4) +
            scale_fill_manual(values = pal) +
            scale_color_manual(values = pal) +
            theme_pca_min() + 
            theme(panel.grid = element_blank(),
                  axis.text.x = element_text(hjust = 0.5),
                  legend.position = "none") +
            labs(title = paste(pc_name, "Distribution by", group_key),
                 x = group_key, y = paste(pc_name, "Score"))
            
        
        if (length(unique(df$group)) > 1) {
            p <- p + stat_summary(fun = median, geom = "crossbar", width = 0.5, 
                                 color = "black", size = 0.4, fatten = 1)
        }
        
        # Statistical testing - normality
        normality_by_group <- by(df$score, df$group, function(x) {
            if (length(x) < 3) return(data.frame(W = NA, p = NA))
            test <- shapiro.test(x)
            data.frame(W = test$statistic, p = test$p.value)
        })
        
        normality_df <- do.call(rbind, lapply(names(normality_by_group), function(g) {
            data.frame(
                PC = pc_name,
                Group = g,
                Test = "Shapiro-Wilk",
                Statistic = normality_by_group[[g]]$W,
                P_Value = normality_by_group[[g]]$p,
                stringsAsFactors = FALSE
            )
        }))
        
        all_normality[[pc_name]] <- normality_df
        
        all_normal <- all(normality_df$P_Value > 0.05, na.rm = TRUE)
        n_groups <- nlevels(df$group)
        
        # Main test
        main_test_df <- NULL
        posthoc_df <- NULL
        
        if (n_groups >= 2) {
            if (all_normal) {
                anova_result <- aov(score ~ group, data = df)
                anova_summary <- summary(anova_result)
                
                main_test_df <- data.frame(
                    PC = pc_name,
                    Grouping = group_key,
                    Test = "ANOVA",
                    N_Groups = n_groups,
                    All_Normal = all_normal,
                    Statistic = anova_summary[[1]][["F value"]][1],
                    P_Value = anova_summary[[1]][["Pr(>F)"]][1],
                    stringsAsFactors = FALSE
                )
                
                if (main_test_df$P_Value < 0.05 && n_groups > 2) {
                    posthoc <- TukeyHSD(anova_result)
                    posthoc_df <- as.data.frame(posthoc$group)
                    posthoc_df$Comparison <- rownames(posthoc_df)
                    posthoc_df$PC <- pc_name
                    posthoc_df$Grouping <- group_key
                    posthoc_df$Test <- "Tukey HSD"
                    posthoc_df <- posthoc_df[, c("PC", "Grouping", "Test", "Comparison", "diff", "lwr", "upr", "p adj")]
                    colnames(posthoc_df) <- c("PC", "Grouping", "Test", "Comparison", "Diff", "Lower_CI", "Upper_CI", "P_Value")
                }
                
            } else {
                kw_result <- kruskal.test(score ~ group, data = df);
                
                main_test_df <- data.frame(
                    PC = pc_name,
                    Grouping = group_key,
                    Test = "Kruskal-Wallis",
                    N_Groups = n_groups,
                    All_Normal = all_normal,
                    Statistic = kw_result$statistic,
                    P_Value = kw_result$p.value,
                    stringsAsFactors = FALSE
                );
                
                if (main_test_df$P_Value < 0.05 && n_groups > 2) {
                    pw_result <- pairwise.wilcox.test(df$score, df$group, 
                                                      p.adjust.method = "bonferroni");
                    
                    pw_matrix <- pw_result$p.value;
                    comparisons <- which(!is.na(pw_matrix), arr.ind = TRUE);
                    
                    if (nrow(comparisons) > 0) {
                        posthoc_df <- data.frame(
                            PC = pc_name,
                            Grouping = group_key,
                            Test = "Pairwise Wilcoxon",
                            Comparison = paste(rownames(pw_matrix)[comparisons[,1]], 
                                             colnames(pw_matrix)[comparisons[,2]], 
                                             sep = " vs "),
                            Diff = NA_real_,
                            Lower_CI = NA_real_,
                            Upper_CI = NA_real_,
                            P_Value = pw_matrix[comparisons],
                            stringsAsFactors = FALSE
                        );
                    }
                }
            }
        }
        
        # Calculate y-axis limits BEFORE adding significance annotations
        y_range <- range(df$score, na.rm = TRUE)
        y_span <- diff(y_range)
        y_min <- y_range[1] - y_span * 0.15;  # Bottom padding
        y_start <- y_range[2] + y_span * 0.15;  # Starting position for brackets
        y_step <- y_span * 0.25;  # Step size between brackets
        
        # Count number of significant comparisons to calculate required space
        n_sig_comparisons <- 0;
        if (!is.null(posthoc_df) && nrow(posthoc_df) > 0) {
            n_sig_comparisons <- sum(posthoc_df$P_Value < 0.05);
        }
        
        # Add significance annotations
        max_y <- y_start;
        if (n_sig_comparisons > 0) {
            posthoc_df$group1 <- sapply(strsplit(posthoc_df$Comparison, " vs "), `[`, 1);
            posthoc_df$group2 <- sapply(strsplit(posthoc_df$Comparison, " vs "), `[`, 2);
            
            sig_comparisons <- posthoc_df[posthoc_df$P_Value < 0.05, ];
            
            for (idx in 1:nrow(sig_comparisons)) {
                y_pos <- y_start + (idx - 1) * y_step;
                max_y <- max(max_y, y_pos);
                
                p_val <- sig_comparisons$P_Value[idx];
                sig_symbol <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else "ns";
                
                x1 <- which(levels(df$group) == sig_comparisons$group1[idx]);
                x2 <- which(levels(df$group) == sig_comparisons$group2[idx]);
                
                p <- p + 
                    annotate("segment", x = x1, xend = x2, y = y_pos, yend = y_pos, color = "black", size = 0.5) +
                    annotate("segment", x = x1, xend = x1, y = y_pos, yend = y_pos - y_span * 0.02, color = "black", size = 0.5) +
                    annotate("segment", x = x2, xend = x2, y = y_pos, yend = y_pos - y_span * 0.02, color = "black", size = 0.5) +
                    annotate("text", x = (x1 + x2) / 2, y = y_pos + y_span * 0.04, label = sig_symbol, size = 5, fontface = "bold");
            }
        }
        
        # Set y-axis limits with generous padding
        # Add extra space above the highest bracket plus text
        y_max <- max_y + y_span * 0.25;  # Generous top padding
        p <- p + coord_cartesian(ylim = c(y_min, y_max), clip = "off");
        
        # Save plot
        save_plot("plots/violin", paste0(pc_name, "_violin_by_", group_key, ".svg"), p, width = 5, height = 5);
        
        if (!is.null(main_test_df)) all_main_tests[[pc_name]] <- main_test_df;
        if (!is.null(posthoc_df)) all_posthocs[[pc_name]] <- posthoc_df;
    }
    
    # Save summary tables
    if (length(all_main_tests) > 0) {
        combined_main_tests <- do.call(rbind, all_main_tests);
        combined_main_tests$Significant <- ifelse(combined_main_tests$P_Value < 0.05, "Yes", "No");
        save_table("tables/statistics", paste0("summary_main_tests_", group_key, ".csv"), combined_main_tests);
    }
    
    if (length(all_posthocs) > 0) {
        combined_posthocs <- do.call(rbind, all_posthocs);
        combined_posthocs$Significant <- ifelse(combined_posthocs$P_Value < 0.05, "Yes", "No");
        save_table("tables/statistics", paste0("summary_posthoc_tests_", group_key, ".csv"), combined_posthocs);
    }
    
    if (length(all_normality) > 0) {
        combined_normality <- do.call(rbind, all_normality);
        save_table("tables/statistics", paste0("summary_normality_tests_", group_key, ".csv"), combined_normality);
    }
}
pc_violin_plots(5, "celltype")
pc_violin_plots(5, "ExpGroup")

# W) Loading contribution per PC (top positive and negative)
loading_contribution_bars <- function(n_pcs = 5, top_n = 10) {
    for (i in 1:min(n_pcs, ncol(pca$rotation))) {
        pc_name <- paste0("PC", i)
        loadings <- pca$rotation[, i]
        
        # Get top positive and negative
        top_pos <- head(sort(loadings, decreasing = TRUE), top_n)
        top_neg <- head(sort(loadings, decreasing = FALSE), top_n)
        
        df <- data.frame(
            protein = c(names(top_pos), names(top_neg)),
            loading = c(top_pos, top_neg),
            direction = c(rep("Positive", top_n), rep("Negative", top_n))
        )
        df$protein <- factor(df$protein, levels = df$protein[order(df$loading)])
        
        p <- ggplot(df, aes(x = protein, y = loading, fill = direction)) +
            geom_col() +
            coord_flip() +
            scale_fill_manual(values = c("Positive" = "#FFB5B5", "Negative" = "#B5D4FF")) +
            theme_pca_min() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(title = paste("Top Loadings for", pc_name),
                 x = "Protein", y = "Loading Value", fill = "Direction")
        
        save_plot("plots/loadings", paste0(pc_name, "_top_loadings_bar.svg"), p, width = 6, height = 4)
    }
}
loading_contribution_bars(5, 5)

# X) Correlation circle plot (biplot style, loadings only)
correlation_circle_plot <- function(pc_x = 1, pc_y = 2, top_n = 30) {
    rot <- as.data.frame(pca$rotation)
    
    # Select top contributors to both PCs
    mag <- sqrt(rot[, pc_x]^2 + rot[, pc_y]^2)
    top_idx <- order(mag, decreasing = TRUE)[1:min(top_n, length(mag))]
    
    df <- data.frame(
        protein = rownames(rot)[top_idx],
        PC1 = rot[top_idx, pc_x],
        PC2 = rot[top_idx, pc_y]
    )
    
    # Create circle
    theta <- seq(0, 2*pi, length.out = 100)
    circle <- data.frame(x = cos(theta), y = sin(theta))
    
    p <- ggplot() +
        geom_path(data = circle, aes(x, y), color = "gray70", linetype = "dashed") +
        geom_segment(data = df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                    arrow = arrow(length = unit(0.2, "cm")), alpha = 0.6, color = "#6B5B95") +
        geom_text_repel(data = df, aes(x = PC1, y = PC2, label = protein),
                       size = 3, max.overlaps = 30) +
        coord_fixed() +
        theme_pca_min() +
        labs(title = paste("Correlation Circle (PC", pc_x, "vs PC", pc_y, ")"),
             x = paste0("PC", pc_x, " Loadings"),
             y = paste0("PC", pc_y, " Loadings"))
    
    save_plot("plots/biplot", paste0("correlation_circle_PC", pc_x, "_PC", pc_y, ".svg"), p)
}
correlation_circle_plot(1, 2, 30)
correlation_circle_plot(2, 3, 30)

# Y) Sample QC metrics vs PC scores
sample_qc_vs_pcs <- function(n_pcs = 2) {
    # Calculate QC metrics
    qc_df <- data.frame(
        sample = colnames(mat),
        n_detected = colSums(!is.na(mat) & mat > 0),
        mean_intensity = colMeans(mat, na.rm = TRUE),
        median_intensity = apply(mat, 2, median, na.rm = TRUE),
        missing_pct = colMeans(is.na(mat)) * 100
    )
    
    # Merge with PC scores
    qc_df <- cbind(qc_df, pca$x[qc_df$sample, 1:min(n_pcs, ncol(pca$x)), drop = FALSE])
    
    # Plot each QC metric vs PCs
    for (qc_col in c("n_detected", "mean_intensity", "missing_pct")) {
        for (i in 1:min(n_pcs, ncol(pca$x))) {
            pc_name <- paste0("PC", i)
            
            p <- ggplot(qc_df, aes_string(x = pc_name, y = qc_col)) +
                geom_point(alpha = 0.7, size = 3, color = "#6B5B95") +
                geom_smooth(method = "lm", se = TRUE, color = "#66C1A4") +
                theme_pca_min() +
                labs(title = paste(qc_col, "vs", pc_name),
                     x = paste(pc_name, "Score"),
                     y = qc_col)
            
            save_plot("plots/qc", paste0(qc_col, "_vs_", pc_name, ".svg"), p)
        }
    }
    
    # Save QC table
    save_table("tables/qc", "sample_qc_metrics.csv", qc_df)
}
sample_qc_vs_pcs(3)

# Z) Hierarchical clustering dendrogram colored by metadata
sample_dendrogram <- function(n_pcs = 10, group_key = "celltype") {
    if (!group_key %in% names(meta)) return(NULL)
    
    pc_mat <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    hc <- hclust(dist(pc_mat), method = "complete")
    
    # Color by group
    groups <- factor(meta[rownames(pc_mat), group_key])
    pal <- make_modern_palette(nlevels(droplevels(groups)))
    colors <- pal[as.numeric(droplevels(groups))]
    
    png(file.path(ensure_dir(subdir("plots/dendrogram")), 
                  paste0("dendrogram_by_", group_key, ".png")),
        width = 12, height = 8, units = "in", res = 150)
    plot(hc, hang = -1, cex = 0.6, main = paste("Sample Dendrogram colored by", group_key))
    rect.hclust(hc, k = nlevels(droplevels(groups)), border = pal)
    dev.off()
}
sample_dendrogram(10, "celltype")

# AA) Scree plot with Kaiser criterion line (eigenvalue = 1)
scree_with_kaiser <- function() {
    eigenvalues <- pca$sdev^2
    df <- data.frame(
        PC = factor(paste0("PC", 1:length(eigenvalues)), 
                   levels = paste0("PC", 1:length(eigenvalues))),
        Eigenvalue = eigenvalues
    )
    df <- head(df, 20)  # Show first 20
    
    p <- ggplot(df, aes(x = PC, y = Eigenvalue)) +
        geom_col(fill = "#6B5B95", alpha = 0.7) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
        annotate("text", x = 2, y = 1.2, label = "Kaiser criterion (λ=1)", 
                color = "red", hjust = 0) +
        theme_pca_min() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Scree Plot with Kaiser Criterion",
             x = "Principal Component", y = "Eigenvalue")
    
    save_plot("plots/variance", "scree_with_kaiser_criterion.svg", p)
}
scree_with_kaiser()

# AB) PC time series (if temporal/sequential metadata available)
pc_timeseries <- function(time_key = "sampleNumber", group_key = "celltype") {
    if (!time_key %in% names(meta) || !group_key %in% names(meta)) {
        message("Required metadata not available for time series")
        return(NULL)
    }
    
    df <- data.frame(
        pca$x[, 1:2],
        time = meta[rownames(pca$x), time_key],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$time) & !is.na(df$group), ]
    df <- df[order(df$time), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(x = time, y = PC1, color = group)) +
        geom_line(alpha = 0.6) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "PC1 Over Sample Order",
             x = "Sample Number", y = "PC1 Score", color = group_key)
    
    save_plot("plots/timeseries", "pc1_timeseries.svg", p)
}
pc_timeseries("sampleNumber", "celltype")

# AC) Batch effect visualization before/after correction
batch_effect_comparison <- function() {
    if (!"plate" %in% names(meta)) return(NULL)
    
    # Before correction
    df_pre <- data.frame(pca$x[, 1:2], plate = meta[rownames(pca$x), "plate"])
    pal <- make_modern_palette(nlevels(droplevels(factor(df_pre$plate))))
    
    p1 <- ggplot(df_pre, aes(PC1, PC2, color = factor(plate))) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "Before Batch Correction", color = "Plate")
    
    save_plot("plots/batch", "batch_effect_before.svg", p1)
}
batch_effect_comparison()

message("All additional supplementary and main figure plots completed!")

# ================== Additional Main Figure Plots ==================

# AD) Side-by-side PCA plots for different groupings (publication quality)
combined_pca_panels <- function() {
    # Create 2x2 panel of key groupings
    groupings <- c("celltype", "ExpGroup", "ReplicateGroup", "plate")
    groupings <- intersect(groupings, names(meta))
    
    if (length(groupings) < 2) {
        message("Not enough grouping variables for panel plot")
        return(NULL)
    }
    
    plot_list <- list()
    for (key in groupings) {
        grp <- build_group(key)
        keep <- !is.na(grp)
        if (sum(keep) < 2) next
        
        grp2 <- droplevels(grp[keep])
        df <- data.frame(pca$x[keep, 1:2, drop=FALSE], group = grp2)
        pal <- make_modern_palette(nlevels(grp2))
        
        varp <- (pca$sdev^2)/sum(pca$sdev^2)
        
        p <- ggplot(df, aes(PC1, PC2, color = group)) +
            geom_point(size = 4, alpha = 0.8, stroke = 0) +
            scale_color_manual(values = pal) +
            theme_pca_min() +
            labs(x = sprintf("PC1 (%.1f%%)", 100*varp[1]),
                 y = sprintf("PC2 (%.1f%%)", 100*varp[2]),
                 color = key,
                 title = key)
        
        plot_list[[key]] <- p
    }
    
    if (length(plot_list) > 0 && requireNamespace("gridExtra", quietly = TRUE)) {
        combined <- gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
        ggsave(file.path(ensure_dir(subdir("plots/main_figure")), "pca_combined_panels.svg"),
               combined, width = 14, height = 12, dpi = 300)
    }
}
combined_pca_panels()

# AE) Publication-ready variance explained with cumulative overlay
variance_main_figure <- function() {
    var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
    cum_var <- cumsum(var_exp)
    n_show <- min(15, length(var_exp))
    
    df <- data.frame(
        PC = 1:n_show,
        PC_label = factor(paste0("PC", 1:n_show), levels = paste0("PC", 1:n_show)),
        Individual = var_exp[1:n_show] * 100,
        Cumulative = cum_var[1:n_show] * 100
    )
    
    # Dual-axis plot
    p <- ggplot(df, aes(x = PC)) +
        geom_col(aes(y = Individual), fill = "#6B5B95", alpha = 0.8, width = 0.7) +
        geom_line(aes(y = Cumulative), color = "#66C1A4", size = 1.5, group = 1) +
        geom_point(aes(y = Cumulative), color = "#66C1A4", size = 3) +
        geom_hline(yintercept = 80, linetype = "dashed", color = "#E89369", size = 0.8) +
        annotate("text", x = n_show * 0.8, y = 82, label = "80% threshold", 
                color = "#E89369", size = 4) +
        scale_x_continuous(breaks = 1:n_show, labels = paste0("PC", 1:n_show)) +
        scale_y_continuous(limits = c(0, 100)) +
        theme_pca_min() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            panel.grid.major.x = element_blank()
        ) +
        labs(title = "Variance Explained by Principal Components",
             x = "Principal Component",
             y = "Variance Explained (%)")
    
    save_plot("plots/main_figure", "variance_explained_main.svg", p, width = 10, height = 6)
}
variance_main_figure()

# AF) Top contributing proteins heatmap (main figure version)
top_proteins_heatmap_main <- function(n_pcs = 3, top_per_pc = 15) {
    if (is.null(pca) || is.null(pca$rotation) || ncol(pca$rotation) == 0 || 
        is.null(mat) || ncol(mat) == 0 || nrow(mat) == 0) {
        message("PCA object or matrix is not properly defined")
        return(NULL)
    }
    
    rot <- as.data.frame(pca$rotation)
    
    # Get top contributors for each PC
    top_proteins <- c()
    for (i in 1:min(n_pcs, ncol(rot))) {
        pc_col <- paste0("PC", i)
        if (!pc_col %in% colnames(rot)) next
        
        # Top positive
        top_pos <- head(names(sort(rot[[pc_col]], decreasing = TRUE)), top_per_pc)
        # Top negative
        top_neg <- head(names(sort(rot[[pc_col]], decreasing = FALSE)), top_per_pc)
        top_proteins <- c(top_proteins, top_pos, top_neg)
    }
    top_proteins <- unique(top_proteins)
    
    # Filter to proteins that exist in mat
    top_proteins <- intersect(top_proteins, rownames(mat))
    
    if (length(top_proteins) < 4) {
        message("Not enough valid proteins found for heatmap")
        return(NULL)
    }
    
    # Create expression matrix for these proteins
    expr_subset <- mat[top_proteins, , drop = FALSE]
    
    # Remove proteins with all NA or no variation
    valid_rows <- rowSums(!is.na(expr_subset)) > 0 & 
                  apply(expr_subset, 1, function(x) var(x, na.rm = TRUE) > 0)
    expr_subset <- expr_subset[valid_rows, , drop = FALSE]
    
    if (nrow(expr_subset) < 4) {
        message("Not enough valid proteins after filtering")
        return(NULL)
    }
    
    # Order samples by PC1
    sample_order <- order(pca$x[, 1])
    expr_subset <- expr_subset[, sample_order, drop = FALSE]
    
    # Annotation
    ann_col <- data.frame(
        celltype = meta[colnames(expr_subset), "celltype"],
        row.names = colnames(expr_subset)
    )
    
    ann_colors <- list(
        celltype = setNames(make_modern_palette(nlevels(factor(ann_col$celltype))),
                           levels(factor(ann_col$celltype)))
    )
    
    pheatmap::pheatmap(
        expr_subset,
        annotation_col = ann_col,
        annotation_colors = ann_colors,
        show_colnames = FALSE,
        scale = "row",
        clustering_method = "complete",
        color = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100),
        filename = file.path(ensure_dir(subdir("plots/main_figure")), 
                           "top_proteins_heatmap_main.png"),
        width = 10,
        height = 12
    )
}
top_proteins_heatmap_main(3, 15)

# AG) PCA with confidence ellipses (main figure version)
pca_with_ellipses_main <- function(group_key = "celltype", conf_level = 0.95) {
    if (!group_key %in% names(meta)) return(NULL)
    
    grp <- build_group(group_key)
    keep <- !is.na(grp)
    grp2 <- droplevels(grp[keep])
    
    if (nlevels(grp2) < 2) {
        message(paste("Not enough groups for", group_key))
        return(NULL)
    }
    
    df <- data.frame(pca$x[keep, 1:2, drop=FALSE], group = grp2)
    pal <- make_modern_palette(nlevels(grp2))
    
    varp <- (pca$sdev^2)/sum(pca$sdev^2)
    
    p <- ggplot(df, aes(PC1, PC2, color = group, fill = group)) +
        stat_ellipse(type = "t", level = conf_level, geom = "polygon", 
                    alpha = 0.15, show.legend = FALSE) +
        stat_ellipse(type = "t", level = conf_level, size = 1, show.legend = FALSE) +
        geom_point(size = 5, alpha = 0.8, stroke = 0, shape = 16) +
        scale_color_manual(values = pal) +
        scale_fill_manual(values = pal) +
        theme_pca_min() +
        theme(
            legend.position = "right",
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 16, face = "bold")
        ) +
        labs(x = sprintf("PC1 (%.1f%%)", 100*varp[1]),
             y = sprintf("PC2 (%.1f%%)", 100*varp[2]),
             color = group_key,
             fill = group_key,
             title = paste("PCA with", conf_level*100, "% Confidence Ellipses"))
    
    save_plot("plots/main_figure", paste0("pca_ellipses_", group_key, "_main.svg"), p, width = 10, height = 8)
}
pca_with_ellipses_main("celltype", 0.95)
pca_with_ellipses_main("ExpGroup", 0.95)

# AH) Loading plot for top contributing proteins (main figure)
loadings_main_figure <- function(pc_x = 1, pc_y = 2, top_n = 25) {
    rot <- as.data.frame(pca$rotation)
    
    # Calculate contribution magnitude
    mag <- sqrt(rot[, pc_x]^2 + rot[, pc_y]^2)
    top_idx <- order(mag, decreasing = TRUE)[1:min(top_n, length(mag))]
    
    df <- data.frame(
        protein = rownames(rot)[top_idx],
        PC1_load = rot[top_idx, pc_x],
        PC2_load = rot[top_idx, pc_y],
        magnitude = mag[top_idx]
    )
    
    # Create arrow plot
    p <- ggplot(df, aes(x = PC1_load, y = PC2_load)) +
        geom_segment(aes(x = 0, y = 0, xend = PC1_load, yend = PC2_load),
                    arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                    color = "#6B5B95", alpha = 0.7, size = 0.8) +
        geom_text_repel(aes(label = protein), size = 3.5, 
                       max.overlaps = 30, box.padding = 0.5) +
        geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
        coord_fixed() +
        theme_pca_min() +
        labs(title = paste("Top", top_n, "Contributing Proteins"),
             x = paste0("PC", pc_x, " Loading"),
             y = paste0("PC", pc_y, " Loading"))
    
    save_plot("plots/main_figure", paste0("loadings_arrows_PC", pc_x, "_PC", pc_y, "_main.svg"), 
             p, width = 10, height = 10)
}
loadings_main_figure(1, 2, 25)

# AI) Statistical summary barplot (ANOVA/Kruskal results)
statistical_summary_barplot <- function(group_key = "celltype") {
    if (!group_key %in% names(meta)) return(NULL)
    
    # Run tests for multiple PCs
    results <- list()
    for (i in 1:min(10, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            score = pca$x[, i],
            group = factor(meta[rownames(pca$x), group_key])
        )
        df <- df[!is.na(df$group), ]
        
        # Check normality
        normality <- by(df$score, df$group, function(x) {
            if (length(x) < 3) return(list(p = NA))
            shapiro.test(x)
        })
        all_normal <- all(sapply(normality, function(x) x$p.value > 0.05), na.rm = TRUE)
        
        # Run appropriate test
        if (all_normal) {
            test_result <- summary(aov(score ~ group, data = df))
            p_val <- test_result[[1]][["Pr(>F)"]][1]
            test_name <- "ANOVA"
        } else {
            test_result <- kruskal.test(score ~ group, data = df)
            p_val <- test_result$p.value
            test_name <- "Kruskal-Wallis"
        }
        
        results[[i]] <- data.frame(
            PC = pc_name,
            P_Value = p_val,
            NegLog10P = -log10(p_val),
            Test = test_name,
            Significant = p_val < 0.05
        )
    }
    
    results_df <- do.call(rbind, results)
    
    p <- ggplot(results_df, aes(x = PC, y = NegLog10P, fill = Significant)) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
                  color = "red", size = 1) +
        annotate("text", x = 2, y = -log10(0.05) + 0.3, 
                label = "p = 0.05", color = "red", size = 4) +
        scale_fill_manual(values = c("FALSE" = "#B2B2B2", "TRUE" = "#6B5B95"),
                         labels = c("ns", "p < 0.05")) +
        theme_pca_min() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("PC Association with", group_key),
             x = "Principal Component",
             y = "-log10(p-value)",
             fill = "Significance")
    
    save_plot("plots/main_figure", paste0("pc_stats_summary_", group_key, ".svg"), 
             p, width = 5, height = 4)
    
    # Save table
    save_table("tables/statistics", 
              paste0("pc_group_associations_", group_key, ".csv"), 
              results_df)
}
statistical_summary_barplot("celltype")
statistical_summary_barplot("ExpGroup")

# AI-UMAP) Statistical summary barplot for UMAP dimensions
statistical_summary_barplot_umap <- function(group_key = "celltype") {
    if (!group_key %in% names(meta)) return(NULL)
    if (is.null(um)) {
        message("UMAP not available; skipping UMAP statistical summary")
        return(NULL)
    }
    
    # Run tests for all UMAP dimensions
    results <- list()
    n_dims <- ncol(um)
    
    for (i in 1:n_dims) {
        dim_name <- paste0("UMAP", i)
        df <- data.frame(
            score = um[, i],
            group = factor(meta[rownames(um), group_key])
        )
        df <- df[!is.na(df$group), ]
        
        # Check normality
        normality <- by(df$score, df$group, function(x) {
            if (length(x) < 3) return(list(p = NA))
            shapiro.test(x)
        })
        all_normal <- all(sapply(normality, function(x) x$p.value > 0.05), na.rm = TRUE)
        
        # Run appropriate test
        if (all_normal) {
            test_result <- summary(aov(score ~ group, data = df))
            p_val <- test_result[[1]][["Pr(>F)"]][1]
            test_name <- "ANOVA"
        } else {
            test_result <- kruskal.test(score ~ group, data = df)
            p_val <- test_result$p.value
            test_name <- "Kruskal-Wallis"
        }
        
        results[[i]] <- data.frame(
            Dimension = dim_name,
            P_Value = p_val,
            NegLog10P = -log10(p_val),
            Test = test_name,
            Significant = p_val < 0.05,
            stringsAsFactors = FALSE
        )
    }
    
    results_df <- do.call(rbind, results)
    
    # Create factor with proper ordering
    results_df$Dimension <- factor(results_df$Dimension, 
                                   levels = paste0("UMAP", 1:n_dims))
    
    p <- ggplot(results_df, aes(x = Dimension, y = NegLog10P, fill = Significant)) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
                  color = "red", size = 1) +
        annotate("text", x = min(2, n_dims/2), y = -log10(0.05) + 0.3, 
                label = "p = 0.05", color = "red", size = 4) +
        scale_fill_manual(values = c("FALSE" = "#B2B2B2", "TRUE" = "#6B5B95"),
                         labels = c("ns", "p < 0.05")) +
        theme_pca_min() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = paste("UMAP Dimension Association with", group_key),
             x = "UMAP Dimension",
             y = "-log10(p-value)",
             fill = "Significance")
    
    save_plot("plots/main_figure", paste0("umap_stats_summary_", group_key, ".svg"), 
             p, width = max(5, n_dims * 0.5), height = 4)
    
    # Save table
    save_table("tables/statistics", 
              paste0("umap_group_associations_", group_key, ".csv"), 
              results_df)
    
    message(paste("Completed UMAP statistical summary for", n_dims, "dimensions by", group_key))
}
statistical_summary_barplot_umap("celltype")
statistical_summary_barplot_umap("ExpGroup")

# AJ) Sample distribution along PC1 (density + rug plot)
pc1_distribution_main <- function(group_key = "celltype") {
    if (!group_key %in% names(meta)) return(NULL)
    
    df <- data.frame(
        PC1 = pca$x[, 1],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(x = PC1, fill = group, color = group)) +
        geom_density(alpha = 0.5, size = 1.2) +
        geom_rug(aes(color = group), alpha = 0.8, size = 1) +
        scale_fill_manual(values = pal) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        theme(legend.position = "top") +
        labs(title = "PC1 Distribution by Group",
             x = "PC1 Score",
             y = "Density",
             fill = group_key,
             color = group_key)
    
    save_plot("plots/main_figure", paste0("pc1_distribution_", group_key, ".svg"), 
             p, width = 10, height = 6)
}
pc1_distribution_main("celltype")

# AK) Multi-panel: PC1 vs PC2, PC2 vs PC3, PC1 vs PC3
multipanel_pca_combinations <- function(group_key = "celltype") {
    if (!group_key %in% names(meta)) return(NULL)
    
    grp <- build_group(group_key)
    keep <- !is.na(grp)
    grp2 <- droplevels(grp[keep])
    pal <- make_modern_palette(nlevels(grp2))
    
    varp <- (pca$sdev^2)/sum(pca$sdev^2) * 100
    
    # PC1 vs PC2
    df12 <- data.frame(pca$x[keep, c(1,2), drop=FALSE], group = grp2)
    p12 <- ggplot(df12, aes(PC1, PC2, color = group)) +
        geom_point(size = 4, alpha = 0.8) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(x = sprintf("PC1 (%.1f%%)", varp[1]),
             y = sprintf("PC2 (%.1f%%)", varp[2]),
             title = "PC1 vs PC2")
    
    # PC2 vs PC3
    if (ncol(pca$x) >= 3) {
        df23 <- data.frame(pca$x[keep, c(2,3), drop=FALSE], group = grp2)
        p23 <- ggplot(df23, aes(PC2, PC3, color = group)) +
            geom_point(size = 4, alpha = 0.8) +
            scale_color_manual(values = pal) +
            theme_pca_min() +
            labs(x = sprintf("PC2 (%.1f%%)", varp[2]),
                 y = sprintf("PC3 (%.1f%%)", varp[3]),
                 title = "PC2 vs PC3")
        
        # PC1 vs PC3
        df13 <- data.frame(pca$x[keep, c(1,3), drop=FALSE], group = grp2)
        p13 <- ggplot(df13, aes(PC1, PC3, color = group)) +
            geom_point(size = 4, alpha = 0.8) +
            scale_color_manual(values = pal) +
            theme_pca_min() +
            labs(x = sprintf("PC1 (%.1f%%)", varp[1]),
                 y = sprintf("PC3 (%.1f%%)", varp[3]),
                 title = "PC1 vs PC3")
        
        if (requireNamespace("gridExtra", quietly = TRUE)) {
            combined <- gridExtra::grid.arrange(p12, p23, p13, ncol = 3)
            ggsave(file.path(ensure_dir(subdir("plots/main_figure")), 
                           paste0("pca_combinations_", group_key, ".svg")),
                   combined, width = 18, height = 6, dpi = 300)
        }
    }
}
multipanel_pca_combinations("celltype")

message("All main figure plots completed!")


# ================== Innovative Main Figure Visualizations ==================

# AL) Arc diagram showing sample relationships based on PC distance
pc_arc_diagram <- function(group_key = "celltype", n_connections = 50) {
    if (!requireNamespace("ggraph", quietly = TRUE) || !requireNamespace("igraph", quietly = TRUE)) {
        message("ggraph/igraph not available; skipping arc diagram")
        return(NULL)
    }
    
    # Calculate distances in PC space (first 5 PCs)
    pc_mat <- pca$x[, 1:min(5, ncol(pca$x)), drop = FALSE]
    dist_mat <- as.matrix(dist(pc_mat))
    
    # Get top N closest pairs
    dist_df <- reshape2::melt(dist_mat)
    dist_df <- dist_df[dist_df$Var1 != dist_df$Var2, ]
    dist_df <- dist_df[order(dist_df$value), ]
    dist_df <- head(dist_df, n_connections)
    
    # Create graph
    g <- igraph::graph_from_data_frame(dist_df[, 1:2], directed = FALSE)
    igraph::V(g)$group <- meta[igraph::V(g)$name, group_key]
    
    pal <- make_modern_palette(length(unique(igraph::V(g)$group)))
    
    p <- ggraph::ggraph(g, layout = "linear") +
        ggraph::geom_edge_arc(aes(alpha = ..index..), strength = 0.3, color = "gray70") +
        ggraph::geom_node_point(aes(color = group), size = 5) +
        scale_color_manual(values = pal) +
        theme_void() +
        theme(legend.position = "bottom") +
        labs(title = "Sample Similarity Network (Arc Diagram)", color = group_key)
    
    save_plot("plots/innovative", "pc_arc_diagram.svg", p, width = 14, height = 6)
}
pc_arc_diagram("celltype", 50)

# AM) Radial/polar PCA plot
radial_pca_plot <- function(group_key = "celltype") {
    df <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    # Convert to polar coordinates
    df$angle <- atan2(df$PC2, df$PC1)
    df$radius <- sqrt(df$PC1^2 + df$PC2^2)
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(x = angle, y = radius, color = group)) +
        geom_point(size = 5, alpha = 0.8) +
        coord_polar(theta = "x") +
        scale_color_manual(values = pal) +
        theme_minimal() +
        theme(
            panel.grid.major = element_line(color = "#ECECEC"),
            axis.text.x = element_blank(),
            axis.title = element_blank()
        ) +
        labs(title = "Radial PCA Projection", color = group_key)
    
    save_plot("plots/innovative", "radial_pca.svg", p, width = 8, height = 8)
}
radial_pca_plot("celltype")

# AN) Alluvial/Sankey diagram showing sample flow across PCA space
alluvial_pca_flow <- function(group_key = "celltype") {
    if (!requireNamespace("ggalluvial", quietly = TRUE)) {
        message("ggalluvial not available; skipping alluvial plot")
        return(NULL)
    }
    
    # Bin PC1 and PC2 into quartiles
    df <- data.frame(
        sample = rownames(pca$x),
        PC1_bin = cut(pca$x[, 1], breaks = 4, labels = c("Low", "Med-Low", "Med-High", "High")),
        PC2_bin = cut(pca$x[, 2], breaks = 4, labels = c("Low", "Med-Low", "Med-High", "High")),
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    # Count flows
    flow_df <- as.data.frame(table(df$PC1_bin, df$PC2_bin, df$group))
    colnames(flow_df) <- c("PC1_bin", "PC2_bin", "group", "Freq")
    flow_df <- flow_df[flow_df$Freq > 0, ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(flow_df, aes(axis1 = PC1_bin, axis2 = PC2_bin, y = Freq, fill = group)) +
        ggalluvial::geom_alluvium(aes(fill = group), alpha = 0.7, width = 1/12) +
        ggalluvial::geom_flow(aes(fill = group), width = 1/12, alpha = 0.5) +
        scale_fill_manual(values = pal) +
        theme_minimal() +
        labs(title = "Sample Distribution Flow (PC1 → PC2)", 
             x = "PC Bins", 
             y = "Number of Samples",
             fill = group_key)
    
    save_plot("plots/innovative", "alluvial_pca_flow.svg", p, width = 10, height = 6)
}
alluvial_pca_flow("celltype")

# AO) Ridgeline plot showing PC score distributions across groups
ridgeline_pc_distributions <- function(n_pcs = 3, group_key = "celltype") {
    if (!requireNamespace("ggridges", quietly = TRUE)) {
        message("ggridges not available; skipping ridgeline plot")
        return(NULL)
    }
    
    # Prepare data for first N PCs
    pc_data_list <- list()
    for (i in 1:min(n_pcs, ncol(pca$x))) {
        pc_name <- paste0("PC", i)
        df <- data.frame(
            PC = pc_name,
            score = pca$x[, i],
            group = factor(meta[rownames(pca$x), group_key])
        )
        pc_data_list[[i]] <- df
    }
    combined_df <- do.call(rbind, pc_data_list)
    combined_df <- combined_df[!is.na(combined_df$group), ]
    
    pal <- make_modern_palette(nlevels(combined_df$group))
    
    p <- ggplot(combined_df, aes(x = score, y = PC, fill = group)) +
        ggridges::geom_density_ridges(alpha = 0.7, scale = 0.9) +
        scale_fill_manual(values = pal) +
        theme_minimal() +
        labs(title = "PC Score Distributions by Group", 
             x = "Score", y = "Principal Component", fill = group_key)
    
    save_plot("plots/innovative", "ridgeline_pc_distributions.svg", p, width = 10, height = 6)
}
ridgeline_pc_distributions(3, "celltype")

# AP) Hexbin density plot with marginal distributions
hexbin_pca_with_marginals <- function(group_key = "celltype") {
    if (!requireNamespace("ggExtra", quietly = TRUE)) {
        message("ggExtra not available; skipping marginal plot")
        return(NULL)
    }
    
    df <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(PC1, PC2, color = group)) +
        geom_point(size = 4, alpha = 0.7) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "PCA with Marginal Densities", color = group_key)
    
    p_marg <- ggExtra::ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
    
    ggsave(file.path(ensure_dir(subdir("plots/innovative")), "pca_marginal_densities.png"),
           p_marg, width = 10, height = 8, dpi = 300)
}
hexbin_pca_with_marginals("celltype")

# AQ) Network graph of PC loadings (protein-protein similarity)
loading_network_graph <- function(pc_x = 1, pc_y = 2, top_n = 50, cor_threshold = 0.7) {
    if (!requireNamespace("ggraph", quietly = TRUE) || !requireNamespace("igraph", quietly = TRUE)) {
        message("ggraph/igraph not available; skipping network graph")
        return(NULL)
    }
    
    # Get top proteins by loading magnitude
    rot <- as.data.frame(pca$rotation)
    mag <- sqrt(rot[, pc_x]^2 + rot[, pc_y]^2)
    top_proteins <- names(head(sort(mag, decreasing = TRUE), top_n))
    
    # Filter to proteins that exist in mat
    top_proteins <- intersect(top_proteins, rownames(mat))
    
    if (length(top_proteins) < 2) {
        message("Not enough valid proteins found for network graph")
        return(NULL)
    }
    
    # Calculate protein-protein correlations in expression space
    expr_subset <- mat[top_proteins, , drop = FALSE]
    cor_mat <- cor(t(expr_subset), use = "pairwise.complete.obs")
    
    # Create edges for high correlations
    cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
    edges <- which(abs(cor_mat) > cor_threshold, arr.ind = TRUE)
    
    if (nrow(edges) > 0) {
        edge_df <- data.frame(
            from = rownames(cor_mat)[edges[, 1]],
            to = colnames(cor_mat)[edges[, 2]],
            correlation = cor_mat[edges]
        )
        
        g <- igraph::graph_from_data_frame(edge_df, directed = FALSE)
        
        # Node attributes: loading values
        igraph::V(g)$PC1_load <- rot[igraph::V(g)$name, pc_x]
        igraph::V(g)$PC2_load <- rot[igraph::V(g)$name, pc_y]
        
        p <- ggraph::ggraph(g, layout = "fr") +
            ggraph::geom_edge_link(aes(alpha = abs(correlation)), color = "gray70") +
            ggraph::geom_node_point(aes(color = PC1_load, size = abs(PC2_load))) +
            ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
            scale_color_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
            theme_void() +
            labs(title = "Protein Co-expression Network (Top Contributors)",
                 color = paste0("PC", pc_x, " Loading"),
                 size = paste0("PC", pc_y, " |Loading|"))
        
        save_plot("plots/innovative", "loading_network_graph.svg", p, width = 12, height = 10)
    }
}
loading_network_graph(1, 2, 50, 0.7)

# AR) Parallel coordinates plot for top PCs
parallel_coordinates_pcs <- function(n_pcs = 5, group_key = "celltype") {
    if (!requireNamespace("GGally", quietly = TRUE)) {
        message("GGally not available; skipping parallel coordinates")
        return(NULL)
    }
    
    pc_subset <- pca$x[, 1:min(n_pcs, ncol(pca$x)), drop = FALSE]
    df <- data.frame(
        pc_subset,
        group = factor(meta[rownames(pc_subset), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- GGally::ggparcoord(df, columns = 1:n_pcs, groupColumn = "group",
                           scale = "std", alphaLines = 0.6) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "Parallel Coordinates: Sample Profiles Across PCs", 
             color = group_key)
    
    save_plot("plots/innovative", "parallel_coordinates_pcs.svg", p, width = 10, height = 6)
}
parallel_coordinates_pcs(5, "celltype")

# AS) Sunburst/treemap of variance contribution
variance_sunburst <- function() {
    if (!requireNamespace("treemap", quietly = TRUE) || !requireNamespace("d3r", quietly = TRUE)) {
        message("treemap/d3r not available; trying basic treemap")
        var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100
        n_show <- min(20, length(var_exp))
        
        df <- data.frame(
            PC = paste0("PC", 1:n_show),
            Category = ifelse(1:n_show <= 5, "Top 5", 
                            ifelse(1:n_show <= 10, "PC 6-10", "PC 11-20")),
            Variance = var_exp[1:n_show]
        )
        
        p <- ggplot(df, aes(area = Variance, fill = Category, label = PC)) +
            treemapify::geom_treemap() +
            treemapify::geom_treemap_text(color = "white", place = "centre") +
            scale_fill_manual(values = c("#6B5B95", "#88B04B", "#B2B2B2")) +
            theme_minimal() +
            labs(title = "Variance Contribution Treemap")
        
        save_plot("plots/innovative", "variance_treemap.svg", p, width = 10, height = 8)
        return(NULL)
    }
}
variance_sunburst()

# AT) Contour plot showing density in PC space
contour_density_pca <- function(group_key = "celltype") {
    df <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = factor(meta[rownames(pca$x), group_key])
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(df$group))
    
    p <- ggplot(df, aes(PC1, PC2, color = group)) +
        geom_density_2d(size = 1) +
        geom_point(alpha = 0.6, size = 3) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        labs(title = "PCA with Density Contours", color = group_key)
    
    save_plot("plots/innovative", "contour_density_pca.svg", p, width = 10, height = 8)
}
contour_density_pca("celltype")

# AU) Chord diagram for group relationships based on PC overlap
chord_diagram_groups <- function(group_key = "celltype", pc_threshold = 2) {
    if (!requireNamespace("circlize", quietly = TRUE)) {
        message("circlize not available; skipping chord diagram")
        return(NULL)
    }
    
    # Calculate group centroids in PC space
    groups <- factor(meta[rownames(pca$x), group_key])
    groups <- groups[!is.na(groups)]
    
    centroids <- aggregate(pca$x[names(groups), 1:5], 
                          by = list(group = groups), 
                          FUN = mean)
    rownames(centroids) <- centroids$group
    centroids$group <- NULL
    
    # Calculate pairwise distances between group centroids
    dist_mat <- as.matrix(dist(centroids))
    
    # Convert to flow matrix (inverse of distance)
    flow_mat <- 1 / (1 + dist_mat)
    diag(flow_mat) <- 0
    
    png(file.path(ensure_dir(subdir("plots/innovative")), "chord_diagram_groups.png"),
        width = 10, height = 10, units = "in", res = 300)
    
    circlize::chordDiagram(flow_mat, 
                          grid.col = make_modern_palette(nrow(flow_mat)),
                          transparency = 0.5)
    title("Group Relationships in PC Space")
    
    dev.off()
    circlize::circos.clear()
}
chord_diagram_groups("celltype")

# AV) Small multiples: individual sample trajectories
sample_trajectory_multiples <- function(n_samples = 16, group_key = "celltype") {
    # Plot all samples, faceted by group
    groups <- factor(meta[rownames(pca$x), group_key])
    
    df <- data.frame(
        sample_id = rownames(pca$x),
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = groups,
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$group), ]
    
    pal <- make_modern_palette(nlevels(droplevels(df$group)))
    
    p <- ggplot(df, aes(x = PC1, y = PC2, color = group)) +
        geom_point(size = 3, alpha = 0.8) +
        facet_wrap(~ group, ncol = 3) +
        scale_color_manual(values = pal) +
        theme_pca_min() +
        theme(strip.text = element_text(size = 10, face = "bold")) +
        labs(title = "Sample Positions in PC Space by Group", color = group_key)
    
    save_plot("plots/innovative", "sample_multiples.svg", p, width = 12, height = 10)
}
sample_trajectory_multiples(16, "celltype")

message("All innovative visualizations completed!")

# hexbin_umap_with_marginals - UMAP version with marginal distributions
hexbin_umap_with_marginals <- function(group_key = "celltype") {
    if (is.null(um)) {
        message("UMAP not available; skipping UMAP marginal plot")
        return(NULL)
    }
    if (!requireNamespace("ggExtra", quietly = TRUE)) {
        message("ggExtra not available; skipping marginal plot")
        return(NULL)
    }
    
    # Recreate um_df to ensure consistency with original UMAP plots
    um_df <- data.frame(UMAP1 = um[,1], UMAP2 = um[,2], meta[rownames(um), , drop = FALSE])
    
    # Filter by group_key
    if (!group_key %in% names(um_df)) {
        message(paste("Group key", group_key, "not found in metadata"))
        return(NULL)
    }
    
    grp <- factor(trim_ws(as.character(um_df[[group_key]])))
    um_df$group <- grp
    um_df <- um_df[!is.na(um_df$group), ]
    
    pal <- make_modern_palette(nlevels(droplevels(um_df$group)))
    
    p <- ggplot(um_df, aes(UMAP1, UMAP2, color = group)) +
        geom_point(size = 6, shape = 16, alpha = 0.8, stroke = 0) +
        scale_color_manual(values = pal, na.translate = FALSE) +
        theme_pca_min() +
        theme(legend.position = "bottom") +
        labs(title = paste("UMAP by", group_key, "with Marginal Densities"), color = group_key)
    
    p_marg <- ggExtra::ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE, size = 10)
    
    ggsave(file.path(ensure_dir(subdir("plots/umap")), 
                     paste0("umap_marginal_densities_by_", group_key, ".svg")),
           p_marg, width = 10, height = 8, dpi = 300)
}
hexbin_umap_with_marginals("celltype")
hexbin_umap_with_marginals("ExpGroup")
hexbin_umap_with_marginals("ReplicateGroup")

