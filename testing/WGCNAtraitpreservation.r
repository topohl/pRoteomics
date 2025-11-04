# ================================
# WGCNA with rich traits + preservation
# ================================

# Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db", ask = FALSE, update = FALSE)

suppressPackageStartupMessages(
  pacman::p_load(WGCNA, flashClust, curl, readxl, ggplot2, svglite, GO.db, install = TRUE)
)
allowWGCNAThreads()

# --------------------------
# Paths and data load
# --------------------------
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

expr_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/male.data.xlsx"
meta_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info.xlsx"

male.data <- read_excel(path = expr_xlsx)
# Keep mouse proteins
male.data <- male.data[grepl("_MOUSE", male.data$gene_symbol), ]
# Expression matrix: rows genes, cols samples
expression.data <- male.data[, -1]
expression.data <- as.data.frame(lapply(expression.data, as.numeric))
expression.data <- as.data.frame(t(expression.data))
names(expression.data) <- male.data$gene_symbol

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
par(cex = 0.6, mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.main = 2); abline(h = 40, col = "red")
svg(file = file.path(output_dir, "sample_clustering_outliers.svg"), width = 8, height = 6)
par(cex = 0.6, mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.main = 2); abline(h = 40, col = "red")
dev.off()

cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
expression.data <- expression.data[cut.sampleTree == 1, ]

# --------------------------
# Soft-threshold selection
# --------------------------
# Soft-threshold selection (improved plotting; opens interactive device when possible)
spt <- pickSoftThreshold(expression.data)

open_plot_device <- function(width = 7, height = 5) {
    if (interactive()) {
        tryCatch({
            if (.Platform$OS.type == "windows") {
                windows(width = width, height = height)
            } else if (capabilities("aqua")) {
                quartz(width = width, height = height)
            } else if (capabilities("X11")) {
                x11(width = width, height = height)
            }
        }, error = function(e) { message("Could not open interactive device: ", e$message) })
    }
}

# 1) Scale independence: open interactive device (if possible) then save SVG
open_plot_device()
par(mar = c(4, 4, 2, 1))
plot(spt$fitIndices[, 1], spt$fitIndices[, 2],
         xlab = "Soft Threshold (power)",
         ylab = "Scale Free Topology Model Fit, signed R^2",
         type = "n", main = "Scale independence")
text(spt$fitIndices[, 1], spt$fitIndices[, 2], labels = spt$fitIndices[, 1], col = "red")
abline(h = 0.80, col = "red")

# save to SVG (always)
out1 <- file.path(output_dir, "soft_threshold_scale_independence.svg")
svglite::svglite(file = out1, width = 7, height = 5)
par(mar = c(4, 4, 2, 1))
plot(spt$fitIndices[, 1], spt$fitIndices[, 2],
         xlab = "Soft Threshold (power)",
         ylab = "Scale Free Topology Model Fit, signed R^2",
         type = "n", main = "Scale independence")
text(spt$fitIndices[, 1], spt$fitIndices[, 2], labels = spt$fitIndices[, 1], col = "red")
abline(h = 0.80, col = "red")
dev.off()

# 2) Mean connectivity: open interactive device (if possible) then save SVG
open_plot_device()
par(mar = c(4, 4, 2, 1))
plot(spt$fitIndices[, 1], spt$fitIndices[, 5],
         xlab = "Soft Threshold (power)",
         ylab = "Mean Connectivity",
         type = "n", main = "Mean connectivity")
text(spt$fitIndices[, 1], spt$fitIndices[, 5], labels = spt$fitIndices[, 1], col = "red")

out2 <- file.path(output_dir, "soft_threshold_mean_connectivity.svg")
svglite::svglite(file = out2, width = 7, height = 5)
par(mar = c(4, 4, 2, 1))
plot(spt$fitIndices[, 1], spt$fitIndices[, 5],
         xlab = "Soft Threshold (power)",
         ylab = "Mean Connectivity",
         type = "n", main = "Mean connectivity")
text(spt$fitIndices[, 1], spt$fitIndices[, 5], labels = spt$fitIndices[, 1], col = "red")
dev.off()

message("Saved SVGs:\n - ", out1, "\n - ", out2, "\nIf you expected an interactive plot in VSCode, ensure the R extension's Plot pane is enabled; the script also opens an interactive device when possible.")

softPower <- 6

# --------------------------
# Network construction
# --------------------------
adjacency <- adjacency(expression.data, power = softPower, type = "signed", corFnc = "bicor", corOptions = list(use="p", maxPOutliers=0.05))
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

svg(file = file.path(output_dir, "gene_dendrogram.svg"), width = 12, height = 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
ModuleColors <- labels2colors(Modules)

svg(file = file.path(output_dir, "gene_dendrogram_module_colors.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, ModuleColors, "Module", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

MElist <- moduleEigengenes(expression.data, colors = ModuleColors)
MEs <- MElist$eigengenes
ME.dissimilarity <- 1 - cor(MEs, use = "p", method = "pearson")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
mergedColors <- merge$colors
mergedMEs <- orderMEs(merge$newMEs)

svg(file = file.path(output_dir, "gene_dendrogram_modules_merged.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), c("Original Module","Merged Module"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

# --------------------------
# Block A: rich trait matrix
# --------------------------
sample_info <- read_excel(path = meta_xlsx)
stopifnot("row.names" %in% names(sample_info))
rownames(sample_info) <- as.character(sample_info$row.names)

# Align metadata to expression samples
Samples <- rownames(expression.data)
sample_info <- sample_info[Samples, , drop = FALSE]

# Build design matrices (no intercept)
X_celltype <- model.matrix(~ 0 + celltype, data = sample_info)
colnames(X_celltype) <- sub("^celltype", "celltype_", colnames(X_celltype))

X_layer <- model.matrix(~ 0 + layer, data = sample_info)
colnames(X_layer) <- sub("^layer", "layer_", colnames(X_layer))

X_region <- model.matrix(~ 0 + region, data = sample_info)
colnames(X_region) <- sub("^region", "region_", colnames(X_region))

X_cond <- model.matrix(~ 0 + ExpGroup, data = sample_info)  # con/res/sus one-hot
colnames(X_cond) <- sub("^ExpGroup", "cond_", colnames(X_cond))

# Optional biologically meaningful interactions
X_CA_neuron <- cbind(
  CA1_neuron = as.numeric(sample_info$region == "CA1" & sample_info$celltype == "neuron"),
  CA2_neuron = as.numeric(sample_info$region == "CA2" & sample_info$celltype == "neuron"),
  CA3_neuron = as.numeric(sample_info$region == "CA3" & sample_info$celltype == "neuron"),
  DG_neuron  = as.numeric(sample_info$region == "DG"  & sample_info$celltype == "neuron")
)

datTraits <- as.data.frame(cbind(X_celltype, X_layer, X_region, X_cond, X_CA_neuron), stringsAsFactors = FALSE)
# Drop constant/NA columns
keep_cols <- vapply(datTraits, function(x) sd(as.numeric(x), na.rm = TRUE) > 0, logical(1))
datTraits <- datTraits[, keep_cols, drop = FALSE]
rownames(datTraits) <- Samples

# Eigengene–trait correlations
nSamples <- nrow(expression.data)
MEcorr <- cor(mergedMEs, datTraits, use = "p", method = "pearson")
MEp <- corPvalueStudent(MEcorr, nSamples)

# Heatmap per trait family (split by prefix to keep readable)
plot_trait_heatmap <- function(matCorr, matP, cols, file) {
  textMatrix <- paste(signif(matCorr[, cols, drop=FALSE], 2), "\n(", signif(matP[, cols, drop=FALSE], 1), ")", sep = "")
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
  layer    = grep("^layer_", trait_names),
  region   = grep("^region_", trait_names),
  cond     = grep("^cond_", trait_names),
  combo    = grep("^(CA[123]_neuron)$", trait_names)
)
for (nm in names(groups)) {
  idx <- groups[[nm]]
  if (length(idx) > 0) {
    plot_trait_heatmap(MEcorr, MEp, idx, file.path(output_dir, paste0("ME_trait_heatmap_", nm, ".svg")))
  }
}

# --------------------------
# kME, GS, and exports
# --------------------------
modNames <- substring(colnames(mergedMEs), 3)
geneModuleMembership <- as.data.frame(cor(expression.data, mergedMEs, use = "p", method = "pearson"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

# Example GS with condition numeric score (con=1,res=2,sus=3)
cond_numeric <- setNames(c(1,2,3), c("con","res","sus"))
ExpGroup_num <- as.numeric(cond_numeric[as.character(sample_info$ExpGroup)])
geneTraitSignificance <- as.data.frame(cor(expression.data, ExpGroup_num, use = "p", method = "pearson"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- "GS.ExpGroup"
names(GSPvalue) <- "p.GS.ExpGroup"

# Export genes per module
modules_of_interest <- unique(mergedColors)
for (module in modules_of_interest) {
  moduleGenes <- mergedColors == module
  if (!any(moduleGenes)) next
  mmcol <- paste0("MM", module)
  if (!mmcol %in% colnames(geneModuleMembership)) next
  gene_info <- data.frame(
    Gene = colnames(expression.data)[moduleGenes],
    Module = module,
    ModuleMembership = geneModuleMembership[moduleGenes, mmcol],
    GeneSignificance = geneTraitSignificance[moduleGenes, "GS.ExpGroup"]
  )
  write.csv(gene_info, file = file.path(output_dir, paste0("genes_in_module_", module, ".csv")), row.names = FALSE)
}

# Top hubs per module
top_hubs <- lapply(sort(unique(mergedColors)), function(mod) {
  sel <- which(mergedColors == mod)
  mmcol <- paste0("MM", mod)
  if (!mmcol %in% colnames(geneModuleMembership)) return(NULL)
  df <- data.frame(Gene = colnames(expression.data)[sel],
                   kME  = geneModuleMembership[sel, mmcol],
                   GS   = geneTraitSignificance[sel, "GS.ExpGroup"],
                   stringsAsFactors = FALSE)
  df[order(-abs(df$kME)), ][1:min(10, nrow(df)), ]
})
names(top_hubs) <- sort(unique(mergedColors))
write.csv(do.call(rbind, lapply(names(top_hubs), function(m) cbind(Module=m, top_hubs[[m]]))),
          file.path(output_dir, "top_hubs_per_module.csv"), row.names = FALSE)

# --------------------------
# Block B (fixed): module preservation across conditions
# --------------------------

# Split indices by ExpGroup
idx_con <- which(sample_info$ExpGroup == "con")
idx_res <- which(sample_info$ExpGroup == "res")
idx_sus <- which(sample_info$ExpGroup == "sus")

# Build raw multi-set list
multiExpr <- list(
  ALL = list(data = expression.data),
  CON = list(data = expression.data[idx_con, , drop = FALSE]),
  RES = list(data = expression.data[idx_res, , drop = FALSE]),
  SUS = list(data = expression.data[idx_sus, , drop = FALSE])
)

# 1) Find genes that are present and good in all sets
good_cols <- lapply(multiExpr, function(e) {
  gsg <- goodSamplesGenes(e$data, verbose = 3)
  which(gsg$goodGenes)
})
common_idx <- Reduce(intersect, good_cols)
common_genes <- colnames(multiExpr[[1]]$data)[common_idx]

# Drop genes with any NA in any set
has_na <- sapply(common_genes, function(g)
  any(vapply(multiExpr, function(e) any(is.na(e$data[, g])), logical(1)))
)
common_genes <- common_genes[!has_na]
stopifnot(length(common_genes) > 0)

# 2) Clean each set: subset common genes, numeric matrix, drop zero-variance samples
multi_expr_clean <- lapply(multiExpr, function(e) {
  dat <- e$data[, common_genes, drop = FALSE]
  dat <- as.data.frame(dat)
  dat[] <- lapply(dat, function(col) as.numeric(col))
  keep_samp <- apply(dat, 1, function(r) sd(r, na.rm = TRUE) > 0)
  dat <- dat[keep_samp, , drop = FALSE]
  list(data = as.matrix(dat))
})
names(multi_expr_clean) <- names(multiExpr)

# Ensure each set has >= 2 samples
if (any(vapply(multi_expr_clean, function(e) nrow(e$data) < 2, logical(1)))) {
  stop("One or more sets have < 2 samples after cleaning; reduce filtering or combine sets.")
}

# 3) Build color list aligned to common genes
ref_colors <- mergedColors[match(common_genes, colnames(expression.data))]
stopifnot(length(ref_colors) == length(common_genes))

multi_color <- list(
  ALL = ref_colors,
  CON = rep("grey", length(common_genes)),
  RES = rep("grey", length(common_genes)),
  SUS = rep("grey", length(common_genes))
)

# 4) Run modulePreservation with documented corOptions as a string
set.seed(12345)
mp <- modulePreservation(
  multi_expr_clean,
  multi_color,
  referenceNetworks = 1,                # 1 = ALL
  nPermutations = 200,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = "use = 'p', maxPOutliers = 0.05",
  verbose = 3
)

# 5) Extract Zsummary
Z_CON <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.CON
Z_RES <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.RES
Z_SUS <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.SUS

write.csv(Z_CON, file.path(output_dir, "module_preservation_CON_vs_ALL_Zsummary.csv"), row.names = TRUE)
write.csv(Z_RES, file.path(output_dir, "module_preservation_RES_vs_ALL_Zsummary.csv"), row.names = TRUE)
write.csv(Z_SUS, file.path(output_dir, "module_preservation_SUS_vs_ALL_Zsummary.csv"), row.names = TRUE)


# Optional: summary plot of Zsummary thresholds
# Interpretation guideline: Z < 2 not preserved, 2–10 weak/moderate, >10 strong. [Langfelder–Horvath]
# (Create custom plots as needed.)
plotZsummary <- function(Z, title, file) {
  modNames <- rownames(Z)
  Zsummary <- Z[, "Zsummary.pres"]
  svg(file = file, width = 7, height = max(4, length(modNames) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  barplot(Zsummary, names.arg = modNames, las = 2,
          col = "lightblue", main = title,
          ylab = "Zsummary preservation statistic")
  abline(h = c(2, 10), col = c("red", "blue"), lty = 2)
  legend("topright", legend = c("Z=2 (not preserved)", "Z=10 (strongly preserved)"),
         col = c("red", "blue"), lty = 2, bty = "n")
  dev.off()
}
plotZsummary(Z_CON, "Module Preservation: CON vs ALL", file.path(output_dir, "module_preservation_CON_vs_ALL_Zsummary.svg"))
plotZsummary(Z_RES, "Module Preservation: RES vs ALL", file.path(output_dir, "module_preservation_RES_vs_ALL_Zsummary.svg"))
plotZsummary(Z_SUS, "Module Preservation: SUS vs ALL", file.path(output_dir, "module_preservation_SUS_vs_ALL_Zsummary.svg"))

# --------------------------
# Optional: project control modules into RES/SUS
# --------------------------
# If desired, build reference modules only in CON and project eigengenes into RES/SUS:
# 1) Re-run network construction on CON only to get mergedColors_CON and mergedMEs_CON (reference).


# 2) For RES/SUS, compute projected eigengenes by correlating each sample's expression with the reference ME loadings.
# 3) Compare kME distributions and eigengene shifts; complements modulePreservation.










# ================================
# WGCNA with spatial traits + preservation + extras
# ================================

# Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db", ask = FALSE, update = FALSE)

suppressPackageStartupMessages(
  pacman::p_load(WGCNA, flashClust, curl, readxl, ggplot2, svglite, GO.db, install = TRUE)
)
allowWGCNAThreads()

# --------------------------
# Paths and data load
# --------------------------
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

expr_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/male.data.xlsx"
meta_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info.xlsx"

male.data <- read_excel(path = expr_xlsx)
male.data <- male.data[grepl("_MOUSE", male.data$gene_symbol), ]  # keep mouse proteins

# Expression matrix: rows = samples, cols = proteins
expression.data <- male.data[, -1]
expression.data <- as.data.frame(lapply(expression.data, as.numeric))
expression.data <- as.data.frame(t(expression.data))
names(expression.data) <- male.data$gene_symbol

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
par(cex = 0.6, mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.main = 2); abline(h = 40, col = "red")
svg(file = file.path(output_dir, "sample_clustering_outliers.svg"), width = 8, height = 6)
par(cex = 0.6, mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.main = 2); abline(h = 40, col = "red")
dev.off()

cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
expression.data <- expression.data[cut.sampleTree == 1, ]

# --------------------------
# Soft-threshold selection
# --------------------------
spt <- pickSoftThreshold(expression.data)

open_plot_device <- function(width = 7, height = 5) {
  if (interactive()) {
    tryCatch({
      if (.Platform$OS.type == "windows") windows(width = width, height = height)
      else if (capabilities("aqua")) quartz(width = width, height = height)
      else if (capabilities("X11")) x11(width = width, height = height)
    }, error = function(e) message("Could not open interactive device: ", e$message))
  }
}

# Scale independence plot
open_plot_device(); par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
text(spt$fitIndices[,1], spt$fitIndices[,2], labels = spt$fitIndices[,1], col = "red"); abline(h = 0.80, col = "red")
svglite::svglite(file = file.path(output_dir, "soft_threshold_scale_independence.svg"), width = 7, height = 5)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
text(spt$fitIndices[,1], spt$fitIndices[,2], labels = spt$fitIndices[,1], col = "red"); abline(h = 0.80, col = "red"); dev.off()

# Mean connectivity plot
open_plot_device(); par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(spt$fitIndices[,1], spt$fitIndices[,5], labels = spt$fitIndices[,1], col = "red")
svglite::svglite(file = file.path(output_dir, "soft_threshold_mean_connectivity.svg"), width = 7, height = 5)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(spt$fitIndices[,1], spt$fitIndices[,5], labels = spt$fitIndices[,1], col = "red"); dev.off()

softPower <- 6

# --------------------------
# Network construction
# --------------------------
adjacency <- adjacency(expression.data, power = softPower, type = "signed", corFnc = "bicor", corOptions = list(use="p", maxPOutliers=0.05))
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

svg(file = file.path(output_dir, "gene_dendrogram.svg"), width = 12, height = 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
ModuleColors <- labels2colors(Modules)

svg(file = file.path(output_dir, "gene_dendrogram_module_colors.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, ModuleColors, "Module", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

MElist <- moduleEigengenes(expression.data, colors = ModuleColors)
MEs <- MElist$eigengenes
ME.dissimilarity <- 1 - cor(MEs, use = "p", method = "pearson")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
mergedColors <- merge$colors
mergedMEs <- orderMEs(merge$newMEs)

svg(file = file.path(output_dir, "gene_dendrogram_modules_merged.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), c("Original Module","Merged Module"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

# --------------------------
# Eigengene network plots (added)
# --------------------------
MET <- orderMEs(mergedMEs)
svg(file = file.path(output_dir, "eigengene_dendrogram.svg"), width = 8, height = 6)
plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE, marDendro = c(0,4,2,0))
dev.off()
svg(file = file.path(output_dir, "eigengene_adjacency_heatmap.svg"), width = 8, height = 6)
par(mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", plotDendrograms = FALSE, marHeatmap = c(5,5,2,2), xLabelsAngle = 90)
dev.off()

# --------------------------
# Block A: rich trait matrix (spatial + condition)
# --------------------------
sample_info <- read_excel(path = meta_xlsx)
stopifnot("row.names" %in% names(sample_info))
rownames(sample_info) <- as.character(sample_info$row.names)

# Align metadata to expression samples
Samples <- rownames(expression.data)
sample_info <- sample_info[Samples, , drop = FALSE]

# One-hot designs
X_celltype <- model.matrix(~ 0 + celltype, data = sample_info); colnames(X_celltype) <- sub("^celltype", "celltype_", colnames(X_celltype))
X_layer    <- model.matrix(~ 0 + layer,    data = sample_info); colnames(X_layer)    <- sub("^layer",    "layer_",    colnames(X_layer))
X_region   <- model.matrix(~ 0 + region,   data = sample_info); colnames(X_region)   <- sub("^region",   "region_",   colnames(X_region))
X_cond     <- model.matrix(~ 0 + ExpGroup, data = sample_info); colnames(X_cond)     <- sub("^ExpGroup", "cond_",     colnames(X_cond))

# Optional biologically meaningful interactions
X_CA_neuron <- cbind(
  CA1_neuron = as.numeric(sample_info$region == "CA1" & sample_info$celltype == "neuron"),
  CA2_neuron = as.numeric(sample_info$region == "CA2" & sample_info$celltype == "neuron"),
  CA3_neuron = as.numeric(sample_info$region == "CA3" & sample_info$celltype == "neuron"),
  DG_neuron  = as.numeric(sample_info$region == "DG"  & sample_info$celltype == "neuron")
)

datTraits <- as.data.frame(cbind(X_celltype, X_layer, X_region, X_cond, X_CA_neuron), stringsAsFactors = FALSE)
keep_cols <- vapply(datTraits, function(x) sd(as.numeric(x), na.rm = TRUE) > 0, logical(1))
datTraits <- datTraits[, keep_cols, drop = FALSE]
rownames(datTraits) <- Samples

# Eigengene–trait correlations for grouped spatial traits
nSamples <- nrow(expression.data)
MEcorr <- cor(mergedMEs, datTraits, use = "p", method = "pearson")
MEp    <- corPvalueStudent(MEcorr, nSamples)

plot_trait_heatmap <- function(matCorr, matP, cols, file) {
  if (length(cols) == 0) return(invisible(NULL))
  textMatrix <- paste(signif(matCorr[, cols, drop=FALSE], 2), "\n(", signif(matP[, cols, drop=FALSE], 1), ")", sep = "")
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
  cond     = grep("^cond_",    trait_names),
  combo    = grep("^(CA[123]_neuron)$", trait_names)
)
for (nm in names(groups)) {
  idx <- groups[[nm]]
  if (length(idx) > 0) {
    plot_trait_heatmap(MEcorr, MEp, idx, file.path(output_dir, paste0("ME_trait_heatmap_", nm, ".svg")))
  }
}

# --------------------------
# Pairwise condition contrast heatmaps (added)
# --------------------------
mk_contrast <- function(vec, a, b){v<-rep(NA_real_,length(vec));v[vec==a]<-0;v[vec==b]<-1;v}
grp <- as.character(sample_info$ExpGroup)
contrasts <- list(con_res = mk_contrast(grp,"con","res"),
                  con_sus = mk_contrast(grp,"con","sus"),
                  res_sus = mk_contrast(grp,"res","sus"))
for (nm in names(contrasts)) {
  v <- contrasts[[nm]]; keep <- !is.na(v)
  cmat <- cor(mergedMEs[keep, , drop=FALSE], v[keep], use="p")
  pmat <- corPvalueStudent(cmat, sum(keep))
  txt <- paste(signif(cmat,2), "\n(", signif(pmat,1), ")", sep = "")
  dim(txt) <- dim(cmat)
  svg(file = file.path(output_dir, paste0("ME_trait_heatmap_", nm, ".svg")), width = 3, height = max(4, ncol(mergedMEs) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = cmat, xLabels = nm, yLabels = colnames(mergedMEs), ySymbols = colnames(mergedMEs),
                 colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt, setStdMargins = FALSE,
                 cex.text = 0.8, zlim = c(-1,1), main = "Module–trait relationships")
  dev.off()
}

# --------------------------
# kME, GS, hub exports
# --------------------------
modNames <- substring(colnames(mergedMEs), 3)
geneModuleMembership <- as.data.frame(cor(expression.data, mergedMEs, use = "p", method = "pearson"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

cond_numeric <- setNames(c(1,2,3), c("con","res","sus"))
ExpGroup_num <- as.numeric(cond_numeric[as.character(sample_info$ExpGroup)])
geneTraitSignificance <- as.data.frame(cor(expression.data, ExpGroup_num, use = "p", method = "pearson"))
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
  write.csv(gene_info, file = file.path(output_dir, paste0("genes_in_module_", module, ".csv")), row.names = FALSE)
}

top_hubs <- lapply(sort(unique(mergedColors)), function(mod) {
  sel <- which(mergedColors == mod)
  mmcol <- paste0("MM", mod); if (!mmcol %in% colnames(geneModuleMembership)) return(NULL)
  df <- data.frame(Gene = colnames(expression.data)[sel],
                   kME  = geneModuleMembership[sel, mmcol],
                   GS   = geneTraitSignificance[sel, "GS.ExpGroup"],
                   stringsAsFactors = FALSE)
  df[order(-abs(df$kME)), ][1:min(10, nrow(df)), ]
})
names(top_hubs) <- sort(unique(mergedColors))
write.csv(do.call(rbind, lapply(names(top_hubs), function(m) cbind(Module=m, top_hubs[[m]]))),
          file.path(output_dir, "top_hubs_per_module.csv"), row.names = FALSE)

# --------------------------
# Per‑stratum eigengene summaries (added)
# --------------------------
agg_ME <- function(by){
  ag <- aggregate(mergedMEs, by = list(group = by), FUN = mean, na.rm = TRUE)
  rownames(ag) <- ag$group; ag$group <- NULL; ag
}
ME_by_region <- agg_ME(sample_info$region)
ME_by_layer  <- agg_ME(sample_info$layer)
ME_by_cell   <- agg_ME(sample_info$celltype)
write.csv(ME_by_region, file.path(output_dir, "ME_by_region_mean.csv"))
write.csv(ME_by_layer,  file.path(output_dir, "ME_by_layer_mean.csv"))
write.csv(ME_by_cell,   file.path(output_dir, "ME_by_celltype_mean.csv"))

# --------------------------
# Module preservation across conditions (fixed)
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
  verbose = 3
)

Z_CON <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.CON
Z_RES <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.RES
Z_SUS <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.SUS

write.csv(Z_CON, file.path(output_dir, "module_preservation_CON_vs_ALL_Zsummary.csv"), row.names = TRUE)
write.csv(Z_RES, file.path(output_dir, "module_preservation_RES_vs_ALL_Zsummary.csv"), row.names = TRUE)
write.csv(Z_SUS, file.path(output_dir, "module_preservation_SUS_vs_ALL_Zsummary.csv"), row.names = TRUE)

plotZsummary <- function(Z, title, file) {
  modNames <- rownames(Z); Zsummary <- Z[, "Zsummary.pres"]
  svg(file = file, width = 7, height = max(4, length(modNames) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  barplot(Zsummary, names.arg = modNames, las = 2, col = "lightblue", main = title, ylab = "Zsummary preservation statistic")
  abline(h = c(2, 10), col = c("red", "blue"), lty = 2)
  legend("topright", legend = c("Z=2 (not preserved)", "Z=10 (strongly preserved)"), col = c("red", "blue"), lty = 2, bty = "n")
  dev.off()
}
plotZsummary(Z_CON, "Module Preservation: CON vs ALL", file.path(output_dir, "module_preservation_CON_vs_ALL_Zsummary.svg"))
plotZsummary(Z_RES, "Module Preservation: RES vs ALL", file.path(output_dir, "module_preservation_RES_vs_ALL_Zsummary.svg"))
plotZsummary(Z_SUS, "Module Preservation: SUS vs ALL", file.path(output_dir, "module_preservation_SUS_vs_ALL_Zsummary.svg"))

# --------------------------
# Optional scaffolds: preservation by region/layer
# --------------------------
preserve_by <- function(fac, label) {
  lev <- levels(factor(fac))
  sets <- lapply(lev, function(lvl){
    idx <- which(fac == lvl)
    if (length(idx) >= 2) list(name = lvl, data = expression.data[idx, , drop = FALSE]) else NULL
  })
  sets <- Filter(Negate(is.null), sets)
  if (length(sets) < 2) return(invisible(NULL))
  multiExpr_reg <- c(list(ALL = list(data = expression.data)),
                     setNames(lapply(sets, `[[`, "data"), vapply(sets, `[[`, "", "name")))
  # Reuse the cleaning and color-alignment steps above, then call modulePreservation(...)
  # Save Zsummary CSVs named with the 'label' and level.
}
# Example usage (uncomment if enough samples per stratum):
# preserve_by(sample_info$region, "region")
# preserve_by(sample_info$layer,  "layer")

# --------------------------
# Optional: MM–GS scatter for selected modules
# --------------------------
plot_MM_GS <- function(module){
  modNames <- sub("^ME","",colnames(mergedMEs)); col_idx <- match(module, modNames); if (is.na(col_idx)) return()
  moduleGenes <- mergedColors == module; if (!any(moduleGenes)) return()
  mm <- abs(as.numeric(geneModuleMembership[moduleGenes, paste0("MM", module)]))
  gs <- abs(as.numeric(geneTraitSignificance[moduleGenes, 1]))
  df <- data.frame(MM=mm, GS=gs)
  p <- ggplot(df, aes(MM, GS)) + geom_point(shape=21, size=3, fill=module, color="white", alpha=0.9, stroke=0.6) +
       geom_smooth(method="lm", se=FALSE, color="black", linetype="dashed") +
       labs(x=paste0("Module Membership (kME) — ", module), y="Gene significance (|GS|)", title=paste("MM vs GS —", module)) +
       theme_classic()
  ggsave(file.path(output_dir, paste0("MM_vs_GS_", module, ".svg")), p, device=svglite::svglite, width=4, height=5, units="in", dpi=300)
}

for (m in unique(mergedColors)) plot_MM_GS(m)

# Example: for (m in unique(mergedColors)) plot_MM_GS(m)









library(ggplot2)

# Plot ME by condition within a chosen stratum (e.g., cell type == "microglia")
plot_ME_by_stratum <- function(module, stratum_vec, stratum_level, cond_vec, MEs, outdir, ylim=NULL) {
  me_col <- module  # module should be like "MEblue"
  keep <- stratum_vec == stratum_level
  df <- data.frame(
    ME  = MEs[keep, me_col, drop=TRUE],
    Cond = factor(cond_vec[keep], levels = c("con","res","sus"))
  )
  p <- ggplot(df, aes(Cond, ME, fill = Cond)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.5) +
    labs(title = paste(me_col, "in", stratum_level), x = NULL, y = "Module eigengene") +
    theme_classic(base_size = 11)
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, paste0(me_col, "_byCond_in_", stratum_level, ".svg")),
         p, device = svglite::svglite, width = 4.5, height = 4.5)
}


# generate ME-by-condition plots for all modules in each cell type
modules <- colnames(mergedMEs)
celltypes <- unique(as.character(sample_info$celltype))
celltypes <- celltypes[!is.na(celltypes)]

outbase <- file.path(output_dir, "plots_ME_by_condition_in_celltypes")
dir.create(outbase, recursive = TRUE, showWarnings = FALSE)

for (ct in celltypes) {
    outdir <- file.path(outbase, ct)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    for (me in modules) {
        tryCatch({
            plot_ME_by_stratum(me, sample_info$celltype, ct, sample_info$ExpGroup, mergedMEs, outdir)
        }, error = function(e) {
            message("Skipping ", me, " in ", ct, " due to error: ", e$message)
        })
    }
}


# Repeat for regions
outdir_regions <- file.path(output_dir, "plots_ME_by_condition_in_regions")
dir.create(outdir_regions, recursive = TRUE, showWarnings = FALSE)

regions <- unique(as.character(sample_info$region))
regions <- regions[!is.na(regions)]

for (reg in regions) {
    safe_reg <- gsub("[^A-Za-z0-9_\\-]", "_", reg)
    od <- file.path(outdir_regions, safe_reg)
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    for (me in colnames(mergedMEs)) {
        tryCatch({
            plot_ME_by_stratum(me, sample_info$region, reg, sample_info$ExpGroup, mergedMEs, od)
        }, error = function(e) {
            message("Skipping ", me, " in ", reg, " due to error: ", e$message)
        })
    }
}


plot_ME_grid_by_stratum <- function(module, stratum_vec, cond_vec, MEs, outdir, facet_name) {
  me_col <- module
  df <- data.frame(
    ME    = MEs[, me_col],
    Cond  = factor(cond_vec, levels = c("con","res","sus")),
    Group = factor(stratum_vec)
  )
  p <- ggplot(df, aes(Cond, ME, fill = Cond)) +
    geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
    geom_boxplot(width = 0.12, outlier.size = 0.3) +
    facet_wrap(~ Group, scales = "free_y") +
    labs(title = paste(me_col, "by condition across", facet_name), x = NULL, y = "Module eigengene") +
    theme_classic(base_size = 10) +
    theme(strip.text = element_text(size = 8))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, paste0(me_col, "_grid_byCond_by_", facet_name, ".svg")),
         p, device = svglite::svglite, width = 10, height = 7)
}


# Generate ME-grid plots for all modules (not just blue/turquoise) across celltype, region, and layer
modules_all <- colnames(mergedMEs)
outbase <- file.path(output_dir, "plots_ME_grids")
dirs <- list(
    celltype = file.path(outbase, "by_celltype"),
    region   = file.path(outbase, "by_region"),
    layer    = file.path(outbase, "by_layer")
)
lapply(dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE))

for (me in modules_all) {
    # by celltype
    tryCatch({
        plot_ME_grid_by_stratum(me, sample_info$celltype, sample_info$ExpGroup, mergedMEs, dirs$celltype, "celltype")
    }, error = function(e) message("Skipping ", me, " for celltype: ", e$message))
    # by region
    tryCatch({
        plot_ME_grid_by_stratum(me, sample_info$region, sample_info$ExpGroup, mergedMEs, dirs$region, "region")
    }, error = function(e) message("Skipping ", me, " for region: ", e$message))
    # by layer
    tryCatch({
        plot_ME_grid_by_stratum(me, sample_info$layer, sample_info$ExpGroup, mergedMEs, dirs$layer, "layer")
    }, error = function(e) message("Skipping ", me, " for layer: ", e$message))
}







library(reshape2)

# Compute p-values for condition and interaction within each stratum
# and add FDR-adjusted p-values (BH) per effect.
test_me_interactions <- function(mes, cond, stratum) {
    mods <- colnames(mes)

    get_pval <- function(tab, pattern) {
        if (is.null(tab)) return(NA_real_)
        ridx <- which(grepl(pattern, rownames(tab), perl = TRUE))
        if (length(ridx) == 0) return(NA_real_)
        pcol <- grep("Pr\\(", colnames(tab), value = TRUE)
        if (length(pcol) == 0) return(NA_real_)
        as.numeric(tab[ridx[1], pcol[1]])
    }

    out <- lapply(mods, function(me) {
        df <- data.frame(
            ME = mes[, me],
            Cond = factor(cond, levels = c("con", "res", "sus")),
            Stratum = factor(stratum)
        )

        tab <- tryCatch({
            anova(lm(ME ~ Cond * Stratum, data = df))
        }, error = function(e) NULL)

        p_cond <- get_pval(tab, "\\bCond\\b")
        p_int  <- get_pval(tab, "Cond:Stratum")

        data.frame(
            module = me,
            p_condition = p_cond,
            p_interaction = p_int,
            stringsAsFactors = FALSE
        )
    })

    do.call(rbind, out)
}

# Build tile plots for a chosen stratum factor and use BH-adjusted p-values
plot_interaction_tiles <- function(mes, cond, stratum, label, outdir) {
    tbl <- test_me_interactions(mes, cond, stratum)

    mt <- reshape2::melt(
        tbl,
        id.vars = "module",
        variable.name = "effect",
        value.name = "p"
    )

    # Treat missing p-values as non-significant (p = 1)
    mt$p[is.na(mt$p)] <- 1

    # Adjust p-values per effect (BH / FDR)
    mt$p_adj <- ave(mt$p, mt$effect, FUN = function(x) p.adjust(x, method = "BH"))

    # Scores: raw and adjusted (use adjusted for plotting)
    mt$score_raw <- -log10(pmax(mt$p, 1e-300))
    mt$score_adj <- -log10(pmax(mt$p_adj, 1e-300))

    # Save table for inspection
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write.csv(mt, file.path(outdir, paste0("ME_condition_interactions_", label, "_pvalues.csv")), row.names = FALSE)

    p <- ggplot(mt, aes(x = effect, y = module, fill = score_adj)) +
        geom_tile(color = "white", size = 0.1) +
        scale_fill_gradient(low = "#f0f9e8", high = "#0868ac") +
        labs(
            title = paste("Condition and interaction effects (BH-adjusted):", label),
            x = NULL, y = NULL, fill = "-log10 FDR"
        ) +
        theme_minimal(base_size = 10) +
        theme(axis.text.y = element_text(size = 7))

    ggsave(
        file.path(outdir, paste0("ME_condition_interactions_", label, ".svg")),
        p,
        device = svglite::svglite,
        width = 6,
        height = 8
    )
}

# Run tiles for celltype / region / layer (uses BH-adjusted p-values)
plot_interaction_tiles(
    mergedMEs,
    sample_info$ExpGroup,
    sample_info$celltype,
    "celltype",
    file.path(output_dir, "plots_interactions")
)
plot_interaction_tiles(
    mergedMEs,
    sample_info$ExpGroup,
    sample_info$region,
    "region",
    file.path(output_dir, "plots_interactions")
)
plot_interaction_tiles(
    mergedMEs,
    sample_info$ExpGroup,
    sample_info$layer,
    "layer",
    file.path(output_dir, "plots_interactions")
)


# Example usage:
# plot_interaction_tiles(mergedMEs, sample_info$ExpGroup, sample_info$celltype,
#                        "celltype", file.path(output_dir, "plots_interactions"))
# plot_interaction_tiles(mergedMEs, sample_info$ExpGroup, sample_info$region,
#                        "region",   file.path(output_dir, "plots_interactions"))
# plot_interaction_tiles(mergedMEs, sample_info$ExpGroup, sample_info$layer,
#                        "layer",    file.path(output_dir, "plots_interactions"))
















library(ggplot2)
library(patchwork)
library(reshape2)
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
  df <- melt(r, varnames = c("module","trait"), value.name = "r")
  p_long <- melt(p, varnames = c("module","trait"), value.name = "p")
  df$p <- p_long$p
  df$stars <- gtools::stars.pval(df$p)
  df$block <- block_name
  # consistent module order
  df$module <- factor(df$module, levels = rownames(r))
  df
}

dfs <- Filter(Negate(is.null),
              mapply(block_cor_df, block_list, names(block_list), SIMPLIFY = FALSE))

panel_plot <- function(dfi, legend = "none") {
  ggplot(dfi, aes(x = trait, y = module, fill = r)) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes(label = stars), size = 3, color = "black") +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = "#2166ac", mid = "white", high = "#b2182b") +
    labs(x = NULL, y = NULL, fill = "Pearson r") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = legend,
          panel.grid = element_blank(),
          plot.margin = margin(2, 2, 2, 2))
}

plots <- lapply(dfs, panel_plot)
# Shared legend from the first panel
p_legend <- panel_plot(dfs[[1]], legend = "right")
legend_only <- cowplot::get_legend(p_legend)

combo <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(title = "Module–trait relationships: cell type, region, layer, condition")

svg(file.path(output_dir, "panel_module_trait_relationships_no_cellclass.svg"), width = 14, height = 4.5)
cowplot::plot_grid(combo, legend_only, rel_widths = c(1, 0.08))
dev.off()


















library(ggplot2)
library(reshape2)
library(gtools)
library(patchwork)
library(cowplot)

# 1) Build combined condition × region × layer × celltype label
region  <- as.character(sample_info$region)
layer   <- as.character(sample_info$layer)
cell    <- as.character(sample_info$celltype)
cond    <- as.character(sample_info$ExpGroup)   # con/res/sus

combo <- paste(cond, region, layer, cell, sep = "_")
combo <- factor(combo)  # e.g., "sus_CA1_slm_neuron", "con_CA1_none_microglia"

# 2) One-hot encode all present combinations
X_combo <- model.matrix(~ 0 + combo)
colnames(X_combo) <- levels(combo)  # no "comb" prefix at all

# 3) Correlate module eigengenes with combined strata
MEcorr_combo <- cor(mergedMEs, X_combo, use = "p")
MEp_combo    <- corPvalueStudent(MEcorr_combo, nrow(expression.data))

# 4) Tidy for ggplot with stars and consistent module order
df_combo <- melt(MEcorr_combo, varnames = c("module","comb"), value.name = "r")
p_combo  <- melt(MEp_combo,    varnames = c("module","comb"), value.name = "p")
df_combo$p <- p_combo$p
df_combo$stars <- gtools::stars.pval(df_combo$p)
df_combo$module <- factor(df_combo$module, levels = rownames(MEcorr_combo))

# 5) Split columns by condition to get grouped panels (con | res | sus) with same module order
split_by_condition <- function(df) {
  # extract condition prefix from combination label
  df$condition <- sub("_.*$", "", as.character(df$comb))
  df$trait     <- sub("^[^_]+_", "", as.character(df$comb))  # region_layer_celltype part
  split(df, df$condition)
}
by_cond <- split_by_condition(df_combo)

# 6) Single-panel tile plot styled like the example (Pearson r with stars, tight tiles)
panel_plot <- function(dfi, panel_title = "") {
  ggplot(dfi, aes(x = trait, y = module, fill = r)) +
    geom_tile(color = "white", size = 0.3) +
    geom_text(aes(label = stars), size = 3, color = "black") +
    scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                         low = "#2166ac", mid = "white", high = "#b2182b") +
    labs(title = panel_title, x = NULL, y = NULL, fill = "Pearson r") +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2),
      legend.position = "none"
    )
}

# 7) Build ordered panels con | res | sus with shared legend
conds <- c("con","res","sus")
plots <- lapply(conds, function(cn) {
  d <- by_cond[[cn]]
  if (is.null(d)) return(NULL)
  # Keep a readable order of trait columns: within a condition, group by region, then layer, then celltype
  # Reorder traits: region alphabetic, then layer order user-defined, then cell type
  parse_parts <- function(tr) {
    parts <- strsplit(tr, "_")[[1]]
    # parts: region, layer, celltype
    if (length(parts) != 3) parts <- c(parts, rep("", 3 - length(parts)))
    parts
  }
  parts <- do.call(rbind, lapply(as.character(d$trait), parse_parts))
  colnames(parts) <- c("region","layer","cell")
  ord <- order(parts[, "region"], parts[, "layer"], parts[, "cell"])
  d$trait <- factor(d$trait, levels = unique(d$trait[ord]))
  panel_plot(d, toupper(cn))
})

# 1) Remove NULL entries and check
plots <- Filter(Negate(is.null), plots)
if (length(plots) == 0) stop("No condition panels were created; check that 'combo' has levels for con/res/sus and that by_cond split returned data.")

# 2) Build a standalone legend using a dummy data frame
legend_df <- data.frame(x = 1, y = 1, r = c(-1, 0, 1))
p_legend <- ggplot(legend_df, aes(x, y, fill = r)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1), oob = scales::squish,
                       low = "#2166ac", mid = "white", high = "#b2182b") +
  theme_void() +
  theme(legend.position = "right") +
  labs(fill = "Pearson r")

legend_only <- cowplot::get_legend(p_legend)

# 3) Combine panels horizontally
combo_row <- wrap_plots(plots, nrow = 1, guides = "collect")

# 4) Safely assemble with legend; if legend is NULL, just save combo_row
out_file <- file.path(output_dir, "panel_ME_vs_condition_region_layer_celltype.svg")
svglite::svglite(out_file, width = 14, height = 4.5)
if (!is.null(legend_only)) {
  cowplot::plot_grid(combo_row, legend_only, rel_widths = c(1, 0.08))
} else {
  print(combo_row)
}
dev.off()













# ================================ currently final
# WGCNA with spatial traits + preservation + condition×region×layer×celltype panels
# ================================ currently final

# Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db", ask = FALSE, update = FALSE)
suppressPackageStartupMessages(
  pacman::p_load(WGCNA, flashClust, curl, readxl, ggplot2, svglite, GO.db,
                 reshape2, gtools, patchwork, cowplot, install = TRUE)
)
allowWGCNAThreads()

# --------------------------
# Paths and data load
# --------------------------
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

expr_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/male.data.xlsx"
meta_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info.xlsx"

male.data <- read_excel(path = expr_xlsx)
male.data <- male.data[grepl("_MOUSE", male.data$gene_symbol), ]

# Expression matrix: rows = samples, cols = proteins
expression.data <- male.data[, -1]
expression.data <- as.data.frame(lapply(expression.data, as.numeric))
expression.data <- as.data.frame(t(expression.data))
names(expression.data) <- male.data$gene_symbol

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
par(cex = 0.6, mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.main = 2); abline(h = 40, col = "red")
svg(file = file.path(output_dir, "sample_clustering_outliers.svg"), width = 8, height = 6)
dev.off()

#cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
expression.data <- expression.data[cut.sampleTree == 1, ]

# --------------------------
# Soft-threshold selection
# --------------------------
spt <- pickSoftThreshold(expression.data)

svglite::svglite(file = file.path(output_dir, "soft_threshold_scale_independence.svg"), width = 7, height = 5)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
text(spt$fitIndices[,1], spt$fitIndices[,2], labels = spt$fitIndices[,1], col = "red"); abline(h = 0.80, col = "red")
dev.off()

svglite::svglite(file = file.path(output_dir, "soft_threshold_mean_connectivity.svg"), width = 7, height = 5)
par(mar = c(4,4,2,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(spt$fitIndices[,1], spt$fitIndices[,5], labels = spt$fitIndices[,1], col = "red")
dev.off()

softPower <- 6

# --------------------------
# Network construction
# --------------------------
adjacency <- adjacency(expression.data, power = softPower, type = "signed", corFnc = "bicor", corOptions = list(use="p", maxPOutliers=0.05))
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

svg(file = file.path(output_dir, "gene_dendrogram.svg"), width = 12, height = 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
ModuleColors <- labels2colors(Modules)

svg(file = file.path(output_dir, "gene_dendrogram_module_colors.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, ModuleColors, "Module", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

MElist <- moduleEigengenes(expression.data, colors = ModuleColors)
MEs <- MElist$eigengenes
ME.dissimilarity <- 1 - cor(MEs, use = "p", method = "pearson")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
mergedColors <- merge$colors
mergedMEs <- orderMEs(merge$newMEs)

svg(file = file.path(output_dir, "gene_dendrogram_modules_merged.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), c("Original Module","Merged Module"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

# --------------------------
# Eigengene network plots
# --------------------------
MET <- orderMEs(mergedMEs)
svg(file = file.path(output_dir, "eigengene_dendrogram.svg"), width = 8, height = 6)
plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE, marDendro = c(0,4,2,0))
dev.off()
svg(file = file.path(output_dir, "eigengene_adjacency_heatmap.svg"), width = 8, height = 6)
par(mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", plotDendrograms = FALSE, marHeatmap = c(5,5,2,2), xLabelsAngle = 90)
dev.off()

# --------------------------
# Spatial + condition trait matrix and heatmaps
# --------------------------
sample_info <- read_excel(path = meta_xlsx)
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
  textMatrix <- paste(signif(matCorr[, cols, drop=FALSE], 2), "\n(", signif(matP[, cols, drop=FALSE], 1), ")", sep = "")
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
    plot_trait_heatmap(MEcorr, MEp, idx, file.path(output_dir, paste0("ME_trait_heatmap_", nm, ".svg")))
  }
}

# --------------------------
# Pairwise condition contrasts (optional)
# --------------------------
mk_contrast <- function(vec, a, b){v<-rep(NA_real_,length(vec));v[vec==a]<-0;v[vec==b]<-1;v}
grp <- as.character(sample_info$ExpGroup)
contrasts <- list(con_res = mk_contrast(grp,"con","res"),
                  con_sus = mk_contrast(grp,"con","sus"),
                  res_sus = mk_contrast(grp,"res","sus"))
for (nm in names(contrasts)) {
  v <- contrasts[[nm]]; keep <- !is.na(v)
  cmat <- cor(mergedMEs[keep, , drop=FALSE], v[keep], use="p")
  pmat <- corPvalueStudent(cmat, sum(keep))
  txt <- paste(signif(cmat,2), "\n(", signif(pmat,1), ")", sep = "")
  dim(txt) <- dim(cmat)
  svg(file = file.path(output_dir, paste0("ME_trait_heatmap_", nm, ".svg")), width = 3, height = max(4, ncol(mergedMEs) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = cmat, xLabels = nm, yLabels = colnames(mergedMEs), ySymbols = colnames(mergedMEs),
                 colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt, setStdMargins = FALSE,
                 cex.text = 0.8, zlim = c(-1,1), main = "Module–trait relationships")
  dev.off()
}

# --------------------------
# kME, GS, hubs
# --------------------------
modNames <- substring(colnames(mergedMEs), 3)
geneModuleMembership <- as.data.frame(cor(expression.data, mergedMEs, use = "p", method = "pearson"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

cond_numeric <- setNames(c(1,2,3), c("con","res","sus"))
ExpGroup_num <- as.numeric(cond_numeric[as.character(sample_info$ExpGroup)])
geneTraitSignificance <- as.data.frame(cor(expression.data, ExpGroup_num, use = "p", method = "pearson"))
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
  write.csv(gene_info, file = file.path(output_dir, paste0("genes_in_module_", module, ".csv")), row.names = FALSE)
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
    verbose = 3
)

Z_CON <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.CON
Z_RES <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.RES
Z_SUS <- mp$preservation$Z$ref.ALL$inColumnsAlsoPresentIn.SUS

write.csv(Z_CON, file.path(output_dir, "module_preservation_CON_vs_ALL_Zsummary.csv"), row.names = TRUE)
write.csv(Z_RES, file.path(output_dir, "module_preservation_RES_vs_ALL_Zsummary.csv"), row.names = TRUE)
write.csv(Z_SUS, file.path(output_dir, "module_preservation_SUS_vs_ALL_Zsummary.csv"), row.names = TRUE)

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
  ggsave(file.path(output_dir, "panel_ME_vs_condition_region_layer_celltype.svg"),
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
  ggsave(file.path(output_dir, "panel_ME_vs_condition_region_layer_celltype.svg"),
         g, device = svglite::svglite, width = 14, height = 4.5)
}

# --------------------------
# End of script
# --------------------------


















































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
                 reshape2, gtools, patchwork, cowplot, install = TRUE)
)

# Parallel setup: prefer physical cores; enable WGCNA threading globally
nCores <- tryCatch({
  pc <- parallel::detectCores(logical = FALSE)
  if (is.na(pc) || pc < 2) 2 else pc
}, error = function(e) 2)
enableWGCNAThreads(nThreads = nCores)  # replaces allowWGCNAThreads(); controls threading globally
WGCNAnThreads()                        # optional: prints current thread count

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

male.data <- readxl::read_excel(path = expr_xlsx)
male.data <- male.data[grepl("_MOUSE", male.data$gene_symbol), ]

# Expression matrix: rows = samples, cols = proteins
expression.data <- male.data[, -1]
expression.data <- as.data.frame(lapply(expression.data, as.numeric))
expression.data <- as.data.frame(t(expression.data))
names(expression.data) <- male.data$gene_symbol

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
par(cex = 0.6, mar = c(0,4,2,0))
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
plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE, marDendro = c(0,4,2,0))
dev.off()
svg(file = fp_net("eigengene_adjacency_heatmap.svg"), width = 8, height = 6)
par(mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", plotDendrograms = FALSE,
                      marHeatmap = c(5,5,2,2), xLabelsAngle = 90)
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
  textMatrix <- paste(signif(matCorr[, cols, drop=FALSE], 2), "\n(", signif(matP[, cols, drop=FALSE], 1), ")", sep = "")
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
contrasts <- list(con_res = mk_contrast(grp,"con","res"),
                  con_sus = mk_contrast(grp,"con","sus"),
                  res_sus = mk_contrast(grp,"res","sus"))
for (nm in names(contrasts)) {
  v <- contrasts[[nm]]; keep <- !is.na(v)
  cmat <- cor(mergedMEs[keep, , drop=FALSE], v[keep], use="p")
  pmat <- corPvalueStudent(cmat, sum(keep))
  txt <- paste(signif(cmat,2), "\n(", signif(pmat,1), ")", sep = "")
  dim(txt) <- dim(cmat)
  svg(file = fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")), width = 3, height = max(4, ncol(mergedMEs) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = cmat, xLabels = nm, yLabels = colnames(mergedMEs), ySymbols = colnames(mergedMEs),
                 colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt, setStdMargins = FALSE,
                 cex.text = 0.8, zlim = c(-1,1), main = "Module–trait relationships")
  dev.off()
}

# --------------------------
# kME, GS, hubs
# --------------------------
modNames <- substring(colnames(mergedMEs), 3)
geneModuleMembership <- as.data.frame(cor(expression.data, mergedMEs, use = "p", method = "pearson"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

cond_numeric <- setNames(c(1,2,3), c("con","res","sus"))
ExpGroup_num <- as.numeric(cond_numeric[as.character(sample_info$ExpGroup)])
geneTraitSignificance <- as.data.frame(cor(expression.data, ExpGroup_num, use = "p", method = "pearson"))
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
# End of script
# --------------------------



# Ensure required objects are defined for the pheatmap block
# Reconstruct df_all / module_levels / traits_in_order / strip_df if they are not present
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
        # mergedColors uses gene-level names; map module names (MEname -> color)
        # module_levels often are like "MEblue" etc. Remove leading "ME" for matching if needed.
        mod_names <- sub("^ME", "", module_levels)
        mod_cols <- mergedColors
        # mergedColors is gene-level mapping; try to use moduleColors vector if present
        if (is.character(mergedColors) && length(names(mergedColors)) == length(mergedColors)) {
            # if mergedColors is named by gene, cannot use here; fall back to simple color palette
            strip_df <- data.frame(module = module_levels, mod_col = rep("grey80", length(module_levels)), stringsAsFactors = FALSE)
        } else {
            strip_df <- data.frame(module = module_levels, mod_col = rep("grey80", length(module_levels)), stringsAsFactors = FALSE)
        }
    } else {
        strip_df <- data.frame(module = module_levels, mod_col = rep("grey80", length(module_levels)), stringsAsFactors = FALSE)
    }
}

# Ensure output_dir exists
if (!exists("output_dir")) output_dir <- getwd()
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Require that the earlier subdirs and helpers exist:
# subdirs$plots_traits and fp_traits <- function(...) file.path(subdirs$plots_traits, ...)

library(pheatmap)
library(reshape2)
library(grid)

# Reconstruct df_all if needed (as in the provided block) ...

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

# 3) Module color mapping (fix)
# Build module-level color mapping from merged module eigengenes
# mergedMEs has columns like "ME<color>" in order; derive color names from those
if (exists("mergedMEs")) {
  me_names <- colnames(mergedMEs)                          # e.g., "MEblue", "MEturquoise"
  mod_colors_levels <- sub("^ME", "", me_names)            # "blue", "turquoise", ...
  # Map rownames(r_mat) (module_levels) to these ME columns
  # Assuming r_mat rownames equal me_names OR can be aligned by removing "ME"
  row_me <- rownames(r_mat)
  row_mod_colors <- sub("^ME", "", row_me)
  # Build a palette lookup from the earlier colorSeq where possible; else use WGCNA standard
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
  # Prefer previously set palette mapping for Modules if available
  # Build a color lookup for module color names
  color_lookup <- colorSeq
  # Fallback: identity map for common WGCNA names if missing
  missing_cols <- setdiff(unique(row_mod_colors), names(color_lookup))
  if (length(missing_cols)) {
    add <- setNames(missing_cols, missing_cols)
    color_lookup <- c(color_lookup, add)
  }
  # Compose annotation row colors
  row_colors_vec <- unname(color_lookup[row_mod_colors])
  # Any NA fallback to "grey80"
  row_colors_vec[is.na(row_colors_vec)] <- "grey80"
  ann_row <- data.frame(ModuleColor = row_colors_vec, row.names = rownames(r_mat))
  ann_colors <- list(ModuleColor = setNames(unique(ann_row$ModuleColor), unique(ann_row$ModuleColor)))
} else {
  # Fallback if mergedMEs not present
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

# 6) Render once
ph <- pheatmap(
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

# 7) Save PNG/PDF/SVG in plots_traits subfolder
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

library(ggplot2)
library(patchwork)
library(reshape2)
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
  df <- melt(r, varnames = c("module","trait"), value.name = "r")
  p_long <- melt(p, varnames = c("module","trait"), value.name = "p")
  df$p <- p_long$p
  df$stars <- gtools::stars.pval(df$p)
  df$block <- block_name
  # consistent module order
  df$module <- factor(df$module, levels = rownames(r))
  df
}

dfs <- Filter(Negate(is.null),
              mapply(block_cor_df, block_list, names(block_list), SIMPLIFY = FALSE))

# panel_plot now returns a combined strip+heat (patchwork) when legend="none",
# or returns the heat ggplot alone when legend != "none" (so the legend can be extracted).
panel_plot <- function(dfi, legend = "none") {
  if (nrow(dfi) == 0) return(NULL)
  # ensure trait is treated as discrete in the given order
  dfi$trait <- factor(dfi$trait, levels = unique(dfi$trait))

  # main heatmap
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

  # if user requested the heatmap with legend, return that ggplot (for legend extraction)
  if (!identical(legend, "none")) return(p_heat)

  # build left color strip: one tile per module in the same order as 'module' factor levels
  modules <- levels(dfi$module)
  # derive color name from module label (e.g. "MEblue" -> "blue")
  mod_cols <- sub("^ME", "", modules)
  strip_df <- data.frame(module = modules, mod_col = mod_cols, stringsAsFactors = FALSE)
  strip_df$module <- factor(strip_df$module, levels = modules)

  p_strip <- ggplot(strip_df, aes(x = 1, y = module, fill = mod_col)) +
    geom_tile() +
    scale_fill_identity() +
    theme_void() +
    theme(plot.margin = margin(2, 0, 2, 2))

  # combine strip + heat using patchwork with narrow strip
  combined <- p_strip + p_heat + plot_layout(widths = c(0.05, 1))
  combined
}

# Build plots list (combined strip+heat)
plots <- lapply(dfs, panel_plot)
# Shared legend from the first panel: build a heat-only ggplot with legend on right and extract
p_legend_heat <- panel_plot(dfs[[1]], legend = "right")
legend_only <- cowplot::get_legend(p_legend_heat)

combo <- wrap_plots(plots, nrow = 1, guides = "collect") +
  plot_annotation(title = "Module–trait relationships: cell type, region, layer, condition")

svg(file.path(output_dir, "panel_module_trait_relationships_no_cellclass.svg"), width = 14, height = 4.5)
cowplot::plot_grid(combo, legend_only, rel_widths = c(1, 0.08))
dev.off()


























































# ================================ parallel-enabled CURRENT FOR SP SG COMPARISONS
# WGCNA with spatial traits + preservation + condition×region×layer×celltype panels
# Single-celltype-safe (celltype one level -> skipped in datTraits)
# Outputs organized into subfolders under output_dir
# ================================ parallel-enabled

# Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("GO.db", quietly = TRUE)) BiocManager::install("GO.db", ask = FALSE, update = FALSE)
suppressPackageStartupMessages(
  pacman::p_load(WGCNA, flashClust, curl, readxl, ggplot2, svglite, GO.db,
                 reshape2, gtools, patchwork, cowplot, pheatmap, install = TRUE)
)

# WGCNA requirements
options(stringsAsFactors = FALSE)
allowWGCNAThreads()  # legacy no-op if using enableWGCNAThreads
nCores <- tryCatch({
  pc <- parallel::detectCores(logical = FALSE)
  if (is.na(pc) || pc < 2) 2 else pc
}, error = function(e) 2)
enableWGCNAThreads(nThreads = nCores)  # controls threading globally
WGCNAnThreads()

# --------------------------
# Paths and data load
# --------------------------
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output_spsg"
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

expr_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/male.data_spsg.xlsx"
meta_xlsx <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info_spsg.xlsx"

# --------------------------
# Load expression and coerce
# --------------------------
male.data <- readxl::read_excel(path = expr_xlsx)
stopifnot("gene_symbol" %in% names(male.data))
male.data <- male.data[grepl("_MOUSE", male.data$gene_symbol), ]

# Expression matrix: rows = samples, cols = proteins
expression.data <- male.data[, -1, drop = FALSE]
expression.data[] <- lapply(expression.data, function(x) as.numeric(x))
expression.data <- as.data.frame(t(expression.data))
colnames(expression.data) <- male.data$gene_symbol

# Sanity checks
stopifnot(nrow(expression.data) >= 2, ncol(expression.data) >= 2)

# --------------------------
# QC and sample clustering
# --------------------------
gsg <- goodSamplesGenes(expression.data, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) printFlush(paste("Removing genes:", paste(colnames(expression.data)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0) printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", ")))
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE, drop = FALSE]
}
stopifnot(nrow(expression.data) >= 2, ncol(expression.data) >= 2)

# Sample dendrogram
sampleTree <- hclust(dist(expression.data), method = "average")
save_svg(fp_qc("sample_clustering_outliers.svg"), 8, 6, {
  par(cex = 0.6, mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.main = 2)
  abline(h = 40, col = "red")
})

# --------------------------
# Soft-threshold selection
# --------------------------
spt <- pickSoftThreshold(expression.data)  # uses Pearson by default
save_svg(fp_qc("soft_threshold_scale_independence.svg"), 7, 5, {
  par(mar = c(4,4,2,1))
  plot(spt$fitIndices[,1], spt$fitIndices[,2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = "Scale independence")
  text(spt$fitIndices[,1], spt$fitIndices[,2], labels = spt$fitIndices[,1], col = "red"); abline(h = 0.80, col = "red")
})
save_svg(fp_qc("soft_threshold_mean_connectivity.svg"), 7, 5, {
  par(mar = c(4,4,2,1))
  plot(spt$fitIndices[,1], spt$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
  text(spt$fitIndices[,1], spt$fitIndices[,5], labels = spt$fitIndices[,1], col = "red")
})

softPower <- 3

# --------------------------
# Network construction
# --------------------------
adjacency <- adjacency(expression.data, power = softPower, type = "signed",
                       corFnc = "bicor", corOptions = list(use="p", maxPOutliers=0.05))
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1 - TOM
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

save_svg(fp_net("gene_dendrogram.svg"), 12, 9, {
  plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
})

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2,
                         pamRespectsDendro = FALSE, minClusterSize = 30)

# Palette for modules
colorSeq <- c(
  "lemon"="lemonchiffon","sage"="darkseagreen","bluegray"="steelblue","mintblue"="lightsteelblue",
  "azure"="deepskyblue","khaki"="khaki","skyblue"="skyblue","babyblue"="lightblue","amber"="goldenrod",
  "tealgreen"="darkcyan","forestgreen"="forestgreen","gold"="gold","violet"="violet","seafoam"="mediumaquamarine",
  "coral"="coral","salmonlight"="lightsalmon","peach"="peachpuff","mint"="palegreen","lime"="limegreen","mauve"="plum",
  "freesia"="lightpink","cocoa"="saddlebrown","lavender"="lavender","magenta"="magenta","salmon"="salmon","rose"="mistyrose",
  "aquamarine"="aquamarine","tomato"="tomato","plum"="plum","hotpink"="hotpink","rust"="sienna"
)

unique_mods <- sort(unique(Modules))
nmods <- length(unique_mods)
palette_vals <- unname(colorSeq)
if (length(palette_vals) < nmods) palette_vals <- rep(palette_vals, length.out = nmods)
mod_colors_map <- setNames(palette_vals[seq_len(nmods)], as.character(unique_mods))
ModuleColors <- as.character(mod_colors_map[as.character(Modules)])
stopifnot(length(ModuleColors) == length(Modules))

save_svg(fp_net("gene_dendrogram_module_colors.svg"), 12, 9, {
  plotDendroAndColors(geneTree, ModuleColors, "Module", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
})

MElist <- moduleEigengenes(expression.data, colors = ModuleColors)
MEs <- MElist$eigengenes
ME.dissimilarity <- 1 - cor(MEs, use = "p", method = "pearson")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)
mergedColors <- merge$colors
mergedMEs <- orderMEs(merge$newMEs)

save_svg(fp_net("gene_dendrogram_modules_merged.svg"), 12, 9, {
  plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors),
                      c("Original Module","Merged Module"),
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors for original and merged modules")
})

# --------------------------
# Eigengene network plots
# --------------------------
MET <- orderMEs(mergedMEs)
save_svg(fp_net("eigengene_dendrogram.svg"), 8, 6, {
  plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE, marDendro = c(0,4,2,0))
})
save_svg(fp_net("eigengene_adjacency_heatmap.svg"), 8, 6, {
  par(mar = c(1,1,1,1))
  plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", plotDendrograms = FALSE,
                        marHeatmap = c(5,5,2,2), xLabelsAngle = 90)
})

# --------------------------
# Traits (single-celltype safe)
# --------------------------
sample_info <- readxl::read_excel(path = meta_xlsx)
stopifnot("row.names" %in% names(sample_info))
rownames(sample_info) <- as.character(sample_info$row.names)

Samples <- rownames(expression.data)
sample_info <- sample_info[Samples, , drop = FALSE]

# Helper: safe one-hot (returns 0 columns if <2 levels present)
one_hot_or_empty <- function(fac, prefix) {
  fac <- droplevels(as.factor(fac))
  if (length(levels(fac)) < 2) {
    m <- matrix(numeric(0), nrow = length(fac), ncol = 0)
    rownames(m) <- names(fac)
    return(m)
  } else {
    mm <- model.matrix(~ 0 + fac)
    colnames(mm) <- paste0(prefix, levels(fac))
    rownames(mm) <- names(fac)
    return(mm)
  }
}

cell_fac  <- setNames(as.character(sample_info$celltype),  Samples)
layer_fac <- setNames(as.character(sample_info$layer),     Samples)
region_fac<- setNames(as.character(sample_info$region),    Samples)
cond_fac  <- setNames(as.character(sample_info$ExpGroup),  Samples)

X_celltype <- one_hot_or_empty(cell_fac,  "celltype_")  # becomes 0-column if single celltype
X_layer    <- one_hot_or_empty(layer_fac, "layer_")
X_region   <- one_hot_or_empty(region_fac,"region_")
X_cond     <- one_hot_or_empty(cond_fac,  "cond_")

datTraits <- as.data.frame(cbind(X_celltype, X_layer, X_region, X_cond))
if (ncol(datTraits) > 0) {
  keep_cols <- vapply(datTraits, function(x) sd(as.numeric(x), na.rm = TRUE) > 0, logical(1))
  datTraits <- datTraits[, keep_cols, drop = FALSE]
  rownames(datTraits) <- Samples
} else {
  message("datTraits has 0 columns after building trait blocks; ME–trait correlations will be skipped.")
}

# ME–trait correlation heatmaps (skipped if no traits)
if (ncol(datTraits) > 0) {
  nSamples <- nrow(expression.data)
  MEcorr <- cor(mergedMEs, datTraits, use = "p", method = "pearson")
  MEp    <- corPvalueStudent(MEcorr, nSamples)

  plot_trait_heatmap <- function(matCorr, matP, cols, file) {
    if (length(cols) == 0) return(invisible(NULL))
    textMatrix <- paste(signif(matCorr[, cols, drop=FALSE], 2), "\n(", signif(matP[, cols, drop=FALSE], 1), ")", sep = "")
    dim(textMatrix) <- dim(matCorr[, cols, drop=FALSE])
    svglite::svglite(file = file, width = 6, height = max(4, nrow(matCorr) * 0.25 + 2))
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
    celltype = if (length(trait_names)) grep("^celltype_", trait_names) else integer(0),
    layer    = if (length(trait_names)) grep("^layer_",    trait_names) else integer(0),
    region   = if (length(trait_names)) grep("^region_",   trait_names) else integer(0),
    cond     = if (length(trait_names)) grep("^cond_",     trait_names) else integer(0)
  )
  for (nm in names(groups)) {
    idx <- groups[[nm]]
    if (length(idx) > 0) {
      plot_trait_heatmap(MEcorr, MEp, idx, fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")))
    }
  }
} else {
  message("Skipping ME–trait heatmaps: datTraits has 0 columns.")
}

# --------------------------
# Pairwise condition contrasts (optional)
# --------------------------
mk_contrast <- function(vec, a, b){v<-rep(NA_real_,length(vec));v[vec==a]<-0;v[vec==b]<-1;v}
grp <- as.character(sample_info$ExpGroup)
contrasts <- list(con_res = mk_contrast(grp,"con","res"),
                  con_sus = mk_contrast(grp,"con","sus"),
                  res_sus = mk_contrast(grp,"res","sus"))
for (nm in names(contrasts)) {
  v <- contrasts[[nm]]; keep <- !is.na(v)
  if (sum(keep) >= 3) {
    cmat <- cor(mergedMEs[keep, , drop=FALSE], v[keep], use="p")
    pmat <- corPvalueStudent(cmat, sum(keep))
    txt <- paste(signif(cmat,2), "\n(", signif(pmat,1), ")", sep = "")
    dim(txt) <- dim(cmat)
    svglite::svglite(file = fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")), width = 3, height = max(4, ncol(mergedMEs) * 0.25 + 2))
    par(mar = c(6, 8.5, 3, 1))
    labeledHeatmap(Matrix = cmat, xLabels = nm, yLabels = colnames(mergedMEs), ySymbols = colnames(mergedMEs),
                   colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt, setStdMargins = FALSE,
                   cex.text = 0.8, zlim = c(-1,1), main = "Module–trait relationships")
    dev.off()
  }
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
    parallelCalculation = TRUE,
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
region <- as.character(sample_info$region)
layer  <- as.character(sample_info$layer)
cell   <- as.character(sample_info$celltype)
cond   <- as.character(sample_info$ExpGroup)

combo <- paste(cond, region, layer, cell, sep = "_")
combo <- factor(combo)

X_combo <- model.matrix(~ 0 + combo)
colnames(X_combo) <- levels(combo)
if (ncol(X_combo) == 0) stop("No combined strata present in X_combo (ncol == 0).")

MEcorr_combo <- cor(mergedMEs, X_combo, use = "p")
MEp_combo    <- corPvalueStudent(MEcorr_combo, nrow(expression.data))

df_combo <- reshape2::melt(MEcorr_combo, varnames = c("module","comb"), value.name = "r")
p_combo  <- reshape2::melt(MEp_combo,    varnames = c("module","comb"), value.name = "p")
df_combo$p <- p_combo$p
df_combo$stars <- gtools::stars.pval(df_combo$p)
df_combo$module <- factor(df_combo$module, levels = rownames(MEcorr_combo))

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
# pheatmap block for ME–trait matrix (robust to missing celltype)
# --------------------------
# Reconstruct df_all / module_levels / traits_in_order / strip_df if they are not present
if (!exists("df_all")) {
  if (exists("df_combo2")) {
    df_all <- df_combo2
    if (!"trait" %in% names(df_all)) {
      df_all$trait <- paste(df_all$condition, df_all$region, df_all$layer, df_all$cell, sep = "_")
    }
  } else if (exists("MEcorr_combo") && exists("MEp_combo")) {
    df_all <- reshape2::melt(MEcorr_combo, varnames = c("module", "comb"), value.var = "r")
    p_tmp  <- reshape2::melt(MEp_combo,    varnames = c("module", "comb"), value.var = "p")
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

if (!exists("module_levels")) {
  if (exists("MEcorr_combo")) {
    module_levels <- rownames(MEcorr_combo)
  } else {
    module_levels <- unique(df_all$module)
  }
}

if (!exists("traits_in_order") || length(traits_in_order) == 0) {
  df_all$trait <- paste(df_all$condition, df_all$region, df_all$layer, df_all$cell, sep = "_")
  ord <- order(df_all$condition, df_all$region, df_all$layer, df_all$cell)
  traits_in_order <- unique(df_all$trait[ord])
}

library(grid)

r_mat  <- reshape2::acast(df_all, module ~ trait, value.var = "r")
p_mat  <- reshape2::acast(df_all, module ~ trait, value.var = "p")
fdr_vec <- p.adjust(as.vector(p_mat), method = "BH")
fdr_mat <- matrix(fdr_vec, nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))

stopifnot(all(module_levels %in% rownames(r_mat)), all(traits_in_order %in% colnames(r_mat)))
r_mat   <- r_mat[module_levels, traits_in_order, drop = FALSE]
p_mat   <- p_mat[module_levels, traits_in_order, drop = FALSE]
fdr_mat <- fdr_mat[module_levels, traits_in_order, drop = FALSE]

hc_cols <- hclust(as.dist(1 - cor(r_mat, use = "pairwise.complete.obs")), method = "average")
k_clusters <- 6
col_grp <- cutree(hc_cols, k = k_clusters)
grp_ord <- col_grp[hc_cols$order]
gap_pos <- which(grp_ord[-1] != head(grp_ord, -1))
xlines <- gap_pos + 0.5

# Module color annotation
if (exists("mergedMEs")) {
  me_names <- colnames(mergedMEs)               # e.g., "MEblue"
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

# Palette and breaks
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", repos = "https://cloud.r-project.org")
r_lim <- 0.8
bk <- seq(-r_lim, r_lim, length.out = 201)
pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(200)

# Significance dots
display_mat <- matrix("", nrow = nrow(r_mat), ncol = ncol(r_mat), dimnames = dimnames(r_mat))
display_mat[fdr_mat < 0.01] <- "•"

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
grid::seekViewport(ph$gtable$layout$name[panel_id])
for (xl in xlines) {
  grid::grid.lines(x = grid::unit(c(xl, xl), "native"),
                   y = grid::unit(c(0, nr), "native"),
                   gp = grid::gpar(col = "white", lwd = 3))
}
grid::upViewport(0); dev.off()

pdf(fp_traits("ME_trait_pheatmap.pdf"), width = 12, height = 7)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
panel_id <- grep("matrix", ph$gtable$layout$name)[1]
grid::seekViewport(ph$gtable$layout$name[panel_id])
for (xl in xlines) {
  grid::grid.lines(x = grid::unit(c(xl, xl), "native"),
                   y = grid::unit(c(0, nr), "native"),
                   gp = grid::gpar(col = "white", lwd = 2))
}
grid::upViewport(0); dev.off()

svglite::svglite(fp_traits("ME_trait_pheatmap.svg"), width = 12, height = 7)
grid::grid.newpage(); grid::grid.draw(ph$gtable)
panel_id <- grep("matrix", ph$gtable$layout$name)[1]
grid::seekViewport(ph$gtable$layout$name[panel_id])
for (xl in xlines) {
  grid::grid.lines(x = grid::unit(c(xl, xl), "native"),
                   y = grid::unit(c(0, nr), "native"),
                   gp = grid::gpar(col = "white", lwd = 2))
}
grid::upViewport(0); dev.off()

# --------------------------
# End of script
# --------------------------



























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
                 reshape2, gtools, patchwork, cowplot, pheatmap, install = TRUE)
)

# Parallel setup: prefer physical cores; enable WGCNA threading globally
nCores <- tryCatch({
  pc <- parallel::detectCores(logical = FALSE)
  if (is.na(pc) || pc < 2) 2 else pc
}, error = function(e) 2)
enableWGCNAThreads(nThreads = nCores)  # replaces allowWGCNAThreads(); controls threading globally
WGCNAnThreads()                        # optional: prints current thread count

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

male.data <- readxl::read_excel(path = expr_xlsx)
male.data <- male.data[grepl("_MOUSE", male.data$gene_symbol), ]

# Expression matrix: rows = samples, cols = proteins
expression.data <- male.data[, -1]
expression.data <- as.data.frame(lapply(expression.data, as.numeric))
expression.data <- as.data.frame(t(expression.data))
names(expression.data) <- male.data$gene_symbol

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
par(cex = 0.6, mar = c(0,4,2,0))
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
plotEigengeneNetworks(MET, "", plotHeatmaps = FALSE, marDendro = c(0,4,2,0))
dev.off()
svg(file = fp_net("eigengene_adjacency_heatmap.svg"), width = 8, height = 6)
par(mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", plotDendrograms = FALSE,
                      marHeatmap = c(5,5,2,2), xLabelsAngle = 90)
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
  textMatrix <- paste(signif(matCorr[, cols, drop=FALSE], 2), "\n(", signif(matP[, cols, drop=FALSE], 1), ")", sep = "")
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
contrasts <- list(con_res = mk_contrast(grp,"con","res"),
                  con_sus = mk_contrast(grp,"con","sus"),
                  res_sus = mk_contrast(grp,"res","sus"))
for (nm in names(contrasts)) {
  v <- contrasts[[nm]]; keep <- !is.na(v)
  cmat <- cor(mergedMEs[keep, , drop=FALSE], v[keep], use="p")
  pmat <- corPvalueStudent(cmat, sum(keep))
  txt <- paste(signif(cmat,2), "\n(", signif(pmat,1), ")", sep = "")
  dim(txt) <- dim(cmat)
  svg(file = fp_traits(paste0("ME_trait_heatmap_", nm, ".svg")), width = 3, height = max(4, ncol(mergedMEs) * 0.25 + 2))
  par(mar = c(6, 8.5, 3, 1))
  labeledHeatmap(Matrix = cmat, xLabels = nm, yLabels = colnames(mergedMEs), ySymbols = colnames(mergedMEs),
                 colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = txt, setStdMargins = FALSE,
                 cex.text = 0.8, zlim = c(-1,1), main = "Module–trait relationships")
  dev.off()
}

# --------------------------
# kME, GS, hubs
# --------------------------
modNames <- substring(colnames(mergedMEs), 3)
geneModuleMembership <- as.data.frame(cor(expression.data, mergedMEs, use = "p", method = "pearson"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

cond_numeric <- setNames(c(1,2,3), c("con","res","sus"))
ExpGroup_num <- as.numeric(cond_numeric[as.character(sample_info$ExpGroup)])
geneTraitSignificance <- as.data.frame(cor(expression.data, ExpGroup_num, use = "p", method = "pearson"))
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
