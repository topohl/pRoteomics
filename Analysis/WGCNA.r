# Use pacman to install/load packages (and BiocManager for Bioconductor packages)
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Ensure GO.db (Bioconductor) is installed
if (!requireNamespace("GO.db", quietly = TRUE)) {
  BiocManager::install("GO.db", ask = FALSE, update = FALSE)
}

# Install (if missing) and load required packages with pacman, suppressing startup messages
suppressPackageStartupMessages(
  pacman::p_load(
    WGCNA,
    flashClust,
    curl,
    readxl,
    ggplot2,
    svglite,
    GO.db,
    install = TRUE
  )
)

# --------------------------------------------------------------------
# Load and prepare data
# --------------------------------------------------------------------

# define output_dir
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# load data
d <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/male.data.xlsx"
male.data <- read_excel(path = d)

# remove rows from male.data where the first column (gene_symbol) doesn't contain "_MOUSE"
male.data <- male.data[grepl("_MOUSE", male.data$gene_symbol), ]

# expression.data <- liver.data[,-c(1:8)] # removing variables not holding expression data
expression.data <- male.data[, -1]

# Make sure expression data is numeric before transposing
expression.data <- as.data.frame(lapply(expression.data, as.numeric))

# transform the data.frame so columns represent genes and rows represent samples
expression.data <- as.data.frame(t(expression.data))
names(expression.data) <- male.data$gene_symbol

# --------------------------------------------------------------------
# Quality check and clustering of samples to detect outliers
# --------------------------------------------------------------------

# Check for genes and samples with too many missing values
gsg <- goodSamplesGenes(expression.data)
summary(gsg)
gsg$allOK

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(
      paste(
        "Removing genes:",
        paste(names(expression.data)[!gsg$goodGenes], collapse = ", ")
      )
    ) # Identifies and prints outlier genes
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(
      paste(
        "Removing samples:",
        paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", ")
      )
    ) # Identifies and prints outlier samples
  }
  expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes offending genes and samples
}

sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 

# Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

# Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)

# Setting the graphical parameters and plot to screen
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2)
abline(h = 40, col = "red")

# Also save the same plot as SVG in output_dir
svg(file = file.path(output_dir, "sample_clustering_outliers.svg"), width = 8, height = 6)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
  cex.axis = 1.5, cex.main = 2)
abline(h = 40, col = "red")
dev.off()

# Determine cluster under the line
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)

# Remove outlier samples
expression.data <- expression.data[cut.sampleTree==1, ]
spt <- pickSoftThreshold(expression.data)

# Plot the results:
svg(file = file.path(output_dir, "soft_threshold_scale_independence.svg"), width = 7, height = 5)
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,2],
  xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
  main = paste("Scale independence"))
text(spt$fitIndices[,1], spt$fitIndices[,2], col="red")
abline(h=0.80, col="red")
dev.off()

svg(file = file.path(output_dir, "soft_threshold_mean_connectivity.svg"), width = 7, height = 5)
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
  xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
  main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1], col="red")
dev.off()

# ---------------------------------------------------------------------
# Constructing a gene co-expression network and identifying modules
# ---------------------------------------------------------------------

# Based on the above plots, we choose a soft-thresholding power of 6
softPower <- 6 
adjacency <- adjacency(expression.data, power = softPower)
TOM <- TOMsimilarity(adjacency)
TOM.dissimilarity <- 1-TOM

# Clustering using TOM-based dissimilarity, creating dendrogram
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

# Also save the same plot as SVG in output_dir
svg_filename <- file.path(output_dir, "gene_dendrogram.svg")
svg(file = svg_filename, width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

# Set the minimum module size, e.g., 30
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
table(Modules)

#ModuleColors <- labels2colors(Modules) #assigns each module number a color
# create ModuleColors list 

table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

svg(file = file.path(output_dir, "gene_dendrogram_module_colors.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, ModuleColors, "Module",
          dendroLabels = FALSE, hang = 0.03,
          addGuide = TRUE, guideHang = 0.05,
          main = "Gene dendrogram and module colors")
dev.off()

# Calculate Eigengenes
MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

# Calculate dissimilarity of module eigengenes
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

# Cluster module eigengenes
METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes
par(mar = c(0,4,2,0)) # setting margin sizes
par(cex = 0.6); # scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75

# Merge modules whose eigengenes are very similar, i.e. their correlation is above 0.75 (height below 0.25)
merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)

# The merged module colors, assigning one color to each module
mergedColors = merge$colors

# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

svg(file = file.path(output_dir, "gene_dendrogram_modules_merged.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors),
          c("Original Module", "Merged Module"),
          dendroLabels = FALSE, hang = 0.03,
          addGuide = TRUE, guideHang = 0.05,
          main = "Gene dendrogram and module colors for original and merged modules")
dev.off()

# ---------------------------------------------------------------------
# Loading trait data (metadata) and associating it with modules
# ---------------------------------------------------------------------

# Load trait data
t <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info.xlsx"
traitData <- read_excel(path = t)
head(traitData)

allTraits <- traitData[, -c(2)] # Removing notes and comments sections
allTraits <- allTraits[, c(1, 5:15)] # Pulling out only continuous traits
Samples <- rownames(expression.data)
traitRows <- match(Samples, allTraits$row.names)
datTraits <- allTraits[traitRows, -1]

# convert to a data.frame (not a tibble) so row names can be set safely
datTraits <- as.data.frame(datTraits, stringsAsFactors = FALSE)
# assign row names as a character vector from the first column of allTraits
rownames(datTraits) <- as.character(allTraits[[1]][traitRows])

# keep only ExpGroup in datTraits for simplicity
datTraits <- datTraits["ExpGroup"]

# ---------------------------------------------------------------------
# Creating module-trait relationships (in this case ExpGroup)
# ---------------------------------------------------------------------

# Helper to create a comparison-ready datTraits by mapping factor levels to numeric values
create_datTraits_comparison <- function(datTraits, column = "ExpGroup", values_map = c(con = 1, res = 2, sus = 3)) {
  if (!column %in% names(datTraits)) stop("Column ", column, " not found in datTraits")
  dt <- datTraits
  colvec <- as.character(dt[[column]])
  newcol <- rep(NA_real_, length(colvec))
  for (nm in names(values_map)) {
    newcol[colvec == nm] <- as.numeric(values_map[[nm]])
  }
  dt[[column]] <- newcol
  dt
}

# Create the three comparisons used downstream
# values_map defines how to convert ExpGroup factor levels to numeric values
# e.g. con = 1, res = 2, sus = 3 for a continuous comparison
# this is done so that the correlation analysis later makes sense
datTraits_consus <- create_datTraits_comparison(datTraits, "ExpGroup", c(con = 1, sus = 3))
datTraits_conres <- create_datTraits_comparison(datTraits, "ExpGroup", c(con = 1, res = 2))
datTraits_susres <- create_datTraits_comparison(datTraits, "ExpGroup", c(sus = 3, res = 2))

# Choose which comparison to use for module-trait correlation analysis and set a stable label
# Update detected_label (used later) and keep label for backward compatibility with existing checks
datTrait <- datTraits_consus
detected_label <- "consus"
label <- detected_label

# Define numbers of genes and samples
nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
module.trait.correlation = cor(mergedMEs, datTrait, use = "p") #p for pearson correlation coefficient
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
           signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
par(mar = c(6, 8.5, 3, 1))

# Display the correlation values within a heatmap plot (to screen)
labeledHeatmap(Matrix = module.trait.correlation,
         xLabels = names(datTraits),
         yLabels = names(mergedMEs),
         ySymbols = names(mergedMEs),
         colorLabels = FALSE,
         colors = blueWhiteRed(50),
         textMatrix = textMatrix,
         setStdMargins = FALSE,
         cex.text = 1,
         zlim = c(-1,1),
         main = paste("Module-trait relationships"))

# Also save the same heatmap as SVG in output_dir
outdir <- file.path(output_dir, "module_trait_relationships")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
label <- if (exists("datTrait") && exists("datTraits_consus") && identical(datTrait, datTraits_consus)) {
  "consus"
} else if (exists("datTrait") && exists("datTraits_conres") && identical(datTrait, datTraits_conres)) {
  "conres"
} else if (exists("datTrait") && exists("datTraits_susres") && identical(datTrait, datTraits_susres)) {
  "susres"
} else {
  # fallback: try to infer from variable name if datTrait was assigned by name
  inferred <- tryCatch({
    nm <- utils::which(sapply(ls(envir = .GlobalEnv), function(x) identical(get(x, envir = .GlobalEnv), datTrait)))
    if (length(nm) && nzchar(ls(envir = .GlobalEnv)[nm[1]])) ls(envir = .GlobalEnv)[nm[1]] else "datTrait"
  }, error = function(e) "datTrait")
  gsub("[^A-Za-z0-9_\\-]", "_", inferred)
}
svg(file = file.path(outdir, paste0("module_trait_relationships_", label, ".svg")), width = 3, height = 8)

par(mar = c(6, 8.5, 3, 1))
labeledHeatmap(Matrix = module.trait.correlation,
         xLabels = names(datTraits),
         yLabels = names(mergedMEs),
         ySymbols = names(mergedMEs),
         colorLabels = FALSE,
         colors = blueWhiteRed(50),
         textMatrix = textMatrix,
         setStdMargins = FALSE,
         cex.text = 0.7,
         zlim = c(-1,1),
         main = paste("Module-trait relationships"))
dev.off()

## experimental
# # Make sure row names match your expression data sample names
# rownames(datTraits_consus) <- $row.names
# 
# # ===== CREATE TRAIT MATRIX FROM YOUR METADATA =====
# # Convert categorical variables to numeric for WGCNA correlation analysis
# 
# # Cell type traits (binary)
# celltype_traits <- model.matrix(~ celltype - 1, data = datTrait)
# colnames(celltype_traits) <- gsub("celltype", "", colnames(celltype_traits))
# 
# # Layer traits (binary)
# layer_traits <- model.matrix(~ layer - 1, data = datTraits)
# colnames(layer_traits) <- paste0("layer_", colnames(layer_traits))
# 
# # Region traits (binary) 
# region_traits <- model.matrix(~ region - 1, data = datTraits)
# colnames(region_traits) <- paste0("region_", colnames(region_traits))
# 
# # ReplicateGroup traits (binary)
# replicate_traits <- model.matrix(~ ReplicateGroup - 1, data = datTraits)
# colnames(replicate_traits) <- paste0("hemisphere_", colnames(replicate_traits))
# 
# # ExpGroup as continuous
# datTraits$ExpGroup_numeric <- as.numeric(as.factor(datTraits$ExpGroup))
# 
# # Create final trait matrix
# datTraits <- data.frame(
#   # Cell type information
#   celltype_traits,
#   
#   # Layer information
#   layer_traits,
#   
#   # Region information  
#   region_traits,
#   
#   # Hemisphere information
#   replicate_traits,
#   
#   # Experimental group (continuous)
#   ExpGroup = sample_info$ExpGroup_numeric,
#   
#   # Biologically meaningful combinations
#   CA1_neuron = as.numeric(sample_info$region == "CA1" & sample_info$celltype == "neuron"),
#   CA2_neuron = as.numeric(sample_info$region == "CA2" & sample_info$celltype == "neuron"), 
#   CA3_neuron = as.numeric(sample_info$region == "CA3" & sample_info$celltype == "neuron"),
#   CA1_microglia = as.numeric(sample_info$region == "CA1" & sample_info$celltype == "microglia"),
#   CA2_microglia = as.numeric(sample_info$region == "CA2" & sample_info$celltype == "microglia"),
#   CA3_microglia = as.numeric(sample_info$region == "CA3" & sample_info$celltype == "microglia"),
#   
#   # Cell layer specificity
#   neuron_soma = as.numeric(sample_info$celltype_layer == "neuron_soma"),
#   neuron_neuropil = as.numeric(sample_info$celltype_layer == "neuron_neuropil"),
#   
#   stringsAsFactors = FALSE
# )
# 
# # remove the columns in expression data where the column name has anything other than "_MOUSE" in it as part of the name, but skip the first column
# expression.data <- expression.data[, grepl("_MOUSE", colnames(expression.data))]
# 
# # Match with expression data samples
# Samples <- rownames(expression.data)
# traitRows <- match(Samples, rownames(sample_info))
# datTraits <- datTraits[traitRows, ]
# rownames(datTraits) <- Samples
# 
# # Remove any columns with all NAs or constant values
# datTraits <- datTraits[, sapply(datTraits, function(x) length(unique(x[!is.na(x)]))) > 1]
# 
# # ===== CONTINUE WITH MODULE-TRAIT CORRELATION ANALYSIS =====
# # Calculate module-trait correlations
# nGenes = ncol(expression.data)
# nSamples = nrow(expression.data)
# module.trait.correlation = cor(mergedMEs, datTraits, use = "p")
# module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples)
# 
# # Create heatmap
# textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
#                    signif(module.trait.Pvalue, 1), ")", sep = "");
# dim(textMatrix) = dim(module.trait.correlation)
# par(mar = c(6, 8.5, 3, 1))
# 
# labeledHeatmap(Matrix = module.trait.correlation,
#                xLabels = names(datTraits),
#                yLabels = names(mergedMEs),
#                ySymbols = names(mergedMEs),
#                colorLabels = FALSE,
#                colors = blueWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.4,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships"))
# 
# # ## experimental

# Define variable ExpGroup from the dataset
ExpGroup = as.data.frame(datTraits_consus$ExpGroup)
names(ExpGroup) = "ExpGroup"

modNames = substring(names(mergedMEs), 3) #extract module names

# Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expression.data, ExpGroup, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(ExpGroup), sep="")
names(GSPvalue) = paste("p.GS.", names(ExpGroup), sep="")
head(GSPvalue)

# Plot module membership vs gene significance for a given module and trait
module <- "grey60"
col_idx <- match(module, modNames)
moduleGenes <- mergedColors == module

if (!any(moduleGenes)) {
  message(sprintf("No genes in module '%s' — skipping plot.", module))
} else {
  mm_vals <- abs(as.numeric(geneModuleMembership[moduleGenes, col_idx]))
  gs_vals <- abs(as.numeric(geneTraitSignificance[moduleGenes, 1]))
  genes <- colnames(expression.data)[moduleGenes]
  df_ms <- data.frame(Gene = genes, MM = mm_vals, GS = gs_vals, stringsAsFactors = FALSE)

  # correlation (pearson). Protect against too few points.
  if (nrow(df_ms) >= 3) {
    ct <- cor.test(df_ms$MM, df_ms$GS, method = "pearson")
    r_val <- as.numeric(ct$estimate)
    p_val <- ct$p.value
    r_txt <- sprintf("r = %.2f", r_val)
    p_txt <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3g", p_val))
  } else {
    r_val <- cor(df_ms$MM, df_ms$GS, use = "complete.obs")
    p_val <- NA
    r_txt <- sprintf("r = %.2f", r_val)
    p_txt <- "p = NA (n < 3)"
  }

  caption_text <- paste(r_txt, p_txt, sep = "   ")

  # Create scatter plot with regression line
  p_mm_gs <- ggplot(df_ms, aes(x = MM, y = GS)) +
    geom_point(shape = 21, size = 6, fill = module, color = "white", stroke = 0.6, alpha = 0.95) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.8) +
    labs(
      x = paste0("Module Membership (kME) — ", module),
      y = "Gene significance (|GS|)",
      title = paste("Module membership vs gene significance —", module),
      caption = caption_text
    ) +
    coord_cartesian(expand = TRUE) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 20),
      axis.title = element_text(size = 14),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(face = "plain", size = 20),
      axis.line = element_line(size = 0.8, colour = "black"),
      axis.ticks = element_line(size = 0.6, colour = "black"),
      plot.caption = element_text(size = 10, hjust = 0, margin = margin(t = 10))
    )

  print(p_mm_gs)

  # save SVG to output_dir
  if (!exists("label")) label <- "datTrait"
  outfile <- file.path(output_dir, paste0("module_membership_vs_gene_significance_", module, "_", label, ".svg"))
  ggsave(filename = outfile, plot = p_mm_gs, device = svglite::svglite, width = 4, height = 5, units = "in", dpi = 300)
}

# Isolate weight from the clinical traits
ExpGroup = as.data.frame(datTraits_consus$ExpGroup);
names(ExpGroup) = "ExpGroup"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, ExpGroup))
# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)

# Plot the dendrogram
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
            plotHeatmaps = FALSE)

# ensure a label exists (fall back to "datTrait" if not)
if (!exists("label")) {
  label <- if (exists("datTrait") && exists("datTraits_consus") && identical(datTrait, datTraits_consus)) {
  "consus"
  } else if (exists("datTrait") && exists("datTraits_conres") && identical(datTrait, datTraits_conres)) {
  "conres"
  } else if (exists("datTrait") && exists("datTraits_susres") && identical(datTrait, datTraits_susres)) {
  "susres"
  } else {
  "datTrait"
  }
}

# Save dendrogram as SVG in output_dir using the label
svg(file = file.path(output_dir, paste0("eigengene_dendrogram_", label, ".svg")), width = 8, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
            plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)

# Ensure a label exists (fall back to "datTrait" if not)
if (!exists("label")) {
  label <- if (exists("datTrait") && exists("datTraits_consus") && identical(datTrait, datTraits_consus)) {
  "consus"
  } else if (exists("datTrait") && exists("datTraits_conres") && identical(datTrait, datTraits_conres)) {
  "conres"
  } else if (exists("datTrait") && exists("datTraits_susres") && identical(datTrait, datTraits_susres)) {
  "susres"
  } else {
  "datTrait"
  }
}

# Save heatmap as SVG in output_dir using the label
svg(file = file.path(output_dir, paste0("eigengene_adjacency_heatmap_", label, ".svg")), width = 8, height = 6)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

# Isolate ExpGroup from the dataset
ExpGroup = as.data.frame(datTraits$ExpGroup);
names(ExpGroup) = "ExpGroup"
layer = as.data.frame(datTraits$layer)
names(layer) = "layer"

# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, ExpGroup, layer))

# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)

# Check for modules of interest
table(mergedColors)

# Extract genes that belong to specific modules and save to csv file in output directory, do that for all available modules in mergedColors
modules_of_interest <- unique(mergedColors)
for (module in modules_of_interest) {
  moduleGenes <- mergedColors == module
  # Skip modules that have no genes to avoid creating data.frames with mismatched row lengths
  if (!any(moduleGenes)) {
    message(sprintf("No genes found in module '%s'; skipping.", module))
    next
  }
  gene_info <- data.frame(
    Gene = colnames(expression.data)[moduleGenes],
    Module = module,
    ModuleMembership = geneModuleMembership[moduleGenes, paste0("MM", module)],
    GeneSignificance = geneTraitSignificance[moduleGenes, "GS.ExpGroup"]
  )
  # ensure a label exists (fall back to "datTrait" if not)
if (!exists("label")) {
  label <- if (exists("datTrait") && exists("datTraits_consus") && identical(datTrait, datTraits_consus)) {
  "consus"
  } else if (exists("datTrait") && exists("datTraits_conres") && identical(datTrait, datTraits_conres)) {
  "conres"
  } else if (exists("datTrait") && exists("datTraits_susres") && identical(datTrait, datTraits_susres)) {
  "susres"
  } else {
  "datTrait"
  }
}
  output_file <- file.path(output_dir, paste0("genes_in_module_", module, "_", label, ".csv"))
  write.csv(gene_info, file = output_file, row.names = FALSE)
}

# Identify top hub genes per module based on kME (Module Membership)
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
if (!exists("label")) {
  label <- if (exists("datTrait") && exists("datTraits_consus") && identical(datTrait, datTraits_consus)) {
  "consus"
  } else if (exists("datTrait") && exists("datTraits_conres") && identical(datTrait, datTraits_conres)) {
  "conres"
  } else if (exists("datTrait") && exists("datTraits_susres") && identical(datTrait, datTraits_susres)) {
  "susres"
  } else {
  "datTrait"
  }
}
# write one summary file
write.csv(do.call(rbind, lapply(names(top_hubs), function(m) cbind(Module=m, top_hubs[[m]]))),
          file.path(output_dir, paste0("top_hubs_per_module_", label, ".csv")), row.names = FALSE)

# explicit module choice or set to NULL to auto-detect
chosen_module <- "grey60"  # e.g. "grey60" or NULL to infer

if (!is.null(chosen_module)) {
  me_col <- paste0("ME", chosen_module)
  if (!me_col %in% colnames(mergedMEs)) stop("Column ", me_col, " not found in mergedMEs")
  moduleEigs <- mergedMEs[[me_col]]
  module_name <- chosen_module
} 

# remove any leftover loop variable to avoid accidental reuse later
if (exists("module")) rm(module)

# prepare plot directory & filename
plot_dir <- file.path(output_dir, "module_eigengene_plots", detected_label)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(plot_dir, paste0("module_eigengene_", module_name, "_", detected_label, ".svg"))

# create group_Factor_valid from datTrait$ExpGroup, ensuring it is a factor with levels in the desired order
if (!"ExpGroup" %in% names(datTrait)) stop("Column 'ExpGroup' not found in datTrait")

# raw group vector (as character)
gf_raw <- as.character(datTrait$ExpGroup)

# attempt to align to expression sample order if possible; if alignment fails, keep original order
if (exists("Samples") && !is.null(rownames(datTrait))) {
  idx_samples <- match(Samples, rownames(datTrait))
  if (all(is.na(idx_samples))) {
    gf_aligned <- gf_raw
  } else {
    gf_aligned <- gf_raw[idx_samples]
  }
} else {
  gf_aligned <- gf_raw
}

# map short labels to readable group names
label_map <- c(con = "Control", res = "Resilient", sus = "Susceptible")
mapped <- unname(label_map[gf_aligned])

# If mapping failed (e.g. ExpGroup contains numeric codes), try numeric-to-label mapping
if (all(is.na(mapped))) {
  numeric_map <- c("1" = "Control", "2" = "Resilient", "3" = "Susceptible")
  mapped <- unname(numeric_map[as.character(gf_aligned)])
}

# fallback to original values if mapping still yields NA (preserve whatever values are present)
mapped[is.na(mapped)] <- as.character(gf_aligned[is.na(mapped)])

# drop positions with NA in either group or module eigengene
moduleEigs_current <- if (exists("moduleEigs")) moduleEigs else stop("moduleEigs not found")
keep_idx <- !is.na(mapped) & !is.na(moduleEigs_current)
# derive mapped vector filtered by keep_idx
mapped_valid <- mapped[keep_idx]

# desired canonical order
desired_order <- c("Control", "Resilient", "Susceptible")

# keep only levels that are actually present, preserving desired order
present_levels <- desired_order[desired_order %in% unique(mapped_valid)]

# create factor restricted to present levels (no extraneous levels)
group_factor_valid <- factor(mapped_valid, levels = present_levels)

# fallback if nothing matched (avoid zero-level factor)
if (length(levels(group_factor_valid)) == 0) group_factor_valid <- factor(mapped_valid)
moduleEigs_valid <- as.numeric(moduleEigs_current[keep_idx])

# final_levels used for plotting (only levels actually present)
final_levels <- levels(droplevels(group_factor_valid))

# assemble plotting dataframe
df_plot <- data.frame(
  ME = moduleEigs_valid,
  Group = factor(as.character(group_factor_valid), levels = final_levels)
)

# colors: ensure mapping contains only the groups present
color_map <- c("Control" = "#457B9D", "Resilient" = "#C6C3BB", "Susceptible" = "#E63946")
plot_colors <- color_map[intersect(names(color_map), levels(df_plot$Group))]
# fallback if none match
if (length(plot_colors) == 0) plot_colors <- rep("#999999", length(levels(df_plot$Group)))

# jitter plot + median and 25-75% intervals + simple stats caption
# compute per-group summary (n, median, Q1, Q3)
stats_df <- do.call(rbind, lapply(split(df_plot$ME, df_plot$Group), function(x) {
  data.frame(
    Group = unique(df_plot$Group[ df_plot$ME %in% x ]),
    n = sum(!is.na(x)),
    median = median(x, na.rm = TRUE),
    Q1 = quantile(x, 0.25, na.rm = TRUE),
    Q3 = quantile(x, 0.75, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
stats_df$Group <- factor(stats_df$Group, levels = levels(df_plot$Group))

# choose test depending on number of groups
gcount <- nlevels(df_plot$Group)
if (gcount == 2) {
  test_res <- wilcox.test(ME ~ Group, data = df_plot)
  test_name <- "Wilcoxon rank-sum"
  pval <- test_res$p.value
} else {
  test_res <- kruskal.test(ME ~ Group, data = df_plot)
  test_name <- "Kruskal-Wallis"
  pval <- test_res$p.value
}
pval_text <- ifelse(pval < 0.001, "<0.001", formatC(pval, format = "f", digits = 3))

# build caption with per-group stats and test result
group_lines <- apply(stats_df, 1, function(r) {
  sprintf("%s: n=%d, median=%.3f, IQR=[%.3f, %.3f]", r["Group"], as.integer(r["n"]),
          as.numeric(r["median"]), as.numeric(r["Q1"]), as.numeric(r["Q3"]))
})
caption_text <- paste(paste(group_lines, collapse = " | "),
                      sprintf("\nTest: %s, p = %s", test_name, pval_text))

# final plot: jittered points, median points, and 25-75% interval errorbars
p <- ggplot(df_plot, aes(x = Group, y = ME, fill = Group)) +
  geom_jitter(width = 0.16, height = 0, size = 7, shape = 21, color = "white", alpha = 0.85, stroke = 0.6) +
  # 25-75% interval as thick errorbar; set inherit.aes = FALSE so the layer does not try to use y = ME from the global mapping
  geom_linerange(data = stats_df, aes(x = Group, ymin = Q1, ymax = Q3),
                 size = 1.1, color = "black", inherit.aes = FALSE) +
  # median as filled point (override global mapping with its own y), being explicit about aesthetics
  geom_point(data = stats_df, aes(x = Group, y = median), shape = 21, size = 3.5, fill = "black", color = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = plot_colors) +
  labs(x = NULL,
       y = paste0("Module Eigengene (", module_name, ")"),
       title = paste0("Module Eigengene by Group — ", detected_label),
       caption = caption_text) +
  coord_cartesian(ylim = c(-0.3, 0.2)) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 20),
    axis.title.y = element_text(margin = margin(r = 20)),
    plot.title = element_text(face = "plain", size = 20),
    axis.line = element_line(size = 0.8, colour = "black"),
    axis.ticks = element_line(size = 0.6, colour = "black"),
    plot.caption = element_text(size = 10, hjust = 0, margin = margin(t = 10))
  )

# save SVG
ggsave(filename = outfile, plot = p, device = svglite::svglite, width = 4, height = 5, units = "in", dpi = 300)