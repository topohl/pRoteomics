# Option 1: Run R as administrator (recommended for system-wide installs)
# Option 2: Install to a user-writable directory:
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
install.packages('BiocManager', lib = user_lib)
library(BiocManager, lib.loc = user_lib)
BiocManager::install('GO.db', lib = user_lib) # Explicitly install GO.db first
BiocManager::install('WGCNA', lib = user_lib)
BiocManager::install('flashClust', lib = user_lib)


library(WGCNA)
library(flashClust)
library(curl)
library(readxl)
library(ggplot2)
library(svglite)


#d <- curl('https://raw.githubusercontent.com/fuzzyatelin/fuzzyatelin.github.io/master/bioanth-stats/module-F21-Group1/FemaleLiver-Data/LiverFemale3600.csv')
#liver.data <- read.csv(file = d, stringsAsFactors = FALSE, header = TRUE)
#head(liver.data)

# define output_dir
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/wgcna/output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# load data
d <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/male.data.xlsx"

male.data <- read_excel(path = d)

# remove rows frome male.data where the first column (gene_symbol) doesnt contain "_MOUSE" anywhere in the name
male.data <- male.data[grepl("_MOUSE", male.data$gene_symbol), ]

#expression.data <- liver.data[,-c(1:8)] #removing variables not holding expression data
expression.data <- male.data[, -1]

# Make sure expression data is numeric before transposing
expression.data <- as.data.frame(lapply(expression.data, as.numeric))

expression.data <- as.data.frame(t(expression.data)) #transforming the data.frame so columns now represent genes and rows represent samples
names(expression.data) <- male.data$gene_symbol
#renaming the columns so we don't lose the gene IDs

gsg <- goodSamplesGenes(expression.data)

summary(gsg)

gsg$allOK

if (!gsg$allOK)
{
if (sum(!gsg$goodGenes)>0) 
printFlush(paste("Removing genes:", paste(names(expression.data)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(expression.data)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
expression.data <- expression.data[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
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

cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10) #returns numeric vector
#Remove outlier
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

# Based on the above plots, we choose a soft-thresholding power of 6
softPower <- 6 
adjacency <- adjacency(expression.data, power = softPower)

TOM <- TOMsimilarity(adjacency)

TOM.dissimilarity <- 1-TOM

#creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

# plotting the dendrogram (display and save to SVG)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

# save as SVG in output_dir
svg_filename <- file.path(output_dir, "gene_dendrogram.svg")
svg(file = svg_filename, width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)

table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module.

ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

svg(file = file.path(output_dir, "gene_dendrogram_module_colors.svg"), width = 12, height = 9)
plotDendroAndColors(geneTree, ModuleColors, "Module",
          dendroLabels = FALSE, hang = 0.03,
          addGuide = TRUE, guideHang = 0.05,
          main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

# Calculate dissimilarity of module eigengenes
ME.dissimilarity = 1-cor(MElist$eigengenes, use="complete") #Calculate eigengene dissimilarity

METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75

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

#t <- curl('https://raw.githubusercontent.com/fuzzyatelin/fuzzyatelin.github.io/master/bioanth-stats/module-F21-Group1/FemaleLiver-Data/ClinicalTraits.csv')
#traitData <- read.csv(t, header = TRUE, stringsAsFactors = FALSE)
#head(traitData)

t <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info.xlsx"
traitData <- read_excel(path = t)
head(traitData)

allTraits <- traitData[, -c(2)] #removing notes and comments sections 
allTraits <- allTraits[, c(1, 5:15) ] #pulling out only continuous traits

Samples <- rownames(expression.data)
traitRows <- match(Samples, allTraits$row.names)
datTraits <- allTraits[traitRows, -1]
# convert to a data.frame (not a tibble) so row names can be set safely
datTraits <- as.data.frame(datTraits, stringsAsFactors = FALSE)
# assign row names as a character vector from the first column of allTraits
rownames(datTraits) <- as.character(allTraits[[1]][traitRows])

# keep only ExpGroup in datTraits for simplicity
datTraits <- datTraits["ExpGroup"]

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
datTraits_consus <- create_datTraits_comparison(datTraits, "ExpGroup", c(con = 1, sus = 3))
datTraits_conres <- create_datTraits_comparison(datTraits, "ExpGroup", c(con = 1, res = 2))
datTraits_susres <- create_datTraits_comparison(datTraits, "ExpGroup", c(sus = 3, res = 2))

datTrait <- datTraits_susres # choose which comparison to use for module-trait correlation analysis

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
         cex.text = 0.4,
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
         cex.text = 0.4,
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

# Define variable weight containing the weight column of datTrait
ExpGroup = as.data.frame(datTraits_consus$ExpGroup)
names(ExpGroup) = "ExpGroup"

modNames = substring(names(mergedMEs), 3) #extract module names

#Calculate the module membership and the associated p-values
geneModuleMembership = as.data.frame(cor(expression.data, mergedMEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#Calculate the gene significance and associated p-values
geneTraitSignificance = as.data.frame(cor(expression.data, ExpGroup, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(ExpGroup), sep="")
names(GSPvalue) = paste("p.GS.", names(ExpGroup), sep="")
head(GSPvalue)

par(mar=c(1,1,1,1))
module = "orange"
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
abs(geneTraitSignificance[moduleGenes,1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for body weight",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

par(mar=c(1,1,1,1))
module = "orange"
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
abs(geneTraitSignificance[moduleGenes,1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for body weight",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

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

# save dendrogram as SVG in output_dir using the label
svg(file = file.path(output_dir, paste0("eigengene_dendrogram_", label, ".svg")), width = 8, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
            plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)

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

# save heatmap as SVG in output_dir using the label
svg(file = file.path(output_dir, paste0("eigengene_adjacency_heatmap_", label, ".svg")), width = 8, height = 6)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)

dev.off()

# Isolate weight and Glucose from the clinical traits
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

# check for modules of interest
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

moduleEigs <- mergedMEs$MEgrey60  # replace with appropriate module name
group <- datTraits$ExpGroup
boxplot(moduleEigs ~ group, names=c("Resilient","Susceptible"))

# Create group factor with defined levels and labels
group_factor <- factor(datTrait, levels = c(2, 3), labels = c("Resilient", "Susceptible"))

# Remove samples with NA in group_factor (control animals)
valid_samples <- !is.na(group_factor)

# Subset module eigengenes and group factor to valid samples only
moduleEigs_valid <- moduleEigs[valid_samples]
group_factor_valid <- group_factor[valid_samples]

# Plot boxplot with filtered data
boxplot(moduleEigs_valid ~ group_factor_valid, 
        main = "Module Eigengene by Group (Susceptible vs Resilient)",
        xlab = "Group", ylab = "Module Eigengene Expression")

table(group_factor)          # Check levels and counts including NAs
table(group_factor_valid)    # Check counts after filtering NA
length(moduleEigs_valid)     # Should match length of group_factor_valid
length(moduleEigs)           # Number of samples in expression data
length(datTrait)             # Should match length of moduleEigs

# auto-detect which datTrait is active and derive a label
detected_label <- "datTrait"
if (exists("datTrait")) {
  if (exists("datTraits_consus") && identical(datTrait, datTraits_consus)) detected_label <- "consus"
  else if (exists("datTraits_conres") && identical(datTrait, datTraits_conres)) detected_label <- "conres"
  else if (exists("datTraits_susres") && identical(datTrait, dataits_susres)) detected_label <- "susres"
}

# get numeric trait vector (support datTrait being a data.frame/single-column df or a vector)
if (exists("datTrait")) {
  trait_vector <- if (is.data.frame(datTrait)) {
    as.numeric(datTrait[[1]])
  } else {
    as.numeric(datTrait)
  }
} else if (exists("datTraits_consus")) {
  # fallback to consus if datTrait wasn't explicitly set
  trait_vector <- as.numeric(datTraits_consus[[1]])
  detected_label <- "consus"
} else {
  stop("No datTrait or datTraits_consus/conres/susres available to infer comparison.")
}

# mapping numeric codes to human labels (matches create_datTraits_comparison mapping: 1=con,2=res,3=sus)
code_to_label <- c("1" = "Control", "2" = "Resilient", "3" = "Susceptible")
present_vals <- sort(unique(na.omit(trait_vector)))
if (length(present_vals) < 2) stop("Trait vector does not contain at least two groups to compare.")
# build ordered group labels for the values present
group_labels <- unname(code_to_label[as.character(present_vals)])
group_labels[is.na(group_labels)] <- paste0("Group", present_vals[is.na(group_labels)])

# build factor for plotting (initial)
group_factor <- factor(trait_vector, levels = present_vals, labels = group_labels)

# REORDER FACTOR LEVELS SO: Control always left, Susceptible always right, others in between
preferred_order <- c("Control", "Resilient", "Susceptible")
# keep only the labels that are present, placing Control first and Sus last if present
middle <- setdiff(levels(group_factor), c("Control", "Susceptible"))
final_levels <- c(intersect(preferred_order[1], levels(group_factor))) # maybe Control
final_levels <- final_levels[!is.na(final_levels)]
# append any middle labels (preserve their current order)
if (length(middle) > 0) final_levels <- c(final_levels, middle)
# append Susceptible if present
if ("Susceptible" %in% levels(group_factor)) final_levels <- c(final_levels, "Susceptible")
# ensure final_levels unique and valid
final_levels <- unique(final_levels)
# if Control wasn't present (unlikely) but to be safe, fallback to existing levels
if (length(final_levels) == 0) final_levels <- levels(group_factor)

# Re-factor with the desired order
group_factor <- factor(as.character(group_factor), levels = final_levels)

# filter valid (non-NA) samples
valid_samples <- !is.na(group_factor)
moduleEigs_valid <- as.numeric(moduleEigs[valid_samples])
group_factor_valid <- droplevels(group_factor[valid_samples])

# infer module name (if module variable present use it, otherwise try to match moduleEigs to mergedMEs)
module_name <- "unknown"
if (exists("module") && is.character(module) && nzchar(module)) {
  module_name <- module
} else if (exists("mergedMEs") && !is.null(dim(mergedMEs)) && length(moduleEigs_valid) > 0) {
  me_sub <- mergedMEs[valid_samples, , drop = FALSE]
  corrs <- sapply(seq_len(ncol(me_sub)), function(i) cor(as.numeric(me_sub[, i]), moduleEigs_valid, use = "complete.obs"))
  nm <- names(corrs)[which.max(abs(corrs))]
  module_name <- if (!is.null(nm)) sub("^ME", "", nm) else "unknown"
}

# prepare plot directory & filename
plot_dir <- file.path(output_dir, "module_eigengene_plots", detected_label)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(plot_dir, paste0("module_eigengene_", module_name, "_boxplot_", detected_label, ".svg"))

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

p <- ggplot(df_plot, aes(x = Group, y = ME, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, colour = "black", fatten = 1, lwd = 1) +
  geom_jitter(color = "white", width = 0.12, size = 6, shape = 21, alpha = 0.8, stroke = 0.8) +
  scale_fill_manual(values = plot_colors) +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5, fill = "black", color = "black") +
  labs(x = NULL, y = paste0("Module Eigengene (", module_name, ")"), title = paste0("Module Eigengene by Group — ", detected_label)) +
  coord_cartesian(ylim = c(-0.3, 0.2)) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 20),
    axis.title.y = element_text(margin = margin(r = 20)),
    plot.title = element_text(face = "plain", size = 20),
    axis.line = element_line(size = 1, colour = "black"),
    axis.ticks = element_line(size = 0.6, colour = "black")
  )

# save SVG (prefer svglite if available)
ggsave(filename = outfile, plot = p, device = svglite::svglite, width = 4, height = 5, units = "in", dpi = 300)