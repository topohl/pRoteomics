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


#d <- curl('https://raw.githubusercontent.com/fuzzyatelin/fuzzyatelin.github.io/master/bioanth-stats/module-F21-Group1/FemaleLiver-Data/LiverFemale3600.csv')
#liver.data <- read.csv(file = d, stringsAsFactors = FALSE, header = TRUE)
#head(liver.data)

d <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/male.data.xlsx"

male.data <- read_excel(path = d)

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


#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
abline(h = 40, col = "red");

cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10) #returns numeric vector
#Remove outlier
expression.data <- expression.data[cut.sampleTree==1, ]

spt <- pickSoftThreshold(expression.data)


par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")


par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")



softPower <- 7
adjacency <- adjacency(expression.data, power = softPower)

TOM <- TOMsimilarity(adjacency)

TOM.dissimilarity <- 1-TOM

#creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
labels = FALSE, hang = 0.04)


Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)

table(Modules) #returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module.

ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")


MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

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



plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
c("Original Module", "Merged Module"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors for original and merged modules")



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

# Define numbers of genes and samples
nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p") #p for pearson correlation coefficient 
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples) #calculate the p-value associated with the correlation

# Will display correlations and their p-values
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
par(mar = c(6, 8.5, 3, 1))
# Display the correlation values within a heatmap plot
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










## experimental

# Make sure row names match your expression data sample names
rownames(sample_info) <- sample_info$row.names

# ===== CREATE TRAIT MATRIX FROM YOUR METADATA =====
# Convert categorical variables to numeric for WGCNA correlation analysis

# Cell type traits (binary)
celltype_traits <- model.matrix(~ celltype - 1, data = sample_info)
colnames(celltype_traits) <- gsub("celltype", "", colnames(celltype_traits))

# Layer traits (binary)
layer_traits <- model.matrix(~ layer - 1, data = sample_info)
colnames(layer_traits) <- paste0("layer_", colnames(layer_traits))

# Region traits (binary) 
region_traits <- model.matrix(~ region - 1, data = sample_info)
colnames(region_traits) <- paste0("region_", colnames(region_traits))

# ReplicateGroup traits (binary)
replicate_traits <- model.matrix(~ ReplicateGroup - 1, data = sample_info)
colnames(replicate_traits) <- paste0("hemisphere_", colnames(replicate_traits))

# ExpGroup as continuous
sample_info$ExpGroup_numeric <- as.numeric(as.factor(sample_info$ExpGroup))

# Create final trait matrix
datTraits <- data.frame(
  # Cell type information
  celltype_traits,
  
  # Layer information
  layer_traits,
  
  # Region information  
  region_traits,
  
  # Hemisphere information
  replicate_traits,
  
  # Experimental group (continuous)
  ExpGroup = sample_info$ExpGroup_numeric,
  
  # Biologically meaningful combinations
  CA1_neuron = as.numeric(sample_info$region == "CA1" & sample_info$celltype == "neuron"),
  CA2_neuron = as.numeric(sample_info$region == "CA2" & sample_info$celltype == "neuron"), 
  CA3_neuron = as.numeric(sample_info$region == "CA3" & sample_info$celltype == "neuron"),
  CA1_microglia = as.numeric(sample_info$region == "CA1" & sample_info$celltype == "microglia"),
  CA2_microglia = as.numeric(sample_info$region == "CA2" & sample_info$celltype == "microglia"),
  CA3_microglia = as.numeric(sample_info$region == "CA3" & sample_info$celltype == "microglia"),
  
  # Cell layer specificity
  neuron_soma = as.numeric(sample_info$celltype_layer == "neuron_soma"),
  neuron_neuropil = as.numeric(sample_info$celltype_layer == "neuron_neuropil"),
  
  stringsAsFactors = FALSE
)

# Match with expression data samples
Samples <- rownames(expression.data)
traitRows <- match(Samples, rownames(sample_info))
datTraits <- datTraits[traitRows, ]
rownames(datTraits) <- Samples

# Remove any columns with all NAs or constant values
datTraits <- datTraits[, sapply(datTraits, function(x) length(unique(x[!is.na(x)]))) > 1]

# ===== CONTINUE WITH MODULE-TRAIT CORRELATION ANALYSIS =====
# Calculate module-trait correlations
nGenes = ncol(expression.data)
nSamples = nrow(expression.data)
module.trait.correlation = cor(mergedMEs, datTraits, use = "p")
module.trait.Pvalue = corPvalueStudent(module.trait.correlation, nSamples)

# Create heatmap
textMatrix = paste(signif(module.trait.correlation, 2), "\n(",
                   signif(module.trait.Pvalue, 1), ")", sep = "");
dim(textMatrix) = dim(module.trait.correlation)
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







## experimental











# Define variable weight containing the weight column of datTrait
ExpGroup = as.data.frame(datTraits$ExpGroup)
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
module = "skyblue"
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
abs(geneTraitSignificance[moduleGenes,1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for body weight",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


par(mar=c(1,1,1,1))
module = "purple"
column = match(module, modNames)
moduleGenes = mergedColors==module
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),
abs(geneTraitSignificance[moduleGenes,1]),
xlab = paste("Module Membership in", module, "module"),
ylab = "Gene significance for body weight",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$ExpGroup);
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


# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0, mar = c(1,1,1,1))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(5,5,2,2),
plotDendrograms = FALSE, xLabelsAngle = 90)







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