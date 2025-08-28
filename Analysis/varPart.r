## ----setup, echo=FALSE, results="hide"------------------------------------------------------------
knitr::opts_chunk$set(
  tidy = FALSE, 
  cache = TRUE,
  fig.align = 'center',
  dev = "png",
  package.startup.message = FALSE,
  message = FALSE, 
  error = FALSE, 
  warning = TRUE
)
options(width = 100)

## ----simResult, cache=TRUE, fig.height=4, fig.width=4---------------------------------------------
# load library
library("variancePartition")
library("readxl")
library("tibble")
library("dplyr")

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
data(varPartData)
# view data varPartData


# load data from two excel files in directory

geneExpr <- read_excel("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/gene_expression.xlsx")
info <- read_excel("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info.xlsx")

# 1. Extract first column name (gene names column)
gene_id_col <- colnames(geneExpr)[1]

# 2. Replace NA or empty gene names with unique placeholders
geneExpr[[gene_id_col]] <- as.character(geneExpr[[gene_id_col]])
na_idx <- which(is.na(geneExpr[[gene_id_col]]) | geneExpr[[gene_id_col]] == "")
geneExpr[[gene_id_col]][na_idx] <- paste0("missing_gene_", seq_along(na_idx))

# 3. Optionally, make all gene names unique across the column
geneExpr[[gene_id_col]] <- make.unique(geneExpr[[gene_id_col]])

# 4. Now safely set as rownames and remove from columns
geneExpr <- geneExpr %>%
  column_to_rownames(var = gene_id_col)

# 5. Convert remaining columns to numeric
geneExpr[] <- lapply(geneExpr, function(x) as.numeric(as.character(x)))

# 6. Convert to matrix
geneExpr <- as.matrix(geneExpr)

# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual and Tissue are both categorical,
# so model them as random effects
# Note the syntax used to specify random effects
form <- ~ (1 | AnimalID) + (1 | celltype) + (1 | region) + (1 | plate) + (1 | layer) + (1 | region:layer) + (1 | celltype_layer) + (1 | celltype_ExpGroup) ExpGroup

# Fit model and extract results
# 1) fit linear mixed model on gene expression
# If categorical variables are specified,
#     a linear mixed model is used
# If all variables are modeled as fixed effects,
# 		a linear model is used
# each entry in results is a regression model fit on a single gene
# 2) extract variance fractions from each model fit
# for each gene, returns fraction of variation attributable
# 		to each variable
# Interpretation: the variance explained by each variables
# after correcting for all other variables
# Note that geneExpr can either be a matrix,
# and EList output by voom() in the limma package,
# or an ExpressionSet
varPart <- fitExtractVarPartModel(geneExpr, form, info)

# sort variables (i.e. columns) by median fraction
# 		of variance explained
vp <- sortCols(varPart)

# Figure 1a
# Bar plot of variance fractions for the first 10 genes
plotPercentBars(vp[1:10, ])

# Figure 1b
# violin plot of contribution of each variable to total variance
plotVarPart(vp)

## ----accessResults, cache=TRUE, warning=FALSE-----------------------------------------------------
# Access first entries
head(varPart)

# Access first entries for Individual
head(varPart$Individual)

# sort genes based on variance explained by Individual
head(varPart[order(varPart$Individual, decreasing = TRUE), ])

## ----savePlot, cache=TRUE, eval=FALSE-------------------------------------------------------------
# fig <- plotVarPart(vp)
# ggsave(file, fig)

## ----plotStratify, cache=TRUE, warning=FALSE, fig.height=4, fig.width=4---------------------------
# get gene with the highest variation across Tissues
# create data.frame with expression of gene i and Tissue
# 		type for each sample
i <- which.max(varPart$ExpGroup)
GE <- data.frame(Expression = geneExpr[i, ], ExpGroup = info$ExpGroup)

# Figure 2a
# plot expression stratified by ExpGroup
plotStratify(Expression ~ ExpGroup, GE, main = rownames(geneExpr)[i])

# get gene with the highest variation across Individuals
# create data.frame with expression of gene i and Tissue
# 		type for each sample
i <- which.max(varPart$AnimalID)
GE <- data.frame(
  Expression = geneExpr[i, ],
  AnimalID = info$AnimalID
)

# Figure 2b
# plot expression stratified by AnimalID
label <- paste("AnimalID:", format(varPart$AnimalID[i] * 100,
  digits = 3
), "%")
main <- rownames(geneExpr)[i]
plotStratify(Expression ~ AnimalID, GE,
  colorBy = NULL,
  text = label, main = main
)



## ----cache=TRUE-----------------------------------------------------------------------------------
library("lme4")

# fit regression model for the first gene
form_test <- geneExpr[1, ] ~ Age + (1 | Individual) + (1 | Tissue)
fit <- lmer(form_test, info, REML = FALSE)

# extract variance statistics
calcVarPart(fit)

## ----canCorPairs, cache=TRUE, results='hide', fig.width=5, fig.height=5---------------------------
form <- ~ Individual + Tissue + Batch + Age + Height

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, info)

# Plot correlation matrix
# between all pairs of variables
plotCorrMatrix(C)