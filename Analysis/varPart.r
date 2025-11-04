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
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
# load CRAN packages
pacman::p_load(readxl, tibble, dplyr, svglite)

# Ensure variancePartition is installed (use BiocManager for Bioconductor package),
# then load it with library() to avoid pacman's argument-parsing ambiguities.
if (!requireNamespace("variancePartition", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install("variancePartition", ask = FALSE, update = FALSE)
}
library(variancePartition)

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
#data(varPartData)
# view data varPartData

# output_dir
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/output"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# load data from two excel files in directory

geneExpr <- read_excel("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/gene_expression.xlsx")
info <- read_excel("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/msdap/variancePartition/data/sample_info.xlsx")
rownames(info) <- info$row.names

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

# make geneExpr a tibble
#geneExpr <- as_tibble(geneExpr)

# Specify variables to consider
# Age is continuous so model it as a fixed effect
# Individual and Tissue are both categorical,
# so model them as random effects
# Note the syntax used to specify random effects
form <- ~ celltype + region + layer + group + ReplicateGroup + (1|AnimalID) + (1|plate)

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

info <- as.data.frame(info)
rownames(info) <- info$row.names

info$ReplicateGroup <- factor(info$ReplicateGroup)
info$celltype       <- factor(info$celltype)
info$region         <- factor(info$region)
info$layer          <- factor(info$layer)
info$ExpGroup       <- factor(info$ExpGroup)

# 3. Define a minimal, non‐collinear formula:
form <- ~ (1 | celltype) + (1|region) + (1|layer) + (1|ReplicateGroup) + (1|ExpGroup) +
         (1|AnimalID) + (1|plate)

# now run
varPart <- fitExtractVarPartModel(geneExpr, form, info)

# sort variables (i.e. columns) by median fraction
# 		of variance explained
vp <- sortCols(varPart)

# Figure 1a
# Bar plot of variance fractions for the first 10 genes
plotPercentBars(vp[1:10, ])

# save in output dir
ggsave(file.path(output_dir, "Figure1a.svg"), plot = last_plot(), device = "svg")

# Figure 1b
# violin plot of contribution of each variable to total variance
plotVarPart(vp)

# save in output dir
ggsave(file.path(output_dir, "Figure1b.svg"), plot = last_plot(), device = "svg")

## ----accessResults, cache=TRUE, warning=FALSE-----------------------------------------------------
# Access first entries
head(varPart)

# Access first entries for Individual
head(varPart$ExpGroup)

# sort genes based on variance explained by Individual
head(varPart[order(varPart$ExpGroup, decreasing = TRUE), ])

## ----savePlot, cache=TRUE, eval=FALSE-------------------------------------------------------------
# fig <- plotVarPart(vp)
# ggsave(file, fig)

## ----plotStratify, cache=TRUE, warning=FALSE, fig.height=4, fig.width=4---------------------------
# get gene with the highest variation across Tissues
# create data.frame with expression of gene i and Tissue
# 		type for each sample

# define the comparison factor
# choose stratifying variable name
strat_var <- "AnimalID"

# find gene with highest variation across that variable
i <- which.max(varPart[[strat_var]])

# construct data.frame; ensure Expression is a numeric vector across samples
GE <- data.frame(Expression = as.numeric(geneExpr[i, ]), factor = info[[strat_var]], stringsAsFactors = FALSE)

# Figure 2a
# plot expression stratified by the chosen variable
main_label <- if (!is.null(rownames(geneExpr))) rownames(geneExpr)[i] else paste0("Gene_", i)
plotStratify(Expression ~ factor, GE, main = main_label)

# save the plot, include the stratifying variable name in the filename
ggsave(file.path(output_dir, paste0("Figure2a_", strat_var, ".svg")), plot = last_plot(), device = "svg")

# get gene with the highest variation across Individuals
# create data.frame with expression of gene i and Tissue
# 		type for each sample
i <- which.max(varPart[[strat_var]])
GE <- data.frame(
    Expression = as.numeric(geneExpr[i, ]),
    stringsAsFactors = FALSE
)
# add the stratifying variable as a real column named by strat_var
GE[[strat_var]] <- info[[strat_var]]

# Figure 2b
# plot expression stratified by ExpGroup
label <- paste(strat_var, ":", format(varPart[[strat_var]][i] * 100,
    digits = 3
), "%")
main <- if (!is.null(rownames(geneExpr))) rownames(geneExpr)[i] else paste0("Gene_", i)
# construct formula from the variable name so the column in GE is used
plotStratify(as.formula(paste("Expression ~", strat_var)), GE,
    colorBy = NULL,
    text = label, main = main
)

# save the plot
ggsave(file.path(output_dir, paste0("Figure2b_", strat_var, ".svg")), plot = last_plot(), device = "svg")

## ----cache=TRUE-----------------------------------------------------------------------------------
library("lme4")

# fit regression model for the first gene
form_test <- geneExpr[1, ] ~ ExpGroup + (1 | AnimalID) + celltype
fit <- lmer(form_test, info, REML = FALSE)

# extract variance statistics
calcVarPart(fit)

## ----canCorPairs, cache=TRUE, results='hide', fig.width=5, fig.height=5---------------------------
form <- ~ (1 | AnimalID) + celltype + (1 | region) + (1 | plate) + (1 | layer) + ExpGroup + ReplicateGroup

# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C <- canCorPairs(form, info)

# Plot correlation matrix
# between all pairs of variables
plotCorrMatrix(C)

# save the plot
ggsave(file.path(output_dir, paste0("Figure3a_CCA.svg")), plot = last_plot(), device = "svg")