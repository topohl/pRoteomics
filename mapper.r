# Install if needed
if (!requireNamespace("HGNChelper", quietly = TRUE)) {
    install.packages("HGNChelper")
}
library(HGNChelper)

# Set directories
working_dir <- "/Users/tobiaspohl/Documents/clusterProfiler"
results_dir <- file.path(working_dir, "Results")
setwd(working_dir)

# Define cell types
cell_types <- c("mcherryG1", "mcherryG2")
file_name <- paste0(paste(cell_types, collapse = "_"), ".csv")
data_path <- file.path(working_dir, "Datasets", file_name)

# Read data
data <- read.csv(data_path, header = TRUE, row.names = 1)


# Suppose your gene symbols are in a vector called gene_symbols
gene_symbols <- c("N42L1", "P33MX", "RIGI", "A1AG1", "FBLN5", "Q792Z1", "Birc5", "Abl1")

# For mouse, set species = "mouse"
results <- checkGeneSymbols(gene_symbols, species = "mouse")

# View results
print(results)