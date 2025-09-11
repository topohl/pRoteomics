# Install and load required packages
packages <- c("readxl", "writexl", "dplyr")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
    install.packages(packages[!installed])
}
lapply(packages, library, character.only = TRUE)

# Define the file path
file_path <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/sample_metadata_R.xlsx"

# Get available sheet names
sheet_names <- excel_sheets(file_path)
print(sheet_names)

# Load the sheet named "samples" into a data frame
df <- read_excel(file_path, sheet = "samples")
head(df)

# extract info from column "sample_id", sample_id have this format: D:\Proteomics\Fabian\Tobi\Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A0003_L_CA1_Microglia_S113_S2-A3_1_13026.d
# in Bluto_20250703_FCo_Evo2_80SPDzoom_Tobias_A0003_L_CA1_so_Neuron_S113_S2-A3_1_13026.d 
# A0003 is AnimalID, _L_ or _R_ indicate left or right replicates and go into "ReplicateGroup".
# _13026 is the run number and goes into "run_order"
df <- df %>%
  mutate(
    AnimalID = sub(".*_(A[0-9]+)_.*", "\\1", sample_id),
    ReplicateGroup = ifelse(grepl("_L_", sample_id), "Left", "Right"),
    run_order = sub(".*_(\\d+)\\.d$", "\\1", sample_id),
    group = sub(".*_(?:A[0-9]+)_(?:L|R)_([^_]+(?:_[^_]+)*)_S\\d+_.*", "\\1", sample_id),
    exclude = ifelse(grepl("CA3_slm", sample_id), TRUE, FALSE),
    plate = ifelse(grepl("20250703", sample_id), "B",
             ifelse(grepl("20250707", sample_id), "C", NA)),
    sample_number = sub(".*_(S\\d+)_.*", "\\1", sample_id),
    sample_location = ifelse(grepl("_[S]\\d+-", sample_id),
                             sub(".*_(S\\d+)-([A-Z]\\d+)_.*", "\\2", sample_id),
                             sub(".*_(S\\d+)_.*", "\\1", sample_id)),
    shortname = paste(plate, sample_number, sample_location, run_order, sep="_"),
    group2 = paste(group, ReplicateGroup, sep = "_"),
    region = sub(".*_(L|R)_((?:CA1|CA2|CA3|DG))_.*", "\\2", sample_id),
    layer = ifelse(grepl(".*_(?:L|R)_[^_]+_((?:so|sp|sr|slm|mo|po|sg))_.*", sample_id),
             sub(".*_(?:L|R)_[^_]+_((?:so|sp|sr|slm|mo|po|sg))_.*", "\\1", sample_id),
             "NA"),
    celltype = ifelse(
      grepl("_(?:L|R)_[^_]+_Microglia", sample_id, perl = TRUE),
      "microglia",
      ifelse(
      grepl("_(?:L|R)_[^_]+_[^_]+_Neuron", sample_id, perl = TRUE),
      "neuron",
      NA
      )
    ),
    celltype_layer = case_when(
      grepl("_(sp|sg)_Neuron", sample_id) ~ "neuron_soma",
      grepl("_Neuron", sample_id) ~ "neuron_neuropil",
      grepl("Microglia", sample_id) ~ "NA",
      TRUE ~ "NA"
    )
  )

# save as excel file using writexl package
write_xlsx(df, "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/proteomics/sample_metadata_R_processed.xlsx")