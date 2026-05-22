# Generate an SDRF-Proteomics-style metadata table.
#
# This script is intentionally conservative. It only fills columns that can be
# inferred from sample_metadata_clean.*. Missing values are left blank and are
# reported by 04_validate_pride_package.R.

source("R/paths.R")
source("R/pride_helpers.R")

ensure_pride_dirs()
meta <- read_sample_metadata()

out_file <- repo_path("PRIDE_package", "00_metadata", "sdrf_proteomics.tsv")

if (nrow(meta) == 0) {
  sdrf <- data.frame(
    `source name` = character(),
    `characteristics[organism]` = character(),
    `characteristics[organism part]` = character(),
    `characteristics[sex]` = character(),
    `characteristics[strain]` = character(),
    `characteristics[brain region]` = character(),
    `characteristics[hippocampal layer]` = character(),
    `factor value[group]` = character(),
    `factor value[sex]` = character(),
    `factor value[region_layer]` = character(),
    `assay name` = character(),
    `comment[data file]` = character(),
    `comment[instrument]` = character(),
    `comment[fraction identifier]` = character(),
    `comment[technical replicate]` = character(),
    `comment[label]` = character(),
    `comment[search engine]` = character(),
    `comment[search engine version]` = character(),
    `comment[protein database]` = character(),
    check.names = FALSE
  )
  write_tsv(sdrf, out_file)
  message("No sample metadata found. Wrote empty SDRF template: ", out_file)
} else {
  sample_id <- coalesce_column(meta, c("sample_id", "sample", "source_name", "source", "animal_id"))
  animal_id <- coalesce_column(meta, c("animal_id", "animal", "mouse_id", "mouse"))
  sex <- coalesce_column(meta, c("sex", "gender"))
  group <- coalesce_column(meta, c("group", "condition", "phenotype", "treatment"))
  region <- coalesce_column(meta, c("region", "brain_region", "hippocampal_region"))
  layer <- coalesce_column(meta, c("layer", "hippocampal_layer", "stratum"))
  region_layer <- coalesce_column(meta, c("region_layer", "regionlayer", "region_layer_id"))
  raw_file <- coalesce_column(meta, c("raw_file_name", "raw_file", "data_file", "filename", "file_name"))
  assay <- coalesce_column(meta, c("assay_name", "assay", "run", "run_id"))
  instrument <- coalesce_column(meta, c("instrument", "ms_instrument", "mass_spectrometer"))
  fraction <- coalesce_column(meta, c("fraction", "fraction_id", "fraction_identifier"))
  tech_rep <- coalesce_column(meta, c("technical_replicate", "technical_rep", "techrep"))
  label <- coalesce_column(meta, c("label", "isobaric_label", "tmt_channel"), default = "label free sample")
  search_engine <- coalesce_column(meta, c("search_engine", "software", "search_software"))
  search_engine_version <- coalesce_column(meta, c("search_engine_version", "software_version"))
  fasta <- coalesce_column(meta, c("fasta_file", "protein_database", "database"))
  strain <- coalesce_column(meta, c("strain", "mouse_strain"))

  source_name <- ifelse(!is.na(sample_id) & sample_id != "", sample_id, animal_id)
  assay_name <- ifelse(!is.na(assay) & assay != "", assay, source_name)

  sdrf <- data.frame(
    `source name` = source_name,
    `characteristics[organism]` = "Mus musculus",
    `characteristics[organism part]` = "brain",
    `characteristics[sex]` = sex,
    `characteristics[strain]` = strain,
    `characteristics[brain region]` = region,
    `characteristics[hippocampal layer]` = layer,
    `factor value[group]` = group,
    `factor value[sex]` = sex,
    `factor value[region_layer]` = ifelse(!is.na(region_layer) & region_layer != "", region_layer, paste(region, layer, sep = "_")),
    `assay name` = assay_name,
    `comment[data file]` = raw_file,
    `comment[instrument]` = instrument,
    `comment[fraction identifier]` = fraction,
    `comment[technical replicate]` = tech_rep,
    `comment[label]` = label,
    `comment[search engine]` = search_engine,
    `comment[search engine version]` = search_engine_version,
    `comment[protein database]` = fasta,
    check.names = FALSE
  )

  write_tsv(sdrf, out_file)
  message("Wrote SDRF table: ", out_file)
}
