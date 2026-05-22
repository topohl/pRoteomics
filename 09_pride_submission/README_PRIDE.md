# PRIDE / ProteomeXchange submission support

This folder adds a reproducibility and packaging layer for preparing a local PRIDE/ProteomeXchange submission package from the pRoteomics analysis outputs.

The repository should contain code, templates, and documentation. Large mass spectrometry files should not be committed to GitHub. Raw/vendor files and full processed search outputs should be assembled locally in `PRIDE_package/` and uploaded through the PRIDE submission workflow.

## Local package structure

Create or populate this local folder:

```text
PRIDE_package/
  00_metadata/
  01_raw/
  02_search_results/
    search_parameters/
    fasta_database/
  03_processed_quantification/
  04_downstream_analysis/
  05_scripts/
```

`PRIDE_package/` is ignored by Git because it may contain large raw data files.

## Expected file classes

Use the following practical file classes in `pride_file_manifest.tsv`:

| File class | Meaning |
| --- | --- |
| `RAW` | Vendor raw files or mzML files from the mass spectrometer |
| `PEAK` | Peak list files such as MGF/MS2, if generated |
| `SEARCH` | Search-engine outputs, e.g. MaxQuant, DIA-NN, Spectronaut, Proteome Discoverer outputs |
| `RESULT` | Processed matrices, differential tables, enrichment tables, EWCE/WGCNA/network outputs |
| `FASTA` | Protein sequence database used for searching |
| `PARAMETERS_FILE` | Search parameters, software settings, `mqpar.xml`, config exports |
| `EXPERIMENTAL_DESIGN` | SDRF, sample metadata, manifest, sample-to-file mapping |
| `OTHER` | Additional documentation or auxiliary files |

## Minimal input metadata

The automation looks for a cleaned sample metadata table at one of these locations:

```text
results/analysis_ready/sample_metadata_clean.tsv
results/analysis_ready/sample_metadata_clean.csv
data/metadata/sample_metadata_clean.tsv
data/metadata/sample_metadata_clean.csv
01_preprocessing/output/sample_metadata_clean.tsv
01_preprocessing/output/sample_metadata_clean.csv
01_preprocessing/sample_metadata_clean.tsv
01_preprocessing/sample_metadata_clean.csv
```

Recommended columns:

```text
sample_id
animal_id
sex
group
region
layer
region_layer
raw_file_name
assay_name
instrument
search_engine
search_engine_version
fasta_file
fraction
technical_replicate
matrix_column_name
```

The critical column is `raw_file_name`. Without it, the SDRF can be generated only as a partial template and validation will fail.

## How to run

From the repository root:

```r
source("09_pride_submission/00_make_pride_package.R")
```

This will create or update:

```text
PRIDE_package/00_metadata/sdrf_proteomics.tsv
PRIDE_package/00_metadata/pride_file_manifest.tsv
PRIDE_package/00_metadata/checksum.md5
PRIDE_package/00_metadata/validation_report.tsv
```

## Interpretation of validation results

`validation_report.tsv` uses three statuses:

| Status | Meaning |
| --- | --- |
| `PASS` | Check passed |
| `WARN` | Package may still be usable, but manual review is needed |
| `FAIL` | Package is not ready for upload or is missing required traceability |

Typical expected early failures are:

- no raw files in `PRIDE_package/01_raw/`
- no search or result files in the package
- no `sample_metadata_clean.tsv/csv`
- no `raw_file_name` column
- raw files listed in metadata but not present in `PRIDE_package/01_raw/`

## What still needs manual review

The scripts cannot infer all PRIDE-relevant experimental metadata. Before upload, manually verify:

- instrument model
- acquisition method
- search engine and version
- FASTA/database release
- enzyme and digestion parameters
- fixed and variable modifications
- precursor and fragment tolerances
- labeling or label-free status
- fractionation information
- biological and technical replicate structure
- whether all raw files are mapped to samples
- whether all processed matrices can be traced back to sample metadata

## Practical workflow

1. Run the normal pRoteomics analysis pipeline.
2. Export or copy the clean sample metadata to `results/analysis_ready/sample_metadata_clean.tsv`.
3. Copy raw files into `PRIDE_package/01_raw/`.
4. Copy search-engine outputs into `PRIDE_package/02_search_results/`.
5. Copy processed protein/peptide/precursor matrices into `PRIDE_package/03_processed_quantification/`.
6. Copy final downstream tables into `PRIDE_package/04_downstream_analysis/`.
7. Run `source("09_pride_submission/00_make_pride_package.R")`.
8. Inspect `validation_report.tsv`.
9. Fix all `FAIL` entries before upload.

## Boundary between GitHub and PRIDE

GitHub should contain:

- analysis scripts
- PRIDE packaging scripts
- metadata templates
- documentation
- small example or template files

PRIDE should contain:

- raw MS data
- full search-engine outputs
- processed quantitative matrices
- SDRF and sample metadata
- checksums
- final downstream result tables needed to reproduce figures and supplements
