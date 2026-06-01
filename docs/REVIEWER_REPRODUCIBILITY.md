# Reviewer Reproducibility

Reviewers can inspect the active pipeline without private raw data:

```bash
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset neuron_neuropil --stage qc --dry-run
Rscript run_dataset_pipeline.R --dataset microglia --stage enrichment --dry-run
Rscript tests/testthat.R
```

The pipeline registry validates that active scripts exist and that scripts mentioned in `RUN_ORDER.md` are either active or explicitly legacy/deprecated. Table contracts live in `inst/schemas/`, and active dataset capabilities live in `R/dataset_config.R`.

For package restoration:

```r
install.packages("renv")
renv::restore()
```

Dry-run failures should be resolved before full analyses. Full scientific runs require the private source matrices, metadata, and raw/vendor files described in the data availability statement.
