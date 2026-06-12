# Maintenance

`pipeline.yml` is the active source of truth for runnable workflow order, dataset support, inputs, outputs, and rerun safety. Documentation tables such as `docs/active_script_io_audit.tsv` are generated/audit snapshots and should not be treated as competing registries.

## Add A New Script

1. Put the script in the numbered module folder that owns the work.
2. Add a header with script purpose, stage, scope, consumed inputs, produced outputs, dataset behavior, and whether it changes core scientific state.
3. Use shared path helpers from `R/paths.R`; avoid machine-local absolute paths.
4. Support `--dry-run` for scripts that will be registered.
5. Write outputs under canonical `results/{tables,figures,source_data,reports,logs}/<module>/<substep>/<dataset>/` or `data/processed/<module>/<substep>/<dataset>/`.

## Register It

Add one entry to the appropriate `pipeline.yml` stage:

- `script`
- `required`
- `datasets`
- `stage`
- `scope`
- `produces`
- `consumes_required`
- `consumes_optional`
- `recomputes_core_state`
- `safe_downstream_rerun`
- `notes`

Use `datasets: ["global"]` and `scope: "global"` only for scripts that are intentionally dataset-agnostic. Preserve existing stage names unless you are doing a coordinated migration.

## Declare Inputs And Outputs

Required inputs should be files or directories that must exist for the script to make valid outputs. Optional inputs should be allowed to be missing and should produce explicit status rows, warnings, or conservative interpretation flags.

Output paths in `pipeline.yml` should be stable contracts, not temporary scratch files. Add manuscript-facing tables to `docs/file_contracts.tsv` when another script or reviewer should depend on them.

## Add Output Validation

Use existing validators in `R/module_contracts.R`, `R/schema_validation.R`, or `R/validation_utils.R` where possible. For new final-facing tables, add a small validator that checks required columns and basic ranges without changing the underlying scientific calculations.

## Run Checks

Use these before committing workflow changes:

```bash
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset all --stage all --dry-run
Rscript tests/smoke_test_active_script_contracts.R
Rscript tests/smoke_test_file_contracts.R
```

For a focused downstream rerun, prefer:

```bash
Rscript run_dataset_pipeline.R --dataset microglia --stage modules_downstream --dry-run
Rscript run_dataset_pipeline.R --dataset microglia --stage modules_downstream
```

## Legacy Material

Do not run scripts from `90_testing/`, `99_deprecated/`, or `09_pride_submission/` for the active workflow unless a migration note explicitly says so. Active export code lives in `09_export_pride_journal/`, and generated deposition payloads belong under `pride_submission/`.
