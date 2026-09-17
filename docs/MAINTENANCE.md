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

Use existing validators in `R/data_contracts/module_contracts.R`, `R/utilities/schema_validation.R`, or `R/utilities/validation_utils.R` where possible. For new final-facing tables, add a small validator that checks required columns and basic ranges without changing the underlying scientific calculations.

## Run Checks

Use these before committing workflow changes:

```bash
Rscript run_dataset_pipeline.R --list-stages
Rscript run_dataset_pipeline.R --dataset all --stage all --dry-run
Rscript tests/smoke_test_active_script_contracts.R
Rscript tests/smoke_test_file_contracts.R
Rscript tools/audit_active_scripts.R
Rscript audits/verify_scientific_contracts.R
```

`audits/verify_scientific_contracts.R` checks the named scientific contracts
(bilateral averaging, CA2-SLM QC, the module identity contract and the rest)
against the code that carries them. It takes no arguments and has no callers by
design: run it by hand after touching anything those contracts describe.

## Regenerate The Architecture Records

These write the derived architecture documents. They read `pipeline.yml`, so
run them after changing the registry, and commit the result:

```bash
Rscript tools/generate_architecture_docs.R    # docs/ANALYSIS_ENTRYPOINTS.md
Rscript tools/generate_results_ownership.R    # config/results_ownership.csv + docs/RESULTS_OWNERSHIP.md
Rscript tools/generate_active_code_tree.R     # audits/phase6e_final_active_tree.csv
Rscript tools/normalize_script_headers.R      # fills missing script-header fields
Rscript tools/generate_legacy_output_registry.R  # config/legacy_output_registry.csv
```

Run the legacy registry generator after adding or removing a registered
writer: a root becomes legacy the moment nothing declares output beneath it,
and `tests/testthat/test-output-layout-contract.R` holds writers to the result.

`tools/build_restructure_migration_map.R` belongs to the repository-split
migration rather than to routine maintenance. It rebuilds
`audits/restructure_migration_map.csv` from the recorded path, boundary and
test-repoint contracts, and is the tool to run if that record needs refreshing.

## Which Document Answers What

One authority per question, so they cannot contradict each other:

| Question | Authority |
| --- | --- |
| What order do scripts run in, and what does each declare? | `pipeline.yml` |
| What command do I type? | `RUN_ORDER.md` |
| Which script is the entry point for an analysis area, and what does it depend on? | `docs/ANALYSIS_ENTRYPOINTS.md` (generated) |
| Which script produces a given manuscript result? | `docs/CANONICAL_ANALYSIS_ENTRYPOINTS.md` |
| Who is the canonical owner of a result family, and who else writes into it? | `config/results_ownership.csv`, rendered as `docs/RESULTS_OWNERSHIP.md` (generated) |
| What role does a given active file play? | `audits/phase6e_final_active_tree.csv` (generated) |
| What was a file called before the Phase 6E renames? | `audits/phase6e_naming_migration.csv` |
| Where does an output go, and what are `work/`, `results/` and `exports/`? | `docs/OUTPUT_LAYOUT.md`, from `config/output_layout.yml` |
| Which output roots are read-only, and why? | `config/legacy_output_registry.csv` (generated) |
| Where was an output written before Phase 6F? | `audits/phase6f_output_inventory.csv` |
| Should the repository be renamed? | `docs/REPOSITORY_NAME_RECOMMENDATION.md` |

For a focused downstream rerun, prefer:

```bash
Rscript run_dataset_pipeline.R --dataset microglia --stage modules_downstream --dry-run
Rscript run_dataset_pipeline.R --dataset microglia --stage modules_downstream
```

## Legacy Material

Do not run scripts from `90_testing/`, `99_deprecated/`, or `09_pride_submission/` for the active workflow unless a migration note explicitly says so. Active export code lives in `09_export_pride_journal/`, and generated deposition payloads belong under `pride_submission/`.
