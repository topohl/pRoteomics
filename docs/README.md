# Documentation Map

`pipeline.yml` is the source of truth for active script order, datasets, required
inputs, optional inputs, and expected outputs. Keep this folder as a small set of
reference documents rather than a second pipeline registry.

## Start Here

- [Workflow](../WORKFLOW.md): plain-language scientific workflow.
- [Run order](../RUN_ORDER.md): human-readable execution guide for the registry.
- [WGCNA workflow](WGCNA_WORKFLOW.md): WGCNA construction, scoring, modelling, annotation, and summary layers.
- [Datasets](DATASETS.md): dataset families, allowed analyses, and guardrails.
- [Input contracts](INPUT_CONTRACTS.md): expected private/local inputs.
- [Output contracts](OUTPUT_CONTRACTS.md): canonical output roots and final tables.
- [Reviewer reproducibility](REVIEWER_REPRODUCIBILITY.md): data-independent checks.

## Manuscript And Interpretation

- [Biological interpretation](BIOLOGICAL_INTERPRETATION.md): claims-table guardrails.
- [Biological interpretation guide](biological_interpretation_guide.md): broader biology/program interpretation notes.
- [Microglia ROI interpretation](MICROGLIA_ROI_INTERPRETATION.md): wording for microglia-enriched ROI results.

## Export And Submission

- [PRIDE export](PRIDE_EXPORT.md): active PRIDE/journal export module.
- [PRIDE submission workflow](PRIDE_submission_workflow.md): practical staging checklist.

## Audit And Maintenance

- [File contracts](file_contracts.tsv): compact producer/consumer table for stable objects.
- [Active script I/O audit](active_script_io_audit.tsv): historical audit snapshot unless explicitly refreshed; useful for cleanup, not the active registry.
- [Naming migration](NAMING_MIGRATION.md): removed aliases and canonical replacements.
- [Refactor backlog](refactor_backlog.md): non-blocking follow-up work.

`current_data_flow.md`, `hardcoded_path_inventory.md`, and
`active_script_io_audit.tsv` are audit/historical snapshots unless refreshed.
Keep them for context, but update `pipeline.yml`, `WORKFLOW.md`,
`RUN_ORDER.md`, and the contract docs first when behavior changes.
