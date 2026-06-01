# Output Contracts

Canonical generated outputs use dataset-aware paths:

```text
data/processed/<module>/<substep>/<dataset>/
results/tables/<module>/<substep>/<dataset>/
results/figures/<module>/<substep>/<dataset>/
results/source_data/<module>/<substep>/<dataset>/
results/reports/<module>/<substep>/<dataset>/
results/logs/<module>/<substep>/<dataset>/
```

Major downstream tables should be validated with `validate_table_schema(df, schema_name, strict = TRUE)` before writing. Active scripts should write:

- `run_manifest.yml`
- `sessionInfo.txt`
- config snapshot where applicable
- input file hashes where applicable

The manuscript export layer collects final outputs under:

```text
results/manuscript/figure_1/
results/manuscript/figure_2/
results/manuscript/extended_data/
results/manuscript/supplementary_tables/
results/manuscript/source_data/
```

These export scripts collect final products only; they do not recompute analyses.
