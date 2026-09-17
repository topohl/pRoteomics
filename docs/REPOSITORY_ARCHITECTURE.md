# Repository architecture

Post Phase 6B/6C. This repository owns scientific inference. Manuscript prose
and journal figure assembly live in the sibling `Exp9_manuscript` repository.

See [RESTRUCTURE_PLAN.md](RESTRUCTURE_PLAN.md) for the migration record and the
equivalence oracle, [ANALYSIS_ENTRYPOINTS.md](ANALYSIS_ENTRYPOINTS.md) for which
script owns which analysis, and [RESULTS_OWNERSHIP.md](RESULTS_OWNERSHIP.md) for
which analysis produces which result family.

## Layers

| Layer | Contains | Does not contain |
| --- | --- | --- |
| `R/` | reusable scientific functions, grouped by domain | runnable analysis workflows |
| `analysis/` | explicit runnable analysis entrypoints | reusable libraries, manuscript rendering |
| `config/` | frozen scientific configuration and contracts | output |
| `data/` | raw, metadata and reference inputs | derived results |
| `results/` | canonical analysis products | code |
| `tests/` | scientific and architectural regression guards | analysis |
| `audits/` | provenance, robustness and migration audits | analysis entrypoints |
| `tools/` | maintenance, export and reference-audit utilities | analysis entrypoints |
| `archive/` | superseded code kept for provenance | anything canonical or executable by the registry |

## R library addressing

`R/` is organised by domain:

```
R/
├── paths.R              bootstrap: repository root, path helpers, the resolver
├── null_coalescing.R    bootstrap dependency of paths.R
├── data_contracts/      dataset, identifier, module and spatial-systems contracts
├── qc/                  QC, compartment-marker and abundance helpers
├── statistics/          module statistics, WGCNA, evidence and audit workbooks
├── spatial/             spatial identity, atlas, CA2-SLM robustness
├── enrichment/          enrichment IO, GO themes, EWCE, clusterProfiler
├── networks/            network construction and position helpers
└── utilities/           registry, validation, export, schema, plotting
```

Libraries are addressed **by bare name**, never by physical path:

```r
source(repo_path("R", "module_stats.R"))   # resolves to R/statistics/module_stats.R
```

`repo_path("R", "<name>.R")` delegates to `r_library_path()`, which looks for a
flat `R/<name>.R` first and then searches one level deep. This is why the domain
layout can be changed again without editing the 811 call sites and 97 test files
that name these libraries. `paths.R` and `null_coalescing.R` stay at the `R/`
root because they define the resolver and cannot use it.

Tests use the same addressing. A test that builds `R/<lib>.R` by hand would
re-introduce the coupling this resolver exists to remove.

## Analysis stage identities

Stage directories describe function. The number is a stable identifier, not an
execution order: `pipeline.yml` remains the sole execution-order authority.

| Identity | Owns |
| --- | --- |
| `analysis/preprocessing` | preprocessing handoff and protein/gene identifier mapping |
| `analysis/qc` | QC, missingness, marker fidelity, confounding |
| `analysis/spatial_validation` | spatial systems, bilateral aggregation, CA2-SLM robustness |
| `analysis/differential_abundance` | differential abundance and enrichment |
| `analysis/wgcna` | WGCNA module and supermodule construction and downstream |
| `analysis/enrichment` | cell-type enrichment |
| `analysis/spatial_networks` | spatial and differential networks |
| `analysis/integration` | biological integration and behaviour/physiology coupling |
| `analysis/publication_source_data` | PRIDE and publication source-data export |

## Output namespaces are deliberately not migrated

Paths under `results/` and `data/` are keyed on **stage identity**, not on
script location, and they carry 129 frozen untracked baseline objects. Moving a
script therefore does not move its outputs, and a stage name appearing inside
`results/tables/06_modules_WGCNA/` is correct rather than stale.

This is the single most load-bearing invariant of the migration. Every
repointing step verified it explicitly.

## Publication boundary

```
analysis/**                      scientific inference
   |
results/source_data/manuscript/  canonical per-identity source data
   |
tools/export_publication_source_data.R
   |
results/publication_source_data/<publication_id>/ + manifest.csv
   |
   v
Exp9_manuscript                  renders from a frozen, hash-verified copy
```

`config/publication_source_data_contract.yml` names the canonical publication
identities and the withheld ones. It is a scientific-side contract, so building
the bundle never requires reading the manuscript repository, and the manuscript
repository never reads a live path here.

Phase 6D closed the last two gaps. Journal figure naming, layout, per-identity
source data, legends and the submission bundle are owned by
`Exp9_manuscript/tools/package_journal_figures.R`, which verifies every
assembled artefact against the publication registry hash before packaging it
and refuses to package one that does not match.

What remains here is `analysis/publication_source_data/08_export_manuscript_figures.R`.
It scans this repository's figure outputs, records for each whether it is an
editable vector, whether it has a PNG companion and whether sibling source data
exists, and stages the candidates. That is a publication-readiness audit of
scientific output, which is a scientific-repository responsibility.

The script keeps its name deliberately. It is listed in
`freeze_protected_export_files()`, whose provenance equivalence check compares
blobs across two historical commits, so renaming it broke that check when
tried. A clearer filename is not worth perturbing a frozen mechanism, and the
ownership boundary is carried by what the code does and by the tests rather
than by the filename.

`tools/audit_cross_repo_boundary.R` measures the result. It strips comments and
requires a mention to sit inside something that actually resolves a path or
loads code before counting it, so provenance text is not mistaken for coupling.
Runtime live cross-repository dependencies in either direction are zero.

`tests/testthat/test-analysis-publication-boundary.R` enforces the boundary:
no `figures/`, no `manuscript/`, no generation-named panel library, no
registered renderer, a complete and hash-exact source-data manifest, and no
active path carrying a `final`/`v7`/`v8`/`v9`/`latest` namespace.

## Version history belongs in git, not in filenames

Active paths do not carry generation names. Historical generations
(`final_truth_v9`, `spatial_v6`, `editorial_v8`, `nature_final_v7`,
`manuscript_figures_v2`, `ED1_FINAL_V9`) survive only in:

- manuscript provenance records, now in `Exp9_manuscript/provenance/`;
- `archive/`, which is excluded from script discovery;
- `docs/NAMING_MIGRATION.md`, whose left column is intentionally historical;
- `tools/restructure_pipeline_folders.sh`, which records an earlier migration.

Two exceptions are deliberate and documented in
[RESTRUCTURE_PLAN.md](RESTRUCTURE_PLAN.md): the untracked output namespace
`results/**/manuscript_candidates/final_truth_v9/` (PB-11) and the
`spatial_v6` fingerprint source tables (PB-12). Both are frozen untracked
objects; renaming them would mean bulk-moving the historical result tree, which
the migration scope forbids.
