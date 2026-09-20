# Phase 6G — output-namespace migration, completion record

Phase 6G moved every analysis domain in this repository off the historical
numbered-stage output layout and onto the canonical namespaces declared in
`config/output_layout.yml`. No scientific value, statistic, figure content or
manuscript claim changed at any point; the migration is an addressing change
with byte-level evidence.

## Canonical namespace policy

```
work/<domain>/<analysis_id>/<scope>/            regenerable, never cited
results/<domain>/<analysis_id>/<scope>/<child>/ canonical scientific results
exports/<bundle>/                               frozen outward-facing copies
```

`child` is one of `tables`, `plots`, `models`, `manifests`, `reports`, with
`tables/source_data` for figure source data.

Three further destinations are declared rather than incidental, and a writer
targeting one is obeying a contract, not evading it:

| destination | declared by | role |
|---|---|---|
| `exports/publication_source_data/` | `output_layout.yml` (`may_be_canonical: false`) | frozen outward bundle; a copy of canonical results, the only thing `Exp9_manuscript` imports |
| `results/manuscript/` | `output_namespaces.yml` `manuscript_export_root`, rule `exporters_write_only_to_manuscript_export_root` | curated manuscript figure and source-data staging |
| `pride_submission/` | allowed non-result root in the output-namespace classifier | PRIDE repository deposit bundle |

## Migrated domains

| domain | writers | legacy write sites at close |
|---|---|---|
| preprocessing | 3 migrated, 1 delegated | 0 |
| qc | 14 migrated, 1 delegated | 0 |
| differential_abundance | 11 migrated | 0 |
| enrichment | 1 migrated | 0 |
| spatial_networks | 6 migrated | 0 |
| spatial_validation | 18 migrated | 0 |
| integration | 17 migrated | 0 |
| wgcna | 22 migrated, 2 delegated | 0 |
| publication_source_data | 3 migrated, 6 delegated | 0 |

**Repo-wide legacy write sites: 0.** 95 writers resolve through the
output-layout API, 10 delegate to a library that does. Split-brain writers: 0.

## Historical read-only compatibility policy

Historical trees are never moved or deleted. They remain as compatibility
fallbacks and provenance, recorded in `config/legacy_output_registry.csv`:

- 80 roots `LEGACY_READ_ONLY`, holding 13,411 files, with **0 active writers**
- 9 roots `ACTIVE_NOT_LEGACY` — the canonical domain roots plus
  `results/manuscript` and `results/publication_source_data`

Every migrated reader resolves normalized-first and falls back to the
historical location, so today's runtime reads exactly the files it read before.

## Evidence

- **Construction enumeration.** A generalized AST scanner
  (`R/utilities/wgcna_construction_scan.R`) plus an independent multiline
  literal cross-check. It detects `path_results`, `file.path`, slash-joined
  literals, alias variables, alias functions, string-constant aliases,
  namespaced calls and `Sys.glob`.
- **Resolution equivalence.** 606 concrete WGCNA instances resolved to
  byte-identical paths; 18 guarded absences preserved their behaviour; 4
  root readers were compared by the downstream file set they select.
- **Scientific immutability.** The authoritative WGCNA baseline is the six
  reproducible stage roots: 5,821 files / 928,087,766 bytes, unchanged. The
  three frozen network states and the failed-run carrier keep their hashes.
- **Runtime witness.** `results/reviewer_audit/input_resolution_audit.csv`
  logged 975 WGCNA resolution events with 0 resolving to a normalized path and
  `file_exists` true for all of them.

### Byte baselines and why only one is a gate

| population | files | bytes | hard gate |
|---|---|---|---|
| authoritative scientific immutability (6 stage roots) | 5,821 | 928,087,766 | yes |
| historical non-reconstructible snapshot (reported 21 roots) | 5,926 | 950,668,976 | no |
| observational full tree (+ `reviewer_audit`) | 6,011 | ~1,005,285,122 | no |

The second is a snapshot of a moment whose root set no audit artifact records;
the third includes an append-only runtime log whose bytes change whenever tests
run. Neither is an erroneous measurement — they answer different population
definitions. Details in `audits/phase6g_wgcna_byte_baseline_reconciliation.csv`.

## Defects found and fixed during the migration

1. **50 figures from the failed WGCNA run `microglia_failed_20260720_133211`**
   were already present in the manuscript figure export and its manifest. The
   exporter never validated the dataset-scope segment. Selection is now guarded
   on both discovery routes by an exact path-segment test, the export was
   regenerated and the freeze re-pinned with the superseded values retained.
2. A publication export glob could admit failed-run artifacts.
3. Two audit-tool globs used `*` in the dataset position.
4. A dry-run reported a different output directory than the run used.
5. `wg_has_files()` judged a directory of sub-directories empty, so a
   non-existent normalized directory could shadow a populated historical one.

## Tests and contracts at close

- full testthat suite: 114 files, 0 failures
- cross-repo runtime boundary: 0 dependencies
- publication freeze: 0 FAIL (3 documented accepted WARNs)
- path-contract rewrites: PASS
- ownership: 107 families, 0 multi-writer outputs, 0 weak owners
- generalized legacy-write audit: repo-wide 0

## Detector corrections, kept separate from migration

The generalized output-namespace classifier replaced an earlier detector that
undercounted. `publication_source_data` went from a reported 3 sites to 29 when
the classifier was generalized: **no new writes were introduced**, 26 had simply
been invisible. The final count is 0. Likewise the legacy-output registry was
stale at the start of Phase 6G.9 and its regeneration records work already
completed in Phase 6G.8, not new migration.

## Debt explicitly excluded from Phase 6G

- `biological_claims_table` schema smoke failure (5 undeclared label-provenance
  columns); both inputs predate Phase 6G
- five unused WGCNA path helpers
- one bare WGCNA stage-root consume in `build_biological_claims_table.R`
- three accepted freeze WARNs
- large `results/manuscript/_superseded_*` and `_failed_*` archive trees
- repository relocation and rename
