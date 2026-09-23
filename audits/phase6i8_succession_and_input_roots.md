# Phase 6I.8 — pipeline succession claims and input-root declarations

One problem in two places: **declaration files asserting repository structure
without proving it.**

No scientific output changed, no archive code was reactivated, nothing was
regenerated.

## Part I — the 19 legacy succession claims

`pipeline.yml` has a top-level `legacy:` key listing 19 scripts excluded from
canonical runnable stages. These are **supersession records, not execution
edges** — Phase 6I.6's finding of zero active→archive runtime edges is
unaffected and was re-verified.

### The field is dual-purpose, by design

`replacement:` sometimes names a successor script and sometimes describes what
the excluded thing *is* (`"repository audit utility"`, `"standalone
pre-promotion DA audit"`). A prose value is **correct usage, not a broken
path**, and an earlier pass that flagged 7 entries as "targets that do not
exist" was measuring the wrong thing.

After the two corrections below: **9 real script paths, 10 prose values.**

### Adjudication

Recorded row by row in `audits/phase6i8_succession_adjudication.csv`.

| verdict | n |
| --- | --- |
| `PROSE_NOT_PATH` | 8 |
| `FUNCTIONAL_SUCCESSOR` | 4 |
| `RELATED_BUT_NOT_SUCCESSOR` | 4 |
| `PARTIAL_SUCCESSOR` | 2 |
| `NO_SUCCESSOR` | 1 |
| **total** | **19** |

Corrected in this phase: **2**.

### The two false claims

Both named `analysis/preprocessing/extract_protigy_contrasts.R`. That target
exists, so every existence check passed. The claim is false anyway, and the
reason is structural — the dependency arrow was inverted:

```
01_impute.r        →  impute/*_missing70pct.xlsx
02_excel_convert.r →  excel_convert/*_with_metadata.{xlsx,gct}
                           ↓  (.gct crosses the external ProTigy boundary)
                      protigy_output/<dataset>/*.gct
                           ↓
extract_protigy_contrasts.R  →  contrast CSVs
```

The named "replacement" **consumes the descendants** of both scripts' outputs.
It is three steps downstream of the first and on the far side of the external
ProTigy boundary from the second.

| | `01_impute.r` | `02_excel_convert.r` |
| --- | --- | --- |
| actually produces | `impute/*_missing70pct.xlsx`, `imputation_qc.csv`, `sessionInfo.txt` | `excel_convert/*_with_metadata.{xlsx,gct}` |
| named replacement produces | contrast CSVs under `results/preprocessing/extract_protigy_contrasts/` | same |
| output overlap | **none** | **none** |
| verdict | `NO_SUCCESSOR` | `RELATED_BUT_NOT_SUCCESSOR` |

The verdicts differ because `02_excel_convert.r` at least sits in the same
ProTigy input/output chain, which is why the wrong claim was plausible;
imputation and contrast extraction have no functional relation at all.

**No successor was invented.** Both entries now read *"no direct successor;
… preserved and still consumed as required inputs"*, with a status recording
where the outputs are, that they remain `REQUIRED` active inputs, and which
resolver line proves it.

This matters beyond tidiness: the false claims are the likely origin of the
belief, carried into Phase 6I.7, that these were *actively regenerable*
intermediates. They are not — they are preserved inputs with archived
producers.

### `replacement:` schema — **P1**

The dual-purpose field is acceptable; document and test it.

Typing the field would **not** have caught either defect, because both false
claims named a real, existing script. The failure was semantic, not
structural, so the remedy is an adjudication that records the *relationship*
and a test that checks it — not a schema migration.

`tests/testthat/test-pipeline-succession-claims.R` enforces: every entry
adjudicated with a stated basis; the recorded `value_kind` still matches the
actual value; prose never carries a positive successor verdict; real paths
never sit at `UNKNOWN`; the two corrected entries cannot silently regain
`extract_protigy_contrasts`; and no archived script appears as a registered
pipeline step.

## Part II — input-root declarations

### The structural gap

Before this phase, `docs/file_contracts.tsv` declared **zero** paths under
`data/raw/`, `data/metadata/` or `data/external/`. The registry described
derived results almost exclusively, so the entire upstream input layer was
invisible to it while being read throughout the pipeline.

### Families found

Nine top-level families exist under the three roots; **all nine have active
readers**. One is excluded as not a scientific input:
`data/metadata/README.md`, read by `R/paths.R:61` through
`rprojroot::has_file("README.md")` as the repository-root sentinel.

The remaining **8 were all undeclared** (`I2_UNDECLARED_REQUIRED_INPUT` /
`I3_UNDECLARED_OPTIONAL_INPUT`), now declared in 7 rows:

| object_id | files | readers | origin |
| --- | --- | --- | --- |
| `raw_protein_group_matrix` | 2 | 7 | instrument/search export |
| `sample_metadata_workbook` | 1 | 6 | manually curated |
| `manual_identifier_mapping` | 1 | 4 | manually curated |
| `uniprot_idmapping_reference` | 2 | 11 | external download |
| `external_reference_marker_sets` | 20 | 2 | external reference data |
| `external_behaviour_data` | 3 | 8 | external assay export |
| `external_published_reference_dataset` | 3 | 1 | external publication |

Registry: 54 → 61 rows, 0 duplicate ids, 0 empty paths.

### No producers were invented

Every one of these is externally sourced, instrument-exported or manually
curated. `created_by` says `"… ; NO repository producer"` rather than naming a
script that does not produce it, and a test asserts that wording stays.

### One correction worth recording

`data/external/behavior` looked declared under a naive check, because the
existing `network_behavior_outputs` row contains the substring "behavior" —
but that row's path is under `results/`. A **source** family and a **derived**
family had been conflated by substring matching. The declaration-status check
now reads the `path` column specifically.

The same trap applies to Phase 6I.7's description prose, which mentions
`data/raw/pg_matrix/quicksearch.pg_matrix.tsv` and
`data/metadata/TPE9_sample_metadata_males.xlsx` inside a *description* field.
Those were never declarations, and are not counted as such.

## Imputed matrices and metadata-joined workbooks

Both families have the same shape, established independently of each other:

| | imputed matrices | metadata-joined workbooks |
| --- | --- | --- |
| location | `01_preprocessing/impute/` | `01_preprocessing/excel_convert/` |
| files | 8 (6 matrices + QC + sessionInfo) | 6 |
| active readers | resolver + ~18 scripts | resolver, module-score path |
| required? | **yes** (`dataset_inputs.R:89`, `:190`) | **yes** (`:121`) |
| active producer | **none** | **none** |
| historical producer | `archive/01_preprocessing/01_impute.r` | `archive/01_preprocessing/02_excel_convert.r` |
| classification | canonical preserved input, historical producer | same |
| regenerable *now*? | not by any active step | not by any active step |
| bytes changed | **0** | **0** |

Family fingerprints, unchanged across this phase:
`a6ffc87ab4c6d1584e772772a7004bc02461075b152ada5701077c876c095441` (impute) and
`80bfab0d364a755c9446cb4e7bb0cd9714dfa52b820fb5b7f4e42e30b5d81f90`
(excel_convert).

Neither archived producer was reactivated. Modernising preprocessing would be
a separate project; this phase documents what is true.

## Archive boundary

Unchanged. Correcting supersession semantics does not make archived scripts
runnable, and the test asserts that no archived script appears as a registered
pipeline step. Phase 6I.6 remains **A1 — CLOSED_BY_ARCHIVAL_BOUNDARY**.
