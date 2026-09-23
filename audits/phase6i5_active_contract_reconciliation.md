# Phase 6I.5 — active contract reconciliation

Two contracts had drifted away from what the repository actually produces. Both
were adjudicated by establishing which side was stale, not by assuming.

No claim, value, direction, gate or label changed. Both repairs are
declaration-only.

## 1. biological_claims_table schema — S1 CONTRACT_STALE

### The failure

```
FAIL smoke_test_active_script_contracts
biological claims table schema: Schema 'biological_claims_table' validation
failed. Unexpected column(s): canonical_display_label, canonical_label_source,
Stage01_ModuleLabel_Final, Stage07_label, Stage06_label
```

Schema-only. The artifact exists at
`results/tables/biological_claims_table.csv`, 1,802 rows × 75 columns, and is
readable; nothing was missing or corrupt. `inst/schemas/biological_claims_table.yml`
declared 49 `required_columns` and 70 entries under `columns:`; there were
**0 declared-but-absent** columns and **5 undeclared** ones.

`validate_table_schema()` computes `allowed <- names(col_defs)` — the `columns:`
block — so a column absent from it fails strict validation regardless of
`required_columns`.

### Which side moved

| date | commit | event |
|---|---|---|
| 2026-07-22 | `170de59` | schema last touched (49 required / 70 declared) |
| 2026-09-11 14:51 | `8be51a6` | producer gains the five columns — "Add WGCNA final label approval table and one canonical display-label resolver". No schema file changed. |
| 2026-09-11 16:11 | — | the artifact is written with 75 columns |

Fifty-one days apart. The producer moved; the contract stood still. An
incomplete migration, and the output is the correct side.

### The five columns

Attached by `attach_canonical_wgcna_display_label()` in
`R/statistics/wgcna_label_activation_utils.R`, whose comment at :406-408 states
the intent directly:

> No historical label is discarded: Stage01_ModuleLabel_Final, Stage06_label and
> Stage07_label are retained alongside so the submission artifact still records
> what each naming stage said.

`tests/testthat/test-wgcna-label-activation.R` already asserts this — the test
is literally named "provenance labels are retained, never discarded".

| column | classification | distinct | NA | note |
|---|---|---|---|---|
| `canonical_display_label` | LABEL | 41 | 0 | one display label per dataset + WGCNA entity |
| `canonical_label_source` | LABEL_PROVENANCE | 3 | 0 | **a live distinction**: 245 `active_reviewed_registry` (microglia, human-adjudicated) vs 834 `stage07_canonical` (automated), 723 blank for non-WGCNA rows |
| `Stage07_label` | LABEL_PROVENANCE | 41 | 0 | superseded stage-07 naming |
| `Stage01_ModuleLabel_Final` | LABEL_PROVENANCE | 1 | **1802** | all-NA in this run |
| `Stage06_label` | LABEL_PROVENANCE | 1 | **1802** | all-NA in this run |

None is a SCIENTIFIC_VALUE.

Two of the five carry no values today, which invites the question of whether
they are inert. They are not: `pick()` returns NA when the resolved label table
has no matching Stage-01/Stage-06 entry, so emptiness records that those naming
stages had nothing to say for this run — which is itself the provenance. The
producer writes them unconditionally and a test requires them, so removing them
would break a passing contract to make a schema shorter.

`canonical_label_source` settles the question of whether this is meaningful
provenance: it distinguishes a reviewer-signed label from an automated one, and
both values actually occur.

### The root cause: the producer validated a frame it does not write

Declaring the columns makes the smoke green, but it would not have stopped this
happening again. The producer had its own copy of the same check, in the wrong
place:

| line | statement |
|---|---|
| 2212 | `validate_table_schema(claims, "biological_claims_table", strict = TRUE)` — on a **70-column** frame |
| 2228 | `claims <- attach_canonical_wgcna_display_label(claims)` — five columns arrive |
| 2234 | `readr::write_csv(claims, csv_out)` — a **75-column** artifact |

So `build_biological_claims_table.R` passed its own contract while writing a
table that failed the identical check in
`smoke_test_active_script_contracts.R`. The producer was structurally incapable
of noticing the drift it was creating, which is why 51 days went by.

The validation is moved to immediately before the write (attach 2227 →
validate 2238 → write 2243). The producer can no longer emit a table its own
contract rejects. Without this, the next column the label resolver adds would
reproduce the failure exactly.

A second `validate_table_schema()` call at :1459 was left alone: it assigns a
`schema_ok` boolean inside a violations summary rather than gating anything,
and it passes either way now that the five columns are optional.

Note that the smoke logic was never wrong. It read the artifact from disk and
truthfully reported five undeclared columns. This is S1, not S4.

### The repair

The five are added to the `columns:` block of
`inst/schemas/biological_claims_table.yml` as `{type: character}`, and
**deliberately not** to `required_columns`:
`attach_canonical_wgcna_display_label()` returns early when a claims table has
no WGCNA entity rows, so requiring them would be a false statement about a
legitimate table.

One authoritative source. No second vocabulary, no hard-coded column list in a
test.

### Claim integrity

The artifact was not touched. Before and after:

| | |
|---|---|
| sha256 | `0b407dc79944b2ae5c45bacfcaad39122cadd343866d02010e1505f9319756b8` |
| rows × columns | 1,802 × 75 |
| mtime | 2026-09-11 16:11:30 |

Claim additions 0, removals 0, numerical changes 0, direction changes 0,
inferential-status changes 0, label changes 0, provenance-field changes 0 —
byte-identity makes each of these trivially true.

## 2. clusterProfiler file contract — F2 WRONG_PATH

### What was declared

`docs/file_contracts.tsv` declared two objects under

```
results/tables/04_differential_expression_enrichment/clusterProfiler/<dataset>/**/protein_group_audits/
```

That directory **exists** — with three empty dataset subdirectories, **0 files**
and **0** `protein_group_audits` directories beneath it.

### Where the artifacts are

| root | exists | files | `protein_group_audits` dirs |
|---|---|---|---|
| declared in the contract | yes | **0** | **0** |
| current canonical namespace `results/differential_abundance/clusterProfiler/` | **no** | – | – |
| `data/processed/04_differential_expression_enrichment/clusterProfiler/` | yes | 1,210 | **54** |

The producer writes through `CANONICAL_PATHS$models`. The executed run used the
older `data/processed/` namespace; the current canonical namespace has never
been written by any run. The contract named a third location that was never
populated.

### Adjudication

**F2 — WRONG_PATH.** Not F1/F5: the family is not retired — it is actively
consumed by `compare_go_enrichment.R` through the canonical manifest reader and
is protected by the Phase 6I.3 retention declaration. Not F3: the producer does
generate it, 54 directories of it. Not F4: zero is not a legitimate cardinality
here, because 1,210 files exist.

Two rows carried the identical wrong root —
`clusterProfiler_protein_group_audits` and `enrichment_term_gene_provenance`.
Both are corrected. Fixing one and knowingly leaving the other false would not
have been a repair.

Nothing was recreated to satisfy a declaration, and no enrichment was rerun.

### Cardinality

The registry has no cardinality vocabulary: its columns are `object_id`, `path`,
`created_by`, `consumed_by`, `required_columns` (CSV columns, not file counts),
`description`, `version`. None of that needed inventing here, because zero was
never the correct answer for this row — the artifacts exist.

### Why nothing caught it

`tests/smoke_test_file_contracts.R` validates the registry's own columns,
rejects duplicate `object_id`s and rejects empty paths. It never checks that a
declared path resolves to anything. A row could name any string and pass.

`tests/testthat/test-active-contract-truthfulness.R` now closes that, in two
strengths. Every one of the 52 rows must resolve to a real location — and
separately, for the two rows adjudicated here, the declared tree must actually
contain `protein_group_audits` directories. The second check is necessary
because the first would **not** have caught this defect: the declared tree
existed, it was merely empty.

The test asserts presence, not a file count: the number of comparisons is
legitimately variable.

## Adjacent observations, not changed

`spatial_edge_validation` and `spatial_network_objects` declare paths under
`results/spatial_networks/`, which exists but holds 0 files. That is a different
situation — those analyses have not been run in this checkout — and is out of
scope for this phase. Both rows resolve to a real root, so they pass the new
test.
