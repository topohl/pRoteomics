# Phase 6I.4 — reviewer-ledger repair and test isolation

Two defects with two different causes. They are reported separately because
forcing one explanation onto both would have been wrong.

## The canonical ledger

| | |
| --- | --- |
| path | `results/reviewer_audit/input_resolution_audit.csv` |
| tracked | no — gitignored, so **no historical copy exists in git** |
| other copies | none. A repository-wide search for `input_resolution_audit*` returns exactly two paths: this file and `inst/schemas/input_resolution_audit.yml`. There is no frozen snapshot, no derived export, no backup and no test fixture. |
| schema | `inst/schemas/input_resolution_audit.yml`, 14 required columns |
| quoting | every field quoted, comma-delimited, CRLF line endings |
| producers | `append_input_resolution_audit()` in `R/paths.R`, reached via `record_input_resolution()` |
| readers | one production consumer: `analysis/publication_source_data/build_biological_claims_table.R` |

### Before and after

| | before | after |
| --- | --- | --- |
| sha256 | `b272aa3004bac50e9145eb25c6c37b0389a83324a7e91baaeac0a8c47f868fd9` | `71261478aa05815c9613f45397794d9493ef57ce7157bac02b7a08d98da4274d` |

The "after" hash is the state at the moment of repair. This ledger is
append-only and grows on every analysis run, so it is a record of the repair,
not a pin — and nothing pins it. The tests below assert structural properties
(parses completely, no unbalanced record, no synthetic row), never a hash or a
row count.

| bytes | 33,304,278 | 33,304,203 (−75) |
| physical lines | 57,551 | 57,548 |
| complete 14-field records | 57,548 | 57,547 |
| rows a plain `read.csv()` recovers | ~1,600 | **57,547** |
| parse coverage | ~2.8% | **100.0000%** |

## Defect 1 — three stray quote characters

`PARTIAL_WRITE`.

Three lines each carried **one** quote character and **two** fields, in a
fourteen-field schema:

| physical line | byte offset | bytes | raw text |
| --- | --- | --- | --- |
| 1595 | 949,653 | 61 | `ntial_expression_enrichment/03_biological_program_summary.r",` |
| 4504 | 2,588,027 | 4 | `lsx"` |
| 4511 | 2,590,892 | 4 | `at",` |

Each is the **tail of a single quoted field** — the remainder of
`...04_differential_expression_enrichment/03_biological_program_summary.r`,
of a `.xlsx` path, and of a `.dat` path.

The decisive measurement is what surrounds them. Lines 1594, 4503 and 4510 are
**complete 14-field records with even quote counts**, and so are 1596, 4505 and
4512. The three fragments are therefore *insertions* between intact records,
not the second halves of split ones. `write.table(append = TRUE)` issues
several writes per call and this ledger is appended to by concurrently running
analyses; a partial write from one process landed inside the stream of another.

Syntactic effect: one unbalanced quote reopens a quoted field, so every
subsequent line is read as part of that field. Three characters made the
50,000+ lines after them unreadable to any CSV parser. 0.005% of the file cost
97% of it.

**Dating.** The surrounding records carry `file_mtime` values from June 2026
(`2026-06-17`, `2026-06-13`, `2026-06-01`), so the damage is months old and
long predates the work in this phase.

### Why the content was not reconstructed

| evidence source | available? |
| --- | --- |
| git history | no — the ledger is gitignored, 0 commits touch it |
| frozen snapshot | none exists |
| backup / derived export | none exists |
| adjacent rows | intact, and therefore contain no part of the lost records |
| the fragment itself | the tail of one field only |

`BYTE_EXACT_HISTORICAL_RECOVERY` is impossible — measured, not assumed.
Inventing the missing fields would be an `INFERRED_REPAIR` that materially
changes field content, which the brief forbids. Nothing was reconstructed.

### The repair

| | |
| --- | --- |
| action | remove exactly the three fragment lines |
| confidence | `STRUCTURALLY_UNAMBIGUOUS_REPAIR` |
| complete records altered | **0** |
| quarantine | `results/reviewer_audit/input_resolution_audit.quarantined_partial_writes.csv` |

Nothing was destroyed: the three lines are preserved verbatim beside the
ledger, with their physical line numbers, byte offsets, byte lengths, quote and
field counts, the pre-repair ledger hash, and an explicit `recoverable = FALSE`.

Method: line-level on the raw byte stream. A dataframe round-trip would have
rewritten all 57,551 CRLF line endings, which is a far larger mutation than the
defect.

**Proof that nothing else changed.** Re-inserting the three quarantined lines at
their original positions reproduces sha256
`b272aa3004bac50e9145eb25c6c37b0389a83324a7e91baaeac0a8c47f868fd9` — the exact
pre-repair hash. The byte delta is 75, which equals `sum(nchar + 2)` over the
three lines. Unintended byte changes: **0**. Every one of the 57,547 records
recoverable before the repair is byte-identical after it.

## Defect 2 — eighty synthetic fixture rows

`TEST_CONTAMINATION`. A different cause, already remediated in Phase 6I.2.

| | |
| --- | --- |
| rows injected | 80 |
| rows remaining | **0** |
| legitimate rows incorrectly removed | **0** |
| source | `tests/testthat/test-addressability-vocabulary-unification.R` calling `read_csv_optional()` with fixture paths |

The deletion pattern was `^"test-vocabulary",`. That could not have matched a
legitimate row: the ledger holds 141 distinct script ids, and `test-vocabulary`
is declared nowhere in the repository except that one test file. The seven
script ids containing the substring "test" are all real analysis scripts named
`test_*.R` (`test_module_behaviour_coupling.R`, `test_enrichment_module_concordance.R`,
`test_network_group_organization.R`) and none begins with `test-vocabulary`.

## The intended row universe

| | |
| --- | --- |
| physical lines (before) | 57,551 |
| header | 1 |
| complete records | 57,547 |
| damaged partial lines | 3 |
| **accounted for** | **57,551 — exact** |

Schema conformance of the recoverable rows: 14 columns, all required present,
0 unexpected, 0 malformed.

Row-count equality is deliberately not used as the integrity criterion. The
ledger is an append-only event log: 57,547 events over 1,738 distinct
`script|stage|input_name|expected_path|file_mtime` keys, because the same input
is resolved again on every run. Duplicates are expected here, not a defect.

10,442 rows carry an empty `script`, from runs where `PROTEOMICS_SCRIPT_ID` was
unset. Pre-existing, unrelated to either defect, and not changed.

## Test isolation

The contamination was possible because a unit test and the production pipeline
wrote to the same path. Isolation is now by **path routing**, using the
repository's existing environment-override convention rather than a new
mechanism:

```
PROTEOMICS_INPUT_RESOLUTION_AUDIT   route the ledger to a disposable file
PROTEOMICS_PROJECT_ROOT             the existing whole-repository sandbox
```

Routing rather than a dry-run guard, deliberately. The appender must stay
unconditional: `tests/testthat/test-preprocessing-writer-namespace.R` asserts it
has **no** dry-run guard, because that is precisely why
`build_module_score_metadata` is classified `PATH_VERIFIED_STRUCTURALLY_ONLY`.
Routing the destination leaves the writer's behaviour untouched and keeps that
classification honest.

The override is read from the environment, so it is inherited by child
processes. That matters: `test-biological-integration-entrypoints.R` spawns
twelve real analysis scripts with `system2()`, and an in-process fixture would
not have reached them.

`local_input_resolution_audit()` wraps the override for in-process use and
restores the prior value on exit.

## Downstream effect of restored readability

The sole production consumer is
`analysis/publication_source_data/build_biological_claims_table.R:1406`, which
lists the ledger as one row of the reviewer-audit index with
`schema_name = NA`, `manuscript_use_allowed = FALSE` and the note
"Provenance only." It records `exists` and `n_rows`; it never parses the
ledger's contents, and `schema_validation_status()` returns early for a row
with no schema name.

Its recorded `n_rows` is 15,921, from a run when the ledger had that many
lines. It is already stale by tens of thousands of rows because the ledger
grows on every run, and it is not regenerated here.

| consumer | classification |
| --- | --- |
| reviewer-audit index | `NO_CHANGE` — not regenerated; on a future run `n_rows` tracks the ledger's growth, of which this repair is −3 |
| ledger readability | `EXPECTED_VISIBILITY_RESTORATION` |
| any scientific output | `NO_CHANGE` — nothing reads the ledger's rows |

No inferential output, figure, publication source-data file, PRIDE payload or
manuscript artifact depends on this file.

## Frozen historical evidence

None exists for this ledger, so the question of rewriting it does not arise.
The quarantine file is new evidence about the repair, not a rewritten snapshot.
