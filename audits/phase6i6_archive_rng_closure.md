# Phase 6I.6 — archive-only RNG audit and closure

**Resolution: A1 — CLOSED_BY_ARCHIVAL_BOUNDARY.**

No archive byte was changed. No code was seeded. No scientific or outward
artifact moved.

## The archive universe

429 tracked R files, classified by tree:

| universe | files | RNG calls | seed-control calls |
|---|---|---|---|
| ACTIVE | 213 | 9 | 21 |
| ARCHIVE | 36 | **2** | 24 |
| AUDIT | 43 | 1 | 2 |
| TEST | 137 | 48 | 70 |

`archive/` holds 36 R files on disk and 36 tracked — **0 untracked**, so the
parse-aware census covered all of it. `docs/REPOSITORY_ARCHITECTURE.md:23`
declares the tree as *"superseded code kept for provenance"*, explicitly not
*"anything canonical or executable by the registry"*.

Scanning was parse-aware throughout: `getParseData()` over
`SYMBOL_FUNCTION_CALL` tokens, so a function name appearing in a comment, a
string, or as a variable is not counted.

## Active reachability — 0 edges

The question that decides whether any of this is archive code at all.

Measured from the execution side, which cannot silently return zero:
every `source`, `sys.source`, `eval`, `evalq` and `debugSource` call in the
393 non-archive R files was located by walking the parse tree, and the full
deparsed call text was checked for `archive`.

| | |
|---|---|
| non-archive R files scanned | 393 |
| execution calls found | **1,303** |
| of those naming `archive/` | **0** |

The non-zero denominator is the point: a scan that found no `source()` calls at
all would report zero offenders while proving nothing.

Corroborating: `R/utilities/pipeline_registry.R` mentions `archive/` on one
line, a comment. `config/output_layout.yml`, `config/output_namespaces.yml` and
`config/legacy_output_registry.csv` mention it **0** times. No archive script is
declared as a pipeline step anywhere.

An earlier string-side scan appeared to find 23 literals naming `archive/` with
no enclosing call, which would have been a vacuous result had it been trusted —
the parent walk had failed. It was replaced by the execution-side measurement
above rather than reported.

## The two archive RNG consumers

### 1. `archive/01_preprocessing/01_impute.r:84` — `rnorm`

| | |
|---|---|
| role | **AR1 — HISTORICAL_EXECUTED_PROVENANCE** |
| seed state | **SEEDED_EXPLICIT** |
| surviving output | **PUBLICATION_SUPPORTING** |
| exactly regenerable | **yes** |

The perseus-style downshifted-normal imputation that produced the matrices
WGCNA consumes. `IMPUTATION_SEED <- 42L` at :20, and each celltype_layer subset
draws with `subset_seed <- IMPUTATION_SEED + idx - 1L` (:152), passed
explicitly at the call site (:153) — the `seed = NULL` default is never taken.

The seed is recorded in **three independent places**:

1. the archived source itself;
2. `results/logs/01_preprocessing/impute/run_manifest.yml` →
   `parameters: imputation_seed: 42`, with the note *"Deterministic
   per-celltype_layer seeds are assigned after sorting celltype_layer labels."*;
3. `data/processed/01_preprocessing/impute/imputation_qc.csv`, which records
   `base_seed` and **`subset_seed` per subset** (42/42 microglia, 42/43
   neuron_neuropil, …) alongside the output path.

Surviving outputs: six imputed matrices across two dated generations
(20260526, 20260601), plus `imputation_qc.csv` and `sessionInfo.txt`.

This is the one archive stochastic call whose output is publication-facing, and
it is fully reproducible. It is the reason this audit is not A5.

### 2. `archive/.../legacy/02_compareGO_superseded_tail.r:2382` — `slice_sample`

| | |
|---|---|
| role | **AR2 — SUPERSEDED_ANALYSIS** |
| seed state | **UNSEEDED_HISTORICAL** |
| surviving output | **SUPERSEDED** / internal provenance |
| publication dependency | none |

The bootstrap-stability loop inside the compareGO tail that Phase 6H.10
archived byte-exactly. It sat below an unconditional `quit()` in the original
file, so it could not execute even before archival.

Already adjudicated in Phase 6H and **not reopened** here, per the brief:
non-inferential, unknown historical seed, orphaned artifact, removed from the
outward PRIDE payload (1,310 → 1,309) and retained internally as provenance.
Nothing found in this audit contradicts that.

### Archive seed-control calls

24 `set.seed` calls across 9 archived files (PCA plot generations, EWCE legacy,
WGCNA trait preservation, control-strata figures). They are recorded for
completeness; none is paired with an RNG consumer that this census found, which
simply means those files' stochastic work happens inside library calls rather
than through the base generators scanned for.

## Why no code was changed

Per the preservation principle: historical unseeded execution is itself a
provenance fact. Adding a seed to `02_compareGO_superseded_tail.r` would make
the archived source describe a run that never happened. The archive does not
claim to be runnable — the architecture document says the opposite — so there
is no reproducibility promise being broken.

Archive code byte changes: **0**, verified against HEAD.

## Active RNG boundary re-check

Not a new remediation sweep; a boundary check that the archive audit did not
move anything into the active tree.

9 active RNG calls across 6 files, **0** with no seed control:

| file | line | call | seeding |
|---|---|---|---|
| `analysis/spatial_validation/validate_spatial_foundations.R` | 161, 165 | `runif` | `set.seed` at :158, directly above |
| `analysis/wgcna/audit_microglia_module_claims.R` | 924 | `sample` | `set.seed` at :919 |
| `R/networks/animal_spatial_network_utils.R` | 277, 300, 301 | `sample` | `set.seed` at :268 / :292 |
| `analysis/spatial_networks/test_differential_network_stability.R` | 177 | `sample` | inside a bootstrap helper; `set.seed(params$seed)` at :426 in the caller |
| `analysis/spatial_networks/test_network_stability.R` | 137 | `sample` | inside `bootstrap_iteration()`; `set.seed(params$seed)` at :249 |
| `R/networks/wgcna_network_position_utils.R` | 131 | `sample` | inside `wnp_draw_null_rows()`; `set.seed(seed)` at :175 in the caller |

The last three are `SEEDED_INHERITED_KNOWN`: the draw sits in a function
definition and the caller seeds with an explicit parameterised value, so the
seed is above the call in execution order even though it is below it in the
file. Consistent with the Phase 6H.8 finding; no contradictory evidence.

Active PRIMARY_INFERENCE unseeded RNG: **0**.

The single AUDIT-tree RNG call, `audits/part29/atlas_lineage.R:402` `sample`, is
seeded at :396.

## An adjacent observation, not acted on

The imputed matrices are not declared in `docs/file_contracts.tsv`. Unlike the
executed-rank record protected in Phase 6I.3, they are genuinely regenerable —
seed 42 is recorded three ways and the producer is byte-preserved — and they
have active readers, so their loss would fail a pipeline loudly rather than
silently. It is a contract-registry gap rather than an RNG or retention defect,
and it is recorded here rather than folded into an RNG phase.

## Closure

`ARCHIVE_ONLY_RNG_DEBT` is closed as **A1**. Every stochastic call under
`archive/` is correctly preserved historical code; no active execution path
reaches it; and no public reproducibility claim depends on rerunning any of it.
The one publication-supporting archive draw is exactly reproducible from a
recorded seed.

`tests/testthat/test-archive-rng-boundary.R` guards the two properties the
closure rests on: that nothing executes archived code, and that the archive
stochastic inventory is still exactly these two calls.
