# Phase 6H.8 — statistical RNG reproducibility audit

Audit and adjudication only. No RNG was seeded, no analysis was rerun, no
artifact was modified.

Phases 6H.6 and 6H.7 closed **rendering** RNG, where a seed can move only
geometry. This phase examines **statistical** RNG, where a seed can move a
sampled observation, a permutation p-value, an FDR status, and therefore a
reported conclusion. The two are kept conceptually and constant-wise separate.

Evidence tables:
- `phase6h_statistical_rng_audit.csv` — 59 rows, the adjudicated inventory
- `phase6h_statistical_rng_sweep_raw.csv` — 334 rows, the raw mechanical sweep
  including `archive/` and `tests/`

## The two findings that changed the conclusion

An adversarial trace overturned two things this audit initially concluded. Both
were verified independently before being accepted, and both are recorded here
because each inverts part of the adjudication.

**1. The `slice_sample()` call cannot execute.** It sits in the script's
disabled legacy tail. `compare_go_enrichment.R:569` is a bare, column-0,
unconditional `quit(status = 0, save = "no")`, followed at :571 by the marker
`# LEGACY_COMPAREGO_TAIL_DISABLED_BY_CANONICAL_EXIT`. Proven from the parse tree
rather than from indentation: top-level expression 103 of 439 is exactly that
`quit()`, and **336 of 439 top-level expressions — lines 572 to 4321 — are
unreachable**. The bootstrap block is top-level expression 311, at lines
2943-2975. The marker was introduced on 2026-09-17 in commit `5919819`.

The script itself remains registered `active_required`; it is the *tail* that is
dead, not the file. But the consequence for this audit is decisive: **there are
no future runs of this call to seed.**

**2. A live copy of its output is staged for PRIDE deposition.** The initial
sweep searched `results/` only and therefore missed `pride_submission/` at the
repository root. There is a copy at

```
pride_submission/supplementary_tables/
  results_tables_04_differential_expression_enrichment_compareGO_neuron_neuropil_
  BP_phenotype_within_unit_08_Bootstrap_Stability_Summary.xlsx
```

5,141 bytes, sha256 `024d3671f2cd6026…`, mtime 2026-05-27. It is row 1564 of
`pride_submission/manifests/pride_file_manifest.tsv` (regenerated 2026-09-21)
with `export_category = pride_staging` and **`intended_for_PRIDE = TRUE`**.

That flag is not a curation decision. `R/utilities/export_helpers.R:881` assigns
`pride_staging` to anything whose path matches `/pride_submission/`, and :933
marks every `pride_staging` file `intended_for_PRIDE = TRUE`. The manifest
builder rescans the directory (`build_pride_manifest.R:36-40`), so the file is
re-endorsed on every build purely for being on disk.

So the artifact is **not** superseded-only, as this audit first recorded. A copy
is staged for public deposition — produced by code that can no longer run, from
an RNG state that was never recorded.

## A. Active RNG inventory

Discovery is `getParseData()` `SYMBOL_FUNCTION_CALL` based, over `analysis/`,
`R/` and `tools/`. Prose and comments that merely mention `sample()` or
`set.seed()` are not calls; a line grep counts them, and that mattered here.

| | calls |
|---|---|
| RNG-related calls in the active tree (excluding seed-control) | 59 |
| — rendering | 38 |
| — **statistical** | **21** |
| of the 21, empirically deterministic despite a stochastic-sounding name | 1 |
| genuine statistical RNG calls | 20 |
| — with a deterministic seed source | **19** |
| — without any seed | **1** (unreachable — see above) |
| **executable statistical RNG calls** | **19** |
| **executable statistical RNG calls lacking a seed** | **0** |
| active `set.seed` / seed-control calls | 21 across 13 files |
| archive-only statistical calls (reported, never modified) | 84 across 16 files |
| archive-only seed-control calls | 24 |
| test-universe statistical calls | 37 across 14 files |

Inferential role of the 20: 4 `PRIMARY_INFERENCE`, 6 `SENSITIVITY_ANALYSIS`,
3 `DESCRIPTIVE_DIAGNOSTIC`, 6 `VISUALIZATION_ONLY`, 1 `UNUSED`.

**Every `PRIMARY_INFERENCE` call is seeded.** The only unseeded call is the
`UNUSED`, unreachable one.

### Hard gate

Unexplained active statistical RNG calls = **0**. All 20 carry a recorded
purpose, seed provenance, downstream artifact and inferential role in the audit
CSV.

### Two corrections the mechanical sweep required

A naive scan produced eight "ungoverned" calls. Six were false positives, for
two reasons that will recur:

1. **A call inside a function definition.** The definition sits above the
   `set.seed` in the main block, so a line-order test reports a governed call as
   ungoverned. Resolved by walking the parse tree for the enclosing function and
   checking whether its *call sites* are seeded. This cleared
   `test_network_stability.R:137`,
   `test_differential_network_stability.R:177`,
   `wgcna_network_position_utils.R:131` and the three
   `animal_spatial_network_utils.R` draws.
2. **A seed applied by a cross-file wrapper.** The three `fgsea` calls in
   `test_microglia_targeted_signatures.R` have no `seed` argument and their file
   contains no `set.seed`, yet all three are wrapped in
   `run_with_stable_gsea_rng(..., gsea_seed = ...)` from
   `R/enrichment/clusterprofiler_reproducibility.R`. Only execution-path
   tracing finds this.

One reclassification: `WGCNA::pickSoftThreshold`
(`build_wgcna_modules.R:2230`) was flagged by a library-name heuristic and
**cleared empirically** — identical output across incoming seeds, and its body
contains no `sample`, `runif`, `rnorm` or `set.seed`. It is deterministic.

### The repository already has a statistical seed convention

`run_with_stable_gsea_rng()` is a well-formed deterministic resampling scope:

```r
withr::with_preserve_seed({
  previous_kind <- RNGkind()
  tryCatch({
    RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
    set.seed(gsea_seed)
    do.call(gsea_fun, dots)
  }, finally = { do.call(RNGkind, as.list(previous_kind)) })
})
```

Its seed is derived from a *semantic identity* rather than a bare constant:
`targeted_enrichment_reproducibility()` builds
`script=…|dataset=…|comparison=…` plus `method=…`, and marks a method stochastic
only when it is (`method %in% c("fgsea", "clusterProfiler_GSEA")`).
`clusterprofiler_fgsea_control_args()` deliberately sets clusterProfiler's own
`seed = FALSE` so its internal flag cannot override that scope, with a comment
saying why. L'Ecuyer-CMRG is the parallel-safe kind, so this is already the
brief's section-12 **option B** pattern.

`wnp_draw_null_rows()`'s caller is a second example: `set.seed(seed)` with an
`on.exit` restore of `.Random.seed`, so its permutation null is deterministic
*and* leaves the global stream untouched.

Verified empirically that the wrapper is load-bearing: `fgsea` genuinely depends
on the RNG — p-values moved by ~0.05 and NES was **not** identical across
incoming RNG states.

## B. `compare_go_enrichment.R:2951` deep trace

The enclosing block is the script's own `--- 12. BOOTSTRAP ENRICHMENT
STABILITY ---` (lines 2943-2975, top-level expression 311). The code's own
terminology is "bootstrap" and "Resample with replacement", and the
implementation does match bootstrap semantics, so that is the term used here.

| | |
|---|---|
| reachability | **unreachable** — 2,374 lines after the unconditional `quit()` at :569 |
| enclosing block | `# --- 12. BOOTSTRAP ENRICHMENT STABILITY ---`, top-level (not in a function) |
| input table | `combined_df` — `bind_rows(enrichment_list)`, built at :1052-1235 from `manifest_filtered$output_table` |
| **one row** | one **(Comparison × GO term) clusterProfiler GSEA result**, carrying pre-computed `ID`, `Description`, `NES`, `p.adjust`, `setSize`, `core_enrichment` |
| sampling unit | **enrichment-result rows** — not genes, not proteins, not animals, not comparisons |
| grouping | `group_by(Comparison)` — stratified within comparison; comparisons themselves are never resampled |
| sample size | `prop = 1`, each stratum redrawn to its own size `n_c` |
| replacement | **yes** |
| repeats | `n_bootstrap <- 100`, a plain sequential `for (b in 1:100)` |
| statistic | `sum(top_terms %in% unique(sampled_df$Description[sampled_df$p.adjust < 0.05]))` |
| summary written | `Mean_Recovery_Rate`, `SD_Recovery_Rate`, `Min_Recovery`, `Max_Recovery`, `Total_TopTerms` |
| output file | `08_Bootstrap_Stability_Summary.xlsx` in `subdirs$tables` |

`top_terms` are significant by construction: `config/compareGO_config.yml` sets
`significant_only: true`, so `candidate_pool` is filtered to `p.adjust < 0.05`
(:1317) before the redundancy-filtered selection that defines `top_terms`
(:1392).

That is config-dependent, and it matters. The bootstrap uses raw `top_terms`,
not the figure-restricted `figure_top_terms <- intersect(top_terms,
significant_term_descriptions)` (:1397). With `significant_only: true` that
intersection is a no-op, so the two agree. Were the flag flipped to `false`,
`candidate_pool` would become the unfiltered `combined_df_best` (:1319) and
`top_terms` could contain terms with no `p.adjust < 0.05` row anywhere — terms
that can **never** be "recovered", deflating `Mean_Recovery_Rate` by a purely
definitional amount unrelated to any resample. A latent defect in already-dead
code, recorded rather than fixed.

### What it does not do

**The bootstrap recomputes no enrichment statistic, no p-value and no
multiple-testing correction.** It re-draws rows that already carry `p.adjust`
and filters on that carried column. Proven on a faithful isolated fixture: the
number of `p.adjust` values appearing in a resample that are **absent from the
source table is 0**, at every seed tested. Nothing is refit, so no new
inferential quantity can come into existence.

### What `Mean_Recovery_Rate` therefore measures

Only the probability that a pre-existing significant row is drawn at least once
in a same-size resample with replacement. That has a closed form: for a top term
`t` with `m_tc` significant rows in comparison stratum `c` of size `n_c`,

```
P(t recovered) = 1 - prod_c (1 - m_tc / n_c)^(n_c)     ->   1 - exp(-M_t)
```

where `M_t` is `t`'s total significant-row count. For a term significant in
exactly one comparison this is `1 - e^-1 = 0.632` — the textbook bootstrap
inclusion constant.

The surviving values are consistent with precisely that:

| artifact | Mean_Recovery_Rate | Total_TopTerms |
|---|---|---|
| `_superseded_20260622/…/08_Bootstrap_Stability_Summary.xlsx` | 0.6535714 | 14 |
| `…compareGO_microglia_BP_phenotype_within_unit_…` | 0.7208 | 25 |
| `…compareGO_neuron_soma_BP_phenotype_within_unit_…` | 0.650625 | 16 |
| `…compareGO_neuron_neuropil_BP_phenotype_within_unit_…` | 0.8846032 | 63 |
| **`pride_submission/…neuron_neuropil…` (staged)** | **0.8738182** | **55** |

Three of the four superseded values sit just above 0.632, i.e. dominated by
singly-significant terms; the neuropil values are higher because more of its top
terms are significant in several comparisons.

**A caution about the two neuropil copies.** They differ — 0.8738182 / 55 top
terms (PRIDE-staged, 2026-05-27) versus 0.8846032 / 63 top terms (superseded,
2026-06-13). But `Total_TopTerms` also differs, so the inputs changed between
those runs. The pair therefore proves the **file** is not reproducible; it is
**not** clean evidence of RNG-induced variation, because an input change is
confounded with it. The clean RNG evidence is the isolated fixture in section D,
where inputs are held fixed.

## C. Seed provenance

| question | answer |
|---|---|
| explicit `seed` argument on the call? | no |
| `set.seed` upstream in the same file? | **no — the 4,321-line script contains none** |
| `withr::with_seed` / wrapper anywhere on the execution path? | no |
| `future.seed` / worker RNG? | n/a — no parallelism in this script |
| current state | **inherited global RNG only** (and, in practice, unreachable) |

`slice_sample` at :2951 is the script's **only** RNG call of any kind.

### Historical reconstruction

**`HISTORICAL_SEED_UNKNOWN`.**

The script writes its own `--- 13. SUPPLEMENTARY: PARAMETER & REPRODUCIBILITY
LOG ---` (`09_Reproducibility_Log.xlsx`) from `analysis_params` at :905-921.
That list records 16 fields — script, version, dataset, timestamp, R version,
platform, ontology, `ensemble_profiling`, condition, base path,
`significant_only`, `target_n_terms`, `redundancy_threshold`, `min_set_size`,
`n_comparisons`, `n_total_proteins` — plus a package-version sheet, and **no
seed and no RNG state**. A log named "reproducibility log" does not record what
would be needed to reproduce the one stochastic number in the script. Confirmed
not only from source but by unzipping all four surviving superseded workbooks
and reading `xl/worksheets/sheet1.xml` and `sheet2.xml` directly: parameters and
package versions only, no seed row, no `RNGkind` row, no `.Random.seed` row.

With no seed set, none recorded, no `.Random.seed` dump, and the producing code
now unreachable, the surviving artifacts cannot be linked to a known seed. Per
the brief, no seed was brute-forced.

## D. Sensitivity to RNG

Measured on a mathematically faithful isolated reimplementation (synthetic data,
same operations, inputs held fixed, 100 iterations, 8 incoming RNG states):

| question | answer |
|---|---|
| does the incoming RNG state change the draw? | **yes** |
| does it change the reported statistic? | **yes** — `Mean_Recovery_Rate` spanned 0.6748–0.6928 (spread 0.018) |
| does it change any p-value or FDR status? | **no** — none is computed in this block |
| does it change a selected term set? | **no** — `top_terms` and `figure_top_terms` are computed before and independently of it |
| does it change a plotted point? | **no** — no figure consumes the statistic |
| does it change a textual claim? | **no** — see section G |

Strongest observed consequence: **the third decimal place of a supplementary
diagnostic number moves.** At `n_bootstrap = 100` the Monte Carlo standard error
of `Mean_Recovery_Rate` is ≈0.009, while the values are recorded to 7 decimal
places — so most of the digits as written are Monte Carlo noise rather than
signal. Raising the estimate to 10,000 iterations reduced that error to ≈0.0009
and converged on the closed-form value.

Closed-form agreement: predicted 0.68317, Monte Carlo mean over 8 seeds
0.68325 — a difference of 8e-05. **The seed moves the Monte Carlo error of the
estimate, never the estimand.**

## E. Inference classification

**`UNUSED`.**

The call is unreachable, and no code reads any copy of its output. Had it been
reachable it would have been `DESCRIPTIVE_DIAGNOSTIC`: it refits nothing and
produces no p-value, so not `PRIMARY_INFERENCE`; it does not perturb an analysis
and re-derive a conclusion, so not `SENSITIVITY_ANALYSIS`; no figure reads it,
so not `VISUALIZATION_ONLY`.

Per the brief's section 10 this classification bounds how aggressive remediation
may be — and `UNUSED` is the least permissive case for code changes.

## F. Parallel and repeated RNG

### The subject call

`compare_go_enrichment.R` uses **no** parallel framework: a case-insensitive
search across all 4,321 lines for `future`, `furrr`, `BiocParallel`,
`mclapply`, `foreach`, `parallel`, `doParallel`, `multisession`, `makeCluster`,
`parLapply` and `bplapply` returns nothing. There are no workers, so no
worker-RNG inheritance question arises. The bootstrap is a plain sequential
`for (b in 1:100)` drawing from one stream, so its 100 replicates would be
independent draws in sequence. There is **no repeated-identical-seed defect** —
the opposite failure mode, where every replicate is identical.

### Parallel fan-out elsewhere in active code

There are two, and both are governed. This section initially covered only the
absence of parallelism in the subject script; the sweep's call-level view missed
these because `future_lapply` is not itself an RNG draw.

**1. `analysis/enrichment/run_ewce_celltype_enrichment.R:1273-1276**

```r
target_results <- future.apply::future_lapply(
  target_grid$TargetRun, run_ewce_target,
  future.seed = analysis_params$seed
)
```

An **explicit numeric** `future.seed` makes `future.apply` derive per-element
L'Ecuyer-CMRG substreams deterministically, so both independence and
reproducibility hold, and neither depends on how work is scheduled across
workers. The script additionally sets `set.seed(42)` at :39.

**2. `analysis/differential_abundance/run_clusterprofiler_enrichment.R** — four
`future_lapply(..., future.seed = TRUE)` sites (:532/547, :550/564, :2130/2142,
:2145/2153).

`future.seed = TRUE` alone would give statistically sound independent streams
but would leave results dependent on the ambient RNG state. The determinism of
the GSEA results does **not** rest on it. Inside each worker,
`analyze_comparison()` reaches (:1273-1292):

```r
seed <- derive_clusterprofiler_gsea_seed(
  analysis_params$gsea_seed_base, comparison_name, analysis_type)
run_seeded_clusterprofiler_gsea(..., gsea_seed = seed, ...)
```

so every GSEA is pinned by a seed derived from
`(gsea_seed_base, comparison, analysis_type)` — **independent of worker
scheduling and of `future.seed`**. That is the correct design for a fan-out:
per-unit deterministic seeds rather than reliance on stream assignment.
`gsea_seed_base` is `20260824L` (:283, and from config at :965).

`R/enrichment/clusterprofiler_reproducibility.R:20-35` documents the measured
limit of this contract, which is worth recording because it is an honest one:
clusterProfiler routes GSEA through `DOSE:::GSEA_fgsea()`, which hardcodes
`nproc = 0`, so the backend falls back to ambient `bpparam()`. Backend choice
was *measured* not to be the source of drift — default SnowParam (30 workers,
unseeded), explicit SnowParam with pinned `RNGseed` at 8 and 4 workers, and a
run preceded by another GSEA call in the same session were all field-identical —
so pinning `BPPARAM` is deliberately not done. Residual cross-context tolerance:
`enrichmentScore` agrees to ~1.6e-15, `setSize`/rank/`leading_edge`/
`core_enrichment` exactly, and that last-bit difference propagates to ≤~2.4e-05
in NES and ≤~2.3e-05 in p-value/FDR.

**No repeated-identical-seed defect was found** at any site: every fan-out
either derives a distinct per-unit seed or uses a sound substream mechanism.

### A note on the shared seed value

`gsea_seed_base = 20260824L` happens to be the same integer as
`NATURE_REPEL_SEED` and `NATURE_JITTER_SEED` from Phases 6H.6 and 6H.7. The
three are **separately declared and independently changeable**, and no code
derives one from another, so the section-13 requirement that statistical and
presentation seeds stay conceptually separate is met. The coincidence of value
is recorded here so that nobody later infers a dependency that does not exist,
or "tidies" them into a single shared constant — which would couple figure
layout to GSEA p-values.

## G. Artifact integrity and the dependency chain

```
slice_sample draw  [UNREACHABLE since 2026-09-17]
  -> sampled_df / overlap_with_top / stability_summary   (transient)
  -> 08_Bootstrap_Stability_Summary.xlsx
       |- 4 copies under results/manuscript/_superseded_20260622/
       '- 1 copy under pride_submission/supplementary_tables/  [intended_for_PRIDE = TRUE]
  -> read by nothing
```

**Code consumers: zero.** Every occurrence of `Mean_Recovery_Rate`,
`SD_Recovery_Rate`, `Num_TopTerms_Recovered`, `stability_results` and
`stability_summary` in the repository is inside the producing block itself
(:2960-2974), plus the Phase 6H.8 test file. No script, figure, export manifest
or freeze entry reads the artifact. The other `stability_summary` /
`bootstrap_*` matches in the tree belong to the **spatial networks** analysis, a
different and **seeded** code path with its own `node_stability_summary()` and
`bootstrap_*_stability_summary.csv` outputs.

**Manuscript consumers: zero.** The sibling repository
(`…/Analysis/Exp9_manuscript`) contains no reference to `Bootstrap_Stability`,
`Mean_Recovery_Rate`, `Recovery_Rate`, `TopTerms_Recovered`, "recovery rate" or
"enrichment stability" — 0 hits for every pattern. **No manuscript claim depends
on this statistic.**

**Deposition consumer: one, by blanket rule.** The PRIDE manifest row described
at the top of this document. Nothing curated it in; nothing will curate it out
while the file remains on disk.

Hashes of every artifact connected to the stochastic call:

| file | bytes | sha256 (first 16) |
|---|---|---|
| `pride_submission/…neuron_neuropil…08_Bootstrap_Stability_Summary.xlsx` | 5,141 | `024d3671f2cd6026` |
| `_superseded…/08_Bootstrap_Stability_Summary.xlsx` | 15,431 | `9fe6b207bfae4129` |
| `_superseded…/…microglia…08_…xlsx` | 5,135 | `93ee406b6715b7b6` |
| `_superseded…/…neuron_soma…08_…xlsx` | 5,136 | `d6d63107bca10c0c` |
| `_superseded…/…neuron_neuropil…08_…xlsx` | 5,143 | `a2b9dcfda9e1ff47` |

Unchanged by this audit:

| artifact | state |
|---|---|
| `figure_export_manifest.csv` | 5,582 rows, `0fd0c9ed…`, 2,560,938 bytes |
| `docs/publication_freeze_manifest.yml` | `b4d37250…`, 35,828 bytes |
| `exports/` | 55 files, 2,545,817 bytes |
| `pride_submission/` | 1,310 files; manifest 1,721 rows |
| scientific files changed | **0** |
| package changed | **0** |
| freeze changed | **0** |

## H. Adversarial assessment

**Could seeding this call change a published numeric result?**
No — and it could not change anything at all, because the call is unreachable.
Seeding it would be a no-op that creates the appearance of remediation. The
published-result question is separately answered no: the statistic is cited in
no manuscript claim, no figure and no other code.

**Could choosing a seed after seeing the current output amount to result
selection?**
Not materially, for two independent reasons. The estimand is closed-form, so
every seed converges on the same quantity and the achievable range is ±0.009 of
Monte Carlo noise around a fixed value — there is nowhere for a motivated choice
to move it. And the call cannot run, so no seed choice has any effect. Were the
tail ever revived, the repository's convention of deriving the seed from a
semantic identity removes the discretion entirely.

**Is the existing output recoverable from a known historical seed?**
No. `HISTORICAL_SEED_UNKNOWN`. Worse than unknown: it is **unregenerable**, since
the producing code has been unreachable since 2026-09-17. The PRIDE-staged copy
cannot be reproduced by running anything in the current repository.

**Is the current output one arbitrary Monte Carlo draw that should instead be
based on many draws?**
It is one arbitrary draw of a 100-replicate estimate, and 100 is too few for the
7 decimal places at which the number is written. But the deeper issue is not the
replicate count: the statistic does not measure what its name implies.
"Bootstrap enrichment stability" suggests robustness of the enrichment result,
whereas the quantity computed is the bootstrap inclusion probability of
already-significant rows, pinned near 0.632 by construction. More replicates
would make a mislabelled quantity more precise. Recorded because section 20 asks
for it; deliberately out of scope for remediation.

**Does the analysis already average over enough resamples that the seed is
immaterial?**
For every *conclusion*, yes — there are no conclusions drawn from it. For
byte-level file reproducibility, no: 100 replicates leave ≈0.009 of Monte Carlo
error, immaterial at 2 decimals and material at the 7 recorded.

**Would a future deterministic seed alter only reproducibility, or also the
estimate?**
Only reproducibility, and only hypothetically, since there is no future run. The
estimand is fixed in closed form by the per-term significant-row multiplicities;
a seed selects which draw of it is written. Confirmed numerically: closed form
0.68317 vs Monte Carlo mean 0.68325.

## I. Recommendation

**R1 — NO CHANGE REQUIRED.**

The brief makes R1 appropriate when the "call is unused/test-only". That is
exactly what it is: unreachable code, downstream of an unconditional `quit()`,
whose output no code, figure, manifest, freeze entry or manuscript claim reads.

Not **R2**, which this audit initially intended to recommend. R2 adds a fixed
seed "for future runs only" — but there are no future runs. Seeding a statement
that cannot execute changes nothing while creating a false record of
remediation, which is worse than leaving it visibly unseeded.

Not **R3**: a master-seed-plus-per-replicate-stream protocol is engineering a
sequential, single-stream, unreachable diagnostic does not need.

Not **R4**, decisively. R4 requires that a frozen scientific output be
unreliable enough to demand reanalysis now. Four independent facts each rule
that out: the block recomputes no inferential quantity (0 new `p.adjust`
values); its estimand is closed-form, so no seed can move it; the statistic is
cited nowhere in the manuscript; and the enrichment results themselves come from
the seeded clusterProfiler path, not from this block.

### The separate, non-RNG defect

R1 is the right answer to the question asked, and it is not the whole picture.
The PRIDE-staged copy is an **irreproducible, unregenerable orphan flagged for
public deposition**:

- it is one unrecorded Monte Carlo draw (`HISTORICAL_SEED_UNKNOWN`);
- its producing code cannot be run to reproduce or supersede it;
- it reports a statistic whose name overstates what it measures;
- it is flagged `intended_for_PRIDE = TRUE` by a path-matching rule, not by a
  curation decision, and is re-endorsed on every manifest rebuild.

That is a payload-curation and dead-code question, not an RNG question, and the
brief's R1–R4 do not cover it. It is recorded here and **not actioned**, in
keeping with sections 19 and 24. It deserves its own decision: whether the file
should be deposited at all, and whether the disabled legacy tail should be
removed rather than left as a source of orphaned outputs.

Two further items are carried as debt rather than actioned:
- `analysis_params` omits RNG state entirely, so the script's "reproducibility
  log" does not capture what its name promises;
- the statistic's name/meaning mismatch in section H.

## Scope notes

- `archive/` holds 84 further statistical RNG calls across 16 files (41 in
  `archive/deprecated/clusterProfiler_newest_Sep16.r` alone) and 24 seed-control
  calls. Reported only; `archive/` is machine-declared non-runnable provenance
  via `pipeline_analysis_script_exclusions()`, and nothing there was modified.
- The test universe holds 37 statistical calls across 14 files with 51
  seed-control calls; 10 of the 14 seed their own fixtures. The remaining four
  (`test-export-canonical-ewce-contract.R`, `test-gsea-wgcna-concordance.R`,
  `test-joint-compartment-qc-publication-figures.R`,
  `test-wgcna-wave1a-semantic-handoff.R`) are classified `TEST_ONLY`; per
  section 3 no seeding is recommended for a fixture whose framework already
  controls its randomness.
- No addressability vocabulary was touched (`integration_utils.R`,
  `evidence_bundle_utils.R`, the vocabulary in `compare_go_enrichment.R`).
- No payload or freeze change; `slice_sample` was neither seeded nor modified.

## Tests

`tests/testthat/test-statistical-rng-audit.R` pins the audit's conclusions so
they degrade loudly:

- the audit table exists and carries the brief's contract columns;
- every active statistical call is explained (the hard gate);
- **only** the adjudicated `compare_go_enrichment.R` call lacks a seed, and every
  `PRIMARY_INFERENCE` call is seeded — a newly added unseeded statistical draw
  fails here;
- no statistical call site uses a rendering seed constant;
- the discovered call set equals the audited call set, so a new draw cannot
  bypass the table;
- the bootstrap block is still unreachable — if the `quit()` at :569 is ever
  removed or made conditional, the reachability assertion fails and the R1
  recommendation must be revisited;
- `pickSoftThreshold` is deterministic (so it is not re-flagged) while `fgsea` is
  not (so its wrapper is not removed as redundant);
- `run_with_stable_gsea_rng` is deterministic, seed-sensitive, and restores both
  `.Random.seed` and `RNGkind`;
- the compareGO statistic tracks its closed form while varying with the seed, and
  creates no `p.adjust` value absent from its source;
- all five surviving artifacts — including the PRIDE-staged copy — plus the
  figure manifest and the freeze are byte-unchanged.
