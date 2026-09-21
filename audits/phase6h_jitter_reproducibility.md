# Phase 6H.7 — deterministic jitter layouts in publication figures

Companion to `phase6h_ggrepel_reproducibility.md`. Same class of defect, a
separate mechanism: seeding one does nothing for the other, which is why this
phase exists as its own repair.

Evidence table: `phase6h_jitter_seed_audit.csv` (26 rows, one per call).

## The defect

`ggplot2::position_jitter()`, `position_jitterdodge()` and `geom_jitter()`
displace points by a random draw. All three default to `seed = NA`, which means
they draw from the **global** RNG stream at render time rather than from a
layer-local stream. A figure rendered twice from identical data therefore got
different point coordinates whenever the incoming RNG state differed — and in a
long analysis script it always differs, because every upstream permutation,
bootstrap, sample or model fit advances the stream before the plot is reached.

The consequence is the same one Phase 6H.6 documented for label placement: a
figure whose bytes move when nothing scientific moved cannot be checked against
a freeze, and two runs of the same script cannot be shown to agree.

Measured in this environment (R 4.5.1, ggplot2 4.0.2, ggrepel 0.9.8,
svglite 2.2.2), identical data rendered under two different incoming RNG states:

| condition | SHA-256 equal | geometry equal |
|---|---|---|
| unseeded jitter | FALSE | FALSE |
| seeded jitter | TRUE | TRUE |

The unseeded row is the defect. Note that it is a genuine render difference, not
a metadata artefact: svglite embeds no timestamp, so byte identity is
achievable and the seeded row attains it.

## Scope

A call is in scope if its script writes into a root the manuscript exporter
scans, because that script's next render becomes a publication candidate. This
phase governs future renders, so eligibility — not the contents of the current
payload — is the right test. Two earlier framings were rejected: matching
literal figure filenames missed every dynamically named figure, and matching
the canonical Phase-6G roots matched nothing, because those roots are still
empty and the payload sits in the pre-6G stage-numbered trees.

| | calls |
|---|---|
| stochastic point-placement calls in `analysis/`, `R/`, `tools/` | 26 |
| publication-facing | 21 |
| — seeded after this phase | 21 |
| — **unseeded** | **0** |
| outside export scope, left unchanged | 5 |

By call type: 21 `position_jitter`, 4 `geom_jitter`, 1 `position_jitterdodge`.

Those three are the complete set of stochastic point-placement mechanisms in
this repository, checked rather than assumed. There are no uses of `geom_sina`,
`stat_sina`, `geom_beeswarm`, `geom_quasirandom`, `position_beeswarm`,
`position_quasirandom`, `geom_dotplot` or `position_dodge2`, and no standalone
`jitter()` calls. `ggforce` appears once, as a package name in a dependency list
in `run_clusterprofiler_enrichment.R`, never as a plotting call. So no fourth
mechanism was left unseeded behind the three that were fixed.

The five out-of-scope calls are in
`analysis/integration/quantify_candidate_network_position.R` (2),
`analysis/integration/test_network_behaviour_coupling.R` (1) and
`analysis/spatial_validation/validate_network_workbook.R` (2). They are recorded
rather than repaired, and `tests/testthat/test-jitter-render-determinism.R`
holds them in an explicit allowlist so that a new unseeded jitter anywhere else
fails rather than passing unnoticed.

The audit count is **26**, not the three call sites named in the brief. The
three named lines (`assess_marker_rank_abundance.R:220`,
`score_module_activity.R:1799` and `:3274`) are all present and all repaired,
but a mechanical sweep found 23 more.

### The swept universe, and what it excludes

The sweep covers `analysis/`, `R/` and `tools/`. It does **not** cover
`archive/`, which holds 15 further unseeded jitter calls. That exclusion is
deliberate and worth stating, because the scope criterion above — "writes into a
root the manuscript exporter scans" — is true of at least one archived script
taken literally: `archive/08_integration/01_compartment_fidelity_summary.R:128`
has an unseeded `position_jitter` and writes under
`results/figures/08_biological_interpretation/`, which the exporter does scan,
and three of its outputs are in the current manifest.

It is excluded because the repository machine-declares `archive/` non-runnable:
`pipeline_analysis_script_exclusions()` in `R/utilities/pipeline_registry.R`
lists `archive` among its excluded roots with the comment "Superseded and
exploratory generations live under `archive/` and are provenance, not runnable
stages", and `docs/active_script_io_audit.tsv` marks that specific script
`registry_excluded` with a named replacement
(`analysis/qc/render_compartment_abundance_figures.R`). Its figures are
historical provenance; seeding it would edit frozen code for no live render.

This is also the established convention rather than a choice made here: the
accepted Phase 6H.6 used the identical three-root universe and left 27 unseeded
`geom_text_repel` calls in `archive/` alone.

## The contract

`NATURE_JITTER_SEED <- 20260824L` in `R/utilities/plotting_nature.R`, beside
`NATURE_REPEL_SEED`. Defined once; not derived from a path, a clock or the
working directory.

Deliberately a **distinct constant** rather than a shared render seed. Jitter
and repel are independent mechanisms, and a single shared constant would mean
that changing one figure's label layout silently moved the point placement of
every jittered panel in the repository. The shared numeric value is incidental —
only stability matters, and nothing scientific depends on which integer it is.

### geom_jitter needs an indirection

`geom_jitter()` has no `seed` argument:

```
mapping, data, stat, position, ..., width, height, na.rm, show.legend, inherit.aes
```

A `seed =` passed to it lands in `...`, raises only
`Ignoring unknown parameters: seed`, and is **not honoured** — measured directly:
`geom_jitter(width = .15, seed = 1L)` still renders differently under two
different incoming RNG states. In a script that already emits warnings, that
looks like a fix and behaves like the defect.

Every `geom_jitter` in scope therefore routes through an explicit
`position = ggplot2::position_jitter(..., seed =)`, which also means its
`width`/`height` move into that call — not as a matter of style but because
ggplot2 **errors** if both are given:

```
Error: Both `position` and `width`/`height` were supplied.
i Choose a single approach to alter the position.
```

The restructuring is otherwise inert. Measured: a bare
`geom_jitter(width = .15)` and an explicit
`geom_point(position = position_jitter(width = .15))` rendered from the same
incoming RNG state are **byte-identical**, because `geom_jitter` forwards an
unsupplied `height` as `NULL`, which is exactly `position_jitter`'s own default.
Both therefore fall back to the same 40%-of-resolution height jitter, so the one
restructured call in this phase did not silently acquire or lose vertical
displacement.

The tests assert the missing formal, the warning, and the non-honouring, so the
reason for the indirection is recorded where someone would otherwise undo it.

## Repaired call sites

14 seed insertions across 10 scripts, covering **16 audit rows**. The two counts
differ because the EWCE stratum plot and the joint-compartment plot are each a
single `geom_jitter(position = position_jitter(..., seed =))` call that the audit
records as two rows — the outer `geom_jitter` (seed delegated) and the inner
`position_jitter` (seed explicit). So 12 single-row calls + 2 double-row calls =
16 rows, which is exactly the number of `NATURE_JITTER_SEED` rows in
`phase6h_jitter_seed_audit.csv` and the sum of the table below. It also
reconciles with the scope table above: 21 publication-facing less the 5
pre-existing `seed = 1` calls left untouched = 16.

| script | layers |
|---|---|
| `analysis/enrichment/run_ewce_celltype_enrichment.R` | 1 call (2 rows) |
| `analysis/qc/assess_joint_compartment_quality.R` | 1 call (2 rows) |
| `analysis/qc/assess_marker_rank_abundance.R` | 1 |
| `analysis/qc/assess_replicate_consistency.R` | 1 |
| `analysis/qc/assess_sample_quality.R` | 2 (one `position_jitterdodge`) |
| `analysis/qc/export_marker_traits.R` | 1 |
| `analysis/qc/summarize_marker_detectability.R` | 3 |
| `analysis/qc/summarize_missingness.R` | 1 |
| `analysis/wgcna/render_module_figures.R` | 1 |
| `analysis/wgcna/score_module_activity.R` | 2 |

Seven of these scripts did not previously source `plotting_nature.R`; the line
was added directly after their `source(paths_file)`.

That addition is the one part of this change that could plausibly alter a
figure — a name collision would silently swap a theme, palette or save helper
underneath the script and so change figure content, which brief constraint 1
forbids. It was checked mechanically rather than assumed:

| check | result |
|---|---|
| names `plotting_nature.R` introduces | 21 |
| top-level expressions that are **not** assignments (side effects) | **0** |
| collisions with `paths.R`, the only thing sourced earlier | **0** |
| names shadowed back by a later source in any of the 7 scripts | **0** |
| scripts where the paths bootstrap precedes the new line, which precedes first use | 10 of 10 |

(The bootstrap line is `source(paths_file)` in nine of the ten;
`assess_sample_quality.R` reaches `repo_path()` through an earlier
`source(early_paths_file)` and sources `paths_file` again later, so the new line
sits after the early one.)

The ordering matters and is deliberate: because the new line runs before every
other `source()` and before the script's own definitions, anything defined later
overwrites `plotting_nature.R` rather than being overwritten by it. Even a future
collision therefore fails safe.

### Pre-existing seeds left alone

Five publication-facing calls already carried `seed = 1`, in
`build_wgcna_modules.R` (3), `render_microglia_module_figures.R` and
`summarize_module_interpretation.R`. They are already deterministic, which is
the property this phase requires. Renaming their literal to the shared constant
would change those figures' layouts on the next render for no correctness gain,
so they were not touched. This was defect repair, not stylistic normalization.

### Aesthetics preserved

Every jitter `width` and `height` is unchanged. The spread of a jittered panel
is an authored aesthetic that affects how the panel reads, so the test pins the
per-script width multiset (`0`, `0.15`, `0.16`, `0.12`, `0.12`, `0.12`,
`0.15`×3, `0.12`, `0.11`, `0.1`×2 — 13 calls carrying a `width`) and fails if a
later edit moves one.

Two calls use their own parameter names and are pinned separately: the
`position_jitterdodge` in `assess_sample_quality.R`
(`jitter.width = 0.12`, `dodge.width = 0.65`) and the EWCE stratum plot, whose
jitter is vertical rather than horizontal (`width = 0`, `height = 0.22`).

## Verified properties

| property | result |
|---|---|
| unseeded layer differs across differing RNG state (the defect reproduces) | confirmed |
| seeded layer byte-identical across differing RNG state | confirmed |
| seeded jitter + seeded repel in one panel, byte-identical | confirmed |
| the same, with the two mechanisms given different seeds (no order dependence) | confirmed |
| `.Random.seed` unchanged by a seeded render | confirmed |
| downstream RNG stream unperturbed — no statistic can shift | confirmed |
| a different jitter seed moves bytes but not axis or category text | confirmed |
| `position_jitterdodge` unseeded differs / seeded identical | confirmed |
| `geom_jitter` has no `seed` formal, warns, and does not honour it | confirmed |
| bare `geom_jitter(width=)` == explicit `position_jitter(width=)`, byte-identical | confirmed |
| ggplot2 errors when `position` and `width` are both supplied | confirmed |

The RNG-isolation rows are the reason this change cannot alter any scientific
result. `position_jitter(seed=)` draws from a layer-local stream via
`withr::with_seed`, restoring the global stream afterwards, so a seeded render
leaves `.Random.seed` exactly as it found it and any downstream permutation test
or bootstrap consumes the same draws it would have consumed before.

## What was not done

- No figure was rerendered. The accepted package is unchanged:
  `figure_export_manifest.csv` still has 5,582 rows and still hashes
  `0fd0c9ed9febbc05ed6928b7fc6bfffdab845daba6255d8b47c607dee5d3c02c`.
- `publication_freeze_manifest.yml` untouched.
- No analysis was rerun; no numerical value, statistic or manuscript claim moved.
- `dplyr::slice_sample()` at `analysis/differential_abundance/compare_go_enrichment.R:2951`
  remains unseeded. It is **statistical** RNG, not render RNG: seeding it would
  change which rows an analysis selects, which is a different decision with a
  different risk profile. Recorded as carried debt, explicitly out of scope here.

## Tests

`tests/testthat/test-jitter-render-determinism.R` — 80 assertions. Determinism
is pinned as "render A == render B" rather than as a coordinate snapshot, so a
legitimate graphics-library upgrade does not produce a spurious failure while a
regression to unseeded placement does.

The fixture deliberately builds its data **once, outside** the render helper.
An earlier version built it inside, which reset the global RNG before each
render and made even the unseeded case deterministic — the defect failed to
reproduce, and the test would have passed against unfixed code. The same trap
was hit in Phase 6H.6.

### The compliance guard is not vacuous

An adversarial review of this change found that the first version of the guard
asked only whether a `seed` argument was *present*. Three ways to defeat that
were closed:

- **`seed = NA` read as compliant.** That is ggplot2's defective default written
  out explicitly, so the guard would have passed the exact defect this phase
  removed. The check now tests the seed's **value** — `NA`, `NULL` and the typed
  `NA_*` variants do not count as pinning the stream.
- **A substring match on "seed".** The `geom_jitter` branch searched the
  deparsed `position` argument for the text "seed", which also matches a symbol
  named `unseeded_pos`. It now parses the position expression and reads its
  `seed` argument.
- **Two calls of one token on a single physical line.** Call text was sliced
  from a regex search, which finds the first occurrence, so the second call
  inherited the first one's seed. Several QC scripts in this repository do put a
  whole plot on one line. Slicing now uses the token's exact parse column.

The guard was then mutation-tested against the live tree: changing one repaired
call to `seed = NA`, and separately deleting its seed altogether, both **fail**
the suite and name the offending file. The mutated file was restored and
verified byte-identical by SHA-256.

Heights are now asserted alongside widths. They were collected by the detector
and never checked, which left half of the aesthetic-preservation contract
unprotected.
