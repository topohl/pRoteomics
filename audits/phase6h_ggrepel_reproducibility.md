# Phase 6H.6 — deterministic ggrepel label placement

## The defect

`ggrepel` positions labels by a randomised physical simulation and its `seed`
argument defaults to `NA`, meaning "draw from whatever the global RNG stream
currently holds". Two renders of identical data therefore produce different
label and leader-line coordinates.

Phase 6H.5C found this the hard way. Two variants of
`Fig_RES_SUS_divergence_publication.svg` differed while:

- the source-data CSV behind them was **byte-identical**
  (`e1ed9fc9…`, 135,885 bytes, 252 rows x 31 columns, same row order);
- the SVG label text was the **same 22-element set**, nothing unique to either;
- all 11 differing lines were `<text>`/`<line>` elements, every one
  coordinate-bearing;
- stripping every number made the two files **identical line-by-line**.

A figure whose bytes move when nothing scientific moved cannot be checked
against a freeze, and it invites exactly the misreading that happened here: the
variant was briefly suspected of containing a withheld scientific correction.

## The contract

`NATURE_REPEL_SEED <- 20260824L`, defined once in `R/utilities/plotting_nature.R`
beside the other `NATURE_*` publication constants, and passed per layer.

Three properties, each measured rather than assumed:

| property | result |
|---|---|
| identical data, identical global RNG state, unseeded | renders agree — so the defect only appears when upstream work has moved the stream |
| identical data, **differing** global RNG state, unseeded | renders **differ** — this is the real-world case |
| identical data, differing global RNG state, **seeded** | renders are **byte-identical** |
| `.Random.seed` before vs after a seeded render | unchanged; the downstream stream is unperturbed |
| different seed values | layout moves, label set unchanged — so the seed governs placement only |

The last two matter. A layer-local seed is the narrowest mechanism available:
`ggrepel` restores the global stream itself, so no `set.seed()` at script scope
is needed and no unrelated stochastic code is perturbed. And because the seed
fully determines placement regardless of incoming RNG state, a canonical run and
a `--validation-only` run of the same producer now agree — which is precisely
the discrepancy Phase 6H.5C documented. The output path never enters the plot,
so two destinations yield byte-identical SVGs.

The value is arbitrary and deliberately so: nothing scientific depends on it and
only its stability matters. It is written in the date-shaped form the repository
already used for seeds elsewhere, and it is not derived from a path or a
timestamp.

## Coverage

12 active `ggrepel` layers across 8 files, enumerated by parsing each call
rather than grepping a line window — a window bleeds into the next layer and
credits it with a neighbour's seed, which is how a first pass mis-attributed two
of them.

| | n |
|---|---|
| active layers | 12 |
| seeded before this phase | 6 |
| newly seeded with `NATURE_REPEL_SEED` | 6 |
| unseeded after | **0** |

Newly seeded: `build_go_program_atlas.R:931`, `compare_go_enrichment.R:2451`,
`assess_marker_rank_abundance.R:190`, and `score_module_activity.R` at `:1727`,
`:1801`, `:1861`. Three of those files gained a
`source(repo_path("R", "plotting_nature.R"))` line; that library is
side-effect-free (4 constant assignments, 15 function definitions, no top-level
calls), so sourcing it more widely is safe. All 6 layers using the constant are
in files that reach it.

The 6 pre-existing seeds were left alone. They are already deterministic, and
changing them would move existing label layouts for no reproducibility gain.
They are inconsistent between themselves — `20260817L`, `20260724L`, `42`, `9` —
which is a tidiness matter, not a correctness one, and not worth churning
published figure layouts to unify.

## Software environment

Determinism is claimed **within the declared environment**, not across versions:

```
R 4.5.1 (2025-06-13 ucrt)
ggplot2 4.0.2
ggrepel 0.9.8
svglite 2.2.2
```

`svglite` embeds no timestamp, so byte identity is achievable and is asserted
rather than weakened to a geometry comparison.

## The accepted payload was not touched

The currently accepted figure is scientifically correct; its arbitrary label
arrangement does not justify churning a frozen publication payload. This phase
governs future renders only.

```
Fig_RES_SUS_divergence_publication.svg   43,462 bytes, unchanged
source_data_RES_SUS_divergence_...csv    e1ed9fc9..., 135,885 B, unchanged
figure_export_manifest.csv               5,582 rows / 0fd0c9ed..., unchanged
publication_freeze_manifest.yml          unchanged
exports/                                 55 files / 2,545,817 bytes
pride_submission/                        1,310 files
```

## Debt this phase found and did not fix

**Three unseeded `position_jitter()` calls in publication figure producers**:
`assess_marker_rank_abundance.R:220` and `score_module_activity.R:1799, :3274`.
`ggplot2::position_jitter()` also defaults to `seed = NA`, and none of those
scripts calls `set.seed()` anywhere, so the *point* positions in those figures
are non-reproducible for the same reason the labels were. Seeding `ggrepel`
does not address it.

Also noted: `compare_go_enrichment.R:2949` uses `dplyr::slice_sample()` for a
bootstrap. That is a statistical computation rather than a rendering detail and
needs separate treatment from figure layout.

Neither was fixed here: this phase was scoped to label placement, and widening
it would have meant changing point geometry in figures during a batch whose
premise is that no accepted figure changes.
