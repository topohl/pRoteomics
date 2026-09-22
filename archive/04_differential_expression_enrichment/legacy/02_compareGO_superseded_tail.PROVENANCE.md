# Provenance — `02_compareGO_superseded_tail.r`

Sidecar for the byte-exact archived payload. Nothing in this record appears
inside the archived file itself: the payload is preserved unmodified, so all
explanatory metadata lives here.

| field | value |
|---|---|
| original file | `analysis/differential_abundance/compare_go_enrichment.R` |
| original full-file SHA-256 | `dc9e41d98a309ba7b582576984a7a7fc0e83712deeb8bc3d079e3cb91aef9a91` |
| original line count | 4,321 |
| archived region | lines **570–4321** (everything after the unconditional exit) |
| archived region lines | 3,752 |
| archived region SHA-256 | `cde1bccf1e4f6381b7a40c9cf66389dfb8d22a58b46bf8f548cc33dcd04403ba` |
| retained active prefix | lines 1–569 |
| retained active prefix SHA-256 | `e5bce156c78365e2f393a78cdb2aac9ea7a48ec0476a8f0ab9ec1f0d7450f5d7` |
| source commit (state archived from) | `c9368b9` |
| commit that disabled the region | `6999f47` (2026-07-16, "Enforce canonical enrichment provenance contracts") |
| commit that introduced the region | `bb317dc` (2026-05-06, "Enhance plotting, data loading and analyses") |
| decision | Phase 6H.9 **C2 — ARCHIVE_UNREACHABLE_TAIL** |
| implemented in | Phase 6H.10 |
| runtime status | **NON_RUNNABLE** |
| active replacement | the canonical manifest-driven path in `analysis/differential_abundance/compare_go_enrichment.R` lines 1–569 |

## Byte preservation

The archived payload is the original bytes of lines 570–4321, unaltered — no
reformatting, linting, seeding, renaming, comment edits or line-ending changes.
Verified by reconstruction rather than by inspection:

```
sha256( active_prefix || archived_tail ) == dc9e41d9…aef9a91 == original file
```

That gate is re-asserted by
`tests/testthat/test-comparego-tail-archival.R`, so the archive cannot drift
from the history it represents.

## Why it was archived

The region was unreachable. `compare_go_enrichment.R:569` is a bare, column-0,
unconditional `quit(status = 0, save = "no")` — top-level expression 103 of 439
— and line 571 carried the marker
`# LEGACY_COMPAREGO_TAIL_DISABLED_BY_CANONICAL_EXIT`. Every top-level expression
after it, 336 of 439, could not execute. `quit` was masked neither in the script
nor in any of the eight libraries it sources before the exit, no `trace()` was
in play, and nothing in the repository sourced, parsed or evaluated the file.

Leaving it in the active tree had measurable cost. Two audit phases were spent
establishing that a statistical RNG call inside it was inert; a ggrepel seed was
applied in Phase 6H.6 to a layer that can never render; and a test in
`test-protein-group-enrichment-utils.R` was passing by matching tokens that
exist only here, giving no assurance about the live path.

## What the region contains

The body of the original standalone `compareGO.r` v2.1 — its own section
numbering restarts at 1, and `analysis_params$script` inside it is
`"compareGO.r"`. Seventeen sections:

summary statistics · term consistency · gene importance ranking · comparison
similarity · enriched-term barplot · NES ridge plot · similarity heatmap ·
redundancy report · term-comparison occurrence matrix · top driver genes · UpSet
plot · bootstrap enrichment stability · parameter & reproducibility log · data
quality summary · Sankey diagram · alluvial diagram · term hierarchy

It carries 56 write calls (16 `write_xlsx`, 14 `svg`, 13 `ggsave`,
10 `write_raw_xlsx`, 2 `writeLines`, 1 `write.csv`) and 157 top-level names, of
which 20 are functions. None of those names is depended on from outside: three
share a name with a definition elsewhere (`jaccard`, `mode_value`,
`optional_read_csv`) and in each case the other file defines its own.

## Runtime and audit status

`archive/` is declared non-runnable provenance by
`pipeline_analysis_script_exclusions()` in `R/utilities/pipeline_registry.R`:
"Superseded and exploratory generations live under `archive/` and are
provenance, not runnable stages." Active-code scanners therefore exclude this
file by repository policy, not by accident, and it is not eligible for outward
deposition.

Do not revive it. Phase 6H.9 recorded that no active script computes an
enrichment-stability or term-recovery statistic today, and that the canonical
path dropped the quantity deliberately rather than losing it. The bootstrap
block's own output was adjudicated
`MISLEADING_OR_UNINTERPRETABLE` — it reports a bootstrap row-inclusion
probability, pinned near `1 − e⁻¹ = 0.632`, under a name that implies enrichment
robustness.
