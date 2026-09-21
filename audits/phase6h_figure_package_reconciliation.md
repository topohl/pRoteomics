# Phase 6H figure package reconciliation

How the manuscript figure package came to disagree with the repository's own
scope policy, and how the two were brought back into agreement. Nothing here
changed a scientific result; every step is packaging, addressing or provenance.

## The three counts, and why they differed

| count | value | what it was |
|---|---|---|
| producer run manifest, legacy address | 5,244 | frozen 2026-09-19, at an address the exporter stopped writing in Phase 6G |
| accepted export, before 6H.5D | 5,598 | the 6H.5 selection, including 16 rows that predated the scope policy |
| accepted export, after 6H.5D | **5,582** | the selection the explicit policy implies |

Two independent problems produced that spread, and neither was a scientific
change:

1. **A stale freeze address.** Phase 6G moved the exporter onto
   `psd_dirs("08_export_manuscript_figures")$manifests`, but
   `freeze_export_payloads()` kept reading the pre-6G
   `results/logs/09_export_pride_journal/manuscript_figures/run_manifest.yml`.
   The exporter was rerun twice during Phase 6H and wrote 5,598 and then 5,582
   inputs canonically, while the freeze compared against a legacy file frozen at
   5,244 and reported `counts_agree: no` for a reason that had nothing to do
   with the export. `freeze_run_manifest_path()` now prefers the canonical
   address and falls back to the legacy one only when the canonical file is
   absent, which is still the case for `09_export_source_data`.

2. **A policy that had never reached the figure exporter.** The proposed-tree
   exclusion existed for `results/tables/` and `results/source_data/` but not
   for `results/figures/`, so 16 `_validation_proposed` figures were selected.

## The 16 removed rows

10 from `compareGO_spatial_atlas_validation_proposed`, 6 from
`microglia_validation_proposed`. Adjudicated in Phase 6H.5C against repository
provenance, not by appearance:

| class | n | basis |
|---|---|---|
| `REMOVE_PROPOSED_DUPLICATE` | 15 | byte-identical to a canonical counterpart that is itself exported |
| `REMOVE_NONDETERMINISTIC_PROPOSED_VARIANT` | 1 | same payload, different unseeded `ggrepel` placement |

The single non-identical file, `Fig_RES_SUS_divergence_publication.svg`, was the
one that needed real evidence. Both variants come from the same script,
`analysis/differential_abundance/build_go_program_atlas.R`, which switches
`SUBSTEP_ID` on a `VALIDATION_ONLY` flag and changes nothing about the
computation. The decisive findings:

- `source_data_RES_SUS_divergence_publication.csv` is **byte-identical** between
  the two scopes (`e1ed9fc9…`, 135,885 bytes, 252 rows x 31 columns, same row
  order), so the upstream scientific payload is the same;
- the SVG label text is the **same 22-element multiset**, with no label unique to
  either side;
- all 11 differing lines are `<text>`/`<line>` elements and every one carries
  coordinates;
- stripping every number from both files makes them **identical line-by-line**;
- the producer calls `ggrepel::geom_text_repel(max.overlaps = 12)` and the script
  contains no `set.seed()`;
- the canonical variant was written **later** (10:40:13 against 10:07:07 on
  2026-08-26), so the proposed one is not a newer correction.

It therefore carries no information the canonical variant lacks, and no
canonical promotion was required before removing it.

An earlier report in this phase described this figure as labelling "different
region/program pairs". That was wrong, and the correction matters: the labels
are the same set, moved.

## Realised delta

```
manifest rows                 5,598 -> 5,582   (-16)
manifest sha256               68850f54... -> 0fd0c9ed...
audit sha256                  f94793ad... -> 6075ed98...

removed selections            16   (15 duplicate + 1 nondeterministic)
additions                      0
unrelated removals             0
failed-run additions           0
canonical counterparts kept   16 / 16

package destinations removed  16   (all physically retired)
package destinations added     0
retained destinations      5,582   with 0 byte changes
```

Source figures were never touched: the 16 proposed sources still exist with
**0 hash changes**, and all **359** figures repaired in Phases 6H.5 and 6H.5B
still hash exactly as recorded. De-selection is not deletion — the proposed
trees remain on disk as validation artifacts.

## Counts after reconciliation

```
figure_export_manifest rows                 5,582
producer run manifest recorded inputs       5,582
freeze manifest_row_count                   5,582
counts_agree                                yes
```

5,244 and 5,598 survive only as `superseded_*` provenance in the freeze record,
alongside the Phase-6G.8 layer beneath them. The legacy run-manifest address is
retained as `superseded_run_manifest_path`.

## Untouched, deliberately

`exports/` 55 files / 2,545,817 bytes · `pride_submission/` 1,310 files · the
three clusterProfiler historical manifests · WGCNA protected states ·
`pipeline.yml` · `results_ownership.csv` · `config/legacy_output_registry.csv` ·
every canonical source-data table.

## Debt this phase exposed and did not fix

`build_go_program_atlas.R` renders with unseeded `ggrepel`, so **its figure bytes
are not reproducible across runs** for identical input data. That is a
reproducibility defect in its own right: it is the reason the two variants of the
divergence figure differ at all, and it means any rerun of that producer yields a
byte-different SVG. Seeding it is deliberately left to a controlled turn, because
doing it here would have required rerendering a figure during a packaging batch.
