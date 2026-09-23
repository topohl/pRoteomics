# Phase 6I.3 — retention of the executed GSEA ranked order

This phase changed **retention semantics only**. No provenance CSV was
modified, no analysis was rerun, and no outward artifact moved.

## What the artifact is

`analysis/differential_abundance/run_clusterprofiler_enrichment.R:1498` writes

```r
write.csv(data.frame(GeneSymbol = names(gene_inputs$sensitivity$median),
                     median_statistic = unname(gene_inputs$sensitivity$median),
                     max_abs_signed_statistic = ...),
          file.path(audit_root, "rank_statistic_sensitivity_audit.csv"))
```

`gene_inputs$sensitivity$median` is assigned from `ranked` in
`R/enrichment/protein_group_enrichment_utils.R:136`, and `ranked` is the object
passed to `gseGO()`. The file therefore stores **the rank order that executed**,
row by row — not an input from which the order could be recomputed.

Phase 6I.2 classified this `BYTE_EXACT_STORED_ORDER` and verified it across all
54 GSEA_GO/BP comparisons: the stored `GeneSymbol` order equals the order
rebuilt from `collapsed_gene_input.csv`, the stored `median_statistic` is
monotone non-increasing and equal to the rebuilt values, and the row count
equals the manifest's `n_genes`.

## Why it read as disposable

Three things, none of which was a decision about this file:

1. **No registry knew it existed.** The clusterProfiler manifest has columns for
   `collapsed_gene_input_file`, `collapsed_gene_provenance_file` and
   `term_gene_provenance_file`, but **no column for this artifact**. It is
   written and never registered.
2. **Nothing reads it.** `audits/phase6h_over_maxpath_reader_check.csv` records
   `files_on_disk = 24`, `active_code_references = 1`, and that single reference
   is `run_clusterprofiler_enrichment.R:1499` — the *write* site. A reader
   census therefore returns zero readers, and an unread file looks deletable.
3. **The frozen over-wall inventory says `REGENERABLE`.** That column takes only
   two values across all 724 rows, `REGENERABLE | RESULTS`, and is derived from
   the root directory: everything under `data/processed/` is `REGENERABLE` and
   everything under `results/` is `RESULTS`. It is a statement about which tree
   a file sits in, not a retention judgement.

The semantic error being corrected:

> computationally reconstructible **≠** disposable provenance

The reconstruction exists only while `collapsed_gene_input.csv` also survives,
and it is this file that records what actually ran.

## What changed

| | before | after |
| --- | --- | --- |
| declared in any registry | no | yes — `clusterprofiler_protected_reference_artifacts()` |
| role | none | `protected_reference_not_consumed` |
| cleanup eligible | unstated, and `REGENERABLE` by tree | `FALSE`, explicitly |
| output contract row | none | `docs/OUTPUT_CONTRACTS.md`, "clusterProfiler executed ranked order" |
| regression protection | none | `tests/testthat/test-executed-rank-retention.R` |

The role token is **not new**. `protected_reference_not_consumed` is already the
repository's class for a protected artifact that is hashed but never read —
see `stress_response_protected_reference_artifacts()` in
`R/statistics/stress_response_biological_audit_utils.R` and the two Stage-11
rows in `docs/OUTPUT_CONTRACTS.md`.

The declared path set is **derived from the manifest**, not hardcoded: the
artifact is a sibling of each comparison's `collapsed_gene_input_file`. A
registry entry that names a path nothing occupies protects nothing, so the test
asserts the declaration resolves onto instances that exist.

## The artifact family, and where its boundary is

The basename occurs **60** times under `data/processed/`. Only **54** are the
executed publication run:

| tree | instances | in scope |
| --- | --- | --- |
| `04_differential_expression_enrichment/clusterProfiler/` | 54 | **yes** — the run behind the published figures |
| `04_differential_expression_enrichment_comparison/animal_level/` | 3 | no |
| `04_differential_expression_enrichment_comparison/legacy_replay/` | 3 | no |

The six excluded instances are all `neuron_soma / DG_sg`, three comparisons in
each of two branches. They are reproducibility and replay runs: they record the
ranked order of *a* GSEA execution, but not of the execution that produced
Figure 3. Protecting them as publication provenance would overstate what they
are, so the declaration resolves from the three canonical clusterProfiler
manifests and reaches exactly the 54. The test asserts both halves — that all
54 are covered, and that the comparison trees are deliberately not.

This is the boundary that matters, and it is a same-name/different-tree
distinction rather than a same-tree/different-name one.

## The in-scope family, measured

| | |
| --- | --- |
| instances | **54** (neuron_neuropil 30, neuron_soma 12, microglia 12) |
| addressability | 36 `present`, 18 `path_over_limit`, **0 `absent`** |
| path length range | 252–274 characters (wall = 260) |
| total rows | 276,054 |
| total bytes | 12,622,903 |
| distinct sha256 | 54 |
| family fingerprint | `a619e6eaca36e47c0e5e2c1b8f299c36b600c1553bc3d1ee7a4e671fa7ee86c0` |

The fingerprint is the sha256 of the sorted list of the 54 per-file hashes. It
is recorded here as the before/after gate for this phase, not as a test pin: a
hash pin would fail the first time the enrichment is legitimately rerun, and
that is not what "protected" means. Byte and semantic stability are tested in
`tests/testthat/test-gsea-rank-provenance.R` against the manifest's own
`n_genes`.

Eighteen of the 54 are past the 260-character wall — more than the nine
`collapsed_gene_input.csv` instances, because this basename is longer. They
resolve only through the Phase 6H.3 staging contract.

## Frozen historical snapshots

Unchanged, deliberately. `audits/phase6h_live_over_maxpath_inventory.csv` and
`audits/phase6h_over_maxpath_reader_check.csv` correctly record what the
classification was when they were taken. A frozen snapshot edited to agree with
the present stops being evidence. The current declaration supersedes them
prospectively; the test asserts the frozen inventory still says what it said.

Frozen Phase 6H snapshots modified: **0**.

## An adjacent defect, recorded and not fixed here

`docs/file_contracts.tsv` declares `clusterProfiler_protein_group_audits` at

```
results/tables/04_differential_expression_enrichment/clusterProfiler/<dataset>/**/protein_group_audits/*.csv
```

That directory holds **zero** files. The producer writes to
`CANONICAL_PATHS$models`, and the 54 real instances live under
`data/processed/04_differential_expression_enrichment/clusterProfiler/...`.
The contract row points somewhere empty, which is part of why this family read
as unowned.

Correcting it is a statement about the current output namespace rather than
about retention, and it touches the writer-namespace contracts, so it is
recorded here rather than changed in a retention-only phase.
