# Output layout

Where analysis output goes, and why there are three places rather than one.

Machine-readable form: [`config/output_layout.yml`](../config/output_layout.yml)
(`output_layout_v2`). It is enforced by
`tests/testthat/test-output-layout-contract.R`, so this document and the code
cannot drift apart.

`config/output_namespaces.yml` is version 1. It is a frozen baseline object and
stays byte-identical. It describes the layout the historical results tree was
written in, and its rules still hold for that tree: historical trees are read
only, and nothing is migrated or deleted automatically.

## Three lifecycles

| Root | Holds | Tracked | Manuscript may cite | Canonical |
| --- | --- | --- | --- | --- |
| `work/` | regenerable intermediates | no | no | no |
| `results/` | canonical scientific results | no | no | yes |
| `exports/` | frozen outward-facing bundles | no, except `manifest.csv` | yes | no |

The distinction is about who is allowed to depend on a file.

**`work/`** can be deleted and rebuilt from raw inputs and config. Nothing may
cite it and no contract may name it. If something under `work/` turns out to be
load-bearing, that is a sign it was misclassified, not a reason to cite it.

**`results/`** is what the analyses own. It is stable enough for another
analysis to consume, and each family has exactly one canonical writer, recorded
in [`config/results_ownership.csv`](../config/results_ownership.csv). The
manuscript does not read this tree.

**`exports/`** is a copy of canonical results, frozen with a manifest carrying
a sha256 per file. This is the only thing outside this repository is allowed to
read. The bundle itself stays untracked because it is derived, but its
`manifest.csv` is tracked: the boundary has to be reviewable in git without
committing megabytes of derived data.

## Canonical result naming

```
results/<domain>/<analysis_id>/<scope>/<child>/<file>
```

- **domain** — the analysis directory under `analysis/`, with no numeric
  prefix: `preprocessing`, `qc`, `spatial_validation`,
  `differential_abundance`, `wgcna`, `enrichment`, `spatial_networks`,
  `integration`, `publication_source_data`.
- **analysis_id** — the canonical owner's script stem. An output is addressed
  by the analysis that owns it, so you can go from a file to the script that
  wrote it without consulting a table.
- **scope** — `global`, or the dataset id for a dataset-specific analysis.
- **child** — `tables/`, `plots/`, `models/`, `manifests/`, `reports/`. Only
  the children an analysis actually writes are created.

So this:

```
results/tables/03_qc_exploration/02_missingness_diagnostics/<dataset>/missingness_diagnostics.xlsx
```

becomes this:

```
results/qc/summarize_missingness/<dataset>/tables/missingness_diagnostics.xlsx
```

Read it left to right: the domain, the analysis, the dataset, the kind of
artefact. Nothing in it requires knowing that missingness diagnostics were once
step 02 of a stage called `03_qc_exploration`.

`plots/`, not `figures/`. In this project `figures/` means an assembled journal
figure, and that belongs to `Exp9_manuscript`. Scientific diagnostic plots
produced by an analysis are plots.

A segment after the scope that carries meaning is preserved. An ontology or a
contrast direction stays in the path:

```
results/differential_abundance/audit_gsea_protein_direction/<dataset>/tables/<ontology>/gsea_contrast_direction_summary.csv
```

Build these paths with the helper rather than by hand, so a typo becomes an
error instead of a new namespace nobody finds:

```r
canonical_result_path("qc", "summarize_missingness.R", dataset, "tables",
                      "missingness_diagnostics.xlsx")
canonical_work_path("preprocessing", "extract_protigy_contrasts.R", dataset)
path_export("publication_source_data", publication_id)
```

An unknown domain or child is rejected against the contract.

## Artifact names

A canonical filename describes the object: `protein_effects.csv`,
`module_membership.csv`, `term_results.csv`, `network_group_tests.csv`.

Avoid `stage*`, `part*`, `final*`, `latest*`, `v[0-9]*`, `revised*`, `new*`,
`old*`. These describe when a file was made, not what is in it, and they stop
being true the next time someone makes a newer one.

Dataset, analysis and scope are already in the path, so repeating them in the
filename adds length without adding information. Prefer one canonical table
with `dataset`, `contrast` and spatial-unit columns over one file per contrast —
but not at the price of reshaping data that is already frozen. Where a
conversion would alter a scientific contract, the existing format stays.

## Legacy outputs

[`config/legacy_output_registry.csv`](../config/legacy_output_registry.csv)
lists the output roots that are read only, with the evidence for each: how many
frozen baseline objects sit beneath it, how many files, and how many registered
writers still target it. A root is legacy when the answer to the last question
is zero. That is a measurement, not a judgement.

The largest group is the manuscript rendering layer. 111 of the 129 frozen
baseline objects under `results/` are produced by renderers that moved to
`Exp9_manuscript` in Phase 6C:

| Root | Frozen objects |
| --- | --- |
| `results/figures/manuscript` | 39 |
| `results/figures/manuscript_candidates` | 35 |
| `results/source_data/manuscript_candidates` | 35 |
| `results/tables/manuscript_candidates` | 2 |

Nothing in this repository writes them any more, and the authoritative copies
are in the manuscript repository, imported and hash-verified in its
`provenance/source_manifests/`. Alongside those sit superseded and comparison
trees kept as provenance: `results/manuscript/_superseded_*`,
`results/manuscript/_failed_*`, the `EWCE_sample_vs_animal_*` comparisons and
the `*_panels` authoring roots.

The policy is deliberately narrow:

- the physical objects are **unchanged**, at the paths the frozen manifests
  already record;
- **reads are permitted**, so results produced before the migration stay
  reachable;
- **writes are forbidden**, and
  `test-output-layout-contract.R` fails if any registered writer declares an
  output beneath a legacy root;
- the registry cannot be quietly narrowed to make that guard pass: the four
  trees above are named in the test, and the frozen-object total is asserted.

This is the forward-facing half of PB-11 and PB-12. The historical artefacts
are not moved or renamed; what changes is that no new canonical run can write
into them, and the manuscript no longer resolves against them.

## The publication source data boundary

```
exports/publication_source_data/<publication_id>/...
exports/publication_source_data/manifest.csv
```

`Exp9_manuscript` imports from here and from nowhere else, and
`tests/testthat/test-import-boundary.R` in that repository enforces it.

The contract for what belongs in the bundle is
[`config/publication_source_data_contract.yml`](../config/publication_source_data_contract.yml):
ten publication identities, each naming the analysis that originates it. This
repository owns that list. It does not own how a journal figure is assembled
from the data.

The bundle carries its own `manifest.csv` with `publication_id`,
`source_analysis`, `source_table`, `exported_file`, row and column counts and a
sha256 per file. Because the Phase 6F relocation was a byte-identical copy,
those hashes validate against the bundle unchanged;
`audits/phase6f_artifact_migration.csv` records the old and new path and hash
for all 55 files.

One transitional detail. That manifest is a byte-identical copy, so its
`exported_file` column still names the pre-6F location,
`results/publication_source_data/...`. The hashes are what the column is for
and they validate against the relocated bundle exactly as they did before.
`tools/export_publication_source_data.R` builds `exported_file` from its bundle
root, which now points at `exports/`, so the next export writes a manifest that
describes itself. Until then, consumers strip either prefix — the manuscript's
`tools/verify_source_bundles.R` does.

The bundle is written by `tools/export_publication_source_data.R`, a tool you
run by hand. No `produces` entry in `pipeline.yml` declares it, which is why an
inventory built from the registry alone reports it as having no writer. That is
a property of the registry, not a missing owner.

One honest note about how this looked before. The bundle the manuscript
imported, `results/publication_source_data/`, sat inside the canonical results
tree and **no registered writer produced it** — it was real, load-bearing, and
ownerless. Separately, the manuscript's 196-row render-input manifest resolved
against arbitrary `results/**` paths. Both are now closed: the bundle has a
declared root and owner, and the render-input manifest is a closed historical
record of the one-off bridge, since this repository produces none of those
files any more.
