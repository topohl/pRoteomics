# Publication Freeze

`docs/publication_freeze_manifest.yml` is the machine-readable record of the
exact state accepted for publication. It is a *record*, not a build product: it
is written by hashing and reading artifacts that already exist, and generating
it never recomputes an analysis, regenerates a figure or source-data table, or
re-exports a deposition payload.

## What the manifest asserts

| Section | Records |
|---|---|
| `freeze_identity` | freeze commit (full SHA), freeze tag, branch, upstream agreement, working-tree cleanliness |
| `software_environment` | R, Bioconductor and the six load-bearing analysis packages, each as freeze version vs observed version vs match |
| `gsea_reproducibility_contract` | seed base, `nPermSimple`, RNG kind, and the config/implementation files they are read from |
| `export_provenance_equivalence` | the historical-vs-freeze commit proof for the export-relevant files |
| `wgcna_protected_identities` | the four protected WGCNA artifacts plus `identity_contract_sha256`, `membership_version`, `frozen_state_sha256` |
| `publication_source_data` | per-set file lists with SHA-256 and a single stable set digest, for Figure 2, Figure 3 and SUS-RES/Stage 11 |
| `export_payloads` | both payloads: manifest path, manifest SHA-256, row count, recorded historical commit, acceptance basis |
| `validation_state` | the full-suite result associated with the freeze |
| `known_gaps` | limitations that are accepted rather than hidden |

## Provenance-equivalence rationale

Both export `run_manifest.yml` files record `git_commit = 0825e42…`. That commit
is no longer reachable from `main`: the merge and push rewrote SHAs after the
exports ran. The payloads were **not** re-exported to change that label, because
re-running could only alter the label, not the content.

Instead the manifest carries a proof. For each export-relevant file it records
the git blob oid and SHA-256 at both the historical export commit and the freeze
commit. A git blob oid is itself a content hash over the exact bytes, so equal
oids are already proof of byte-identity; the SHA-256 values are recorded so the
claim is checkable without relying on git internals.

The conclusion is deliberately narrow:

> The existing export payloads are accepted for the publication freeze because
> the export-relevant implementation is byte-identical across the historical
> export commit and the freeze commit; therefore the stale `git_commit` field is
> a provenance-label drift, not evidence of payload drift.

This is **not** a claim that the two commits have identical repository trees.
The scope is exactly the file list in `freeze_protected_export_files()`, and a
test fails if that scope is silently widened or narrowed.

## Reproducibility contract

Frozen: the ranked input, the per-comparison derived seed, `RNGkind`,
`nPermSimple`, `by = "fgsea"`, and clusterProfiler's logical `seed` flag (kept
`FALSE`). Numerical reproducibility is expected within one execution context.

Bit-exact equality across *different* execution contexts is **not** guaranteed.
`R/enrichment/clusterprofiler_reproducibility.R` records a measured tolerance instead
(`enrichmentScore` to ~1.6e-15, propagating to ≤ ~2.4e-05 in NES and ≤ ~2.3e-05
in p-value/FDR), with no FDR-0.05 crossings and an unchanged Figure-2f display
selection in the audited comparison. The manifest restates that contract and
does not strengthen it.

## Known gaps

These are recorded as `warn`, surface as WARN in the validator, and are not
fatal:

1. **`renv.lock` incomplete** — pins only `renv`, `yaml` and `testthat`, so the
   freeze is not reproducible from the lockfile alone. See
   `docs/RENV_LOCK_STATUS.md` for the safe refresh procedure. Untouched here.
2. **Export dry-run semantics** — only `08_export_manuscript_figures.R` and
   `09_export_source_data.R` implement a side-effect-free dry-run guard.
   `RUN_EXPORT.R` restricts the step list under `--dry-run`, so the orchestrator
   is safe, but the other seven steps remain individually unguarded.
3. **Wildcard glob filter** — `R/utilities/export_helpers.R:658`, in
   `processed_files_for_dataset()`, filters globs with
   `grepl("\\*", g, fixed = TRUE)`, which searches for a literal backslash-star
   and therefore never matches an ordinary `*`. The filter degenerates to a
   dataset-name test, keeping 1 of 6 configured globs for `microglia` and 0 of 6
   for the other datasets. It affects the PRIDE processed-data selection and is
   masked because both callers pass
   `include_derived = has("--include-derived-results")`, which is `FALSE` by
   default.

## Validate or regenerate

```powershell
# read-only check of the current checkout against the freeze
Rscript tools/validate_publication_freeze.R

# rewrite the manifest from current state (fails if the environment drifted)
Rscript tools/generate_publication_freeze_manifest.R --test-log <suite.log>
```

The validator distinguishes **PASS**, **WARN** (documented accepted gap) and
**FAIL** (changed frozen source-data, changed WGCNA state, drifted package
version, changed export manifest), and exits non-zero only on FAIL.

The generator refuses to write if the active R installation disagrees with the
asserted publication environment or if a protected WGCNA state has changed, so a
drifted state cannot be silently recorded as a new freeze.

## Determinism

Lists are built in fixed sorted order, paths are project-relative, and the only
volatile field — `metadata.generated_at` — is excluded from every section that
identifies scientific content. Regenerating on an unchanged checkout changes
nothing but that field.
