# renv Lockfile Status

## Current status

`renv.lock` records **327 packages**: the 77 direct dependencies of the active
project plus their transitive `Depends`/`Imports`/`LinkingTo` closure, with
base-priority packages excluded. It pins **R 4.5.1** and **Bioconductor 3.22**,
and every recorded version equals the version installed in the analysis library.
The scientific stack is present at its frozen versions:

| Package | Version | Source |
|---|---|---|
| clusterProfiler | 4.18.4 | Bioconductor |
| fgsea | 1.36.2 | Bioconductor |
| DOSE | 4.4.0 | Bioconductor |
| BiocParallel | 1.44.0 | Bioconductor |
| limma | 3.66.0 | Bioconductor |
| WGCNA | 1.74 | CRAN |

Validate completeness with:

```powershell
Rscript tools/audit_renv_lock.R
```

## Why the lockfile was previously incomplete

The earlier lockfile held three records (`renv`, `yaml`, `testthat`). The cause
was not a misconfigured snapshot: **this project has no renv infrastructure at
all**. There is no `.Rprofile`, no `renv/activate.R`, no `renv/settings.json`,
no `DESCRIPTION`, no `.renvignore`, and no project library — `.libPaths()` has a
single entry, the user R library. The `renv` package itself is not installed.

So the old lockfile could not have been produced by `renv::snapshot()`; it was
hand-seeded to document the lightweight contract-test bootstrap, exactly as its
own status note said. It had also drifted: it recorded `yaml 2.3.10` and
`testthat 3.2.3` while the library actually holds 2.3.12 and 3.3.2.

## What changed

`renv.lock` is now generated from the installed library by
`tools/generate_renv_lockfile.R`, backed by `R/utilities/renv_lock_audit.R`. Nothing is
installed, updated or loaded; versions and source metadata are read from
installed `DESCRIPTION` files only.

- `renv` was **removed** from the lockfile. It is not installed, so no version
  could be recorded faithfully, and it appears in active code only inside
  error-message text suggesting `renv::restore()`. `renv::restore()` bootstraps
  renv itself, so nothing is lost.
- `yaml` and `testthat` records were corrected to the installed versions. No
  package was upgraded; the record was corrected to what is already present.

## How dependencies are discovered

Two passes, because neither alone is sufficient:

1. **Active-surface scan.** 231 files: the 78 registry-declared pipeline scripts
   from `pipeline.yml`, `R/` helpers, `tools/`, and `tests/`. Archived trees
   (`99_deprecated/`, `90_testing/`, `legacy/`, `_scratchpad/`) are excluded, so
   a package used only by deprecated code is not treated as a dependency.
   Patterns matched: `library`, `require`, `requireNamespace`, `loadNamespace`,
   `getNamespace`, `pkg::`, `pkg:::`.
2. **Closure.** `tools::package_dependencies(which = c("Depends","Imports","LinkingTo"), recursive = TRUE)`
   against `installed.packages()`.

Nine scanner hits are excluded as verified false positives: seven string
literals that contain `::` and are used as composite keys (`"Neuropil::CA1"` and
similar), the loop parameter `pkg`, and `renv`. The two places that load
packages dynamically —
`analysis/differential_abundance/run_clusterprofiler_enrichment.R` (`master_packages`)
and `analysis/enrichment/run_ewce_celltype_enrichment.R` (`cran_packages` /
`bioc_packages`) — resolve to literal vectors that the static scan already
covers.

## What the lockfile represents

The library that produced the publication outputs:
`C:/Users/topohl/AppData/Local/Programs/R/R-4.5.1/library`, under R 4.5.1 and
Bioconductor 3.22. Recorded provenance, taken verbatim from installed metadata:

- 279 packages with `Repository: CRAN`
- 40 Bioconductor packages recording `https://bioc-release.r-universe.dev`
  (the R-universe mirror of the Bioconductor release), so that repository is
  declared in the lockfile
- 5 recording `Bioconductor 3.22`
- 3 annotation packages (`GO.db`, `org.Hs.eg.db`, `org.Mm.eg.db`) whose
  `DESCRIPTION` carries no `Repository` field; the field is omitted rather than
  guessed

`janitor` and `snakecase` carry `RemoteType: standard` with
`RemoteRef: <package>` and `RemoteSha: <version>`. That is the signature of a
`pak`/`remotes` install from CRAN, not a git remote, so they are recorded as
CRAN and no git metadata is invented.

## Known limitations

- **No per-package `Hash`.** renv computes it with an internal algorithm that
  cannot be reproduced without renv installed. The field is omitted rather than
  fabricated, so renv cache lookups will fall back to reinstalling.
- **No renv project infrastructure.** `renv::restore()` from a clean checkout
  still requires installing renv and activating the project first; the lockfile
  describes the environment but does not bootstrap it.
- **Exact historical versions.** 213 of 327 recorded versions match current
  repository heads; the remaining 114 are older and would need the CRAN archive
  or a Bioconductor release archive to reinstall. This is expected — it is the
  consequence of pinning the environment that was actually used rather than
  today's latest.

## Restore validation performed

| Level | Result |
|---|---|
| A — structural | **PASS**: parses, 327 records, 0 malformed, 0 duplicates, 0 unresolved `Requirements` |
| B — availability | **PASS**: all 327 package names resolvable from CRAN + Bioconductor 3.22 repository indexes |
| C — isolated restore | **not attempted** |

Level C was not attempted deliberately. It would require installing renv and
building a project library, which would materially alter the live analysis
environment that the freeze is meant to preserve.

## Refreshing the lockfile

```powershell
Rscript tools/generate_renv_lockfile.R --dry-run   # report only
Rscript tools/generate_renv_lockfile.R            # rewrite renv.lock
Rscript tools/audit_renv_lock.R                   # sentinel check
Rscript tools/validate_publication_freeze.R       # freeze contract
```

Run this only on the analysis machine whose library produced the canonical
results. The generator records what is installed; it never resolves versions
from a repository, so it cannot silently upgrade the recorded environment.
