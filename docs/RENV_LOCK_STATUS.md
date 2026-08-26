# renv Lockfile Status

## Current status

As audited on 2026-08-26, `renv.lock` contains three package records (`renv`,
`testthat`, and `yaml`). It is sufficient to document the lightweight contract
test bootstrap, but it is **not** a complete description of the scientific
analysis environment. In particular, it does not record the packages used for
differential abundance, enrichment, WGCNA, mixed models, or publication
figures.

Run the dependency-free audit with:

```powershell
Rscript tools/audit_renv_lock.R
```

The audit deliberately exits nonzero while expected scientific sentinel
packages are absent. Sentinel presence is only a plausibility check; it does
not prove that every dependency is present or that the environment reproduces
the analysis.

## Safe refresh workflow

Refresh the lockfile only on the intended data-containing analysis machine,
using the R and Bioconductor libraries that actually produced the canonical
results. Do not infer or pin package versions from static imports.

From the repository root in PowerShell/R:

```r
install.packages("renv", repos = "https://cloud.r-project.org")
renv::activate()
deps <- renv::dependencies(".")
stopifnot(!anyNA(deps$Package), all(nzchar(deps$Package)))
required <- sort(unique(deps$Package))
missing <- setdiff(required, rownames(installed.packages()))
stopifnot(length(missing) == 0L)
renv::snapshot(packages = required, prompt = FALSE)
renv::status()
```

Then review the lockfile diff, run `Rscript tools/audit_renv_lock.R`, execute
the full repository validation suite, and verify `renv::restore()` in a clean
checkout before treating the lockfile as a reproducible scientific snapshot.
Commit `renv.lock` together with any `renv` activation files intentionally
created by that workflow. The present repository does not include committed
activation files, so `renv::restore()` alone is not currently a supported
full-analysis bootstrap.
