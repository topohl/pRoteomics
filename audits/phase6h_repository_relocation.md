# Phase 6H — repository relocation, preflight and deferral

```
relocation_status    = DEFERRED_BY_USER
relocation_executed  = FALSE
filesystem_changes   = 0
```

Phase 6H evaluated relocating the two sibling repositories to a shorter parent
directory and renaming the scientific checkout. The preflight completed and
established that the relocation would have been technically safe. The user then
deferred it as not currently necessary. **No repository was moved or renamed and
no directory was created.** Both checkouts remain at:

```
S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Analysis\proteomics
S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Analysis\Exp9_manuscript
```

This record exists so the preflight does not have to be repeated if relocation
is revisited.

## Phase-6G marker

Annotated tag `phase6g-output-namespace-migration-complete` →
`58169de7049a2a49f2d35e7a635d61659a99cd9c`, carrying the Phase-6G closure
evidence in its message. Note a naming divergence worth knowing about: the six
pre-existing tags all carry a `-YYYY-MM` suffix; this one does not, because the
name was specified without it.

## What the preflight established

### Runtime independence from the current absolute path

212 old-absolute-path literals across 76 files, classified in
`phase6h_old_path_literal_scan.csv`:

| class | proteomics | Exp9_manuscript |
|---|---|---|
| HISTORICAL_ARCHIVED_SCRIPT | 145 | 0 |
| FROZEN_AUDIT_RECORD | 56 | 0 |
| TEST_FIXTURE | 4 | 0 |
| COMMENT_OR_DOC | 2 | 0 |
| GENERATED_PROVENANCE | 0 | 2 |
| ACTIVE_RUNTIME_DEPENDENCY | 2 | 1 |

All three ACTIVE_RUNTIME candidates cleared on inspection, so the true count is
**0**:

- `R/enrichment/enrichment_io.R:80` and its vendored twin
  `Exp9_manuscript/R/vendor/enrichment_io.R:76` — the string occurs inside
  user-facing guidance text ("Run the repository from a shorter project root
  such as `P:\` before launching workers"), not in a path expression.
- `tools/audit_active_scripts.R:40` — `"S:/Lab_Member/Tobi"` is an entry in a
  detection allowlist, alongside `"C:/"`, `"/Users/"` and `setwd(`. It is a
  pattern the tool *searches for*, not a path the tool depends on.

`repo_root()` resolves `PROTEOMICS_PROJECT_ROOT` → `rprojroot` git-root → an
upward marker walk → error. There is no absolute default anywhere in the chain.
Every data root (`path_raw`, `path_metadata`, `path_external`, `path_processed`,
`path_results`, `path_work`, `path_export`) is defined as `repo_path(...)`, so
nothing reaches outside the checkout — in particular nothing reads the unrelated
raw mass-spec tree at `Exp9_Social-Stress\proteomics`.

### Sibling discovery, both directions

Four tools discover the manuscript repository, all using the desired contract —
explicit override first, sibling-relative default second, no absolute default:

| site | mechanism |
|---|---|
| `tools/audit_cross_repo_boundary.R:27` | `EXP9_MANUSCRIPT_ROOT`, else `<repo_root>/../Exp9_manuscript` |
| `tools/build_restructure_migration_map.R:49` | same |
| `tools/verify_path_contract_rewrites.R:140` | same |
| `tools/enumerate_output_consumers.R:68` | override only, no default |

The manuscript repository never reaches back at runtime: its `..` traversals are
all within-repo (`figures/` → root, `tests/testthat/` → root), and
`tools/import_render_inputs.R` requires `PROTEOMICS_ROOT` explicitly and errors
if it is unset, documenting itself as a one-off import.

Two name dependencies constrain any future rename:

- Three sites in the manuscript repository use `sub(".*proteomics[/\\]", "", f)`
  to build display-relative paths
  (`figures/final_truth_v9_semantics.R:676,881`,
  `figures/final_truth_v9_heatmap_scale_audit.R:75`). The scientific checkout's
  directory name must **end in `proteomics`** or these silently emit absolute
  paths instead of failing.
- `tools/verify_restructure_equivalence.R:45` holds
  `ALLOWED_REPOS <- c("pRoteomics", "Exp9_manuscript")`, but these are values of
  a `destination_repo` data column, not directory names. A directory rename does
  not affect it.

### Payload baseline

| repository | files | bytes |
|---|---|---|
| proteomics | 464,009 | 40,368,577,700 |
| Exp9_manuscript | 2,060 | 502,702,954 |
| total | 466,069 | 40,871,280,654 |

Content fingerprint over the 1,184 scientifically meaningful files
(`exports/`, `audits/`, `config/`, `docs/`, WGCNA tables, manuscript source data
and provenance, plus five named carriers):
rollup SHA-256 `a3e405a27b01c72c111e035a0a5307d01e76cc4eacaf281d581cb770c6f770a2`,
0 unreadable. `results/manuscript/figure_export_manifest.csv` hashes to
`1bf4fe2e16f57ccc3bd684a048e367bae9cdc5c2c92f8d643ad5fcaa1d7b22e9`, matching the
value pinned in `docs/publication_freeze_manifest.yml`.

### Move mechanics

Source and destination would have shared the ancestor `S:\Lab_Member\Tobi`, with
no reparse points anywhere along the path, on one SMB share
(`\\mdc-berlin.net\fs\AG_Hoernberg`, 3.1 TB free). The move would therefore have
been a server-side rename — metadata only, not a 40 GB copy.

## The MAX_PATH finding, which outlives the deferral

`LongPathsEnabled = 0`. R 4.5.1 **enumerates** paths at or beyond 260 characters
through `list.files()` but **cannot open** them. The cliff is exact:

| absolute length | files | R `file.exists()` on a 60-file sample |
|---|---|---|
| ≤ 259 | 300,399 | 60/60 |
| ≥ 260 | 165,670 | 0/60 |

PowerShell 7 reads the same files, because .NET Core opts into long paths through
its manifest regardless of the registry setting. The limitation is specific to
the R runtime — which is the runtime that runs the science.

Excluding archive, `.git` and scratch trees, **724 live files sit at or beyond
the wall** (680 under `data/processed/`, 44 under `results/`), inventoried in
`phase6h_live_over_maxpath_inventory.csv`. Of the 165,670 over-wall files
overall, 164,946 are the two frozen archive trees
`results/manuscript/_superseded_20260622` and
`results/manuscript/_failed_20260901_maxpath` — the latter named for this exact
problem.

**No Phase-6G gated population is affected.** Every gate sits below the cliff:
WGCNA tables max 219, WGCNA figures 259, `exports/` 152, `pride_submission/`
259, `audits/` 146, `config/` 133. No Phase-6G baseline, byte count or hash was
ever computed over a file R could not open.

## Relocation options that were measured

From `phase6h_relocation_candidate_matrix.csv`. The longest repo-relative path is
217 characters, so a prefix of 42 or fewer clears the wall outright.

| parent | scientific repo dir | prefix | files ≥ 260 | longest | headroom |
|---|---|---|---|---|---|
| `…\Exp9_Social-Stress\Analysis` (current) | `proteomics` | 69 | 165,670 | 287 | −27 |
| `S:\Lab_Member\Tobi\exp9` | `proteomics` | 34 | 0 | 252 | 8 |
| `S:\Lab_Member\Tobi\exp9` | `exp9-proteomics` | 39 | 0 | 257 | 3 |
| `S:\Lab_Member\Tobi\exp9` | `exp9-spatial-proteomics` | 47 | 234 | 265 | −5 |
| `S:\Lab_Member\Tobi\e9` | `exp9-spatial-proteomics` | 45 | 6 | 263 | −3 |

The preferred rename `proteomics` → `exp9-spatial-proteomics` costs 13
characters, which is what prevented it from clearing the wall at the natural
parent. `S:\Lab_Member\Tobi\exp9`, `S:\Lab_Member\Tobi\e9` and
`S:\Lab_Member\Tobi\Experiments\exp9` were all free of collisions; the three
neighbouring `proteomics`-like directories
(`Exp9_Social-Stress\proteomics` raw mass-spec data, `Experiments\Proteomics`
documents, `Exp9_Social-Stress_backup`) are none of them git checkouts.

## Verification not performed

Because nothing moved, the post-move gates were not run and were not required:
git identity re-verification, filesystem equivalence diff, old-path independence
proof, multi-CWD bootstrap, and the full test/smoke/boundary/freeze re-run. The
repository state is unchanged from the Phase-6G close at `58169de`, whose gates
already pass.

## A verification caveat found while writing this record

This repository sets `status.showUntrackedFiles = no` in `.git/config`. Plain
`git status --porcelain` therefore reports **nothing** for untracked files, and
a worktree holding stray new files still looks clean. Any cleanliness assertion
must use `git status --porcelain --untracked-files=all`.

Checked with `-uall` at the time of writing, the only untracked entries in the
repository were the six Phase-6H artifacts listed below, so earlier "clean"
claims were correct in substance — but they were established by a weaker test
than they implied. This sits alongside the existing rule that git cleanliness is
not evidence of payload integrity, since `results/` and `exports/**` are
gitignored: the two limitations are different, and both apply.

## Preserved artifacts

| file | rows | contents |
|---|---|---|
| `phase6h_old_path_literal_scan.csv` | 212 | every old-absolute-path literal, classified |
| `phase6h_path_length_distribution.csv` | 21 | lifecycle × guard band, files and bytes |
| `phase6h_live_over_maxpath_inventory.csv` | 724 | every live file at or beyond the wall |
| `phase6h_relocation_candidate_matrix.csv` | 16 | parent × repo-name projections |
| `phase6h_preflight_content_fingerprint.csv` | 1,184 | pre-existing SHA-256 content baseline |

## Follow-on

The MAX_PATH exposure documented above is a property of the current checkout and
does not depend on relocation. It is carried forward as a separate, targeted
audit of active path lengths against the repository's own guards
(`validate_clusterprofiler_output_path_lengths(safe_limit = 240L)` in
`R/enrichment/enrichment_io.R`, and the 260-character budget in
`R/utilities/export_helpers.R`).
