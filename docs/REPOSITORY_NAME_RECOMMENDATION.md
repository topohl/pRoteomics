# Repository name: assessment and recommendation

Phase 6F section 14 asks whether `pRoteomics` should eventually become
`exp9-proteomics` or `exp9-spatial-proteomics`, and says not to rename the Git
remote automatically. Nothing has been renamed. This document is the
assessment.

## Current state

| | |
| --- | --- |
| Git remote | `https://github.com/topohl/pRoteomics.git` |
| Local working directory | `.../Exp9_Social-Stress/Analysis/proteomics` |
| Sibling repository | `Exp9_manuscript` (no remote; local only) |

The remote and the working directory already disagree in case and spelling:
`pRoteomics` against `proteomics`. Documentation and paths in this repository
say `proteomics`, while provenance records say `topohl/pRoteomics`.

## Recommendation

**Rename to `exp9-spatial-proteomics`**, as an administrative step after the
`repository-architecture-migration` branch is reviewed and merged.

Three reasons.

**It says which experiment it belongs to.** The sibling repository is already
`Exp9_manuscript`. A reader who finds `pRoteomics` cannot tell whether it is
Exp9's analysis repository or a general-purpose proteomics toolkit, and the
lab has more than one experiment. `exp9-` removes that ambiguity in the place
people actually see it: the clone URL.

**It says what kind of proteomics.** The science here is spatial: 18 spatial
units, bilateral aggregation to animal level, CA2-SLM robustness, spatial
networks and a spatial atlas. `exp9-proteomics` is accurate but loses the one
word that most distinguishes this dataset from an ordinary bulk proteomics
project.

**The intercapped R has a cost and no remaining benefit.** `pRoteomics` is a
pun on R, and it has produced a real inconsistency: the remote is `pRoteomics`
and the checkout is `proteomics`. Case-sensitive tooling and case-insensitive
Windows filesystems disagree about whether those are the same word, which is
exactly the class of problem this migration has been removing everywhere else.

`exp9-proteomics` is the acceptable second choice if a shorter name is
preferred. Either is better than the status quo.

## Why not now

A remote rename is not a code change and should not ride along with one.
GitHub leaves a redirect, but the following would need a deliberate pass, and
doing it mid-migration would mix two kinds of risk:

- the clone URL in any CI configuration, and `.github/` workflows;
- `renv` and any lockfile or cache path that embeds the repository name;
- the sibling repository's own references, twelve of which are in its test
  suite;
- documentation and README links in both repositories.

It should happen once, after the migration branch is merged, so there is a
single before-and-after rather than a rename layered on top of an unreviewed
restructure.

## What must not be rewritten

Provenance records name the repository as it was at the commit they describe.
The frozen manifests and the manuscript's import records carry
`source_repo: topohl/pRoteomics`, and that is correct for those commits:

- `exports/publication_source_data/manifest.csv` (`source_repo` column);
- `Exp9_manuscript/source_data/pRoteomics/manifest.csv`;
- `Exp9_manuscript/provenance/source_manifests/render_inputs_manifest.csv`.

A rename must not retroactively edit them. The same rule the rest of this
migration has followed applies here: records of history keep the historical
name, and only live interfaces move. On the manuscript side the local bundle
directory `source_data/pRoteomics/` is a path in a frozen import, so it stays
as well.

## Checklist, for when it happens

1. Merge `repository-architecture-migration`.
2. Rename the GitHub repository; keep the redirect.
3. Update the clone URL in `.github/` workflows and any CI configuration.
4. Update README links in both repositories.
5. Rename the local checkout from `proteomics` to match, so the remote and the
   working directory finally agree.
6. Leave every `source_repo`, manifest and provenance record untouched.
7. Do not rename `Exp9_manuscript/source_data/pRoteomics/`.
