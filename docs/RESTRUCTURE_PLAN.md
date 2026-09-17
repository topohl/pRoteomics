# Phase 6B/6C — Repository Architecture Migration Plan

Structural migration only. No scientific result, numerical value, statistical
analysis, figure content or manuscript claim changes in this phase.

## 1. Authoritative baseline

| Item | Value |
| --- | --- |
| Baseline commit | `6801edbce8a5d222f4af46e06b6db4e99f6a9761` |
| Baseline tag | `pre-restructure-scientific-freeze-2026-09` |
| Baseline manifest | `manuscript/prerestructure_freeze_manifest.csv` |
| Manifest rows | 239 |
| Distinct filesystem objects | 200 |
| Migration branch | `repository-architecture-migration` |

### 1.1 The manifest is 239 assertions over 200 objects

The freeze manifest records one row per *(object, publication role)* pair, not
one row per file. 14 paths appear more than once because the same canonical
table feeds several panels:

- `results/tables/10_biological_integration/gsea_wgcna_concordance/global/ontology_aware_gsea_theme_assignments_all_contrasts.csv` — 11 rows (5 Extended Data 6 panels, 5 Figure 3 panels, 1 protected-state assertion)
- `manuscript/figure1_bridge_mmmsociability/behavior_sex_effect_contract.csv` — 5 rows
- `manuscript/figure1_bridge_mmmsociability/source_data/figure1_panel_statistics.csv` — 5 rows
- `manuscript/figure1_bridge_mmmsociability/source_data/figure1a_timeline_source.csv` — 5 rows
- 10 further paths at 2–3 rows each

239 rows − 200 distinct paths = 39 repeat references. Every duplicated path
carries an identical `sha256` across its rows (verified: 0 inconsistent).

Consequence for the equivalence oracle: `unique_baseline_paths = 239` as stated
in the phase brief is arithmetically unreachable. The oracle is therefore
enforced at two levels, both of which must hold:

| Level | Invariant |
| --- | --- |
| Assertion | 239 manifest rows each resolve to an existing destination whose sha256 matches |
| Object | 200 distinct objects each map to exactly one destination; no object maps twice |

This is the limitation the phase brief anticipated in §17 and is implemented in
`tools/verify_restructure_equivalence.R`.

### 1.2 130 of 200 frozen objects are untracked

| | tracked | untracked |
| --- | --- | --- |
| `results/` | 0 | 129 |
| `data/` | 0 | 1 |
| `manuscript/` | 50 | 0 |
| `tests/` | 10 | 0 |
| `config/` | 4 | 0 |
| `docs/` | 3 | 0 |
| `figures/` | 2 | 0 |
| `pipeline.yml` | 1 | 0 |
| **total** | **70** | **130** |

The 130 untracked objects are gitignored regenerable canonical outputs. Per §4
they are inventoried and classified but **not** physically relocated; they are
recorded in the migration map with `migration_class = RETAINED_LEGACY_PATH`,
which §16 explicitly permits. This is also why the `results/` reorganisation in
§6 is expressed as a forward interface rather than a bulk move.

## 2. The three-way constraint and how it was resolved

The phase brief requires all three of the following, and they cannot all hold:

1. §14/§28 — the manuscript layer is physically removed from pRoteomics;
2. §40 — the pRoteomics suite passes at HEAD;
3. §16/§34 — all 239 frozen assertions are byte-identical (`hash_mismatch = 0`).

Evidence chain:

- `pipeline.yml` is a frozen `configuration_contract` object. Its
  `manuscript_candidates` stage spans 883 lines and registers 51
  `figures/*.R` renderer scripts.
- `tests/testthat/test-pipeline-registry.R` is also a frozen object, and
  asserts `validate_pipeline_scripts_exist(registry)` — every registered
  script must exist on disk.
- Extracting the renderers makes those 51 paths dangle, so the suite fails.
  Repairing it requires editing `pipeline.yml` or the test, both frozen.

**Resolution adopted: full split with re-freeze of path contracts.** The
manuscript layer is removed, the architecture is re-treed, and the small set of
*non-scientific path-contract* objects that encode the old layout is rewritten
and re-frozen under a new tag. Every scientific and publication object stays
byte-identical.

Objects permitted to change hash in this phase are confined to layout
contracts, and each is logged in `audits/restructure_migration_map.csv` with
`migration_class = REWRITTEN_PATH_CONTRACT`:

| Object | Why it must change |
| --- | --- |
| `pipeline.yml` | stage paths + removal of the `manuscript_candidates` renderer stage |
| `tests/testthat/test-pipeline-registry.R` | asserts registry script existence and stage-relative paths |
| `tests/testthat/test-output-namespace-contract.R` | asserts namespace/stage paths |
| `config/output_namespaces.yml` | output namespace roots |
| `docs/publication_freeze_manifest.yml` | records renderer and stage paths |
| `docs/MANUSCRIPT_STATISTICAL_CONTRACT.md` | cites stage-relative script paths |
| `config/clusterProfiler_config.yml`, `config/manuscript_spatial_order.yml` | stage-relative input roots |

No table of numbers, no panel, no figure, no prose and no statistic is touched
by any of these edits. The re-freeze tag is
`post-restructure-architecture-freeze-2026-09`.

## 3. Destination split

From `docs/restructure_inventory.csv` (786 rows, 0 `UNKNOWN`):

| destination | objects | frozen objects |
| --- | --- | --- |
| pRoteomics | 634 | 148 |
| Exp9_manuscript | 152 | 52 |

Exp9_manuscript is created at
`S:\Lab_Member\Tobi\Experiments\Exp9_Social-Stress\Analysis\Exp9_manuscript`
as a sibling of this repository. It is local/private; no public remote is
created or pushed.

## 4. Analysis / publication boundary

Verified before extraction: the five canonical renderers named by
`manuscript/canonical_publication_registry.csv` contain **zero** occurrences of
`lm(`, `glm(`, `lmer(`, `bam(`, `cor.test(`, `p.adjust(`, `gseGO(`, `fgsea` or
`WGCNA`. The rendering/inference boundary already holds in practice; this phase
makes it structural.

The renderers transitively source 28 `R/` libraries. These split as:

- **Panel implementation libraries** (generation-named: `final_truth_v9_*`,
  `editorial_v8_*`, `nature_final_v7_*`, `nature_v2_*`, `story_v3/v4/v5_*`,
  `spatial_v6_*panels*`, `candidate_figure_*`, `manuscript_figure*_utils`) →
  move to `Exp9_manuscript/R/panels/`.
- **Scientific and infrastructure libraries** (`enrichment_io.R`,
  `module_contracts.R`, `spatial_grammar_utils.R`,
  `sus_res_spatial_dap_atlas_utils.R`, `dataset_config.R`, `paths.R`,
  `validation_utils.R`, …) → stay in pRoteomics.

The canonical v9 generation inherits panel implementations from every earlier
generation, so those earlier generations are **`ACTIVE_SUPPORT`, not
archivable** (§19). They move with the manuscript layer rather than to
`archive/`.

## 5. Known deferrals

| Item | Status | Reason |
| --- | --- | --- |
| PB-11 `manuscript_candidates/final_truth_v9` output namespace | DEFERRED | Resolution requires relocating 129 untracked frozen result objects, which §4 forbids bulk-moving; the forward interface avoids inheriting the name |
| PB-12 `spatial_v6` fingerprint source table | DEFERRED, documented | Phase 5 established this is the only generation that ever existed and its content is current; rehoming is optional and would move a frozen untracked object |
| Renderer repoint onto frozen `source_data/` | Phase 6D | §15 forbids editing renderers during extraction; they move byte-identical and are not executed in the manuscript repo this phase |

## 6. Commit sequence

| Commit | Content |
| --- | --- |
| A | Target architecture, inventory, this plan. No moves. |
| B | Reorganise reusable R functions and analysis entrypoints |
| C | Normalise canonical results and publication source-data paths |
| D | Separate audits, tools and superseded provenance |
| E | Extract the manuscript publication layer |
| F | Remove manuscript rendering from the scientific repository |
| G | Post-migration equivalence guards |
