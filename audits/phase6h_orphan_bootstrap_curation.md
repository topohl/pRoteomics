# Phase 6H.9 — curation of the orphaned bootstrap artifact and the unreachable code tail

Adjudication only. Nothing was moved, deleted, seeded, rerun or regenerated. The
workbook remains in `pride_submission/`; the code tail remains in place.

Two decisions are taken independently, as the brief requires, and they are
reached on different grounds:

- **Code:** `C2 — ARCHIVE_UNREACHABLE_TAIL`
- **Artifact:** `P2 — RETAIN_AS_INTERNAL_PROVENANCE_ONLY`

Evidence tables:
- `phase6h_orphan_bootstrap_artifact_inventory.csv` — 5 rows, one per surviving copy
- `phase6h_unreachable_comparego_tail.csv` — 14 rows, the structural proof and tail metrics
- `phase6h_unreachable_comparego_tail_outputs.csv` — 56 rows, the tail's historical output surface

## A. Artifact identity

| | |
|---|---|
| path | `pride_submission/supplementary_tables/results_tables_04_differential_expression_enrichment_compareGO_neuron_neuropil_BP_phenotype_within_unit_08_Bootstrap_Stability_Summary.xlsx` |
| sha256 | `024d3671f2cd6026e8bfd7feb6d9839c7ff23eeb41c98d8f66235d854185baba` |
| size | 5,141 bytes |
| mtime | 2026-05-27 16:05:13 |
| sheets | **1** (`Sheet1`) |
| dimensions | **1 row × 5 columns** |
| columns | `Mean_Recovery_Rate`, `SD_Recovery_Rate`, `Min_Recovery`, `Max_Recovery`, `Total_TopTerms` |
| values | 0.8738182, 0.04100906, 43, 53, 55 |
| metadata / provenance sheet | **none** |
| caption or methods note | **none** |

Five copies survive in total, all 1×5, all with distinct SHA-256, **none** with a
metadata or provenance sheet:

| copy | Mean_Recovery_Rate | Total_TopTerms | bytes |
|---|---|---|---|
| `pride_submission/…neuron_neuropil…` (2026-05-27) | 0.8738182 | 55 | 5,141 |
| `_superseded_20260622/…neuron_neuropil…` (2026-06-13) | 0.8846032 | 63 | 5,143 |
| `_superseded_20260622/…microglia…` | 0.7208 | 25 | 5,135 |
| `_superseded_20260622/…neuron_soma…` | 0.650625 | 16 | 5,136 |
| `_superseded_20260622/08_Bootstrap_Stability_Summary.xlsx` | 0.6535714 | 14 | 15,431 |

The artifact is uniquely identified by its SHA-256; no two copies share one.

Note that the two neuropil copies differ in `Total_TopTerms` (55 vs 63) as well
as in the statistic, so their difference reflects an input change as well as a
different RNG draw. It demonstrates that the **file** is not reproducible; it is
not clean evidence of RNG-induced variation on its own.

## B. Scientific meaning

### What is actually computed

```r
sampled_df  <- combined_df |> group_by(Comparison) |> slice_sample(prop = 1, replace = TRUE) |> ungroup()
sampled_sig <- unique(sampled_df$Description[sampled_df$p.adjust < 0.05])
overlap     <- sum(top_terms %in% sampled_sig)
```

Only the **boolean** `p.adjust < 0.05` is consulted. `NES` is never referenced,
and the magnitude of a p-value below the threshold is irrelevant. No enrichment
statistic, p-value or multiple-testing correction is recomputed.

Tested rather than asserted. Holding the significance *pattern* fixed and
varying everything else about the biology:

| variant | rate | sd |
|---|---|---|
| `p.adjust = 0.0499`, `NES = +0.9` (barely significant, weak effect) | 0.6677500 | 0.1023465 |
| `p.adjust = 1e-30`, `NES = +4.8` (overwhelming, huge effect) | 0.6677500 | 0.1023465 |
| `p.adjust = 1e-12`, `NES = −3.3` (**opposite direction**) | 0.6677500 | 0.1023465 |

Identical to seven decimal places. Effect size, p-value magnitude and even the
**direction** of regulation are invisible to the statistic.

Varying only the multiplicity of the significance pattern:

| each top term significant in | observed rate | 1 − e^−M |
|---|---|---|
| 1 comparison | 0.6505 | 0.6321 |
| 2 comparisons | 0.8705 | 0.8647 |
| 3 comparisons | 0.9535 | 0.9502 |
| 4 comparisons | 0.9855 | 0.9817 |

(Observed sits slightly above the asymptote because `(1 − m/n)^n > e^−m` at
finite stratum size — the correct direction.)

### Claimed versus actual interpretation

The names `Bootstrap_Stability_Summary`, "BOOTSTRAP ENRICHMENT STABILITY" and
`Mean_Recovery_Rate` invite the reading *"87% of our top enriched GO terms are
robust"*. What the number actually reports is the probability that a row already
flagged significant appears at least once in a same-size resample with
replacement — a **bootstrap row-inclusion probability**, determined entirely by
how many comparisons each top term was already significant in.

So the labels `bootstrap stability`, `enrichment stability` and `robustness` do
**not** accurately describe the calculation. This is not a marginal
mischaracterisation: a term significant in exactly one comparison scores
0.632 regardless of how overwhelming its enrichment is, and a term significant
in four scores 0.98 regardless of how marginal each one is.

### Closed-form comparison

Every reported Monte Carlo quantity has a closed-form expectation. The exact
per-term expectation needs the 2026-05-27 input state, which no longer exists,
so the comparison is made self-contained from the workbook's own numbers.

| copy | reported rate | implied mean multiplicity M | SD (count) observed | homogeneous-P bound | observed range | predicted mean ± 2.5 SD |
|---|---|---|---|---|---|---|
| PRIDE-staged | 0.8738182 | 2.070 | 2.255 | 2.463 | [43, 53] | [42.4, 53.7] ✓ |
| superseded neuropil | 0.8846032 | 2.159 | 2.034 | 2.536 | [51, 60] | [50.6, 60.8] ✓ |
| superseded microglia | 0.7208 | 1.276 | 2.252 | 2.243 | [12, 23] | [12.4, 23.6] ✓ |
| superseded soma | 0.650625 | 1.052 | 1.995 | 1.907 | [5, 15] | [5.4, 15.4] ✓ |
| superseded bare | 0.6535714 | 1.060 | 1.888 | 1.780 | [3, 13] | [4.4, 13.9] ✗ (min off by 1) |

Every reported value inverts to an implied mean multiplicity between **1.05 and
2.16**, i.e. each is fully explained by "the top terms were significant in one
or two comparisons on average". The two lowest values (1.052, 1.060) are
essentially exactly M = 1, i.e. the bare `1 − e^−1 = 0.632` constant.

Two honest qualifications:

- The homogeneous-P variance bound is *exceeded* by 0.4–6% in three copies. That
  is expected and informative rather than contradictory: top terms share
  comparison strata, so an unlucky resample of one stratum causes several terms
  to be missed together. That positive within-stratum correlation inflates
  variance above the independent-Bernoulli bound. The Bernoulli-sum model is
  therefore a good first-order description, not an exact one.
- The range prediction holds for 4 of 5 copies. The miss is one unit on a count
  of 14, where the normal ±2.5σ approximation is crude because of discreteness.

### Monte Carlo discrepancy versus recorded precision

`se(Mean_Recovery_Rate) = SD_Recovery_Rate / sqrt(100)`:

| copy | value | se | decimal places that are signal | decimal places recorded |
|---|---|---|---|---|
| PRIDE-staged | 0.8738182 | 0.00410 | ~2 | 7 |
| superseded neuropil | 0.8846032 | 0.00323 | ~2 | 7 |
| superseded microglia | 0.7208000 | 0.00901 | ~2 | 4 |
| superseded soma | 0.6506250 | 0.01247 | ~1 | 6 |
| superseded bare | 0.6535714 | 0.01348 | ~1 | 7 |

Most recorded digits are Monte Carlo noise.

### The workbook contains no additional quantities

`Mean_Recovery_Rate`, `SD_Recovery_Rate`, `Min_Recovery` and `Max_Recovery` are
four summaries of one length-100 vector; `Total_TopTerms` is `length(top_terms)`
and is deterministic. So there is nothing else in the workbook that might carry
independent information — the brief's caution against overstating does not
change the conclusion here.

## C. Downstream dependency

Searching `08_Bootstrap_Stability_Summary` across the repository:

| scope | files scanned | references |
|---|---|---|
| manifests / indexes | 5 | **1** |
| code (`analysis/`, `R/`, `tools/`, `tests/`) | 340 | 6 |
| docs | 37 | **0** |
| manuscript repository | 770 | **0** |
| `Mean_Recovery_Rate` / `Recovery_Rate` in manuscript | 770 | **0** |

The one manifest reference is `pride_submission/manifests/pride_file_manifest.tsv`
line 1564. The six code references are the dead producer's own `write_xlsx` at
`compare_go_enrichment.R:2972` plus five assertions in the Phase 6H.8 audit test
that exist to pin this artifact.

| dependency type | count |
|---|---|
| manuscript citations | **0** |
| table / figure dependencies | **0** |
| textual claim dependencies | **0** |
| source-data dependencies | **0** |
| deposit-only references | **1** (rule-generated) |

One clarification so the zero is not misread. `manuscript_draft.md` *does*
discuss bootstrapping — "5,000 non-parametric bootstrap resamples over animals"
(:801), leave-one-animal-out resampling (:633), a bootstrap interval including
zero (:195). Every one of those is an **animal-level behavioural or abundance**
resample from a different, seeded analysis. None refers to GO term recovery, and
the manuscript cites no supplementary table by number at all. So the paper uses
bootstraps; it just never uses this one.

## D. Deposition role

**`DEPOSIT_UNDECLARED_EXTRA`.**

Physical presence and declared deposition are not equivalent here, and the
distinction is the whole point. The single manifest row is generated by a
path-matching rule, not by curation:

```r
# R/utilities/export_helpers.R:881
if (grepl("/pride_submission/", p)) return("pride_staging")
# :933
intended_for_PRIDE = cats %in% c(..., "pride_staging", ...)
```

and `analysis/publication_source_data/build_pride_manifest.R:36-40` rescans the
directory, so the file is re-endorsed on every manifest build purely for being
on disk. The manifest was last regenerated 2026-09-21 — four months after the
producing code became unreachable.

Curated references: **0**. Rule-based references: **1**. No submission document,
README, checklist, SDRF/IDF file, supplementary-table index or source-data
export names it.

## E. Contribution classification

**`E — MISLEADING_OR_UNINTERPRETABLE`.**

`D — REDUNDANT_COMPUTATIONAL_DIAGNOSTIC` also describes its *computational*
content accurately: the quantity has a closed form, so 100 Monte Carlo
iterations estimate something that could be written down exactly from
multiplicity counts, and those counts are already derivable from the deposited
enrichment tables.

`E` is chosen because it is the property that governs a **deposition** decision.
As it would be deposited — a 1×5 sheet named `Bootstrap_Stability_Summary` with
a column `Mean_Recovery_Rate`, no caption, no methods note, no provenance sheet
— a reader has nothing available to prevent the reading "87% of our enriched
terms are robust", which is not what the number means. Misleading in the
deposited form, and nothing accompanies it to correct that.

## F. The unreachable tail

| | |
|---|---|
| unreachable region | lines **572–4321** |
| lines | **3,750** |
| top-level expressions | **336 of 439** |
| exit boundary | line **569**, `quit(status = 0, save = "no")` |
| section banners in the tail | **17** |
| write calls in the tail | **56** |
| top-level names defined | **157** (20 functions) |

The tail's own section numbering restarts at 1, and `analysis_params$script` is
`"compareGO.r"` with `version = "2.1 (enhanced)"` — the tail is the body of the
original standalone `compareGO.r`, with a canonical reimplementation prepended.

Its 17 sections: summary statistics, term consistency, gene importance ranking,
comparison similarity, enriched-term barplot, NES ridge plot, similarity
heatmap, redundancy report, term-comparison occurrence matrix, top driver genes,
UpSet plot, **bootstrap enrichment stability**, parameter & reproducibility log,
data quality summary, Sankey diagram, alluvial diagram, term hierarchy.

Its 56 dead write calls: 16 `write_xlsx`, 14 `svg`, 13 `ggsave`,
10 `write_raw_xlsx`, 2 `writeLines`, 1 `write.csv`.

### Structural proof of unreachability

R evaluates top-level expressions in order and has no goto, so the tail runs
only if the `quit()` is evaluated and fails to exit. Each bypass was checked:

| check | result |
|---|---|
| `quit()` is its own top-level expression, not nested in `if`/`for`/`function` | **yes** (expression 103 of 439) |
| its arguments are constants (no conditional status) | **yes** |
| `quit` redefined anywhere in the file | **no** |
| `quit` redefined in any of the 8 libraries sourced before the exit | **no** |
| `trace()` / `untrace()` used in the file | **no** |
| any file `source()`s this script | **no** — all 8 `source(` matches are the script loading its own libraries |
| any file `parse()`s or `eval()`s it | **no** — 19 external references: 17 path/comment mentions, 2 `readLines` |
| tail-defined functions imported elsewhere | **no** |

On the last point: a name scan flagged `jaccard`, `mode_value` and
`optional_read_csv` as appearing in other files. These are **name collisions,
not dependencies** — `annotate_neuropil_reference.R:268` and
`summarize_biological_programs.R:264,269` each define their own copies, and since
nothing sources `compare_go_enrichment.R`, no tail definition can be imported by
anything.

That conclusion took three passes, and the reason is worth recording because the
same trap will recur in any dead-code audit. Name-based dependency detection
needs **token-type discrimination**, not name matching:

| criterion | result for `jaccard` | why it is wrong |
|---|---|---|
| name appears anywhere in the file | 4 files — looks like a dependency | matches comments, strings, column names |
| `SYMBOL` or `SYMBOL_SUB` tokens | 4 files — still looks like one | `jaccard = …` as a named argument, `.data$jaccard` as a column, `nature_palette("jaccard")` as a role |
| **`SYMBOL_FUNCTION_CALL` only, minus files with a local definition** | **0 files** | the only caller, `annotate_neuropil_reference.R`, defines its own at :268 |

Outside `compare_go_enrichment.R`, `jaccard` is a column name, a named argument
and a palette role — and a function call in none of them.

**Conclusion: the tail is unreachable by every ordinary entry path, and nothing
outside depends on its definitions.**

### Four tests depend on its literal presence

This is the real consequence of moving or deleting it, and it is the complete
bookkeeping list for C2:

1. **`tests/testthat/test-comparego-canonical-contract.R:120-128`** —
   "canonical compareGO path terminates before legacy UniProt logic". It does
   `marker <- grep("LEGACY_COMPAREGO_TAIL_DISABLED_BY_CANONICAL_EXIT", script)`
   then `expect_length(marker, 1L)`. Removing the tail removes the marker and
   fails this test.

2. **`tests/testthat/test-protein-group-enrichment-utils.R:53-57`** — "compareGO
   prefers manifest-provided collapsed gene inputs". It asserts the script text
   matches `comparison_input_file` and `GeneSymbol`. **Both tokens occur only in
   the tail** — `comparison_input_file` 3 times, `GeneSymbol` once, and **zero
   times** in the canonical head.

3. **`tests/testthat/test-ggrepel-render-determinism.R:108`** —
   `expect_gte(total, 12L)` over every `geom_text_repel` / `geom_label_repel`
   layer in `analysis/`, `R/` and `tools/`. The tree total is **exactly 12**,
   with no slack, and one of the twelve is
   `compare_go_enrichment.R:2452` — **inside the dead tail**. Archiving drops
   the count to 11 and fails the assertion.

   This also corrects a framing in Phase 6H.6, which recorded 12
   "publication-facing" ggrepel layers and seeded all 12. One of those twelve
   can never render. Seeding it was harmless and the other eleven are live, so
   no 6H.6 conclusion changes — but the live count is 11, not 12.

4. **`tests/testthat/test-statistical-rng-audit.R:155-183`** — the Phase 6H.8
   unreachability tripwire, which asserts exactly one `slice_sample` token in
   the file and that it sits after the exit. Archiving removes the token and
   fails both assertions. This one is *designed* to fail on C2; it exists so the
   R1 recommendation cannot be silently inherited after the code moves.

The second item is a latent test-validity defect worth stating plainly: a test
whose name claims to verify live compareGO behaviour currently passes by
matching text in code that cannot execute. It provides no assurance about the
canonical path. Moving the tail would expose that and force the test to be
rewritten against the head — which is an argument *for* C2, not against it.

### A weakness found in this phase's own guard

The Phase 6H.8 test that pins the PRIDE-staged workbook gated its whole block on
`skip_if_not(file.exists(f))`, so deleting the workbook would have made every
assertion below it — including the manifest-row guard — vanish into a **silent
skip** rather than a failure. A tripwire that disappears along with the thing it
watches is not a tripwire. Fixed in this phase: file presence is now asserted
rather than skipped, and the presence-versus-manifest-endorsement invariant is
checked with no dependence on the file existing, so it holds under either
decision.

### Why the tail exists

**Superseded legacy analysis.** Established from history, not from its position:

| date | commit | event |
|---|---|---|
| 2026-05-06 | `bb317dc` | "Enhance plotting, data loading and analyses" — introduces the bootstrap block; the tail is **live** |
| 2026-05-27 | — | PRIDE-staged workbook produced (tail live) |
| 2026-06-13 | — | the four superseded copies produced (tail live) |
| **2026-07-16** | **`6999f47`** | **"Enforce canonical enrichment provenance contracts" — introduces the canonical exit and the disable marker; the tail becomes dead** |
| 2026-09-17 | `5919819` | repository-wide rename touches the file |
| 2026-09-21 | — | PRIDE manifest regenerated, re-endorsing the now-orphaned file |

This corrects Phase 6H.8, which dated the marker to `5919819` (2026-09-17);
`git log -S` without `--follow` could not see past that rename commit. The
6H.8 record has been amended.

**No active script computes an enrichment-stability or term-recovery statistic
today.** The canonical path deliberately dropped the quantity rather than
reimplementing it — so this is not a capability the repository lost by accident.

## G. Code curation — C1 / C2 / C3

**Chosen: `C2 — ARCHIVE_UNREACHABLE_TAIL`.**

**C1 (keep in place)** is rejected. Its stated precondition is provenance value
that *cannot be preserved more cleanly elsewhere*, and this repository has an
established archive convention that preserves it more cleanly. Its listed risks
are not hypothetical here: two consecutive audit phases (6H.8 and 6H.9) were
spent discovering that a dead RNG call and 56 dead output writers were dead, and
a test is currently drawing false assurance from the region.

**C3 (delete)** is rejected. Nothing blocks it technically — git history holds
the content and no external code depends on it. But 3,750 lines comprising a
complete superseded analysis with 17 sections and 56 output writers is
substantial provenance, and repository policy explicitly retains this class of
content: `pipeline_analysis_script_exclusions()` declares `archive/` as
"Superseded and exploratory generations … provenance, not runnable stages", and
`archive/deprecated/clusterProfiler_newest_Sep16.r` is already exactly this kind
of artifact. Deleting when a documented retention convention exists would
discard provenance the repository's own policy says to keep. The brief also warns
against choosing C3 simply because the code is dead.

**C2** satisfies every criterion in section 21: the active script ends at its
real execution boundary (:569), the dead RNG and output writers leave the active
audit surface, content is preserved under the convention already used for
superseded generations, and nothing is regenerated. It also forces the two
literal-presence tests to be rewritten against the canonical head, converting a
hollow assurance into a real one.

Bookkeeping C2 requires (not performed):

1. relocate lines 572–4321 under `archive/`, leaving the active script ending at
   its real boundary;
2. rewrite `test-comparego-canonical-contract.R:120-128` to assert the head ends
   at the canonical exit without relying on the marker travelling with the tail;
3. rewrite `test-protein-group-enrichment-utils.R:53-57` against the canonical
   head or against `R/enrichment/protein_group_enrichment_utils.R` — this is the
   change that converts a hollow assurance into a real one;
4. lower `test-ggrepel-render-determinism.R:108` from `expect_gte(total, 12L)` to
   11, since one of the twelve seeded repel layers lives in the tail;
5. update the `test-statistical-rng-audit.R` unreachability tripwire, which by
   design fails once the tail is gone;
6. refresh `docs/active_script_io_audit.tsv` and any legacy registry entry.

Items 4 and 5 are why this phase adjudicates rather than implements: the change
is small in the source file and touches four test files, and each of those
touches deserves to be a deliberate edit rather than a hurried fix to a red
suite.

## H. Artifact curation — P1 / P2 / P3 / P4

**Chosen: `P2 — RETAIN_AS_INTERNAL_PROVENANCE_ONLY`.**

**P1 (keep in public deposit)** is rejected. Its precondition is that the
artifact conveys scientifically useful, interpretable information despite its
provenance limits. Section 7 asks whether its presence improves transparency,
reproducibility, reuse or interpretation. It improves none: there is nothing to
be transparent *about*, since no claim rests on it; it is not reproducible from
any active code path and its seed is unknown; it cannot be reused because it is
a single closed-form-reducible number with no metadata; and it actively harms
interpretation through its naming. P1 would additionally require documenting
stochastic generation, unknown seed and non-inferential status — three
disclaimers longer than the five numbers they qualify.

**P4 (replace with a deterministic analytical artifact)** is rejected. Its
precondition is a genuine reason the information should be public. There is
none: no manuscript claim uses it, and any reader wanting the multiplicity
structure can compute it directly from the deposited enrichment tables, exactly
and without Monte Carlo. Publishing a deterministic replacement would add a
public artifact nobody needs.

**P3 (remove from the deposition payload)** is close and was seriously
considered. It is rejected only because it preserves evidence "through Git/audit
record", and `pride_submission/` is gitignored — so P3's preservation route does
not actually apply to this file's bytes.

**P2** is the accurate label for the right action: take it out of the outward
deposition selection while keeping the bytes locally for provenance, alongside
the four sibling copies already under
`results/manuscript/_superseded_20260622/supplementary_tables/`. The provenance
limitations are documented here rather than shipped to a public repository.

Because `intended_for_PRIDE` is assigned by path, relocating the file out of
`pride_submission/` and rebuilding the manifest drops the row automatically — no
manual manifest edit, and no change to the classification rule that would affect
other files.

## I. Simulated consequences of the recommended actions

Computed read-only; nothing was written.

| quantity | now | after | delta |
|---|---|---|---|
| `pride_submission/` physical files | 1,310 | 1,309 | **−1** |
| supplementary tables in the validator count | 496 | 495 | **−1** |
| `pride_file_manifest.tsv` rows | 1,721 | 1,720 | **−1** |
| manuscript citations affected | 0 | 0 | 0 |
| figures / tables affected | 0 | 0 | 0 |
| freeze entries affected | 0 | 0 | 0 |

Files requiring a bookkeeping update:

| file | why |
|---|---|
| `pride_submission/manifests/pride_file_manifest.tsv` | row 1564 (the only row naming it; drops automatically on rebuild) |
| `pride_submission/validation/validation_report.tsv` | line 11 hard-records `supplementary_tables_present PASS 496 supplementary table(s)` → 495 |
| `pride_submission/validation/validation_summary.md` | line 28 records the same 496 count |
| `tests/testthat/test-statistical-rng-audit.R` | its Phase 6H.8 assertions pin this path, size and sha |

The two validator outputs were **missed by this audit's first pass**, which
checked manifests, indexes and the freeze but not
`pride_submission/validation/`. They are written by
`analysis/publication_source_data/validate_pride_submission.R`, so both
regenerate when the validator is rerun — but they do embed a hard count, and any
removal that skipped them would leave the package self-inconsistent.

Confirmed **not** affected — each contains zero references to the workbook, its
SHA-256 or its size: `results/manuscript/figure_export_manifest.csv`,
`figure_publication_audit.csv`, `source_data_export_manifest.csv`,
`docs/publication_freeze_manifest.yml`. The freeze records no
`pride_submission` file count or hash, and **no test asserts a
`pride_submission` file count (1,310) or manifest row count (1,721)** — the
literals 1310, 1721 and 1564 appear nowhere in `tests/` or `docs/`.

### One piece of evidence that does *not* support the decision

`pride_submission/supplementary_tables/_supplementary_table_staging_manifest.tsv`
(written by `build_supplementary_tables.R`) has 60 data rows and **no row for the
workbook** — which looks at first like proof that the curated staging process
never selected it. It is not. The directory holds 496 supplementary tables, so
roughly 436 other present files are equally absent from that manifest, which
evidently records one staging run rather than the whole directory. Its silence
about the workbook is therefore uninformative.

The load-bearing evidence for `DEPOSIT_UNDECLARED_EXTRA` remains the
path-matching rule in `export_helpers.R:881`/`:933` plus the directory rescan in
`build_pride_manifest.R:36-40` — not this absence.

### Reproducibility, asked in both directions

Section 17 asks two questions that can point opposite ways, and here they do
not.

*Does removing it hide information needed to reproduce a reported result?* **No.**
No reported result depends on it. Every input needed to recompute the quantity
exactly — the enrichment tables with their `p.adjust` columns — is itself
deposited, so a reader who wants it can derive it in closed form.

*Does keeping it falsely suggest the artifact is reproducible?* **Yes.** A
workbook in a deposition payload carries an implicit claim that it was produced
by the deposited analysis. This one was produced by code that has been
unreachable since 2026-07-16, from an RNG state that was never recorded, and
cannot be regenerated by anything in the repository. Keeping it creates
misleading provenance; removing it removes nothing a reader could use.

Both answers therefore favour P2.

## J. Integrity

| | |
|---|---|
| payload changed | **no** |
| package changed | **no** |
| freeze changed | **no** |
| `figure_export_manifest.csv` | 5,582 rows, `0fd0c9ed…` |
| `docs/publication_freeze_manifest.yml` | `b4d37250…` |
| `exports/` | 55 files, 2,545,817 bytes |
| `pride_submission/` | 1,310 files, manifest 1,721 rows |
| the workbook | still in place, `024d3671…`, 5,141 bytes |
| the code tail | still in place, lines 572–4321 |

Neither decision was implemented. `quit()` was not removed, `slice_sample()` was
not seeded, the bootstrap was not modernised, rerun or regenerated, and nothing
was removed from `pride_submission/`.

### A transient integrity breach during this phase, and what caught it

Two `grep.exe.stackdump` files — msys crash dumps produced by this phase's own
search tooling, not repository content — were written at 2026-09-22 11:42 to the
repository root and to
`pride_submission/supplementary_tables/`. The second took the payload from 1,310
to **1,311 files**. Both were verified to be stack traces
(`msys-2.0.dll+0x215D7` frames) and removed; the count is restored to 1,310 with
zero stackdumps anywhere in the tree.

Two things are worth recording rather than quietly fixing:

1. **It was caught by a pre-existing repository test, not by this phase's
   integrity checks.** `test-supplementary-filename-budget.R:184` flagged
   `grep.exe.stackdump` as a malformed staged filename. That test — added in
   Phase 6H.4/6H.5 for an entirely different reason — did exactly its job.
2. **This phase's own payload verification ran too early to see it.** The
   `pride_submission/ files = 1310` check was performed during the removal
   simulation, *before* the crash occurred. An integrity check that runs before
   the work is not an integrity check. The final verification re-measures rather
   than recalling the earlier value, which is why the count above is a fresh
   measurement.

No audit conclusion depended on the payload file count, and the two decisions
rest on the reference graph, the unreachability proof and the closed-form
analysis, none of which touched that directory.

## Scope notes

Untouched, per section 18: the addressability vocabularies in
`integration_utils.R`, `evidence_bundle_utils.R` and `compare_go_enrichment.R`;
the `biological_claims_table` schema; archive rendering RNG; repository
relocation.

One pre-existing record was corrected: the marker-introduction date in
`phase6h_statistical_rng_reproducibility.md`, which Phase 6H.8 attributed to the
wrong commit.
