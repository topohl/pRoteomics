# Manuscript statistical contract

One row per analysis that reaches the manuscript. This file answers, once:
what exactly was tested, on how many animals, against which multiple-testing
family, and what may therefore be claimed.

Generated from the frozen v9 figure contract and the Part A audit tables by
`audits/publication_hardening/02_statistical_audits.R` and rendered by
`audits/publication_hardening/03_statistical_contract_doc.R`.

## 1. Analyses

### preprocessing_bilateral

- **Question:** How are hemispheres combined into one animal-level value?
- **Biological n:** 9 animals (unit of analysis: protein group x animal x spatial unit)
- **Model / test:** bilateral aggregation to animal level
- **Design:** 3 groups x 18 spatial units; contrast: not applicable
- **Repeated structure:** two hemispheres per animal, aggregated
- **Multiple-testing scope:** not applicable
- **Software:** in-house R (R 4.5.1)
- **Claim status:** descriptive
- **Major limitation:** hemispheres are repeated tissue, not independent replicates
- **Methods sentence:** Hemispheric measurements were aggregated to one value per animal and spatial unit before any group comparison, so that the biological replicate is the animal.

### differential_abundance

- **Question:** Which proteins differ between outcome groups within a spatial unit?
- **Biological n:** 3 per group (unit of analysis: protein group)
- **Model / test:** moderated linear model (Protigy, animal_level_protigy_da_v1)
- **Design:** 3 groups x 18 spatial units; contrast: SUS-RES, SUS-CON, RES-CON
- **Repeated structure:** one value per animal; no within-animal repeats at this stage
- **Multiple-testing scope:** BH within each comparison
- **Software:** Protigy (animal_level_protigy_da_v1)
- **Claim status:** inferential
- **Major limitation:** n = 3 per group limits power; absence is not evidence of absence
- **Methods sentence:** Protein-level differential abundance was assessed with a moderated linear model on animal-level values within each spatial unit, with Benjamini-Hochberg correction within each comparison.

### gsea

- **Question:** Which biological programs are coordinately shifted?
- **Biological n:** 3 per group (unit of analysis: gene (median moderated t per official symbol))
- **Model / test:** clusterProfiler::gseGO over GO-BP
- **Design:** ranked list per comparison; contrast: SUS-RES, SUS-CON, RES-CON
- **Repeated structure:** collapsed to one statistic per gene
- **Multiple-testing scope:** BH over every GO-BP set returned in that comparison
- **Software:** clusterProfiler / fgsea (clusterProfiler 4.18.4)
- **Claim status:** inferential
- **Major limitation:** p-values are floored at eps = 1e-10; a term at the floor has an unresolved true p
- **Methods sentence:** Genes were ranked by the median moderated t statistic per official gene symbol and tested against GO biological process with clusterProfiler::gseGO (minGSSize 10, maxGSSize 800, BH within each comparison, eps = 1e-10).

### camera_sensitivity

- **Question:** Does the program-level direction survive a correlation-aware competitive test?
- **Biological n:** 3 per group (unit of analysis: gene)
- **Model / test:** limma::cameraPR
- **Design:** preranked, inter.gene.cor = 0.01; contrast: all three pairwise contrasts
- **Repeated structure:** none
- **Multiple-testing scope:** BH over the full comparable GO-BP family within each contrast
- **Software:** limma (limma 3.66.0)
- **Claim status:** sensitivity only (LEVEL 4)
- **Major limitation:** preranked CAMERA is weaker than a full expression-matrix CAMERA
- **Methods sentence:** As a sensitivity analysis, the same ranked statistics were tested with limma::cameraPR (inter-gene correlation fixed at 0.01) over the comparable GO-BP family; this is a concordance check and not independent validation.

### go_program_atlas

- **Question:** How are enrichment results summarised for display?
- **Biological n:** 3 per group (unit of analysis: GO term grouped into ontology-defined families)
- **Model / test:** median NES per family
- **Design:** 7 primary families x 18 spatial units; contrast: all three pairwise contrasts
- **Repeated structure:** none
- **Multiple-testing scope:** NONE - descriptive aggregation
- **Software:** in-house R + GO.db (R 4.5.1)
- **Claim status:** descriptive
- **Major limitation:** no theme-level p-value or FDR exists or is implied
- **Methods sentence:** Canonical GO-BP results were mapped to seven ontology-defined program families (registry manuscript_go_themes_v3) and summarised as the median NES of their constituent terms; this summary is descriptive and carries no theme-level statistical test.

### wgcna

- **Question:** Are co-abundance modules associated with outcome?
- **Biological n:** 9 animals (unit of analysis: module eigengene)
- **Model / test:** lmerTest eigengene ~ StressGroup + SpatialUnit + (1|AnimalID)
- **Design:** modules x spatial units; contrast: all three pairwise contrasts
- **Repeated structure:** repeated spatial units within animal, modelled as a random intercept
- **Multiple-testing scope:** BH within the module-phenotype family
- **Software:** WGCNA / lmerTest (WGCNA 1.74)
- **Claim status:** descriptive - no association survived correction
- **Major limitation:** 0 of 45 module-phenotype tests survived correction
- **Methods sentence:** Module eigengenes were related to outcome group with a linear mixed model including a random intercept for animal; no module-phenotype association survived Benjamini-Hochberg correction.

### ewce

- **Question:** Which external cell types are enriched among module members?
- **Biological n:** not applicable - external reference (unit of analysis: gene set)
- **Model / test:** EWCE bootstrap enrichment
- **Design:** module vs reference panel; contrast: not applicable
- **Repeated structure:** none
- **Multiple-testing scope:** BH within the EWCE family
- **Software:** EWCE (see manifest)
- **Claim status:** external context only
- **Major limitation:** affinity is context, never cell-intrinsic identity
- **Methods sentence:** Module cell-type affinity was assessed against external single-cell reference panels; these results provide cell-type context and do not establish the cellular origin of an enriched-ROI measurement.

### spatial_identity

- **Question:** Is the spatial molecular architecture reproducible?
- **Biological n:** 3 CON animals (unit of analysis: spatial unit profile)
- **Model / test:** bilateral correlation and ICC
- **Design:** CON only; contrast: not a group contrast
- **Repeated structure:** two hemispheres per animal
- **Multiple-testing scope:** none - descriptive
- **Software:** in-house R (R 4.5.1)
- **Claim status:** descriptive
- **Major limitation:** CON only; says nothing about stress
- **Methods sentence:** Spatial reproducibility was quantified as the left-right concordance of each prespecified anatomical contrast in control animals.

### external_validation

- **Question:** Do internal anatomical contrasts recover published signatures?
- **Biological n:** 3 CON animals (unit of analysis: gene set)
- **Model / test:** GSEA against external signatures
- **Design:** CON-only anatomical contrasts; contrast: anatomical, not phenotypic
- **Repeated structure:** none
- **Multiple-testing scope:** BH within the external-validation inventory
- **Software:** clusterProfiler (clusterProfiler 4.18.4)
- **Claim status:** inferential - the only external validation
- **Major limitation:** validates anatomy, not the stress result
- **Methods sentence:** Internal control-only anatomical contrasts were tested against independently published hippocampal signatures; this is the only externally anchored validation in the study.

### network

- **Question:** Does the spatial molecular network differ by outcome?
- **Biological n:** 9 animals (27 animal x dataset instances) (unit of analysis: animal x dataset network)
- **Model / test:** distance from a leave-one-CON-animal-out consensus
- **Design:** 3 groups x 3 datasets; contrast: SUS/RES vs CON consensus
- **Repeated structure:** three dataset instances per animal
- **Multiple-testing scope:** BH within the network family
- **Software:** in-house R (R 4.5.1)
- **Claim status:** descriptive - no detectable difference
- **Major limitation:** 27 instances arise from 9 animals; they are not 27 independent replicates
- **Methods sentence:** Animal-level spatial networks were compared with a leave-one-control-animal-out consensus; no whole-network group difference and no edge-behaviour association survived correction at this sample size.

### behaviour_correlation

- **Question:** Do network edges track behavioural outcome?
- **Biological n:** 9 animals (unit of analysis: network edge)
- **Model / test:** correlation with behavioural score
- **Design:** 8 tested neuropil edges; contrast: edge x behaviour
- **Repeated structure:** one value per animal
- **Multiple-testing scope:** BH within the edge-coupling family
- **Software:** in-house R (R 4.5.1)
- **Claim status:** descriptive - none survived correction
- **Major limitation:** only 8 neuropil edges were tested, not every spatial-unit pair
- **Methods sentence:** Edge-behaviour associations were tested for eight neuropil spatial-unit pairs; none survived correction, and with nine animals a single correlation has very little resolution.

## 2. Effect and sign conventions

| Analysis | Contrast | Formal effect definition | Positive means |
|---|---|---|---|
| differential abundance | SUS - RES | group2 - group1 on animal-level log2 abundance within a spatial unit; Protigy contract animal_level_protigy_da_v1 | higher abundance in susceptible animals |
| differential abundance | SUS - CON | group3 - group1 | higher abundance in susceptible than control |
| differential abundance | RES - CON | group2 - group1 | higher abundance in resilient than control |
| ranked GSEA | RES - CON / SUS - CON / SUS - RES | gseGO on genes ranked by the median moderated t per official gene symbol | coordinated higher abundance of the gene set in the first-named group |
| GO-program atlas | RES - CON / SUS - CON / SUS - RES | median of constituent canonical GO-term NES per dataset x spatial unit x theme | coordinated higher abundance of the program in the first-named group |
| CAMERA sensitivity | RES - CON / SUS - CON / SUS - RES | cameraPR on the exact canonical ranked statistic, inter.gene.cor = 0.01 | gene set shifted upward relative to the rest of the ranking |
| WGCNA module phenotype | RES - CON / SUS - CON / SUS - RES | lmerTest eigengene ~ StressGroup + SpatialUnit + (1/AnimalID), named contrast estimate | higher module eigengene in the first-named group |
| spatial network | SUS/RES vs CON consensus | Euclidean distance from a LEAVE-ONE-CON-ANIMAL-OUT consensus | greater divergence from the control network |
| bilateral reproducibility | not a group contrast | correlation between hemispheres of the same animal | higher left-right agreement |
| external anatomical validation | CON-only anatomical contrasts | internal anatomical contrast tested against an external signature | agreement with the published anatomical signature |

## 3. Statistic identity

Each displayed quantity has exactly one correct name and a list of names that
must never be substituted for it.

| Quantity | Correct name | Never call it | Where displayed |
|---|---|---|---|
| protein differential abundance | log2 fold change | expression; effect size; abundance change significance | F3 g/h/i leading-edge dot plots |
| gene ranking statistic | moderated t | log2FC; fold change | GSEA input (not displayed) |
| gene-set enrichment | normalised enrichment score (NES) | significance; effect size | F3 d/e/f strips, ED6 c/d/e |
| theme aggregation | median NES across a theme's GO terms | mean NES; theme significance; theme FDR | F3b, ED6 a/b |
| module member abundance | mean module-member CON z-score | expression; eigengene | ED WGCNA a |
| module eigengene effect | module eigengene difference | expression; abundance | ED WGCNA b |
| profile similarity | median Spearman profile correlation (CON) | connectivity; anatomical connection | ED8 a |
| bilateral precision | intraclass correlation (ICC) | reliability significance | ED1 b |
| cell-type affinity | EWCE z-score / reference-panel overlap | cell-type identity; purity | ED WGCNA c |
| baseline abundance | CON z-score | expression | F2 d, F2 e, ED2 a |

## 4. Dataset-specific interpretation rules

| Dataset | Resolution | Legal | Illegal |
|---|---|---|---|
| neuron_neuropil | region x layer (10 units) | layer-resolved and laminar wording; region x layer comparisons | cell-type attribution of a neuropil measurement |
| neuron_soma | region only (4 units) | region-level statements; neuronal-soma compartment | any layer or laminar statement about this compartment |
| microglia_enriched | region only (4 units) | microglia-enriched ROI; local microenvironment | layer or laminar wording; microglial proteome; cell-intrinsic or cell-autonomous |

## 5. Standing constraints

- The biological replicate is the **animal**: n = 3 per group. Acquisitions,
  spatial units, hemispheres and animal x dataset network instances are never
  biological replicates, and are labelled as such wherever they are counted.
- **No theme-level p-value or FDR exists.** Atlas themes are descriptive
  aggregations of canonical GO terms and carry no multiple-testing family.
- **GSEA p-values are floored at eps = 1e-10.** This affects 90 of the 851
  displayed FDR-supported occurrences, including all three exemplars. A term
  at the floor has a true p the method does not resolve, so its FDR bounds
  the evidence rather than measuring it.
- **Absence of FDR support is never absence of effect.** At three animals per
  group, write 'did not survive correction' or 'no detectable ... at the
  present sample size', never 'no difference', 'unchanged' or 'equivalent'.
- **CAMERA is LEVEL 4 sensitivity**, never independent validation. The only
  genuinely external validation is the CON-only anatomical signature test
  (F2g, ED2b).
- **Specificity and selectivity claims require a heterogeneity test.** The
  only one executed is the WGCNA stress x spatial-unit omnibus (0 of 35
  FDR-supported, smallest FDR 0.27). No equivalent test exists at GO-program
  level, so program-level results are reported as the count of spatial units
  in which the contrast was FDR-supported, not as selectivity.

## 6. Spatial wording contract

One phrase, used consistently, for what the spatial design delivers.

- **`spatially resolved` is the default descriptive wording.** It says that
  the design and analysis resolve effects across defined hippocampal spatial
  contexts - 18 prespecified spatial units, laminar in the neuropil and
  region-level in the neuronal soma and microglia-enriched ROI. It asserts
  nothing about where effects are or are not present.
- **`restricted`, `selective` and `specific` require either a formal
  corresponding test or an explicitly factual description of observed
  support.** 'FDR-supported in 6 of 18 spatial units' is factual and allowed.
  'Spatially restricted program differences' is an inferential claim and is
  not, because the only heterogeneity test in the package is at WGCNA level
  and is FDR-negative.
- **Absence of support in some contexts does not by itself establish
  specificity.** At three animals per group, a unit without FDR support is a
  unit where nothing was detected, not a unit where nothing is happening.
  What is restricted is the detection, not the effect.

Enforced by the `spatially restricted / specific` rule in the S9 scan of
`figures/final_truth_v9_semantics.R`, which is anchored to the adverb
`spatially` so that a restricted *interpretation* or a specificity
*inventory* is not flagged.
