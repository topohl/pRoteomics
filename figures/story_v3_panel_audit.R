#!/usr/bin/env Rscript

# Panel-by-panel visual/scientific audit of the whole candidate library:
# canonical Figure 2/3 panels, the Part-16 candidate layer and the Part-17
# Nature-v2 layer, treated as ONE pool of ideas rather than as generations
# where the newest wins.
#
# Every row is a judgement recorded after inspecting the rendered figure at
# print size. The judgements are versioned here in the script so they can be
# reviewed and argued with; the CSV is the emitted artifact.
#
# This script renders nothing and promotes nothing.

source(file.path("R", "paths.R"))
source(repo_path("R", "integration_utils.R"))
suppressPackageStartupMessages({ library(readr) })
Sys.setenv(PROTEOMICS_SCRIPT_ID = "figures/story_v3_panel_audit.R")

r <- function(panel, source_version, scientific_question, two_second_takeaway,
              current_visual_encoding, encoding_works, problem, better_encoding,
              narrative_role, dataset_scope, recommended_destination) {
  data.frame(panel, source_version, scientific_question, two_second_takeaway,
             current_visual_encoding, encoding_works, problem, better_encoding,
             narrative_role, dataset_scope, recommended_destination,
             stringsAsFactors = FALSE)
}

audit <- rbind(

# ---------------------------------------------------------------- FIGURE 2
r("2b depth", "canonical", "How deep is the proteome per compartment?",
  "We identify thousands of proteins in every compartment.",
  "stacked ranked area bars per compartment", "partly",
  "Opens the figure with a technical QC statement and takes a full cell.",
  "small boxplot, one per compartment, <=40 mm wide",
  "technical context", "all three", "main_small_or_extended_data"),

r("2c PCA", "canonical", "Do measured proteomes occupy structured molecular spaces?",
  "Compartments and regions separate; the measurement captures real structure.",
  "sample-level PCA scatter with many point labels", "yes",
  "Over-labelled; dozens of sample labels compete with the geometry.",
  "same scatter, compartment-coloured, region-shaped, at most 4 anchor labels",
  "global structure", "all three", "main_small"),

r("2d compartment identity", "canonical",
  "Do soma, neuropil and microglia-enriched ROI show the expected markers?",
  "Intended markers are higher in the compartment they should mark.",
  "intended-minus-comparator dots grouped by marker family", "yes",
  "Dot cloud is fine but the grouping headers are large relative to the data.",
  "keep the logic; compact boxplot or effect forest by marker family",
  "compartment validation", "all three", "main"),

r("2e Kaulich external validation", "canonical",
  "Does an independent external atlas recover the same spatial identity?",
  "Our internal spatial contrasts match external signatures.",
  "purple sequential dot plot, signature x contrast", "partly",
  "Sequential purple reads as magnitude-only; direction and support are unclear.",
  "diverging NES tile with an explicit FDR keyline, one shared scale",
  "external validation", "neuropil + soma", "main"),

r("2f internal anatomical GO", "canonical",
  "Do internal anatomical contrasts recover known regional programs?",
  "Known anatomical programs come back from our own data.",
  "grouped GO dot plot across region and CA1 layers", "partly",
  "Large and sparse; many rows carry no supported term.",
  "top supported terms only, diverging NES tiles",
  "internal validation", "neuropil", "main_small_or_extended_data"),

r("2x_bilateral", "part16",
  "Does the spatial architecture reproduce across paired hemispheres?",
  "Left and right give the same anatomical effect.",
  "three L-vs-R scatterplots with four annotation lines each", "yes",
  "Four stacked statistics per facet is more text than the point needs.",
  "same scatter, identity line, ONE statistic, shared axis range",
  "reproducibility", "all three", "main_large"),

r("2x_precision", "part16",
  "How much reliability does bilateral averaging add?",
  "Averaging hemispheres raises reliability.",
  "per-endpoint lines plus median crossbar, faceted 3x3", "partly",
  "Nine facets for one message; the improvement is not immediate.",
  "dumbbell: single-side to bilateral-mean, one row per endpoint class",
  "measurement quality", "all three", "main_small_or_extended_data"),

r("2x_wgcna_bilateral", "part16",
  "Do modules keep level and spatial pattern across hemispheres?",
  "Module structure is bilaterally reproducible.",
  "scatter of absolute r vs spatial-profile r", "partly",
  "Two reproducibility axes at once is abstract for a main panel.",
  "fold into the WGCNA structure panel as one annotation track",
  "reproducibility detail", "all three", "extended_data"),

r("2x_compartment_bilateral", "part16",
  "Is compartment identity supported by markers AND by hemisphere concordance?",
  "Compartment identity holds and reproduces.",
  "marker boxplots stacked over a concordance bar chart", "no",
  "Two questions in one panel; neither reads in two seconds.",
  "split: keep marker panel in main, move concordance to the bilateral panel",
  "compartment validation", "all three", "extended_data"),

r("n2_anchor", "part17", "What was sampled, where, and what is the replicate?",
  "18 spatial units across 4 subfields and 3 compartments, bilateral, 9 animals.",
  "data-backed sampling map (region x layer/compartment tiles)", "yes",
  "Honest and legible, but it is a matrix where a reader expects anatomy.",
  "keep the map as the specification; pair it with a labelled hippocampal schematic",
  "experimental anchor", "all three", "main_large"),

r("n2_depth", "part17", "How deep is the proteome per compartment?",
  "Thousands of proteins per sample in all three compartments.",
  "compact boxplot per compartment", "yes",
  "None; this is the right size for the message.",
  "keep", "technical context", "all three", "main_small"),

r("n2_bilateral", "part17", "Does the architecture reproduce across hemispheres?",
  "Left effect equals right effect.",
  "three L-vs-R scatters, identity line, r only", "yes",
  "None; this is the cleanest form of the panel so far.",
  "keep", "reproducibility", "all three", "main_large"),

r("n2_precision", "part17", "Does bilateral averaging improve reliability?",
  "Reliability goes up.",
  "two-column line plot with median crossbars", "partly",
  "Reads as a generic before/after; endpoint classes are invisible.",
  "dumbbell by endpoint class with faint individual endpoint pairs",
  "measurement quality", "all three", "main_small"),

r("n2_compartment", "part17", "Do compartments carry their expected markers?",
  "Intended markers are enriched where they should be.",
  "horizontal marker-family boxplots coloured by compartment", "yes",
  "None.", "keep", "compartment validation", "all three", "main"),

r("n2_external", "part17", "Does an external atlas agree with our spatial identity?",
  "External signatures match our contrasts.",
  "diverging NES tiles with FDR dots", "yes",
  "Axis labels are long; otherwise good.",
  "keep, shorten signature labels", "external validation", "neuropil + soma", "main"),

r("n2_internal", "part17", "Do internal contrasts recover known anatomical programs?",
  "Known regional programs are recovered internally.",
  "top-8 supported GO terms as diverging tiles", "yes",
  "Sparse when given full width.", "keep but size to content",
  "internal validation", "neuropil", "main"),

# ---------------------------------------------------------------- FIGURE 3
r("3a DAP counts", "canonical",
  "How many proteins differ between SUS and RES, and where?",
  "One spatial unit has almost all the hits.",
  "signed horizontal count bars per spatial unit", "no",
  "One giant CA2-SLM bar implies a single dominant finding; 22 of those 28 are not claimable.",
  "compact spatial-unit overview with quiet robustness shading",
  "sets up sparsity", "all three", "main_small"),

r("3b circular supermodule atlas", "canonical",
  "How are module programs organised across space?",
  "There is a whole coherent system here.",
  "multi-ring circular heatmap with rotated labels", "partly",
  "Unreadable at 183 mm and its phenotype interpretation no longer holds.",
  "SIMPLIFIED circle for module STRUCTURE only: identity, spatial peak, cell type, reproducibility",
  "systems context", "neuropil", "main_or_extended_data_test_both"),

r("3c supermodule GO heatmap", "canonical",
  "What GO biology defines each supermodule?",
  "Modules have distinct GO identities.",
  "sparse diagonal GO x supermodule dot matrix", "partly",
  "Large, sparse, and subsumed by the ranked-GSEA theme atlas.",
  "drop from main; the theme atlas covers it",
  "module annotation", "neuropil", "extended_data"),

r("3d WGCNA 15-module heatmap", "canonical",
  "Do module scores differ between groups?",
  "Some modules look different between groups.",
  "module x contrast diverging heatmap", "no",
  "Zero of 45 cells reach tier-specific FDR (min 0.245); the colour field implies inference that does not exist.",
  "descriptive strip in Extended Data with the status stated",
  "phenotype claim that is not supported", "neuropil", "extended_data"),

r("3e m12 protein heatmap", "canonical",
  "Which individual proteins carry the signal?",
  "Real molecules, real spatial pattern.",
  "protein x spatial-unit log2FC heatmap", "yes",
  "The DESIGN is good; the biology was m12/CA2-SLM which is QC-compromised.",
  "reuse the design with QC-clean leading-edge proteins from supported programs",
  "molecular tangibility", "neuropil", "main_redesigned_biology"),

r("3x_dap_status", "part16", "How many hits, and how many are claimable?",
  "Most CA2-SLM hits are not claimable.",
  "stacked status bars per spatial unit", "partly",
  "Correct, but the QC decomposition is a second story inside panel a.",
  "quiet status shading in main; full 6/10/12 decomposition to Extended Data",
  "sparsity plus honesty", "all three", "main_small"),

r("3x_gsea_atlas_susres", "part16",
  "Where does coordinated program signal sit across compartments?",
  "Programs differ even where single proteins do not.",
  "theme x spatial-unit NES heatmap with FDR dots", "yes",
  "Long theme labels eat width; the reference-contrast strip doubles the height.",
  "keep as centrepiece; short theme labels; drop the reference strip to Extended Data",
  "quantitative centrepiece", "all three", "main_largest"),

r("3x_gsea_atlas_blocks", "part16", "Same as the atlas, three contrast blocks.",
  "Contrast geometry is comparable.",
  "three aligned contrast blocks", "partly",
  "Duplicates the atlas; only one layout can ship.",
  "keep the SUS-RES layout in main, blocks to Extended Data",
  "alternative layout", "all three", "extended_data"),

r("3x_wgcna_annotated", "part16", "Do module effects align with module identity?",
  "Modules differ by group and have annotation.",
  "module x contrast heatmap plus annotation tracks", "no",
  "Same unsupported phenotype claim as canonical 3d, now with more ink.",
  "keep only the annotation tracks; drop the phenotype colour field from main",
  "unsupported phenotype", "neuropil", "extended_data"),

r("3x_convergence", "part16", "Which programs recur across methods?",
  "Some programs appear in several analyses.",
  "green evidence table, 5 columns", "no",
  "An audit table, not a biological panel; invites reading five columns as independent validations.",
  "compact per-compartment support summary, or Extended Data",
  "internal synthesis", "cross-dataset", "extended_data"),

r("3x_program_examples", "part16", "What do the strongest programs look like per unit?",
  "A few programs are supported in several places.",
  "median-NES bars per spatial unit, 3 programs", "no",
  "Another summary of the atlas; adds no new evidence type.",
  "replace with DIRECT ranked evidence per compartment",
  "should be direct evidence", "all three", "replace"),

r("3x_stress_identity", "part16",
  "Do stress effects sit outside a protein's dominant baseline niche?",
  "Effects avoid the protein's own peak unit.",
  "subset bar chart plus rank histogram", "partly",
  "Introduces a second conceptual question mid-figure.",
  "keep the result, move to Extended Data unless the story needs it",
  "second question", "all three", "extended_data"),

r("n3_da", "part17", "How sparse are protein-level effects and where?",
  "Effects are sparse and spatially restricted.",
  "horizontal stacked status bars, faceted by compartment", "yes",
  "Good; status could be quieter still.",
  "keep, mute the status palette", "sets up sparsity", "all three", "main_small"),

r("n3_atlas", "part17", "Where does coordinated program signal sit?",
  "Programs differ across compartments and space.",
  "theme x unit NES heatmap, compartment-blocked, FDR dots", "yes",
  "None; this is the strongest quantitative panel in the library.",
  "keep as the largest panel", "quantitative centrepiece", "all three", "main_largest"),

r("n3_ex_neuropil / soma / microglia", "part17",
  "What does one supported program look like in the ranked proteome?",
  "The program's genes concentrate at one end of the ranking.",
  "ranked statistic curve with leading-edge rug", "yes",
  "Rug sits at the axis and is easy to miss; annotation crowds the corner.",
  "emphasise leading-edge marks, move NES/FDR to a clean corner block",
  "direct evidence triptych", "one per compartment", "main"),

r("n3_identity", "part17", "Where does the effect sit in a protein's own ranking?",
  "Effects sit away from the protein's peak unit.",
  "relative-rank dot strip by compartment", "partly",
  "Sparse at 26 mm and still a second question.",
  "Extended Data unless the narrative uses it",
  "second question", "all three", "extended_data"),

r("n3_synthesis", "part17", "Which programs recur across compartments?",
  "A few programs recur everywhere.",
  "horizontal bars split by compartment", "partly",
  "Restates the atlas rather than advancing the argument.",
  "Extended Data, or fold into the atlas margin",
  "synthesis", "all three", "extended_data"),

r("n3_wgcna_small", "part17", "Do module effects show a consistent direction?",
  "Module effects are small and unsupported.",
  "module x contrast diverging strip, no significance marks", "no",
  "Honest but it spends main-figure area on a null.",
  "replace with WGCNA STRUCTURE, or drop from main",
  "unsupported phenotype", "neuropil", "extended_data"),

r("ned_modules", "part17", "How spatially specific and reproducible is each module?",
  "Modules are spatially specific and bilaterally reproducible.",
  "tau lollipop coloured by reproducibility, faceted", "yes",
  "Y labels cramped at 88 mm.", "keep in Extended Data, or feed a main structure panel",
  "module structure", "all three", "extended_data_or_main_structure"),

r("ned_celltype", "part17", "Which cell types are modules associated with?",
  "Modules map onto known cell types.",
  "compartment x cell-type count tiles", "yes",
  "None.", "keep", "module annotation", "all three", "extended_data"),

r("ned_wgcna_heatmap", "part17", "What are module-level group effects?",
  "Descriptive only; nothing survives FDR.",
  "module x contrast heatmap", "partly",
  "Correct home, but wastes width at 175 mm for 3 columns.",
  "narrow it and pair with the annotation tracks",
  "descriptive phenotype", "neuropil", "extended_data")
)

out <- path_results("tables", "manuscript_candidates", "story_v3")
dir_create(out)
write_csv_safe(audit, file.path(out, "panel_visual_scientific_audit.csv"))

cat("\n===== Panel visual/scientific audit =====\n")
cat("panels audited:", nrow(audit), "\n\n")
print(as.data.frame(table(audit$source_version)), row.names = FALSE)
cat("\nrecommended destination:\n")
print(as.data.frame(table(audit$recommended_destination)), row.names = FALSE)
cat("\nencoding works?\n")
print(as.data.frame(table(audit$encoding_works)), row.names = FALSE)
cat("\nwritten:", relative_to(file.path(out, "panel_visual_scientific_audit.csv")), "\n")
