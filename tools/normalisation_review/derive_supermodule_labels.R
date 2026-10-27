#!/usr/bin/env Rscript
# Supermodule labels generated from their members' labels (normalisation review, 2026-10-23), so that every supermodule
# label follows one stated template and can be checked from the module labels alone:
#   member part  the member's own label without its suffix ('-enriched', '-associated'); a two-signal member contributes
#                'X and Y'; a several-signal member its first two named signals, 'X and Y'; '/' only ever joins terms of
#                one signal. The first letter is lower case unless the first word is an abbreviation (ER, RNA) or a name
#                (Golgi).
#   singleton    the member's label with its ID: 'X-enriched (m05)'; weak: 'unannotated (weak X association; m05)';
#                unannotated: 'unannotated (m05)'
#   concordant   the largest member's wording with all IDs, its suffix set by the supermodule's confidence:
#                'X-associated (m09 + m16)' (or '-enriched' when high)
#   compound     'Compound: X (mA) + Y (mB)', members in ID order
#   mixed        the labelled members' parts, then '+ unannotated member (mC)': 'Mixed: X (mA) + unannotated member (mC)'
#   unannotated  several unlabelled members: 'unannotated (weak X (mA) and Y (mB) associations)'; plain unannotated
#                members are listed by ID after a semicolon; without any weak member 'unannotated (mA, mB)'
# Also applies the module corrections listed in MODULE_FIXES, then those of REGENERATION_FIXES (both networks rebuilt
# with the official gene mapping and GO term sizes 10-800, 2026-10-24/26). Rewrites annotation_review/rule_based/labels.tsv and the
# change log annotation_review/label_changes_2026-10-23.csv (all changes against the labels preserved in the review
# archive of 2026-10-22, with a short reason each; supermodules are matched by their members, since the cut-height rule of
# 2026-10-26 merged and renumbered some). Run check_rule_labels.R afterwards.
# USAGE  Rscript tools/normalisation_review/derive_supermodule_labels.R      (working directory: repo root)

suppressPackageStartupMessages({ library(dplyr); library(readr) })
REV <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21"); ANN <- file.path(REV, "annotation_review"); RB <- file.path(ANN, "rule_based")
PREV <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/_archive_proteomics_review_2026-10-22/results/reviewer_audit/normalisation_choice_2026-10-21/annotation_review/rule_based/labels.tsv"
STYLE <- "Spelling unified across labels: GO term words without added hyphens, and GO's spelling 'signaling'."
MODULE_FIXES <- tribble(~version, ~dataset, ~entity_id, ~label, ~summary,
  "new", "neuron_soma", "WGCNA_m08", "Synapse organization/actin cytoskeleton organization-enriched",
  "The signal's leading and largest term is synapse organization (19 proteins); actin cytoskeleton organization (17) is tied with it on q in the same signal. The previous label named other tied terms.",
  "new", "neuron_soma", "WGCNA_m05", "Actin cytoskeleton organization-enriched", NA,   # carried from original soma m05 (rule 8): keeps its words
  "original", "microglia", "WGCNA_m01", "Several signals: organic-acid metabolism, mitochondrial respiration, mitochondrial organization/gene expression",
  "Signal 1 is led by carboxylic acid metabolic process (245 proteins), with organic acid metabolic process as its largest term (247); the previous 'organic-acid catabolism' named a smaller tied term (89).",
  "new", "microglia", "WGCNA_m02", "Organic-acid metabolism- and mitochondrial organization/gene expression-associated",
  "Carried from original microglia m01 (98% of its proteins) and takes that module's corrected words for the shared signal (rule 8); its own largest term is organic acid metabolic process (126 proteins).",
  "new", "neuron_soma", "WGCNA_m02", "Organic-acid metabolism-associated",
  "The signal's leading and largest term is carboxylic acid metabolic process (31 proteins); 'catabolism' named a smaller tied term.",
  "original", "neuron_neuropil", "WGCNA_m05", "Vesicle fusion/membrane organization-enriched",
  "'Membrane trafficking' is not a GO term: the signal's largest term is membrane organization (91 proteins), with organelle membrane fusion (leading term, 37) and vesicle fusion in the same signal.",
  "original", "neuron_neuropil", "WGCNA_m15", "unannotated (weak cornified envelope association)", STYLE,
  "original", "neuron_soma", "WGCNA_m03", "unannotated (weak anchoring junction association)", STYLE,
  "original", "microglia", "WGCNA_m05", "unannotated (weak actin filament/synaptic transmission association)", NA,
  "original", "microglia", "WGCNA_m09", "unannotated (weak cornified envelope association)", STYLE,
  "original", "microglia", "WGCNA_m10", "unannotated (weak actin filament association)", STYLE,
  "original", "microglia", "WGCNA_m13", "unannotated (weak ER protein complex association)", STYLE,
  "new", "microglia", "WGCNA_m08", "unannotated (weak COPII coat association)", STYLE,
  "new", "microglia", "WGCNA_m11", "unannotated (weak ER membrane association)", STYLE,
  "new", "microglia", "WGCNA_m12", "unannotated (weak synaptic vesicle cycle association)", STYLE,
  "new", "microglia", "WGCNA_m15", "unannotated (weak proteasome core association)", STYLE,
  "new", "neuron_neuropil", "WGCNA_m10", "Synaptic transmission- and immune response-regulating signaling-associated", STYLE,
  "new", "microglia", "WGCNA_m03", "Synaptic signaling-enriched", STYLE)
# reader-facing reasons for the corrections of the 2026-10-22 review (applied 2026-10-23)
SUMMARY <- tribble(~version, ~dataset, ~entity_id, ~summary,
  "original", "neuron_neuropil", "WGCNA_m03", "The dominant signal holds ER-to-Golgi vesicle transport (15 proteins; the module's strongest terms are coat complexes) as well as endosome-to-lysosome transport; the label names both.",
  "original", "neuron_neuropil", "WGCNA_m12", "About a quarter of the module is ribosomal: cytoplasmic translation (q = 2.4e-45, 71 proteins) is in the same signal as RNA processing.",
  "original", "microglia", "WGCNA_m05", paste("Four terms tie at q = 0.026; the rule lets the term with most proteins lead (actin filament-based process, 11), and synaptic-transmission terms carry the rest of the signal. The previous wording named the smallest tied term.", STYLE),
  "new", "neuron_neuropil", "WGCNA_m01", "The third signal's leading and largest term is intracellular protein transport (49 proteins), not its 15-protein sub-term.",
  "new", "neuron_neuropil", "WGCNA_m03", "As in original neuropil m12: cytoplasmic translation (71 proteins) is in the same signal as RNA processing.")
SUMMARY <- bind_rows(SUMMARY, MODULE_FIXES |> filter(!is.na(summary)) |> select(version, dataset, entity_id, summary))
# module labels after the regeneration of both networks on the corrected code (2026-10-24/26; applied after MODULE_FIXES)
REGENERATION_FIXES <- tribble(~version, ~dataset, ~entity_id, ~label, ~summary,
  "original", "neuron_neuropil", "WGCNA_m02", "Cellular respiration/organic acid metabolism-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): respiration and organic acid metabolism form one signal whose leading and largest term is the broad parent small molecule metabolic process (405 proteins); named after its specific terms cellular respiration and organic acid metabolic process.",
  "original", "neuron_neuropil", "WGCNA_m03", "Endosome to lysosome/ER to Golgi transport- and chromosome segregation-associated", "Wording (rule 7, label check 2026-10-26): named with GO's own words: leading term endosome to lysosome transport (12 proteins) with endoplasmic reticulum to Golgi vesicle-mediated transport (15) in the same signal; 'endolysosomal' is not a GO word and no hyphens are added.",
  "original", "neuron_neuropil", "WGCNA_m04", "Synaptic vesicle cycle/lysosomal protein localization- and protein ubiquitination-associated", "Regeneration (official gene mapping, GO term sizes 10-800): synaptic vesicle cycle and protein localization to lysosome now fall into one signal (largest term protein transport, 75 proteins, a broad parent); protein ubiquitination is the second signal.",
  "original", "neuron_neuropil", "WGCNA_m05", "Vesicle fusion/lysosomal transport-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): the signal's leading and largest term is now the broad parent intracellular transport (154 proteins); named after its specific terms vesicle fusion (36, q = 4.8e-08) and lysosomal transport (26); lysosome is the module's top CC term.",
  "original", "neuron_neuropil", "WGCNA_m06", "unannotated (weak proteasome complex association)", "Regeneration (official gene mapping, GO term sizes 10-800): one sub-threshold term group over BP, CC and MF (proteasome complex).",
  "original", "neuron_neuropil", "WGCNA_m10", "Several signals: intracellular signaling, positive regulation of cellular component biogenesis, proteasomal protein catabolism", "Regeneration (official gene mapping, GO term sizes 10-800): three comparable signals: intracellular signaling cassette (109 proteins); a signal led by positive regulation of cellular component biogenesis (47; largest term cell projection organization, 104, a broad parent), into which the synaptic transmission terms merged; protein modification with proteasomal protein catabolism (85).",
  "original", "neuron_neuropil", "WGCNA_m12", "RNA processing/splicing/cytoplasmic translation-enriched", "About a quarter of the module is ribosomal: cytoplasmic translation (69 proteins, q = 1.3e-43) is in the same signal as RNA processing.",
  "original", "neuron_neuropil", "WGCNA_m14", "Cytoplasmic translational initiation-associated", "Wording (rule 7, label check 2026-10-26): GO's spelling of the leading term, cytoplasmic translational initiation (14 proteins).",
  "original", "neuron_soma", "WGCNA_m02", "RNA processing/chromatin organization/transcription regulation-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): RNA processing and chromatin organization/transcription regulation now form one signal (largest term regulation of RNA metabolic process, 455 proteins), without a nuclear cellular component.",
  "original", "neuron_soma", "WGCNA_m03", "Cell-cell signaling-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): cell-cell signaling (43 proteins; with synaptic signaling) now passes the threshold; single signal, high confidence.",
  "original", "neuron_soma", "WGCNA_m04", "Several signals: cytoplasmic translation, protein catabolism, ER organization/Golgi vesicle transport", "Regeneration (official gene mapping, GO term sizes 10-800): the third signal now joins ER organization (27 proteins) and Golgi vesicle transport (61) in the CC term organelle subcompartment (leading term cytosolic transport, largest endomembrane system organization, 87).",
  "original", "neuron_soma", "WGCNA_m06", "Ion transport- and mitochondrial ATP synthesis-associated", "Regeneration (official gene mapping, GO term sizes 10-800): the dominant signal's leading and largest term is now the broad parent transmembrane transport (191 proteins), named after monoatomic ion transport; most proton transmembrane transport proteins (44 of 55) lie in the ATP synthesis signal.",
  "original", "neuron_soma", "WGCNA_m07", "Synaptic signaling/synapse organization-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): synaptic signaling (leading term, 135 proteins) and synapse organization now form one signal; single signal, high confidence.",
  "original", "microglia", "WGCNA_m01", "Oxidative phosphorylation/organic acid metabolism-associated", "Regeneration (official gene mapping, GO term sizes 10-800): oxidative phosphorylation (leading term, 92 proteins) and organic acid metabolism now form one signal (largest term small molecule metabolic process, 459, a broad parent); single signal, moderate confidence.",
  "original", "microglia", "WGCNA_m03", "Several signals: ER, microtubule, main axon", "Regeneration (official gene mapping, GO term sizes 10-800): the three largest comparable signals are the ER (largest term endoplasmic reticulum, 155 proteins; leading term organelle subcompartment, a broad parent), microtubule (largest term microtubule cytoskeleton, 122) and main axon (leading term; largest term axon, 120); secretory vesicle (81) is fourth.",
  "original", "microglia", "WGCNA_m05", "unannotated (weak actin filament-based process/synaptic transmission association)", "Wording (rule 7, label check 2026-10-26): the weak group is led by the process actin filament-based process (11 proteins; four terms tie at q = 0.028) with chemical synaptic transmission terms; 'actin filament' named a structure.",
  "original", "microglia", "WGCNA_m06", "Cell-substrate adhesion/extracellular matrix organization-enriched", "Wording (rule 7, label check 2026-10-26): the signal is led by cell adhesion (18 proteins, its largest term), with cell-substrate adhesion (12), extracellular matrix organization (9) and basement membrane organization (6); 'basement membrane' named a structure for these processes.",
  "original", "microglia", "WGCNA_m08", "Synapse organization-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): protein ubiquitination is no longer a comparable signal; single signal, high confidence.",
  "original", "microglia", "WGCNA_m10", "unannotated (weak actin filament-based process association)", "Wording (rule 7, label check 2026-10-26): the weak group's leading and largest term is the process actin filament-based process (10 proteins); 'actin filament' named a structure.",
  "original", "microglia", "WGCNA_m12", "Endomembrane organization/cation transport regulation-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): endomembrane organization and cation transport regulation now form one signal (leading and largest term protein transport, 41 proteins, a broad parent); single signal, high confidence.",
  "new", "neuron_neuropil", "WGCNA_m01", "Vacuolar transport- and cytoplasmic translational initiation/protein catabolism-associated", "Regeneration (official gene mapping, GO term sizes 10-800): two comparable signals: one led by vacuolar transport (28 proteins; largest term protein transport, 92) and one led by cytoplasmic translational initiation (17, q = 2.7e-06; largest term macromolecule catabolic process, 73) that also holds protein catabolic process (59).",
  "new", "neuron_neuropil", "WGCNA_m02", "Cellular respiration-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): the signal no longer has a mitochondrial cellular component, so it is named from its leading term; continues primary neuropil m02, whose organic acid part is not in this module's signal.",
  "new", "neuron_neuropil", "WGCNA_m03", "RNA processing/splicing/cytoplasmic translation-enriched", "As in primary neuropil m12: cytoplasmic translation (69 proteins, q = 2.0e-46) is in the same signal as RNA processing.",
  "new", "neuron_neuropil", "WGCNA_m06", "Organic acid metabolism-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): mitochondrial gene expression is no longer a comparable signal; uses primary neuropil m02's words for the shared organic acid signal (rule 8).",
  "new", "neuron_neuropil", "WGCNA_m09", "Golgi vesicle transport- and protein K48-linked ubiquitination-associated", "Regeneration (official gene mapping, GO term sizes 10-800): two comparable signals: Golgi vesicle transport (largest term protein transport, 200 proteins) and protein K48-linked ubiquitination (19; largest term macromolecule catabolic process, 131).",
  "new", "neuron_neuropil", "WGCNA_m10", "Intracellular signaling/neuron projection development-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): one signal led by intracellular signaling cassette (115 proteins; words of primary neuropil m10, rule 8); its largest term, cell projection organization (119), is a broad parent and is named by its specific descendant neuron projection development (94).",
  "new", "neuron_soma", "WGCNA_m01", "RNA processing/chromatin organization/transcription regulation-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): uses primary soma m02's words (rule 8): RNA processing and chromatin organization/transcription regulation now form one signal.",
  "new", "neuron_soma", "WGCNA_m02", "Organic acid catabolism-associated", "Regeneration (official gene mapping, GO term sizes 10-800): the leading term is now organic acid catabolic process (17 proteins; largest term carboxylic acid metabolic process, 32); no primary module's words apply (lineage revised).",
  "new", "neuron_soma", "WGCNA_m03", "Synaptic signaling/synapse organization-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): the vesicle docking/exocytosis group, previously a second comparable signal, now lies in the synaptic signaling signal; uses primary soma m07's words (rule 8).",
  "new", "neuron_soma", "WGCNA_m06", "Mitochondrial ATP synthesis/ion transport-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): ATP synthesis and ion transport now form one signal in the mitochondrial inner membrane (leading term proton motive force-driven ATP synthesis, 57 proteins; largest transmembrane transport, 166, with monoatomic ion transport, 144); uses primary soma m06's words (rule 8).",
  "new", "neuron_soma", "WGCNA_m07", "Mitochondrial translation- and early to late endosome transport-associated", "Wording (rule 7, label check 2026-10-26): GO's spelling without added hyphens (early endosome to late endosome transport).",
  "new", "neuron_soma", "WGCNA_m08", "Cell morphogenesis/actin cytoskeleton organization-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): the signal's leading term is now cell morphogenesis (23 proteins), with actin cytoskeleton organization; synapse organization no longer leads it.",
  "new", "microglia", "WGCNA_m02", "Organic acid catabolism-associated", "Regeneration (official gene mapping, GO term sizes 10-800): the leading term is organic acid catabolic process (58 proteins; largest term small molecule metabolic process, 205, a broad parent); mitochondrial organization/gene expression is no longer comparable; no primary module's words apply (lineage revised).",
  "new", "microglia", "WGCNA_m03", "Several signals: neuron projection morphogenesis, synaptic signaling, intracellular signaling", "Regeneration (official gene mapping, GO term sizes 10-800): four comparable signals: cell projection organization with neuron projection morphogenesis (185 proteins), synaptic signaling (141; largest term cell-cell signaling, 162), intracellular signaling cassette (160), protein modification (138). The module continues primary microglia m08, but synapse organization is not significant here, so its words are not carried.",
  "new", "microglia", "WGCNA_m06", "Cell-substrate adhesion/extracellular matrix organization-enriched", "Wording (rule 7, label check 2026-10-26): uses primary microglia m06's words (rule 8): cell-substrate adhesion leads the signal (13 proteins) with extracellular matrix organization (9) and basement membrane organization (6).",
  "new", "microglia", "WGCNA_m07", "Synapse organization-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): actin cytoskeleton organization is now in the same signal as synapse organization; uses primary microglia m08's words (rule 8).",
  "new", "microglia", "WGCNA_m08", "unannotated", "Regeneration (official gene mapping, GO term sizes 10-800): the single COPII coat term is no longer significant.",
  "new", "microglia", "WGCNA_m09", "Oxidative phosphorylation-enriched", "Regeneration (official gene mapping, GO term sizes 10-800): uses primary microglia m01's words for the shared respiration signal (rule 8); oxidative phosphorylation is a term of its signal (mitochondrial inner membrane).",
  "new", "microglia", "WGCNA_m11", "System process-associated", "Regeneration (official gene mapping, GO term sizes 10-800): system process (40 proteins) now passes the threshold (term sizes up to 800); its signal holds this single top-level term.",
  "new", "microglia", "WGCNA_m12", "Vesicle-mediated transport in synapse-associated", "Regeneration (official gene mapping, GO term sizes 10-800): vesicle-mediated transport in synapse (12 proteins, with synaptic vesicle cycle) now passes the threshold; single signal, moderate confidence.",
  "new", "microglia", "WGCNA_m14", "Small molecule metabolism-associated", "Regeneration (official gene mapping, GO term sizes 10-800): small molecule metabolic process (35 proteins) now passes the threshold (term sizes up to 800); its signal holds this single top-level term.",
  "new", "microglia", "WGCNA_m16", "Oxidative phosphorylation- and transmembrane transport-associated", "Regeneration (official gene mapping, GO term sizes 10-800): transmembrane transport is now a second comparable signal; uses primary microglia m01's words for the shared respiration signal (rule 8).")
SUMMARY <- bind_rows(SUMMARY |> anti_join(REGENERATION_FIXES, by = c("version", "dataset", "entity_id")),
                     REGENERATION_FIXES |> select(version, dataset, entity_id, summary))

lab <- read_tsv(file.path(RB, "labels.tsv"), show_col_types = FALSE, quote = "", col_types = cols(.default = "c"))
prev <- read_tsv(PREV, show_col_types = FALSE, quote = "", col_types = cols(.default = "c"))
mods <- read_csv(file.path(RB, "module_rule_classes.csv"), show_col_types = FALSE) |> select(version, dataset, ModuleID, signal_class, n_proteins)
sms <- read_csv(file.path(RB, "supermodule_rule_classes.csv"), show_col_types = FALSE)
# module corrections
for (i in seq_len(nrow(MODULE_FIXES))) {
  k <- lab$version == MODULE_FIXES$version[i] & lab$dataset == MODULE_FIXES$dataset[i] & lab$level == "module" & lab$entity_id == MODULE_FIXES$entity_id[i]
  stopifnot(sum(k) == 1); lab$label[k] <- MODULE_FIXES$label[i]
}
for (i in seq_len(nrow(REGENERATION_FIXES))) {
  k <- lab$version == REGENERATION_FIXES$version[i] & lab$dataset == REGENERATION_FIXES$dataset[i] & lab$level == "module" & lab$entity_id == REGENERATION_FIXES$entity_id[i]
  stopifnot(sum(k) == 1); lab$label[k] <- REGENERATION_FIXES$label[i]
}
mlab <- lab |> filter(level == "module") |> select(version, dataset, ModuleID = entity_id, label) |> left_join(mods, by = c("version", "dataset", "ModuleID"))
short <- function(m) sub("^WGCNA_", "", m)
lower_first <- function(x) {
  w <- regmatches(x, regexpr("^[A-Za-z]+", x))
  if (!length(w) || w %in% c("Golgi") || grepl("^[A-Z0-9]{2,}$", w)) x else paste0(tolower(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}
part_of <- function(label, cls) lower_first(switch(cls,
  "single signal" = sub("-(enriched|associated)$", "", label),
  "two signals" = sub("- and ", " and ", sub("-associated$", "", label), fixed = TRUE),
  "several signals" = { n <- strsplit(sub("^Several signals: ", "", label), ", ", fixed = TRUE)[[1]]; paste(n[1:2], collapse = " and ") },
  stop("no part for class ", cls)))
weak_of <- function(label) sub("^unannotated [(]weak (.*) association[)]$", "\\1", label)
and_list <- function(x) if (length(x) == 1) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
sm_label <- function(r) {
  m <- mlab |> filter(version == r$version, dataset == r$dataset, ModuleID %in% strsplit(r$member_modules, ";", fixed = TRUE)[[1]]) |> arrange(ModuleID)
  id <- short(m$ModuleID); labelled <- !grepl("^unannotated", m$signal_class); weak <- m$signal_class == "unannotated (weak)"
  switch(r$rule_class,
    "singleton" = sprintf("%s (%s)", m$label, id),
    "unannotated (weak)" = sprintf("unannotated (weak %s association; %s)", weak_of(m$label), id),
    "concordant" = sprintf("%s (%s)", if (identical(r$rule_confidence, "high")) m$label[which.max(m$n_proteins)] else
      sub("-enriched$", "-associated", m$label[which.max(m$n_proteins)]), paste(id, collapse = " + ")),   # suffix from the supermodule's confidence
    "compound" = paste0("Compound: ", paste(sprintf("%s (%s)", mapply(part_of, m$label, m$signal_class), id), collapse = " + ")),
    "mixed" = paste0("Mixed: ", paste(c(sprintf("%s (%s)", mapply(part_of, m$label[labelled], m$signal_class[labelled]), id[labelled]),
                                         sprintf("unannotated member (%s)", id[!labelled])), collapse = " + ")),
    "unannotated" = if (nrow(m) == 1) sprintf("unannotated (%s)", id) else if (any(weak)) {
      w <- sprintf("%s (%s)", vapply(m$label[weak], weak_of, ""), id[weak])
      paste0("unannotated (weak ", and_list(w), if (sum(weak) > 1) " associations" else " association",
             if (any(!weak)) paste0("; ", paste(id[!weak], collapse = ", ")) else "", ")")
    } else sprintf("unannotated (%s)", paste(id, collapse = ", ")),
    stop("unknown class ", r$rule_class))
}
gen <- sms |> rowwise() |> mutate(generated = sm_label(pick(everything()))) |> ungroup() |> select(version, dataset, entity_id = SupermoduleID, generated)
# supermodule rows are the current supermodules: their IDs follow the dendrogram leaf order, so a new cut height renumbers them
lab <- bind_rows(lab |> filter(level == "module"), gen |> transmute(version, dataset, level = "supermodule", entity_id, label = generated)) |>
  arrange(match(version, unique(prev$version)), match(level, c("module", "supermodule")), match(dataset, unique(prev$dataset)), entity_id)   # order of 2026-10-22
key <- function(d) sort(paste(d$version, d$dataset, d$entity_id)[d$level == "module"])
stopifnot(!anyNA(lab$label), identical(key(lab), key(prev)))
write_tsv(lab, file.path(RB, "labels.tsv"), na = "", escape = "none")

# change log: everything that differs from the labels preserved on 2026-10-22, with a reason. Modules are matched by ID,
# supermodules by their members: the cut-height rule of 2026-10-26 merged some supermodules and renumbered others.
CUT_CHANGE <- tibble(dataset = c("neuron_soma", "microglia"), cut_before = c(0.35, 0.40), cut_now = c(0.55, 0.45))
members <- sms |> transmute(version, dataset, entity_id = SupermoduleID, rule_class, member_modules)
prev_sm <- read_csv(file.path(dirname(PREV), "supermodule_rule_classes.csv"), show_col_types = FALSE) |>
  transmute(version, dataset, old_id = SupermoduleID, old_members = member_modules) |>
  left_join(prev |> filter(level == "supermodule") |> select(version, dataset, old_id = entity_id, old_label = label), by = c("version", "dataset", "old_id"))
stopifnot(!anyNA(prev_sm$old_label))   # every supermodule of 2026-10-22 had a label; a missing one would drop its change-log row
regen <- paste(REGENERATION_FIXES$version, REGENERATION_FIXES$dataset, REGENERATION_FIXES$entity_id)   # dated 2026-10-26 in the log
sm_before <- function(v, d, mm) {   # the supermodule with the same members before, or the earlier supermodules it joins
  p <- prev_sm[prev_sm$version == v & prev_sm$dataset == d, ]
  if (sum(p$old_members == mm) == 1) return(c(p$old_id[p$old_members == mm], p$old_label[p$old_members == mm], "same members"))
  ids <- strsplit(mm, ";", fixed = TRUE)[[1]]; om <- strsplit(p$old_members, ";", fixed = TRUE)
  p <- p[vapply(om, function(x) any(x %in% ids), logical(1)), ]
  stopifnot(nrow(p) > 1, all(unlist(strsplit(p$old_members, ";", fixed = TRUE)) %in% ids))   # a higher cut only merges
  c(paste(p$old_id, collapse = " + "), paste(sprintf("%s: %s", p$old_id, p$old_label), collapse = " | "), "merged")
}
b <- t(mapply(sm_before, members$version, members$dataset, members$member_modules, USE.NAMES = FALSE))
cur <- lab |> rename(new_label = label)
chg_mod <- prev |> filter(level == "module") |> rename(previous_label = label) |> inner_join(cur, by = c("version", "dataset", "level", "entity_id")) |>
  mutate(previous_entity_id = entity_id, match = "same ID")
chg_sm <- members |> mutate(level = "supermodule", previous_entity_id = b[, 1], previous_label = b[, 2], match = b[, 3]) |>
  inner_join(cur, by = c("version", "dataset", "level", "entity_id"))
chg <- bind_rows(chg_mod, chg_sm) |> filter(previous_label != new_label | previous_entity_id != entity_id) |>
  left_join(SUMMARY, by = c("version", "dataset", "entity_id")) |> select(-any_of(c("rule_class", "member_modules"))) |>
  left_join(members, by = c("version", "dataset", "entity_id")) |> left_join(CUT_CHANGE, by = "dataset") |>
  mutate(member_changed = mapply(function(v, d, mm) { if (is.na(mm)) return(FALSE); ids <- strsplit(mm, ";", fixed = TRUE)[[1]]
           any(prev$label[prev$version == v & prev$dataset == d & prev$entity_id %in% ids] != lab$label[lab$version == v & lab$dataset == d & lab$entity_id %in% ids]) },
           version, dataset, member_modules),
         member_regen = mapply(function(v, d, mm) !is.na(mm) && any(paste(v, d, strsplit(mm, ";", fixed = TRUE)[[1]]) %in% regen), version, dataset, member_modules),
         label_reason = case_when(!is.na(summary) ~ summary,
                                  previous_label == new_label ~ NA_character_,
                                  level == "supermodule" & !grepl(";", member_modules, fixed = TRUE) ~ "Follows its member's corrected label.",
                                  level == "supermodule" & member_changed ~ "Generated from the members' current labels with the supermodule template (rule 9, sheet 'How labels are formed'); a member's label was corrected.",
                                  level == "supermodule" ~ "Generated from the members' current labels with the supermodule template (rule 9, sheet 'How labels are formed'): members' own wording, 'and' between a member's two signals, no shortening.",
                                  TRUE ~ NA_character_),
         cut_reason = case_when(match == "merged" ~ sprintf("Formed when the supermodule cut height moved from %.2f to %.2f (stated rule, 2026-10-26): joins the members of %s; labelled from its members' current labels with the supermodule template (rule 9).",
                                                            cut_before, cut_now, vapply(strsplit(previous_entity_id, " + ", fixed = TRUE), and_list, "")),
                                previous_entity_id != entity_id ~ sprintf("Same members as %s at cut height %.2f; renumbered at %.2f (stated rule, 2026-10-26), as supermodule IDs follow the dendrogram order.",
                                                                          previous_entity_id, cut_before, cut_now),
                                TRUE ~ NA_character_),
         reason = ifelse(match == "merged", cut_reason, trimws(paste(coalesce(label_reason, ""), coalesce(cut_reason, "")))),
         changed_on = ifelse(!is.na(cut_reason) | paste(version, dataset, entity_id) %in% regen | member_regen, "2026-10-26", "2026-10-23")) |>
  select(version, dataset, level, entity_id, previous_entity_id, previous_label, new_label, reason, changed_on)
stopifnot(!anyNA(chg$reason), all(nzchar(chg$reason)), !any(grepl(" NA ", chg$reason, fixed = TRUE)))
write_csv(chg, file.path(ANN, "label_changes_2026-10-23.csv"), na = "")
cat(sprintf("labels: %d changed against 2026-10-22 (%d modules, %d supermodules)\n", nrow(chg), sum(chg$level == "module"), sum(chg$level == "supermodule")))
print(as.data.frame(chg |> filter(level == "supermodule") |> transmute(version, dataset = substr(dataset, 1, 9), entity_id, previous_entity_id, new_label = substr(new_label, 1, 120))), right = FALSE)
