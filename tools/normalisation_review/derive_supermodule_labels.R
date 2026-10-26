#!/usr/bin/env Rscript
# Supermodule labels generated from their members' labels (normalisation review, 2026-10-23), so that every supermodule
# label follows one stated template and can be checked from the module labels alone:
#   member part  the member's own label without its suffix ('-enriched', '-associated'); a two-signal member contributes
#                'X and Y'; a several-signal member its first two named signals, 'X and Y'; '/' only ever joins terms of
#                one signal. The first letter is lower case unless the first word is an abbreviation (ER, RNA) or a name
#                (Golgi).
#   singleton    the member's label with its ID: 'X-enriched (m05)'; weak: 'unannotated (weak X association; m05)';
#                unannotated: 'unannotated (m05)'
#   concordant   the members' shared label with all IDs: 'X-enriched (m09 + m16)'
#   compound     'Compound: X (mA) + Y (mB)', members in ID order
#   mixed        the labelled members' parts, then '+ unannotated member (mC)': 'Mixed: X (mA) + unannotated member (mC)'
#   unannotated  several unlabelled members: 'unannotated (weak X (mA) and Y (mB) associations)'; plain unannotated
#                members are listed by ID after a semicolon; without any weak member 'unannotated (mA, mB)'
# Also applies the module corrections listed in MODULE_FIXES, then those of REGENERATION_FIXES (both networks rebuilt
# with the official gene mapping and GO term sizes 10-800, 2026-10-25/26). Rewrites annotation_review/rule_based/labels.tsv and the
# change log annotation_review/label_changes_2026-10-23.csv (all changes against the labels preserved in the review
# archive of 2026-10-22, with a short reason each). Run check_rule_labels.R afterwards.
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
# module labels after the regeneration of both networks on the corrected code (2026-10-25/26; applied after MODULE_FIXES)
REGENERATION_FIXES <- tribble(~version, ~dataset, ~entity_id, ~label, ~summary,
  "original", "neuron_neuropil", "WGCNA_m02", "Cellular respiration/organic-acid metabolism-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): respiration and organic-acid metabolism now form one signal (leading and largest term: small molecule metabolic process, 405 proteins); single signal, high confidence.",
  "original", "neuron_neuropil", "WGCNA_m04", "Synaptic vesicle cycle/lysosomal protein localization- and protein ubiquitination-associated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): synaptic vesicle cycle and lysosomal protein localization now fall into one signal (largest term: protein transport, 75 proteins); protein ubiquitination is the second signal.",
  "original", "neuron_neuropil", "WGCNA_m06", "unannotated (weak proteasome complex association)", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): one sub-threshold term group over BP, CC and MF (proteasome complex).",
  "original", "neuron_neuropil", "WGCNA_m10", "Several signals: intracellular signaling, cellular component biogenesis regulation, proteasomal protein catabolism", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): the three largest comparable signals are now intracellular signaling (109 proteins), cell projection organization with cellular component biogenesis regulation (104) and protein modification with proteasomal protein catabolism (85); synaptic transmission is no longer comparable.",
  "original", "neuron_soma", "WGCNA_m02", "RNA processing/chromatin organization/transcription regulation-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): RNA processing and chromatin organization/transcription regulation now form one signal (largest term: regulation of RNA metabolic process, 455 proteins), without a nuclear cellular component.",
  "original", "neuron_soma", "WGCNA_m03", "Cell-cell signaling-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): cell-cell signaling (43 proteins; with synaptic signaling) now passes the threshold; single signal, high confidence.",
  "original", "neuron_soma", "WGCNA_m07", "Synaptic signaling/synapse organization-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): synaptic signaling (leading term, 135 proteins) and synapse organization now form one signal; single signal, high confidence.",
  "original", "microglia", "WGCNA_m01", "Oxidative phosphorylation/organic-acid metabolism-associated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): oxidative phosphorylation (leading term, 92 proteins) and organic-acid metabolism now form one signal (largest term: small molecule metabolic process, 459); single signal, moderate confidence.",
  "original", "microglia", "WGCNA_m03", "Several signals: ER membrane, microtubule, main axon", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): the three largest comparable signals are now ER (155 proteins), microtubule cytoskeleton (122) and axon (120; leading term main axon); secretory vesicle (81) is fourth.",
  "original", "microglia", "WGCNA_m06", "Cell-substrate adhesion/extracellular matrix organization-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): basement membrane is no longer a term of the signal; extracellular matrix organization and cell-substrate adhesion are.",
  "original", "microglia", "WGCNA_m08", "Synapse organization-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): protein ubiquitination is no longer a comparable signal; single signal, high confidence.",
  "original", "microglia", "WGCNA_m12", "Endomembrane organization/cation transport regulation-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): endomembrane organization and cation transport regulation now form one signal (leading and largest term: protein transport, 41 proteins); single signal, high confidence.",
  "new", "neuron_neuropil", "WGCNA_m01", "Lysosomal transport- and protein catabolism-associated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): two comparable signals: lysosomal/vacuolar transport (largest term: protein transport, 92 proteins) and protein catabolism (73).",
  "new", "neuron_neuropil", "WGCNA_m02", "Cellular respiration-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): the signal no longer has a mitochondrial cellular component, so it is named from its leading term; carried from original neuropil m02 (rule 8), whose organic-acid part is not in this module's signal.",
  "new", "neuron_neuropil", "WGCNA_m06", "Organic-acid metabolism-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): mitochondrial gene expression is no longer a comparable signal; carries original neuropil m02's words for the shared organic-acid signal (rule 8).",
  "new", "neuron_neuropil", "WGCNA_m09", "Golgi vesicle transport- and protein K48-linked ubiquitination-associated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): two comparable signals: Golgi vesicle transport (largest term: protein transport, 200 proteins) and protein K48-linked ubiquitination (largest term: macromolecule catabolic process, 131).",
  "new", "neuron_neuropil", "WGCNA_m10", "Intracellular signaling/neuron development-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): one signal (leading term: intracellular signaling cassette, 115 proteins; with neuron development and synaptic transmission); carries original neuropil m10's words for the shared signal (rule 8).",
  "new", "neuron_soma", "WGCNA_m01", "RNA processing/chromatin organization/transcription regulation-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): carried from original soma m02 (rule 8): RNA processing and chromatin organization/transcription regulation now form one signal.",
  "new", "neuron_soma", "WGCNA_m03", "Synaptic signaling/synapse organization-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): carried from original soma m07 (rule 8): synaptic signaling and synapse organization now form one signal.",
  "new", "neuron_soma", "WGCNA_m06", "Mitochondrial ATP synthesis/proton transport-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): ATP synthesis and proton/ion transport now form one signal in the mitochondrial inner membrane; carries original soma m06's words (rule 8).",
  "new", "neuron_soma", "WGCNA_m08", "Cell morphogenesis/actin cytoskeleton organization-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): the signal's leading term is now cell morphogenesis (23 proteins), with actin cytoskeleton organization; synapse organization no longer leads it.",
  "new", "microglia", "WGCNA_m02", "Organic-acid metabolism-associated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): mitochondrial organization/gene expression is no longer a comparable signal; carries original microglia m01's words for the shared organic-acid signal (rule 8).",
  "new", "microglia", "WGCNA_m03", "Several signals: neuron projection morphogenesis, synaptic signaling, intracellular signaling", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): four comparable signals now: cell projection organization with neuron projection morphogenesis (185 proteins), synaptic signaling (162), intracellular signaling (160), protein modification (138); synapse organization is not a term of the dominant signal.",
  "new", "microglia", "WGCNA_m07", "Synapse organization-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): actin cytoskeleton organization is now in the same signal as synapse organization; carried from original microglia m08 (rule 8).",
  "new", "microglia", "WGCNA_m08", "unannotated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): the single COPII coat term is no longer significant.",
  "new", "microglia", "WGCNA_m09", "Oxidative phosphorylation-enriched", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): carries original microglia m01's words for the shared respiration signal (rule 8); oxidative phosphorylation is a term of its signal (mitochondrial inner membrane).",
  "new", "microglia", "WGCNA_m11", "System process-associated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): system process (40 proteins) now passes the threshold (term sizes up to 800); single signal, moderate confidence.",
  "new", "microglia", "WGCNA_m12", "Vesicle-mediated transport in synapse-associated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): vesicle-mediated transport in synapse (12 proteins, with synaptic vesicle cycle) now passes the threshold; single signal, moderate confidence.",
  "new", "microglia", "WGCNA_m14", "Small molecule metabolism-associated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): small molecule metabolic process (35 proteins) now passes the threshold (term sizes up to 800); single signal, moderate confidence.",
  "new", "microglia", "WGCNA_m16", "Oxidative phosphorylation- and transmembrane transport-associated", "Regenerated 2026-10-25 (official gene mapping, GO term sizes 10-800): transmembrane transport is now a second comparable signal; carries original microglia m01's words for the shared respiration signal (rule 8).")
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
    "concordant" = sprintf("%s (%s)", m$label[which.max(m$n_proteins)], paste(id, collapse = " + ")),
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
lab <- lab |> left_join(gen, by = c("version", "dataset", "entity_id")) |>
  mutate(label = ifelse(level == "supermodule" & !is.na(generated), generated, label)) |> select(-generated)
stopifnot(!anyNA(lab$label), nrow(lab) == nrow(prev))
write_tsv(lab, file.path(RB, "labels.tsv"), na = "", escape = "none")

# change log: everything that differs from the labels preserved on 2026-10-22, with a reason
members <- sms |> transmute(version, dataset, entity_id = SupermoduleID, rule_class, member_modules)
chg <- prev |> rename(previous_label = label) |> inner_join(lab |> rename(new_label = label), by = c("version", "dataset", "level", "entity_id")) |>
  filter(previous_label != new_label) |> left_join(SUMMARY, by = c("version", "dataset", "entity_id")) |> left_join(members, by = c("version", "dataset", "entity_id")) |>
  mutate(member_changed = mapply(function(v, d, mm) { if (is.na(mm)) return(FALSE); ids <- strsplit(mm, ";", fixed = TRUE)[[1]]
           any(prev$label[prev$version == v & prev$dataset == d & prev$entity_id %in% ids] != lab$label[lab$version == v & lab$dataset == d & lab$entity_id %in% ids]) },
           version, dataset, member_modules),
         reason = case_when(!is.na(summary) ~ summary,
                            level == "supermodule" & rule_class == "singleton" ~ "Follows its member's corrected label.",
                            level == "supermodule" & member_changed ~ "Generated from the members' current labels with the supermodule template (rule 9, sheet 'How labels are formed'); a member's label was corrected.",
                            level == "supermodule" ~ "Generated from the members' current labels with the supermodule template (rule 9, sheet 'How labels are formed'): members' own wording, 'and' between a member's two signals, no shortening.",
                            TRUE ~ NA_character_),
         changed_on = "2026-10-23") |>
  select(version, dataset, level, entity_id, previous_label, new_label, reason, changed_on)
stopifnot(!anyNA(chg$reason))
write_csv(chg, file.path(ANN, "label_changes_2026-10-23.csv"), na = "")
cat(sprintf("labels: %d changed against 2026-10-22 (%d modules, %d supermodules)\n", nrow(chg), sum(chg$level == "module"), sum(chg$level == "supermodule")))
print(as.data.frame(chg |> filter(level == "supermodule") |> transmute(version, dataset = substr(dataset, 1, 9), entity_id, new_label = substr(new_label, 1, 150))), right = FALSE)
