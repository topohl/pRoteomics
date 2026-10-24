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
# Also applies the module corrections listed in MODULE_FIXES. Rewrites annotation_review/rule_based/labels.tsv and the
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

lab <- read_tsv(file.path(RB, "labels.tsv"), show_col_types = FALSE, quote = "", col_types = cols(.default = "c"))
prev <- read_tsv(PREV, show_col_types = FALSE, quote = "", col_types = cols(.default = "c"))
mods <- read_csv(file.path(RB, "module_rule_classes.csv"), show_col_types = FALSE) |> select(version, dataset, ModuleID, signal_class, n_proteins)
sms <- read_csv(file.path(RB, "supermodule_rule_classes.csv"), show_col_types = FALSE)
# module corrections
for (i in seq_len(nrow(MODULE_FIXES))) {
  k <- lab$version == MODULE_FIXES$version[i] & lab$dataset == MODULE_FIXES$dataset[i] & lab$level == "module" & lab$entity_id == MODULE_FIXES$entity_id[i]
  stopifnot(sum(k) == 1); lab$label[k] <- MODULE_FIXES$label[i]
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
