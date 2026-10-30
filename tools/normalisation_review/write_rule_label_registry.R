#!/usr/bin/env Rscript
# Write the rule-based module and supermodule labels (tools/normalisation_review/label_modules_by_rule.R and
# derive_supermodule_labels.R; rule_based/labels.tsv, version "original") into the reviewed label registries of the
# primary tree, so that the pipeline displays them (Stage 07, analysis/wgcna/summarize_module_interpretation.R).
# Decision of 2026-10-30 (Tobias Pohl): every neuropil and soma module and supermodule takes its rule-based label; in
# microglia only the five July rows that contradicted the rule-based labels are replaced; the other microglia rows keep
# their reviewed July labels.
#   neuron_neuropil, neuron_soma: partial registries (with member_ModuleIDs). Every row is written from the rule-based
#     labels on each run, except the soma rows reviewed on 2026-10-27 (m01, m03, SM04: the rule-based labels in reviewed
#     wording), which are kept as they are. The September neuropil m11 row is replaced (REPLACED).
#   microglia: full registry; the REPLACED rows are rewritten in place, every other row is kept. A singleton supermodule
#     repeats its member's label (the full registry's singleton contract). Microglia SM04 (kept, reviewed 2026-10-27)
#     names m09 in its text; that reference was updated by hand on 2026-10-30.
# Labels drop the rule's -enriched/-associated suffix, as the reviewed soma rows do; microglia labels start with a capital,
# as its other rows do ("unannotated" stays lower-case). Rule confidence "none" is recorded as low (the registries allow
# high, moderate and low).
# subcellular_context names the module's three strongest GO cellular-component terms (BH-adjusted p < 0.05, all module
# proteins, as the rule tool reads them); a multi-module supermodule's short label joins each member's first named part.
# The automatic label a row replaces is quoted from Stage 06 (annotate_module_microenvironment), which this does not change.
# Usage (primary tree root): Rscript tools/normalisation_review/write_rule_label_registry.R <rule_based folder> [--out DIR]
#   --out DIR writes the registries to DIR instead of config/wgcna_labels (for review); nothing else is written.

suppressPackageStartupMessages({ library(dplyr); library(readr) })
args <- commandArgs(trailingOnly = TRUE)
out_i <- match("--out", args)
out_dir <- if (!is.na(out_i)) args[[out_i + 1L]] else file.path("config", "wgcna_labels")
args <- if (!is.na(out_i)) args[-c(out_i, out_i + 1L)] else args
if (length(args) != 1L || !file.exists(file.path(args[[1]], "labels.tsv")))
  stop("usage: Rscript tools/normalisation_review/write_rule_label_registry.R <rule_based folder> [--out DIR]", call. = FALSE)
RB <- args[[1]]
if (!file.exists("pipeline.yml")) stop("Run from the primary tree root.", call. = FALSE)
norm <- yaml::read_yaml(file.path("config", "normalisation.yml"))
if (!identical(norm$sample_centring, norm$approved_sample_centring))
  stop("Not the primary tree (config/normalisation.yml); the registries belong to the original network.", call. = FALSE)

TOOL <- "tools/normalisation_review/write_rule_label_registry.R"
REVIEWER <- "Tobias Pohl"
REVIEW_DATE <- "2026-10-30"
KEPT_REVIEWED <- list(neuron_soma = c("WGCNA_m01", "WGCNA_m03", "SM04"), neuron_neuropil = character(0))
# reviewed labels the rule-based labels replace (decision of 2026-10-30), quoted in the new rows' rationale
REPLACED <- tribble(~dataset, ~entity_id, ~label, ~reviewer, ~review_date, ~note,
  "neuron_neuropil", "WGCNA_m11", "Enriched for myelin-associated proteins", "Part-29 registry adjudication", "2026-09-14", "",
  "microglia", "WGCNA_m02", "Membrane organization and ion handling", "Tobias Pohl", "2026-07-21", "",
  "microglia", "WGCNA_m04", "Synaptic vesicle endocytosis and membrane trafficking", "Tobias Pohl", "2026-07-21", "",
  "microglia", "WGCNA_m09", "Mixed metabolic and neuronal-structural proteins", "Tobias Pohl", "2026-07-21", "",
  "microglia", "SM03", "Mixed: vesicle/exocytic trafficking; microtubule-organizing structures", "Tobias Pohl", "2026-07-21", "",
  "microglia", "SM06", "Membrane organization and ion handling", "Tobias Pohl", "2026-07-21",
  " [Same members as SM07 at cut height 0.40; renumbered at 0.45 on 2026-10-26.]")
SHORT_LABELS <- c("microglia/SM03" = "Cell cycle / synapse / endomembrane")   # plot labels wrap at 18-34 characters
ROI <- c(neuron_neuropil = "Neuropil ROI", neuron_soma = "Neuronal soma ROI", microglia = "Shared microglia-neuropil ROI")
REG_COLS <- c("dataset", "level", "entity_id", "reviewed_biological_label", "reviewed_short_label", "subcellular_context",
              "roi_context", "label_confidence", "manual_review_required", "aggregation_evidence_class", "structural_status",
              "rationale", "reviewer", "proposal_prepared_by", "review_date", "adjudication_status")

rd <- function(f) read_csv(file.path(RB, f), show_col_types = FALSE, col_types = cols(.default = "c"))
newer <- function(f) if (file.mtime(file.path(RB, f)) < file.mtime(file.path(RB, "labels.tsv")))
  stop(f, " is older than labels.tsv; rerun check_rule_labels.R first.", call. = FALSE)
lab <- read_tsv(file.path(RB, "labels.tsv"), show_col_types = FALSE, quote = "", col_types = cols(.default = "c")) |>
  filter(version == "original")
form <- rd("label_form_check.csv") |> filter(version == "original")
if (any(nzchar(coalesce(form$problem, "")))) stop("label_form_check.csv reports problems; rerun check_rule_labels.R first.", call. = FALSE)
newer("label_form_check.csv"); newer("label_terms.csv")
mods <- rd("module_rule_classes.csv") |> filter(version == "original")
sigs <- rd("module_signals.csv") |> filter(version == "original")
sms <- rd("supermodule_rule_classes.csv") |> filter(version == "original")
terms <- rd("label_terms.csv") |> filter(version == "original")

num <- function(x) suppressWarnings(as.numeric(x))
fmt_q <- function(q) formatC(num(q), format = "g", digits = 2)
# the member-label form used inside supermodule labels: no -enriched/-associated suffix, lower-case first letter
part_of <- function(x) {
  x <- sub("-(enriched|associated)$", "", x)
  x <- gsub("- and ", " and ", x, fixed = TRUE)
  ifelse(grepl("^[A-Z][a-z]", x), paste0(tolower(substr(x, 1, 1)), substr(x, 2, nchar(x))), x)
}
capital <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
short_id <- function(x) sub("^WGCNA_", "", x)
replaces <- function(auto) if (is.na(auto)) "" else sprintf(" Replaces the automatic Stage 06 label '%s'.", auto)
replaced_note <- function(ds, id) {
  r <- REPLACED[REPLACED$dataset == ds & REPLACED$entity_id == id, ]
  if (!nrow(r)) return("")
  sprintf(" Replaces the reviewed label '%s' (%s, %s), which the rule-based labels %s.%s", r$label, r$reviewer,
          r$review_date, if (ds == "microglia") "contradict" else "name in GO's own words", r$note)
}

go_cc <- function(ds) {   # cellular-component terms of each module (all module proteins), as the rule tool reads them
  s <- readRDS(file.path("results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"))
  as.data.frame(s$GO_enrichment) |> filter(ModuleProteinSetType == "all", Ontology == "CC", p.adjust < 0.05) |>
    arrange(p.adjust, desc(Count))
}
GO_CC <- list()
cc_context <- function(ds, id) {
  if (is.null(GO_CC[[ds]])) GO_CC[[ds]] <<- go_cc(ds)
  x <- head(GO_CC[[ds]][GO_CC[[ds]]$ModuleID == id, , drop = FALSE], 3)
  if (!nrow(x)) return("")
  sprintf("%s (GO CC; strongest q = %s)", paste(x$Description, collapse = "; "), fmt_q(x$p.adjust[1]))
}
network <- function(ds) {   # current modules and supermodule members of the built network
  a <- read_csv(file.path("results", "wgcna", "build_wgcna_modules", ds, "tables", "supermodules",
                          "wgcna_module_supermodule_annotation.csv"), show_col_types = FALSE, col_types = cols(.default = "c"))
  a <- distinct(a, ModuleID, SupermoduleID)
  members <- tapply(a$ModuleID, a$SupermoduleID, function(v) paste(sort(v), collapse = ";"))
  list(modules = sort(unique(a$ModuleID)), members = members)
}
stage06 <- function(ds) {   # automatic labels the registry rows replace
  t6 <- file.path("results", "wgcna", "annotate_module_microenvironment", ds, "tables")
  m <- read_csv(file.path(t6, "WGCNA_module_biological_annotation.csv"), show_col_types = FALSE, col_types = cols(.default = "c"))
  s <- read_csv(file.path(t6, "WGCNA_supermodule_biological_annotation.csv"), show_col_types = FALSE, col_types = cols(.default = "c"))
  s_lab <- sub("^SM[0-9]+ \\S+ ", "", s$Supermodule_DisplayLabel)
  s_lab <- ifelse(s_lab == "singleton", paste0("singleton: ", s$DominantMemberTheme), s_lab)
  c(setNames(m$cleaned_biological_label, m$ModuleID), setNames(s_lab, s$SupermoduleID))
}

module_row <- function(ds, id, auto) {
  r <- lab |> filter(dataset == ds, level == "module", entity_id == id)
  f <- form |> filter(dataset == ds, level == "module", entity_id == id)
  m <- mods |> filter(dataset == ds, ModuleID == id)
  s <- sigs |> filter(dataset == ds, ModuleID == id, ontology_set == m$signal_basis) |> arrange(num(signal))
  stopifnot(nrow(r) == 1L, nrow(f) == 1L, nrow(m) == 1L)
  words <- terms |> filter(dataset == ds, ModuleID == id)
  named <- s |> filter(passing == "TRUE", signal %in% words$named_signal)   # the signals the label's words come from
  other <- s |> filter(passing == "TRUE", !signal %in% words$named_signal)
  sig_text <- if (nrow(named)) {
    paste0(if (nrow(words)) paste0("Words traced to GO terms: ", paste(sprintf("'%s' = %s (%s)", words$part, words$term, words$term_id),
                                                                       collapse = "; "), ". ") else "",
           "Named signal", if (nrow(named) > 1L) "s" else "", " (", m$signal_basis, "): ",
           paste(sprintf("%s, leading term %s (q = %s, %s proteins)", named$signal, named$leading_term, fmt_q(named$leading_q),
                         named$leading_count), collapse = "; "), ".",
           if (nrow(other)) paste0(" Other passing signal", if (nrow(other) > 1L) "s" else "", ", not named: ",
                                   paste(sprintf("%s (q = %s, %s proteins)", other$leading_term, fmt_q(other$leading_q),
                                                 other$leading_count), collapse = "; "), ".") else "")
  } else "No GO signal passes the rule thresholds."
  weak <- if (nzchar(coalesce(m$weak_association, ""))) sprintf(" Weak, sub-threshold association: %s (%s).", m$weak_association, m$weak_ontologies) else ""
  n_hub <- round(10 * num(m$hub_support))
  hubs <- if (!is.na(n_hub) && nrow(named))
    sprintf(" %d of the ten highest-kME hubs %s in the dominant signal.", n_hub, if (n_hub == 1) "is" else "are") else ""
  stab <- if (!is.na(num(m$stability_n_rebuilds)))
    sprintf(" Recovered in %s of %s leave-one-animal-out rebuilds (median Jaccard %.2f).", m$stability_n_recovered,
            m$stability_n_rebuilds, num(m$stability_median_jaccard)) else ""
  ctx <- cc_context(ds, id)
  label <- part_of(r$label)
  if (ds == "microglia" && !startsWith(label, "unannotated")) label <- capital(label)
  data.frame(dataset = ds, level = "module", entity_id = id,
    reviewed_biological_label = label, reviewed_short_label = label,
    subcellular_context = if (nzchar(ctx)) ctx else "no GO cellular-component term at q < 0.05",
    roi_context = paste0(ROI[[ds]], " co-abundance module; ", if (startsWith(label, "unannotated")) "no process label (unannotated)" else
      "a process label", ", not a cell-type identity"),
    label_confidence = if (f$confidence %in% c("high", "moderate", "low")) f$confidence else "low",
    manual_review_required = FALSE, aggregation_evidence_class = "not_applicable_module", structural_status = "module",
    rationale = paste0(sprintf("Rule-based label (class %s, rule confidence %s; %d proteins). ", f$class, f$confidence, as.integer(m$n_proteins)),
                       sig_text, weak, hubs, stab, replaces(auto), replaced_note(ds, id)),
    reviewer = REVIEWER,
    proposal_prepared_by = paste0("rule-based label (tools/normalisation_review/label_modules_by_rule.R, regenerated 2026-10-26); row written by ", TOOL),
    review_date = REVIEW_DATE, adjudication_status = "reviewed", stringsAsFactors = FALSE)
}

supermodule_row <- function(ds, id, members, auto) {
  r <- lab |> filter(dataset == ds, level == "supermodule", entity_id == id)
  f <- form |> filter(dataset == ds, level == "supermodule", entity_id == id)
  sm <- sms |> filter(dataset == ds, SupermoduleID == id)
  stopifnot(nrow(r) == 1L, nrow(f) == 1L, nrow(sm) == 1L)
  ids <- strsplit(members, ";", fixed = TRUE)[[1]]
  if (!setequal(ids, strsplit(sm$member_modules, ";", fixed = TRUE)[[1]]))
    stop(ds, " ", id, ": the rule-based label was made for other members (", sm$member_modules, "); rerun the label tools.", call. = FALSE)
  single <- length(ids) == 1L
  ctx <- vapply(ids, function(m) cc_context(ds, m), "")
  ctx <- if (single) ctx else ifelse(nzchar(ctx), paste0(short_id(ids), ": ", ctx), "")
  ctx <- paste(Filter(nzchar, ctx), collapse = " | ")
  member_label <- function(m) lab$label[lab$dataset == ds & lab$level == "module" & lab$entity_id == m]
  annotated <- vapply(ids, function(m) !grepl("^unannotated", member_label(m)), TRUE)
  # singleton: the member's label without the rule's suffix, with its module ID; else the rule label
  bio <- if (single && grepl(" \\(m[0-9]+\\)$", r$label)) paste0(part_of(sub(" \\(m[0-9]+\\)$", "", r$label)), " (", short_id(ids), ")") else r$label
  # short label: each annotated member's first named part, then the unannotated members
  first_part <- function(m) { x <- sub("^several signals: ", "", part_of(member_label(m))); strsplit(x, "/| and |, ")[[1]][1] }
  short <- if (single) bio else paste(c(sprintf("%s (%s)", vapply(ids[annotated], first_part, ""), short_id(ids[annotated])),
                                         if (any(!annotated)) sprintf("unannotated (%s)", paste(short_id(ids[!annotated]), collapse = ", "))),
                                       collapse = " + ")
  if (!is.na(SHORT_LABELS[paste0(ds, "/", id)])) short <- unname(SHORT_LABELS[paste0(ds, "/", id)])
  data.frame(dataset = ds, level = "supermodule", entity_id = id,
    reviewed_biological_label = bio, reviewed_short_label = short,
    subcellular_context = if (nzchar(ctx)) ctx else "no GO cellular-component term at q < 0.05 in its members",
    roi_context = paste0(ROI[[ds]], "; see the member modules"),
    label_confidence = if (f$confidence %in% c("high", "moderate", "low")) f$confidence else "low",
    manual_review_required = FALSE,
    aggregation_evidence_class = if (single) "not_applicable_singleton" else if (any(annotated)) "mixed_supported_themes" else "mixed_unresolved",
    structural_status = if (single) "singleton" else "multi_module",
    rationale = paste0(sprintf("Rule-based supermodule label (class %s, confidence %s), generated from the member labels with the supermodule template (rule 9): members %s.",
                               f$class, f$confidence, paste(short_id(ids), collapse = ", ")),
                       if (!single && !any(annotated)) " No member carries a passing GO signal." else "",
                       replaces(auto), replaced_note(ds, id)),
    reviewer = REVIEWER,
    proposal_prepared_by = paste0("rule-based supermodule label (tools/normalisation_review/derive_supermodule_labels.R); row written by ", TOOL),
    review_date = REVIEW_DATE, adjudication_status = "reviewed", stringsAsFactors = FALSE)
}

# the full registry's singleton: the member's label, short label, contexts, confidence and manual flag
singleton_row <- function(ds, id, member, member_row) {
  row <- member_row
  row$level <- "supermodule"; row$entity_id <- id
  row$aggregation_evidence_class <- "not_applicable_singleton"; row$structural_status <- "singleton"
  row$rationale <- paste0("Singleton inheriting current member ", member, ": ", sub(" Replaces the reviewed label .*$", "", member_row$rationale),
                          replaced_note(ds, id))
  row$proposal_prepared_by <- paste0(member_row$proposal_prepared_by, "; the member's label, inherited by its singleton supermodule")
  row
}

eol_of <- function(path) if (file.exists(path) && grepl("\r\n", readChar(path, min(file.size(path), 4096L), useBytes = TRUE), fixed = TRUE)) "\r\n" else "\n"
write_registry <- function(x, path, eol) {
  x$manual_review_required <- as.logical(x$manual_review_required)
  con <- file(path, "wb"); on.exit(close(con))
  utils::write.csv(x, con, row.names = FALSE, na = "", eol = eol)
}
read_registry <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character",
                                                 na.strings = character(0), encoding = "UTF-8")
entity_order <- function(x) order(x$level != "module", x$entity_id)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (ds in c("neuron_neuropil", "neuron_soma")) {
  src <- file.path("config", "wgcna_labels", paste0(ds, ".csv"))
  old <- if (file.exists(src)) read_registry(src) else NULL
  if (!is.null(old) && !"member_ModuleIDs" %in% names(old)) old$member_ModuleIDs <- ""   # older 16-column file
  kept <- if (is.null(old)) NULL else old[old$entity_id %in% KEPT_REVIEWED[[ds]], c(REG_COLS, "member_ModuleIDs"), drop = FALSE]
  if (NROW(kept) != length(KEPT_REVIEWED[[ds]])) stop(ds, ": a kept reviewed row is missing from ", src, call. = FALSE)
  net <- network(ds); auto <- stage06(ds)
  rows <- list()
  for (id in setdiff(net$modules, kept$entity_id))
    rows[[length(rows) + 1L]] <- cbind(module_row(ds, id, auto[[id]]), member_ModuleIDs = "")
  for (id in setdiff(names(net$members), kept$entity_id))
    rows[[length(rows) + 1L]] <- cbind(supermodule_row(ds, id, net$members[[id]], auto[[id]]), member_ModuleIDs = net$members[[id]])
  new <- rbind(kept, do.call(rbind, rows))
  new <- new[entity_order(new), , drop = FALSE]
  n_ent <- length(net$modules) + length(net$members)
  if (nrow(new) != n_ent || anyDuplicated(paste(new$level, new$entity_id)))
    stop(ds, ": expected one row per module and supermodule (", n_ent, "), got ", nrow(new), call. = FALSE)
  write_registry(new, file.path(out_dir, paste0(ds, ".csv")), eol_of(src))
  cat(sprintf("%s: %d rows (%d kept as reviewed, %d written from the rule-based labels)\n", ds, nrow(new), NROW(kept), length(rows)))
}

# microglia: rewrite the REPLACED rows in place (full registry, 16 columns); modules first, as singletons copy them
src <- file.path("config", "wgcna_labels", "microglia.csv")
mg <- read_registry(src)
net <- network("microglia")
ids <- REPLACED$entity_id[REPLACED$dataset == "microglia"]
for (id in ids[order(!startsWith(ids, "WGCNA_"))]) {
  k <- which(mg$entity_id == id)
  stopifnot(length(k) == 1L)
  members <- strsplit(net$members[id], ";", fixed = TRUE)[[1]]
  row <- if (mg$level[k] == "module") module_row("microglia", id, NA) else if (length(members) == 1L)
    singleton_row("microglia", id, members, mg[mg$level == "module" & mg$entity_id == members, REG_COLS, drop = FALSE]) else
    supermodule_row("microglia", id, net$members[[id]], NA)   # microglia displayed its July label, not Stage 06
  if (mg$level[k] == "supermodule" && length(members) > 1L) {
    kept_members <- mg[mg$level == "module" & mg$entity_id %in% members & !grepl(TOOL, mg$proposal_prepared_by, fixed = TRUE), , drop = FALSE]
    if (nrow(kept_members)) row$rationale <- paste0(row$rationale, " Its member parts are the rule-based member labels; ",
      paste(sprintf("%s displays its reviewed label '%s'", short_id(kept_members$entity_id), kept_members$reviewed_biological_label),
            collapse = " and "), ".")
  }
  mg[k, REG_COLS] <- row[, REG_COLS]
}
write_registry(mg[, REG_COLS], file.path(out_dir, "microglia.csv"), eol_of(src))
cat(sprintf("microglia: %d rows (%d written from the rule-based labels: %s)\n", nrow(mg), length(ids), paste(ids, collapse = ", ")))
