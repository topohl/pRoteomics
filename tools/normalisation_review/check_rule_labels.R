# Mechanical check of the hand-written label wording against the scripted classes (normalisation review 2026-10-22):
# the form of every label must follow from its class and confidence, and every supermodule label must name all members.
#   module   single signal: '...-enriched' when confidence is high, '...-associated' otherwise; two signals: 'A- and
#            B-associated'; several signals: 'Several signals: ...'; unannotated (weak): 'unannotated (weak ... association)';
#            unannotated: 'unannotated'.
#   supermodule  singleton: the member's label with its ID; compound: 'Compound: ... (mXX) + ... (mYY)'; concordant: the
#            shared label with all IDs; mixed: 'Mixed: ...' with all IDs; unannotated: 'unannotated ...' with the IDs.
# Writes annotation_review/rule_based/label_form_check.csv. USAGE  Rscript tools/normalisation_review/check_rule_labels.R
suppressPackageStartupMessages({ library(dplyr); library(readr) })
RB <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "annotation_review", "rule_based")
lab <- read_tsv(file.path(RB, "labels.tsv"), show_col_types = FALSE, quote = "")
mods <- read_csv(file.path(RB, "module_rule_classes.csv"), show_col_types = FALSE)
sms <- read_csv(file.path(RB, "supermodule_rule_classes.csv"), show_col_types = FALSE)
m_ok <- function(cls, conf, l) switch(cls,
  "single signal" = if (conf == "high") grepl("-enriched$", l) else grepl("-associated$", l) && !grepl("- and ", l),
  "two signals" = grepl("- and .*-associated$", l),
  "several signals" = grepl("^Several signals: ", l),
  "unannotated (weak)" = grepl("^unannotated [(]weak .*association[)]$", l),
  "unannotated" = l == "unannotated", FALSE)
mc <- mods |> select(version, dataset, entity_id = ModuleID, signal_class, rule_confidence) |>
  full_join(lab |> filter(level == "module"), by = c("version", "dataset", "entity_id")) |>
  mutate(level = "module", ok = mapply(function(c, f, l) !is.na(l) && !is.na(c) && m_ok(c, f, l), signal_class, rule_confidence, label),
         problem = ifelse(ok, "", ifelse(is.na(label), "no label", ifelse(is.na(signal_class), "no class", "form does not match class"))))
sc <- sms |> select(version, dataset, entity_id = SupermoduleID, member_modules, rule_class, rule_confidence) |>
  full_join(lab |> filter(level == "supermodule"), by = c("version", "dataset", "entity_id")) |>
  rowwise() |> mutate(level = "supermodule",
    ids = list(sub("WGCNA_", "", strsplit(member_modules, ";", fixed = TRUE)[[1]])),
    all_ids = all(vapply(ids, function(i) grepl(paste0("\\b", i, "\\b"), label), logical(1))),
    form = switch(rule_class, compound = grepl("^Compound: ", label), mixed = grepl("^Mixed: ", label),
                  unannotated = grepl("^unannotated", label), "unannotated (weak)" = grepl("^unannotated [(]weak ", label),
                  concordant = !grepl("^(Compound|Mixed)", label), singleton = !grepl("^(Compound|Mixed)", label), FALSE),
    ok = !is.na(label) && all_ids && form,
    problem = if (ok) "" else paste(c(if (is.na(label)) "no label", if (!isTRUE(all_ids)) "member ID missing", if (!isTRUE(form)) "form does not match class"), collapse = "; ")) |>
  ungroup()
# singleton supermodules carry the member's label
single <- sc |> filter(rule_class %in% c("singleton", "unannotated (weak)", "unannotated"), lengths(ids) == 1) |>
  mutate(member = paste0("WGCNA_", unlist(ids))) |>
  left_join(lab |> filter(level == "module") |> select(version, dataset, member = entity_id, member_label = label), by = c("version", "dataset", "member")) |>
  mutate(same = startsWith(label, sub(" [(]?weak .*$", "", member_label)) | grepl("^unannotated", label)) |> filter(!same)
out <- bind_rows(mc |> select(version, dataset, level, entity_id, class = signal_class, confidence = rule_confidence, label, problem),
                 sc |> select(version, dataset, level, entity_id, class = rule_class, confidence = rule_confidence, label, problem))
if (nrow(single)) out$problem[match(paste(single$version, single$dataset, single$entity_id), paste(out$version, out$dataset, out$entity_id))] <- "singleton label differs from its member's"
write_csv(out, file.path(RB, "label_form_check.csv"), na = "")
cat(sum(out$problem == ""), "of", nrow(out), "labels pass\n")
print(out |> filter(problem != "") |> select(version, dataset, entity_id, class, confidence, problem, label), width = 220)
