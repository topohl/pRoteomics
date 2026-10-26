# Mechanical check of the hand-written label wording against the scripted classes and the GO evidence:
#   form     the form of every label follows from its class and confidence, and every supermodule label names all members.
#     module   single signal: '...-enriched' when confidence is high, '...-associated' otherwise; two signals: 'A- and
#              B-associated'; several signals: 'Several signals: ...'; unannotated (weak): 'unannotated (weak ... association)';
#              unannotated: 'unannotated'.
#     supermodule  singleton: the member's label with its ID; compound: 'Compound: ... (mXX) + ... (mYY)'; concordant: the
#              largest member's wording with all IDs; mixed: 'Mixed: ...' with all IDs; unannotated: 'unannotated ...' with the IDs.
#   terms    every term named in a module label (the parts between ',', '/' and the templates) is a GO term of the signal it
#            names: its words are the term's words, shortened only as rule 7 allows ('metabolic process' -> 'metabolism',
#            'catabolic process' -> 'catabolism', 'biosynthetic process' -> 'biosynthesis', 'endoplasmic reticulum' -> 'ER',
#            an adjective for its noun, e.g. 'lysosomal' for 'lysosome'), optionally qualified by the signal's cellular
#            component; no hyphen is added to GO's spelling. Two signals name the dominant signal first; several signals name
#            the three largest comparable signals, largest first; a weak label names the single sub-threshold term group.
#   unique   no two module labels of one network are the same, except 'unannotated'.
# Writes annotation_review/rule_based/label_form_check.csv (one row per label) and label_terms.csv (one row per named term,
# with its GO ID and signal). USAGE  Rscript tools/normalisation_review/check_rule_labels.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })
RB <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "annotation_review", "rule_based")
lab <- read_tsv(file.path(RB, "labels.tsv"), show_col_types = FALSE, quote = "")
mods <- read_csv(file.path(RB, "module_rule_classes.csv"), show_col_types = FALSE)
sms <- read_csv(file.path(RB, "supermodule_rule_classes.csv"), show_col_types = FALSE)
rsig <- bind_rows(read_csv(file.path(RB, "module_signals.csv"), show_col_types = FALSE),
                  read_csv(file.path(RB, "module_signals_all_ontologies.csv"), show_col_types = FALSE))
rterm <- read_csv(file.path(RB, "signal_terms.csv"), show_col_types = FALSE)
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
flag <- function(out, keys, msg) {   # append a problem to the rows named by keys (version, dataset, level, entity_id)
  if (!nrow(keys)) return(out)
  i <- match(paste(keys$version, keys$dataset, keys$level, keys$entity_id), paste(out$version, out$dataset, out$level, out$entity_id))
  out$problem[i] <- ifelse(nzchar(out$problem[i]), paste(out$problem[i], msg, sep = "; "), msg); out
}
out <- flag(out, single |> transmute(version, dataset, level = "supermodule", entity_id), "singleton label differs from its member's")

# ----------------------------------------------------------------------------------- named terms
# a module label as its named signals, each a list of terms ('/' joins terms of one signal)
parts_of <- function(cls, l) {
  body <- switch(cls,
    "single signal" = sub("-(enriched|associated)$", "", l),
    "two signals" = strsplit(sub("-associated$", "", l), "- and ", fixed = TRUE)[[1]],
    "several signals" = strsplit(sub("^Several signals: ", "", l), ", ", fixed = TRUE)[[1]],
    "unannotated (weak)" = sub("^unannotated [(]weak (.*) association[)]$", "\\1", l), character())
  tibble(named_signal = seq_along(body), part = body) |> separate_rows(part, sep = "/") |> mutate(part = trimws(part))
}
phrases <- mc |> filter(signal_class != "unannotated", !is.na(label)) |> rowwise() |>
  reframe(version, dataset, ModuleID = entity_id, signal_class, parts_of(signal_class, label))
go_names <- local({ ids <- unique(rterm$term_id)
  t <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = ids, columns = "TERM", keytype = "GOID")); setNames(t$TERM, t$GOID) })
SHORT <- c("metabolic process" = "metabolism", "catabolic process" = "catabolism", "biosynthetic process" = "biosynthesis", "endoplasmic reticulum" = "er")
ADJ <- c(lysosomal = "lysosome", mitochondrial = "mitochondrion", vacuolar = "vacuole", proteasomal = "proteasome", endosomal = "endosome",
         synaptic = "synapse", nuclear = "nucleus", ribosomal = "ribosome")
norm <- function(x) { x <- gsub(",", "", tolower(x)); for (k in names(SHORT)) x <- gsub(k, SHORT[[k]], x, fixed = TRUE); gsub("\\s+", " ", trimws(x)) }
words <- function(x) unique(strsplit(norm(x), "[ -]+")[[1]])
fits <- function(part, term, qualifier) {   # words of the named part that the term (or the signal's component) does not supply
  w <- words(part); tw <- words(term); qw <- if (is.na(qualifier) || !nzchar(qualifier)) character() else words(qualifier)
  miss <- setdiff(w, c(tw, qw)); miss <- miss[!(miss %in% names(ADJ) & ADJ[miss] %in% tw)]; length(miss) == 0
}
hyphen_ok <- function(part, term) {   # every hyphenated word of the part is spelled so in the term; no added hyphens
  h <- grep("-", strsplit(tolower(part), " ", fixed = TRUE)[[1]], value = TRUE); all(vapply(h, function(x) grepl(x, tolower(term), fixed = TRUE), logical(1)))
}
grouping_of <- function(cls, basis) ifelse(cls == "unannotated (weak)", "all ontologies (decides weak or unannotated)",
                                           ifelse(basis == "BP", "BP", "all ontologies (no passing BP signal)"))
sig_terms <- rterm |> mutate(term = unname(go_names[term_id])) |>
  left_join(rsig |> select(version, dataset, ModuleID, ontology_set, signal, comparable, passing, compartment), by = c("version", "dataset", "ModuleID", "ontology_set", "signal"))
cand <- phrases |> left_join(mods |> select(version, dataset, ModuleID, signal_basis), by = c("version", "dataset", "ModuleID")) |>
  mutate(ontology_set = grouping_of(signal_class, signal_basis), row = row_number()) |>
  inner_join(sig_terms, by = c("version", "dataset", "ModuleID", "ontology_set"), relationship = "many-to-many") |>
  filter(signal_class == "unannotated (weak)" | comparable %in% TRUE) |>
  rowwise() |> filter(fits(part, term, compartment)) |> ungroup() |>
  mutate(extra = mapply(function(p, t) length(setdiff(words(t), words(p))), part, term), exact = norm(part) == norm(term),
         qualified = mapply(function(p, t) !fits(p, t, NA), part, term))   # needs the signal's component for a word
# each named signal is the one signal whose terms supply all its parts; within it, each part names the term closest to its
# words, or the term chosen here where the closest words name another process (ATP synthesis, not electron transport)
TERM_CHOICE <- tribble(~version, ~dataset, ~ModuleID, ~part, ~choice,
  "original", "neuron_soma", "WGCNA_m06", "mitochondrial ATP synthesis", "GO:0042776",
  "new", "neuron_soma", "WGCNA_m06", "Mitochondrial ATP synthesis", "GO:0042776")
cand <- cand |> left_join(TERM_CHOICE, by = c("version", "dataset", "ModuleID", "part")) |> mutate(chosen = !is.na(choice) & term_id == choice)
stopifnot(all(TERM_CHOICE$choice %in% cand$term_id))
sig_pick <- cand |> distinct(row, version, dataset, ModuleID, named_signal, signal) |> count(version, dataset, ModuleID, named_signal, signal, name = "n_parts") |>
  inner_join(phrases |> count(version, dataset, ModuleID, named_signal, name = "parts"), by = c("version", "dataset", "ModuleID", "named_signal")) |>
  filter(n_parts == parts) |> group_by(version, dataset, ModuleID, named_signal) |> slice_min(signal, n = 1) |> ungroup() |> select(-n_parts, -parts)
terms <- phrases |> mutate(row = row_number()) |>
  left_join(cand |> semi_join(sig_pick, by = c("version", "dataset", "ModuleID", "named_signal", "signal")) |> group_by(row) |>
              arrange(desc(chosen), desc(exact), qualified, extra, .by_group = TRUE) |> slice_head(n = 1) |> ungroup() |> select(row, ontology_set, signal, term_id, term, exact, qualified),
            by = "row") |>
  mutate(hyphens_ok = ifelse(is.na(term), NA, mapply(hyphen_ok, part, coalesce(term, ""))))
# problems per module: an unmatched part, an added hyphen, signals named out of order or not the largest comparable ones
order_bad <- terms |> filter(!is.na(signal), signal_class %in% c("two signals", "several signals")) |> distinct(version, dataset, ModuleID, named_signal, signal) |>
  group_by(version, dataset, ModuleID) |> summarise(bad = !identical(as.integer(signal), seq_along(signal)), .groups = "drop") |> filter(bad)
prob_mod <- bind_rows(
  terms |> filter(is.na(term_id)) |> transmute(version, dataset, entity_id = ModuleID, msg = sprintf("'%s' is not a term of the named signal", part)),
  terms |> filter(hyphens_ok %in% FALSE) |> transmute(version, dataset, entity_id = ModuleID, msg = sprintf("'%s' adds a hyphen to GO's '%s'", part, term)),
  order_bad |> transmute(version, dataset, entity_id = ModuleID, msg = "signals not named largest first from the dominant one"),
  mc |> filter(!is.na(label), label != "unannotated") |> group_by(version, dataset, label) |> filter(n() > 1) |> ungroup() |>
    transmute(version, dataset, entity_id, msg = "label not unique in its network")) |>
  group_by(version, dataset, entity_id) |> summarise(msg = paste(unique(msg), collapse = "; "), .groups = "drop") |> mutate(level = "module")
for (i in seq_len(nrow(prob_mod))) out <- flag(out, prob_mod[i, ], prob_mod$msg[i])
# supermodule labels repeat their members' words; no hyphen may be added there either
sm_h <- sc |> filter(grepl("[a-z]-(acid|to)\\b", label)) |> transmute(version, dataset, level = "supermodule", entity_id)
out <- flag(out, sm_h, "hyphen added to a GO term")
write_csv(out, file.path(RB, "label_form_check.csv"), na = "")
write_csv(terms |> select(version, dataset, ModuleID, signal_class, named_signal, part, grouping = ontology_set, signal, term_id, term, exact_wording = exact, qualified),
          file.path(RB, "label_terms.csv"), na = "")
cat(sum(out$problem == ""), "of", nrow(out), "labels pass;", sum(!is.na(terms$term_id)), "of", nrow(terms), "named terms traced to a GO term of their signal\n")
print(out |> filter(problem != "") |> select(version, dataset, entity_id, class, confidence, problem, label), width = 220)
