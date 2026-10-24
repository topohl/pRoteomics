# Mechanical check of every GO count quoted in annotation_review/proposed_annotations.tsv against the evidence tables:
# "term a/b module, a/c term" must match Count = a, GeneRatio denominator = b, BgRatio numerator = c of that term in
# that module (module rows) or in the direct supermodule enrichment ("direct: term a/b", supermodule rows).
# Writes annotation_review/proposal_number_check.csv. USAGE  Rscript tools/normalisation_review/check_proposal_numbers.R
suppressPackageStartupMessages({ library(dplyr); library(readr) })
ANN <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "annotation_review")
prop <- read_tsv(file.path(ANN, "proposed_annotations.tsv"), show_col_types = FALSE, quote = "")
mg <- read_csv(file.path(ANN, "module_go_top_terms.csv"), show_col_types = FALSE)
# the module table keeps the top 15 terms per ontology; look terms up in the full build tables as well
full <- bind_rows(lapply(unique(prop$dataset), function(ds) {
  s <- readRDS(file.path("results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"))
  as.data.frame(s$GO_enrichment) |> filter(ModuleProteinSetType == "all") |>
    transmute(dataset = ds, ModuleID, Ontology, Description, GeneRatio, BgRatio, Count, p.adjust)
}))
sg <- read_csv(file.path(ANN, "supermodule_direct_go_top_terms.csv"), show_col_types = FALSE)
num <- function(x, part) as.integer(if (part == 1) sub("/.*", "", x) else sub(".*/", "", x))
out <- list()
for (i in seq_len(nrow(prop))) {
  r <- prop[i, ]
  txt <- r$key_evidence
  if (is.na(txt)) next
  txt <- sub("^direct: ", "", txt)
  parts <- trimws(strsplit(txt, ";", fixed = TRUE)[[1]])
  for (p in parts) {
    m <- regmatches(p, regexec("^(.+?) ([0-9]+)/([0-9]+)( module)?(, ([0-9]+)/([0-9]+) term)?", p))[[1]]
    if (!length(m)) next
    term <- trimws(m[2]); a <- as.integer(m[3]); b <- as.integer(m[4]); c_term <- if (nzchar(m[8])) as.integer(m[8]) else NA_integer_
    if (r$level == "module") {
      hit <- full |> filter(dataset == r$dataset, ModuleID == r$entity_id, tolower(Description) == tolower(term))
    } else {
      hit <- sg |> filter(dataset == r$dataset, SupermoduleID == r$entity_id, tolower(Description) == tolower(term))
    }
    status <- if (!nrow(hit)) "term not found" else {
      ok_a <- any(hit$Count == a); ok_b <- any(num(hit$GeneRatio, 2) == b)
      ok_c <- is.na(c_term) || any(num(hit$BgRatio, 1) == c_term)
      if (ok_a && ok_b && ok_c) "ok" else sprintf("mismatch: file has %s", paste(unique(sprintf("%s/%s (term %s)", hit$Count, num(hit$GeneRatio, 2), num(hit$BgRatio, 1))), collapse = " | "))
    }
    out[[length(out) + 1L]] <- tibble(dataset = r$dataset, level = r$level, entity_id = r$entity_id, quoted = p, status = status)
  }
}
res <- bind_rows(out)
write_csv(res, file.path(ANN, "proposal_number_check.csv"))
print(count(res, status != "ok"))
print(as.data.frame(filter(res, status != "ok")), row.names = FALSE)
