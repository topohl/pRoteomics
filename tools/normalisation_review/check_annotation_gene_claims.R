# Check the gene-level statements in annotation_review/proposed_annotations.tsv ("hub X is in no significant term",
# "hub X carries term Y") against every significant GO term of the module (full build table, not the top-15 table),
# and list hubs without an Entrez ID (they can never be in a GO term).
# Writes annotation_review/gene_claim_check.csv. USAGE  Rscript tools/normalisation_review/check_annotation_gene_claims.R
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })
ANN <- file.path("results", "reviewer_audit", "normalisation_choice_2026-10-21", "annotation_review")
go <- bind_rows(lapply(c("neuron_neuropil", "neuron_soma", "microglia"), function(ds) {
  s <- readRDS(file.path("results", "wgcna", "build_wgcna_modules", ds, "models", "wgcna_final_model_state.rds"))
  as.data.frame(s$GO_enrichment) |> filter(ModuleProteinSetType == "all", p.adjust < 0.05) |>
    transmute(dataset = ds, ModuleID, Ontology, Description, geneID)
})) |> separate_rows(geneID, sep = "/")
mem <- bind_rows(lapply(c("neuron_neuropil", "neuron_soma", "microglia"), function(ds)
  read_csv(file.path("results", "wgcna", "build_wgcna_modules", ds, "tables", "modules", "WGCNA_modules_long.csv"),
           show_col_types = FALSE, col_select = c(ModuleID, GeneSymbol, EntrezID)) |> mutate(dataset = ds)))
# claim: dataset, module, gene, expected = "none" (in no significant term) or a regular expression the terms must match
claims <- tribble(
  ~dataset, ~module, ~genes, ~expected,
  "neuron_neuropil", "WGCNA_m01", "SCG2,SYT17", "none",
  "neuron_neuropil", "WGCNA_m01", "PTPRN,PDYN", "not gene-mapped",
  "neuron_neuropil", "WGCNA_m01", "PAM", "^(protein catabolic process|cytoplasmic vesicle membrane|perinuclear region of cytoplasm)$",
  "neuron_neuropil", "WGCNA_m01", "EHD1,EHD3", "lysosom|vacuol|endosom|endocyt",
  "neuron_neuropil", "WGCNA_m07", "RASGRF1,DIAPH1", "^(?!.*(P-body|CCR4-NOT)).*$",
  "neuron_neuropil", "WGCNA_m07", "SPTB", "not gene-mapped",
  "neuron_neuropil", "WGCNA_m07", "CNOT1", "P-body|CCR4-NOT",
  "neuron_neuropil", "WGCNA_m11", "ATP2A2,DYNC1I1,MAPRE2,GRK2", "^(?!.*(folding|chaperon|unfolded)).*$",
  "neuron_soma", "WGCNA_m03", "KHDRBS3,PDE1A,RYR3,RASAL1", "^(?!.*(synap|exocyt|vesicle tethering)).*$",
  "neuron_soma", "WGCNA_m06", "CAPN2,NECAB2", "^(?!.*(mitochond|respirat|ATP synthesis|oxidative)).*$",
  "microglia", "WGCNA_m02", "MPP2,ATP1A1,LIMCH1,PDE1B", "^(?!.*mitochond).*$",
  "microglia", "WGCNA_m02", "MPST,FAHD1", "mitochond",
  "microglia", "WGCNA_m02", "GCDH", "not gene-mapped",
  "microglia", "WGCNA_m07", "PEX5L,ASAP2,NLRX1", "^(?!.*(synap|postsynaptic|glutamate|actin)).*$",
  "microglia", "WGCNA_m07", "GRIA1,GRIA2,SHISA7,GIT1", "synap|postsynaptic|glutamate",
  "microglia", "WGCNA_m16", "SLC1A2,GPR37L1", "^(?!.*(respirat|respiratory|oxidative|mitochond)).*$",
  "neuron_soma", "WGCNA_m08", "ALDH2", "none",
  "neuron_neuropil", "WGCNA_m10", "CAMK1,CAMK1D,NRXN3,PLXNA1", "synap|kinase|phosphat",
  "microglia", "WGCNA_m10", "MAG,CNP,MOG,NEFH,INA", "none",
  "microglia", "WGCNA_m10", "SIRT2", "glycosyltransferase",
  "microglia", "WGCNA_m12", "CHORDC1,HMGB3,GNG2,RHOB,MYL3,MGLL", "^(?!.*(synap|vesicle)).*$",
  "microglia", "WGCNA_m17", "SLC4A4", "transporter|symporter",
  "microglia", "WGCNA_m17", "ALDH1L1", "none",
  "microglia", "WGCNA_m14", "WDR1,INF2,MYH10,WIPF3,DSTN,CAP2,KALRN", "^(?!.*actin).*$") |>
  separate_rows(genes, sep = ",") |> rename(gene = genes)
out <- claims |> rowwise() |> mutate(
  entrez = paste(unique(na.omit(mem$EntrezID[mem$dataset == dataset & mem$ModuleID == module & toupper(mem$GeneSymbol) == gene])), collapse = ";"),
  terms = paste(unique(go$Description[go$dataset == dataset & go$ModuleID == module & go$geneID %in% strsplit(entrez, ";")[[1]]]), collapse = " | "),
  in_module = any(mem$dataset == dataset & mem$ModuleID == module & toupper(mem$GeneSymbol) == gene),
  status = if (expected == "not gene-mapped") { if (in_module && !nzchar(entrez)) "ok" else "WRONG: gene-mapped or not in the module" }
  else if (!in_module) "WRONG: hub not in the module" else if (!nzchar(entrez)) "WRONG: not gene-mapped (claim untestable)" else if (expected == "none") { if (nzchar(terms)) "WRONG: in significant terms" else "ok" } else {
    t <- strsplit(terms, " | ", fixed = TRUE)[[1]]
    if (startsWith(expected, "^(?!")) { if (all(grepl(expected, t, perl = TRUE))) "ok" else "WRONG: in a signal term" }
    else if (startsWith(expected, "^(")) { if (length(t) && all(grepl(expected, t, perl = TRUE))) "ok" else "WRONG: other terms" }
    else { if (any(grepl(expected, t, perl = TRUE))) "ok" else "WRONG: not in the signal terms" } }) |> ungroup()
# module-level: no actin term at all in microglia m14
m14 <- go |> filter(dataset == "microglia", ModuleID == "WGCNA_m14") |> distinct(Description)
out <- bind_rows(out, tibble(dataset = "microglia", module = "WGCNA_m14", gene = "(module)", expected = "no significant actin term",
                             entrez = "", terms = paste(m14$Description, collapse = " | "),
                             status = if (any(grepl("actin", m14$Description))) "WRONG" else "ok"))
# hubs cited as GO support that have no Entrez ID
hub_unmapped <- tribble(~dataset, ~module, ~gene, "neuron_neuropil", "WGCNA_m03", "DDX5", "neuron_neuropil", "WGCNA_m06", "GCDH", "microglia", "WGCNA_m06", "HSPG2") |>
  rowwise() |> mutate(expected = "has Entrez ID", entrez = paste(unique(na.omit(mem$EntrezID[mem$dataset == dataset & mem$ModuleID == module & toupper(mem$GeneSymbol) == gene])), collapse = ";"),
                      terms = "", status = if (nzchar(entrez)) "mapped" else "no Entrez ID (cannot be in a GO term)") |> ungroup()
out <- bind_rows(out, hub_unmapped)
write_csv(out, file.path(ANN, "gene_claim_check.csv"))
print(out |> count(status))
print(out |> filter(!status %in% c("ok", "mapped")) |> select(dataset, module, gene, status, terms), width = 200)
