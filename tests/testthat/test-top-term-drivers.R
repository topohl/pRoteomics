# summarize_top_term_drivers.R and the three consumers repointed with it
# (2026-11-16): compare_recurrent_module_proteins.R,
# build_curated_overlap_programs.R and render_network_chord_diagram.R.

source(testthat::test_path("..", "..", "R", "paths.R"))
source(repo_path("R", "enrichment_io.R"))
source(repo_path("R", "top_term_driver_utils.R"))

TTD_SCRIPT <- "analysis/differential_abundance/summarize_top_term_drivers.R"
TTD_UTILS <- "R/enrichment/top_term_driver_utils.R"
TTD_RECURRENT <- "analysis/wgcna/compare_recurrent_module_proteins.R"
TTD_CURATED <- "analysis/wgcna/build_curated_overlap_programs.R"
TTD_CHORD <- "analysis/spatial_networks/render_network_chord_diagram.R"

ttd_code <- function(rel) {
  l <- readLines(repo_path(rel), warn = FALSE)
  paste(l[!grepl("^\\s*#", l)], collapse = "\n")
}

# Three neuron_neuropil comparisons in three spatial units. C3 has no
# significant term; T2 is redundant with T1 (Jaccard 9/11); T5 is below the
# set-size floor; T8 is never significant. Rpl3 is carried by two protein
# groups in every comparison.
TTD_COMPARISONS <- c(C1 = "CA1slmsus_CA1slmcon", C2 = "CA2slmsus_CA2slmcon", C3 = "CA3sosus_CA3socon")
TTD_UNITS <- c(C1 = "CA1_slm", C2 = "CA2_slm", C3 = "CA3_so")

ttd_fixture <- function() {
  genes <- c("Abcd1", "Ndufa1", "Rpl3", sprintf("G%02d", 4:60))
  g <- function(i) genes[i]
  sets <- list(
    T1 = g(1:10), T2 = g(c(1:9, 11)), T3 = g(12:21), T4 = g(22:31),
    T5 = g(32:35), T6 = g(36:45), T7 = g(46:55), T8 = g(56:60)
  )
  ids <- sprintf("GO:%07d", 1:8)
  desc <- c("alpha process", "beta process", "gamma process", "delta process",
            "small process", "epsilon process", "zeta process", "eta process")
  set_size <- c(12, 12, 15, 15, 8, 14, 14, 11)
  nes <- list(
    C1 = c(2.5, 2.4, 2.3, 2.2, 2.6, -2.0, -1.8, 1.0),
    C2 = c(1.2, 1.1, 2.0, 1.5, 1.3, -2.2, -0.8, 0.9),
    C3 = c(0.9, 0.8, 1.1, 1.0, 0.7, -1.0, -0.9, 0.5)
  )
  padj <- list(
    C1 = c(0.001, 0.002, 0.01, 0.02, 0.001, 0.01, 0.04, 0.3),
    C2 = c(0.4, 0.45, 0.01, 0.2, 0.3, 0.005, 0.5, 0.6),
    C3 = c(0.3, 0.35, 0.2, 0.25, 0.5, 0.2, 0.3, 0.7)
  )
  # leading edge per comparison: T1 shrinks in C2 and C3, the others are fixed
  core_of <- function(cid, t) {
    if (t == 1 && cid == "C2") return(sets$T1[1:5])
    if (t == 1 && cid == "C3") return(sets$T1[1:3])
    sets[[t]]
  }
  pg_of <- function(gene) {
    if (gene == "Rpl3") return(data.frame(pg = c("PG:t:rpl3a", "PG:t:rpl3b"),
                                          acc = c("Q00001", "Q00003;P00002")))
    if (gene == "Abcd1") return(data.frame(pg = "PG:t:abcd1", acc = "P03888"))
    data.frame(pg = paste0("PG:t:", tolower(gene)), acc = paste0("A", toupper(gene)))
  }
  terms <- list(); prov <- list()
  for (cid in names(TTD_COMPARISONS)) {
    comp <- TTD_COMPARISONS[[cid]]
    for (t in seq_along(ids)) {
      core <- core_of(cid, t)
      terms[[length(terms) + 1L]] <- data.frame(
        ID = ids[t], Description = desc[t], setSize = set_size[t],
        NES = nes[[cid]][t], p.adjust = padj[[cid]][t],
        core_enrichment = paste(core, collapse = "/"),
        term_id = ids[t], term_description = desc[t], dataset = "neuron_neuropil",
        comparison = comp, result_type = "GSEA_GO", ontology = "BP",
        stringsAsFactors = FALSE)
      for (gene in core) {
        p <- pg_of(gene)
        prov[[length(prov) + 1L]] <- data.frame(
          dataset = "neuron_neuropil", comparison = comp, result_type = "GSEA_GO",
          ontology = "BP", term_id = ids[t], term_description = desc[t],
          official_gene_symbol = gene, official_entrez_id = "1000",
          ProteinGroupID = p$pg, member_accessions = p$acc,
          protein_group_gene_annotation_status = "concordant_official_gene",
          gene_level_claim_allowed = TRUE, rank_statistic = 1.5,
          core_enrichment_member = TRUE,
          enrichment_contract_version = "protein_group_enrichment_v3_term_gene_provenance",
          gene_annotation_contract_version = "mouse_gene_annotation_v1",
          stringsAsFactors = FALSE)
      }
    }
    # a ranked but non-core gene of T1, which must never count
    prov[[length(prov) + 1L]] <- data.frame(
      dataset = "neuron_neuropil", comparison = comp, result_type = "GSEA_GO",
      ontology = "BP", term_id = ids[1], term_description = desc[1],
      official_gene_symbol = "Notcore1", official_entrez_id = "1001",
      ProteinGroupID = "PG:t:notcore1", member_accessions = "ZNOT01",
      protein_group_gene_annotation_status = "concordant_official_gene",
      gene_level_claim_allowed = TRUE, rank_statistic = -0.2,
      core_enrichment_member = FALSE,
      enrichment_contract_version = "protein_group_enrichment_v3_term_gene_provenance",
      gene_annotation_contract_version = "mouse_gene_annotation_v1",
      stringsAsFactors = FALSE)
  }
  list(terms = do.call(rbind, terms), provenance = do.call(rbind, prov))
}

ttd_run <- function(fx, params = top_term_driver_parameters()) {
  term_table <- top_term_table(fx$terms)
  core <- top_term_core_rows(fx$provenance)
  best <- top_term_best_instances(term_table)
  selected <- top_term_selection(best, top_term_leading_edges(core), params)
  drivers <- top_term_driver_genes(term_table, core, unique(selected$Description), params)
  list(selected = selected, drivers = drivers)
}

# ---- selection ------------------------------------------------------------

testthat::test_that("parameters come from the compareGO config with the former defaults", {
  p <- top_term_driver_parameters()
  testthat::expect_identical(p$target_n_terms, 5L)
  testthat::expect_identical(p$redundancy_threshold, 0.7)
  testthat::expect_identical(p$min_set_size, 10L)
  testthat::expect_identical(p$top_driver_genes_per_term, 100L)
  testthat::expect_true(p$significant_only)
  cfg <- yaml::read_yaml(repo_path("config", "compareGO_config.yml"))
  q <- top_term_driver_parameters(cfg)
  testthat::expect_identical(q$target_n_terms, as.integer(cfg$target_n_terms))
  testthat::expect_identical(q$redundancy_threshold, as.numeric(cfg$redundancy_threshold))
  testthat::expect_identical(q$min_set_size, as.integer(cfg$min_set_size))
  testthat::expect_identical(top_term_driver_parameters(list(target_n_terms = 2))$target_n_terms, 2L)
  testthat::expect_error(top_term_driver_parameters(list(min_set_size = NA)), "Invalid")
})

testthat::test_that("the Jaccard helper is local and treats an empty set as no overlap", {
  testthat::expect_identical(top_term_jaccard(character(), c("a", "b")), 0)
  testthat::expect_equal(top_term_jaccard(c("a", "b", "c"), c("b", "c", "d")), 0.5)
  testthat::expect_equal(top_term_jaccard(c("a", "a", "b"), c("a", "b")), 1)
})

testthat::test_that("the best instance per term and comparison has the largest |NES|", {
  tt <- data.frame(comparison = c("c", "c", "c"), term_id = c("GO:2", "GO:1", "GO:3"),
                   Description = c("x", "x", "y"), NES = c(-2.1, 1.5, 1), p.adjust = 0.01,
                   setSize = 20L, stringsAsFactors = FALSE)
  b <- top_term_best_instances(tt)
  testthat::expect_identical(nrow(b), 2L)
  testthat::expect_identical(b$term_id[b$Description == "x"], "GO:2")
})

testthat::test_that("redundant and small terms are skipped, up and down are kept separately", {
  out <- ttd_run(ttd_fixture())$selected
  c1 <- out[out$comparison == TTD_COMPARISONS[["C1"]], ]
  testthat::expect_identical(c1$Description[c1$direction == "up"],
                             c("alpha process", "gamma process", "delta process"))
  testthat::expect_identical(c1$Description[c1$direction == "down"],
                             c("epsilon process", "zeta process"))
  testthat::expect_identical(c1$selection_rank[c1$direction == "up"], 1:3)
  # beta is redundant with alpha, small is below min_set_size, eta is never significant
  testthat::expect_false(any(c("beta process", "small process", "eta process") %in% out$Description))
  testthat::expect_true(all(c1$max_jaccard_with_kept[-1] <= 0.7, na.rm = TRUE))
  # only significant best instances form the pool; C3 has none
  c2 <- out[out$comparison == TTD_COMPARISONS[["C2"]], ]
  testthat::expect_setequal(c2$Description, c("gamma process", "epsilon process"))
  testthat::expect_false(TTD_COMPARISONS[["C3"]] %in% out$comparison)
})

testthat::test_that("target_n_terms caps each direction, NES order decides", {
  out <- ttd_run(ttd_fixture(), top_term_driver_parameters(list(target_n_terms = 2)))$selected
  c1 <- out[out$comparison == TTD_COMPARISONS[["C1"]], ]
  testthat::expect_identical(c1$Description[c1$direction == "up"], c("alpha process", "gamma process"))
  testthat::expect_identical(c1$Description[c1$direction == "down"], c("epsilon process", "zeta process"))
  # a looser threshold lets the redundant term back in
  loose <- ttd_run(ttd_fixture(), top_term_driver_parameters(list(target_n_terms = 2, redundancy_threshold = 0.9)))$selected
  l1 <- loose[loose$comparison == TTD_COMPARISONS[["C1"]] & loose$direction == "up", ]
  testthat::expect_identical(l1$Description, c("alpha process", "beta process"))
})

# ---- driver table -----------------------------------------------------------

testthat::test_that("Freq, Mean_NES, n_significant_comparisons and the per-term cap", {
  d <- ttd_run(ttd_fixture(), top_term_driver_parameters(list(top_driver_genes_per_term = 4)))$drivers
  testthat::expect_identical(names(d), c("Description", "term_id", "Gene", "ProteinGroupIDs",
    "member_accessions", "Freq", "n_significant_comparisons", "Mean_NES", "Comparisons"))
  a <- d[d$Description == "alpha process", ]
  # all comparisons count, significant or not; ties go to the gene symbol
  testthat::expect_identical(a$Gene, c("Abcd1", "Ndufa1", "Rpl3", "G04"))
  testthat::expect_identical(a$Freq, c(3L, 3L, 3L, 2L))
  testthat::expect_identical(a$n_significant_comparisons, c(1L, 1L, 1L, 1L))
  testthat::expect_equal(a$Mean_NES[1], mean(c(2.5, 1.2, 0.9)))
  testthat::expect_equal(a$Mean_NES[4], mean(c(2.5, 1.2)))
  testthat::expect_identical(a$Comparisons[1], paste(sort(unname(TTD_COMPARISONS), method = "radix"), collapse = "; "))
  # negative terms use |NES|; epsilon is significant in C1 and C2
  e <- d[d$Description == "epsilon process", ]
  testthat::expect_true(all(e$Mean_NES > 0))
  testthat::expect_true(all(e$n_significant_comparisons == 2L))
  testthat::expect_true(all(table(d$Description) <= 4L))
  testthat::expect_false("Notcore1" %in% d$Gene)
  testthat::expect_setequal(unique(d$Description),
    c("alpha process", "gamma process", "delta process", "epsilon process", "zeta process"))
})

testthat::test_that("a gene with two protein groups counts once per comparison and keeps both groups", {
  d <- ttd_run(ttd_fixture())$drivers
  r <- d[d$Description == "alpha process" & d$Gene == "Rpl3", ]
  testthat::expect_identical(nrow(r), 1L)
  testthat::expect_identical(r$Freq, 3L)
  testthat::expect_identical(r$ProteinGroupIDs, "PG:t:rpl3a;PG:t:rpl3b")
  testthat::expect_identical(r$member_accessions, "Q00001 | Q00003;P00002")
  pairs <- top_term_protein_group_pairs(r$ProteinGroupIDs, r$member_accessions)
  testthat::expect_identical(pairs$ProteinGroupID, c("PG:t:rpl3a", "PG:t:rpl3b"))
  testthat::expect_identical(top_term_representative_accession(pairs$member_accessions), c("Q00001", "P00002"))
  merged <- top_term_merge_protein_groups(c(r$ProteinGroupIDs, "PG:t:rpl3a"), c(r$member_accessions, "Q00001"))
  testthat::expect_identical(unname(merged[["ProteinGroupIDs"]]), r$ProteinGroupIDs)
  testthat::expect_error(top_term_protein_group_pairs("PG:1", "A | B"), "more protein groups")
})

testthat::test_that("the selection and the driver table do not depend on input order", {
  fx <- ttd_fixture()
  a <- ttd_run(fx)
  set.seed(20261116)
  fx$terms <- fx$terms[sample(nrow(fx$terms)), ]
  fx$provenance <- fx$provenance[sample(nrow(fx$provenance)), ]
  b <- ttd_run(fx)
  rownames(a$selected) <- rownames(b$selected) <- NULL
  testthat::expect_identical(a$selected[order(a$selected$comparison, a$selected$direction, a$selected$selection_rank), ],
                             b$selected[order(b$selected$comparison, b$selected$direction, b$selected$selection_rank), ])
  testthat::expect_identical(a$drivers, b$drivers)
})

# ---- canonical bundle -------------------------------------------------------

ttd_write_bundle <- function(root, fx = ttd_fixture()) {
  fdir <- file.path(root, "fixture")
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  term_file <- file.path(fdir, "compareGO_term_comparison.csv")
  prov_file <- file.path(fdir, "compareGO_term_gene_provenance.csv")
  status_file <- file.path(fdir, "compareGO_analysis_status_summary.csv")
  cp_file <- file.path(fdir, "clusterProfiler_manifest.csv")
  utils::write.csv(fx$terms, term_file, row.names = FALSE)
  utils::write.csv(fx$provenance, prov_file, row.names = FALSE)
  utils::write.csv(data.frame(dataset = "neuron_neuropil", comparison = unname(TTD_COMPARISONS),
    result_type = "GSEA_GO", ontology = "BP", analysis_status = "success_with_terms",
    n_terms = 8L, comparego_action = "included"), status_file, row.names = FALSE)
  writeLines("x", cp_file)
  manifest <- data.frame(
    dataset = "neuron_neuropil", comparison = unname(TTD_COMPARISONS), result_type = "GSEA_GO",
    ontology = "BP", route_category = "phenotype_within_unit", route_unit = unname(TTD_UNITS),
    analysis_status = "success_with_terms", comparego_analysis_status = "included",
    input_manifest = cp_file, term_comparison_file = term_file,
    term_gene_provenance_output_file = prov_file, analysis_status_summary_file = status_file,
    enrichment_contract_version = canonical_clusterprofiler_manifest_contract_version(),
    comparego_contract_version = canonical_comparego_manifest_contract_version(),
    stringsAsFactors = FALSE)
  mdir <- file.path(root, "results", "differential_abundance", "compare_go_enrichment",
                    "neuron_neuropil", "models")
  dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(manifest, file.path(mdir, "compareGO_input_manifest.csv"), row.names = FALSE)
  file.path(mdir, "compareGO_input_manifest.csv")
}

testthat::test_that("the bundle reader can read a provenance column subset under the same contract", {
  root <- withr::local_tempdir("ttd_bundle_")
  mf <- ttd_write_bundle(root)
  full <- read_canonical_comparego_bundle(mf, "neuron_neuropil", repository_root = root)
  part <- read_canonical_comparego_bundle(mf, "neuron_neuropil", repository_root = root,
    provenance_columns = c("member_accessions", "core_enrichment_member"))
  testthat::expect_setequal(names(full$provenance), term_gene_provenance_columns())
  testthat::expect_setequal(names(part$provenance),
    c(term_gene_provenance_key_columns(), "member_accessions", "core_enrichment_member"))
  testthat::expect_identical(part$provenance, full$provenance[names(part$provenance)])
  testthat::expect_identical(part$terms, full$terms)
  # the strict checks still run on what was read
  fx <- ttd_fixture()
  fx$provenance$gene_level_claim_allowed[1] <- FALSE
  bad <- ttd_write_bundle(withr::local_tempdir("ttd_bad_"), fx)
  testthat::expect_error(read_canonical_comparego_bundle(bad, "neuron_neuropil",
    provenance_columns = "core_enrichment_member"), "Ineligible")
})

# ---- static -----------------------------------------------------------------

testthat::test_that("the producer reads the canonical bundle and owns its rules", {
  script <- ttd_code(TTD_SCRIPT)
  utils_txt <- ttd_code(TTD_UTILS)
  both <- paste(script, utils_txt)
  testthat::expect_match(script, "canonical_comparego_manifest_path(DATASET)", fixed = TRUE)
  testthat::expect_match(script, "read_canonical_comparego_bundle(", fixed = TRUE)
  testthat::expect_match(script, "current_dataset_from_cli(", fixed = TRUE)
  testthat::expect_match(both, "gene_level_claim_allowed", fixed = TRUE)
  testthat::expect_match(both, "core_enrichment_member", fixed = TRUE)
  testthat::expect_match(utils_txt, "top_term_jaccard <- function", fixed = TRUE)
  for (txt in c(paste(readLines(repo_path(TTD_SCRIPT), warn = FALSE), collapse = "\n"),
                paste(readLines(repo_path(TTD_UTILS), warn = FALSE), collapse = "\n"))) {
    testthat::expect_false(grepl("07_Top_Genes_Driving_TopTerms", txt, fixed = TRUE))
  }
  # nothing sources or parses the archive
  for (rel in c(TTD_SCRIPT, TTD_UTILS)) {
    pd <- utils::getParseData(parse(repo_path(rel), keep.source = TRUE))
    calls <- pd$text[pd$token == "SYMBOL_FUNCTION_CALL"]
    testthat::expect_false(any(calls %in% c("sys.source", "eval", "evalq")), info = rel)
    src <- pd$text[pd$token == "STR_CONST"]
    testthat::expect_false(any(grepl("archive/|\"archive\"", src)), info = rel)
  }
})

testthat::test_that("the three consumers no longer name their legacy inputs", {
  rec <- ttd_code(TTD_RECURRENT)
  testthat::expect_false(grepl("07_Top_Genes_Driving_TopTerms", rec, fixed = TRUE))
  testthat::expect_false(grepl("04_differential_expression_enrichment", rec, fixed = TRUE))
  testthat::expect_match(rec, "top_term_driver_genes.csv", fixed = TRUE)
  testthat::expect_false(grepl("readxl::read_excel", rec, fixed = TRUE))

  cur <- ttd_code(TTD_CURATED)
  testthat::expect_false(grepl("overlapping_proteins_min3_datasets", cur, fixed = TRUE))
  testthat::expect_false(grepl("04_differential_expression_enrichment", cur, fixed = TRUE))
  testthat::expect_match(cur, "overlapping_proteins_min3_spatial_units.xlsx", fixed = TRUE)
  testthat::expect_match(cur, "Overlap_summary", fixed = TRUE)

  cho <- ttd_code(TTD_CHORD)
  testthat::expect_false(grepl("neuron-phenotypeWithinUnit", cho, fixed = TRUE))
  testthat::expect_false(grepl("02_id_mapping", cho, fixed = TRUE))
  testthat::expect_match(cho, "preprocessing_mapped_contrast_dir(", fixed = TRUE)
  testthat::expect_false(grepl('c("gene_symbol"', cho, fixed = TRUE))

  testthat::skip_if_not_installed("yaml")
  reg <- yaml::read_yaml(repo_path("pipeline.yml"))
  rows <- unlist(lapply(reg$stages, `[[`, "scripts"), recursive = FALSE)
  row_of <- function(s) Filter(function(x) identical(x$script, s), rows)[[1]]
  legacy <- c("07_Top_Genes_Driving_TopTerms", "overlapping_proteins_min3_datasets",
              "neuron-phenotypeWithinUnit", "results/tables/07_spatial_networks")
  for (s in c(TTD_RECURRENT, TTD_CURATED, TTD_CHORD)) {
    consumed <- c(unlist(row_of(s)$consumes_required), unlist(row_of(s)$consumes_optional))
    for (l in legacy) testthat::expect_false(any(grepl(l, consumed, fixed = TRUE)), info = paste(s, l))
  }
  testthat::expect_true("data/external/MOUSE_10090_idmapping.dat" %in% unlist(row_of(TTD_CURATED)$consumes_required))
  new_row <- row_of(TTD_SCRIPT)
  testthat::expect_false(isTRUE(new_row$required))
  testthat::expect_false(isTRUE(new_row$recomputes_core_state))
  scripts <- vapply(reg$stages$enrichment$scripts, `[[`, character(1), "script")
  testthat::expect_identical(match(TTD_SCRIPT, scripts),
    match("analysis/differential_abundance/compare_go_enrichment.R", scripts) + 1L)
})

# ---- end to end, in an isolated project root --------------------------------

ttd_project_root <- function(prefix) {
  root <- withr::local_tempdir(prefix, .local_envir = parent.frame())
  file.copy(repo_path("R"), root, recursive = TRUE)
  file.copy(repo_path("config"), root, recursive = TRUE)
  root
}

ttd_rscript <- function(root, rel, args = character()) {
  ## resolved here, before the project root is redirected
  script <- shQuote(repo_path(rel))
  out <- withr::with_dir(root, withr::with_envvar(
    c(PROTEOMICS_PROJECT_ROOT = root,
      PROTEOMICS_INPUT_RESOLUTION_AUDIT = file.path(root, "ledger.csv")),
    suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      c(script, args), stdout = TRUE, stderr = TRUE))))
  list(status = as.integer(attr(out, "status") %||% 0L), output = paste(out, collapse = "\n"))
}

testthat::test_that("producer, recurrent proteins and curated programs run end to end", {
  for (p in c("yaml", "readr", "readxl", "openxlsx", "dplyr", "tidyr", "purrr", "stringr", "tibble", "withr"))
    testthat::skip_if_not_installed(p)
  root <- ttd_project_root("ttd_chain_")
  ttd_write_bundle(root)
  ext <- file.path(root, "data", "external")
  dir.create(ext, recursive = TRUE)
  writeLines(c("P03888\tUniProtKB-ID\tNU1M_MOUSE", "P03888\tGO\tGO:0005739",
               "Q00001\tUniProtKB-ID\tRL3_MOUSE", "P00002\tUniProtKB-ID\tRL3B_MOUSE"),
             file.path(ext, "MOUSE_10090_idmapping.dat"))

  snapshot <- function() sort(list.files(root, recursive = TRUE, all.files = TRUE, include.dirs = TRUE))
  before <- snapshot()
  dry <- ttd_rscript(root, TTD_SCRIPT, c("--dataset", "neuron_neuropil", "--dry-run"))
  testthat::expect_identical(dry$status, 0L, info = dry$output)
  testthat::expect_identical(setdiff(snapshot(), before), character(0))

  run <- ttd_rscript(root, TTD_SCRIPT, c("--dataset", "neuron_neuropil"))
  testthat::expect_identical(run$status, 0L, info = run$output)
  tables <- file.path(root, "results", "differential_abundance", "summarize_top_term_drivers",
                      "neuron_neuropil", "tables")
  drivers <- utils::read.csv(file.path(tables, "top_term_driver_genes.csv"), stringsAsFactors = FALSE)
  selected <- utils::read.csv(file.path(tables, "selected_top_terms.csv"), stringsAsFactors = FALSE)
  expected <- ttd_run(ttd_fixture())
  testthat::expect_identical(nrow(drivers), nrow(expected$drivers))
  testthat::expect_identical(drivers$Gene, expected$drivers$Gene)
  testthat::expect_identical(drivers$Freq, expected$drivers$Freq)
  testthat::expect_identical(nrow(selected), nrow(expected$selected))
  testthat::expect_true(all(selected$route_unit %in% TTD_UNITS))
  testthat::expect_true(file.exists(file.path(root, "results", "differential_abundance",
    "summarize_top_term_drivers", "neuron_neuropil", "manifests", "run_manifest.yml")))
  testthat::expect_identical(list.files(tables, pattern = "[.]tmp-"), character(0))

  rec <- ttd_rscript(root, TTD_RECURRENT, c("--dataset", "neuron_neuropil"))
  testthat::expect_identical(rec$status, 0L, info = rec$output)
  wb <- file.path(root, "results", "wgcna", "compare_recurrent_module_proteins", "neuron_neuropil",
                  "tables", "recurrent_top_gene_drivers", "overlapping_proteins_min3_spatial_units.xlsx")
  ov <- readxl::read_excel(wb, sheet = "Overlap_summary")
  # each comparison is its own unit, so a gene is recurrent when it drives a
  # term in all three comparisons
  testthat::expect_setequal(ov$Gene, unique(expected$drivers$Gene[expected$drivers$Freq == 3L]))
  testthat::expect_true(all(c("Abcd1", "Ndufa1", "Rpl3") %in% ov$Gene))
  testthat::expect_false("G04" %in% ov$Gene)
  testthat::expect_true(all(ov$N_datasets == 3))
  testthat::expect_identical(ov$ProteinGroupIDs[ov$Gene == "Rpl3"], "PG:t:rpl3a;PG:t:rpl3b")
  testthat::expect_identical(ov$member_accessions[ov$Gene == "Rpl3"], "Q00001 | Q00003;P00002")

  cur <- ttd_rscript(root, TTD_CURATED)
  testthat::expect_identical(cur$status, 0L, info = cur$output)
  ml <- readxl::read_excel(file.path(root, "results", "wgcna", "build_curated_overlap_programs",
                                     "global", "tables", "curated_overlap_programs.xlsx"),
                           sheet = "Modules_long")
  ribo <- ml[ml$ModuleID == "Neuropil_ribosome_translation", ]
  testthat::expect_setequal(ribo$UniProt, c("Q00001", "P00002"))
  testthat::expect_setequal(ribo$ProteinGroupID, c("PG:t:rpl3a", "PG:t:rpl3b"))
  testthat::expect_identical(ribo$member_accessions[ribo$UniProt == "P00002"], "Q00003;P00002")
  mito <- ml[ml$ModuleID == "Neuropil_mito_bioenergetics", ]
  # Abcd1 is mitochondrial only through its member accession
  testthat::expect_setequal(mito$GeneSymbol, c("Abcd1", "Ndufa1"))
  testthat::expect_true(all(c("ModuleSet", "ModuleID", "Module", "UniProt", "GeneSymbol", "Source") %in% names(ml)))
  testthat::expect_false("ProteinID" %in% names(ml))
})

testthat::test_that("the chord renderer joins on ProteinGroupID and writes a note for an empty direction", {
  for (p in c("circlize", "magick", "RColorBrewer", "openxlsx", "readxl", "digest", "dplyr", "tidyr", "withr"))
    testthat::skip_if_not_installed(p)
  root <- ttd_project_root("ttd_chord_")
  pf <- file.path(root, "results", "preprocessing", "map_protein_identifiers", "neuron_neuropil",
                  "tables", "mapped", "forward", "per_file")
  dir.create(pf, recursive = TRUE)
  # PG1 and PG2 are two protein groups of one gene, and gene_symbol (display
  # only) repeats on a third group; PG5 may not carry a protein-level claim
  contrast <- function(padj, log2fc) data.frame(
    ProteinGroupID = c("PG1", "PG2", "PG3", "PG4", "PG5"),
    official_gene_symbol = c("Shared1", "Shared1", "Other3", "Diff4", "Amb5"),
    gene_symbol = c("Shared1", "Shared1", "Other3", "Shared1", "Amb5"),
    protein_level_claim_allowed = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    padj = padj, log2fc = log2fc, stringsAsFactors = FALSE)
  utils::write.csv(contrast(c(0.01, 0.02, 0.01, 0.5, 0.01), c(1, 1.2, -1, 0.3, 1)),
                   file.path(pf, "CA1slmsus_CA1slmcon.csv"), row.names = FALSE)
  utils::write.csv(contrast(c(0.01, 0.03, 0.6, 0.01, 0.01), c(0.8, 0.9, -0.2, 1, 1)),
                   file.path(pf, "CA2slmsus_CA2slmcon.csv"), row.names = FALSE)
  # a contrast between two units is not phenotype-within-unit and is skipped
  utils::write.csv(contrast(rep(0.001, 5), rep(2, 5)),
                   file.path(pf, "CA1slmsus_CA2slmcon.csv"), row.names = FALSE)

  dry <- ttd_rscript(root, TTD_CHORD, c("--dataset", "neuron_neuropil", "--dry-run"))
  testthat::expect_identical(dry$status, 0L, info = dry$output)
  testthat::expect_match(dry$output, "Phenotype-within-unit CSV count: 2", fixed = TRUE)
  testthat::expect_false(dir.exists(file.path(root, "results", "spatial_networks")))

  run <- ttd_rscript(root, TTD_CHORD, c("--dataset", "neuron_neuropil"))
  testthat::expect_identical(run$status, 0L, info = run$output)
  base <- file.path(root, "results", "spatial_networks", "render_network_chord_diagram", "neuron_neuropil")
  up <- readxl::read_excel(file.path(base, "tables", "shared_proteins_up_summary.xlsx"), sheet = "summary")
  details <- readxl::read_excel(file.path(base, "tables", "shared_proteins_up_summary.xlsx"), sheet = "details")
  testthat::expect_identical(sort(up$ProteinGroupID), c("PG1", "PG2"))
  testthat::expect_identical(unique(up$official_gene_symbol), "Shared1")
  testthat::expect_identical(nrow(details), 2L)
  testthat::expect_setequal(c(details$comparison_a, details$comparison_b),
                            c("CA1slmsus_CA1slmcon", "CA2slmsus_CA2slmcon"))
  testthat::expect_true(file.exists(file.path(base, "plots", "chord_upregulated.svg")))
  # nothing shared downwards: a note, no figure, no side-by-side panel
  testthat::expect_true(file.exists(file.path(base, "plots", "chord_downregulated_note.txt")))
  testthat::expect_false(file.exists(file.path(base, "plots", "chord_downregulated.svg")))
  testthat::expect_false(file.exists(file.path(base, "plots", "chord_up_down_side_by_side.png")))
  down <- readxl::read_excel(file.path(base, "tables", "shared_proteins_down_summary.xlsx"), sheet = "summary")
  testthat::expect_identical(nrow(down), 0L)
  testthat::expect_identical(list.files(file.path(base, "plots"), pattern = "[.]xlsx$"), character(0))
})
