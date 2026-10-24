testthat::local_edition(3)

source(testthat::test_path("..", "..", "R", "paths.R"))

# The script's helpers, without running it (its entrypoint is guarded by sys.nframe()).
load_overlap_helpers <- function() {
  env <- new.env(parent = globalenv())
  old_wd <- setwd(repo_path())
  on.exit(setwd(old_wd), add = TRUE)
  sys.source(repo_path("analysis/wgcna", "compare_module_enrichment_overlap.R"), envir = env)
  env
}

pg <- function(x, dataset = "neuron_soma") sprintf("PG:%s:%s", dataset, x)

universe_fixture <- function() {
  data.frame(
    ProteinGroupID = pg(c("A1", "B2", "C3", "D4", "E5")),
    included_in_wgcna = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    ModuleID = c("WGCNA_m01", "WGCNA_m01", "WGCNA_m02", "WGCNA_m02", "WGCNA_m02"),
    member_accessions = c("P11111", "P22222", "P33333", "P44444", "P55555"),
    official_gene_symbol = c("Syp", "Snap25", "Gfap", NA, "Aqp4"),
    stringsAsFactors = FALSE
  )
}

write_fixture <- function(df, dir, name) {
  path <- file.path(dir, name)
  utils::write.csv(df, path, row.names = FALSE)
  path
}

testthat::test_that("DE set is the significant ProteinGroupIDs of the DA table", {
  helpers <- load_overlap_helpers()
  dir <- withr::local_tempdir()
  da <- write_fixture(data.frame(
    ProteinGroupID = pg(c("A1", "B2", "C3")),
    official_gene_symbol = c("Syp", "Snap25", "Gfap"),
    padj = c(0.01, 0.5, 0.04),
    significant = c(TRUE, FALSE, TRUE)
  ), dir, "da.csv")
  testthat::expect_setequal(helpers$read_de_set(da), pg(c("A1", "C3")))
  testthat::expect_identical(helpers$read_de_set(NA_character_), character())

  no_flag <- write_fixture(data.frame(ProteinGroupID = pg("A1"), padj = 0.01), dir, "no_flag.csv")
  testthat::expect_error(helpers$read_de_set(no_flag), "missing column")
  testthat::expect_error(helpers$read_de_set(file.path(dir, "absent.csv")), "not readable")
})

testthat::test_that("leading edge keeps core members of one comparison's significant terms", {
  helpers <- load_overlap_helpers()
  dir <- withr::local_tempdir()
  provenance <- write_fixture(data.frame(
    comparison = c("X", "X", "X", "X", "Y"),
    term_id = c("GO:1", "GO:1", "GO:2", "GO:1", "GO:1"),
    official_gene_symbol = c("Syp", "Snap25", "Gfap", "Aqp4", "Gfap"),
    ProteinGroupID = pg(c("A1", "B2", "C3", "E5", "C3")),
    core_enrichment_member = c(TRUE, TRUE, TRUE, FALSE, TRUE),
    gene_level_claim_allowed = TRUE
  ), dir, "provenance.csv")

  # clusterProfiler: one term table per comparison, no comparison column.
  per_run <- write_fixture(data.frame(ID = c("GO:1", "GO:2"), p.adjust = c(0.01, 0.2)), dir, "per_run.csv")
  testthat::expect_setequal(helpers$read_leading_edge_set(per_run, provenance, "X"), pg(c("A1", "B2")))

  # compareGO: pooled table; GO:2 is significant only in another comparison.
  pooled <- write_fixture(data.frame(
    comparison = c("X", "X", "Y", "Y"), ID = c("GO:1", "GO:2", "GO:1", "GO:2"),
    p.adjust = c(0.01, 0.2, 0.5, 0.001)
  ), dir, "pooled.csv")
  testthat::expect_setequal(helpers$read_leading_edge_set(pooled, provenance, "X"), pg(c("A1", "B2")))
  testthat::expect_identical(helpers$read_leading_edge_set(pooled, provenance, "Y"), character())
})

testthat::test_that("sets map to the WGCNA universe by exact ProteinGroupID and fail closed across datasets", {
  helpers <- load_overlap_helpers()
  bundle <- gww_build_universe_bundle(universe_fixture(), "neuron_soma")
  testthat::expect_setequal(
    helpers$map_to_wgcna_universe(pg(c("A1", "E5", "Z9")), bundle, "neuron_soma"),
    pg("A1")
  )
  testthat::expect_identical(helpers$map_to_wgcna_universe(pg(c("E5", "Z9")), bundle, "neuron_soma"), character())
  testthat::expect_identical(helpers$map_to_wgcna_universe(character(), bundle, "neuron_soma"), character())
  testthat::expect_error(
    helpers$map_to_wgcna_universe(c(pg("A1"), pg("A1", "microglia")), bundle, "neuron_soma"),
    "outside the PG:neuron_soma: namespace"
  )
  # A gene symbol is not an identity here.
  testthat::expect_error(helpers$map_to_wgcna_universe(c(pg("A1"), "Syp"), bundle, "neuron_soma"), "namespace")
})

testthat::test_that("overlap statistics label overlap proteins with their gene symbols", {
  helpers <- load_overlap_helpers()
  universe <- pg(c("A1", "B2", "C3", "D4"))
  labels <- stats::setNames(c("Syp", "Snap25", "Gfap", NA), universe)
  stats <- helpers$overlap_stats(pg(c("C3", "D4")), pg(c("C3", "D4", "A1")), universe, labels = labels)
  testthat::expect_identical(stats$n_overlap, 2L)
  testthat::expect_identical(stats$overlap_proteins, paste(pg(c("C3", "D4")), collapse = ";"))
  testthat::expect_identical(stats$overlap_genes, paste(c("Gfap", pg("D4")), collapse = ";"))
})

testthat::test_that("a run on fixtures overlaps on ProteinGroupID and skips compareGO re-collations", {
  helpers <- load_overlap_helpers()
  dir <- withr::local_tempdir()
  # Every path the run reads or writes resolves inside the temp directory.
  helpers$wgcna_dirs <- function(...) list(tables = dir, logs = dir)
  helpers$wgcna_modules_dir <- function(...) dir
  helpers$wgcna_modules_artifact <- function(filename, ...) file.path(dir, filename)
  helpers$path_results <- function(...) file.path(dir, utils::tail(c(...), 1L))
  helpers$validate_wgcna_module_definitions <- function(...) invisible(TRUE)
  # the fixture manifest carries only the fields this step reads, not the full clusterProfiler contract
  helpers$read_canonical_clusterprofiler_manifest <- function(path, dataset) readr::read_csv(path, show_col_types = FALSE)

  write_fixture(universe_fixture(), dir, "WGCNA_feature_universe.csv")
  write_fixture(data.frame(
    ModuleID = rep(c("WGCNA_m01", "WGCNA_m02"), each = 2), ModuleColor = rep(c("#486A8A", "#7FA6C1"), each = 2),
    ProteinGroupID = pg(c("A1", "B2", "C3", "D4"))
  ), dir, "WGCNA_module_definitions_for_downstream.csv")
  da <- write_fixture(data.frame(ProteinGroupID = pg(c("A1", "B2", "C3", "D4")), significant = c(TRUE, FALSE, FALSE, FALSE)), dir, "da_x.csv")
  da_none <- write_fixture(data.frame(ProteinGroupID = pg(c("A1", "B2")), significant = FALSE), dir, "da_z.csv")
  go <- write_fixture(data.frame(ID = "GO:1", p.adjust = 0.01), dir, "go.csv")
  kegg <- write_fixture(data.frame(ID = "mmu1", p.adjust = 0.5), dir, "kegg.csv")
  go_prov <- write_fixture(data.frame(
    comparison = "X", term_id = "GO:1", ProteinGroupID = pg(c("C3", "D4")),
    core_enrichment_member = TRUE, gene_level_claim_allowed = TRUE
  ), dir, "go_prov.csv")
  kegg_prov <- write_fixture(data.frame(
    comparison = "X", term_id = "mmu1", ProteinGroupID = pg("A1"),
    core_enrichment_member = TRUE, gene_level_claim_allowed = TRUE
  ), dir, "kegg_prov.csv")
  manifest <- data.frame(
    analysis_id = c("go_X", "kegg_X", "go_Z"), dataset = "neuron_soma",
    contrast = c("X", "X", "Z"), comparison = c("X", "X", "Z"), route_category = "phenotype_within_unit",
    result_type = c("GSEA_GO", "GSEA_KEGG", "GSEA_GO"), analysis_status = "success_with_terms",
    input_gene_file = c(da, da, da_none), output_table = c(go, kegg, go),
    term_gene_provenance_file = c(go_prov, kegg_prov, go_prov)
  )
  write_fixture(manifest, dir, "clusterProfiler_manifest.csv")
  write_fixture(manifest[1, ], dir, "compareGO_input_manifest.csv")

  helpers$run_wgcna_de_gsea_overlap("neuron_soma", dry_run = FALSE)
  rows <- readr::read_csv(file.path(dir, "WGCNA_vs_DE_GSEA_overlap.csv"), show_col_types = FALSE)
  testthat::expect_identical(nrow(rows), 4L)
  testthat::expect_true(all(rows$route == "clusterProfiler"))
  testthat::expect_identical(unique(rows$n_universe), 4)
  go_m02 <- rows[rows$ModuleID == "WGCNA_m02" & rows$result_type == "GSEA_GO", ]
  testthat::expect_identical(go_m02$n_leading_edge_overlap, 2)
  testthat::expect_identical(go_m02$top_overlap_genes, paste(c("Gfap", pg("D4")), collapse = ";"))
  kegg_m01 <- rows[rows$ModuleID == "WGCNA_m01" & rows$result_type == "GSEA_KEGG", ]
  testthat::expect_identical(c(kegg_m01$n_DE_overlap, kegg_m01$n_leading_edge_overlap), c(1, 0))
  testthat::expect_identical(kegg_m01$input_status, "DE_only_available")
  status <- readr::read_csv(file.path(dir, "WGCNA_vs_DE_GSEA_overlap_status.csv"), show_col_types = FALSE)
  testthat::expect_match(status$input_availability, "no_significant_set=1", fixed = TRUE)
})

testthat::test_that("the script keys modules on ProteinGroupID, not a coalesced accession or symbol", {
  txt <- paste(readLines(repo_path("analysis/wgcna", "compare_module_enrichment_overlap.R"), warn = FALSE), collapse = "\n")
  testthat::expect_false(grepl("coalesce(.data$UniProt", txt, fixed = TRUE))
  testthat::expect_false(grepl("tokenize_proteins", txt, fixed = TRUE))
  testthat::expect_true(grepl("gww_build_universe_bundle(universe_raw, dataset)", txt, fixed = TRUE))
})
