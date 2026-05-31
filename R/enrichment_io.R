# Shared enrichment IO and conservative biological-program mapping helpers.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}
if (!exists("safe_name", mode = "function")) {
  source(repo_path("R", "validation_utils.R"))
}

biological_program_patterns <- function() {
  data.frame(
    biological_program = c(
      "RNA_RNP_processing",
      "Translation_Ribosome",
      "Mitochondria_OXPHOS_Energy",
      "Synapse_Vesicle",
      "Proteostasis_Ubiquitin_Folding",
      "Cytoskeleton_Motility",
      "Immune_Microglia",
      "Myelin_Glial",
      "Stress_Hormone_Response"
    ),
    pattern = c(
      "rna|ribonucleoprotein|splice|splicing|mrna|transcription|nucleolus|ribosome biogenesis",
      "translation|ribosom|peptide biosynthetic|trna|elongation|initiation factor",
      "mitochond|oxidative phosphorylation|electron transport|respiratory chain|atp synthesis|oxidoreduct|energy",
      "synap|vesicle|neurotransmitter|axon|dendrit|postsynap|presynap|exocytosis|endocytosis",
      "proteas|ubiquitin|folding|chaperone|heat shock|proteostasis|protein quality",
      "cytoskeleton|actin|tubulin|microtubule|motility|adhesion|migration|extracellular matrix",
      "immune|microglia|inflamm|complement|cytokine|phagocyt|antigen|interferon|leukocyte",
      "myelin|oligodendro|glial|axon ensheath|ensheathment",
      "stress|glucocorticoid|corticosterone|hormone|response to stimulus|oxidative stress|apoptotic"
    ),
    stringsAsFactors = FALSE
  )
}

map_terms_to_programs <- function(df, description_col = "Description") {
  if (!description_col %in% names(df)) {
    df$biological_program <- NA_character_
    return(df)
  }
  patterns <- biological_program_patterns()
  desc <- tolower(as.character(df[[description_col]]))
  hits <- lapply(seq_len(nrow(patterns)), function(i) grepl(patterns$pattern[[i]], desc, ignore.case = TRUE))
  hit_mat <- do.call(cbind, hits)
  program <- rep(NA_character_, length(desc))
  first_hit <- max.col(hit_mat, ties.method = "first")
  any_hit <- rowSums(hit_mat, na.rm = TRUE) > 0
  program[any_hit] <- patterns$biological_program[first_hit[any_hit]]
  df$biological_program <- program
  df
}

read_csv_if_exists <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

first_existing_path <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (!length(paths)) return(NA_character_)
  paths <- unique(normalizePath(paths, winslash = "/", mustWork = FALSE))
  hit <- paths[file.exists(paths)]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

latest_file <- function(root, pattern) {
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  if (!dir.exists(root)) return(NA_character_)
  files <- list.files(root, pattern = pattern, full.names = TRUE, recursive = TRUE)
  files <- files[file.exists(files)]
  if (!length(files)) return(NA_character_)
  info <- file.info(files)
  normalizePath(rownames(info)[order(info$mtime, decreasing = TRUE)[1]], winslash = "/", mustWork = FALSE)
}
