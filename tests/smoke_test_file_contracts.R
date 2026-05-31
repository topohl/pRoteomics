#!/usr/bin/env Rscript

paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
source(paths_file)

contract_file <- repo_path("docs", "file_contracts.tsv")
required_cols <- c("object_id", "path", "created_by", "consumed_by", "required_columns", "description", "version")

fail <- character()
if (!file.exists(contract_file)) {
  fail <- c(fail, paste("Missing", contract_file))
} else {
  contracts <- utils::read.delim(contract_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  missing_cols <- setdiff(required_cols, names(contracts))
  if (length(missing_cols)) fail <- c(fail, paste("docs/file_contracts.tsv missing columns:", paste(missing_cols, collapse = ", ")))
  if ("object_id" %in% names(contracts) && any(duplicated(contracts$object_id))) {
    fail <- c(fail, paste("Duplicate contract object keys:", paste(unique(contracts$object_id[duplicated(contracts$object_id)]), collapse = ", ")))
  }
  if ("path" %in% names(contracts) && any(!nzchar(contracts$path))) {
    fail <- c(fail, "At least one file contract has an empty path.")
  }
}

if (length(fail)) {
  message("FAIL smoke_test_file_contracts")
  message(paste(fail, collapse = "\n"))
  quit(status = 1, save = "no")
}

message("PASS smoke_test_file_contracts")
