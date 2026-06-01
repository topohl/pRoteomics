# Lightweight table-schema validation for publication-facing contracts.

if (!exists("repo_path", mode = "function")) {
  paths_file <- if (file.exists(file.path("R", "paths.R"))) file.path("R", "paths.R") else file.path("..", "R", "paths.R")
  source(paths_file)
}

read_table_schema <- function(schema_name) {
  schema_path <- repo_path("inst", "schemas", paste0(schema_name, ".yml"))
  if (!file.exists(schema_path)) {
    stop("Unknown schema '", schema_name, "'. Missing file: ", schema_path, call. = FALSE)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for schema validation. Install it with install.packages('yaml') or renv::restore().", call. = FALSE)
  }
  yaml::read_yaml(schema_path)
}

validate_table_schema <- function(df, schema_name, strict = TRUE) {
  schema <- read_table_schema(schema_name)
  required <- unlist(schema$required_columns %||% character(), use.names = FALSE)
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop(
      "Schema '", schema_name, "' validation failed. Missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  col_defs <- schema$columns %||% list()
  type_errors <- character()
  value_errors <- character()
  for (col in intersect(names(col_defs), names(df))) {
    expected <- col_defs[[col]]$type %||% "any"
    ok <- switch(
      expected,
      any = TRUE,
      character = is.character(df[[col]]) || is.factor(df[[col]]),
      numeric = is.numeric(df[[col]]) || is.integer(df[[col]]) || all(is.na(df[[col]]) | suppressWarnings(!is.na(as.numeric(df[[col]])))),
      integer = is.integer(df[[col]]) || all(is.na(df[[col]]) | suppressWarnings(!is.na(as.integer(df[[col]])))),
      logical = is.logical(df[[col]]) || all(is.na(df[[col]]) | tolower(as.character(df[[col]])) %in% c("true", "false", "0", "1")),
      TRUE
    )
    if (!isTRUE(ok)) {
      type_errors <- c(type_errors, paste0(col, " expected ", expected, " got ", paste(class(df[[col]]), collapse = "/")))
    }

    allowed_values <- unlist(col_defs[[col]]$allowed_values %||% character(), use.names = FALSE)
    if (length(allowed_values)) {
      observed <- unique(as.character(df[[col]][!is.na(df[[col]])]))
      invalid <- setdiff(observed, as.character(allowed_values))
      if (length(invalid)) {
        value_errors <- c(
          value_errors,
          paste0(col, " has invalid value(s): ", paste(invalid, collapse = ", "),
                 ". Allowed: ", paste(allowed_values, collapse = ", "))
        )
      }
    }

    range <- col_defs[[col]]$range %||% NULL
    if (!is.null(range)) {
      numeric_values <- suppressWarnings(as.numeric(df[[col]]))
      non_missing <- !is.na(df[[col]])
      range_bad <- is.na(numeric_values) & non_missing
      min_value <- range$min %||% -Inf
      max_value <- range$max %||% Inf
      range_bad <- range_bad | (!is.na(numeric_values) & (numeric_values < min_value | numeric_values > max_value))
      if (any(range_bad)) {
        bad_values <- unique(as.character(df[[col]][range_bad]))
        value_errors <- c(
          value_errors,
          paste0(col, " outside allowed range [", min_value, ", ", max_value, "]: ",
                 paste(bad_values, collapse = ", "))
        )
      }
    }
  }

  if (length(type_errors)) {
    stop(
      "Schema '", schema_name, "' validation failed. Type mismatch(es): ",
      paste(type_errors, collapse = "; "),
      call. = FALSE
    )
  }

  if (length(value_errors)) {
    stop(
      "Schema '", schema_name, "' validation failed. Value constraint error(s): ",
      paste(value_errors, collapse = "; "),
      call. = FALSE
    )
  }

  if (isTRUE(strict)) {
    allowed <- names(col_defs)
    extra <- setdiff(names(df), allowed)
    if (length(allowed) && length(extra)) {
      stop(
        "Schema '", schema_name, "' validation failed. Unexpected column(s): ",
        paste(extra, collapse = ", "),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
