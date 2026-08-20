# Downstream-only helpers for the manuscript Figure 3 renderer.

figure3_gene_symbol_columns <- function() {
  c("official_gene_symbol", "representative_gene_symbol", "gene_symbol")
}

figure3_display_gene_symbol <- function(x) {
  if (!"ProteinGroupID" %in% names(x)) {
    stop("Figure 3 DA input requires ProteinGroupID for the display fallback.", call. = FALSE)
  }
  out <- rep(NA_character_, nrow(x))
  available <- intersect(figure3_gene_symbol_columns(), names(x))
  for (column in available) {
    value <- trimws(as.character(x[[column]]))
    value[is.na(value) | !nzchar(value)] <- NA_character_
    fill <- is.na(out) & !is.na(value)
    out[fill] <- value[fill]
  }
  fallback <- trimws(as.character(x$ProteinGroupID))
  fallback[is.na(fallback) | !nzchar(fallback)] <- NA_character_
  out[is.na(out)] <- fallback[is.na(out)]
  out
}

figure3_validation_record <- function(
    check, check_type, observed, expected, passed = NA,
    enforcement = "report_only", source_artifact = NA_character_,
    source_sha256 = NA_character_, note = NA_character_) {
  status <- if (is.na(passed)) {
    "observed"
  } else if (isTRUE(passed)) {
    "pass"
  } else {
    "fail"
  }
  data.frame(
    check = as.character(check),
    check_type = as.character(check_type),
    observed = paste(as.character(observed), collapse = ";"),
    expected = paste(as.character(expected), collapse = ";"),
    status = status,
    enforcement = as.character(enforcement),
    source_artifact = as.character(source_artifact),
    source_sha256 = as.character(source_sha256),
    note = as.character(note),
    stringsAsFactors = FALSE
  )
}

figure3_assert_structural_checks <- function(validation) {
  failed <- validation$check_type == "structural_contract" &
    validation$enforcement == "stop" & validation$status != "pass"
  if (any(failed)) {
    stop(
      "Figure 3 structural contract failed: ",
      paste(validation$check[failed], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(validation)
}

figure3_warn_snapshot_mismatches <- function(validation) {
  failed <- validation$check_type == "frozen_manuscript_result_regression" &
    validation$status == "fail"
  if (any(failed)) {
    warning(
      "Figure 3 frozen manuscript-result snapshot differs from the current source artifact: ",
      paste(validation$check[failed], collapse = ", "),
      ". This is not a generic analytical-validity failure; review manuscript text and provenance.",
      call. = FALSE
    )
  }
  invisible(validation)
}
