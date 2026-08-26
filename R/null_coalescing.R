# Canonical project null/empty/scalar-missing coalescing semantics.
#
# The repository does not currently declare R >= 4.4 as its minimum supported
# version, and base `%||%` only falls back for NULL. Historical project callers
# also rely on empty values and scalar NA values falling back. Vectors are
# always retained, including partially missing and all-missing vectors.

proteomics_null_coalesce <- function(x, fallback) {
  if (
    is.null(x) ||
      length(x) == 0L ||
      (length(x) == 1L && is.na(x))
  ) {
    fallback
  } else {
    x
  }
}

`%||%` <- proteomics_null_coalesce
