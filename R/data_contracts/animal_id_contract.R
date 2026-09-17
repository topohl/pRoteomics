# =====================================================================
# Canonical AnimalID contract.
#
# WHY THIS EXISTS
#
# Three ad-hoc animal-identifier normalisers grew up independently in this
# repository and disagreed with each other:
#
#   analysis/08_integration/test_network_behaviour_coupling.R
#       str_extract(x, "A[0-9]{3,4}|[0-9]{3,4}") then "A" + str_pad(x, 4)
#   analysis/05_wgcna/score_module_activity.R
#       str_extract(x, "\\d+$") then sprintf("%04d")
#   analysis/08_integration/test_behaviour_proteomics_associations.R
#       gsub("[^0-9]", "") then sub("^0+", "")
#
# The first one is not merely lossy, it MANUFACTURES WRONG ANIMALS, in three
# independent ways:
#
#   1. TRUNCATION. "[0-9]{3,4}" is capped at four digits and is anchored
#      nowhere, so it keeps only the FIRST four digits of a longer run.
#      "13856" and "13857" both become "1385"; "00690".."00696" all become
#      "0069". Across the 117 physiology rows this collapses 117 distinct
#      animals to 101 through 7 collision groups covering 23 animals, and the
#      downstream distinct(AnimalID, .keep_all = TRUE) then SILENTLY DROPS the
#      16 losers while misattributing the survivor's phenotype to the key.
#   2. BRANCH ASYMMETRY. Values that already start with "A" are returned with
#      their original digit width, while bare numerals are zero-padded to four.
#      So "A111" stays "A111" but "OR111" becomes "A0111". The function is not
#      idempotent, and the two forms of the SAME animal never compare equal.
#      This is what silently reduced the network-behaviour join to one animal.
#   3. MINIMUM WIDTH. "[0-9]{3,4}" needs at least three digits, so the
#      canonical AnimalID "3" resolves to NA.
#
# THE CANONICAL FORM
#
# The authoritative value is the one stored in the canonical merged metadata:
# a BARE, UNPADDED decimal string. For Exp9 that is
#   "3", "111", "127", "129", "135", "139", "755", "764", "765".
# Note "3", not "0003" and not "A0003".
#
# RESOLUTION ORDER (deliberate, and each step is allowed to fail closed)
#   1. exact canonical match
#   2. explicit alias lookup from config/animal_id_aliases.csv
#   3. formatting-only normalisation: case, whitespace, a documented leading
#      cohort prefix, and leading zeros
#   4. FAIL CLOSED - unresolved or ambiguous ids are NA in audit mode and an
#      error in strict mode
#
# Formatting normalisation may never delete terminal digits to force a match.
# aid_format_normalize() asserts that the value it returns is a digit-suffix
# of the input's own digit run, so a truncating change cannot pass review.
# =====================================================================

# The nine Exp9 animals, frozen as a contract constant. A test asserts that
# the canonical metadata still agrees with this list; it is NOT a substitute
# for reading the metadata.
aid_expected_exp9_animals <- function() {
  c("3", "111", "127", "129", "135", "139", "755", "764", "765")
}

aid_expected_exp9_groups <- function() {
  c("3" = "RES", "111" = "SUS", "127" = "CON", "129" = "CON", "135" = "RES",
    "139" = "RES", "755" = "SUS", "764" = "SUS", "765" = "CON")
}

aid_alias_path <- function() repo_path("config", "animal_id_aliases.csv")

# ------------------------------------------------------------------ aliases

aid_alias_table <- function(path = aid_alias_path()) {
  if (!file.exists(path)) {
    return(data.frame(source_system = character(), raw_id = character(),
                      canonical_AnimalID = character(), evidence = character(),
                      status = character(), stringsAsFactors = FALSE))
  }
  a <- utils::read.csv(path, stringsAsFactors = FALSE, colClasses = "character")
  need <- c("source_system", "raw_id", "canonical_AnimalID", "evidence", "status")
  miss <- setdiff(need, names(a))
  if (length(miss)) {
    stop("animal_id_aliases.csv is missing column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  a$source_system <- trimws(a$source_system)
  a$raw_id <- trimws(a$raw_id)
  a$canonical_AnimalID <- trimws(a$canonical_AnimalID)
  a <- a[a$status %in% "active", , drop = FALSE]

  # A raw id may not resolve to two different animals within one source system.
  key <- paste(a$source_system, a$raw_id, sep = "\r")
  amb <- unique(key[duplicated(key)])
  amb <- amb[vapply(amb, function(k)
    length(unique(a$canonical_AnimalID[key == k])) > 1L, logical(1))]
  if (length(amb)) {
    stop("animal_id_aliases.csv maps one raw id to more than one AnimalID: ",
         paste(sub("\r", "/", amb), collapse = ", "), call. = FALSE)
  }
  a[!duplicated(key), , drop = FALSE]
}

# ------------------------------------------------- formatting normalisation

# Case, whitespace, a documented leading alphabetic cohort prefix (OQ/OR/A/...)
# and leading zeros. Nothing else. In particular no terminal digit is ever
# dropped: the returned digits must be a suffix of the input's digit run.
aid_format_normalize <- function(x) {
  raw <- as.character(x)
  u <- toupper(trimws(raw))
  u[!nzchar(u)] <- NA_character_
  u <- sub("[.]0+$", "", u)              # 111.0 written by a spreadsheet

  # A leading alphabetic cohort prefix is dropped ONLY when what remains is
  # entirely digits. "OQ754" -> "754"; "A1B2" is left alone and fails closed.
  stripped <- sub("^[A-Z]+", "", u)
  ok <- !is.na(stripped) & grepl("^[0-9]+$", stripped)
  out <- rep(NA_character_, length(u))
  out[ok] <- sub("^0+(?=[0-9])", "", stripped[ok], perl = TRUE)

  # Guard: the result must be a digit-suffix of the input's own digit run.
  # A truncating implementation (the historical defect) fails this.
  chk <- !is.na(out)
  if (any(chk)) {
    digits_in <- gsub("[^0-9]", "", u[chk])
    bad <- substring(digits_in, nchar(digits_in) - nchar(out[chk]) + 1L) != out[chk]
    if (any(bad)) {
      stop("aid_format_normalize dropped non-leading digits for: ",
           paste(utils::head(raw[chk][bad], 5), collapse = ", "), call. = FALSE)
    }
  }
  out
}

# --------------------------------------------------------------- resolution

# Resolve raw identifiers to canonical AnimalID.
#   strict = TRUE  -> stop() on any unresolved id
#   strict = FALSE -> NA for unresolved ids (audit mode)
# `canonical` restricts the accepted output vocabulary; when supplied, a value
# that normalises to something outside it is treated as unresolved rather than
# invented.
aid_resolve <- function(x, source_system, strict = TRUE, canonical = NULL,
                        aliases = aid_alias_table()) {
  raw <- as.character(x)
  out <- rep(NA_character_, length(raw))
  route <- rep(NA_character_, length(raw))

  # 1. exact canonical match
  if (!is.null(canonical)) {
    hit <- trimws(raw) %in% canonical
    out[hit] <- trimws(raw)[hit]
    route[hit] <- "exact_canonical"
  }

  # 2. explicit alias lookup
  todo <- is.na(out)
  if (any(todo) && nrow(aliases)) {
    al <- aliases[aliases$source_system == source_system, , drop = FALSE]
    if (nrow(al)) {
      m <- match(trimws(raw[todo]), al$raw_id)
      v <- al$canonical_AnimalID[m]
      idx <- which(todo)[!is.na(v)]
      out[idx] <- v[!is.na(v)]
      route[idx] <- "alias"
    }
  }

  # 3. formatting-only normalisation
  todo <- is.na(out)
  if (any(todo)) {
    f <- aid_format_normalize(raw[todo])
    if (!is.null(canonical)) f[!is.na(f) & !(f %in% canonical)] <- NA_character_
    idx <- which(todo)[!is.na(f)]
    out[idx] <- f[!is.na(f)]
    route[idx] <- "format_normalized"
  }

  # 4. fail closed
  if (isTRUE(strict) && any(is.na(out) & !is.na(raw) & nzchar(trimws(raw)))) {
    bad <- unique(raw[is.na(out) & !is.na(raw) & nzchar(trimws(raw))])
    stop("unresolved AnimalID in source system '", source_system, "': ",
         paste(utils::head(bad, 10), collapse = ", "),
         if (length(bad) > 10) paste0(" (+", length(bad) - 10, " more)") else "",
         ". Add an explicit row to config/animal_id_aliases.csv; ",
         "do not widen the formatting rule.", call. = FALSE)
  }
  attr(out, "route") <- route
  out
}

# Hard-fail when normalisation itself merges two distinct raw identifiers.
# This is the check the historical defect would not have survived.
aid_assert_no_collision <- function(raw, resolved, context = "") {
  raw <- as.character(raw); resolved <- as.character(resolved)
  keep <- !is.na(resolved)
  if (!any(keep)) return(invisible(TRUE))
  sp <- split(unique(data.frame(raw = raw[keep], res = resolved[keep],
                                stringsAsFactors = FALSE)), resolved[keep][!duplicated(
    paste(raw[keep], resolved[keep]))])
  tab <- tapply(raw[keep], resolved[keep], function(z) length(unique(z)))
  bad <- names(tab)[tab > 1L]
  if (length(bad)) {
    detail <- vapply(bad, function(k)
      paste0(k, " <- {", paste(sort(unique(raw[keep][resolved[keep] == k]),
                                    decreasing = FALSE), collapse = ", "), "}"),
      character(1))
    stop("AnimalID normalisation collision", if (nzchar(context))
      paste0(" in ", context) else "", ": ",
      paste(detail, collapse = "; "),
      ". Two distinct source identifiers resolved to one AnimalID.",
      call. = FALSE)
  }
  invisible(TRUE)
}

# Audit-mode resolution report: one row per distinct raw id.
aid_resolution_report <- function(x, source_system, canonical = NULL,
                                  aliases = aid_alias_table()) {
  raw <- unique(as.character(x))
  res <- aid_resolve(raw, source_system, strict = FALSE, canonical = canonical,
                     aliases = aliases)
  route <- attr(res, "route")
  tab <- tapply(raw, res, function(z) length(unique(z)))
  n_share <- ifelse(is.na(res), NA_integer_, as.integer(tab[res]))
  data.frame(
    source_system = source_system,
    raw_id = raw,
    canonical_AnimalID = as.character(res),
    resolution_route = route,
    resolved = !is.na(res),
    n_raw_ids_sharing_this_AnimalID = n_share,
    collision = !is.na(n_share) & n_share > 1L,
    stringsAsFactors = FALSE)
}

# The historical broken normaliser, kept ONLY so the audit and the regression
# tests can demonstrate what it did. Never call it to produce analysis values.
aid_legacy_broken_normalizer <- function(x) {
  x <- trimws(as.character(x))
  m <- regmatches(x, regexpr("A[0-9]{3,4}|[0-9]{3,4}", x))
  out <- rep(NA_character_, length(x))
  has <- regexpr("A[0-9]{3,4}|[0-9]{3,4}", x) > 0
  out[has] <- m
  # str_pad(x, 4, pad = "0") on a character vector, reproduced without stringr.
  # formatC(flag = "0") pads character input with SPACES, not zeros.
  pad0 <- function(z) ifelse(nchar(z) >= 4L, z,
                             paste0(strrep("0", 4L - nchar(z)), z))
  ifelse(is.na(out), NA_character_,
         ifelse(substr(out, 1, 1) == "A", out, paste0("A", pad0(out))))
}
