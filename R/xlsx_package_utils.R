# OOXML package validation and repair for workbooks written with openxlsx.
#
# WHY THIS EXISTS
#   openxlsx (4.2.8.1) emits, for EVERY worksheet, a `drawing` relationship in
#   xl/worksheets/_rels/sheetN.xml.rels and a matching Override in
#   [Content_Types].xml pointing at /xl/drawings/drawingN.xml, but it only
#   serialises that drawing part when the sheet actually carries drawing
#   content.  A workbook with no images, shapes or comments therefore ships
#   relationships and content-type overrides that reference parts which are not
#   in the archive.
#
#   This reproduces with the smallest possible workbook
#   (createWorkbook + addWorksheet + writeData), so it is a library defect and
#   not a consequence of styles, freeze panes, tables or comments.
#
#   Lenient readers (readxl, LibreOffice) ignore the dangling references.
#   Strict OOXML consumers resolve every relationship and fail with errors of
#   the PartDoesNotExist family.
#
# WHAT THIS DOES
#   `xlsx_repair_package()` rewrites the archive, dropping relationships whose
#   target part is absent and Content_Types overrides for parts that were never
#   written.  It is deliberately generic: it removes any dangling relationship,
#   not just drawing ones, and it never adds or edits cell content.

# `%||%` comes from the canonical R/null_coalescing.R, loaded via R/paths.R.
# It is deliberately not redefined here: the repository permits exactly one
# definition so that source order cannot change coalescing semantics.

# Collapse "a/b/../c" and "./c" into a package-root-relative part name.
.xlsx_normalize_part <- function(base, target) {
  target <- gsub("\\\\", "/", as.character(target))
  if (startsWith(target, "/")) {
    parts <- sub("^/", "", target)
  } else {
    parts <- if (nzchar(base)) paste0(base, "/", target) else target
  }
  segments <- strsplit(parts, "/", fixed = TRUE)[[1]]
  out <- character(0)
  for (segment in segments) {
    if (!nzchar(segment) || identical(segment, ".")) next
    if (identical(segment, "..")) {
      if (length(out)) out <- out[-length(out)]
      next
    }
    out <- c(out, segment)
  }
  paste(out, collapse = "/")
}

# Every relationship declared anywhere in the package, with a resolved target
# and whether that target is present in the archive.
xlsx_package_relationships <- function(path) {
  entries <- utils::unzip(path, list = TRUE)$Name
  rels_files <- grep("_rels/[^/]*\\.rels$", entries, value = TRUE)
  rows <- lapply(rels_files, function(rels) {
    con <- unz(path, rels)
    on.exit(close(con), add = TRUE)
    txt <- paste(readLines(con, warn = FALSE), collapse = "")
    nodes <- regmatches(txt, gregexpr("<Relationship\\b[^>]*/>", txt))[[1]]
    if (!length(nodes)) return(NULL)
    base <- dirname(dirname(rels))
    if (identical(base, ".")) base <- ""
    attr_of <- function(node, name) {
      hit <- regmatches(node, regexpr(paste0(name, '="[^"]*"'), node))
      if (!length(hit)) return(NA_character_)
      sub('"$', "", sub(paste0(name, '="'), "", hit))
    }
    target <- vapply(nodes, attr_of, character(1), name = "Target", USE.NAMES = FALSE)
    mode <- vapply(nodes, attr_of, character(1), name = "TargetMode", USE.NAMES = FALSE)
    external <- !is.na(mode) & mode == "External"
    resolved <- ifelse(
      external, NA_character_,
      vapply(target, function(t) .xlsx_normalize_part(base, t), character(1),
             USE.NAMES = FALSE)
    )
    data.frame(
      rels_part = rels,
      id = vapply(nodes, attr_of, character(1), name = "Id", USE.NAMES = FALSE),
      type = vapply(nodes, attr_of, character(1), name = "Type", USE.NAMES = FALSE),
      target = target,
      external = external,
      resolved_part = resolved,
      exists = external | resolved %in% entries,
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    return(data.frame(
      rels_part = character(), id = character(), type = character(),
      target = character(), external = logical(), resolved_part = character(),
      exists = logical(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

# Relationships whose target part is missing from the archive.
xlsx_package_dangling_relationships <- function(path) {
  rel <- xlsx_package_relationships(path)
  rel[!rel$exists, , drop = FALSE]
}

# Content_Types overrides naming parts that are not in the archive.
xlsx_package_orphan_overrides <- function(path) {
  entries <- utils::unzip(path, list = TRUE)$Name
  con <- unz(path, "[Content_Types].xml")
  on.exit(close(con), add = TRUE)
  txt <- paste(readLines(con, warn = FALSE), collapse = "")
  nodes <- regmatches(txt, gregexpr("<Override\\b[^>]*/>", txt))[[1]]
  if (!length(nodes)) return(character(0))
  part <- sub('"$', "", sub('PartName="', "", regmatches(
    nodes, regexpr('PartName="[^"]*"', nodes)
  )))
  part <- sub("^/", "", part)
  part[!(part %in% entries)]
}

# TRUE when every internal relationship resolves and no override is orphaned.
xlsx_package_is_valid <- function(path) {
  nrow(xlsx_package_dangling_relationships(path)) == 0L &&
    length(xlsx_package_orphan_overrides(path)) == 0L
}

xlsx_assert_package_valid <- function(path, label = basename(path)) {
  dangling <- xlsx_package_dangling_relationships(path)
  orphans <- xlsx_package_orphan_overrides(path)
  if (nrow(dangling) || length(orphans)) {
    stop(
      "Malformed OOXML package in ", label, ": ",
      nrow(dangling), " dangling relationship(s), ",
      length(orphans), " orphaned content-type override(s).",
      if (nrow(dangling)) paste0(
        " First: ", dangling$rels_part[[1]], " -> ", dangling$target[[1]], "."
      ) else "",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Rewrite the archive without dangling relationships or orphaned overrides.
#
# Cell content, styles, tables and sheet order are untouched: only package
# plumbing that points at absent parts is removed.  Sheet XML elements whose
# r:id no longer resolves (<drawing>, <legacyDrawing>, <picture>) are dropped
# too, so the sheet never references a relationship that is gone.
xlsx_repair_package <- function(path) {
  if (!requireNamespace("zip", quietly = TRUE)) {
    stop("Repairing an xlsx package requires the 'zip' package.", call. = FALSE)
  }
  dangling <- xlsx_package_dangling_relationships(path)
  orphans <- xlsx_package_orphan_overrides(path)
  if (!nrow(dangling) && !length(orphans)) {
    return(invisible(list(path = path, removed_relationships = 0L,
                          removed_overrides = 0L)))
  }

  workdir <- file.path(tempdir(), paste0("xlsxfix_", basename(tempfile(""))))
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(path, exdir = workdir)

  # 1. strip dangling <Relationship> nodes, remembering the ids per rels part
  removed_ids <- list()
  for (rels in unique(dangling$rels_part)) {
    file <- file.path(workdir, rels)
    txt <- paste(readLines(file, warn = FALSE), collapse = "")
    ids <- dangling$id[dangling$rels_part == rels]
    for (id in ids) {
      pattern <- paste0('<Relationship[^>]*Id="', id, '"[^>]*/>')
      txt <- sub(pattern, "", txt)
    }
    removed_ids[[rels]] <- ids
    writeLines(txt, file, useBytes = TRUE)
  }

  # 2. drop sheet elements that referenced those now-removed relationship ids
  for (rels in names(removed_ids)) {
    sheet <- file.path(workdir, sub("_rels/", "", sub("\\.rels$", "", rels)))
    if (!file.exists(sheet)) next
    txt <- paste(readLines(sheet, warn = FALSE), collapse = "")
    changed <- FALSE
    for (id in removed_ids[[rels]]) {
      for (tag in c("drawing", "legacyDrawing", "legacyDrawingHF", "picture")) {
        pattern <- paste0("<", tag, "[^>]*r:id=\"", id, "\"[^>]*/>")
        if (grepl(pattern, txt)) {
          txt <- gsub(pattern, "", txt)
          changed <- TRUE
        }
      }
    }
    if (changed) writeLines(txt, sheet, useBytes = TRUE)
  }

  # 3. strip <Override> nodes for parts that were never written
  if (length(orphans)) {
    ct <- file.path(workdir, "[Content_Types].xml")
    txt <- paste(readLines(ct, warn = FALSE), collapse = "")
    for (part in orphans) {
      pattern <- paste0('<Override PartName="/', part, '"[^>]*/>')
      txt <- sub(pattern, "", txt)
    }
    writeLines(txt, ct, useBytes = TRUE)
  }

  # 4. repack, keeping [Content_Types].xml first as readers conventionally expect
  files <- list.files(workdir, recursive = TRUE, all.files = TRUE,
                      no.. = TRUE)
  files <- files[!dir.exists(file.path(workdir, files))]
  files <- c("[Content_Types].xml", setdiff(files, "[Content_Types].xml"))
  # "mirror" keeps each entry's path relative to `root`; "cherry-pick" would
  # flatten the package to basenames and destroy the OOXML directory layout.
  tmp_zip <- paste0(path, ".repair")
  if (file.exists(tmp_zip)) unlink(tmp_zip)
  zip::zip(tmp_zip, files = files, root = workdir, mode = "mirror")
  stored <- zip::zip_list(tmp_zip)$filename
  missing <- setdiff(files, stored)
  if (length(missing)) {
    unlink(tmp_zip)
    stop("Refusing to write a repaired package that lost ", length(missing),
         " part(s), first: ", missing[[1]], ".", call. = FALSE)
  }
  file.copy(tmp_zip, path, overwrite = TRUE)
  unlink(tmp_zip)

  invisible(list(path = path,
                 removed_relationships = nrow(dangling),
                 removed_overrides = length(orphans)))
}

# Save an openxlsx workbook and guarantee a well-formed package on disk.
xlsx_save_valid_workbook <- function(wb, path, overwrite = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, path, overwrite = overwrite)
  repaired <- xlsx_repair_package(path)
  xlsx_assert_package_valid(path)
  invisible(repaired)
}
