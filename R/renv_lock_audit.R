renv_lock_package_names <- function(lockfile) {
  if (!file.exists(lockfile)) {
    stop("renv lockfile does not exist: ", lockfile, call. = FALSE)
  }

  lines <- readLines(lockfile, warn = FALSE)
  packages_start <- grep('^  "Packages": \\{$', lines)
  if (length(packages_start) != 1L) {
    stop("renv lockfile must contain exactly one Packages object", call. = FALSE)
  }

  package_lines <- lines[seq.int(packages_start + 1L, length(lines))]
  package_keys <- regexec('^    "([^"]+)": \\{$', package_lines)
  matches <- regmatches(package_lines, package_keys)
  sort(unique(vapply(matches[lengths(matches) == 2L], `[[`, character(1), 2L)))
}

renv_lock_scientific_sentinels <- function() {
  c(
    "clusterProfiler", "dplyr", "fgsea", "ggplot2", "limma", "lme4",
    "readr", "WGCNA"
  )
}

audit_renv_lock <- function(
    lockfile,
    scientific_sentinels = renv_lock_scientific_sentinels()) {
  packages <- renv_lock_package_names(lockfile)
  missing <- setdiff(scientific_sentinels, packages)

  list(
    lockfile = normalizePath(lockfile, winslash = "/", mustWork = TRUE),
    package_count = length(packages),
    packages = packages,
    scientific_sentinels = scientific_sentinels,
    missing_scientific_sentinels = missing,
    plausibly_full_scientific_lock = length(missing) == 0L
  )
}

# ===========================================================================
# Lockfile construction from the installed analysis library.
#
# renv is NOT installed in this environment and there is no renv project
# infrastructure (no .Rprofile, no renv/activate.R, no renv/settings.json, no
# DESCRIPTION). The previous three-record lockfile could therefore not have
# been produced by renv::snapshot(); it was hand-seeded for the contract-test
# bootstrap. Rather than install renv and create a project library -- which
# would materially alter the live analysis environment -- the lockfile is
# written directly from the installed library, which is exactly what "record
# the environment that exists" requires.
#
# Nothing here installs, updates or loads a package. Versions and source
# metadata are read from installed DESCRIPTION files only, and no field is
# invented: where provenance is absent it is omitted rather than guessed.
# ===========================================================================

# Active project surface. Archived trees are excluded: a package referenced
# only by deprecated or scratch code is not a dependency of the current
# publication analysis.
renv_lock_excluded_trees <- function() {
  c("99_deprecated/", "90_testing/", "legacy/", "_scratchpad/")
}

renv_lock_active_source_files <- function(root = ".") {
  registry <- file.path(root, "pipeline.yml")
  scripts <- character()
  if (file.exists(registry) && requireNamespace("yaml", quietly = TRUE)) {
    reg <- yaml::read_yaml(registry)
    scripts <- unique(unlist(lapply(reg$stages, function(st)
      vapply(st$scripts, function(s) s$script, character(1)))))
    scripts <- file.path(root, scripts)
  }
  helpers <- list.files(file.path(root, "R"), pattern = "[.][Rr]$", full.names = TRUE)
  tools <- list.files(file.path(root, "tools"), pattern = "[.][Rr]$", full.names = TRUE)
  tests <- c(
    list.files(file.path(root, "tests"), pattern = "[.][Rr]$", full.names = TRUE),
    list.files(file.path(root, "tests", "testthat"), pattern = "[.][Rr]$", full.names = TRUE)
  )
  files <- unique(c(scripts, helpers, tools, tests))
  files <- files[file.exists(files)]
  excluded <- renv_lock_excluded_trees()
  keep <- !vapply(files, function(p) {
    any(vapply(excluded, function(x) grepl(x, p, fixed = TRUE), logical(1)))
  }, logical(1))
  sort(files[keep], method = "radix")
}

# Tokens that the scanner picks up but which are not packages: string literals
# containing "::" (used as composite keys such as "Neuropil::CA1"), the loop
# parameter `pkg`, and `renv` itself, which appears only inside error-message
# text advising the reader to run renv::restore().
renv_lock_token_false_positives <- function() {
  c("A1", "A111", "BP", "GSEA_GO", "Neuropil", "ROI", "Soma", "pkg", "renv")
}

renv_lock_dependency_patterns <- function() {
  c(
    library   = "library\\s*\\(\\s*[\"']?([A-Za-z][A-Za-z0-9._]*)[\"']?\\s*[,)]",
    require   = "(?<!requireNamespace)\\brequire\\s*\\(\\s*[\"']?([A-Za-z][A-Za-z0-9._]*)[\"']?\\s*[,)]",
    requireNS = "requireNamespace\\s*\\(\\s*[\"']([A-Za-z][A-Za-z0-9._]*)[\"']",
    loadNS    = "loadNamespace\\s*\\(\\s*[\"']([A-Za-z][A-Za-z0-9._]*)[\"']",
    getNS     = "getNamespace\\s*\\(\\s*[\"']([A-Za-z][A-Za-z0-9._]*)[\"']",
    colon2    = "([A-Za-z][A-Za-z0-9._]*)::",
    colon3    = "([A-Za-z][A-Za-z0-9._]*):::"
  )
}

renv_lock_discover_direct_dependencies <- function(
    files = renv_lock_active_source_files(),
    installed = utils::installed.packages(),
    false_positives = renv_lock_token_false_positives()) {
  patterns <- renv_lock_dependency_patterns()
  found <- character()
  for (f in files) {
    txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
    if (!length(txt)) next
    code <- sub("#.*$", "", txt)
    for (p in patterns) {
      m <- regmatches(code, gregexpr(p, code, perl = TRUE))
      for (grp in m) {
        if (!length(grp)) next
        pkg <- gsub("[^A-Za-z0-9._]", "", sub(p, "\\1", grp, perl = TRUE))
        found <- c(found, pkg[nzchar(pkg)])
      }
    }
  }
  found <- setdiff(unique(found), false_positives)
  base_priority <- rownames(installed)[
    !is.na(installed[, "Priority"]) & installed[, "Priority"] == "base"
  ]
  sort(intersect(setdiff(found, base_priority), rownames(installed)), method = "radix")
}

renv_lock_dependency_closure <- function(direct,
                                         installed = utils::installed.packages()) {
  base_priority <- rownames(installed)[
    !is.na(installed[, "Priority"]) & installed[, "Priority"] == "base"
  ]
  deps <- tools::package_dependencies(
    direct, db = installed,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE
  )
  all <- sort(unique(c(direct, unlist(deps, use.names = FALSE))), method = "radix")
  all <- setdiff(all, base_priority)
  list(
    packages = all[all %in% rownames(installed)],
    unresolved = sort(setdiff(unlist(deps, use.names = FALSE),
                              c(rownames(installed), base_priority)), method = "radix")
  )
}

# Repositories that packages in this library actually came from, plus the
# Bioconductor 3.22 repository set. Nothing speculative is added.
renv_lock_repositories <- function() {
  repos <- list(
    c(Name = "CRAN", URL = "https://cloud.r-project.org"),
    c(Name = "BioCsoft", URL = "https://bioconductor.org/packages/3.22/bioc"),
    c(Name = "BioCann", URL = "https://bioconductor.org/packages/3.22/data/annotation"),
    c(Name = "BioCexp", URL = "https://bioconductor.org/packages/3.22/data/experiment"),
    c(Name = "BioCworkflows", URL = "https://bioconductor.org/packages/3.22/workflows"),
    c(Name = "BioCbooks", URL = "https://bioconductor.org/packages/3.22/books"),
    # 40 Bioconductor packages in this library record this mirror as their
    # install source, so it must be declared for the lockfile to be resolvable.
    c(Name = "bioc-release-r-universe", URL = "https://bioc-release.r-universe.dev")
  )
  lapply(repos, function(r) list(Name = unname(r[["Name"]]), URL = unname(r[["URL"]])))
}

# One faithful record per package. Bioconductor membership is decided by the
# presence of biocViews in the installed DESCRIPTION, which is the field
# Bioconductor itself requires; the Repository string is copied verbatim when
# present and omitted when the installed package does not carry one (true for
# the AnnotationDbi-built annotation packages).
renv_lock_package_record <- function(pkg, installed = utils::installed.packages()) {
  dsc <- utils::packageDescription(pkg)
  field <- function(name) {
    v <- dsc[[name]]
    if (is.null(v) || !length(v) || is.na(v[[1]])) NA_character_ else as.character(v)[[1]]
  }
  repository <- field("Repository")
  is_bioc <- !is.na(field("biocViews"))

  base_priority <- rownames(installed)[
    !is.na(installed[, "Priority"]) & installed[, "Priority"] == "base"
  ]
  reqs <- unique(unlist(tools::package_dependencies(
    pkg, db = installed, which = c("Depends", "Imports", "LinkingTo"),
    recursive = FALSE
  ), use.names = FALSE))
  reqs <- sort(setdiff(reqs, base_priority), method = "radix")

  record <- list(
    Package = pkg,
    Version = as.character(utils::packageVersion(pkg)),
    Source = if (is_bioc) "Bioconductor" else "Repository"
  )
  if (!is.na(repository)) {
    record$Repository <- repository
  }
  if (length(reqs)) record$Requirements <- reqs
  record
}

# Deterministic renv-schema lockfile. Hash is deliberately absent: renv
# computes it with an internal algorithm that cannot be reproduced without renv
# installed, and inventing a value would misrepresent provenance.
build_renv_lockfile <- function(root = ".",
                                installed = utils::installed.packages(),
                                bioc_version = NULL) {
  direct <- renv_lock_discover_direct_dependencies(
    renv_lock_active_source_files(root), installed = installed
  )
  closure <- renv_lock_dependency_closure(direct, installed = installed)
  if (length(closure$unresolved)) {
    stop("Dependency closure references packages absent from the library: ",
         paste(closure$unresolved, collapse = ", "), call. = FALSE)
  }
  if (is.null(bioc_version)) {
    bioc_version <- tryCatch(as.character(BiocManager::version()),
                             error = function(e) NA_character_)
  }
  records <- lapply(closure$packages, renv_lock_package_record, installed = installed)
  names(records) <- closure$packages

  lock <- list(
    R = list(
      Version = paste(R.version$major, R.version$minor, sep = "."),
      Repositories = renv_lock_repositories()
    )
  )
  if (!is.na(bioc_version)) {
    lock$Bioconductor <- list(Version = bioc_version)
  }
  lock$Packages <- records
  attr(lock, "direct_dependencies") <- direct
  lock
}

write_renv_lockfile <- function(lock, path = "renv.lock") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to write the lockfile.", call. = FALSE)
  }
  json <- jsonlite::toJSON(lock, auto_unbox = TRUE, pretty = 2, null = "null")
  # renv terminates the lockfile with a single trailing newline.
  cat(json, "\n", sep = "", file = path)
  invisible(path)
}

# Completeness audit: every discovered direct dependency must be recorded, and
# every Requirements reference must resolve inside the lockfile.
audit_renv_lock_completeness <- function(lockfile = "renv.lock", root = ".",
                                         installed = utils::installed.packages()) {
  recorded <- renv_lock_package_names(lockfile)
  direct <- renv_lock_discover_direct_dependencies(
    renv_lock_active_source_files(root), installed = installed
  )
  parsed <- if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::fromJSON(lockfile, simplifyVector = FALSE)
  } else NULL

  unresolved <- character()
  duplicates <- character()
  if (!is.null(parsed) && !is.null(parsed$Packages)) {
    nms <- names(parsed$Packages)
    duplicates <- unique(nms[duplicated(nms)])
    for (nm in nms) {
      reqs <- parsed$Packages[[nm]]$Requirements
      if (is.null(reqs)) next
      missing <- setdiff(unlist(reqs, use.names = FALSE), nms)
      if (length(missing)) {
        unresolved <- c(unresolved, paste0(nm, " -> ", paste(missing, collapse = ",")))
      }
    }
  }

  list(
    lockfile = normalizePath(lockfile, winslash = "/", mustWork = TRUE),
    recorded_count = length(recorded),
    direct_dependencies = direct,
    missing_direct = setdiff(direct, recorded),
    unresolved_requirements = unresolved,
    duplicate_records = duplicates,
    r_version = if (!is.null(parsed)) parsed$R$Version else NA_character_,
    bioconductor_version = if (!is.null(parsed)) parsed$Bioconductor$Version else NULL,
    complete = length(setdiff(direct, recorded)) == 0L &&
      length(unresolved) == 0L && length(duplicates) == 0L
  )
}
