#!/usr/bin/env Rscript

# Part B, sections 13-16 and 19-22: what the repository actually contains, how
# it is wired, and where its structure contradicts its own rules.
#
# AUDIT ONLY. Nothing is moved, renamed or deleted. The file list comes from
# git rather than from a directory walk so that untracked scratch output cannot
# inflate the inventory, and the active/inactive split is taken from the
# pipeline registry's own exclusion rules rather than restated here.

setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/proteomics")
source("audits/publication_hardening/00_checkpoint.R")
source("R/paths.R")
source(repo_path("R", "pipeline_registry.R"))
suppressMessages(library(yaml))

dir.create(PH_TAB, recursive = TRUE, showWarnings = FALSE)
REG <- read_pipeline_registry()
ENTRIES <- pipeline_registry_entries(REG)
tracked <- system2("git", c("ls-files"), stdout = TRUE)
scripts <- grep("[.][Rr]$", tracked, value = TRUE)

# ============================================== B13 architecture inventory
# 99_audits, 99_deprecated and 90_testing also carry a two-digit prefix, so
# they must be classified before the generic numbered-stage rule or they are
# silently absorbed into it.
LAYER <- function(p) {
  top <- sub("/.*", "", p)
  if (top == "99_audits") "audit layer"
  else if (top == "99_deprecated") "deprecated"
  else if (top == "90_testing") "testing scaffold"
  else if (grepl("^[0-9]{2}_", top)) "numbered analysis stage"
  else if (top == "R") "shared helper library"
  else if (top == "figures") "figure / manuscript layer"
  else if (top == "tests") "test suite"
  else if (top == "tools") "developer tooling"
  else if (top == "config") "configuration"
  else if (top == "inst") "packaged assets"
  else if (!grepl("/", p)) "repository root"
  else "other"
}
inv <- data.frame(
  path = scripts, top_level = sub("/.*", "", scripts),
  layer = vapply(scripts, LAYER, character(1)),
  registered_in_pipeline = scripts %in% ENTRIES$script,
  stringsAsFactors = FALSE)
inv$registered_stage <- ENTRIES$stage[match(inv$path, ENTRIES$script)]
inv$lines <- vapply(inv$path, function(p)
  tryCatch(length(readLines(p, warn = FALSE)), error = function(e) NA_integer_),
  integer(1))
# the registry's own view of which scripts are supposed to be pipeline steps
act <- tryCatch(active_analysis_scripts(), error = function(e) character(0))
inv$registry_considers_active <- inv$path %in% act
rownames(inv) <- NULL
utils::write.csv(inv, file.path(PH_TAB, "repository_architecture_inventory.csv"),
                 row.names = FALSE)

# ================================================= B14 active dependency graph
#
# An edge is a source()/sourcing of one repository file by another. Paths are
# written through repo_path(), so the second argument list is reassembled
# rather than matched as a literal path.
# the lookbehind keeps helper names that merely end in "source" - fp_source(),
# path_source() - from being read as a source() call
SRC <- paste0("(?<![A-Za-z0-9_.])(source|sys[.]source)\\s*\\(\\s*",
              "(repo_path\\s*\\(([^)]*)\\)|[\"']([^\"']+)[\"'])")
edges <- do.call(rbind, lapply(scripts, function(p) {
  ln <- readLines(p, warn = FALSE)
  m <- regmatches(ln, gregexpr(SRC, ln, perl = TRUE))
  m <- unlist(m)
  if (!length(m)) return(NULL)
  tgt <- vapply(m, function(x) {
    if (grepl("repo_path", x)) {
      a <- sub(".*repo_path\\s*\\(", "", x)
      paste(gsub("[\"' )]", "", strsplit(a, ",")[[1]]), collapse = "/")
    } else sub(".*[\"']([^\"']+)[\"'].*", "\\1", x)
  }, character(1))
  data.frame(from = p, to = gsub("^[.]{2}/|^[.]/", "", tgt),
             stringsAsFactors = FALSE)
}))
edges <- unique(edges[nzchar(edges$to), , drop = FALSE])
edges$to_exists <- file.exists(edges$to)
edges$from_layer <- vapply(edges$from, LAYER, character(1))
edges$to_layer <- vapply(edges$to, function(p)
  if (file.exists(p)) LAYER(p) else "unresolved", character(1))
rownames(edges) <- NULL
utils::write.csv(edges, file.path(PH_TAB, "repository_dependency_edges.csv"),
                 row.names = FALSE)

fan_in <- sort(table(edges$to[edges$to_exists]), decreasing = TRUE)

# ================================================= B15 structural anti-patterns
AP <- function(pattern, severity, detail, evidence)
  data.frame(anti_pattern = pattern, severity = severity, detail = detail,
             evidence = evidence, stringsAsFactors = FALSE)
ap <- list()

# B15A near-duplicate stage directories
stage_dirs <- sort(unique(inv$top_level[inv$layer == "numbered analysis stage"]))
theme <- sub("^[0-9]{2}_", "", stage_dirs)
dupe <- stage_dirs[duplicated(substr(theme, 1, 10)) |
                   duplicated(substr(theme, 1, 10), fromLast = TRUE)]
if (length(dupe)) ap[[length(ap) + 1L]] <- AP(
  "near-duplicate stage directories", "P2 clarity",
  "two numbered stages share a theme prefix, so a reader must guess which is current",
  paste(sprintf("%s (%d scripts)", dupe,
                vapply(dupe, function(d) sum(inv$top_level == d), integer(1))),
        collapse = "; "))

# B15B analysis scripts sitting at the repository root
root_scripts <- inv$path[inv$layer == "repository root"]
if (length(root_scripts)) ap[[length(ap) + 1L]] <- AP(
  "analysis scripts at the repository root", "P2 clarity",
  "root scripts have no stage, so their run order is not expressible in the registry",
  paste(root_scripts, collapse = "; "))

# B15C unregistered scripts inside numbered stages
unreg <- inv$path[inv$layer == "numbered analysis stage" &
                  !inv$registered_in_pipeline]
if (length(unreg)) ap[[length(ap) + 1L]] <- AP(
  "unregistered script inside a numbered stage", "P1 reachability",
  "lives where a pipeline step lives but is not a registry step, so it can be run by hand and is never validated",
  paste(unreg, collapse = "; "))

# B15D registry steps whose file is gone
missing <- ENTRIES$script[!file.exists(ENTRIES$script)]
if (length(missing)) ap[[length(ap) + 1L]] <- AP(
  "registry step with no file", "P0 correctness",
  "the registry names a script that does not exist",
  paste(missing, collapse = "; "))

# B15E unresolved source() targets. The test suite quotes source paths inside
# assertion strings, so a literal there is test data rather than a real edge.
bad_edge <- edges[!edges$to_exists & edges$from_layer != "test suite",
                  , drop = FALSE]
if (nrow(bad_edge)) ap[[length(ap) + 1L]] <- AP(
  "unresolved source() target", "P1 reachability",
  "a script sources a path that does not resolve from the repository root",
  paste(sprintf("%s -> %s", bad_edge$from, bad_edge$to), collapse = "; "))

# B15F deprecated or testing code reachable from active code
leak <- edges[edges$to_exists &
              edges$to_layer %in% c("deprecated", "testing scaffold") &
              !edges$from_layer %in% c("deprecated", "testing scaffold",
                                       "audit layer", "test suite"),
              , drop = FALSE]
ap[[length(ap) + 1L]] <- AP(
  "active code sourcing deprecated or testing code",
  if (nrow(leak)) "P0 correctness" else "none - clean",
  "a deprecated or scaffold file must never be reachable from a publication path",
  if (nrow(leak)) paste(sprintf("%s -> %s", leak$from, leak$to),
                        collapse = "; ") else "0 edges")

ap <- do.call(rbind, ap)
utils::write.csv(ap, file.path(PH_TAB, "repository_anti_patterns.csv"),
                 row.names = FALSE)

# =============================================== B19 script naming audit
#
# The numbered stages use NN_name.r; the figure layer uses a layer prefix. A
# name that does not announce its stage or layer is a navigation cost, not a
# defect, so everything here is P2.
nm <- data.frame(path = scripts, file = basename(scripts),
                 layer = inv$layer, stringsAsFactors = FALSE)
nm$has_numeric_prefix <- grepl("^[0-9]{2}_", nm$file)
nm$extension <- ifelse(grepl("[.]R$", nm$file), "R", "r")
nm$convention <- ifelse(
  nm$layer == "numbered analysis stage",
  ifelse(nm$has_numeric_prefix, "OK - NN_name", "OFF - no step number"),
  ifelse(nm$layer %in% c("shared helper library", "figure / manuscript layer",
                         "test suite", "audit layer"),
         "OK - layer convention", "n/a"))
utils::write.csv(nm, file.path(PH_TAB, "repository_script_naming_audit.csv"),
                 row.names = FALSE)
ext_mix <- tapply(nm$extension, nm$layer, function(x) length(unique(x)))

# ================================ B21 one authoritative script per output
#
# Two scripts that write the same path are the concrete form of "more than one
# plausible authoritative script". Writers are detected from the literal path
# arguments each script passes to a write function.
W <- "(write[._]csv[a-z_]*|write[._]tsv|writeLines|ggsave|write[._]yaml|saveRDS)"
writers <- do.call(rbind, lapply(scripts, function(p) {
  ln <- readLines(p, warn = FALSE)
  i <- grep(W, ln, perl = TRUE)
  if (!length(i)) return(NULL)
  f <- unlist(regmatches(ln[i], gregexpr(
    "[\"'][A-Za-z0-9_./-]+[.](csv|tsv|svg|pdf|png|md|yml|rds)[\"']", ln[i])))
  if (!length(f)) return(NULL)
  data.frame(script = p, output_file = gsub("[\"']", "", f),
             stringsAsFactors = FALSE)
}))
writers <- unique(writers)
multi <- aggregate(script ~ output_file, writers, function(x)
  paste(sort(unique(x)), collapse = " | "))
multi$n_writers <- vapply(strsplit(multi$script, " [|] "), length, integer(1))
multi <- multi[multi$n_writers > 1, , drop = FALSE]
names(multi)[names(multi) == "script"] <- "writing_scripts"
# A literal without a separator is only the basename; its directory is supplied
# at runtime, so two scripts sharing it may well write to different places. Only
# a full relative path proves a contested output.
multi$evidence_class <- ifelse(grepl("/", multi$output_file),
                               "contested - identical relative path",
                               "inconclusive - shared basename only")
multi <- multi[order(multi$evidence_class, -multi$n_writers), , drop = FALSE]
rownames(multi) <- NULL
utils::write.csv(multi, file.path(PH_TAB,
  "repository_contested_outputs.csv"), row.names = FALSE)

cat("\n===== PART B INVENTORY =====\n")
cat("B13 tracked R scripts:", nrow(inv),
    "| registered pipeline steps:", sum(inv$registered_in_pipeline),
    "| registry-active but unregistered:",
    sum(inv$registry_considers_active & !inv$registered_in_pipeline), "\n")
print(table(inv$layer))
cat("\nB14 source() edges:", nrow(edges),
    "| unresolved:", sum(!edges$to_exists), "\n")
cat("most-depended-on files:\n"); print(head(fan_in, 8))
cat("\nB15 anti-patterns:\n")
print(ap[, c("anti_pattern", "severity")])
cat("\nB19 extension mixing per layer (1 = consistent):\n"); print(ext_mix)
cat("off-convention names in numbered stages:",
    sum(nm$convention == "OFF - no step number"), "\n")
cat("\nB21 output files written by more than one script:", nrow(multi), "\n")
print(table(multi$evidence_class))
cont <- multi[multi$evidence_class == "contested - identical relative path", ]
if (nrow(cont)) print(cont[, c("output_file", "n_writers", "writing_scripts")])
