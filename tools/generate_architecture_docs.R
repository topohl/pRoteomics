#!/usr/bin/env Rscript

# Generate docs/ANALYSIS_ENTRYPOINTS.md: one row per canonical analysis area,
# answering "which script owns this analysis, what does it need, and who reads
# what it produces".
#
# Derived from pipeline.yml so it cannot drift from the registry.
#
# Result ownership is a different contract with its own generator:
# tools/generate_results_ownership.R writes docs/RESULTS_OWNERSHIP.md and
# config/results_ownership.csv. One generator per document.

source(file.path("R", "paths.R"))
source(repo_path("R", "dataset_config.R"))
source(repo_path("R", "pipeline_registry.R"))

registry <- read_pipeline_registry(repo_path("pipeline.yml"))
steps <- pipeline_steps(registry, pipeline_stage_names(registry),
                        dataset = "all", include_unsupported = TRUE)
steps <- steps[!duplicated(steps$script), , drop = FALSE]

## Analysis identity comes from the analysis/ directory a script lives in.
area_of <- function(script) {
  m <- regmatches(script, regexpr("^analysis/[0-9]{2}[a-z]?_[A-Za-z_]+", script))
  if (!length(m) || !nzchar(m)) return(NA_character_)
  sub("^analysis/", "", m)
}
steps$area <- vapply(steps$script, area_of, character(1))

AREA_TITLE <- c(
  "01_preprocessing"        = "Preprocessing and identifier mapping",
  "02_qc"                   = "Quality control and marker fidelity",
  "03_spatial_validation"   = "Spatial systems validation, bilateral QC and CA2-SLM robustness",
  "04_differential_abundance" = "Differential abundance and enrichment",
  "05_wgcna"                = "WGCNA modules and supermodules",
  "06_gsea"                 = "Cell-type enrichment (EWCE/GSEA)",
  "07_spatial_networks"     = "Spatial and differential networks",
  "08_integration"          = "Biological integration and behaviour coupling",
  "09_publication_exports"  = "Publication source data and PRIDE export"
)

## pipeline_steps() collapses list-valued fields with "|"; a few hand-written
## registry entries use ";". Split on both, or paths get glued together.
split_paths <- function(x) {
  p <- unlist(strsplit(paste(x, collapse = "|"), "[|;]"))
  p <- trimws(p)
  unique(p[nzchar(p)])
}

fmt <- function(x, n = 3L) {
  if (!length(x)) return("-")
  out <- paste0("`", utils::head(x, n), "`", collapse = ", ")
  if (length(x) > n) out <- paste0(out, ", +", length(x) - n, " more")
  out
}

## ---- path coverage -------------------------------------------------------
## Registry paths carry <dataset> placeholders and globs. A produced entry
## covers a consumed entry when the consumed path matches it exactly under
## wildcard expansion, or when the produced entry is a directory the consumed
## path sits beneath. Nothing looser: a loose rule invents dependencies.
as_rx <- function(p) {
  q <- gsub("([.^$+?(){}\\[\\]\\\\])", "\\\\\\1", p)
  q <- gsub("<[A-Za-z_]+>", "[^/]*", q)
  q <- gsub("\\*\\*", "\001", q)
  q <- gsub("\\*", "[^/]*", q)
  q <- gsub("\001", ".*", q)
  paste0("^", q, "$")
}
literal_stem <- function(p) sub("(<|\\*).*$", "", p)

covers <- function(produced, consumed) {
  if (grepl(as_rx(produced), consumed)) return(TRUE)
  if (grepl("/$", produced)) return(startsWith(consumed, literal_stem(produced)))
  FALSE
}

area_names <- names(AREA_TITLE)
area_names <- area_names[area_names %in% steps$area]

produced_by <- lapply(area_names, function(a)
  split_paths(steps$produces[!is.na(steps$area) & steps$area == a]))
## Dependency edges come from required inputs only. Including optional inputs
## makes the graph cyclic (WGCNA optionally reads behaviour-coupling output,
## which requires WGCNA), which is true but useless as a dependency contract.
consumed_by <- lapply(area_names, function(a)
  split_paths(steps$consumes_required[!is.na(steps$area) & steps$area == a]))
names(produced_by) <- names(consumed_by) <- area_names

## upstream[a] = areas producing something area a consumes
upstream <- lapply(area_names, function(a) {
  cons <- consumed_by[[a]]
  hits <- vapply(area_names, function(b) {
    if (identical(a, b)) return(FALSE)
    any(vapply(produced_by[[b]], function(p)
      any(vapply(cons, function(c) covers(p, c), logical(1))), logical(1)))
  }, logical(1))
  area_names[hits]
})
names(upstream) <- area_names
downstream <- lapply(area_names, function(a)
  area_names[vapply(area_names, function(b) a %in% upstream[[b]], logical(1))])
names(downstream) <- area_names

## ---- config a script actually reads --------------------------------------
config_of <- function(scripts) {
  out <- character(0)
  for (s in scripts) {
    if (!file.exists(s)) next
    txt <- paste(readLines(s, warn = FALSE), collapse = "\n")
    m <- unlist(regmatches(txt, gregexpr('config/[A-Za-z0-9_./-]+', txt)))
    m2 <- unlist(regmatches(txt, gregexpr('repo_path\\(\\s*"config"\\s*,\\s*"[^"]+"', txt)))
    m2 <- paste0("config/", sub('.*"([^"]+)"$', "\\1", m2))
    out <- c(out, m, m2)
  }
  out <- unique(out[grepl("[.](ya?ml|csv|tsv|json)$", out)])
  sort(out)
}

PUB_RX <- "^results/(publication_)?source_data/"

## ------------------------------------------------------ ANALYSIS_ENTRYPOINTS
lines <- c(
  "# Canonical analysis entrypoints",
  "",
  "Generated by `tools/generate_architecture_docs.R` from `pipeline.yml`.",
  "Do not hand-edit.",
  "",
  "`pipeline.yml` remains the sole execution-order authority. This document",
  "answers a different question: which script owns a given analysis, what that",
  "analysis needs, and who consumes what it produces.",
  "",
  paste("Result-file ownership is a separate contract, in",
        "`docs/RESULTS_OWNERSHIP.md` and `config/results_ownership.csv`."),
  "",
  paste("Manuscript rendering is intentionally not listed here. It lives in the",
        "Exp9_manuscript repository and consumes the frozen bundles under",
        "`results/publication_source_data/`; it never reads any other path in",
        "this repository."),
  "",
  paste("`upstream_dependencies` and `downstream_consumers` are computed by",
        "matching produced paths against consumed paths across areas, so they",
        "reflect the registry rather than a hand-maintained diagram."),
  ""
)

for (a in area_names) {
  s <- steps[!is.na(steps$area) & steps$area == a, , drop = FALSE]
  req <- s[s$required %in% TRUE, , drop = FALSE]
  lead <- if (nrow(req)) req$script[1] else s$script[1]
  outs <- produced_by[[a]]
  pub <- grep(PUB_RX, outs, value = TRUE)
  cfg <- config_of(s$script)
  envs <- split_paths(s$env)

  lines <- c(lines,
    paste0("## ", a),
    "",
    paste0("**", AREA_TITLE[[a]], "**"),
    "",
    "| field | value |",
    "| --- | --- |",
    paste0("| analysis_id | `", a, "` |"),
    paste0("| human_name | ", AREA_TITLE[[a]], " |"),
    paste0("| entrypoint | `", lead, "` |"),
    paste0("| scripts_in_area | ", nrow(s), " (", sum(s$required %in% TRUE), " required) |"),
    paste0("| inputs | ", fmt(split_paths(s$consumes_required)), " |"),
    paste0("| optional_inputs | ", fmt(split_paths(s$consumes_optional)), " |"),
    paste0("| outputs | ", fmt(outs, 4L), " |"),
    paste0("| required_config | ", fmt(c(cfg, envs), 4L), " |"),
    paste0("| upstream_dependencies | ",
           if (length(upstream[[a]])) paste0("`", upstream[[a]], "`", collapse = ", ") else "-", " |"),
    paste0("| downstream_consumers | ",
           paste(c(if (length(downstream[[a]])) paste0("`", downstream[[a]], "`", collapse = ", "),
                   if (length(pub)) "Exp9_manuscript (via frozen source data)"),
                 collapse = ", "),
           if (!length(downstream[[a]]) && !length(pub)) "-" else "", " |"),
    paste0("| publication_source_outputs | ", fmt(pub, 3L), " |"),
    paste0("| dependency_stages | ", paste0("`", unique(s$stage), "`", collapse = ", "), " |"),
    "",
    "Scripts:",
    "",
    paste0("- `", s$script, "` - stage `", s$stage, "`, scope `", s$scope, "`",
           ifelse(s$required %in% TRUE, ", required", "")),
    ""
  )
}
## ---- area-level cycles ---------------------------------------------------
## Two areas exchange files in both directions. That is real, not a matching
## artefact: each direction rests on one exact file. The pipeline is a DAG at
## script level, and pipeline.yml's stage order is what defines execution;
## areas are a naming grouping, so they can legitimately interleave.
mutual <- list()
for (a in area_names) for (b in area_names) {
  if (a >= b) next
  if (a %in% upstream[[b]] && b %in% upstream[[a]]) mutual[[length(mutual) + 1L]] <- c(a, b)
}
links_between <- function(from, to) {
  prod <- produced_by[[from]]
  cons <- consumed_by[[to]]
  keep <- character(0)
  for (p in prod) for (c in cons) if (covers(p, c)) keep <- c(keep, c)
  unique(keep)
}
lines <- c(lines, "## Area-level cycles", "")
if (!length(mutual)) {
  lines <- c(lines, "None: the area graph is acyclic.", "")
} else {
  lines <- c(lines,
    paste("These area pairs exchange files in both directions. Execution order",
          "is defined by `pipeline.yml` stage order, not by area; an area is a",
          "naming grouping and scripts within two areas can legitimately",
          "interleave. Each direction below is created by a specific file."),
    "")
  for (m in mutual) {
    lines <- c(lines, paste0("### `", m[1], "` and `", m[2], "`"), "")
    for (d in list(m, rev(m))) {
      l <- links_between(d[1], d[2])
      lines <- c(lines,
        paste0("- `", d[2], "` requires from `", d[1], "`: ",
               paste0("`", utils::head(l, 3), "`", collapse = ", "),
               if (length(l) > 3) paste0(", +", length(l) - 3, " more") else ""))
    }
    lines <- c(lines, "")
  }
}

writeLines(lines, repo_path("docs", "ANALYSIS_ENTRYPOINTS.md"), useBytes = TRUE)

cat("wrote docs/ANALYSIS_ENTRYPOINTS.md\n")
cat("analysis areas       :", length(area_names), "\n")
cat("registered scripts   :", nrow(steps), "\n")
cat("areas with no upstream:",
    sum(vapply(upstream, function(x) length(x) == 0L, logical(1))), "\n")
cat("areas exporting publication source data:",
    sum(vapply(area_names, function(a)
      any(grepl(PUB_RX, produced_by[[a]])), logical(1))), "\n")
cat("result ownership is generated separately by tools/generate_results_ownership.R\n")
