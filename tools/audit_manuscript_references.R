#!/usr/bin/env Rscript

# Manuscript figure-reference integrity audit.
#
# Every figure and Extended Data reference in the manuscript must resolve to a
# panel that is declared in the manuscript contract, belongs to a CANONICAL
# publication, and exists as a rendered artefact. The reverse must also hold: a
# canonical numbered figure that nothing cites is a submission defect too, and
# one that is easy to create by promoting a figure the prose never mentions.
#
# This replaces an earlier scratch audit that matched only "Fig. Nx" and
# therefore audited ZERO Extended Data citations while reporting zero unresolved.
# A detector that cannot see the class of reference you are worried about is
# worse than none, because it produces a green result.

source(file.path("R", "paths.R"))

OUT <- repo_path("manuscript", "figure_panel_reference_audit.csv")
contract <- yaml::read_yaml(repo_path("figures", "figure_contract.yml"))
registry <- utils::read.csv(repo_path("manuscript", "canonical_publication_registry.csv"),
                            stringsAsFactors = FALSE)

# panel id -> publication identity, for every numbered publication
decl <- do.call(rbind, lapply(names(contract$figures), function(k) {
  f <- contract$figures[[k]]
  if (!isTRUE(f$is_numbered_manuscript_figure)) return(NULL)
  pid <- as.character(f$canonical_publication_id)
  data.frame(contract_key = k, publication_id = pid,
             is_ed = !is.null(f$extended_data_number),
             number = if (is.null(f$extended_data_number)) as.integer(k) else as.integer(f$extended_data_number),
             panel_id = vapply(f$panels, function(p) as.character(p$id), character(1)),
             letter = sub("^[0-9]+", "", vapply(f$panels, function(p) as.character(p$id), character(1))),
             stringsAsFactors = FALSE)
}))

status_of <- function(pid) {
  i <- match(pid, registry$publication_id)
  if (is.na(i)) "NOT_IN_REGISTRY" else registry$status[i]
}

files <- c(repo_path("manuscript", "manuscript_draft.md"),
           list.files(repo_path("manuscript"), pattern = "_legend[.]md$", full.names = TRUE),
           list.files(repo_path("manuscript"), pattern = "[.]csv$", full.names = TRUE))
# The audit must not read its own output, and the records that EXIST in order to
# describe defective or withheld references - the promotion audit, the readiness
# register, the blocker list - are not citation surfaces. Including them makes
# the audit report a defect for correctly recording a defect.
NOT_CITATION_SURFACES <- c("figure_panel_reference_audit.csv",
                           "extended_data_promotion_audit.csv",
                           "submission_readiness.csv",
                           "figure_promotion_blockers.csv",
                           "canonical_publication_registry.csv",
                           "figure_generation_inventory.csv",
                           "figure_claim_support_matrix.csv",
                           "prerestructure_freeze_manifest.csv")
files <- files[file.exists(files) & !(basename(files) %in% NOT_CITATION_SURFACES)]

expand <- function(tail_) {
  ls_ <- regmatches(tail_, gregexpr("[a-i]", tail_))[[1]]
  idx <- match(ls_, letters)
  if (grepl("[–-]", tail_) && length(idx) >= 2L && !anyNA(idx)) {
    ls_ <- letters[seq(min(idx), max(idx))]
  }
  ls_
}

rows <- list()
for (f in files) {
  ln <- readLines(f, warn = FALSE)
  for (i in seq_along(ln)) {
    # Extended Data first, so its "Fig. N" is not mistaken for a main figure
    for (pat in c(ed = "Extended Data Fig[.][ ]?[0-9]+[a-i]?(([,–-])[a-i])*",
                  main = "(?<!Extended Data )Fig[.][ ]?[123][a-i](([,–-])[a-i])*")) {
      hits <- regmatches(ln[i], gregexpr(pat, ln[i], perl = TRUE))[[1]]
      is_ed <- identical(pat, unname(pat["ed"])) || grepl("^Extended", pat)
      for (hh in hits) {
        is_ed_hit <- grepl("^Extended Data", hh)
        num <- as.integer(sub("^(Extended Data )?Fig[.][ ]?([0-9]+).*$", "\\2", hh))
        tail_ <- sub("^(Extended Data )?Fig[.][ ]?[0-9]+", "", hh)
        ls_ <- expand(tail_)
        if (!length(ls_)) ls_ <- NA_character_          # whole-figure citation
        for (L in ls_) {
          d <- decl[decl$is_ed == is_ed_hit & decl$number == num, , drop = FALSE]
          pid <- if (nrow(d)) unique(d$publication_id)[1] else NA_character_
          ok_panel <- if (is.na(L)) nrow(d) > 0L else any(d$letter == L)
          art <- if (!is.na(pid)) repo_path("results", "figures", "manuscript", pid,
                                            "assembled", paste0(pid, ".svg")) else NA_character_
          rows[[length(rows) + 1L]] <- data.frame(
            source_file = basename(f), citation = hh,
            kind = if (is_ed_hit) "extended_data" else "main_figure",
            publication_id = pid, panel_letter = L,
            declared_in_contract = ok_panel,
            publication_status = if (is.na(pid)) "NO_SUCH_FIGURE" else status_of(pid),
            exported_artifact_exists = !is.na(art) && file.exists(art),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
}
aud <- unique(do.call(rbind, rows))
aud <- aud[order(aud$kind, aud$publication_id, aud$panel_letter, aud$source_file), ]

aud$resolved <- aud$declared_in_contract & aud$publication_status == "CANONICAL" &
  aud$exported_artifact_exists

# the reverse defect: a canonical publication nothing cites
cited <- unique(aud$publication_id[aud$resolved])
canon <- registry$publication_id[registry$status == "CANONICAL"]
uncited <- setdiff(canon, cited)

con <- file(OUT, open = "wb")
utils::write.csv(aud, con, row.names = FALSE, na = "", eol = "\n")
close(con)

cat("manuscript reference audit ->", OUT, "\n")
cat("  citation instances :", nrow(aud), "\n")
cat("    main figure      :", sum(aud$kind == "main_figure"), "\n")
cat("    extended data    :", sum(aud$kind == "extended_data"), "\n")
cat("  UNRESOLVED         :", sum(!aud$resolved), "\n")
if (any(!aud$resolved)) print(unique(aud[!aud$resolved,
  c("source_file", "citation", "publication_status", "declared_in_contract")]), row.names = FALSE)
cat("  canonical but never cited:", if (length(uncited)) paste(uncited, collapse = ", ") else "none", "\n")
