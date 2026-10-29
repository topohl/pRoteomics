# Shared by the neuropil-lamina analyses of Figure 2j-l (map_input_territories, map_afferent_sources,
# compare_lamina_synapses; 2026-10-28). Each analysis applies the same carry-over rule to the joint
# bundle itself, so none depends on another script's table.

# Somatic carry-over into the DG polymorphic layer. DG po borders the granule cell layer. Soma-only
# proteins (observed in >= 90% of the acquisitions of every soma unit and < 50% of every other unit
# except DG po itself) index carry-over; a DG po acquisition with more soma-only detections than any
# acquisition of the other nine neuropil laminae is excluded.
#   info  data.frame of every acquisition of the bundle: Sample, AnimalID, group, hemisphere, unit
#         (region_layer)
#   obs   logical matrix, proteins x samples (observed and finite), columns in info$Sample order
# Returns one row per neuropil acquisition, with soma_only_detected, the threshold and excluded.
dg_po_carryover <- function(info, obs, laminae, soma_units) {
  stopifnot(identical(colnames(obs), info$Sample), "DG_po" %in% laminae, all(soma_units %in% info$unit))
  units_all <- sort(unique(info$unit))
  det <- vapply(units_all, function(k) rowMeans(obs[, info$unit == k, drop = FALSE]), numeric(nrow(obs)))
  soma_only <- rownames(det)[apply(det[, soma_units, drop = FALSE] >= 0.9, 1, all) &
                             apply(det[, setdiff(units_all, c(soma_units, "DG_po")), drop = FALSE] < 0.5, 1, all)]
  np <- info[info$unit %in% laminae, , drop = FALSE]
  np$soma_only_detected <- colSums(obs[soma_only, np$Sample, drop = FALSE])
  threshold <- max(np$soma_only_detected[np$unit != "DG_po"])
  np$carryover_threshold <- threshold
  np$excluded <- np$unit == "DG_po" & np$soma_only_detected > threshold
  np$n_soma_only_proteins <- length(soma_only)
  np
}

# What a run manifest needs so a reader can tell that a table is the one this run wrote with committed
# code: the SHA-256 of every output and whether the script or these helpers had uncommitted changes.
run_output_provenance <- function(outputs, script_id) {
  files <- unlist(outputs)
  # untracked files count too (this repository hides them from plain git status)
  status <- tryCatch(system2("git", c("-C", repo_root(), "status", "--porcelain", "--untracked-files=all", "--", script_id, "R/spatial"),
                             stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
  list(output_sha256 = as.list(stats::setNames(unname(tools::sha256sum(files)), basename(files))),
       code_uncommitted = if (length(status) == 1L && is.na(status)) NA else length(status) > 0L)
}
