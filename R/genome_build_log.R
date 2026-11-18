# Requested versus built genome features, checked before and after building.
#
# simGenome() places uORFs by random sampling inside each leader, so a leader can
# end up with fewer uORFs than asked for. Two things make that visible instead of
# silent: a feasibility check that stops before any sequence is made when a gene's
# request cannot be met at all, and a per-gene build log written next to the genome.

# Necessary conditions for a gene's uORF request; a gene that fails one of them can
# never get its uORFs, whatever the random draws. Passing them does not guarantee
# success (placement is random), which is why the build is checked again afterwards.
uorf_feasibility <- function(n_uorfs, leader_length, uorf_max_length,
                             uorfs_can_overlap, uorfs_can_overlap_cds, tx_names) {
  k <- as.integer(n_uorfs)
  L <- as.integer(leader_length)
  reason <- rep(NA_character_, length(k))
  # A uORF needs its start codon inside the leader and at least 6 nt (start + stop).
  shortest_leader_only <- 6L
  # uORFs that may not overlap need room side by side; one of them may end in the CDS
  # (then only its start codon and one more codon must fit in the leader).
  into_cds <- if (uorfs_can_overlap_cds >= 1 && uorf_max_length >= 12) 1L else 0L
  if (uorfs_can_overlap_cds == 2) {
    # Every uORF must run into the CDS; without overlap only one can.
    bad <- k > 1L & !uorfs_can_overlap
    reason[bad] <- sprintf("%d uORFs that must all run into the CDS cannot be placed without overlapping each other (uorfs_can_overlap = FALSE)", k[bad])
    need <- ifelse(k > 0L, 4L + pmax(k - 1L, 0L), 0L)
  } else if (!uorfs_can_overlap) {
    leader_only <- pmax(k - into_cds, 0L)
    need <- leader_only * shortest_leader_only + ifelse(k > leader_only, 4L, 0L)
  } else {
    # Overlapping uORFs still need their own start codon position each.
    need <- ifelse(k > 0L, pmax(k + 3L, ifelse(into_cds == 1L, 4L, shortest_leader_only)), 0L)
  }
  short <- is.na(reason) & k > 0L & L < need
  reason[short] <- sprintf("%d uORF(s) need a leader of at least %d nt, got %d",
                           k[short], need[short], L[short])
  data.table::data.table(transcript_id = tx_names, requested_uorfs = k,
                         leader_length = L, min_leader_needed = need,
                         feasible = is.na(reason), reason = reason)
}

# Stop (or warn) before building when some genes cannot get their uORFs at all.
stop_if_infeasible <- function(feasibility, on_infeasible = c("error", "warning")) {
  on_infeasible <- match.arg(on_infeasible)
  bad <- feasibility[feasible == FALSE]
  if (!nrow(bad)) return(invisible(TRUE))
  shown <- head(bad, 10)
  text <- paste0(nrow(bad), " gene(s) cannot get the requested uORFs:\n",
                 paste0("  ", shown$transcript_id, ": ", shown$reason, collapse = "\n"),
                 if (nrow(bad) > 10) paste0("\n  ... and ", nrow(bad) - 10, " more") else "",
                 "\nUse longer leaders, fewer uORFs, or allow overlaps.")
  if (on_infeasible == "error") stop(text, call. = FALSE)
  warning(text, call. = FALSE)
  invisible(FALSE)
}

# Count built uORFs per transcript and add them to the feasibility table.
uorf_build_log <- function(feasibility, uorf_ranges, dropped_invalid = 0L) {
  built <- if (length(uorf_ranges)) table(names(uorf_ranges)) else integer(0)
  log <- data.table::copy(feasibility)
  log[, built_uorfs := as.integer(ifelse(transcript_id %in% names(built),
                                         built[transcript_id], 0L))]
  log[, status := data.table::fifelse(built_uorfs == requested_uorfs, "ok",
                    data.table::fifelse(feasible, "fewer_than_requested", "infeasible"))]
  log[status == "fewer_than_requested", reason := paste0(
    "placed ", built_uorfs, " of ", requested_uorfs,
    " (random placement found no room, or an overlap broke a uORF's start/stop and it was dropped)")]
  log[]
}

write_genome_build_log <- function(log, genome_path) {
  path <- sub("\\.fasta$", "_build_log.tsv", genome_path)
  data.table::fwrite(log, path, sep = "\t")
  path
}
