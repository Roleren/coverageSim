# Constructive gene architectures (RDG plan step 2, 2026-11-18).
#
# A gene architecture lists the translons a transcript should contain besides
# its main CDS: uORFs, upstream ORFs overlapping the CDS (uoORF), N-terminal
# extensions (NTE), internal out-of-frame ORFs (intORF), downstream ORFs
# overlapping the CDS stop (doORF) and downstream ORFs in the trailer (dORF).
# Each feature can be "translated" (gets reads later) or a "decoy" (exists in
# the sequence with a proper start and stop, but gets no reads).
#
# Building happens in four steps, all in transcript coordinates:
#   1. place: every feature gets a start/end inside the window its class allows
#      (checked before any sequence exists; an empty window is an error);
#   2. write: start and stop codons are fixed, and every codon inside every
#      feature (and the main CDS) is kept free of stop codons in that feature's
#      frame, so overlapping frames are satisfied together;
#   3. verify: each feature is read back from the final sequence;
#   4. label: every other ORF found in the transcripts is a decoy, classified
#      with the same rules, so ORF callers can be scored for false positives.

feature_classes <- c("uORF", "uoORF", "NTE", "intORF", "doORF", "dORF")

#' Describe the extra translons of a gene
#'
#' One row per feature. Several genes can share one architecture, or
#' \code{simGenome(architecture = )} takes a list with one per gene.
#' @param class feature class: "uORF" (start and stop in the leader), "uoORF"
#'  (start in the leader, stop inside the CDS, other frame), "NTE" (in-frame
#'  upstream start in the leader, ends at the CDS stop), "intORF" (start and stop
#'  inside the CDS, other frame), "doORF" (start inside the CDS, other frame,
#'  stop in the trailer) or "dORF" (start and stop in the trailer).
#' @param length length in nt including start and stop codon, a multiple of 3.
#'  For "NTE" it is the length of the extension (nt upstream of the CDS start).
#' @param start_codon start codon of the feature, e.g. "ATG", "CTG".
#' @param status "translated" (gets reads) or "decoy" (built, but no reads).
#' @param frame for out-of-frame classes, 1 or 2 nt after the CDS frame, or NA
#'  for either.
#' @param may_overlap FALSE keeps this feature from overlapping any other extra
#'  feature of the gene (the main CDS is not counted).
#' @return a data.table, one row per feature
#' @export
#' @examples
#' gene_architecture(c("uORF", "uoORF"), length = c(30, 90))
gene_architecture <- function(class, length, start_codon = "ATG",
                              status = "translated", frame = NA_integer_,
                              may_overlap = TRUE) {
  if (!is.numeric(length) || anyNA(length) || any(length != floor(length)))
    stop("Feature lengths must be whole numbers")
  n <- max(lengths(list(class, length)))
  arch <- data.table::data.table(class = rep_len(as.character(class), n),
                                 length = rep_len(as.integer(length), n),
                                 start_codon = rep_len(toupper(start_codon), n),
                                 status = rep_len(status, n),
                                 frame = rep_len(as.integer(frame), n),
                                 may_overlap = rep_len(as.logical(may_overlap), n))
  bad <- setdiff(arch$class, feature_classes)
  if (length(bad)) stop("Unknown feature class: ", paste(bad, collapse = ", "),
                        "; use one of ", paste(feature_classes, collapse = ", "))
  if (any(!arch$status %in% c("translated", "decoy")))
    stop("status must be 'translated' or 'decoy'")
  if (any(nchar(arch$start_codon) != 3 | grepl("[^ACGT]", arch$start_codon)))
    stop("start_codon must be 3 letters of A, C, G, T")
  if (any(arch$start_codon %in% c("TAA", "TAG", "TGA")))
    stop("start_codon cannot be a stop codon")
  if (any(!is.na(arch$frame) & !arch$frame %in% 1:2))
    stop("frame must be 1, 2 or NA")
  orf <- arch$class != "NTE"
  if (any(arch$length[orf] %% 3 != 0 | arch$length[orf] < 6))
    stop("Feature lengths must be multiples of 3 and at least 6 (start + stop)")
  if (any(arch$length[!orf] %% 3 != 0 | arch$length[!orf] < 3))
    stop("NTE lengths must be multiples of 3 and at least 3")
  arch[]
}

#' Ready-made gene architectures
#'
#' @param name one of "uORF", "two_uORFs", "uoORF", "ATF4_like" (short uORF,
#'  then a uoORF across the CDS start), "NTE_CUG" (near-cognate in-frame
#'  extension), "intORF", "doORF", "dORF", "uORF_decoy" (a translated uORF
#'  plus an untranslated decoy uORF), "all" (one of each class).
#' @return a data.table from \code{\link{gene_architecture}}
#' @export
architecture_preset <- function(name) {
  switch(name,
    uORF = gene_architecture("uORF", 30),
    two_uORFs = gene_architecture(c("uORF", "uORF"), c(30, 45), may_overlap = FALSE),
    uoORF = gene_architecture("uoORF", 90),
    ATF4_like = gene_architecture(c("uORF", "uoORF"), c(12, 180), may_overlap = FALSE),
    NTE_CUG = gene_architecture("NTE", 60, start_codon = "CTG"),
    intORF = gene_architecture("intORF", 90),
    doORF = gene_architecture("doORF", 60),
    dORF = gene_architecture("dORF", 45),
    uORF_decoy = gene_architecture(c("uORF", "uORF"), c(30, 30),
                                   status = c("translated", "decoy"), may_overlap = FALSE),
    all = gene_architecture(c("uORF", "uoORF", "intORF", "doORF", "dORF"),
                            c(24, 90, 60, 60, 45)),
    stop("Unknown preset '", name, "'"))
}

# A table built or edited by hand goes through the same checks as
# gene_architecture(), with its defaults for missing columns.
revalidate_architecture <- function(a) {
  a <- data.table::as.data.table(a)
  if (!all(c("class", "length") %in% names(a)))
    stop("An architecture table needs at least the columns class and length")
  if (!nrow(a)) return(gene_architecture(character(), integer()))
  get <- function(col, default) if (col %in% names(a)) a[[col]] else default
  gene_architecture(a$class, a$length, start_codon = get("start_codon", "ATG"),
                    status = get("status", "translated"), frame = get("frame", NA_integer_),
                    may_overlap = get("may_overlap", TRUE))
}

# Recycle the architecture argument of simGenome to one table per gene.
architecture_per_gene <- function(architecture, n, tx_names) {
  if (data.table::is.data.table(architecture) || is.data.frame(architecture)) {
    architecture <- rep(list(revalidate_architecture(architecture)), n)
  } else if (is.character(architecture)) {
    # Presets are immutable here: rbindlist creates the per-gene result.
    # Construct each distinct preset once instead of once per transcript.
    preset_names <- rep_len(architecture, n)
    unique_names <- unique(preset_names)
    presets <- lapply(unique_names, architecture_preset)
    architecture <- presets[match(preset_names, unique_names)]
  } else if (is.list(architecture)) {
    if (length(architecture) != n)
      stop("architecture as a list needs one entry per gene (", n, "), got ",
           length(architecture))
    architecture <- lapply(architecture, function(a) {
      if (is.null(a) || (is.character(a) && !length(a))) return(gene_architecture(character(), integer()))
      if (is.character(a)) architecture_preset(a) else revalidate_architecture(a)
    })
  } else stop("architecture must be a gene_architecture() table, a preset name, or a list of them")
  arch <- data.table::rbindlist(architecture, idcol = "gene", fill = TRUE)
  if (!nrow(arch)) return(arch)
  arch[, transcript_id := tx_names[gene]]
  arch[, feature_id := paste0(transcript_id, "_", class, "_", seq_len(.N)), by = gene]
  arch[]
}

# Where each class may start, in transcript coordinates (1-based, the first
# base of the start codon), given leader length L and CDS end C (last base of the
# CDS stop codon) and transcript length Tn. Out-of-frame classes keep their
# start and stop codons off the CDS start and stop codons.
feature_start_window <- function(class, len, L, C, Tn) {
  switch(class,
    uORF   = c(1L, L - len + 1L),
    uoORF  = c(max(1L, L + 7L - len), min(L - 2L, C - 2L - len)),
    NTE    = if (len <= L) rep(L + 1L - len, 2L) else c(1L, 0L),
    intORF = c(L + 4L, C - 2L - len),
    doORF  = c(L + 4L, min(C - 5L, Tn - len + 1L), C - len + 2L),
    dORF   = c(C + 1L, Tn - len + 1L))
}

# Frame of a start position relative to the CDS (0 = in frame).
cds_frame <- function(start, L) as.integer((start - (L + 1L)) %% 3L)

# One gene, one attempt (plain vectors: this runs once per gene). Returns
# starts/ends (NA when not placed), the problem text, and whether the failure
# is definite (an empty class window).
place_gene_features <- function(class, len, start_codon, want_frame, may_overlap,
                                L, C, Tn, max_tries) {
  k <- length(class)
  tx_start <- tx_end <- rep(NA_integer_, k); problem <- rep(NA_character_, k)
  definite <- FALSE
  fixed <- c(L + 1:3, C - 2:0)
  # Spans already placed; the CDS first.
  ps <- L + 1L; pe <- C; pf <- 0L; pok <- TRUE; pcds <- TRUE
  # Longest first: they have the fewest positions to choose from.
  for (r in order(-len)) {
    cl <- class[r]
    w <- feature_start_window(cl, len[r], L, C, Tn)
    cand <- if (w[1] <= w[2]) seq.int(w[1], w[2]) else integer()
    if (cl == "doORF") cand <- cand[cand >= w[3]] # stop must reach the trailer
    if (cl %in% c("uoORF", "intORF", "doORF")) {
      f <- cds_frame(cand, L)
      cand <- cand[f != 0L & (is.na(want_frame[r]) | f == want_frame[r])]
    }
    is_nte <- cl == "NTE"
    if (is_nte && start_codon[r] %in% c("TAA", "TAG", "TGA")) cand <- integer()
    if (!length(cand)) {
      problem[r] <- sprintf("%s of %d nt does not fit (leader %d, CDS %d, trailer %d nt)",
                            cl, len[r], L, C - L, Tn - C)
      definite <- TRUE
      next
    }
    ok <- FALSE
    for (try in seq_len(min(max_tries, length(cand) * 3L))) {
      s <- cand[sample.int(length(cand), 1L)]
      e <- if (is_nte) C else s + len[r] - 1L
      codons <- if (is_nte) s + 0:2 else c(s + 0:2, e - 2:0)
      if (any(codons %in% fixed)) next
      fr <- cds_frame(s, L)
      ov <- ps <= e & pe >= s
      # Same frame overlap would put one stop inside the other (NTE shares the CDS).
      if (any(ov & pf == fr & !(is_nte & pcds))) next
      if (!may_overlap[r] && any(ov & !pcds)) next
      if (any(ov & !pok & !pcds)) next
      ok <- TRUE; break
    }
    if (!ok) {
      problem[r] <- sprintf("no room for %s of %d nt next to the other features", cl, len[r])
      next
    }
    tx_start[r] <- s; tx_end[r] <- e
    fixed <- c(fixed, codons)
    ps <- c(ps, s); pe <- c(pe, e); pf <- c(pf, fr); pok <- c(pok, may_overlap[r])
    pcds <- c(pcds, FALSE)
  }
  list(tx_start = tx_start, tx_end = tx_end, problem = problem, definite = definite)
}

# Place all features of all genes; returns the architecture with tx_start/tx_end,
# or stops listing every gene whose request cannot be met.
place_features <- function(arch, leader_length, cds_total_length, trailer_length,
                           max_tries = 500L, restarts = 50L) {
  if (!nrow(arch)) return(arch)
  arch <- data.table::copy(arch)
  arch[, `:=`(tx_start = NA_integer_, tx_end = NA_integer_, problem = NA_character_)]
  rows_by_gene <- split(seq_len(nrow(arch)), arch$gene)
  cls <- arch$class; len <- arch$length; sc <- arch$start_codon
  fr <- arch$frame; mo <- arch$may_overlap
  new_start <- new_end <- rep(NA_integer_, nrow(arch)); new_problem <- rep(NA_character_, nrow(arch))
  for (gene_rows in rows_by_gene) {
    g <- arch$gene[gene_rows[1]]
    L <- as.integer(leader_length[g]); C <- L + as.integer(cds_total_length[g])
    Tn <- C + as.integer(trailer_length[g])
    # Greedy random placement can block itself (a long feature taking the room
    # a short one needed), so a stuck gene is placed again from scratch before
    # it counts as impossible. A class window that is empty never retries.
    for (attempt in seq_len(restarts)) {
      res <- place_gene_features(cls[gene_rows], len[gene_rows], sc[gene_rows],
                                 fr[gene_rows], mo[gene_rows], L, C, Tn, max_tries)
      if (!anyNA(res$tx_start) || res$definite) break
    }
    new_start[gene_rows] <- res$tx_start; new_end[gene_rows] <- res$tx_end
    new_problem[gene_rows] <- res$problem
  }
  arch[, `:=`(tx_start = new_start, tx_end = new_end, problem = new_problem)]
  bad <- arch[!is.na(problem)]
  if (nrow(bad)) {
    shown <- head(bad, 10)
    stop(nrow(bad), " feature(s) cannot be built:\n",
         paste0("  ", shown$transcript_id, ": ", shown$problem, collapse = "\n"),
         if (nrow(bad) > 10) paste0("\n  ... and ", nrow(bad) - 10, " more") else "",
         "\nUse longer regions, shorter features, or allow overlaps.", call. = FALSE)
  }
  arch[, problem := NULL]
  arch[, frame := cds_frame(tx_start, leader_length[gene])]
  arch[]
}

# Codon start positions (global, 1-based) inside features: from the codon after
# the start codon to the codon before the stop codon.
internal_codon_starts <- function(start, end) {
  n <- (end - start + 1L) %/% 3L - 2L
  keep <- n > 0L
  if (!any(keep)) return(integer())
  sequence(n[keep], from = start[keep] + 3L, by = 3L)
}

# Write start/stop codons of all features into the transcripts and remove stop
# codons from every internal codon of every feature and CDS, in every frame at
# once. Only bases that are not part of a fixed start/stop codon are redrawn.
write_feature_sequences <- function(tx_seqs, arch, leader_length, cds_total_length,
                                    stop_codons, max_rounds = 200L) {
  if (!identical(getOption("coverageSim.engine"), "R") && max_rounds > 0L) {
    stops <- if (nrow(arch)) ifelse(arch$class == "NTE", NA_character_,
                    sample(stop_codons, nrow(arch), replace = TRUE)) else character()
    if (nrow(arch)) arch[, stop_codon := stops]
    result <- write_feature_sequences_cpp(tx_seqs, as.integer(leader_length),
      as.integer(cds_total_length), as.integer(arch$gene), as.integer(arch$tx_start),
      as.integer(arch$tx_end), arch$class, arch$start_codon, stops,
      stop_codons, as.integer(max_rounds))
    return(list(tx_seqs = result$tx_seqs, arch = arch, rounds = result$rounds))
  }
  widths <- nchar(tx_seqs)
  offset <- c(0L, cumsum(widths))[seq_along(widths)]
  bases <- strsplit(paste(tx_seqs, collapse = ""), "", fixed = TRUE)[[1]]
  fixed <- logical(length(bases))
  cds_s <- offset + leader_length + 1L
  cds_e <- offset + leader_length + cds_total_length
  # The CDS keeps its own start and stop codons.
  fixed[c(cds_s, cds_s + 1L, cds_s + 2L, cds_e - 2L, cds_e - 1L, cds_e)] <- TRUE
  if (nrow(arch)) {
    g_s <- offset[arch$gene] + arch$tx_start
    g_e <- offset[arch$gene] + arch$tx_end
    stops <- ifelse(arch$class == "NTE", NA_character_,
                    sample(stop_codons, nrow(arch), replace = TRUE))
    arch[, stop_codon := stops]
    for (k in 0:2) {
      bases[g_s + k] <- substr(arch$start_codon, k + 1L, k + 1L)
      fixed[g_s + k] <- TRUE
      orf <- arch$class != "NTE"
      bases[g_e[orf] - 2L + k] <- substr(stops[orf], k + 1L, k + 1L)
      fixed[g_e[orf] - 2L + k] <- TRUE
    }
  } else g_s <- g_e <- integer()
  codon_pos <- internal_codon_starts(c(cds_s, g_s), c(cds_e, g_e))
  codon_pos <- unique(codon_pos)
  for (round in seq_len(max_rounds)) {
    bad <- feature_stop_positions(bases, codon_pos, stop_codons)
    if (!length(bad)) break
    pos <- c(bad, bad + 1L, bad + 2L)
    pos <- unique(pos[!fixed[pos]])
    if (!length(pos)) stop("A stop codon inside a feature is made only of fixed start/stop bases; ",
                           "this placement is impossible (please report)")
    bases[pos] <- sample(c("A", "C", "G", "T"), length(pos), replace = TRUE)
  }
  if (length(bad)) stop("Could not remove internal stop codons in ", max_rounds, " rounds")
  big <- paste(bases, collapse = "")
  list(tx_seqs = substring(big, offset + 1L, offset + widths), arch = arch, rounds = round)
}

# Read each feature back: start codon, stop codon, and no stop in between.
verify_features <- function(tx_seqs, arch, stop_codons, cds_end) {
  if (!nrow(arch)) return(arch)
  seq_of <- substring(tx_seqs[arch$gene], arch$tx_start, arch$tx_end)
  start_ok <- substr(seq_of, 1, 3) == arch$start_codon
  stop_ok <- substring(seq_of, nchar(seq_of) - 2L, nchar(seq_of)) %in% stop_codons
  internal_ok <- vapply(seq_along(seq_of), function(i) {
    x <- seq_of[i]; n <- nchar(x) %/% 3L
    if (n <= 2L) return(TRUE)
    !any(substring(x, seq(4L, by = 3L, length.out = n - 2L),
                   seq(6L, by = 3L, length.out = n - 2L)) %in% stop_codons)
  }, logical(1))
  in_frame <- (arch$tx_end - arch$tx_start + 1L) %% 3L == 0L
  arch[, verified := start_ok & stop_ok & internal_ok & in_frame &
         (class != "NTE" | tx_end == cds_end[gene])]
  if (any(!arch$verified)) {
    bad <- arch[verified == FALSE]
    stop("Built features failed verification: ",
         paste(head(bad$feature_id, 10), collapse = ", "), call. = FALSE)
  }
  arch[]
}

# Class of any ORF by where its start and stop fall (same rules as the builder).
classify_orfs <- function(start, end, L, C) {
  fr <- cds_frame(start, L)
  data.table::fcase(
    end <= L, "uORF",
    start <= L & end == C & fr == 0L, "NTE",
    start <= L & end > C, "spanning",
    start <= L, "uoORF",
    start == L + 1L & end == C, "CDS",
    start <= C & end <= C & fr == 0L, "N_truncation",
    start <= C & end <= C, "intORF",
    start <= C, "doORF",
    default = "dORF")
}

# Every ORF in the final transcripts that is not a built feature or the CDS.
label_decoy_orfs <- function(tx_seqs, arch, tx_names, leader_length, cds_total_length,
                             start_codons, stop_codons) {
  seqs <- Biostrings::DNAStringSet(tx_seqs)
  names(seqs) <- seq_along(seqs)
  orfs <- ORFik::findORFs(seqs, startCodon = paste(start_codons, collapse = "|"),
                          stopCodon = paste(stop_codons, collapse = "|"),
                          longestORF = FALSE, minimumLength = 0)
  if (!length(unlist(orfs))) return(data.table::data.table())
  gene <- as.integer(rep(names(orfs), lengths(orfs)))
  orfs <- unlist(orfs, use.names = FALSE)
  dt <- data.table::data.table(gene = gene, tx_start = start(orfs), tx_end = end(orfs))
  L <- leader_length[dt$gene]; C <- L + cds_total_length[dt$gene]
  dt[, class := classify_orfs(tx_start, tx_end, L, C)]
  dt <- dt[class != "CDS"]
  if (nrow(arch)) {
    built <- paste(arch$gene, arch$tx_start, arch$tx_end)
    dt <- dt[!paste(gene, tx_start, tx_end) %in% built]
  }
  dt[, `:=`(transcript_id = tx_names[gene],
            start_codon = substr(tx_seqs[gene], tx_start, tx_start + 2L),
            frame = cds_frame(tx_start, leader_length[gene]),
            status = "decoy", origin = "incidental")]
  dt[, feature_id := paste0(transcript_id, "_decoy_", seq_len(.N)), by = gene]
  # A decoy that ends at the stop of a built feature or the CDS is a nested
  # in-frame truncation of it, not an independent ORF.
  owners <- data.table::data.table(gene = seq_along(tx_seqs), tx_end = leader_length + cds_total_length,
                                   owner = "CDS")
  if (nrow(arch)) owners <- rbind(owners, arch[, .(gene, tx_end, owner = feature_id)])
  dt[owners, nested_in := i.owner, on = .(gene, tx_end)]
  dt[]
}

# Replace the spliced CDS inside the genomic CDS string (exons + introns), using
# the exon layout of create_cds_ranges(): cds_exons - 1 exons of
# floor(width / cds_exons) nt, the last exon takes the rest.
splice_cds_back <- function(cds_string_genomic, cds_string, cds_exons, intron_length) {
  if (cds_exons <= 1) return(Biostrings::DNAStringSet(cds_string))
  genomic <- as.character(cds_string_genomic)
  intron_length <- rep_len(intron_length, length(genomic))
  out <- vapply(seq_along(genomic), function(i) {
    cds <- cds_string[i]; w <- nchar(cds); split <- floor(w / cds_exons)
    pieces <- character(2L * cds_exons - 1L)
    exon_starts <- (seq_len(cds_exons) - 1L) * split + 1L
    exon_ends <- c(exon_starts[-1] - 1L, w)
    g_pos <- 1L
    for (k in seq_len(cds_exons)) {
      pieces[2L * k - 1L] <- substr(cds, exon_starts[k], exon_ends[k])
      g_pos <- g_pos + exon_ends[k] - exon_starts[k] + 1L
      if (k < cds_exons) {
        pieces[2L * k] <- substr(genomic[i], g_pos, g_pos + intron_length[i] - 1L)
        g_pos <- g_pos + intron_length[i]
      }
    }
    paste(pieces, collapse = "")
  }, character(1))
  Biostrings::DNAStringSet(out)
}

# Map transcript-coordinate features to the genome through the exons.
features_to_genome <- function(features, leader_ranges, exon_cds_ranges, trailer_ranges,
                               tx_names) {
  if (!nrow(features)) return(GenomicRanges::GRangesList())
  strip <- function(x) {
    S4Vectors::mcols(x) <- S4Vectors::DataFrame(transcript_id = x$transcript_id); x
  }
  exons <- c(strip(leader_ranges), strip(exon_cds_ranges), strip(trailer_ranges))
  grl <- split(exons, exons$transcript_id)
  grl <- ORFik::sortPerGroup(GenomicRanges::reduce(grl))
  grl <- grl[features$transcript_id]
  out <- ORFik::pmapFromTranscriptF(IRanges::IRanges(features$tx_start, features$tx_end),
                                    grl, removeEmpty = TRUE)
  names(out) <- features$feature_id
  out
}

# Build the architectures into the transcript sequences made by simGenome.
# Returns the new leader, CDS and trailer strings and the feature table.
build_architectures <- function(arch, leader_string, cds_string, trailer_string,
                                start_codons, stop_codons, tx_names) {
  L <- Biostrings::width(leader_string); Lc <- Biostrings::width(cds_string)
  Tl <- Biostrings::width(trailer_string)
  tx <- paste0(as.character(leader_string), as.character(cds_string),
               as.character(trailer_string))
  written <- write_feature_sequences(tx, arch, L, Lc, stop_codons)
  tx <- written$tx_seqs; arch <- written$arch
  if (nrow(arch)) {
    arch <- verify_features(tx, arch, stop_codons, cds_end = L + Lc)
    arch[, `:=`(origin = "requested")]
  }
  decoys <- label_decoy_orfs(tx, arch, tx_names, L, Lc,
                             unique(c(start_codons, arch$start_codon)), stop_codons)
  if (nrow(arch)) arch[, nested_in := NA_character_]
  cols <- c("feature_id", "transcript_id", "gene", "class", "status", "origin",
            "tx_start", "tx_end", "frame", "start_codon", "nested_in")
  features <- data.table::rbindlist(list(
    if (nrow(arch)) arch[, ..cols], if (nrow(decoys)) decoys[, ..cols]), fill = TRUE)
  list(leader_string = Biostrings::DNAStringSet(substr(tx, 1L, L)),
       cds_string = Biostrings::DNAStringSet(substr(tx, L + 1L, L + Lc)),
       trailer_string = Biostrings::DNAStringSet(substr(tx, L + Lc + 1L, L + Lc + Tl)),
       features = features, rounds = written$rounds)
}

# Save the features (genomic ranges + table); translated uORF/uoORF also as
# true_uORFs.rds, the format the read simulation reads (named by transcript).
export_architecture_features <- function(features, leader_ranges, exon_cds_ranges,
                                         trailer_ranges, tx_names, out_dir, genome_path) {
  gr <- features_to_genome(features, leader_ranges, exon_cds_ranges, trailer_ranges,
                           tx_names)
  out <- c(features = file.path(out_dir, "true_features.rds"),
           features_table = sub("\\.fasta$", "_features.tsv", genome_path))
  saveRDS(gr, out[["features"]])
  # Every translated feature, of any class, gets reads through the extra-ORF
  # ("uorf") region of the read simulation, which samples each ORF in its own
  # frame. Its reads are labelled "uorf:<k>" (k = position in true_uORFs.rds);
  # read_group links that label back to the feature (see annotate_ground_truth()).
  translated <- features$status == "translated" & features$origin == "requested"
  features[, read_group := NA_character_]
  if (any(translated)) {
    features[translated, read_group := paste0("uorf:", seq_len(sum(translated)))]
    uorf_ranges <- gr[translated]
    names(uorf_ranges) <- features$transcript_id[translated]
    out <- c(out, uorfs = file.path(out_dir, "true_uORFs.rds"))
    saveRDS(uorf_ranges, out[["uorfs"]])
  }
  data.table::fwrite(features, out[["features_table"]], sep = "\t")
  out
}

# Per-gene build log for an architecture run (same uORF columns as the legacy
# log, plus all feature classes). Every requested feature was verified, or the
# build stopped before this point.
architecture_build_log <- function(features, tx_names) {
  req <- if (nrow(features)) features[origin == "requested"] else features
  log <- data.table::data.table(transcript_id = tx_names)
  if (nrow(req)) {
    counts <- req[, .(requested_uorfs = sum(class %in% c("uORF", "uoORF")),
                      requested_features = .N,
                      classes = paste(paste0(class, ifelse(status == "decoy", "(decoy)", "")),
                                      collapse = ",")),
                  by = transcript_id]
    log <- merge(log, counts, by = "transcript_id", all.x = TRUE, sort = FALSE)
  } else log[, `:=`(requested_uorfs = 0L, requested_features = 0L, classes = "")]
  log[is.na(requested_features), `:=`(requested_uorfs = 0L, requested_features = 0L, classes = "")]
  log[, `:=`(built_uorfs = requested_uorfs, built_features = requested_features,
             status = "ok", reason = NA_character_)]
  log[match(tx_names, transcript_id)]
}

# The features table as written by simGenome(architecture = ), with empty
# fields (written for NA) read back as NA.
read_features_table <- function(simGenome) {
  if (is.na(simGenome["features_table"]))
    stop("simGenome has no features table: build it with simGenome(architecture = )")
  x <- data.table::fread(simGenome[["features_table"]], na.strings = c("", "NA"))
  for (col in c("read_group", "nested_in"))
    if (col %in% names(x)) data.table::set(x, j = col, value = as.character(x[[col]]))
  x[]
}

#' Add the feature behind every simulated read to a ground-truth table
#'
#' Reads of the extra translons of a genome built with
#' \code{simGenome(architecture = )} are labelled \code{"uorf:<k>"} in the ground
#' truth. This adds \code{feature_id}, \code{class} (uORF, uoORF, NTE, intORF,
#' doORF, dORF), \code{feature_frame} (its frame relative to the CDS) and the
#' feature's transcript coordinates \code{feature_start}/\code{feature_end} for those reads, and
#' \code{class = "CDS"/"leader"/"trailer"} for the others.
#' @param truth the ground-truth table (data.table or a path to it, as written
#'  by \code{simNGScoverage(ground_truth = TRUE)})
#' @param simGenome the output of \code{simGenome()}
#' @return the truth table with the added columns
#' @export
annotate_ground_truth <- function(truth, simGenome) {
  if (is.character(truth)) truth <- data.table::fread(truth)
  truth <- data.table::as.data.table(truth)
  if (is.na(simGenome["features_table"]))
    stop("simGenome has no features table (it was not built with architecture =)")
  features <- read_features_table(simGenome)[!is.na(read_group)]
  region <- sub(":.*", "", truth$sampling_group)
  truth[, `:=`(feature_id = NA_character_, class = data.table::fifelse(
    region == "cds", "CDS", region), feature_frame = NA_integer_,
    feature_start = NA_integer_, feature_end = NA_integer_)]
  truth[features, `:=`(feature_id = i.feature_id, class = i.class, feature_frame = i.frame,
                       feature_start = i.tx_start, feature_end = i.tx_end),
        on = .(sampling_group = read_group)]
  # Transcript coordinate of each read's ribosome site, and its 0-based position
  # inside its feature (so feature_position %% 3 is the frame within the feature).
  exons <- ORFik::loadRegion(ORFik::loadTxdb(simGenome[["txdb"]]), "mrna")
  exons <- exons[truth$transcript_id]
  chrom <- as.character(GenomicRanges::seqnames(unlist(exons, use.names = FALSE)))[
    cumsum(c(1L, head(lengths(exons), -1L)))]
  site <- GenomicRanges::GRanges(chrom, IRanges::IRanges(truth$signal_position, width = 1L),
                                 strand = truth$strand)
  tx <- ORFik::pmapToTranscriptF(site, exons)
  truth[, tx_position := as.integer(BiocGenerics::start(tx))]
  truth[, feature_position := data.table::fifelse(
    !is.na(feature_start) & tx_position >= feature_start & tx_position <= feature_end,
    tx_position - feature_start, NA_integer_)]
  truth[]
}
