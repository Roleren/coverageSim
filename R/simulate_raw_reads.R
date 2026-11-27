# Raw read simulation: turn footprints (inserts) into the FASTQ a sequencer would
# give for a chosen library protocol, with a per-read truth table.
#
# A read is laid out 5' to 3' as
#
#   [5' UMI N_k][5' barcode][untemplated 5' bases][insert][untemplated 3' base]
#   [poly(A) tail][3' UMI N_m][3' barcode][linker][3' adapter ... flow cell][filler]
#
# and cut to the read length. Everything in this file is new and self-contained:
# it only reads the shipped end tables through method_end_bias().

RAW_READ_TRUSEQ_R1 <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC[index]ATCTCGTATGCCGTCTTCTGCTTG"
RAW_READ_SMALLRNA_R1 <- "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC[index]ATCTCGTATGCCGTCTTCTGCTTG"
RAW_READ_TRUSEQ_R2 <- "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
RAW_READ_SMALLRNA_R2 <- "GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"

raw_read_presets <- function() {
  list(
    ingolia2012 = list(
      method = "circularization", linker = "CTGTAGGCACCATCAAT",
      adapter3 = RAW_READ_TRUSEQ_R1, adapter5 = RAW_READ_TRUSEQ_R2,
      source = "Ingolia et al. 2012 Nat Protoc: 3' linker, circularization, no UMI"),
    mcglincy2017 = list(
      method = "circularization", umi5 = 2L, umi3 = 5L, barcode3 = "ATCGT",
      adapter3 = RAW_READ_TRUSEQ_R1, adapter5 = RAW_READ_TRUSEQ_R2,
      source = paste("McGlincy & Ingolia 2017 Methods: linker NNNNN + 5 nt sample barcode +",
                     "AGATCGGAAGAGCAC, RT primer adds NN at the read start")),
    nextflex_small_rna = list(
      method = "dual_ligation", umi5 = 4L, umi3 = 4L,
      adapter3 = RAW_READ_SMALLRNA_R1, adapter5 = RAW_READ_SMALLRNA_R2,
      source = "NEXTflex Small RNA-seq v3: 4 random nt on each adapter"),
    truseq_small_rna = list(
      method = "dual_ligation",
      adapter3 = RAW_READ_SMALLRNA_R1, adapter5 = RAW_READ_SMALLRNA_R2,
      source = "Illumina TruSeq Small RNA: 3' adapter TGGAATTCTCGG"),
    illumina_truseq = list(
      method = "dual_ligation",
      adapter3 = RAW_READ_TRUSEQ_R1, adapter5 = RAW_READ_TRUSEQ_R2,
      source = "Illumina TruSeq: adapter AGATCGGAAGAG"),
    smarter = list(
      method = "template_switching", add5 = list(base = "G", length_prob = c("3" = 1)),
      polya = c(10L, 30L),
      adapter3 = RAW_READ_TRUSEQ_R1, adapter5 = RAW_READ_TRUSEQ_R2,
      source = paste("SMARTer smRNA-seq (Takara): poly(A) tailing and template switching;",
                     "first 3 nt of read 1 come from the switch")),
    ottr = list(
      method = "template_switching", umi5 = 7L,
      add3 = list(bases = c(A = 0.5, G = 0.5), rate = 1),
      adapter3 = substring(RAW_READ_TRUSEQ_R1, 2), adapter5 = RAW_READ_TRUSEQ_R2,
      source = paste("OTTR (Ferguson et al. 2023): 7 nt UMI at the read start, one untemplated",
                     "3' purine before GATCGGAAGAGC; untemplated 5' bases off by default (set add5)")),
    polya_cats = list(
      method = "template_switching", add5 = list(base = "G", length_prob = c("3" = 1)),
      polya = c(10L, 40L),
      adapter3 = RAW_READ_TRUSEQ_R1, adapter5 = RAW_READ_TRUSEQ_R2,
      source = "CATS-like: poly(A) tailing, template switching (3 nt), TruSeq adapter"),
    ribolace = list(
      method = "ribolace", linker = "TCTCCTTGCATAATCACCAACC",
      adapter3 = RAW_READ_TRUSEQ_R1, adapter5 = RAW_READ_TRUSEQ_R2,
      source = paste("RiboLace / LACE-seq (Immagina): linker TCTCCTTGCATAATCACCAACC as in the",
                     "massiveNGSpipe adapter list; UMIs of the kit are not modelled (set umi5/umi3)")),
    custom = list(method = NULL, adapter3 = RAW_READ_TRUSEQ_R1, adapter5 = RAW_READ_TRUSEQ_R2,
                  source = "user defined")
  )
}

#' Library layouts for raw read simulation
#'
#' Describes what surrounds the footprint in a raw read: UMIs (random flanks),
#' fixed barcodes, untemplated bases, poly(A) tail, linker and adapter. Start
#' from a preset and override any field, or use \code{"custom"}.
#'
#' Read layout, 5' to 3': \code{umi5} random nt, \code{barcode5},
#' \code{add5}, the insert, \code{add3}, poly(A) tail, \code{umi3} random nt,
#' \code{barcode3}, \code{linker}, \code{adapter3}.
#'
#' Presets: \code{ingolia2012} (linker CTGTAGGCACCATCAAT), \code{mcglincy2017}
#' (5' NN, 3' NNNNN + 5 nt sample barcode, then AGATCGGAAGAGCAC),
#' \code{nextflex_small_rna} (4N / 4N, TGGAATTCTCGG), \code{truseq_small_rna}
#' (TGGAATTCTCGG), \code{illumina_truseq} (AGATCGGAAGAG), \code{smarter}
#' (template switching: 3 untemplated 5' G, poly(A) tail), \code{ottr} (7 nt 5'
#' UMI, one untemplated 3' purine), \code{polya_cats} (poly(A) tailing, CATS
#' like), \code{ribolace} (linker TCTCCTTGCATAATCACCAACC) and \code{custom}
#' (nothing but a TruSeq adapter). The circularization presets get the
#' untemplated 5' base of that chemistry from the shipped table (see
#' \code{\link{method_end_bias}}; T on about 10 percent of reads) unless
#' \code{add5} is given. \code{raw_read_protocols()} lists the presets.
#'
#' @param preset name of a preset, see Details.
#' @param name label of the layout, default the preset name.
#' @param method library chemistry for the terminal end tables of
#'   \code{\link{method_end_bias}}: "circularization", "dual_ligation",
#'   "template_switching", "ribolace" or NULL (none).
#' @param umi5,umi3 number of random nucleotides before / after the insert.
#' @param barcode5,barcode3 fixed barcode(s) before / after the insert: a
#'   character vector of equally long sequences (reads draw one each, a pooled
#'   lane), optionally named with the sampling probabilities as values of
#'   \code{barcode5_prob} / \code{barcode3_prob}. "" for none.
#' @param barcode5_prob,barcode3_prob sampling weights of the barcodes, default
#'   equal.
#' @param add5 untemplated 5' bases: \code{list(base = "G", length_prob =
#'   c("0" = 0.2, "1" = 0.5, "2" = 0.3))} (names are lengths), or NULL.
#' @param add3 one untemplated 3' base: \code{list(bases = c(A = 0.5, G = 0.5),
#'   rate = 1)}, or NULL.
#' @param polya poly(A) tail length range \code{c(min, max)} (uniform), or NULL.
#' @param linker sequence ligated to the insert 3' end ahead of the sequencing
#'   adapter ("" for none). Read 2 of a pair starts in it.
#' @param adapter3 sequence read 1 runs into after the linker; "[index]" is
#'   replaced by \code{index}.
#' @param adapter5 sequence read 2 runs into after the insert (already reverse
#'   complemented, as seen in read 2).
#' @param index the i7 index written into \code{adapter3}.
#' @return a list of class \code{rawReadProtocol}.
#' @export
#' @examples
#' raw_read_protocol("mcglincy2017")
#' raw_read_protocol("custom", umi5 = 6, barcode5 = c("ACGT", "TGCA"),
#'                   adapter3 = "CTGTCTCTTATACACATCT")
#' raw_read_protocols()
raw_read_protocol <- function(preset = "custom", name = preset, method, umi5, umi3,
                              barcode5, barcode3, barcode5_prob = NULL,
                              barcode3_prob = NULL, add5, add3, polya, linker,
                              adapter3, adapter5, index = "ATCACG") {
  presets <- raw_read_presets()
  if (length(preset) != 1L || !preset %in% names(presets)) {
    stop("preset must be one of: ", paste(names(presets), collapse = ", "))
  }
  protocol <- utils::modifyList(
    list(method = NULL, umi5 = 0L, umi3 = 0L, barcode5 = "", barcode3 = "",
         add5 = NULL, add3 = NULL, polya = NULL, linker = ""),
    presets[[preset]])
  given_add5 <- !missing(add5)
  if (!missing(method)) protocol["method"] <- list(method)
  if (!missing(umi5)) protocol$umi5 <- umi5
  if (!missing(umi3)) protocol$umi3 <- umi3
  if (!missing(barcode5)) protocol$barcode5 <- barcode5
  if (!missing(barcode3)) protocol$barcode3 <- barcode3
  if (given_add5) protocol["add5"] <- list(add5)
  if (!missing(add3)) protocol["add3"] <- list(add3)
  if (!missing(polya)) protocol["polya"] <- list(polya)
  if (!missing(linker)) protocol$linker <- linker
  if (!missing(adapter3)) protocol$adapter3 <- adapter3
  if (!missing(adapter5)) protocol$adapter5 <- adapter5
  if (!given_add5 && is.null(protocol$add5) &&
      (identical(protocol$method, "circularization") || identical(protocol$method, "ribolace"))) {
    added <- method_end_bias(protocol$method, nuclease = NULL, added_base = TRUE)$added_base
    if (added$rate > 0) {
      protocol$add5 <- list(base = added$base,
                            length_prob = c("0" = 1 - added$rate, "1" = added$rate))
    }
  }
  protocol$name <- name
  protocol$preset <- preset
  protocol$index <- index
  protocol$barcode5_prob <- barcode5_prob
  protocol$barcode3_prob <- barcode3_prob
  protocol$adapter3 <- sub("[index]", index, protocol$adapter3, fixed = TRUE)
  protocol$umi5 <- as.integer(protocol$umi5)
  protocol$umi3 <- as.integer(protocol$umi3)
  validate_raw_read_protocol(protocol)
  structure(protocol, class = "rawReadProtocol")
}

validate_raw_read_protocol <- function(protocol) {
  dna <- function(x) is.character(x) && !anyNA(x) && all(grepl("^[ACGTN]*$", x))
  for (field in c("umi5", "umi3")) {
    value <- protocol[[field]]
    if (length(value) != 1L || is.na(value) || value < 0L) {
      stop(field, " must be one non-negative whole number")
    }
  }
  for (field in c("barcode5", "barcode3")) {
    value <- protocol[[field]]
    if (!dna(value) || length(value) < 1L || length(unique(nchar(value))) != 1L) {
      stop(field, " must be DNA sequences of one common length")
    }
    prob <- protocol[[paste0(field, "_prob")]]
    if (!is.null(prob) && (length(prob) != length(value) || any(prob < 0) || sum(prob) <= 0)) {
      stop(field, "_prob must give one non-negative weight per barcode")
    }
  }
  for (field in c("linker", "adapter3", "adapter5")) {
    if (!dna(protocol[[field]]) || length(protocol[[field]]) != 1L) {
      stop(field, " must be one DNA sequence")
    }
  }
  add5 <- protocol$add5
  if (!is.null(add5)) {
    lengths <- suppressWarnings(as.integer(names(add5$length_prob)))
    if (!dna(add5$base) || nchar(add5$base) != 1L || length(lengths) == 0L ||
        anyNA(lengths) || any(lengths < 0L) || any(add5$length_prob < 0) ||
        sum(add5$length_prob) <= 0) {
      stop("add5 must be list(base = one nucleotide, length_prob = weights named by length)")
    }
  }
  add3 <- protocol$add3
  if (!is.null(add3)) {
    if (!dna(names(add3$bases)) || any(nchar(names(add3$bases)) != 1L) ||
        length(add3$rate) != 1L || add3$rate < 0 || add3$rate > 1) {
      stop("add3 must be list(bases = weights named by nucleotide, rate = share of reads)")
    }
  }
  polya <- protocol$polya
  if (!is.null(polya) && (length(polya) != 2L || anyNA(polya) || polya[1] < 0 || polya[2] < polya[1])) {
    stop("polya must be c(min, max) tail lengths")
  }
  invisible(TRUE)
}

#' @export
print.rawReadProtocol <- function(x, ...) {
  part <- function(label, value) if (nzchar(value)) paste0("[", label, " ", value, "]") else ""
  add5 <- if (is.null(x$add5)) "" else paste0(x$add5$base, "x", paste(names(x$add5$length_prob), collapse = "/"))
  add3 <- if (is.null(x$add3)) "" else paste(names(x$add3$bases), collapse = "/")
  polya <- if (is.null(x$polya)) "" else paste0(x$polya[1], "-", x$polya[2])
  cat("Raw read layout '", x$name, "' (", if (is.null(x$method)) "no chemistry table" else x$method, ")\n  ",
      part("UMI", if (x$umi5 > 0L) paste0("N", x$umi5) else ""),
      part("barcode", paste(x$barcode5, collapse = "|")), part("+5'", add5),
      "[insert]", part("+3'", add3), part("polyA", polya),
      part("UMI", if (x$umi3 > 0L) paste0("N", x$umi3) else ""),
      part("barcode", paste(x$barcode3, collapse = "|")), part("linker", x$linker),
      "[adapter ", substr(x$adapter3, 1, 20), "...]\n", sep = "")
  invisible(x)
}

#' @rdname raw_read_protocol
#' @export
raw_read_protocols <- function() {
  presets <- raw_read_presets()
  data.table::rbindlist(lapply(names(presets), function(preset) {
    protocol <- raw_read_protocol(preset)
    trim <- raw_read_trim_settings(protocol)
    data.table::data.table(
      preset = preset, method = if (is.null(protocol$method)) NA_character_ else protocol$method,
      umi5 = protocol$umi5, barcode5 = nchar(protocol$barcode5[1]),
      umi3 = protocol$umi3, barcode3 = nchar(protocol$barcode3[1]),
      polya = !is.null(protocol$polya), adapter = trim$adapter,
      trim_front = trim$trim_front, trim_tail = trim$trim_tail,
      source = presets[[preset]]$source)
  }))
}

#' Trimming settings that undo a raw read layout
#'
#' The adapter to remove and the number of bases to cut from each end after
#' adapter removal, for a \code{\link{raw_read_protocol}}. Untemplated 5' bases
#' are only counted when every read has the same number of them; a variable
#' number cannot be removed by a fixed cut and stays in the truth table.
#'
#' @param protocol a \code{\link{raw_read_protocol}} or a preset name.
#' @param adapter_length number of adapter bases to report, default 20.
#' @return a list: \code{adapter} (first bases of linker + adapter, or A's for a
#'   poly(A) layout), \code{trim_front}, \code{trim_tail} (cut after adapter
#'   removal), \code{polya} (TRUE when the adapter is a poly(A) tail) and
#'   \code{fastp}, two argument strings: pass 1 removes the adapter and the 5'
#'   bases, pass 2 the 3' bases.
#' @export
#' @examples
#' raw_read_trim_settings("mcglincy2017")
raw_read_trim_settings <- function(protocol, adapter_length = 20L) {
  if (is.character(protocol)) protocol <- raw_read_protocol(protocol)
  fixed_add5 <- 0L
  if (!is.null(protocol$add5)) {
    used <- protocol$add5$length_prob[protocol$add5$length_prob > 0]
    if (length(used) == 1L) fixed_add5 <- as.integer(names(used))
  }
  fixed_add3 <- if (!is.null(protocol$add3) && protocol$add3$rate == 1) 1L else 0L
  polya <- !is.null(protocol$polya) && protocol$polya[1] > 0
  adapter <- if (polya) strrep("A", max(6L, min(10L, protocol$polya[1]))) else {
    substr(paste0(protocol$linker, protocol$adapter3), 1L, adapter_length)
  }
  front <- protocol$umi5 + nchar(protocol$barcode5[1]) + fixed_add5
  tail <- protocol$umi3 + nchar(protocol$barcode3[1]) + fixed_add3
  common <- "--disable_quality_filtering --disable_length_filtering --disable_trim_poly_g"
  list(adapter = adapter, trim_front = front, trim_tail = tail, polya = polya,
       fastp = c(pass1 = paste("--adapter_sequence", adapter, "--trim_front1", front, common),
                 pass2 = paste("--disable_adapter_trimming --trim_tail1", tail, common)))
}

# ---------------------------------------------------------------------------
# Sequence sources and fragment drawing

raw_read_source <- function(sequences, encoded = FALSE) {
  if (is.character(sequences) && length(sequences) == 1L && file.exists(sequences)) {
    sequences <- Biostrings::readDNAStringSet(sequences)
    names(sequences) <- sub("\\s.*$", "", names(sequences))
  } else if (is.character(sequences)) {
    sequences <- Biostrings::DNAStringSet(sequences)
  }
  if (!methods::is(sequences, "DNAStringSet") || length(sequences) == 0L) {
    stop("sequences must be a DNAStringSet, a named character vector or a FASTA path")
  }
  if (is.null(names(sequences))) names(sequences) <- paste0("seq", seq_along(sequences))
  width <- Biostrings::width(sequences)
  big <- unlist(sequences, use.names = FALSE)
  # Decode just the portion of a footprint needed by each FASTQ record.
  if (encoded) big <- as.raw(big)
  list(big = big, names = names(sequences),
       width = width, offset = cumsum(c(0, width[-length(width)])))
}

raw_read_default_lengths <- function(enzyme) {
  shape <- switch(enzyme,
    "RNase I" = c(mean = 29.5, sd = 1.6, min = 24, max = 35),
    "RNase A" = c(mean = 28, sd = 2.5, min = 21, max = 36),
    "RNase T1" = c(mean = 31, sd = 3.5, min = 22, max = 42),
    "RNase T1 + S7" = c(mean = 31, sd = 3, min = 22, max = 40),
    c(mean = 32, sd = 3.5, min = 22, max = 42))
  lengths <- shape[["min"]]:shape[["max"]]
  prob <- stats::dnorm(lengths, shape[["mean"]], shape[["sd"]])
  stats::setNames(prob / sum(prob), lengths)
}

raw_read_length_table <- function(insert_lengths, enzyme) {
  if (is.null(insert_lengths)) insert_lengths <- raw_read_default_lengths(enzyme)
  if (is.null(names(insert_lengths))) {
    insert_lengths <- stats::setNames(rep(1, length(insert_lengths)), insert_lengths)
  }
  lengths <- suppressWarnings(as.integer(names(insert_lengths)))
  if (anyNA(lengths) || any(lengths < 1L) || any(insert_lengths < 0) || sum(insert_lengths) <= 0) {
    stop("insert_lengths must be lengths, or weights named by length")
  }
  list(length = lengths, prob = as.numeric(insert_lengths) / sum(insert_lengths))
}

# Draw n footprints from a source: sequences weighted by length times a
# log-normal abundance; `hotspots` restricts them to that many distinct
# fragments with Zipf weights (rRNA / tRNA like).
raw_fragment_sampler <- function(source, lengths, expression_sd, pad = 3L) {
  usable <- which(source$width >= max(lengths$length) + 2L * pad + 1L)
  if (length(usable) == 0L) {
    usable <- which(source$width >= min(lengths$length) + 2L * pad + 1L)
    if (length(usable) == 0L) stop("no sequence is long enough for the insert lengths")
  }
  weight <- source$width[usable] * stats::rlnorm(length(usable), 0, expression_sd)
  cdf <- function(weight) {
    cumulative <- cumsum(weight)
    if (!is.finite(tail(cumulative, 1L)) || tail(cumulative, 1L) <= 0) {
      stop("footprint sampling weights must have a finite positive sum")
    }
    cumulative / tail(cumulative, 1L)
  }
  source_cdf <- cdf(weight)
  length_cdf <- cdf(lengths$prob)
  source_width <- source$width
  function(n) {
    raw_draw_fragments_cpp(n, usable, source_width, lengths$length, source_cdf, length_cdf, pad)
  }
}

raw_read_weight_vector <- function(profile) {
  weights <- rep(1, 16L)
  if (is.null(profile)) return(weights)
  table <- profile$table
  strength <- if (is.null(profile$strength)) 1 else profile$strength
  weights[c(1L, 2L, 4L, 8L)] <- table$weight[match(c("A", "C", "G", "T"), table$kmer)]^strength
  weights[is.na(weights)] <- 1
  weights
}

# Nuclease and chemistry preferences as a small shift of each read end: every
# end may move up to max_shift nt, with probability proportional to the weight
# of the nucleotide just outside the read (cut site) times the weight of the
# terminal nucleotide inside it.
raw_end_bias_sampler <- function(source, bias, max_shift, fixed_shifts = FALSE, min_length = 12L) {
  max_shift <- min(max_shift, max(source$width) - 1L)
  if (max_shift < 1L) return(identity)
  inner5 <- raw_read_weight_vector(bias$five_prime_bias)
  outer5 <- raw_read_weight_vector(bias$five_prime_flank_bias)
  inner3 <- raw_read_weight_vector(bias$three_prime_bias)
  outer3 <- raw_read_weight_vector(bias$three_prime_flank_bias)
  key <- if (fixed_shifts) sample.int(.Machine$integer.max, 2L, replace = TRUE)
  function(fragments) {
    shifted <- raw_shift_ends_cpp(source$big, source$offset, source$width, fragments$seq_idx,
                                  fragments$start, fragments$end, inner5, outer5,
                                  inner3, outer3, max_shift, min_length, fragments$copy_rank, key)
    fragments$start <- shifted$start
    fragments$end <- shifted$end
    fragments
  }
}

raw_source_inserts <- function(source, fragments) {
  if (nrow(fragments) == 0L) return(character())
  global <- source$offset[fragments$seq_idx]
  as.character(Biostrings::extractAt(
    source$big, IRanges::IRanges(global + fragments$start, global + fragments$end)))
}

raw_reverse_complement <- function(x) {
  if (length(x) == 0L) return(character())
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}

# Footprints given as a table on the sequences themselves.
raw_table_fragments <- function(fragments, source) {
  fragments <- data.table::as.data.table(fragments)
  missing_columns <- setdiff(c("seqname", "start", "end"), names(fragments))
  if (length(missing_columns) > 0L) {
    stop("fragments needs the columns seqname, start, end (and optionally count)")
  }
  seq_idx <- match(as.character(fragments$seqname), source$names)
  if (anyNA(seq_idx)) stop("fragments has seqnames that are not in sequences")
  start <- as.integer(fragments$start)
  end <- as.integer(fragments$end)
  if (anyNA(start) || anyNA(end) || any(start < 1L) || any(end < start) ||
      any(end > source$width[seq_idx])) {
    stop("fragments must lie inside their sequence with start <= end")
  }
  rows <- data.table::data.table(seq_idx = seq_idx, start = start, end = end)
  list(n = nrow(rows), take = function(index) as.list(rows[index]), count = fragments$count)
}

# Footprints given as genomic ranges or alignments (for example a coverageSim
# library read with ORFik::fimport()): the insert is the sequence of the range
# on its strand, spliced reads are joined over their exons.
raw_genomic_source <- function(fragments, sequences) {
  seqname <- as.character(unique(GenomicRanges::seqnames(fragments)))
  if (is.character(sequences) && length(sequences) == 1L && file.exists(sequences)) {
    sequences <- Biostrings::readDNAStringSet(sequences)
    names(sequences) <- sub("\\s.*$", "", names(sequences))
  }
  if (!all(seqname %in% names(sequences))) {
    stop("fragments has seqnames that are not in sequences")
  }
  raw_read_source(sequences[names(sequences) %in% seqname])
}

raw_genomic_fragments <- function(fragments, source) {
  spliced <- methods::is(fragments, "GAlignments")
  blocks <- if (spliced) GenomicAlignments::grglist(fragments) else NULL
  ranges <- if (spliced) GenomicRanges::granges(fragments) else fragments
  seqname <- as.character(GenomicRanges::seqnames(ranges))
  piece <- function(x) {
    raw_source_inserts(source, data.table::data.table(
      seq_idx = match(as.character(GenomicRanges::seqnames(x)), source$names),
      start = GenomicRanges::start(x), end = GenomicRanges::end(x)))
  }
  if (spliced) {
    flat <- unlist(blocks, use.names = FALSE)
    exon <- data.table::data.table(read = rep(seq_along(blocks), lengths(blocks)),
                                   start = GenomicRanges::start(flat), sequence = piece(flat))
    data.table::setorder(exon, read, start)
    insert <- exon[, list(sequence = paste(sequence, collapse = "")), by = read]$sequence
  } else {
    insert <- piece(ranges)
  }
  strand <- as.character(GenomicRanges::strand(ranges))
  minus <- strand == "-"
  insert[minus] <- raw_reverse_complement(insert[minus])
  data.table::data.table(seqname = seqname, strand = strand,
                         insert_start = GenomicRanges::start(ranges),
                         insert_end = GenomicRanges::end(ranges), insert = insert)
}

raw_fragment_table <- function(source, fragments, label, reference) {
  list(source = label, seqname = source$names[fragments$seq_idx],
       strand = "+", insert_start = fragments$start, insert_end = fragments$end,
       insert = "", reference = reference,
       offset = source$offset[fragments$seq_idx] + fragments$start - 1)
}

# Index copies with cumulative counts, avoiding a row per copy in memory.
raw_supplied_sampler <- function(fragments, replace, keep_rank) {
  count <- fragments$count
  cumulative <- if (!is.null(count)) {
    count <- trunc(as.numeric(count))
    if (anyNA(count) || any(count < 0) || any(!is.finite(count))) {
      stop("fragment counts must be finite and non-negative")
    }
    cumsum(count)
  }
  total <- if (is.null(cumulative)) fragments$n else sum(tail(cumulative, 1L))
  if (total > .Machine$integer.max) stop("too many supplied footprint copies")
  total <- as.integer(total)
  key <- if (!replace && total > 1L) sample.int(.Machine$integer.max, 4L, replace = TRUE)
  list(total = total, draw = function(n, first) {
    if (total == 0L) stop("cannot sample from zero supplied footprints")
    rank <- if (replace) sample.int(total, n, replace = TRUE) else if (total == 1L) 1L else {
      raw_permuted_indices_cpp(first, n, total, key)
    }
    index <- if (is.null(cumulative)) rank else findInterval(rank - 1L, cumulative) + 1L
    rows <- fragments$take(index)
    if (keep_rank) rows$copy_rank <- rank
    rows
  })
}

# Hypergeometric chunk counts followed by a shuffle give a random order with
# fixed library-wide source totals and no library-wide permutation.
raw_chunk_counts <- function(remaining, n) {
  if (length(remaining) == 1L) return(n)
  if (sum(remaining) == n) return(remaining)
  counts <- integer(length(remaining))
  available <- sum(remaining)
  for (i in seq_len(length(remaining) - 1L)) {
    counts[i] <- stats::rhyper(1L, remaining[i], available - remaining[i], n)
    available <- available - remaining[i]
    n <- n - counts[i]
  }
  counts[length(counts)] <- n
  counts
}

# ---------------------------------------------------------------------------
# Read assembly, qualities, errors, writing

raw_random_kmers <- function(n, k) {
  if (k == 0L) return("")
  if (n == 0L) return(character())
  if (k > 6L) return(paste0(raw_random_kmers(n, 6L), raw_random_kmers(n, k - 6L)))
  kmers <- do.call(paste0, rev(expand.grid(rep(list(c("A", "C", "G", "T")), k),
                                           stringsAsFactors = FALSE)))
  kmers[sample.int(length(kmers), n, replace = TRUE)]
}

raw_pick_barcode <- function(n, barcodes, prob) {
  if (length(barcodes) == 1L) return(barcodes)
  barcodes[sample.int(length(barcodes), n, replace = TRUE, prob = prob)]
}

raw_added_five <- function(n, add5) {
  if (is.null(add5)) return("")
  lengths <- as.integer(names(add5$length_prob))
  bases <- strrep(add5$base, lengths)
  if (length(bases) == 1L) return(bases)
  bases[sample.int(length(bases), n, replace = TRUE, prob = add5$length_prob)]
}

raw_added_three <- function(n, add3) {
  if (is.null(add3)) return("")
  bases <- names(add3$bases)
  drawn <- bases[sample.int(length(bases), n, replace = TRUE, prob = add3$bases)]
  drawn[stats::runif(n) >= add3$rate] <- ""
  drawn
}

raw_polya_lengths <- function(n, polya) {
  if (is.null(polya)) return(0L)
  as.integer(polya[1] + floor(stats::runif(n) * (polya[2] - polya[1] + 1)))
}

#' Simulate raw sequencing reads (FASTQ) from footprints
#'
#' Wraps footprints in everything a library protocol and a sequencer add, and
#' writes the reads as FASTQ with a truth table: random flanks (UMIs), fixed
#' barcodes, untemplated bases, poly(A) tail, linker and adapter, cut to the
#' read length, with Illumina-like qualities and substitution errors, adapter
#' dimers, contaminants (rRNA, tRNA) and reads of a second organism. Made for
#' testing trimming and layout detection: every read has a known answer.
#'
#' Footprints come from one of three inputs:
#' \itemize{
#'   \item \code{fragments = NULL}: \code{n_reads} footprints are drawn from
#'     \code{sequences} (transcripts or any other sequences; forward strand).
#'     Sequences are used in proportion to their length times a log-normal
#'     abundance, lengths follow \code{insert_lengths} and both ends are moved
#'     by up to \code{end_shift} nt according to the nuclease cut-site and
#'     chemistry end tables of \code{\link{method_end_bias}}.
#'   \item a table with \code{seqname}, \code{start}, \code{end} (1-based, on
#'     \code{sequences}) and optionally \code{count}. Used as given;
#'     \code{end_bias = TRUE} also moves these ends.
#'   \item a \code{GRanges} or \code{GAlignments} of reads on a genome, for
#'     example a library written by \code{\link{simNGScoverage}} and read with
#'     \code{ORFik::fimport()}, with \code{sequences} the genome. Used as given
#'     (coverageSim already placed the ends); a \code{score} column counts
#'     copies and spliced reads are joined over their exons.
#' }
#' With \code{fragments} given, \code{n_reads = NULL} uses each footprint once
#' and \code{n_reads} samples that many with replacement.
#'
#' Read 1 is \code{[umi5][barcode5][add5][insert][add3][poly(A)][umi3][barcode3][linker][adapter3]}
#' followed by filler (G on two-colour instruments, A otherwise), cut to
#' \code{read_length}. Read 2 of a pair is the reverse complement of everything
#' up to and including the linker, followed by \code{adapter5}. Adapter dimers
#' have an empty insert but keep the UMIs and barcodes.
#'
#' Files written: \code{<out_prefix>.fastq.gz} (or \code{_1} / \code{_2}),
#' \code{<out_prefix>_truth.csv.gz} with one row per read (\code{read_id},
#' \code{source}, \code{seqname}, \code{strand}, \code{insert_start},
#' \code{insert_end}, \code{insert_len}, \code{umi5}, \code{barcode5},
#' \code{add5}, \code{add3}, \code{polya_len}, \code{umi3}, \code{barcode3},
#' \code{insert_start_read}, \code{insert_end_read}, \code{adapter_start_read},
#' \code{insert_complete}, \code{n_errors}, \code{n_errors_insert} (wrong or N
#' bases in read 1 and in its insert part),
#' \code{adapter}, \code{enzyme}, \code{protocol}) and
#' \code{<out_prefix>_summary.csv} with one row for the library. Positions in
#' the read are 1-based and not cut to the read length:
#' \code{adapter_start_read} is the first base after the 3' barcode and
#' \code{insert_complete} says whether the whole insert is inside read 1.
#'
#' Footprints and read noise are generated in chunks. The same seed and chunk
#' size reproduce the simulation within this implementation. Seeded reads and
#' truth rows differ from earlier versions, which drew all footprints before
#' assembly and used R's per-base random number stream.
#'
#' @param sequences the sequences footprints come from: a \code{DNAStringSet},
#'   a named character vector or a FASTA path.
#' @param fragments NULL, a table or genomic ranges; see Details.
#' @param n_reads total number of reads to write, including dimers and
#'   contaminants. With \code{fragments} given and NULL, dimers and contaminants
#'   are added on top so that their shares hold.
#' @param protocol a \code{\link{raw_read_protocol}} or a preset name.
#' @param out_prefix output path without extension.
#' @param read_length length of the reads, for example 50, 75, 100 or 150.
#' @param paired TRUE to also write read 2.
#' @param enzyme nuclease of the cut-site table and the default insert lengths:
#'   "RNase I" (default), "RNase A", "RNase T1", "RNase S7", "RNase T1 + S7",
#'   "MNase", or NULL for no cut-site preference.
#' @param insert_lengths footprint lengths to draw: weights named by length, or
#'   NULL for a default per enzyme (RNase I around 29-30 nt; the others are
#'   wider, rough approximations).
#' @param end_bias apply the end tables when footprints are drawn (default
#'   TRUE) or given as a table (default FALSE there).
#' @param end_shift largest shift of a read end by the end tables, default 2.
#' @param dimer_rate share of reads that are adapter dimers (empty insert).
#' @param contaminants named list of other read sources, each
#'   \code{list(sequences = , fraction = , hotspots = NULL, insert_lengths =
#'   NULL)}. \code{hotspots = 200} draws all its reads from 200 distinct
#'   fragments with Zipf weights, as rRNA and tRNA fragments are; leave NULL
#'   for a second organism. The name is the \code{source} in the truth table.
#' @param organism \code{source} label of the reads from \code{sequences}.
#' @param quality "novaseq" (four binned values, G filler), "hiseq" (full
#'   range, A filler) or "perfect" (Q40, no errors).
#' @param error_rate substitution rate added on top of what the qualities imply.
#' @param expression_sd log-normal spread of the abundance per sequence.
#' @param library name used in the read ids, default the file name.
#' @param seed random seed, or NULL to leave the random state alone.
#' @param chunk_size reads built in memory at a time.
#' @param compress write gzip files. On Unix-like systems, the build uses ISA-L
#'   or libdeflate when detected by pkg-config, with zlib as the fallback.
#' @param threads number of compression workers. The decompressed FASTQ and
#'   truth files are identical across worker counts for the same seed and
#'   chunk size.
#' @return invisibly a list: \code{files}, \code{summary} (the one-row library
#'   table) and \code{trim} (\code{\link{raw_read_trim_settings}}).
#' @export
#' @examples
#' set.seed(1)
#' transcripts <- Biostrings::DNAStringSet(vapply(1:20, function(i)
#'   paste(sample(c("A", "C", "G", "T"), 600, replace = TRUE), collapse = ""), ""))
#' names(transcripts) <- paste0("tx", 1:20)
#' out <- simulate_raw_reads(transcripts, n_reads = 2000, protocol = "mcglincy2017",
#'                           out_prefix = file.path(tempdir(), "example"),
#'                           read_length = 75, dimer_rate = 0.05, seed = 1)
#' out$summary
#' out$trim
simulate_raw_reads <- function(sequences, fragments = NULL, n_reads = NULL,
                               protocol = "illumina_truseq", out_prefix,
                               read_length = 50L, paired = FALSE, enzyme = "RNase I",
                               insert_lengths = NULL, end_bias = is.null(fragments),
                               end_shift = 2L, dimer_rate = 0, contaminants = NULL,
                               organism = "main",
                               quality = c("novaseq", "hiseq", "perfect"),
                               error_rate = 0, expression_sd = 1.5,
                               library = basename(out_prefix), seed = NULL,
                               chunk_size = 250000L, compress = TRUE, threads = 2L) {
  quality <- match.arg(quality)
  if (is.character(protocol)) protocol <- raw_read_protocol(protocol)
  if (!inherits(protocol, "rawReadProtocol")) stop("protocol must be a raw_read_protocol()")
  read_length <- as.integer(read_length)
  if (length(read_length) != 1L || is.na(read_length) || read_length < 10L) {
    stop("read_length must be one whole number of at least 10")
  }
  if (!is.numeric(dimer_rate) || length(dimer_rate) != 1L || dimer_rate < 0 || dimer_rate >= 1) {
    stop("dimer_rate must be in [0, 1)")
  }
  if (!is.numeric(error_rate) || length(error_rate) != 1L || error_rate < 0 || error_rate > 1) {
    stop("error_rate must be in [0, 1]")
  }
  fractions <- vapply(contaminants, function(x) as.numeric(x$fraction), numeric(1))
  if (length(contaminants) > 0L && (is.null(names(contaminants)) || any(fractions < 0))) {
    stop("contaminants must be a named list with non-negative fractions")
  }
  other <- dimer_rate + sum(fractions)
  if (other >= 1) stop("dimer_rate and contaminant fractions must sum to less than 1")
  if (is.null(fragments) && is.null(n_reads)) stop("n_reads is needed when fragments is NULL")
  if (!is.null(n_reads) && (!is.numeric(n_reads) || length(n_reads) != 1L ||
      !is.finite(n_reads) || n_reads < 0 || n_reads > .Machine$integer.max)) {
    stop("n_reads must be one finite non-negative read count")
  }
  if (!is.null(seed)) set.seed(seed)
  enzyme_label <- if (is.null(enzyme)) NA_character_ else as.character(enzyme)
  protocol_label <- as.character(protocol$name)
  organism <- as.character(organism)
  if (length(organism) != 1L) stop("organism must be one source label")
  if (length(protocol_label) != 1L) stop("protocol must have one name")
  bias <- method_end_bias(protocol$method, nuclease = enzyme)
  lengths <- raw_read_length_table(insert_lengths, if (is.null(enzyme)) "RNase I" else enzyme)

  chunk_size <- as.integer(chunk_size)
  threads <- as.integer(threads)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("chunk_size must be one positive integer")
  }
  if (length(threads) != 1L || is.na(threads) || threads < 1L) {
    stop("threads must be one positive integer")
  }
  end_shift <- as.integer(end_shift)
  if (length(end_shift) != 1L || is.na(end_shift)) stop("end_shift must be one integer")
  make_end_bias <- function(source, fixed_shifts = FALSE) {
    raw_end_bias_sampler(source, bias, end_shift, fixed_shifts)
  }
  # Each source provides chunk-sized coordinates and materializes inserts on demand.
  genomic <- methods::is(fragments, "GRanges") || methods::is(fragments, "GAlignments")
  if (genomic) {
    source <- if (length(fragments)) raw_genomic_source(fragments, sequences)
    supplied <- list(n = length(fragments), count = S4Vectors::mcols(fragments)$score,
                     take = function(index) raw_genomic_fragments(fragments[index], source))
  } else {
    source <- raw_read_source(sequences, encoded = TRUE)
    supplied <- if (!is.null(fragments)) raw_table_fragments(fragments, source)
  }
  fixed_shifts <- !is.null(supplied) && !genomic && isTRUE(end_bias) && !is.null(n_reads)
  main <- if (!is.null(supplied)) raw_supplied_sampler(supplied, !is.null(n_reads), fixed_shifts)
  total <- if (is.null(n_reads)) round(main$total / (1 - other)) else as.integer(n_reads)
  n_main <- if (is.null(n_reads)) main$total else round(n_reads * (1 - other))
  take_main <- if (!is.null(main)) main$draw else if (n_main > 0L) {
    draw_main <- raw_fragment_sampler(source, lengths, expression_sd)
    function(n, first) draw_main(n)
  }
  samplers <- list(take_main)
  shift_ends <- list(if (!genomic && n_main > 0L && isTRUE(end_bias)) {
    make_end_bias(source, fixed_shifts)
  } else identity)
  materialize <- list(if (genomic) {
    function(rows) cbind(source = organism, rows, reference = 0L, offset = 0)
  } else {
    function(rows) raw_fragment_table(source, rows, organism, 1L)
  })

  n_dimer <- round(total * dimer_rate)
  labels <- c(organism, names(contaminants), "adapter_dimer")
  references <- vector("list", length(labels))
  if (!genomic) references[[1L]] <- source$big
  counts <- as.integer(c(n_main, round(total * fractions), n_dimer))
  n <- sum(counts)
  if (!n) stop("no reads to write")
  if (n > .Machine$integer.max) stop("too many reads to write")
  active <- which(counts > 0L)
  for (i in seq_along(contaminants) + 1L) {
    if (counts[i] == 0L) next
    entry <- contaminants[[i - 1L]]
    entry_source <- raw_read_source(entry$sequences, encoded = TRUE)
    references[[i]] <- entry_source$big
    entry_lengths <- if (is.null(entry$insert_lengths)) lengths else {
      raw_read_length_table(entry$insert_lengths, "RNase I")
    }
    sampler <- raw_fragment_sampler(entry_source, entry_lengths, expression_sd)
    if (is.null(entry$hotspots)) {
      samplers[[i]] <- local({
        draw <- sampler
        function(n, first) draw(n)
      })
      shift_ends[[i]] <- make_end_bias(entry_source)
    } else {
      samplers[[i]] <- local({
        hotspots <- as.integer(entry$hotspots)
        pool <- data.table::setDT(sampler(hotspots))
        prob <- 1 / seq_len(hotspots)
        function(n, first) as.list(pool[sample.int(hotspots, n, replace = TRUE, prob = prob)])
      })
      shift_ends[[i]] <- identity
    }
    materialize[[i]] <- local({
      src <- entry_source
      label <- labels[i]
      reference <- i
      function(rows) raw_fragment_table(src, rows, label, reference)
    })
  }
  dimer <- length(labels)
  samplers[[dimer]] <- function(n, first) n
  shift_ends[[dimer]] <- identity
  materialize[[dimer]] <- function(n) list(
    source = "adapter_dimer", seqname = NA_character_, strand = NA_character_,
    insert_start = rep.int(NA_integer_, n), insert_end = NA_integer_,
    insert = "", reference = 0L, offset = 0)

  # Output
  separator <- getOption("datatable.fwrite.sep", ",")
  if (length(separator) != 1L || is.na(separator) || nchar(separator, type = "bytes") != 1L) {
    stop("the truth table separator must be one single-byte character")
  }
  logical01 <- isTRUE(getOption("datatable.logical01", FALSE))

  if (length(library) != 1L) stop("library must be one name")
  id_prefix <- paste0(library, ".")
  id_width <- max(7L, nchar(n))

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  extension <- if (compress) ".fastq.gz" else ".fastq"
  fastq <- if (paired) paste0(out_prefix, c("_1", "_2"), extension) else paste0(out_prefix, extension)
  truth_file <- paste0(out_prefix, "_truth.csv", if (compress) ".gz" else "")
  summary_file <- paste0(out_prefix, "_summary.csv")
  if (file.exists(truth_file)) file.remove(truth_file)
  writers <- list()
  on.exit(for (writer in writers) try(raw_close_writer_cpp(writer), silent = TRUE), add = TRUE)
  for (path in c(fastq, truth_file)) {
    writers[[length(writers) + 1L]] <- raw_open_writer_cpp(path.expand(path), compress, threads)
  }
  truth_writer <- writers[[length(writers)]]
  profile <- match(quality, c("perfect", "novaseq", "hiseq")) - 1L
  model <- raw_read_model_cpp(read_length, profile, error_rate)
  buffer <- raw_read_buffer_cpp()
  on.exit(raw_release_buffer_cpp(buffer), add = TRUE)
  adapter <- substr(paste0(protocol$linker, protocol$adapter3), 1L, 20L)
  totals <- list(errors = 0, complete = 0, adapter_in_read = 0, insert_bases = 0)
  remaining <- counts[active]
  composition_order <- character()
  last <- 0L
  while (last < n) {
    first <- last + 1L
    size <- min(chunk_size, n - last)
    last <- last + size
    take <- raw_chunk_counts(remaining, size)
    part <- which(take > 0L)
    groups <- active[part]
    drawn <- lapply(part, function(i) samplers[[active[i]]](take[i], counts[active[i]] - remaining[i] + 1L))
    shifted <- Map(function(draw, i) shift_ends[[i]](draw), drawn, groups)
    parts <- Map(function(draw, i) materialize[[i]](draw), shifted, groups)
    reads <- if (length(parts) == 1L) parts[[1L]] else {
      data.table::rbindlist(lapply(parts, data.table::as.data.table), use.names = TRUE)[sample.int(size)]
    }
    if (length(composition_order) < length(active)) {
      composition_order <- union(composition_order, unique(reads$source))
    }
    remaining <- remaining - take

    umi5 <- raw_random_kmers(size, protocol$umi5)
    barcode5 <- raw_pick_barcode(size, protocol$barcode5, protocol$barcode5_prob)
    add5 <- raw_added_five(size, protocol$add5)
    add3 <- raw_added_three(size, protocol$add3)
    polya_len <- raw_polya_lengths(size, protocol$polya)
    umi3 <- raw_random_kmers(size, protocol$umi3)
    barcode3 <- raw_pick_barcode(size, protocol$barcode3, protocol$barcode3_prob)
    pieces <- list(umi5, barcode5, add5, add3, umi3, barcode3)
    inserts <- list(reads$insert, reads$reference, reads$offset, reads$insert_start, reads$insert_end)
    made <- raw_read_records_cpp(model, buffer, id_prefix, id_width, first, size, pieces, polya_len,
                                 protocol$linker, protocol$adapter3, references, inserts, FALSE)
    raw_write_fastq_cpp(writers[[1L]], buffer)
    if (paired) {
      raw_read_records_cpp(model, buffer, id_prefix, id_width, first, size, pieces, polya_len,
                           protocol$linker, protocol$adapter5, references, inserts, TRUE)
      raw_write_fastq_cpp(writers[[2L]], buffer)
    }
    truth <- list(
      source = reads$source, seqname = reads$seqname,
      strand = reads$strand, insert_start = reads$insert_start,
      insert_end = reads$insert_end, insert_len = made$insert_len,
      umi5 = umi5, barcode5 = barcode5, add5 = add5, add3 = add3,
      polya_len = polya_len, umi3 = umi3, barcode3 = barcode3,
      insert_start_read = made$insert_start_read, insert_end_read = made$insert_end_read,
      adapter_start_read = made$adapter_start_read,
      insert_complete = made$insert_end_read <= read_length,
      n_errors = made$n_errors, n_errors_insert = made$n_errors_insert,
      adapter = adapter,
      enzyme = enzyme_label, protocol = protocol_label)
    raw_write_truth_cpp(truth_writer, truth, id_prefix, id_width, first, size, separator, logical01)
    totals$errors <- totals$errors + made$errors
    totals$complete <- totals$complete + sum(truth$insert_complete)
    totals$adapter_in_read <- totals$adapter_in_read + sum(made$adapter_start_read <= read_length)
    totals$insert_bases <- totals$insert_bases + sum(made$insert_len)
  }
  for (writer in writers) raw_close_writer_cpp(writer)

  trim <- raw_read_trim_settings(protocol)
  composition <- data.table::data.table(source = labels[active], reads = counts[active])[
    , list(reads = sum(reads)), by = source][match(composition_order, source)]
  describe_add5 <- if (is.null(protocol$add5)) "" else {
    paste0(protocol$add5$base, ":", paste0(names(protocol$add5$length_prob), "=",
                                          signif(protocol$add5$length_prob / sum(protocol$add5$length_prob), 3),
                                          collapse = ";"))
  }
  describe_add3 <- if (is.null(protocol$add3)) "" else {
    paste0(paste(names(protocol$add3$bases), collapse = "/"), ":rate=", protocol$add3$rate)
  }
  summary <- data.table::data.table(
    library = library, protocol = protocol$name, preset = protocol$preset,
    method = if (is.null(protocol$method)) NA_character_ else protocol$method,
    enzyme = enzyme_label, organism = organism, n_reads = n, read_length = read_length,
    paired = paired, umi5 = protocol$umi5, barcode5 = paste(protocol$barcode5, collapse = ";"),
    add5 = describe_add5, add3 = describe_add3,
    polya = if (is.null(protocol$polya)) "" else paste(protocol$polya, collapse = "-"),
    umi3 = protocol$umi3, barcode3 = paste(protocol$barcode3, collapse = ";"),
    linker = protocol$linker, adapter3 = protocol$adapter3,
    trim_adapter = trim$adapter, trim_front = trim$trim_front, trim_tail = trim$trim_tail,
    dimer_rate = n_dimer / n,
    composition = paste0(composition$source, "=", signif(composition$reads / n, 4), collapse = ";"),
    mean_insert_len = totals$insert_bases / max(1, n - n_dimer),
    insert_complete = totals$complete / n, adapter_in_read = totals$adapter_in_read / n,
    error_rate = totals$errors / (as.numeric(n) * read_length), quality = quality,
    seed = if (is.null(seed)) NA_integer_ else as.integer(seed),
    fastq = paste(fastq, collapse = ";"), truth = truth_file)
  data.table::fwrite(summary, summary_file)
  invisible(list(files = c(fastq = fastq, truth = truth_file, summary = summary_file),
                 summary = summary, trim = trim))
}
