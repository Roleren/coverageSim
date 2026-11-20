normalize_technical_artifacts <- function(technical_artifacts) {
  defaults <- list(
    duplication_rate = 0,
    duplicate_copies = 1L,
    multimapping_rate = 0,
    secondary_alignments = 1L,
    spike_sites = 0L,
    spike_library_fraction = 0.25,
    spike_multiplicity = c(1000, 100000),
    spike_lengths = c(20L, 25L, 32L),
    spike_seed = 1L
  )
  if (is.null(technical_artifacts)) return(defaults)
  if (!is.list(technical_artifacts)) stop("technical_artifacts must be a list")
  unknown <- setdiff(names(technical_artifacts), names(defaults))
  if (length(unknown)) {
    stop("Unknown technical_artifacts field(s): ", paste(unknown, collapse = ", "))
  }
  result <- utils::modifyList(defaults, technical_artifacts)
  for (name in c("duplication_rate", "multimapping_rate")) {
    value <- result[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 0 || value > 1) {
      stop(name, " must be one finite probability from zero to one")
    }
  }
  for (name in c("duplicate_copies", "secondary_alignments")) {
    value <- result[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != as.integer(value)) {
      stop(name, " must be one positive integer")
    }
    result[[name]] <- as.integer(value)
  }
  validate_spike_settings(result)
}

# Oligo spikes (carina-claude, E17-S11): one read sequence at one position,
# repeated 10^3-10^5 times in a few libraries, with a length outside the
# footprint range. Measured in backclaude E16 (genes12: LSM4 20 nt, TEKTL1
# 25 nt) and E6 (PRJNA807729 H1-2, 32 nt).
validate_spike_settings <- function(result) {
  n <- result$spike_sites
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 ||
      n != as.integer(n)) {
    stop("spike_sites must be one non-negative integer")
  }
  result$spike_sites <- as.integer(n)
  f <- result$spike_library_fraction
  if (!is.numeric(f) || length(f) != 1L || !is.finite(f) || f <= 0 || f > 1) {
    stop("spike_library_fraction must be one number in (0, 1]")
  }
  m <- result$spike_multiplicity
  if (!is.numeric(m) || length(m) != 2L || any(!is.finite(m)) || m[1] < 1 ||
      m[2] < m[1]) {
    stop("spike_multiplicity must be two numbers c(min, max) with 1 <= min <= max")
  }
  l <- result$spike_lengths
  if (!is.numeric(l) || !length(l) || any(!is.finite(l)) || any(l < 1) ||
      any(l != as.integer(l))) {
    stop("spike_lengths must be positive integers")
  }
  result$spike_lengths <- as.integer(l)
  s <- result$spike_seed
  if (!is.numeric(s) || length(s) != 1L || !is.finite(s)) {
    stop("spike_seed must be one number")
  }
  result
}

has_active_technical_artifacts <- function(technical_artifacts) {
  technical_artifacts$duplication_rate > 0 ||
    technical_artifacts$multimapping_rate > 0 ||
    isTRUE(technical_artifacts$spike_sites > 0L)
}

# Draw the spike sites of one simulation, shared by all its libraries: which
# transcript, where, which length, and in which libraries with how many copies.
# Uses its own seed (spike_seed) and puts the random state back, so the other
# draws of the run do not change. A site lies inside the CDS when `cds_start`
# and `cds_length` (transcript coordinates, named by transcript) are given,
# else anywhere on the transcript. Copies are log-uniform in spike_multiplicity;
# every site is in at least one library.
draw_spike_sites <- function(models, n_libraries, settings,
                             cds_start = NULL, cds_length = NULL) {
  empty <- data.table::data.table(
    spike_id = character(0), transcript_id = character(0),
    transcript_start = integer(0), fragment_length = integer(0),
    sequence = character(0), library = integer(0), copies = integer(0)
  )
  if (settings$spike_sites < 1L || !length(models) || n_libraries < 1L) return(empty)
  with_seed(settings$spike_seed, {
    ids <- names(models)
    lens <- settings$spike_lengths
    sites <- lapply(seq_len(settings$spike_sites), function(k) {
      len <- lens[sample.int(length(lens), 1L)]
      if (!is.null(cds_start)) {
        room <- cds_length[ids] - len + 1L
        ok <- ids[!is.na(room) & room >= 1L]
        if (!length(ok)) return(NULL)
        tx <- ok[sample.int(length(ok), 1L)]
        start <- cds_start[[tx]] + sample.int(room[[tx]], 1L) - 1L
      } else {
        room <- vapply(models, function(m) m$length, numeric(1)) - len + 1L
        ok <- ids[room >= 1L]
        if (!length(ok)) return(NULL)
        tx <- ok[sample.int(length(ok), 1L)]
        start <- sample.int(room[[tx]], 1L)
      }
      libs <- which(stats::runif(n_libraries) < settings$spike_library_fraction)
      if (!length(libs)) libs <- sample.int(n_libraries, 1L)
      lm <- log(settings$spike_multiplicity)
      copies <- as.integer(round(exp(stats::runif(length(libs), lm[1], lm[2]))))
      data.table::data.table(
        spike_id = sprintf("spike_%d", k), transcript_id = tx,
        transcript_start = as.integer(start), fragment_length = len,
        sequence = substr(models[[tx]]$sequence, start, start + len - 1L),
        library = libs, copies = copies
      )
    })
    out <- data.table::rbindlist(sites)
    if (!nrow(out)) empty else out
  })
}

# SAM records for the spike sites of one library: `copies` identical primary
# reads per site (same QNAME stem, same sequence, same position).
spike_records <- function(spikes, models) {
  if (is.null(spikes) || !nrow(spikes)) return(NULL)
  index <- match(spikes$transcript_id, names(models))
  aln <- fragment_alignments(models[unique(index)],
                             match(index, unique(index)),
                             spikes$transcript_start, spikes$fragment_length)
  strand <- vapply(models[spikes$transcript_id], function(m) m$strand, character(1))
  seqnames <- vapply(models[spikes$transcript_id], function(m)
    as.character(GenomicRanges::seqnames(m$exons)[1]), character(1))
  sites <- data.table::data.table(
    seqnames = seqnames, start = aln$pos, end = aln$end, strand = unname(strand),
    cigar = aln$cigar, fragment_length = spikes$fragment_length,
    sequence = spikes$sequence, transcript_id = spikes$transcript_id,
    fragment_id = spikes$spike_id, copies = spikes$copies
  )
  records <- sites[rep.int(seq_len(.N), copies)]
  records[, `:=`(
    qname = paste0(fragment_id, "_", base::sequence(rle(fragment_id)$lengths)),
    reference_sequence = as_reference_orientation(sequence, strand),
    score = 1L,
    flag = ifelse(strand == "+", 0L, 16L),
    mapq = 30L, nh = 1L, is_duplicate = FALSE, is_secondary = FALSE,
    is_spike = TRUE, source_molecule = NA_integer_, copies = NULL
  )]
  records
}

expand_fragment_records <- function(fragment_table) {
  index <- rep.int(seq_len(nrow(fragment_table)), fragment_table$score)
  records <- data.table::copy(fragment_table[index])
  records[, `:=`(
    qname = sprintf("molecule_%d", seq_len(.N)),
    flag = ifelse(strand == "+", 0L, 16L),
    mapq = 30L, nh = 1L, is_duplicate = FALSE,
    is_secondary = FALSE, is_spike = FALSE, source_molecule = seq_len(.N)
  )]
  records
}

# Draw from a vector of values. Unlike sample(), a single numeric value is never
# treated as the range 1:value.
resample_values <- function(values, size, replace = FALSE) {
  values[sample.int(length(values), size, replace = replace)]
}

# SAM reports SEQ relative to the genome's plus strand, while `sequence` is kept
# in the direction the fragment was read off the transcript. A record placed on
# the minus strand therefore stores the reverse complement. This mirrors how
# `reference_sequence` is first built in fragment_geometry.R, and is needed again
# whenever a record is moved to a location whose strand may differ.
as_reference_orientation <- function(sequence, strand) {
  flip <- strand == "-"
  if (any(flip)) {
    sequence[flip] <- as.character(Biostrings::reverseComplement(
      Biostrings::DNAStringSet(sequence[flip])
    ))
  }
  sequence
}

simulate_alignment_artifacts <- function(fragment_table, technical_artifacts,
                                         spikes = NULL, models = NULL) {
  settings <- normalize_technical_artifacts(technical_artifacts)
  primary <- expand_fragment_records(fragment_table)
  original_count <- nrow(primary)

  duplicate_count <- stats::rbinom(
    1L, original_count, settings$duplication_rate
  )
  duplicates <- primary[0]
  if (duplicate_count > 0L) {
    selected <- sample.int(original_count, duplicate_count, replace = FALSE)
    duplicate_index <- rep(selected, each = settings$duplicate_copies)
    duplicates <- data.table::copy(primary[duplicate_index])
    duplicates[, `:=`(
      qname = paste0(qname, "_dup", sequence(.N), "_", seq_len(.N)),
      flag = bitwOr(flag, 1024L), is_duplicate = TRUE
    )]
  }

  alignment_key <- paste(primary$seqnames, primary$start, primary$cigar,
                         primary$strand, sep = ":")
  alternatives <- split(seq_len(original_count), primary$fragment_length)
  eligible <- vapply(seq_len(original_count), function(i) {
    any(alignment_key[alternatives[[as.character(primary$fragment_length[i])]]] !=
          alignment_key[i])
  }, logical(1))
  requested_multi <- stats::rbinom(
    1L, original_count, settings$multimapping_rate
  )
  selected_multi <- if (requested_multi > 0L && any(eligible)) {
    resample_values(which(eligible), min(requested_multi, sum(eligible)))
  } else integer()
  secondary <- primary[0]
  if (length(selected_multi)) {
    primary[selected_multi, nh := 1L + settings$secondary_alignments]
    secondary <- data.table::rbindlist(lapply(selected_multi, function(i) {
      choices <- alternatives[[as.character(primary$fragment_length[i])]]
      choices <- choices[alignment_key[choices] != alignment_key[i]]
      chosen <- resample_values(choices, settings$secondary_alignments,
                                replace = length(choices) < settings$secondary_alignments)
      result <- data.table::copy(primary[chosen])
      # A multimapper is one molecule reported at several locations, so every
      # record has to carry that molecule's own sequence rather than the
      # sequence of the fragment whose coordinates were borrowed. Both columns
      # are set: `sequence` in read direction, and `reference_sequence` turned
      # to the strand of this location, because that is the one written to SAM.
      result[, `:=`(
        qname = primary$qname[i],
        sequence = primary$sequence[i],
        reference_sequence = as_reference_orientation(
          rep(primary$sequence[i], .N), strand
        ),
        flag = bitwOr(ifelse(strand == "+", 0L, 16L), 256L),
        nh = 1L + settings$secondary_alignments,
        is_secondary = TRUE,
        source_molecule = primary$source_molecule[i]
      )]
      result
    }))
  }
  spiked <- spike_records(spikes, models)
  records <- data.table::rbindlist(
    list(primary, duplicates, secondary, spiked), use.names = TRUE, fill = TRUE
  )
  records[is.na(is_spike), is_spike := FALSE]
  attr(records, "diagnostics") <- c(
    biological_molecules = original_count,
    pcr_duplicate_records = nrow(duplicates),
    requested_multimapping_molecules = requested_multi,
    multimapping_molecules = length(selected_multi),
    secondary_alignment_records = nrow(secondary),
    spike_records = if (is.null(spiked)) 0L else nrow(spiked),
    output_records = nrow(records)
  )
  records
}

write_artifact_sam <- function(records, path, seqinfo) {
  header <- c(paste(
    "@SQ", paste0("SN:", GenomeInfoDb::seqnames(seqinfo)),
    paste0("LN:", GenomeInfoDb::seqlengths(seqinfo)), sep = "\t"
  ), simulated_pg_header())
  lines <- paste(
    records$qname, records$flag, records$seqnames, records$start,
    records$mapq, records$cigar, "*", 0L, 0L, records$reference_sequence, "*",
    paste0("NH:i:", records$nh), sep = "\t"
  )
  writeLines(c(header, lines), path)
  path
}

write_artifact_library <- function(records, file_base, format, seqinfo) {
  sam_path <- paste0(file_base, ".sam")
  write_artifact_sam(records, sam_path, seqinfo)
  if (format == "sam") return(c(default = sam_path, sam = sam_path))
  bam_path <- Rsamtools::asBam(sam_path, destination = file_base, overwrite = TRUE)
  unlink(sam_path)
  c(default = bam_path, bam = bam_path)
}

# Spike truth of one library: site, position, length, sequence and copies.
write_spike_ground_truth <- function(spikes, file_base, ground_truth) {
  if (identical(ground_truth, FALSE) || is.null(ground_truth) ||
      is.null(spikes) || !nrow(spikes)) return(NULL)
  path <- if (identical(ground_truth, TRUE)) {
    paste0(file_base, "_spike_truth.tsv")
  } else {
    dir.create(ground_truth, recursive = TRUE, showWarnings = FALSE)
    file.path(ground_truth, paste0(basename(file_base), "_spike_truth.tsv"))
  }
  data.table::fwrite(spikes, path, sep = "\t")
  path
}

# CDS start (first CDS nt, transcript coordinates) and CDS length per transcript,
# named by transcript, for placing spike sites inside the CDS.
cds_transcript_coordinates <- function(cds_ranges, models) {
  ids <- intersect(names(cds_ranges), names(models))
  starts <- vapply(ids, function(tx) {
    ex <- cds_ranges[[tx]]
    g <- if (as.character(GenomicRanges::strand(ex)[1]) == "-") {
      max(GenomicRanges::end(ex))
    } else min(GenomicRanges::start(ex))
    genomic_site_to_transcript(models[[tx]], g)
  }, integer(1))
  lens <- vapply(ids, function(tx) as.integer(sum(GenomicRanges::width(cds_ranges[[tx]]))),
                 integer(1))
  keep <- !is.na(starts)
  list(start = starts[keep], length = lens[keep])
}

write_artifact_ground_truth <- function(records, file_base, ground_truth) {
  if (identical(ground_truth, FALSE) || is.null(ground_truth)) return(NULL)
  path <- if (identical(ground_truth, TRUE)) {
    paste0(file_base, "_artifact_truth.tsv")
  } else {
    dir.create(ground_truth, recursive = TRUE, showWarnings = FALSE)
    file.path(ground_truth, paste0(basename(file_base), "_artifact_truth.tsv"))
  }
  columns <- c(
    "qname", "source_molecule", "fragment_id", "seqnames", "start",
    "strand", "cigar", "fragment_length", "flag", "nh",
    "is_duplicate", "is_secondary", "is_spike"
  )
  data.table::fwrite(records[, ..columns], path, sep = "\t")
  path
}
