# Fill in the defaults and reject settings that cannot mean anything.
#
# Callers supply only the fields they care about, so this completes the rest and
# checks the result once, up front, rather than letting a bad value surface deep
# inside a loop. Unknown field names are rejected rather than ignored: a
# misspelled field would otherwise be silently dropped and the default used, so
# the simulation would quietly not do what was asked.
normalize_fragment_geometry <- function(fragment_geometry) {
  defaults <- list(
    source = "default",
    site_reference = "a_site",
    distribution = NULL,
    site_offset = NULL,
    boundary_action = "renormalize",
    five_prime_bias = list(source = "none"),
    three_prime_bias = list(source = "none"),
    five_prime_flank_bias = list(source = "none"),
    three_prime_flank_bias = list(source = "none"),
    codon_bias = list(source = "none"),
    frame_bias = list(source = "none"),
    added_base = list(base = "T", rate = 0)
  )
  if (is.null(fragment_geometry)) return(defaults)
  stopifnot(is.list(fragment_geometry))
  unknown <- setdiff(names(fragment_geometry), names(defaults))
  if (length(unknown)) {
    stop("Unknown fragment_geometry field(s): ", paste(unknown, collapse = ", "))
  }
  geometry <- utils::modifyList(defaults, fragment_geometry)
  if (!geometry$source %in% c("none", "default", "user", "learned")) {
    stop("fragment_geometry$source must be one of: none, default, user, learned")
  }
  if (!geometry$site_reference %in% c("p_site", "a_site")) {
    stop("fragment_geometry$site_reference must be 'p_site' or 'a_site'")
  }
  if (!is.null(geometry$site_offset) &&
      (!is.numeric(geometry$site_offset) || !length(geometry$site_offset) ||
       anyNA(geometry$site_offset) || any(geometry$site_offset < 0) ||
       any(geometry$site_offset != as.integer(geometry$site_offset)))) {
    stop("fragment_geometry$site_offset must be NULL or non-negative integers")
  }
  if (identical(geometry$boundary_action, "adjust_offset")) {
    stop("fragment_geometry$boundary_action 'adjust_offset' no longer exists: ",
         "offsets are not shifted at transcript ends any more. Use 'renormalize' ",
         "(positions without room for a fragment get no reads) or 'error'")
  }
  if (!geometry$boundary_action %in% c("error", "renormalize")) {
    stop("fragment_geometry$boundary_action must be 'error' or 'renormalize'")
  }
  if (!is.null(geometry$site_offset)) {
    geometry$site_offset <- as.integer(geometry$site_offset)
  }
  for (bias_name in c("five_prime_bias", "three_prime_bias",
                      "five_prime_flank_bias", "three_prime_flank_bias")) {
    bias <- geometry[[bias_name]]
    geometry[[bias_name]] <- normalize_end_bias(bias, bias_name)
  }
  geometry$added_base <- normalize_added_base(geometry$added_base)
  geometry$codon_bias <- normalize_fragment_feature_bias(
    geometry$codon_bias, "codon_bias", "codon"
  )
  geometry$frame_bias <- normalize_fragment_feature_bias(
    geometry$frame_bias, "frame_bias", "frame"
  )
  geometry
}

normalize_fragment_feature_bias <- function(bias, bias_name, feature) {
  if (!is.list(bias) || is.null(bias$source) || length(bias$source) != 1L ||
      !bias$source %in% c("none", "default", "user", "learned")) {
    stop("fragment_geometry$", bias_name,
         " must specify source none, default, user, or learned")
  }
  strength <- if (is.null(bias$strength)) 1 else bias$strength
  if (!is.numeric(strength) || length(strength) != 1L ||
      !is.finite(strength) || strength < 0) {
    stop(bias_name, " strength must be one finite non-negative number")
  }
  if (bias$source %in% c("none", "default")) {
    return(list(source = bias$source, table = NULL, strength = strength))
  }
  table <- data.table::copy(data.table::as.data.table(bias$table))
  required <- c(feature, "fragment_length", "weight")
  if (!all(required %in% names(table)) || !nrow(table)) {
    stop(bias_name, " table requires columns: ", paste(required, collapse = ", "))
  }
  if (feature == "codon") {
    table[, codon := toupper(as.character(codon))]
    if (anyNA(table$codon) || any(!grepl("^[ACGT]{3}$", table$codon))) {
      stop(bias_name, " codons must be three-base DNA strings")
    }
  } else {
    if (anyNA(table$frame) || any(!table$frame %in% 0:2)) {
      stop(bias_name, " frames must be 0, 1, or 2")
    }
    table[, frame := as.integer(frame)]
  }
  if (anyNA(table$fragment_length) ||
      any(table$fragment_length <= 0 | !is.finite(table$fragment_length) |
          table$fragment_length != as.integer(table$fragment_length)) ||
      anyNA(table$weight) || any(!is.finite(table$weight) | table$weight <= 0)) {
    stop(bias_name, " lengths and weights must be finite and positive")
  }
  table[, fragment_length := as.integer(fragment_length)]
  if (anyDuplicated(table[, c(feature, "fragment_length"), with = FALSE])) {
    stop(bias_name, " table contains duplicate feature/length rows")
  }
  list(source = bias$source, table = table[, c(feature, "fragment_length", "weight"), with = FALSE],
       strength = strength)
}

fragment_feature_weight <- function(value, fragment_length, bias, feature) {
  if (is.null(bias$table) || identical(bias$strength, 0)) return(1)
  index <- bias$table[[feature]] == value &
    bias$table$fragment_length == fragment_length
  hit <- bias$table$weight[index]
  if (!length(hit)) 1 else hit[[1L]] ^ bias$strength
}

validate_end_bias_table <- function(table, bias_name = "end bias") {
  table <- data.table::copy(data.table::as.data.table(table))
  if (!all(c("kmer", "weight") %in% names(table))) {
    stop(bias_name, " table requires columns: kmer, weight")
  }
  if (anyNA(table$kmer) || any(!grepl("^[ACGT]+$", toupper(table$kmer))) ||
      any(!is.finite(table$weight)) || any(table$weight <= 0)) {
    stop(bias_name, " weights must be finite and strictly positive, with DNA k-mers")
  }
  table[, kmer := toupper(as.character(kmer))]
  if ("fragment_length" %in% names(table)) {
    if (any(!is.na(table$fragment_length) &
            (!is.finite(table$fragment_length) | table$fragment_length <= 0 |
             table$fragment_length != as.integer(table$fragment_length)))) {
      stop(bias_name, " fragment_length must contain positive integers")
    }
    table[, fragment_length := as.integer(fragment_length)]
  } else table[, fragment_length := NA_integer_]
  if (anyDuplicated(table[, .(kmer, fragment_length)])) {
    stop(bias_name, " table contains duplicate kmer/fragment_length rows")
  }
  table[, .(kmer, weight, fragment_length)]
}

normalize_end_bias <- function(bias, bias_name) {
  if (!is.list(bias) || is.null(bias$source) || length(bias$source) != 1L ||
      !bias$source %in% c("none", "default", "user", "learned")) {
    stop("fragment_geometry$", bias_name,
         " must specify source none, default, user, or learned")
  }
  strength <- if (is.null(bias$strength)) 1 else bias$strength
  if (!is.numeric(strength) || length(strength) != 1L ||
      !is.finite(strength) || strength < 0) {
    stop(bias_name, " strength must be one finite non-negative number")
  }
  if (!is.null(bias$flank)) {
    stop(bias_name, "$flank is not supported: model the nuclease cut site with ",
         "fragment_geometry$five_prime_flank_bias / three_prime_flank_bias (learn_end_bias(flank = TRUE))")
  }
  if (bias$source %in% c("none", "default")) {
    return(list(source = bias$source, table = NULL, strength = strength))
  }
  if (is.null(bias$table)) {
    stop("fragment_geometry$", bias_name, " requires a kmer/weight table")
  }
  list(source = bias$source, table = validate_end_bias_table(bias$table, bias_name),
       strength = strength)
}

# The non-templated 5' base (opt-in, rate 0 = off): library preparation adds `base` in front of the footprint 5'
# end of a share `rate` of the fragments (reverse transcriptase in circularization libraries, mostly T; the TSO rG
# in template switching). See add_nontemplated_5prime_base().
normalize_added_base <- function(added_base) {
  if (is.null(added_base)) added_base <- list()
  if (!is.list(added_base)) stop("fragment_geometry$added_base must be a list(base, rate)")
  unknown <- setdiff(names(added_base), c("base", "rate"))
  if (length(unknown)) stop("Unknown fragment_geometry$added_base field(s): ", paste(unknown, collapse = ", "))
  added_base <- utils::modifyList(list(base = "T", rate = 0), added_base)
  if (!is.character(added_base$base) || length(added_base$base) != 1L ||
      !toupper(added_base$base) %in% c("A", "C", "G", "T")) {
    stop("fragment_geometry$added_base$base must be one of A, C, G, T")
  }
  rate <- added_base$rate
  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate) || rate < 0 || rate > 1) {
    stop("fragment_geometry$added_base$rate must be one number from 0 to 1")
  }
  list(base = toupper(added_base$base), rate = as.numeric(rate))
}

has_added_base <- function(geometry) {
  !is.null(geometry$added_base) && geometry$added_base$rate > 0
}

# Cut-site (flank) k-mers: the bases just outside the fragment, read 5' to 3' along the transcript. The 5' flank
# is the k bases upstream of the 5' end (its last base touches the end), the 3' flank the k bases downstream of the
# 3' end (its first base touches the end). They carry the nuclease cut-site preference, while five_prime_bias and
# three_prime_bias (the bases inside the fragment) carry the ligation / RT chemistry. Where the flank runs past the
# transcript the k-mer is "" and the fragment gets weight 1 for that end.
end_flank_kmer <- function(transcript_sequence, tx_start, tx_end, bias, three_prime = FALSE) {
  if (is.null(bias$table)) return(NULL)
  k <- unique(nchar(bias$table$kmer))
  if (length(k) != 1L) stop("Each end-bias table must use one k-mer length")
  from <- if (three_prime) tx_end + 1L else tx_start - k
  to <- from + k - 1L
  kmer <- substr(transcript_sequence, from, to)
  kmer[from < 1L | to > nchar(transcript_sequence)] <- ""
  kmer
}

# How much a fragment's own end sequences change its chance of being sequenced.
#
# Ligation and digestion are not indifferent to sequence, so two fragments from
# the same position can survive library preparation at different rates purely
# because of the bases at their ends. The two ends are looked up independently
# and their weights multiplied, since each end is handled by its own enzymatic
# step.
#
# A table may give weights per fragment length or once for all lengths; the
# length-specific entry wins where both exist. A k-mer the table does not mention
# is treated as neutral rather than as impossible. `strength` scales the whole
# effect as an exponent, so 0 disables the bias and 2 doubles its reach on the
# log scale.
#
# Takes any number of fragments at once and returns one weight each. The lookup
# is a plain vector match rather than a data.table filter per fragment, because
# it runs for every candidate fragment of every read: filtering the table once
# per fragment cost about 4.5 ms, which made a 300-gene sample with an end bias
# take more than twice as long as one without. Each distinct k-mer and length is
# looked up only once, since millions of candidate fragments share a few hundred
# of them.
end_bias_weight <- function(sequence, fragment_length,
                            five_bias, three_bias,
                            five_flank_bias = NULL, three_flank_bias = NULL,
                            five_flank_kmer = NULL, three_flank_kmer = NULL) {
  lookup <- function(kmer, bias) {
    if (is.null(bias$table) || identical(bias$strength, 0)) {
      return(rep(1, length(sequence)))
    }
    distinct_kmer <- unique(kmer)
    distinct_length <- unique(fragment_length)
    # Every combination of a distinct k-mer and a distinct length.
    value <- rep(toupper(distinct_kmer), times = length(distinct_length))
    value_length <- rep(distinct_length, each = length(distinct_kmer))
    table_kmer <- bias$table$kmer
    table_length <- bias$table$fragment_length
    table_weight <- bias$table$weight
    # A weight given for this fragment length wins over one given for all
    # lengths; a k-mer the table does not mention stays neutral.
    by_length <- !is.na(table_length)
    weight <- table_weight[by_length][match(
      paste(value, value_length),
      paste(table_kmer[by_length], table_length[by_length])
    )]
    general <- is.na(weight)
    weight[general] <- table_weight[!by_length][
      match(value[general], table_kmer[!by_length])
    ]
    weight[is.na(weight)] <- 1
    weight <- weight ^ bias$strength
    weight[match(kmer, distinct_kmer) +
             (match(fragment_length, distinct_length) - 1L) * length(distinct_kmer)]
  }
  extract_kmer <- function(bias, from_end = FALSE) {
    if (is.null(bias$table)) return(NULL)
    k <- unique(nchar(bias$table$kmer))
    if (length(k) != 1L) stop("Each end-bias table must use one k-mer length")
    if (from_end) substr(sequence, nchar(sequence) - k + 1L, nchar(sequence))
    else substr(sequence, 1L, k)
  }
  weight <- lookup(extract_kmer(five_bias), five_bias) *
    lookup(extract_kmer(three_bias, TRUE), three_bias)
  # Cut-site terms (bases outside the fragment), multiplied in like a third and fourth end; a flank k-mer the
  # table does not list (also "" past the transcript end) is neutral.
  if (!is.null(five_flank_bias$table)) weight <- weight * lookup(five_flank_kmer, five_flank_bias)
  if (!is.null(three_flank_bias$table)) weight <- weight * lookup(three_flank_kmer, three_flank_bias)
  weight
}

end_bias_kmer <- function(sequence, bias, from_end = FALSE) {
  if (is.null(bias$table)) return(NA_character_)
  k <- unique(nchar(bias$table$kmer))
  if (length(k) != 1L) stop("Each end-bias table must use one k-mer length")
  if (from_end) substr(sequence, nchar(sequence) - k + 1L, nchar(sequence))
  else substr(sequence, 1L, k)
}

#' Create a simple synthetic sequence-end bias profile.
#'
#' Library preparation is not indifferent to sequence: ligation and digestion
#' favour some bases at a fragment's ends over others, so the observed ends are
#' not a faithful record of where the enzyme cut. This builds a deliberately
#' simple version of that effect, for testing whether an analysis notices it.
#'
#' The result is passed as `table` in `five_prime_bias` or `three_prime_bias`.
#' Weights are relative: 1 is neutral, above 1 makes a fragment ending in that
#' k-mer more likely to be kept, below 1 less likely.
#' @param k integer, how many bases at the end the bias depends on. 1 makes it
#'   depend on the final base alone, 2 on the final two, and so on.
#' @param enriched_kmer optional k-mer to favour; every other k-mer stays
#'   neutral.
#' @param enriched_weight how strongly to favour it. 4 means fragments ending
#'   in that k-mer are four times as likely to survive as any other.
#' @param depleted_kmer optional k-mer to disfavour.
#' @param depleted_weight how strongly to disfavour it; below 1.
#' @return a data.table with `kmer` and `weight` columns.
#' @export
make_synthetic_end_bias <- function(k = 1L, enriched_kmer = "G",
                                    enriched_weight = 2,
                                    depleted_kmer = NULL,
                                    depleted_weight = 0.5) {
  if (length(k) != 1L || is.na(k) || k < 1 || k != as.integer(k)) {
    stop("k must be one positive integer", call. = FALSE)
  }
  k <- as.integer(k)
  kmers <- apply(expand.grid(rep(list(c("A", "C", "G", "T")), k)), 1L,
                 paste0, collapse = "")
  result <- data.table::data.table(kmer = kmers, weight = 1)
  update <- function(kmer, new_weight) {
    if (is.null(kmer)) return(invisible(NULL))
    if (length(kmer) != 1L || nchar(kmer) != k ||
        !grepl("^[ACGT]+$", toupper(kmer)) || !is.finite(new_weight) || new_weight <= 0) {
      stop("Synthetic end-bias k-mers and weights are invalid", call. = FALSE)
    }
    value <- toupper(kmer)
    result[kmer == value, weight := new_weight]
  }
  update(enriched_kmer, enriched_weight)
  update(depleted_kmer, depleted_weight)
  result[]
}

# The biologically plausible default P-site (ribosome peptidyl-site) offset
# from a fragment's 5' end, by fragment length, plus 3 nt more for the
# A-site (one codon further into the ribosome). This is the simulator's own
# default geometry (default_fragment_distribution() below); elsewhere it is
# only an initial guess used to center a search window when learning real
# offsets (learn_fragment_geometry(), and the rust_helpers analysis scripts
# via coverageSim:::default_ribosome_site_offset(), kept in sync with this
# rather than an independently duplicated copy of the same formula).
default_ribosome_site_offset <- function(fragment_length, site_reference) {
  p_offset <- ifelse(fragment_length <= 27L, 11L,
                     ifelse(fragment_length <= 30L, 12L, 13L))
  as.integer(p_offset + ifelse(site_reference == "a_site", 3L, 0L))
}

# Pooled RFP read-length distribution of 702 human runs (mode 26-35 nt, each run
# normalised over 25-36 nt, then averaged; carina-claude C2, 2026-11-18).
# Lengths outside 25-36 get a normal(30.5, 2.3) weight on the same scale.
measured_rfp_length_weights <- function(lengths) {
  table_p <- c(`25` = .046, `26` = .062, `27` = .073, `28` = .090, `29` = .123,
               `30` = .157, `31` = .152, `32` = .105, `33` = .071, `34` = .054,
               `35` = .045, `36` = .022)
  w <- unname(table_p[as.character(lengths)])
  outside <- is.na(w)
  w[outside] <- stats::dnorm(lengths[outside], mean = 30.5, sd = 2.3) *
    (.157 / stats::dnorm(30, mean = 30.5, sd = 2.3))
  w
}

default_fragment_distribution <- function(fragment_lengths,
                                          site_reference = "a_site") {
  length_counts <- table(as.integer(fragment_lengths))
  lengths <- as.integer(names(length_counts))
  length_weight <- measured_rfp_length_weights(lengths) *
    as.numeric(length_counts)
  if (!any(length_weight > 0)) length_weight <- as.numeric(length_counts)
  offset <- default_ribosome_site_offset(lengths, site_reference)
  offset <- pmin(offset, lengths - 1L)
  data.table::data.table(
    fragment_length = lengths,
    site_offset = as.integer(offset),
    probability = length_weight / sum(length_weight)
  )
}

validate_fragment_distribution <- function(distribution) {
  distribution <- data.table::as.data.table(distribution)
  required <- c("fragment_length", "site_offset", "probability")
  if (!all(required %in% names(distribution))) {
    stop("Fragment distribution requires columns: ", paste(required, collapse = ", "))
  }
  distribution <- distribution[, ..required]
  if (!nrow(distribution) || anyNA(distribution) ||
      any(distribution$fragment_length <= 0) ||
      any(distribution$fragment_length != as.integer(distribution$fragment_length)) ||
      any(distribution$site_offset < 0) ||
      any(distribution$site_offset != as.integer(distribution$site_offset)) ||
      any(distribution$site_offset >= distribution$fragment_length) ||
      any(!is.finite(distribution$probability)) ||
      any(distribution$probability < 0) || sum(distribution$probability) <= 0) {
    stop("Fragment distribution contains invalid lengths, offsets, or probabilities")
  }
  distribution[, `:=`(
    fragment_length = as.integer(fragment_length),
    site_offset = as.integer(site_offset)
  )]
  distribution <- distribution[, .(probability = sum(probability)),
                               by = .(fragment_length, site_offset)]
  distribution[, probability := probability / sum(probability)]
  distribution
}

resolve_fragment_distribution <- function(fragment_lengths, geometry) {
  if (!is.null(geometry$distribution)) {
    if (geometry$source == "none") {
      stop("source = 'none' cannot be combined with a fragment distribution")
    }
    return(validate_fragment_distribution(geometry$distribution))
  }
  if (geometry$source %in% c("user", "learned")) {
    if (is.null(geometry$site_offset)) {
      stop("A user/learned source requires distribution or site_offset")
    }
  }
  if (!is.null(geometry$site_offset)) {
    distribution <- expand.grid(
      fragment_length = as.integer(fragment_lengths),
      site_offset = geometry$site_offset,
      KEEP.OUT.ATTRS = FALSE
    )
    distribution$probability <- 1
    return(validate_fragment_distribution(distribution))
  }
  if (geometry$source == "none") {
    stop("Simulated RPFs require a default, user, or learned geometry source")
  }
  default_fragment_distribution(fragment_lengths, geometry$site_reference)
}

# Shared core for building the per-transcript position/strand model used by
# fragment-geometry projection and placement, learned end-bias, and
# region-proportion/fragment-geometry learning: exon coordinates, each
# exon's cumulative transcript-position start, transcript length, and
# strand. fasta_file, when given, also fetches and attaches the transcript
# sequence (needed to look up sequence/codon motifs during simulation and
# end-bias learning; the geometry-learning callers don't need it and don't
# pass it). check_chromosome additionally rejects a transcript whose exons
# span more than one chromosome, which only the geometry-learning callers
# have historically required.
# Precompute, once per transcript, everything the fragment code keeps needing.
#
# Placing fragments means repeatedly converting between genome and transcript
# coordinates, which needs the exons in transcript order, where each one starts
# along the transcript, and the total length. Working that out per fragment would
# repeat it millions of times, so it is done once here and carried along.
#
# Exons are sorted by rank rather than by position, because on the minus strand
# the first exon of the transcript has the highest genomic coordinate.
build_transcript_models <- function(ranges, fasta_file = NULL,
                                    check_chromosome = FALSE) {
  sequences <- if (!is.null(fasta_file)) ORFik::txSeqsFromFa(ranges, fasta_file)
  models <- lapply(seq_along(ranges), function(i) {
    exons <- ranges[[i]]
    exon_rank <- S4Vectors::mcols(exons)$exon_rank
    if (!is.null(exon_rank)) {
      ord <- order(exon_rank)
      if (!identical(ord, seq_along(exon_rank))) exons <- exons[ord]
    }
    strand_value <- as.character(unique(GenomicRanges::strand(exons)))
    strand_ok <- length(strand_value) == 1L && strand_value %in% c("+", "-")
    if (check_chromosome) {
      chromosome <- unique(as.character(GenomicRanges::seqnames(exons)))
      if (!strand_ok || length(chromosome) != 1L) {
        stop("Transcript has inconsistent chromosome or strand: ", names(ranges)[i])
      }
    } else if (!strand_ok) {
      stop("Transcript has missing or inconsistent strand: ", names(ranges)[i])
    }
    model <- list(
      transcript_id = names(ranges)[i],
      exons = exons,
      cumulative_start = cumsum(c(1L, head(GenomicRanges::width(exons), -1L))),
      length = sum(GenomicRanges::width(exons)),
      strand = strand_value
    )
    if (!is.null(fasta_file)) model$sequence <- as.character(sequences[[i]])
    model
  })
  names(models) <- names(ranges)
  models
}

transcript_models <- function(mrna_ranges, fasta_file) {
  # fasta_file is required here (unlike build_transcript_models()'s optional
  # default): every caller needs the attached sequence, so fail fast on a
  # missing one instead of silently building sequence-less models.
  if (is.null(fasta_file)) stop("fasta_file is required by transcript_models()")
  build_transcript_models(mrna_ranges, fasta_file = fasta_file)
}

# How many nucleotides on each side of a footprint the RNase kernel reaches:
# the RFP kernel's half-width, rounded down. Shared by both fragment_mode's
# RNase-extension helpers (this file's append_rnase_to_simulated_rpf_table()
# and Sim_reads_from_counts_helpers.R's append_rnase_to_dt()) so the two
# fragment modes always extend by the same amount for the same kernel.
rnase_kernel_reach <- function(rnase_bias) {
  floor(length(rnase_bias[["RFP"]]) / 2L)
}

genomic_site_to_transcript <- function(model, genomic_position) {
  exon_index <- which(
    GenomicRanges::start(model$exons) <= genomic_position &
      GenomicRanges::end(model$exons) >= genomic_position
  )
  if (length(exon_index) != 1L) return(NA_integer_)
  i <- exon_index[[1]]
  within_exon <- if (model$strand == "+") {
    genomic_position - GenomicRanges::start(model$exons)[i]
  } else {
    GenomicRanges::end(model$exons)[i] - genomic_position
  }
  as.integer(model$cumulative_start[i] + within_exon)
}

transcript_position_to_genomic <- function(model, transcript_position) {
  transcript_position <- pmin(pmax(as.integer(transcript_position), 1L), model$length)
  exon_index <- findInterval(transcript_position, model$cumulative_start)
  within_exon <- transcript_position - model$cumulative_start[exon_index]
  if (model$strand == "+") {
    GenomicRanges::start(model$exons)[exon_index] + within_exon
  } else {
    GenomicRanges::end(model$exons)[exon_index] - within_exon
  }
}

append_rnase_to_simulated_rpf_table <- function(dt_range, rnase_bias,
                                                transcript_models) {
  reach <- rnase_kernel_reach(rnase_bias)
  if (reach == 0L) return(dt_range)
  # Preserve split/rbindlist semantics for unusual grouping inputs.
  if (nrow(dt_range) && !anyNA(dt_range$genes) && !is.factor(dt_range$genes) &&
      is.null(data.table::key(dt_range)) && is.null(attr(dt_range, "index"))) {
    groups <- split(seq_len(nrow(dt_range)), dt_range$genes)
    indices <- starts <- ends <- positions <- vector("list", length(groups))
    has_pos <- "region_position" %in% names(dt_range)
    for (i in seq_along(groups)) {
      rows <- groups[[i]]
      first <- rows[1L]
      last <- rows[length(rows)]
      model <- transcript_models[[dt_range$transcript_id[first]]]
      first_tx <- genomic_site_to_transcript(model, dt_range$start[first])
      last_tx <- genomic_site_to_transcript(model, dt_range$start[last])
      indices[[i]] <- c(rep(first, reach), rows, rep(last, reach))
      starts[[i]] <- c(transcript_position_to_genomic(
        model, seq.int(first_tx - reach, first_tx - 1L)), dt_range$start[rows],
        transcript_position_to_genomic(model, seq.int(last_tx + 1L, last_tx + reach)))
      ends[[i]] <- c(head(starts[[i]], reach), dt_range$end[rows],
                     tail(starts[[i]], reach))
      if (has_pos) positions[[i]] <- c(
        seq.int(dt_range$region_position[first] - reach, dt_range$region_position[first] - 1L),
        dt_range$region_position[rows],
        seq.int(dt_range$region_position[last] + 1L, dt_range$region_position[last] + reach))
    }
    out <- dt_range[unlist(indices, use.names = FALSE)]
    data.table::set(out, j = "start", value = unlist(starts, use.names = FALSE))
    data.table::set(out, j = "end", value = unlist(ends, use.names = FALSE))
    if (has_pos) data.table::set(out, j = "region_position",
                                  value = unlist(positions, use.names = FALSE))
    return(out)
  }
  groups <- split(dt_range, dt_range$genes, keep.by = TRUE)
  data.table::rbindlist(lapply(groups, function(group) {
    model <- transcript_models[[group$transcript_id[1]]]
    first_tx <- genomic_site_to_transcript(model, group$start[1])
    last_tx <- genomic_site_to_transcript(model, group$start[nrow(group)])
    prefix <- group[rep(1L, reach)]
    suffix <- group[rep(nrow(group), reach)]
    prefix[, start := transcript_position_to_genomic(
      model, seq.int(first_tx - reach, first_tx - 1L)
    )]
    suffix[, start := transcript_position_to_genomic(
      model, seq.int(last_tx + 1L, last_tx + reach)
    )]
    prefix[, end := start]
    suffix[, end := start]
    if ("region_position" %in% names(group)) {
      prefix[, region_position := seq.int(
        group$region_position[1L] - reach, group$region_position[1L] - 1L
      )]
      suffix[, region_position := seq.int(
        group$region_position[nrow(group)] + 1L,
        group$region_position[nrow(group)] + reach
      )]
    }
    data.table::rbindlist(list(prefix, group, suffix))
  }))
}

# Translate a fragment's position along the transcript into genome coordinates.
#
# A fragment is defined by where it sits in the spliced transcript, but an
# alignment file describes the genome, where the same fragment may be split
# across exons with introns in between. This walks the exons the fragment covers
# and writes the result as a CIGAR string: stretches that match the genome as M,
# the skipped introns as N.
#
# Minus-strand transcripts run against the genome, so their coordinates are
# measured from the exon's far end; the pieces are then sorted into ascending
# genomic order, because that is the order an alignment file expects regardless
# of which strand the gene is on.
transcript_fragment_alignment <- function(model, transcript_start,
                                          fragment_length) {
  transcript_end <- transcript_start + fragment_length - 1L
  exon_tx_start <- model$cumulative_start
  exon_tx_end <- exon_tx_start + GenomicRanges::width(model$exons) - 1L
  used <- which(exon_tx_start <= transcript_end & exon_tx_end >= transcript_start)

  pieces <- lapply(used, function(i) {
    overlap_start <- max(transcript_start, exon_tx_start[i])
    overlap_end <- min(transcript_end, exon_tx_end[i])
    offset_start <- overlap_start - exon_tx_start[i]
    offset_end <- overlap_end - exon_tx_start[i]
    if (model$strand == "+") {
      genomic_start <- GenomicRanges::start(model$exons)[i] + offset_start
      genomic_end <- GenomicRanges::start(model$exons)[i] + offset_end
    } else {
      genomic_start <- GenomicRanges::end(model$exons)[i] - offset_end
      genomic_end <- GenomicRanges::end(model$exons)[i] - offset_start
    }
    c(start = genomic_start, end = genomic_end)
  })
  pieces <- do.call(rbind, pieces)
  pieces <- pieces[order(pieces[, "start"]), , drop = FALSE]
  widths <- pieces[, "end"] - pieces[, "start"] + 1L
  if (length(widths) == 1L) {
    cigar <- paste0(widths, "M")
  } else {
    introns <- pieces[-1L, "start"] -
      pieces[-nrow(pieces), "end"] - 1L
    cigar <- paste0(
      paste0(widths, "M"),
      c(paste0(introns, "N"), ""),
      collapse = ""
    )
  }
  list(
    pos = as.integer(min(pieces[, "start"])),
    end = as.integer(max(pieces[, "end"])),
    cigar = cigar
  )
}

# Genome coordinates and CIGAR of many fragments at once: the same result as
# transcript_fragment_alignment() for each fragment, without a loop over them.
# `models` are transcript models, `model_index` says which one each fragment
# belongs to, and the fragments are given in transcript coordinates.
fragment_alignments <- function(models, model_index, transcript_start,
                                fragment_length) {
  n <- length(transcript_start)
  if (!n) return(list(pos = integer(0), end = integer(0), cigar = character(0)))
  # All exons of all models in one table, in transcript order. `base` moves
  # each model's transcript coordinates past those of the models before it,
  # so one sorted vector of exon starts serves every fragment.
  exon_count <- vapply(models, function(model) length(model$exons), integer(1))
  model_length <- vapply(models, function(model) as.numeric(model$length), numeric(1))
  base <- c(0, cumsum(model_length))[seq_along(models)]
  exon_model <- rep(seq_along(models), exon_count)
  exon_tx_start <- unlist(lapply(models, function(model) model$cumulative_start))
  exon_width <- unlist(lapply(models, function(model) GenomicRanges::width(model$exons)))
  exon_genome_start <- unlist(lapply(models, function(model) GenomicRanges::start(model$exons)))
  exon_genome_end <- unlist(lapply(models, function(model) GenomicRanges::end(model$exons)))
  exon_minus <- rep(vapply(models, function(model) model$strand != "+", logical(1)),
                    exon_count)
  exon_key <- base[exon_model] + exon_tx_start

  # The first and last exon each fragment touches; it covers the ones between.
  transcript_end <- transcript_start + fragment_length - 1L
  first <- findInterval(base[model_index] + transcript_start, exon_key)
  last <- findInterval(base[model_index] + transcript_end, exon_key)
  piece_count <- last - first + 1L
  fragment <- rep(seq_len(n), piece_count)
  exon <- first[fragment] + sequence(piece_count) - 1L

  # The part of each exon the fragment covers, in genome coordinates.
  offset_start <- pmax(transcript_start[fragment], exon_tx_start[exon]) - exon_tx_start[exon]
  offset_end <- pmin(transcript_end[fragment], exon_tx_start[exon] + exon_width[exon] - 1L) -
    exon_tx_start[exon]
  minus <- exon_minus[exon]
  piece_start <- ifelse(minus, exon_genome_end[exon] - offset_end,
                        exon_genome_start[exon] + offset_start)
  piece_end <- ifelse(minus, exon_genome_end[exon] - offset_start,
                      exon_genome_start[exon] + offset_end)

  # Pieces left to right on the genome; a gap between two is an intron (N).
  order_pieces <- order(fragment, piece_start)
  fragment <- fragment[order_pieces]
  piece_start <- piece_start[order_pieces]
  piece_end <- piece_end[order_pieces]
  first_piece <- !duplicated(fragment)
  last_piece <- !duplicated(fragment, fromLast = TRUE)
  intron <- as.integer(c(piece_start[-1L], NA) - piece_end - 1L)
  piece_text <- paste0(as.integer(piece_end - piece_start + 1L), "M",
                       ifelse(last_piece, "", paste0(intron, "N")))
  cigar <- piece_text[first_piece]
  end <- piece_end[last_piece]
  spliced <- which(piece_count > 1L)
  if (length(spliced)) {
    in_spliced <- fragment %in% spliced
    cigar[spliced] <- vapply(split(piece_text[in_spliced], fragment[in_spliced]),
                             paste, character(1), collapse = "")
    # The rightmost end, also if one exon lies inside another.
    end[spliced] <- vapply(split(piece_end[in_spliced], fragment[in_spliced]),
                           max, numeric(1))
  }
  list(
    pos = as.integer(piece_start[first_piece]),
    end = as.integer(end),
    cigar = cigar
  )
}

# Whether candidate weights depend on the fragment sequence: an end-bias table
# (k-mers at the fragment ends), a cut-site flank table (bases just outside the
# fragment) or a codon-bias table (codon at the site) with
# non-zero strength. Mirrors the checks in end_bias_weight() and
# fragment_feature_weights(), which return weight 1 otherwise.
needs_candidate_sequence <- function(geometry) {
  active <- function(bias) !is.null(bias$table) && !identical(bias$strength, 0)
  active(geometry$five_prime_bias) || active(geometry$three_prime_bias) ||
    active(geometry$five_prime_flank_bias) || active(geometry$three_prime_flank_bias) ||
    active(geometry$codon_bias)
}

# Which joint fragment draw to use: "joint" (default, one multinomial over all
# candidates of a group) or "two_stage" (same distribution, positions first, then
# length and offset; faster, but different reads for the same seed). The
# two-stage draw is used only for joint draws without sequence-dependent biases
# and without boundary_action = "error"; otherwise the joint draw runs.
fragment_draw_option <- function() {
  draw <- getOption("coverageSim.fragment_draw", "joint")
  if (!is.character(draw) || length(draw) != 1L || !draw %in% c("joint", "two_stage")) {
    stop('options(coverageSim.fragment_draw) must be "joint" or "two_stage"')
  }
  draw
}

has_active_fragment_bias <- function(geometry) {
  bias_names <- c("five_prime_bias", "three_prime_bias", "five_prime_flank_bias",
                  "three_prime_flank_bias", "codon_bias", "frame_bias")
  any(vapply(geometry[bias_names], function(bias) {
    !is.null(bias$table) && bias$strength > 0
  }, logical(1)))
}

# Transcript positions of several genome positions on one transcript; NA where
# a position is not inside exactly one exon. The same arithmetic as
# genomic_site_to_transcript(), for all positions at once.
genomic_sites_to_transcript <- function(model, genomic_positions) {
  exon_start <- GenomicRanges::start(model$exons)
  exon_end <- GenomicRanges::end(model$exons)
  inside <- outer(genomic_positions, exon_start, ">=") &
    outer(genomic_positions, exon_end, "<=")
  exon_index <- max.col(inside, ties.method = "first")
  within_exon <- if (model$strand == "+") {
    genomic_positions - exon_start[exon_index]
  } else {
    exon_end[exon_index] - genomic_positions
  }
  site <- as.integer(model$cumulative_start[exon_index] + within_exon)
  site[rowSums(inside) != 1L] <- NA_integer_
  site
}

# Codon or frame weights for many fragments: each distinct value and length is
# looked up once with fragment_feature_weight().
fragment_feature_weights <- function(values, fragment_lengths, bias, feature) {
  if (is.null(bias$table) || identical(bias$strength, 0)) {
    return(rep(1, length(values)))
  }
  # One number per distinct (value, length) pair; cheaper than pasting the
  # two together for millions of candidate fragments.
  length_id <- match(fragment_lengths, unique(fragment_lengths))
  key <- (match(values, unique(values)) - 1) * max(length_id, 0L) + length_id
  first <- !duplicated(key)
  weights <- vapply(which(first), function(j) {
    fragment_feature_weight(values[j], fragment_lengths[j], bias, feature)
  }, numeric(1))
  weights[match(key, key[first])]
}

# Every fragment each signal position could give rise to, with its weights.
#
# Each position is combined with every length-and-offset pair of `candidates`,
# and the pairs whose fragment fits inside the transcript are kept. Positions
# near a transcript end where nothing fits are handled by `boundary_action`
# (see make_simulated_rpf_fragments()). Rows come in signal-row order, then in
# the order of `candidates`.
fragment_candidates <- function(signal_table, transcript_models, candidates,
                                geometry) {
  n_signal <- nrow(signal_table)
  transcript_id <- signal_table$transcript_id
  used_ids <- unique(transcript_id)
  used_models <- lapply(seq_along(used_ids), function(k) {
    transcript_models[[used_ids[k]]]
  })
  model_index <- match(transcript_id, used_ids)
  has_model <- !vapply(used_models, is.null, logical(1))[model_index]

  site_tx <- rep(NA_integer_, n_signal)
  rows_per_model <- split(seq_len(n_signal), model_index)
  for (k in as.integer(names(rows_per_model))) {
    if (is.null(used_models[[k]])) next
    rows <- rows_per_model[[as.character(k)]]
    site_tx[rows] <- genomic_sites_to_transcript(
      used_models[[k]], signal_table$signal_position[rows]
    )
  }
  transcript_length <- vapply(used_models, function(model) {
    if (is.null(model)) NA_real_ else as.numeric(model$length)
  }, numeric(1))[model_index]

  # Every position against every length-and-offset pair; keep what fits.
  n_pairs <- nrow(candidates)
  signal_row <- rep(seq_len(n_signal), each = n_pairs)
  candidate_index <- rep(seq_len(n_pairs), times = n_signal)
  tx_start <- site_tx[signal_row] - candidates$site_offset[candidate_index]
  fits <- !is.na(tx_start) & tx_start >= 1L &
    tx_start + candidates$fragment_length[candidate_index] - 1L <=
      transcript_length[signal_row]
  pairs <- data.table::data.table(
    signal_row = signal_row[fits],
    fragment_length = candidates$fragment_length[candidate_index[fits]],
    site_offset = candidates$site_offset[candidate_index[fits]],
    probability = candidates$probability[candidate_index[fits]]
  )

  # With "renormalize", a position where no fragment fits simply has no
  # candidates; make_simulated_rpf_fragments() moves or drops its reads.
  nothing_fits <- has_model & !is.na(site_tx) &
    tabulate(pairs$signal_row, nbins = n_signal) == 0L &
    geometry$boundary_action == "error"

  # Stop at the first position that cannot be placed, with the same message as
  # when positions were placed one at a time.
  problem <- which(!has_model | is.na(site_tx) | nothing_fits)
  if (length(problem)) {
    i <- problem[1L]
    if (!has_model[i]) stop("No transcript model for: ", transcript_id[i])
    if (is.na(site_tx[i])) {
      stop("Signal position is not exonic in transcript ", transcript_id[i])
    }
    stop(
      "Cannot create a complete fragment for transcript ", transcript_id[i],
      " at signal position ", signal_table$signal_position[i],
      "; affected count: ", signal_table$score[i]
    )
  }

  pairs[, `:=`(
    tx_start = site_tx[signal_row] - site_offset,
    tx_end = site_tx[signal_row] - site_offset + fragment_length - 1L
  )]
  # Cutting out the sequence of every candidate is the costly part (tens of
  # millions of candidates, of which only a few become reads). It is needed here
  # only when an end-bias or codon-bias table weights the candidates; otherwise
  # the columns stay NA and make_simulated_rpf_fragments() fills them for the
  # drawn fragments only, with the same values.
  if (needs_candidate_sequence(geometry)) {
    transcript_sequence <- vapply(used_models, function(model) model$sequence,
                                  character(1))[model_index[pairs$signal_row]]
    pairs[, sequence := substr(transcript_sequence, tx_start, tx_end)]
    pairs[, codon := substr(transcript_sequence, site_tx[signal_row],
                            site_tx[signal_row] + 2L)]
  } else {
    pairs[, `:=`(sequence = NA_character_, codon = NA_character_)]
  }
  pairs[, frame := if ("region_position" %in% names(signal_table)) {
    as.integer((signal_table$region_position[signal_row] - 1L) %% 3L)
  } else 0L]
  five_flank <- three_flank <- NULL
  if (needs_candidate_sequence(geometry)) {
    five_flank <- end_flank_kmer(transcript_sequence, pairs$tx_start, pairs$tx_end,
                                 geometry$five_prime_flank_bias)
    three_flank <- end_flank_kmer(transcript_sequence, pairs$tx_start, pairs$tx_end,
                                  geometry$three_prime_flank_bias, three_prime = TRUE)
  }
  pairs[, end_bias_weight := end_bias_weight(
    sequence, fragment_length,
    geometry$five_prime_bias, geometry$three_prime_bias,
    geometry$five_prime_flank_bias, geometry$three_prime_flank_bias,
    five_flank, three_flank
  )]
  pairs[, codon_bias_weight := fragment_feature_weights(
    codon, fragment_length, geometry$codon_bias, "codon"
  )]
  pairs[, frame_bias_weight := fragment_feature_weights(
    frame, fragment_length, geometry$frame_bias, "frame"
  )]
  pairs[, bias_weight := end_bias_weight * codon_bias_weight * frame_bias_weight]
  pairs[, probability := probability / sum(probability), by = signal_row]
  pairs[, geometry_probability := probability]
  pairs
}

# Preserve region budgets when present; direct callers supply counts per transcript.
# Decide which of the possible fragments actually get sequenced.
#
# Each signal position could give rise to fragments of several lengths and
# several offsets. Every such candidate is weighted by three things multiplied
# together: how likely that length-and-offset combination is at all, how much the
# sequence at its ends helps or hinders it surviving library preparation, and how
# much signal the position carries. The region's read budget is then drawn over
# the candidates in those proportions, so the total comes out exactly right
# rather than approximately.
#
# `joint` decides what competes with what. Drawing per transcript lets fragments
# compete across the whole transcript, which is what happens in a real library;
# drawing per position keeps each position's count fixed, which is useful when
# the point is to reproduce a given profile exactly.
#
# A group none of whose positions has room for a complete fragment cannot
# receive its reads; they are left out and listed in the "dropped_reads"
# attribute of the result.
allocate_fragment_candidates <- function(candidates, signal_table, joint) {
  # `labels` is the ground-truth region of each position (e.g. "uorf:3");
  # `groups` decides what is drawn together. They differ when drawing per
  # position: the draw then uses the row, but every fragment must still carry
  # its region label, or reads cannot be attributed to their feature.
  labels <- fragment_sampling_groups(signal_table)
  groups <- if (joint) labels else seq_len(nrow(signal_table))
  # The read budget of each group, from all of its positions.
  budgets <- data.table::data.table(draw_group = groups,
                                    sampling_group = labels,
                                    transcript_id = signal_table$transcript_id)
  if ("read_count" %in% names(signal_table)) {
    budgets[, read_count := signal_table$read_count]
    budgets <- budgets[, .(sampling_group = sampling_group[1],
                           transcript_id = transcript_id[1],
                           total = if (data.table::uniqueN(read_count) == 1L) {
                             read_count[1]
                           } else NA_real_), by = draw_group]
  } else {
    budgets[, score := signal_table$score]
    budgets <- budgets[, .(sampling_group = sampling_group[1],
                           transcript_id = transcript_id[1],
                           total = round(sum(score))), by = draw_group]
  }
  if (anyNA(budgets$total) || any(budgets$total < 0) ||
      any(budgets$total != floor(budgets$total))) {
    stop("Invalid fragment sampling budget")
  }
  # Groups are handled as integer ids of the signal rows: the candidate table
  # has tens of millions of rows, and matching or grouping it by character
  # labels cost far more than the draws themselves.
  group_id <- match(groups, budgets$draw_group)
  present <- unique(group_id[candidates$signal_row])
  dropped <- budgets[total > 0 & !seq_len(.N) %in% present,
                     .(reads = as.numeric(sum(total))),
                     by = .(sampling_group = as.character(sampling_group),
                            transcript_id = as.character(transcript_id))]

  candidates[, sampling_group := labels[signal_row]]
  candidates[, probability := geometry_probability * bias_weight *
               signal_table$score[signal_row]]
  if (any(!is.finite(candidates$probability))) {
    stop("End-bias strength produces non-finite weights; reduce strength")
  }
  # No signal positions at all: there is nothing to draw, and the table stays
  # empty instead of reaching the draw with no budget.
  if (!nrow(candidates)) {
    candidates[, score := integer(0)]
  } else {
    # One multinomial draw per group, in the order in which the groups first
    # appear in the candidate table (the order data.table's `by` used, so the
    # random draws and every number are unchanged).
    candidate_group <- group_id[candidates$signal_row]
    group_rank <- match(candidate_group, unique(candidate_group))
    rows <- split(seq_along(group_rank), group_rank)
    probability <- candidates$probability
    score <- integer(length(probability))
    for (k in seq_along(rows)) {
      idx <- rows[[k]]
      w <- probability[idx]
      if (sum(w) <= 0) stop("End-bias weights exclude every feasible fragment")
      w <- w / max(w)
      p <- w / sum(w)
      probability[idx] <- p
      score[idx] <- as.integer(stats::rmultinom(
        1L, budgets$total[candidate_group[idx[1L]]], p))
    }
    data.table::set(candidates, j = "probability", value = probability)
    data.table::set(candidates, j = "score", value = score)
  }
  data.table::setattr(candidates, "dropped_reads", dropped)
  candidates
}

# The groups whose read budgets are kept: the region of a transcript when the
# simulation labelled it (sampling_group), otherwise the whole transcript.
fragment_sampling_groups <- function(signal_table) {
  groups <- if ("sampling_group" %in% names(signal_table)) {
    signal_table$sampling_group
  } else signal_table$transcript_id
  if (anyNA(groups)) stop("sampling_group and transcript_id must not be missing")
  groups
}

# Reads at a position where none of the configured fragment lengths and
# offsets fits inside the transcript.
#
# This happens close to either end of an mRNA, where a ribosome has no room for
# the footprints the distribution allows. Such reads go to the other positions
# of the same region of the same transcript, in proportion to the reads those
# positions already have, so the region keeps its read count. If no other
# position of the region has both reads and room (positions without reads were
# removed before this step), the region loses these reads; they are returned in
# `dropped`. Used when each position keeps its own count (its `read_count` if
# given, else its `score`); when reads are drawn per group (end bias),
# allocate_fragment_candidates() does this by itself.
move_unplaceable_reads <- function(signal_table, placeable) {
  column <- if ("read_count" %in% names(signal_table)) "read_count" else "score"
  count <- signal_table[[column]]
  stuck <- which(!placeable & count > 0)
  if (!length(stuck)) return(list(signal_table = signal_table, dropped = NULL))
  signal_table <- data.table::copy(signal_table)
  groups <- fragment_sampling_groups(signal_table)
  involved <- which(groups %in% groups[stuck])
  rows_per_group <- split(involved, factor(groups[involved],
                                           levels = unique(groups[stuck])))
  dropped <- list()
  for (rows in rows_per_group) {
    from <- rows[!placeable[rows] & count[rows] > 0]
    to <- rows[placeable[rows] & count[rows] > 0]
    reads <- sum(count[from])
    count[from] <- count[from] * 0L
    if (length(to)) {
      count[to] <- count[to] + as.vector(stats::rmultinom(1L, reads, count[to]))
    } else {
      dropped[[length(dropped) + 1L]] <- data.table::data.table(
        sampling_group = as.character(groups[from[1]]),
        transcript_id = as.character(signal_table$transcript_id[from[1]]),
        reads = as.numeric(reads)
      )
    }
  }
  data.table::set(signal_table, j = column, value = count)
  list(signal_table = signal_table, dropped = data.table::rbindlist(dropped))
}

# Opt-in faster joint draw: options(coverageSim.fragment_draw = "two_stage")
# (default "joint").
#
# The joint draw spreads a group's read budget over every (position, length,
# offset) candidate with weight w = geometry probability x bias weight x
# position score. A multinomial over all candidates has the same distribution
# as two nested ones: first the budget over positions, with each position's
# summed weight W_i, then each position's reads over its own candidates. Only
# the positions that receive reads need their candidates built. Without
# sequence-dependent biases (end-bias / codon-bias tables) a position's summed
# weight depends only on its frame, unless it sits so close to a transcript
# end that some candidates do not fit; only those edge positions are expanded
# to compute W_i. Same distribution as allocate_fragment_candidates(), but a
# different use of the random numbers, so outputs are not identical to the
# default path.
allocate_fragments_two_stage <- function(signal_table, transcript_models,
                                         candidates, geometry) {
  n_signal <- nrow(signal_table)
  labels <- fragment_sampling_groups(signal_table)
  budgets <- data.table::data.table(draw_group = labels, sampling_group = labels,
                                    transcript_id = signal_table$transcript_id)
  if ("read_count" %in% names(signal_table)) {
    budgets[, read_count := signal_table$read_count]
    budgets <- budgets[, .(sampling_group = sampling_group[1],
                           transcript_id = transcript_id[1],
                           total = if (data.table::uniqueN(read_count) == 1L) {
                             read_count[1]
                           } else NA_real_), by = draw_group]
  } else {
    budgets[, score := signal_table$score]
    budgets <- budgets[, .(sampling_group = sampling_group[1],
                           transcript_id = transcript_id[1],
                           total = round(sum(score))), by = draw_group]
  }
  if (anyNA(budgets$total) || any(budgets$total < 0) ||
      any(budgets$total != floor(budgets$total))) {
    stop("Invalid fragment sampling budget")
  }
  group_id <- match(labels, budgets$draw_group)

  # Transcript coordinate of every position and its transcript length.
  transcript_id <- signal_table$transcript_id
  used_ids <- unique(transcript_id)
  used_models <- lapply(used_ids, function(id) transcript_models[[id]])
  model_index <- match(transcript_id, used_ids)
  if (any(vapply(used_models, is.null, logical(1)))) {
    stop("No transcript model for: ", used_ids[vapply(used_models, is.null, logical(1))][1])
  }
  site_tx <- rep(NA_integer_, n_signal)
  rows_per_model <- split(seq_len(n_signal), model_index)
  for (k in as.integer(names(rows_per_model))) {
    rows <- rows_per_model[[as.character(k)]]
    site_tx[rows] <- genomic_sites_to_transcript(used_models[[k]],
                                                 signal_table$signal_position[rows])
  }
  if (anyNA(site_tx)) {
    stop("Signal position is not exonic in transcript ", transcript_id[which(is.na(site_tx))[1]])
  }
  tx_length <- vapply(used_models, function(m) as.numeric(m$length), numeric(1))[model_index]

  # Summed candidate weight per position (before the position score).
  prior <- candidates$probability / sum(candidates$probability)
  frame <- if ("region_position" %in% names(signal_table)) {
    as.integer((signal_table$region_position - 1L) %% 3L)
  } else integer(n_signal)
  per_frame <- vapply(0:2, function(f) {
    sum(prior * fragment_feature_weights(rep(f, nrow(candidates)),
                                         candidates$fragment_length,
                                         geometry$frame_bias, "frame"))
  }, numeric(1))
  interior <- site_tx >= 1L + max(candidates$site_offset) &
    site_tx <= tx_length - max(candidates$fragment_length - 1L - candidates$site_offset)
  site_weight <- per_frame[frame + 1L]
  edge <- which(!interior)
  if (length(edge)) {
    edge_candidates <- fragment_candidates(signal_table[edge], transcript_models,
                                           candidates, geometry)
    edge_sum <- edge_candidates[, .(w = sum(geometry_probability * bias_weight)),
                                by = signal_row]
    site_weight[edge] <- 0
    site_weight[edge[edge_sum$signal_row]] <- edge_sum$w
  }
  fits_any <- site_weight > 0
  weight <- site_weight * signal_table$score
  if (any(!is.finite(weight))) {
    stop("End-bias strength produces non-finite weights; reduce strength")
  }

  # Groups whose positions have no room for any fragment lose their reads.
  has_room <- unique(group_id[fits_any])
  dropped <- budgets[total > 0 & !seq_len(.N) %in% has_room,
                     .(reads = as.numeric(sum(total))),
                     by = .(sampling_group = as.character(sampling_group),
                            transcript_id = as.character(transcript_id))]

  # Stage 1: each group's budget over its positions, groups in order of first
  # appearance among positions that have room.
  position_reads <- integer(n_signal)
  group_total_weight <- numeric(nrow(budgets))
  room_rows <- which(fits_any)
  room_group <- group_id[room_rows]
  rows_by_group <- split(room_rows, match(room_group, unique(room_group)))
  for (rows in rows_by_group) {
    g <- group_id[rows[1L]]
    w <- weight[rows]
    if (sum(w) <= 0) stop("End-bias weights exclude every feasible fragment")
    group_total_weight[g] <- sum(w)
    position_reads[rows] <- as.integer(stats::rmultinom(1L, budgets$total[g], w / sum(w)))
  }

  # Stage 2: each position's reads over its own candidates.
  hit <- which(position_reads > 0L)
  allocated <- fragment_candidates(signal_table[hit], transcript_models,
                                   candidates, geometry)
  allocated[, signal_row := hit[signal_row]]
  allocated[, sampling_group := labels[signal_row]]
  w <- allocated$geometry_probability * allocated$bias_weight
  score <- integer(length(w))
  rows_by_position <- split(seq_along(w), allocated$signal_row)
  for (idx in rows_by_position) {
    i <- allocated$signal_row[idx[1L]]
    score[idx] <- as.integer(stats::rmultinom(1L, position_reads[i], w[idx] / sum(w[idx])))
  }
  # Final probability of each candidate within its group, as in the joint draw.
  probability <- w * signal_table$score[allocated$signal_row] /
    group_total_weight[group_id[allocated$signal_row]]
  data.table::set(allocated, j = "probability", value = probability)
  data.table::set(allocated, j = "score", value = score)
  data.table::setattr(allocated, "dropped_reads", dropped)
  allocated
}

# Turn signal positions into whole fragments with real coordinates.
#
# Up to here a read is just a position: where the ribosome was. A sequenced
# fragment is a stretch of RNA around that position, and how far it extends in
# each direction depends on its length and on where inside it the ribosome site
# falls. This picks a length and an offset for every read, checks the fragment
# fits inside its transcript, and hands back fragments ready to be written out.
#
# Transcript ends are the awkward case: near them, some length-and-offset
# combinations would run past the edge, and letting a fragment hang over the
# edge would invent sequence. `boundary_action` chooses: "error" stops, and
# "renormalize" keeps only the combinations that do fit and redistributes the
# probability among them, so positions near an end draw from a narrower set of
# geometries. Where none of the configured combinations fits (a ribosome too
# close to the end for the footprints the distribution allows), the position
# gets no reads: they go to other positions of the same region, or are left out
# if no other position of the region has both signal and room.
# Reads left out are listed in the "dropped_reads" attribute of the result.
make_simulated_rpf_fragments <- function(signal_table, transcript_models,
                                         fragment_lengths, fragment_geometry) {
  geometry <- normalize_fragment_geometry(fragment_geometry)
  candidates <- resolve_fragment_distribution(fragment_lengths, geometry)
  if (any(!is.finite(signal_table$score) | signal_table$score < 0) ||
      (!"read_count" %in% names(signal_table) &&
       any(signal_table$score != as.integer(signal_table$score)))) {
    stop("Simulated RPF scores must be non-negative integers")
  }

  joint <- has_active_fragment_bias(geometry)
  draw <- fragment_draw_option()
  if (draw == "two_stage" && joint && !needs_candidate_sequence(geometry) &&
      geometry$boundary_action != "error") {
    allocated <- allocate_fragments_two_stage(signal_table, transcript_models,
                                              candidates, geometry)
    return(finish_simulated_rpf_fragments(allocated, attr(allocated, "dropped_reads"),
                                          signal_table, transcript_models, geometry))
  }
  # Positions without reads cannot produce fragments. When each position keeps
  # its own count (no joint draw), their candidates are skipped: a multinomial
  # draw of size 0 uses no random numbers, so the result is unchanged.
  build_rows <- seq_len(nrow(signal_table))
  if (!joint && geometry$boundary_action != "error") {
    count_column <- if ("read_count" %in% names(signal_table)) "read_count" else "score"
    build_rows <- which(signal_table[[count_column]] > 0)
  }
  if (length(build_rows) == nrow(signal_table)) {
    candidate_table <- fragment_candidates(
      signal_table, transcript_models, candidates, geometry
    )
  } else {
    candidate_table <- fragment_candidates(
      signal_table[build_rows], transcript_models, candidates, geometry
    )
    candidate_table[, signal_row := build_rows[signal_row]]
  }
  dropped <- NULL
  if (!joint) {
    placeable <- tabulate(candidate_table$signal_row,
                          nbins = nrow(signal_table)) > 0L
    moved <- move_unplaceable_reads(signal_table, placeable)
    signal_table <- moved$signal_table
    dropped <- moved$dropped
  }
  allocated <- allocate_fragment_candidates(candidate_table, signal_table, joint)
  if (joint) dropped <- attr(allocated, "dropped_reads")
  finish_simulated_rpf_fragments(allocated, dropped, signal_table,
                                 transcript_models, geometry)
}

# Drawn candidates -> fragment table with genome coordinates (the second half of
# make_simulated_rpf_fragments()).
finish_simulated_rpf_fragments <- function(allocated, dropped, signal_table,
                                           transcript_models, geometry) {
  # Only the candidates that were drawn become fragments, often a small part of
  # all candidates. They stay in the order of the candidate table (signal row,
  # then candidate order), which sets the fragment ids.
  drawn <- allocated[score > 0L]
  transcript_id <- signal_table$transcript_id[drawn$signal_row]
  signal_position <- signal_table$signal_position[drawn$signal_row]
  # Each transcript's model is looked up once, the same way as above.
  used_ids <- unique(transcript_id)
  used_models <- lapply(seq_along(used_ids), function(k) {
    transcript_models[[used_ids[k]]]
  })
  model_index <- match(transcript_id, used_ids)
  # Sequence and codon of the drawn fragments, when fragment_candidates() left
  # them out (no sequence-dependent bias): the same substrings it would cut.
  if (nrow(drawn) && anyNA(drawn$sequence)) {
    tx_sequence <- vapply(used_models, function(model) model$sequence,
                          character(1))[model_index]
    site_tx <- drawn$tx_start + drawn$site_offset
    drawn[, `:=`(sequence = substr(tx_sequence, tx_start, tx_end),
                 codon = substr(tx_sequence, site_tx, site_tx + 2L))]
  }
  strand <- vapply(used_models, function(model) model$strand,
                   character(1))[model_index]
  chromosome <- vapply(used_models, function(model) {
    as.character(GenomicRanges::seqnames(model$exons)[1])
  }, character(1))[model_index]
  alignments <- fragment_alignments(used_models, model_index, drawn$tx_start,
                                    drawn$fragment_length)
  start <- alignments$pos
  end <- alignments$end
  # `sequence` is transcript-oriented (5'->3' of the mRNA), matching the
  # biological fragment/read direction; it is what codon/end-motif biases and
  # the ground-truth table use. SAM/BAM SEQ must instead be reported relative
  # to the reference (genome plus-strand), so for minus-strand fragments we
  # store the reverse complement separately rather than reusing `sequence`
  # when writing alignments.
  reference_sequence <- drawn$sequence
  minus <- strand != "+"
  reference_sequence[minus] <- as.character(Biostrings::reverseComplement(
    Biostrings::DNAStringSet(drawn$sequence[minus])
  ))
  result <- data.table::data.table(
    seqnames = chromosome,
    start = start,
    end = end,
    strand = strand,
    cigar = alignments$cigar,
    score = drawn$score,
    sampling_group = drawn$sampling_group,
    transcript_id = transcript_id,
    ribosome_site = signal_position,
    signal_position = signal_position,
    fragment_start = start,
    fragment_end = end,
    five_prime_end = ifelse(minus, end, start),
    three_prime_end = ifelse(minus, start, end),
    fragment_length = drawn$fragment_length,
    site_offset = drawn$site_offset,
    site_reference = geometry$site_reference,
    geometry_probability = drawn$geometry_probability,
    sequence = drawn$sequence,
    reference_sequence = reference_sequence,
    codon = drawn$codon,
    frame = drawn$frame,
    five_prime_kmer = if (!is.null(drawn$five_prime_kmer)) drawn$five_prime_kmer
                      else end_bias_kmer(drawn$sequence, geometry$five_prime_bias),
    three_prime_kmer = if (!is.null(drawn$three_prime_kmer)) drawn$three_prime_kmer
                       else end_bias_kmer(drawn$sequence, geometry$three_prime_bias, TRUE),
    end_bias_weight = drawn$end_bias_weight,
    codon_bias_weight = drawn$codon_bias_weight,
    frame_bias_weight = drawn$frame_bias_weight,
    fragment_bias_weight = drawn$bias_weight,
    final_probability = drawn$probability
  )
  if (has_added_base(geometry)) {
    transcript_sequence <- vapply(used_models, function(model) model$sequence,
                                  character(1))[model_index]
    exon_start <- lapply(used_models, function(model) model$cumulative_start)[model_index]
    result <- add_nontemplated_5prime_base(result, geometry$added_base, transcript_sequence,
                                           drawn$tx_start, exon_start)
  }
  result[, fragment_id := sprintf("fragment_%d", seq_len(.N))]
  if (is.null(dropped) || !nrow(dropped)) {
    dropped <- data.table::data.table(sampling_group = character(0),
                                      transcript_id = character(0),
                                      reads = numeric(0))
  }
  data.table::setattr(result, "dropped_reads", dropped)
  result
}

# The non-templated 5' base (fragment_geometry$added_base, opt-in).
#
# Reverse transcriptase (circularization libraries) or the template-switching oligo adds a base in front of the
# footprint 5' end of some reads. Each fragment's reads get it with probability `rate` (binomial split of the
# score). The read is then one base longer than the footprint. An aligner soft clips the added base when it does
# not match the genome base upstream of the footprint ("clipped": cigar 1S at the read 5' end), and aligns it when
# it does ("hidden": the alignment starts one base early and is one base longer, so the read looks like a footprint
# of length L + 1 whose 5' end starts with the added base). The upstream base is taken from the transcript; at an
# exon start (or the transcript start) the genome base upstream is intronic or unknown here, and the added base is
# then treated as clipped.
#
# Only the alignment changes (start, end, cigar, reference_sequence); the ground-truth columns (fragment_start,
# fragment_end, five_prime_end, fragment_length, sequence) keep describing the footprint. The column
# `added_base` says "none", "clipped" or "hidden" for every row. Rows that get the base are split off with their
# own score, right after the row they came from.
add_nontemplated_5prime_base <- function(result, added_base, transcript_sequence, tx_start, exon_start) {
  n <- nrow(result)
  added_n <- stats::rbinom(n, result$score, added_base$rate)
  base <- added_base$base
  upstream <- substr(transcript_sequence, tx_start - 1L, tx_start - 1L)
  at_exon_start <- vapply(seq_len(n), function(i) tx_start[i] %in% exon_start[[i]], logical(1))
  hidden <- upstream == base & !at_exon_start
  keep <- result[, `:=`(added_base = "none", source_row = seq_len(n))]
  keep[, score := score - added_n]
  add <- data.table::copy(keep[added_n > 0])
  if (nrow(add)) {
    add[, score := added_n[added_n > 0]]
    h <- hidden[added_n > 0]
    plus <- add$strand == "+"
    grow_first <- function(cigar) {
      m <- as.integer(sub("^([0-9]+)M.*$", "\\1", cigar))
      paste0(m + 1L, sub("^[0-9]+M", "M", cigar))
    }
    grow_last <- function(cigar) {
      m <- as.integer(sub("^(.*[A-Z])?([0-9]+)M$", "\\2", cigar))
      paste0(sub("[0-9]+M$", "", cigar), m + 1L, "M")
    }
    complement <- chartr("ACGT", "TGCA", base)
    add[plus, reference_sequence := paste0(base, reference_sequence)]
    add[!plus, reference_sequence := paste0(reference_sequence, complement)]
    add[h & plus, `:=`(start = start - 1L, cigar = grow_first(cigar), added_base = "hidden")]
    add[h & !plus, `:=`(end = end + 1L, cigar = grow_last(cigar), added_base = "hidden")]
    add[!h & plus, `:=`(cigar = paste0("1S", cigar), added_base = "clipped")]
    add[!h & !plus, `:=`(cigar = paste0(cigar, "1S"), added_base = "clipped")]
  }
  out <- data.table::rbindlist(list(keep[score > 0L], add))
  out <- out[order(source_row, added_base != "none")]
  out[, source_row := NULL]
  out[]
}

simulated_rpf_alignments <- function(fragment_table, seqinfo) {
  GenomicAlignments::GAlignments(
    seqnames = fragment_table$seqnames,
    pos = fragment_table$start,
    cigar = fragment_table$cigar,
    strand = fragment_table$strand,
    seqinfo = seqinfo,
    score = fragment_table$score,
    fragment_id = fragment_table$fragment_id,
    transcript_id = fragment_table$transcript_id,
    ribosome_site = fragment_table$ribosome_site,
    signal_position = fragment_table$signal_position,
    fragment_start = fragment_table$fragment_start,
    fragment_end = fragment_table$fragment_end,
    five_prime_end = fragment_table$five_prime_end,
    three_prime_end = fragment_table$three_prime_end,
    fragment_length = fragment_table$fragment_length,
    site_offset = fragment_table$site_offset,
    site_reference = fragment_table$site_reference,
    geometry_probability = fragment_table$geometry_probability,
    sequence = fragment_table$sequence,
    reference_sequence = fragment_table$reference_sequence,
    codon = fragment_table$codon,
    frame = fragment_table$frame,
    end_bias_weight = fragment_table$end_bias_weight,
    codon_bias_weight = fragment_table$codon_bias_weight,
    frame_bias_weight = fragment_table$frame_bias_weight,
    fragment_bias_weight = fragment_table$fragment_bias_weight
  )
}

# Compatibility aliases for scripts written before the simulated-RPF naming.
make_physical_fragments <- function(...) {
  warning("make_physical_fragments() is deprecated; use make_simulated_rpf_fragments()")
  make_simulated_rpf_fragments(...)
}

physical_fragment_alignments <- function(...) {
  warning("physical_fragment_alignments() is deprecated; use simulated_rpf_alignments()")
  simulated_rpf_alignments(...)
}

append_rnase_to_physical_table <- function(...) {
  warning(paste0(
    "append_rnase_to_physical_table() is deprecated; use ",
    "append_rnase_to_simulated_rpf_table()"
  ))
  append_rnase_to_simulated_rpf_table(...)
}

write_fragment_ground_truth <- function(fragment_table, file_base,
                                        ground_truth) {
  if (identical(ground_truth, FALSE) || is.null(ground_truth)) return(NULL)
  path <- if (identical(ground_truth, TRUE)) {
    paste0(file_base, "_ground_truth.tsv")
  } else {
    if (length(ground_truth) != 1L || !is.character(ground_truth)) {
      stop("ground_truth must be FALSE, TRUE, or a directory path")
    }
    dir.create(ground_truth, recursive = TRUE, showWarnings = FALSE)
    file.path(ground_truth, paste0(basename(file_base), "_ground_truth.tsv"))
  }
  columns <- c(
    "fragment_id", "transcript_id", "ribosome_site", "signal_position",
    "fragment_start", "fragment_end", "five_prime_end", "three_prime_end",
    "fragment_length", "site_offset", "site_reference", "geometry_probability",
    "strand", "cigar", "sequence", "score", "sampling_group",
    "codon", "frame", "five_prime_kmer", "three_prime_kmer",
    "end_bias_weight", "codon_bias_weight", "frame_bias_weight",
    "fragment_bias_weight", "final_probability"
  )
  if ("added_base" %in% names(fragment_table)) columns <- c(columns, "added_base")
  data.table::fwrite(fragment_table[, ..columns], path, sep = "\t")
  path
}

# Next to the library file: the reads of each transcript region that were left
# out because no position of the region with reads has room for a complete
# fragment (see make_simulated_rpf_fragments()). Written only when something
# was left out.
write_dropped_reads <- function(dropped, file_base) {
  path <- paste0(file_base, "_dropped_reads.tsv")
  if (!nrow(dropped)) {
    # Not left over from an earlier run into the same folder.
    if (file.exists(path)) file.remove(path)
    return(NULL)
  }
  group <- as.character(dropped$sampling_group)
  data.table::fwrite(data.table::data.table(
    transcript_id = dropped$transcript_id,
    region = ifelse(grepl(":", group, fixed = TRUE), sub(":.*", "", group), NA_character_),
    reads = dropped$reads,
    reason = ifelse(group == "uorf:no_uorf", "no uORF in the genome",
                    "no position with reads and room for a whole fragment")
  ), path, sep = "\t")
  path
}
