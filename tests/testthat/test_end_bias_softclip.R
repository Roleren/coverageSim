# Soft clips and the hidden non-templated 5' base in learn_end_bias().

# Add a non-templated base in front of the 5' end of a share p of the fragments. Where the genome base
# upstream of the fragment is that base, the read aligns one nt longer (hidden); otherwise it is a 1S clip.
add_nontemplated_base <- function(fragments, genome, p, base = "T") {
  fragments <- data.table::copy(data.table::as.data.table(fragments))
  added_n <- stats::rbinom(nrow(fragments), fragments$score, p)
  added <- fragments[added_n > 0]
  added[, score := added_n[added_n > 0]]
  fragments[, score := score - added_n]
  fragments <- fragments[score > 0]
  plus <- added$strand == "+"
  genomic_end <- added$start + suppressWarnings(GenomicAlignments::cigarWidthAlongReferenceSpace(added$cigar)) - 1L
  upstream <- ifelse(plus, substring(genome, added$start - 1L, added$start - 1L),
                     chartr("ACGT", "TGCA", substring(genome, genomic_end + 1L, genomic_end + 1L)))
  hidden <- upstream == base
  grow_first <- function(cigar) vapply(cigar, function(z) {
    sub("^[0-9]+M", paste0(as.integer(sub("M.*$", "", z)) + 1L, "M"), z) }, character(1), USE.NAMES = FALSE)
  grow_last <- function(cigar) vapply(cigar, function(z) {
    n <- as.integer(sub("^(.*[A-Z])?([0-9]+)M$", "\\2", z))
    sub("[0-9]+M$", paste0(n + 1L, "M"), z) }, character(1), USE.NAMES = FALSE)
  added[plus, reference_sequence := paste0(base, reference_sequence)]
  added[!plus, reference_sequence := paste0(reference_sequence, chartr("ACGT", "TGCA", base))]
  added[hidden & plus, `:=`(start = start - 1L, cigar = grow_first(cigar))]
  added[hidden & !plus, cigar := grow_last(cigar)]
  added[!hidden & plus, cigar := paste0("1S", cigar)]
  added[!hidden & !plus, cigar := paste0(cigar, "1S")]
  out <- data.table::rbindlist(list(fragments, added), fill = TRUE)
  attr(out, "hidden") <- sum(added$score[hidden]); attr(out, "clipped") <- sum(added$score[!hidden])
  out
}

added_base_bam <- function(fixture, p, seed = 7L) {
  set.seed(seed)
  geometry <- fixture$geometry
  geometry$five_prime_bias <- list(source = "user", table =
    make_synthetic_end_bias(enriched_kmer = "G", enriched_weight = 4))
  geometry$three_prime_bias <- list(source = "user", table =
    make_synthetic_end_bias(enriched_kmer = "C", enriched_weight = 3))
  signal <- data.table::copy(fixture$signal)[, score := score * 4L]
  fragments <- make_simulated_rpf_fragments(signal, fixture$models, 28:29, geometry)
  genome <- as.character(Biostrings::readDNAStringSet(fixture$fasta)[[1]])
  reads <- add_nontemplated_base(fragments, genome, p)
  bam <- tempfile(fileext = ".bam")
  write_bam_library(simulated_rpf_alignments(reads, GenomeInfoDb::seqinfo(fixture$transcripts)), bam)
  list(bam = unname(bam[1]), hidden = attr(reads, "hidden"), clipped = attr(reads, "clipped"))
}

test_that("terminal soft clips are split from the aligned core", {
  clips <- terminal_soft_clips(c("28M", "1S28M", "28M2S", "3S12M60N16M1S"))
  expect_equal(clips$leading, c(0L, 1L, 0L, 3L))
  expect_equal(clips$trailing, c(0L, 0L, 2L, 1L))
  expect_equal(clips$core, c("28M", "28M", "28M", "12M60N16M"))
  expect_error(validate_added_base("U"), "added_base")
  expect_silent(validate_added_base(NULL))
})

test_that("a hidden added 5' T inflates T and added_base = 'T' removes it", {
  fixture <- make_end_learning_fixture()  # site offsets 15/16: hidden 29-nt reads land in frame
  sim <- added_base_bam(fixture, p = 0.5)
  expect_gt(sim$hidden, 0); expect_gt(sim$clipped, sim$hidden)
  drop <- fit_learning_fixture(fixture, sim$bam)
  fixed <- fit_learning_fixture(fixture, sim$bam, added_base = "T")
  auto <- fit_learning_fixture(fixture, sim$bam, added_base = "auto")
  t_drop <- relative_end_weight(drop$five_prime_bias, "T")
  t_fixed <- relative_end_weight(fixed$five_prime_bias, "T")
  # truth: T/A = 1, G/A = 4
  expect_gt(t_drop, 1.25)
  expect_lt(abs(log(t_fixed)), abs(log(t_drop)) / 2)
  expect_equal(relative_end_weight(fixed$five_prime_bias, "G"), 4, tolerance = .15)
  expect_equal(relative_end_weight(fixed$three_prime_bias, "C"), 3, tolerance = .15)
  expect_equal(auto$diagnostics$added_base$base, "T")
  expect_equal(relative_end_weight(auto$five_prime_bias, "T"), t_fixed, tolerance = 1e-6)
  # the clip report sees the 1-nt 5' clips, mostly T
  sc <- drop$diagnostics$soft_clips
  expect_equal(sc$ends[end == "5'" & clip == "1", reads], sim$clipped)
  expect_equal(fixed$diagnostics$added_base$trimmed_reads > 0, TRUE)
  expect_gt(fixed$diagnostics$added_base$hidden_reads_removed, 0)
  # truth for the estimate: hidden reads of length 29 that reached the fit
  expect_true(is.finite(fixed$diagnostics$added_base$upstream_share))
})

test_that("no correction without an added base, and drop/trim agree on clean data", {
  fixture <- make_end_learning_fixture()
  bam <- make_learning_bam(fixture)
  drop <- fit_learning_fixture(fixture, bam)
  trim <- fit_learning_fixture(fixture, bam, soft_clips = "trim")
  auto <- fit_learning_fixture(fixture, bam, added_base = "auto")
  expect_equal(drop$five_prime_bias$table$weight, trim$five_prime_bias$table$weight)
  expect_equal(auto$five_prime_bias$table$weight, drop$five_prime_bias$table$weight)
  expect_true(is.na(auto$diagnostics$added_base$base))
  expect_match(auto$diagnostics$added_base$reason, "auto")
  expect_equal(drop$diagnostics$soft_clips$unclipped_share, 1)
})

test_that("hidden reads off frame (equal site offsets) need no 5' correction", {
  fixture <- make_end_learning_fixture()
  fixture$geometry$distribution$site_offset <- c(15L, 15L)
  sim <- added_base_bam(fixture, p = 0.5)
  fixed <- fit_learning_fixture(fixture, sim$bam, added_base = "T")
  expect_false(any(fixed$diagnostics$added_base$per_length$in_frame))
  expect_equal(fixed$diagnostics$added_base$hidden_reads_removed, 0)
})
