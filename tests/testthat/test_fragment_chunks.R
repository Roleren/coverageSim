# make_simulated_rpf_fragments() builds and draws candidates in chunks of whole draw groups
# (options(coverageSim.fragment_chunk_rows)); the result must equal one table for all positions (Inf).

chunk_fixture <- function() {
  set.seed(11)
  chromosome <- paste(sample(c("A", "C", "G", "T"), 3000L, TRUE), collapse = "")
  fasta <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(chr1 = chromosome)), filepath = fasta)
  Rsamtools::indexFa(fasta)
  transcripts <- GenomicRanges::GRangesList(
    tx_plus = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(101L, 301L), width = c(120L, 150L)), "+",
                                     exon_rank = 1:2),
    tx_minus = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1401L, 1101L), width = c(100L, 130L)), "-",
                                      exon_rank = 1:2))
  GenomeInfoDb::seqlengths(transcripts) <- c(chr1 = 3000L)
  models <- transcript_models(transcripts, fasta)
  # Positions across both transcripts, including the first and last exonic bases (no room for a
  # fragment), with interleaved regions on tx_plus and a region that has no room at all.
  plus <- c(101L, 102L, 110:140, 301L, 330:360, 450L)
  minus <- c(1500L, 1499L, 1450:1470, 1200:1215, 1101L)
  signals <- data.table::data.table(
    seqnames = "chr1", strand = rep(c("+", "-"), c(length(plus), length(minus))),
    transcript_id = rep(c("tx_plus", "tx_minus"), c(length(plus), length(minus))),
    signal_position = c(plus, minus))
  signals[, `:=`(start = signal_position, end = signal_position)]
  signals[, sampling_group := paste0(transcript_id, ":", c("a", "b")[(seq_len(.N) %% 2L) + 1L]), by = transcript_id]
  signals[signal_position %in% c(101L, 1101L), sampling_group := paste0(transcript_id, ":edge")]
  signals[, region_position := seq_len(.N), by = sampling_group]
  signals[, score := sample(1:6, .N, TRUE)]
  signals[signal_position %in% c(101L, 1101L), score := 4L]
  list(models = models, signals = signals)
}

chunk_geometry <- function(bias = TRUE, boundary = "renormalize") {
  list(source = "default", site_reference = "a_site", distribution = NULL, boundary_action = boundary,
       five_prime_bias = if (bias) list(source = "user", table = make_synthetic_end_bias(
         k = 2L, enriched_kmer = "GC", enriched_weight = 3)) else list(source = "none"),
       three_prime_bias = list(source = "none"), codon_bias = list(source = "none"),
       frame_bias = list(source = "none"))
}

run_chunked <- function(chunk_rows, signals, models, geometry, seed = 4) {
  old <- options(coverageSim.fragment_chunk_rows = chunk_rows); on.exit(options(old))
  set.seed(seed)
  result <- make_simulated_rpf_fragments(signals, models, 26:32, geometry)
  list(result = result, seed = .Random.seed)
}

test_that("chunked fragment draw equals one candidate table (joint and per position)", {
  fx <- chunk_fixture()
  setups <- list(
    joint = list(signals = fx$signals, geometry = chunk_geometry(TRUE)),
    joint_read_count = list(signals = data.table::copy(fx$signals)[, read_count := sum(score), by = sampling_group],
                            geometry = chunk_geometry(TRUE)),
    per_position = list(signals = data.table::copy(fx$signals)[seq(1L, .N, 4L), score := 0L],
                        geometry = chunk_geometry(FALSE)),
    per_position_read_count = list(signals = data.table::copy(fx$signals)[, read_count := score + 1L],
                                   geometry = chunk_geometry(FALSE)))
  for (nm in names(setups)) {
    s <- setups[[nm]]
    ref <- run_chunked(Inf, s$signals, fx$models, s$geometry)
    expect_gt(nrow(ref$result), 10)
    for (chunk in c(1, 40, 700, 2e6)) {
      got <- run_chunked(chunk, s$signals, fx$models, s$geometry)
      expect_identical(got$result, ref$result, info = paste(nm, chunk))
      expect_identical(attr(got$result, "dropped_reads"), attr(ref$result, "dropped_reads"), info = paste(nm, chunk))
      expect_identical(got$seed, ref$seed, info = paste(nm, chunk))
    }
  }
  # The edge-only regions have no room: their reads are reported as dropped.
  joint <- run_chunked(40, fx$signals, fx$models, chunk_geometry(TRUE))$result
  expect_setequal(attr(joint, "dropped_reads")$sampling_group, c("tx_plus:edge", "tx_minus:edge"))
})

test_that("chunked fragment draw refuses the same positions with the same message", {
  fx <- chunk_fixture()
  for (bias in c(TRUE, FALSE)) {
    geometry <- chunk_geometry(bias, boundary = "error")
    msg <- function(chunk) {
      old <- options(coverageSim.fragment_chunk_rows = chunk); on.exit(options(old))
      tryCatch(make_simulated_rpf_fragments(fx$signals, fx$models, 26:32, geometry),
               error = conditionMessage)
    }
    expect_match(msg(Inf), "Cannot create a complete fragment")
    expect_identical(msg(40), msg(Inf))
  }
  missing <- data.table::copy(fx$signals)[5L, transcript_id := "tx_none"]
  for (chunk in c(Inf, 40)) {
    old <- options(coverageSim.fragment_chunk_rows = chunk)
    expect_error(make_simulated_rpf_fragments(missing, fx$models, 26:32, chunk_geometry(TRUE)),
                 "No transcript model for: tx_none")
    options(old)
  }
  expect_error(withr::with_options(list(coverageSim.fragment_chunk_rows = 0), fragment_chunk_option()),
               "must be a number")
})

test_that("fragment_row_chunks never splits a draw group", {
  group_id <- c(1L, 1L, 2L, 1L, 3L, 3L, 4L, 2L, 5L, 6L, 6L)
  for (target in 1:12) {
    chunks <- fragment_row_chunks(group_id, target)
    expect_identical(unlist(chunks), seq_along(group_id))
    owner <- rep(seq_along(chunks), lengths(chunks))
    expect_true(all(tapply(owner, group_id, function(x) length(unique(x))) == 1L))
  }
  expect_identical(lengths(fragment_row_chunks(group_id, 1L)), c(8L, 1L, 2L))
  expect_identical(fragment_row_chunks(integer(0), 5L), list())
})
