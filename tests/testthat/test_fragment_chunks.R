# make_simulated_rpf_fragments() builds and draws candidates in chunks of whole draw groups
# (options(coverageSim.fragment_chunk_rows)), and without sequence biases draws in compiled code
# (options(coverageSim.fragment_compiled_draw)); the result must equal one table for all positions (Inf).

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

chunk_geometry <- function(bias = TRUE, boundary = "renormalize", frame = FALSE) {
  list(source = "default", site_reference = "a_site", distribution = NULL, boundary_action = boundary,
       five_prime_bias = if (bias) list(source = "user", table = make_synthetic_end_bias(
         k = 2L, enriched_kmer = "GC", enriched_weight = 3)) else list(source = "none"),
       three_prime_bias = list(source = "none"), codon_bias = list(source = "none"),
       frame_bias = if (frame) list(source = "user", strength = 0.25, table = data.table::data.table(
         frame = rep(0:2, 4), fragment_length = rep(27:30, each = 3),
         weight = c(2, 0.5, 0.4, 1.8, 0.6, 0.3, 1.5, 0.7, 0.9, 3, 0.2, 0.1))) else list(source = "none"))
}

run_chunked <- function(chunk_rows, signals, models, geometry, seed = 4, compiled = TRUE) {
  old <- options(coverageSim.fragment_chunk_rows = chunk_rows,
                 coverageSim.fragment_compiled_draw = compiled); on.exit(options(old))
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
                                   geometry = chunk_geometry(FALSE)),
    joint_frame = list(signals = fx$signals, geometry = chunk_geometry(FALSE, frame = TRUE)),
    joint_frame_read_count = list(signals = data.table::copy(fx$signals)[, read_count := sum(score), by = sampling_group],
                                  geometry = chunk_geometry(FALSE, frame = TRUE)),
    joint_frame_no_region_position = list(signals = data.table::copy(fx$signals)[, region_position := NULL],
                                          geometry = chunk_geometry(FALSE, frame = TRUE)))
  for (nm in names(setups)) {
    s <- setups[[nm]]
    ref <- run_chunked(Inf, s$signals, fx$models, s$geometry)
    expect_gt(nrow(ref$result), 10)
    for (compiled in c(FALSE, TRUE)) for (chunk in c(1, 40, 700, 2e6)) {
      got <- run_chunked(chunk, s$signals, fx$models, s$geometry, compiled = compiled)
      info <- paste(nm, chunk, if (compiled) "compiled" else "R")
      expect_identical(got$result, ref$result, info = info)
      expect_identical(attr(got$result, "dropped_reads"), attr(ref$result, "dropped_reads"), info = info)
      expect_identical(got$seed, ref$seed, info = info)
    }
  }
  # The edge-only regions have no room: their reads are reported as dropped.
  joint <- run_chunked(40, fx$signals, fx$models, chunk_geometry(TRUE))$result
  expect_setequal(attr(joint, "dropped_reads")$sampling_group, c("tx_plus:edge", "tx_minus:edge"))
})

test_that("chunked fragment draw refuses the same positions with the same message", {
  fx <- chunk_fixture()
  for (bias in c(TRUE, FALSE)) for (frame in c(FALSE, TRUE)) {
    geometry <- chunk_geometry(bias, boundary = "error", frame = frame)
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
  expect_error(withr::with_options(list(coverageSim.fragment_compiled_draw = NA), fragment_cpp_option()),
               "must be TRUE or FALSE")
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

test_that("compiled draw stops like the R draw on impossible weights and keeps empty inputs", {
  fx <- chunk_fixture()
  geometry <- chunk_geometry(FALSE, frame = TRUE)
  zero <- data.table::copy(fx$signals)[, read_count := sum(score), by = sampling_group][, score := 0L]
  for (compiled in c(FALSE, TRUE)) {
    expect_error(run_chunked(2e6, zero, fx$models, geometry, compiled = compiled),
                 "exclude every feasible fragment")
  }
  # A missing region position (NA frame) gives a non-finite weight on every path.
  na_frame <- data.table::copy(fx$signals)[c(3L, 50L), region_position := NA_integer_]
  for (compiled in c(FALSE, TRUE)) for (chunk in c(Inf, 40)) {
    expect_error(run_chunked(chunk, na_frame, fx$models, geometry, compiled = compiled), "non-finite weights")
  }
  # Only edge positions: no candidates, every read dropped, no random numbers used.
  edge <- fx$signals[signal_position %in% c(101L, 1101L)]
  ref <- run_chunked(Inf, edge, fx$models, geometry)
  got <- run_chunked(2e6, edge, fx$models, geometry)
  expect_equal(nrow(got$result), 0L)
  expect_identical(got$result, ref$result)
  expect_identical(attr(got$result, "dropped_reads"), attr(ref$result, "dropped_reads"))
  expect_identical(got$seed, ref$seed)
})

test_that("fragment_fits_cpp matches the R fit rule", {
  site <- c(1L, 5L, 20L, NA, 40L, 3L)
  len <- c(30, 30, 30, 30, NA, 3)
  pairs <- data.table::data.table(fragment_length = c(4L, 6L), site_offset = c(1L, 3L))
  r <- vapply(seq_along(site), function(i) {
    ok <- !is.na(site[i]) & site[i] - pairs$site_offset >= 1L &
      site[i] - pairs$site_offset + pairs$fragment_length - 1L <= len[i]
    isTRUE(any(ok))
  }, logical(1))
  expect_identical(fragment_fits_cpp(site, len, pairs$fragment_length, pairs$site_offset), r)
})
