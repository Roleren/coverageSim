# Cut-site (flank) term, per-method default tables and the simulated non-templated 5' base.

flank_geometry <- function(fixture, five = 4, three = 3, five_flank = NULL, three_flank = NULL,
                           added = NULL) {
  geometry <- fixture$geometry
  geometry$five_prime_bias <- list(source = "user", table =
    make_synthetic_end_bias(enriched_kmer = "G", enriched_weight = five))
  geometry$three_prime_bias <- list(source = "user", table =
    make_synthetic_end_bias(enriched_kmer = "C", enriched_weight = three))
  if (!is.null(five_flank)) geometry$five_prime_flank_bias <- list(source = "user", table = five_flank)
  if (!is.null(three_flank)) geometry$three_prime_flank_bias <- list(source = "user", table = three_flank)
  if (!is.null(added)) geometry$added_base <- added
  geometry
}

simulate_flank_bam <- function(fixture, geometry, scale = 4L, seed = 11L) {
  set.seed(seed)
  signal <- data.table::copy(fixture$signal)[, score := score * scale]
  fragments <- make_simulated_rpf_fragments(signal, fixture$models, 28:29, geometry)
  bam <- tempfile(fileext = ".bam")
  write_bam_library(simulated_rpf_alignments(fragments, GenomeInfoDb::seqinfo(fixture$transcripts)), bam)
  list(bam = unname(bam[1]), fragments = fragments)
}

test_that("flank k-mers are the bases just outside the fragment, neutral past the transcript", {
  tx <- c("ACGTACGTAC", "ACGTACGTAC")
  five <- list(table = data.table::data.table(kmer = c("A", "C", "G", "T"), weight = 1))
  expect_equal(end_flank_kmer(tx, c(3L, 1L), c(6L, 10L), five), c("C", ""))
  expect_equal(end_flank_kmer(tx, c(3L, 1L), c(6L, 10L), five, three_prime = TRUE), c("G", ""))
  two <- list(table = data.table::data.table(kmer = "AC", weight = 2))
  expect_equal(end_flank_kmer(tx[1], 4L, 6L, two), "CG")
  expect_equal(end_flank_kmer(tx[1], 4L, 6L, two, three_prime = TRUE), "GT")
  flank <- normalize_end_bias(list(source = "user", table = data.table::data.table(
    kmer = c("A", "C", "G", "T"), weight = c(1, 1, 1, 3))), "five_prime_flank_bias")
  none <- list(source = "none", table = NULL, strength = 1)
  w <- end_bias_weight(c("AAAA", "AAAA", "AAAA"), 4L, none, none, flank, NULL, c("T", "A", ""), NULL)
  expect_equal(w, c(3, 1, 1))
  g <- normalize_fragment_geometry(NULL)
  expect_equal(g$added_base, list(base = "T", rate = 0))
  expect_null(g$five_prime_flank_bias$table)
  expect_error(normalize_fragment_geometry(list(added_base = list(base = "U", rate = .1))), "base")
  expect_error(normalize_fragment_geometry(list(added_base = list(base = "T", rate = 2))), "rate")
  expect_error(normalize_fragment_geometry(list(added_base = list(prob = .1))), "Unknown")
})

test_that("the simulator weights the cut site and learn_end_bias(flank = TRUE) recovers both tables", {
  fixture <- make_end_learning_fixture()
  five_flank <- data.table::data.table(kmer = c("A", "C", "G", "T"), weight = c(1, 1, 1, 3))
  three_flank <- data.table::data.table(kmer = c("A", "C", "G", "T"), weight = c(2, 1, 1, 1))
  geometry <- flank_geometry(fixture, five_flank = five_flank, three_flank = three_flank)
  sim <- simulate_flank_bam(fixture, geometry)
  # the drawn fragments carry the flank preference: T upstream is enriched over the transcript composition
  f <- sim$fragments
  seqs <- vapply(fixture$models, `[[`, character(1), "sequence")[f$transcript_id]
  tx_start <- mapply(function(id, five_end) {
    model <- fixture$models[[id]]; genomic_site_to_transcript(model, five_end)
  }, f$transcript_id, f$five_prime_end)
  upstream <- substr(seqs, tx_start - 1L, tx_start - 1L)
  t_share <- sum(f$score[upstream == "T"]) / sum(f$score[nzchar(upstream)])
  expect_gt(t_share, 0.4)

  fit <- fit_learning_fixture(fixture, sim$bam, flank = TRUE)
  expect_equal(relative_end_weight(fit$five_prime_bias, "G"), 4, tolerance = .15)
  expect_equal(relative_end_weight(fit$three_prime_bias, "C"), 3, tolerance = .15)
  expect_equal(relative_end_weight(fit$five_prime_flank_bias, "T"), 3, tolerance = .15)
  expect_equal(relative_end_weight(fit$five_prime_flank_bias, "C"), 1, tolerance = .15)
  expect_equal(relative_end_weight(fit$three_prime_flank_bias, "C"), 0.5, tolerance = .15)
  expect_equal(relative_end_weight(fit$three_prime_flank_bias, "G", "T"), 1, tolerance = .15)
  expect_true(fit$provenance$flank)
  # the default (no flank term) is unchanged and returns no flank tables
  plain <- fit_learning_fixture(fixture, sim$bam)
  expect_null(plain$five_prime_flank_bias)
  expect_false(plain$provenance$flank)
  # by_length gives one flank table per length
  by_len <- fit_learning_fixture(fixture, sim$bam, flank = TRUE, by_length = TRUE)
  expect_setequal(unique(by_len$five_prime_flank_bias$table$fragment_length), c(28L, 29L))
  # the learned flank tables drive a new simulation directly
  geometry2 <- fixture$geometry
  geometry2[c("five_prime_bias", "three_prime_bias", "five_prime_flank_bias", "three_prime_flank_bias")] <-
    fit[c("five_prime_bias", "three_prime_bias", "five_prime_flank_bias", "three_prime_flank_bias")]
  expect_gt(nrow(make_simulated_rpf_fragments(fixture$signal, fixture$models, 28:29, geometry2)), 0)
})

test_that("the simulated added 5' base is clipped or hidden as the upstream base decides", {
  fixture <- make_end_learning_fixture()
  geometry <- flank_geometry(fixture, added = list(base = "T", rate = 0.5))
  set.seed(3)
  f <- make_simulated_rpf_fragments(fixture$signal, fixture$models, 28:29, geometry)
  expect_equal(sum(f$score), sum(fixture$signal$score))
  expect_setequal(unique(f$added_base), c("none", "clipped", "hidden"))
  expect_equal(sum(f$score[f$added_base != "none"]) / sum(f$score), 0.5, tolerance = .03)
  plus <- f$strand == "+"
  expect_true(all(grepl("^1S", f$cigar[f$added_base == "clipped" & plus])))
  expect_true(all(grepl("1S$", f$cigar[f$added_base == "clipped" & !plus])))
  read_width <- vapply(regmatches(f$cigar, gregexpr("[0-9]+[MS]", f$cigar)), function(x)
    sum(as.integer(sub("[MS]$", "", x))), numeric(1))
  expect_equal(read_width, nchar(f$reference_sequence))
  expect_true(all(read_width[f$added_base != "none"] == f$fragment_length[f$added_base != "none"] + 1L))
  # hidden: the transcript base upstream of the footprint is the added base; clipped: it is not
  seqs <- vapply(fixture$models, `[[`, character(1), "sequence")[f$transcript_id]
  tx_start <- mapply(function(id, five_end) genomic_site_to_transcript(fixture$models[[id]], five_end),
                     f$transcript_id, f$five_prime_end)
  upstream <- substr(seqs, tx_start - 1L, tx_start - 1L)
  expect_true(all(upstream[f$added_base == "hidden"] == "T"))
  expect_false(any(upstream[f$added_base == "clipped"] == "T" &
                     !mapply(function(id, s) s %in% fixture$models[[id]]$cumulative_start,
                             f$transcript_id[f$added_base == "clipped"], tx_start[f$added_base == "clipped"])))
  # the hidden alignment starts one base before the footprint at the 5' end
  h <- f$added_base == "hidden"
  expect_true(all(f$start[h & plus] == f$fragment_start[h & plus] - 1L))
  expect_true(all(f$end[h & !plus] == f$fragment_end[h & !plus] + 1L))
  # off by default: no column, same draws
  set.seed(3)
  plain <- make_simulated_rpf_fragments(fixture$signal, fixture$models, 28:29, flank_geometry(fixture))
  expect_false("added_base" %in% names(plain))
  # multimapping artefacts do not support it yet
  expect_error(simulate_alignment_artifacts(f, list(multimapping_rate = .1)), "added_base")
})

test_that("truth check: learn_end_bias(added_base = 'auto') recovers the simulated truth", {
  fixture <- make_end_learning_fixture()  # offsets 15/16: hidden 29-nt reads land in frame
  geometry <- flank_geometry(fixture, added = list(base = "T", rate = 0.4))
  sim <- simulate_flank_bam(fixture, geometry)
  expect_gt(sum(sim$fragments[added_base == "hidden", score]), 0)
  drop <- fit_learning_fixture(fixture, sim$bam)
  auto <- fit_learning_fixture(fixture, sim$bam, added_base = "auto")
  expect_equal(auto$diagnostics$added_base$base, "T")
  # truth: 5' T/A = 1, G/A = 4, 3' C/A = 3
  expect_gt(relative_end_weight(drop$five_prime_bias, "T"), 1.2)
  expect_equal(relative_end_weight(auto$five_prime_bias, "T"), 1, tolerance = .1)
  expect_equal(relative_end_weight(auto$five_prime_bias, "G"), 4, tolerance = .15)
  expect_equal(relative_end_weight(auto$three_prime_bias, "C"), 3, tolerance = .15)
  # the estimated added rate matches the simulated one
  expect_equal(auto$diagnostics$added_base$added_rate, 0.4, tolerance = .15)
  # with the flank term the added base is not learned as a depleted cut site (truth: flank neutral)
  auto_flank <- fit_learning_fixture(fixture, sim$bam, added_base = "auto", flank = TRUE)
  expect_equal(relative_end_weight(auto_flank$five_prime_flank_bias, "T"), 1, tolerance = .15)
  expect_equal(relative_end_weight(auto_flank$five_prime_bias, "T"), 1, tolerance = .1)
  expect_gt(auto_flank$diagnostics$added_base$upstream_reads_restored, 0)
})

test_that("method_end_bias() gives chemistry, nuclease and added-base defaults", {
  circ <- method_end_bias("circularization", added_base = TRUE)
  expect_named(circ, c("five_prime_bias", "three_prime_bias", "five_prime_flank_bias",
                       "three_prime_flank_bias", "added_base"))
  for (p in circ[1:4]) {
    expect_setequal(p$table$kmer, c("A", "C", "G", "T"))
    expect_equal(exp(mean(log(p$table$weight))), 1, tolerance = .01)
  }
  expect_equal(circ$added_base$base, "T")
  expect_gt(circ$added_base$rate, 0.05)
  expect_lt(circ$five_prime_bias$table[kmer == "G", weight], 1)   # circularization 5' G depleted
  t1 <- method_end_bias(NULL, nuclease = "RNase T1")
  expect_named(t1, c("five_prime_flank_bias", "three_prime_flank_bias"))
  expect_gt(t1$five_prime_flank_bias$table[kmer == "G", weight], 4)   # T1 cuts 3' of G
  expect_equal(method_end_bias("dual_ligation", added_base = TRUE)$added_base$rate, 0)
  expect_error(method_end_bias("circularization", nuclease = "DNase"), "nuclease")
  expect_true(all(attr(circ, "evidence")$studies >= 1))
  expect_equal(method_end_bias("template_switching", strength = 0)$five_prime_bias$strength, 0)
  # usable as fragment_geometry
  fixture <- make_end_learning_fixture()
  geometry <- utils::modifyList(fixture$geometry, method_end_bias("circularization", added_base = TRUE))
  set.seed(1)
  f <- make_simulated_rpf_fragments(fixture$signal, fixture$models, 28:29, geometry)
  expect_equal(sum(f$score), sum(fixture$signal$score))
  expect_true(any(f$added_base == "clipped"))
})

test_that("footprints starting an internal exon are not restored by the added-base correction", {
  counted <- data.table::data.table(transcript_id = "tx", fragment_length = 28L, site_tx = 1:4,
    five = c("A", "A", "C", "C"), upstream = c("T", "T", "G", "G"),
    exon_start = c(FALSE, TRUE, FALSE, FALSE), count = c(60, 60, 100, 100))
  trimmed <- data.table::copy(counted)[, count := c(0, 0, 40, 0)]
  dist <- data.table::data.table(fragment_length = 28L, site_offset = 15L, probability = 1)
  out <- hidden_added_base_correction(counted, trimmed, dist, "T")
  # p = 40 / (40 + 200) for this length: the exon-start row keeps its count, the other T-upstream row is restored
  expect_equal(out$data[site_tx == 2, count], 60)
  expect_equal(out$data[site_tx == 1, count], 60 / (1 - 40 / 240))
  expect_equal(out$diagnostics$upstream_reads_restored, 60 / (1 - 40 / 240) - 60)
})
