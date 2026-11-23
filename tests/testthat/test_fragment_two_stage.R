two_stage_fixture <- function() {
  genome_dir <- tempfile("two-stage-genome-"); dir.create(genome_dir)
  set.seed(31)
  g <- suppressWarnings(simGenome(n = 8, out_dir = genome_dir, genome_name = "two_stage",
                                  max_uorfs = 1, length_model = "human"))
  counts <- simCountTables(ORFik::loadRegion(g["txdb"], "cds"), libtypes = "RFP",
                           replicates = 1, conditions = "WT", print_statistics = FALSE,
                           plot_PCA = FALSE, interceptMean = 8)
  list(genome = g, regions = simCountTablesRegions(count_table = counts))
}

run_fragment_draw <- function(fixture, draw, seed = 5) {
  old <- options(coverageSim.codon_bias = codon_bias_preset("realistic"),
                 coverageSim.fragment_draw = draw)
  on.exit(options(old))
  truth <- tempfile("two-stage-truth-")
  set.seed(seed)
  do.call(simNGScoverage, c(list(fixture$genome, fixture$regions,
                                 exp_name = basename(tempfile("ts")),
                                 exp_save_dir = tempfile("two-stage-exp-"),
                                 libFormats = list(RFP = "ofst"), validate = FALSE,
                                 ground_truth = truth), realistic_rfp_settings()))
  data.table::fread(list.files(truth, pattern = "ground_truth", full.names = TRUE)[1])
}

test_that("two-stage fragment draw keeps region budgets and the target distribution", {
  fixture <- two_stage_fixture()
  joint <- run_fragment_draw(fixture, "joint")
  two <- run_fragment_draw(fixture, "two_stage")
  # Same region budgets.
  per_group <- function(x) x[, .(n = sum(score)), by = sampling_group][order(sampling_group)]
  expect_equal(per_group(two), per_group(joint))
  expect_identical(names(two), names(joint))
  expect_true(all(nchar(two$sequence) == two$fragment_length))
  expect_true(all(nchar(two$codon) == 3L))
  # Same target distribution: a fragment drawn by both has the same final
  # probability (its share of the region's total candidate weight).
  key <- c("transcript_id", "signal_position", "fragment_length", "site_offset")
  both <- merge(unique(joint[, c(key, "final_probability"), with = FALSE]),
                unique(two[, c(key, "final_probability"), with = FALSE]),
                by = key, suffixes = c("_joint", "_two"))
  expect_gt(nrow(both), 20)
  expect_equal(both$final_probability_two, both$final_probability_joint, tolerance = 1e-8)
})

test_that("two-stage draw matches the joint draw in repeated sampling", {
  skip_on_cran()
  fixture <- two_stage_fixture()
  per_position <- function(draw) {
    data.table::rbindlist(lapply(1:8, function(s) {
      run_fragment_draw(fixture, draw, seed = 100 + s)[, .(n = sum(score)), by = .(transcript_id, signal_position)]
    }))[, .(n = sum(n)), by = .(transcript_id, signal_position)]
  }
  m <- merge(per_position("joint"), per_position("two_stage"),
             by = c("transcript_id", "signal_position"), all = TRUE, suffixes = c("_joint", "_two"))
  m[is.na(n_joint), n_joint := 0]; m[is.na(n_two), n_two := 0]
  m <- m[(n_joint + n_two) >= 10]
  # Two-sample chi-square per position: about 1 per position if both draws
  # share one distribution.
  ratio <- sum((m$n_joint - m$n_two)^2 / (m$n_joint + m$n_two)) / nrow(m)
  expect_gt(nrow(m), 50)
  expect_lt(abs(ratio - 1), 0.25)
})

test_that("fragment_draw option is validated", {
  old <- options(coverageSim.fragment_draw = "fast"); on.exit(options(old))
  expect_error(fragment_draw_option(), "joint")
})
