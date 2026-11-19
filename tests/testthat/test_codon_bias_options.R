test_that("codon bias defaults reproduce the original and options validate", {
  withr::local_options(coverageSim.codon_bias = NULL)
  cb <- codon_bias_options()
  expect_equal(cb$sharpen, "per_gene"); expect_equal(cb$pause_sd, 0); expect_false(cb$per_genome_draws)
  withr::local_options(coverageSim.codon_bias = list(oops = 1))
  expect_error(codon_bias_options(), "unknown")
})

test_that("t57 preset: identical alphas across libraries, pauses on, session RNG untouched", {
  alpha <- list(g1 = c(0.03, stats::runif(80, 0.5, 2)), g2 = stats::runif(60, 0.5, 2))
  withr::local_options(coverageSim.codon_bias = codon_bias_preset("t57"))
  set.seed(5); a <- sim_sequence_bias(NULL, lengths(alpha) * 3, alpha, shapes(9)); r1 <- runif(1)
  set.seed(6); b <- sim_sequence_bias(NULL, lengths(alpha) * 3, alpha, shapes(9))
  expect_equal(a, b)
  set.seed(5); sim_sequence_bias(NULL, lengths(alpha) * 3, alpha, shapes(9)); expect_equal(runif(1), r1)
  withr::local_options(coverageSim.codon_bias = modifyList(codon_bias_preset("t57"), list(pause_sd = 0)))
  c0 <- sim_sequence_bias(NULL, lengths(alpha) * 3, alpha, shapes(9))
  expect_false(isTRUE(all.equal(a, c0)))
})
