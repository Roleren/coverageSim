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

test_that("zero_prob masks the same codon sites in all libraries, at the requested fraction; default unchanged", {
  alpha <- list(g1 = stats::runif(3000, 0.5, 2), g2 = stats::runif(3000, 0.5, 2))
  run <- function(zp, seed) {
    withr::local_options(coverageSim.codon_bias = modifyList(codon_bias_preset("realistic"), list(zero_prob = zp)))
    set.seed(seed); sim_sequence_bias(NULL, lengths(alpha) * 3, alpha, shapes(9))
  }
  a <- run(0.3, 5); b <- run(0.3, 6)
  expect_equal(a, b)                                  # identical across libraries (different session seeds)
  site <- function(o) { x <- o[[1]]; x <- x[4:(length(x) - 3)]; colSums(matrix(x[seq_len(3 * (length(x) %/% 3))], 3)) < 1e-10 }  # codon sites (margin 3 nt)
  expect_equal(mean(site(a)), 0.3, tolerance = 0.05)
  expect_equal(mean(site(run(0.5, 5))), 0.5, tolerance = 0.05)
  # default 0 is identical to the option not being set at all
  withr::local_options(coverageSim.codon_bias = codon_bias_preset("realistic"))
  set.seed(7); u <- sim_sequence_bias(NULL, lengths(alpha) * 3, alpha, shapes(9))
  expect_equal(u, run(0, 7))
  expect_equal(codon_bias_options()$zero_prob, 0)
  withr::local_options(coverageSim.codon_bias = list(zero_prob = 1))
  expect_error(codon_bias_options())
})
