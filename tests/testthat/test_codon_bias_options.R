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

test_that("tau2: replicate-shared window variance, per region, recovered by the E7/A cross-replicate estimator; default 0 identical", {
  alpha <- stats::setNames(lapply(1:60, function(i) stats::runif(300, 0.5, 2)), paste0("g", 1:60))
  nt <- lengths(alpha) * 3
  run <- function(tau2, seed, region = "cds", nu = alpha) {
    withr::local_options(coverageSim.codon_bias = modifyList(codon_bias_preset("realistic"), list(tau2 = tau2)))
    set.seed(seed); sim_sequence_bias(NULL, nt, nu, shapes(9), region = region)
  }
  # default (0 / not set) is identical to before, and consumes no session random numbers
  withr::local_options(coverageSim.codon_bias = codon_bias_preset("realistic"))
  set.seed(7); u <- sim_sequence_bias(NULL, nt, alpha, shapes(9)); r1 <- stats::runif(1)
  expect_equal(u, run(0, 7))
  expect_equal(codon_bias_options()$tau2, c(cds = 0, utr5 = 0, utr3 = 0))
  # identical across libraries (different session seeds), different from tau2 = 0
  a <- run(0.1, 5); b <- run(0.1, 6)
  expect_equal(a, b); expect_false(isTRUE(all.equal(a, u)))
  # per-region scale: cds-only tau2 leaves the UTR regions untouched; scalar applies to all
  expect_equal(run(c(cds = 0.1), 5, region = "leader"), run(0, 5, region = "leader"))
  expect_false(isTRUE(all.equal(run(c(utr5 = 0.1), 5, region = "leader"), run(0, 5, region = "leader"))))
  expect_equal(normalize_tau2(0.3), c(cds = 0.3, utr5 = 0.3, utr3 = 0.3))
  withr::local_options(coverageSim.codon_bias = list(tau2 = c(cds = 0.1, foo = 1))); expect_error(codon_bias_options(), "unknown tau2")
  # block factors: window sums (90 nt, aligned to the block start) have mean 1 and relative variance tau2
  tau <- 0.2; w0 <- run(0, 5); w1 <- run(tau, 5)
  rel <- unlist(lapply(seq_along(w0), function(i) {
    m <- (length(w0[[i]]) - nt[i]) %/% 2L; idx <- m + seq_len(nt[i] %/% 90L * 90L)
    s0 <- colSums(matrix(w0[[i]][idx], 90)); s1 <- colSums(matrix(w1[[i]][idx], 90)); s1 / s0 * sum(s0) / sum(s1)
  }))
  expect_equal(stats::var(rel), tau, tolerance = 0.2)
})

test_that("tau2 estimator (E7/A: replicate cross-product) recovers the configured tau2 within 20 %", {
  alpha <- stats::setNames(lapply(1:150, function(i) stats::runif(450, 0.5, 2)), paste0("g", 1:150))
  nt <- lengths(alpha) * 3
  w <- function(tau2) { withr::local_options(coverageSim.codon_bias = modifyList(codon_bias_preset("realistic"), list(tau2 = tau2)))
    set.seed(1); sim_sequence_bias(NULL, nt, alpha, shapes(9)) }
  w0 <- w(0); w1 <- w(0.1)
  win <- function(v, i) { m <- (length(v) - nt[i]) %/% 2L; colSums(matrix(v[m + seq_len(nt[i] %/% 90L * 90L)], 90)) }
  set.seed(11); depth <- 60000
  d <- data.table::rbindlist(lapply(seq_along(w1), function(i) {
    s0 <- win(w0[[i]], i); s1 <- win(w1[[i]], i); p <- w1[[i]] / sum(w1[[i]])
    A <- win(stats::rmultinom(1, depth, p)[, 1], i); B <- win(stats::rmultinom(1, depth, p)[, 1], i)
    E <- s0 / sum(s0) * depth
    data.table::data.table(gene = i, E = E, OA = A, OB = B)
  }))
  d[, `:=`(EA = E, EB = E)]
  # E7/A estimator, method of moments: tau2 = sum c_g (OA - EA)(OB - EB) / sum EA EB, c_g = nw / (nw - 1)
  d[, cg := .N / (.N - 1), by = gene]
  est <- d[, sum(cg * (OA - EA) * (OB - EB)) / sum(EA * EB)]
  expect_equal(est, 0.1, tolerance = 0.2)
  # zero tau2: the same estimator is ~0 (here: Poisson/half-specific noise cancels)
  d0 <- data.table::rbindlist(lapply(seq_along(w0), function(i) {
    s0 <- win(w0[[i]], i); p <- w0[[i]] / sum(w0[[i]])
    data.table::data.table(gene = i, E = s0 / sum(s0) * depth, OA = win(stats::rmultinom(1, depth, p)[, 1], i),
                           OB = win(stats::rmultinom(1, depth, p)[, 1], i))
  }))
  d0[, cg := .N / (.N - 1), by = gene]
  expect_lt(abs(d0[, sum(cg * (OA - E) * (OB - E)) / sum(E * E)]), 0.02)
})
