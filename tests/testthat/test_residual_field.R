# Port S14: residual log-field, P-site smear, per-gene contrast exponent, study-shared field, gene dispersion
# heterogeneity and footprint-GC bias (claude-megafst explainability loop, Fable review 2026-11-22, NC-1, M22).

s14_genes <- function(n, len = 400, seed = 3) {
  tab <- load_seq_bias(type = "codon", shift = "a-site")
  sense <- setdiff(tab$seqs, c("###", "%%%", "&&&", "***", "TAA", "TAG", "TGA"))
  withr::local_seed(seed)
  genes <- lapply(seq_len(n), function(i) sample(sense, len, TRUE))
  list(genes = genes, alpha = lapply(genes, function(g) tab$alpha[match(g, tab$seqs)]))
}

test_that("S14 options are off by default, validate and leave the simulation unchanged", {
  withr::local_options(coverageSim.codon_bias = NULL)
  o <- codon_bias_options()
  expect_false(o$field$active); expect_true(o$smooth); expect_equal(o$sharpen, "per_gene")
  expect_false(s13_needs_codons(o))
  expect_error(residual_field_options(field_var = c(0.1, 0.1)))            # length differs from field_scales
  expect_error(residual_field_options(shared = 1.5))
  expect_error(residual_field_options(smear = c(0, 1)))
  expect_error(residual_field_options(gc_knots = list(gc = c(0.6, 0.4), log = c(0, 0))))
  withr::local_options(coverageSim.codon_bias = list(sharpen = "contrast", contrast_sd = -1)); expect_error(codon_bias_options())
  withr::local_options(coverageSim.codon_bias = list(field = list(bogus = 1))); expect_error(codon_bias_options())
  # a residual_field_options() result (with its 'active' flag) is accepted as is
  withr::local_options(coverageSim.codon_bias = list(field = residual_field_options(gc_bias = TRUE)))
  expect_true(codon_bias_options()$field$active); expect_true(s13_needs_codons())
  # an empty / all-off field list gives the same weights as no option at all
  g <- s14_genes(5, 200)
  withr::local_options(coverageSim.codon_bias = codon_bias_preset("realistic"))
  ref <- sim_sequence_bias(NULL, 600, g$alpha, shapes(9))
  withr::local_options(coverageSim.codon_bias = c(codon_bias_preset("realistic"), list(field = list())))
  expect_equal(sim_sequence_bias(NULL, 600, g$alpha, shapes(9), codons = g$genes, a_shift = 0L, study = "x"), ref)
})

test_that("AR(1) field components have the requested variance and correlation exp(-1 / L)", {
  withr::local_seed(1)
  x <- field_ar1(2e5, 3, 0.2)
  expect_equal(var(x), 0.2, tolerance = 0.03)
  expect_equal(stats::cor(x[-1], x[-length(x)]), exp(-1 / 3), tolerance = 0.01)
  expect_equal(stats::cor(x[-(1:3)], x[seq_len(length(x) - 3)]), exp(-1), tolerance = 0.05)
  y <- field_ar1(2e5, 300, 0.1)
  expect_equal(stats::cor(y[-(1:30)], y[seq_len(length(y) - 30)]), exp(-0.1), tolerance = 0.03)
  expect_equal(field_ar1(10, 3, 0), numeric(10))
})

test_that("the field is shared between studies by the requested fraction; gene amplitudes are lognormal", {
  f <- residual_field_options(field_var = c(0.1, 0.05, 0.05), shared = 0.7, gene_sd = 0)
  n <- rep(500, 400)
  a <- unlist(draw_residual_field(n, f, "A")); b <- unlist(draw_residual_field(n, f, "B"))
  expect_equal(unlist(draw_residual_field(n, f, "A")), a)                   # same study, same field
  expect_equal(stats::cor(a, b), 0.7, tolerance = 0.03)                      # two studies share 70%
  expect_equal(var(a), 0.2, tolerance = 0.1)
  # fully shared: every study gets the same field
  f1 <- residual_field_options(field_var = c(0.1, 0, 0), shared = 1)
  expect_equal(draw_residual_field(n[1:3], f1, "A"), draw_residual_field(n[1:3], f1, "B"))
  # gene_sd: per-gene amplitude s_g = exp(gene_sd z - gene_sd^2), so E[s^2] = 1 and sd(log s) = gene_sd
  fg <- residual_field_options(field_var = c(0.1, 0, 0), shared = 1, gene_sd = 0.5)
  x <- draw_residual_field(rep(400, 2000), fg, "A")
  s <- vapply(x, stats::sd, numeric(1)) / sqrt(0.1)
  expect_equal(stats::sd(log(s)), 0.5, tolerance = 0.06)
  expect_equal(mean(s^2), 1, tolerance = 0.1)
  expect_gt(stats::quantile(s, 0.9) / stats::quantile(s, 0.5), stats::quantile(s, 0.5) / stats::quantile(s, 0.1) - 0.05)
  # draws leave the session RNG alone
  withr::local_seed(9); r1 <- stats::runif(1)
  withr::local_seed(9); invisible(draw_residual_field(n[1:5], f, "C")); expect_equal(stats::runif(1), r1)
})

test_that("the smear kernel spreads a codon symmetrically and keeps the interior total", {
  x <- c(rep(1, 10), 11, rep(1, 10))
  y <- apply_codon_smear(x, c(1, 0.5, 0.25))
  k <- c(0.25, 0.5, 1, 0.5, 0.25) / 2.5
  expect_equal(y[11], 1 + 10 * k[3]); expect_equal(y[c(10, 12)], 1 + 10 * k[c(2, 4)])
  expect_equal(y[c(9, 13)], 1 + 10 * k[c(1, 5)]); expect_equal(y[c(8, 14)], c(1, 1))
  expect_equal(sum(y), sum(x))
  expect_equal(apply_codon_smear(rep(2, 6), c(1, 0.5, 0.25)), rep(2, 6))    # edges renormalised: flat stays flat
  expect_identical(apply_codon_smear(x, 1), x); expect_identical(apply_codon_smear(x, NULL), x)
})

test_that("contrast mode raises each gene's codon weights to its own exponent N(mean, sd)", {
  withr::local_options(coverageSim.codon_bias = list(sharpen = "contrast", contrast_mean = 1.03, contrast_sd = 0.15))
  e <- draw_contrast_exponents(5000, codon_bias_options())
  expect_equal(mean(e), 1.03, tolerance = 0.01); expect_equal(stats::sd(e), 0.15, tolerance = 0.05)
  expect_equal(draw_contrast_exponents(5000, codon_bias_options()), e)       # per genome
  # no smoothing and no RNase: the log weights are an exact linear function of the log table, slope = exponent
  g <- s14_genes(20, 300)
  withr::local_options(coverageSim.codon_bias = list(sharpen = "contrast", smooth = FALSE, per_genome_draws = TRUE))
  out <- sim_sequence_bias(NULL, 900, g$alpha, shapes(9), rnase_acf = NULL)
  ex <- draw_contrast_exponents(20, codon_bias_options())
  slope <- vapply(seq_along(out), function(i) {
    w <- out[[i]][seq(1, length(out[[i]]), 3)]
    stats::coef(stats::lm(log(w) ~ log(g$alpha[[i]])))[[2]] }, numeric(1))
  expect_equal(slope, ex, tolerance = 1e-6)
})

test_that("footprint GC: window GC per P-site codon and the inverted-U factor of M22", {
  cod <- c("GCG", "GCG", "AAA", "AAA", "GCA", "TTT")
  # window -1..+1 around each P-site codon, cut at the ends; weights indexed by P site (a_shift = 1)
  expect_equal(footprint_gc(cod, 1L, c(-1L, 1L)), c(1, 2 / 3, 1 / 3, 2 / 9, 2 / 9, 1 / 3))
  # A-site indexed weights (a_shift = 0): the codon whose P site is p sits at p + 1; position 1 gets the gene mean
  expect_equal(footprint_gc(cod, 0L, c(-1L, 1L)), c(mean(c(1, 1, 0, 0, 2 / 3, 0)), 1, 2 / 3, 1 / 3, 2 / 9, 2 / 9))
  f <- residual_field_options(gc_bias = TRUE)
  h <- gc_log_factor(c(0.2, 0.37, 0.41, 0.5, 0.6, 0.65, 0.7, 0.9), f)
  expect_equal(h, c(-0.13, -0.13, -0.035, 0.06, 0.06, -0.035, -0.13, -0.13))
  expect_equal(gc_log_factor(0.5, residual_field_options(gc_strength = 2)), 0.12)
  # in sim_sequence_bias (no other layer): log weight change = h(footprint GC)
  g <- s14_genes(30, 300)
  base <- list(sharpen = "fixed", sharpen_power = 1, smooth = FALSE, per_genome_draws = TRUE)
  withr::local_options(coverageSim.codon_bias = base)
  off <- sim_sequence_bias(NULL, 900, g$alpha, shapes(9), rnase_acf = NULL, codons = g$genes, a_shift = 0L)
  withr::local_options(coverageSim.codon_bias = c(base, list(field = list(gc_bias = TRUE))))
  on <- sim_sequence_bias(NULL, 900, g$alpha, shapes(9), rnase_acf = NULL, codons = g$genes, a_shift = 0L)
  d <- unlist(lapply(seq_along(on), function(i) { r <- log(on[[i]] / off[[i]])[seq(1, 900, 3)]; r - mean(r) }))
  h <- unlist(lapply(g$genes, function(cc) { x <- gc_log_factor(footprint_gc(cc, 0L), f); x - mean(x) }))
  expect_equal(d, h, tolerance = 1e-8)
  # low/high GC footprints lose about 0.19 log relative to 45-60% GC
  gc <- unlist(lapply(g$genes, footprint_gc, a_shift = 0L))
  expect_lt(abs(mean(d[gc >= 0.45 & gc <= 0.6]) - mean(d[gc <= 0.37 | gc >= 0.7]) - 0.19), 0.03)
})

test_that("simNGScoverage hands the library's study to the field", {
  ct <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1:4, 2)),
                                                   colData = S4Vectors::DataFrame(study = c("PRJ1", "PRJ2")))
  expect_equal(field_study(ct[, 2]), "PRJ2")
  ct2 <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(1:4, 2)))
  expect_equal(field_study(ct2[, 1]), "study1")
  colnames(ct) <- c("RFP_WT_1", "RFP_WT_2")
  expect_equal(field_study(ct[, 2], unit = "library"), "RFP_WT_2")
  expect_error(residual_field_options(study_unit = "lab"))
  g <- s14_genes(10, 300)
  withr::local_options(coverageSim.codon_bias = list(sharpen = "fixed", sharpen_power = 1, per_genome_draws = TRUE,
                                                     field = list(field_var = c(0.1, 0.1, 0.1), shared = 0.5)))
  a <- sim_sequence_bias(NULL, 900, g$alpha, shapes(9), study = "PRJ1")
  expect_equal(sim_sequence_bias(NULL, 900, g$alpha, shapes(9), study = "PRJ1"), a)
  expect_false(isTRUE(all.equal(sim_sequence_bias(NULL, 900, g$alpha, shapes(9), study = "PRJ2"), a)))
})

# Two merged half-pools ("A", "B": the two study-disjoint halves of the merged real track, study_unit = "library") of
# 500 random-codon genes (300-1,500 codons, lognormal depth around 40 reads per codon), each the sum of 4
# Dirichlet-multinomial libraries (dmn scale 6). Residuals vs the A-site codon table (the codon-table-only explainer).
s14_pool_sim <- function(cb, n = 500, seed = 21) {
  tab <- load_seq_bias(type = "codon", shift = "a-site")
  sense <- setdiff(tab$seqs, c("###", "%%%", "&&&", "***", "TAA", "TAG", "TGA"))
  withr::local_seed(seed)
  len <- round(exp(stats::runif(n, log(300), log(1500))))
  genes <- lapply(len, function(m) sample(sense, m, TRUE)); alpha <- lapply(genes, function(g) tab$alpha[match(g, tab$seqs)])
  dep <- exp(stats::rnorm(n, log(40), 1.4))
  withr::local_options(coverageSim.codon_bias = cb)
  pool <- lapply(c("A", "B"), function(st) {
    res <- scale_dmn_alpha(sim_sequence_bias(NULL, lengths(alpha) * 3, alpha, shapes(9), codons = genes, a_shift = 0L,
                                             study = st), 6)
    lapply(seq_along(res), function(i) {
      nc <- length(genes[[i]]); m <- (length(res[[i]]) - 3 * nc) %/% 2; s <- 0
      for (k in 1:4) s <- s + stats::rmultinom(1, round(dep[i] * nc / 8), draw_site_probabilities(res[[i]], TRUE))[, 1]
      colSums(matrix(s[m + seq_len(3 * nc)], 3))
    })
  })
  data.table::rbindlist(lapply(seq_len(n), function(i) data.table::data.table(t = i, k = seq_along(alpha[[i]]),
    OA = pool[[1]][[i]], OB = pool[[2]][[i]], e = alpha[[i]])))
}

test_that("human_v2 preset: block variance, window ACF, peak flanks, gene dispersion, split halves, table-E flags", {
  C <- s14_pool_sim(codon_bias_preset("human_v2"))
  C[, O := OA + OB]; C[, E := sum(O) * e / sum(e), by = t]; C[, `:=`(n = .N, g = mean(O)), by = t]
  # 1. within-gene variance of log O/E by block (genes >= 200 reads/codon; real 0.34 / 0.26 / 0.17 / 0.10 / 0.055)
  D <- C[g >= 200 & n >= 300]
  bv <- vapply(c(1, 3, 10, 30, 100), function(bs) {
    X <- D[, .(O = sum(O), E = sum(E), m = .N), by = .(t, b = (k - 1) %/% bs)][m == bs]
    X[, r := log((O + 1) / (E + 1))]; X[, rc := r - mean(r), by = t]; X[, stats::var(rc)] }, numeric(1))
  expect_true(all(abs(log(bv / c(0.34, 0.26, 0.17, 0.10, 0.055))) < log(1.5)), info = paste(round(bv, 3), collapse = " "))
  expect_gt(bv[4] / bv[1], 3 * 0.34 / 30 / 0.34)                             # far redder than white noise (1/30)
  # 2. gene-centred window ACF (30 codons, E > 1000, genes >= 20 windows; real 0.51 / 0.41 / 0.35 / 0.29)
  W <- C[, .(O = sum(O), E = sum(E), OA = sum(OA), OB = sum(OB), m = .N), by = .(t, w = (k - 1) %/% 30)][m >= 15]
  W[, lr := log(pmax(O, 0.5) / pmax(E, 0.5))]
  V <- W[E > 1000]; V[, lc := lr - mean(lr), by = t]; V[, nw := .N, by = t]; V <- V[nw >= 20]
  acf <- vapply(1:4, function(L) { V[, sh := data.table::shift(lc, L), by = t]; V[, stats::cor(lc, sh, use = "c")] }, numeric(1))
  expect_true(all(abs(acf - c(0.51, 0.41, 0.35, 0.29)) < 0.15), info = paste(round(acf, 3), collapse = " "))
  # 3. mean gene-centred log residual around strong peaks at 0 / +-1 / +-2 / +-5 (real 2.15 / 1.36 / 0.96 / 0.70)
  P <- C[g >= 100]; P[, r := log((O + 1) / (E + 1))]; P[, r := r - mean(r), by = t]; P[, i := .I]
  pk <- P[r > log(4) & O > 5 * g & k > 60, i]
  expect_gt(length(pk), 100)
  fl <- vapply(c(0, 1, 2, 5), function(d) { j <- c(pk - d, pk + d); s <- rep(pk, 2)
    ok <- j >= 1 & j <= nrow(P); ok[ok] <- P$t[j[ok]] == P$t[s[ok]]; mean(P$r[j[ok]]) }, numeric(1))
  expect_true(all(abs(fl - c(2.15, 1.36, 0.96, 0.70)) < 0.35), info = paste(round(fl, 2), collapse = " "))
  expect_true(all(diff(fl) < 0))                                              # symmetric decay, no queue
  # 4. split halves: window residual r between the two pools (real: 20-study pools 0.6 at E > 50, NC-1 0.78 deep)
  W[, `:=`(EA = E * sum(OA) / sum(E), EB = E * sum(OB) / sum(E)), by = t]
  X <- W[EA > 50 & EB > 50]
  r <- X[, stats::cor(log(pmax(OA, 0.5) / EA), log(pmax(OB, 0.5) / EB))]
  expect_gt(r, 0.6); expect_lt(r, 0.9)
  # 5. gene residual sd (windows E > 300, genes >= 5 windows; real 10-90% 0.16-0.47, heavy right tail)
  gs <- W[E > 300, .(nw = .N, s = stats::sd(lr)), by = t][nw >= 5]$s
  q <- stats::quantile(gs, c(0.1, 0.5, 0.9))
  expect_true(abs(q[[1]] - 0.16) < 0.05 && abs(q[[3]] - 0.47) < 0.1, info = paste(round(q, 3), collapse = " "))
  expect_gt(q[[3]] - q[[2]], q[[2]] - q[[1]])
  # 6. table-E density flags with the real explainer band (NP + frozen tau2; real density-fail 8.7%): 7-9% target
  vC <- function(E) ifelse(E < 682, 2 * pmax(E, 0.5) + 0.0162 * E^2, 0.13^2 * E^2)
  f <- W[OA >= 10 & OB >= 10, mean((O - E)^2 / (vC(E) + 0.0386 * E^2) > 3.84 & abs(O / E - 1) >= 0.25)]
  expect_gt(f, 0.04); expect_lt(f, 0.13)
})
