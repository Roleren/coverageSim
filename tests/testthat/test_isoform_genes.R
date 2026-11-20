test_that("isoform genes follow the real distribution of n_iso and reference weight", {
  x <- simIsoformGenes(4000, seed = 1)
  g <- x$genes
  expect_equal(mean(g$w_ref < 0.5), 3708 / 17935, tolerance = 0.12)   # relative: 20.7 % of genes MANE-minor
  expect_equal(mean(g$w_ref < 0.1), 1065 / 17935, tolerance = 0.2)
  expect_equal(mean(g$n_iso == 1), 4918 / 17935, tolerance = 0.1)
  expect_true(all(g$w_ref[g$n_iso == 1] == 1))
  # weights of every gene sum to one and the reference carries w_ref
  s <- x$iso[, .(s = sum(w_rna), r = w_rna[is_ref]), by = gene]
  expect_equal(s$s, rep(1, nrow(s)), tolerance = 1e-8)
  expect_equal(s$r, g$w_ref, tolerance = 1e-8)
  # exactly one reference per gene; segments tile the reference CDS
  expect_equal(x$iso[, sum(is_ref), by = gene]$V1, rep(1L, nrow(g)))
  tile <- x$seg[, .(n = sum(end - start + 1L)), by = gene]
  expect_equal(tile$n, g$cds_codons)
})

test_that("ref_minor_frac sets the fraction of MANE-minor multi-isoform genes; default NULL keeps the real mix", {
  x <- simIsoformGenes(3000, ref_minor_frac = 0.6, seed = 2)$genes
  expect_equal(mean(x$w_ref[x$n_iso > 1] < 0.5), 0.6, tolerance = 0.08)
  x0 <- simIsoformGenes(3000, ref_minor_frac = 0, seed = 2)$genes
  expect_equal(sum(x0$w_ref[x0$n_iso > 1] < 0.5), 0)
  y <- simIsoformGenes(300, n_iso = 4, seed = 3)$genes
  expect_true(all(y$n_iso == 4))
})

test_that("expected profiles: truth mixture differs from the flat reference only where isoforms skip segments", {
  x <- simIsoformGenes(300, ref_minor_frac = 1, seed = 4)
  ex <- isoform_expected(x)
  expect_equal(unname(vapply(ex$truth_ribo, sum, 0)), rep(1, 300), tolerance = 1e-8)
  expect_equal(unname(vapply(ex$flat, sum, 0)), rep(1, 300), tolerance = 1e-8)
  one <- which(x$genes$n_iso == 1)
  if (length(one)) expect_equal(ex$truth_ribo[one], ex$flat[one])
  multi <- which(x$genes$n_iso > 1 & x$genes$frac_specific > 0.1)
  expect_gt(length(multi), 20)
  dev <- vapply(multi, function(i) max(abs(ex$truth_ribo[[i]] / ex$flat[[i]] - 1)), 0)
  expect_gt(median(dev), 0.1)
  # without TE and without RNA mismatch, ribosomes and RNA see the same isoform usage
  u <- vapply(multi, function(i) cor(ex$truth_ribo[[i]] / ex$flat[[i]], ex$truth_rna[[i]]), 0)
  expect_gt(min(u, na.rm = TRUE), 0.999)
  # with TE the ribosome mixture departs from the RNA mixture
  xt <- simIsoformGenes(300, ref_minor_frac = 1, te_sd = 1, seed = 4)
  et <- isoform_expected(xt); m <- which(xt$genes$n_iso > 2 & xt$genes$frac_specific > 0.2)
  ut <- vapply(m, function(i) cor(et$truth_ribo[[i]] / et$flat[[i]], et$truth_rna[[i]]), 0)
  expect_lt(median(ut, na.rm = TRUE), 0.99)
})

test_that("simIsoformReads: ribo and RNA tracks come from the same weights, counts match the expectation, RNG restored", {
  x <- simIsoformGenes(200, ref_minor_frac = 1, seed = 5)
  set.seed(9); r1 <- stats::runif(1)
  set.seed(9); r <- simIsoformReads(x, ribo_depth = 20000, rna_depth = 20000, replicates = 2, seed = 6); expect_equal(stats::runif(1), r1)
  expect_equal(sum(r$O1) / sum(r$E_truth), 1, tolerance = 0.02)
  expect_equal(sum(r$R) / sum(r$E_rna), 1, tolerance = 0.02)
  expect_equal(r[, sum(ref_transcript), by = gene]$V1, r[, sum(O1 + O2) / 2, by = gene]$V1, tolerance = 1e-8)
  # the truth E explains the reads (mean Poisson deviance per codon ~ 1) far better than the flat reference in minor genes
  dv <- function(O, E) mean(2 * (ifelse(O > 0, O * log(O / E), 0) - (O - E)))
  sel <- r$w_ref < 0.3
  expect_gt(sum(sel), 5000)
  expect_lt(dv(r$O1[sel], r$E_truth[sel]), 1.1)
  expect_gt(dv(r$O1[sel], r$E_flat[sel]), 1.5)
})

test_that("isoform_mixture_alpha feeds sim_sequence_bias()", {
  x <- simIsoformGenes(5, cds_codons = 120, seed = 7)
  a <- isoform_mixture_alpha(x)
  withr::local_options(coverageSim.codon_bias = NULL)
  out <- sim_sequence_bias(NULL, lengths(a) * 3, a, shapes(9))
  expect_length(out, 5)
  expect_true(all(is.finite(unlist(out))))
})
