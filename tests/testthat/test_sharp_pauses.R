# Port S13: sharp pauses, Pro-Pro stalls and Pro-run depletion (claude-megafst M16-acf, 2026-11-21).

test_that("S13 options are off by default, validate and leave the simulation unchanged", {
  withr::local_options(coverageSim.codon_bias = NULL)
  o <- codon_bias_options()
  expect_false(o$pause_sharp); expect_false(o$pp_stall); expect_false(o$pro_run)
  expect_false(s13_needs_codons(o))
  withr::local_options(coverageSim.codon_bias = list(pp_stall_take = 0.2)); expect_error(codon_bias_options())
  withr::local_options(coverageSim.codon_bias = list(pause_flank_take = 0)); expect_error(codon_bias_options())
  withr::local_options(coverageSim.codon_bias = list(pro_run_factor = 0.5)); expect_error(codon_bias_options())
  p <- codon_bias_preset("human")
  expect_true(p$pp_stall && p$pro_run)
  withr::local_options(coverageSim.codon_bias = p); expect_true(s13_needs_codons())
  # codons passed without an S13 layer change nothing
  withr::local_options(coverageSim.codon_bias = codon_bias_preset("realistic"))
  tab <- load_seq_bias(type = "codon", shift = "a-site")
  set.seed(1); g <- list(sample(c("CCA", "CCT", "GAA", "AAA"), 200, TRUE))
  a <- list(tab$alpha[match(g[[1]], tab$seqs)])
  expect_equal(sim_sequence_bias(NULL, 600, a, shapes(9), codons = g, a_shift = 0L),
               sim_sequence_bias(NULL, 600, a, shapes(9)))
})

test_that("sharp pauses take their gain from the flanks; dips are only multiplied", {
  x <- c(1, 2, 3, 2, 1, 4, 1)
  m <- c(1, 1, 2, 1, 1, 1, 1)
  y <- sharp_pause_redistribute(x, m, width = 1L, take = 0.9)
  expect_equal(sum(y), sum(x))
  expect_equal(y[3], 6)                       # 3 + (2 - 1) * 3
  expect_equal(y[c(2, 4)], c(2, 2) * (1 - 3 / 4))
  expect_equal(y[-(2:4)], x[-(2:4)])
  # the gain is capped at `take` of the flank total
  expect_equal(sharp_pause_redistribute(x, m, 1L, take = 0.25)[3], 3 + 0.25 * 4)
  # width 2 spreads the loss over four codons
  y2 <- sharp_pause_redistribute(x, m, width = 2L, take = 0.9)
  expect_equal(sum(y2), sum(x)); expect_equal(y2[3], 6)
  # a dip is multiplicative and leaves its neighbours alone
  expect_equal(sharp_pause_redistribute(x, c(1, 1, 0.5, 1, 1, 1, 1), 1L, 0.5), x * c(1, 1, 0.5, 1, 1, 1, 1))
  expect_true(all(sharp_pause_redistribute(x, c(9, 9, 9, 0.1, 9, 9, 9), 2L, 1) > 0))
})

test_that("the Pro-Pro stall kernel and the Pro-run depletion use P-site coordinates", {
  # P sites: codon 3 has Pro at E (2) and P (3) and a non-Pro A site (4): a stall. In the run PPP (6-8) only codon 8
  # (Pro at E and P, Lys at A) is one.
  cod <- c("ATG", "CCA", "CCT", "GAA", "AAA", "CCC", "CCG", "CCA", "AAA", "AAA")
  expect_equal(pp_stall_sites(cod), c(3L, 8L))
  x <- rep(1, 10)
  # P-site indexed weights (a_shift = 1): the stall codon is position 3
  y <- apply_pp_stall(x, cod, a_shift = 1L, take = c(0.2, 0.1), peak = 1)
  expect_equal(y[2:4], c(0.8, 1.3, 0.9)); expect_equal(sum(y), 10)
  expect_equal(apply_pp_stall(x, cod, 1L, c(0.2, 0.1), peak = 0.9)[3], 1.3 * 0.9)
  # A-site indexed weights (a_shift = 0): the codon whose P site is 3 sits at position 4
  y0 <- apply_pp_stall(x, cod, a_shift = 0L, take = c(0.2, 0.1), peak = 1)
  expect_equal(y0[3:5], c(0.8, 1.3, 0.9))
  # width 2
  y2 <- apply_pp_stall(x, cod, 1L, take = c(0.1, 0.2, 0.1, 0.05), peak = 1)
  expect_equal(y2[1:5], c(0.9, 0.8, 1.45, 0.9, 0.95))
  # Pro run: Pro count in P-site +-2 per codon
  pro <- as.numeric(codon_amino_acid(cod) == "P")
  nP <- pro + flank_sum(pro, 2L)
  expect_equal(nP, c(2, 2, 2, 3, 3, 3, 3, 3, 2, 1))
  f <- c(`3` = 0.5)
  expect_equal(apply_pro_run(x, cod, 1L, f), ifelse(nP == 3, 0.5, 1))
  expect_equal(apply_pro_run(x, cod, 0L, f), c(1, ifelse(nP == 3, 0.5, 1)[-10]))
  expect_error(apply_s13_codon_layers(x, cod[-1], 1L, codon_bias_options()), "do not match")
})

# Deep simulated libraries under codon_bias_preset("human"): 1,000 genes x 400 random sense codons, 20 libraries of
# 1,500 reads per gene merged (30,000 reads per gene, ~75 per codon), like the merged real track.
s13_deep_sim <- function(cb, genes, alpha, depth = 30000, libs = 20) {
  withr::local_options(coverageSim.codon_bias = cb)
  res <- scale_dmn_alpha(sim_sequence_bias(NULL, lengths(alpha) * 3, alpha, shapes(9), codons = genes, a_shift = 0L), 6)
  withr::local_seed(11)
  lapply(seq_along(res), function(i) {
    s <- 0
    for (k in seq_len(libs)) s <- s + stats::rmultinom(1, depth / libs, draw_site_probabilities(res[[i]], TRUE))[, 1]
    m <- (length(s) - 3 * length(genes[[i]])) %/% 2
    colSums(matrix(s[m + seq_len(3 * length(genes[[i]]))], 3))
  })
}

test_that("human preset: background ACF, Pro-Pro stall profile and Pro-run depletion match real data", {
  tab <- load_seq_bias(type = "codon", shift = "a-site")
  sense <- setdiff(tab$seqs, c("###", "%%%", "&&&", "***", "TAA", "TAG", "TGA"))
  withr::local_seed(3)
  genes <- lapply(1:1000, function(i) sample(sense, 400, TRUE))
  alpha <- lapply(genes, function(g) tab$alpha[match(g, tab$seqs)])
  cb <- codon_bias_preset("human")
  on <- s13_deep_sim(cb, genes, alpha)
  off <- s13_deep_sim(modifyList(cb, list(pp_stall = FALSE, pro_run = FALSE)), genes, alpha)
  cen <- lapply(on, function(o) { x <- log(o + 1); x - mean(x) })
  cen_off <- lapply(off, function(o) { x <- log(o + 1); x - mean(x) })
  # 1. per-codon ACF of gene-centred log(O + 1), lags 1-5 (real, M16-acf: 0.58 / 0.42 / 0.37 / 0.33 / 0.28)
  acf <- vapply(1:5, function(l) stats::cor(unlist(lapply(cen, function(x) utils::head(x, -l))),
                                            unlist(lapply(cen, function(x) utils::tail(x, -l)))), numeric(1))
  expect_true(all(abs(acf - c(0.58, 0.42, 0.37, 0.33, 0.28)) <= 0.08), info = paste(round(acf, 3), collapse = " "))
  # 2. residual (S13 on vs off) around Pro-Pro stalls at P-site offsets -1 / 0 / +1 (real: -0.14 / +0.10 / -0.09).
  # Weights are A-site indexed here (a_shift = 0), so the codon of P-site p is at p + 1.
  st <- lapply(genes, function(g) pp_stall_sites(g) + 1L)
  expect_gt(sum(lengths(st)), 1000)
  resid <- vapply(-1:1, function(d) mean(unlist(mapply(function(x, y, s) {
    j <- s + d; j <- j[j >= 1 & j <= length(x)]; x[j] - y[j] }, cen, cen_off, st))), numeric(1))
  real <- c(-0.14, 0.10, -0.09)
  expect_equal(sign(resid), sign(real))
  expect_true(all(abs(resid - real) <= 0.5 * abs(real)), info = paste(round(resid, 3), collapse = " "))
  # 3. Pro-run depletion: on minus off by Pro count in P-site +-2 (target exp factors: 3 -0.27, 4-5 -0.80)
  nP <- lapply(genes, function(g) { p <- as.numeric(codon_amino_acid(g) == "P"); n <- p + flank_sum(p, 2L)
    c(NA, utils::head(pmin(n, 4), -1)) })
  d <- tapply(unlist(cen) - unlist(cen_off), unlist(nP), mean)
  expect_lt(d[["3"]] - d[["0"]], -0.27 * 0.5); expect_gt(d[["3"]] - d[["0"]], -0.27 * 1.5)
  expect_lt(d[["4"]] - d[["0"]], -0.80 * 0.5); expect_gt(d[["4"]] - d[["0"]], -0.80 * 1.5)
})

test_that("simNGScoverage passes codons to the S13 layers", {
  set.seed(77)
  genome_dir <- tempfile("s13-genome-"); exp_dir <- tempfile("s13-exp-"); out_dir <- tempfile("s13-out-")
  dir.create(genome_dir); dir.create(exp_dir); dir.create(out_dir)
  sg <- suppressWarnings(simGenome(n = 4, out_dir = genome_dir, max_uorfs = 0, cds_length = rep(1500, 4),
                                   export_txdb = TRUE))
  cds <- ORFik::loadRegion(sg["txdb"], "cds")
  ct <- simCountTables(cds, libtypes = "RFP", print_statistics = FALSE, plot_PCA = FALSE, interceptMean = 11,
                       replicates = 1)
  rt <- simCountTablesRegions(ct, regionsToSample = "cds")
  pro_share <- function(cb, tag) {
    withr::local_options(coverageSim.codon_bias = cb)
    od <- file.path(out_dir, tag); dir.create(od)
    suppressMessages(simNGScoverage(sg, rt, out_dir = od, exp_name = tag, exp_save_dir = exp_dir,
                                    validate = FALSE, ground_truth = TRUE, rnase_bias = list(RFP = NULL)))
    t <- data.table::fread(list.files(od, "_ground_truth.tsv$", full.names = TRUE)[1])
    t[codon %in% c("CCA", "CCC", "CCG", "CCT"), .N] / nrow(t)
  }
  base <- codon_bias_preset("realistic")
  # every codon with a Pro in its P-site +-2 window (which includes the site codon itself) nearly vanishes
  depleted <- modifyList(base, list(pro_run = TRUE, pro_run_factor = c(`1` = 1e-3, `2` = 1e-3, `3` = 1e-3,
                                                                       `4` = 1e-3, `5` = 1e-3)))
  before <- pro_share(base, "off"); after <- pro_share(depleted, "on")
  expect_gt(before, 0.02)
  expect_lt(after, 0.2 * before)
})
