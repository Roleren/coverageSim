test_that("realistic preset keeps t57 and only changes the pause settings", {
  t57 <- codon_bias_preset("t57"); real <- codon_bias_preset("realistic")
  expect_equal(t57$pause_sd, 1.2); expect_equal(t57$pause_prob, 0.45)
  expect_equal(real$pause_sd, 1.0); expect_equal(real$pause_prob, 0.7)
  expect_equal(real[setdiff(names(real), c("pause_sd", "pause_prob"))],
               t57[setdiff(names(t57), c("pause_sd", "pause_prob"))])
  withr::local_options(coverageSim.codon_bias = real)
  expect_equal(codon_bias_options()$pause_sd, 1.0)
})

test_that("realistic_rfp_settings returns the validated t69 arguments", {
  s <- realistic_rfp_settings()
  expect_equal(names(s), c("dmn_alpha_scale", "fragment_geometry"))
  expect_equal(s$dmn_alpha_scale, 6)
  g <- s$fragment_geometry
  expect_equal(g$source, "learned"); expect_equal(g$site_reference, "a_site")
  expect_equal(g$frame_bias$strength, 0.25)
  expect_equal(sort(unique(g$distribution$fragment_length)), 25:36)
  expect_equal(nrow(g$distribution), 58L)
  expect_equal(sum(g$distribution$probability), 1, tolerance = 0.01)
  expect_equal(nrow(g$frame_bias$table), 36L)
  expect_equal(nrow(coverageSim:::normalize_fragment_geometry(g)$distribution), 58L)
})

test_that("a tiny simulation runs with the realistic settings", {
  set.seed(9)
  genome <- suppressWarnings(suppressMessages(simGenome(
    n = 4, out_dir = tempfile("genome-"), cds_length = rep(600, 4), max_uorfs = 0, debug_on = FALSE)))
  cds <- ORFik::loadRegion(genome["txdb"], "cds")
  gene <- matrix(2000L, length(cds), 1L, dimnames = list(names(cds), "RFP_x"))
  counts <- SummarizedExperiment::SummarizedExperiment(
    assays = list(gene = gene, leader = gene * 0L, cds = gene), rowRanges = cds,
    colData = S4Vectors::DataFrame(libtype = factor("RFP"), condition = factor("x"),
                                   replicate = "1", row.names = colnames(gene)))
  withr::local_options(coverageSim.codon_bias = codon_bias_preset("realistic"))
  out_dir <- tempfile("reads-"); exp_dir <- tempfile("exp-"); dir.create(out_dir); dir.create(exp_dir)
  on.exit(unlink(c(out_dir, exp_dir), recursive = TRUE))
  experiment <- suppressMessages(do.call(simNGScoverage, c(list(
    simGenome = genome, count_table = counts, out_dir = out_dir, exp_name = "realistic",
    exp_save_dir = exp_dir, libFormats = list(RFP = "ofst"), validate = FALSE),
    realistic_rfp_settings())))
  total <- sum(S4Vectors::mcols(ORFik::fimport(ORFik::filepath(experiment, "default")[1]))$score)
  expect_gt(total, 0)
})
