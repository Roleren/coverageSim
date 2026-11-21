test_that("utr3 defaults switch everything off and options validate", {
  withr::local_options(coverageSim.utr3 = NULL)
  o <- utr3_options()
  expect_false(utr3_active(o))
  expect_null(utr3_run_features(c(a = 900L), "a"))
  expect_equal(nrow(draw_utr3_features(c(a = 900L), opts = o)), 0L)
  withr::local_options(coverageSim.utr3 = list(oops = 1))
  expect_error(utr3_options(), "unknown")
  withr::local_options(coverageSim.utr3 = list(plateau_fold = 1))
  expect_error(utr3_options())
  withr::local_options(coverageSim.utr3 = utr3_preset("human"))
  expect_true(utr3_active())
  expect_equal(utr3_options()$apa_ratio, 0.51)
})

test_that("3'UTR features: fixed per genome, measured rates and positions, session RNG untouched", {
  withr::local_options(coverageSim.utr3 = utr3_preset("human"))
  lengths <- stats::setNames(rep(c(450L, 1800L, 3600L), 400), paste0("tx", 1:1200))
  set.seed(5); a <- draw_utr3_features(lengths); r1 <- stats::runif(1)
  set.seed(6); b <- draw_utr3_features(lengths)
  expect_equal(a, b)
  set.seed(5); draw_utr3_features(lengths); expect_equal(stats::runif(1), r1)
  pl <- a[type == "plateau"]
  # never in the first or the last window, always one 90-nt window
  nw <- lengths[pl$gene] %/% 90L
  expect_true(all(pl$start > 90L & pl$end <= (nw - 1L) * 90L & pl$end - pl$start + 1L == 90L))
  eligible <- sum(pmax(lengths %/% 90L - 2L, 0L))
  expect_equal(nrow(pl) / eligible, 0.045, tolerance = 0.15)
  expect_equal(mean((pl$start + pl$end) / 2 <= lengths[pl$gene] / 2), 0.85, tolerance = 0.06)
  expect_equal(stats::median(pl$factor), 3, tolerance = 0.1)
  apa <- a[type == "apa"]
  expect_equal(nrow(apa) / length(lengths), 0.54, tolerance = 0.08)
  expect_true(all(apa$end == lengths[apa$gene] & apa$start > 0.1 * lengths[apa$gene]))
  expect_equal(nrow(a[type == "end"]), length(lengths))
})

test_that("3'UTR multipliers: plateaus and the depleted end only in ribo libraries, APA in all", {
  withr::local_options(coverageSim.utr3 = list(end_depletion = 0.1, end_width = 10L))
  o <- utr3_options()
  f <- data.table::data.table(transcript_id = "t", gene = 1L, type = c("plateau", "apa", "end"),
                              start = c(91L, 271L, 351L), end = c(180L, 360L, 360L), factor = c(3, 0.5, 0.1))
  rfp <- utr3_factor_vector(360L, f, "RFP", o)
  rna <- utr3_factor_vector(360L, f, "RNA", o)
  expect_equal(rfp[c(1, 91, 180, 181, 271, 351)], c(1, 3, 3, 1, 0.5, 0.05))
  expect_equal(rna[c(1, 91, 271, 351)], c(1, 1, 0.5, 0.5))
  # a margin around the 3'UTR (RNase padding of DMN sampling) is left alone
  w <- apply_utr3_factors(rep(1, 366), 360L, 1L, "RFP", list(opts = o, by_gene = list(f)))
  expect_equal(w[c(1:3, 364:366)], rep(1, 6)); expect_equal(w[3 + 91], 3)
})

test_that("simNGScoverage plants 3'UTR plateaus and APA steps in the reads and writes the truth table", {
  set.seed(41)
  genome_dir <- tempfile("utr3-genome-"); exp_dir <- tempfile("utr3-exp-"); out_dir <- tempfile("utr3-out-")
  dir.create(genome_dir); dir.create(exp_dir); dir.create(out_dir)
  sg <- suppressWarnings(simGenome(n = 4, out_dir = genome_dir, max_uorfs = 0, cds_length = rep(300, 4),
                                   leader_length = rep(30, 4), trailer_length = rep(1800, 4), strand = rep("+", 4),
                                   export_txdb = TRUE))
  cds <- ORFik::loadRegion(sg["txdb"], "cds")
  ct <- simCountTables(cds, libtypes = "RFP", print_statistics = FALSE, plot_PCA = FALSE, interceptMean = 11,
                       replicates = 1)
  rt <- simCountTablesRegions(ct, regionsToSample = c("cds", "trailer"),
                              region_proportion = list(cds = list(RFP = 0.2), trailer = list(RFP = 0.8)))
  withr::local_options(coverageSim.utr3 = list(plateau_prob = 0.25, plateau_fold = 4, plateau_fold_sdlog = 0,
                                               apa_prob = 1, apa_ratio = 0.5, seed = 3L))
  experiment <- simNGScoverage(sg, rt[, 1], out_dir = out_dir, exp_name = basename(out_dir), exp_save_dir = exp_dir,
                               validate = FALSE, ground_truth = TRUE, seq_bias = NULL,
                               rnase_bias = list(RFP = NULL), auto_correlation = NULL)
  feats <- data.table::fread(file.path(out_dir, "utr3_features.tsv"))
  expect_equal(feats, draw_utr3_features(stats::setNames(rep(1800L, 4), rownames(rt)), rownames(rt)),
               ignore_attr = TRUE)
  truth <- data.table::fread(list.files(out_dir, "_ground_truth.tsv$", full.names = TRUE)[1])
  utr <- truth[startsWith(sampling_group, "trailer")]
  # ribosome_site is genomic; all genes are on the plus strand, so 3'UTR position 1 follows the last CDS codon
  first_site <- truth[startsWith(sampling_group, "cds"), max(ribosome_site), by = transcript_id]
  utr[first_site, pos := ribosome_site - (i.V1 + 2L), on = "transcript_id"]
  utr[, gene := match(transcript_id, rownames(rt))]
  sites <- data.table::rbindlist(lapply(seq_len(nrow(rt)), function(g) {
    data.table::data.table(gene = g, pos = 121:1679,
                           f = utr3_factor_vector(1800L, feats[gene == g], "RFP", utr3_options())[121:1679])
  }))
  sites[utr[, .(reads = sum(score)), by = .(gene, pos)], reads := i.reads, on = c("gene", "pos")]
  sites[is.na(reads), reads := 0]
  # density relative to the gene's own ordinary positions (genes differ in depth and in how much is planted)
  sites[, rel := reads / mean(reads[f == 1]), by = gene]
  dens <- sites[, .(d = mean(rel)), by = f]
  # plateau windows hold ~4x, APA-downstream ~0.5x the reads of an ordinary window
  expect_equal(dens[f == 4, d], 4, tolerance = 0.2)
  expect_equal(dens[f == 0.5, d], 0.5, tolerance = 0.2)
})
