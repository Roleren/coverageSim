# Default RFP region shares and read lengths are the measured human medians
# (757 libraries / 74 studies, 702 runs for lengths; carina-claude C2, 2026-11-18).

test_that("default RFP region shares fall within the measured human p25-p75", {
  set.seed(1)
  dir <- tempfile("defaults_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = 40, out_dir = dir, max_uorfs = 1, debug_on = FALSE))
  counts <- suppressMessages(simCountTables(ORFik::loadRegion(genome["txdb"], "cds"),
    libtypes = "RFP", conditions = "WT", replicates = 1, plot_PCA = FALSE,
    print_statistics = FALSE))
  SummarizedExperiment::assay(counts)[] <- 2000
  regions <- suppressMessages(simCountTablesRegions(counts))
  total <- sum(SummarizedExperiment::assay(regions, "gene"))
  share <- function(a) sum(SummarizedExperiment::assay(regions, a)) / total
  # In real data uORF reads are part of the leader.
  expect_gt(share("leader") + share("uorf"), 0.030); expect_lt(share("leader") + share("uorf"), 0.070)
  expect_gt(share("cds"), 0.903); expect_lt(share("cds"), 0.951)
  expect_gt(share("trailer"), 0.014); expect_lt(share("trailer"), 0.023)
})

test_that("default RFP read lengths peak at 30-31 nt", {
  dist <- default_fragment_distribution(25:36)
  weight_col <- setdiff(names(dist), c("fragment_length", "site_offset"))[1]
  w <- dist[[weight_col]]
  expect_true(dist$fragment_length[which.max(w)] %in% 30:31)
  expect_equal(measured_rfp_length_weights(c(30L, 36L)), c(.157, .022))
  expect_true(all(measured_rfp_length_weights(c(24L, 40L)) > 0))
})
