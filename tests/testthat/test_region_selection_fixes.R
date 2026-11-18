# Fixes from the 2026-11-18 review: region subsets must keep every read, and
# simulations without uORFs or without a CDS region must run.

small_rfp_counts <- function(genome, reads = 1000) {
  counts <- suppressMessages(simCountTables(
    ORFik::loadRegion(genome["txdb"], "cds"), libtypes = "RFP",
    conditions = "WT", replicates = 1, plot_PCA = FALSE, print_statistics = FALSE))
  SummarizedExperiment::assay(counts)[] <- reads
  counts
}

test_that("a region subset keeps every read of the gene", {
  set.seed(2026)
  dir <- tempfile("regions_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = 3, out_dir = dir, max_uorfs = 0,
                                       debug_on = FALSE))
  counts <- small_rfp_counts(genome)
  regions <- suppressMessages(simCountTablesRegions(
    counts, regionsToSample = c("leader", "cds", "trailer")))
  per_region <- Reduce(`+`, lapply(c("leader", "cds", "trailer"), function(r)
    SummarizedExperiment::assay(regions, r)))
  expect_equal(per_region, SummarizedExperiment::assay(regions, "gene"))
  expect_message(simCountTablesRegions(counts, regionsToSample = c("leader", "cds", "trailer")),
                 "rescaled")
})

test_that("selected regions with proportion 0 for a library type stop with an error", {
  set.seed(1)
  dir <- tempfile("regions_zero_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = 3, out_dir = dir, max_uorfs = 0,
                                       debug_on = FALSE))
  counts <- small_rfp_counts(genome)
  expect_error(simCountTablesRegions(counts, regionsToSample = "leader",
    region_proportion = list(leader = list(RFP = 0), cds = list(RFP = 1))),
    "nowhere to go")
  expect_error(simCountTablesRegions(counts, regionsToSample = "trailer",
    region_proportion = list(cds = list(RFP = 1))), "without a region_proportion")
})

test_that("a genome without uORFs works with the default region table", {
  set.seed(123)
  dir <- tempfile("no_uorf_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = 4, out_dir = dir, max_uorfs = 0,
                                       debug_on = FALSE))
  counts <- small_rfp_counts(genome)
  regions <- suppressMessages(simCountTablesRegions(counts))
  expect_message(experiment <- simNGScoverage(genome, regions, out_dir = dir,
    exp_save_dir = dir, fragment_mode = "legacy_point",
    libFormats = list(RFP = "ofst"), validate = FALSE), "No uORF annotation")
  reads <- ORFik::fimport(ORFik::filepath(experiment, "default")[[1]])
  uorf_reads <- sum(SummarizedExperiment::assay(regions, "uorf"))
  expect_equal(sum(S4Vectors::mcols(reads)$score), 4 * 1000 - uorf_reads)
})

test_that("a leader-only simulation runs and keeps all reads", {
  set.seed(123)
  dir <- tempfile("leader_only_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = 4, out_dir = dir, max_uorfs = 0,
                                       debug_on = FALSE))
  counts <- small_rfp_counts(genome)
  regions <- suppressMessages(simCountTablesRegions(counts, regionsToSample = "leader",
    region_proportion = list(leader = list(RFP = 1))))
  experiment <- suppressMessages(simNGScoverage(genome, regions, out_dir = dir,
    exp_save_dir = dir, fragment_mode = "legacy_point",
    libFormats = list(RFP = "ofst"), validate = FALSE))
  reads <- ORFik::fimport(ORFik::filepath(experiment, "default")[[1]])
  expect_equal(sum(S4Vectors::mcols(reads)$score), 4 * 1000)
})
