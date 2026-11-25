# Ribosome flow on the RDG (plan step 4): analytic cases and the coupled read split.

# Analytic cases of the scanning model with the codon table for every start (p_init_translon = NULL, uoorf_cds_factor = 1).
flat <- function(...) rdg_parameters(conditions = list(), p_init_translon = NULL, uoorf_cds_factor = 1, ...)

test_that("a uORF that captures every ribosome leaves none for the CDS", {
  p <- flat(p_init = c(ATG = 1), reinit_max = 0)
  x <- rdg_gene_loads(starts = 11L, ends = 40L, classes = "uORF", start_codons = "ATG",
                      L = 100L, C = 400L, params = p)
  expect_equal(x$initiations, 1)
  expect_equal(x$cds_initiations, 0)
})

test_that("leaky scanning past a uORF gives the CDS (1 - p) * p", {
  p <- flat(p_init = c(ATG = 0.8), reinit_max = 0)
  x <- rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L, params = p)
  expect_equal(x$cds_initiations, 0.2 * 0.8)
  expect_equal(x$load, 0.8 * 30)
  expect_equal(x$cds_load, 0.16 * 300)
})

test_that("near-cognate starts initiate less, and reinitiation needs distance", {
  p <- flat(p_init = c(ATG = 0.8, CTG = 0.1), reinit_max = 0.5, reinit_scale_codons = 1e9,
            reinit_distance = 50)
  x <- rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 100L, C = 400L, params = p)
  expect_equal(x$initiations, 0.1)
  # 0.9 leak + 0.1 * 0.5 reinitiating (stop at 40 + 50 = 90 < CDS start 101).
  expect_equal(x$cds_initiations, 0.8 * (0.9 + 0.05))
  far <- rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 60L, C = 360L, params = p)
  expect_equal(far$cds_initiations, 0.8 * 0.9) # rejoins at 90, after the CDS start at 61
})

test_that("ATF4-like: under ISR ribosomes skip the uoORF and reach the CDS", {
  # uORF1 1-12, uoORF 100-489 (ends inside the CDS), CDS 400-1000.
  args <- list(starts = c(1L, 100L), ends = c(12L, 489L), classes = c("uORF", "uoORF"),
               start_codons = c("ATG", "ATG"), L = 399L, C = 1000L)
  params <- rdg_parameters(reinit_distance = 50, p_init_translon = NULL, uoorf_cds_factor = 1,
                           conditions = list(ISR = list(reinit_distance = 200)))
  wt <- do.call(rdg_gene_loads, c(args, list(params = rdg_condition_parameters(params, "WT"))))
  isr <- do.call(rdg_gene_loads, c(args, list(params = rdg_condition_parameters(params, "ISR"))))
  expect_gt(wt$initiations[2], isr$initiations[2])
  expect_gt(isr$cds_initiations, 2 * wt$cds_initiations)
})

test_that("simRDGflow keeps each gene's reads and the read simulation follows its shares", {
  set.seed(4)
  n <- 6
  dir <- tempfile("rdg_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = n, out_dir = dir, architecture = "all",
    leader_length = rep(200, n), cds_length = rep(600, n), trailer_length = rep(300, n),
    debug_on = FALSE))
  counts <- suppressMessages(simCountTables(
    ORFik::loadRegion(genome["txdb"], "cds"), libtypes = "RFP",
    conditions = c("WT", "ISR"), replicates = 1, plot_PCA = FALSE, print_statistics = FALSE))
  SummarizedExperiment::assay(counts)[] <- 4000
  regions <- suppressMessages(simCountTablesRegions(counts))
  flowed <- simRDGflow(regions, genome)
  for (a in c("gene", "leader", "trailer"))
    expect_equal(SummarizedExperiment::assay(flowed, a), SummarizedExperiment::assay(regions, a))
  translated <- function(x) SummarizedExperiment::assay(x, "cds") + SummarizedExperiment::assay(x, "uorf")
  expect_equal(translated(flowed), translated(regions))
  share <- S4Vectors::metadata(flowed)$rdg_feature_share
  expect_equal(dim(share), c(n * 5L, 2L))
  expect_true(all(abs(tapply(share[, 1], rep(seq_len(n), each = 5), sum) - 1) < 1e-9))
  suppressMessages(simNGScoverage(genome, flowed, out_dir = dir, exp_save_dir = dir,
    ground_truth = TRUE, libFormats = list(RFP = "ofst"), validate = FALSE))
  truth_file <- list.files(dir, "RFP_WT_1_ground_truth", recursive = TRUE, full.names = TRUE)
  truth <- annotate_ground_truth(truth_file, genome)
  reads <- truth[!is.na(feature_id), .(reads = sum(score)), by = sampling_group]
  expected <- share[reads$sampling_group, 1]
  # Translons with a larger share get more reads.
  expect_gt(cor(reads$reads, expected), 0.5)
})

test_that("frameshift and readthrough loads come from CDS ribosomes", {
  p <- flat(p_init = c(ATG = 0.5), p_readthrough = 0.1, p_frameshift = 0.2)
  x <- rdg_gene_loads(starts = c(1001L, 301L), ends = c(1060L, 450L), classes = c("CTE", "FS"),
                      start_codons = c(NA, NA), L = 100L, C = 1000L, params = p,
                      efficiency = c(NA, 0.3))
  expect_equal(x$cds_initiations, 0.5)
  expect_equal(x$initiations, c(0.5 * 0.1, 0.5 * 0.3))
  expect_equal(x$load, c(0.05 * 60, 0.15 * 150))
  # Shifted ribosomes leave the CDS frame at the shift site.
  expect_equal(x$cds_load, 0.5 * 900 - 0.15 * (1000 - 301 + 1))
})

test_that("simRDGflow and simNGScoverage work on gene and library subsets", {
  set.seed(8)
  n <- 6
  dir <- tempfile("rdg_subset_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = n, out_dir = dir, architecture = "uORF",
    leader_length = rep(150, n), cds_length = rep(300, n), trailer_length = rep(150, n),
    debug_on = FALSE))
  counts <- suppressMessages(simCountTables(ORFik::loadRegion(genome["txdb"], "cds"),
    libtypes = "RFP", conditions = c("WT", "ISR"), replicates = 1, plot_PCA = FALSE,
    print_statistics = FALSE))
  SummarizedExperiment::assay(counts)[] <- 2000
  regions <- suppressMessages(simCountTablesRegions(counts))
  sub_genes <- simRDGflow(regions[1:4, ], genome)
  share <- S4Vectors::metadata(sub_genes)$rdg_feature_share
  expect_true(all(is.na(share[!features_in_genes(genome, rownames(regions)[1:4]), ])))
  flowed <- simRDGflow(regions, genome)
  expect_no_error(suppressMessages(simNGScoverage(genome, flowed[, 1], out_dir = dir,
    exp_save_dir = dir, fragment_mode = "legacy_point", libFormats = list(RFP = "ofst"),
    validate = FALSE)))
})

test_that("the default parameters treat every condition alike", {
  expect_length(rdg_parameters()$conditions, 0)
})
