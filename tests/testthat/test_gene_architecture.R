# Constructive gene architectures (RDG plan step 2): exact placement, sequence
# written for every frame at once, verification, decoy labelling, genome mapping.

stops <- c("TAA", "TAG", "TGA")

has_internal_stop <- function(x) vapply(x, function(z) {
  n <- nchar(z) %/% 3L
  if (n <= 2L) return(FALSE)
  any(substring(z, seq(4L, by = 3L, length.out = n - 2L),
                seq(6L, by = 3L, length.out = n - 2L)) %in% stops)
}, logical(1), USE.NAMES = FALSE)

build_genome <- function(architecture, n = 6, seed = 7, ...) {
  set.seed(seed)
  dir <- tempfile("arch_"); dir.create(dir)
  suppressMessages(simGenome(n = n, out_dir = dir, architecture = architecture,
                             leader_length = rep(150, n), cds_length = rep(300, n),
                             trailer_length = rep(150, n), debug_on = FALSE, ...))
}

feature_seqs <- function(genome) {
  ranges <- readRDS(genome[["features"]])
  as.character(ORFik::txSeqsFromFa(ranges, Rsamtools::FaFile(genome[["genome"]]),
                                   is.sorted = TRUE))
}

test_that("every requested feature exists exactly in the genome, on both strands and across introns", {
  genome <- build_genome("all", n = 8, cds_exons = 2)
  features <- data.table::fread(genome[["features_table"]])
  requested <- features$origin == "requested"
  expect_equal(sum(requested), 8 * 5)
  expect_setequal(unique(features$class[requested]),
                  c("uORF", "uoORF", "intORF", "doORF", "dORF"))
  seqs <- feature_seqs(genome)
  expect_equal(unname(nchar(seqs)), features$tx_end - features$tx_start + 1L)
  expect_true(all(substr(seqs, 1, 3)[requested] == features$start_codon[requested]))
  expect_true(all(substring(seqs, nchar(seqs) - 2L)[requested] %in% stops))
  expect_false(any(has_internal_stop(seqs[requested])))
  # Out-of-frame classes really are out of the CDS frame.
  oof <- requested & features$class %in% c("uoORF", "intORF", "doORF")
  expect_true(all(features$frame[oof] %in% 1:2))
  # The main CDS is still a valid ORF after all other frames were written.
  fa <- Rsamtools::FaFile(genome[["genome"]])
  cds <- as.character(ORFik::txSeqsFromFa(
    ORFik::loadRegion(ORFik::loadTxdb(genome[["txdb"]]), "cds"), fa, is.sorted = TRUE))
  expect_true(all(substr(cds, 1, 3) == "ATG"))
  expect_false(any(has_internal_stop(cds)))
  expect_true(all(substring(cds, nchar(cds) - 2L) %in% stops))
})

test_that("near-cognate N-terminal extensions run in frame to the CDS stop", {
  genome <- build_genome("NTE_CUG", n = 4)
  features <- data.table::fread(genome[["features_table"]])
  nte <- features[origin == "requested"]
  expect_true(all(nte$class == "NTE" & nte$frame == 0L & nte$start_codon == "CTG"))
  seqs <- feature_seqs(genome)[features$origin == "requested"]
  expect_true(all(startsWith(seqs, "CTG")))
  expect_false(any(has_internal_stop(seqs)))
  expect_equal(unname(nchar(seqs)), rep(60L + 306L, 4))
})

test_that("decoys are built untranslated and incidental ORFs are labelled", {
  genome <- build_genome("uORF_decoy", n = 4)
  features <- data.table::fread(genome[["features_table"]])
  requested <- features[origin == "requested"]
  expect_equal(as.vector(table(requested$status)), c(4L, 4L)) # decoy, translated
  incidental <- features[origin == "incidental"]
  expect_true(all(incidental$status == "decoy"))
  expect_true(all(incidental$class %in% c(feature_classes, "N_truncation")))
  # Only the translated uORFs reach the read simulation.
  uorfs <- readRDS(genome[["uorfs"]])
  expect_equal(length(uorfs), 4L)
  expect_setequal(names(uorfs), unique(requested$transcript_id))
})

test_that("a different architecture per gene, including none, is built as asked", {
  arch <- list(architecture_preset("two_uORFs"), NULL, "dORF")
  genome <- build_genome(arch, n = 3)
  features <- data.table::fread(genome[["features_table"]])[origin == "requested"]
  expect_equal(features[, .N, by = transcript_id][order(transcript_id)]$N, c(2L, 1L))
  expect_equal(sort(features$class), c("dORF", "uORF", "uORF"))
  two <- features[class == "uORF"]
  expect_true(two$tx_end[1] < two$tx_start[2] || two$tx_end[2] < two$tx_start[1])
})

test_that("impossible architectures stop before any sequence is made", {
  dir <- tempfile("arch_bad_"); dir.create(dir)
  expect_error(simGenome(n = 2, out_dir = dir, leader_length = c(150, 20),
                         architecture = gene_architecture("uORF", 30), debug_on = FALSE),
               "ENSTTEST10002: uORF of 30 nt does not fit")
  expect_false(file.exists(file.path(dir, "Homo_sapiens_dummy.fasta")))
  expect_error(simGenome(n = 1, out_dir = dir, leader_length = 40,
                         architecture = gene_architecture(c("uORF", "uORF"), c(30, 30),
                                                          may_overlap = FALSE),
                         debug_on = FALSE),
               "no room for uORF of 30 nt")
  expect_error(simGenome(n = 1, out_dir = dir, max_uorfs = 2, architecture = "uORF",
                         debug_on = FALSE), "not both")
  expect_error(gene_architecture("uORF", 31), "multiples of 3")
  expect_error(gene_architecture("xORF", 30), "Unknown feature class")
})

test_that("reads can be simulated on a genome built from an architecture", {
  genome <- build_genome("ATF4_like", n = 4)
  counts <- suppressMessages(simCountTables(
    ORFik::loadRegion(genome["txdb"], "cds"), libtypes = "RFP",
    conditions = "WT", replicates = 1, plot_PCA = FALSE, print_statistics = FALSE))
  SummarizedExperiment::assay(counts)[] <- 1000
  regions <- suppressMessages(simCountTablesRegions(counts))
  dir <- dirname(genome[["genome"]])
  experiment <- suppressMessages(simNGScoverage(genome, regions, out_dir = dir,
    exp_save_dir = dir, fragment_mode = "legacy_point",
    libFormats = list(RFP = "ofst"), validate = FALSE))
  reads <- ORFik::fimport(ORFik::filepath(experiment, "default")[[1]])
  expect_equal(sum(S4Vectors::mcols(reads)$score), 4 * 1000)
})
