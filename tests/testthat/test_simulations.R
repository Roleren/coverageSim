context("Full simulations")
library(coverageSim)

test_that("Default genome simulation works", {
  simGenome6 <- simGenome(n = 6, max_uorfs = 0, cds_length = c(rep.int(330, 3), rep.int(300, 3)),
                          export_txdb = FALSE)
})

test_that("uORF (n=1) genome simulation works", {
  simGenome_uORFs_1 <- simGenome(n = 6, max_uorfs = 1, cds_length = c(rep.int(330, 3), rep.int(300, 3)),
                                 export_txdb = FALSE)
})

test_that("uORF (n=2) genome simulation works", {
  simGenome_uORFs_2 <- simGenome(n = 6, max_uorfs = 2, cds_length = c(rep.int(330, 3), rep.int(300, 3)),
                                 export_txdb = FALSE)
})

test_that("Default gene count table simulation (RFP) works", {
  gene_count_table <-simCountTables(loadRegion(simGenome6["txdb"], "cds"),
                                    libtypes = "RFP", print_statistics = FALSE, interceptMean = 10)
})

test_that("Default gene region count table simulation (RFP) works", {
  region_count_table <- simCountTablesRegions(gene_count_table,
                                              regionsToSample = c("cds"))
})

test_that("Default read coverage simulation (RFP) works", {
  simNGScoverage(simGenome6, region_count_table[,1])
})


