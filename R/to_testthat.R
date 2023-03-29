# debug(create_cds_seq)
# # Simulate Ribo-seq only (only cds)
# simGenome6 <- simGenome(n = 6, max_uorfs = 0, cds_length = c(rep.int(330, 3), rep.int(300, 3)))
# gene_count_table <-simCountTables(loadRegion(simGenome6["txdb"], "cds"),
#                                   libtypes = "RFP", print_statistics = FALSE, interceptMean = 10)
# region_count_table <- simCountTablesRegions(gene_count_table,
#                                             regionsToSample = c("cds"))
# simNGScoverage(simGenome6, region_count_table[,1])
#
# # Simulate Ribo-seq only with (cds, leader, trailers)
# simGenome6 <- simGenome(n = 6, max_uorfs = 0, cds_length = c(rep.int(330, 3), rep.int(300, 3)))
# gene_count_table <-simCountTables(loadRegion(simGenome6["txdb"], "cds"),
#                                   libtypes = "RFP", print_statistics = FALSE, interceptMean = 9)
# region_count_table <- simCountTablesRegions(gene_count_table,
#                                             regionsToSample = c("leader", "cds", "trailer"))
# simNGScoverage(simGenome6, region_count_table[,1])
# # Simulate Ribo-seq only with (cds, leader, trailers, uorfs)
# simGenome6 <- simGenome(n = 6, max_uorfs = 2, cds_length = c(rep.int(330, 3), rep.int(300, 3)))
# gene_count_table <-simCountTables(loadRegion(simGenome6["txdb"], "cds"),
#                                   libtypes = c("RFP", "RNA"), print_statistics = FALSE, interceptMean = 9); assay(gene_count_table)
# region_count_table <- simCountTablesRegions(gene_count_table)
# simNGScoverage(simGenome6, region_count_table[,c(1)], debug_coverage = F)
#
# #simNGScoverage(simGenome6, region_count_table[,1], debug_coverage = T)
# df <- read.experiment("simulated_data")
# uorfs <- readRDS(simGenome6["uorfs"])
