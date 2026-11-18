test_that("compiled feature writer reproduces R draws and overlapping frames", {
  arch <- data.table::data.table(gene = c(1L, 1L, 2L), tx_start = c(4L, 13L, 4L),
    tx_end = c(30L, 39L, 30L), class = c("uORF", "uORF", "uORF"),
    start_codon = c("ATG", "CTG", "ATG"))
  tx <- c(paste(rep("C", 120), collapse = ""), paste(rep("A", 120), collapse = ""))
  compare <- function(engine, seed) {
    old <- options(coverageSim.engine = engine); on.exit(options(old))
    set.seed(seed)
    tryCatch(list(value = write_feature_sequences(tx, data.table::copy(arch), c(40L, 40L), c(60L, 60L), c("TAA", "TAG", "TGA")), seed = .Random.seed), error = function(e) conditionMessage(e))
  }
  for (seed in 1:20) expect_identical(compare("R", seed), compare("cpp", seed))
  arch <- arch[0]
  expect_identical(compare("R", 42), compare("cpp", 42))
})
