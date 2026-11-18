test_that("compiled stop matching preserves R results and order", {
  set.seed(147)
  bases <- sample(c("A", "C", "G", "T"), 10000, replace = TRUE)
  positions <- sample(seq_len(length(bases) - 2L), 3000, replace = TRUE)
  for (stops in list(c("TAA", "TAG", "TGA"), c("AAA", "CGT"), character())) {
    expected <- positions[paste0(bases[positions], bases[positions + 1L],
                                 bases[positions + 2L]) %in% stops]
    expect_identical(feature_stop_positions(bases, positions, stops), expected)
  }
  expect_identical(feature_stop_positions(bases, integer(), "TAA"), integer())
})
