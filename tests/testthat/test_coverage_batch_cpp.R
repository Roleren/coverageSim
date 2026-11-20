test_that("batched uniform multinomial draws preserve R draws and padding", {
  widths <- c(1L, 3L, 100L, 4L)
  budgets <- c(5, 0, 103, 21)
  for (padding in c(0L, 3L)) for (seed in 1:5) {
    set.seed(seed)
    actual <- coverageSim:::sample_uniform_regions_cpp(widths, budgets, padding)
    actual_seed <- .Random.seed
    set.seed(seed)
    expected <- unlist(lapply(seq_along(widths), function(i) {
      c(rep(0L, padding), as.vector(rmultinom(1, budgets[i], rep(1, widths[i]))),
        rep(0L, padding))
    }), use.names = FALSE)
    expect_identical(actual, expected)
    expect_identical(actual_seed, .Random.seed)
  }
})
