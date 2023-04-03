context("Unit tests")
library(coverageSim)

test_that("Loading sequence bias tables work", {
  dt <- load_seq_bias()
  expect_is(dt, "data.table")
  expect_equal(unique(dt$variable), "R2")
})

