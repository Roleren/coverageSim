# Graded ternary-complex re-acquisition (rdg_parameters(reinit_d0 = )).
flat <- function(...) rdg_parameters(conditions = list(), p_init_translon = NULL, uoorf_cds_factor = 1, ...)

test_that("reinit_d0 = NULL keeps the hard step: results frozen from the code before reinit_d0", {
  a <- rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 100L, C = 400L,
                      params = flat(p_init = c(ATG = 0.8, CTG = 0.1), reinit_max = 0.5, reinit_scale_codons = 1e9))
  expect_equal(a$cds_initiations, 0.8 * (0.9 + 0.05))
  expect_null(rdg_parameters()$reinit_d0)
  # values computed with rdg_flow.R of feature/rdg-architecture 9179e3d (before this option)
  cases <- list(list(c(5L, 60L), c(40L, 140L), c("uORF", "uoORF"), c("ATG", "CTG"), 100L, 700L),
                list(11L, 40L, "uORF", "CTG", 100L, 400L),
                list(c(1L, 100L), c(12L, 489L), c("uORF", "uoORF"), c("ATG", "ATG"), 399L, 1000L))
  frozen_default <- list(
    list(initiations = c(0.3, 0.21), load = c(10.8, 17.01), cds_initiations = 0.390820869303306, cds_load = 234.492521581983),
    list(initiations = 0.3, load = 9, cds_initiations = 0.663180508722626, cds_load = 198.954152616788),
    list(initiations = c(0.3, 0.257259359228319), load = c(3.6, 100.331150099044), cds_initiations = 0.384173976447623,
         cds_load = 230.888559845022))
  frozen_flat <- list(
    list(initiations = c(0.5, 0.0140219202762138), load = c(18, 1.13577554237332), cds_initiations = 0.343537046767239,
         cds_load = 206.122228060343),
    list(initiations = 0.02, load = 0.6, cds_initiations = 0.494299187863443, cds_load = 148.289756359033),
    list(initiations = c(0.5, 0.381275997856442), load = c(6, 148.697639164012), cds_initiations = 0.190637998928221,
         cds_load = 114.573437355861))
  pf <- rdg_parameters(p_init = c(ATG = 0.5), p_init_translon = NULL, uoorf_cds_factor = 1, reinit_distance = 20)
  for (i in seq_along(cases)) {
    expect_equal(do.call(rdg_gene_loads, c(cases[[i]], list(rdg_parameters()))), frozen_default[[i]], tolerance = 1e-12)
    expect_equal(do.call(rdg_gene_loads, c(cases[[i]], list(pf))), frozen_flat[[i]], tolerance = 1e-12)
  }
  # positional calls keep their meaning (the new argument is last)
  pp <- rdg_parameters(c(ATG = 0.5), 0.01, 0.4, 20, 30, 0.1)
  expect_equal(pp$reinit_after_cds, 0.1); expect_equal(pp$reinit_distance, 30); expect_null(pp$reinit_d0)
})

test_that("at exactly reinit_distance the step gives all, the graded reload none", {
  base <- list(p_init = c(ATG = 0.8, CTG = 0.1), reinit_max = 0.5, reinit_scale_codons = 1e9, reinit_distance = 61)
  step <- rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 100L, C = 400L, params = do.call(flat, base))  # CDS start 101 = 40 + 61
  expect_equal(step$cds_initiations, 0.8 * (0.9 + 0.05))
  for (d0 in c(1e-3, 1, 100)) {
    g <- rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 100L, C = 400L, params = do.call(flat, c(base, reinit_d0 = d0)))
    expect_equal(g$cds_initiations, 0.8 * 0.9)
  }
})

test_that("graded reload gives 1 - exp(-(d - reinit_distance) / d0) of the reinitiating ribosomes", {
  for (rd in c(0, 30)) for (d0 in c(20, 100)) {
    p <- flat(p_init = c(ATG = 0.8, CTG = 0.1), reinit_max = 0.5, reinit_scale_codons = 1e9,
              reinit_distance = rd, reinit_d0 = d0)
    x <- rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 100L, C = 400L, params = p)
    d <- 101 - 40                                             # CDS start 101, uORF stop ends at 40
    expect_equal(x$cds_initiations, 0.8 * (0.9 + 0.1 * 0.5 * -expm1(-(d - rd) / d0)))
  }
  # a start closer than reinit_distance gets no reinitiating ribosome
  p <- flat(p_init = c(ATG = 0.8, CTG = 0.1), reinit_max = 0.5, reinit_scale_codons = 1e9,
            reinit_distance = 50, reinit_d0 = 10)
  expect_equal(rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 60L, C = 360L, params = p)$cds_initiations, 0.8 * 0.9)
})

test_that("CDS initiation rises smoothly with distance and approaches the step for small d0", {
  cds_at <- function(L, d0) rdg_gene_loads(11L, 40L, "uORF", "ATG", L = L, C = L + 300L,
    params = flat(p_init = c(ATG = 0.5), reinit_max = 0.6, reinit_scale_codons = 1e9, reinit_distance = 0,
                  reinit_d0 = d0))$cds_initiations
  v <- vapply(c(45L, 60L, 100L, 200L, 400L), cds_at, 0, d0 = 100)
  expect_true(all(diff(v) > 0))
  expect_equal(cds_at(400L, 1e-6), rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 400L, C = 700L,
    params = flat(p_init = c(ATG = 0.5), reinit_max = 0.6, reinit_scale_codons = 1e9, reinit_distance = 0))$cds_initiations)
})

test_that("unloaded ribosomes skip a start and can load before the next one (memoryless)", {
  # uORF 1-30, a second uORF at 61 (d = 31), CDS at 201 (d = 171). Ribosomes not yet loaded at 61 keep loading.
  p <- flat(p_init = c(ATG = 1), reinit_max = 1, reinit_scale_codons = 1e9, reinit_distance = 0, reinit_d0 = 50)
  x <- rdg_gene_loads(c(1L, 61L), c(30L, 90L), c("uORF", "uORF"), c("ATG", "ATG"), L = 200L, C = 500L, params = p)
  f1 <- -expm1(-31 / 50)                                       # loaded by the 2nd uORF: all of them initiate there
  expect_equal(x$initiations, c(1, f1))
  # CDS: the rest of uORF1's ribosomes (memoryless, loaded between 61 and 201) + uORF2's reinitiating ones (d = 111)
  expect_equal(x$cds_initiations, (1 - f1) * -expm1(-140 / 50) + f1 * -expm1(-111 / 50))
})

test_that("loaded ribosomes that do not initiate at a start stay in the scanning flux", {
  # p = 0.4 at the second uORF: loaded ribosomes that leak past it reach the CDS together with those still loading.
  p <- flat(p_init = c(ATG = 1, CTG = 0.4), reinit_max = 1, reinit_scale_codons = 1e9, reinit_distance = 0,
            reinit_d0 = 50)
  x <- rdg_gene_loads(c(1L, 61L), c(30L, 90L), c("uORF", "uORF"), c("ATG", "CTG"), L = 200L, C = 500L, params = p)
  f1 <- -expm1(-31 / 50)
  ini2 <- 0.4 * f1
  expect_equal(x$initiations, c(1, ini2))
  expect_equal(x$cds_initiations, f1 - ini2 + (1 - f1) * -expm1(-140 / 50) + ini2 * -expm1(-111 / 50))
})

test_that("reinit_d0 is validated, also per condition, and a condition change counts as ISR-type", {
  expect_error(rdg_parameters(reinit_d0 = 0), "reinit_d0")
  expect_error(rdg_parameters(reinit_d0 = c(1, 2)), "reinit_d0")
  expect_error(rdg_parameters(reinit_d0 = NA_real_), "reinit_d0")
  expect_error(rdg_parameters(conditions = list(ISR = list(reinit_d0 = -1))), "reinit_d0")
  p <- rdg_parameters(conditions = list(ISR = list(reinit_d0 = 300)))
  expect_true(rdg_reinit_regulated(p))
  expect_equal(rdg_condition_parameters(p, "ISR")$reinit_d0, 300)
  expect_false(rdg_reinit_regulated(rdg_parameters(reinit_d0 = 100, conditions = list(B = list(reinit_d0 = 100)))))
  expect_false(rdg_reinit_regulated(rdg_parameters()))
  expect_false(rdg_reinit_regulated(rdg_parameters(reinit_d0 = 100, conditions = list(B = list(reinit_d0 = 100L)))))
  expect_true(rdg_reinit_regulated(rdg_parameters(reinit_d0 = 100, conditions = list(B = list(reinit_d0 = NULL)))))
})
