# Difficulty presets (Carina item 1) and the cheap known-limitation fixes (item 4): nuclease cut-site context of the
# end bias (flank), uORF / CDS scanning coupling.

codon_bias_options_from <- function(p) withr::with_options(p$options, codon_bias_options())

test_that("the three presets exist, validate and order their difficulty knobs", {
  ps <- lapply(c(clean = "clean", realistic = "realistic", noisy = "noisy"), simulation_preset)
  for (p in ps) {
    expect_s3_class(p, "simulation_preset")
    expect_true(all(c("options", "genome", "count_tables", "regions", "coverage") %in% names(p)))
    withr::with_options(p$options, {
      expect_silent(o <- codon_bias_options())
      expect_true(is.list(utr3_options()))
    })
    r <- p$regions$region_proportion
    expect_equal(sum(vapply(r, function(x) x$RFP, numeric(1))), 1, tolerance = 1e-9)
    d <- normalize_fragment_geometry(p$coverage$fragment_geometry)$distribution
    expect_equal(sum(d$probability), 1, tolerance = 1e-6)
  }
  expect_identical(ps$realistic$options$coverageSim.codon_bias, codon_bias_preset("human_v2"))
  # CDS share, gene dispersion, DMN concentration, multimapping: clean > realistic > noisy (or reverse)
  cds <- vapply(ps, function(p) p$regions$region_proportion$cds$RFP, numeric(1))
  expect_true(cds[["clean"]] > cds[["realistic"]] && cds[["realistic"]] > cds[["noisy"]])
  disp <- vapply(ps, function(p) p$count_tables$dispMeanRel(1e6), numeric(1))
  expect_true(disp[["clean"]] < disp[["realistic"]] && disp[["realistic"]] < disp[["noisy"]])
  alpha <- vapply(ps, function(p) p$coverage$dmn_alpha_scale, numeric(1))
  expect_true(alpha[["clean"]] > alpha[["realistic"]] && alpha[["realistic"]] > alpha[["noisy"]])
  mm <- vapply(ps, function(p) normalize_technical_artifacts(p$coverage$technical_artifacts)$multimapping_rate, numeric(1))
  expect_true(mm[["clean"]] < mm[["realistic"]] && mm[["realistic"]] < mm[["noisy"]])
  expect_equal(normalize_technical_artifacts(ps$clean$coverage$technical_artifacts)$spike_sites, 0L)
  expect_gt(normalize_technical_artifacts(ps$noisy$coverage$technical_artifacts)$spike_sites, 0L)
  expect_false(codon_bias_options_from(ps$clean)$field$active)
  expect_true(codon_bias_options_from(ps$noisy)$field$gc_bias)
  expect_output(print(ps$noisy), "noisy")
  expect_error(simulation_preset("bogus"))
})

test_that("clean geometry has one offset per length, noisy geometry spreads offsets and adds short reads", {
  geo <- realistic_rfp_settings()$fragment_geometry$distribution
  one <- single_offset_geometry(geo)
  expect_equal(nrow(one), data.table::uniqueN(geo$fragment_length))
  expect_equal(sum(one$probability), 1, tolerance = 1e-9)
  # length shares are kept
  lp <- data.table::as.data.table(geo)[, .(p = sum(probability)), by = fragment_length]
  expect_equal(one$probability, lp$p / sum(lp$p), tolerance = 1e-9)
  wide <- spread_offset_geometry(geo, floor = 0.6, short_share = 0.12)
  expect_equal(sum(wide$probability), 1, tolerance = 1e-9)
  expect_equal(sum(wide[fragment_length < 25, probability]), 0.12, tolerance = 1e-9)
  r <- wide[fragment_length >= 25, .(ratio = min(probability) / max(probability), n = .N), by = fragment_length]
  expect_true(all(r$ratio >= 0.6 - 1e-9))
  expect_true(all(r$n >= 3))
})

test_that("noisy codon preset is human_v2 plus documented extra noise; clean has no pauses", {
  v2 <- codon_bias_preset("human_v2"); nz <- codon_bias_preset("noisy"); cl <- codon_bias_preset("clean")
  expect_equal(nz$field$field_var, 2 * v2$field$field_var)
  expect_true(nz$field$gc_bias); expect_equal(nz$zero_prob, 0.3)
  expect_equal(cl$pause_sd, 0); expect_null(cl$field)
  withr::with_options(list(coverageSim.codon_bias = nz), expect_true(codon_bias_options()$field$active))
})

test_that("run_simulation_preset runs a tiny clean simulation and restores options", {
  skip_on_cran()
  withr::local_options(coverageSim.codon_bias = NULL)
  out <- tempfile("preset-"); on.exit(unlink(out, recursive = TRUE))
  e <- suppressWarnings(suppressMessages(run_simulation_preset("clean", n = 6, out_dir = out, replicates = 1,
                                                               interceptMean = 2, seed = 3)))
  expect_null(getOption("coverageSim.codon_bias"))
  f <- ORFik::filepath(e, "default")
  expect_length(f, 1L)
  expect_gt(sum(S4Vectors::mcols(ORFik::fimport(f[1]))$score), 0)
})

test_that("end bias flank = 0 is the old behaviour; flank > 0 reads the bases across the cut", {
  geometry <- function(flank) list(source = "user",
    distribution = data.frame(fragment_length = 8L, site_offset = c(1L, 2L), probability = c(0.5, 0.5)),
    five_prime_bias = list(source = "user", flank = flank,
                           table = data.table::data.table(kmer = "GA", weight = 4)))
  # fixture: poly-A with G at 28; signal at 30. offset 1 starts at 29 (base before the cut = G),
  # offset 2 starts at 28 (first fragment base = G).
  fx <- end_selection_fixture(); fx$signal <- fx$signal[1]
  share_offset1 <- function(flank) {
    set.seed(5)
    r <- make_simulated_rpf_fragments(fx$signal, fx$models, 8L, geometry(flank))
    r[site_offset == 1L, sum(score)] / sum(r$score)
  }
  expect_equal(share_offset1(0L), 0.2, tolerance = 0.02)   # inside k-mer "GA" at offset 2 is favoured
  expect_equal(share_offset1(1L), 0.8, tolerance = 0.02)   # cut-site k-mer "G|A" at offset 1 is favoured
  r <- make_simulated_rpf_fragments(fx$signal, fx$models, 8L, geometry(1L))
  expect_true(all(r[site_offset == 1L, five_prime_kmer] == "GA"))
  expect_true(all(r[site_offset == 2L, five_prime_kmer] == "AG"))
  # three prime end: flank bases come from downstream of the fragment; beyond the transcript is neutral ("N")
  g <- normalize_fragment_geometry(list(three_prime_bias = list(source = "user", flank = 2L,
                                        table = data.table::data.table(kmer = "AAA", weight = 2))))
  expect_equal(end_bias_context_kmer("ACGTACGTAC", 3L, 6L, g$three_prime_bias, TRUE), "CGT")
  expect_equal(end_bias_context_kmer("ACGTACGTAC", 3L, 9L, g$three_prime_bias, TRUE), "ACN")
  expect_error(normalize_fragment_geometry(list(five_prime_bias = list(source = "user", flank = 2L,
                 table = data.table::data.table(kmer = "GA", weight = 2)))), "flank")
  expect_error(normalize_fragment_geometry(list(five_prime_bias = list(source = "user", flank = -1,
                 table = data.table::data.table(kmer = "GA", weight = 2)))), "flank")
})

test_that("scan_coupling = 0 makes uORF translation independent of the CDS; default 1 keeps leaky scanning", {
  p1 <- rdg_parameters(p_init = c(ATG = 0.8), reinit_max = 0)
  p0 <- rdg_parameters(p_init = c(ATG = 0.8), reinit_max = 0.5, scan_coupling = 0)
  ph <- rdg_parameters(p_init = c(ATG = 0.8), reinit_max = 0, scan_coupling = 0.5)
  x1 <- rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L, params = p1)
  x0 <- rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L, params = p0)
  xh <- rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L, params = ph)
  expect_equal(x1$cds_initiations, 0.2 * 0.8)
  expect_equal(x0$cds_initiations, 0.8)                 # no suppression, no extra reinitiation flux
  expect_equal(x0$load, 0.8 * 30)                       # uORF load unchanged
  expect_equal(xh$cds_initiations, (1 - 0.5 * 0.8) * 0.8)
  expect_error(rdg_parameters(scan_coupling = 1.5))
})
