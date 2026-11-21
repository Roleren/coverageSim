test_that("codon dwell defaults are off and options validate", {
  withr::local_options(coverageSim.codon_dwell = NULL)
  o <- codon_dwell_options()
  expect_equal(o$global, "none"); expect_false(o$study_dwell)
  expect_false(codon_dwell_active(o, S4Vectors::DataFrame(libtype = "RFP", study = "S1")))
  expect_true(codon_dwell_active(o, S4Vectors::DataFrame(tissue_group = "striated")))
  expect_null(codon_dwell_run_setup(S4Vectors::DataFrame(libtype = c("RFP", "RFP")), c("a", "b")))
  expect_equal(library_dwell_for_region(NULL, NULL, "a"), list(factors = NULL, power = 1, var_ratio = 1))
  withr::local_options(coverageSim.codon_dwell = list(oops = 1)); expect_error(codon_dwell_options(), "unknown")
  withr::local_options(coverageSim.codon_dwell = list(between_r = 0.9)); expect_error(codon_dwell_options())
  tab <- codon_dwell_table()
  expect_equal(nrow(tab), 16L * 61L)
  expect_equal(sort(unique(tab$offset)), -10:5)
  # H1 reference codon is 0 at every offset; A-site acidic codons (E, D) are slow in H1 (E13-E)
  expect_true(all(tab[codon == "AAA", h1_beta] == 0))
  expect_gt(tab[offset == 0 & codon %in% c("GAA", "GAC", "GAT"), min(h1_beta)], 0.3)
  # striated group: all four Pro A-site codons slower by 0.10-0.24, CAG faster
  expect_true(all(tab[offset == 0 & codon %in% c("CCA", "CCC", "CCG", "CCT"), striated_delta] >= 0.09))
  expect_lt(tab[offset == 0 & codon == "CAG", striated_delta], -0.1)
})

test_that("study and library deviations reproduce the human within / between-study profile correlation", {
  withr::local_options(coverageSim.codon_dwell = list(global = "h1", study_dwell = TRUE))
  o <- codon_dwell_options()
  studies <- rep(paste0("S", 1:40), each = 4)
  P <- vapply(seq_along(studies), function(i) library_dwell_profile(studies[i], paste0("L", i), NA, o)$log_effect,
              numeric(976))
  C <- stats::cor(P); same <- outer(studies, studies, "=="); ut <- upper.tri(C)
  expect_equal(stats::median(C[same & ut]), 0.88, tolerance = 0.03)
  expect_equal(stats::median(C[!same & ut]), 0.37, tolerance = 0.12)
  # draws depend on the names only (not on order) and leave the session RNG alone
  set.seed(3); a <- library_dwell_profile("S7", "L9", NA, o); r1 <- stats::runif(1)
  set.seed(3); library_dwell_profile("S1", "L1", NA, o); b <- library_dwell_profile("S7", "L9", NA, o)
  expect_equal(a, b); expect_equal({set.seed(3); library_dwell_profile("S7", "L9", NA, o); stats::runif(1)}, r1)
  # tissue delta is added on top
  s <- library_dwell_profile("S7", "L9", "striated", o)
  tab <- codon_dwell_table()
  expect_equal(s$log_effect - a$log_effect, tab$striated_delta)
  expect_error(library_dwell_profile("S7", "L9", "liver", o), "unknown tissue_group")
})

test_that("dwell factors use offsets relative to the A site, shifted for P-site placement", {
  profile <- data.table::data.table(offset = c(0L, 1L), codon = c("GAA", "TTT"), log_effect = log(c(2, 3)))
  codons <- list(c("ATG", "GAA", "TTT", "AAA"))
  # A-site placement: position i is the A site, offset +1 is the next codon
  expect_equal(dwell_codon_factors(codons, profile, 0L)[[1]], c(1, 2 * 3, 1, 1))
  # P-site placement: the A site of position i is codon i + 1
  expect_equal(dwell_codon_factors(codons, profile, 1L)[[1]], c(2 * 3, 1, 1, 1))
})

test_that("codon factors and the peaky power act on the codon weights; peaky libraries get 1.5x top-10 share", {
  withr::local_options(coverageSim.codon_bias = codon_bias_preset("realistic"))
  tab <- load_seq_bias(type = "codon", shift = "a-site")
  sense <- setdiff(tab$seqs, c("###", "%%%", "&&&", "***", "TAA", "TAG", "TGA"))
  set.seed(3)
  genes <- lapply(1:12, function(i) sample(sense, 400, TRUE))
  alpha <- lapply(genes, function(g) tab$alpha[match(g, tab$seqs)])
  nt <- lengths(alpha) * 3
  base <- sim_sequence_bias(NULL, nt, alpha, shapes(9))
  expect_equal(sim_sequence_bias(NULL, nt, alpha, shapes(9), codon_factors = lapply(alpha, function(a) rep(1, length(a)))),
               base)
  doubled <- sim_sequence_bias(NULL, nt, alpha, shapes(9),
                               codon_factors = lapply(genes, function(g) ifelse(g == "GAA", 4, 1)))
  site <- function(v, g) { m <- (length(v) - 3 * length(g)) %/% 2; colSums(matrix(v[m + seq_len(3 * length(g))], 3)) }
  ratio <- mapply(function(d, b, g) { x <- site(d, g) / site(b, g); mean(x[g == "GAA"]) / mean(x[g != "GAA"]) },
                  doubled, base, genes)
  expect_gt(stats::median(ratio), 2.5)
  expect_error(sim_sequence_bias(NULL, nt, alpha, shapes(9), codon_factors = list(1)), "do not match")
  top10 <- function(power, ratio) {
    set.seed(9)
    res <- scale_dmn_alpha(sim_sequence_bias(NULL, nt, alpha, shapes(9), peaky_power = power), 6 / ratio)
    mean(vapply(seq_along(res), function(i) mean(replicate(10, {
      cc <- site(stats::rmultinom(1, 5000, draw_site_probabilities(res[[i]], dirichlet = TRUE))[, 1], genes[[i]])
      sum(sort(cc, TRUE)[1:10]) / sum(cc)
    })), numeric(1)))
  }
  p <- codon_dwell_preset("human")
  expect_equal(top10(p$peaky_power, p$peaky_var_ratio) / top10(1, 1), 1.5, tolerance = 0.1)
})

test_that("simNGScoverage applies the library dwell profile per library and writes it", {
  set.seed(77)
  genome_dir <- tempfile("dwell-genome-"); exp_dir <- tempfile("dwell-exp-"); out_dir <- tempfile("dwell-out-")
  dir.create(genome_dir); dir.create(exp_dir); dir.create(out_dir)
  sg <- suppressWarnings(simGenome(n = 4, out_dir = genome_dir, max_uorfs = 0, cds_length = rep(1500, 4),
                                   export_txdb = TRUE))
  cds <- ORFik::loadRegion(sg["txdb"], "cds")
  ct <- simCountTables(cds, libtypes = "RFP", print_statistics = FALSE, plot_PCA = FALSE, interceptMean = 11,
                       replicates = 1)
  rt <- simCountTablesRegions(ct, regionsToSample = "cds")
  rt <- assign_dwell_groups(rt, study = "S1", tissue_group = c(NA, "test"))
  withr::local_options(coverageSim.codon_bias = codon_bias_preset("realistic"),
                       coverageSim.codon_dwell = list(tissue_deltas = list(
                         test = data.frame(offset = 0L, codon = "GAA", delta = log(8)))))
  experiment <- simNGScoverage(sg, rt, out_dir = out_dir, exp_name = basename(out_dir), exp_save_dir = exp_dir,
                               validate = FALSE, ground_truth = TRUE, rnase_bias = list(RFP = NULL))
  prof <- data.table::fread(file.path(out_dir, "codon_dwell_profiles.tsv"))
  expect_equal(nrow(prof), 2L * 976L)
  expect_equal(prof[offset == 0 & codon == "GAA", log_effect], c(0, log(8)))
  truth <- lapply(colnames(rt), function(l) data.table::fread(file.path(out_dir, paste0(l, "_ground_truth.tsv"))))
  share <- vapply(truth, function(t) t[codon == "GAA", sum(score)] / sum(t$score), numeric(1))
  odds <- share / (1 - share)
  expect_gt(odds[2] / odds[1], 4)
  expect_lt(odds[2] / odds[1], 12)
})
