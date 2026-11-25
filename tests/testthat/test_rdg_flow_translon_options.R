# uORF -> CDS coupling options fitted on 62 human marker genes with detected uORFs (carina-claude, 2026-11-25):
# p_init_translon (default 0.3), translon_classes, uoorf_cds_factor (default 0.8), translon_gene_sd (default 0, opt-in),
# translon_exempt (default "auto": genes with uORF + uoORF in ISR-type scenarios keep the codon table).
# The tests pin both the new defaults and the old behaviour (p_init_translon = NULL, uoorf_cds_factor = 1).
old_params <- function(...) rdg_parameters(p_init_translon = NULL, uoorf_cds_factor = 1, ...)

test_that("defaults: p_init_translon 0.3 and uoorf_cds_factor 0.8 on, translon_gene_sd off; inputs validated", {
  p <- rdg_parameters()
  expect_equal(p$p_init_translon, 0.3); expect_equal(p$uoorf_cds_factor, 0.8); expect_equal(p$translon_gene_sd, 0)
  expect_equal(p$translon_classes, c("uORF", "uoORF", "NTE"))
  o <- old_params()
  expect_null(o$p_init_translon); expect_equal(o$uoorf_cds_factor, 1)
  expect_error(rdg_parameters(p_init_translon = 1.2), "within \\[0, 1\\]")
  expect_error(rdg_parameters(p_init_translon = c(0.2, 0.3)), "one probability")
  expect_error(rdg_parameters(translon_gene_sd = -1), ">= 0")
  expect_error(rdg_parameters(uoorf_cds_factor = 0), "\\(0, 1\\]")
})

test_that("old behaviour: leader translons use the codon table, no extra CDS repression", {
  p <- old_params(p_init = c(ATG = 0.8, CTG = 0.1), reinit_max = 0)
  x <- rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 100L, C = 400L, params = p)
  expect_equal(x$initiations, 0.1)
  expect_equal(x$cds_initiations, 0.8 * 0.9)
  args <- list(starts = 60L, ends = 150L, classes = "uoORF", start_codons = "ATG", L = 100L, C = 400L)
  y <- do.call(rdg_gene_loads, c(args, list(params = old_params(p_init = c(ATG = 0.8), reinit_max = 0))))
  expect_equal(y$initiations, 0.8)
  expect_equal(y$cds_initiations, 0.8 * 0.2)
})

test_that("new defaults: leader translons start at 0.3 whatever their codon; the CDS keeps its ATG value", {
  p <- rdg_parameters(p_init = c(ATG = 0.8, CTG = 0.1), reinit_max = 0)
  x <- rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 100L, C = 400L, params = p)
  expect_equal(x$initiations, 0.3)                  # not the CTG value 0.1
  expect_equal(x$cds_initiations, 0.8 * 0.7)        # no uoORF, so no extra repression
  args <- list(starts = 60L, ends = 150L, classes = "uoORF", start_codons = "ATG", L = 100L, C = 400L)
  y <- do.call(rdg_gene_loads, c(args, list(params = rdg_parameters(p_init = c(ATG = 0.8), reinit_max = 0))))
  expect_equal(y$initiations, 0.3)
  expect_equal(y$cds_initiations, 0.8 * 0.7 * 0.8)  # one overlapping uoORF: x 0.8
})

test_that("gene_factor scales the extra features' initiation, capped at 1", {
  p <- rdg_parameters(p_init = c(ATG = 0.8), reinit_max = 0, p_init_translon = 0.3)
  x <- rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L, params = p, gene_factor = 2)
  expect_equal(x$initiations, 0.6)
  y <- rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L, params = p, gene_factor = 10)
  expect_equal(y$initiations, 1)
})

test_that("uoorf_cds_factor represses the CDS once per overlapping uoORF", {
  args <- list(starts = c(11L, 60L), ends = c(40L, 150L), classes = c("uORF", "uoORF"),
               start_codons = c("ATG", "ATG"), L = 100L, C = 400L)
  base <- do.call(rdg_gene_loads, c(args, list(params = old_params(reinit_max = 0))))
  half <- rdg_parameters(p_init_translon = NULL, uoorf_cds_factor = 0.5, reinit_max = 0)
  rep <- do.call(rdg_gene_loads, c(args, list(params = half)))
  expect_equal(rep$cds_initiations, base$cds_initiations * 0.5)
  expect_equal(rep$initiations, base$initiations)  # the translons themselves are unchanged
  # a gene with only a leader uORF is not affected
  u <- rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L, params = half)
  expect_equal(u$cds_initiations, rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L,
                                                 params = old_params(reinit_max = 0))$cds_initiations)
})

test_that("p_init_translon applies only to translon_classes (default leader classes)", {
  d <- rdg_gene_loads(500L, 560L, "dORF", "ATG", L = 100L, C = 400L, params = rdg_parameters(p_init = c(ATG = 0.8), reinit_max = 0))
  d_old <- rdg_gene_loads(500L, 560L, "dORF", "ATG", L = 100L, C = 400L, params = old_params(p_init = c(ATG = 0.8), reinit_max = 0))
  expect_equal(d$initiations, d_old$initiations)    # dORF keeps the codon table under the new defaults
  p2 <- rdg_parameters(p_init = c(ATG = 0.8), reinit_max = 0, translon_classes = "dORF")
  u <- rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L, params = p2)
  expect_equal(u$initiations, 0.8)                  # uORF not in translon_classes keeps the codon table
})

test_that("simRDGflow: old and new defaults draw no extra random numbers; translon_gene_sd draws one factor per gene", {
  set.seed(4)
  n <- 6
  dir <- tempfile("rdg_opt_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = n, out_dir = dir, architecture = "all",
    leader_length = rep(200, n), cds_length = rep(600, n), trailer_length = rep(300, n), debug_on = FALSE))
  counts <- suppressMessages(simCountTables(ORFik::loadRegion(genome["txdb"], "cds"), libtypes = "RFP",
    conditions = c("WT", "ISR"), replicates = 1, plot_PCA = FALSE, print_statistics = FALSE))
  SummarizedExperiment::assay(counts)[] <- 4000
  regions <- suppressMessages(simCountTablesRegions(counts))
  cds <- function(x) SummarizedExperiment::assay(x, "cds")
  set.seed(1); a <- simRDGflow(regions, genome); r_a <- runif(1)
  set.seed(1); o <- simRDGflow(regions, genome, old_params()); r_o <- runif(1)
  expect_identical(r_a, r_o)                        # same random stream after the flow
  expect_null(S4Vectors::metadata(a)$rdg_gene_factor)
  expect_true(sum(cds(a)) > sum(cds(o)))            # the new defaults give the CDS more reads
  set.seed(1); g <- simRDGflow(regions, genome, rdg_parameters(translon_gene_sd = 0.8))
  f <- S4Vectors::metadata(g)$rdg_gene_factor
  expect_length(f, n)
  expect_true(all(f > 0) && length(unique(f)) == n)
  translated <- function(x) cds(x) + SummarizedExperiment::assay(x, "uorf")
  expect_equal(translated(g), translated(regions))  # reads are only moved between CDS and translons
  expect_equal(translated(a), translated(regions))
})

test_that("translon_exempt: auto only in ISR-type scenarios; ids, functions and NULL", {
  cl <- list(g1 = "uORF", g2 = "dORF", g3 = c("uORF", "uoORF")); en <- list(30L, 900L, c(30L, 450L)); L <- list(100L, 100L, 399L)
  ex <- function(p) rdg_translon_exempt(p, names(cl), cl, en, L)
  expect_false(any(ex(rdg_parameters())))                                            # no condition changes reinit_distance
  expect_false(any(ex(rdg_parameters(conditions = list(ISR = list(p_init_scale = 0.8))))))
  expect_false(any(ex(rdg_parameters(reinit_distance = 50, conditions = list(B = list(reinit_distance = 50))))))
  isr <- rdg_parameters(conditions = list(ISR = list(reinit_distance = 150)))
  expect_equal(ex(isr), c(g1 = FALSE, g2 = FALSE, g3 = TRUE))                         # only uORF + uoORF (ATF4 layout)
  expect_false(any(ex(rdg_parameters(p_init_translon = NULL, conditions = list(ISR = list(reinit_distance = 150))))))
  expect_false(any(ex(rdg_parameters(translon_exempt = NULL, conditions = list(ISR = list(reinit_distance = 150))))))
  expect_equal(ex(rdg_parameters(translon_exempt = "g2")), c(g1 = FALSE, g2 = TRUE, g3 = FALSE))
  f <- function(gene, classes, ends, L) "uoORF" %in% classes
  expect_equal(ex(rdg_parameters(translon_exempt = f)), c(g1 = FALSE, g2 = FALSE, g3 = TRUE))
  expect_error(rdg_parameters(translon_exempt = 3), "translon_exempt")
  # under ISR a gene with a single uORF and no uoORF is not exempt and keeps 0.3
  p_isr <- rdg_condition_parameters(isr, "ISR")
  single <- rdg_gene_loads(11L, 40L, "uORF", "ATG", L = 100L, C = 400L, params = p_isr, exempt = ex(isr)[["g1"]])
  expect_equal(single$initiations, 0.3)
  # an exempt gene keeps the codon table and ignores the per-gene factor
  p <- rdg_parameters(p_init = c(ATG = 0.8, CTG = 0.1), reinit_max = 0)
  x <- rdg_gene_loads(11L, 40L, "uORF", "CTG", L = 100L, C = 400L, params = p, gene_factor = 2, exempt = TRUE)
  expect_equal(x$initiations, 0.1)
})

test_that("simRDGflow: an ISR-type scenario exempts reinitiation-type genes and records them", {
  set.seed(4)
  n <- 4
  dir <- tempfile("rdg_ex_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = n, out_dir = dir, architecture = "all",
    leader_length = rep(200, n), cds_length = rep(600, n), trailer_length = rep(300, n), debug_on = FALSE))
  counts <- suppressMessages(simCountTables(ORFik::loadRegion(genome["txdb"], "cds"), libtypes = "RFP",
    conditions = c("WT", "ISR"), replicates = 1, plot_PCA = FALSE, print_statistics = FALSE))
  SummarizedExperiment::assay(counts)[] <- 4000
  regions <- suppressMessages(simCountTablesRegions(counts))
  cds <- function(x) SummarizedExperiment::assay(x, "cds")
  isr <- list(ISR = list(reinit_distance = 200))
  set.seed(1); a <- simRDGflow(regions, genome, rdg_parameters(conditions = isr))
  ids <- S4Vectors::metadata(a)$rdg_translon_exempt
  expect_setequal(ids, rownames(regions))                     # architecture "all": every gene has uORF + uoORF
  set.seed(1); b <- simRDGflow(regions, genome, rdg_parameters(conditions = isr, translon_exempt = ids))
  expect_equal(cds(a), cds(b))
  set.seed(1); o <- simRDGflow(regions, genome, rdg_parameters(conditions = isr, translon_exempt = NULL))
  expect_null(S4Vectors::metadata(o)$rdg_translon_exempt)
  expect_false(isTRUE(all.equal(cds(a), cds(o))))
  set.seed(1); w <- simRDGflow(regions, genome)                # no ISR-type scenario: nothing exempt
  expect_null(S4Vectors::metadata(w)$rdg_translon_exempt)
})
