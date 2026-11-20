# Exact uORF counts per gene, the feasibility check before building, and the
# per-gene build log (2026-11-18).

build_log_of <- function(genome) data.table::fread(genome[["build_log"]])

test_that("an exact uORF request is met and logged per gene", {
  set.seed(1)
  dir <- tempfile("exact_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = 6, out_dir = dir, n_uorfs = 2,
    leader_length = rep(140, 6), uorfs_can_overlap = FALSE, debug_on = FALSE))
  log <- build_log_of(genome)
  expect_equal(nrow(log), 6)
  expect_true(all(log$requested_uorfs == 2 & log$built_uorfs == 2 & log$status == "ok"))
  uorfs <- readRDS(genome[["uorfs"]])
  expect_equal(as.vector(table(names(uorfs))), rep(2L, 6))
})

test_that("different counts per gene, including 0, are built exactly", {
  set.seed(3)
  dir <- tempfile("per_gene_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = 3, out_dir = dir, n_uorfs = c(0, 1, 2),
    leader_length = rep(150, 3), uorfs_can_overlap = FALSE, debug_on = FALSE))
  log <- build_log_of(genome)
  expect_equal(log$built_uorfs, c(0L, 1L, 2L))
  expect_true(all(log$status == "ok"))
})

test_that("an impossible request stops before any sequence is made", {
  dir <- tempfile("infeasible_"); dir.create(dir)
  expect_error(simGenome(n = 2, out_dir = dir, n_uorfs = 2, leader_length = c(140, 8),
                         uorfs_can_overlap = FALSE, uorfs_can_overlap_cds = 0,
                         debug_on = FALSE),
               "ENSTTEST10002: 2 uORF\\(s\\) need a leader of at least 12 nt, got 8")
  expect_false(file.exists(file.path(dir, "Homo_sapiens_dummy.fasta")))
  expect_error(simGenome(n = 1, out_dir = dir, n_uorfs = 2, leader_length = 140,
                         uorfs_can_overlap = FALSE, uorfs_can_overlap_cds = 2,
                         debug_on = FALSE),
               "must all run into the CDS")
})

test_that("on_infeasible = 'warning' builds what it can and logs the rest", {
  set.seed(1)
  dir <- tempfile("warn_"); dir.create(dir)
  expect_warning(genome <- suppressMessages(simGenome(n = 2, out_dir = dir, n_uorfs = 2,
    leader_length = c(140, 8), uorfs_can_overlap = FALSE, uorfs_can_overlap_cds = 0,
    on_infeasible = "warning", debug_on = FALSE)), "cannot get the requested uORFs")
  log <- build_log_of(genome)
  expect_equal(log$status, c("ok", "infeasible"))
  expect_match(log$reason[2], "need a leader of at least 12 nt")
})

test_that("max_uorfs keeps its 'at most' meaning and is logged", {
  set.seed(2)
  dir <- tempfile("max_"); dir.create(dir)
  genome <- suppressMessages(simGenome(n = 3, out_dir = dir, max_uorfs = 0,
                                       leader_length = rep(140, 3), debug_on = FALSE))
  log <- build_log_of(genome)
  expect_true(all(log$requested_uorfs == 0 & log$built_uorfs == 0 & log$status == "ok"))
  expect_false("uorfs" %in% names(genome))
})
