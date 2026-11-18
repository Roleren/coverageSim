# Regression tests from Codex's independent review of steps 1-2 (2026-11-18): each failed on 546b9b7.
local_dir <- function() { d <- tempfile('rev_'); dir.create(d); d }
test_that('short architecture CDS is independent of legacy uORF maximum', {
  expect_silent(suppressMessages(simGenome(n=1, out_dir=local_dir(), architecture='uORF',
    leader_length=30, cds_length=9, trailer_length=10, export_txdb=FALSE, debug_on=FALSE)))
})
test_that('fractional feature lengths are rejected instead of truncated', {
  expect_error(gene_architecture('uORF', 12.9), 'length')
})
test_that('feature starts cannot be a configured stop codon', {
  expect_error(suppressMessages(simGenome(n=1, out_dir=local_dir(),
    architecture=gene_architecture('uORF', 12, start_codon='TAA'),
    leader_length=90, cds_length=120, trailer_length=90, export_txdb=FALSE,
    debug_on=FALSE)), 'start|stop')
})
test_that('leader-to-trailer incidental ORFs are not classified as uoORFs', {
  expect_false(coverageSim:::classify_orfs(10L, 151L, 30L, 120L) == 'uoORF')
})
test_that('rebuilding without uORFs cannot reuse a stale uORF file', {
  d <- local_dir(); set.seed(2)
  suppressMessages(simGenome(n=1, out_dir=d, architecture='uORF',
    leader_length=90, cds_length=120, trailer_length=90, export_txdb=FALSE, debug_on=FALSE))
  g <- suppressMessages(simGenome(n=1, out_dir=d, architecture='dORF',
    leader_length=90, cds_length=120, trailer_length=90, export_txdb=FALSE, debug_on=FALSE))
  expect_false(file.exists(file.path(dirname(g[['genome']]), 'true_uORFs.rds')))
})
test_that('caller supplied architecture tables are validated', {
  a <- gene_architecture('uORF', 12); a$length <- 13L
  expect_error(suppressMessages(simGenome(n=1,out_dir=local_dir(),architecture=a,
    leader_length=90,cds_length=120,trailer_length=90,export_txdb=FALSE,
    debug_on=FALSE)), 'length|multiple')
})
test_that('build log reports architecture uORFs actually requested and built', {
  d <- local_dir(); set.seed(1)
  g <- suppressMessages(simGenome(n=1,out_dir=d,architecture='uORF',leader_length=90,
    cds_length=120,trailer_length=90,export_txdb=FALSE,debug_on=FALSE))
  b <- fread(g[['build_log']])
  expect_equal(b$requested_uorfs, 1L)
  expect_equal(b$built_uorfs, 1L)
})
test_that('NTE extending beyond the leader fails for the right gene before sequence writing', {
  a <- coverageSim:::architecture_per_gene(list(NULL, gene_architecture('NTE',12)),
                                           2L, c('g1','g2'))
  expect_error(coverageSim:::place_features(a, c(30L,6L), c(126L,126L), c(30L,30L)), 'g2')
})
