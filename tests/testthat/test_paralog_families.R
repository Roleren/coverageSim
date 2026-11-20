test_that("paralog families: all-mapper counts every multimapper in full at each member, unique track keeps member-specific reads", {
  x <- simParalogFamilies(30, k = 3, divergence = 0.25, seed = 1, depth = 3000)
  tr <- data.table::copy(x$tracks)
  expect_equal(mean(x$members$q_pos), 0.25, tolerance = 0.01)
  tr[, `:=`(nG = sum(!specific), m = sum(Y * !specific)), by = .(fam, condition, replicate, codon)]
  # member-specific codon: all-mapper = unique = own reads; shared codon: all-mapper = reads of all members that codon maps to
  expect_true(all(tr[specific == TRUE, A == Y & U == Y]))
  expect_true(all(tr[specific == FALSE, A == m]))
  # unique reads only where the codon is specific or maps to one member only
  expect_true(all(tr[specific == FALSE & nG >= 2L, U == 0L]))
  expect_true(all(tr[specific == FALSE & nG == 1L, U == Y]))
  expect_true(all(tr$A >= tr$U))
  # identical members (divergence 0): no unique reads, every member shows the family total
  z <- simParalogFamilies(10, k = 2, divergence = 0, seed = 2)
  expect_equal(sum(z$tracks$U), 0)
  tt <- paralog_totals(z)
  expect_equal(tt$A, tt$M)
})

test_that("shares, planted abundance change, flank mode, leakage and RNG handling", {
  set.seed(9); r1 <- stats::runif(1); set.seed(9)
  x <- simParalogFamilies(200, k = 2, shares = c(0.85, 0.15), planted = list(member = 1, factor = 1 / 3), seed = 3, depth = 4000)
  expect_equal(stats::runif(1), r1)                       # session RNG restored
  tt <- paralog_totals(x)
  s <- tt[condition == "A", .(sh = sum(Y[member == 1]) / sum(Y)), by = replicate]
  expect_equal(s$sh, rep(0.85, 2), tolerance = 0.02)
  fc <- tt[member == 1, .(r = sum(Y[condition == "B"]) / sum(Y[condition == "A"]))]$r
  expect_equal(fc, 1 / 3, tolerance = 0.03)
  fc2 <- tt[member == 2, .(r = sum(Y[condition == "B"]) / sum(Y[condition == "A"]))]$r
  expect_equal(fc2, 1, tolerance = 0.03)
  expect_true(all(x$members[member == 1]$planted)); expect_false(any(x$members[member == 2]$planted))
  expect_identical(simParalogFamilies(5, seed = 4)$tracks, simParalogFamilies(5, seed = 4)$tracks)
  # flank mode: one contiguous block of specific codons per member
  fl <- simParalogFamilies(10, k = 2, divergence = 0.2, specific_mode = "flank", seed = 5)$tracks
  runs <- fl[condition == "A" & replicate == 1, sum(diff(c(FALSE, specific)) == 1L), by = .(fam, member)]$V1
  expect_true(all(runs == 1L))
  # leakage: ~20 % of the unique reads at specific codons come from other members
  lk <- simParalogFamilies(150, k = 2, divergence = 0.2, leak = 0.2, seed = 6, depth = 6000)$tracks
  lk <- lk[specific == TRUE]
  expect_equal(sum(lk$U - lk$Y) / sum(lk$U), 0.2, tolerance = 0.1)
  expect_error(simParalogFamilies(2, k = 1))
})

test_that("EM on unique + multimapper counts recovers member totals where the all-mapper track overcounts", {
  x <- simParalogFamilies(150, k = 2, divergence = 0.2, shares = c(0.7, 0.3), seed = 7, depth = 6000)
  tt <- paralog_totals(x)[condition == "A" & replicate == 1]
  em <- function(U, M, q, it = 200) { T <- U / q + 1e-9; for (i in seq_len(it)) { w <- T * (1 - q); T <- U + M * w / sum(w) }; T }
  d <- tt[, { T <- em(U, M[1], q_unique); .(est = T, Y = Y, A = A) }, by = fam]
  expect_gt(median(d$A / d$Y), 1.5)
  m <- x$members; expect_true(all(m$q_unique >= m$q_pos))                         # all-mapper track: multimappers counted at every member
  expect_lt(abs(median(d$est / d$Y) - 1), 0.1)              # EM: family total conserved, shares recovered
})
