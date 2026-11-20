# Paralog families / multimappers (backclaude E17-S8, 2026-11-20).
#
# Near-identical genes (RPL/RPS pseudogene families, actins, histones) share most of their CDS: a read from such a position
# aligns equally well to every member (a multimapper). Aligners then count it at each member in an all-mapper track
# and drop it from a unique-mapper track; only reads over member-specific positions (the unique flank, a fraction q of the
# member's reads) are unique (E5/B, E6/C, E7/B, E8/B, E9/C). This builds such families constructively on a codon grid
# (all members have the same CDS length and the same dwell profile) and returns the true member reads plus both tracks.

#' Simulate paralog families with all-mapper and unique-mapper tracks
#'
#' Per family, \code{k} members with identical dwell profile. Member j has a set S_j of member-specific codon positions
#' (fraction \code{divergence}, scattered or one contiguous flank). A read originating at a codon that is specific to its
#' member is unique to that member; a read at a codon shared with other members is a multimapper among all members for which
#' that codon is not specific (unique if only its own member qualifies). All-mapper track: every multimapper read is counted
#' in full at each member it maps to. Unique track: only unique reads. Optional \code{leak}: the fraction of the unique
#' reads at a member's specific positions that really come from other members (mismatch-tolerant alignment, E9/C: 19.8\% in
#' human pseudogene regions); those reads are added to the member's unique and all-mapper track.
#' @param n_fam number of families.
#' @param k family size (2 or more; one number or one per family).
#' @param shares member expression shares: \code{"equal"} (default), \code{"unequal"} (Dirichlet(1), seeded) or a numeric
#'  vector of length k (recycled over families; normalised).
#' @param divergence fraction q of member-specific codon positions per member: one number or a vector of length k
#'  (default 0.2; 0 = identical members, no unique reads).
#' @param specific_mode \code{"scattered"} (default: random codons per member) or \code{"flank"} (one contiguous block).
#' @param cds_codons CDS length in codons, default 300.
#' @param leak fraction of unique reads at specific positions that come from other members (default 0).
#' @param dwell_sd sd of the log-normal per-codon dwell shared by all members of a family (default 0.5).
#' @param conditions,replicates libraries: names of the conditions and number of replicates per condition.
#' @param depth median family total reads per library; \code{depth_sd} sd of its log-normal spread between families.
#' @param planted \code{NULL} (default) or a list \code{member}, \code{factor}, \code{condition} (default the last one) and
#'  \code{fraction} (share of families, default 1) that multiplies the expression of one member in one condition (planted
#'  abundance change, e.g. factor 1/3).
#' @param seed integer or NULL (session RNG restored afterwards).
#' @return an object of class \code{paralogFamilies}: \code{members} (fam, member, share, planted, q_pos = fraction of
#'  member-specific codons, q_unique = fraction of codons where the member's reads are unique (specific, or the only member
#'  that fits: the mappability an analyst would compute from k-mers), q_reads = expected fraction of the member's reads that
#'  are unique), and \code{tracks}
#'  (fam, member, condition, replicate, codon, specific, Y = true member reads, A = all-mapper track, U = unique track).
#' @export
#' @examples
#' f <- simParalogFamilies(5, k = 3, seed = 1)
#' f$members
simParalogFamilies <- function(n_fam, k = 2L, shares = "equal", divergence = 0.2, specific_mode = c("scattered", "flank"),
                               cds_codons = 300L, leak = 0, dwell_sd = 0.5, conditions = c("A", "B"), replicates = 2L,
                               depth = 5000, depth_sd = 0.5, planted = NULL, seed = NULL) {
  specific_mode <- match.arg(specific_mode)
  stopifnot(n_fam >= 1, all(k >= 2), all(divergence >= 0), all(divergence < 1), leak >= 0, leak < 1, cds_codons >= 10)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, globalenv()))
    set.seed(seed)
  }
  k <- rep_len(as.integer(k), n_fam); L <- as.integer(cds_codons)
  if (!is.null(planted)) {
    planted <- modifyList(list(member = 1L, factor = 1 / 3, condition = utils::tail(conditions, 1L), fraction = 1), planted)
    stopifnot(planted$condition %in% conditions, planted$factor > 0)
  }
  libs <- data.table::CJ(condition = factor(conditions, levels = conditions), replicate = seq_len(replicates), sorted = FALSE)
  mem <- list(); trk <- list()
  for (f in seq_len(n_fam)) {
    kk <- k[f]; fid <- sprintf("F%04d", f)
    p <- if (is.character(shares)) { if (shares[1] == "unequal") { d <- stats::rgamma(kk, 1); d / sum(d) } else rep(1 / kk, kk) } else {
      v <- rep_len(as.numeric(shares), kk); v / sum(v) }
    q <- rep_len(divergence, kk)
    dw <- if (dwell_sd > 0) exp(stats::rnorm(L, 0, dwell_sd)) else rep(1, L); dw <- dw / sum(dw)
    S <- matrix(FALSE, L, kk)
    for (j in seq_len(kk)) {
      nsp <- as.integer(round(q[j] * L))
      if (nsp > 0L) {
        idx <- if (specific_mode == "scattered") sample.int(L, nsp) else { s0 <- sample.int(L - nsp + 1L, 1L); s0 + seq_len(nsp) - 1L }
        S[idx, j] <- TRUE
      }
    }
    fac <- rep(1, kk); pl <- rep(FALSE, kk)
    do_plant <- !is.null(planted) && planted$member <= kk && stats::runif(1) < planted$fraction
    gf <- exp(stats::rnorm(1, 0, depth_sd))
    nG <- rowSums(!S)
    Uq <- S | nG == 1L        # codons where a member's reads are unique (specific to it, or the only member that fits)
    mem[[f]] <- data.table::data.table(fam = fid, member = seq_len(kk), share = p, planted = FALSE,
                                       q_pos = colMeans(S), q_unique = colMeans(Uq), q_reads = as.vector(crossprod(dw, Uq)))
    if (do_plant) mem[[f]]$planted[planted$member] <- TRUE
    for (li in seq_len(nrow(libs))) {
      fac <- rep(1, kk)
      if (do_plant && as.character(libs$condition[li]) == planted$condition) fac[planted$member] <- planted$factor
      mu <- depth * gf * outer(dw, p * fac)
      Y <- matrix(stats::rpois(L * kk, mu), L, kk)
      X <- Y * (!S); m <- rowSums(X)                              # reads at codons shared with other members
      A <- ifelse(S, Y, m); U <- ifelse(S | nG == 1L, Y, 0L)     # multimappers counted in full at every member they map to
      if (leak > 0) for (j in seq_len(kk)) {
        sj <- S[, j]; den <- sum(m[sj])
        if (any(sj) && den > 0) {
          lam <- min(1, leak / (1 - leak) * sum(Y[sj, j]) / den)
          lk <- stats::rbinom(sum(sj), m[sj], lam)
          A[sj, j] <- A[sj, j] + lk; U[sj, j] <- U[sj, j] + lk
        }
      }
      trk[[length(trk) + 1L]] <- data.table::data.table(
        fam = fid, member = rep(seq_len(kk), each = L), condition = libs$condition[li], replicate = libs$replicate[li],
        codon = rep(seq_len(L), kk), specific = as.vector(S), Y = as.vector(Y), A = as.vector(A), U = as.vector(U))
    }
  }
  structure(list(members = data.table::rbindlist(mem), tracks = data.table::rbindlist(trk),
                 options = list(k = k, divergence = divergence, specific_mode = specific_mode, leak = leak, planted = planted)),
            class = "paralogFamilies")
}

#' @export
print.paralogFamilies <- function(x, ...) {
  m <- x$members
  cat("paralogFamilies:", data.table::uniqueN(m$fam), "families,", nrow(m), "members; mean unique-read fraction q",
      round(mean(m$q_reads), 3), "\n")
  invisible(x)
}

#' Member-level totals of a paralog-family simulation
#'
#' Sums the tracks per library and member: true reads (\code{Y}), all-mapper reads (\code{A}), unique reads (\code{U}), and
#' the multimapper reads of the family counted once (\code{M}: reads at codons shared by two or more members).
#' @param x a \code{paralogFamilies} object.
#' @return a data.table, one row per family, member, condition and replicate.
#' @export
paralog_totals <- function(x) {
  tr <- x$tracks
  tr[, nG := sum(!specific), by = .(fam, condition, replicate, codon)]
  tr[, Ysh := ifelse(!specific, Y, 0L)]
  mm <- tr[nG >= 2L & !specific, .(M = sum(Y)), by = .(fam, condition, replicate)]   # each multimapper read once
  out <- tr[, .(Y = sum(Y), A = sum(A), U = sum(U)), by = .(fam, member, condition, replicate)]
  out <- merge(out, mm, by = c("fam", "condition", "replicate"), all.x = TRUE)
  out[is.na(M), M := 0L]
  tr[, c("nG", "Ysh") := NULL]
  merge(out, x$members, by = c("fam", "member"))[order(fam, member, condition, replicate)]
}
