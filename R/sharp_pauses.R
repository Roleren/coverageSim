# Sharp pauses, Pro-Pro stalls and Pro-run depletion (opt-in, explainability port S13).
#
# Real human Ribo-seq (claude-megafst M16-acf, 2026-11-21, 7.8 M codons of 14,310 MANE genes) has a smooth, positively
# autocorrelated background (lag 1-5 ACF 0.58 / 0.42 / 0.37 / 0.33 / 0.28), but a pause is a SHARP single-codon peak
# whose neighbours are depleted relative to the codon dwell model: around 42,195 Pro-Pro stalls (Pro in the E and P
# sites, non-Pro A site) the residual log density at P-site offsets -1 / 0 / +1 is -0.14 / +0.10 / -0.09. The peak
# gains what the flanks lose. Multiplying a codon after the smoothing (S2) gives the peak but not the depleted flanks;
# these helpers move the extra reads of a pause out of its +-width neighbouring codons, conserving the local total.
# All act on the codon-level weights after the smoothing step of sim_sequence_bias(); see codon_bias_options().

# Does any S13 layer need the codon identities of the genes?
s13_needs_codons <- function(cb = codon_bias_options()) {
  isTRUE(cb$pp_stall) || isTRUE(cb$pro_run) || isTRUE(cb$field$gc_bias)
}

# Sum of v over the flanks -width..-1, 1..width of every position (positions outside the gene count as 0).
flank_sum <- function(v, width) {
  n <- length(v); out <- numeric(n)
  for (d in c(-seq_len(width), seq_len(width))) {
    j <- seq_len(n) + d; ok <- j >= 1L & j <= n
    out[ok] <- out[ok] + v[j[ok]]
  }
  out
}

# Sharp multiplicative pauses that conserve the local total.
# x: codon weights, m: pause factor per codon (1 = no pause). A pause (m > 1) at codon c gains delta = (m - 1) * x[c],
# which its flanks c - width .. c + width (c excluded) give in proportion to their weights. A dip (m < 1) is only
# multiplied: a codon read fast does not push ribosomes onto its neighbours (giving its loss to the flanks would act
# like a second smoothing and raise the autocorrelation). A pause takes at most `take` of its flanks' total, and a codon
# that borders several pauses keeps at least 5% of its weight (where these caps bind, the local total changes a little).
sharp_pause_redistribute <- function(x, m, width = 1L, take = 0.5) {
  xd <- x * pmin(m, 1)
  if (width < 1L || all(m <= 1)) return(xd)
  delta <- pmax(m - 1, 0) * x
  flank <- flank_sum(xd, width)
  delta <- pmin(delta, take * flank)
  r <- ifelse(flank > 0, delta / flank, 0)
  pmax(xd + delta - xd * pmin(flank_sum(r, width), 0.95), 0.05 * xd)
}

# Codon index (in the weight vector) of P-site codon p is p + 1 - a_shift (a_shift = 1: the vector is indexed by the
# P site; a_shift = 0: by the A site). Returns the weight-vector index of each P-site position (NA outside the gene).
psite_to_index <- function(n, a_shift) {
  i <- seq_len(n) + 1L - a_shift
  ifelse(i >= 1L & i <= n, i, NA_integer_)
}

# Pro-Pro stall kernel: P-site codon p with Pro at p - 1 (E site) and p (P site) and a non-Pro A site (p + 1).
# The stall codon gains take[k] of the weight of its neighbour at offset k (take ordered -w..-1, +1..+w), the
# neighbour loses it; then the stall codon is multiplied by `peak` (1 = pure redistribution).
pp_stall_sites <- function(codons) {
  aa <- codon_amino_acid(codons)
  n <- length(aa)
  if (n < 3L) return(integer())
  p <- seq.int(2L, n - 1L)
  p[aa[p - 1L] == "P" & aa[p] == "P" & aa[p + 1L] != "P"]
}

apply_pp_stall <- function(x, codons, a_shift = 1L, take = c(0.20, 0.10), peak = 0.9) {
  sites <- pp_stall_sites(codons)
  if (!length(sites)) return(x)
  n <- length(x); w <- length(take) %/% 2L
  idx <- sites + 1L - a_shift
  idx <- idx[idx >= 1L & idx <= n]
  offsets <- c(-rev(seq_len(w)), seq_len(w))
  y <- x
  for (k in seq_along(offsets)) {
    j <- idx + offsets[k]; ok <- j >= 1L & j <= n
    moved <- take[k] * x[j[ok]]
    y[j[ok]] <- y[j[ok]] - moved
    y[idx[ok]] <- y[idx[ok]] + moved
  }
  y[idx] <- y[idx] * peak
  pmax(y, 0)
}

# Pro-run depletion: the codon whose P site is p gets factor[as.character(nP)], nP = number of Pro codons in p-2..p+2
# (counts without a factor are left at 1).
apply_pro_run <- function(x, codons, a_shift = 1L, factor = c(`3` = 0.766, `4` = 0.45, `5` = 0.45)) {
  pro <- as.numeric(codon_amino_acid(codons) == "P")
  n <- length(pro)
  nP <- pro + flank_sum(pro, 2L)
  f <- unname(factor[as.character(nP)]); f[is.na(f)] <- 1
  idx <- psite_to_index(n, a_shift)
  ok <- !is.na(idx)
  x[idx[ok]] <- x[idx[ok]] * f[ok]
  x
}

codon_amino_acid <- function(codons) {
  aa <- unname(Biostrings::GENETIC_CODE[toupper(as.character(codons))])
  aa[is.na(aa)] <- "X"
  aa
}

# Step 1c3 of sim_sequence_bias(): the codon-identity layers of S13 for one gene.
apply_s13_codon_layers <- function(x, codons, a_shift, cb) {
  if (is.null(codons)) return(x)
  if (length(codons) != length(x)) stop("codons do not match the codon weights of a gene")
  if (isTRUE(cb$pro_run)) x <- apply_pro_run(x, codons, a_shift, cb$pro_run_factor)
  if (isTRUE(cb$pp_stall)) x <- apply_pp_stall(x, codons, a_shift, cb$pp_stall_take, cb$pp_stall_peak)
  x
}
