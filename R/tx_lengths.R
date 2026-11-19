# Realistic transcript region lengths (t73, 2026-11-19).
#
# Log-normal fits to the 19,099 human Ensembl 108 canonical protein-coding
# transcripts with a 5'UTR, a 3'UTR and a stop codon (GTF
# Homo_sapiens.GRCh38.108_ensembl.gtf, tag Ensembl_canonical; exon lengths of
# five_prime_utr / CDS + stop_codon / three_prime_utr features). The three
# regions are correlated on the log scale (r = 0.21 to 0.23), reproduced with
# one shared latent normal (rho = 0.22). Draws are cut at the real 0.5 % and
# 99.5 % quantiles (CDS lower bound 192 nt incl. start and stop, above the smallest cds_length simGenome accepts).
human_length_model <- list(
  meanlog = c(leader = 4.8535, cds = 7.1710, trailer = 6.8250),  # CDS incl. start and stop codon
  sdlog   = c(leader = 1.0107, cds = 0.7568, trailer = 1.2816),
  rho     = 0.22,
  lower   = c(leader = 8, cds = 192, trailer = 20),
  upper   = c(leader = 1266, cds = 11126, trailer = 11844)
)

#' Draw realistic transcript region lengths
#'
#' Per gene log-normal draws for the leader, the coding sequence and the trailer,
#' fitted to human Ensembl 108 canonical protein-coding transcripts
#' (medians 136 / 1284 / 1018 nt, see \code{human_length_model}).
#' @param n number of genes.
#' @param min_leader smallest leader to return, default 20 (a leader shorter
#'   than the longest uORF you ask for is infeasible; raise it then).
#' @return a list with \code{leader_length}, \code{cds_length} (without start and
#'   stop codon, divisible by 3, as \code{simGenome} wants it) and \code{trailer_length}.
#' @export
human_tx_lengths <- function(n, min_leader = 20L) {
  m <- human_length_model
  shared <- stats::rnorm(n)
  draw <- function(i) {
    z <- sqrt(m$rho) * shared + sqrt(1 - m$rho) * stats::rnorm(n)
    pmin(pmax(exp(m$meanlog[i] + m$sdlog[i] * z), m$lower[i]), m$upper[i])
  }
  leader <- pmax(as.integer(round(draw(1))), as.integer(min_leader))
  cds_full <- 3L * as.integer(round(draw(2) / 3))
  trailer <- as.integer(round(draw(3)))
  list(leader_length = leader, cds_length = cds_full - 6L, trailer_length = trailer)
}
