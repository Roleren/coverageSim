#' 3'UTR background and ribosome-specific 3'UTR plateaus (opt-in, explainability port S10)
#'
#' Set with \code{options(coverageSim.utr3 = list(...))}; unset fields keep the defaults, which switch every
#' feature off (the original flat 3'UTR). \code{utr3_preset("human")} gives the measured human values.
#' Three features act on the 3'UTR (\code{trailer}) region of \code{simNGScoverage()}; their positions are drawn once
#' per genome (seed \code{seed + 5000}, session random state untouched), so they are shared by all libraries, as in
#' real data. The planted sites are written to \code{utr3_features.tsv} in \code{out_dir} for truth checks.
#'
#' \describe{
#'  \item{Ribosome-specific plateaus}{Each 3'UTR window of \code{plateau_width} nt (counted from the stop codon) except
#'   the first and the last gets a plateau with probability \code{plateau_prob}, \code{plateau_proximal} of them in the
#'   proximal half. A plateau multiplies the expected coverage of its window by a fold \code{1 + exp(N(log(plateau_fold
#'   - 1), plateau_fold_sdlog))} (median \code{plateau_fold}), in library types \code{ribo_libtypes} only.
#'   Human (backclaude E10-C, E11-B, E12-A/B): 3.8-5.2\% of 3'UTR windows keep an excess after the UTR tau2 band
#'   (\code{E10/C/csv/utr_tau2_pooled_Emin.csv}, \code{utr3_remaining_*.csv}), median O/E 3, 46-56\% of excess reads
#'   at O/E-1 in (1, 2], 6-11\% above 5; 99\% in middle windows, 83-90\% in the proximal half; unframed; broad
#'   (18-40 distinct 5' ends per window); present in 36 of 40 studies; RNA of the same samples does not carry them
#'   (fold 1.02 vs ribo 1.61, E12-B).}
#'  \item{Alternative polyadenylation (APA) step}{With probability \code{apa_prob} a 3'UTR gets one internal poly(A)
#'   site, placed uniformly between 10\% and 90\% of its length (3'UTRs of at least 100 nt); coverage downstream of the
#'   site is multiplied by \code{apa_ratio} in all library types (\code{apa_libtypes = NULL}), since the shorter isoform
#'   lacks that part for every assay. Human: 54\% of 1,231 MANE 3'UTRs carry at least one internal GENCODE v44 poly(A)
#'   site (backclaude \code{E11/B/csv/tx_features.csv}); windows downstream of an internal site have median O/E 0.57
#'   vs 1.11 upstream (claude-megafst A1 screen, explainability_loop/ITERATIONS.md), ratio 0.51.}
#'  \item{Depleted transcript end}{The last \code{end_width} nt of the 3'UTR get \code{end_depletion} times the
#'   coverage in \code{ribo_libtypes}. Human (E12-A): last 25 nt 0.04-0.09x the 3'UTR mean (controls 0.07-0.13x).}
#' }
#' The window-level shared 3'UTR variance (tau2 0.25-0.45) is a separate option:
#' \code{codon_bias_options(tau2 = c(utr3 = ...))} (S9).
#' @param plateau_prob numeric in [0, 1], probability per eligible 3'UTR window (0 = off).
#' @param plateau_width integer, plateau and window length in nt (90 = the 30-codon window of the density test).
#' @param plateau_fold numeric > 1, median fold (O/E) of a plateau window.
#' @param plateau_fold_sdlog numeric >= 0, sd of log(fold - 1) between plateaus.
#' @param plateau_proximal numeric in [0, 1], share of plateaus in the proximal half of the 3'UTR.
#' @param ribo_libtypes character, library types that carry plateaus and the depleted end (default "RFP").
#' @param apa_prob numeric in [0, 1], probability that a 3'UTR has an internal poly(A) site (0 = off).
#' @param apa_ratio numeric > 0, coverage ratio downstream / upstream of the internal site.
#' @param apa_libtypes character or NULL (all library types) that see the APA step.
#' @param end_depletion numeric > 0, factor on the last \code{end_width} nt (1 = off).
#' @param end_width integer, length in nt of the depleted 3' end.
#' @param seed integer seed of the per-genome draws.
#' @return a list of options
#' @examples
#' options(coverageSim.utr3 = utr3_preset("human"))
#' utr3_options()$plateau_prob
#' options(coverageSim.utr3 = NULL)
#' @export
utr3_options <- function(plateau_prob = 0, plateau_width = 90L, plateau_fold = 3, plateau_fold_sdlog = 0.65,
                         plateau_proximal = 0.85, ribo_libtypes = "RFP", apa_prob = 0, apa_ratio = 0.51,
                         apa_libtypes = NULL, end_depletion = 1, end_width = 25L, seed = 1L) {
  def <- list(plateau_prob = plateau_prob, plateau_width = plateau_width, plateau_fold = plateau_fold,
              plateau_fold_sdlog = plateau_fold_sdlog, plateau_proximal = plateau_proximal,
              ribo_libtypes = ribo_libtypes, apa_prob = apa_prob, apa_ratio = apa_ratio,
              apa_libtypes = apa_libtypes, end_depletion = end_depletion, end_width = end_width, seed = seed)
  set <- getOption("coverageSim.utr3", list())
  bad <- setdiff(names(set), names(def)); if (length(bad)) stop("unknown coverageSim.utr3 fields: ", toString(bad))
  def[names(set)] <- set
  stopifnot(def$plateau_prob >= 0, def$plateau_prob <= 1, def$plateau_width >= 3, def$plateau_fold > 1,
            def$plateau_fold_sdlog >= 0, def$plateau_proximal >= 0, def$plateau_proximal <= 1,
            def$apa_prob >= 0, def$apa_prob <= 1, def$apa_ratio > 0, def$end_depletion > 0, def$end_width >= 1,
            is.character(def$ribo_libtypes), is.null(def$apa_libtypes) || is.character(def$apa_libtypes))
  def
}

#' Named 3'UTR settings
#' @param name "original" (all off) or "human" (plateaus 4.5\% of windows, median O/E 3; APA in 54\% of 3'UTRs with
#'  ratio 0.51; transcript end 0.07x; see \code{\link{utr3_options}} for the sources).
#' @return a list for \code{options(coverageSim.utr3 = ...)}
#' @export
utr3_preset <- function(name = c("original", "human")) {
  switch(match.arg(name),
         original = list(),
         human = list(plateau_prob = 0.045, plateau_fold = 3, plateau_fold_sdlog = 0.65, plateau_proximal = 0.85,
                      apa_prob = 0.54, apa_ratio = 0.51, end_depletion = 0.07, end_width = 25L))
}

utr3_active <- function(opts = utr3_options()) {
  opts$plateau_prob > 0 || opts$apa_prob > 0 || opts$end_depletion != 1
}

#' Draw the 3'UTR features of a genome
#'
#' The sites \code{simNGScoverage()} plants when \code{\link{utr3_options}} are active, drawn with the fixed seed, so
#' the same 3'UTR lengths (in the same order) always give the same sites.
#' @param lengths integer, 3'UTR lengths in nt, one per transcript.
#' @param transcript_id character, transcript names (default \code{names(lengths)}).
#' @param opts options, default \code{utr3_options()}.
#' @return a data.table with one row per feature: transcript_id, gene (index), type ("plateau", "apa", "end"),
#'  start, end (1-based 3'UTR coordinates) and factor (coverage multiplier).
#' @export
draw_utr3_features <- function(lengths, transcript_id = names(lengths), opts = utr3_options()) {
  lengths <- as.integer(lengths)
  if (is.null(transcript_id)) transcript_id <- as.character(seq_along(lengths))
  empty <- data.table::data.table(transcript_id = character(0), gene = integer(0), type = character(0),
                                  start = integer(0), end = integer(0), factor = numeric(0))
  cbs <- list(per_genome_draws = TRUE, seed = opts$seed + 5000L)
  out <- with_gene_seed(cbs, length(lengths), function(n) {
    w <- opts$plateau_width
    plateaus <- if (opts$plateau_prob > 0) data.table::rbindlist(lapply(seq_len(n), function(i) {
      nw <- lengths[i] %/% w
      if (nw < 3L) return(NULL)
      win <- seq.int(2L, nw - 1L)
      proximal <- (win - 0.5) * w <= lengths[i] / 2
      p <- pmin(1, opts$plateau_prob * ifelse(proximal, 2 * opts$plateau_proximal, 2 * (1 - opts$plateau_proximal)))
      hit <- win[stats::runif(length(win)) < p]
      if (!length(hit)) return(NULL)
      fold <- 1 + exp(stats::rnorm(length(hit), log(opts$plateau_fold - 1), opts$plateau_fold_sdlog))
      data.table::data.table(gene = i, type = "plateau", start = (hit - 1L) * w + 1L, end = hit * w, factor = fold)
    }))
    apa <- if (opts$apa_prob > 0) data.table::rbindlist(lapply(seq_len(n), function(i) {
      hit <- stats::runif(1) < opts$apa_prob
      site <- stats::runif(1, 0.1, 0.9)
      if (!hit || lengths[i] < 100L) return(NULL)
      s <- as.integer(ceiling(site * lengths[i]))
      data.table::data.table(gene = i, type = "apa", start = s + 1L, end = lengths[i], factor = opts$apa_ratio)
    }))
    list(plateaus, apa)
  })
  ends <- if (opts$end_depletion != 1) {
    keep <- which(lengths > opts$end_width)
    data.table::data.table(gene = keep, type = "end", start = lengths[keep] - opts$end_width + 1L,
                           end = lengths[keep], factor = opts$end_depletion)
  }
  res <- data.table::rbindlist(c(out, list(ends)), use.names = TRUE)
  if (!nrow(res)) return(empty)
  res[, transcript_id := as.character(transcript_id[gene])]
  data.table::setcolorder(res, "transcript_id")
  res[order(gene, start)]
}

# Per-nucleotide multiplier of one 3'UTR of length n for one library type.
utr3_factor_vector <- function(n, features, libtype, opts = utr3_options()) {
  f <- rep(1, n)
  if (is.null(features) || !nrow(features)) return(f)
  ribo <- libtype %in% opts$ribo_libtypes
  apa <- is.null(opts$apa_libtypes) || libtype %in% opts$apa_libtypes
  for (j in seq_len(nrow(features))) {
    type <- features$type[j]
    if ((type %in% c("plateau", "end") && !ribo) || (type == "apa" && !apa)) next
    idx <- seq.int(features$start[j], min(features$end[j], n))
    f[idx] <- f[idx] * features$factor[j]
  }
  f
}

# The 3'UTR multipliers of a run, or NULL when off: drawn once per simNGScoverage() call.
utr3_run_features <- function(lengths, transcript_id) {
  opts <- utr3_options()
  if (!utr3_active(opts)) return(NULL)
  feats <- draw_utr3_features(lengths, transcript_id, opts)
  list(opts = opts, features = feats, by_gene = split(feats, factor(feats$gene, levels = seq_along(lengths))))
}

# Multiply the weights of 3'UTR gene i (vector with an optional symmetric margin around the 3'UTR) in place.
apply_utr3_factors <- function(weights, n, i, libtype, utr3) {
  margin <- max((length(weights) - n) %/% 2L, 0L)
  f <- utr3_factor_vector(n, utr3$by_gene[[i]], libtype, utr3$opts)
  idx <- margin + seq_len(n)
  weights[idx] <- weights[idx] * f
  weights
}
