# Residual field, P-site smear, per-gene contrast and footprint-GC bias (opt-in, explainability port S14).
#
# Real human Ribo-seq (claude-megafst explainability loop, Fable review 2026-11-22, NC-1 / NC-1b / M22) deviates from
# a codon-level dwell model in a smooth, multi-scale, mostly study-shared way: within-gene variance of log O/E by
# block of 1 / 3 / 10 / 30 / 100 codons 0.34 / 0.26 / 0.17 / 0.10 / 0.055 (white noise would give 0.011 at 30 codons),
# window (30-codon) lag-1..4 autocorrelation 0.51 / 0.41 / 0.35 / 0.29, 70-75% of the reproducible part shared between
# study-disjoint halves, a lognormal spread of gene residual sd (0.16-0.47, 10-90%), and a per-gene slope of log O on
# log E around 1.03 (sd ~0.15). Peaks are smeared symmetrically onto their neighbours (offset heterogeneity between
# libraries; 2.15 / 1.36 / 0.96 / 0.70 log at 0 / +-1 / +-2 / +-5 codons around strong peaks), not queued at -10.
# The helpers below act on the codon-level weights of sim_sequence_bias() after the pauses and the codon-identity
# layers (S12, S13); see residual_field_options() and codon_bias_options(field = ...).

#' Residual field options (port S14)
#'
#' Opt-in layers that give simulated coverage the multi-scale, study-shared residual structure of real Ribo-seq
#' relative to a codon dwell model. Pass the result as \code{codon_bias_options(field = ...)}, i.e. inside
#' \code{options(coverageSim.codon_bias = list(..., field = residual_field_options(...)))}; with the defaults every layer
#' is off and the simulation is unchanged. \code{codon_bias_preset("human_v2")} turns on the calibrated set.
#'
#' Order inside \code{sim_sequence_bias()} (codon-level weights of CDS / uORF genes): contrast exponent and smoothing
#' (\code{codon_bias_options(sharpen, smooth)}), pauses, library dwell, Pro-Pro stall / Pro run, then
#' (1) the footprint-GC factor, (2) the log-field \code{exp(f)}, (3) the P-site smear kernel; then the nucleotide
#' expansion and the RNase blur as before.
#'
#' \describe{
#'  \item{log-field}{All field draws are per genome (fixed seeds, whatever \code{per_genome_draws}), so they depend
#'   on gene order and count, like the pauses. f = s_g (sqrt(shared) F_shared + sqrt(1 - shared) F_study), where each F is a sum of independent
#'   stationary AR(1) (discrete Ornstein-Uhlenbeck) processes over codons with correlation exp(-1 / L_k) per codon,
#'   length scales \code{field_scales} (codons) and variances \code{field_var}. F_shared is drawn per genome
#'   (seed + 6000), identical in all libraries; F_study per study (seed + 6500 + a hash of the study name), so two
#'   libraries of different studies share a fraction \code{shared} of the field variance. The study is
#'   \code{colData(count_table)$study} (\code{study_unit = "study"}; without that column all libraries are one study
#'   and the field is fully shared) or the library itself (\code{study_unit = "library"}). A library sees the full
#'   field variance; a merged pool of k equally deep studies sees sum(field_var) (shared + (1 - shared) / k). s_g = exp(gene_sd z - gene_sd^2), z ~ N(0, 1) per gene and genome, so E[s_g^2] = 1: genes differ in
#'   residual amplitude (lognormal, heavy right tail) without changing the mean field variance.}
#'  \item{smear}{A symmetric codon kernel (weights for offsets 0, 1, 2, ...; the same weight at -k and +k) applied after
#'   the pauses and the field, so a pause spreads onto its neighbours as P-site offset heterogeneity between merged
#'   libraries does. Edge codons are renormalised by the kernel weight inside the gene.}
#'  \item{footprint GC}{A technical library bias (M22): factor exp(h(GC)) on the codon whose P site is p, with GC the
#'   GC fraction of the footprint (P-site codons \code{gc_window}, default -4..+5 = 30 nt) and h piecewise linear
#'   through \code{gc_knots} (default: -0.13 at <= 37\% and >= 70\%, +0.06 at 45-60\%, the M22 inverted U), times
#'   \code{gc_strength}. Needs codon identities (codon-level seq_bias).}
#' }
#' @param field_var numeric, variance (log scale) of each AR(1) component; all 0 = no field (default).
#' @param field_scales numeric, length scale in codons of each component (same length as \code{field_var}).
#' @param shared numeric in [0, 1], fraction of the field variance shared by all studies.
#' @param gene_sd numeric >= 0, sd of the log per-gene field amplitude (0 = same amplitude in every gene).
#' @param smear NULL (off) or numeric weights for codon offsets 0, 1, ..., K (symmetric kernel).
#' @param gc_bias logical, apply the footprint-GC factor.
#' @param gc_knots data.frame-like list with \code{gc} (increasing, in [0, 1]) and \code{log} (log factor at gc).
#' @param gc_window integer(2), P-site codon offsets spanned by the footprint.
#' @param gc_strength numeric multiplier of the log factor.
#' @param study_unit "study" (default: the non-shared part is drawn per \code{colData(count_table)$study}, one study for
#'  all libraries when the column is missing) or "library" (drawn per library: every simulated library is its own
#'  study pool, as the two study-disjoint halves of the merged real track are).
#' @param seed integer seed of the field draws.
#' @return a validated list
#' @examples
#' residual_field_options()$active
#' f <- residual_field_options(field_var = c(0.1, 0.05, 0.03), smear = c(1, 0.3))
#' f$active
#' @export
residual_field_options <- function(field_var = c(0, 0, 0), field_scales = c(3, 30, 300), shared = 0.72,
                                   gene_sd = 0, smear = NULL, gc_bias = FALSE,
                                   gc_knots = list(gc = c(0.37, 0.45, 0.60, 0.70), log = c(-0.13, 0.06, 0.06, -0.13)),
                                   gc_window = c(-4L, 5L), gc_strength = 1, study_unit = "study", seed = 1L) {
  o <- list(field_var = as.numeric(field_var), field_scales = as.numeric(field_scales), shared = shared,
            study_unit = study_unit,
            gene_sd = gene_sd, smear = if (!is.null(smear)) as.numeric(smear), gc_bias = gc_bias,
            gc_knots = list(gc = as.numeric(gc_knots$gc), log = as.numeric(gc_knots$log)),
            gc_window = as.integer(gc_window), gc_strength = gc_strength, seed = as.integer(seed))
  stopifnot(length(o$field_var) == length(o$field_scales), all(o$field_var >= 0), all(o$field_scales > 0),
            o$shared >= 0, o$shared <= 1, o$gene_sd >= 0, is.logical(o$gc_bias),
            length(o$gc_knots$gc) == length(o$gc_knots$log), length(o$gc_knots$gc) >= 1L,
            !is.unsorted(o$gc_knots$gc), length(o$gc_window) == 2L, o$gc_window[1] <= o$gc_window[2],
            is.numeric(o$gc_strength), o$study_unit %in% c("study", "library"))
  if (!is.null(o$smear)) stopifnot(length(o$smear) >= 1L, all(o$smear >= 0), o$smear[1] > 0)
  o$active <- any(o$field_var > 0) || !is.null(o$smear) || isTRUE(o$gc_bias)
  o
}

# One stationary AR(1) path of length n with variance v and correlation exp(-1 / L) per step.
field_ar1 <- function(n, L, v) {
  if (n < 1L || v <= 0) return(numeric(n))
  phi <- exp(-1 / L)
  e <- stats::rnorm(n, 0, sqrt(v * (1 - phi^2)))
  e[1] <- stats::rnorm(1, 0, sqrt(v))
  as.numeric(stats::filter(e, phi, method = "recursive"))
}

# Sum of the AR(1) components for genes of n codons (one vector per gene), drawn with a fixed seed.
field_paths <- function(n, f, seed) {
  with_gene_seed(list(per_genome_draws = TRUE, seed = seed), length(n), function(k) lapply(n, function(m) {
    out <- numeric(m)
    for (j in seq_along(f$field_var)) out <- out + field_ar1(m, f$field_scales[j], f$field_var[j])
    out
  }))
}

# Log-field per gene for one study: s_g (sqrt(shared) F_shared + sqrt(1 - shared) F_study).
draw_residual_field <- function(n, f, study = "study1") {
  if (!any(f$field_var > 0)) return(lapply(n, numeric))
  sh <- if (f$shared > 0) field_paths(n, f, f$seed + 6000L)
  st <- if (f$shared < 1) field_paths(n, f, f$seed + 6500L + name_seed(study))
  amp <- with_gene_seed(list(per_genome_draws = TRUE, seed = f$seed + 6900L), length(n), function(k)
    exp(f$gene_sd * stats::rnorm(k) - f$gene_sd^2))
  lapply(seq_along(n), function(i) {
    x <- numeric(n[i])
    if (!is.null(sh)) x <- x + sqrt(f$shared) * sh[[i]]
    if (!is.null(st)) x <- x + sqrt(1 - f$shared) * st[[i]]
    amp[i] * x
  })
}

# Symmetric codon smear: y_i = sum_k w_|k| x_{i+k} / (sum of the kernel weights inside the gene).
apply_codon_smear <- function(x, w) {
  if (is.null(w) || length(w) == 1L || length(x) < 2L) return(x)
  kern <- c(rev(w[-1]), w)
  pad <- rep(0, length(w) - 1L)
  num <- rolling_weighted_sum(c(pad, x, pad), kern)
  den <- rolling_weighted_sum(c(pad, rep(1, length(x)), pad), kern)
  keep <- length(pad) + seq_along(x)
  num[keep] / den[keep]
}

# GC fraction of the footprint of every P-site codon p (codons p + window[1] .. p + window[2], cut at the gene ends),
# returned in weight-vector order (index p + 1 - a_shift); weights without a P-site codon get the gene mean.
footprint_gc <- function(codons, a_shift = 1L, window = c(-4L, 5L)) {
  cc <- toupper(as.character(codons))
  gc <- (nchar(gsub("[^GC]", "", cc))) / pmax(nchar(cc), 1L)
  n <- length(gc); cs <- c(0, cumsum(gc))
  lo <- pmax(seq_len(n) + window[1], 1L); hi <- pmin(seq_len(n) + window[2], n)
  fp <- (cs[hi + 1L] - cs[lo]) / (hi - lo + 1L)
  out <- rep(mean(gc), n)
  idx <- psite_to_index(n, a_shift); ok <- !is.na(idx)
  out[idx[ok]] <- fp[ok]
  out
}

# Log factor of a footprint GC fraction (piecewise linear through the knots, constant outside).
gc_log_factor <- function(gc, f) {
  k <- f$gc_knots
  if (length(k$gc) == 1L) return(rep(k$log * f$gc_strength, length(gc)))
  f$gc_strength * stats::approx(k$gc, k$log, xout = gc, rule = 2)$y
}

# Steps 1f-1h of sim_sequence_bias(): footprint GC, log-field and smear for all genes of one library.
apply_residual_field <- function(res, f, codons = NULL, a_shift = 1L, study = "study1") {
  if (is.null(f) || !isTRUE(f$active)) return(res)
  if (isTRUE(f$gc_bias) && !is.null(codons)) {
    if (length(codons) != length(res)) stop("codons do not match the genes of the codon weights")
    res <- lapply(seq_along(res), function(i) {
      if (length(codons[[i]]) != length(res[[i]])) stop("codons do not match the codon weights of a gene")
      res[[i]] * exp(gc_log_factor(footprint_gc(codons[[i]], a_shift, f$gc_window), f))
    })
  }
  if (any(f$field_var > 0)) {
    lf <- draw_residual_field(lengths(res), f, if (is.null(study) || is.na(study)) "study1" else study)
    res <- lapply(seq_along(res), function(i) res[[i]] * exp(lf[[i]]))
  }
  if (!is.null(f$smear)) res <- lapply(res, apply_codon_smear, w = f$smear)
  res
}

# Per-gene contrast exponents N(mean, sd), floored at 0.05, drawn per genome (seed + 5000).
draw_contrast_exponents <- function(n, cb) {
  with_gene_seed(list(per_genome_draws = TRUE, seed = cb$seed + 5000L), n, function(k)
    pmax(stats::rnorm(k, cb$contrast_mean, cb$contrast_sd), 0.05))
}

# Study of the (single) library of a count table, for the study part of the residual field: the library name
# (study_unit "library"), else colData$study, else one study for all libraries.
field_study <- function(count_table, unit = codon_bias_options()$field$study_unit) {
  if (identical(unit, "library")) return(as.character(colnames(count_table)[1]))
  st <- SummarizedExperiment::colData(count_table)$study
  if (is.null(st) || !length(st) || is.na(st[1])) "study1" else as.character(st[1])
}
