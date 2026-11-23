#' Default read-end tables per library method and nuclease
#'
#' Real Ribo-seq read ends carry two separate sequence preferences: where the
#' nuclease cut (the cut site, i.e. the nucleotide just outside the read) and how
#' the library chemistry treated the footprint ends (ligase, reverse
#' transcriptase, template-switching oligo: the first and last nucleotide inside
#' the read). This returns both as \code{fragment_geometry} entries, so a
#' simulation can look like a library made with a given method and nuclease
#' without learning from a BAM first.
#'
#' The tables are medians over public studies (each study counted once) of
#' per-library fits with ORFik's \code{codon_dwell_end_bias(end_nt = 2, flank =
#' TRUE)} on unclipped reads, from libraries whose read ends were not damaged by
#' trimming (class D of the 2026-11-22 trimming audit). The 5' terminal term of
#' the added non-templated base was corrected for hidden added bases first, so
#' that it is not counted twice when \code{added_base = TRUE} simulates the base
#' itself. Weights are relative (geometric mean 1 per end). The number of studies
#' behind each table is in \code{attr(result, "evidence")}; the non-RNase-I
#' nuclease tables rest on one study each and are a first approximation. Template
#' switching: one of its two studies poly(A)-tails the footprints, so its 3'
#' terminal A is strongly depleted (trimming the tail removes genomic A's).
#'
#' @param method library method: "circularization", "dual_ligation",
#'   "template_switching" or "ribolace" (RiboLace pull-down, circularization
#'   library), or NULL for no chemistry table.
#' @param nuclease "RNase I" (default; the same in all methods), "RNase T1",
#'   "RNase A", "RNase S7", "RNase T1 + S7", "MNase", or NULL for no cut-site
#'   table.
#' @param added_base TRUE to also return the non-templated 5' base of the method
#'   (\code{fragment_geometry$added_base}: base and rate, the median over
#'   studies of the share of reads that got it). Circularization T 10 percent,
#'   RiboLace T 16 percent; dual ligation and template switching 0 (no base
#'   dominates their 1-nt 5' clips). FALSE (default) leaves it out.
#' @param strength exponent applied to every table (1 = as measured, 0 = off).
#' @return a list with any of five_prime_bias, three_prime_bias (chemistry),
#'   five_prime_flank_bias, three_prime_flank_bias (nuclease) and added_base, to
#'   merge into \code{fragment_geometry}, e.g.
#'   \code{utils::modifyList(list(source = "default"), method_end_bias("circularization"))}.
#' @export
method_end_bias <- function(method = c("circularization", "dual_ligation",
                                       "template_switching", "ribolace"),
                            nuclease = "RNase I", added_base = FALSE,
                            strength = 1) {
  tables <- data.table::fread(system.file("extdata", "end_bias_method_tables.csv",
                                          package = "coverageSim", mustWork = TRUE))
  if (!is.null(method)) method <- match.arg(method)
  if (!is.null(nuclease)) {
    available <- tables[table == "nuclease", unique(group)]
    if (length(nuclease) != 1L || !nuclease %in% available) {
      stop("nuclease must be NULL or one of: ", paste(available, collapse = ", "))
    }
  }
  if (!is.numeric(strength) || length(strength) != 1L || !is.finite(strength) || strength < 0) {
    stop("strength must be one finite non-negative number")
  }
  profile <- function(what, group_name, end_name) {
    rows <- tables[table == what & group == group_name & end == end_name]
    list(source = "user", table = rows[, .(kmer, weight)], strength = strength)
  }
  result <- list()
  evidence <- NULL
  if (!is.null(method)) {
    result$five_prime_bias <- profile("chemistry", method, "five")
    result$three_prime_bias <- profile("chemistry", method, "three")
    evidence <- rbind(evidence, tables[table == "chemistry" & group == method,
                                       .(studies = studies[1], libraries = libraries[1]), by = .(table, group, end)])
  }
  if (!is.null(nuclease)) {
    result$five_prime_flank_bias <- profile("nuclease", nuclease, "five")
    result$three_prime_flank_bias <- profile("nuclease", nuclease, "three")
    evidence <- rbind(evidence, tables[table == "nuclease" & group == nuclease,
                                       .(studies = studies[1], libraries = libraries[1]), by = .(table, group, end)])
  }
  if (isTRUE(added_base)) {
    if (is.null(method)) stop("added_base = TRUE needs a method")
    added <- data.table::fread(system.file("extdata", "end_bias_added_base.csv",
                                           package = "coverageSim", mustWork = TRUE))
    chosen <- method
    row <- added[method == chosen]
    result$added_base <- list(base = row$base, rate = row$rate)
  }
  attr(result, "evidence") <- evidence
  result
}
