#' Realistic RFP settings (opt-in) learned from a real human Ribo-seq library
#'
#' Returns the extra \code{simNGScoverage()} arguments validated in backclaude t69 (null-gene false-flag test, 3
#' seeds, mid/high depth). Use them together with \code{options(coverageSim.codon_bias = codon_bias_preset("realistic"))}.
#' Nothing changes unless you pass them: package defaults are untouched.
#'
#' The fragment geometry (\code{source = "learned"}, \code{site_reference = "a_site"}) uses the length/offset
#' distribution learned from SRR12651187 (6 M read subsample, 2264 canonical transcripts,
#' \code{\link{learn_fragment_geometry}}), keeping every offset with at least 5\% of the length's best in-frame
#' fraction (weight = learned length probability x in-frame fraction), plus the by-length frame bias
#' (\code{learn_end_bias(by_length = TRUE)$frame_bias}) at \code{strength = 0.25} (full strength over-sharpens frame 0).
#' The data ship in \code{inst/extdata/realistic_rfp_geometry_SRR12651187.rds}.
#'
#' Caveat: with these settings the simulated replicate-to-replicate window band is narrower than in real data
#' (false density flags 0.2-0.3\% against a 5\% nominal level), so a false-flag test on simulated null genes is
#' conservative. Clumping is 3.4x Poisson (real 2.7x). See \code{docs/REALISTIC_RFP.md}.
#' @param frame_strength numeric, strength of the learned frame bias (default 0.25 as validated).
#' @return a list with \code{dmn_alpha_scale} (6) and \code{fragment_geometry}, to pass with
#'  \code{do.call(simNGScoverage, c(list(...), realistic_rfp_settings()))}. \code{dmn_gdm_scale} stays at its default 1.
#' @examples
#' options(coverageSim.codon_bias = codon_bias_preset("realistic"))
#' s <- realistic_rfp_settings()
#' names(s)
#' @export
realistic_rfp_settings <- function(frame_strength = 0.25) {
  path <- system.file("extdata", "realistic_rfp_geometry_SRR12651187.rds", package = "coverageSim")
  if (!nzchar(path)) stop("realistic_rfp_geometry_SRR12651187.rds not found in the coverageSim installation", call. = FALSE)
  learned <- readRDS(path)
  list(
    dmn_alpha_scale = 6,
    fragment_geometry = list(
      source = "learned",
      site_reference = "a_site",
      distribution = data.table::as.data.table(learned$distribution),
      frame_bias = list(source = "learned", table = data.table::as.data.table(learned$frame_bias_table),
                        strength = frame_strength)
    )
  )
}
