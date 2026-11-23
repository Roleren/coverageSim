# Difficulty presets (Carina's item 1, 2026-11-23): "clean" (easier than real data), "realistic" (human_v2, the
# calibrated human setting) and "noisy" (harder than real data, with edge cases), to test ORF callers and metrics at
# known difficulty levels. Each preset is one documented parameter set over all simulator layers. Validation against
# real human libraries: docs/PRESETS.md.

# Fragment geometry with every read length's offsets collapsed onto its most likely offset (clean preset): after a
# single per-length P-site offset every read lands in its true frame, as in ultraperiodic libraries.
single_offset_geometry <- function(distribution) {
  d <- data.table::as.data.table(distribution)
  best <- d[, .(site_offset = site_offset[which.max(probability)], probability = sum(probability)),
            by = fragment_length]
  best[, probability := probability / sum(probability)]
  best[]
}

# Fragment geometry with a flat offset spread (noisy preset): every length keeps its offsets and gets the two
# neighbouring offsets, with the off-peak offsets lifted to `floor` of the length's best offset, and an extra
# short-read tail (20-24 nt, offsets 8-12) carrying `short_share` of the reads.
spread_offset_geometry <- function(distribution, floor = 0.6, short_share = 0.12) {
  d <- data.table::as.data.table(distribution)
  ext <- d[, .(site_offset = seq(min(site_offset) - 1L, max(site_offset) + 1L)), by = fragment_length]
  d <- d[ext, on = c("fragment_length", "site_offset")]
  d[is.na(probability), probability := 0]
  d[, probability := pmax(probability, floor * max(probability)), by = fragment_length]
  d[, probability := probability / sum(probability) * (1 - short_share)]
  short <- data.table::CJ(fragment_length = 20:24, site_offset = 8:12)
  short[, probability := short_share / .N]
  out <- rbind(short, d)
  out[, probability := probability / sum(probability)]
  data.table::setorder(out, fragment_length, site_offset)
  out[]
}

#' Simulation difficulty presets: clean, realistic and noisy
#'
#' One parameter set per difficulty level over every simulator layer (genome, gene counts, region shares, codon
#' structure, fragment geometry and frame, technical artefacts). "realistic" is the calibrated human setting
#' (\code{codon_bias_preset("human_v2")}, \code{realistic_rfp_settings()}, \code{utr3_preset("human")}, B1 region
#' shares). "clean" is easier than real data: high periodicity (one P-site offset per read length, as in the
#' ultraperiodic human libraries with ~90\% frame 0), little codon-level noise, low gene dispersion, few UTR reads,
#' no spikes, no multimapping. "noisy" is harder than real data: low periodicity (flat offset spread, short-read
#' tail, no frame bias), strong codon-level noise (doubled residual field, footprint GC bias, zero mask), high gene
#' dispersion, UTR background (contamination-like reads outside the CDS), spikes, multimapping and PCR duplicates,
#' and gene architectures with more uORFs. Where each preset falls relative to real human libraries:
#' \code{docs/PRESETS.md}.
#'
#' @param name "clean", "realistic" or "noisy".
#' @return a list of class "simulation_preset": \code{name}, \code{description}, \code{options} (session options to
#'  set: \code{coverageSim.codon_bias}, \code{coverageSim.utr3}), \code{genome} (extra \code{simGenome()} arguments),
#'  \code{count_tables} (\code{simCountTables()}), \code{regions} (\code{simCountTablesRegions()}) and
#'  \code{coverage} (\code{simNGScoverage()}). Use \code{\link{run_simulation_preset}} or the parts by hand:
#'  \preformatted{p <- simulation_preset("noisy")
#' withr::with_options(p$options, {
#'   g <- do.call(simGenome, c(list(n = 500, out_dir = dir), p$genome))
#'   ct <- do.call(simCountTables, c(list(n = ORFik::loadRegion(g["txdb"], "cds")), p$count_tables))
#'   rc <- do.call(simCountTablesRegions, c(list(count_table = ct), p$regions))
#'   do.call(simNGScoverage, c(list(simGenome = g, count_table = rc), p$coverage))
#' })}
#' @export
#' @examples
#' p <- simulation_preset("clean")
#' p$description
#' names(p$coverage)
simulation_preset <- function(name = c("realistic", "clean", "noisy")) {
  name <- match.arg(name)
  real <- realistic_rfp_settings()
  geo <- real$fragment_geometry
  v2 <- codon_bias_preset("human_v2")
  b1 <- list(leader = list(RFP = 0.025, RNA = 0.3, CAGE = 1, PAS = 0),
             cds = list(RFP = 0.938, RNA = 0.4, CAGE = 0, PAS = 0),
             trailer = list(RFP = 0.017, RNA = 0.3, CAGE = 0, PAS = 1),
             uorf = list(RFP = 0.02, RNA = 0, CAGE = 0, PAS = 0))
  # the learned geometry with another offset distribution and no frame bias
  with_geometry <- function(g, distribution) {
    g$distribution <- distribution
    g$frame_bias <- list(source = "none")
    g
  }
  regions_with <- function(leader, cds, trailer, uorf) {
    r <- b1
    r$leader$RFP <- leader; r$cds$RFP <- cds; r$trailer$RFP <- trailer; r$uorf$RFP <- uorf
    r
  }
  p <- switch(name,
    clean = list(
      description = paste("Easier than real data: one P-site offset per read length and a narrow RNase kernel",
                          "(frame 0 ~0.9 like the ultraperiodic human libraries), mild codon dwell without pauses or residual field,",
                          "low gene dispersion, 98% CDS reads, no spikes, no multimapping."),
      options = list(coverageSim.codon_bias = codon_bias_preset("clean"),
                     coverageSim.utr3 = list()),
      genome = list(length_model = "human", max_uorfs = 1, on_infeasible = "warning"),
      count_tables = list(dispMeanRel = function(x) 1 / x + 0.001),
      regions = list(region_proportion = regions_with(0.006, 0.98, 0.006, 0.008), region_dmn_concentration = 50),
      coverage = list(dmn_alpha_scale = 30,
                      # RNase kernel with little smear (default c(0.5, 2, 1, 10, 2, 1, 0.5) puts 41% off the site)
                      rnase_bias = list(RFP = c(0.05, 0.4, 0.05, 10, 0.4, 0.05, 0.05), RNA = 1, CAGE = 1, PAS = 1),
                      fragment_geometry = with_geometry(geo, single_offset_geometry(geo$distribution)),
                      technical_artifacts = NULL)),
    realistic = list(
      description = paste("Calibrated human setting: codon_bias_preset('human_v2'), realistic_rfp_settings()",
                          "(learned geometry and frame bias of SRR12651187), utr3_preset('human'), B1 region",
                          "shares (median of 757 human libraries), 15% multimapped reads."),
      options = list(coverageSim.codon_bias = v2, coverageSim.utr3 = utr3_preset("human")),
      genome = list(length_model = "human", max_uorfs = 1, on_infeasible = "warning"),
      count_tables = list(dispMeanRel = function(x) 4 / x + 0.01),
      regions = list(region_proportion = b1, region_dmn_concentration = 1),
      coverage = c(real, list(technical_artifacts = list(multimapping_rate = 0.15, secondary_alignments = 1L)))),
    noisy = list(
      description = paste("Harder than real data: flat P-site offset spread, a wide RNase kernel and a 20-24 nt",
                          "short-read tail (low periodicity), doubled residual field with footprint GC bias and a 30% zero mask,",
                          "high gene dispersion, 25% of reads outside the CDS (contamination-like background),",
                          "spikes in a third of the libraries, 35% multimapped reads, 20% PCR duplicates,",
                          "up to 3 uORFs per leader."),
      options = list(coverageSim.codon_bias = codon_bias_preset("noisy"),
                     coverageSim.utr3 = utils::modifyList(utr3_preset("human"), list(plateau_prob = 0.12))),
      genome = list(length_model = "human", max_uorfs = 3, on_infeasible = "warning"),
      count_tables = list(dispMeanRel = function(x) 8 / x + 0.4),
      regions = list(region_proportion = regions_with(0.09, 0.75, 0.13, 0.03), region_dmn_concentration = 0.5),
      coverage = list(dmn_alpha_scale = 1.5,
                      # wider RNase smear than the default kernel
                      rnase_bias = list(RFP = c(1, 3, 2, 10, 3, 2, 1), RNA = 1, CAGE = 1, PAS = 1),
                      fragment_geometry = with_geometry(geo, spread_offset_geometry(geo$distribution)),
                      technical_artifacts = list(multimapping_rate = 0.35, secondary_alignments = 2L,
                                                 duplication_rate = 0.2, duplicate_copies = 1L,
                                                 spike_sites = 4L, spike_library_fraction = 0.34,
                                                 spike_multiplicity = c(300, 3000)))))
  structure(c(list(name = name), p), class = "simulation_preset")
}

#' @export
print.simulation_preset <- function(x, ...) {
  cat("coverageSim simulation preset '", x$name, "'\n", strwrap(x$description, prefix = "  "), sep = "\n")
  cat("  layers: options (", toString(names(x$options)), "), genome, count_tables, regions, coverage (",
      toString(names(x$coverage)), ")\n", sep = "")
  invisible(x)
}

#' Run a whole RFP simulation with a difficulty preset
#'
#' Builds a genome, gene and region count tables and the reads with one \code{\link{simulation_preset}}. The
#' preset's options are set only for the duration of the call. RFP output is BAM when the preset has technical
#' artefacts (they need SAM/BAM), else ofst.
#' @param preset a \code{simulation_preset()} result or its name.
#' @param n number of genes.
#' @param out_dir directory for genome, reads and experiment (created).
#' @param conditions,replicates libraries: \code{length(conditions) * replicates} RFP libraries.
#' @param interceptMean log mean count per gene (\code{simCountTables()}).
#' @param seed RNG seed.
#' @param ... further \code{simNGScoverage()} arguments (override the preset).
#' @return the ORFik experiment of the simulated libraries (invisibly), with the genome in attribute "genome".
#' @export
run_simulation_preset <- function(preset = "realistic", n = 500, out_dir = tempfile("preset-"),
                                  conditions = "WT", replicates = 2, interceptMean = 6, seed = 1, ...) {
  if (is.character(preset)) preset <- simulation_preset(preset)
  stopifnot(inherits(preset, "simulation_preset"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  old <- options(preset$options)
  on.exit(options(old), add = TRUE)
  set.seed(seed)
  g <- do.call(simGenome, c(list(n = n, out_dir = file.path(out_dir, "genome"), debug_on = FALSE), preset$genome))
  cds <- ORFik::loadRegion(g["txdb"], "cds")
  ct <- do.call(simCountTables, c(list(n = cds, libtypes = "RFP", conditions = conditions, replicates = replicates,
                                       interceptMean = interceptMean, plot_PCA = FALSE, print_statistics = FALSE),
                                  preset$count_tables))
  rc <- do.call(simCountTablesRegions, c(list(count_table = ct), preset$regions))
  args <- preset$coverage
  if (has_active_technical_artifacts(normalize_technical_artifacts(args$technical_artifacts)))
    args$libFormats <- list(RFP = "bam")
  dots <- list(...)
  args[names(dots)] <- dots
  reads <- file.path(out_dir, "reads"); exps <- file.path(out_dir, "experiment")
  dir.create(reads, showWarnings = FALSE); dir.create(exps, showWarnings = FALSE)
  e <- do.call(simNGScoverage, c(list(simGenome = g, count_table = rc, out_dir = reads,
                                      exp_name = paste0("preset_", preset$name), exp_save_dir = exps,
                                      validate = FALSE), args))
  attr(e, "genome") <- g
  invisible(e)
}
