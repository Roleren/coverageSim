#' Codon dwell by context, study, tissue and library (opt-in, explainability port S12)
#'
#' Set with \code{options(coverageSim.codon_dwell = list(...))}; unset fields keep the defaults, which leave the
#' simulation unchanged. The layers multiply the codon-level weights of the codon-biased regions (CDS, uORFs) after the
#' smoothing step of the sequence bias (where the replicate-shared pauses of \code{\link{codon_bias_options}} also act),
#' as \code{exp(sum_k M[codon at i + k, k])} for the decoding-site codon i and offsets k = -10..+5 relative to the A site
#' (for P-site placements the A site is i + 1). They need a codon-level \code{seq_bias} (the default "AUTO" table).
#' Values ship in \code{inst/extdata/codon_dwell_human.csv} (built by \code{analysis/R/build_codon_dwell_human.R}).
#'
#' \describe{
#'  \item{global}{"h1": the H1 codon-dwell GLM of the explainer (16 offsets x 60 codon contrasts, reference AAA; Codex
#'   task explainer_h1_codon_dwell_20261118/model.rds), "atlas": the consensus of the human dwell atlas (median of 367
#'   study means, 2,960 libraries, claude-megafst m6_study_fingerprint/beta_matrix_w10_5.rds), "none" (default).
#'   The bundled A-site \code{seq_bias} tables already carry the A-site effect, including the acidic one (D 2.0x,
#'   E 1.47x of the median amino acid; E13-E: observed E/D density 1.39x vs 1.35x expected by H1), so with them use
#'   \code{global_offsets = setdiff(-10:5, 0)} (the preset does) to add only the context offsets.}
#'  \item{study_dwell}{Each study (\code{colData(count_table)$study}; one study for all libraries when the column is
#'   missing) gets its own codon x offset deviation, N(0, (study_scale * study_sd_k)^2), and each library a smaller
#'   one, N(0, (library_scale * library_sd_k)^2), with the per-offset sds of the atlas. By default the two scales are
#'   solved from the target correlations of library dwell profiles: \code{between_r} for libraries of different
#'   studies and \code{within_r} for libraries of one study (human: same study 0.85-0.89, other study 0.30-0.35,
#'   study explains 0.81 of the library variance; M6, m6_study_fingerprint/pair_correlation_classes.csv,
#'   variance_decomposition.csv). The real within-study value includes estimation noise, which a simulated library
#'   adds again through its finite depth, so measured simulated profiles correlate a little less.}
#'  \item{tissue groups}{Libraries with \code{colData(count_table)$tissue_group} get that group's delta added.
#'   Built in: "striated" (muscle + heart, M1-a / M1-c, dwell_by_group/muscle_detail.csv: E/P/A deltas, all four Pro
#'   A-site codons +0.10..+0.24, CTT/CTA slower, CAG/CAC faster). Further groups via \code{tissue_deltas}, a named list
#'   of data.frames with columns offset, codon, delta.}
#'  \item{peaky libraries}{Libraries with \code{colData(count_table)$peaky == TRUE} raise their codon weights to
#'   \code{peaky_power} (more reads on the top codons) and divide their Dirichlet concentration by
#'   \code{peaky_var_ratio} (more window-level extra-Poisson variance). Human (M1-d): striated libraries have a 1.5x
#'   top-10-codon share on 16 housekeeping genes at matched depth (0.21-0.30 vs 0.14-0.18) and 1.59x window-level
#'   extra-Poisson variance (dwell_by_group/peakiness_housekeeping.csv, phi_housekeeping.csv).}
#' }
#' The per-library profiles used (offset x codon log effects) are written to \code{codon_dwell_profiles.tsv} in
#' \code{out_dir}. Study and library draws use seeds derived from \code{seed} and the study / library names, so they do
#' not depend on library order and leave the session random state untouched.
#' @param global "none", "h1" or "atlas": the global context dwell layer.
#' @param global_offsets integer offsets (relative to the A site) of the global layer that are used.
#' @param study_dwell logical, add study and library deviations.
#' @param between_r,within_r target correlations of library dwell profiles between and within studies.
#' @param study_scale,library_scale numeric multipliers of the atlas per-offset sds; NULL = solved from
#'  \code{between_r} / \code{within_r} (relative to the variance of the global profile, the atlas consensus when
#'  \code{global = "none"}).
#' @param tissue_deltas named list of extra tissue-group deltas (data.frames: offset, codon, delta).
#' @param peaky_power numeric > 0, power on the codon weights of peaky libraries.
#' @param peaky_var_ratio numeric >= 1, extra-Poisson variance ratio of peaky libraries.
#' @param seed integer seed of the study and library draws.
#' @return a list of options (with the solved \code{study_scale} and \code{library_scale})
#' @examples
#' options(coverageSim.codon_dwell = codon_dwell_preset("human"))
#' o <- codon_dwell_options()
#' c(o$study_scale, o$library_scale)
#' options(coverageSim.codon_dwell = NULL)
#' @export
codon_dwell_options <- function(global = "none", global_offsets = -10:5, study_dwell = FALSE, between_r = 0.37,
                                within_r = 0.88, study_scale = NULL, library_scale = NULL, tissue_deltas = NULL,
                                peaky_power = 1, peaky_var_ratio = 1, seed = 1L) {
  def <- list(global = global, global_offsets = global_offsets, study_dwell = study_dwell, between_r = between_r,
              within_r = within_r, study_scale = study_scale, library_scale = library_scale,
              tissue_deltas = tissue_deltas, peaky_power = peaky_power, peaky_var_ratio = peaky_var_ratio, seed = seed)
  set <- getOption("coverageSim.codon_dwell", list())
  bad <- setdiff(names(set), names(def)); if (length(bad)) stop("unknown coverageSim.codon_dwell fields: ", toString(bad))
  def[names(set)] <- set
  stopifnot(def$global %in% c("none", "h1", "atlas"), all(def$global_offsets %in% -10:5),
            is.logical(def$study_dwell), def$between_r > 0, def$within_r > def$between_r, def$within_r <= 1,
            def$peaky_power > 0, def$peaky_var_ratio >= 1, is.null(def$tissue_deltas) || is.list(def$tissue_deltas))
  if (is.null(def$study_scale) || is.null(def$library_scale)) {
    tab <- codon_dwell_table()
    ref <- if (def$global == "h1") tab$h1_beta else tab$atlas_consensus
    # pooled over all codon x offset cells, as the correlation of two profile vectors pools them
    vg <- mean((ref - mean(ref))^2)
    total <- vg / def$between_r
    if (is.null(def$study_scale))
      def$study_scale <- sqrt((def$within_r - def$between_r) * total / mean(tab$study_sd^2))
    if (is.null(def$library_scale))
      def$library_scale <- sqrt((1 - def$within_r) * total / mean(tab$library_sd^2))
  }
  def
}

#' Named codon dwell settings
#' @param name "original" (all off) or "human": H1 context offsets (A site left to the seq_bias table), study and
#'  library deviations with between / within-study profile correlation 0.37 / 0.88, peaky libraries with power 1.3
#'  and variance ratio 1.6 (see \code{\link{codon_dwell_options}}).
#' @return a list for \code{options(coverageSim.codon_dwell = ...)}
#' @export
codon_dwell_preset <- function(name = c("original", "human")) {
  switch(match.arg(name),
         original = list(),
         human = list(global = "h1", global_offsets = setdiff(-10:5, 0L), study_dwell = TRUE, between_r = 0.37,
                      within_r = 0.88, peaky_power = 1.3, peaky_var_ratio = 1.6))
}

#' Mark the study, tissue group and peakiness of simulated libraries
#'
#' Adds the colData columns read by the codon dwell layers (\code{\link{codon_dwell_options}}).
#' @param count_table a SummarizedExperiment from \code{simCountTables()} / \code{simCountTablesRegions()}.
#' @param study character, study (library group) per library, recycled.
#' @param tissue_group character, tissue group per library (NA = none), recycled.
#' @param peaky logical, peaky library flag, recycled.
#' @return count_table with columns study, tissue_group and peaky.
#' @examples
#' \dontrun{
#' ct <- assign_dwell_groups(ct, study = rep(c("S1", "S2"), each = 3), tissue_group = c(rep(NA, 3), rep("striated", 3)))
#' }
#' @export
assign_dwell_groups <- function(count_table, study, tissue_group = NA_character_, peaky = FALSE) {
  n <- ncol(count_table)
  SummarizedExperiment::colData(count_table)$study <- rep_len(as.character(study), n)
  SummarizedExperiment::colData(count_table)$tissue_group <- rep_len(as.character(tissue_group), n)
  SummarizedExperiment::colData(count_table)$peaky <- rep_len(as.logical(peaky), n)
  count_table
}

# The bundled human dwell table (offset, codon, h1_beta, atlas_consensus, study_sd, library_sd, striated_delta).
codon_dwell_table <- function() {
  path <- system.file("extdata", "codon_dwell_human.csv", package = "coverageSim")
  if (!nzchar(path)) stop("codon_dwell_human.csv not found in the coverageSim installation", call. = FALSE)
  data.table::fread(path)
}

# Deterministic integer from a name, for per-study / per-library seeds independent of library order.
name_seed <- function(x) {
  v <- utf8ToInt(enc2utf8(as.character(x)))
  as.integer(sum(v * seq_along(v) * 31L) %% 1000003L)
}

# Is any dwell layer active for these libraries?
codon_dwell_active <- function(opts, coldata) {
  opts$global != "none" || isTRUE(opts$study_dwell) ||
    (!is.null(coldata$tissue_group) && any(!is.na(coldata$tissue_group))) ||
    (!is.null(coldata$peaky) && any(coldata$peaky %in% TRUE))
}

#' Codon dwell profile of one simulated library
#'
#' The offset x codon log effect that \code{simNGScoverage()} applies for a library under
#' \code{\link{codon_dwell_options}}: global layer + study deviation + library deviation + tissue delta.
#' @param study character, the library's study.
#' @param library character, the library name.
#' @param tissue_group character or NA.
#' @param opts options, default \code{codon_dwell_options()}.
#' @return a data.table with offset, codon and log_effect (61 sense codons x 16 offsets).
#' @export
library_dwell_profile <- function(study, library, tissue_group = NA_character_, opts = codon_dwell_options()) {
  tab <- codon_dwell_table()
  eff <- numeric(nrow(tab))
  if (opts$global != "none") {
    g <- if (opts$global == "h1") tab$h1_beta else tab$atlas_consensus
    eff <- eff + ifelse(tab$offset %in% opts$global_offsets, g, 0)
  }
  if (isTRUE(opts$study_dwell)) {
    draw <- function(seed, sd) with_gene_seed(list(per_genome_draws = TRUE, seed = seed), length(sd),
                                              function(n) stats::rnorm(n, 0, sd))
    eff <- eff + draw(opts$seed + 7000L + name_seed(study), opts$study_scale * tab$study_sd)
    eff <- eff + draw(opts$seed + 8000L + name_seed(library), opts$library_scale * tab$library_sd)
  }
  if (!is.na(tissue_group)) {
    delta <- if (identical(tissue_group, "striated")) {
      tab$striated_delta
    } else if (!is.null(opts$tissue_deltas[[tissue_group]])) {
      d <- data.table::as.data.table(opts$tissue_deltas[[tissue_group]])
      stopifnot(all(c("offset", "codon", "delta") %in% names(d)))
      v <- d$delta[match(paste(tab$offset, tab$codon), paste(d$offset, d$codon))]
      ifelse(is.na(v), 0, v)
    } else stop("unknown tissue_group '", tissue_group, "': give its delta in codon_dwell_options(tissue_deltas)",
                call. = FALSE)
    eff <- eff + delta
  }
  data.table::data.table(offset = tab$offset, codon = tab$codon, log_effect = eff)
}

# Per-gene multiplicative factors for the decoding-site codons of one library.
# codons: list of codon vectors (one per gene); a_shift: 0 when position i is the A site, 1 when it is the P site.
dwell_codon_factors <- function(codons, profile, a_shift = 0L) {
  offsets <- sort(unique(profile$offset))
  lev <- sort(unique(profile$codon))
  M <- matrix(0, length(lev), length(offsets), dimnames = list(lev, offsets))
  M[cbind(match(profile$codon, lev), match(profile$offset, offsets))] <- profile$log_effect
  used <- which(colSums(M != 0) > 0)
  lapply(codons, function(cc) {
    n <- length(cc); ci <- match(cc, lev); eta <- numeric(n)
    for (j in used) {
      pos <- seq_len(n) + a_shift + offsets[j]
      ok <- pos >= 1L & pos <= n
      v <- M[ci[pos[ok]], j]; v[is.na(v)] <- 0
      eta[ok] <- eta[ok] + v
    }
    exp(eta)
  })
}

# Library dwell settings of one run: profiles per library (or NULL when off). Called once by simNGScoverage().
codon_dwell_run_setup <- function(coldata, libnames) {
  opts <- codon_dwell_options()
  if (!codon_dwell_active(opts, coldata)) return(NULL)
  study <- if (!is.null(coldata$study)) as.character(coldata$study) else rep("study1", length(libnames))
  tissue <- if (!is.null(coldata$tissue_group)) as.character(coldata$tissue_group) else rep(NA_character_, length(libnames))
  peaky <- if (!is.null(coldata$peaky)) coldata$peaky %in% TRUE else rep(FALSE, length(libnames))
  profiles <- lapply(seq_along(libnames), function(s) library_dwell_profile(study[s], libnames[s], tissue[s], opts))
  names(profiles) <- libnames
  list(opts = opts, study = stats::setNames(study, libnames), tissue_group = stats::setNames(tissue, libnames),
       peaky = stats::setNames(peaky, libnames), profiles = profiles)
}

# The dwell inputs of sim_sequence_bias() for one library and region.
library_dwell_for_region <- function(dwell, codons, library, a_shift = 0L) {
  none <- list(factors = NULL, power = 1, var_ratio = 1)
  if (is.null(dwell) || !library %in% names(dwell$profiles)) return(none)
  peaky <- isTRUE(dwell$peaky[[library]])
  list(factors = if (!is.null(codons)) dwell_codon_factors(codons, dwell$profiles[[library]], a_shift),
       power = if (peaky) dwell$opts$peaky_power else 1,
       var_ratio = if (peaky) dwell$opts$peaky_var_ratio else 1)
}

write_codon_dwell_profiles <- function(dwell, out_dir) {
  if (is.null(dwell)) return(invisible(NULL))
  tab <- data.table::rbindlist(lapply(names(dwell$profiles), function(l) {
    data.table::data.table(library = l, study = dwell$study[[l]], tissue_group = dwell$tissue_group[[l]],
                           peaky = dwell$peaky[[l]], dwell$profiles[[l]])
  }))
  data.table::fwrite(tab, file.path(out_dir, "codon_dwell_profiles.tsv"), sep = "\t")
  invisible(tab)
}
