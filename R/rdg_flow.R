# Ribosome flow on the ribosome decision graph (RDG plan step 4, 2026-11-18).
#
# Each transcript gets a scanning flux of 1 at the cap. Scanning ribosomes meet
# the start codons of the translated features and the CDS in 5' -> 3' order. At
# each start a fraction initiates (by start codon), the rest scan on (leaky
# scanning). Ribosomes that finish a short ORF can resume scanning
# (reinitiation; less likely after longer ORFs) but can only initiate again
# after re-acquiring an initiator ternary complex, modelled as a minimum
# distance downstream of the stop. A small fraction reinitiates after the CDS.
# A feature's ribosome load is its initiations times its length; loads set how
# a gene's translated reads are split between the CDS and its translons.
# Conditions change the parameters per library (e.g. the integrated stress
# response lengthens the distance needed to reinitiate, so ribosomes skip an
# ATF4-like uoORF and reach the CDS).

#' Parameters of the ribosome flow model
#'
#' @param p_init initiation probability per start codon for a scanning
#'  ribosome (named numeric); codons not listed use \code{p_init_other}.
#' @param p_init_other initiation probability of start codons not in \code{p_init}.
#' @param reinit_max reinitiation probability after a very short ORF.
#' @param reinit_scale_codons length scale: reinitiation after an ORF of
#'  \code{k} codons is \code{reinit_max * exp(-k / reinit_scale_codons)}.
#' @param reinit_distance nt after a stop codon before reinitiating ribosomes
#'  can initiate again (ternary complex re-acquisition).
#' @param reinit_after_cds reinitiation probability after the CDS stop.
#' @param p_frameshift default fraction of CDS ribosomes that shift frame at an
#'  "FS" feature (a feature's own \code{efficiency} overrides it).
#' @param p_readthrough default fraction of CDS ribosomes that read through the
#'  CDS stop into a "CTE" feature.
#' @param p_init_translon one initiation probability (default 0.3) for the start of the extra translated features
#'  whose class is in \code{translon_classes}, used instead of \code{p_init}; the CDS start and other classes keep
#'  \code{p_init} (also a per-condition \code{p_init}; \code{p_init_scale} still applies). NULL gives the behaviour
#'  before 2026-11 (every start uses \code{p_init} by its start codon). Fitted on 62 human marker genes with
#'  Ribo-seq-detected uORFs (clean-CDS usage, 5-fold cross-validation by gene): a codon-specific fit drove
#'  near-cognate starts to ~1 (selection effect: detected translons are used whatever their start codon), while one
#'  shared value of 0.3 cut the held-out RMSE on the logit scale from 1.67 (NULL) to 1.00.
#'  With 0.3 a translated uORF captures few ribosomes, so delayed reinitiation (ATF4-like ISR via
#'  \code{reinit_distance}) moves little to the CDS; for that regulation use NULL and a strong uORF start.
#' @param translon_classes feature classes that \code{p_init_translon} applies to, default the classes that start in
#'  the leader and were in the fit: \code{c("uORF", "uoORF", "NTE")}. Downstream and internal ORFs (dORF, doORF,
#'  intORF) were not in the fitted panel and keep \code{p_init}.
#' @param translon_gene_sd numeric >= 0, default 0 (off). sd of a per-gene log-normal factor (mean 1) on the initiation
#'  probabilities of the gene's extra features in \code{translon_classes}, drawn once per \code{simRDGflow()} call (shared by all libraries and
#'  conditions; stored in \code{metadata(result)$rdg_gene_factor}). Real genes vary much more than uORF count and
#'  overlap predict (Kozak context, uORF peptide, structure); 0.5 matches the spread of real clean-CDS usage in the
#'  same 62 genes (variance of logit usage 1.66 vs 1.69 real).
#' @param uoorf_cds_factor numeric in (0, 1], default 0.8 (1 = off, the behaviour before 2026-11). Extra factor on
#'  the CDS initiations for every overlapping uoORF of the gene, on top of the ribosomes the uoORF already takes from
#'  the scanning flux. Overlapping uoORFs repress the clean CDS more than the flow model alone predicts: with
#'  \code{p_init_translon = 0.3}, 0.8 lowers the held-out RMSE (logit) to 0.93 and moves the uoORF effect (logit
#'  slope) from -0.65 to -1.01 (real -1.15).
#' @param conditions named list: per condition (as in the count table's
#'  \code{condition} column) a list of parameters that replace the defaults, or
#'  \code{p_init_scale}, a factor on all initiation probabilities. Default none:
#'  every condition uses the defaults. Example for an integrated stress response
#'  condition: \code{list(ISR = list(reinit_distance = 150, p_init_scale = 0.8))}.
#' @return a list of class "rdg_parameters"
#' @export
#' @examples
#' rdg_parameters(conditions = list(ISR = list(reinit_distance = 150)))
rdg_parameters <- function(p_init = c(ATG = 0.8, CTG = 0.12, GTG = 0.07, TTG = 0.05,
                                      ACG = 0.04),
                           p_init_other = 0.02, reinit_max = 0.6,
                           reinit_scale_codons = 30, reinit_distance = 50,
                           reinit_after_cds = 0.05, p_frameshift = 0.05,
                           p_readthrough = 0.02,
                           p_init_translon = 0.3, translon_classes = c("uORF", "uoORF", "NTE"),
                           translon_gene_sd = 0, uoorf_cds_factor = 0.8,
                           conditions = list()) {
  params <- list(p_init = p_init, p_init_other = p_init_other, reinit_max = reinit_max,
                 reinit_scale_codons = reinit_scale_codons,
                 reinit_distance = reinit_distance, reinit_after_cds = reinit_after_cds,
                 p_frameshift = p_frameshift, p_readthrough = p_readthrough,
                 p_init_translon = p_init_translon,
                 translon_classes = translon_classes, translon_gene_sd = translon_gene_sd, uoorf_cds_factor = uoorf_cds_factor,
                 p_init_scale = 1, conditions = conditions)
  probs <- c(p_init, p_init_other, reinit_max, reinit_after_cds, p_frameshift, p_readthrough,
             p_init_translon)
  if (any(!is.finite(probs) | probs < 0 | probs > 1))
    stop("Probabilities in rdg_parameters must be within [0, 1]")
  if (!is.null(p_init_translon) && length(p_init_translon) != 1L)
    stop("p_init_translon must be NULL or one probability")
  if (length(translon_gene_sd) != 1L || !is.finite(translon_gene_sd) || translon_gene_sd < 0)
    stop("translon_gene_sd must be one number >= 0")
  if (length(uoorf_cds_factor) != 1L || !is.finite(uoorf_cds_factor) || uoorf_cds_factor <= 0 || uoorf_cds_factor > 1)
    stop("uoorf_cds_factor must be one number in (0, 1]")
  if (reinit_scale_codons <= 0 || reinit_distance < 0)
    stop("reinit_scale_codons must be > 0 and reinit_distance >= 0")
  structure(params, class = "rdg_parameters")
}

# Parameters for one condition.
rdg_condition_parameters <- function(params, condition) {
  over <- params$conditions[[as.character(condition)]]
  if (is.null(over)) return(params)
  for (nm in names(over)) params[[nm]] <- over[[nm]]
  params
}

# Ribosome load of the CDS and each translated feature of one transcript.
# starts/ends/classes/start_codons describe the translated extra features
# (transcript coordinates); L = leader length, C = CDS end (last stop base).
# gene_factor: the gene's factor on the initiation probabilities of its extra features (translon_gene_sd).
rdg_gene_loads <- function(starts, ends, classes, start_codons, L, C, params,
                           cds_start_codon = "ATG", efficiency = NULL, gene_factor = 1) {
  # Frameshift segments and readthrough extensions are fed by ribosomes that
  # elongate the CDS, not by scanning: handled after the scanning pass.
  continued <- classes %in% c("FS", "CTE")
  if (any(continued)) {
    inner <- rdg_gene_loads(starts[!continued], ends[!continued], classes[!continued],
                            start_codons[!continued], L, C, params, cds_start_codon,
                            gene_factor = gene_factor)
    eff <- if (is.null(efficiency)) rep(NA_real_, length(starts)) else efficiency
    default <- ifelse(classes == "FS", params$p_frameshift, params$p_readthrough)
    eff <- ifelse(is.na(eff), default, eff)
    initiations <- load <- numeric(length(starts))
    initiations[!continued] <- inner$initiations
    load[!continued] <- inner$load
    initiations[continued] <- inner$cds_initiations * eff[continued]
    load[continued] <- initiations[continued] * (ends[continued] - starts[continued] + 1L)
    # Ribosomes that shift frame leave the CDS frame at the shift site.
    fs <- classes == "FS"
    cds_load <- inner$cds_load - sum(initiations[fs] * pmax(0, C - starts[fs] + 1L))
    return(list(initiations = initiations, load = load,
                cds_initiations = inner$cds_initiations, cds_load = max(cds_load, 0)))
  }
  ids <- c("CDS", seq_along(starts))
  s <- c(L + 1L, starts)
  e <- c(C, ends)
  codon <- c(cds_start_codon, start_codons)
  p <- params$p_init[codon]
  p[is.na(p)] <- params$p_init_other
  # One initiation probability for the translon_classes features (position 1 is the CDS), and the per-gene factor.
  sel <- c(FALSE, classes %in% (if (is.null(params$translon_classes)) c("uORF", "uoORF", "NTE") else params$translon_classes))
  if (!is.null(params$p_init_translon)) p[sel] <- params$p_init_translon
  if (gene_factor != 1) p[sel] <- pmin(1, p[sel] * gene_factor)
  p <- pmin(1, p * params$p_init_scale)
  len <- e - s + 1L
  # Reinitiation after an ORF depends on its length; the CDS uses its own value.
  # An NTE ends at the CDS stop like the CDS itself.
  reinit <- params$reinit_max * exp(-(len / 3) / params$reinit_scale_codons)
  reinit[ids == "CDS" | c(FALSE, classes == "NTE")] <- params$reinit_after_cds
  order_by_start <- order(s)
  initiations <- numeric(length(s))
  scanning <- 1
  pending_pos <- numeric(0); pending_amount <- numeric(0)
  for (i in order_by_start) {
    joined <- pending_pos <= s[i]
    if (any(joined)) {
      scanning <- scanning + sum(pending_amount[joined])
      pending_pos <- pending_pos[!joined]; pending_amount <- pending_amount[!joined]
    }
    ini <- p[i] * scanning
    initiations[i] <- ini
    scanning <- scanning - ini
    # Terminating ribosomes that resume scanning rejoin after re-acquiring a
    # ternary complex downstream of the stop.
    pending_pos <- c(pending_pos, e[i] + params$reinit_distance)
    pending_amount <- c(pending_amount, ini * reinit[i])
  }
  # Extra repression of the CDS by every overlapping uoORF (uoorf_cds_factor).
  n_uo <- sum(classes == "uoORF")
  if (n_uo > 0L && !is.null(params$uoorf_cds_factor) && params$uoorf_cds_factor != 1)
    initiations[1] <- initiations[1] * params$uoorf_cds_factor^n_uo
  load <- initiations * len
  list(initiations = initiations[-1], load = load[-1],
       cds_initiations = initiations[1], cds_load = load[1])
}

#' Couple translated reads of every gene to its RDG with the ribosome flow model
#'
#' For genes built with \code{simGenome(architecture = )}, each library's
#' translated reads of a gene (its "cds" plus "uorf" counts) are split between
#' the CDS and its extra translons by their ribosome loads
#' (\code{\link{rdg_parameters}}), using the parameters of the library's
#' condition. Leader and trailer background are not changed. The per-library
#' share of every translon within its gene is stored in
#' \code{metadata(result)$rdg_feature_share} (translons x libraries, in the order
#' of \code{true_uORFs.rds}); \code{simNGScoverage()} uses it automatically.
#' Only RFP libraries are changed.
#' @param count_table_regions output of \code{simCountTablesRegions()}, with
#'  the regions "cds" and "uorf".
#' @param simGenome output of \code{simGenome(architecture = )}.
#' @param params \code{rdg_parameters()}.
#' @return the count table with new "cds" and "uorf" assays for RFP libraries
#' @export
simRDGflow <- function(count_table_regions, simGenome, params = rdg_parameters()) {
  stopifnot(is(count_table_regions, "SummarizedExperiment"), inherits(params, "rdg_parameters"))
  if (!all(c("cds", "uorf") %in% SummarizedExperiment::assayNames(count_table_regions)))
    stop("count_table_regions needs the regions 'cds' and 'uorf'")
  if (is.na(simGenome["features_table"]))
    stop("simGenome has no features table: build it with simGenome(architecture = )")
  features <- read_features_table(simGenome)[!is.na(read_group)]
  txdb <- ORFik::loadTxdb(simGenome[["txdb"]])
  leaders <- ORFik::loadRegion(txdb, "leaders")
  cds <- ORFik::loadRegion(txdb, "cds")
  genes <- rownames(count_table_regions)
  L <- setNames(rep(0L, length(genes)), genes)
  L[names(leaders)[names(leaders) %in% genes]] <- as.integer(sum(width(leaders[names(leaders) %in% genes])))
  C <- L + setNames(as.integer(sum(width(cds[genes]))), genes)
  cds_codon <- setNames(rep("ATG", length(genes)), genes)

  col <- SummarizedExperiment::colData(count_table_regions)
  libtype <- as.character(col$libtype)
  condition <- if ("condition" %in% names(col)) as.character(col$condition) else rep("", nrow(col))
  cds_counts <- SummarizedExperiment::assay(count_table_regions, "cds")
  orf_counts <- SummarizedExperiment::assay(count_table_regions, "uorf")
  share <- matrix(NA_real_, nrow(features), ncol(cds_counts),
                  dimnames = list(features$read_group, colnames(cds_counts)))
  by_gene <- split(seq_len(nrow(features)), features$transcript_id)
  # Only genes in this count table get shares; rows of other genes stay NA and
  # are dropped by simNGScoverage() together with their ranges.
  by_gene <- by_gene[names(by_gene) %in% genes]
  # Opt-in per-gene factor on the extra features' initiation (mean 1, log-normal), drawn once for all libraries.
  # No random numbers are drawn when translon_gene_sd is 0, so default runs keep their random stream.
  gene_factor <- setNames(rep(1, length(by_gene)), names(by_gene))
  sd_g <- if (is.null(params$translon_gene_sd)) 0 else params$translon_gene_sd
  if (sd_g > 0 && length(by_gene))
    gene_factor[] <- exp(stats::rnorm(length(by_gene), -sd_g^2 / 2, sd_g))
  for (cond in unique(condition)) {
    p <- rdg_condition_parameters(params, cond)
    libs <- which(condition == cond)
    gene_load <- lapply(names(by_gene), function(g) {
      rows <- by_gene[[g]]
      rdg_gene_loads(features$tx_start[rows], features$tx_end[rows], features$class[rows],
                     features$start_codon[rows], L[[g]], C[[g]], p, cds_codon[[g]],
                     efficiency = if ("efficiency" %in% names(features)) features$efficiency[rows],
                     gene_factor = gene_factor[[g]])
    })
    names(gene_load) <- names(by_gene)
    for (g in names(by_gene)) {
      rows <- by_gene[[g]]
      loads <- gene_load[[g]]
      total <- loads$cds_load + sum(loads$load)
      share[rows, libs] <- if (sum(loads$load) > 0) loads$load / sum(loads$load) else 1 / length(rows)
      for (j in libs[libtype[libs] == "RFP"]) {
        translated <- cds_counts[g, j] + orf_counts[g, j]
        if (translated <= 0 || total <= 0) next
        to_cds <- stats::rbinom(1L, translated, loads$cds_load / total)
        cds_counts[g, j] <- to_cds
        orf_counts[g, j] <- translated - to_cds
      }
    }
  }
  # Genes without translated extra features translate only their CDS.
  no_features <- setdiff(genes, names(by_gene))
  rfp <- which(libtype == "RFP")
  if (length(no_features) && length(rfp)) {
    cds_counts[no_features, rfp] <- cds_counts[no_features, rfp] + orf_counts[no_features, rfp]
    orf_counts[no_features, rfp] <- 0
  }
  SummarizedExperiment::assay(count_table_regions, "cds") <- cds_counts
  SummarizedExperiment::assay(count_table_regions, "uorf") <- orf_counts
  S4Vectors::metadata(count_table_regions)$rdg_feature_share <- share
  S4Vectors::metadata(count_table_regions)$rdg_parameters <- params
  if (sd_g > 0) S4Vectors::metadata(count_table_regions)$rdg_gene_factor <- gene_factor
  count_table_regions
}

# Which translated features (rows of the features table, read_group order) belong
# to the given transcripts; used in tests and by callers that subset genes.
features_in_genes <- function(simGenome, genes) {
  f <- read_features_table(simGenome)[!is.na(read_group)]
  f$transcript_id %in% genes
}
