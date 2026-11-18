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
#' @param conditions named list: per condition (as in the count table's
#'  \code{condition} column) a list of parameters that replace the defaults, or
#'  \code{p_init_scale}, a factor on all initiation probabilities.
#' @return a list of class "rdg_parameters"
#' @export
#' @examples
#' rdg_parameters(conditions = list(ISR = list(reinit_distance = 150)))
rdg_parameters <- function(p_init = c(ATG = 0.8, CTG = 0.12, GTG = 0.07, TTG = 0.05,
                                      ACG = 0.04),
                           p_init_other = 0.02, reinit_max = 0.6,
                           reinit_scale_codons = 30, reinit_distance = 50,
                           reinit_after_cds = 0.05,
                           conditions = list(ISR = list(reinit_distance = 150,
                                                        p_init_scale = 0.8))) {
  params <- list(p_init = p_init, p_init_other = p_init_other, reinit_max = reinit_max,
                 reinit_scale_codons = reinit_scale_codons,
                 reinit_distance = reinit_distance, reinit_after_cds = reinit_after_cds,
                 p_init_scale = 1, conditions = conditions)
  probs <- c(p_init, p_init_other, reinit_max, reinit_after_cds)
  if (any(!is.finite(probs) | probs < 0 | probs > 1))
    stop("Probabilities in rdg_parameters must be within [0, 1]")
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
rdg_gene_loads <- function(starts, ends, classes, start_codons, L, C, params,
                           cds_start_codon = "ATG") {
  ids <- c("CDS", seq_along(starts))
  s <- c(L + 1L, starts)
  e <- c(C, ends)
  codon <- c(cds_start_codon, start_codons)
  p <- params$p_init[codon]
  p[is.na(p)] <- params$p_init_other
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
  for (cond in unique(condition)) {
    p <- rdg_condition_parameters(params, cond)
    libs <- which(condition == cond)
    gene_load <- lapply(names(by_gene), function(g) {
      rows <- by_gene[[g]]
      rdg_gene_loads(features$tx_start[rows], features$tx_end[rows], features$class[rows],
                     features$start_codon[rows], L[[g]], C[[g]], p, cds_codon[[g]])
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
  count_table_regions
}
