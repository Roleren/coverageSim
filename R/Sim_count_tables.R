#' Create count tables of simulated data
#'
#' Follows the rules defined in DESeq paper, using negative binomial distributions.
#' The true values statistics are found in the meta columns.
#'
#' Info about params:\cr
#' - dispMeanRel: see DESeq2, section: Estimation of dispersions
#' @inheritParams DESeq2::makeExampleDESeqDataSet
#' @param n either transcripts as GRangesList or GRanges, or an integer number of transcripts (rows in count table).
#' Default 500. If integer, will create random ranges assigned to each transcript.
#' @param replicates = 2
#' @param conditions = c("WT", "Mutant")
#' @param libtypes = c("RFP", "RNA", "CAGE", "PAS")
#' @param betaSD integer, default: 1 (log2 value). The standard deviation of dispersion between the conditions.
#' The higher it is, the lower the correlation between conditions. Sane values: Between 1-3
#' @param betaLibSD named numeric vector, default: \code{c("RFP" = 1, "RNA" = 1, "CAGE" = 1, "PAS" = 1)}
#' Standard deviation of dispersion, relative to the anchor distribution (log2 value).
#' The higher these values the lower correlation between library types.
#' @param interceptMean numeric, default 4 (log2 value). The mean value for the anchoring distribution.
#' The higher it is, the more counts all libraries will get. Sane values: Between 1-12
#' @param interceptSD numeric, default 2 (log2 value). The standard deviation for the anchoring distribution.
#' The higher it is, the more variance counts of the genes can have. Sane values: Between 1-3
#' @param m number of total samples, \code{m = length(libtypes) * length(conditions) * replicates}
#' @references https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8
#' @return a ranged summarized experiment object
#' @export
#' @examples
#' # Sim 500 genes with: 4 library types, 2 conditions, 2 replicates each
#' simCountTables()
#' #' # Sim 500 genes with: Only Ribo-seq, 2 conditions, 2 replicates each
#' simCountTables(libtypes = "RFP")
simCountTables <- function (n = 500, libtypes = c("RFP", "RNA", "CAGE", "PAS"),
                            conditions = c("WT", "Mutant"), replicates = 2,
                            m = length(libtypes) * length(conditions) * replicates,
                            betaSD = 1, betaLibSD = c("RFP" = 1, "RNA" = 1, "CAGE" = 1, "PAS" = 1),
                            interceptMean = 4, interceptSD = 2,
                            dispMeanRel = function(x) 4/x + 0.1, sizeFactors = rep(1, m),
                            plot_PCA = TRUE, print_statistics = TRUE)
{
  message("Simulating count tables --")
  stopifnot(m > 0)
  stopifnot(!is.null(names(betaLibSD)))
  betaLibSD <- betaLibSD[libtypes]
  if (!is(n, "numeric")) {
    rowRanges <- n
    n <- length(n)
  } else {
    rowRanges <- GRanges("1", IRanges(start = (1:n - 1) * 100 +
                                        1, width = 100))
    names(rowRanges) <- paste0("gene", 1:n)
  }


  # The anchor normal sample (Lib types are based on this)
  beta_anchor <- rnorm(n, interceptMean, interceptSD)
  # The GLM dispersion between library types
  beta_per_libtype <- matrix(beta_anchor + unlist(lapply(betaLibSD, function(x) rnorm(n, 0, x))), ncol = length(betaLibSD))
  colnames(beta_per_libtype) <- libtypes

  if (length(conditions) > 1) { # The GLM dispersion between conditions
    beta_between_cond <- rnorm(n, 0, betaSD)
  } else { # Add 0 betaSDs for metacolumns
    beta_between_cond <- NULL
  }
  dispersion <- dispMeanRel(2^(beta_per_libtype))
  beta <- cbind(beta_per_libtype, beta_between_cond)

  colData <- DataFrame(libtype = factor(rep(libtypes, each = replicates*length(conditions)), levels = unique(libtypes)),
                       condition = factor(rep(conditions, each = replicates), levels = unique(conditions)),
                       replicate = as.character(seq(replicates)))
  # The design matrix
  x <- if (length(conditions) > 1 | length(libtypes) > 1) {
    if (length(conditions) > 1 & length(libtypes) > 1){
      stats::model.matrix.default(~ 0 + libtype + condition, data = colData) #+ colData$libtype:colData$condition
    } else if (length(libtypes) > 1){
      stats::model.matrix.default(~ 0 + colData$libtype)
    } else if (length(conditions) > 1) {
      stats::model.matrix.default(~ 1 + colData$condition)
    }
  }
  else {
    if (m == 1) {
      rep(1, m)
    } else cbind(rep(1, m), rep(0, m))
  }

  mu <- t(2^(x %*% t(beta)) * sizeFactors)
  # Size is here a dispersion shape parameter (non integer)
  countData <- matrix(rnbinom(m * n, mu = mu, size = 1/dispersion),
                      ncol = m)

  mode(countData) <- "integer"
  colnames(countData) <- paste(colData$libtype, colData$condition, colData$replicate, sep = "_")

  object <- SummarizedExperiment(assays = countData, colData = colData, rowRanges = rowRanges)
  if (m == 1) beta <- cbind(beta, rep(0, nrow(beta)))
  trueVals <- DataFrame(trueIntercept = beta[, 1], trueBeta = beta[,2], trueDisp = dispersion)
  colnames(trueVals)[2 + seq_along(betaLibSD)] <- paste0("trueDisp_", names(betaLibSD))
  mcols(trueVals) <- DataFrame(type = rep("input", ncol(trueVals)),
                               description = c("simulated intercept values", "simulated beta values",
                                               paste0("simulated dispersion values: ", names(betaLibSD))))
  mcols(object) <- cbind(mcols(object), trueVals)
  message("Count table statistics --")

  if (max(colSums(assay(object))) > 1e4) {
    message("Reads per library (in millions):")
    print(round(colSums(assay(object)) / 1e6, 2))
  } else {
    message("Reads per library (in thousands):")
    print(round(colSums(assay(object)) / 1e3, 2))
  }

  # Statistics
  print <- (length(conditions) > 1 | length(libtypes) > 1) & print_statistics
  if (print) {
    message("PCA plot:")
    if (length(conditions) > 1 & length(libtypes) > 1) {
      d_mat <- DESeq2::DESeqDataSet(object, design = ~ libtype + condition)
    } else if (length(conditions) > 1) {
      d_mat <- DESeq2::DESeqDataSet(object, design = ~ condition)
    } else if (length(libtypes) > 1) d_mat <- DESeq2::DESeqDataSet(object, design = ~ libtype)
    if (plot_PCA) {
      rd <- DESeq2::rlog(d_mat)
      plot(DESeq2::plotPCA(rd, c("libtype", "condition")))
      plot(DESeq2::plotPCA(rd, "libtype"))
    }
    if (length(libtypes) > 1) {
      message("Correlation between libtypes:")
      print(cor(assay(DESeq2::collapseReplicates(d_mat, groupby = d_mat$libtype))))
    }
    if (length(conditions) > 1) {
      message("Correlation between condition:")
      print(cor(assay(DESeq2::collapseReplicates(d_mat, groupby = d_mat$condition))))
    }
  }
  return(object)
}

#' Info about params:\cr
#' - dispMeanRel: see DESeq2, section: Estimation of dispersions
#' @param count_table a
#' @param count_table a
#' @param count_table a
#' @param count_table a
#' @param count_table a
#' @return a rangedSummarizedExperiment
#' @export
#' @examples
#' # Sim 500 genes with: 4 library types, 2 conditions, 2 replicates each
#' # With default region spread
#' simCountTablesRegions()
#' # Sim 500 genes with: Only Ribo-seq, 2 conditions, 2 replicates each
#' # With default region spread
#' simCountTablesRegions(simCountTables(libtypes = "RFP"))
#' # Ribo-seq with all uORFs set to 0 TRUE coverage
#'
simCountTablesRegions <- function(count_table = simCountTables(),
                                  regionsToSample = c("leader", "cds", "trailer", "uorf"),
                                  region_proportion =
                                    list(leader =  list(RFP = 0.05, RNA = 0.3, CAGE = 1, PAS = 0),
                                         cds =     list(RFP = 0.75, RNA = 0.4, CAGE = 0, PAS = 0),
                                         trailer = list(RFP = 0.1, RNA = 0.3, CAGE = 0, PAS = 1),
                                         uorf =    list(RFP = 0.1, RNA = 0,   CAGE = 0, PAS = 0)),
                                  sampling = c(RFP = "MN", RNA = "MN", CAGE = "MN", PAS = "MN")) {
  stopifnot(is(count_table, "SummarizedExperiment"))
  message("Simulating region count tables -")
  stopifnot(all(sampling %in% c("MN", "DMN")))
  prop_regions <- names(region_proportion)
  all_allowed_regions <- c("leader", "cds", "trailer", "uorf")
  stopifnot(all(regionsToSample %in% all_allowed_regions))
  stopifnot(all(prop_regions %in% all_allowed_regions))
  prop_mat <- t(matrix(unlist(region_proportion), ncol = length(region_proportion),
                       dimnames = list(names(region_proportion[[1]]), names(region_proportion))))
  stopifnot(all(colSums(prop_mat) == 1))

  assay <- assay(count_table)
  colnames <- colnames(assay)
  rownames <- rownames(assay)
  libtypes <- as.character(colData(count_table)$libtype)
  stopifnot(all(libtypes %in% names(sampling)))
  mat <- matrix(nrow = nrow(assay), ncol = 0)
  for (s in seq_along(colnames)) {
    libtype <- libtypes[s]
    if (sampling[libtype] == "MN") {
      dt <- t(as.data.table(lapply(assay[,s], function(x)
        rmultinom(1, x, prob = prop_mat[,libtype]))))
    } else {
      dt <- t(as.data.table(lapply(assay[,s], function(x)
        extraDistr::rdirmnom(n = 1, size = x,
                             alpha = prop_mat[,libtype]))))
    }
    colnames(dt) <- names(region_proportion)
    mat <-cbind(mat, dt)

  }
  split_mat <- lapply(regionsToSample, function(region) {
    mat_region <- as.matrix(mat[, colnames(mat) == region])
    colnames(mat_region) <- colnames
    rownames(mat_region) <- rownames
    mat_region
  })
  assays(count_table) <- c(assays(count_table), split_mat)
  assayNames(count_table) <- c("gene", regionsToSample)
  message("Done")
  return(count_table)
}
