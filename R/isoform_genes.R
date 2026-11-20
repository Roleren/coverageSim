# Multi-isoform genes (backclaude E17-S7, 2026-11-20).
#
# Real genes make several transcript isoforms. Ribo-seq and RNA-seq reads of a gene are then a MIXTURE over its
# isoforms (weights = isoform RNA abundance, times a translation efficiency for ribosomes), while the analysis
# usually works on one reference transcript (MANE). For 3,708 of 17,935 human genes the MANE isoform has an RNA
# weight below 0.5 (1,065 below 0.1): their reads look "unexplained" against a flat MANE expectation (E10/A, E14/G,
# E16/domiso). This builds such genes constructively, in the codon coordinates of the reference isoform:
#   - n_iso, reference weight w_ref and the isoform-specific CDS fraction are resampled jointly from the real
#     distribution (inst/extdata/isoform_human_summary.csv, 17,935 genes of E11/C iso_gene.csv);
#   - the reference CDS is cut into shared segments and isoform-specific segments; each specific segment is absent
#     (exon skipping / alternative start or end) from a random non-empty set of the non-reference isoforms;
#   - ribosome reads follow  dwell(codon) * sum over isoforms containing the codon of (RNA weight * TE), the RNA-seq
#     track the same without dwell and TE, so the explainer's RNA-mixture expectation can be tested against truth.
# Nothing here changes simGenome() / simNGScoverage(): the object is a separate, opt-in layer on the reference
# transcript, and isoform_mixture_alpha() hands the mixture to sim_sequence_bias() as a gene alpha matrix.

human_isoform_table <- function() {
  f <- system.file("extdata", "isoform_human_summary.csv", package = "coverageSim")
  if (!nzchar(f)) f <- file.path("inst", "extdata", "isoform_human_summary.csv")
  data.table::fread(f)
}

#' Simulate genes with several isoforms
#'
#' @param n number of genes.
#' @param n_iso \code{"human"} (default): number of isoforms, reference weight and specific-CDS structure are
#'  resampled jointly from the 17,935 real human genes. A whole number (or vector) fixes the isoform count; the
#'  other quantities are then drawn from real genes with that count (counts above 10 are treated as 10).
#' @param ref_minor_frac \code{NULL} (default) keeps the real mix (20.7\% of all genes, 28.5\% of multi-isoform genes,
#'  have a reference weight below 0.5). A number in [0, 1] sets the fraction of multi-isoform genes whose reference
#'  (MANE-like) isoform is a minor isoform (weight < 0.5); the reference weight is then redrawn from the real genes
#'  of that class.
#' @param cds_codons reference CDS length in codons: one number, a vector, or \code{"human"} (drawn with
#'  \code{human_tx_lengths()}). Default 300.
#' @param dir_alpha Dirichlet concentration splitting the non-reference weight over the other isoforms
#'  (small = one dominant alternative).
#' @param te_sd sd of a log-normal per-isoform translation efficiency (default 0 = ribosomes follow RNA weights).
#'  Above 0 the ribosome mixture differs from the RNA mixture, as in real data (slope of log O/E on the RNA ratio 0.07-0.37).
#' @param rna_weight_sd sd of a log-normal perturbation of the isoform weights seen by the RNA-seq track (default 0);
#'  models an RNA library from another tissue than the ribosome library.
#' @param dwell_sd sd of the log-normal per-codon dwell (shared by all isoforms of a gene), default 0.5.
#' @param seed integer or NULL; when given, the session RNG state is restored afterwards.
#' @return an object of class \code{isoformGenes}: a list with \code{genes} (one row per gene: n_iso, ref isoform,
#'  w_ref, cds_codons, w_dom), \code{iso} (gene, iso, is_ref, w_rna, te, w_ribo, w_rna_seen), \code{seg}
#'  (gene, seg, start, end codon in reference coordinates, specific) and \code{member} (gene, seg, iso, present);
#'  \code{dwell} is a list of per-codon dwell vectors.
#' @export
#' @examples
#' x <- simIsoformGenes(20, ref_minor_frac = 0.5, seed = 1)
#' x$genes
simIsoformGenes <- function(n, n_iso = "human", ref_minor_frac = NULL, cds_codons = 300, dir_alpha = 0.5,
                            te_sd = 0, rna_weight_sd = 0, dwell_sd = 0.5, seed = NULL) {
  stopifnot(n >= 1, dir_alpha > 0, te_sd >= 0, rna_weight_sd >= 0, dwell_sd >= 0)
  if (!is.null(ref_minor_frac)) stopifnot(ref_minor_frac >= 0, ref_minor_frac <= 1)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, globalenv()))
    set.seed(seed)
  }
  tab <- human_isoform_table()
  tab[, n_iso := pmin(n_iso, 10L)]
  multi <- which(tab$n_iso > 1L)
  minor_rows <- multi[tab$w_mane[multi] < 0.5]; major_rows <- multi[tab$w_mane[multi] >= 0.5]
  if (identical(n_iso, "human")) {
    if (is.null(ref_minor_frac)) rows <- sample.int(nrow(tab), n, TRUE)
    else {
      single <- which(tab$n_iso == 1L)
      is_multi <- stats::runif(n) < length(multi) / nrow(tab)
      minor <- stats::runif(n) < ref_minor_frac
      rows <- ifelse(!is_multi, sample(single, n, TRUE),
                     ifelse(minor, sample(minor_rows, n, TRUE), sample(major_rows, n, TRUE)))
    }
  } else {
    k <- pmin(rep_len(as.integer(n_iso), n), 10L); stopifnot(all(k >= 1L))
    rows <- vapply(k, function(kk) {
      pool <- if (kk == 1L) which(tab$n_iso == 1L) else which(tab$n_iso == kk)
      if (kk > 1L && !is.null(ref_minor_frac)) {
        pool <- pool[(tab$w_mane[pool] < 0.5) == (stats::runif(1) < ref_minor_frac)]
        if (!length(pool)) pool <- which(tab$n_iso == kk)
      }
      if (!length(pool)) pool <- multi
      pool[sample.int(length(pool), 1L)]
    }, integer(1))
  }
  r <- tab[rows]
  k <- ifelse(is.na(r$n_iso), 1L, r$n_iso)
  if (!identical(n_iso, "human")) k <- pmin(rep_len(as.integer(n_iso), n), 10L)
  L <- if (identical(cds_codons, "human")) pmax(64L, as.integer(human_tx_lengths(n)$cds_length / 3)) else
    as.integer(rep_len(cds_codons, n))
  stopifnot(all(L >= 10L))
  gid <- sprintf("G%05d", seq_len(n))
  genes <- list(); iso <- list(); seg <- list(); mem <- list(); dwell <- vector("list", n)
  for (g in seq_len(n)) {
    kk <- k[g]; wref <- if (kk == 1L) 1 else min(max(r$w_mane[g], 0), 1)
    ref <- sample.int(kk, 1L)       # position of the reference isoform among the k (irrelevant to truth, kept for the table)
    w <- rep(0, kk)
    if (kk > 1L) {
      d <- stats::rgamma(kk - 1L, dir_alpha); d <- d / sum(d)
      w[-ref] <- (1 - wref) * d
    }
    w[ref] <- wref
    w <- if (sum(w) > 0) w / sum(w) else rep(1 / kk, kk)
    te <- if (te_sd > 0) exp(stats::rnorm(kk, 0, te_sd) - te_sd^2 / 2) else rep(1, kk)
    wr <- w * te; wr <- wr / sum(wr)
    ws <- if (rna_weight_sd > 0) { v <- w * exp(stats::rnorm(kk, 0, rna_weight_sd)); v / sum(v) } else w
    iso[[g]] <- data.table::data.table(gene = gid[g], iso = seq_len(kk), is_ref = seq_len(kk) == ref,
                                       w_rna = w, te = te, w_ribo = wr, w_rna_seen = ws)
    # reference CDS: alternating shared / specific segments
    ns <- if (kk == 1L) 0L else as.integer(r$n_specific_seg[g])
    fs <- if (kk == 1L) 0 else r$frac_cds_specific[g]
    ns <- min(ns, max(L[g] %/% 2L, 0L), 40L)
    if (ns == 0L || fs <= 0) { ns <- 0L; lens_spec <- integer(0) } else {
      tot <- min(max(round(fs * L[g]), ns), L[g] - 1L)
      lens_spec <- 1L + stats::rmultinom(1, tot - ns, rep(1, ns))[, 1]
    }
    n_shared <- L[g] - sum(lens_spec)
    lens_shared <- if (ns == 0L) n_shared else stats::rmultinom(1, n_shared, rep(1, ns + 1L))[, 1]
    pieces <- integer(0); spec <- logical(0)
    for (j in seq_len(ns + 1L)) {
      pieces <- c(pieces, lens_shared[j]); spec <- c(spec, FALSE)
      if (j <= ns) { pieces <- c(pieces, lens_spec[j]); spec <- c(spec, TRUE) }
    }
    keep <- pieces > 0L; pieces <- pieces[keep]; spec <- spec[keep]
    end <- cumsum(pieces); start <- end - pieces + 1L
    seg[[g]] <- data.table::data.table(gene = gid[g], seg = seq_along(pieces), start = start, end = end, specific = spec)
    others <- setdiff(seq_len(kk), ref)
    m <- do.call(rbind, lapply(seq_along(pieces), function(s) {
      pres <- rep(TRUE, kk)
      if (spec[s] && length(others)) {
        absent <- others[stats::runif(length(others)) < 0.5]
        if (!length(absent)) absent <- others[sample.int(length(others), 1L)]
        pres[absent] <- FALSE
      }
      data.table::data.table(gene = gid[g], seg = s, iso = seq_len(kk), present = pres)
    }))
    mem[[g]] <- m
    dwell[[g]] <- if (dwell_sd > 0) exp(stats::rnorm(L[g], 0, dwell_sd) - dwell_sd^2 / 2) else rep(1, L[g])
    genes[[g]] <- data.table::data.table(gene = gid[g], n_iso = kk, ref_iso = ref, w_ref = wref, cds_codons = L[g],
                                         w_dom = max(w[-ref], 0), frac_specific = sum(lens_spec) / L[g])
  }
  names(dwell) <- gid
  structure(list(genes = data.table::rbindlist(genes), iso = data.table::rbindlist(iso),
                 seg = data.table::rbindlist(seg), member = data.table::rbindlist(mem), dwell = dwell,
                 options = list(n_iso = n_iso, ref_minor_frac = ref_minor_frac, te_sd = te_sd,
                                rna_weight_sd = rna_weight_sd, dwell_sd = dwell_sd)),
            class = "isoformGenes")
}

#' @export
print.isoformGenes <- function(x, ...) {
  g <- x$genes
  cat("isoformGenes:", nrow(g), "genes,", sum(g$n_iso), "isoforms; reference weight < 0.5 in",
      sum(g$w_ref < 0.5), "genes, < 0.1 in", sum(g$w_ref < 0.1), "\n")
  invisible(x)
}

# usage per reference codon: the summed weight of the isoforms that contain it (column = "w_ribo", "w_rna" or "w_rna_seen")
isoform_codon_usage <- function(x, col = "w_ribo") {
  cov <- merge(x$member[present == TRUE], x$iso[, c("gene", "iso", col), with = FALSE], by = c("gene", "iso"))
  su <- cov[, .(u = sum(get(col))), by = .(gene, seg)]
  sg <- merge(x$seg, su, by = c("gene", "seg"), all.x = TRUE)
  sg[is.na(u), u := 0]
  data.table::setorder(sg, gene, seg)
  split(rep(sg$u, sg$end - sg$start + 1L), rep(sg$gene, sg$end - sg$start + 1L))[x$genes$gene]
}

#' Expected codon profiles of simulated isoform genes
#'
#' Per gene, in the codon coordinates of the reference (MANE-like) isoform, each profile scaled to a total of 1:
#' \code{truth_ribo} = dwell x usage of the codon by the ribosome-weighted isoform mixture (the simulator's own mixture
#' expectation E); \code{flat} = dwell only, i.e. the expectation if the reference isoform were the only one;
#' \code{truth_rna} = the RNA-seq mixture (uniform coverage, with the RNA weights as seen by the RNA library).
#' @param x an \code{isoformGenes} object.
#' @return a list with the three profile lists named by gene.
#' @export
isoform_expected <- function(x) {
  norm1 <- function(l) lapply(l, function(v) v / sum(v))
  ur <- isoform_codon_usage(x, "w_ribo"); us <- isoform_codon_usage(x, "w_rna_seen")
  list(truth_ribo = norm1(Map(function(d, u) d * u, x$dwell, ur)),
       flat = norm1(x$dwell),
       truth_rna = norm1(us))
}

#' Mixture alpha for sim_sequence_bias()
#'
#' The ribosome mixture as a gene alpha matrix (one codon-level weight vector per gene, mean preserved), so the
#' isoform mixture can be pushed through the normal codon-level bias chain of \code{simNGScoverage()}.
#' @param x an \code{isoformGenes} object.
#' @return a named list of numeric vectors, one per gene.
#' @export
isoform_mixture_alpha <- function(x) {
  lapply(isoform_expected(x)$truth_ribo, function(v) pmax(v / mean(v), 1e-6))
}

#' Simulate ribosome and RNA reads of isoform genes
#'
#' Poisson reads per codon from the isoform mixture. Gene totals are log-normal around \code{ribo_depth} (reads over the
#' reference CDS in the truth mixture).
#' @param x an \code{isoformGenes} object.
#' @param ribo_depth,rna_depth median reads per gene (ribo: per replicate).
#' @param depth_sd sd of the log gene totals (ribo and RNA share the gene factor).
#' @param replicates number of ribosome replicates (columns O1..On); RNA is one merged library.
#' @param seed integer or NULL (session RNG restored afterwards).
#' @return a data.table, one row per reference codon: gene, codon, w_ref, O1.., R (RNA reads), E_truth, E_flat, E_rna
#'  (expected counts of replicate 1 / RNA), plus \code{ref_transcript}: the reference isoform's own output (E_flat
#'  scaled to the observed gene total).
#' @export
simIsoformReads <- function(x, ribo_depth = 3000, rna_depth = 3000, depth_sd = 0.5, replicates = 2, seed = NULL) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, globalenv()))
    set.seed(seed)
  }
  ex <- isoform_expected(x); n <- nrow(x$genes)
  gf <- exp(stats::rnorm(n, 0, depth_sd))
  out <- lapply(seq_len(n), function(g) {
    L <- x$genes$cds_codons[g]; id <- x$genes$gene[g]
    et <- ex$truth_ribo[[g]] * ribo_depth * gf[g]
    d <- data.table::data.table(gene = id, codon = seq_len(L), w_ref = x$genes$w_ref[g])
    for (j in seq_len(replicates)) data.table::set(d, j = paste0("O", j), value = stats::rpois(L, et))
    er <- ex$truth_rna[[g]] * rna_depth * gf[g]
    d[, `:=`(R = stats::rpois(L, er), E_truth = et, E_flat = ex$flat[[g]] * ribo_depth * gf[g], E_rna = er)]
    d
  })
  res <- data.table::rbindlist(out)
  res[, ref_transcript := E_flat / sum(E_flat) * sum(Reduce(`+`, .SD)) / replicates, by = gene,
      .SDcols = paste0("O", seq_len(replicates))]
  res[]
}
