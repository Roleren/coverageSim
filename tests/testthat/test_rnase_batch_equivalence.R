test_that("batched RNase extension equals the original for reordered groups", {
  original <- function(dt_range, rnase_bias, transcript_models) {
    reach <- rnase_kernel_reach(rnase_bias)
    if (reach == 0L) return(dt_range)
    groups <- split(dt_range, dt_range$genes, keep.by = TRUE)
    data.table::rbindlist(lapply(groups, function(group) {
      model <- transcript_models[[group$transcript_id[1]]]
      first_tx <- genomic_site_to_transcript(model, group$start[1])
      last_tx <- genomic_site_to_transcript(model, group$start[nrow(group)])
      prefix <- group[rep(1L, reach)]; suffix <- group[rep(nrow(group), reach)]
      prefix[, start := transcript_position_to_genomic(model, seq.int(first_tx - reach, first_tx - 1L))]
      suffix[, start := transcript_position_to_genomic(model, seq.int(last_tx + 1L, last_tx + reach))]
      prefix[, end := start]; suffix[, end := start]
      if ("region_position" %in% names(group)) {
        prefix[, region_position := seq.int(group$region_position[1L] - reach, group$region_position[1L] - 1L)]
        suffix[, region_position := seq.int(group$region_position[nrow(group)] + 1L, group$region_position[nrow(group)] + reach)]
      }
      data.table::rbindlist(list(prefix, group, suffix))
    }))
  }
  models <- list(a = list(exons = GenomicRanges::GRanges('chr1', IRanges::IRanges(c(1, 21), c(10, 30))), cumulative_start = c(1L, 11L), length = 20L, strand = '+'),
                 b = list(exons = GenomicRanges::GRanges('chr1', IRanges::IRanges(c(21, 1), c(30, 10))), cumulative_start = c(1L, 11L), length = 20L, strand = '-'))
  dt <- data.table::data.table(genes = c('b','a','b','a'), transcript_id = c('b','a','b','a'), start = c(25L, 7L, 8L, 22L), end = c(25L, 7L, 8L, 22L), region_position = c(1L, 1L, 2L, 2L), score = 1:4)
  keyed <- data.table::copy(dt); data.table::setkey(keyed, genes)
  indexed <- data.table::copy(dt); data.table::setindex(indexed, genes)
  factor_dt <- data.table::copy(dt); factor_dt[, genes := factor(genes)]
  for (x in list(dt, dt[, !'region_position'], dt[0], factor_dt, keyed, indexed)) {
    for (bias in list(list(RFP = 1), list(RFP = c(1,2,1)))) {
      expect_equal(append_rnase_to_simulated_rpf_table(data.table::copy(x), bias, models), original(data.table::copy(x), bias, models), tolerance = 0)
    }
  }
})
