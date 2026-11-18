test_that("recycled presets retain the per-gene table and random state", {
  for (presets in list("all", c("uORF", "dORF", "uORF"))) {
    ids <- paste0("tx", 1:7)
    set.seed(711)
    expected <- data.table::rbindlist(lapply(rep_len(presets, 7), architecture_preset), idcol = "gene", fill = TRUE)
    expected[, transcript_id := ids[gene]]
    expected[, feature_id := paste0(transcript_id, "_", class, "_", seq_len(.N)), by = gene]
    state <- .Random.seed
    set.seed(711)
    expect_equal(architecture_per_gene(presets, 7, ids), expected, tolerance = 0)
    expect_identical(.Random.seed, state)
  }
})
