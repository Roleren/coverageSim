raw_test_transcripts <- function(n = 40, width = 900, seed = 11) {
  set.seed(seed)
  sequences <- vapply(seq_len(n), function(i) {
    paste(sample(c("A", "C", "G", "T"), width, replace = TRUE), collapse = "")
  }, character(1))
  Biostrings::DNAStringSet(stats::setNames(sequences, paste0("tx", seq_len(n))))
}

raw_read_fastq <- function(path) {
  lines <- readLines(path)
  index <- seq(1L, length(lines), by = 4L)
  data.table::data.table(read_id = sub("^@", "", lines[index]), sequence = lines[index + 1L],
                         quality = lines[index + 3L])
}

# The trimmer a pipeline would run with known settings: cut at the adapter (or
# the poly(A) tail), then remove fixed numbers of bases from both ends.
raw_trim_known <- function(sequence, trim, seed_length = 10L) {
  seed <- substr(trim$adapter, 1L, seed_length)
  position <- regexpr(seed, sequence, fixed = TRUE)
  found <- position > 0L
  kept <- ifelse(found, substr(sequence, 1L, position - 1L), sequence)
  kept <- substr(kept, trim$trim_front + 1L, nchar(kept) - trim$trim_tail)
  list(insert = kept, found = found)
}

raw_simulate_test <- function(protocol, ..., n_reads = 3000, read_length = 100) {
  prefix <- tempfile("rawreads_")
  out <- simulate_raw_reads(raw_test_transcripts(), n_reads = n_reads, protocol = protocol,
                            out_prefix = prefix, read_length = read_length, seed = 5, ...)
  list(out = out, truth = data.table::fread(out$files[["truth"]], na.strings = "NA",
                                            colClasses = list(character = c(
                                              "umi5", "barcode5", "add5", "add3", "umi3",
                                              "barcode3", "seqname"))),
       reads = raw_read_fastq(out$files[grep("^fastq", names(out$files))][1]))
}

test_that("every preset builds and lists its trim settings", {
  presets <- raw_read_protocols()
  expect_true(all(c("ingolia2012", "mcglincy2017", "nextflex_small_rna", "truseq_small_rna",
                    "illumina_truseq", "smarter", "ottr", "polya_cats", "ribolace",
                    "custom") %in% presets$preset))
  expect_identical(presets[preset == "ingolia2012", adapter], "CTGTAGGCACCATCAATAGA")
  expect_identical(presets[preset == "truseq_small_rna", substr(adapter, 1, 12)], "TGGAATTCTCGG")
  expect_identical(presets[preset == "illumina_truseq", substr(adapter, 1, 12)], "AGATCGGAAGAG")
  expect_identical(presets[preset == "nextflex_small_rna", c(trim_front, trim_tail)], c(4L, 4L))
  expect_identical(presets[preset == "mcglincy2017", c(trim_front, trim_tail)], c(2L, 10L))
  expect_identical(presets[preset == "ottr", c(trim_front, trim_tail)], c(7L, 1L))
  expect_identical(presets[preset == "smarter", c(trim_front, adapter)], c("3", "AAAAAAAAAA"))
  expect_error(raw_read_protocol("nope"), "preset must be one of")
  expect_error(raw_read_protocol("custom", barcode5 = c("ACG", "AC")), "one common length")
})

test_that("known trimmer settings recover the inserts of every ligation preset", {
  transcripts <- raw_test_transcripts()
  for (preset in c("ingolia2012", "mcglincy2017", "nextflex_small_rna", "truseq_small_rna",
                   "illumina_truseq", "ottr", "ribolace")) {
    # add5 off: a variable untemplated base cannot be cut by a fixed setting
    protocol <- raw_read_protocol(preset, add5 = NULL)
    sim <- raw_simulate_test(protocol, quality = "perfect", dimer_rate = 0.1)
    truth <- sim$truth
    expect_identical(sim$reads$read_id, truth$read_id)
    expect_true(all(nchar(sim$reads$sequence) == 100L))
    trimmed <- raw_trim_known(sim$reads$sequence, sim$out$trim)
    expect_true(all(trimmed$found), info = preset)
    expected <- as.character(Biostrings::subseq(
      transcripts[truth$seqname[truth$source == "main"]],
      truth$insert_start[truth$source == "main"], truth$insert_end[truth$source == "main"]))
    expect_identical(trimmed$insert[truth$source == "main"], unname(expected), info = preset)
    expect_true(all(trimmed$insert[truth$source == "adapter_dimer"] == ""), info = preset)
    expect_equal(mean(truth$source == "adapter_dimer"), 0.1, tolerance = 0.01)
  }
})

test_that("truth columns describe the read", {
  protocol <- raw_read_protocol("custom", umi5 = 3, barcode5 = c("ACGTA", "TTGCA"),
                                barcode5_prob = c(3, 1), umi3 = 2, barcode3 = "GGC",
                                add5 = list(base = "G", length_prob = c("0" = 1, "1" = 1, "2" = 1)))
  sim <- raw_simulate_test(protocol, quality = "perfect")
  truth <- sim$truth
  read <- sim$reads$sequence
  five <- paste0(truth$umi5, truth$barcode5, truth$add5)
  expect_identical(substr(read, 1L, nchar(five)), five)
  expect_identical(truth$insert_start_read, nchar(five) + 1L)
  expect_identical(substr(read, truth$insert_end_read + 1L, truth$adapter_start_read - 1L),
                   paste0(truth$umi3, truth$barcode3))
  expect_identical(substr(read, truth$adapter_start_read, truth$adapter_start_read + 11L),
                   rep("AGATCGGAAGAG", nrow(truth)))
  expect_setequal(truth$barcode5, c("ACGTA", "TTGCA"))
  expect_gt(mean(truth$barcode5 == "ACGTA"), 0.7)
  expect_setequal(nchar(truth$add5), 0:2)
  expect_true(all(nchar(truth$umi5) == 3L))
  expect_gt(length(unique(truth$umi5)), 50L)
})

test_that("poly(A) and template switching presets are recovered up to genomic A's", {
  transcripts <- raw_test_transcripts()
  for (preset in c("smarter", "polya_cats")) {
    sim <- raw_simulate_test(preset, quality = "perfect", read_length = 100)
    truth <- sim$truth
    expect_true(all(truth$add5 == "GGG"))
    expect_true(all(truth$polya_len >= 10L))
    insert <- as.character(Biostrings::subseq(transcripts[truth$seqname], truth$insert_start,
                                              truth$insert_end))
    trimmed <- raw_trim_known(sim$reads$sequence, sim$out$trim)
    # a tail hides A's at the insert end, and 10 A's inside an insert end it early
    no_run <- !grepl("AAAAAAAAAA", insert, fixed = TRUE)
    expect_identical(trimmed$insert[no_run], sub("A+$", "", unname(insert))[no_run], info = preset)
  }
})

test_that("reads are cut to the read length and short reads lose the adapter", {
  sim <- raw_simulate_test("nextflex_small_rna", quality = "perfect", read_length = 30)
  expect_true(all(nchar(sim$reads$sequence) == 30L))
  expect_identical(sim$truth$insert_complete, sim$truth$insert_end_read <= 30L)
  expect_lt(sim$out$summary$insert_complete, 0.5)
  long <- raw_simulate_test("nextflex_small_rna", quality = "perfect", read_length = 150)
  expect_equal(long$out$summary$adapter_in_read, 1)
  expect_true(all(grepl("A{20}$", long$reads$sequence)))
})

test_that("qualities and errors follow the profile", {
  perfect <- raw_simulate_test("illumina_truseq", quality = "perfect")
  expect_true(all(perfect$truth$n_errors == 0L))
  expect_true(all(perfect$reads$quality == strrep("I", 100)))
  nova <- raw_simulate_test("illumina_truseq", quality = "novaseq", read_length = 150)
  expect_true(all(grepl("G{10}$", nova$reads$sequence[nova$truth$n_errors == 0L])))
  expect_true(all(unlist(strsplit(nova$reads$quality[1:200], "")) %in% c("#", "-", "8", "F")))
  expect_gt(nova$out$summary$error_rate, 0.0005)
  expect_lt(nova$out$summary$error_rate, 0.03)
  trimmed <- raw_trim_known(nova$reads$sequence, nova$out$trim)
  exact <- nova$truth$n_errors == 0L
  expect_gt(mean(exact), 0.4)
  expect_true(all(nchar(trimmed$insert[exact]) == nova$truth$insert_len[exact]))
  hiseq <- raw_simulate_test("illumina_truseq", quality = "hiseq", error_rate = 0.02,
                             read_length = 150)
  expect_gt(length(unique(unlist(strsplit(hiseq$reads$quality[1:200], "")))), 10L)
  expect_gt(hiseq$out$summary$error_rate, 0.02)
  expect_true(all(grepl("A{10}$", hiseq$reads$sequence[hiseq$truth$n_errors == 0L])))
  expect_true(all(hiseq$truth$n_errors_insert <= hiseq$truth$n_errors))
})

test_that("contaminants and a second organism keep their shares and labels", {
  rrna <- raw_test_transcripts(2, 1800, seed = 3)
  other <- raw_test_transcripts(10, 900, seed = 4)
  names(other) <- paste0("yeast", 1:10)
  sim <- raw_simulate_test("truseq_small_rna", quality = "perfect", n_reads = 4000,
                           dimer_rate = 0.05, organism = "human",
                           contaminants = list(
                             rRNA = list(sequences = rrna, fraction = 0.25, hotspots = 20),
                             yeast = list(sequences = other, fraction = 0.1)))
  truth <- sim$truth
  expect_equal(nrow(truth), 4000L)
  shares <- table(truth$source) / 4000
  expect_equal(as.numeric(shares[c("rRNA", "yeast", "adapter_dimer", "human")]),
               c(0.25, 0.1, 0.05, 0.6), tolerance = 0.005)
  expect_true(all(grepl("^yeast", truth$seqname[truth$source == "yeast"])))
  expect_lte(nrow(unique(truth[source == "rRNA", .(seqname, insert_start, insert_end)])), 20L)
  expect_gt(nrow(unique(truth[source == "human", .(seqname, insert_start, insert_end)])), 2000L)
  expect_error(simulate_raw_reads(rrna, n_reads = 10, out_prefix = tempfile(), dimer_rate = 0.6,
                                  contaminants = list(x = list(sequences = rrna, fraction = 0.5))),
               "less than 1")
})

test_that("the nuclease cut-site table shapes the read ends", {
  sequences <- raw_test_transcripts(30, 1500)
  ends <- function(enzyme) {
    prefix <- tempfile("rawreads_")
    out <- simulate_raw_reads(sequences, n_reads = 20000, protocol = raw_read_protocol("custom"),
                              out_prefix = prefix, enzyme = enzyme, quality = "perfect", seed = 9)
    truth <- data.table::fread(out$files[["truth"]])
    upstream <- as.character(Biostrings::subseq(sequences[truth$seqname],
                                                truth$insert_start - 1L, width = 1L))
    prop.table(table(upstream))
  }
  none <- ends(NULL)
  expect_equal(as.numeric(none), rep(0.25, 4), tolerance = 0.08)
  table_t1 <- method_end_bias(NULL, "RNase T1")$five_prime_flank_bias$table
  favoured <- table_t1$kmer[which.max(table_t1$weight)]
  t1 <- ends("RNase T1")
  expect_identical(names(which.max(t1)), favoured)
  expect_gt(t1[[favoured]], none[[favoured]] + 0.05)
})

test_that("footprints given as a table are used as given", {
  transcripts <- raw_test_transcripts()
  fragments <- data.table::data.table(seqname = c("tx1", "tx2", "tx40"), start = c(1L, 101L, 871L),
                                      end = c(30L, 128L, 900L), count = c(5L, 1L, 2L))
  prefix <- tempfile("rawreads_")
  out <- simulate_raw_reads(transcripts, fragments = fragments, protocol = "ingolia2012",
                            out_prefix = prefix, quality = "perfect", compress = FALSE, seed = 1)
  truth <- data.table::fread(out$files[["truth"]])
  expect_identical(truth[, .N, keyby = .(seqname, insert_start, insert_end)]$N, c(5L, 1L, 2L))
  reads <- raw_read_fastq(out$files[["fastq"]])
  first <- truth[seqname == "tx40"][1]
  expect_identical(
    substr(reads$sequence[match(first$read_id, reads$read_id)], first$insert_start_read,
           first$insert_end_read),
    as.character(Biostrings::subseq(transcripts[["tx40"]], 871, 900)))
  sampled <- simulate_raw_reads(transcripts, fragments = fragments, n_reads = 500,
                                protocol = "ingolia2012", out_prefix = tempfile("rawreads_"),
                                dimer_rate = 0.2, quality = "perfect", seed = 1)
  expect_equal(sampled$summary$n_reads, 500L)
  expect_equal(sampled$summary$dimer_rate, 0.2)
  expect_error(simulate_raw_reads(transcripts, fragments = data.table::data.table(
    seqname = "tx1", start = 890L, end = 905L), out_prefix = tempfile()), "inside their sequence")
})

test_that("genomic ranges and spliced alignments give the insert on its strand", {
  genome <- raw_test_transcripts(2, 500, seed = 21)
  names(genome) <- c("chr1", "chr2")
  ranges <- GenomicRanges::GRanges(c("chr1", "chr2"), IRanges::IRanges(c(11, 201), width = 28),
                                   strand = c("+", "-"), score = c(2L, 3L))
  out <- simulate_raw_reads(genome, fragments = ranges, protocol = "truseq_small_rna",
                            out_prefix = tempfile("rawreads_"), quality = "perfect", seed = 1)
  truth <- data.table::fread(out$files[["truth"]])
  reads <- raw_read_fastq(out$files[["fastq"]])
  inserts <- substr(reads$sequence, truth$insert_start_read, truth$insert_end_read)
  plus <- as.character(Biostrings::subseq(genome[["chr1"]], 11, 38))
  minus <- as.character(Biostrings::reverseComplement(Biostrings::subseq(genome[["chr2"]], 201, 228)))
  expect_identical(sort(inserts), sort(c(rep(plus, 2), rep(minus, 3))))
  expect_identical(sort(truth$strand), c("+", "+", "-", "-", "-"))

  spliced <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle(factor(c("chr1", "chr1"))), pos = c(50L, 300L),
    cigar = c("10M100N18M", "28M"), strand = GenomicRanges::strand(c("-", "+")))
  out <- simulate_raw_reads(genome, fragments = spliced, protocol = "truseq_small_rna",
                            out_prefix = tempfile("rawreads_"), quality = "perfect", seed = 1)
  truth <- data.table::fread(out$files[["truth"]])
  reads <- raw_read_fastq(out$files[["fastq"]])
  inserts <- substr(reads$sequence, truth$insert_start_read, truth$insert_end_read)
  joined <- Biostrings::xscat(Biostrings::subseq(genome[["chr1"]], 50, 59),
                              Biostrings::subseq(genome[["chr1"]], 160, 177))
  expect_identical(inserts[truth$strand == "-"],
                   as.character(Biostrings::reverseComplement(joined)))
  expect_identical(inserts[truth$strand == "+"],
                   as.character(Biostrings::subseq(genome[["chr1"]], 300, 327)))
})

test_that("read 2 is the reverse complement of the construct", {
  prefix <- tempfile("rawreads_")
  out <- simulate_raw_reads(raw_test_transcripts(), n_reads = 500, protocol = "nextflex_small_rna",
                            out_prefix = prefix, read_length = 60, paired = TRUE,
                            quality = "perfect", seed = 3)
  truth <- data.table::fread(out$files[["truth"]])
  read1 <- raw_read_fastq(out$files[["fastq1"]])
  read2 <- raw_read_fastq(out$files[["fastq2"]])
  expect_identical(read1$read_id, read2$read_id)
  construct <- substr(read1$sequence, 1L, truth$adapter_start_read - 1L)
  mate <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(construct)))
  expect_identical(substr(read2$sequence, 1L, nchar(mate)), mate)
  expect_identical(substr(read2$sequence, nchar(mate) + 1L, nchar(mate) + 12L),
                   rep("GATCGTCGGACT", 500))
})

test_that("a seed reproduces the library and fastp recovers the inserts", {
  first <- raw_simulate_test("mcglincy2017", dimer_rate = 0.05)
  second <- raw_simulate_test("mcglincy2017", dimer_rate = 0.05)
  expect_identical(first$reads$sequence, second$reads$sequence)
  expect_identical(first$truth$umi3, second$truth$umi3)

  skip_if(!nzchar(Sys.which("fastp")), "fastp not installed")
  protocol <- raw_read_protocol("mcglincy2017", add5 = NULL)
  sim <- raw_simulate_test(protocol, quality = "novaseq", read_length = 75)
  trim <- sim$out$trim
  pass1 <- tempfile(fileext = ".fastq")
  pass2 <- tempfile(fileext = ".fastq")
  report <- tempfile()
  run <- function(input, output, arguments) {
    system2("fastp", c("-i", input, "-o", output, strsplit(arguments, " ")[[1]],
                       "--json", paste0(report, ".json"), "--html", paste0(report, ".html")),
            stdout = FALSE, stderr = FALSE)
  }
  expect_identical(run(sim$out$files[["fastq"]], pass1, trim$fastp[["pass1"]]), 0L)
  expect_identical(run(pass1, pass2, trim$fastp[["pass2"]]), 0L)
  trimmed <- raw_read_fastq(pass2)
  truth <- sim$truth[match(sub(" .*$", "", trimmed$read_id), read_id)]
  clean <- truth$source == "main" & truth$n_errors == 0L
  expect_gt(mean(nchar(trimmed$sequence[clean]) == truth$insert_len[clean]), 0.97)
})

test_that("chunked paired FASTQ and truth contents reproduce across compression workers", {
  transcripts <- raw_test_transcripts(4L)
  fragments <- data.table::data.table(seqname = c("tx1", "tx2", "tx3"), start = c(21L, 77L, 91L),
                                      end = c(48L, 109L, 118L), count = c(3001L, 0L, 3002L))
  run <- function(workers) {
    simulate_raw_reads(transcripts, fragments = fragments, read_length = 150L,
                       protocol = "polya_cats", paired = TRUE, library = "seeded",
                       out_prefix = tempfile(), seed = 17L, chunk_size = 4500L,
                       dimer_rate = 0.05, threads = workers)
  }
  first <- run(1L)
  second <- run(2L)
  paths <- c("fastq1", "fastq2", "truth")
  for (name in paths) {
    expect_identical(readLines(first$files[[name]]), readLines(second$files[[name]]))
  }
  truth <- data.table::fread(first$files[["truth"]])
  expect_identical(truth[source == "main", .N, keyby = seqname]$N, c(3001L, 3002L))
  expect_identical(truth$read_id, raw_read_fastq(first$files[["fastq2"]])$read_id)
})

test_that("sampled copies of one supplied footprint retain their end shifts", {
  result <- simulate_raw_reads(raw_test_transcripts(1L), fragments = data.frame(
    seqname = "tx1", start = 20L, end = 48L), n_reads = 41L, end_bias = TRUE,
    out_prefix = tempfile(), chunk_size = 13L, quality = "perfect", seed = 9L,
    compress = FALSE)
  truth <- data.table::fread(result$files[["truth"]])
  expect_equal(nrow(unique(truth[, .(insert_start, insert_end)])), 1L)
})

test_that("raw reads retain single-worker throughput at one million reads", {
  skip_on_cran()
  transcripts <- raw_test_transcripts(256L)
  elapsed <- system.time({
    out <- simulate_raw_reads(transcripts, n_reads = 1000000L, read_length = 150L,
                              protocol = "mcglincy2017", out_prefix = tempfile(),
                              seed = 42L, threads = 1L)
  })
  on.exit(unlink(out$files))
  expect_equal(out$summary$n_reads, 1000000L)
  # CPU time excludes contention on shared CI runners. Leave room for zlib
  # builds and slower machines while catching the previous implementation.
  expect_lt(unname(sum(elapsed[c("user.self", "sys.self", "user.child", "sys.child")])), 8)
})

test_that("native truth output preserves CSV quoting and ambiguous DNA", {
  dna <- Biostrings::DNAString("ACGTMRSVWYHKDBN-+.")
  label <- 'tx,"quoted"'
  fragments <- data.frame(seqname = label, start = 1L, end = length(dna))
  out <- simulate_raw_reads(stats::setNames(as.character(dna), label), fragments = fragments,
                            library = 'lib,"one"', organism = "a,b", paired = TRUE,
                            out_prefix = tempfile(), quality = "perfect", compress = FALSE)
  truth <- utils::read.csv(out$files[["truth"]], colClasses = "character")
  expect_identical(truth$seqname, label)
  expect_identical(truth$source, "a,b")
  expect_identical(truth$umi5, "")
  expect_identical(truth$read_id, 'lib,"one".0000001')
  expect_identical(substr(raw_read_fastq(out$files[["fastq1"]])$sequence, 1L, length(dna)),
                   as.character(dna))
  expect_identical(substr(raw_read_fastq(out$files[["fastq2"]])$sequence, 1L, length(dna)),
                   as.character(Biostrings::reverseComplement(dna)))
})
