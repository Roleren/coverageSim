# Raw reads with a known answer: simulate_raw_reads()
#
# simNGScoverage() ends at aligned footprints. simulate_raw_reads() goes the
# other way from there: it wraps footprints in what a library protocol and a
# sequencer add (UMIs, barcodes, untemplated bases, poly(A) tail, linker,
# adapter, read length, qualities, errors) and writes FASTQ plus a truth table.
# Use it to test trimming and layout detection against a known answer.
#
# Run from the repository root:
#   env -u LC_ALL Rscript analysis/workflows/simulate_raw_reads_example.R [out_dir]

devtools::load_all(".", quiet = TRUE)
arguments <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(arguments) > 0L) arguments[1] else file.path(tempdir(), "raw_reads_example")

## 1. Layouts -----------------------------------------------------------------
# Presets for published protocols, each with the trimming that undoes it.
print(raw_read_protocols()[, .(preset, method, umi5, umi3, barcode3, polya, adapter,
                               trim_front, trim_tail)])
print(raw_read_protocol("mcglincy2017"))
# Any field can be overridden, or build a layout from nothing.
pooled <- raw_read_protocol("custom", name = "pooled_lane", umi5 = 6,
                            barcode5 = c("ACGTAC", "TGCATG", "GATCGA"),
                            adapter3 = "CTGTAGGCACCATCAAT")

## 2. Footprints drawn from sequences -----------------------------------------
# Any sequences work: transcripts, CDSs, a small genome. Here a toy transcriptome.
set.seed(1)
random_sequences <- function(n, width, prefix) {
  sequences <- vapply(seq_len(n), function(i) {
    paste(sample(c("A", "C", "G", "T"), width, replace = TRUE), collapse = "")
  }, character(1))
  Biostrings::DNAStringSet(stats::setNames(sequences, paste0(prefix, seq_len(n))))
}
transcripts <- random_sequences(300, 1500, "tx")
rrna <- random_sequences(4, 2000, "rRNA")
yeast <- random_sequences(50, 1200, "yeast_tx")

library_one <- simulate_raw_reads(
  transcripts, n_reads = 200000, protocol = "mcglincy2017",
  out_prefix = file.path(out_dir, "mcglincy_rnase1"),
  read_length = 75, enzyme = "RNase I", dimer_rate = 0.05, organism = "human",
  contaminants = list(
    rRNA = list(sequences = rrna, fraction = 0.2, hotspots = 200),  # few, deep fragments
    yeast = list(sequences = yeast, fraction = 0.05)),              # second organism
  quality = "novaseq", seed = 1)
print(t(library_one$summary))
truth <- data.table::fread(library_one$files[["truth"]])
print(truth[1:5])
print(truth[, .(reads = .N, mean_insert = mean(insert_len)), by = source])

## 3. Round trip: the known settings recover the inserts -----------------------
trim <- library_one$trim
if (nzchar(Sys.which("fastp"))) {
  pass1 <- file.path(out_dir, "pass1.fastq.gz")
  pass2 <- file.path(out_dir, "trimmed.fastq.gz")
  fastp <- function(input, output, arguments) {
    system2("fastp", c("-i", input, "-o", output, strsplit(arguments, " ")[[1]],
                       "--json", file.path(out_dir, "fastp.json"),
                       "--html", file.path(out_dir, "fastp.html")), stdout = FALSE, stderr = FALSE)
  }
  fastp(library_one$files[["fastq"]], pass1, trim$fastp[["pass1"]])  # adapter + 5' bases
  fastp(pass1, pass2, trim$fastp[["pass2"]])                         # 3' UMI + barcode
  trimmed <- Biostrings::readDNAStringSet(pass2, format = "fastq")
  recovered <- data.table::data.table(read_id = sub(" .*$", "", names(trimmed)),
                                      trimmed_len = Biostrings::width(trimmed))
  check <- merge(truth, recovered, by = "read_id")
  # The untemplated 5' T of circularization stays on about 10 percent of reads.
  print(check[source != "adapter_dimer",
              .(exact = mean(trimmed_len == insert_len),
                exact_or_added_base = mean(trimmed_len == insert_len + nchar(add5)))])
}

## 4. Footprints from coverageSim or from a table ------------------------------
# A table of positions on the sequences (count = copies) is used as given:
fragments <- data.table::data.table(seqname = c("tx1", "tx2"), start = c(100L, 250L),
                                    end = c(128L, 279L), count = c(500L, 300L))
simulate_raw_reads(transcripts, fragments = fragments, protocol = pooled,
                   out_prefix = file.path(out_dir, "pooled_lane"), read_length = 50, seed = 2)
# A coverageSim library goes in as genomic ranges with the genome as sequences:
#   genome <- simGenome(n = 6, max_uorfs = 0)
#   experiment <- simNGScoverage(genome, simCountTablesRegions(simCountTables(
#     loadRegion(genome["txdb"], "cds"), libtypes = "RFP", print_statistics = FALSE)))
#   footprints <- ORFik::fimport(ORFik::filepath(experiment, "default")[1])
#   simulate_raw_reads(genome[["genome"]], fragments = footprints, protocol = "ottr",
#                      out_prefix = file.path(out_dir, "from_coverageSim"), read_length = 100)

## 5. Other layouts ------------------------------------------------------------
for (preset in c("ingolia2012", "nextflex_small_rna", "smarter", "ottr")) {
  out <- simulate_raw_reads(transcripts, n_reads = 20000, protocol = preset,
                            out_prefix = file.path(out_dir, preset), read_length = 100,
                            paired = preset == "nextflex_small_rna", seed = 3)
  cat(preset, ": adapter", out$trim$adapter, "front", out$trim$trim_front, "tail",
      out$trim$trim_tail, "\n")
}
cat("Output in", out_dir, "\n")
