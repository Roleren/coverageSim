#' Simulate genome
#'
#' Step 1 of 4 in the package. Create a genome of 'n' genes distributed on
#' 'seqnames' number of chromosomes.
#' You can specify all details, 5' UTR length, CDS length, CDS exons,
#' cds intron length, number of uORFs, how they overlap etc.
#'
#' Theory:\cr
#' GC content per region:\cr
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC403801/
#' @param out_dir = "~/Desktop/benchmark_ORFPrediction",
#' @param genome_name = "Homo_sapiens_dummy",
#' @param n = 500, integer, number of genes to create
#' @param cds_length = rep.int(300, n)
#' @param cds_exons = 2,
#' @param cds_intron_length = sample(c(30, 30), size = n, replace = TRUE),
#' only used any cds_exons > 1. Is the size excluding the '5 GT motif.
#' @param leader_length = rep.int(140, n),
#' @param trailer_length = rep.int(105, n),
#' @param gene_names = paste0("ENSGTEST1000", seq(n)),
#' @param tx_names = paste0("ENSTTEST1000", seq(n)),
#' @param seqnames = \code{list(paste0("chr", seq(n)), paste0("chr", seq(23)))[[1]]},
#' @param flank_length = 305,
#' @param strand = c("+", sample(c("+", "-"), n-1, replace = T)),
#' @param start_codons = "ATG",
#' @param stop_codons = c("TAG", "TGA", "TAA"),
#' @param max_uorfs = 1,
#' @param uorfs_can_overlap = TRUE,
#' @param uorfs_can_overlap_cds = TRUE, logical/integer, 0/FALSE is no overlap, 1/TRUE is can overlap
#' , 2 is must overlap.
#' @param export_txdb logical, default TRUE. Export TxDb object of gtf. Highly adviced,
#' but for testing you can turn it of to speed up raw output.
#' @param uorf_max_length = 60
#' @import Biostrings ORFik data.table
#' @return a named character vector with paths to output files.
#' genome (fasta), transcriptome (gtf), transcript Database (TxDb) and uorfs (rds)
#' @export
#' @examples
#' # 6 genes on 6 chromosomes'
#' # Default stores to tempdir (lost after restart, so save to location
#' #  if you need it!)
#' simGenome(n = 6)
simGenome <- function(n = 500,
                      out_dir = tempdir(),
                      genome_name = "Homo_sapiens_dummy",
                      cds_length = rep.int(300, n),
                      cds_exons = 2,
                      cds_intron_length = sample(c(30, 30), size = n, replace = TRUE),
                      leader_length = rep.int(140, n),
                      trailer_length = rep.int(105, n),
                      gene_names = paste0("ENSGTEST1000", seq(n)),
                      tx_names = paste0("ENSTTEST1000", seq(n)),
                      seqnames = list(paste0("chr", seq(n)), paste0("chr", seq(23)))[[1]],
                      flank_length = 305,
                      strand = c("+", sample(c("+", "-"), n-1, replace = TRUE)),
                      start_codons = "ATG",
                      stop_codons = c("TAG", "TGA", "TAA"),
                      max_uorfs = 1,
                      uorfs_can_overlap = TRUE,
                      uorfs_can_overlap_cds = TRUE,
                      uorf_max_length = 60,
                      export_txdb = TRUE,
                      debug_on = TRUE) {
  genome_input_test_controller()
  genome_exon_flanks_controller()

  ## Create annotation
  # Transcriptome
  make_transcriptome_sequences_controller()
  # Chromosomes
  make_genome_sequences_controller()
  ## Make exon ranges
  genome_exon_ranges_controller()
  ## Sub-annotation ranges
  # uORFS
  if (max_uorfs > 0) {
    uorf_list <- create_uORFs(leader_string, chromosome_seqs, n,
                              max_uorfs, start_codons, stop_codons,
                              leader_ranges, cds_ranges, cds_grl, seqnames, strand,
                              uorf_max_length, uorfs_can_overlap_cds, uorfs_can_overlap,
                              tx_names, debug_on)
    uorf_ranges <- uorf_list[[1]]; chromosome_seqs <- uorf_list[[2]]; rm(uorf_list)
  }
  # Validate simulation
  genome_output_test_controller()
  # Make gtf as ranges object with metadata
  gtf_ranges <- make_gtf_from_ranges(leader_ranges, cds_ranges, trailer_ranges,
                                     exon_cds_ranges, gene_ranges,
                                     transcript_ranges, tx_names, gene_names)
  export_genome_output_controller()
  return(simGenome)
}
