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
#' @param n = 500. Number of genes to create, one whole number of at least 1.
#' @param length_model "fixed" (default, the old 300 / 140 / 105 nt for every gene) or
#'   "human": for each of \code{cds_length}, \code{leader_length} and
#'   \code{trailer_length} that is not given, draw per-gene log-normal lengths fitted
#'   to human Ensembl 108 canonical transcripts (\code{\link{human_tx_lengths}}; medians
#'   136 / 1284 / 1018 nt incl. start and stop codon). Mean transcript is about 3,800 nt
#'   instead of 551 nt, so genome and read simulation take longer.
#' @param cds_length = NULL (300 for every gene with length_model "fixed"). Coding sequence length per gene, in
#'   nucleotides, \strong{without} the start and stop codon, which are added to it
#'   (300 gives a 306 nt CDS). One whole number per gene, at least 9 and at
#'   least \code{uorf_max_length} + 3, divisible by 3 since a coding sequence is
#'   read in whole codons.
#' @param cds_exons = 2. How many exons the coding sequence is split into, one
#'   whole number of at least 1 and smaller than every \code{cds_length}.
#'   More than one means reads can cross a splice junction, which is what
#'   exercises the spliced-alignment handling downstream.
#' @param cds_intron_length = sample(c(30, 30), size = n, replace = TRUE).
#'   Intron length, including the GT motif at its 5' end; only used (and
#'   checked) when \code{cds_exons} > 1. Whole numbers of at least 3, divisible
#'   by 3, one per gene or one value for all genes.
#' @param leader_length = NULL (140 with "fixed"). Length of the 5' untranslated region
#'   ahead of the start codon, one whole number of at least 1 per gene. This is
#'   where uORFs are placed, so a short leader leaves little room for them.
#'   The length does not change how many reads the leader gets (that is
#'   \code{region_proportion} in \code{simCountTablesRegions()}), only how
#'   densely they lie: the same reads spread over a longer leader.
#' @param trailer_length = NULL (105 with "fixed"). Length of the 3' untranslated region
#'   after the stop codon, one whole number of at least 0 per gene (0 = the
#'   transcript ends with the stop codon). As for the leader, the length sets
#'   the read density there, not the number of reads.
#' @param gene_names = paste0("ENSGTEST1000", seq(n)),
#' @param tx_names = paste0("ENSTTEST1000", seq(n)),
#' @param seqnames Character vector of chromosome names. When NULL,
#'   the default compact layout uses up to \code{chromosome_count} chromosomes.
#' @param chromosome_layout Either \code{"compact"} (the default), which places
#'   multiple genes on chromosomes, or
#'   \code{"legacy_one_gene_per_contig"}, which recreates the historical
#'   \code{chr1}, ..., \code{chrN} layout. Only used when \code{seqnames} is NULL.
#' @param chromosome_count = 24. Maximum number of chromosomes in the default
#'   compact layout, one whole number of at least 1.
#' @param chromosome_weights Optional sampling weights controlling how genes are
#'   distributed among chromosomes. Must have one positive value per chromosome.
#' @param flank_length = 305. Genomic sequence added on either side of each
#'   gene: whole numbers of at least 0, one per gene or one value for all
#'   genes. Reads never fall here, but the space means a gene is not flush
#'   against a chromosome edge, where coordinate handling behaves differently.
#' @param strand = c("+", sample(c("+", "-"), n-1, replace = T)). Which strand
#'   each gene sits on. The default deliberately mixes both: minus-strand genes
#'   run against the genome's coordinates and are a common source of
#'   off-by-one and reverse-complement errors, so a test genome should contain
#'   them.
#' @param start_codons = "ATG". Codon(s) that open a coding sequence: three
#'   letters from A, C, G, T each.
#' @param stop_codons = c("TAG", "TGA", "TAA"), a non-empty subset (no
#'   duplicates) of the standard genetic code's three stop codons. Use a
#'   smaller subset to simulate stop-codon usage bias (e.g. \code{"TGA"}
#'   alone). Codons outside \code{c("TAG", "TGA", "TAA")} are not supported:
#'   the CDS/uORF integrity checks rely on \code{Biostrings::translate()}'s
#'   default (standard) genetic code to recognize stops, so simulating an
#'   organism/organelle with a different genetic code (e.g. mitochondrial,
#'   where TGA is not a stop) is out of scope here.
#' @param max_uorfs = 1. Upper limit on how many upstream ORFs a gene may get,
#'   one whole number of at least 0 (0 = no uORFs). A gene can get fewer when
#'   its leader has no room left; a leader shorter than 4 nt gets none. If no
#'   gene has room for a uORF, simGenome() stops with an error.
#'   These are short reading frames in the leader that ribosomes can translate
#'   before reaching the main coding sequence, so they carry real signal and are
#'   not merely decoration.
#' @param n_uorfs NULL (default) or the exact number of uORFs per gene: one whole
#'   number for all genes or one per gene. Unlike \code{max_uorfs} (an upper limit;
#'   a gene silently gets fewer when its leader has no room), every gene must then get
#'   exactly this many: genes that cannot (leader too short for the request under the
#'   overlap settings) stop the run before any sequence is made, and genes that end up
#'   with fewer after the random placement are reported (see \code{on_infeasible}).
#'   Overrides \code{max_uorfs}.
#' @param on_infeasible "error" (default) or "warning": what to do when the
#'   \code{n_uorfs} request cannot be met, before building (impossible) or after it
#'   (placement left a gene short). Either way the details are in the build log.
#'   A per-gene build log (requested and built uORFs, leader length, minimum leader
#'   needed, status and reason) is always written next to the genome as
#'   \code{<genome_name>_build_log.tsv} and returned as \code{"build_log"}.
#' @param architecture NULL (default) or the extra translons each gene should
#'   contain, built exactly and verified: a \code{\link{gene_architecture}} table
#'   (all genes), a preset name from \code{\link{architecture_preset}}, or a list
#'   with one of these (or NULL) per gene. Classes: uORF, uoORF, NTE, intORF, doORF,
#'   dORF; each translated or a decoy. Positions are chosen before any sequence is
#'   made, so impossible requests stop early with the gene and reason. Every other
#'   ORF found in the final transcripts is labelled as an incidental decoy. Writes
#'   \code{true_features.rds} (genomic GRangesList, one element per feature) and
#'   \code{<genome_name>_features.tsv}; translated uORFs and uoORFs are also written
#'   to \code{true_uORFs.rds} so the read simulation uses them. All translated
#'   extra features of any class are read through the single "uorf" region of
#'   \code{simCountTablesRegions()}: its share (\code{region_proportion$uorf}) is the
#'   TOTAL for all extra translons of a gene, split between them uniformly, by length,
#'   or by \code{simRDGflow()}; \code{annotate_ground_truth()} maps reads back to features. Cannot be combined
#'   with \code{max_uorfs} or \code{n_uorfs}.
#' @param uorfs_can_overlap = TRUE. Whether two uORFs of the same gene may share
#'   sequence: TRUE or FALSE.
#' @param uorfs_can_overlap_cds = TRUE. Whether a uORF may run into the CDS:
#'   0 (FALSE) = no overlap, 1 (TRUE) = may overlap, 2 = must overlap.
#'   Such a uORF starts in the leader and ends after the CDS start codon, so it
#'   needs \code{uorf_max_length} of at least 12: with 2 a smaller value
#'   stops with an error, with 1 all uORFs then stay in the leader.
#' @param export_txdb TRUE or FALSE, default TRUE. Export TxDb object of gtf.
#'   Highly advised, but for testing you can turn it off to speed up raw output.
#' @param uorf_max_length = 60. Longest a uORF may be, in nucleotides
#'   including its start and stop codon: one whole number of at least 6 (a
#'   start codon directly followed by a stop), divisible by 3. Real
#'   uORFs are short, which is part of what distinguishes them from a second
#'   coding sequence.
#' @param debug_on TRUE or FALSE, default TRUE. Activate debug on suspected error. This can happen
#' on uORF sampling, with corner case created transcripts. Useful to debug to understand why it failed.
#' If you are not interested and you believe settings are valid, just rerun until it works.
#' @import Biostrings ORFik data.table
#' @return a named character vector with paths to output files.
#' genome (fasta), transcriptome (gtf), transcript Database (TxDb) and uorfs (rds)
#' The settings this run used, and the random state it started from, are saved
#' next to the genome as \code{<genome_name>_simulation_record.rds};
#' \code{simNGScoverage()} adds them to its \code{simulation_report.tsv}.
#' @export
#' @examples
#' # 6 genes on 6 chromosomes'
#' # Default stores to tempdir (lost after restart, so save to location
#' #  if you need it!)
#' simGenome(n = 6)
simGenome <- function(n = 500,
                      out_dir = tempdir(),
                      genome_name = "Homo_sapiens_dummy",
                      cds_length = NULL,
                      cds_exons = 2,
                      cds_intron_length = sample(c(30, 30), size = n, replace = TRUE),
                      leader_length = NULL,
                      trailer_length = NULL,
                      length_model = c("fixed", "human"),
                      gene_names = paste0("ENSGTEST1000", seq(n)),
                      tx_names = paste0("ENSTTEST1000", seq(n)),
                      seqnames = NULL,
                      chromosome_layout = c("compact", "legacy_one_gene_per_contig"),
                      chromosome_count = 24L,
                      chromosome_weights = NULL,
                      flank_length = 305,
                      strand = c("+", sample(c("+", "-"), n-1, replace = TRUE)),
                      start_codons = "ATG",
                      stop_codons = c("TAG", "TGA", "TAA"),
                      max_uorfs = 1,
                      n_uorfs = NULL,
                      on_infeasible = c("error", "warning"),
                      architecture = NULL,
                      uorfs_can_overlap = TRUE,
                      uorfs_can_overlap_cds = TRUE,
                      uorf_max_length = 60,
                      export_txdb = TRUE,
                      debug_on = TRUE) {
  # Read before any setting is used, so recording cannot change the draws.
  rng_state <- current_rng_state()
  supplied <- names(match.call())[-1]
  chromosome_layout <- match.arg(chromosome_layout)
  check_whole_numbers(n, "n", min = 1)
  # Region lengths not given by the caller: the old fixed 300 / 140 / 105 nt, or
  # per-gene log-normal draws fitted to human canonical transcripts.
  length_model <- match.arg(length_model)
  if (is.null(cds_length) || is.null(leader_length) || is.null(trailer_length)) {
    lens <- if (length_model == "human") human_tx_lengths(n) else
      list(leader_length = rep.int(140, n), cds_length = rep.int(300, n),
           trailer_length = rep.int(105, n))
    if (is.null(cds_length)) cds_length <- lens$cds_length
    if (is.null(leader_length)) leader_length <- lens$leader_length
    if (is.null(trailer_length)) trailer_length <- lens$trailer_length
  }
  # chromosome_count only counts when the chromosomes are not named directly.
  chromosome_count_used <- is.null(seqnames) &&
    chromosome_layout != "legacy_one_gene_per_contig"
  seqnames <- resolve_simulated_seqnames(
    n = n,
    seqnames = seqnames,
    chromosome_layout = chromosome_layout,
    chromosome_count = chromosome_count
  )
  on_infeasible <- match.arg(on_infeasible)
  genome_input_test_controller()
  # Exact uORF counts per gene, checked against each leader before anything is built.
  exact_uorfs <- !is.null(n_uorfs)
  if (!is.null(architecture)) {
    if (exact_uorfs || ("max_uorfs" %in% supplied && any(max_uorfs > 0)))
      stop("Use either architecture or max_uorfs/n_uorfs, not both", call. = FALSE)
    max_uorfs <- 0L
  }
  if (exact_uorfs) {
    check_whole_numbers(n_uorfs, "n_uorfs", min = 0, per_gene = n, allow_single = TRUE)
    uorfs_per_gene <- rep_len(as.integer(n_uorfs), n)
  } else uorfs_per_gene <- rep_len(as.integer(max_uorfs), n)
  feasibility <- uorf_feasibility(uorfs_per_gene, rep_len(leader_length, n),
                                  uorf_max_length, uorfs_can_overlap,
                                  uorfs_can_overlap_cds, tx_names)
  if (exact_uorfs) stop_if_infeasible(feasibility, on_infeasible)
  max_uorfs <- uorfs_per_gene
  # Place every architecture feature now, so impossible requests stop here.
  if (!is.null(architecture)) {
    arch <- place_features(architecture_per_gene(architecture, n, tx_names),
                           rep_len(leader_length, n), rep_len(cds_length, n) + 6L,
                           rep_len(trailer_length, n))
  }
  # The next step turns these two into per-gene coordinates, so keep them as
  # they were given for the record.
  seqnames_given <- seqnames
  flank_length_given <- flank_length
  genome_exon_flanks_controller()

  ## Create annotation
  # Transcriptome
  make_transcriptome_sequences_controller()
  if (!is.null(architecture)) {
    built <- build_architectures(arch, leader_string, cds_string, trailer_string,
                                 start_codons, stop_codons, tx_names)
    leader_string <- built$leader_string; trailer_string <- built$trailer_string
    cds_string <- built$cds_string
    cds_string_genomic <- splice_cds_back(cds_string_genomic, as.character(cds_string),
                                          cds_exons, cds_intron_length)
    gene_string <- DNAStringSet(paste0(leader_string, cds_string_genomic, trailer_string))
    names(gene_string) <- tx_names
    gene_string[names(strand)[strand == "-"]] <-
      reverseComplement(gene_string[names(strand)[strand == "-"]])
    features <- built$features; rm(built)
  }
  # Chromosomes
  make_genome_sequences_controller()
  ## Make exon ranges
  genome_exon_ranges_controller()
  ## Sub-annotation ranges
  # uORFS
  if (any(max_uorfs > 0)) {
    uorf_list <- create_uORFs(leader_string, chromosome_seqs, n,
                              max_uorfs, start_codons, stop_codons,
                              leader_ranges, cds_ranges, cds_grl, seqnames, strand,
                              uorf_max_length, uorfs_can_overlap_cds, uorfs_can_overlap,
                              tx_names, debug_on)
    uorf_ranges <- uorf_list[[1]]; chromosome_seqs <- uorf_list[[2]]; rm(uorf_list)
  }
  build_log <- uorf_build_log(feasibility,
                              if (exists("uorf_ranges", inherits = FALSE)) uorf_ranges else NULL)
  if (!is.null(architecture)) build_log <- architecture_build_log(features, tx_names)
  if (exact_uorfs && any(build_log$status != "ok")) {
    short <- build_log[status != "ok"]
    text <- paste0(nrow(short), " gene(s) got fewer uORFs than requested: ",
                   paste(head(paste0(short$transcript_id, " (", short$built_uorfs, "/",
                                     short$requested_uorfs, ")"), 10), collapse = ", "),
                   ". See the build log; another seed or longer leaders may help.")
    if (on_infeasible == "error") stop(text, call. = FALSE) else warning(text, call. = FALSE)
  }
  # Validate simulation
  genome_output_test_controller()
  # Make gtf as ranges object with metadata
  gtf_ranges <- make_gtf_from_ranges(leader_ranges, cds_ranges, trailer_ranges,
                                     exon_cds_ranges, gene_ranges,
                                     transcript_ranges, tx_names, gene_names)
  export_genome_output_controller()
  simGenome <- c(simGenome, build_log = unname(write_genome_build_log(build_log, simGenome[["genome"]])))
  if (!is.null(architecture)) {
    simGenome <- c(simGenome, export_architecture_features(
      features, leader_ranges, exon_cds_ranges, trailer_ranges, tx_names,
      out_dir, simGenome[["genome"]]))
  }
  # Recorded last, once every setting has been used in its own order.
  step_record <- make_step_record(
    "simGenome", sys.function(), environment(), supplied, rng_state,
    not_used = if (!chromosome_count_used) "chromosome_count"
  )
  step_record$settings$seqnames <- seqnames_given
  step_record$settings$flank_length <- flank_length_given
  saveRDS(step_record, genome_record_path(simGenome["genome"]))
  return(simGenome)
}
