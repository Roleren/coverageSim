# coverageSim: plan for RDG-structured simulation (draft, claude-megafst 2026-11-18)

Goal (Håkon): coverageSim should be able to test tools on very complicated Ribosome Decision Graph (RDG) structures:
overlapping translons, uORF/uoORF/dORF/doORF, internal out-of-frame ORFs, N/C-terminal extensions, leaky scanning,
reinitiation, frameshifts and readthrough, with ground truth precise enough to score ORF callers.
Inputs: code review of branch review/carina-dev-claude-20261118, the manuscript (summary by carina-claude), the
literature table (no published read-level simulator has RDG coupling; TASEP/ICIER models have the biology but no reads;
RiboBA 2026 is the closest competitor and resamples real reads). Pending: Codex test/bug/performance report, realism
benchmark, Carina's design notes. This plan will be updated with them.

## 1. What is there and what is missing

| Area | Now | Gap for RDG testing |
|---|---|---|
| Features | leader uORFs, may overlap each other or run out of frame into the CDS start | no dORF/doORF, internal out-of-frame ORF, NTE (in-frame upstream start), CTE/readthrough, frameshift, near-cognate main starts, isoforms sharing ORFs |
| Construction | random sequence, then rejection-sampled uORF positions, codons overwritten on the chromosome, internal stops mutated away, CDS repaired, broken uORFs dropped | `max_uorfs` is only "at most": a request for 2 uORFs can silently give 0-1; no feasibility check before building; no per-gene record of requested vs built; random leaders contain unlabelled ATG ORFs; `browser()` in package code |
| Reads per feature | gene count split into leader/CDS/trailer/uORF by fixed proportions (default 5/75/10/10 %), uORF share split uniformly or by length | no coupling: uORF translation never takes ribosomes from the CDS; no leaky scanning/reinitiation; conditions change only gene totals, so "uORF up, CDS down" (ISR) cannot be simulated |
| Truth | true_uORFs.rds; per read `sampling_group` = region:index | no translated/untranslated status for every ORF in the sequence (no decoys with labels); no per-read feature ID across all feature types |
| Performance | per-gene R loops; per-uORF `replaceLetterAt` copies the whole chromosome | manuscript: 36 s for 5,000 genes without uORFs vs 420 s with 2 uORFs/gene |

## 2. Gene architecture: declare, check, build, verify, log

Replace generate-then-patch with a constructive builder.

1. **Declare.** A per-gene architecture spec, as a data.table (one row per feature) or via helper constructors:
   `feature` (uORF, uoORF, CDS, intORF, dORF, doORF, NTE, CTE/readthrough, frameshift), `length` or range,
   `frame` relative to the CDS, `start_codon` (ATG/CTG/GTG/...), `kozak` (strong/weak), `anchor` (leader/CDS/trailer,
   position or "any"), `overlaps` (ids it must or must not overlap), `status` (translated / decoy), plus gene-level
   lengths. Presets cover common designs ("2 uORF + CDS", "uoORF + CDS", "CDS + dORF", "ATF4-like",
   "frameshift antizyme-like", "decoy-rich leader").
2. **Check (fail early).** Before any sequence exists, compute each gene's minimal leader/CDS/trailer length from its
   features and their spacing (start/stop codons, required overlaps, frames, the rule that uORF/CDS overlaps must be
   out of frame, etc.). Impossible genes stop with one error that lists the gene, the request and the reason
   (e.g. "gene 7: 2 non-overlapping uORFs of >= 30 nt need a leader of >= 66 nt, got 40"), or with
   `on_infeasible = "extend"` the leader/trailer is lengthened and logged. No silent fewer features.
3. **Build.** Place feature coordinates by sampling inside the feasible windows. Then write sequence position by position
   from a constraint table: each nucleotide knows the codon position it has in every overlapping translated feature, so
   internal codons are drawn to avoid stops in all of their frames, starts and stops are fixed, and (optionally)
   unwanted ATG/near-cognate starts are suppressed in leaders. Small backtracking windows resolve conflicts; no global
   mutate-and-repair pass. Genes are built as transcripts first and concatenated into chromosomes once (removes the
   per-uORF whole-chromosome copy).
4. **Verify.** Re-detect all ORFs on the final transcripts (ORFik::findORFs for all start codons used). Every requested
   feature must exist exactly (start, stop, frame); everything else is an annotated decoy. Verification failure is an
   error, not a debugger.
5. **Log.** Per run: `genome_build_log.tsv` (gene, feature, requested, built, adjusted, reason),
   `features.gtf`/`features.tsv` with class, frame, start codon, status (translated / decoy) and ids used by the read
   simulator, and a short summary in `simulation_report.tsv`. No `browser()` in package code; `debug_on` only adds detail
   to the log.

## 3. Ribosome flow on the RDG (coupled reads per feature)

A per-transcript scanning/initiation model (ICIER/TASEP-inspired, kept simple and explicit):
- Scanning flux enters at the 5' cap. At every start of a translated feature the scanning ribosome initiates with
  probability `p_init` (by start codon and Kozak strength), otherwise leaks past (leaky scanning).
- After a uORF stop, it reinitiates downstream with probability `p_reinit` (rising with distance, decreasing with uORF
  length), else releases. Optional elongating-ribosome collisions and stalls (queue at a pause site; emit disome
  reads). Frameshift: fraction `p_fs` continues in the new frame at the shift site. Readthrough: fraction `p_rt`
  continues past the stop into the CTE.
- Output: ribosome load per feature (initiations x length / elongation rate), which sets each feature's share of the
  gene's Ribo count. Leader/trailer background stays as now.
- Conditions: the count-table conditions can change `p_init`, `p_reinit`, `p_fs` and `p_rt` per gene, e.g. ISR mode:
  lower eIF2-TC means more leaky scanning and delayed reinitiation, so ATF4-like genes gain CDS and lose uORF reads.
  This makes differential-translation and uORF-regulation tools testable.
- Backwards compatible: the current fixed region proportions stay available as `region_model = "proportions"`; the new
  model is `region_model = "rdg"`.

## 4. Reads and ground truth

- Each feature is sampled in its own frame with the existing codon/autocorrelation/RNase/fragment machinery, so
  overlapping features give mixed-frame signal by construction (as uoORF/CDS overlaps do now).
- Ground truth per read: `feature_id`, `feature_class`, `frame`, `ribosome_site`, plus the existing fragment columns.
  Per feature: status (translated/decoy), simulated load, reads. Decoys get 0 reads by construction but sit in realistic
  sequence context, so callers can be scored for FDR.
- Scoring helpers: given a caller's ORF table, compute precision/recall/MCC per feature class, overlap type, depth and
  start codon (the manuscript already has MCC functions to extend).

## 5. Performance

Profile first (Codex report pending), then replace hot R loops with C++ (Rcpp) where they dominate:
`apply_autocorrelation_kernel` (per-position vapply, O(L x k) in R), `rolling_weighted_sum`, per-gene
`lapply(eval(ideal_cov))`, sequence construction, and fragment candidate generation. Target: 5,000 genes with
complex architectures in roughly the time of the current no-uORF run. Results must stay identical to the R version
(seeded tests comparing both).

## 6. Order of work (branch feature/rdg-architecture)

1. Architecture spec + feasibility check + log (no sequence change yet), with tests for impossible cases.
2. Constructive builder + verifier + decoy annotation; tests that every requested feature exists exactly, on both
   strands and with introns.
3. Feature-aware read sampling + per-read ground truth for all feature types.
4. RDG flow model (`region_model = "rdg"`) + condition effects + tests on analytic expectations (e.g. one uORF with
   p_init = 1 and p_reinit = 0 leaves no CDS reads).
5. C++ hot spots, identical-output tests, benchmark table.
6. Scoring helpers and an example benchmark of 2-3 ORF callers (ORFik-based scoring, RiboCode, PRICE if available).

Each step: tests first (AGENTS.md), Codex independent review of the diff, no change to default behaviour until the
step is complete.
