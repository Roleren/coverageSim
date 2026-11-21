# Explainability loop -> coverageSim ports

The explainability loop (claude-megafst for Håkon Tjeldnes, from 2026-11-18) explains real human Ribo-seq coverage
window by window, checks the explainer on coverageSim data with known truth, and puts every adopted finding back into
the simulator. Rule (Håkon, 2026-11-20): an iteration is complete only when its finding is in coverageSim as an opt-in
setting with a test, or this file says why it cannot be simulated.

This file lists every iteration of the ledger with: the finding, the evidence (report and numbers), the simulator
setting that reproduces it (or why there is none), its default, and how to switch it on. Paths are on the shared
server: `OW` = `/home/rstudio/shared_workspace/open_workspace`, `BC` = `OW/backclaude`, `MF` = `OW/claude-megafst`,
ledger = `MF/explainability_loop/ITERATIONS.md` (+ `README.md`, `metrics.csv`). Real numbers are from the unique-mapper
track of the merged human collection unless noted. "NU" = not understood.

All settings below are opt-in: with none of them set, every simulation is identical to the one before the port.

## Switches at a glance

| Port | Setting | Default | Human setting |
| --- | --- | --- | --- |
| S2 | `codon_bias_options(pause_sd, pause_prob, per_genome_draws, sharpen)` | off (`codon_bias_preset("original")`) | `codon_bias_preset("realistic")` |
| S4 | `realistic_rfp_settings()`: `dmn_alpha_scale = 6`, learned fragment geometry + frame bias | package defaults | `do.call(simNGScoverage, c(list(...), realistic_rfp_settings()))` |
| S5 | `simGenome(length_model = "human")` | `"fixed"` | `"human"` |
| S6 | `codon_bias_options(zero_prob)` | 0 | 0.3-0.4 |
| S7 | `simIsoformGenes()`, `simIsoformReads()`, `isoform_expected()`, `isoform_mixture_alpha()` | separate layer | `simIsoformGenes(n, n_iso = "human")` |
| S8 | `simParalogFamilies()`, `paralog_totals()` | separate layer | see below |
| S9 | `codon_bias_options(tau2 = c(cds, utr5, utr3), tau2_window)` | 0 | about 2x the human tau2 (CDS 0.026-0.039, 5'UTR 0.15-0.19, 3'UTR 0.25-0.45) |
| S10 | `options(coverageSim.utr3 = ...)`: `utr3_options()`, `utr3_preset()`, `draw_utr3_features()` | off | `utr3_preset("human")` |
| S11 | `simNGScoverage(technical_artifacts = list(spike_sites, ...))` | `spike_sites = 0` | e.g. `spike_sites = 2` |
| S12 | `options(coverageSim.codon_dwell = ...)`: `codon_dwell_options()`, `codon_dwell_preset()`, `assign_dwell_groups()`, `library_dwell_profile()` | off | `codon_dwell_preset("human")` + `assign_dwell_groups()` |
| S13 | `codon_bias_options(pp_stall, pp_stall_take, pp_stall_peak, pro_run, pro_run_factor, pause_sharp, pause_flank_width, pause_flank_take)` | off | `codon_bias_preset("human")` |

Truth tables written next to the libraries: `<library>_ground_truth.tsv` (`ground_truth = TRUE`),
`<library>_spike_truth.tsv` (S11), `utr3_features.tsv` (S10), `codon_dwell_profiles.tsv` (S12).

## Simulator ports (S rows)

**S2 (11-19) shared single-codon pauses.** Real replicates share single-codon structure the simulator lacked (replicate
codon r 0.37 sim vs 0.54 real). Evidence: BC/t57/REPORT.md, BC/t63/REPORT.md (pause_sd 1.2 / pause_prob 0.45: codon r
0.52-0.53 at 600-1200 P-sites vs real 0.45-0.54; false flags vs the true alphas 4.4 / 2.8%). Setting:
`codon_bias_options(pause_sd, pause_prob)` with per-genome draws, `codon_bias_preset("t57")`. Default off. Lesson: truth
checks must use the alphas the simulator used, not the global A-site table.

**S3 (11-19) frame noise.** The simulator had 3x the real replicate frame noise (P-site clumping) and uniform gene frame
profiles; the human frame band was confirmed (BC/t64/REPORT.md). Fixed by S4; gene-level frame variation (real sd
0.06-0.15) is still not simulated (sim 0.005; docs/REALISTIC_RFP.md).

**S4 (11-19) learned geometry, concentration.** `dmn_alpha_scale = 6`, fragment geometry and by-length frame bias learned
from SRR12651187 (BC/t69/REPORT.md, BC/t74/REPORT.md): false density flags 4.0 / 2.6 -> 0.3 / 0.2%, frame false flags
4.7 / 2.9 -> 0 / 0%, replicate frame sd 0.164 -> 0.064 (real 0.052). Setting: `realistic_rfp_settings()`.

**S5 (11-19) realistic preset and human lengths.** pause 1.0 / 0.7 (`codon_bias_preset("realistic")`, BC/t75/REPORT.md:
peak 9.0 vs real 9-12), human transcript region lengths (`simGenome(length_model = "human")`, BC/t73). False flags
0.5-0.6%. Open: p95 peak, skew, clumping 4.25 vs real 2.7.

**S6 (11-19) shared codon zero mask.** Zero share 0.37 -> 0.53 at `zero_prob = 0.4` (real 0.57-0.60), frame sd 0.067,
false flags < 1% (BC/t77/REPORT.md). Cost: clumping 4.25 -> 6.6 (real 2.7). Setting: `codon_bias_options(zero_prob)`,
default 0, 0.3 milder.

**S7 (11-20) multi-isoform genes, incl. MANE-minor genes** (from E4, E10-A, GW1.1, E16-B). `simIsoformGenes()` resamples
isoform count, MANE weight and specific-CDS fraction jointly from the 17,935 real genes
(`inst/extdata/isoform_human_summary.csv`): weight of the reference < 0.5 in 20.9% of simulated genes (real 20.7%).
Truth check (BC/E17/simports/REPORT.md, scripts/s7_truth_check.R): the gated gene-fitted mixture cuts false flags in
minor genes from 34% to 1.6% (matched RNA), 6.5-7.1% with isoform TE / tissue mismatch. Not inside `simGenome()` / GTF.

**S8 (11-20) paralog families / multimappers** (from E5-B, E7-B, E8-B, GW1.3, M1-b). `simParalogFamilies(k, shares,
divergence, leak, planted, ...)` gives true member reads, the all-mapper track (a multimapper counted at every member)
and the unique track; `paralog_totals()` recovers member totals by EM on unique flanks (median 1.00 of truth vs naive
all-mapper 1.64 / 2.28 / 2.92 for k = 2 / 3 / 4; BC/E17/simports/REPORT.md). Grid level, not inside `simGenome()`.

**S9 (11-20) shared per-window variance tau2, per region** (from E7-A, E8-A, E10-C). `codon_bias_options(tau2 = c(cds =,
utr5 =, utr3 =), tau2_window = 90)`: one log-normal factor per 90-nt block shared by all libraries. The E7/A estimator
recovers 0.92-1.06 of the configured value on aligned windows, about half on windows straddling two blocks.

**S10 (11-21) 3'UTR background, ribosome-specific 3'UTR plateaus, APA step** (from E10-C, E11-B, E12-A, E12-B, A1-screen;
this port). `options(coverageSim.utr3 = utr3_preset("human"))`, see `?utr3_options`:
- plateaus: per eligible 90-nt 3'UTR window (not the first or last) with probability `plateau_prob` (human 0.045:
  3.8-5.2% of 3'UTR windows keep an excess after the tau2 band, BC/E10/C/csv/utr_tau2_pooled_Emin.csv and
  utr3_remaining_*.csv), 85% in the proximal half (real 83-90%), fold `1 + exp(N(log(2), 0.65))` (median O/E 3; real
  46-56% of excess reads at O/E-1 in (1, 2], 6-11% above 5), unframed, flat over the window (broad: 18-40 distinct 5'
  ends per window, BC/E12/A, BC/E12/B), in RFP only (RNA of the same samples: fold 1.02 vs ribo 1.61, E12-B).
- APA step: an internal poly(A) site in 54% of 3'UTRs (`apa_prob`; 1,231 MANE 3'UTRs with >= 1 internal GENCODE v44
  site, BC/E11/B/csv/tx_features.csv), uniform between 10 and 90% of the 3'UTR; downstream coverage x `apa_ratio` 0.51
  (A1 screen: median O/E 0.57 downstream vs 1.11 upstream), in all library types.
- depleted transcript end: last 25 nt x 0.07 in RFP (E12-A: 0.04-0.09x the 3'UTR mean, controls 0.07-0.13x).
Sites are drawn once per genome (seed + 5000) and written to `utr3_features.tsv`. Test (test_utr3_features.R): planted
plateau windows carry 4.0x and APA-downstream positions 0.5x the reads of ordinary positions in simulated RPF reads
(tolerance 20%); rates, positions and folds of the draws match the settings. Explainer truth check (backclaude: tau2
band and 3'UTR excess flags on the new simulation) still to run. Not simulated: the shorter plateau reads (25-27 nt 23%
vs 15% in the CDS, E12-B) and the extra duplicate-rich component of PRJNA521751 (study-specific amplitude).

**S11 (11-20) technical spikes** (from E6-B, E7-C, E16, E16-R2; carina-claude, commit 7e27f0a). `technical_artifacts =
list(spike_sites = n, spike_library_fraction = 0.25, spike_multiplicity = c(1e3, 1e5), spike_lengths = c(20, 25, 32),
spike_seed = 1)`: one read sequence at one CDS position in a fraction of the libraries, log-uniform multiplicity
1e3-1e5 (real: LSM4 20 nt, TEKTL1 25 nt, 1,000-173,000 copies; PRJNA807729 H1-2 32 nt ~19 M reads). Needs SAM/BAM
output and `fragment_mode = "simulated_rpf"`. Truth in `<library>_spike_truth.tsv`. E16-R2 later showed these are
contamination (28S rRNA, Mycoplasma 23S), not oligos; the simulated signature (one sequence, one position, few libraries)
is the same.

**S12 (11-21) study / tissue / library codon dwell, H1 context dwell** (from H1, E13-E, M1-a, M1-c, M1-d, M1-e, M5-a, M6;
this port). `options(coverageSim.codon_dwell = codon_dwell_preset("human"))` plus `assign_dwell_groups(count_table,
study, tissue_group, peaky)`, see `?codon_dwell_options`. Data: `inst/extdata/codon_dwell_human.csv` (built by
`analysis/R/build_codon_dwell_human.R`). Layers, multiplied onto the codon weights after the smoothing step:
- global context dwell: `global = "h1"` (H1 GLM, 16 offsets -10..+5 relative to the A site x 60 codon contrasts,
  `hakon_playground/agent_sync/codex_tasks/explainer_h1_codon_dwell_20261118/model.rds`) or `"atlas"` (consensus of 367
  study means, MF/m6_study_fingerprint/beta_matrix_w10_5.rds; H1 and atlas agree, per-offset r 0.62-0.98). The
  preset uses offsets != 0, since the bundled A-site table already has the A-site effect.
- study deviation per codon x offset and a smaller library deviation, scaled so that true library profiles correlate
  `within_r = 0.88` inside a study and `between_r = 0.37` between studies (M6, MF/m6_study_fingerprint/REPORT.md: same
  study 0.85-0.89, other study 0.30-0.35, study explains 0.81 of library dwell variance); per-offset sds from the atlas.
  Check (test_library_dwell.R, 160 libraries / 40 studies): within 0.88, between 0.38.
- tissue group "striated" (M1-a, MF/dwell_by_group/muscle_detail.csv: E/P/A deltas, Pro A-site +0.10..+0.24, CTT/CTA
  slower, CAG/CAC faster; M1-c: predicts muscle-gene codon profiles out of sample); own groups via `tissue_deltas`.
- peaky libraries: codon weights ^ 1.3 and Dirichlet concentration / 1.6 (M1-d: striated libraries have 1.5x the
  top-10-codon share, 0.21-0.30 vs 0.14-0.18, and 1.59x window extra-Poisson variance). Check: top-10 share 0.137 ->
  0.20 (1.47-1.51x) at 3,000-15,000 reads per gene.
Profiles are written to `codon_dwell_profiles.tsv`. Test: a planted +log(8) A-site GAA delta raises the GAA read odds 4-12x
in the planted library only. Explainer truth check (study-mixture dwell of v1.6 on the new simulation) still to run.

**S13 (11-21) sharp pauses, Pro-Pro stalls, Pro-run depletion, autocorrelated background** (from M15, M16, M16-acf;
this port, requested by Håkon). Real data (MF/m16_repeats/acf_real.R and raw_pro.R; 7.8 M codons of 14,310
well-covered MANE genes, merged unique track, v1.8b explainer expectation, median ~1,200 P-sites per codon):
- background per-codon ACF of gene-centred log(O+1) at lags 1-5: 0.58 / 0.42 / 0.37 / 0.33 / 0.28 (lag 10 0.15); the
  codon dwell model alone 0.45 / 0.32 / 0.31 / 0.26 / 0.17; residual log O/E 0.64 at lag 1, 0.47 at lag 2, 0.22 at lag
  10 (smooth red noise, no periodicity, M15).
- 42,195 Pro-Pro stalls (Pro at E and P, non-Pro A; offset 0 = P-site codon): raw log density at offsets -3..+3
  -0.28 / -0.31 / -0.16 / +0.20 / -0.12 / -0.40 / -0.36, residual vs the model -0.13 / -0.19 / -0.14 / +0.10 / -0.09 /
  -0.06 / -0.04: a sharp peak whose neighbours are depleted relative to the model.
- Pro runs: raw log density by Pro count in P-site +-2: 0 -0.28, 1 -0.32, 2 -0.44, 3 -0.64, 4-5 -1.18; the additive
  model predicts -0.09 .. -0.20 (observed minus model relative to 0 Pro: 3 -0.27, 4-5 -0.80). PPP at E-P-A -0.99.

The original simulator multiplies a pause onto one codon after the smoothing (S2): a peak without depleted flanks, and
with the "realistic" preset (fixed sharpening, pause_sd 1.0 on 70% of codons) almost no autocorrelation beyond lag 1
(0.39 / 0.03 / 0 ...; the original per-gene sharpening has the opposite problem, 0.84 / 0.72 / 0.61 / 0.52 / 0.44).
Settings (`?codon_bias_options`), all after the smoothing step:
- `pp_stall = TRUE`: the codon whose P site is a stall takes `pp_stall_take` of each neighbour's weight (offsets
  -1, +1; a vector of length 2w for -w..-1, +1..+w) and is then multiplied by `pp_stall_peak`. Defaults c(0.20, 0.10)
  and 0.9, calibrated on simulated reads: the direct translation of the residuals, 1 - exp(c(-0.14, -0.09)) =
  c(0.13, 0.09), gives only -0.08 / -0.08 on the flanks because the simulator's RNase blur spreads part of the loss;
  the real peak gains less than the flanks lose, hence the peak factor < 1. Applied after the S12 dwell layers, so the
  stall takes relative to the model weights. The decoding-site frame is respected (`a_shift`: A- or P-site weights).
- `pro_run = TRUE`: the codon whose P site has nP Pro codons in P-site +-2 gets `pro_run_factor[nP]` (default nP 3
  0.766 = exp(-0.266), 4-5 0.45 = exp(-0.799); add `"2" = 0.90` for the milder nP = 2).
- `pause_sharp = TRUE`: the S2 pauses (m > 1) take their extra weight from their `pause_flank_width` neighbours (at most
  `pause_flank_take` of the flank total), conserving the local total; dips stay multiplicative. Not in the preset: with
  pause_prob 0.5-0.9 the cap shrinks the codon-level variance, so the smooth component dominates and lags 2-10 rise
  (0.68 / 0.59 / 0.54 / 0.48 / 0.43 at pause_sd 1.1, pause_prob 0.9); use it with sparse pauses.
- `codon_bias_preset("human")`: the original per-gene sharpening + smoothing (the smooth background), per-genome draws,
  pause_sd 0.9 / pause_prob 0.8, CDS tau2 0.06 (S9, about 2x human), `pp_stall`, `pro_run`. A larger tau2 (0.35) also
  fits the ACF but adds about 10x the real window-level variance.

Test (test_sharp_pauses.R; 1,000 genes x 400 random sense codons, 20 libraries of 1,500 reads per gene merged =
30,000 reads per gene): simulated ACF at lags 1-5 0.64 / 0.42 / 0.35 / 0.31 / 0.26 (lag 10 0.10; tolerance +-0.08 of
real); residual (S13 on vs off) at Pro-Pro stall offsets -1 / 0 / +1 -0.139 / +0.097 / -0.087 (1,558 stalls, se 0.01;
real -0.14 / +0.10 / -0.09, tolerance 50%), -2 -0.11 (real -0.19, from the Pro-run layer: windows at -2..+1 hold both
stall Pro); Pro-run on minus off by Pro count: 3 -0.24, 4-5 -0.79 (target -0.27 / -0.80). At 3,000 reads per gene the
ACF is lower (0.52 / 0.34 / 0.29 / 0.24 / 0.21; Poisson noise), the stall residual similar (-0.11 / +0.10 / -0.08).
Unit tests: local conservation of sharp pauses, flank takes and P/A-site indexing of the stall kernel, Pro counts, and
that `simNGScoverage()` hands the codons to the layers. Not simulated: the real raw stall profile in absolute terms
(simulated raw at -1 / 0 / +1 -0.40 / -0.04 / -0.15 vs real -0.16 / +0.20 / -0.12: the bundled A-site table and the
random sequences differ from the human dwell around Pro), the residual at -3 / +2 / +3, Gly runs (M16 fits a Pro x Gly
table; only Pro is ported). Explainer truth check (v1.8b on the new simulation) still to run.

## Explainer iterations (E rows)

**E1 (11-18) baseline explainer.** Windows vs ORF scanner, between-half noise band, mappability mask; CDS u/n/k 46/23/31%
(BC/t42/REPORT.md; caller benchmark on coverageSim BC/t40, t47). Method, nothing to simulate.

**E2 (11-18) v2 terms.** H1 codon-dwell baseline, isoform dip flag, pause allowance, frame floor, NTE periodicity, culprit
libraries (13 lengths reshifted): CDS understood 46 -> 69% (MF/region_explainer/ITERATION_v2.md). Simulator: the H1
dwell is S12 (`global = "h1"`); pauses S2; frame S4; isoforms S7. The v2 ORF scanner invented ORFs on sim nulls (v1 kept).

**E3 (11-18) v3.** Per-gene frame profiles, replicated pause clusters, depth-specific band: CDS deep 69.4 -> 72.2%
(BC/t56/REPORT.md). Simulator: replicated pauses = S2; per-gene frame profiles are not simulated (gene frame sd 0.005 vs
real 0.06-0.15, docs/REALISTIC_RFP.md).

**E3b (11-19) disome support for bumps.** Bumps are not collision-driven (chance level; BC/t66/REPORT.md). Rejected:
nothing to simulate.

**E4 (11-19) isoform mixture from RNA (rule 3).** No net gain (NU 11.2 -> 12.0-12.4%; BC/t65/REPORT.md). Rejected as a
rule; the underlying biology (MANE-minor genes) became S7.

**E5 (11-19) inhibitor profiles, all-mapper families, truth check on the realistic sim.** Family-level status from the
all-mapper track (CDS deep 72.2 -> 87.0% incl. family level); inhibitor profiles rejected; v3 fails the realistic sim
(false density 15.4 / 9.8%) because the human dwell shape differs from the sim codon table (BC/E5/REPORT.md). Simulator:
families = S8; the dwell mismatch led to truth checks against the simulator's own alphas (S2) and to S12.

**E6 (11-19) self-fitted dwell, library-artefact flag, family test truth check.** Sim expectation carries shared per-site
variance (sd 0.21 log) that no codon model explains -> tau2 (S9); artefact flag (<= 3 libraries hold >= 70% of a window
excess) -> S11 spikes; family test blind to abundance (BC/E6/REPORT.md).

**E7 (11-19) tau2 band, family-total test, culprit QC.** First passing truth check: sim false density 15.4 / 9.8 ->
3.9 / 2.7%; CDS deep NU 11.2 -> 8.2%; MCF7 PRJNA807729 H1-2 = one 32-nt read sequence at one P-site in 3 replicates
(BC/E7/REPORT.md). Simulator: tau2 = S9 (`codon_bias_options(tau2)`); the single-sequence site = S11; family totals = S8.

**E8 (11-19) disjoint tau2, RNA reference for families, frozen settings on 1,400 new genes.** Gain holds on new genes
(deep 88.9% understood); null-sim family abundance test 0.0% false flags (BC/E8/REPORT.md). Simulator: S9 values (CDS
0.026 deep / 0.039 mid), S8.

**E9 (11-19) ER-targeting term, gene class, pseudogene leakage.** ER term rejected (OR 1.05, p 0.82); the pilot gap is
gene class + low mappability; leakage real for 3 of 20 pseudogene pairs (BC/E9/REPORT.md). Simulator: mappability /
leakage = S8 (`leak`); ER term not simulated (rejected; see M2-a).

**E10 (11-19) isoform table, shared-noise decomposition, UTR tau2.** 3'UTR NU deep 17.5 -> 6.0%, 5'UTR 5.3 -> 0.9% with UTR
tau2; shared noise = real cross-study signal not explained by sequence features (BC/E10/REPORT.md, E10/C/REPORT_C.md).
Simulator: UTR tau2 = S9 (`tau2 = c(utr5 = 0.15-0.19, utr3 = 0.25-0.45)`); the sim UTR was a Poisson background
(sim tau2 0.004), now structured by S10; isoforms S7.

**E11 (11-19) protocol attributes, 3'UTR excess.** No protocol choice explains the shared excess (footprint quality only,
rho ~0.45); the unframed 3'UTR excess: 37% follows the RNA shape (APA / isoform-like), 53% local hot spots, near annotated
poly(A) sites 77% vs 38% (BC/E11/REPORT.md, E11/B/REPORT_B.md). Simulator: S10 (plateaus + APA step); no protocol term
(not supported).

**E12-A (11-19) poly(A)/PABP complex footprints.** Not supported (local-null enrichment 0.94x, PAS 0.99x, frame-0 0.333,
no 25-27 nt excess; present in 36 / 40 studies; BC/E12/A/REPORT_A.md). Simulator: the unframed, broad, multi-study
plateau of S10; the transcript-end depletion (last 25 nt 0.04-0.09x) is S10 `end_depletion`.

**E12-B (11-19) library-level check of the 3'UTR plateaus.** Not a two-study artefact (71% of windows stay without the
two dominant studies), not mapping / rRNA / low complexity, RNA of the same samples has no plateau (BC/E12/B/REPORT_B.md).
Simulator: S10 plateaus in `ribo_libtypes = "RFP"` only. Not simulated: shorter plateau reads, study-specific amplitude.

**E12-C (11-19) MGP deep dive.** Reproducible in-frame peak at codons 49-57 (R-E-R, Gla domain start, O/E 2.9) ending at
exon junction 3/4 in 3 tissue studies + 2 brain lines, absent in 6 others; mirrored dip (BC/E12/C/REPORT_C.md, follow-ups
E13-A/B/D). Not simulated as such: cause open (TruSeq / ART-seq prep, E13-D); the study-specific part is covered in kind
by S12 study dwell, single-gene pauses by S2.

**E13-E (11-20) A-site acidic effect.** H1 already captures it: A-coefficient vs observed density per amino acid r 0.98;
acidic density 1.39x observed vs 1.35x H1; cross-fit correction moves 0.31% of bump reads (BC/E13/E/REPORT_E.md).
Simulator: already in the bundled A-site `seq_bias` tables (median of 10 libraries: D 2.0x, E 1.47x of the median amino
acid; codon table GAC 2.4x, GAT 1.9x); the H1 betas ship with S12 (`global = "h1"`, `global_offsets = -10:5` to use the
H1 A-site betas instead, e.g. with a flat `seq_bias`). Default: the table effect is on, H1 off.

**E13-A/B/C/D, E14-F/G, E15-H (11-19/20) leading to GW1.1.** MGP literature, exon-junction test, dip-vs-RNA calibration,
MGP protocol (TruSeq 4/5 peak studies vs 0/6), top-50 NU genes, RNA-mixture expectation, v1.1 five-status split
(reports in BC/E13-E15). Simulator: isoform rule and mixture = S7; tissue / protocol effects = S12.

**E16 (11-20) top-13 NU genes.** Artefacts (LSM4, TEKTL1 single sequences), isoforms (BSG, VCAN, SCP2, TTN, CD74),
tissue-protocol (MYH7, CKM, NEB, MYBPC1), biology candidates (CALR, TMBIM6) (BC/E16/genes12/REPORT.md). Simulator:
spikes S11, isoforms S7, tissue / study dwell S12 (striated group, study deviations, peaky libraries).

**E16-B (11-20) dominant isoform as reference.** Worse than flat MANE (NU 7.00 vs 6.42%); mixture kept
(BC/E16/domiso/REPORT.md). Simulator: S7 (`ref_minor_frac`, mixtures).

**E16-R2 (11-20) origin of spikes.** LSM4 = 28S rRNA fragments (28S/18S missing from the pipeline rRNA.fasta), TEKTL1 =
Mycoplasma 23S rRNA; 23 confirmed spike windows; 197 multi-sequence single-codon peaks = pause candidates
(BC/E16/spikes/REPORT.md). Simulator: S11 reproduces the signature (one sequence, one site, a few libraries, high
multiplicity); the contamination source itself (reads from a missing reference) is not modelled.

**E19-R1 (11-20) tissue-matched RNA for MANE-minor genes.** No genome gain from tissue matching; isoform-mixture shape with
ridge weights helps (NU 5.69 -> 5.59) (BC/E19/R1/REPORT.md). Simulator: S7 `rna_weight_sd` (RNA groups with distinct
true weights) and `te_sd`.

**A1-screen (11-20) new data screens.** Windows downstream of an internal GENCODE v44 poly(A) site: median O/E 0.57 vs
1.11 upstream; a plain RNA-shaped 3'UTR expectation fits worse; 3'UTR NU reads overlap annotated ORFs 2x (ledger).
Simulator: S10 APA step (`apa_prob = 0.54`, `apa_ratio = 0.51`); 3'UTR ORFs (dORFs) are the existing
`simGenome()` architecture / translon features, not changed here.

## Genome map runs (GW rows)

The GW rows run the frozen explainer on all 17,935 MANE genes; they add rules, not new biology, so their simulator
counterparts are those of the rules they use.

**GW1 (11-19) map v1.** CDS reads understood / family / artefact / NU / unknown 87.7 / 5.7 / 0.1 / 6.4 / 0.1; top 1% of
genes hold 45% of NU (BC/E11/C/REPORT_C.md). Simulator: inherited from E7 / E8 (S8, S9).

**GW1.1 / GW1.1-final (11-20) map v1.1.** Isoform dip rule, gene-fitted RNA mixture for MANE-minor genes, UTR tau2: CDS NU
6.42 -> 5.69%, 3'UTR NU 17.7 -> 9.2%; 3,708 MANE-minor genes hold 38% of genome NU (BC/E15/H, BC/E16/v11/REPORT.md).
Simulator: S7, S9, S10.

**GW1.2 (11-20) map v1.2.** Isoform-mixture shape with ridge weights, spike genes, contamination status: CDS NU 5.69 ->
5.38% (BC/E22/v12/REPORT.md). Simulator: S7, S11.

**GW1.3 (11-20) exact 28-mer mappability.** Exact routing: CDS understood 90.29 / family 4.03 / NU 5.41%; conservative
(<= 2 mismatch) routing 76.53 / 17.82 / 5.39%; unique > all at edge read lengths = P-shift offset mismatch
(BC/E23/v13/REPORT.md). Simulator: unique vs all-mapper tracks = S8. Offset mismatches per read length are not simulated
(the fragment geometry of S4 uses one learned offset table).

**GW1.4 (11-20) striated dwell layer.** CDS NU 5.41 -> 5.41%: dwell is real (M1-c) but acts at codon scale
(MF/v14). Simulator: S12 tissue group "striated".

**GW1.6a/b/c, GW1.6 adopted (11-20/21) study-mixture dwell.** Gated log-space study-composition correction for
study-dominated genes: CDS NU 5.426 -> 5.295% (-38 M reads); study-dominated genes |log O/E| 0.346 -> 0.309
(MF/v16; ledger). Simulator: S12 study deviations (every library of a study shares its codon x offset deviation), so a
gene expressed in one dominant study gets that study's codon profile.

## Mechanism iterations (M rows)

**M1-a (11-20) tissue dwell, formal test.** Striated muscle + heart (5 studies) replicate their study deltas (r 0.21,
p 0.023); tissue explains ~11% of between-study variance; Pro A-site +0.10..+0.24, CTT/CTA slower, CAG/CAC faster
(MF/dwell_by_group/tissue_test.csv, muscle_detail.csv, LIT_MUSCLE_TRANSLATION.md). Simulator: S12 `tissue_group =
"striated"` (deltas shipped in codon_dwell_human.csv, column striated_delta).

**M1-b (11-20) NEB near-identical windows.** Windows with lowmap > 0.9 fail 22% (vs 8%); hybrid routing gives no gain;
these windows also fail on the all-mapper track (ledger). Simulator: S8 (`divergence` near 0). Not adopted.

**M1-c (11-20) striated dwell out of sample.** Striated betas raise the explained codon-profile deviance of CKM 0.112 ->
0.195, MYH7 0.072 -> 0.159, DES 0.148 -> 0.322 and beat 30 random study sets for 7 / 9 genes
(MF/dwell_by_group/m1c_*.csv). Simulator: S12 striated group.

**M1-d (11-20) why muscle windows fail.** Peakiness (top-3 codons 30-90% of window reads), a replication blind spot of
study halves, and peaky libraries: 1.5x top-10-codon share, 1.59x window extra-Poisson variance
(MF/dwell_by_group/peakiness_housekeeping.csv, phi_housekeeping.csv). Simulator: S12 `peaky = TRUE` libraries
(`peaky_power = 1.3`, `peaky_var_ratio = 1.6`).

**M1-e (11-20) library-dispersion noise multiplier.** Real but not the explanation (symmetric phi: CDS NU 5.41 -> 5.49%;
MF/v15). Simulator: S12 peaky libraries give per-library dispersion; per-gene phi is not a separate setting.

**M2-a (11-20) ER targeting.** Signal-peptide genes are NU more often along the whole CDS (10-15% vs 6-8%), density over
the signal itself overpredicted (log O/E -0.19) (MF/m2_er_targeting/sp_profile.csv). Not simulated: a gene-class effect
of unknown mechanism, not adopted by the explainer; the simulated genome has no signal peptides.

**M3 / M3-b (11-20) library-split replication; cross-study replication of muscle profiles.** Study halves cover only 31%
of reads; muscle window patterns are mostly study-specific (study-vs-study r 0.14-0.60) (MF/m3_replication/REPORT.md).
Simulator: S12 study deviations make window patterns study-specific in the simulation too.

**M4 (11-20) study-level artefact test.** 100.5 M NU reads (6.4% of genome NU) are study artefacts in the strict sense,
87.6 M from one skeletal-muscle study (MF/m4_study_artefact/REPORT.md). Simulator: S12 (a study with a dominant share of a
gene's reads and its own dwell deviation); study composition per gene comes from the count table design.

**M5-a (11-20) study-specific dwell.** The dominant study's own betas explain its CKM profile 0.41 vs consensus 0.09; study
betas 2x more extreme (sd 0.20 vs 0.10) (MF/dwell_by_group/motif_test.csv). Simulator: S12 study deviations (an extreme
study can be given via `study_scale` or a custom `tissue_deltas` group).

**M6 (11-20) study dwell fingerprint.** Study explains 0.81 of library dwell-beta variance; cell line / tissue / inhibitor
add +0.027; same study pair r 0.85-0.89, other study 0.30-0.35; technical covariates predict 0.16
(MF/m6_study_fingerprint/REPORT.md). Simulator: S12 `study_dwell = TRUE`, `between_r` / `within_r`, per-offset sds of the
atlas.

**M7 (11-20) codon demand / tRNA depletion.** Muscle depletion not supported (cor -0.25); weak, partly circular general
coupling (0.17-0.24) (MF/dwell_by_group/demand_test.csv). Not simulated: no supported mechanism.

**M8 (11-20) nascent-chain terms.** Coiled-coil onset dip ~4%, exit-tunnel charge slope +0.0032 log O/E per unit charge,
R2 0.004 (MF/m8_nascent/nascent.R; MF/lit_residual_nu/REVIEW.md). Not simulated: small, not adopted; the simulated genome
has no protein features.

**M9 (11-21) co-translational assembly, direct test.** Mechanism confirmed in the Bertolini disome data, but no dip in the
merged monosome windows; a fitted distance term worsens held-out fit (MF/m8_nascent/disome.R). Not simulated: not adopted.

**M15 (11-21) UniProt motifs and periodicity.** Motifs / sites / PTMs |delta| <= 0.05; poly-Pro -0.39, poly-Gly -0.35
(additive dwell saturates on repeats); the residual is red noise without 3.6 / 2 / 7 codon periods
(MF/m15_uniprot_motifs/REPORT.md). Simulator: the red-noise background and the Pro runs are S13; no motif term.

**M16 / M16-acf (11-21) Pro / Gly repeat saturation, per-codon autocorrelation.** Pro x Gly run factor cross-fitted
(per-codon R2 +0.65%, not adopted for the window map); raw coverage falls with Pro count; Pro-Pro stalls are sharp
peaks with depleted flanks; real ACF 0.58 / 0.42 / 0.37 / 0.33 / 0.28 (MF/m16_repeats/). Simulator: S13
(`codon_bias_preset("human")`).

## Simulating realistic human Ribo-seq

The ports combine. A study design with three studies, one of them striated muscle with peaky libraries:

```r
library(coverageSim)
options(coverageSim.codon_bias  = modifyList(codon_bias_preset("realistic"),            # S2, S5 (or "human": S13)
                                             list(zero_prob = 0.3,                       # S6
                                                  tau2 = c(cds = 0.06, utr5 = 0.3, utr3 = 0.6))),  # S9 (~2x human)
        coverageSim.utr3        = utr3_preset("human"),                                  # S10
        coverageSim.codon_dwell = codon_dwell_preset("human"))                           # S12
genome <- simGenome(n = 500, length_model = "human", out_dir = tempfile())                # S5
counts <- simCountTables(ORFik::loadRegion(genome["txdb"], "cds"), libtypes = "RFP",
                         conditions = c("A", "B", "C"), replicates = 2)
regions <- simCountTablesRegions(counts, regionsToSample = c("leader", "cds", "trailer"))
regions <- assign_dwell_groups(regions, study = rep(c("S1", "S2", "S3"), each = 2),     # S12
                               tissue_group = rep(c(NA, NA, "striated"), each = 2),
                               peaky = rep(c(FALSE, FALSE, TRUE), each = 2))
exp <- do.call(simNGScoverage, c(list(simGenome = genome, count_table = regions,
                                      libFormats = list(RFP = "bam"), ground_truth = TRUE,
                                      technical_artifacts = list(spike_sites = 2)),   # S11
                                 realistic_rfp_settings()))                            # S4
```

Isoform mixtures (S7) and paralog families (S8) are separate grid-level simulators (`simIsoformGenes()` +
`simIsoformReads()`, `simParalogFamilies()`); use them for the isoform and multimapper truth checks.

For the codon-level structure (smooth autocorrelated background, Pro-Pro stalls with depleted flanks, Pro-run
depletion) use `codon_bias_preset("human")` (S13) in place of "realistic" and add the UTR tau2 values of S9; S13 has not
yet been checked against the S4-S6 clumping and peak metrics, which were tuned on "realistic".

Known gaps (real vs simulated): per-gene frame profiles (S3), clumping (S5-S6), shorter 3'UTR plateau reads and
study-specific plateau amplitude (S10), read-length-specific P-site offset errors (GW1.3), ER / nascent-chain classes
(M2, M8, M9), contamination from missing references (E16-R2, signature only via S11).
