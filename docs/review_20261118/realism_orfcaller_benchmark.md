# coverageSim benchmark (branch review/carina-dev-claude-20261118 @ 8e284eb, copy in ./repo, unmodified)
Sim: simGenome(n=300, max_uorfs=3, uorfs_can_overlap_cds=TRUE), 2 RFP reps, simulated_rpf, ground_truth=TRUE, defaults (seed 1).
Default depth is only ~10k reads/lib (0.36 P-sites/codon), so a deep run on the same genome was added (interceptMean=10, 0.85 M).
Real: PRJNA556704 HEK293 CHX WT reps 1-2, P-shifted ofst, 10,823 longest tx with >=100 CDS P-sites. Sim P-site = truth
A-site offset - 3; ORFik detectRibosomeShifts on sim finds the same offsets (27:-11, 28:-12, 29:-12).

## Part A: realism (partA_metrics.csv, partA_figure.png)
| metric | real | sim default | sim deep |
|---|---|---|---|
| read length (mode; central 90 %) | 33; 31-34 | 29; 27-29 | 29; 27-29 |
| frame 0/1/2 | .66/.14/.19 | .64/.12/.24 | .65/.11/.24 |
| start peak / stop peak (x CDS mean) | 6.8 / 2.3 | 7.2 / 1.8 | 6.3 / 1.2 |
| ramp (codons 2-20 / CDS interior) | 0.70 | 1.06 | 0.95 |
| leader/CDS, trailer/CDS density | 0.17, 0.020 | 0.41, 0.36 | 0.43, 0.38 |
| share of reads leader/CDS/trailer | 2/96/2 % | 14/76/9 % | 15/75/10 % |
| Gini, zero codons @0.25 reads/codon | .87, 83 % | .85, 81 % | .84, 81 % |
| Gini, CV, zero codons @1 read/codon | .76, 1.82, 60 % | .69, 1.45, 53 % (20 genes) | .70, 1.54, 54 % |
| replicate r, codons @1 / @0.25 reads/codon | 0.70 / 0.40 | 0.27 / 0.13 | 0.29 / 0.12 |
| replicate r, genes (log) | 0.998 | 0.76 | 0.95 |
Realistic: frame shares (within 2-5 points), start peak, per-gene roughness at matched depth (Gini .84 vs .87).
Not realistic:
- UTRs: trailer 18x and leader 2.5x too dense (default region_proportion leader .05, trailer .10, uorf .10).
- Roughness is noise, not signal: codon-level replicate r is 0.28 vs 0.70 real at equal depth (no reproducible pausing).
- No 5' ramp, weak or no stop peak, one narrow 27-29 length set (the real library is 31-34; learn it per library).
- Default depth (~33 reads/gene) leaves most uORFs with 0-2 reads.

## Part B: ORF callers (sim300*/partB_auc.csv, partB_orf_features.csv)
findUORFs(ATG, longestORF=FALSE, cds) recovers all 888 true uORFs plus 439 ATG-stop decoys. True: 458 leader-only,
178 uoORF, 252 overlapping another uORF. Decoys: 57 leader-only, 84 uoORF, 124 overlapping a true uORF, and 174 in-frame
start-variants of a true uORF (same stop). AUC against decoys of the same class:
| class | depth | ORFscore | reads/nt | reads | frac codons covered |
|---|---|---|---|---|---|
| leader-only | default / deep | .66 / .85 | .63 / .67 | .62 / .68 | .63 / .49 |
| uoORF | default / deep | .59 / .71 | .47 / .46 | .37 / .35 | .49 / .56 |
| overlapping uORFs | default / deep | .66 / .84 | .62 / .65 | .61 / .64 | .62 / .58 |
| true vs nested in-frame decoys | deep | .60 | .54 | .48 | .50 |
- uoORF: coverage features have AUC < 0.5 because out-of-frame decoys running into the CDS collect CDS reads (median
  464 vs 199). Only frame (ORFscore) works, which is the intended hard case.
- In-frame start-variants cannot be separated (no initiation peak on uORFs); score them as "wrong start", not FP.
- Overlap frames (labelled run, 178 uoORF-CDS overlaps, 76k reads): 90 % of reads come from the CDS, 10 % from the uORF.
  Pooled, 65 % are in CDS frame 0 and 20 % in uORF frame 0: a CDS-dominated mix. Each source keeps its own frame
  (CDS 70 %, uORF 59 % in their own frame 0).
- Coupling: only compositional. uORF vs CDS counts rho 0.99 (shared gene total); uORF share vs CDS share -0.63
  (multinomial split). The CDS share is 0.750 with a uoORF and 0.752 without. The uORF share does not depend on the
  number of uORFs (rho -0.03) or their length (0.11).

## Part C: run time (1 sample, 1 core; peak RSS of the R process)
| n | wall | simGenome | simNGScoverage | peak RSS |
|---|---|---|---|---|
| 300 | 126 s | 74 s | 30 s | 2.17 GB |
| 3000 | 589 s | 422 s | 142 s | 2.53 GB |
simGenome dominates and scales about linearly. n=300 with 2 reps, deep: 158 s, 2.34 GB.

## Missing for RDG / overlapping-translon tests
1. No regulatory coupling: no uORF-dependent CDS repression, leaky scanning or reinitiation. The region split is a fixed
   multinomial, and uORF reads are shared uniformly regardless of length, Kozak context or position.
2. ATG starts only; no non-AUG uORFs, N-terminal extensions, internal out-of-frame ORFs, dORFs or isoforms (1 tx/gene).
3. No initiation or termination peaks on uORFs; no queuing or collisions.
4. uORF number barely varies (288/300 genes have 3, 12 have 2), so there are no 0- or 1-uORF controls.
5. Bug or gap: by default the ground-truth sampling_group is a per-position row index, not a region. Region labels
   (uorf:k/cds:k) appear only when a fragment bias is active. I used a neutral 5' end bias (all weights 1).
6. No explicit negative (untranslated decoy) set is exported, and nested in-frame start variants need a truth rule.

## Verdict
Usable for frame-based ORF-caller benchmarks (ORFscore AUC .84-.85 for leader-only and overlapping uORFs) once depth is
raised. The default depth gives about 0.6 for every feature. uoORFs reproduce the realistic failure of coverage-based
callers. UTR densities and pause reproducibility are unrealistic, and RDG coupling is absent (items 1-4), so it cannot
yet test decision-graph inference. Fit region_proportion, geometry and end bias from a real BAM before benchmarking.
Files: scripts/*.R, partA_*.csv, partA_figure.png, sim300{,deep,lab}/partB_*.csv, time_*.json; extraDistr in ./Rlib only.
