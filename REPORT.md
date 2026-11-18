# coverageSim speed core

- B3k: 497.0 → 287.2 s total (single runs, 42% less); B20k: 1,919.383 → 1,515.129 s (21% less).
- Step elapsed seconds (baseline → optimized):

| Genes | Genome | Counts | Regions | Coverage | Total |
|---:|---:|---:|---:|---:|---:|
| 3,000 | 61.018 → 57.906 | 0.547 → 0.390 | 2.135 → 1.591 | 433.234 → 227.267 | 496.933 → 287.154 |
| 20,000 | 283.346 → 218.323 | 0.999 → 0.469 | 13.687 → 7.006 | 1621.351 → 1289.331 | 1919.383 → 1515.129 |

- B3k baseline coverage Rprof inclusive seconds: nucleotide coverage 209.14,
  sequence table 101.70, RNase extension 76.80, sequence bias 63.48,
  fragment generation 52.70, transcript models 42.14. Optimized: 128.00,
  29.00, 6.84, 32.56, 26.70, 28.16 respectively. These overlap.
- B20k coverage Rprof inclusive seconds (baseline → optimized): nucleotide
  coverage 789.22 → 790.82, fragment generation 216.38 → 224.02, sequence
  bias 212.34 → 213.48, models 203.30 → 110.84, sequence table 329.14 →
  87.16, RNase extension 269.28 → 27.30. Inclusive times overlap.
- B3k genome Rprof: feature sequence writer 12.84 s, feature placement 4.44 s;
  B20k baseline: 97.00 and 35.92 s; preset construction within placement
  19.10 s. Optimized B20k writer 91.44, placement 17.02, presets 0.50 s.
- B3k baseline counts/regions Rprof self leaders were garbage collection and
  data.table conversion; B20k region conversion had 9.00 s inclusive.
  Profiles include `system.time()` setup GC; elapsed times exclude that setup.
- One thread, seeds 6618, leader 50–400, CDS 300–3000 in multiples of 3,
  trailer 100–1500, all architectures, 3 CDS exons, 2 conditions × 3
  replicates for each of RFP/RNA, simulated RPF, OFST, ground truth, validation
  disabled. Counts use default interceptMean=4. Ground truth is RFP-only.
- Measurements overlap other jobs; single samples, no confidence intervals.
  Rprof at 20 ms; profiles and full logs are in `../bench/`, benchmark driver
  `analysis/benchmarks/speed_core.R`; command sets `R_LIBS` to task Rlib.
- Exact: all 24 corresponding OFST, truth and dropped-read files at each scale
  matched SHA-256; B3k/B20k FASTA and feature tables matched. This does not establish
  equivalence for all settings. Changes intentionally preserve random draw
  order; no output differences are intended against 130374f.
- Batched RNase extension preserves unusual factor, missing, keyed and indexed
  groups via the original path. Repeated architecture presets are constructed
  once. The C++17/Rcpp matcher preserves position order and random draws.
- New focused tests: RNase 12, preset 4, compiled matcher 4 expectations.
  Full `tests/testthat`: 1,599 passed, 0 failed, 0 errors, 14 warnings,
  10 skipped (suite log: `../bench/suite_counts.log`).
- Incomplete: the requested flat exon table / vectorized model rearchitecture,
  and optimizations of sequence-bias matching, fragment allocation and per-gene
  coverage were not implemented; measured B20k took 25.3 minutes.
- Next candidates: flat transcript models (at most ~12% B3k coverage if entirely
  eliminated), sequence-bias processing (~14% inclusive), and nucleotide
  coverage generation (~56% inclusive); actual gains will be smaller and overlap.
