# coverageSim speed core pass 2

- B3k, 1 thread: 287.154 → 245.787 s; B20k, 1 thread: 1515.129 → 1552.559 s.
- B3k, 8 threads: 109.554 s; B20k, 8 threads: stopped at 50.78 GiB PSS; no total.
- B20k target <300 s is not achieved by the changes implemented here; coverage,
  fragment allocation, sequence matching and flat transcript models are unchanged.

| Genes | Threads | Genome before → after | Counts | Regions | Coverage | Total |
|---:|---:|---:|---:|---:|---:|---:|
| 3,000 | 1 | 57.906 → 35.351 | 0.390 → 0.299 | 1.591 → 0.961 | 227.267 → 209.176 | 287.154 → 245.787 |
| 20,000 | 1 | 218.323 → 171.063 | 0.469 → 0.800 | 7.006 → 10.812 | 1289.331 → 1369.885 | 1515.129 → 1552.559 |
| 3,000 | 8 | — → 33.878 | — → 0.272 | — → 0.913 | — → 74.491 | — → 109.554 |
| 20,000 | 8 | — → 127.006 | — → 0.454 | — → 6.469 | stopped | unavailable |

- B20k genome Rprof inclusive: writer 91.44 → 52.80 s; placement 17.02 → 10.04 s.
- Single runs; different system load may affect timings. Benchmark driver:
  `analysis/benchmarks/speed_core.R`; pass 1 comparison: `../coveragesim_speed_core_20261118/bench/compiled{3,20}`.
- Compiled writer preserves random draw order; `options(coverageSim.engine = "R")`
  selects the original writer. Twenty overlapping-frame seeds and an empty feature
  table matched in output and final RNG state; focused test: 41 passes.
- B3k FASTA, feature table, GTF, and all 24 coverage, truth and dropped-read files
  match pass 1 SHA-256. B20k FASTA, feature table, GTF and all 24 coverage, truth and dropped-read
  files also match. All 24 B3k outputs match between 1 and 8 threads. No output
  changes are intended.
- Full test suite: no test failures; 14 warnings, 10 skips. Three R finalizer error
  messages were printed outside the test-failure summary; see `../suite.log`.
- Not implemented: coverage C++ pass, fragment/allocation and sequence-bias kernels,
  flat model mapping, and placement rewrite. Five-seed distribution comparisons
  are not done. The complete requested performance target remains open.
- Eight-thread B20k stopped at 50.78 GiB proportional memory to keep headroom
  under the requested 60 GB reservation; log: `../memory_watch.log`.
