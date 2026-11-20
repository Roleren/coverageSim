# coverageSim review

Reviewed branch review/carina-dev-claude-20261118, commit 8e284eb94b0c.

## Findings (severity ranked)

1. **High — silently lost reads when selecting regions.** `repo/R/Sim_count_tables.R:250` samples all default proportions, but lines 255–270 retain only `regionsToSample`; `Sim_reads_from_counts.R:333` treats the retained sum as the budget. `repro.R`, seed 2026, sets 3 genes to 1,000 RFP reads each and uses the README's leader/CDS/trailer selection: 3,000 gene reads become 2,682 region reads, with no dropped-read record. Renormalize to selected categories or explicitly report the omitted counts.
2. **Medium — a genome with no uORFs cannot use the default region table.** `repo/R/Sim_reads_from_counts_helpers.R:415` unconditionally reads `true_uORFs.rds`, while `repo/R/Sim_genome_helpers.R:452` does not write it for `max_uorfs=0`. `repro_no_cds.R`: default `simCountTablesRegions()` followed by `simNGScoverage()` errors `cannot open the connection`; warning names the missing RDS. Documented accounting for absent uORFs is never reached.
3. **Medium — non-CDS region sets fail.** `repo/R/Sim_reads_from_counts.R:338` uses `cds_ranges` even when only requested regions were loaded at lines 304–305. `repro_no_cds.R` selects only leader with probability 1 and gets `object 'cds_ranges' not found`.

## Checks

- The initial mandated test run stopped at load: missing `extraDistr`; installed CRAN extraDistr 1.10.0.5 into task-local `Rlib`, then reran under the requested 30 GB claim: **1,510 pass, 0 fail, 0 test warnings, 10 skip**; testthat duration 1,263.4 s, wall time including startup 1,283.8 s. No failure messages. Skip reasons are in `tests-local.log`; a separate dbDisconnect warning and intermittent database finalizer messages appeared outside the testthat warning count.
- `mapping.R` checks transcript positions 1, 10, 11, 20 against a two-exon plus/minus model and round-trips all positions; a six-base junction fragment produces `3M90N3M` at position 108 on both strands.
- `integration.R` and `verify_bam.py`: 429 and 516 BAM records; `samtools quickcheck` and `samtools view` pass. Positions, strands, CIGARs and sequences exactly match expanded ground truth, including 68 splice-junction and 857 minus-strand records. Seeded truth tables are identical for 1 and 4 requested cores (2 samples, hence 2 actual workers); all four OFST files at both benchmark sizes are byte-identical across 1 and 4 cores.
- `repro.R` checks stop counts across 10 seeds and 30 genes per seed with up to five overlapping uORFs; all CDSs and annotated uORFs had exactly one stop; the drop path removed 15 invalid uORFs across 8 of those runs. See `repro.log`. R printed a database finalizer error during repeated runs but continued; its origin is unverified. These checks do not exhaust overlap, dropped-uORF, boundary, or mapping cases.
- Results and profiles for the README-scale performance runs: `timing.tsv`, `bench/*/times.tsv`, `bench/*/profile_summary.txt`. Peak RSS is sampled process-tree RSS (shared fork pages counted repeatedly), not peak unique memory; R startup is included in wall time. BLAS/OpenMP threads limited to one per process.

## Performance

| Genes | Cores | Wall s | Genome s | Coverage s | Peak tree RSS GiB |
|---:|---:|---:|---:|---:|---:|
| 500 | 1 | 89.81 | 44.13 | 24.94 | 1.96 |
| 500 | 4 | 72.57 | 37.62 | 13.37 | 9.73 |
| 3,000 | 1 | 224.75 | 126.83 | 75.92 | 2.49 |
| 3,000 | 4 | 177.38 | 108.45 | 53.47 | 11.66 |

All four runs succeeded. Counts/regions steps and exact values are in `bench/*/times.tsv`; wall includes R startup and report/profile work. One-core runs had Rprof enabled; four-core runs did not, so profiling overhead and ordinary shared-host variability affect comparisons. Process-tree RSS is sampled once per second and double-counts shared pages. `timing.tsv` also gives the maximum single-process RSS. Four OFST files match byte-for-byte between core counts at each size.

The slowest measured step was 3,000-gene genome creation (126.83 s, 1 core). The sampled R profile attributes 96.96 s to `create_uORFs`, including 37.82 s under `[[<-` and 15.56 s under `create_uorf_ranges`; see `profile_uorfs.txt` and `bench/3000_1/profile_summary.txt`. These are inclusive sampled times, not direct wall-time measurements of isolated functions.
