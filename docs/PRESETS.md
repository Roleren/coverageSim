# Difficulty presets: clean, realistic, noisy

`simulation_preset("clean" | "realistic" | "noisy")` returns one documented parameter set over every simulator layer;
`run_simulation_preset()` runs genome -> gene counts -> region counts -> reads with it (options are set only for the call).
The codon layer is also available alone as `codon_bias_preset("clean")` / `codon_bias_preset("noisy")` next to
`codon_bias_preset("human_v2")`. Purpose (Carina, 2026-11-23): test ORF callers and metrics at known difficulty levels.

## Parameter sets

| Layer | clean (easier than real) | realistic (calibrated human) | noisy (harder than real) |
| --- | --- | --- | --- |
| Genome | human lengths, <= 1 uORF | human lengths, <= 1 uORF | human lengths, <= 3 uORFs |
| Gene dispersion `dispMeanRel` | 1/x + 0.001 | 4/x + 0.01 | 8/x + 0.4 |
| RFP regions leader / CDS / trailer / uORF | 0.006 / 0.98 / 0.006 / 0.008 | 0.025 / 0.938 / 0.017 / 0.02 (B1) | 0.09 / 0.75 / 0.13 / 0.03 |
| Region DMN concentration | 50 | 1 | 0.5 |
| Codon layer | `codon_bias_preset("clean")`: per-genome dwell, fixed sharpening, no pauses, no field | `"human_v2"` | `"noisy"`: human_v2 + pause_sd 1.0, zero mask 0.3, field variances x2 (50% shared, gene sd 0.8), wider smear, footprint GC bias |
| 3'UTR layer | off | `utr3_preset("human")` | human, plateau_prob 0.12 |
| DMN `dmn_alpha_scale` | 30 | 6 | 1.5 |
| RNase kernel (RFP) | 0.05, 0.4, 0.05, **10**, 0.4, 0.05, 0.05 | default 0.5, 2, 1, **10**, 2, 1, 0.5 | 1, 3, 2, **10**, 3, 2, 1 |
| Fragment geometry | learned (SRR12651187) collapsed to the best offset per length, no frame bias | `realistic_rfp_settings()` | learned offsets +-1 nt, off-peak offsets >= 0.6 x best, 12% 20-24 nt reads (offsets 8-12), no frame bias |
| Technical artefacts | none (ofst output) | 15% multimapped molecules (1 secondary) | 35% multimapped (2 secondaries), 20% PCR duplicates, 4 spike sites in 34% of libraries (300-3000 copies) |

## Validation against real libraries

One metric code for real and simulated data (`open_workspace/claude-megafst/coverageSim_presets/scripts/metrics_lib.R`):
unit = P-sites on one canonical transcript per gene (ORFik `filterTranscripts(30, 150, 30, longestPerGene = TRUE)`;
the simulated genome has one transcript per gene), P-site offset per read length = the offset in 9-17 nt that maximises
the CDS frame-0 share (the same rule for real and simulated reads), depth = every library thinned to 5e5 mRNA P-sites.
Real: one replicate pair from each of 40 human studies (80 libraries, >= 2e6 unique P-sites, not flagged by the re-trim
audit) and 25 yeast libraries, unique mappers (all mappers for the multimapper share), plus the ultraperiodic merged
tracks. Simulated: 1,000 genes, 2 replicates per preset.

See `presets_and_limitations.md` in the same folder for the table and figures. In short: "realistic" lies inside the
human real range for frame 0, region shares, read lengths, multimapping and replicate dispersion; "clean" sits at the
ultraperiodic reference for frame 0 and above every real library for CDS share; "noisy" sits beyond the real range
for frame 0, UTR share, short reads, multimapping and dispersion. Known gap for all presets: codon-level
heterogeneity (excess CV^2 of codon counts, max/mean codon count) is below every real library, including the
ultraperiodic one (the human_v2 calibration targets windows, not single codons).

## Yeast

The shipped codon, RNase and end-bias tables are human (`*_human.csv`, `realistic_rfp_geometry_SRR12651187.rds`), and
`human_tx_lengths()` / B1 shares are human. A yeast preset needs: yeast transcript lengths (Ensembl yeast has almost
no UTR annotation; use TIF-seq / CAGE UTRs), yeast region shares (real median here: CDS 0.96, 5'UTR 0.019, 3'UTR 0.011),
a yeast codon dwell table (the codon dwell atlas has species fits in `codon_dwell/atlas_<species>/`), a learned yeast
geometry (`learn_fragment_geometry()` on one yeast library, e.g. the ultraperiodic set, modal length 28) and a
multimapper rate (median 0.49 of mRNA reads here, much higher than human because of duplicated ribosomal-protein genes).
The real yeast distribution is measured with the same code (figure), so a yeast preset can be validated the same way.
