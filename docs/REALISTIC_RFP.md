# Realistic RFP simulation settings (opt-in)

```r
options(coverageSim.codon_bias = codon_bias_preset("realistic"))   # t57 base, pause_sd 0.8, pause_prob 1.0
extra <- realistic_rfp_settings()   # dmn_alpha_scale = 6, learned fragment_geometry (A-site, frame_bias strength 0.25)
exp <- do.call(simNGScoverage, c(list(simGenome = g, count_table = counts, ...), extra))
```

`codon_bias_preset("t57")` and all defaults are unchanged. Source of the learned geometry: SRR12651187, 6 M read
subsample, 2264 canonical transcripts, `learn_fragment_geometry()` (lengths 25-36, offsets with >= 5% of the length's
best in-frame fraction) and `learn_end_bias(by_length = TRUE)$frame_bias`.

## Narrower window band than real (conservative false-flag test)
With these settings the simulated replicate window noise is below real at E >= 50 (rel sd .09/.08/.07 vs real .16/.14/.18)
because real data has extra library-specific density noise that the knobs cannot add without also moving the frame.
False density flags against the true simulated alphas fall to 0.2-0.3% (nominal 5%): a false-flag test on simulated
null genes is conservative, i.e. it understates the real false-positive rate. Per-nt clumping is 3.4x Poisson (real 2.7).
Gene-level frame variation (real sd 0.06-0.15) is not simulated (sim 0.005).

## t69 before / after (3 seeds; high / mid depth; t57 = before)
| metric | before (t57) | after (realistic) | real / target |
|---|---|---|---|
| false density flags vs truth E, real band | 4.0 / 2.6% | 0.3 / 0.2% | <=5% |
| same vs A-site table E (t57 def.) | 20 / 15% | 5 / 4% | |
| codon r 600-1200 | 0.52 / 0.53 | 0.52 / 0.51 | 0.54 / 0.45 |
| codon r 1200-2500 | 0.56 / 0.53 | 0.63 / 0.62 | 0.56 / 0.57 |
| rel sd window totals, E 15/28/56/112/212/445/806 | .30/.24/.20/.19/.16/.17/.19 | .25/.19/.15/.12/.09/.08/.07 | .31/.25/.18/.17/.16/.14/.18 |
| shared misfit sd vs truth E | 0.007 / 0.045 | 0.015 / 0.000 | 0.18-0.25 |
| replicate window frame sd (psi0 ratio) | 0.164 (3.1) | 0.064 (1.3) | 0.052 |
| false frame flags fail / not understood, human band | 4.7 / 2.9% | 0 / 0% | 3.4 / 3.2 (thinned real) |
| per-nt count clumping (x Poisson) | 6.2 | 3.4 | 2.7 |
| pooled frame 0/1/2 after real-style single offset per length | .65/.12/.24 | .62/.15/.23 | .62/.17/.20 |
| gene-level shared frame sd | 0.016 | 0.005 | 0.06-0.15 |

Caveats: n = 3 seeds; geometry and frame strength chosen on one seed (101, mid depth). Frame 1 vs 2 asymmetry remains.
