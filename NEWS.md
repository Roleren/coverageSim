# coverageSim (development version)

* `simulation_preset()` / `run_simulation_preset()`: clean / realistic / noisy
  difficulty presets, plus `codon_bias_preset("clean")` and `"noisy"`; validation
  in docs/PRESETS.md. `rdg_parameters(scan_coupling = 1)`: 0 decouples uORF
  initiation from the CDS. Defaults unchanged.
* `learn_end_bias()`: `soft_clips = "drop"` (default) / `"trim"` and
  `added_base = NULL / "auto" / A / C / G / T` correct for the non-templated 5'
  base; `flank = TRUE` (default FALSE) also learns the nuclease cut site, the
  base just outside each read end, as `five_prime_flank_bias` /
  `three_prime_flank_bias`.
* `fragment_geometry`: the cut-site tables `five_prime_flank_bias` /
  `three_prime_flank_bias` (same names as the learner output, multiplied with the
  terminal `five_prime_bias` / `three_prime_bias` chemistry weights) and
  `added_base = list(base, rate)` (rate 0 = off). Both opt-in; defaults unchanged.
  `method_end_bias(method, nuclease, added_base)` gives default chemistry and
  cut-site tables per library method and nuclease. This is the one cut-site API: the experimental `five_prime_bias$flank` field
  (a terminal k-mer spanning the cut, never on master) was dropped and is
  rejected with a pointer to the flank tables.
* `technical_artifacts`: linear multimapper eligibility, skipped when unused;
  faster fragment allocation (outputs unchanged).
* New opt-in faster fragment draw: `options(coverageSim.fragment_draw = "two_stage")`
  (default `"joint"`, unchanged). Reads are drawn first over positions, then over
  fragment lengths and offsets. The distribution and the per-region read budgets
  are the same as with the default draw, but for a given seed the simulated reads
  are **not bit-identical** to the default draw. It applies to joint draws without
  end-bias or codon-bias tables; otherwise the default draw runs.
  When to use it: for deep or large simulations where speed matters (fragment
  step about 10x faster with `realistic_rfp_settings()`). Keep the default when
  outputs must match earlier runs exactly. With an end-bias or codon-bias table
  there is no speed-up (the default draw runs, and the end bias still works).
