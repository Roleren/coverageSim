# coverageSim (development version)

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
