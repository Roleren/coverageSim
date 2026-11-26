# coverageSim (development version)

* `rdg_parameters(reinit_d0 = )` (opt-in, default NULL = the hard step, results
  unchanged): graded ternary-complex re-acquisition after a uORF. Reinitiating
  ribosomes load gradually while they scan, so the share reaching a start d nt past
  the stop is `1 - exp(-(d - reinit_distance) / reinit_d0)`; unloaded ones scan past
  starts without initiating. Motivation: in the genome-wide human translon panel
  (691 single-uORF genes) CDS density relative to uORF reads rises +0.33 log2 per
  doubling of the uORF-to-CDS distance (se 0.12), a graded rise the step cannot
  give. Caveat: confounded with leader length (not separable from cap distance +
  uORF length + distance), weaker with the classifier P (+0.18), correlational; it
  supports the shape, not a value, so `reinit_d0` is not fitted. A condition that
  changes `reinit_d0` counts as ISR-type for `translon_exempt = "auto"`.

* `rdg_parameters(translon_exempt = "auto")`: genes that keep the start-codon table for
  their leader translons although `p_init_translon` is set. "auto" (default): when a
  condition changes `reinit_distance` (ISR-type scenario), genes with the
  delayed-reinitiation layout (a translated uORF followed by a uoORF over the CDS start:
  ATF4 / ATF5 / PPP1R15A-like) are exempt, so the ATF4-like CDS gain under ISR stays
  3.7x (1.1x without the exemption); genes with only a uORF keep the fitted 0.3.
  Without such a scenario no gene is exempt and results are unchanged. Also NULL, transcript ids or a function; `simRDGflow()` records
  the exempt genes in `metadata(result)$rdg_translon_exempt`.

* `rdg_parameters()`: new defaults for the uORF -> CDS coupling (a deliberate realism
  change; every `simRDGflow()` result changes). Fitted on 62 human marker genes with
  Ribo-seq-detected uORFs (clean-CDS usage, 5-fold cross-validation by gene, held-out
  RMSE on the logit scale 1.67 with the old defaults):
  - `p_init_translon = 0.3`: one initiation probability for the leader translons
    instead of the codon table (RMSE 1.00); `translon_classes` (default uORF,
    uoORF, NTE: the leader classes of the fit) sets which classes; dORF, doORF and
    intORF keep the codon table.
  - `uoorf_cds_factor = 0.8`: extra CDS repression per overlapping uoORF (RMSE 0.93
    together with `p_init_translon = 0.3`; 2000 gene bootstraps: gain positive in
    100%, 0.3 chosen in 61%, 0.6-0.9 for the factor in 91%).
  - `translon_gene_sd` (opt-in, default 0): per-gene log-normal factor on the
    translons' initiation, drawn once per `simRDGflow()` call; 0.5 matches the real
    spread of clean-CDS usage. It changes the random stream, so it stays off.
  Impact (500 genes, `architecture = "all"`, same seeds): CDS share of translated
  reads 0.44 -> 0.73 (per-gene median 0.30 -> 0.73; 349 of 500 genes change by
  > 0.1); uORF read share 0.167 -> 0.045, uoORF 0.251 -> 0.134. With leader
  translons at 0.3 the ATF4-like response to a longer `reinit_distance` is much
  weaker (CDS gain under ISR 3.7x -> 1.1x in the test gene). The old behaviour:
  `rdg_parameters(p_init_translon = NULL, uoorf_cds_factor = 1)`.
  `simulation_preset()` does not call `simRDGflow()` and is unchanged.

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
