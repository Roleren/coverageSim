# coverageSim speed core pass 3

- Branch `speed/core3`, implementation commit `7a345a7`: batch uniform region multinomial sampling in C++, using R's RNG. Arbitrary profiles and deferred probabilities retain the R path; `options(coverageSim.engine = "R")` selects it throughout.
- B3k, one thread, same driver: baseline genome/counts/regions/coverage 36.414/0.273/0.903/249.406 s (total 286.996 s); optimized 34.960/0.278/0.884/235.304 s (total 271.426 s), 5.4% faster. Single runs; timings are sensitive to system load. Logs: `../baseline_one_thread.log`, `../optimized.log`.
- Focused test: 20 passes, zero failures (`../focused.log`); full suite: no reported failures, 14 warnings and 10 skips, with three R finalizer error messages outside the test summary (`../suite.log`).
- Five paired B3k seeds, 12 libraries each: zero region-budget or truth-accounting mismatches and zero truth-file differences. Position and fragment histograms were identical in these pairs (`../five_seed_checks.log`). Seeded equivalence is evidence for tested cases, not a proof for all inputs.
- Final B20k, one thread, 32 GB reservation: running; results pending (`../final20.log`). The first launch failed before starting simulation because the dependency library was absent from the environment; the final run uses the existing read-only library.
- Independent Codex review was not run: delegation is disabled for this task. B20k speed target below 300 s remains unverified.
