# Prime sieve measurement suite

This developer suite checks the experimental inspector and compares segment
budgets, preparation costs, and effective worker counts. It does not change the
production sizing policy. Start with [the findings report](REPORT.md) and
[the inspector guide](../cache-inspector.md).

First create an isolated build with the
[developer toolkit](../prime-sieve-inspector/README.md). Run the commands below
from the repository root with that build's `library/` and `package/` paths. R, testthat, Python 3, and a C++
compiler are required; the suite adds no dependencies to the package.

```sh
Rscript .github/development/prime-sieve-benchmarks/run.R /path/to/build/library /path/to/build/package /path/to/new-results 7
Rscript .github/development/prime-sieve-benchmarks/controls.R /path/to/build/library /path/to/new-results
python3 .github/development/prime-sieve-benchmarks/confirm-native.py /path/to/new-results --source-root /path/to/build/package --max-cores 10 --max-threads 10
Rscript .github/development/prime-sieve-benchmarks/analyze.R /path/to/new-results
```

Run these sequentially on an otherwise quiet machine. Use real host values for
the native script's core/thread arguments. In a restricted execution sandbox,
R's physical-core detection can fall back to one and change algorithm selection;
check the saved metadata before interpreting results. Native confirmation uses
fixed workloads near `4e15` and `8e15` and rejects a core policy under which the
first is not naturally Big. The R driver adapts those bases where possible, but
the fixed native confirmation then needs an explicitly revised workload grid.

- `validate.R` runs the existing sieve/inspector test files and 30 independent
  primality-reference checks, including prime squares, segment boundaries, and
  the supported numeric limit. Reference work is outside measurement.
- `run.R` screens seven serial budgets, nominates two using equal-weight
  geometric means, tests one/two/four requested workers, and confirms the
  finalists with fresh shifted intervals. Every timed output is checked exactly.
  Seven screening repeats and fifteen confirmation repeats produced 879 measured
  calls on the recorded ten-core host. The main metric sums coarse C++ elapsed
  phases, including R-result conversion but excluding diagnostic-report assembly.
- `controls.R` compares full ordinary/inspected calls across Small, Medium, and
  Big, using fifteen randomized paired rounds of ten calls. Its R elapsed timer
  has roughly millisecond resolution; interpret very short cases cautiously.
- `confirm-native.py` creates a temporary source copy with precisely one
  asserted dispatch substitution: select uninstrumented workers even when
  supplying the override. It retains coarse setup/timing, compiles C++11 `-O3`,
  checks complete outputs outside timing, and saves compiler/source metadata.
  Native elapsed time excludes result concatenation and R conversion. Production
  source files are never edited. The temporary build is retained for inspection.
- `analyze.R` writes medians, quartiles, ratios, worker/phase tables, and plots.
  Raw observations remain alongside summaries. Results are descriptive; this
  suite does not infer statistical significance or hardware cache behavior.

`results-2026-11-01/` contains the recorded run. It has no prime vectors or
compiled binaries. Source hashes identify the uncommitted inspector build in
addition to its base Git revision. Build/install and sanitizer validation of
that implementation are documented in the inspector guide; sanitizer timings
must not be mixed with release benchmarks.

The original measurement scripts are preserved under
`results-2026-11-01/scripts/`. Their original source paths assumed the inspector
was in the checkout. Use the current scripts above with the isolated build for
new runs; the historical source and artifact hashes remain unchanged.
