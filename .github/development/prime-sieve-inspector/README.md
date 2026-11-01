# Developer-only prime sieve inspector

This toolkit preserves the inspector without adding code to the production
package. Everything here is excluded from package distributions by the existing
`.Rbuildignore` rule for `.github/`. The production sieve remains the sole
maintained algorithm implementation.

`instrumentation.patch` contains the changes to the sieve, its interface/header,
and the generated cpp11 registration. The additional C++, header, R wrapper, and
test files live here in their package-relative paths. `build.R` assembles them
only in an external build directory, installs into its own R library, and runs
the two focused test files in a fresh R process. It never edits, stages, or
commits anything in the checkout, or installs into your usual R library.

## Build from current sources

From the repository root:

```sh
MAKEFLAGS=-j4 Rscript .github/development/prime-sieve-inspector/build.R /tmp/rcppalgos-inspector
```

Choose a new output directory with an existing parent. It must be outside the
repository; an existing directory is never overwritten. `MAKEFLAGS=-j4` is
optional and controls compilation parallelism, not sieve worker counts. On
Windows, set that environment variable separately or omit it.

The current-source mode includes tracked working-tree files (including local
edits), plus non-ignored untracked files under package directories such as
`src/`, `R/`, `inst/`, and `tests/`. It excludes `.Rbuildignore` matches, compiled
objects/libraries, and unrelated untracked files at the repository root. Keep
any new required root-level package files tracked. Source paths containing
newlines are not supported by this developer script.

The output contains:

- `package/`: the instrumented source copy, available for inspection.
- `library/`: the separately installed experimental RcppAlgos package.
- `install.log` and `tests.log`: build and focused-test results.
- `build.rds`, `base-files.csv`, and `instrumented-files.csv`: provenance,
  validation status, and source hashes.

Git, R, the package's build dependencies, and testthat must already be available.
No dependencies are downloaded by this script. A patch conflict stops the build
before compilation; resolve it deliberately instead of forcing an application.
Successful patch application alone is not evidence that new algorithm changes
were instrumented appropriately—review phase boundaries when the sieve changes.

## Reconstruct the recorded study

`base-revision.txt` records `0364f49bc64476829c6feae63ec55148521e5c79`.
To reconstruct its exact instrumented sources, without switching branches:

```sh
MAKEFLAGS=-j4 Rscript .github/development/prime-sieve-inspector/build.R /tmp/rcppalgos-inspector-study 0364f49bc64476829c6feae63ec55148521e5c79
```

This mode archives that local Git revision, applies the preserved patch and
helpers, and verifies all eight inspector-related source files against the
recorded hashes before installation. The archived `base/` and `base.tar` also
remain in the output directory. The revision must exist locally; no fetch or
checkout occurs. Exact source reconstruction does not guarantee identical
binaries or timings across toolchains or hosts.

## Use the isolated build

Start a fresh R process so an already-loaded production DLL cannot be reused:

```r
library(RcppAlgos, lib.loc = "/tmp/rcppalgos-inspector/library")
result <- RcppAlgos:::primeSieveInspect(4e15, 4e15 + 1e6,
                                      nThreads = 1, segmentBytes = 1024^2)
result$summary
result$phases
result$kernels
```

The override requires an interval naturally assigned to Big; that cutoff depends
on detected cores and requested threads. See the
[inspector guide](../cache-inspector.md) for interpretation and limits, and the
[benchmark README](../prime-sieve-benchmarks/README.md) for the complete suite.
For example, from the repository root:

```sh
Rscript .github/development/prime-sieve-benchmarks/run.R /tmp/rcppalgos-inspector/library /tmp/rcppalgos-inspector/package /tmp/rcppalgos-inspector-results 7
python3 .github/development/prime-sieve-benchmarks/confirm-native.py /tmp/rcppalgos-inspector-results --source-root /tmp/rcppalgos-inspector/package --max-cores 10 --max-threads 10
```

Use your actual policy counts for the native command and run timing experiments
sequentially. Removing an external build directory later removes that build;
the toolkit, source patch, and recorded results remain in the repository.

## Maintaining the toolkit

When production changes conflict, revise the instrumentation in a disposable
copy and update the patch and helper files here. Avoid maintaining a second full
`Eratosthenes.cpp`. Regenerate native registration with `cpp11::cpp_register()`
in that disposable copy when interfaces change; the patch currently preserves
the previously generated registration, not hand-written registration code.
Review the generated R wrapper against the project's hand-written wrappers.

The base revision and expected hashes describe a specific reconstruction. Update
them together only when deliberately replacing that baseline; preserve a dated
study's original reconstruction materials if doing so. Historical measurements
and their hashes should not be rewritten to imply they came from a newer build.

The November 1 report and results remain in
[`../prime-sieve-benchmarks/`](../prime-sieve-benchmarks/REPORT.md).

## Separation validation

On November 1, 2026, both the pinned-revision and current-source workflows
compiled and passed `testCacheInspector.R` and `testPrimeSieve.R`. The pinned
copy matched all eight recorded source hashes. The relocated benchmark driver
completed 531 measured output comparisons and 30 independent boundary checks;
its native confirmation completed 36 further measured comparisons. Those were
integration checks in temporary directories, not replacements for the recorded
study. The analyzer also ran successfully on the temporary results.

A normal `R CMD build --no-build-vignettes --no-manual --no-resave-data` archive
was inspected: it contained neither this toolkit nor inspector helpers/hooks.
That archive installed successfully, exposed no inspector R function or native
registration, and passed the prime-sieve tests. Existing-build protection and
the benchmark driver's library/source pairing checks were also exercised.
Validation used the local macOS/R toolchain; a full `R CMD check` and other
platforms were not run for this relocation.
