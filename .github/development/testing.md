# Testing workflows

Run commands from the repository root. Use an isolated source copy and library
so validation neither replaces the installed package nor reuses stale objects.
This guide describes current tools. Architecture, mathematical invariants,
memory estimates, and retained performance decisions are in
[number-theory.md](number-theory.md); study entry points are in
[README.md](README.md).

## Package installation and focused tests

Package tests must use at most two threads, including automatic parallel paths.
Wider concurrency checks belong in the CRAN-excluded development fixtures.

The following exports the **Git index**, including staged changes and excluding
unstaged edits. To test other contents, prepare that source snapshot explicitly
and use its path as `package_source`. Neither building the working tree nor
loading an already installed package verifies the index.

```sh
validation_root=$(mktemp -d)
package_source="$validation_root/source"
package_library="$validation_root/library"
mkdir "$package_source" "$package_library"
git checkout-index --all --prefix="$package_source/"
MAKEFLAGS=-j4 R CMD INSTALL --preclean --no-multiarch \
    --library="$package_library" "$package_source"
Rscript --vanilla - "$package_library" "$package_source" <<'RS'
args <- commandArgs(TRUE)
package_library <- normalizePath(args[1])
package_source <- normalizePath(args[2])
library(RcppAlgos, lib.loc = package_library)
dll <- normalizePath(getLoadedDLLs()[["RcppAlgos"]][["path"]])
stopifnot(startsWith(dll, paste0(package_library, "/")))
files <- c("testPrimeCount.R", "testNthPrime.R", "testPartialSieve.R",
           "testPrimeSieve.R", "testPrimeFactorize.R", "testPrimeFactorizeSieve.R",
           "testEulerPhiSieve.R", "testNumDivisorSieve.R", "testIsPrimeRcpp.R",
           "testDivisorsRcpp.R", "testDivisorsSieve.R", "testErrors.R",
           "testCombPermSample.R", "testPartitionsSampleRank.R",
           "testPermuteGeneral.R", "testExpandGrid.R", "testComboGroups.R")
for (file in files) {
    testthat::test_file(file.path(package_source, "tests/testthat", file),
                       reporter = "summary", stop_on_failure = TRUE,
                       package = "RcppAlgos", load_package = "installed")
}
registered <- getDLLRegisteredRoutines("RcppAlgos")$.Call
stopifnot(registered$`_RcppAlgos_PrimeCountCpp`$numParameters == 4L,
          registered$`_RcppAlgos_PrimeSieveCpp`$numParameters == 6L,
          registered$`_RcppAlgos_NthPrimeCpp`$numParameters == 4L,
          !"run_testthat_tests" %in% names(registered))
for (name in c("primeCount", "primePhi", "semiPrimeCount", "triPrimeCount",
               "primeSieve", "nthPrime")) {
    stopifnot(length(tools::checkRd(file.path(package_source, "man",
                                             paste0(name, ".Rd")))) == 0L)
    example(name, package = "RcppAlgos", lib.loc = package_library,
            ask = FALSE, echo = FALSE)
}
stopifnot(length(tools::checkRd(file.path(package_source, "inst/NEWS.Rd"))) == 0L,
          length(tools::codoc(package = "RcppAlgos",
                             lib.loc = package_library)) == 0L)
RS
```

`MAKEFLAGS=-j4` controls compilation, not algorithm workers. Dependencies must
already be installed. Run tests from the chosen snapshot, and retain its revision
or source hashes with the logs. Prime-count tests run first in a fresh R process
because earlier calls can populate the shared Phi cache. Use separate processes
when a particular test requires an initially empty cache.

These focused tests do not constitute a full `R CMD check` or cross-platform
validation. Run broader package checks when the change's scope warrants them;
record their actual results separately. Examples and `codoc` check interface and
documentation consistency, not performance or race freedom.

## Full-suite timing and CRAN test workloads

Use the isolated library above to time the complete test suite. For test-only
changes, the same compiled package can serve both runs; pass the original and
modified source directories separately. Run each in a fresh R process, with no
other builds or benchmarks running, and keep the seed and reporter consistent.

```sh
Rscript --vanilla - "$package_library" "$package_source" <<'RS'
args <- commandArgs(TRUE)
package_library <- normalizePath(args[1])
.libPaths(c(package_library, .libPaths()))
library(RcppAlgos)
stopifnot(startsWith(normalizePath(getLoadedDLLs()[["RcppAlgos"]][["path"]]),
                    paste0(package_library, "/")))
set.seed(104729)
elapsed <- system.time(results <- testthat::test_dir(
    file.path(args[2], "tests/testthat"), reporter = "silent",
    package = "RcppAlgos", load_package = "installed", stop_on_failure = TRUE
))
results <- as.data.frame(results)
timings <- aggregate(real ~ file, results, sum)
print(timings[order(-timings$real), ], row.names = FALSE)
print(elapsed)
print(colSums(results[c("failed", "error", "warning", "skipped")]))
RS
```

Keep expensive checks tied to the behavior they exercise. The exhaustive
multiset oracle retains every frequency, width, target, and rank, but compares
named batches to reduce testthat overhead. Partition iterator fixtures retain
full native traversal and state checks, with bounded endpoint/interior random
access and short windows at both ends of GMP-sized result spaces. Parallel
iterator batches must still reach the production worker thresholds: 20,000
rows generally and 40,000 for the uncapped partition/composition fixtures in
`testParallelIterators.R`. Capped partition paths elsewhere can require 150,000;
verify the dispatch before reducing their sizes. Small data alone does not
establish parallel coverage.

On November 19, 2026, the test-only reductions were measured against
`f2359a5d`, using the same isolated optimized build on macOS arm64 / R 4.6.1.
Fresh processes used seed 104729 and the same ListReporter timing harness:

| Workload | Original | Reduced |
| --- | ---: | ---: |
| Complete suite, wall time | 115.39 s | 59.70 s |
| Partition iterator tests | 45.60 s | 12.51 s |
| Exhaustive multiset oracle | 19.91 s | 0.97 s |
| Parallel iterator tests | 3.23 s | 0.16 s |

Both full runs exercised the same 247 test blocks in 42 files, with no failures
or warnings and one existing multiset-composition specification skip.
These are local test timings, not Windows measurements or an overall
`R CMD check` time.
The partition iterator timing combines reduced random-access sampling, shorter
GMP traversal windows, and removal of forced garbage collections; it does not
isolate the benefit of shortening sequential traversal.

A separate covr/gcov comparison ran both complete suites in fresh processes
against one instrumented build, resetting native counters between runs. Both
exercised exactly the same 18,995 of 19,847 instrumented production lines:
17,466 C++ lines, 900 header lines, and 629 R lines. Generated registration and
native test sources were excluded from those totals. No previously covered
lines were lost, and both instrumented runs passed. Line coverage complements
the fixture/dispatch review; it does not establish identical branch coverage
or constitute a correctness proof. Instrumented runtimes are not used in the
timing comparison above.

That coverage comparison predates moving the direct native regressions out of
the package on November 20. Those internal checks now run separately using the
optional development workflow below; an ordinary package check does not run
them.

## Optional native regressions outside the package build

The direct C++ regressions live in
[native-regressions/](native-regressions/README.md), excluded from package
archives by `.Rbuildignore`. They preserve the raw GMP conversion, poisoned
buffer/storage reuse, forced-GMP counting/ranking, and native primality guard
checks. Run them explicitly when changing these internals and before release:

```sh
Rscript --vanilla .github/development/native-regressions/run.R
```

The runner compiles a fresh temporary bridge from the current production
sources using cpp11 and GMP, and retains its source hashes and results. It uses
ordinary exceptions for assertions, with no Catch/testthat C++ dependency.
Normal package installation no longer compiles these direct native regressions
or Catch, or exports `run_testthat_tests`; `testthat` remains in `Suggests` for
the R tests. `xml2` is no longer needed for native test report parsing or
coverage CI. It is declared under `Config/Needs/website` for `dev/build_site.R`.

On November 20, a clean `R CMD INSTALL --preclean --no-multiarch --install-tests`
comparison on the same macOS/R setup with `MAKEFLAGS=-j4` took 23.43 s before
this removal and 19.04 s after. The installed shared library decreased from
2,788,736 to 2,143,184 bytes. These single-run measurements include package
compilation, linking, and installation checks, with dependencies already
installed; they do not measure dependency provisioning or CRAN's total check
time. The development bridge separately passed all seven retained checks.

## Package archive and broader check

From the temporary validation directory, build the selected source snapshot
with rendered vignettes, then check the resulting archive:

```sh
MAKEFLAGS=-j4 R CMD build --no-manual "$package_source"
MAKEFLAGS=-j4 R CMD check --no-manual RcppAlgos_2.10.3.tar.gz
```

Use the filename printed by `R CMD build` when the version changes. These
commands require Pandoc for the R Markdown vignettes. On this macOS host it is
bundled with RStudio; command-line R can find it with the temporary setting
`RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64`.
Skipping vignette creation produces missing-`inst/doc` warnings, so that
shortcut cannot establish a clean archive check. `--no-manual` omits the PDF
reference manual; it does not omit the R tests, examples, or vignette rebuilds.

On November 20, the final archive passed `R CMD check --no-manual` with
`Status: OK`, including examples and vignette rebuilding. The R suite reported
5,669 passing expectations, no failures or warnings, and its one existing skip.
Remote CRAN/Bioconductor index lookups were unavailable in the sandbox;
dependency checks against the installed libraries passed. The PDF manual and
Windows checks were not run.

## Native correctness and sanitizer checks

### Package floating-conversion and mapping checks

The numeric regressions in `tests/testthat/testNumericConstraints.R` should also
run against a package built with UBSan, including `float-cast-overflow`. Ordinary
R output does not expose an unsafe cast or signed overflow reliably. Use the
isolated source snapshot and library described above; preserve any existing user
Makevars settings in the temporary file before appending these flags:

```make
PKG_CXXFLAGS += -fsanitize=undefined,float-cast-overflow -fno-sanitize-recover=all
PKG_LIBS += -fsanitize=undefined,float-cast-overflow -fno-sanitize-recover=all
```

```sh
R_MAKEVARS_USER="$validation_root/Makevars.ubsan" MAKEFLAGS=-j4 \
    R CMD INSTALL --preclean --no-multiarch \
    --library="$package_library" "$package_source"
UBSAN_OPTIONS=halt_on_error=1 Rscript --vanilla - \
    "$package_library" "$package_source" <<'RS'
args <- commandArgs(TRUE)
library(RcppAlgos, lib.loc = args[1])
stopifnot(startsWith(normalizePath(getLoadedDLLs()[["RcppAlgos"]][["path"]]),
                    paste0(normalizePath(args[1]), "/")))
files <- c("testGmpConvert.R", "testNumericConstraints.R", "testZeroPartitions.R",
           "testIsPrimeRcpp.R", "testPrimeFactorize.R", "testDivisorsRcpp.R",
           "testPartitionsDesign.R", "testPartitionTypeRefactor.R")
for (file in files) {
    testthat::test_file(file.path(args[2], "tests/testthat", file),
                       reporter = "summary", stop_on_failure = TRUE,
                       package = "RcppAlgos", load_package = "installed")
}
RS
```

This configuration was verified with Apple Clang 17 and R 4.6.1 on macOS arm64.
The sanitizer runtime linked and loaded without preloading. Do not substitute
an ordinary installed package for the instrumented library. Changes to shared
input preparation also warrant the full package test suite.

The numeric regressions cover native mapping and its range errors, including
widths 1024/1025 with source values below `2^53`, `INT64_MIN` source values,
and positive/negative mapping-product overflow. The original native search
runs only when a conservative bound guarantees safe arithmetic; it has no
GMP fallback. Existing GMP counting, ranking, and sampling remain separate.
When changing the search arithmetic, also run `testPartitions*.R`,
`testConstraintsClass.R`, `testCombPermResults.R`, and `testCombPermRank.R`
against the same instrumented library to cover counting, generation,
sampling/ranking, iterators, and parallel dispatch.

The zero-only regressions cover repeated partitions and compositions (including
weak compositions), retaining the requested width and integer/numeric result
type across generation, counting, sampling, ranking, and iterator exhaustion
and random access. They also check the constrained combination/permutation
routes and nearby frequency-based and positive-target cases.

The primality and factorization regressions cover the native/GMP dispatch
boundaries in serial and with two requested workers. The optional
`native-regressions/primality.cpp` checks also exercise the `IsPrime()`
upper-bound guard directly; R-facing calls route large inputs to GMP before
reaching that guard. Native
modular products rely on nonnegative residues below a modulus at most `INT_MAX`,
so their product is below `2^62` and needs no floating-point chunking fallback.

### Standalone native fixtures

The standalone fixtures need a C++11 compiler and thread support, but no R
runtime. Compile in a temporary directory. ASan/UBSan and TSan require separate
binaries; do not combine their flags. The following shell helper keeps flags
consistent without changing any production sources:

```sh
native_checks=$(mktemp -d)
check_sanitizers=address,undefined
build_check() {
    check_name="$1"
    shift
    clang++ -std=c++11 -O1 -g -pthread -Wall -Wextra -Wpedantic \
        -fsanitize="$check_sanitizers" -fno-omit-frame-pointer \
        -Iinst/include "$@" -o "$native_checks/$check_name"
}
```

Choose the relevant builds below. All paths refer to the source tree in which
the commands run; for an exact snapshot, first change to that snapshot's root.

```sh
build_check phi .github/development/phi-check.cpp \
    src/Phi.cpp src/Eratosthenes.cpp
build_check phi-parallel .github/development/phi-parallel-check.cpp \
    src/PhiParallel.cpp src/Phi.cpp src/Eratosthenes.cpp
build_check cache .github/development/prime-count-cache-check.cpp \
    src/PiLegendre.cpp src/PhiParallel.cpp src/PiPrime.cpp src/Eratosthenes.cpp
build_check partial-sieve .github/development/partial-sieve-check.cpp \
    src/P2.cpp src/P3.cpp src/Phi.cpp src/Eratosthenes.cpp
build_check methods .github/development/prime-count-methods-check.cpp \
    src/PiMeissel.cpp src/PiLehmer.cpp src/PiLegendre.cpp src/PiAuto.cpp \
    src/PhiParallel.cpp src/P2.cpp src/P3.cpp src/Phi.cpp src/Eratosthenes.cpp
build_check pi-prime .github/development/pi-prime-check.cpp src/PiPrime.cpp
build_check nth-prime .github/development/nth-prime-check.cpp \
    src/PiAuto.cpp src/PiLegendre.cpp src/PiMeissel.cpp src/PiLehmer.cpp \
    src/PiPrime.cpp src/Phi.cpp src/PhiParallel.cpp src/P2.cpp src/P3.cpp \
    src/Eratosthenes.cpp
```

The cache fixture includes `Phi.cpp` directly, the PiPrime fixture includes
`Eratosthenes.cpp`, and the nthPrime fixture includes `NthPrime.cpp`. **Do not
link those included sources a second time.** This exposes private boundaries
without adding test hooks to the package.

Run each selected binary with halt-on-error enabled, including the alternate
Phi/cache modes when relevant:

```sh
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$native_checks/phi"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$native_checks/phi" --no-cache
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$native_checks/cache" 4
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$native_checks/cache" 4 --allocation-failure
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$native_checks/nth-prime"
```

Other binaries take no required arguments. For race checks, set
`native_checks` to a **new** temporary directory, set `check_sanitizers=thread`,
repeat the selected builds, and run with `TSAN_OPTIONS=halt_on_error=1`.
On the recorded macOS toolchain, `:symbolize=0` suppresses external-symbolizer
warnings without disabling detection. Cache checks accept a worker count;
exercise more than one count when changing concurrency.

| Fixture | Independent evidence and special coverage |
| --- | --- |
| `phi` | Exhaustive partial sieves, inclusion-exclusion, tiny-table boundaries, invalid inputs, cache-disabled behavior, concurrent evaluations |
| `phi-parallel` | General cutoffs versus independent references and serial Phi; shared-table callers; basic `ThreadBatch` lifetimes |
| `cache` | Every populated cell versus independent partial sieves; page boundaries, allocation failure/recovery, concurrent Legendre calls |
| `partial-sieve` | P2/P3 enumeration including repeated factors, roots, table coverage, invalid budgets, concurrent callers and counting identities |
| `methods` | Independent prime sieve, known counts, square/cube/fourth-power neighbors, prepared-table reuse and insufficient coverage |
| `pi-prime` | Inclusive intervals, all sieve kernels, output types/order, synthetic hardware budgets, supported-domain boundaries |
| `nth-prime` | Exact interval updates, forward/backward selection, composite plateaus, multiple bounded windows, fixed budgets and table reuse |

`pi-prime --wide` additionally crosses the billion-position count-window limit.
Use a separate optimized build for that expensive integration check; the normal
fixture deliberately avoids it. Passing known-value tests does not detect data
races; use TSan for that purpose.

`phi-parallel --batch-only` checks callable transfers and exception-safe batch
lifetimes with at most two workers, without running the broader Phi cases.

Validation on November 18, 2026 (macOS arm64, `cd688747`): all 34 native runs
passed under ASan/UBSan and TSan without diagnostics. Temporary fixture-only
budget changes limited execution to two workers. Wider concurrency,
six-worker automatic Lehmer execution, and `pi-prime --wide` were not run;
direct Lehmer and policy-selection checks remained covered.

## Scheduler bounds and worker failures

The preparers verify expected production structure, then create instrumented
copies under the supplied temporary directory. They never modify production
files. A rejected pattern means the fixture must be reviewed against the new
source, not that an unchecked replacement should be forced.

### Phi scheduler and worker exception handling

```sh
phi_checks=$(mktemp -d)
python3 .github/development/prepare-prime-count-schedule-check.py "$phi_checks"
python3 .github/development/prepare-prime-count-thread-check.py "$phi_checks"
build_check schedule -I"$phi_checks" \
    .github/development/prime-count-schedule-check.cpp src/Phi.cpp src/Eratosthenes.cpp
build_check phi-failures -I"$phi_checks" \
    .github/development/prime-count-thread-check.cpp src/Phi.cpp src/Eratosthenes.cpp
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$native_checks/schedule"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$native_checks/phi-failures" --threads 2
Rscript --vanilla .github/development/prime-count-thread-check.R "$phi_checks/r-bridge"
```

These commands assume the ASan/UBSan helper configuration above. Build a
separate TSan worker-failure binary too. The scheduler substitutes synchronous
workers to verify once-only interval coverage and extreme budgets without
launching hundreds of OS threads; it is not a concurrency test. The failure
fixture uses real workers and gated joins, testing worker exceptions, later
launch failures, exception precedence, scratch lifetimes, and recovery. Both
preparers follow current `PhiParallel.cpp`, `PiLegendre.cpp`, and `ThreadBatch.h`.

### Prime sieve worker exception handling

```sh
sieve_checks=$(mktemp -d)
python3 .github/development/prepare-prime-sieve-thread-check.py "$sieve_checks"
build_check sieve-failures -I"$sieve_checks" \
    .github/development/prime-sieve-thread-check.cpp
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$native_checks/sieve-failures" --threads 2
Rscript --vanilla .github/development/prime-sieve-thread-check.R "$sieve_checks/r-bridge"
```

Also build/run this fixture with TSan separately. Its generated source already
contains the sieve and modified helper; no separate production `.cpp` is linked.
Coverage includes generation and count-only worker failures, partial creation,
join-before-propagation, recovery, and synthetic hardware-sizing guards.

The R bridges compile temporary cpp11 wrappers and check error translation plus
successful subsequent calls. They require installed cpp11/build dependencies.
Native exceptions are propagated only after joining; R APIs must still run on
the calling R thread.

### Package-wide launchers

```sh
package_threads=$(mktemp -d)
python3 .github/development/prepare-package-thread-check.py "$package_threads"
build_check package-failures -Wno-unused-parameter -I"$package_threads" \
    .github/development/package-thread-check.cpp
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$native_checks/package-failures" --threads 2
```

Repeat with the TSan build configuration above. The fixture retains the actual
16 launchers and gates real threads to check failure propagation and scratch
lifetimes; mathematical kernels are stand-ins. The family-specific R tests
above check real native/GMP results and attributes.
Recorded outcomes are in [DEFERRED_ISSUES.md](DEFERRED_ISSUES.md).

## P2 private policy and stream checks

```sh
python3 .github/development/p2-hybrid/production-integration/validate.py \
    /absolute/path/to/new-validation-output
```

The current-source runner snapshots the required sources/headers and runs
ASan/UBSan and TSan checks for independent P2/P3/count identities, the dynamic
queue, selector agreement, actual worker budgets, and private streaming windows.
It uses small retained reference fixtures, not an installed historical package
or the discarded variant-replay framework. See the
[validation README](p2-hybrid/production-integration/README.md) for the output
contract and [P2 findings](p2-hybrid/CONSOLIDATED.md) for policy rationale.

## R method routing and result conversion

```sh
Rscript --vanilla .github/development/prime-count-dispatch-check.R \
    /absolute/path/to/new-dispatch-output
```

This copies the real R wrapper, adapter and conversion sources, but substitutes
sentinel counters. It tests routing, hardware-capped budgets, validation order,
cpp11 translation and integer/numeric conversion. It does not test counting
mathematics or worker safety. The temporary bridge requires cpp11 and the
package's existing GMP build dependency.

## Native registration

When a registered signature changes, run `cpp11::cpp_register()` on an isolated
source copy. Review and retain the generated `src/cpp11.cpp`. The generator also
creates `R/cpp11.R`;
this package uses handwritten R wrappers and intentionally does not track that
file. Preserve any pre-existing local file or edits instead of deleting them
blindly. Update the handwritten wrapper, help, exports and arity checks together.

```sh
Rscript --vanilla - "$package_source" <<'RS'
cpp11::cpp_register(commandArgs(TRUE)[1], quiet = FALSE)
RS
```

Generation in the isolated copy is a verification step, not permission to
replace unrelated working-tree edits. Install and smoke-test the selected
wrapper/registration contents together.

## Performance measurements and the sieve inspector

Use optimized, unsanitized binaries. Keep setup and kernel timing distinct,
verify results, separate first and warmed calls using fresh processes, alternate
variant order, and finish builds/tests before timing. Requested workers may be
reduced; report effective budgets where available. Payload/capacity and process
peak RSS answer different memory questions.

Three small native drivers remain useful:

| Driver | Arguments | Linked sources |
| --- | --- | --- |
| `prime-count-check.cpp` | `n threads repetitions expected_count` | `PiLegendre.cpp PhiParallel.cpp Phi.cpp PiPrime.cpp Eratosthenes.cpp` |
| `prime-count-methods-benchmark.cpp` | `legendre\|meissel\|lehmer n threads repetitions expected_count` | `PiMeissel.cpp PiLehmer.cpp PiLegendre.cpp PhiParallel.cpp P2.cpp P3.cpp Phi.cpp Eratosthenes.cpp` |
| `partial-sieve-benchmark.cpp` | `p2\|p3\|reference\|verify x root_power threads repetitions expected` | `P2.cpp P3.cpp Phi.cpp Eratosthenes.cpp` |

Compile the selected driver from `.github/development/` with those `src/` inputs,
`-std=c++11 -O3 -pthread -Iinst/include`. The first two report whole-call timing;
the partial-sieve driver reports table setup separately. `reference` is a
straightforward serial-Phi P2 baseline, not an optimized external library.

The optional [sieve inspector](prime-sieve-inspector/README.md) creates an
isolated diagnostic package. It measures phases and storage estimates, not
hardware cache misses. Follow its current build recipe and interpret results
using [cache-inspector.md](cache-inspector.md). Retained studies and their limits
are indexed in [README.md](README.md); historical successful checks are evidence
for their recorded snapshots, not a substitute for rerunning the relevant tools.
