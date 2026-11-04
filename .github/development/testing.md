# Testing workflows

Run these commands from the repository root. Use a temporary library so an
investigation does not replace the normally installed package.

## Prime-counting source layout

- `src/PrimeCount.cpp` owns the registered R entry point, input conversion,
  small-input dispatch, thread arguments, and R result conversion.
- `src/PiPrime.cpp` implements the sieve counter used for small inputs.
- `src/PiLegendre.cpp` implements Legendre's formula and its serial or
  batched parallel summation. Each call prepares a local, immutable `PhiTables`
  object, and the exception-safe batches keep it alive until all workers join.
- `src/Phi.cpp` implements the general partial-sieve evaluator, periodic tiny
  tables, and shared atomic cache. `inst/include/NumbersUtils/Phi.h` declares
  `Phi` and the reusable tables, including a read-only prime sequence.
- `inst/include/NumbersUtils/PrimeCount.h` declares the two native entry points
  without R dependencies. It does not expose cache internals or phi helpers.

This split preserves the existing R interface and dispatch thresholds. It does
not introduce a method selector or a Lehmer implementation. Future algorithms
can use the shared phi infrastructure after their mathematical and architectural
design is agreed upon. The legacy `PhiTinyLookup.h` is retained but is no longer
included by production code; the shared evaluator owns the tiny tables.

### Source-split validation

The November 3, 2026 split was checked against baseline `a44c760`. Exact text
comparisons confirmed that the sieve, cache, phi recurrence, scheduler, and
setup bodies were preserved. The R adapter changed only its native callee name;
registration, R wrappers, dispatch thresholds, and result conversion did not
change. Default native arguments moved into the declaration header.

A clean package build completed without compiler warnings. `testPrimeCount.R`
(15 expectations), `testErrors.R`, and extra R dispatch-boundary/return-type
checks passed. The native checker passed ThreadSanitizer at four workers and
AddressSanitizer/UndefinedBehaviorSanitizer at one and four workers. Each
default run checked 1,172 small sieve counts, 662 Legendre counts, 76 signed
phi results, and 1,394,803 cache entries against independent references. Both
sanitizer builds also passed cache-allocation failure and recovery checks.

Native timing checks used `-std=c++11 -O3 -pthread`, five fresh processes per
version/input/worker combination, and one cold plus two warm calls per process.
Execution order was shuffled with seed 20261103 after other checks finished.
All 270 results matched. Median milliseconds were:

| Input | Workers | Before cold | Split cold | Before warm | Split warm |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `1e8` | 1 | 0.19150 | 0.17688 | 0.05567 | 0.04904 |
| `1e8` | 4 | 0.19737 | 0.19179 | 0.06575 | 0.06142 |
| `1e8` | 8 | 0.20913 | 0.20638 | 0.08631 | 0.08206 |
| `1e12` | 1 | 77.320 | 69.993 | 47.317 | 43.226 |
| `1e12` | 4 | 25.973 | 23.568 | 16.140 | 15.279 |
| `1e12` | 8 | 23.117 | 21.106 | 15.422 | 14.138 |
| `1e13` | 1 | 425.334 | 383.244 | 370.642 | 335.266 |
| `1e13` | 4 | 140.947 | 126.890 | 118.577 | 107.575 |
| `1e13` | 8 | 121.513 | 111.865 | 108.701 | 98.564 |

No regression appeared in this timing check. Lower medians do not identify the
cause: source organization and helper linkage changed together, and compiler
optimization effects were not isolated. Validation used arm64 macOS; a full
package check and other platforms were not run for this refactor.

## Shared Phi integration

Legendre prepares one `PhiTables` instance per call and shares it by const
reference with its workers. Its starting value is `Phi(n, 6, tables)`, and each
original foreman interval subtracts `Phi(n / p[i + 1], i, tables)`. This is the
same telescoping recurrence as the original negative-sign worker. For the
Legendre domain `n >= 100000`, the prepared root table contains at least the
first six primes. Its prime count supplies `a`; root preparation itself retains
the previous floating-point conversion. The general Phi evaluator still uses
exactly corrected roots internally.

The old Legendre-specific recurrence, tiny lookup use, mutable prime/pi globals,
and duplicate cache were removed. The batch loops, interval formulas, thread
thresholds/clamping, result order, and exception handling were preserved.
Temporary scheduler probes retained the actual setup and threads, replacing
only foreman computation with logging. Before/after logs matched exactly in
65 cases: 3,389 batches, 12,998 launches and worker receipts, 23 serial
intervals, and 17 partial batches. Logs included batch region/capacity, launch
order, destination result slot, bounds, and initial arithmetic. Inputs covered
threshold neighbors and values through `1e15`, requested workers 1/2/3/4/8/10,
and clamping/serial controls. This verifies scheduling, not the counting math.

The integrated native checker passed ASan/UBSan at 1/4/8 workers and TSan at
4 workers, plus both cache-allocation failure modes. Each default run verified
662 independent Legendre counts, 1,172 small sieve counts, 76 signed phi
results, two concurrent native Legendre calls with different tables, and
1,425,075 populated cache cells. The standalone general-Phi checker also passed
both sanitizers with and without caching: 2,621,310 exhaustive cases, 165
targeted cases, 220 large-input cases, 3,328 concurrent evaluations, 57 invalid
inputs, and known prime counts through `1e12` (through `1e8` without caching).

The existing failure harness passed again with shared Phi: 28 native scenarios,
24 blocked-join checks, 32 successful counts, and 20 R errors followed by 20
successful recoveries. A clean package installation had no compiler warnings;
`testPrimeCount.R` and `testErrors.R` passed. Native builds explicitly used
C++11 on Apple clang/arm64 macOS. A full package check and other platforms were
not run.

### Integration performance investigation

Direct integration has a material performance regression. A diagnostic timing
comparison used `-std=c++11 -O3 -pthread`, three fresh processes per
input/worker/version, one cold and two warm calls per process, and shuffled
process pairs/version order with seed 33. It compared the baseline with its
own Legendre recurrence, the direct integration, and three temporary variants.
All 405 counts matched. Timing included table preparation and allocation, and
ran after other builds/checks completed. Temporary variants were not applied
to production. Before/direct-integration medians in milliseconds were:

| Input | Workers | Before cold | Shared cold | Before warm | Shared warm |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `1e8` | 1 | 0.15883 | 0.21017 | 0.04685 | 0.06348 |
| `1e8` | 4 | 0.19092 | 0.23175 | 0.06598 | 0.08019 |
| `1e8` | 8 | 0.20483 | 0.24883 | 0.08577 | 0.10008 |
| `1e10` | 1 | 2.80229 | 3.58738 | 1.06065 | 1.35804 |
| `1e10` | 4 | 2.07304 | 2.70346 | 0.73527 | 0.93921 |
| `1e10` | 8 | 2.00571 | 2.55829 | 0.66467 | 0.86629 |
| `1e12` | 1 | 68.56250 | 90.65433 | 43.16560 | 55.59440 |
| `1e12` | 4 | 23.77713 | 29.19937 | 15.98650 | 18.67583 |
| `1e12` | 8 | 21.94583 | 26.76554 | 14.01560 | 17.44877 |

Warm direct-integration calls were roughly 17-35% slower in this comparison.
Per-call tables allocate/initialize pi storage afresh, whereas the previous
global pi vector reused capacity. General Phi also has extra validation,
initialization guards, and evaluator calls. Temporary product-based exact-root
checks removed two integer divisions per recursive expansion but recovered
little of the large-input gap; restoring the old style of pi-table filling
helped modestly. At `1e12`/one worker, warm medians were 43.166 ms for baseline,
55.594 ms for direct integration, and 54.975 ms with both diagnostic changes.
These measurements do not isolate the dominant cause. Further hot-path
optimization requires discussion before expanding the integration scope.

## Legendre scheduling bounds

The excessive-worker fix caps the effective worker count at the number of
terms in `[6, a)`, where `a` is the prepared root table's prime count. This
keeps the small-input chunk positive. Previously, `PiLegendre(10000000,
442, 442)` included index 5 and returned 504739 instead of 664579; larger
effective worker requests had the same problem. The R hardware-thread cap
can mask this defect on smaller machines.

Large-input first batches now stop at `a`, including a partial final batch,
and return the accumulated result if they have completed the whole sum.
Every launched worker is joined before that return. Middle batches also
bound their endpoints explicitly rather than relying on `pow` and repeated
multiplication rounding identically. The original chunk formulas, worker
exception handling, and ordinary schedules are retained.

The standalone scheduler checker instruments a temporary production-source
copy and substitutes synchronous workers. It verifies contiguous, once-only
coverage of `[6, a)`, table bounds, worker limits, and completed joins without
requiring hundreds of OS threads. Known-count cases execute the actual Phi
evaluator; schedule-only cases skip that computation. This fixture checks
partitioning and lifetime bookkeeping, not concurrent execution; the existing
cache and worker-failure checks exercise real threads.

```sh
schedule_check=$(mktemp -d)
python3 .github/development/prepare-prime-count-schedule-check.py "$schedule_check"
clang++ -std=c++11 -O1 -g -pthread -fsanitize=address,undefined \
    -fno-omit-frame-pointer -Iinst/include -I"$schedule_check" \
    .github/development/prime-count-schedule-check.cpp \
    src/Phi.cpp src/Eratosthenes.cpp -o "$schedule_check/check"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$schedule_check/check"
```

The final checker passed ASan/UBSan with 82 cases, 14 real counts, and 63,946
bounded contiguous intervals, including effective-worker boundary neighbors
and requests through `INT_MAX`. Before-change probes rejected `[6, 5)` at
`1e7`/1000 workers and `[9589, 9596)` beyond `a = 9592` at `1e10 + 1`/2000
workers. Dumps from 53 ordinary or hardware-clamped cases matched before and
after: 708 records including batch, result slot, and interval boundaries.

The production cache checker and all 28 worker/launch-failure scenarios also
passed ASan/UBSan with real workers. A clean temporary package build, followed
by rebuilding the final middle-batch guard, completed without compiler warnings;
`testPrimeCount.R` and `testErrors.R` passed against the final build. Validation
used Apple clang on arm64 macOS, with C++11 for the native checks. A full
package check, other platforms, and performance benchmarks were not run for
this focused correctness change.

## Package installation and focused tests

The prime-count cache fix was built from a temporary copy of the package source,
without pre-existing object files. `package_source` below is that copy; using
the working tree instead will update its ignored build artifacts.

```sh
package_library=$(mktemp -d)
MAKEFLAGS=-j4 R CMD INSTALL --preclean --no-multiarch \
    --library="$package_library" "$package_source"
Rscript --vanilla - "$package_library" <<'RS'
library(RcppAlgos, lib.loc = commandArgs(TRUE)[1])
testthat::test_file("tests/testthat/testPrimeCount.R",
                    reporter = "summary", stop_on_failure = TRUE)
testthat::test_file("tests/testthat/testErrors.R",
                    reporter = "summary", stop_on_failure = TRUE)
RS
```

The parallel cases precede the serial cases in `testPrimeCount.R`. A standalone
run in a new R process therefore tests an initially empty phi cache. Running the
whole suite may populate it earlier. Count comparisons do not detect data races;
use ThreadSanitizer separately.

## Native prime-count sanitizer checks

This checks the production sieve counter and Legendre implementation, together
with `Eratosthenes.cpp`, without the R adapter. No source extraction or rewriting
is needed. The standalone checks call the native entry points directly; use
supported inputs and realistic thread counts.

```sh
prime_count_check=$(mktemp -d)
clang++ -std=c++11 -O1 -g -pthread -fsanitize=thread \
    -fno-omit-frame-pointer -Iinst/include \
    .github/development/prime-count-cache-check.cpp \
    src/PiLegendre.cpp src/PiPrime.cpp src/Eratosthenes.cpp \
    -o "$prime_count_check/prime-count-tsan"
TSAN_OPTIONS=halt_on_error=1 \
    "$prime_count_check/prime-count-tsan" 4
TSAN_OPTIONS=halt_on_error=1 \
    "$prime_count_check/prime-count-tsan" 4 --allocation-failure
```

The checker directly includes `src/Phi.cpp` in its translation
unit, so it can inspect every populated cache entry without exposing production
test hooks. Do not also list that source as a separate compiler input for this
checker; it would define the shared evaluator twice. Link the Legendre source
separately so the same checker exercises its production integration.
Its optional second argument replaces page allocations with failures in the
standalone executable. The package itself does not replace allocation functions.

The default mode compares `PiPrime` with an independent sieve at small inputs,
prime-square neighbors, segment boundaries, and the upper end of its R dispatch
range. It also checks both signs of the phi evaluator against independent partial
sieves around its direct-table transition, with cache publication disabled.
It then runs known counts in increasing/decreasing order, compares small counts
and prime-square neighbors with an independent Eratosthenes sieve, and compares
every populated phi cell with independent partial-sieve counts.
Two concurrent native Legendre calls also exercise separate prepared tables
and shared cache use. The failure mode checks the result with no successful
page allocations, verifies that every bounded cache key remains a miss, then
checks recovery, signs, zero entries, page boundaries, and the final valid key.

AddressSanitizer and UndefinedBehaviorSanitizer use a separate executable:

```sh
clang++ -std=c++11 -O1 -g -pthread -fsanitize=address,undefined \
    -fno-omit-frame-pointer -Iinst/include \
    .github/development/prime-count-cache-check.cpp \
    src/PiLegendre.cpp src/PiPrime.cpp src/Eratosthenes.cpp \
    -o "$prime_count_check/prime-count-asan-ubsan"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$prime_count_check/prime-count-asan-ubsan" 4
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$prime_count_check/prime-count-asan-ubsan" 4 --allocation-failure
```

The original implementation produced a ThreadSanitizer data-race report during
vector growth on a cold `10^10`, four-thread call. On the macOS toolchain that
report included
external-symbolizer warnings; `TSAN_OPTIONS=halt_on_error=1:symbolize=0` disables
symbolization without disabling race detection.

The integrated atomic cache passed both modes with ThreadSanitizer at 1, 2, 3,
4, and 10 threads, and with AddressSanitizer/UndefinedBehaviorSanitizer at four
threads. Each default-mode run checked 662 independent prime counts and
1,394,803 populated phi cells. Each failure-mode run rejected 646 page
allocations while retaining the correct result. There were no sanitizer reports.

## Legendre worker exception handling

Each existing launch region now uses a local `PhiThreadBatch`. It allocates one
exception slot per possible worker and reserves thread storage before starting
any worker. Each worker catches native C++ exceptions into its own slot;
`finish()` joins every worker before inspecting those slots or allowing partial
sums to be combined. It rethrows the first stored exception in launch order on
the calling thread. When called through the registered R interface, cpp11 then
translates that exception into an R error.

The batch destructor also joins any still-joinable workers. If a later launch
fails, earlier workers therefore finish before the launch exception propagates
and before their result storage is destroyed. A launch exception takes
precedence over any worker exceptions in that incomplete batch. The helper
does not call R, replace cpp11's error handling, reuse threads, or change the
existing intervals, batch boundaries, recurrence, or result order. General
The batch now holds a reference to its call's `PhiTables`; those tables and the
result vectors outlive the batch and all its workers.

The standalone launcher checker injects worker and thread-creation failures in
a temporary source copy. Its preparer verifies the expected source structure,
then substitutes a test thread wrapper, adds a worker-entry hook, and labels
the four launch regions. Production code has no fault-injection hooks.

```sh
thread_check=$(mktemp -d)
python3 .github/development/prepare-prime-count-thread-check.py "$thread_check"
clang++ -std=c++11 -O1 -g -pthread -fsanitize=address,undefined \
    -fno-omit-frame-pointer -Iinst/include -I"$thread_check" \
    .github/development/prime-count-thread-check.cpp src/Phi.cpp src/Eratosthenes.cpp \
    -o "$thread_check/prime-count-thread-asan-ubsan"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$thread_check/prime-count-thread-asan-ubsan"
clang++ -std=c++11 -O1 -g -pthread -fsanitize=thread \
    -fno-omit-frame-pointer -Iinst/include -I"$thread_check" \
    .github/development/prime-count-thread-check.cpp src/Phi.cpp src/Eratosthenes.cpp \
    -o "$thread_check/prime-count-thread-tsan"
TSAN_OPTIONS=halt_on_error=1:symbolize=0 \
    "$thread_check/prime-count-thread-tsan"
Rscript --vanilla .github/development/prime-count-thread-check.R \
    "$thread_check/r-bridge"
```

Both sanitizer builds passed 28 injected-failure scenarios across all four
launch regions, with 24 gated checks proving the caller waits for workers to
finish and join before propagation. Cases include standard and nonstandard
worker exceptions, `std::bad_alloc`, failures before any launch and after
partial creation, simultaneous worker/launch failures, and deterministic
selection by launch slot even when exceptions occur in reverse order. Every
scenario is followed by a successful known-count recovery; there are 32
successful counts including baseline and serial controls.

The small, first, and middle regions use `1e10` or `1e10 + 1`; the final region
uses `1e12 - 1` to exercise a two-worker tail in a four-worker configuration.
An initial fixture assertion correctly rejected the assumption that
`1e10 + 1` had a two-worker tail: it actually has one. Only the test fixture
changed. The checker asserts its observed region sizes before injecting faults.

The R companion builds a temporary cpp11 wrapper around the same instrumented
Legendre source and checks 20 errors (five categories in each launch region).
It verifies error messages, completed joins, no later batch after failure, and
a correct subsequent count after every error. All 20 error/recovery pairs
passed in a fresh build without compiler warnings. The first bridge build
exposed a missing-return warning when the native checker's `main` was renamed;
an explicit successful return fixed that test-only issue. R selected
`-std=gnu++20` in this environment despite the `cpp_source` C++11 request; the
standalone sanitizer builds explicitly used `-std=c++11`. No package build
configuration or required language standard was changed.

The November 3, 2026 change passed the existing native checker with
AddressSanitizer/UndefinedBehaviorSanitizer at 1, 4, and 8 workers and
ThreadSanitizer at 4 workers. Each default run checked 1,172 small sieve counts,
662 Legendre counts, 76 signed phi results, and 1,394,803 populated cache cells.
Both sanitizer builds also passed cache-allocation failure and recovery checks.
A clean temporary package installation produced no compiler warnings, and
`testPrimeCount.R` (15 expectations) and `testErrors.R` passed. These checks
used C++11 for native builds and Apple clang on arm64 macOS. A full package
check and other platforms were not run.

### Exception-handling timing check

Before/after binaries used `-std=c++11 -O3 -pthread`. The baseline was the
staged source split immediately before adding `PhiThreadBatch`. Each
input/worker/version combination used five fresh processes with one cold and
two warm calls. Process pairs and version order were shuffled with seed
20261103; builds and sanitizer checks had finished before timing began. All
450 results matched known counts. Median times in milliseconds were:

| Input | Workers | Before cold | After cold | Before warm | After warm |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `1e8` | 1 | 0.16462 | 0.17133 | 0.04860 | 0.04846 |
| `1e8` | 4 | 0.20037 | 0.19346 | 0.06327 | 0.06054 |
| `1e8` | 8 | 0.22596 | 0.22387 | 0.08052 | 0.08150 |
| `1e10` | 1 | 3.13379 | 3.07254 | 1.17931 | 1.13777 |
| `1e10` | 4 | 2.27188 | 2.26737 | 0.78863 | 0.80537 |
| `1e10` | 8 | 2.06121 | 2.07488 | 0.66738 | 0.66142 |
| `1e10 + 1` | 1 | 2.88183 | 2.89913 | 1.07110 | 1.07258 |
| `1e10 + 1` | 4 | 1.74571 | 1.80625 | 0.84027 | 0.85102 |
| `1e10 + 1` | 8 | 1.46412 | 1.50871 | 0.76106 | 0.75625 |
| `1e12` | 1 | 71.34300 | 69.49175 | 43.59158 | 44.00892 |
| `1e12` | 4 | 23.74208 | 23.62504 | 15.84662 | 15.75675 |
| `1e12` | 8 | 21.75829 | 21.67425 | 14.39752 | 14.14800 |
| `1e13` | 1 | 384.43608 | 385.43658 | 335.98429 | 337.95840 |
| `1e13` | 4 | 126.61146 | 127.52312 | 108.88590 | 109.51781 |
| `1e13` | 8 | 112.97242 | 114.63958 | 99.09154 | 99.51102 |

Warm parallel changes ranged from about -4.3% to +2.1%, while serial controls
also varied. These measurements show no consistent performance regression;
they do not establish a speedup or isolate the small helper overhead.

## Atomic-cache contract

The shared cache preserves the original bounds (`a < 100`, `x <= 65535`) and
reuse across worker batches and calls. Each row owns up to sixteen fixed blocks
of 4,096 atomic cells, allocated on demand. Cells and pointers are explicitly
initialized for C++11. Published blocks remain alive and stationary until cache
destruction, which requires all workers to have joined.

The row mutex protects allocation only. Release/acquire publication makes block
initialization visible to readers. Cell loads/stores are relaxed: zero is a
miss, and each positive value is the complete phi answer. Duplicate computations
store the same value; a missed cache value only causes recomputation. Failed
page allocation skips the cache write and leaves the calculated answer intact.

Following shared-Phi integration, Legendre's prepared prime and pi tables belong
to each call; the atomic cache is the shared mutable state. Concurrent native
Legendre calls with different inputs are covered by the sanitizer checker.
This does not authorize R API calls from worker threads. Excessive-worker
interval bounds are addressed separately under "Legendre scheduling bounds".

On the tested arm64 macOS environment, both pointer and 16-bit atomics were
lock-free and a cell occupied two bytes. Correctness does not require lock-free
atomics; their size and performance may differ elsewhere.

## Native timing checks

Compile the production native sources with the separate timing driver:

```sh
clang++ -std=c++11 -O3 -pthread -Iinst/include \
    .github/development/prime-count-check.cpp \
    src/PiLegendre.cpp src/Phi.cpp src/PiPrime.cpp src/Eratosthenes.cpp \
    -o "$prime_count_check/prime-count-timing"
"$prime_count_check/prime-count-timing" 1000000000000 10 4 37607912018
```

Arguments are `n threads repetitions expected_count`. Each invocation starts a
fresh process; repetition zero is cold, and later repetitions reuse the cache.
Output columns are `n,threads,repetition,count,milliseconds`. Timings include
native setup and cache allocation, but exclude R and process startup. Each
result is checked against the expected count. Sanitizer timings must not be
used as performance results. Keep before/after binaries, alternate their
execution order, and compare cold and warm calls separately.

The November 2, 2026 integration check used seven fresh processes per version
and thread count, each with one cold and three warm calls at `10^12`. Version
and thread-count order were shuffled with seed 20261102. All 168 results
matched. Median native times in milliseconds were:

| Workers | State | Original | Atomic cache |
| ---: | --- | ---: | ---: |
| 1 | cold | 85.135 | 85.724 |
| 1 | warm | 51.137 | 52.967 |
| 4 | cold | 102.453 | 28.635 |
| 4 | warm | 18.037 | 17.961 |
| 10 | cold | 161.335 | 24.433 |
| 10 | warm | 15.525 | 16.296 |

The original implementation is an unsafe timing baseline, not a correctness
reference. Cold parallel calls improve substantially; warmed single/ten-worker
calls retain a small overhead. Earlier prototype measurements also found a
small first-call allocation cost at `10^5` (about 19 to 41 microseconds).

The package rebuilt successfully under R 4.6.1, and both focused R test files
above passed. A full package check was not run for this localized internal
change. Validation used Apple clang 17 on arm64 macOS; other platforms have not
been checked in this investigation.

## Worker-reuse experiment (not retained)

On November 3, 2026, a per-call `PhiBatchPool` prototype reused workers for
parallel inputs above `1e10`, preserving the existing intervals and batch
barriers. Each worker received one interval per batch and slept between
batches using condition variables. A mutex protected dispatch and completion;
the counting work ran outside that mutex. Workers were joined at call exit.
Serial counting and the single-batch path at or below `1e10` were unchanged.

Decision: retain the original create/join-per-batch implementation. The small,
inconsistent timing differences below did not justify the extra synchronization
and lifecycle code. The prototype was removed; the atomic cache, 32-bit pi
table, setup improvements, and batch-boundary regression test were retained.
A fresh package build and both focused R test files passed again after
restoring the pre-pool source.

The prototype passed the existing native checker with
ThreadSanitizer at 2, 3, 4, 8, and 10 workers and AddressSanitizer/
UndefinedBehaviorSanitizer at four workers. Each run checked 662 independent
prime counts and 1,394,803 populated cache entries. A clean package build had
no compiler warnings; `testPrimeCount.R` (15 expectations) and `testErrors.R`
passed. The retained R regression crosses `1e10 - 1`, `1e10`, and `1e10 + 1` and
repeats the multiple-batch call with available worker counts from 2 through 4.

Two temporary diagnostics provided additional evidence:

- A compiled comparison replaced `phiForeman` with interval logging and
  compared the old and new schedules, sorting within each batch and excluding
  empty slots. All 3,398 batches and 13,027 nonempty intervals matched in 40
  cases: threshold neighbors and powers of ten through `1e15`, with 2, 3, 4,
  8, and 10 workers. These cases also had contiguous interval coverage.
- A harness extracting the pool class tested 804 batches, empty slots,
  once-only slot execution, worker exceptions followed by reuse, and idle
  shutdown. A modified diagnostic injected 100 constructor failures after
  zero through three successful launches. ThreadSanitizer reported no races.
  The fault injection was confined to the diagnostic, not production code.

These checks used Apple clang on arm64 macOS and C++11 for the native builds.
The experiment did not include a full package check or testing on other
platforms. It did not address the previously identified bounds problems at
excessive effective worker counts.

### Worker-reuse timings

Before/after native binaries used the same compiler and `-std=c++11 -O3
-pthread`. The baseline (`6cff583`) already contained the atomic cache, 32-bit pi table,
prime-vector move, and reuse of the setup count. Execution order was shuffled
with seed 20261103, with no other investigation builds or diagnostics running.
Calls included setup and pool destruction, but excluded process/R startup.

At `1e10`, `1e10 + 1`, `1e12`, and `1e13`, each version/worker combination
used seven fresh processes with one cold and three warm calls. The `1e14`
cases used three processes with one cold and two warm calls. The serial
`1e12` control used five processes with one cold and three warm calls.
All 524 counts matched. Median times in milliseconds were:

| Input | Workers | Before cold | Pool cold | Before warm | Pool warm |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `1e10` (control) | 4 | 2.524 | 2.349 | 0.944 | 0.889 |
| `1e10` (control) | 8 | 2.333 | 2.394 | 0.784 | 0.770 |
| `1e10 + 1` | 4 | 1.922 | 1.907 | 0.997 | 0.993 |
| `1e10 + 1` | 8 | 1.557 | 1.570 | 0.843 | 0.868 |
| `1e12` (serial control) | 1 | 90.300 | 88.509 | 55.134 | 54.359 |
| `1e12` | 4 | 28.685 | 27.915 | 18.804 | 18.385 |
| `1e12` | 8 | 25.131 | 25.194 | 16.900 | 16.552 |
| `1e13` | 4 | 169.418 | 165.382 | 145.663 | 144.065 |
| `1e13` | 8 | 141.419 | 138.838 | 126.845 | 127.900 |
| `1e14` | 4 | 1142.360 | 1103.988 | 1130.194 | 1095.499 |
| `1e14` | 8 | 935.681 | 941.031 | 934.904 | 944.292 |

The launch reduction was substantial: for example, `1e14` with eight workers
used eight threads instead of 875 launches across 110 batches. Elapsed-time
benefits were small and inconsistent, however. Four-worker `1e14` calls
improved about 3%, while eight-worker calls were about 1% slower. Unchanged
control paths also varied, so small differences should not be attributed
confidently to the pool. The prototype passed correctness checks, but did not
demonstrate a compelling performance advantage on this machine.

## Monotonic direct-table suffix

The November 3, 2026 loop change splits `phiWorker` into a recursive prefix
and a direct-table suffix. The prefix stops at the first successful `isPix`
check; the suffix includes that same term and all remaining terms. No terms
are reordered or omitted, and initialization, signs, cache behavior, and
thread scheduling are unchanged. The simple two-loop implementation repeats
one division at the transition rather than carrying extra state between loops.

For a fixed frame, let `q = floor(x / p[i + 1])` and let `L` be the pi-table
length. As `i` increases, `q` cannot increase and `p[i + 1]^2` increases.
Thus once `q < L` and `q < p[i + 1]^2` hold, both continue to hold. This
justifies using the direct formula unconditionally for the remaining terms.
The equivalent thresholds are `x < p^3` and `x < p * L`; production code
retains the quotient comparisons to avoid introducing larger products.

The native checker adds 76 signed phi comparisons using independently built
tables for parent input `1e8`. Cases cover `p^3 - 1`, `p^3`, and `p^3 + 1`
for primes 41, 43, and 47, with cutoffs before, at, and after each transition.
It also checks `127 * L` and its neighbors, where table membership controls
the transition, plus an entirely direct loop and an early-direct degenerate
case. Cache publication is disabled for these checks so neither sign can
bypass the changed loop through a previously computed cache value.

A clean package build had no compiler warnings, and `testPrimeCount.R`
(15 expectations) and `testErrors.R` passed. The extended native checker
passed ThreadSanitizer with four workers and AddressSanitizer/
UndefinedBehaviorSanitizer with 1, 4, and 8 workers. Each default run checked
662 independent prime counts, 76 signed phi results, and 1,394,803 populated
cache entries. Both sanitizer builds also passed cache-allocation failure and
recovery checks. These checks used C++11 and Apple clang on arm64 macOS;
a full package check and other platforms have not been tested for this change.

The native timing comparison used baseline `4c487bb` and the two-loop version,
both built with `-std=c++11 -O3 -pthread`. Execution order was shuffled with
seed 20261103 after builds and diagnostics finished. Inputs through `1e13`
used seven fresh processes per version/worker combination, each with one cold
and three warm calls. The `1e14` cases used three fresh processes with one cold
and two warm calls. Timings included native setup and excluded R/process
startup. All 764 counts matched. Median milliseconds were:

| Input | Workers | Before cold | Suffix cold | Before warm | Suffix warm |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `1e5` | 1 | 0.03658 | 0.03908 | 0.00083 | 0.00083 |
| `1e8` | 1 | 0.17821 | 0.18804 | 0.05717 | 0.05433 |
| `1e8` | 4 | 0.20354 | 0.20567 | 0.06458 | 0.06175 |
| `1e8` | 8 | 0.22117 | 0.23613 | 0.08533 | 0.08950 |
| `1e10` | 1 | 3.137 | 3.311 | 1.256 | 1.244 |
| `1e10` | 4 | 2.379 | 2.334 | 0.897 | 0.834 |
| `1e10` | 8 | 2.281 | 2.245 | 0.777 | 0.725 |
| `1e12` | 1 | 84.695 | 78.300 | 52.756 | 47.207 |
| `1e12` | 4 | 28.112 | 26.165 | 18.263 | 16.475 |
| `1e12` | 8 | 25.230 | 23.288 | 16.702 | 15.177 |
| `1e13` | 1 | 475.864 | 427.472 | 415.801 | 370.456 |
| `1e13` | 4 | 157.516 | 139.549 | 135.088 | 121.667 |
| `1e13` | 8 | 135.592 | 124.833 | 120.819 | 108.320 |
| `1e14` | 4 | 1060.480 | 921.348 | 1054.790 | 922.452 |
| `1e14` | 8 | 944.286 | 832.732 | 926.667 | 832.038 |

Warmed calls from `1e12` through `1e14` improved by approximately 9–13% in
these measurements; cold calls improved by approximately 7–13%. Smaller
workloads were mixed: some cold medians increased by 5–7%, including serial
`1e10` (3.137 to 3.311 ms). These results do not establish a universal speedup
or isolate the cause of the smaller regressions.

## Native registration

No registration changes are needed for internal prime-counting fixes. When a
registered interface changes, the verified workflow is
`cpp11::cpp_register(".")`; see [cache-inspector.md](cache-inspector.md)
for this repository's handling of the generated R wrappers.
