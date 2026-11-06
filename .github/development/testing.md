# Testing workflows

Run these commands from the repository root. Use a temporary library so an
investigation does not replace the normally installed package.

## Prime-counting source layout

- `R/NumberTheory.R` exposes `primeCount(n, nThreads = NULL, method = "auto")`
  and matches the method argument. `src/PrimeCount.cpp` owns the registered R
  entry point, input conversion, method dispatch, thread arguments, and result
  conversion.
- `src/PiPrime.cpp` implements the sieve counter used for small inputs.
- `src/PiLegendre.cpp` prepares a local, immutable `PhiTables` object, clamps
  the requested workers, and applies Legendre's formula using `PhiParallel`.
- `src/PiMeissel.cpp` and `src/PiLehmer.cpp` prepare square-root tables once,
  evaluate Phi at the cube-root or fourth-root cutoff, and subtract their
  P2 or P2/P3 corrections. They reuse the existing helpers and worker policies.
- `src/PhiParallel.cpp` evaluates the outer Phi sum over the requested cutoff
  `a`, independently of the table extent. It retains the existing Legendre
  interval formulas, batch barriers, and serial/parallel thresholds; their
  performance with other cutoffs has not yet been tuned.
- `inst/include/NumbersUtils/ThreadBatch.h` owns one batch's threads and
  exception slots. Callables supply the work; captured tables and result
  storage must outlive the batch, which joins before propagating failures.
- `src/Phi.cpp` implements the general partial-sieve evaluator, periodic tiny
  tables, and shared atomic cache. `inst/include/NumbersUtils/Phi.h` declares
  `Phi`, `PhiParallel`, and the reusable tables, including a read-only prime
  sequence. Worker evaluation through `Phi` remains serial.
- `src/P2.cpp` counts two-prime products using direct pi lookups and shared
  serial Phi for larger quotients. `src/P3.cpp` counts three-prime
  products with direct lookups. `NumbersUtils/PartialSieve.h` declares these
  internal kernels; `PartialSieveUtils.h` shares validation and corrected roots.
  Neither kernel changes the R driver or selects a prime-counting method.
- `inst/include/NumbersUtils/PrimeCount.h` declares the four native counters
  without R dependencies. It does not expose cache internals or phi helpers.

The R method selector exposes Legendre, Meissel and Lehmer. The default `auto`
policy preserves the historical small-input dispatch and Legendre cutoff;
broader performance and memory measurements are needed before tuning that
policy. The legacy `PhiTinyLookup.h` is retained but is no longer included by
production code; the shared evaluator owns the tiny tables.

## R prime-counting method selection

`primeCount(n, nThreads = NULL, method = "auto")` appends the new argument so
existing positional thread arguments keep their meaning. `match.arg` accepts
`auto`, `legendre`, `meissel`, and `lehmer`, with case-sensitive unique prefixes;
`le` is ambiguous and errors, while `leg` and `leh` distinguish the two methods.
As with standard `match.arg`, `NULL` selects the first choice, `auto`.

The current automatic policy is the previous implementation: a tiny lookup
below 10, the sieve counter from 10 through 99,999, and Legendre from 100,000
upward. Explicit choices call the named native algorithm for every accepted
positive input, including small values. Legendre's native positive-input
contract is now documented as `1 <= n <= 2^53 - 1`; its algorithm is unchanged.
The automatic policy may be tuned after broader measurement and does not
promise the fastest method for every input or machine.

Input validation and result types retain the existing contract: `n` must be a
positive whole number at most `2^53 - 1`, and the result is integer when the
count fits in R's integer range and double otherwise. `nThreads = NULL` supplies
one requested worker. All methods ignore `nThreads` below 100,000 and execute
serially, preserving the old validation threshold. At larger inputs, positive
whole thread requests are converted as before, capped by the detected maximum,
and can be reduced further by individual kernels. Result allocation now uses
`cpp11::as_sexp` so R allocation errors unwind C++ objects, including the new
method string, safely.

The registered `PrimeCountCpp` signature has a fourth argument containing the
matched method string. Its native validation also rejects unknown names.
Registration was regenerated with `cpp11::cpp_register()` in a temporary
package copy; only the generated `src/cpp11.cpp` was copied back. The repository
continues to use its handwritten R wrappers, with no generated `R/cpp11.R`.

`tests/testthat/testPrimeCount.R` includes an independent base-R sieve, cutoff
neighbors, known counts and result types, positional compatibility, method
matching, parallel calls, and validation/recovery. Existing parallel cache
tests remain first so focused runs retain their cold-cache coverage.

To verify dispatch independently of the algorithms agreeing on a count, run:

```sh
Rscript --vanilla .github/development/prime-count-dispatch-check.R
```

An optional output-directory argument preserves the generated build in a
chosen temporary location. This fixture compiles unchanged copies of the real
adapter and converters against sentinel counters, then evaluates the unchanged
R wrapper through a local bridge to the generated cpp11 entry point. It checks
which counter was called, forwarded thread arguments, integer/double conversion
at `INT_MAX`, invalid inputs, and C++ exception translation/recovery. It does
not validate counting mathematics, actual worker scheduling, or worker failures;
the independent counting and thread-failure checks cover those separately.

### API validation

A clean isolated installation passed 540 expectations across the ten focused
number-theory/error test files: 115 in `testPrimeCount.R`, 207 in `testErrors.R`,
and 218 across prime sieving, factorization, primality, Euler phi, and divisor
tests. No failures, errors, warnings, or skips were reported. The sentinel
dispatch fixture passed all 138 checks. Installed formals match the documented
signature, native registration reports four arguments, all five help examples
run, and `tools::codoc()` reports no package-wide documentation mismatches.
Both updated Rd files passed parsing and `tools::checkRd`.

Validation used arm64 macOS and R 4.6.1. The R build selected its configured
compiler standard; the dispatch fixture is an R/cpp11 integration check, not
a separate C++11 portability or sanitizer run. Existing numerical kernels,
cache operations and worker scheduling were not modified, and their earlier
sanitizer evidence remains separate. No full `R CMD check`, other-platform
run, or new timing/peak-memory study was performed for this API integration.

Commands, build/test logs, the unchanged source copies used by the sentinel
fixture, and source hashes are retained under the directory recorded by
`/tmp/rcppalgos-prime-api-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-prime-api-loxkonkr`.

## Source-split validation

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

## Meissel and Lehmer native drivers

`PiMeissel(n, nThreads, maxThreads)` and `PiLehmer(n, nThreads, maxThreads)`
accept integer `0 <= n <= 2^53 - 1`. Invalid inputs throw `std::out_of_range`
before preparation; zero and one return zero immediately. Both thread controls
default to one. Their minimum supplies the phase budget, and budgets at most
one select serial execution. The native integration left Legendre unchanged;
the R selector described above now exposes all three counters.

Each call prepares one immutable `PhiTables(floor(sqrt(n)))` and borrows it
across sequential phases. Meissel chooses `a = pi(floor(cbrt(n)))` and returns
`PhiParallel(n,a) + a - 1 - P2(n,a)`. Lehmer chooses
`a = pi(floor(n^(1/4)))` and additionally subtracts `P3(n,a)`. Corrected integer
roots supply both cutoffs: `floor(sqrt(floor(sqrt(n))))` is exactly the fourth
root floor, so no new root helper or floating fourth-root estimate is needed.

The identities follow by partitioning the survivors counted by Phi. They
include one, the `pi(n) - a` surviving primes, and composite survivors. After
the cube-root cutoff every surviving factor exceeds `cbrt(n)`, so a composite
has exactly two prime factors. After the fourth-root cutoff it has two or
three. Repeated factors are included by P2/P3. At the latter cutoff the
smallest surviving prime `p` satisfies `p^2 > sqrt(n)`; therefore
`floor(n/p^2) <= floor(sqrt(n))`, covering all P3 lookups in the shared table.
The first cutoff cannot exceed the prepared square-root limit either.

`PhiParallel`, P2 and P3 finish their workers before the next phase starts.
Their serial worker kernels avoid nested parallel launches, and their existing
`ThreadBatch` lifetimes propagate exceptions only after joining. All phases
share prepared tables; Phi and P2 also reuse the bounded atomic Phi cache.
Dense table memory remains `O(sqrt(n))`. The native drivers added no helper,
dependency, or new scheduling policy. Their later R integration extends the
existing registered entry point with a method argument.

### Memory requirements

Legendre, Meissel and Lehmer each prepare one `PhiTables` context per call.
For a prepared table limit `L`, the dense `int32_t` pi array contains `L + 1`
elements, requiring `4 * (L + 1)` bytes of element storage. The
`int64_t` prime vector, shared Phi cache, and other runtime allocations add
to this `O(sqrt(n))` storage requirement.

At `n = 2^53 - 1`, `L = 94,906,265`. Calculated sizes for the current
implementation are:

| Component | Bytes | MiB (`2^20` bytes) |
| --- | ---: | ---: |
| Dense pi element payload | 379,625,064 | 362.039 |
| Prime-vector reservation request: 5,534,701 `int64_t` elements | 44,277,608 | 42.226 |
| Sum of those two sizes | 423,902,672 | 404.265 |

The reservation comes from the current `EstimatePiPrime` call in prime-table
preparation. These are arithmetic storage estimates, not measured peak resident
memory or a bound on total allocation. Actual vector capacities, allocator
overhead, lazy Phi cache pages, tiny Phi tables, P3 row metadata, and worker
resources also contribute. Prime-generation scratch is destroyed before the
dense pi array is allocated, though the allocator may retain freed memory.
Full-range runtime and peak memory have not yet been measured.

Sequential phases and their workers share the same prepared table, so an
eight-worker call does not allocate eight copies. Concurrent top-level calls
each allocate their own table; the bounded Phi cache is shared process-wide.
This storage design is also used by the existing Legendre implementation.
The supported numeric range does not guarantee that every input is practical
for the available memory.

### Reproducing native validation and timings

```sh
prime_methods_check=$(mktemp -d)
clang++ -std=c++11 -O1 -g -pthread -Wall -Wextra \
    -fsanitize=address,undefined -fno-omit-frame-pointer -Iinst/include \
    .github/development/prime-count-methods-check.cpp \
    src/PiMeissel.cpp src/PiLehmer.cpp src/PiLegendre.cpp \
    src/PhiParallel.cpp src/P2.cpp src/P3.cpp src/Phi.cpp src/Eratosthenes.cpp \
    -o "$prime_methods_check/check"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$prime_methods_check/check"
```

Build a separate executable with `-fsanitize=thread` and use
`TSAN_OPTIONS=halt_on_error=1:symbolize=0` for concurrency checks. The fixture
uses an independent segmented Eratosthenes sieve for exhaustive inputs through
3,000, deterministic samples, and prime square/cube/fourth-power neighbors.
It also checks known counts through `1e12`, thread budgets, invalid-input
recovery, simultaneous different methods, and supplemental Legendre results.

```sh
clang++ -std=c++11 -O3 -pthread -Wall -Wextra -Wpedantic -Iinst/include \
    .github/development/prime-count-methods-benchmark.cpp \
    src/PiMeissel.cpp src/PiLehmer.cpp src/PiLegendre.cpp \
    src/PhiParallel.cpp src/P2.cpp src/P3.cpp src/Phi.cpp src/Eratosthenes.cpp \
    -o "$prime_methods_check/benchmark"
"$prime_methods_check/benchmark" lehmer 1000000000000 4 5 37607912018
```

The first argument selects `legendre`, `meissel`, or `lehmer`. Output columns
are `method,n,threads,repeat,count,total_ms`. Every count is verified, and total
time includes table construction and destruction. Repetition zero is cold;
later calls reuse the process-wide cache but reconstruct the tables. Use fresh
processes for each configuration, since running one method can warm another's
cache. The threads column is a budget, not a promise that each phase launches
that many workers. Existing P2/P3 kernel timings exclude preparation and have
different cache histories, so they cannot simply be added to predict totals.

### Driver validation and timings (November 6, 2026)

Both ASan/UBSan and TSan runs passed all 6,272 checks: 6,002 exhaustive counts,
140 sampled counts, 72 cutoff boundaries, 14 known counts, 20 thread-budget
cases, six supplemental Legendre results, eight invalid inputs, eight recovery
counts, and two concurrent results. The independent sieve streamed through 88,529,282
to cover the largest fourth-power neighbors while retaining a 64 KiB block.
No compiler warnings or sanitizer diagnostics were reported.

A clean isolated package installation compiled both new sources and passed
`testPrimeCount.R` (15 expectations) and `testErrors.R` (207 expectations),
without warnings, failures or skips. Source review checked the identities,
table coverage, signed accumulation, small inputs, and exception lifetimes.
Existing production helpers and R dispatch were unchanged, as verified against
the initial working-tree snapshot; the Git index was also unchanged.

Complete native timings used Apple clang, C++11, `-O3 -pthread`, and arm64
macOS. Three fresh processes per method/input/thread configuration each made
one cold and four warm calls. Configuration order and the three method orders
were shuffled with seed 202611061 after other builds and checks finished.
All 540 timed results matched known counts, extending driver validation through
`1e13`. Each cell below is **cold / warm milliseconds**: median first call and
median of per-process warm medians.

| Input | Threads | Legendre | Meissel | Lehmer |
| ---: | ---: | ---: | ---: | ---: |
| `1e8` | 1 | 0.242 / 0.065 | 0.250 / 0.047 | 0.271 / 0.052 |
| `1e8` | 4 | 0.235 / 0.071 | 0.292 / 0.083 | 0.310 / 0.089 |
| `1e8` | 8 | 0.258 / 0.100 | 0.317 / 0.106 | 0.337 / 0.117 |
| `1e10` | 1 | 3.510 / 1.226 | 3.462 / 1.249 | 3.585 / 1.301 |
| `1e10` | 4 | 2.596 / 0.867 | 1.479 / 0.661 | 1.461 / 0.673 |
| `1e10` | 8 | 2.434 / 0.774 | 1.442 / 0.625 | 1.502 / 0.692 |
| `1e12` | 1 | 83.097 / 49.760 | 82.711 / 49.386 | 84.757 / 51.188 |
| `1e12` | 4 | 26.818 / 16.710 | 25.902 / 15.902 | 25.728 / 16.141 |
| `1e12` | 8 | 23.817 / 15.829 | 23.538 / 14.544 | 22.002 / 14.060 |
| `1e13` | 1 | 449.170 / 388.530 | 453.270 / 392.303 | 454.806 / 394.427 |
| `1e13` | 4 | 142.698 / 121.994 | 136.561 / 117.832 | 134.650 / 116.550 |
| `1e13` | 8 | 128.483 / 108.476 | 112.804 / 95.584 | 109.401 / 94.147 |

At `1e13` with eight requested threads, Meissel and Lehmer warm medians were
about 12% and 13% below Legendre, respectively. Their single-thread medians
were close to Legendre's, and small-input/cold behavior was mixed. These first
drivers reuse Phi-backed P2 and Legendre's Phi scheduling; the measurements do
not establish a generally fastest method or justify an automatic dispatch rule.
No scheduling or kernel tuning was included in this integration.

No full package check, other-platform validation, or nonempty driver count at
the `2^53 - 1` ceiling was run. Prior helper tests cover exact-root boundaries
there, but full-range runtime and peak memory remain unmeasured. The independent
oracle and known results are computational evidence alongside the identities
above, not proof from exhaustive testing of the entire domain.

The initial source snapshot, test/build logs, benchmark executable, raw timing
JSON and shuffled timing script are retained under the directory recorded by
`/tmp/rcppalgos-prime-methods-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-prime-methods-a8o1v9kw`.

## Shared parallel Phi orchestration

`PhiParallel(x, a, nThreads, tables)` uses the same input domain and error
messages as serial `Phi`. It borrows the caller's immutable tables and returns
`phi(x, a)`; it does not choose a prime-counting method or prepare another table.
The upper summation index is the requested `a`, not `tables.primeCount()`.

For `a >= 6`, telescoping the Phi recurrence gives
`phi(x,a) = phi(x,6) - sum(i=6..a-1, phi(floor(x/p[i+1]),i))`.
Disjoint half-open intervals partition exactly these terms; worker results are
combined only after joining. This identity does not assume Legendre's choice
of `a`. `PiLegendre` still chooses its prepared square-root prime count and adds
`a - 1` after the shared evaluation. That extraction preserved root preparation
and the then-existing public R interface.

Validation precedes scheduling, including invalid `a` when `x == 0`. Tiny
cutoffs (`a <= 6`) and inputs within the prepared table use serial `Phi`, which
retains its boundary and direct pi-table shortcuts. Other inputs retain the
existing outer sum even on the serial path. Parallel work requires `x >= 1e7`
and more than one requested worker, and workers are capped at `a - 6` before
allocation. Recursive worker calls use serial `Phi`; callers choose the outer
parallel level to avoid nested worker launches.

`ThreadBatch` accepts copyable C++11 callables and owns their threads and fixed
exception slots. Bounds are captured by value; tables and distinct result slots
are captured by reference and must outlive the batch. All workers are joined
before a worker exception is rethrown in launch order. If starting a later
thread fails, the batch destructor joins earlier workers before the launch
exception propagates. Exceeding the declared capacity throws `std::length_error`.
The helper does not call R APIs or reuse workers between batches.

The interval policy remains the one tuned for Legendre, including its first,
middle, final, and small-input batches. General cutoff correctness does not
establish that these chunk sizes are optimal for Meissel or Lehmer. The new
drivers use this policy unchanged; tuning it for their cutoffs remains open.

The direct parallel-Phi checker uses independent partial sieves and
inclusion-exclusion, differential checks against serial Phi, concurrent callers,
and focused generic batch-lifetime checks:

```sh
phi_parallel_check=$(mktemp -d)
clang++ -std=c++11 -O1 -g -pthread -fsanitize=address,undefined \
    -fno-omit-frame-pointer -Iinst/include \
    .github/development/phi-parallel-check.cpp \
    src/PhiParallel.cpp src/Phi.cpp src/Eratosthenes.cpp \
    -o "$phi_parallel_check/check"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$phi_parallel_check/check"
```

Use a separate executable with `-fsanitize=thread` and
`TSAN_OPTIONS=halt_on_error=1:symbolize=0` for race detection. The scheduler and
worker-failure preparers below now instrument the extracted sources together,
without adding production test hooks.

### Extraction validation (November 5, 2026)

The before baseline is the working source immediately before this extraction,
including the approved private-evaluator check optimization and the maintainer's
staged formatting/timing comments. The Git index was not changed.

- Before/after scheduler dumps matched all 708 records from 53 existing cases.
  The extended synchronous scheduler checker passed ASan/UBSan with 136 cases:
  14 known prime counts, 54 general Phi counts, and 65,451 bounded contiguous
  intervals, including oversized worker requests and `a` below the table count.
- The new direct `PhiParallel` checker passed ASan/UBSan and TSan: 6,966 small
  independent cases, 10 near the parallel threshold, 8 exact large inputs
  (including `2^53 - 1`), 14 serial/parallel differential cases, two concurrent
  callers, 18 invalid inputs, and four generic batch checks.
- The real-worker failure harness passed both sanitizer builds: 28 scenarios,
  24 blocked-join checks, and 32 successful counts. Its cpp11/R bridge passed
  20 error/recovery pairs, verifying errors reach R only after workers join.
- The production cache/count checker passed both sanitizer builds at four
  workers, including allocation failure and recovery. Each normal run checked
  662 independent Legendre counts, 1,172 small counts, 84 signed Phi values,
  two concurrent Legendre calls, and 1,425,075 populated cache cells.
- A clean isolated package installation passed `testPrimeCount.R` and
  `testErrors.R`. The native builds explicitly used C++11 with Apple clang on
  arm64 macOS. No compiler warnings or sanitizer diagnostics were reported.

A source review also compared the extracted foreman, chunk-factor function,
and scheduler against the snapshot: interval arithmetic is unchanged after the
intended name/call substitutions and removal of the cancelling `a - a` term.
No full package check or validation on other platforms was run for this focused
internal extraction.

### Extraction timing check

The before/after native binaries used `-std=c++11 -O3 -pthread`. Five fresh
processes per version/input/worker combination each made one cold and three
warm calls. Paired version order and case order were shuffled with seed
202611054, after all other builds and checks finished. All 560 counts matched.
Times below are milliseconds: cold medians of first calls, warm medians of
per-process warm medians.

| Input | Workers | Before cold | After cold | Before warm | After warm | Warm change |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e8` | 1 | 0.228 | 0.228 | 0.065 | 0.065 | -0.38% |
| `1e8` | 4 | 0.246 | 0.241 | 0.070 | 0.073 | +4.28% |
| `1e8` | 8 | 0.261 | 0.245 | 0.101 | 0.104 | +3.03% |
| `1e10` | 1 | 3.450 | 3.503 | 1.229 | 1.222 | -0.59% |
| `1e10` | 4 | 2.561 | 2.575 | 0.869 | 0.865 | -0.55% |
| `1e10` | 8 | 2.418 | 2.457 | 0.784 | 0.761 | -2.91% |
| `1e10 + 1` | 4 | 2.006 | 2.019 | 0.965 | 0.970 | +0.52% |
| `1e10 + 1` | 8 | 1.620 | 1.706 | 0.834 | 0.863 | +3.51% |
| `1e12` | 1 | 81.881 | 81.792 | 48.313 | 48.335 | +0.05% |
| `1e12` | 4 | 26.637 | 26.590 | 16.541 | 16.526 | -0.09% |
| `1e12` | 8 | 23.620 | 23.720 | 15.038 | 15.067 | +0.20% |
| `1e13` | 1 | 445.507 | 448.454 | 383.914 | 385.695 | +0.46% |
| `1e13` | 4 | 139.254 | 138.902 | 119.287 | 118.936 | -0.29% |
| `1e13` | 8 | 125.481 | 124.974 | 108.775 | 108.161 | -0.56% |

The larger warmed cases (`1e12` and `1e13`) differed by less than 0.6% in
either direction. Small warmed cases ranged from -2.9% to +4.3%; these small
absolute differences and mixed directions do not establish a consistent
regression or improvement. This comparison measures the extraction with
Legendre inputs, not the suitability of its scheduling policy for future
Meissel/Lehmer cutoffs.

Sources, fixture hashes, build/test logs, raw timing JSON, and the timing script
are retained under the directory recorded by `/tmp/rcppalgos-phi-parallel-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-phi-parallel-fnn6bu6s`.

## P2 and P3 counting kernels

The initial [P2 profiling study](p2-profile/README.md) separates complete-driver
phase timing, isolated cold/warm P2, OS stack sampling, and diagnostic work
counters. It retains reproducible harnesses and measured results without
changing the production kernels. At `1e13`, warm P2 accounts for about 67%/87%
of serial Meissel/Lehmer time; the dominant sampled work is inside Phi.
The same note now includes exact cache-ineligible state traces and serial
memoization replays. Reuse concentrates in small descendants; a prospective
doubling of the cache's argument limit removes 9.66%/14.34% of evaluator calls
at `1e13`, with 7.20 MiB of additional occupied-page cell payload. These are
work and storage estimates, not measured speedups or an implemented change.
The subsequent doubled-cache benchmark in that note measures all
three complete drivers. Warm `1e13` calls improve about 24–28%, while first
`1e10` calls regress about 14–20%; observed large-input cache storage increases
by 7.21 MiB. The candidate passed ASan/UBSan and TSan. The approved change is
now applied to production, with `x <= 131071` and 32 pages per row; the original
comparison patch and data remain as historical records.

`P2(x, a, nThreads, tables)` and `P3(x, a, nThreads, tables)` count integers
through `x` with exactly two or three prime factors, counted with multiplicity,
all greater than `p[a]`. `a == 0` permits all primes. They borrow immutable
`PhiTables`, whose new `pi()` accessor exposes only a const reference. Both
validate `0 <= x <= 2^53 - 1` and `0 <= a <= tables.primeCount()` before any
empty-sum return. `nThreads <= 1` requests serial evaluation. The APIs are
internal C++ interfaces; the R method selector calls the complete counting
drivers rather than exposing these partial counts directly.

### Preparation and counting invariants

P2 requires a table through `floor(sqrt(x))` for `x >= 4`. Let
`b = pi(floor(sqrt(x)))`. Its sum over `a < i <= b` is
`pi(floor(x/p[i])) - (i - 1)`. Each pair is nondecreasing, `p[i] <= p[j]`, so
squares count once and swapped factors are not duplicated. If `a >= b`, the
sum is zero without accessing `p[a+1]`.

P2 answers quotients within the prepared table directly. The exact split is
`floor(x/p) <= limit` iff `p > floor(x/(limit+1))`; it separates this suffix once
instead of testing table membership in every expensive iteration. Larger
quotients use `pi(q) = Phi(q, pi(floor(sqrt(q))), tables) +
pi(floor(sqrt(q))) - 1`. Since `q >= p[i] >= 2` and `sqrt(q) <= sqrt(x)`, both
the identity and the prepared-table coverage hold. The existing serial Phi
implementation and its bounded atomic cache are reused. No additional dense
pi table, segmented sieve, or vector of quotient results is allocated.

P2 workers obtain contiguous chunks of 32 indices from a relaxed atomic
counter. Each claim owns disjoint work; the counter does not publish other
data. Thread launch makes immutable tables visible, and join makes independent
subtotals visible. Dynamic assignment accommodates the uneven cost of pi
queries. Worker selection initially requires 2,048 terms per worker, and
small jobs remain serial. These are measured implementation choices, not
mathematical limits or universally optimal settings for arbitrary cutoffs.

An initial segmented P2 prototype passed correctness checks and scaled across
threads, but the simpler Phi-backed formula was substantially faster for some
warmed workloads. At `1e10` with the fourth-root cutoff, one-worker warm medians
were 8.457 ms for streaming versus 0.951 ms for the straightforward Phi-backed
formula. At the cube-root cutoff they were 1.159 versus 0.670 ms. Cold behavior
was mixed: at `1e10`/cube root streaming took 1.200 ms versus 2.088 ms, while at
the fourth-root cutoff it took 8.480 ms versus 2.883 ms. This comparison used
three shuffled fresh processes and one cold plus two warm calls each, with
prepared-table setup excluded. All 234 timed partial counts in that initial
screen matched their expected values.

The maintainer chose the simpler Phi-backed implementation after this result.
The segmented source, executable, and full timing results are retained in the
`segmented-experiment` directory under the artifact path recorded below. P2
still evaluates its explicit pair-count sum; neither kernel delegates the
whole requested prime count to another method. Performance at the 53-bit
ceiling has not been established.

P3 requires table coverage through `floor(cbrt(x))` for `x >= 8`. If triples
exist, it additionally requires coverage through
`floor(floor(x/p[a+1])/p[a+1])`. Insufficient preparation throws
`std::out_of_range`; the kernel does not silently allocate a larger table.
At Lehmer's cutoff a square-root table suffices. A smaller arbitrary `a` may
need more coverage than `PhiTables` currently supports.

For nondecreasing triples `i <= j <= k`, the outer index ends at
`pi(floor(cbrt(x)))`; the inner index runs from `i` through
`pi(floor(sqrt(floor(x/p[i]))))`. Each pair contributes
`pi(floor(floor(x/p[i])/p[j])) - (j - 1)`. The largest required lookup also
covers every inner square-root lookup and prime index. Equal factors are
included. With multiple requested workers, cumulative inner-iteration counts
partition contiguous outer rows. At least 131,072 total inner iterations are
required before launching workers, with a 32,768-iteration grain per worker.
A cheap first-row bound skips scheduling allocation when the sum is certainly
small. The exact work sum handles other small cases. Serial evaluation avoids
this scheduling storage.

Both kernels use `ThreadBatch` for join-before-propagation and partial-launch
cleanup. Referenced tables, rows and results outlive their batches. Roots are
corrected by exact integer comparisons; successive divisions avoid forming
potentially overflowing products. Contributions are nonnegative and the total
count is at most `x`, so signed 64-bit accumulation is sufficient.

### Native checks and timings

```sh
partial_sieve_check=$(mktemp -d)
clang++ -std=c++11 -O1 -g -pthread -fsanitize=address,undefined \
    -fno-omit-frame-pointer -Iinst/include \
    .github/development/partial-sieve-check.cpp \
    src/P2.cpp src/P3.cpp src/Phi.cpp src/Eratosthenes.cpp \
    -o "$partial_sieve_check/check"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$partial_sieve_check/check"
```

Build a separate executable with `-fsanitize=thread` and use
`TSAN_OPTIONS=halt_on_error=1:symbolize=0` for concurrency checks.

```sh
clang++ -std=c++11 -O3 -pthread -Iinst/include \
    .github/development/partial-sieve-benchmark.cpp \
    src/P2.cpp src/P3.cpp src/Phi.cpp src/Eratosthenes.cpp \
    -o "$partial_sieve_check/benchmark"
"$partial_sieve_check/benchmark" verify 10000000000 4 4 1 455052511
```

The verification mode checks the Meissel/Lehmer identity against a supplied
known prime count and reports `P2,P3,Phi,a`; it evaluates the expression in the
fixture rather than calling the production drivers added afterward. Timing
modes `p2`, `p3`, and `reference`
accept the same arguments but take the expected partial count as the last
argument. The `reference` mode independently evaluates every P2 quotient using
serial Phi/Legendre with shared tables. It is a straightforward formula
baseline, not a comparison with an optimized external prime-counting library.
Output is `method,x,root_power,threads,repetition,count,setup_ms,kernel_ms`;
prepared-table setup is measured separately from the kernel. Use fresh
processes, shuffle order, and finish other builds/checks before timing.

### Final validation and performance (November 6, 2026)

The final Phi-backed P2 and tuned P3 passed both ASan/UBSan and TSan:

- 100,548 independent small counts from a smallest-prime-factor/Omega sieve.
- 72 repeated-factor boundaries, 18 larger/quotient/worker-boundary P2 cases,
  seven independent P3 checks (including explicit triple enumeration), and
  two concurrent calls.
- 196 integer-root checks against independent binary-search roots, including
  large square/cube neighbors and the 53-bit boundary; 22 preparation/empty
  cases, including valid maximal-x empty P3; 44 invalid-input/coverage cases.
- Twelve known-prime-count identities: Meissel through `1e12` and Lehmer
  through `1e10`. An additional optimized driver verified both cutoffs through
  `1e12` and the fourth-root identity at `1e13`.

An initial clean isolated package build and a rebuild of the final P2/P3
sources both passed `testPrimeCount.R` and `testErrors.R`. Native builds used
Apple clang, C++11, and `-Wall -Wextra -Wpedantic`, without warnings or sanitizer
diagnostics. Review checked the counting identities, table coverage, exact-once
chunk/row ownership, signed arithmetic, and failure-path object lifetimes.
No full package check, other-platform validation, or full-range P2 performance
run was performed.

The final timing matrix used `-std=c++11 -O3 -pthread`, three fresh processes
per case, one cold and two warm calls, shuffled with seed 202611055 after all
builds and tests finished. All 234 timed counts matched. The following are
kernel milliseconds, excluding shared-table preparation. Cold values are
medians of first calls; warm values are medians of per-process warm medians.
Threads are requested budgets: small jobs can use fewer workers.

| Kernel | Input | Cutoff root | Threads | Cold | Warm |
| --- | ---: | ---: | ---: | ---: | ---: |
| P2 | `1e8` | 3 | 1 | 0.111 | 0.018 |
| P2 | `1e8` | 3 | 4 | 0.115 | 0.018 |
| P2 | `1e10` | 3 | 1 | 2.085 | 0.641 |
| P2 | `1e10` | 3 | 4 | 0.789 | 0.233 |
| P2 | `1e10` | 3 | 8 | 0.728 | 0.235 |
| P2 | `1e10` | 4 | 1 | 2.750 | 0.902 |
| P2 | `1e10` | 4 | 4 | 0.941 | 0.308 |
| P2 | `1e10` | 4 | 8 | 0.960 | 0.308 |
| P2 | `1e12` | 3 | 1 | 55.749 | 31.438 |
| P2 | `1e12` | 3 | 4 | 15.915 | 8.585 |
| P2 | `1e12` | 3 | 8 | 13.264 | 7.148 |
| P2 | `1e12` | 4 | 1 | 69.963 | 40.868 |
| P2 | `1e12` | 4 | 4 | 19.821 | 11.416 |
| P2 | `1e12` | 4 | 8 | 16.089 | 9.245 |
| P3 | `1e10` | 4 | 1 | 0.032 | 0.031 |
| P3 | `1e10` | 4 | 4 | 0.033 | 0.031 |
| P3 | `1e10` | 4 | 8 | 0.033 | 0.031 |
| P3 | `1e12` | 4 | 1 | 0.673 | 0.585 |
| P3 | `1e12` | 4 | 4 | 0.277 | 0.260 |
| P3 | `1e12` | 4 | 8 | 0.249 | 0.234 |
| P3 | `1e13` | 4 | 1 | 2.858 | 2.643 |
| P3 | `1e13` | 4 | 4 | 1.236 | 1.081 |
| P3 | `1e13` | 4 | 8 | 0.923 | 0.806 |

Table preparation was approximately 0.034 ms at `1e8`, 0.26 ms at `1e10`,
2.5 ms at `1e12`, and 7.7 ms at `1e13`; the new drivers share that setup
with Phi and the applicable kernels. At `1e12`/fourth-root cutoff, warmed P2 improved
from 40.868 ms with one worker to 9.245 ms with eight (about 4.4x). P3 at that
input improved from 0.585 to 0.234 ms; at `1e13` it improved from 2.643 to
0.806 ms. P2 at `1e10` caps requested eight workers at four. P3 at `1e10`
stays serial, avoiding the launch overhead seen in the initial version.

The simpler P2 is not uniformly faster than the rejected segmented prototype.
At `1e12`/cube-root cutoff, the prototype had one-worker cold/warm medians of
29.429/29.378 ms versus 55.749/31.438 ms for the final implementation; with
eight workers they were 6.034/5.888 versus 13.265/7.148 ms. Conversely, at the
fourth-root cutoff, prototype one-worker cold/warm times were 363.199/358.478
ms versus 69.963/40.868 ms, and eight-worker times were 67.991/71.293 versus
16.089/9.245 ms. These are separate shuffled runs on the same machine, not a
claim of a universal crossover. The selected design favors simplicity and
the shared Phi infrastructure; a hybrid was explicitly deferred.

A focused P3 tuning screen also rejected doubling its per-worker grain: that
slowed a `5e10` four-worker midpoint from about 0.059 to 0.071 ms. The retained
change only adds the minimum total-work cutoff, preserving useful parallelism
on larger workloads. Its timing artifacts are under `p3-small-work`.

Sources, old and final timing JSON, scripts, and validation logs are retained
under the directory recorded by `/tmp/rcppalgos-partial-sieve-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-partial-sieve-k6jpym97`.

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

### November 5 performance investigation (diagnostic only)

No production changes were applied. Temporary variants used the current
working `Phi.cpp` (including explanatory comments), `PiLegendre.cpp` from
`034712d`, and the retained pre-integration Legendre snapshot. All native
builds used Apple clang 17, arm64 macOS, and `-std=c++11 -O3 -pthread`.
Clock readings inside the counter separated table preparation from the
summation; the external driver also measured total call time, including table
destruction. Instrumentation was identical across compared implementations.

The first three screens varied translation-unit visibility, evaluator reuse
per foreman (both signs), tiny-table reference reuse, constant six-prime
arithmetic, exact-root correction by products, prepared-table representation,
and redundant private checks. Each screen shuffled fresh processes, with
three processes per input/worker/variant and one cold plus four warm calls.
Combining source files, reusing the evaluator, and tiny-table changes did not
recover the main regression. Product-based root correction also recovered
little. These experiments do not prove those costs are zero or predict other
compilers' results.

Separate serial work counters at `1e12` showed identical warmed counts of
1,382,174 recursive expansions, 9,514,107 cache hits, 33,501,450 direct suffix
terms, and 2,764,350 tiny-table evaluations before/after integration. The
public Phi fast path accounts for one fewer private evaluator entry. These
counts support a higher cost per operation, rather than extra recursive work,
for this input; instrumented counter builds were not used for timings.

The most promising diagnostic combination (`lean_raw` in the artifacts):

- Captured prepared-table data pointers and sizes once per evaluator, while
  retaining the same immutable owning tables and lifetimes.
- Retained public zero/range validation, but omitted the private evaluator's
  zero-x/zero-a checks: public callers handle these cases, and recursive calls
  use positive quotients and positive sieve indices.
- Omitted only the recursive loop's redundant lookahead-bound test. Since
  `i < b <= a`, the prepared prime at `i + 1 <= a` exists. The general entry
  shortcut retained its missing-lookahead guard.

A second, independently shuffled confirmation used five fresh processes per
input/worker/version, one cold and four warm calls per process. The table shows
medians of the five per-process warm medians, in milliseconds. No other builds
or tests ran concurrently with these benchmarks.

| Input | Workers | Before integration | Current | Diagnostic combination | Reduction vs current |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `1e10` | 1 | 1.071 | 1.372 | 1.285 | 6.3% |
| `1e10` | 4 | 0.753 | 0.939 | 0.904 | 3.7% |
| `1e10` | 8 | 0.657 | 0.832 | 0.807 | 3.0% |
| `1e12` | 1 | 43.089 | 55.762 | 51.701 | 7.3% |
| `1e12` | 4 | 15.506 | 18.734 | 17.406 | 7.1% |
| `1e12` | 8 | 14.112 | 17.043 | 16.072 | 5.7% |
| `1e13` | 1 | 334.463 | 448.815 | 413.201 | 7.9% |
| `1e13` | 4 | 107.485 | 138.361 | 128.175 | 7.4% |
| `1e13` | 8 | 97.076 | 122.289 | 113.417 | 7.3% |

At `1e12`/one worker, median preparation times were about 1.92 ms before
integration and 2.56 ms currently, versus total warm times of 43.09 and
55.76 ms. Preparation explains only a small part of this serial gap. The
combination recovers about 4.06 ms, but still trails the old implementation
by about 8.61 ms. The remaining dominant machine-code cause has not been
isolated; do not attribute the entire regression to validation or object
construction. Cold-call timings exclude process startup, including the old
implementation's eager tiny-table initialization.

All 2,475 timed counts across the screens and confirmation matched known
answers. The diagnostic combination also passed the standalone general-Phi
checker and native cache checker under ASan/UBSan: exhaustive/general-domain
and invalid-input cases, independent counts, signed phi values, concurrent
calls, and cache-cell validation. These checks support further investigation;
they are not approval to apply the prototype. No R package build, TSan run,
full package check, or other-platform benchmark was performed for these
unapplied variants.

Reproduction scripts, exact variant sources, compiler commands, source hashes,
raw timings, work counts, and sanitizer logs are retained in the temporary
artifact directory recorded by `/tmp/rcppalgos-phi-performance-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-phi-performance-kbl3iz9q`.
The benchmark seeds are 20261105 for the screens and 11052026 for confirmation.

### Follow-up: private evaluator invariants (diagnostic only)

The follow-up found a smaller candidate that retains the existing vector
references, cache, table ownership, exact-root correction, public API, and
scheduling. These variants were temporary during the investigation; the
subsequent approved implementation is recorded below.

The private evaluator has stronger preconditions than the public Phi API:

- Public entry already handles zero inputs, `x <= p[a]`, and `a <= 6`.
- A recursive loop exists only if `b > c = min(b, 6)`. Therefore `c = 6`
  and its index `i >= 6`.
- For `j = i + 1 <= b`, `p[j] <= floor(sqrt(x))`, including when the
  prepared table ends below that root. Thus the child quotient satisfies
  `floor(x / p[j]) >= p[j] > p[i]`.

By induction, every private call has `a >= 6` and `x > p[a]`. Its zero-input
and `x <= p[a]` branches are unreachable; its tiny case can explicitly use
`a == 6` and the fixed six-prime formula. The general expansion must still
retain `c = min(b, 6)`: `Phi(18, 7, PhiTables(17))` needs `b = c = 2`.

A further duplicate check can be moved out of private entry. The parent tests
`isPi(quotient, i)` before recursing, and a successful test breaks into the
direct suffix. Therefore every recursive child has already failed that test.
Moving the initial direct-table shortcut into public Phi permits omitting it
at private entry. Public entry retains its table and lookahead bounds;
inside the loop, `i + 1 <= b <= a` guarantees the lookahead prime exists.
Immutable tables make these conditions stable across calls and threads.

The resulting `checks_hoist` candidate makes only these invariant-based
changes; `lean_hoist` also retains the preceding raw-pointer/size experiment.
Additional screens tested a directly addressed global cache, outlining cold
allocation code, const pointer members, separated expansion helpers, and
controlled recursive inlining. These did not supply a compelling additional
improvement. Compiler-specific noinline attributes were diagnostic only.
Static assembly counts included cold initialization/unwinding code and were
not treated as dynamic instruction costs.

Confirmation used the same Apple clang 17/arm64/C++11/-O3 setup, five fresh
processes per input/worker/variant, one cold plus four warm calls per process,
and shuffled order (seed 202611052). No builds or validation ran concurrently.
Warm numbers below are medians of the five per-process medians, in ms.
"Earlier candidate" is `lean_raw` from the preceding investigation.

| Input | Workers | Before integration | Current | Earlier candidate | Checks only | Checks + pointers | Reduction, checks only |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e8` | 1 | 0.050 | 0.069 | 0.065 | 0.064 | 0.066 | 7.5% |
| `1e10` | 1 | 1.056 | 1.340 | 1.286 | 1.221 | 1.245 | 8.9% |
| `1e12` | 1 | 42.967 | 55.828 | 51.796 | 48.412 | 48.150 | 13.3% |
| `1e12` | 4 | 15.923 | 19.115 | 17.562 | 16.635 | 16.532 | 13.0% |
| `1e12` | 8 | 14.591 | 17.030 | 16.288 | 15.477 | 15.138 | 9.1% |
| `1e13` | 1 | 336.469 | 449.383 | 414.700 | 384.973 | 384.240 | 14.3% |
| `1e13` | 4 | 111.081 | 139.122 | 130.993 | 124.666 | 119.374 | 10.4% |
| `1e13` | 8 | 98.056 | 123.767 | 116.433 | 110.185 | 108.271 | 11.0% |

The checks-only candidate recovers about 57% of the serial integration gap
at both `1e12` and `1e13`, while still trailing the old code by about 13-14%.
Raw pointers add only about 0.5% and 0.2% over that candidate in these serial
cases. Parallel medians sometimes favor pointers, but process-median ranges
overlap substantially: at `1e13`/four workers, 117.75-125.90 ms for vector
references versus 118.30-124.48 ms for pointers. Pointer conversion is therefore
a separate choice, not a prerequisite for the stronger optimization. Exact
root correction by products also added little and is not part of the smaller
candidate. The remaining regression has not been fully explained.

All 2,055 timed counts across two screens and confirmation matched known
answers. Both `checks_hoist` and `lean_hoist` passed the general-Phi and native
cache checkers under ASan/UBSan and TSan, including cached/uncached modes,
allocation failure, and concurrent calls. Each candidate also checked the
explicit `Phi(18,7,PhiTables(17)) = 1` boundary by independent enumeration.
There were no compiler warnings or sanitizer findings.

Validation linked the unchanged production PiLegendre launcher: its diagnostic
copy has global phase timers appropriate only to the timing driver's sequential
calls, not concurrent-call tests. The private signed checker still includes
some direct-shortcut inputs (e.g. `(169,7)`) outside the narrowed private-entry
contract; their existing expectations were retained and passed. Public Phi
continues to support these inputs normally. The invariant proof, rather than
these private tests alone, justifies the narrower recursive contract.

This investigation supported the subsequent checks-only implementation. The
diagnostic variants did not restore pre-integration performance completely and
were not tested through an R package build or on other platforms. Artifacts (sources, compiler
commands, scripts, timing JSON, proof, sanitizer logs) are under the directory
recorded by `/tmp/rcppalgos-phi-deeper-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-phi-deeper-1jrq1mc8`.

### Implemented: private evaluator checks

The approved checks-only change was applied on November 5, 2026. Public Phi
retains input validation and all general base cases, and now handles the initial
direct pi-table shortcut. Private recursive entry relies on the documented
`a >= 6`, `x > p[a]` invariants and the parent's failed shortcut check. Its tiny
case explicitly uses six primes; the expansion retains `c = min(b, 6)`.
The unused private shortcut helper was removed, and the loop checks table
membership and the prime square directly, using its proven lookahead bound.
Vector references, cache synchronization, table ownership, root correction,
public signatures, and Legendre scheduling are unchanged.

Regression coverage adds the short-table cases at 17/18/19 and 24/25/26 with
`a = 7`, and both signs of the private six-prime terminal case at 289 and the
30030-period neighbors. All existing cases and expectations remain. The general
checker now has 171 targeted cases; the cache checker has 84 signed phi checks.

Validation of the actual production source passed without compiler warnings
or sanitizer findings:

- General Phi, cached and uncached, under ASan/UBSan and TSan: 2,621,310
  exhaustive cases, 171 targeted cases, 220 large-input cases, 3,328 concurrent
  evaluations, 57 invalid inputs, and known counts through `1e12` (through
  `1e8` without caching).
- Native cache checks at four workers under both sanitizer builds, including
  allocation failure/recovery: 662 independent Legendre counts, 1,172 small
  counts, 84 signed phi checks, two concurrent Legendre calls, and 1,425,075
  populated cache cells.
- The worker-failure harness under ASan/UBSan: 28 scenarios, 24 blocked-join
  checks, and 32 successful recovery counts.
- A clean temporary R package installation, followed by `testPrimeCount.R`
  and `testErrors.R` against that installation.

A final native timing comparison compiled the before/after production sources
with Apple clang 17 and `-std=c++11 -O3 -pthread`, using the unchanged timing
driver and scheduler. No internal phase instrumentation was used. Three fresh
processes per version/input/worker combination each made one cold and three
warm calls. Process order was shuffled with seed 202611053 after all builds
and tests finished. All 192 results matched known answers. Cold values are
medians of first calls; warm values are medians of per-process warm medians,
in milliseconds.

| Input | Workers | Before cold | After cold | Before warm | After warm | Warm reduction |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e8` | 1 | 0.246 | 0.233 | 0.067 | 0.065 | 2.9% |
| `1e10` | 1 | 3.710 | 3.422 | 1.339 | 1.225 | 8.5% |
| `1e12` | 1 | 90.594 | 81.644 | 54.588 | 48.299 | 11.5% |
| `1e12` | 4 | 29.125 | 26.656 | 18.609 | 16.567 | 11.0% |
| `1e12` | 8 | 25.886 | 23.673 | 16.958 | 15.102 | 10.9% |
| `1e13` | 1 | 507.436 | 444.734 | 445.575 | 384.194 | 13.8% |
| `1e13` | 4 | 157.022 | 139.370 | 137.705 | 118.235 | 14.1% |
| `1e13` | 8 | 139.517 | 119.535 | 121.755 | 106.496 | 12.5% |

The larger warmed cases improved about 11-14% in this run. This does not remove
the entire regression against the pre-integration implementation. No full
package check or other-platform validation was run for this focused change.
The existing worker harness was rerun under ASan/UBSan; its earlier TSan results
are not claimed as a new run here. Native checks explicitly used C++11.

Before-source snapshots, build/test logs, benchmark script, and raw timing JSON
are retained under the directory recorded by `/tmp/rcppalgos-phi-checks-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-phi-checks-phr1r05m`.

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
    src/PiLegendre.cpp src/PhiParallel.cpp src/PiPrime.cpp src/Eratosthenes.cpp \
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
    src/PiLegendre.cpp src/PhiParallel.cpp src/PiPrime.cpp src/Eratosthenes.cpp \
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

The shared cache uses `a < 100`, `x <= 131071` and retains reuse across worker
batches and calls. Each row owns up to thirty-two fixed blocks
of 4,096 atomic cells, allocated on demand. Cells and pointers are explicitly
initialized for C++11. Published blocks remain alive and stationary until cache
destruction, which requires all workers to have joined.

The larger argument limit retains 16-bit answer cells. Only positive arguments
with `a >= 7` reach the cache, so monotonicity gives
`1 <= Phi(x,a) <= Phi(131071,6) = 25142 < 65536`. The page count derives from
the explicit argument bound, independently of the answer type. The permanent
cache checker covers both the former boundary and the new endpoint/page edges.

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
    src/PiLegendre.cpp src/PhiParallel.cpp src/Phi.cpp src/PiPrime.cpp src/Eratosthenes.cpp \
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
