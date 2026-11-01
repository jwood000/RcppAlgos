# Prime sieve cache experiments

`RcppAlgos:::primeSieveInspect(bound1, bound2 = NULL, nThreads = NULL,
segmentBytes = NULL)` is an internal developer helper. It returns the primes
along with phase timings and estimates of storage used by the sieve. It is
available only in the isolated developer build created by the
[inspector toolkit](prime-sieve-inspector/README.md). The production package
contains no inspector hooks, wrapper, registration, or extra kernel variants.
Ordinary `primeSieve()` retains its interface and default sizing policy.

The inspector measures elapsed time and records data structure sizes. It does
**not** measure cache hits, misses, cache residency, memory bandwidth, or hardware
cache capacity, and it cannot place a vector in a particular cache. The storage
estimates help form hypotheses about cache pressure; timings help test them.

See the [measurement suite and findings](prime-sieve-benchmarks/REPORT.md) for
the recorded segment-size sweep, parallel comparisons, and recommendations.

In the experimental build, the guarded normal/profiled kernel variants increase
generated code. In a local
Apple Clang 17 C++11 `-O3` build, `Eratosthenes.o`'s text section grew from
85,633 to 142,797 bytes versus commit `0364f49` (about 56 KiB). A single timed
compilation took 0.75 versus 1.16 seconds; those compile times are indicative,
not a repeated benchmark. This is the cost of retaining one algorithm source
while keeping phase-clock checks out of ordinary marking/extraction loops.

## First inspection

```r
library(RcppAlgos, lib.loc = "/tmp/rcppalgos-inspector/library")
inspect <- RcppAlgos:::primeSieveInspect
lower <- 4e15
upper <- lower + 1e6

result <- inspect(lower, upper, nThreads = 1)
result$summary
result$phases
result$kernels
result$shared
result$host
stopifnot(identical(result$primes,
                    RcppAlgos::primeSieve(lower, upper, nThreads = 1)))
```

Use the same bounds and output requirements when comparing configurations.
Returning the complete prime vector deliberately includes normal result costs;
the inspector is not a prime-count benchmark.

## Segment-size override

`segmentBytes` specifies a **per-worker sieve payload budget in bytes**, not a
hardware cache size or a total memory limit. The Big sieve stores one logical
bit per integer position. The budget is converted to a multiple of the
30030-position wheel, rounding down:

```text
segment span = floor(8 * segmentBytes / 30030) * 30030
```

The minimum accepted budget is 3754 bytes, enough for one wheel. The argument
must be a finite, positive, whole numeric scalar within the implementation's
checked integer range. Actual `vector<bool>` capacity can include allocation
rounding; the reported segment span is the authoritative algorithmic size.

An override is accepted only when the entire interval would naturally use Big
under the existing dispatch policy. It does not force a kernel, and an
inapplicable override reports an error mentioning Big. The cutoff depends on
host core counts and requested threads, so an example's bounds may need to be
increased on another host. Keep bounds within the package's supported domain.
Small and Medium are still inspectable with `segmentBytes = NULL`.

```r
budgets <- c(64 * 1024, 256 * 1024, 1e6)
results <- lapply(budgets, function(bytes) {
    inspect(lower, upper, nThreads = 1, segmentBytes = bytes)
})
stopifnot(all(vapply(results, function(x) {
    identical(x$primes, result$primes)
}, logical(1))))
lapply(results, function(x) x$summary)
```

Changing the span also changes the dense/bucket prime split and the number of
segments available to workers. Consequently, this experiment varies more than
the flag array's size. Start with one requested worker; use the recorded kernel
and worker counts when later comparing parallel runs.

## Reading the results

- `primes`: the normal ordered result, including its integer or double type.
- `summary`: requested and effective worker counts, requested byte budget,
  effective override span, and elapsed seconds for base-prime preparation,
  policy/setup work, sieve execution, and conversion to the R result. An
  effective override span of zero means automatic sizing; consult `kernels`
  for the actual span of each kernel invocation.
- `phases`: elapsed seconds and call counts by worker, kernel, and phase.
  Worker identifiers start at one. Calls within one phase are accumulated.
- `kernels`: invocation bounds, segment span, segments visited, primes emitted,
  and storage estimates for flags, offsets, buckets, extraction buffer, and
  output capacity. One worker can invoke several kernels when bounds cross
  dispatch cutoffs.
- `shared`: storage estimates for base primes, the dense/bucket prime arrays,
  and the shared wheel lookup.
- `host`: operating system, architecture, R version, and the core/thread counts
  available to the package. This is context, not a cache-topology report.

The phase labels are `setup` (allocation and local preparation), `offsets`
(dense-prime starting positions), `buckets` (scheduling sparse strikes),
`dense_mark`, `bucket_mark`, `extract`, and `reset`. Not every kernel uses
every phase. The Small lookup-table path has no segmented flag storage.
`bucket_entries` counts scheduled strike offsets, not distinct sieving primes;
`bucket_capacity_bytes` sums the bucket vectors' allocated payload capacities
after preparation.

Times use a monotonic clock. Instrumentation itself adds work, especially for
small segments and short phases. The sieve summary is elapsed wall time across
execution; adding worker times does not reproduce wall time when workers run
concurrently. Neither the C++ phase timings nor their sum is the full R call
time: R dispatch, diagnostic-result assembly, and host metadata also cost time.

Within the experimental build, normal and inspected calls compile separate
variants of the same kernel source.
Phase-clock checks are compiled out of the normal loops; a final timer-destructor
check remains. This keeps timing work out of production loops at the cost of
additional compiled code, while maintaining one source implementation of each
algorithm.

Storage fields estimate payload or capacity, not process RSS or peak live
memory. Allocator metadata, container bookkeeping, thread stacks, cache-line
rounding, and the R result are not fully represented. Summing per-invocation
capacities can overcount reused output storage and does not establish how much
memory was simultaneously resident. Bucket storage and counts should be read
as implementation measurements, not cache-miss counters.

## A reproducible comparison

1. Fix a workload grid: starting magnitude, interval width, and alignment.
   Include short interior intervals and intervals crossing several segments.
2. Record the revision, compiler/build flags, host context, and whether this
   is a release or diagnostic build. Keep unrelated machine activity low.
3. Warm up each configuration, then repeat in randomized or alternating order.
   Save every observation; compare medians and variability rather than a
   single fastest run.
4. Verify complete prime outputs separately from the timed region. Use an
   independent primality check on manageable intervals as well as comparison
   with the unchanged production entry point.
5. Inspect preparation, marking, extraction, and total costs together. Confirm
   any promising result with ordinary `primeSieve()` or a minimally instrumented
   experiment before changing the production sizing policy.

Do not automatically interpret a timing discontinuity as an L1/L2/L3 boundary.
The number of segments, bucket scheduling, allocation, and worker eligibility
can also change at the same point. Hardware performance counters would require
a separate platform-specific profiling experiment.

## Developer validation

These are historical checks of the instrumented implementation. For the current
build-and-test command, use the [toolkit README](prime-sieve-inspector/README.md).
The production source files have since been restored to the uninstrumented base.

Normal calls were compared with commit `0364f49` using a standalone Apple
Clang 17 ARM64 C++11 `-O3` build. After warming each version, 15 runs alternated
execution order; these are median native milliseconds, excluding R conversion:

| Interval | Committed baseline | Inspector disabled |
|---|---:|---:|
| `1 ... 1e7` | 2.034 | 2.037 |
| `1 ... 1e8` | 22.808 | 22.647 |
| `1e12 ... 1e12 + 1e6` | 1.112 | 1.124 |
| `1e15 ... 1e15 + 1e6` | 13.590 | 13.661 |

These native checks used a fixed policy of one physical core and ten available
threads, with one requested worker. That makes the final workload Big; it is
not the normal ten-physical-core policy used by the later measurement suite.
These local measurements check for instrumentation regressions, not a portable
performance guarantee. Twenty-eight additional native comparisons matched the
committed implementation under ASan/UBSan, including segment overrides and
parallel execution.

The native registration was regenerated with `cpp11::cpp_register(".")`.
This updates `src/cpp11.cpp` and also creates `R/cpp11.R`. This repository uses
hand-written R wrappers and does not track that generated R file. During this
change, `R/cpp11.R` was absent before generation, and the newly generated file
was removed afterward. Preserve any pre-existing file or local edits when
repeating the workflow; review generated changes before retaining them.

The installation command verified for this change was:

```sh
MAKEFLAGS=-j4 R CMD INSTALL --preclean --no-multiarch \
    --library=/tmp/rcppalgos-inspector-library /tmp/rcppalgos-inspector-package
```

Here the source directory was a fresh copy of tracked package sources plus the
new helper files, with no cached object files, and the development library
directory already existed. Use your own source and library paths when repeating
the command.

After installation, load that specific library in a fresh R process and run
both focused test files:

```r
library(RcppAlgos, lib.loc = "/path/to/development-library")
testthat::test_file("/path/to/build/package/tests/testthat/testCacheInspector.R",
                    reporter = "summary", stop_on_failure = TRUE)
testthat::test_file("/path/to/build/package/tests/testthat/testPrimeSieve.R",
                    reporter = "summary", stop_on_failure = TRUE)
```

These commands were verified against the modified package. The inspector tests
cover default and overridden sizes, independent prime outputs, kernel changes,
parallel reporting, serial fallback, and invalid or unsafe byte budgets.
