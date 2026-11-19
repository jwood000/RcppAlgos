# Deferred Investigations

## Unsafe floating-point to integer conversions

Status: Resolved; reviewed by the maintainer and merged into `main`
Priority: High
Discovered during: PartitionType refactoring

### Observation

`SetTolerance()` and `CheckPartition()` cast unbounded numeric source values
to `std::int64_t` while testing integrality. `IsDecimal()` does the same for
scalar numeric input before normal validation. `CheckIsInteger()` uses the
same pattern, although its current callers bound the values before that cast.

### Confirmed failures, November 18, 2026

An isolated build of `3da5ad74` with
`-fsanitize=undefined,float-cast-overflow -fno-sanitize-recover=all` confirmed
nine invalid-conversion probes and one signed-overflow probe among 14 small
public-call checks. These included infinite/large source-vector elements,
non-finite scalar inputs, and a mean target scaled beyond `int64_t` at width
1025 (the neighboring width 1024 passed).

The investigation also confirmed two related partition-mapping failures:

- Individually representable values can have an unrepresentable spacing or
  overflowing intermediate sums. For example, repeated pairs from
  `c(2^62, 2^62 + 1024)` overflowed while preparing a search with target 1.
- The exact spacing of `c(-(2^53 - 1), 2)` is `2^53 + 1`. Computing it in
  double rounded it and caused a valid width-two partition with target
  `-(2^53 - 3)` to fail with "Strange mapping!!!".

### Implemented fix

The maintainer approved addressing both conversion and mapping arithmetic.
Whole-number tests now use finite/trunc checks; actual conversions require
the half-open range `[-2^63, 2^63)`. Non-finite numeric scalars reach existing
validation, while large numeric vectors retain general constraint handling.

`CheckPartition()` computes integer spacing exactly. Partition mapping keeps
its original `int64_t` search. An upfront check requires
`m * max(abs(v)) + abs(target) <= INT64_MAX`, evaluated by division without
signed absolute-value overflow. The initial width bounds later searches if
the width decreases. Cases outside this conservative bound produce a clear
range error; some mathematically valid cases are deliberately outside the
supported native arithmetic range.

Mapped-target storage remains limited to `int`. The shift calculation checks
its native multiplication and subtraction before evaluating them. No GMP
source vector, templated search, or wide-intermediate fallback is added.
Existing GMP counting/ranking functionality is unchanged. The target's
existing `abs(target) < 2^53` validation also remains unchanged; this is a
native arithmetic limit, not a new source-vector magnitude restriction.

The maintainer also approved retiring the obsolete large-product chunking
fallback in `PollardRhoUtils.cpp` and its `my63Max`/`Sqrt63Max` constants.
Native primality multiplication now uses exact integer arithmetic: residues
below a modulus at most `INT_MAX` have products below `2^62`. `IsPrime()`
enforces that upper bound with a native exception, and its installed header
documents the existing small-factor trial-division precondition. R-facing
native/GMP dispatch remains unchanged.

Validation uses macOS arm64, Apple Clang 17 / R 4.6.1. Earlier revisions
passed the full UBSan package suite (43,591 expectations) and independent
mapping/primality comparisons; those results apply to the recorded earlier
snapshots, not the final guard-only mapping revision.

Final guard-only validation:

- A fresh UBSan build passed 2,064 R expectations across 16 relevant files,
  including all 88 numeric regressions, and all 20 native assertions. There
  were no failures, warnings, or sanitizer diagnostics; the existing
  unsupported-multiset-composition test was skipped. Existing GMP counting,
  ranking, sampling, iterator, and two-thread cases remain covered.
- C++11 ASan/UBSan checks passed 4,862 guard-boundary cases against exact GMP
  reference calculations, including 870 diagnostic reconstructions. GMP is
  used as a test oracle here, not in partition preparation.
- The two native search helpers and all of `PartitionsDesign.cpp` match
  baseline `3da5ad74` exactly. Source/test hashes for that guard-only revision
  matched its tested package snapshot; later zero-only validation is recorded
  below.

Permanent public regressions are in `tests/testthat/testNumericConstraints.R`.
The package sanitizer recipe is in
[testing.md](testing.md#package-floating-conversion-and-mapping-checks).
Baseline failures and earlier validation logs are retained under
`/tmp/rcppalgos-float-probe-w90bjvot`,
`/tmp/rcppalgos-float-fixed-_ss_s7et`, and
`/tmp/rcppalgos-mapping-dispatch-1rp58mtj`.
The bounded-input overflow example from the domain review is recorded under
`/tmp/rcppalgos-supported-range-oupxzbda`.
Final build/test logs, source hashes, and an index of native guard checks are
under `/tmp/rcppalgos-native-mapping-y6co4ddt`. Primality-retirement checks are
retained separately under `/tmp/rcppalgos-native-primality-l49ftkjg`.
No cross-platform check or performance benchmark was run for this change.

### Decision

The maintainer reviewed the fix and validation evidence. The fix was merged
into `main` in [PR #90](https://github.com/jwood000/RcppAlgos/pull/90), addressing
[Issue #89](https://github.com/jwood000/RcppAlgos/issues/89).

## Zero-only repeated partitions can divide by zero during setup

Status: Resolved; reviewed by the maintainer and merged into `main`
Discovered during: Review of partition-mapping arithmetic, November 18, 2026

`partitionsGeneral(0, 2, repetition = TRUE)` reaches integer division by zero.
`StandardDesign()` in `src/PartitionsUtils.cpp` reduces the requested width
to the zero target before calculating the mapping shift. The original shift
expression and the guarded `MappingShift()` both require a positive divisor.

Isolated UBSan processes reproduced this on baseline `3da5ad74` and the prior
floating-conversion fix, with verified loaded-library paths. Both stopped on
division by zero; this predates the upfront native/GMP dispatch change.
Reproduction scripts and logs are under
`/tmp/rcppalgos-zero-width-probe-yzvvc534`.

The maintainer approved a focused fix. Repeated zero-only inputs with target
zero now use the existing positive mapping path instead of `StandardDesign()`.
Mapping each zero to one makes the mapped target equal to the positive width:
there is exactly one mapped result, consisting entirely of ones. Mapping back
gives one all-zero result at the requested width, consistent with the existing
finite-frequency zero-only cases. No search or enumeration algorithm changes.
Other domains retain their existing routing.

A fresh UBSan build passed 1,417 R expectations across seven related files and
all 20 native assertions, with no failures, warnings, skips, or sanitizer
diagnostics. The 59 new expectations in `tests/testthat/testZeroPartitions.R`
cover partitions, compositions with both weak settings, counting, sampling,
ranking, iterator exhaustion/random access, integer/numeric results, and
constrained combinations/permutations. The original target-omitted reproducer
also returned one row of two zeros under UBSan. Build/test logs and snapshot
hashes are under `/tmp/rcppalgos-zero-fix-1b0mpvpj`.

Validation used macOS arm64, Apple Clang 17 / R 4.6.1. Full `R CMD check` and
cross-platform checks were not rerun for this routing-only change. The local
deferred records were checked for a duplicate; GitHub was not consulted and no
issue was published. The fix was reviewed by the maintainer and merged into
`main` in [PR #90](https://github.com/jwood000/RcppAlgos/pull/90).

## Prime sieve worker and thread-launch exception handling

Status: Resolved; reviewed by the maintainer and merged into `main`
Discovered during: Considering native prime-sieve reuse in nthPrime

### Observation

Before the fix, `PrimeSieveMain()` in `src/Eratosthenes.cpp` launched workers
directly with `std::thread` and joined them after the launch loop. Unlike the
prime-counting paths, it did not use `ThreadBatch` or catch worker exceptions. Worker sieve
and result-vector allocations can throw. An uncaught worker exception would
terminate the process; a later launch failure could also unwind through
earlier joinable threads.

### Focused review, November 12, 2026

The November 1 sieve review covered boundary correctness, serial/parallel
equivalence, sanitizers, and benchmarks. The recorded 28 worker/launch-failure
scenarios were added for Legendre in commit `9832e53`; that commit did not
change the sieve launcher.

Eight isolated native cases were run against an instrumented temporary copy
of `Eratosthenes.cpp` (SHA256
`98904804fb7c2f79d31cb1161d4313e6c2c7e7f99089730694bfb019c67275c7`):

- Serial and two-thread controls both matched an independent sieve through
  `1e6` (78,498 primes).
- Injected `std::bad_alloc` and `std::runtime_error` at worker entry were
  caught on the calling thread in serial mode; both reached `std::terminate`
  in parallel mode.
- An injected `std::system_error` before the first launch was caught normally.
  The same simulated failure before the second launch, after a barrier proved
  the first worker was active, reached `std::terminate` during unwinding.

The terminate handler exited the child process without a core dump. No actual
memory/thread-resource exhaustion was attempted. These tests establish the
exception-control-flow behavior, not the probability of resource failures.
That initial review did not test a fixed variant. Its temporary source, harness, build log and results
are under the directory recorded by `/tmp/rcppalgos-sieve-launch-probe-path`.

### Implemented fix, November 12, 2026

The maintainer approved replacing this launch/join block with the existing
`ThreadBatch`. Interval boundaries, worker selection, output slots, and sieve
kernels are unchanged. Mutable interval bounds and loop indices are captured
by value; shared data outlives the batch on normal and exceptional exits.

The durable sieve failure harness passed ASan/UBSan and TSan: 16 injected
failure scenarios, 14 gated join-before-propagation checks, and 19 independent
reference/recovery sieves. The cpp11 bridge translated 19 errors into R errors
and completed a subsequent sieve after each one. All 1,003 relevant R test
expectations passed. A bounded before/after comparison found no consistent
parallel slowdown. See [testing.md](testing.md#prime-sieve-worker-exception-handling)
for commands and measurement details.

### Decision

The maintainer reviewed the fix and validation evidence; the fix is merged
into `main`. No separate GitHub issue was opened for this finding.
The known sieve issues listed below concern separate endpoint and coverage
defects, not exception propagation.

## Package thread launchers and worker-side R errors

Status: Resolved; reviewed by the maintainer and merged into `main`

Sixteen raw launch regions could terminate R on worker or partial-launch
failures. They now use `RcppAlgos::ThreadBatch`, declared after borrowed storage
so unwinding joins workers before that storage is destroyed. Starting states
are captured by value, then moved into workers. Worker-reachable R errors throw
native exceptions with the original messages; cpp11 translates them on the calling thread.

Validation, November 17, 2026 (macOS arm64):

- 43,519 R expectations passed; a temporary runtime guard verified that no
  test batch exceeded two workers.
- ASan/UBSan and TSan passed 282 launcher failure scenarios (18 preparation)
  at two workers, plus callable-transfer, Phi and sieve failure checks.
- Earlier `R CMD check` reported no errors and two warnings for deliberately
  unbuilt vignettes. Tests ran separately; manual/vignette rebuilds were omitted.

Retained checks: `package-thread-check.cpp` plus its preparer, and focused
regressions in the existing family-specific R tests. Commands are in
[testing.md](testing.md#package-wide-launchers). The native fixture uses real
launchers with stand-in mathematical kernels; R tests exercise real native/GMP
results. This is not cross-platform validation or proof of kernel race freedom.
One-off guard/bridge diagnostics and detailed logs are kept outside the repo.
Tracked in [GitHub issue #87](https://github.com/jwood000/RcppAlgos/issues/87).

## Prime sieve hardware-sizing arithmetic on high-core-count hosts

Status: Resolved; reviewed by the maintainer and merged into `main` with the
PiPrime refactor. Ordinary hardware sizing is preserved.
Discovered during: Production integration of the refined nthPrime cost model,
November 13, 2026.

At discovery, `PrimeSieveMain()` in `src/Eratosthenes.cpp` computed
`maxCores * (1024 * 1024) * 6` using `int` intermediates before assigning it
to `std::int_fast64_t`. With a 32-bit `int`, 341 cores produces 2,145,386,496,
but 342 produces 2,151,677,952, beyond `INT_MAX`. Signed overflow therefore
occurs before the wider destination can help. The later Big-sieve expression
`(maxThreads / nThreads) * 2 * 1024 * 1024` similarly overflows when that
integer ratio reaches 1024, if that path executes.

These expressions predate the nthPrime integration. That integration reached
the launcher through a native sieve finish; it was subsequently rolled back
to the original local sieve. `primeSieve()` still reaches this launcher with
R's detected hardware limits. Ordinary ten-core validation does not exercise
these values.
The evidence here is the operand types and exact arithmetic, not a failed
end-to-end run on high-core hardware. Wider intermediates and subsequent size
and square calculations needed to be checked together when developing the fix.

The deferred records and known sieve issues #82–#84 were checked. A GitHub
issue search for `sieve overflow` returned only the separate bucket-boundary
issue #83; no matching hardware-sizing record was found. No issue was opened.
On November 15, 2026, the maintainer approved the narrow arithmetic fix while
sharing this launcher with PiPrime. The shared `SieveMain()` now widens both
hardware products, avoids squaring above the supported numeric range, and
caps Big segments so their next marking offsets fit `int`. Bucket construction
uses wide positions and a checked bucket-count conversion. The ordinary sizing
formulas are unchanged.

Review also confirmed an out-of-bounds access with independently supplied
native budgets: `PiPrimeRange(2, 50000, INT_MAX, INT_MAX, 1)` put the Medium
cutoff below the Small cutoff and ran Small beyond its prepared base primes.
The launcher now enforces Small <= Medium, preserving the transition order.
Regression checks cover these synthetic budgets without allocating huge sieve
segments; this remains synthetic validation, not a run on high-core hardware.

## Findings transferred to GitHub

The following sieve fixes and regression coverage were reviewed by the
maintainer and merged into `main` in
[PR #85](https://github.com/jwood000/RcppAlgos/pull/85). The issues preserve the
original reproductions and sanitizer evidence; follow GitHub for issue-closure
status.

- [Issue #82: Precomputed prime-table endpoint read](https://github.com/jwood000/RcppAlgos/issues/82).
  Commit [`204b6f9`](https://github.com/jwood000/RcppAlgos/commit/204b6f9)
  bounds table traversal in `PrimeSieveSmall()`, including the final prime,
  31627.
- [Issue #83: Scheduling past the final bucket](https://github.com/jwood000/RcppAlgos/issues/83).
  Commit [`3b8e620`](https://github.com/jwood000/RcppAlgos/commit/3b8e620)
  restricts `PrimeSieveBig()` scheduling to `idx < myRange`, preserving the
  inclusive upper endpoint without admitting offsets beyond it.
- [Issue #84: Insufficient base primes for an extended range](https://github.com/jwood000/RcppAlgos/issues/84).
  Commit [`11d464c`](https://github.com/jwood000/RcppAlgos/commit/11d464c)
  makes `sqrtBigPrimes()` prepare base primes using `myUpper`, including the
  extension needed by `bAddExtraPrime`.

The implementations are in `src/Eratosthenes.cpp`; permanent regression tests
are in `tests/testthat/testPrimeSieve.R`.
