# Deferred Investigations

## Unsafe floating-point to integer conversions

Status: Needs investigation
Priority: High
Discovered during: PartitionType refactoring

### Observation

SetTolerance(), CheckIsInteger(), and CheckPartition()
may cast out-of-range, NaN, or infinite double values
to std::int64_t, resulting in undefined behavior.

### Proposed investigation

- Reproduce the issue with targeted tests.
- Review all affected conversion paths.
- Evaluate std::isfinite() and std::round() approaches.
- Agree on the fix before implementation.

### Decision

Defer to a dedicated investigation and PR.

## Prime sieve worker and thread-launch exception handling

Status: Fix implemented and validated in the working tree; pending maintainer review
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

Implementation and evidence are ready for maintainer review; no GitHub issue
has been opened or closed for this finding.
The known sieve issues listed below concern separate endpoint and coverage
defects, not exception propagation.

## Package thread launchers and worker-side R errors

Status: Fix implemented and validated; pending maintainer review

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

Status: Narrow fix implemented with the approved PiPrime refactor; pending
maintainer review. Ordinary hardware sizing is preserved.
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
and square calculations should be checked together before proposing a fix.

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

The following sieve findings have fixes and regression coverage in
[PR #85](https://github.com/jwood000/RcppAlgos/pull/85). As of November 1, 2026,
the PR and issues remain open; the fixes are implemented on
`investigate-prime-sieve` and await merge. Follow GitHub for current review,
merge, and issue-closure status. The issues preserve the original reproductions
and sanitizer evidence.

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
