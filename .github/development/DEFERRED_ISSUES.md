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
