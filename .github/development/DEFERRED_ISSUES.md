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

## Prime sieve lookup-table endpoint reads past the table

Status: Confirmed defect; fix pending discussion
Discovered during: Prime sieve performance investigation

### Observation and evidence

`PrimeSieveSmall()` in `src/Eratosthenes.cpp` takes the precomputed-table
path when `maxNum <= smallPrimeBase.back()`. Its output loop checks
`smallPrimeBase[ind] <= maxNum` without checking the array boundary.
At `maxNum == 31627`, it emits the final entry and then reads element 3402,
one past the 3402-element table in `inst/include/NumbersUtils/PrimesSegSieve.h`.
This is undefined behavior, even if an ordinary build returns expected results.

A standalone build of the unchanged sieve source with Apple Clang 17,
`-std=c++11 -O1 -g -fsanitize=address -fno-omit-frame-pointer`, confirmed a
`global-buffer-overflow` (4-byte read immediately after `smallPrimeBase`).
The reproduction called `PrimeSieve::PrimeSieveMain<int>` with bounds
`1, 31627`, one thread, and initially empty output vectors. The R call
`primeSieve(31627)` reaches the same path; the sanitizer reproduction used
the native entry point rather than an installed R package.

### Follow-up

Review the table loop's termination and add coverage for the exact final
table entry. The existing `primeSieve(31625, 31630)` test takes a different
path and does not cover this boundary. No implementation fix or GitHub
issue was created during the performance investigation.
