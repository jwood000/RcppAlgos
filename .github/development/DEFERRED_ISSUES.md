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

## Prime sieve bucket scheduling can exceed the final bucket

Status: Confirmed defect; fix pending discussion
Discovered during: Prime sieve optimization investigation

`PrimeSieveBig()` in `src/Eratosthenes.cpp` sets `maxIndex = myRange + 1`,
although the final requested offset is `myRange - 1`. An admitted strike
past the upper bound can therefore index beyond `myBuckets`.

A native reproduction calls `PrimeSieve::PrimeSieveMain<double>` with bounds
`1037517761278999, 1037517761279999`, initially empty output vectors,
`Parallel = false`, and `nThreads = 1, maxThreads = 10, maxCores = 1`.
The segment size is 19,219,200 and only one bucket is allocated.
The sieving prime 10,000,019 schedules its multiple
`10000019 * 103751579 = 1037517761280001` at offset `segSize + 1`,
leading to `myBuckets[1].push_back(...)`.

Apple Clang 17 with AddressSanitizer and UndefinedBehaviorSanitizer confirmed
a `container-overflow` in the unchanged native implementation. This was a
native reproduction; public R dispatch depends on detected core/thread settings.

Review the inclusive scheduling limit and add regression coverage around
aligned final-segment boundaries. Keep this correctness fix separate from
performance experiments. No implementation fix or GitHub issue was created.

## Prime sieve base-prime preparation ignores the extended output bound

Status: Confirmed defect; source fix pending discussion
Discovered during: Maintainer's buffer optimization experiments

`sqrtBigPrimes()` in `src/Eratosthenes.cpp` prepares `sqrtSievePrimes` from
`floor(sqrt(sqrtBound))`, but with `bAddExtraPrime` it then asks
`PrimeSieveSmall()` to sieve through `sqrtBound + 225`.

For `primeSieve(1.5e9, 1.5e9 + 1e5)`, these bounds are 38,731 and 38,956.
The base-prime list ends at its sentinel 197. Since 197^2 = 38,809 is within
the extended interval, the final-segment loop consumes that sentinel and
reads past the vector while looking for the next prime.

A bounded standalone ASan/UBSan reproduction confirmed a container-overflow
in both the working source and committed source at `40d06ef`, before entering
`PrimeSieveMedium()`. Thus the buffer changes did not introduce the defect.
The maintainer observed a hang in R; the native sanitizer build aborts.

A temporary copy computing the base-prime bound from `myUpper` instead of
`sqrtBound` passed both reported ranges under sanitizers, returning 4,703
primes for the 1.5e9 interval and 4,745 for the 2e9 interval, matching counts
from `isPrimeRcpp()`. Review this bound correction and add regression coverage
before applying it to the repository. No source fix or GitHub issue was made.
