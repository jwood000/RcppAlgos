# Number-theory implementation notes

These notes preserve the reasoning behind the current design. Public argument
definitions and supported ranges live in the R help pages; runnable checks live
in [testing.md](testing.md). Measurements below are historical evidence, not
current-checkout or cross-platform certification.

## Source map

| Responsibility | Main source |
| --- | --- |
| R adapters and conversion | `R/NumberTheory.R`, `PrimeCount.cpp`, `PartialSieve.cpp`, `NthPrimeR.cpp`, `PrimeSieve.cpp` |
| Prime-count formulas and automatic policy | `PiLegendre.cpp`, `PiMeissel.cpp`, `PiLehmer.cpp`, `PiAuto.cpp`, `NumbersUtils/PrimeCount.h` |
| Prepared tables, serial Phi, shared cache | `Phi.cpp`, `NumbersUtils/Phi.h` |
| Outer Phi scheduling | `PhiParallel.cpp` |
| Two- and three-prime products | `P2.cpp`, `P3.cpp`, `NumbersUtils/PartialSieveUtils.h` |
| Inclusive interval counting | `PiPrime.cpp` |
| Shared Small/Medium/Big sieve | `Eratosthenes.cpp` |
| Exact one-based prime selection | `NthPrime.cpp` |
| Thread ownership and failure propagation | `ThreadBatch.h` |

Source filenames are under `src/`; headers are under `inst/include/`.
`CppConvert/Constants.h` owns `Significand53`, its integer counterpart
`Significand53Int`, and `Sqrt53Max`. Algorithm-specific cache and scheduling
constants stay with their implementations.

## Mathematical contracts

Prime and pi tables have an inclusive limit. The prime vector is one-based,
with a sentinel at zero; `pi[t]` is the number of primes at most `t`.
`PhiTables` owns these vectors and can be shared read-only by concurrent
evaluations. A borrowed table must cover the requested calculation, even when
it was prepared for a different input. Each algorithm derives its cutoff from
the current input, not from the full extent of a reused table.

`Phi(x, a)` counts integers in `[1, x]` with no prime factor among the first
`a` primes. It is not Euler's totient and is not a count of surviving primes:
1 survives, and composites can survive. Its recurrence is

```text
Phi(x, 0) = x
Phi(x, a) = Phi(x, a - 1) - Phi(floor(x / p[a]), a - 1).
```

For `a <= 6`, divisibility is periodic modulo the product of those primes.
The tiny tables hold cumulative survivor counts for a partial period, not
individual survivor flags. Full periods contribute their totient. One immutable
function-local object is initialized safely under C++11 and shared by all calls.

P2 counts `p[i] * p[j] <= x` with `a < i <= j`; P3 counts
`p[i] * p[j] * p[k] <= x` with `a < i <= j <= k`. Factors are nondecreasing,
repetition is allowed, and each integer is counted once. For `x >= 1`, the
counting identities are:

| Method | Cutoff | Identity |
| --- | --- | --- |
| Legendre | `a = pi(floor(sqrt(x)))` | `pi(x) = Phi(x,a) + a - 1` |
| Meissel | `a = pi(floor(cbrt(x)))` | `pi(x) = Phi(x,a) + a - 1 - P2(x,a)` |
| Lehmer | `a = pi(floor(x^(1/4)))` | `pi(x) = Phi(x,a) + a - 1 - P2(x,a) - P3(x,a)` |

The native kernels handle small and empty cases before their main sums. Integer
roots correct floating-point estimates using integer comparisons. Quotient
comparisons avoid unnecessary products that could overflow.

The direct-table suffix of a Phi sum follows from monotonicity: as its prime
increases, `q = floor(x/p)` cannot increase while `p*p` increases. Once both
table coverage and `q < p*p` hold, they continue to hold for subsequent terms.
The suffix includes the transition term; it must not skip or duplicate it.

P2's exact split between Phi, streaming, and table lookup is documented with
its current performance policy in [CONSOLIDATED.md](p2-hybrid/CONSOLIDATED.md).

## Cache and worker lifetime

The shared Phi cache accepts positive `x <= 131071` and `a < 100`; `a <= 6`
returns through the tiny evaluator first. Monotonicity gives
`1 <= Phi(x,a) <= Phi(131071,6) = 25142`, so every stored answer fits in 16 bits.
Zero is the miss sentinel. Each row has 32 lazily allocated pages of 4,096
explicitly zero-initialized atomic cells.

Row mutexes protect allocation only. Publishing a page pointer with release
and reading it with acquire makes construction visible to readers. Each cell
contains a complete independent answer, so relaxed atomic loads/stores suffice.
Two workers may compute and store the same answer; equal writes would still
be a data race without atomics. A missed value only causes recomputation.
Failed optional page allocation does not invalidate the calculated result.

Published pages stay alive and stationary until the shared cache is destroyed;
all evaluations must have ended before destruction. Atomic pointers alone do
not solve reclamation. Atomics were lock-free on the measured arm64 host, but
correctness does not depend on that property or on a particular atomic size.

`RcppAlgos::ThreadBatch` owns threads and exception slots, not borrowed tables
or result buffers. Referenced storage must outlive the batch. Worker failures and later
launch failures are propagated only after all started workers join. Declare
worker-referenced storage before the batch so stack unwinding preserves that
ordering. Workers do native C++ work; they must not call R error APIs or use
`cpp11::stop` as a worker error boundary.

The current Phi scheduler retains create/join batches; P2 uses a dynamic work
queue for its own Phi prefix. The same helper now protects package-wide
launchers, including caller-side preparation of the next worker's starting
state. That migration and its review status are recorded in
[DEFERRED_ISSUES.md](DEFERRED_ISSUES.md).

## Memory

The dense pi payload is exactly `4 * (limit + 1)` bytes. At
`limit = Sqrt53Max = 94906265`, this is about 362 MiB before the prime vector,
cache, allocator overhead, and workers. The three prime-count methods currently
share this square-root table design. The often-quoted cube-root memory bound
for a different Meissel implementation does not describe this implementation.

The cache's configured cell payload is 25 MiB when atomic cells occupy two
bytes. Allocation is lazy. The doubled-cache study observed about 11.39 MiB
including page payload and the cache object, about 7.21 MiB above the old
configuration. These are payload observations, not process peak RSS; see the
[cache study](p2-profile/README.md) for input sizes and first-call costs.

Count-only sieving avoids retained prime output, but still needs base primes,
offsets, segment flags, and Big-sieve buckets. `PiPrimeRange` limits each counting
window to one billion integer positions; long intervals can rebuild preparation
between windows. The final nthPrime generation cap likewise limits output, not
total scratch memory. Releasing tables before sieving ends their ownership but
does not guarantee that the allocator returns memory to the operating system.

## nthPrime and exact interval updates

Ranks through 25 use a lookup. Larger ranks start from
`n * (log(n) + log(log(n)) - 1 + (log(log(n)) - 2) / log(n))`.
At most two whole-prefix counts establish an exact pair `(anchor, pi(anchor))`.
The first count hands off to sieving when the estimated displacement
`abs((n - pi(anchor)) * log(anchor))` fits the following allowance:

| Search location, not rank | Allowance in integer positions |
| --- | ---: |
| Below `1e10` | `1e6` |
| Below `1e12` | `1e7` |
| Below `1e14` | `1e8` |
| Otherwise | `1e9` |

These broad bands are work limits, not bounds on the estimate's error.
Prepared tables are reused for prefix counting, then released. Up to two
interval-count corrections can improve the anchor without another prefix count:

```text
next > anchor: pi(next) = pi(anchor) + PiPrimeRange(anchor + 1, next)
next < anchor: pi(next) = pi(anchor) - PiPrimeRange(next + 1, anchor).
```

Exact counts maintain a bracket containing the requested prime; the density
estimate only chooses a location in that bracket. Forward selection takes the
`(n-count)`th prime strictly above the anchor. Backward selection takes the
`(count-n+1)`th prime at or below it. Equality still requires selection because
the anchor may be composite.

`PrimeSieveMain<int64_t>` generates adjacent windows of at most one million
positions. Worker rows are numerically ordered and can be walked directly in
either direction. Empty outer output storage is valid: the launcher resizes
it before workers start. Exhausted rows/windows reduce the exact remaining
rank; a poor width estimate changes runtime, not the answer.

## Sieving and nthPrime measurements

The November 16 nthPrime comparison used baseline commit `a3db472`, Apple
clang 17, C++11 `-O3`, and macOS arm64. Eight ranks and one/eight requested
threads gave 192 checked calls across 64 fresh processes: two processes per
configuration/variant, three calls per process, alternating variant order.
First calls were separated from the two warmed calls. Actual R hardware budgets
were `maxThreads = 10`, `maxCores = 1`. No builds or tests overlapped timing.
All samples remain in [nth-prime-native-timings.csv](nth-prime-native-timings.csv).

Representative medians, earlier implementation to native sieve:

| Rank | Threads | First ms | Warmed ms | Peak RSS MiB |
| --- | ---: | ---: | ---: | ---: |
| `1e8` | 8 | 1.691 → 1.666 | 0.509 → 0.586 | 9.54 → 10.05 |
| `1e9` | 1 | 8.390 → 6.250 | 2.663 → 2.862 | 13.85 → 14.54 |
| `1e12` | 8 | 358.165 → 253.904 | 314.382 → 214.479 | 42.73 → 44.88 |
| `29844570422669` | 1 | 31902.093 → 19715.180 | 29231.637 → 19497.233 | 150.77 → 183.02 |
| `29844570422669` | 8 | 6836.402 → 4610.024 | 7040.813 → 4713.181 | 152.93 → 188.30 |

The final rank returns `999999999999989`, just below `1e15`. Large warmed
calls improved about 32–34%, with roughly 32–35 MiB more peak RSS near that
input. Smaller warmed regressions include the two rows above and `1e9` with
eight threads (1.630 → 1.766 ms); their first calls improved. Rank `1e11` was
roughly unchanged serially. Two process pairs on one machine establish neither
portable crossover points nor behavior at the maximum supported rank.

The shared sieve's separate November 15 study compared generation with
count-only output using three fresh processes, four calls each, fixed hardware
budgets of eight, rotated order, and separate processes for RSS measurements:

| Inclusive interval | Threads | Generate/count warmed ms | Generate/count peak MiB |
| --- | ---: | ---: | ---: |
| `[1e12, 1e12 + 1e7]` | 1 | 5.50 / 5.25 | 5.20 / 2.41 |
| `[1e15, 1e15 + 1e9]` | 1 | 741 / 735 | 276 / 89.4 |
| `[1e15, 1e15 + 1e9]` | 8 | 193 / 199 | 713 / 506 |

Count-only primarily saves output memory; it is not invariably faster. The last
interval contains `1e9 + 1` integers, so its final endpoint starts a second
counting window. The first-call ordering differs from the warmed ordering.
These figures are preserved from the original report; raw artifacts were in
`/private/tmp/rcppalgos-piprime-refactor-ohfp5ctm/` and are not a durable repository
dependency. Earlier sieve optimization data remain in
[prime-sieve-benchmarks](prime-sieve-benchmarks/README.md).

## Decisions not to repeat without new evidence

- The original shared cache had unsynchronized reads concurrent with writes.
  A global lock removed the race but harmed parallel work. Atomic cells and
  safely published pages preserve sharing without locking every lookup.
- A per-call Phi worker pool was rejected. At `1e14`/eight workers it reduced
  875 launches to eight, but measured time was about 1% worse; four-worker
  calls improved about 3%. Benefits were small and inconsistent relative to
  the added lifecycle/synchronization code. The full November 3 study is in
  `36303eb:.github/development/testing.md`.
- Shared-Phi integration initially slowed warmed calls. The retained improvement
  validates contracts at the public boundary and uses established invariants
  inside recursion, rather than repeatedly validating recursive subproblems.
- Enlarging the Phi cache improves many warmed calls at the cost of memory
  and some first calls. It is an explicit tradeoff, not free memoization.
- Wider P2 streams and more complicated selectors are not uniformly faster.
  First-query scheduling changed the comparisons; the final policy keeps
  conservative fallbacks instead of fitting every measured input.
- nthPrime hardware calibration and a fitted runtime cost model were removed
  in favor of broad magnitude bands and exact sieve correction. Earlier studies
  are outside the repository at
  `/Users/owacasa/RcppAlgos-archives/nth-prime-studies-2026-11-15.tar.gz`;
  the recorded SHA-256 is
  `e4fc37744974a3bf4fce2d4cf7578ad002d9c0fc9528f2ff47c5e7500a59cb34`.
  This is a historical archive reference, not a required build input.

Future work should start from the recorded limitations and a concrete question.
Passing regression tests does not prove a formula, and a benchmark win does not
establish a universal dispatch rule.
