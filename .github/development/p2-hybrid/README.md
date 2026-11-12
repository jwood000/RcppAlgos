# Bounded streaming in P2 — November 6, 2026

This experiment compares the committed Phi-backed P2 with temporary hybrids
that stream a bounded interval of prime-count queries. Production code is
unchanged. Every build uses the adopted `Phi(x,a)` cache endpoint of 131,071.
The source snapshot is commit `8b43ed790435b22eb9ec02627379082fe7deb103`.

This page records the original `2L`/`4L`/`8L` study. The separate
[16L/32L follow-up](WIDE.md) extends the range while retaining these results.
The later [64L/128L follow-up](EXTENDED.md) looks for diminishing returns.
The [frozen dynamic-selector study](SELECTOR.md) evaluates intermediate inputs,
thread counts and arbitrary cutoffs, with first and warmed calls kept separate.
The [interval-aware v2 study](SELECTOR_V2.md) compares a short-range exception
and worker-budget rule with v1 and the fixed variants, including new inputs.
The [stream-density v3 study](SELECTOR_V3.md) targets v2's warmed regressions
with enlarged tables and tests exact fractional table sizes around its new guard.
The [matched-cutoff profile](MATCHED_PROFILE.md) explains the opposite choices
using worker timelines, cache histories and a temporary split-Phi control.
The [smaller-initial-chunks study](INITIAL_CHUNKS.md) tests first-chunk sizes
1/2/4/8 while retaining 32-query later chunks and the original single batch.

The 8L endpoint gives the strongest measured gains at `1e12`–`1e14`: warmed
complete-call paired reductions are approximately 13–33%, with first calls
also improving. At `1e10` with eight requested threads, warmed calls instead
regress about 23–24%. The results support reviewing 8L for larger workloads,
with further crossover measurements before choosing a production policy.

## Experimental implementation

Let `L = tables.limit()` (normally `floor(sqrt(x))`) and
`q[i] = floor(x / p[i])`. P2 still evaluates

```text
sum(i = a+1, ..., pi(sqrt(x))) [pi(q[i]) - (i - 1)].
```

The variants select `Q = min(k * L, largest requested quotient)`, with
`k = 2, 4, 8`, and divide the terms into three disjoint parts:

- `q > Q`: existing serial Phi evaluations inside dynamic 32-query chunks.
- `L < q <= Q`: bounded segmented sieving, visiting prime indices in reverse
  so the queried endpoints increase.
- `q <= L`: direct lookup in the existing prime-count table.

The split uses `floor(x / p) <= Q` exactly when
`p > floor(x / (Q + 1))`. Factors remain nondecreasing and prime squares count
once. Each streaming interval starts with a known `pi(start)` and counts
additional primes through each inclusive requested endpoint. The first
interval uses `pi(L)`; subsequent parallel intervals use one serial Phi seed.
This shares interval counting without assuming that the distinct outer
quotients repeat.

`P2Hybrid.cpp` adapts the retained earlier streaming prototype: one byte per
odd number, persistent next-crossing positions, and a buffer containing
`max(32768, floor(sqrt(finish)))` odd positions. It is a counting-only sieve;
it does not materialize every prime or a dense pi array through Q.

The Phi prefix completes before streaming workers launch. Its worker budget
is calculated from the original full non-table range, retaining production's
2,048-term grain, and capped by the remaining prefix terms. Removing the
cheap tail therefore does not automatically reduce workers on the expensive
prefix. Streaming uses disjoint approximately equal quotient spans, with a
1,048,576-number minimum span per worker. This scheduling choice is part of
the measured experiment; it is not established as optimal.

## Refreshed profile

`refresh.py` instruments only temporary copies. Forty fresh processes produced
120 checked P2 results, covering complete-driver and isolated-P2 contexts,
inputs `1e12`/`1e13`, and one/eight requested threads. Counter builds and timing
probes are separate. The extra timing build partitions queries by their
quotients at `2L`, `4L`, and `8L` rather than assuming those boundaries coincide
with equal index quantiles.

At warm `1e13`, deterministic work counts are:

| Measure | Meissel | Lehmer |
| --- | ---: | ---: |
| Outer Phi queries | 225,230 | 227,372 |
| Phi evaluator entries | 48,010,204 | 64,336,606 |
| Uncached nodes expanded | 4,408,008 | 6,050,405 |
| Direct pi-lookup suffix terms | 152,453,268 | 183,855,264 |
| Eligible cache lookups, all hits | 39,194,188 | 52,235,796 |
| Cache misses / mutex attempts / stores | 0 / 0 / 0 | 0 / 0 / 0 |

These counts overlap and must not be added as disjoint work categories. They
agree between warmed repeats, one/eight workers, and isolated/complete-driver
contexts. Relative to the [old 65,535-cache profile](../p2-profile/README.md),
evaluator entries fall 29.73%/31.39%, expanded nodes 41.24%/43.40%, and suffix
terms 18.71%/18.88%. These are work reductions, not runtime percentages.

Serial diagnostic Phi-time shares in the proposed streaming tails:

| Endpoint | Meissel queries | Meissel Phi time | Lehmer queries | Lehmer Phi time |
| --- | ---: | ---: | ---: | ---: |
| `2L` | 47.86% | 12.57% | 47.41% | 9.68% |
| `4L` | 72.99% | 25.89% | 72.30% | 19.89% |
| `8L` | 86.20% | 41.06% | 85.39% | 31.50% |

The first index quarter accounts for about 72.42%/78.63% of summed serial Phi
time. Eight-worker duration spreads are approximately 0.073/0.060 ms; these
are durations, not absolute completion timestamps. The tail shares include
probe overhead and are not predicted speedups: streaming and its launch/seed
costs still have to be paid. Uninstrumented timings determine the result.

Isolated first P2 calls allocate 1,454 cache pages. Following the complete
driver's Phi phase, first P2 adds 284 pages at `1e12` for both methods, but only
zero/one pages at `1e13`. This is why complete-driver phase order matters when
comparing first calls.

The uninstrumented refreshed baseline at warm `1e13` still attributes most
runtime to P2:

| Method | Requested threads | P2 ms | Complete ms | P2 share |
| --- | ---: | ---: | ---: | ---: |
| Meissel | 1 | 191.833 | 283.502 | 67.67% |
| Meissel | 8 | 42.807 | 72.463 | 59.07% |
| Lehmer | 1 | 250.563 | 288.545 | 86.84% |
| Lehmer | 8 | 58.772 | 76.180 | 77.15% |

Shares are ratios of medians. The phase attribution is consistent with the
earlier bottleneck diagnosis, despite the substantial reduction in Phi work.

## Timing methodology

`bench.py` builds four uninstrumented executables from the same frozen
production sources and common object files. Only P2 differs. Compiler flags
include C++11, `-O3`, `-pthread`, debug symbols and frame pointers. The run used
Apple Clang 17 on arm64 macOS 15.7.4.
Compiler commands, platform details and source/binary hashes are retained.

Complete-call cases use Meissel and Lehmer, `1e10`, `1e12`, `1e13`, `1e14`,
and thread budgets one/eight. Each configuration has three matched blocks;
each block runs all four variants in randomized order, in separate fresh
processes. Each process performs one first and three warmed calls. Case order
is also randomized. The harness follows the production driver order: local
tables, PhiParallel, P2, optional P3, then table destruction. Timings exclude
R conversion and process startup.

The isolated-P2 matrix uses `1e12`, `1e13`, `1e14` with the same methods,
budgets, blocks and repetitions. These processes start with an empty Phi
cache. Use their `p2_ms` column: isolated `total_ms` includes table setup only
on repetition zero and destruction only on the final repetition. Complete
calls rebuild tables every time and retain only the process-wide Phi cache.

Warm summaries take a median within each process, then a median across the
three processes. Percentage changes take candidate/baseline ratios within
each matched block before taking their median. Displayed ratios need not equal
the ratio of separately displayed medians. Ranges of paired changes describe
observed variability, not confidence intervals. Builds, diagnostics, memory
probes and sanitizer runs finish before timing starts.

## Complete-call results

All 768 timed complete counts matched their known values, and all P2 values
agreed across variants, repeats and worker budgets. `8L` gave the largest
median reduction of the tested variants in every `1e12`–`1e14` configuration,
for both first and warmed calls. Every matched `8L` comparison at these inputs
improved; the smallest individual improvement was approximately 10.6%.

Times below are milliseconds, shown as **baseline / hybrid 8L**. The final
column uses matched process ratios, as described above.

| Input | Method | Threads | First ms | Warm ms | Warm paired change |
| --- | --- | ---: | ---: | ---: | ---: |
| `1e10` | Meissel | 1 | 4.039 / 3.258 | 0.999 / 0.816 | -18.8% |
| `1e10` | Meissel | 8 | 1.614 / 1.686 | 0.587 / 0.718 | +24.3% |
| `1e10` | Lehmer | 1 | 4.241 / 3.329 | 1.086 / 0.909 | -17.6% |
| `1e10` | Lehmer | 8 | 1.724 / 1.935 | 0.639 / 0.788 | +22.9% |
| `1e12` | Meissel | 1 | 87.494 / 64.964 | 38.476 / 30.803 | -19.8% |
| `1e12` | Meissel | 8 | 24.364 / 18.868 | 12.020 / 10.557 | -13.4% |
| `1e12` | Lehmer | 1 | 86.590 / 64.031 | 39.026 / 31.273 | -20.1% |
| `1e12` | Lehmer | 8 | 22.996 / 17.570 | 11.928 / 10.183 | -13.7% |
| `1e13` | Meissel | 1 | 422.520 / 311.010 | 283.502 / 211.740 | -25.3% |
| `1e13` | Meissel | 8 | 108.007 / 83.705 | 72.463 / 57.614 | -21.1% |
| `1e13` | Lehmer | 1 | 428.863 / 316.013 | 288.545 / 218.615 | -24.1% |
| `1e13` | Lehmer | 8 | 109.580 / 82.345 | 76.180 / 58.995 | -22.1% |
| `1e14` | Meissel | 1 | 2500.188 / 1738.933 | 2285.870 / 1539.333 | -32.7% |
| `1e14` | Meissel | 8 | 594.751 / 437.738 | 560.094 / 385.460 | -28.5% |
| `1e14` | Lehmer | 1 | 2580.613 / 1840.982 | 2325.568 / 1623.594 | -31.6% |
| `1e14` | Lehmer | 8 | 633.039 / 443.807 | 588.671 / 398.917 | -31.2% |

The smaller endpoints also improve the larger cases, but leave more queries
with Phi. Across `1e12`–`1e14`, warmed complete-call paired changes range from
about -2% to -12% for `2L`, -7% to -21% for `4L`, and -13% to -33% for `8L`.
The raw and summary CSV files retain every configuration and the observed
minimum/maximum paired change.

The small parallel case is a real tradeoff. At `1e10`, eight requested threads
use four Phi workers and one streaming worker. All three endpoints regress
in warmed complete-call median comparisons; `8L` adds about 0.13–0.15 ms.
The first-call comparisons at that size are more variable and do not establish
a reliable improvement. Serial calls at `1e10` improve. A universal endpoint
or unconditional switch is not supported by this matrix.

## Isolated P2 and a separate confirmation

All 576 isolated timed P2 results matched the reference values derived from
the checked complete identities. Unlike the complete calls above, isolated
first P2 has no preceding Phi phase to warm its cache. Eight-thread P2 times
are shown below as **baseline / hybrid 8L**, in milliseconds:

| Input | Cutoff | First P2 ms | Warm P2 ms | Warm paired change |
| --- | --- | ---: | ---: | ---: |
| `1e12` | Meissel | 13.474 / 7.220 | 5.368 / 3.865 | -28.4% |
| `1e12` | Lehmer | 18.052 / 11.835 | 7.589 / 5.827 | -23.2% |
| `1e13` | Meissel | 67.962 / 36.838 | 43.241 / 28.794 | -33.4% |
| `1e13` | Lehmer | 83.834 / 58.893 | 55.985 / 40.990 | -26.8% |
| `1e14` | Meissel | 393.848 / 194.985 | 356.688 / 178.028 | -50.1% |
| `1e14` | Lehmer | 529.137 / 328.284 | 490.616 / 298.893 | -38.5% |

One original serial Lehmer `1e12` hybrid process was slower than its peers:
its warm median was 39.700 ms versus its paired baseline's 32.437 ms (+22.4%).
The other two pairs improved approximately 20.7% and 24.0%. No cause was
identified. The slow process is retained in the primary CSV and summaries.

A separate five-pair confirmation (`confirm.py`) reused the same baseline and
8L executables, randomized their order in each pair, and checked another 40
P2 results. Every pair improved. First-call paired changes range from -34.5%
to -30.2%, median -32.2%; warm changes range from -25.8% to -21.2%, median
-24.2%. The confirmation supports the direction of the gain without replacing
the original variability or establishing a universal runtime guarantee.

## Interpretation and next decision

Bounded streaming is a useful candidate for the larger measured workloads.
At `1e13`, it removes roughly 86% of outer Phi queries with the 8L endpoint,
but leaves the costly leading queries and the complete driver's initial Phi
phase. The measured benefit therefore differs between isolated P2 and complete
prime counting. The experiment establishes gains for this implementation and
matrix, not a speedup from the algorithm name alone.

The 8L endpoint deserves review for larger calls. The small parallel regression
supports retaining a Phi-only path where streaming does not amortize its work.
The matrix is too coarse to choose a final crossover threshold, thread policy,
or behavior for arbitrary `semiPrimeCount(n, a)`. A targeted follow-up around
the crossover would be appropriate before adopting a dispatch rule. Production
P2, Phi and the R method defaults remain unchanged pending that discussion.

## Working storage

`memory.cpp` and `build_memory.py` use separately compiled probe variants;
timing builds contain no memory-probe writes. All 96 first/warm complete-count
rows and query/worker invariants passed. Buffer capacities are identical on
first and warmed calls for the same configuration.

The tracked scratch upper bound is 33,464–41,672 bytes serial and at most
318,568 bytes (about 311 KiB) across the eight-thread cases. It is the maximum
of prefix subtotal capacity and streaming subtotal capacity plus all stream
sieve/next-crossing capacities. Worker allocations need not all overlap, so
this is a conservative payload bound. It excludes ThreadBatch vectors,
thread stacks, allocator overhead, the diagnostic snapshot, the shared Phi
cache, existing tables and their construction scratch. It is not peak RSS.

For `2L`/`4L`/`8L`, actual streaming workers at `1e12` are 1/2/6; at `1e13`
they are 3/8/8; at `1e14` all use eight. The Phi prefix retains eight workers
in those cases. Serial streaming uses the table seed; parallel streaming uses
one table seed and `workers - 1` Phi seeds in this matrix.

## Correctness and scope

Baseline and each hybrid passed 60,106 checks under both ASan/UBSan and TSan,
without compiler warnings or sanitizer diagnostics. Checks cover independent
small factorization counts, exact and oversized contexts, quotient split and
segment boundaries, prime powers, nonpositive and positive worker budgets,
independently checked multiworker streaming, concurrent calls sharing tables,
and known complete Meissel/Lehmer counts through `1e12`. Larger timed calls
are also checked against known complete counts and matching P2 results.

That total belongs to the original fixture snapshot. The 16L/32L follow-up
added 108 boundary checks; the 64L/128L follow-up adds another 108, bringing
the current checker to 60,322 checks per binary. The optional large-coordinate
stream checks are described in that later report. The original runner defaults
still select baseline, `2L`, `4L`, and `8L`.

The existing ThreadBatch failure-injection coverage remains separate; no new
failure injection was added. These experiments do not establish performance
for arbitrary user-supplied `a`, the full 53-bit range, other CPUs, or a final
dispatch policy. No production code, R interfaces, or package tests change,
and no full R package check was run for these standalone experimental copies.

## Reproduction and retained artifacts

Run from the repository root. Use a fresh temporary directory and finish
validation/probes before starting either timing command:

```sh
p2_hybrid_dir=$(mktemp -d)
mkdir "$p2_hybrid_dir/source"
git archive 8b43ed790435b22eb9ec02627379082fe7deb103 src inst/include | \
  tar -x -C "$p2_hybrid_dir/source"
git rev-parse 8b43ed790435b22eb9ec02627379082fe7deb103 > "$p2_hybrid_dir/commit.txt"
python3 .github/development/p2-hybrid/refresh.py build "$p2_hybrid_dir"
python3 .github/development/p2-hybrid/refresh.py run "$p2_hybrid_dir"
python3 .github/development/p2-hybrid/bench.py build "$p2_hybrid_dir"
python3 .github/development/p2-hybrid/validate.py "$p2_hybrid_dir"
python3 .github/development/p2-hybrid/build_memory.py build "$p2_hybrid_dir"
python3 .github/development/p2-hybrid/build_memory.py run "$p2_hybrid_dir"
python3 .github/development/p2-hybrid/bench.py phases "$p2_hybrid_dir"
python3 .github/development/p2-hybrid/bench.py p2 "$p2_hybrid_dir"
# Separate confirmation; it does not overwrite the primary matrix.
python3 .github/development/p2-hybrid/confirm.py "$p2_hybrid_dir"
```

The [results directory](results/) retains raw complete/isolated timings,
summaries with paired ranges, the separate confirmation, diagnostic summaries,
quotient bands, memory/work observations, validation totals and build metadata.
The production-preservation audit verified 392 source, header, R, manual, test
and package metadata files against their pre-experiment hashes without changes.

The original full artifact directory is recorded by
`/tmp/rcppalgos-p2-hybrid-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-p2-hybrid-fd4j7ai5`.
It contains the immutable sources, executables, raw diagnostic logs, build
commands, compiler logs and sanitizer logs, in addition to result tables.
