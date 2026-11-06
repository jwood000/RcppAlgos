# P2 profiling studies — November 6, 2026

These studies profiled Phi-backed P2 with the original `x <= 65535` cache and
then tested a doubled bound in isolated copies. The approved `x <= 131071`
bound has now been applied to production. All measurements below retain their
original baseline and meaning; the R method-selection policy is unchanged.

The dominant cost is evaluating Phi for the outer quotients.
The existing dynamic scheduling distributes that uneven work effectively in
the measured cases. The initial investigations made no production changes.
The follow-up below measures reuse among states excluded by the original cache,
including subtree work and possible storage costs.

For historical reproduction, prepare the pre-adoption sources once from the
repository root before using the command blocks below:

```sh
p2_baseline_source=$(mktemp -d)
git archive 0bd47cd07ae0c747b51f4696b3949994553f1e65 src inst/include | \
  tar -x -C "$p2_baseline_source"
```

This commit preserves the computational sources used in the studies; one
header differs from the original snapshot only in API documentation comments.
The original artifact snapshots remain available for literal file-hash
reproduction. The reuse and cache-double preparers deliberately reject the
new production cache shape instead of silently changing historical experiments.

## Complete-call attribution

These are **uninstrumented warm medians in milliseconds** at `x = 1e13`,
including table preparation and destruction in the complete call. Each
configuration used three fresh processes with one initial and three warm calls.
P2's share is the ratio of the two medians.

| Method | Requested threads | P2 | Complete count | P2 share |
| --- | ---: | ---: | ---: | ---: |
| Meissel | 1 | 261.202 | 388.632 | 67.2% |
| Meissel | 8 | 57.469 | 96.562 | 59.5% |
| Lehmer | 1 | 345.011 | 394.354 | 87.5% |
| Lehmer | 8 | 74.787 | 94.964 | 78.8% |

P2 itself speeds up about 4.55x/4.61x with eight requested threads in these
Meissel/Lehmer cases. The full matrix also covers `1e10` and `1e12`, and four
requested threads. At `1e10`, P2 caps an eight-thread request at four workers.
At `1e12`, warm P2 accounts for 66.5%/84.8% of serial Meissel/Lehmer time and
50.2%/68.9% at eight requested threads.

The harness follows each driver's actual phase order: tables, PhiParallel,
P2, then P3 for Lehmer, followed by table destruction. It does not estimate
driver time by adding independently warmed kernel benchmarks.

## Where P2 spends its time

Five-second macOS stack samples of the uninstrumented warm P2 loop at `1e13`
placed 4,300/4,309 Meissel and 4,283/4,292 Lehmer samples (about 99.8%) in
`PhiEvaluator::evaluate`. This includes inlined cache operations and lookup
shortcuts; it is not a measurement of arithmetic alone.

In the eight-worker Lehmer sample, 33,284/33,396 worker samples (about 99.7%)
were in the evaluator. The calling thread was waiting in `std::thread::join`.
That expected wait is not evidence that workers are blocked on a cache mutex.
Sampling is statistical, and optimized inline source-line attribution should
not be treated as exact instruction-cost accounting.

Separate diagnostic counters give the following work totals for one warm
P2 call at `1e13`. They agree between isolated P2 and the complete phase
sequence, between one and eight workers, and between the two warm repeats.

| Work counter | Meissel cutoff | Lehmer cutoff |
| --- | ---: | ---: |
| Outer Phi queries | 225,230 | 227,372 |
| Phi evaluator entries | 68,318,691 | 93,769,798 |
| Uncached nodes expanded | 7,501,776 | 10,689,403 |
| Direct pi-lookup terms inside Phi's suffix loop | 187,553,851 | 226,651,637 |
| Eligible descendant cache lookups, all hits | 53,315,139 | 72,390,992 |
| Eligible cache misses / cache stores / mutex attempts | 0 / 0 / 0 | 0 / 0 / 0 |

These counters overlap and should not be added as disjoint work categories.
Evaluator entries include recursive calls that return cached or tiny results;
expanded nodes exclude those returns. The direct lookup terms are inside Phi,
not the separate direct-table suffix in P2.

This makes cache **misses and mutex contention an unlikely explanation for
the warm bottleneck in these cases**. The remaining work includes millions
of cache-ineligible nodes and hundreds of millions of loop terms. Atomic
read costs, memory bandwidth, and cache locality have not been isolated with
hardware counters; the observed speedup alone does not identify their impact.

## Why the outer queries are costly

For the drivers' minimal table limit `L = floor(sqrt(x))`, P2's direct-table
suffix has at most one term: only `p[i] = L`, with `L` prime and
`x < L * (L + 1)`. For all power-of-ten inputs in this matrix it is empty.
Consequently every outer term invokes Phi. When `L >= 65535`, its quotient
`q > L` cannot fit the cache's `q <= 65535` range. The counters confirm zero
eligible top-level cache queries; the cache benefit comes from descendants.

The outer quotients are distinct: for integer denominators at most `sqrt(x)`,
increasing the denominator strictly decreases `floor(x/p)`. Memoizing entire
outer answers therefore cannot remove repeated queries within one P2 call.
There may still be reusable descendant states outside the current cache bounds;
this profile does not measure their unique-key distribution.

The cost is concentrated near the beginning. Separate timing probes at warm
`1e13` attribute approximately 45.5%/57.7% of serial summed Phi elapsed time
to the first eighth of the Meissel/Lehmer index range, and 65.1%/73.0% to the
first quarter. Lehmer adds only 2,142 outer queries (less than 1% more) compared
with Meissel, but generates about 37% more evaluator entries. Both use the same
table limit and identical queries on their common index suffix.

P2 computes `sqrt(q)` to choose Phi's cutoff, and an uncached top-level Phi
evaluation computes that root again. The outer quotient/root/pi-lookup probes
sum to only about 3.5–3.7 ms in the serial `1e13` cases; clock overhead biases
such small components upward. Removing that duplicate work is a possible small
improvement, but the measurements point to the recursive work as the larger
opportunity.

## Scheduling and cold-cache behavior

The first chunks are expensive, supporting dynamic assignment rather than
equal-size static partitions. Eight-worker warm runtime durations are closely
matched: the per-call duration spreads are about 0.05–0.06 ms for Meissel and
0.06–0.09 ms for Lehmer at `1e13`. These are worker **durations**, not absolute
completion timestamps. Summed per-query/worker elapsed time includes parallel
overlap and descheduling; it is not process CPU time.

The largest measured 32-query chunks are about 1.9–2.0 ms for Meissel and
11.8–13.1 ms for Lehmer. The serial path processes its entire range as one
chunk, so its maximum-chunk measurement is not comparable. There is no evidence
here that a scheduling redesign would address the main warm cost.

Cache starting state matters. At `1e13` with one worker:

| P2 context | Meissel first / warm | Lehmer first / warm |
| --- | ---: | ---: |
| Isolated P2, initially empty cache | 316.512 / 262.564 | 407.508 / 346.683 |
| P2 after the complete driver's Phi phase | 271.072 / 261.202 | 386.232 / 345.011 |

Isolated initial calls allocate 533 cache pages in the diagnostic matrix
(about 4.16 MiB of cells on this platform). The preceding Phi phase reduces
P2's new-page count to zero for Meissel and one for Lehmer at `1e13`; at `1e12`
both add three pages. Warm P2 adds none. The immutable table vectors themselves
have approximately 13.81 MiB of capacity at `1e13`, shared by all workers.
These figures do not establish maximum-range or R-process peak memory.

## Initial-study measurement limits

Count instrumentation makes warm P2 about 64–81% slower; use those runs only
for work counts. Timing probes compile recursive counters out. Relative to
the baseline, their warm timings differ by +3.4–11.3% at `1e12` and −1.6% to
+1.8% at `1e13`. Negative differences reflect run/compiler-layout variability,
not negative probe overhead. Only uninstrumented timings are used for the
complete-call and cold/warm tables above. Diagnostic quantile timings are
approximate and should guide investigation rather than predict a speedup.

The follow-up below measures potential reuse among **cache-ineligible Phi
subproblems** generated by the expensive leading queries. The current scheduler
and cache remain unchanged.

## Cache-ineligible state reuse

The larger recursive subtrees do not repeat in these measurements. Reuse is
concentrated in relatively small descendants close to the existing cache's
argument limit. This supports testing a modest extension to the paged cache
before considering a general key/value memo table.

Here `x` is P2's input and `(y, a)` is an individual Phi state. An exact serial
trace records every private evaluator state with `a > 6` and either `y > 65535`
or `a >= 100`. Each fresh process performs one unrecorded P2 warmup, followed
by one recorded P2 call. Thus the existing cache is warm; each hypothetical
additional memo starts empty. These are work counts, not benchmark timings.

| Input | Cutoff | Ineligible visits | Distinct keys | Repeated visits | Repeated share |
| --- | --- | ---: | ---: | ---: | ---: |
| `1e12` | Meissel | 879,801 | 861,244 | 18,557 | 2.11% |
| `1e12` | Lehmer | 1,186,596 | 1,142,674 | 43,922 | 3.70% |
| `1e13` | Meissel | 7,501,776 | 6,031,345 | 1,470,431 | 19.60% |
| `1e13` | Lehmer | 10,689,403 | 7,621,470 | 3,067,933 | 28.70% |

Every repeated visit crosses an outer-query boundary. There are no repeated
cache-ineligible keys within an individual outer Phi query, and no repeated
top-level keys.
At `1e13`, repeated states have `65536 <= y <= 509409` and `7 <= a <= 67`
(combining the two cutoffs). Every repeated state's subtree contains exactly
one cache-ineligible node: itself. Its children already reach cache/tiny-table
returns or direct pi shortcuts. The largest repeated subtree contains 13
evaluator entries, including itself, and at most 52 direct suffix terms.
These are measured properties of these cases, not general Phi theorems.

The expensive prefix is predominantly distinct work. At `1e13`, the first
eighth of the Meissel/Lehmer index ranges contains 3,076,238/6,156,785 ineligible
visits, of which 379,041/1,377,620 repeat an earlier key. The additional 2,142
Lehmer queries preceding the Meissel range generate 3,187,627 visits and
2,549,559 distinct keys. An unlimited added memo would remove 12.05% of that
prefix's evaluator entries and 3.13% of its suffix terms. That prefix also
populates 959,434 keys subsequently seen in the common Meissel suffix.

### Work that memoization could remove

Replay follows the original preorder traversal, publishes a value only when
its subtree completes, and skips the entire subtree on a subsequent hit.
It never credits both a hit and its skipped descendants. A hit still enters
the evaluator, so saved evaluator calls are the inclusive subtree count minus
one per hit. Saved suffix terms and expanded nodes are reported separately;
they are not summed into an invented time estimate.

For a warm `1e13` call, with an initially empty additional memo and no eviction:

| Added domain | Evaluator calls avoided, Meissel / Lehmer | Suffix terms avoided, Meissel / Lehmer | Additional paged cell payload |
| --- | ---: | ---: | ---: |
| All ineligible keys | 11.36% / 17.16% | 3.28% / 4.72% | Not estimated as a dense domain |
| `y <= 131071, a < 100` | 9.66% / 14.34% | 3.16% / 4.55% | 7.20 MiB |
| `y <= 262143, a < 64` | 11.30% / 17.05% | 3.27% / 4.72% | 21.16 MiB |
| `y <= 262143, a < 100` | 11.30% / 17.05% | 3.28% / 4.72% | 26.83 MiB |

The last domain captures 98.77%/98.58% of all repeat hits. Extending to
`y <= 1048575` captures the remaining hits, but gains little work removal
while sharply increasing the candidate population. Increasing only the row
limit while retaining `y <= 65535` captures no states in any of these runs.
For these P2 calls, when `y <= 65535` (hence within the prepared table),
the failed direct shortcut implies `y >= p[a+1]^2`; `a >= 100` would
require `y >= 547^2 = 299209`.
Do not apply that argument to every general Phi call: public calls can lack
the lookahead prime needed for the shortcut.

At `1e12`, even memoizing every excluded key would remove only 0.60%/1.09%
of evaluator calls and 0.019%/0.030% of suffix terms. A cache extension's
initialization, allocation and lookup costs could outweigh those savings.

### Memory interpretation and value bounds

A straightforward sparse representation using a 64-bit argument, 32-bit
cutoff and 64-bit answer has 20 bytes of field payload, or 24 bytes with
ordinary 8-byte alignment. Retaining all ineligible keys at `1e13` would
therefore consume about 138.05/174.44 MiB at 24 bytes per key, before hash
buckets, allocator and synchronization overhead. Even the `262143, a < 100`
domain would use 84.12/100.76 MiB with that representation. These are
illustrative layout costs, not lower bounds on all cache designs.

The table instead counts distinct `(a, floor(y/4096))` pages occupied by
candidate keys. The three bounded domains use respectively 921, 2,708 and
3,434 **additional** pages at `1e13`, for both cutoffs. At 4,096 two-byte
cells per page, these require 7.20, 21.16 and 26.83 MiB of cell payload.
They exclude existing-cache storage (about 4.16 MiB in the initial study),
page-pointer metadata and allocator overhead. Two-byte atomic cells were
verified on this platform. These are prospective occupied-page counts,
not measured allocations by an implemented larger cache or process peak RSS.

The value width can remain 16 bits for these domains, with a stronger bound
than `Phi(y,a) <= y`. Every state that reaches this cache has positive `y`
and `a >= 7`, hence for `y <= 262143`:

```
1 <= Phi(y,a) <= Phi(262143,6) = 50282 < 65536.
```

The value 50,282 was independently checked by inclusion–exclusion over the
first six primes and by a periodic survivor prefix. This preserves zero as
the missing sentinel. It does not justify arbitrary further extension with
`uint16_t`, or establish the correctness of unimplemented page-addressing
changes. Expanding the current 100-row layout fully to this endpoint would
raise its maximum cell payload from 12.5 MiB to 50 MiB across all workloads;
the occupied-page figures above describe these P2 traces only.

The next recommended experiment is a temporary doubling to `y <= 131071`,
preserving the existing paged atomic design, followed by cold/warm and
serial/parallel benchmarks. Because Phi's cache is shared, that experiment
must also check Legendre and complete Meissel/Lehmer calls, including smaller
inputs. No cache extension or optimization has been implemented in this study.

### Reuse-study validation and reproduction

All four traces preserve the expected P2 answer and exactly match the prior
independent instrumentation's query, evaluator, expanded-node and suffix
totals. Each query's recorded root is checked against its actual `(y,a)`;
root subtree costs cover the entire trace, and intervals are checked for
proper nesting. Two synthetic replay checks exercise overlapping savings
and completion-before-publication. A second fresh-process run of the largest
case reproduces every reported statistic exactly. The C++11 build has no
warnings. Capacity rejection and output-alias guards were also exercised.

The trace records occupy 32 bytes each; sorting uses a separate 4-byte index
per record, released before replay. The largest measured trace plus sorting
index uses 366.99 MiB. The driver caps traces at 12 million records, rejects
incomplete traces, and releases the dense prime tables before sorting. These
are diagnostic-buffer figures, not production memory requirements or peak
RSS. No sanitizer rerun, package check, cold-cache replay, parallel memo
simulation, eviction study, maximum-range run or actual cache-speed benchmark
was performed. In particular, serial reuse opportunities do not establish
parallel cache performance or simultaneous-worker hit availability.

[Reuse results](results/reuse-results.json) retain exact distributions, nine
domain-policy replays, work savings and page counts;
[metadata](results/reuse-metadata.json) retain compiler and source hashes.
Temporary sources, executables, logs, guard checks and per-case output are
under the directory recorded by `/tmp/rcppalgos-p2-reuse-path`.

To reproduce independently from the repository root:

```sh
p2_reuse_dir=$(mktemp -d)
cp -R "$p2_baseline_source" "$p2_reuse_dir/source"
python3 .github/development/p2-profile/run_reuse.py "$p2_reuse_dir"
```

The runner reuses these historical sources, instruments only copies, and runs
all cases sequentially. It rejects repository output directories and
pre-existing artifact symlinks. Run it separately from timed benchmarks.

## Temporary doubled-cache benchmark

The approved experiment doubles the argument endpoint from 65,535 to 131,071
in **temporary source copies only**. It retains 100 rows, 4,096-cell pages,
16-bit atomic answers, lazy allocation and the existing synchronization.
The new page-pointer array has 32 slots per row instead of 16. An independent
partial sieve and inclusion–exclusion both give `Phi(131071,6) = 25142`,
bounding every cached answer because cached evaluations have `a >= 7`.
The [original comparison patch](results/cache-double-candidate.patch) records
this experiment. Its computational change is now applied to production, with
a clarified comment that all `a <= 6` cases return before cache lookup.

### Complete-call timings

The main matrix covers all three native drivers, inputs `1e6`, `1e8`, `1e10`,
`1e12`, `1e13`, and thread budgets 1 and 8. Every configuration has three
paired baseline/candidate processes, with randomly ordered configurations and
randomized adjacent variant order. Each fresh process makes one first call
and three warmed calls. Warmed calls rebuild the prime tables but retain the
process-wide Phi cache. Native timing includes table construction/destruction
and excludes process startup and the R adapter. Timing binaries contain no
probe code. No other builds, sanitizer runs or measurements ran concurrently.

At `1e13`, milliseconds for the complete native call:

| Method | Threads | First: current / doubled | Warm: current / doubled | Warm paired runtime change |
| --- | ---: | ---: | ---: | ---: |
| Legendre | 1 | 447.02 / 419.14 | 386.07 / 278.26 | -27.93% |
| Meissel | 1 | 445.82 / 421.64 | 384.99 / 278.31 | -27.71% |
| Lehmer | 1 | 456.48 / 430.58 | 394.69 / 286.63 | -27.34% |
| Legendre | 8 | 122.85 / 116.77 | 108.00 / 80.70 | -25.08% |
| Meissel | 8 | 110.11 / 104.01 | 95.02 / 71.14 | -24.98% |
| Lehmer | 8 | 109.36 / 104.04 | 94.30 / 71.44 | -24.24% |

Displayed times are medians of per-process medians. Percentage changes are
medians of the **paired** candidate/baseline ratios; they need not equal the
ratio of the two separately displayed medians. Negative means lower runtime.
Raw data and summaries retain the range of paired changes; those ranges are
variability summaries, not confidence intervals.

Warm improvements at `1e13` hold in every process pair. First-call medians
improve by 4.94–6.25%, although one Meissel serial pair and one Legendre
eight-thread pair regress. The tradeoff differs at smaller inputs:

- At `1e10`, first calls are consistently **14.03–20.44% slower**, while warm
  calls are 9.75–20.39% faster across method/thread configurations.
- At `1e12`, serial first calls are 2.26–4.66% slower; eight-thread first
  results are mixed. Warm calls improve 16.46–23.06%.
- At `1e6`–`1e8`, absolute times are very small. Several `1e8` first-call
  regressions recur, while some parallel results vary substantially. The
  full matrix is retained; these cases do not establish a universal gain.

These warm gains are larger than the preceding single-call replay predicted.
That replay started the additional memo empty, whereas here later calls reuse
all entries populated in earlier calls, including states seen only once per
call. Complete-driver Phi phases can also warm entries before P2 begins.

A larger confirmation at `1e14` uses three process pairs per method/thread
configuration, one first and two warmed calls per process. Every pair improves
both first and warm calls. Individual first-call pairs improve 18.8–31.7%,
with configuration medians 20.12–24.23%; warm
configuration medians improve 24.66–27.85%.

At `1e15`, a short screen uses **one** process pair per method at eight
threads, with one first and one warm call. Seconds:

| Method | First: current / doubled | Warm: current / doubled |
| --- | ---: | ---: |
| Legendre | 6.594 / 4.926 | 6.608 / 5.147 |
| Meissel | 6.183 / 5.056 | 6.176 / 5.011 |
| Lehmer | 5.813 / 4.992 | 5.930 / 4.746 |

That screen supports a larger-input benefit but does not have the repeated
process coverage of the `1e14` confirmation. No maximum-range run was made.

### Measured cache storage

Separate memory-only copies scan allocated page pointers after each complete
driver has returned and joined its workers. At every tested input `1e12`,
`1e13`, and `1e14`, all methods and both thread budgets report:

| Cache measure | Current | Doubled |
| --- | ---: | ---: |
| Allocated pages | 533 | 1,454 |
| Allocated page payload | 4.1641 MiB | 11.3594 MiB |
| Cache object, including pointer slots/mutexes | 19,200 bytes | 32,000 bytes |
| Page payload plus cache object | 4.1824 MiB | 11.3899 MiB |
| Maximum configured page payload | 12.5 MiB | 25 MiB |

The observed increase is **7.2075 MiB** including the larger cache object.
Observed occupancy grows by about 2.72x even though configured capacity
doubles. At `1e8`, occupancy is lower: Legendre page payload grows from 1.23
to 1.42 MiB, and Meissel/Lehmer from 1.55 to 2.01 MiB. Every memory-probe
configuration has the same page count on its first and second calls.
These figures exclude tiny tables, allocator
bookkeeping and per-call prime tables; they are not process peak RSS.

The maintainer accepted the memory and first-call tradeoff and approved
production adoption after reviewing these results. No input-dependent cache
policy or public API change was added.

### The primecount observation

A separate diagnostic run of the installed **primecount 8.0** used
`primecount 1e15 --threads=8 --METHOD --status --time`. Each method ran in a
fresh process after our benchmarks. Its own rounded phase times, seconds:

| Method | Phi | P2 | P3 | Total |
| --- | ---: | ---: | ---: | ---: |
| Legendre | 1.466 | — | — | 1.466 |
| Meissel | 0.383 | 0.178 | — | 0.561 |
| Lehmer | 0.093 | 3.595 | 0.020 | 3.708 |

Here P2 accounts for about **97%** of Lehmer's reported total; P3 is small.
This reproduces the observed ordering and identifies its expensive phase.
These are single status-enabled diagnostic runs, not a new statistical
comparison between libraries; status printing and rounding affect timings.

Current upstream [P2.cpp](https://github.com/kimwalisch/primecount/blob/master/src/P2.cpp)
uses a prime-count seed per worker interval, then advances counts with
`primesieve` iterators. Our P2 evaluates Phi separately for each non-table
quotient. The installed run reports cutoff `y=100000` for Meissel and `y=5623`
for Lehmer; the latter gives an approximately 17.8x larger `x/y` endpoint.
In the upstream sieve-based implementation, that wider range is a plausible
explanation for the expensive Lehmer P2 phase. Tagged v8.0 source was not
retrievable during this study, so the source-mechanism comparison describes
current upstream rather than a verified reading of the installed binary's
exact source. Phase dominance alone does not establish our Phi-cache issue
inside primecount.

### Validation and reproduction of the doubled-cache experiment

All **840 uninstrumented timing results** and **96 memory-probe results**
matched known prime counts. The doubled candidate passed the focused boundary
checker, expanded cache checker, cache allocation-failure/recovery mode,
general Phi checker with/without caching, and native-method checker under
both ASan/UBSan and TSan. Coverage includes 2,025,842 populated cache cells,
6,272 native-method checks, old/new endpoints and page boundaries,
concurrent publication, signed results, and the independent value bound.
Builds emitted no warnings and sanitizers reported no errors. No R package
build/check was run because this experiment links only temporary native copies.

Retained files include [main timings](results/cache-double.csv),
[main summary](results/cache-double-summary.csv),
[larger-input summary](results/cache-double-large-summary.csv),
[memory observations](results/cache-double-memory.json),
[primecount phase output](results/primecount-status.json), and
[build/validation metadata](results/cache-double-build-metadata.json).
Complete sources, executables, sanitizer logs and executed validation commands
are under the directory recorded by `/tmp/rcppalgos-phi-cache-double-path`.

From the repository root, using a fresh artifact directory:

```sh
cache_double_dir=$(mktemp -d)
cp -R "$p2_baseline_source" "$cache_double_dir/source"
python3 .github/development/p2-profile/cache_double_bench.py build "$cache_double_dir"
# Complete sanitizer checks before starting any timings.
python3 .github/development/p2-profile/cache_double_bench.py bench "$cache_double_dir"
python3 .github/development/p2-profile/cache_double_bench.py large "$cache_double_dir"
python3 .github/development/p2-profile/cache_double_bench.py memory "$cache_double_dir"
```

The new focused checker is `cache_double_check.cpp`; define `PHI_SOURCE` as
the quoted path to the generated `variants/double/Phi.cpp` and link only
`Eratosthenes.cpp` alongside it. The general Phi and native-method fixtures
are unchanged; the older white-box cache fixture was adapted only in the
temporary validation directory to cover the new upper bound. Actual build
commands are retained in `validation/cache-double-commands.txt`. Main and
larger timing seeds are 202611081 and 202611082 respectively.

### Production adoption

After approval, `src/Phi.cpp` now contains the tested 131,071 endpoint and
32-page layout. Its computational code matches the benchmarked candidate;
only the explanatory base-case comment was clarified. The permanent cache
checker now scans the expanded range, rejects 131,072, and retains checks at
the former boundary alongside the new page edges.

Fresh builds against production passed the cache checker, allocation-failure
mode, focused range/concurrency checker and all 6,272 native-method checks
under both ASan/UBSan and TSan, without warnings or sanitizer diagnostics.
The cache scan verified 2,025,842 populated entries. Commands, logs and the
candidate-equivalence/preservation audit are under the directory recorded by
`/tmp/rcppalgos-phi-cache-adopt-path`. Performance was not remeasured because
the computational change is identical to the benchmarked candidate. A full
R package build/check was not rerun for this native-only adoption.

## Initial-study validation, retained data, and reproduction

All 432 uninstrumented rows matched known complete counts or the corresponding
baseline P2 result. All 96 diagnostic rows preserved those P2 results, with
consistent query, recursion, cache, and aggregation counters. Sampling loops
checked every P2 result. The expected large P2 values are obtained from the
uninstrumented complete identities, not a new independent P2 oracle; earlier
independent partial-sieve checks remain relevant. Three C++11 builds with
Apple clang 17, `-O3 -g -fno-omit-frame-pointer -pthread -Wall -Wextra
-Wpedantic` emitted no warnings. Source review checked the diagnostic lifetime
and counter semantics. No package changes, full-range run, sanitizer rerun,
or hardware-counter study was performed.

The [results](results/) directory retains the raw baseline CSV, median summary,
compact diagnostic counters/timers, expected P2 values, and stack-sample totals.
Full diagnostic output, complete stack traces, executables, immutable production
source copies, compiler commands, and initial Git state are under the directory
recorded by `/tmp/rcppalgos-p2-profile-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-p2-profile-tvi2xx3b`.

From the repository root, create a fresh temporary directory and build the
normal driver before any timed runs:

```sh
p2_profile_dir=$(mktemp -d)
clang++ -std=c++11 -O3 -g -fno-omit-frame-pointer -pthread \
  -Wall -Wextra -Wpedantic -I"$p2_baseline_source/inst/include" \
  .github/development/p2-profile/profile.cpp \
  "$p2_baseline_source/src/P2.cpp" "$p2_baseline_source/src/P3.cpp" \
  "$p2_baseline_source/src/PhiParallel.cpp" "$p2_baseline_source/src/Phi.cpp" \
  "$p2_baseline_source/src/Eratosthenes.cpp" \
  -o "$p2_profile_dir/profile"
python3 .github/development/p2-profile/run.py "$p2_profile_dir"
```

`run.py` uses seed 202611071, three fresh processes per configuration and four
calls per process; configurations and paired method order are shuffled. Cold
means the first call in a fresh process, not cold hardware caches. Warm medians
are medians of per-process warm medians. Finish other builds/tests before timing.

Create diagnostic copies and build **separate** executables:

```sh
python3 .github/development/p2-profile/instrument.py "$p2_profile_dir/instrumented" \
  --source-root "$p2_baseline_source"
clang++ -std=c++11 -O3 -g -fno-omit-frame-pointer -pthread \
  -Wall -Wextra -Wpedantic -DP2_PROFILE_METRICS \
  -I"$p2_baseline_source/inst/include" -I"$p2_profile_dir/instrumented" \
  .github/development/p2-profile/profile.cpp \
  "$p2_profile_dir/instrumented/P2.cpp" "$p2_baseline_source/src/P3.cpp" \
  "$p2_baseline_source/src/PhiParallel.cpp" \
  "$p2_profile_dir/instrumented/Phi.cpp" "$p2_baseline_source/src/Eratosthenes.cpp" \
  -o "$p2_profile_dir/profile-counts"
```

Repeat that compile command with `-DP2_PROFILE_TIMING` added and output
`"$p2_profile_dir/profile-timing"`. Then run:

```sh
python3 .github/development/p2-profile/diagnostics.py "$p2_profile_dir"
python3 .github/development/p2-profile/sample.py "$p2_profile_dir" lehmer 1
python3 .github/development/p2-profile/sample.py "$p2_profile_dir" meissel 1
python3 .github/development/p2-profile/sample.py "$p2_profile_dir" lehmer 8
python3 .github/development/p2-profile/analyze.py "$p2_profile_dir" "$p2_profile_dir/results"
```

`sample.py` is macOS-specific and samples only its own temporary child. The
sandbox blocked the first attempt in this study; sampling succeeded after
permission to run outside it. Sampling runs are separate from baseline timings.
Other systems can use their native stack profiler on the driver's `sample` mode.
The remaining harness uses portable C++11 and Python's standard library.

Instrumentation performs checked substitutions into copied `P2.cpp` and
`Phi.cpp`. Counters are thread-local; each worker publishes once into a disjoint
slot and the caller reports after joining. The harness enables metrics only
around P2, so the driver's preceding Phi phase does not pollute its counters.
