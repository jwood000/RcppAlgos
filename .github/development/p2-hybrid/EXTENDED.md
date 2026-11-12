# P2 at 64L and 128L — November 6, 2026

This follow-up extends the temporary [16L/32L experiment](WIDE.md) to `64L`
and `128L`, with fresh baseline and `32L` controls. Production remains
unchanged at commit `8b43ed790435b22eb9ec02627379082fe7deb103`. The earlier
[2L/4L/8L](README.md) and [16L/32L](WIDE.md) results remain separate.

The improvements have reached a size-dependent limit in this matrix:

- At `1e13`, warmed complete calls favor `32L`. `64L` is about 6–8% slower
  serially and 1–2% slower with eight threads.
- At `1e14`, `64L` still helps: approximately 1–2% serially and 6–7% with
  eight threads, relative to `32L`.
- At every measured input from `1e12` through `1e14`, every matched first
  and warmed complete-call `128L` comparison is slower than `64L`.

These are the best of the endpoints measured here, not a proof that the
optimal endpoint is exactly `32L` or `64L`. The input grid and thread budgets
are too coarse to choose a universal multiplier or dispatch rule.

## Experimental scope and method

`L` is the prepared table limit, normally `floor(sqrt(x))`. The endpoint is
`Q = min(multiplier * L, largest requested quotient)`, capped before forming
the product. P2 uses Phi for quotients above Q and bounded streaming for
quotients in `(L, Q]`, retaining direct lookup for quotients at most L.
The dense pi table is not enlarged.

Only the experimental multiplier allowlist/comments changed in `P2Hybrid.cpp`;
its entire namespace/kernel body is byte-for-byte unchanged. Cache paging,
atomics, dynamic 32-query chunks, stream intervals, seed calculations and
buffer policy are identical. The wider endpoint removes outer Phi queries
while increasing the streamed range and changing its seed arguments. This
experiment measures their combined cost, without separately profiling each
component of the streaming phase.

The four uninstrumented executables are **baseline / 32L / 64L / 128L**.
They use the same frozen production snapshot, common objects, Apple Clang 17,
C++11, `-O3 -pthread`, debug symbols and frame pointers on arm64 macOS 15.7.4.
Builds, sanitizer checks and separate memory probes finished before timing.

Complete calls cover both Meissel and Lehmer, inputs `1e10`, `1e12`, `1e13`,
`1e14`, and thread budgets one/eight. Isolated P2 covers the same method
cutoffs and budgets at `1e12`–`1e14`. Each configuration has three matched
blocks, with variants in randomized order and a separate fresh process per
variant. Each process performs one first and three warmed calls. Seeds are
`202611111` (complete) and `202611112` (isolated).

The complete harness follows driver order: table construction, Phi, P2,
optional P3 and table destruction. It excludes the R adapter and process
startup. Warmed complete calls rebuild tables but retain the process-wide
Phi cache. Isolated first P2 has no preceding Phi phase.

Times are medians of per-process medians. Percentage changes pair process
medians within the same block before taking the median ratio; they need not
equal ratios of the displayed times. The two percentage columns below measure
the successive steps `32L -> 64L` and `64L -> 128L`. Negative means faster.
Paired minimum/maximum changes are retained in the CSVs; these are observed
ranges, not confidence intervals. Comparisons use contemporaneous controls,
not timings from the previous run.

## Complete calls

All **768 complete counts** matched known values. P2 subtotals also agree
across variants, repetitions and worker budgets.

### Warmed calls

| Input | Method | Threads | Baseline ms | 32L ms | 64L ms | 128L ms | 64L vs 32L | 128L vs 64L |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e10` | Meissel | 1 | 1.000 | 1.332 | 1.690 | 1.702 | +28.3% | +0.7% |
| `1e10` | Meissel | 8 | 0.571 | 0.914 | 0.770 | 0.746 | -15.8% | -5.0% |
| `1e10` | Lehmer | 1 | 1.079 | 1.389 | 2.260 | 4.060 | +63.2% | +79.6% |
| `1e10` | Lehmer | 8 | 0.629 | 0.960 | 0.961 | 1.150 | -0.9% | +19.7% |
| `1e12` | Meissel | 1 | 37.227 | 29.914 | 36.052 | 45.692 | +20.8% | +26.7% |
| `1e12` | Meissel | 8 | 11.948 | 10.074 | 10.907 | 12.335 | +7.4% | +13.6% |
| `1e12` | Lehmer | 1 | 39.152 | 31.498 | 37.448 | 55.627 | +18.8% | +48.4% |
| `1e12` | Lehmer | 8 | 11.555 | 9.704 | 10.499 | 13.641 | +9.6% | +29.5% |
| `1e13` | Meissel | 1 | 276.912 | 169.609 | 180.854 | 237.270 | +6.4% | +31.0% |
| `1e13` | Meissel | 8 | 72.111 | 46.711 | 47.153 | 55.745 | +1.4% | +18.2% |
| `1e13` | Lehmer | 1 | 285.279 | 175.524 | 188.823 | 244.902 | +7.5% | +30.1% |
| `1e13` | Lehmer | 8 | 71.781 | 46.215 | 47.132 | 56.509 | +2.0% | +19.9% |
| `1e14` | Meissel | 1 | 2266.158 | 1217.252 | 1200.295 | 1348.046 | -1.6% | +12.0% |
| `1e14` | Meissel | 8 | 576.329 | 317.405 | 296.553 | 313.766 | -5.9% | +5.8% |
| `1e14` | Lehmer | 1 | 2297.922 | 1257.972 | 1238.265 | 1385.382 | -1.3% | +11.9% |
| `1e14` | Lehmer | 8 | 584.094 | 309.108 | 287.183 | 311.214 | -7.1% | +8.4% |

The `1e14` serial improvements are small but consistent across these three
blocks. The paired `64L` changes at that input range from -2.6% to -1.3%
for Meissel and -1.8% to -1.1% for Lehmer. With eight threads, median gains
are 5.9% and 7.1%; observed paired changes range from -7.9% to -5.1% and
-8.0% to -1.2%, respectively. Larger input or a different machine could
move this tradeoff.

At `1e12`, widening beyond `32L` loses ground for both methods and budgets.
At `1e13`, `64L` gives a small parallel regression and a clearer serial one;
`128L` is substantially slower. At `1e14`, moving to `128L` loses ground
relative to `64L`, with warmed regressions of about 12% serially and
6–8% with eight threads relative to `64L`.

At `1e10`, every warmed candidate is slower than production. Meissel is a
special case: `64L` already reaches the largest required quotient, so `128L`
performs identical work. Their measured difference is timing variation, not
an algorithmic benefit from increasing the multiplier.

### First calls

| Input | Method | Threads | Baseline ms | 32L ms | 64L ms | 128L ms | 64L vs 32L | 128L vs 64L |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e10` | Meissel | 1 | 4.067 | 3.236 | 3.465 | 3.474 | +7.1% | +0.2% |
| `1e10` | Meissel | 8 | 1.590 | 1.680 | 1.477 | 1.507 | -9.9% | +2.0% |
| `1e10` | Lehmer | 1 | 4.135 | 3.289 | 3.805 | 5.192 | +14.9% | +37.6% |
| `1e10` | Lehmer | 8 | 1.709 | 1.783 | 1.645 | 1.870 | -4.8% | +16.4% |
| `1e12` | Meissel | 1 | 85.109 | 54.874 | 59.940 | 69.181 | +9.9% | +15.3% |
| `1e12` | Meissel | 8 | 22.968 | 16.641 | 17.349 | 18.556 | +2.3% | +10.3% |
| `1e12` | Lehmer | 1 | 88.395 | 56.104 | 61.930 | 78.587 | +10.4% | +26.9% |
| `1e12` | Lehmer | 8 | 22.120 | 15.613 | 16.193 | 19.563 | +4.8% | +19.5% |
| `1e13` | Meissel | 1 | 420.758 | 257.659 | 271.207 | 323.817 | +5.7% | +19.5% |
| `1e13` | Meissel | 8 | 104.709 | 70.755 | 69.782 | 78.655 | -0.5% | +12.7% |
| `1e13` | Lehmer | 1 | 422.201 | 263.796 | 276.025 | 330.591 | +4.6% | +19.8% |
| `1e13` | Lehmer | 8 | 104.492 | 68.403 | 68.709 | 75.204 | +1.9% | +9.9% |
| `1e14` | Meissel | 1 | 2460.285 | 1414.196 | 1390.050 | 1538.919 | -1.2% | +10.4% |
| `1e14` | Meissel | 8 | 606.217 | 366.137 | 341.185 | 361.720 | -3.5% | +5.9% |
| `1e14` | Lehmer | 1 | 2504.303 | 1455.932 | 1438.374 | 1572.421 | -1.8% | +9.6% |
| `1e14` | Lehmer | 8 | 630.999 | 351.649 | 329.356 | 352.360 | -6.3% | +6.8% |

At `1e14`, first-call medians also favor `64L`, but the eight-thread Lehmer
comparison has one +0.7% pair alongside two improvements. At `1e13`, the
parallel `32L`/`64L` first-call comparisons are close and mixed. All first-call
`128L` comparisons at `1e12`–`1e14` regress against `64L`. The smaller cases
show why first and warmed behavior should both inform a future policy.

## Isolated P2

All **576 isolated P2 results** matched the checked reference subtotals.
Warmed P2 times and successive paired changes are:

| Input | Cutoff | Threads | Baseline ms | 32L ms | 64L ms | 128L ms | 64L vs 32L | 128L vs 64L |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e12` | Meissel | 1 | 23.876 | 16.600 | 22.752 | 32.353 | +36.9% | +42.3% |
| `1e12` | Meissel | 8 | 5.443 | 3.544 | 4.342 | 5.960 | +21.7% | +36.7% |
| `1e12` | Lehmer | 1 | 33.085 | 25.306 | 31.476 | 49.626 | +24.4% | +57.2% |
| `1e12` | Lehmer | 8 | 7.782 | 5.896 | 6.593 | 9.634 | +13.1% | +47.0% |
| `1e13` | Meissel | 1 | 190.230 | 82.311 | 94.094 | 151.168 | +13.8% | +60.7% |
| `1e13` | Meissel | 8 | 43.806 | 18.636 | 18.647 | 28.235 | +0.9% | +49.3% |
| `1e13` | Lehmer | 1 | 246.640 | 140.317 | 152.173 | 209.492 | +8.1% | +37.0% |
| `1e13` | Lehmer | 8 | 57.176 | 31.180 | 31.566 | 40.231 | +1.2% | +27.4% |
| `1e14` | Meissel | 1 | 1518.962 | 489.180 | 468.002 | 608.921 | -4.9% | +30.1% |
| `1e14` | Meissel | 8 | 378.883 | 108.137 | 93.218 | 113.399 | -12.6% | +23.1% |
| `1e14` | Lehmer | 1 | 2024.340 | 1031.126 | 1048.430 | 1128.528 | +0.0% | +7.9% |
| `1e14` | Lehmer | 8 | 498.112 | 229.340 | 213.950 | 232.738 | -5.9% | +8.9% |

At `1e14`, warmed parallel P2 improves about 12.6%/5.9% from `32L` to `64L`
for Meissel/Lehmer. Serial Meissel improves about 4.9%, but serial Lehmer
shows no clear isolated gain: its paired changes span -1.9% to +5.1%, with
a median near zero. The complete serial gain at that input is also small;
these experiments do not support treating it as a large, universal benefit.
At `1e13`, the warmed parallel `32L`/`64L` difference is small, while serial
calls regress. Every matched warmed `128L` comparison is slower than `64L`.

Isolated first P2 starts without the driver's preceding Phi phase and has
more mixed comparisons than complete first calls. In particular, one
Lehmer `1e14` eight-thread `64L` first call takes 325.442 ms, versus
247.965/252.374 ms in the other blocks. Its paired `128L` call takes
270.732 ms, giving a -16.8% change; the other two first-call pairs regress
7.5% and 10.3%, and every warmed pair regresses. One Meissel `1e12`
eight-thread first-call pair also favors `128L` by 0.9%. These observations
remain in the raw data and summaries, without assigning a cause.

One warmed Meissel `1e14` eight-thread `32L` process is slower than its peers
(127.447 ms versus 105.383/108.137 ms). Its `64L` pair improves 26.9%; the
other two improve 12.2% and 12.6%. The retained median paired gain is 12.6%,
and all observations remain included. No separate confirmation replaces the
primary matrix.

Use `p2_ms`: isolated `total_ms` includes table setup only on repetition zero
and destruction only on the last repetition. First-call summaries and all
paired ranges are retained in the CSVs.

## Work and memory

All **96 independently instrumented memory/work rows** passed complete-count
and query/worker checks. Every `32L` row exactly matches the earlier study;
first/warmed work and capacity fields are identical.

At `1e13`, the outer Phi-prefix queries continue to decline:

| Cutoff | 32L | 64L | 128L |
| --- | ---: | ---: | ---: |
| Meissel | 7,068 | 2,659 | 318 |
| Lehmer | 9,210 | 4,801 | 2,460 |

Both cutoffs stream 218,162 / 222,571 / 224,912 queries respectively. These
outer counts exclude the stream's Phi seeds and do not measure recursive
work. The runtime regressions show that saving these remaining prefix
queries is insufficient to pay for the expanded streaming work in this case.

At `1e10`, Meissel's `64L` and `128L` both cap Q at 4,627,487, stream all
9,267 queries, launch no Phi-prefix workers and use four stream workers with
an eight-thread budget. Their work and memory measurements match exactly.
At `1e12`, Meissel's `128L` cap is 99,930,048: it streams all 77,269 queries,
whereas `64L` retains 592 Phi-prefix queries. A larger multiplier cannot
extend the stream beyond the largest query actually needed.

| Variant | Serial scratch range, bytes | Maximum across matrix, bytes | Maximum, KiB |
| --- | ---: | ---: | ---: |
| 32L | 34,984–49,168 | 361,640 | 353.2 |
| 64L | 35,368–55,080 | 396,344 | 387.1 |
| 128L | 35,368–66,190 | 447,607 | 437.1 |

The maximum tracked payload grows by **85,967 bytes (about 84 KiB)** from
`32L` to `128L`. The existing buffer rule uses
`max(32768, floor(sqrt(finish)))` bytes per stream, so `128L` at `1e14`
exercises larger buffers. Its aggregate sieve capacity is 265,863 bytes;
next-crossing positions account for 181,680 bytes, plus subtotal storage.

These are conservative vector-capacity bounds, not peak RSS. They take the
maximum of prefix subtotal storage and streaming subtotals plus all worker
sieve/next-crossing capacities. They exclude existing tables, the Phi cache,
ThreadBatch storage, stacks, allocator overhead and diagnostic metadata;
worker capacities need not all be live simultaneously.

## Correctness and numerical boundaries

Baseline, `32L`, `64L` and `128L` each passed **60,322 checks under both
ASan/UBSan and ThreadSanitizer**, with no compiler warnings or sanitizer
reports. This adds 108 independent split-boundary cases for `64L`/`128L`,
retaining all earlier cases. The prepared limits satisfy `x <= L^2`, and
all required oracle quotients fit the independent dense table through 9e6.

A separate optional private-stream fixture passed **88 checks per sanitizer**.
It tests short windows around the transition from a 32,768-byte buffer to a
larger buffer, across `INT32_MAX`, and near the maximum canonical `128L`
endpoint **12,148,001,920**. It crosses fill boundaries and repeats endpoints,
using independent trial division for interval counts. A synthetic initial
count of `INT32_MAX` checks 64-bit accumulation; this is an increment test,
not a claim that the seed equals the actual prime count there. Test-only
inclusion exposes the private class without adding hooks to the kernel.

The integer bounds remain unchanged. For a nonempty sum, `p[first] >= 2` and
`finish <= Q <= floor(x / p[first])`, hence
`floor(sqrt(finish)) <= floor(sqrt(x)) <= L`. Every required base prime and
pi-table index is therefore covered. Multiplication is performed only when
`L <= floor(largestQuotient / multiplier)`, which proves the product fits.
Also `Q + 1 <= 2^52` in the supported domain. Coordinates and crossing
positions use `int64_t`; the largest canonical buffer index is bounded by
110,217. No 32-bit coordinate narrowing is introduced.

The performance matrix covers one machine, the standard Meissel/Lehmer
cutoffs and selected inputs through `1e14`; it does not establish behavior
for arbitrary public `semiPrimeCount(n, a)` workloads, the full range, or
other machines. No full R package check or new exception-failure injection
was run for these standalone copies. Production source, public interfaces,
cache settings and dispatch are unchanged.

## Reproduction and retained artifacts

Run from the repository root with a fresh artifact directory. Finish all
builds, checks and probes before timing:

```sh
p2_extended_dir=$(mktemp -d)
mkdir "$p2_extended_dir/source"
git archive 8b43ed790435b22eb9ec02627379082fe7deb103 src inst/include | \
  tar -x -C "$p2_extended_dir/source"
git rev-parse 8b43ed790435b22eb9ec02627379082fe7deb103 > "$p2_extended_dir/commit.txt"
python3 .github/development/p2-hybrid/bench.py build "$p2_extended_dir" --variants baseline hybrid32 hybrid64 hybrid128
python3 .github/development/p2-hybrid/validate.py "$p2_extended_dir" --variants baseline 32 64 128 --stream-window
python3 .github/development/p2-hybrid/build_memory.py build "$p2_extended_dir" --multipliers 32 64 128
python3 .github/development/p2-hybrid/build_memory.py run "$p2_extended_dir" --multipliers 32 64 128
python3 .github/development/p2-hybrid/bench.py phases "$p2_extended_dir" --variants baseline hybrid32 hybrid64 hybrid128 --seed 202611111
python3 .github/development/p2-hybrid/bench.py p2 "$p2_extended_dir" --variants baseline hybrid32 hybrid64 hybrid128 --seed 202611112
```

The [extended results directory](results-extended/) retains raw timings,
paired summaries, memory/work rows, validation totals, exact build commands
and source/binary hashes. The preservation audit verifies 392 production
files unchanged, along with all 24 original result files and all 18 files
from the 16L/32L study. Nothing is staged or committed.
Runner defaults still select the original variant sets and seeds. The current
checker includes the later boundary extensions; the extra private-stream
fixture is enabled explicitly with `--stream-window`.

The full artifact directory is recorded by `/tmp/rcppalgos-p2-hybrid-extended-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-p2-hybrid-extended-a4qm1p4b`.
It contains the frozen sources, pre-extension experiment snapshot, executables,
build/sanitizer logs and raw per-process timing outputs. Nothing is staged
or committed.
