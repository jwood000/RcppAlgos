# Wider P2 streaming endpoints — November 6, 2026

This follow-up adds temporary `16L` and `32L` variants to the
[original bounded-streaming experiment](README.md). Here `L = tables.limit()`,
normally `floor(sqrt(x))`, and the multiplier bounds the quotient range streamed
by P2. It does not enlarge the dense pi table. Production code and dispatch
remain unchanged at commit `8b43ed790435b22eb9ec02627379082fe7deb103`.

The subsequent [64L/128L experiment](EXTENDED.md) extends this study. The
results below retain the original 16L/32L measurements and conclusions.

The gains continue at the larger measured inputs. At `1e13`–`1e14`, warmed
complete calls improve approximately **18–22% beyond 8L** with `32L`; `16L`
improves approximately 10–14%. First calls also benefit from `32L`. The result
is size-dependent: `1e12` warmed serial calls favor `16L`, and both wider
endpoints make warmed `1e10` calls slower than `8L`. These data do not support
a universal multiplier or establish an optimal dispatch policy.

## Scope and methodology

Only the experimental multiplier guard changed in `P2Hybrid.cpp`; its kernel
body is identical to the original experiment. Paging, atomics, Phi evaluation,
32-query scheduling chunks, and streaming worker/buffer policies are unchanged.
The endpoint remains `Q = min(multiplier * L, largest requested quotient)`,
using the existing overflow-safe cap. The supported new multipliers change
how much work goes to each path, without changing the counting identity.

Four fresh executables compare **baseline / 8L / 16L / 32L**, using the same
frozen production snapshot and common objects. This run remeasures baseline
and `8L`; it does not compare new observations directly with the earlier run.
Builds use Apple Clang 17, C++11 and `-O3 -pthread` on arm64 macOS 15.7.4.
All builds, sanitizer checks and capacity probes finished before timing.

Complete calls cover Meissel/Lehmer, inputs `1e10`, `1e12`, `1e13`, `1e14`, and
one/eight requested threads. Isolated P2 covers the same cutoffs and budgets
at `1e12`–`1e14`. Each configuration has three randomized matched blocks;
each variant runs in a fresh process with one first and three warmed calls.
The complete harness follows driver order (tables, Phi, P2, optional P3,
destruction), excluding the R adapter and process startup. Complete warmed
calls rebuild tables and reuse the process-wide Phi cache. Isolated first P2
has no preceding Phi phase.

Tables report medians of per-process medians. Percentage changes first pair
candidate and `8L` within each block, then take the median ratio. Consequently,
percentages need not equal ratios of the separately displayed times. The CSVs
also retain paired minimum/maximum changes; these ranges are observations,
not confidence intervals. Seeds are `202611101` (complete) and `202611102`
(isolated). The earlier study and its raw data remain separate and unchanged.

## Complete-call results

All **768 complete counts** matched known values, and P2 subtotals agreed
across variants, repetitions and worker budgets. Negative percentage changes
mean faster than the contemporaneous `8L` variant.

### Warmed calls

| Input | Method | Threads | Baseline ms | 8L ms | 16L ms | 32L ms | 16L vs 8L | 32L vs 8L |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e10` | Meissel | 1 | 1.009 | 0.816 | 0.946 | 1.322 | +15.8% | +61.5% |
| `1e10` | Meissel | 8 | 0.566 | 0.754 | 0.952 | 0.981 | +20.8% | +29.6% |
| `1e10` | Lehmer | 1 | 1.119 | 0.913 | 1.055 | 1.450 | +14.2% | +57.1% |
| `1e10` | Lehmer | 8 | 0.677 | 0.783 | 1.042 | 0.994 | +33.1% | +21.9% |
| `1e12` | Meissel | 1 | 37.410 | 29.649 | 28.801 | 29.964 | -3.0% | +1.0% |
| `1e12` | Meissel | 8 | 12.182 | 10.673 | 10.262 | 10.126 | -3.4% | -5.0% |
| `1e12` | Lehmer | 1 | 38.930 | 31.348 | 30.403 | 31.390 | -2.2% | -0.2% |
| `1e12` | Lehmer | 8 | 11.745 | 10.292 | 9.842 | 9.764 | -4.4% | -5.0% |
| `1e13` | Meissel | 1 | 277.810 | 208.887 | 183.755 | 169.459 | -12.0% | -18.9% |
| `1e13` | Meissel | 8 | 73.958 | 58.322 | 51.903 | 47.143 | -11.0% | -19.2% |
| `1e13` | Lehmer | 1 | 285.393 | 216.373 | 190.610 | 175.616 | -12.2% | -18.4% |
| `1e13` | Lehmer | 8 | 72.677 | 57.311 | 51.805 | 46.967 | -10.4% | -18.0% |
| `1e14` | Meissel | 1 | 2258.826 | 1531.366 | 1347.166 | 1232.159 | -11.8% | -19.7% |
| `1e14` | Meissel | 8 | 573.260 | 407.463 | 348.986 | 324.324 | -14.4% | -20.1% |
| `1e14` | Lehmer | 1 | 2298.672 | 1565.636 | 1397.170 | 1260.485 | -11.4% | -19.5% |
| `1e14` | Lehmer | 8 | 566.661 | 381.979 | 327.943 | 298.125 | -13.2% | -22.0% |


At `1e12`, `16L` improves warmed serial calls roughly 2–3% over `8L`.
Moving to `32L` gives no consistent additional serial benefit: Meissel is
about 1% slower than `8L`, and Lehmer is approximately tied. Eight-thread
`32L` calls improve about 5% over `8L`. At `1e13` and `1e14`, both wider
endpoints improve every matched warmed comparison with `8L`, and `32L`
has the lowest median in all eight configurations.

At `1e10`, widened endpoints cost more than `8L` in every matched warmed
comparison. The absolute cost is small, but meaningful for these short calls;
`32L` also exceeds the production baseline in all four warmed configurations.

### First calls

| Input | Method | Threads | Baseline ms | 8L ms | 16L ms | 32L ms | 16L vs 8L | 32L vs 8L |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e10` | Meissel | 1 | 4.032 | 3.177 | 3.089 | 3.267 | -2.3% | +2.8% |
| `1e10` | Meissel | 8 | 1.691 | 1.794 | 1.770 | 1.968 | +6.0% | +1.3% |
| `1e10` | Lehmer | 1 | 4.309 | 3.318 | 3.229 | 3.500 | -3.5% | +1.8% |
| `1e10` | Lehmer | 8 | 2.170 | 1.796 | 2.058 | 1.860 | +8.5% | +3.6% |
| `1e12` | Meissel | 1 | 84.889 | 62.790 | 56.149 | 54.622 | -10.6% | -13.0% |
| `1e12` | Meissel | 8 | 23.374 | 18.750 | 17.741 | 16.985 | -5.4% | -11.3% |
| `1e12` | Lehmer | 1 | 86.356 | 64.724 | 57.916 | 56.116 | -9.8% | -13.4% |
| `1e12` | Lehmer | 8 | 22.404 | 17.858 | 16.467 | 15.939 | -7.3% | -10.4% |
| `1e13` | Meissel | 1 | 433.271 | 303.828 | 272.615 | 257.179 | -10.7% | -15.4% |
| `1e13` | Meissel | 8 | 108.003 | 83.107 | 74.595 | 72.825 | -9.8% | -12.4% |
| `1e13` | Lehmer | 1 | 422.395 | 310.356 | 278.986 | 264.288 | -10.2% | -14.8% |
| `1e13` | Lehmer | 8 | 107.053 | 79.429 | 73.768 | 68.781 | -7.1% | -13.9% |
| `1e14` | Meissel | 1 | 2467.467 | 1715.838 | 1526.742 | 1415.993 | -11.2% | -17.6% |
| `1e14` | Meissel | 8 | 608.616 | 451.903 | 393.730 | 366.250 | -12.9% | -18.2% |
| `1e14` | Lehmer | 1 | 2507.930 | 1763.699 | 1654.069 | 1452.075 | -5.5% | -17.7% |
| `1e14` | Lehmer | 8 | 608.485 | 426.356 | 387.458 | 340.286 | -10.2% | -18.6% |

Every matched first-call `32L` comparison at `1e12`–`1e14` improves on `8L`.
One first-call serial Lehmer `1e14` comparison with `16L` regresses 4.7%;
the other two improve, and all three `32L` comparisons improve about 17–18%.
That observation is retained without attributing a cause. First-call behavior
at `1e10` is mixed, especially in the short parallel runs.

## Isolated P2 results

All **576 isolated P2 results** matched the checked reference values. The
table reports warmed P2 milliseconds, with paired changes against `8L`.

| Input | Cutoff | Threads | Baseline ms | 8L ms | 16L ms | 32L ms | 16L vs 8L | 32L vs 8L |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `1e12` | Meissel | 1 | 23.655 | 15.893 | 15.182 | 16.634 | -4.6% | +4.7% |
| `1e12` | Meissel | 8 | 5.890 | 4.166 | 3.710 | 3.635 | -10.9% | -11.7% |
| `1e12` | Lehmer | 1 | 32.853 | 25.143 | 24.058 | 25.324 | -4.5% | +0.9% |
| `1e12` | Lehmer | 8 | 7.588 | 6.085 | 5.928 | 5.698 | -5.8% | -7.6% |
| `1e13` | Meissel | 1 | 189.639 | 121.707 | 95.294 | 81.893 | -21.4% | -32.1% |
| `1e13` | Meissel | 8 | 42.521 | 28.584 | 21.640 | 18.003 | -22.0% | -38.0% |
| `1e13` | Lehmer | 1 | 247.376 | 179.208 | 152.484 | 138.988 | -15.0% | -22.1% |
| `1e13` | Lehmer | 8 | 54.744 | 40.428 | 35.893 | 29.864 | -15.0% | -26.1% |
| `1e14` | Meissel | 1 | 1519.646 | 795.460 | 598.698 | 488.657 | -24.8% | -38.8% |
| `1e14` | Meissel | 8 | 362.978 | 186.838 | 133.284 | 107.296 | -28.7% | -40.4% |
| `1e14` | Lehmer | 1 | 2032.282 | 1311.286 | 1098.612 | 999.602 | -15.3% | -24.0% |
| `1e14` | Lehmer | 8 | 496.899 | 312.537 | 260.265 | 225.693 | -16.7% | -27.8% |

At `1e13`–`1e14`, the paired median warmed P2 reductions from `8L` to `32L`
range from 22% to 40%, larger than the complete-call reductions because the
other driver phases remain. At `1e12`, warmed serial P2 favors `16L`, while
`32L` slightly regresses against `8L`; the parallel medians favor `32L`.
All matched first-call `32L` comparisons improve on `8L`; their median
reductions range from approximately 13% to 38% across this matrix.

One isolated serial Meissel `1e13` `8L` process is slower than its peers:
first/warm P2 are 244.008/203.575 ms, versus approximately 156/120–122 ms in
the other two processes. Its paired `32L` warm change is -60.3%, whereas the
other two are -32.1% and -31.6%. All observations remain in the raw data and
summaries, with no cause assigned and no samples discarded. The median paired
gain is -32.1%; the slow observation does not establish a 60% typical gain.

Use `p2_ms`, because isolated `total_ms` includes table setup only on the
first repetition and destruction only on the last repetition.

## Work and memory

The independently compiled memory probes checked **96 complete results** and
query/worker invariants. First and warmed calls have identical tracked buffer
capacities. At `1e13`, widening the stream removes more outer Phi queries:

| Method | Outer Phi queries at 8L | At 16L | At 32L |
| --- | ---: | ---: | ---: |
| Meissel | 31,082 | 15,373 | 7,068 |
| Lehmer | 33,224 | 17,515 | 9,210 |

Streaming answers 194,148 / 209,857 / 218,162 queries respectively for either
cutoff. These outer-query counts exclude the parallel stream's Phi seed
queries and say nothing by themselves about recursive work or elapsed time.
The wider interval costs more sieve work; measured timings determine whether
that cost is worthwhile.

| Variant | Serial tracked scratch range, bytes | Maximum across matrix, bytes | Maximum, KiB |
| --- | ---: | ---: | ---: |
| 8L | 34,008–41,672 | 318,568 | 311.1 |
| 16L | 34,416–44,864 | 336,536 | 328.6 |
| 32L | 34,984–49,168 | 361,640 | 353.2 |

The largest tracked payload rises by **43,072 bytes (about 42 KiB)** between
`8L` and `32L`. It remains bounded streaming storage, rather than a pi array
through `32L`. These are conservative vector-capacity bounds: the maximum of
prefix subtotals and streaming subtotals plus all stream sieve/next-crossing
buffers. Worker allocations need not all overlap. They exclude existing
tables, the Phi cache, ThreadBatch storage, stacks, allocator overhead and
diagnostic metadata; they are **not peak RSS measurements**. All `8L` work
and capacity rows exactly match the earlier study.

With eight requested threads, actual streaming workers at `1e10` are 1/1/2
for `8L`/`16L`/`32L`; at `1e12` they are 6/8/8; at `1e13` and `1e14` all use
eight. This follows the unchanged minimum-span policy, so the timings include
both the wider interval and any resulting worker-budget change.

## Validation and limits

Baseline, `8L`, `16L` and `32L` each passed **60,214 checks under both
ASan/UBSan and ThreadSanitizer**, without compiler warnings or sanitizer
reports. The checker retains the earlier independent factorization/pi oracle,
concurrency and differential checks, and adds 108 cases around the `16L` and
`32L` split boundaries, including both ends of floor-division plateaus.
The prepared limits and independent oracle cover every added case.

The existing arithmetic remains bounded: a nonempty sum has
`Q <= x / p[first] <= x / 2`, so `Q + 1` is representable in the supported
range, and sieve primes satisfy `floor(sqrt(Q)) <= floor(sqrt(x)) <= L`. Multiplication
occurs only after the cap proves the product does not exceed the largest
requested quotient. Tests provide additional evidence, not a replacement
for these bounds.

This matrix covers one platform, selected inputs through `1e14`, and the
standard Meissel/Lehmer cutoffs. It does not establish performance for arbitrary
`semiPrimeCount(n, a)`, the full supported range, other machines, or an optimal
crossover. Full R package checks and new exception-failure injection were not
run: these are standalone temporary kernels, with existing ThreadBatch
failure coverage unchanged. No production algorithm, R interface or default
has been changed.

The useful next decision is which size/thread regions deserve further
crossover measurements. `32L` is the strongest endpoint measured here at
`1e13`–`1e14`; the smaller cases show why that is not yet a universal choice.

## Reproduction and retained artifacts

Use a fresh directory. Finish all builds and checks before starting timing:

```sh
p2_wide_dir=$(mktemp -d)
mkdir "$p2_wide_dir/source"
git archive 8b43ed790435b22eb9ec02627379082fe7deb103 src inst/include | \
  tar -x -C "$p2_wide_dir/source"
git rev-parse 8b43ed790435b22eb9ec02627379082fe7deb103 > "$p2_wide_dir/commit.txt"
python3 .github/development/p2-hybrid/bench.py build "$p2_wide_dir" --variants baseline hybrid8 hybrid16 hybrid32
python3 .github/development/p2-hybrid/validate.py "$p2_wide_dir" --variants baseline 8 16 32
python3 .github/development/p2-hybrid/build_memory.py build "$p2_wide_dir" --multipliers 8 16 32
python3 .github/development/p2-hybrid/build_memory.py run "$p2_wide_dir" --multipliers 8 16 32
python3 .github/development/p2-hybrid/bench.py phases "$p2_wide_dir" --variants baseline hybrid8 hybrid16 hybrid32 --seed 202611101
python3 .github/development/p2-hybrid/bench.py p2 "$p2_wide_dir" --variants baseline hybrid8 hybrid16 hybrid32 --seed 202611102
```

The [wide results directory](results-wide/) keeps this follow-up separate from
[the original results](results/). It retains raw timings, summaries, memory
observations, validation totals, exact build commands and source/binary hashes.
The preservation audit verifies all 392 production files against their
pre-experiment hashes and all 24 original result files against the saved
pre-extension copy. No changes were staged or committed.
The original default variant sets and seeds still reproduce the earlier
experiment's design. This study's checker added 108 boundary cases; subsequent
extensions retain those checks and are documented in their own report.

The full artifact directory, recorded by `/tmp/rcppalgos-p2-hybrid-wide-path`, is
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-p2-hybrid-wide-2qzxd8mn`.
It includes immutable production sources, the pre-extension experiment copy,
executables, build/sanitizer logs and all raw timing process outputs.
