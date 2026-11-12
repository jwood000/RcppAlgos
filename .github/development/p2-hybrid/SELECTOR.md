# Frozen P2 selector v1 — November 6, 2026

This study adds a temporary runtime selector to the fixed streaming variants.
Production code, public interfaces, cache settings and method dispatch remain
unchanged at commit `8b43ed790435b22eb9ec02627379082fe7deb103`.
The policy was frozen before collecting the new timings; it is not retuned
using this evaluation matrix. Its exact header hash is recorded in `suite.json`.

## Frozen policy

The policy uses `x`, the prepared limit `L`, the actual non-table prime-index
range determined by `a`, and the original Phi worker budget. It runs once
inside P2, outside the query loops. No clocks, cache scans or allocations are
used to choose the boundary.

Apply the following rules in order:

1. Fewer than 4,096 non-table queries: use Phi only.
2. `x < 1e11` with more than one usable Phi worker: use Phi only.
3. Otherwise choose the base multiplier from this table.
4. Promote `16L` to `32L` when `x >= 1e12` and the original Phi worker budget
   is at least four.
5. Cap the endpoint at the largest requested quotient, using the existing
   overflow-safe cap. If fewer than 4,096 queries would actually stream,
   fall back to Phi only.

| Input range | Base multiplier |
| --- | ---: |
| `x < 1e11` | 8L |
| `1e11 <= x < 3e12` | 16L |
| `3e12 <= x < 3e13` | 32L |
| `x >= 3e13` | 64L |

The Phi budget is computed from the entire original non-table range using
its existing 2,048-term grain and the normalized requested thread count.
These constants are empirical trial choices, not mathematical boundaries.
The query guards make arbitrary cutoffs and larger prepared tables affect
selection. `a` affects both the available query count and largest quotient.

A Phi-only choice has `Q=L` and an empty streaming range; it uses the existing
Phi prefix and launches no streaming workers. Empty sums and fully table-
covered sums return before selection. A positive multiplier is a cap:
`Q=min(k*L, floor(x/p[a+1]))` can be smaller than `k*L`.

The exact identity `floor(x/p) <= Q` iff `p > floor(x/(Q+1))` partitions the
same P2 sum into disjoint exhaustive ranges. Changing Q changes which exact
prime-counting path answers each query, not the integer result. The counting
loops, scheduling, paging and atomics are unchanged. A preprocessed fixed32
build matches the preceding experimental source after whitespace normalization.

## Evaluation design

Six fresh executables compare baseline Phi-only, fixed `8L`, `16L`, `32L`,
`64L`, and selector v1. All use a frozen production snapshot and common
objects, Apple Clang 17, C++11, `-O3 -pthread`, debug symbols and frame pointers
on arm64 macOS 15.7.4. Instrumented memory builds are separate from timing.
All compilation, references, sanitizer checks and probes finish before timing.

There are 80 complete-call configurations and 122 isolated P2 configurations.
Each has three randomized matched blocks, six fresh variant processes per
block, and one first plus three warmed calls per process. The complete seed
is `202611121`; the isolated seed is `202611122`. Exact cases and frozen policy
hash are retained in [suite.json](results-selector/suite.json).

Complete calls cover both Meissel and Lehmer at `3e10`, `1e11`, `3e11`,
`1e12`, `3e12`, `1e13`, `3e13`, `5e13`, `1e14`, with 1/2/4/8 requested
threads. Additional `2e12` and `2e13` cases use one/four threads.

Isolated cases cover inputs `1e8`, `2e9`, `1e10`, `3e12`, `3e13`, selected
budgets 1/2/4/8, and larger `2*floor(sqrt(x))` tables at `1e10`. The cutoff
labels include zero, one, the standard Meissel/Lehmer choices, half/double
the Meissel **index**, half the square-root prime count, and sums retaining
up to 32, 4,095, 4,096, one or zero terms. Definitions and actual numeric `a`
values appear in every row. This is a targeted matrix, not a full Cartesian
product of every input, cutoff and budget.

The selector trace is computed after timed P2, without evaluating Phi or
warming data before the first isolated call. Actual selection remains inside
P2 timing. Complete total time is the sum of setup, Phi, P2, optional P3 and
teardown clocks, excluding diagnostic collection and the R adapter. This
phase-sum definition is shared by all six variants; use contemporaneous
controls rather than comparing absolute times directly with earlier studies.
Isolated `p2_ms` excludes table setup/destruction; its first call has no
preceding Phi phase. Complete warmed calls rebuild tables and retain Phi cache.

Summaries first take a median within each process, then across processes.
Ratios and differences pair process medians within matched blocks. They compare
against baseline, the fixed candidate selected by the policy, and the fastest
fixed median for that case/state. The last is an optimistic hindsight
comparison among Phi-only and fixed8/16/32/64, not an implementable oracle.
Matching-choice comparisons measure the entire compiled variant, including
possible code-generation and timing variation; they do not isolate the few
instructions in the selector itself.

Every sample is retained. Before isolated timing, the analysis convention was
set to emphasize absolute differences for cases whose fastest fixed median
is below 0.01 ms, where percentage changes are unstable. Relative summaries
for larger cases use 5%, 10% and 20% bands. Candidates with identical capped
partitions can have different observed timings without different work.

## Evaluation results

All **14,544 timed results** match their expected values: 5,760 complete-call
results and 8,784 isolated P2 results. P2 subtotals also agree across methods,
variants, budgets and repeats wherever `(x, a, table limit)` agrees. The
policy remains exactly as frozen above.

The selector captures most of the complete-call benefit, but its two fallback
rules are too coarse for arbitrary cutoffs. This is evidence for retaining a
dynamic choice, not for promoting v1 unchanged.

### Complete Meissel and Lehmer calls

These are counts of the 80 configurations, using paired comparisons against
the fastest fixed candidate for each configuration and cache state:

| Comparison with fastest fixed candidate | First | Warmed |
| --- | ---: | ---: |
| Within 5% (or faster) | 64 / 80 | 74 / 80 |
| Within 10% (or faster) | 71 / 80 | 77 / 80 |
| Within 20% (or faster) | 78 / 80 | 80 / 80 |
| Median paired time change | +0.72% | +0.51% |

Against production Phi-backed P2, median paired time changes across the cases
are -34.15% first and -26.78% warmed; 76 of 80 configurations improve in each
state. These are unweighted distribution summaries over this matrix, not an
overall speedup for an assumed user workload. The following table shows the
median paired change against production across the methods/thread budgets
tested at each input. Negative values mean less elapsed time.

| Input | Configurations | First change | Warmed change |
| --- | ---: | ---: | ---: |
| `3e10` | 8 | +0.0% | -0.9% |
| `1e11` | 8 | -19.7% | -5.6% |
| `3e11` | 8 | -26.8% | -13.2% |
| `1e12` | 8 | -30.6% | -18.7% |
| `2e12` | 4 | -29.6% | -22.0% |
| `3e12` | 8 | -32.8% | -26.8% |
| `1e13` | 8 | -36.5% | -37.0% |
| `2e13` | 4 | -38.4% | -41.2% |
| `3e13` | 8 | -38.3% | -43.0% |
| `5e13` | 8 | -40.6% | -44.9% |
| `1e14` | 8 | -43.5% | -47.5% |

The clearest input-threshold miss is `3e10` with two threads. The selector
retains Phi, but fixed16 improves first calls in every matched block. The
selector takes 30.27% longer for Meissel and 30.63% longer for Lehmer than
fixed16, with paired differences of about 1.32 and 1.22 ms. Warmed gaps are
smaller, 9.60% and 6.25%. The blanket parallel fallback below `1e11` therefore
extends too far for these complete-call cases.

The median selector difference from its **own chosen fixed candidate** is
-0.013% first and +0.033% warmed. However, individual same-choice comparisons
vary materially: warmed Lehmer at `1e11`/two threads is 19.77% slower than
fixed16, and at `5e13`/eight threads is 12.51% slower than fixed64. In both
cases the selector chooses the best fixed label. Setup and Phi timings also
vary, so these are not evidence that changing the multiplier would fix the
gap, nor a measurement of selection-instruction overhead. Three process
blocks on one host do not fully separate system variability and code-generation
effects. No samples have been removed.

### Arbitrary cutoffs and table coverage

For cases whose fastest fixed median is at least 0.01 ms, the selector is
within 5% of that candidate in 35/90 first-call configurations and 65/86 warmed
configurations. The median paired gaps are +9.74% first and +0.11% warmed, but
these medians hide severe misses at the query guard. The remaining 32 first
and 36 warmed configurations are retained and assessed using absolute time.

For serial `x=3e13`, `L=5,477,225`, two adjacent cutoffs expose an artificial
performance cliff. Times below are medians of process medians in milliseconds:

| Cutoff `a` | Non-table queries | Selector choice | Selector first | Selector warm | Fixed64 first | Fixed64 warm |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 375238 | 4095 | Phi only | 1.989 | 1.078 | 0.077 | 0.031 |
| 375237 | 4096 | Stream, capped 64L | 0.072 | 0.034 | 0.070 | 0.032 |

The 4,095-query stream ends at 5,541,380: only **64,155 integers beyond L**.
Every tested fixed multiplier caps at that same endpoint and streams all
queries with one worker. Thus the lesson is about permitting this short
stream, not preferring 64 over another multiplier. Production Phi takes
1.078 ms warmed in the 4,095-query case, essentially the selector's time.
The large missed saving is not a new regression against production.

The same cliff occurs at `3e12`: serial 4,095-term warmed calls take 0.512 ms
with the selector versus about 0.030 ms with a fixed stream. Larger prepared
tables can trigger it too. At `1e10`, table multiplier two and the half-root
cutoff leave just 337 non-table queries: first calls take 0.082 ms with the
selector versus 0.012 ms with streaming; warmed differences are about 15 us.
Even a 32-term near-root tail at `3e13` saves roughly 70–77 us on first calls
when streamed, although its fastest fixed median is below the relative-summary
floor. Query count alone is not a good proxy for the cost of either path.

The small-parallel fallback also misses streaming savings for restricted
cutoffs. At `1e10`, half-root and 4,096-term tails with 2/4/8 requested threads
can roughly halve warmed P2 time when streamed. Conversely, first and warmed
preferences genuinely differ for some ordinary cutoffs: serial isolated
Meissel at `2e9` favors fixed64 on the first call (0.431 versus 0.729 ms for
the selector), while warmed calls favor Phi (production 0.147 ms). A rule
based only on input, range and threads cannot always minimize both states.

### Interpretation and next discussion

Keep the size-dependent 8/16/32/64 strategy as a useful starting point. Before
promoting a selector, review a second temporary policy that considers the
actual interval length and worker budgets when deciding whether to stream a
short tail. It should also revisit the broad `x < 1e11` parallel fallback.
The present evidence does not specify an optimal replacement threshold or
justify runtime timing/cache inspection. Retain v1 and evaluate any revised
policy separately, with additional cases rather than silently fitting and
relabeling this run.

Raw samples are retained losslessly as [phases.csv.gz](results-selector/phases.csv.gz)
and [p2.csv.gz](results-selector/p2.csv.gz), alongside uncompressed
[complete summaries](results-selector/phases-summary.csv),
[isolated summaries](results-selector/p2-summary.csv), and
[aggregate definitions/results](results-selector/evaluation-summary.json).
The [retention manifest](results-selector/retained-timings.json) records raw
and compressed hashes and row counts. Full per-process CSVs remain in the
artifact directory.

## References and validation

Complete expected counts at all 11 input values agree between the installed
primecount 8.0 Gourdon and Deleglise–Rivat implementations, using one thread.
Known power-of-ten anchors are also checked. CLI commands, version and binary
hash are retained. The isolated suite has 49 input/label/table reference contexts;
12 cutoff labels at `1e8` also agree with a separately computed dense
pi oracle, covering 24 input/thread configurations. Larger isolated expectations
are explicitly labeled differential references, not independent proofs.

Six variants each pass 60,322 common checks under both ASan/UBSan and TSan.
Each sanitizer also passes 88 bounded stream-window checks and a separate
selector fixture with 78 decision cases and 54 kernel integration cases.
Coverage includes every input threshold, 4,095/4,096 query boundaries,
8,191/8,192 effective-worker boundaries, nonpositive budgets, oversized
contexts and concurrent differing choices. Kernel checks use independent
bounded pi oracles where feasible, with larger differential cases labeled.

The harness passed 628 result/CSV checks and seven expected invalid-input
rejections. Small arbitrary cutoffs agree with independent multiplicity and
dense-pi oracles. Final strict builds and sanitizer logs have no diagnostics.
An initial runner build omitted `selector_cases.h` from its snapshot copy;
that copy-list error was corrected before reference/probe/timing execution.

## Work and memory

All 404 first/warmed memory/work results pass expected-count and actual
kernel/trace checks. Choice and capacity fields are identical between repeats.
The policy exercises every decision reason; the largest tracked scratch bound
is 396,344 bytes (about 387 KiB), the same maximum as the earlier fixed64 study.
It excludes tables, Phi cache, stacks, ThreadBatch storage, allocator overhead
and diagnostic metadata. This is a conservative capacity bound, not peak RSS.
The selector itself does not allocate.

## Scope and reproduction

This is a temporary policy experiment on one machine. The matrix does not
establish an optimal rule, performance throughout the full supported range,
or behavior for every public `semiPrimeCount(n, a)` workload. In particular,
the isolated matrix focuses on Phi-only/8L/32L/64L choices; 16L is exercised
by complete first/warmed calls. No production promotion, R package check,
new dependency, or public API change is included.

From the repository root, create a fresh artifact directory and finish all
builds and checks before timing:

```sh
p2_selector_dir=$(mktemp -d)
mkdir "$p2_selector_dir/source"
git archive 8b43ed790435b22eb9ec02627379082fe7deb103 src inst/include | \
  tar -x -C "$p2_selector_dir/source"
git rev-parse 8b43ed790435b22eb9ec02627379082fe7deb103 > "$p2_selector_dir/commit.txt"
cp .github/development/p2-hybrid/results-selector/suite.json "$p2_selector_dir/suite.json"
python3 .github/development/p2-hybrid/selector_bench.py build "$p2_selector_dir"
python3 .github/development/p2-hybrid/validate.py "$p2_selector_dir" --variants baseline selector 8 16 32 64 --stream-window --selector-check
python3 .github/development/p2-hybrid/selector_bench.py references "$p2_selector_dir"
python3 .github/development/p2-hybrid/selector_bench.py memory "$p2_selector_dir"
python3 .github/development/p2-hybrid/selector_bench.py phases "$p2_selector_dir"
python3 .github/development/p2-hybrid/selector_bench.py p2 "$p2_selector_dir"
```

The reference step uses an already installed `primecount` CLI. The runner
checks source/binary/suite hashes and refuses to overwrite an existing timing
matrix. Results and metadata are retained under [results-selector](results-selector/).
The full artifact directory is recorded by `/tmp/rcppalgos-p2-selector-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-p2-selector-5kevyr5k`.
It contains frozen sources, the pre-selector snapshot, executables, logs,
reference/oracle outputs and every per-process timing CSV.

The final [preservation audit](results-selector/production-preservation.json)
confirms all 392 production files and all 61 prior-result files are unchanged.
All 282 build inputs and seven binaries retain their recorded hashes. The
policy and runner still match the versions frozen before timing. No changes
were staged or committed.
