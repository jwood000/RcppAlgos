# Interval-aware P2 selector v2

This is a second temporary selector, compared with frozen v1, production
Phi-backed P2, and fixed 8L/16L/32L/64L streaming variants. Production remains
at commit `8b43ed790435b22eb9ec02627379082fe7deb103`. The kernel, scheduling,
Phi cache, public interfaces and method dispatch are unchanged. Only the
experimental choice of where to stream differs.

The [v1 study](SELECTOR.md) identified two missed opportunities: narrow complete
ranges rejected by the 4,096-query guard, and small parallel calls rejected
without considering how much parallel work the stream could supply. This
policy is frozen before collecting its evaluation timings. It is not retuned
using the new results. Both header hashes and the complete matrix are recorded
in [suite.json](results-selector-v2/suite.json).

## Frozen policy

`P2SelectorV2.h` evaluates one candidate once per call:

1. Retain Phi for at most one non-table query.
2. Retain v1's input-based 8/16/32/64 cap and its promotion from 16 to 32 when
   `x >= 1e12` and the original Phi worker budget is at least four.
3. For `x < 1e11` with more than one usable Phi worker, try 16L instead of 8L.
4. Compute the exact capped split and the actual interval from L to the largest
   quotient assigned to streaming. The actual endpoint can be below the cap.
5. Accept a candidate covering the **entire non-table range** when the interval
   is at most 131,072 integers. This bypasses the remaining guards.
6. Otherwise retain Phi if there are fewer than 4,096 non-table queries, or
   fewer than 4,096 queries would stream.
7. For the remaining small parallel cases, accept the 16L candidate only if
   its stream worker budget is at least the original Phi worker budget.
8. Accept the remaining candidates.

The 131,072-integer threshold corresponds to two standard 32,768-odd-position
segments. Requiring the whole range avoids adding a stream alongside a Phi
prefix for this short-range exception. The one-query guard acknowledges that
stream initialization can exceed a warmed Phi lookup; it does not claim to
minimize first-call time.

Both worker budgets use the existing launch rules: Phi has a 2,048-query grain,
and streaming has a 1,048,576-integer grain. Requested thread counts are
normalized to at least one and bounded by available work. A stream budget is
not a promise that every quotient interval contains work; traces also record
the actual launches. Short complete ranges can profitably use one stream
worker even when the original Phi budget was larger, so that exception
precedes the budget comparison.

The policy uses no cache inspection, clocks or allocations. Diagnostics retain
the rejected candidate's span and worker budget as well as the actual choice.
Every positive selected multiplier matches an existing fixed variant; a
Phi-only choice has Q=L and no stream. The exact quotient partition and all
counting operations remain those of the earlier experiment. Preprocessed v1
and fixed32 builds, with and without instrumentation, are unchanged after
whitespace normalization.

## Evaluation design

Seven timing builds use the same production snapshot and native harness:
baseline, fixed8/16/32/64, v1 and v2. Separate instrumented v1/v2 executables
measure work and vector capacities. Strict C++11 builds use Apple Clang 17,
`-O3 -pthread`, debug symbols and frame pointers. All compilation, correctness
checks, references and memory probes finish before measured runs begin.

There are 88 complete-call configurations and 266 isolated P2 configurations.
Of these, 170 repeat selected v1 configurations and 184 add previously unmeasured
input/cutoff/table combinations. The additional inputs include `3e9`, `7e10`,
`9e10`, `7e11`, `7e12` and `7e13`. Cutoff labels add tails of 2/4/8/16/31/64/256/
1,024/2,048/4,094/4,097/8,192 terms and quarter/three-quarter square-root prime
indices. The actual integer cutoffs and definitions are retained in every row.
Labels can coincide numerically, and doubled tables can make a tail entirely
table-covered. This is a targeted matrix, not every possible combination.

Each configuration has three randomized matched blocks. Every variant starts
a fresh process and makes one first call followed by three warmed calls.
This gives 7,392 complete and 22,344 isolated timed results. The seeds are
`202611131` and `202611132`; the frozen suite records all cases.

Complete timing sums setup, Phi, P2, optional P3 and teardown, excluding the R
adapter and diagnostic collection. Isolated timing uses P2 alone with prepared
tables retained between calls; its first call has no preceding Phi phase.
Actual selection is timed inside P2, while trace collection follows the timed
kernel. Complete warmed calls rebuild tables but retain the process-wide Phi
cache. These definitions match the v1 study.

Summaries take the median of each process's warmed calls, then the median
across processes. Changes pair process medians within matched blocks. Controls
include production, the fastest fixed median, each selector's own chosen fixed
variant, and the other selector. The fastest fixed reference excludes both
selectors and is an optimistic hindsight comparison. Same-choice differences
can include code generation and system variation; they do not isolate selector
instruction overhead. Every sample is retained. For isolated cases whose
fastest fixed median is below 0.01 ms, emphasize absolute differences instead
of unstable percentages. The evaluation also separates repeated configurations
from the new combinations.

## Results

All **29,736 timed results** match their expected values: 22,344 isolated P2
results and 7,392 complete-call results. P2 subtotals also agree across modes,
variants, thread counts and repeats wherever the numeric input, cutoff and
table limit agree. The policy has not been adjusted in response to its results.

### Isolated P2

V2 removes the severe narrow-tail misses and improves coverage of the fastest
fixed candidates. However, it is not a strict improvement over v1. The worker
budget test admits some candidates that substantially regress warmed calls
with enlarged tables, and the short-range exception trades first-call savings
for microsecond-scale warmed costs on some small inputs.

The table counts configurations finishing within 5% of the fastest fixed
candidate, using matched-block comparisons. Relative aggregates include only
cases whose fastest fixed median is at least 0.01 ms. No raw cases are dropped.

| Configurations | First v1 | First v2 | Warmed v1 | Warmed v2 |
| --- | ---: | ---: | ---: | ---: |
| All eligible | 83 / 180 | 124 / 180 | 126 / 166 | 147 / 166 |
| Repeated from v1 | 39 / 90 | 57 / 90 | 70 / 86 | 81 / 86 |
| New combinations | 44 / 90 | 67 / 90 | 56 / 80 | 66 / 80 |

The remaining 86 first and 100 warmed configurations are assessed with
absolute differences. V2 changes the actual partition in 106 of 266 isolated
configurations: 102 short-complete-range choices and four worker-coverage
choices. First calls improve in all three matched blocks in 104 of those 106
cases. Warmed calls improve in 68, regress in 36 and have zero paired median
difference in two, classifying each case by its paired median. The tiny-call
tradeoffs should not be hidden by a single
percentage summary.

Among changed cases above the relative floor, 46/48 improve first and 28/34
improve warmed, with median paired changes of -84.71% and -70.34% against v1.
For the 132 unchanged-partition cases above the floor in each state, median
paired differences are +0.120% first and +0.034% warmed. These are unweighted
distribution summaries over selected cases, not an overall user-workload
speedup or an isolated measurement of selector overhead.

Representative serial narrow-tail times (milliseconds, medians of process
medians) show both a repeated v1 configuration and an additional input:

| Input and remaining terms | First v1 | First v2 | Warmed v1 | Warmed v2 |
| --- | ---: | ---: | ---: | ---: |
| `3e13`, 4,095 terms | 1.877 | 0.076 | 1.083 | 0.033 |
| `7e12`, 4,094 terms | 1.386 | 0.076 | 0.658 | 0.031 |

The new `7e12` case streams 62,044 integers with one worker, improving about
95% in both states by paired comparison. At these cutoffs every fixed cap
covers the same entire range, so this supports allowing the stream rather
than selecting a particular multiplier.

The main new regressions occur when a doubled table changes the worker
balance. These isolated cases use the Meissel cutoff and four requested
threads; they are **not complete Meissel timings**:

| Input, table multiplier 2 | First v1 | First v2 | Warmed v1 | Warmed v2 |
| --- | ---: | ---: | ---: | ---: |
| `1e10` | 1.006 | 0.707 | 0.256 | 0.470 |
| `7e10` | 2.158 | 1.302 | 0.511 | 0.754 |

First calls improve 29.81% and 39.44%, but warmed calls regress **83.54% and
47.99%**. V2 chooses 16L with two/two and four/four Phi/stream worker budgets,
respectively. It tracks its corresponding fixed16 control closely, so these
are poor warmed choices rather than unexplained selector overhead. Matching
worker budgets does not establish matching costs: interval size, the remaining
Phi prefix, seed work and multiple launch phases can matter. Their individual
contributions have not been profiled here.

The zero-cutoff doubled-table cases also regress versus v1: first calls by
11.73% at `1e10` and 11.03% at `7e10`, and warmed calls by 33.43% and 14.37%.
Thus the worker-coverage exception has regressions beyond a first/warmed
preference reversal.

Short complete streams also have limits. At `3e9`, the serial 2,048-term case
improves from 0.186 to 0.034 ms first but regresses from 0.0074 to 0.0157 ms
warmed: a paired increase of about 8.3 us. Its largest quotient is 91,271,
within the bounded Phi cache's argument range. Absolute quotient size can
therefore matter alongside interval length. The production Phi control is
0.0174 ms warmed in this case despite using the same mathematical path as
v1; this implementation/timing difference prevents attributing every gap to
the selector or making a blanket claim that Phi is faster there. The other
small warmed regressions are retained in the raw data, including two-term
tails where the added time is a fraction of a microsecond.

The rule remains too conservative elsewhere. At `7e10`, an 8,192-term tail
with four requested threads is rejected for insufficient stream workers,
yet a single fixed stream is substantially faster. This is a missed saving
already present in v1, distinct from the new regressions above.

### Complete calls

V2 changes the partition in eight of the 88 complete-call configurations,
always replacing the small-parallel Phi fallback with 16L. All eight improve
in **every matched block**, both first and warmed. The following values are
paired median time changes against v1; negative means less time:

| Method | Input | Threads | First change | Warmed change |
| --- | --- | ---: | ---: | ---: |
| Meissel | `3e10` | 2 | -20.05% | -6.82% |
| Lehmer | `3e10` | 2 | -17.35% | -3.22% |
| Meissel | `7e10` | 2 | -25.46% | -10.59% |
| Lehmer | `7e10` | 2 | -25.20% | -9.40% |
| Meissel | `9e10` | 2 | -24.66% | -11.36% |
| Lehmer | `9e10` | 2 | -25.25% | -10.07% |
| Meissel | `9e10` | 4 | -19.87% | -3.54% |
| Lehmer | `9e10` | 4 | -19.23% | -1.40% |

Across all 88 configurations, proximity to the fastest fixed candidate improves:

| Within threshold of fastest fixed | First v1 | First v2 | Warmed v1 | Warmed v2 |
| --- | ---: | ---: | ---: | ---: |
| 5% | 54 / 88 | 63 / 88 | 79 / 88 | 82 / 88 |
| 10% | 67 / 88 | 74 / 88 | 85 / 88 | 87 / 88 |
| 20% | 74 / 88 | 81 / 88 | 88 / 88 | 88 / 88 |

Among the 80 unchanged-partition cases, median v2/v1 differences are +0.105%
first and +0.054% warmed. Individual differences are larger: first Lehmer at
`3e10`/eight threads is 7.23% slower, and warmed Lehmer at `9e10`/one thread
is 9.80% slower, despite identical choices. These comparisons include possible
compiler and system variation and do not isolate selection-instruction cost.
Across all cases, v2 versus its own chosen fixed control has median differences
of +0.033% first and +0.168% warmed.

The complete matrix still exposes missed first-call savings: eight-thread
Lehmer at `9e10` retains Phi and is 26.67% slower than fixed32. The budget
guard is therefore too conservative for some inputs as well as too permissive
for some enlarged-table isolated cases. Improved coverage does not establish
that worker parity is a sufficient cost model.

### Interpretation

The interval exception removes important misses, and the revised parallel
choice helps the standard complete calls where it changes their work. Keep
these as separate useful findings. A production decision should also consider
the short-range warmed tradeoffs and enlarged-table regressions; frozen v2
should not be described as uniformly superior to v1.

The next discussion can focus on whether to account for absolute quotient
size/cache eligibility in the short-range rule and how to compare the cost of
a Phi prefix plus streaming against Phi alone. This experiment does not
identify an optimal replacement rule or isolate launch, seed and cache costs.
No further selector has been implemented or production behavior changed.

## Validation and memory

Complete expected counts at all 13 inputs agree between the installed
primecount Gourdon and Deleglise–Rivat implementations. There are 121 distinct
input/cutoff-label/table P2 reference contexts. The 24 configurations at `1e8`
also agree with the independent dense-pi oracle; larger isolated references are
explicitly production-P2 differential checks. All 170 repeated configurations
retain the same expected results as the v1 study.

All 1,416 first/warmed memory/work results pass expected-count and actual
kernel/trace checks. V2 changes the actual Phi/stream partition in 114 of 354
configurations and leaves it unchanged in 240. Both versions reach the same
maximum tracked scratch capacity, 396,344 bytes (about 387 KiB). This excludes
prepared tables, the Phi cache, stacks, ThreadBatch storage, allocator overhead
and diagnostic metadata; it is not peak RSS. Selection itself allocates no
storage. Every v2 policy reason is exercised, in addition to empty/direct-only
bypasses, and decisions/capacities remain identical between first and warmed
probes.

All seven variants pass 60,322 common checks under both ASan/UBSan and TSan.
Each sanitizer also passes the 88 bounded stream-window checks, the unchanged
v1 fixture (78 choices and 54 kernel cases), and the new v2 fixture (201 choices
and 173 kernel cases). V2 coverage includes exact 131,071/131,072/131,073 spans,
both ends of floor-division plateaus, full versus partial ranges, 1/2 and
4,095/4,096/4,097-query boundaries, two/four-worker parity boundaries,
nonpositive budgets, oversized contexts, MAX_X choices and concurrent opposing
decisions. Of the 173 v2 kernel cases, 171 use independent bounded sieve
comparisons and two are explicitly differential; additional concurrent checks
share immutable contexts with thread-local diagnostics.
All 56 build/run logs are free of compiler warnings and sanitizer diagnostics.

The separate harness verification passes 1,325 result/schema/cutoff/trace rows
and eight intentional argument rejections. All 26 labels are checked against
independent trial-factor multiplicities and dense-pi results at small inputs.
An independent linear partition checks diagnostic indices, spans and budgets;
instrumented snapshots agree. Six complete-formula checks at `1e6` also pass.
No experimental source fixes were required after freezing the policies.

## Reproduction and scope

This is a native experiment on one machine. It does not modify production or
the public R interface, and does not include a full R package check or establish
performance throughout the supported range. First/warmed tradeoffs remain
possible for a policy that deliberately ignores cache state.

From the repository root, prepare a fresh artifact directory:

```sh
p2_selector_v2_dir=$(mktemp -d)
mkdir "$p2_selector_v2_dir/source"
git archive 8b43ed790435b22eb9ec02627379082fe7deb103 src inst/include | \
  tar -x -C "$p2_selector_v2_dir/source"
git rev-parse 8b43ed790435b22eb9ec02627379082fe7deb103 > "$p2_selector_v2_dir/commit.txt"
cp .github/development/p2-hybrid/results-selector-v2/suite.json "$p2_selector_v2_dir/suite.json"
python3 .github/development/p2-hybrid/selector_v2_bench.py build "$p2_selector_v2_dir"
python3 .github/development/p2-hybrid/validate.py "$p2_selector_v2_dir" --variants baseline selector selector-v2 8 16 32 64 --stream-window --selector-check --selector-v2-check
python3 .github/development/p2-hybrid/selector_v2_bench.py references "$p2_selector_v2_dir"
python3 .github/development/p2-hybrid/selector_v2_bench.py memory "$p2_selector_v2_dir"
python3 .github/development/p2-hybrid/selector_v2_bench.py p2 "$p2_selector_v2_dir"
python3 .github/development/p2-hybrid/selector_v2_bench.py phases "$p2_selector_v2_dir"
python3 .github/development/p2-hybrid/selector_v2_analyze.py "$p2_selector_v2_dir"
```

Reference generation uses the already installed primecount CLI. The runner
verifies source, executable and policy hashes and refuses to overwrite a timing
matrix. The full artifact path is recorded in `/tmp/rcppalgos-p2-selector-v2-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-p2-selector-v2-pkbdf4cx`.

All raw samples are retained losslessly in [p2.csv.gz](results-selector-v2/p2.csv.gz)
and [phases.csv.gz](results-selector-v2/phases.csv.gz). Uncompressed
[isolated summaries](results-selector-v2/p2-summary.csv),
[complete summaries](results-selector-v2/phases-summary.csv), and
[aggregate definitions/results](results-selector-v2/evaluation-summary.json)
support the tables above. The [retention manifest](results-selector-v2/retained-timings.json)
records raw/compressed hashes and row counts; full per-process CSVs remain in
the artifact directory. No sample was discarded.

The [preservation audit](results-selector-v2/production-preservation.json)
checks all 392 production files, all 90 earlier result files, both frozen
policies, and the timing build inputs/executables. No changes were staged or
committed. The experiment's source snapshot and hashes are retained separately
from v1 in [experiment-sources.json](results-selector-v2/experiment-sources.json).
