# Stream-density refinement: temporary P2 selector v3

V3 refines one decision in [frozen v2](SELECTOR_V2.md): whether a small-input
parallel candidate should stream after its worker budget matches the original
Phi budget. Production remains unchanged at commit
`8b43ed790435b22eb9ec02627379082fe7deb103`. V1, v2, their harnesses and all
earlier result files are preserved. Counting, scheduling, cache settings and
public interfaces are unchanged.

## Frozen rule and rationale

V2's four doubled-table admissions scanned approximately 637 or 701 integers
per streamed query. Its eight successful complete-call admissions scanned
approximately 180–189. Matching worker budgets therefore did not distinguish
these different amounts of sieve work per answer.

V3 adds this requirement **only** after a candidate passes v2's small-parallel
worker-coverage check:

```text
actual_stream_span <= 256 * streamed_query_count
```

The implementation compares `span / 256 + (span % 256 != 0)` with the query
count, avoiding either an overflowing multiplication or `span + 255`.
It uses the largest actual streamed quotient minus the prepared table limit,
not the possibly larger candidate cap. A rejected candidate returns Phi-only
with reason `stream-too-sparse`, retaining candidate diagnostics for inspection.

256 is a round experimental threshold, not a proven crossover or a fitted
timing model. It is frozen before collecting this matrix. V2's short-full-range
exception retains precedence, including its known small warmed-call tradeoffs.
Serial calls, existing query/worker fallbacks, and the rules for `x >= 1e11`
are unchanged. In particular, this trial does not claim to address every
enlarged-table workload or v2's other missed opportunities.

All selectors use immutable tables and exact integer partitions. Selection
does not inspect cache state, time execution, allocate storage, or change any
counting operation. Earlier v1/v2/fixed32 preprocessed kernels, with and without
instrumentation, remain identical after whitespace normalization.

## Evaluation design

Eight timing builds compare production Phi-only, fixed8/16/32/64, v1, v2 and
v3. Instrumented v1/v2/v3 builds separately measure work and vector capacities.
Strict C++11 builds use Apple Clang 17, `-O3 -pthread`, debug symbols and frame
pointers. All builds, sanitizer checks, reference generation, scope checks and
memory probes finish before measured runs begin.

The frozen [suite](results-selector-v3/suite.json) has 228 isolated P2 and
30 complete-call configurations. The isolated grid emphasizes zero, Meissel
and Lehmer cutoffs at `1e10`, `5e10`, `7e10` and `8e10`, with 2/4-thread grids
and selected serial/eight-thread controls. Additional cases cover the short
exception, empty/direct results, the `1e11` policy boundary, and small independent
oracles. Complete calls retain v2's eight changed configurations, additional
small inputs and controls through `1e14`.

Prepared limits use exact ratios 1, 9/8, 6/5, 5/4, 3/2, 7/4, 2 and 3:
`L = floor(isqrt(x) * numerator / denominator)`. This means enlarging the
integer-root limit, not scaling the real square root before flooring. The
product and native table ceiling are checked before allocation. Both ratio
components and actual L appear in every row. The legacy `table_multiplier`
column means the numerator; use the explicit `table_numerator` and
`table_denominator` together for fractional contexts. Cutoffs still use
`b = pi(isqrt(x))`, independently of table enlargement.

49 configurations repeat exact labeled v2 configurations; 209 add combinations
not present in that matrix. Some labels or ratios can describe identical numeric
work. These are evaluation cohorts, not statistically independent populations.
The suite is targeted rather than a full Cartesian product.

Each configuration uses three randomized matched blocks, eight fresh variant
processes per block, and one first plus three warmed calls per process. There
are 21,888 isolated and 2,880 complete timed results (24,768 total). Seeds are
`202611141` for isolated P2 and `202611142` for complete calls.

Timing definitions match v2: isolated `p2_ms` excludes retained-table setup;
complete totals sum setup, Phi, P2, P3 when applicable and teardown, excluding
diagnostic collection and the R adapter. First isolated calls have no preceding
Phi phase. Complete warmed calls rebuild tables but retain the process-wide
Phi cache. Actual selection occurs inside the P2 timer; trace collection follows
it. No competing builds or validation workloads run during timing.

Comparisons pair process medians within matched blocks, reporting first and
warmed calls separately. The fixed reference is the fastest median among
production and fixed8/16/32/64, excluding all selectors; it is an optimistic
hindsight comparison. Own-choice controls help distinguish changed work from
compiler/system variation and do not isolate the selector's instruction cost.
The predeclared 0.01 ms fastest-fixed floor governs relative isolated aggregates;
smaller cases use absolute differences. No samples are discarded or policy
constants retuned after timing.

## Results

All **24,768 timed results** match their expected values: 21,888 isolated P2
and 2,880 complete-call results. The rule remains frozen at 256. Results below
use contemporaneous controls, not timings copied from the earlier studies.

### The four doubled-table regressions

V3 rejects all four problematic v2 admissions and restores the Phi path.
Warmed calls improve by 12–46% versus v2. The zero-cutoff cases also improve
first calls, while the Meissel cutoffs give back v2's first-call savings:

| Input | Cutoff, four threads and 2x table | First v3 vs v2 | Warmed v3 vs v2 |
| --- | --- | ---: | ---: |
| `1e10` | Zero | -11.81% | -26.97% |
| `1e10` | Meissel | +58.26% | -45.53% |
| `7e10` | Zero | -6.27% | -12.12% |
| `7e10` | Meissel | +66.54% | -31.04% |

These are paired median changes in isolated P2 time; negative means less time.
For example, the `1e10` Meissel-cutoff warmed medians are 0.468 ms for v2 and
0.256 ms for v3, close to v1's 0.255 ms. First medians are 0.701 and 1.078 ms.
The new choice expresses a first/warmed tradeoff, not a uniform improvement.

### New rational-table cases

V3 changes 99 isolated partitions relative to v2. Warmed paired medians improve
in 89 and regress in ten; first medians improve in 44 and regress in 55.
Across those changed cases, the median paired change is **-14.78% warmed** and
**+45.51% first**. Among the 95 newly changed combinations, 85 improve warmed.
These are unweighted summaries of this focused matrix, not a predicted average
speedup for a user's workload.

A new positive example is the Meissel cutoff at `5e10`, four threads and a
3/2-sized table. Its candidate density is about 398.79 integers per query.
V3 improves warmed P2 from 0.579 to 0.463 ms (-20.11% paired), while first calls
regress 71.03%. All three matched blocks agree on both directions.

The boundary cases show why 256 is not a universal crossover. At `8e10`, two
threads and a 6/5-sized table, both cutoffs below propose the **same** streaming
tail: 5,089,103 integers for 19,161 queries, or about 265.60 integers per query.
Both have the same Phi and streaming worker budgets. V3 rejects both:

| Cutoff | First v2 ms | First v3 ms | Warmed v2 ms | Warmed v3 ms | Warmed paired change |
| --- | ---: | ---: | ---: | ---: | ---: |
| Zero | 22.111 | 21.103 | 8.422 | 7.634 | -10.62% |
| Meissel | 2.575 | 5.270 | 1.258 | 1.374 | +8.86% |

Times are medians of process medians; percentages use paired process ratios
and therefore need not equal ratios of the displayed medians. Every warmed
matched block favors Phi for the zero cutoff and streaming for the Meissel
cutoff. The Meissel first call also regresses in every pair, by 102.71% at the
paired median. V3 remains close to its selected Phi control in this example,
so this is a poor choice rather than an unexplained selector-time outlier.

The shared tail does not make the whole calls equivalent: their preceding Phi
work differs. This suggests investigating interactions among prefix work, cache
reuse, scheduling and stream setup. This study has not profiled their individual
costs. It shows empirically that density and worker parity alone do not predict
the preferred choice for these contexts.

In an exploratory, post-hoc grouping, all ten new warmed regressions fall
between 256 and 300 integers per query. All 76 changed cases above 300 improve
at their warmed paired medians. This observation is not a validated 300 rule;
the same-density counterexample above would remain a problem for any rule that
uses only density. A retained-stream limitation also remains below 256: the
zero cutoff at `8e10`, two threads, 9/8 table and density 234.66 is 11.26% slower
warmed than production Phi, despite matching its fixed16 control closely.

Some large observed differences should not be treated as clean policy effects.
At `8e10`, Lehmer cutoff, four threads and 6/5 table, warmed v3 is 48.26% slower
than v2, but it is also 47.95% slower than its own Phi control. V3 process
medians vary from 1.086 to 1.880 ms. All samples are retained; compiler/system
variation is not isolated here. The matched-cutoff example above is a stronger
basis for judging the decision rule.

### Comparison with the fixed variants

Counts below use matched comparisons against the fastest fixed candidate. The
relative isolated aggregate includes only cases with a fastest fixed median
of at least 0.01 ms; 31 first and 37 warmed configurations are retained and
assessed using absolute differences instead.

| Within 5% of fastest fixed | v1 | v2 | v3 |
| --- | ---: | ---: | ---: |
| Isolated first | 80 / 197 | 108 / 197 | 91 / 197 |
| Isolated warmed | 150 / 191 | 80 / 191 | 149 / 191 |
| Complete first | 13 / 30 | 25 / 30 | 27 / 30 |
| Complete warmed | 15 / 30 | 26 / 30 | 27 / 30 |

This matrix deliberately emphasizes enlarged tables, so these counts should
not be compared directly with differently composed v1/v2 study totals.

All 30 complete-call partitions are identical between v2 and v3, including
v2's eight successful admissions. Median v3/v2 differences are +0.085% first
and -0.155% warmed; individual case changes range from -5.83% to +4.12% first
and -3.40% to +3.58% warmed. Their proximity-count changes are same-work timing
variation, not demonstrated gains from the new rule. Median changes against
v3's own selected fixed control are -0.126% first and -0.006% warmed.

### Interpretation

The extra guard recovers the targeted warmed performance and helps many new
enlarged-table cases while preserving v2's complete-call choices. Its cost is
substantial first-call regressions for many standard cutoffs, and it introduces
real warmed misses near the threshold. It remains a temporary experiment.

The most informative next investigation is the pair of cutoffs with identical
stream work but opposite warmed preferences. Understanding their prefix/cache/
scheduling costs would give a stronger basis for another policy. No threshold
was retuned, further selector implemented, or production change made here.

## Validation and memory

Complete reference counts at nine inputs agree between the installed primecount
Gourdon and Deleglise–Rivat implementations. There are 125 isolated input/label/
ratio reference contexts; 30 configurations at `1e8` also agree with an
independent dense-pi oracle. Larger isolated expectations are explicitly
production-P2 differential references. All 49 repeated v2 configurations retain
their original expected results. Both reference generation and timed checks
require P2 results to agree for the same numeric `(x, a)` across different valid
prepared-table limits, in addition to variant/thread comparisons.

All 1,548 first/warmed memory/work probes pass. V3 changes the partition in 99
of the 228 isolated configurations: the four repeated doubled-table admissions
and 95 new combinations. All 30 complete-call partitions agree with v2. The
largest tracked scratch capacity is 396,344 bytes (about 387 KiB) for each
selector. This excludes tables, the Phi cache, stacks, ThreadBatch storage,
allocator overhead and diagnostic metadata; it is not peak RSS. Choice and
capacity fields are identical between first/warmed probes.

A separate [scope replay](results-selector-v3/scope-check/results.json) checks
actual frozen v2/v3 choices on all 354 prior contexts. One independently sieved
prime vector covers their limits; the replay calls neither Phi nor P2. All
332 evaluated candidate descriptions match, and 22 empty/direct-table bypasses
are preserved. Exactly the four known doubled-table admissions change; the
other 350 contexts remain unchanged. Source, input, commands and hashes are
retained with the replay. This is an executable scope check, not a proof of
performance or correctness over every input.

Both ASan/UBSan and TSan pass 60,322 common checks for each of eight variants,
88 stream-window checks, and the unchanged v1/v2 policy fixtures. The new v3
fixture passes 163 explicit choices and 147 kernel cases against an independent
sieve, plus repeated concurrent comparisons. Its 4,724 assertions include exact
`256N-1`, `256N`, `256N+1` spans, a cap-versus-actual-endpoint counterexample,
table ratios, guard precedence, worker boundaries, serial/large-input exceptions
and pure MAX_X choices. All 64 build/run logs are free of compiler warnings
and sanitizer diagnostics.

The separate harness check passes 2,870 rows, including 2,582 independently
computed linear partition/worker traces, all 26 cutoff labels and eight table
ratios. Small results agree with independent multiplicity/dense-pi oracles;
29 larger contexts are explicitly differential. It also checks 16 complete
formula rows, 14 intentional argument rejections, and rational/MAX_X table-limit
validation without giant allocations. No frozen policy or harness corrections
were needed during preflight.

## Reproduction and scope

This is a native experiment on one machine; no full R package check or other
platform run is included. Production and the R API are unchanged. From the
repository root, prepare a fresh artifact directory:

```sh
p2_selector_v3_dir=$(mktemp -d)
mkdir "$p2_selector_v3_dir/source"
git archive 8b43ed790435b22eb9ec02627379082fe7deb103 src inst/include | \
  tar -x -C "$p2_selector_v3_dir/source"
git rev-parse 8b43ed790435b22eb9ec02627379082fe7deb103 > "$p2_selector_v3_dir/commit.txt"
cp .github/development/p2-hybrid/results-selector-v3/suite.json "$p2_selector_v3_dir/suite.json"
python3 .github/development/p2-hybrid/selector_v3_bench.py build "$p2_selector_v3_dir"
python3 .github/development/p2-hybrid/validate.py "$p2_selector_v3_dir" --variants baseline selector selector-v2 selector-v3 8 16 32 64 --stream-window --selector-check --selector-v2-check --selector-v3-check
python3 .github/development/p2-hybrid/selector_v3_bench.py references "$p2_selector_v3_dir"
python3 .github/development/p2-hybrid/selector_v3_bench.py memory "$p2_selector_v3_dir"
python3 .github/development/p2-hybrid/selector_v3_bench.py p2 "$p2_selector_v3_dir"
python3 .github/development/p2-hybrid/selector_v3_bench.py phases "$p2_selector_v3_dir"
python3 .github/development/p2-hybrid/selector_v3_analyze.py "$p2_selector_v3_dir"
```

References use the already installed primecount CLI. The runner verifies policy,
source, executable and suite hashes and refuses to overwrite a timing matrix.
The full artifact path is recorded in `/tmp/rcppalgos-p2-selector-v3-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-p2-selector-v3-ag0tbj0m`.

All samples are retained losslessly in [p2.csv.gz](results-selector-v3/p2.csv.gz)
and [phases.csv.gz](results-selector-v3/phases.csv.gz). Uncompressed
[isolated summaries](results-selector-v3/p2-summary.csv),
[complete summaries](results-selector-v3/phases-summary.csv) and
[aggregate definitions/results](results-selector-v3/evaluation-summary.json)
support the report. The [retention manifest](results-selector-v3/retained-timings.json)
records raw/compressed hashes and row counts. Full per-process CSVs remain in
the artifact directory. All 48 distinct numeric `(x, a)` counts agree across
137 `(x, a, L)` contexts, including different rational table sizes.

The [final preservation audit](results-selector-v3/production-preservation.json)
confirms all 392 production files and all 125 earlier result files are unchanged.
All 284 build inputs and 11 binaries retain their recorded hashes, as do the
three frozen policies and timing runner. Experiment sources and hashes are
retained in [experiment-sources.json](results-selector-v3/experiment-sources.json).
Nothing was staged or committed.
