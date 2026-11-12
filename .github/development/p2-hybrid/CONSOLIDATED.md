# Consolidated P2 candidate: common scheduling and conservative streaming

The temporary candidate is ready for review as a first production policy.
It combines the useful initial scheduling change with the established large-input
streaming caps, and limits small parallel admissions to a conservative region.
The fresh-input results support this direction. It does not select the fastest
measured variant for every workload, and the recurring small first-call regression
remains visible below.

**Production is unchanged.** All work is under `.github/development/p2-hybrid/`.
The source baseline is commit `8b43ed790435b22eb9ec02627379082fe7deb103`.
No public API, Legendre implementation, Phi cache, thread budget, counting identity,
or package requirement changes in this experiment. No commits were made.

## What the candidate does

Both Phi-only P2 and the hybrid's Phi prefix now use the same queue:
the first 32 non-table queries are individual jobs, followed by 32-query jobs.
Workers claim from one queue; there is no barrier after the first 32 queries.
The existing barrier between the entire Phi prefix and the segmented stream
remains. Serial evaluation and the full-range worker-budget calculation are
unchanged. The queue allocates no storage and makes at most 31 additional
successful claims compared with the original 32-query scheduling.

The [candidate selector](consolidated/P2SelectorCandidate.h) is one deterministic
policy. It does not switch among version names, inspect the cache, calibrate
itself at runtime, or identify the calling method.

Let `L` be the prepared pi-table limit. Its proposed stream cap retains the
earlier schedule:

| Input | Proposed cap |
| --- | --- |
| `x < 1e11` | `8L` serial; `16L` when the Phi budget exceeds one |
| `1e11 <= x < 3e12` | `16L`; promoted to `32L` for `x >= 1e12` and Phi budget >= 4 |
| `3e12 <= x < 3e13` | `32L` |
| `x >= 3e13` | `64L` |

The cap cannot exceed the largest actual quotient. Zero or one non-table query
retains Phi before any streaming exception is considered. As in v2, a complete
non-table range spanning at most 131,072 integers can stream without the usual
query-count guards. Otherwise, both the original non-table range and proposed
streamed tail must contain at least 4,096 queries.

For **small parallel inputs** (`x < 1e11`, Phi budget > 1), streaming additionally
requires all three conditions:

1. The stream worker budget covers the original Phi worker budget.
2. `L == floor(sqrt(x))`.
3. The first remaining prime exceeds `floor(x^(1/4))`.

The last two conditions are the new conservative guards. They admit the usual
Lehmer/Meissel cutoffs and larger numeric cutoffs with root-sized tables; earlier
cutoffs and enlarged tables retain Phi with K=1 initial scheduling. These are
empirical limits on admission, **not mathematical performance boundaries**.
They intentionally avoid another fitted stream-density threshold. The short
complete-range exception takes precedence over these guards.

The [generated hybrid](results-consolidated/generated/P2HybridInitial.cpp) retains
the existing counting and streaming code. The exact partition remains
`floor(x/p) <= Q` iff `p > floor(x/(Q+1))`; each query belongs to one of the Phi,
stream, or direct-table ranges. The queue's disjoint ticket ranges cover its
input exactly once. These invariants explain why changing dispatch/scheduling
preserves the P2 sum; the tests independently check its counts.

## Prospective comparison

The [plan](consolidated/plan.json), candidate and both cohorts were frozen before
timing. The [suite](results-consolidated/suite.json) records the actual seeds and
source hashes. The candidate was not retuned after either cohort.

The focused cohort contains 45 configurations and ten variants: original
production, Phi-only K=1, v1/v2/v3, fixed 8/16/32/64, and the candidate. The old
selectors' **decision logic** is frozen, but all hybrid controls use the common
K=1 prefix scheduling. Their timings therefore are new comparisons, not copies
of the earlier studies.

The held-out cohort contains 28 configurations at eight previously untimed inputs
in these retained studies: `4e9`, `6e10`, `8.5e10`, `1.1e11`, `8e11`, `1.2e12`,
`6e12`, and `8e13`. Its controls are production, Phi-only K=1, and fixed 16/32/64.
It includes three/six requested threads, arbitrary cutoffs, root and enlarged
tables, short complete ranges, and complete Meissel/Lehmer computations.
The [focused decision record](results-consolidated/focus-decision.json) was saved
before running this cohort.

Each fresh process supplies one first call and three warmed calls. Cases below
`1e13` use five randomized matched blocks; larger cases use three. First calls
are separate from the median warmed call within each process. Percentage changes
are medians of **paired block ratios**, so they need not equal ratios of the
displayed median times. No samples were discarded. Builds, sanitizer runs,
references and allocation probes finished before timing began.
Timing builds use Apple Clang 17 with `-std=c++11 -O3 -pthread`, debug symbols
and frame pointers, consistently across variants.

Isolated P2 timings exclude prepared-table setup and have no preceding Phi call.
Complete timings reproduce the native setup/Phi/P2/P3/teardown sequence and sum
those phases; they exclude the R adapter and diagnostics. They rebuild tables
on each call while retaining the process-wide Phi cache. Actual dispatch runs
inside the P2 timer; diagnostic selection descriptions run afterward.

The predeclared regression screen requires a paired median increase exceeding
both 5% and 0.05 ms against production or K=1 Phi, with every block slower.
Short-full-range warmed calls instead have a 0.02 ms absolute allowance.
These are review criteria, not statistical guarantees or permission to ignore
smaller systematic costs.

## Results

All **12,104 timed results** matched their expected counts: 8,840 focused and
3,264 held-out results, from 3,026 fresh processes. There was one focused
regression flag and no held-out flags.

### Complete computations

Representative warmed native totals, in milliseconds:

| Cohort | Input | Method | Requested threads | Production | Candidate | Paired change |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| Focus | `1e12` | Meissel | 8 | 12.237 | 10.175 | -18.21% |
| Focus | `1e12` | Lehmer | 8 | 11.994 | 9.988 | -16.27% |
| Focus | `1e14` | Meissel | 8 | 618.226 | 325.777 | -47.07% |
| Focus | `1e14` | Lehmer | 8 | 602.999 | 309.025 | -48.36% |
| Held out | `8e11` | Meissel | 3 | 14.519 | 12.052 | -16.99% |
| Held out | `8e11` | Lehmer | 3 | 14.429 | 11.866 | -17.62% |
| Held out | `8e13` | Meissel | 6 | 568.552 | 308.136 | -45.80% |
| Held out | `8e13` | Lehmer | 6 | 559.999 | 291.502 | -46.14% |

Every matched block improves in these examples. First calls improve too:
approximately 30–31% at `1e12`/eight threads, 44–46% at `1e14`/eight threads,
and 42–43% at held-out `8e13`/six threads. Comparisons against K=1 Phi retain
the large gains, so these are not merely the initial scheduling improvement.

Small complete calls need not improve. At `1e10`, eight-thread warmed Lehmer is
1.80% slower than production, about 0.011 ms, in every block. At held-out
`8.5e10`/six threads the worker-budget guard retains Phi; complete times are
near the controls, with mixed first-call changes. At `6e10`/three threads,
admitted streaming improves complete first calls by about 20%, while warmed
improvements are only about 2–4%.

### Initial scheduling and conservative fallbacks

At the matched `8e10`, zero-cutoff, `6/5`-table context, the candidate retains
Phi. Warmed P2 changes from 7.552 to 2.750 ms with eight requested threads
(-63.50% paired), and from 7.599 to 6.298 ms with two (-17.62%). It remains
close to the corresponding K=1 Phi control. Fresh zero-cutoff inputs also
benefit: warmed reductions against production are 40.35% at `6e10`/three
threads and 61.08% at `8.5e10`/six threads.

The four previously problematic doubled-table cases improve warmed P2 by about
8–47% against v2 on the common scheduler. This protection has a first-call
cost: the Meissel-cutoff cases give back v2's first-call savings, becoming
44.83% and 61.83% slower than v2 at `1e10` and `7e10`, respectively.

The policy also deliberately misses useful streams. At matched `8e10`,
Meissel cutoff, two threads and `6/5` table, its Phi fallback is 111.90% slower
first and 12.02% slower warmed than v2. With K=1 shared by both paths, even
the matched zero-cutoff case favors v2 modestly warmed: the candidate is
2.74% slower. The original opposite-preference example therefore cannot be
used unchanged to justify a new dispatch boundary after changing scheduling.

Against the fastest hindsight control (production, K=1 Phi, or fixed cap),
the focused candidate is within 5% in 25/27 warmed isolated cases, versus
22/27 for each of v2 and v3 and 17/27 for v1. First calls reverse some of that
advantage: candidate 16/27, v2 22/27, v3 18/27, v1 7/27. In held-out cases,
the candidate is within 5% in 14/18 warmed isolated cases and all 10 warmed
complete cases. These targeted, unweighted counts describe a tradeoff, not
an expected speedup for a user's workload. The best-control floor is 0.01 ms;
every case in these counts exceeds it. These proximity counts use ratios of
aggregate median times; the reported percentage changes elsewhere use paired
block ratios.

Held-out examples expose the same conservatism: warmed zero-cutoff P2 at
`6e10` is 4.011 ms versus fixed16's 3.603 ms; the slightly enlarged
`1001/1000` table with a Meissel cutoff is 0.846 versus 0.787 ms. These are
missed opportunities against a hindsight control, not regressions against
the retained Phi path. We have not added exceptions to capture them.
That 0.1% table enlargement exposes the strict guard's discontinuity: first
calls are 89.92% slower than fixed16 in every block, while the warmed 9.92%
paired increase has mixed block directions.

### Costs that remain visible

The sole material flag is isolated P2 at `1e10`, Lehmer cutoff and eight
requested threads (four actual Phi workers): first calls are **6.76% slower
than production**, a paired increase of **0.078916 ms**, in all five blocks.
Median times are 1.246 versus 1.160 ms. This call selects Phi-only, not streaming.
K=1 Phi is also slower at its paired median, but candidate-versus-K=1 changes
are mixed. This is consistent with the earlier scheduling tradeoff; it does
not isolate all of the cost's cause. Warmed calls do not repeat that increase.

Smaller systematic observations are retained too. Held-out isolated Lehmer
P2 at `8.5e10`/six threads is 4.30% / 0.164041 ms slower first than K=1 Phi
in all five blocks. Both choose the same work, and the cause remains unresolved;
it falls below the relative flag threshold. Focused complete first-call
increases against K=1 Phi also occur at `1e10` Meissel/two threads
(2.68% / 0.07025 ms) and `7e10` Lehmer/eight threads (2.31% / 0.099377 ms).
They should not be dismissed merely because the screen does not flag them.

The four short-full-range configurations all improve at their first and warmed
paired medians in this comparison. Their warmed absolute allowance was not
needed. This does not erase the microsecond-scale warmed costs in earlier
short-range experiments.

## Correctness, memory and reproducibility

All ten variants pass the unchanged native checker under both ASan/UBSan and
ThreadSanitizer: 60,322 checks per variant, covering independent small
factorization counts, independent larger pi-table counts, boundaries, production
differential cases, concurrent calls, invalid inputs and complete counter values.
The shared queue passes 174,907 checks per sanitizer. The added independent
selector checker passes 4,674,256 checks per sanitizer, including the exact
fourth-power admission boundary, `L == root` versus `root + 1`, and precedence
of the short-complete-range exception. Five source-drift/overwrite checks pass.
There were no compiler warnings or sanitizer findings.

The timing harness additionally passes 240 independent dense-pi checks.
Complete benchmark expectations come from matching installed primecount Gourdon
and Deleglise–Rivat results. The matched `8e10` isolated counts use the retained
independent sieve/primecount terms; other large isolated expectations are
explicitly production differential comparisons. Passing these tests is evidence
alongside the partition/coverage invariants, not a proof by benchmarking.

All 73 separate allocation/selection probes match the actual kernel partition
to the diagnostic description. Maximum recorded stream/subtotal vector capacity
is **396,344 bytes (about 387 KiB)**; maximum prepared-table vector capacity is
45,316,764 bytes (about 43.2 MiB), at the largest tested input. These are component
capacities, not process peak memory. They exclude shared Phi cache pages,
allocator overhead and thread stacks/state. No new heap allocation was added
by the scheduling change. The existing large-table memory limitation remains.

An [independent audit](results-consolidated/independent-focus-holdout-audit.json)
recomputes all 7,416 summary and 9,036 paired distributions, both flag sets,
and the optimistic best-control comparisons from the raw data. Frozen build
inputs and binaries are unchanged. The separate
[preservation check](results-consolidated/preservation.json) confirms that all
532 inventoried production files and 349 earlier experiment files are unchanged.
The only post-validation generator correction clarified that *diagnostic*
decisions are outside timers; the validation record verifies every compiled
kernel remained byte-identical.

This is a targeted, single-machine Apple arm64/C++11 study, through `1e14`.
Other platforms, full-range peak RSS, R-adapter timings and a full R package
check were not run. Production was not edited or installed. Adoption should
include focused R integration tests after the reviewed native implementation
is moved into the package.

The retained [results directory](results-consolidated/) includes raw compressed
timings, process medians, comparisons, references, memory probes, source hashes,
validation logs, generated sources and the decision record. The code lives in
[consolidated/](consolidated/). To reproduce with a fresh artifact directory,
extract the baseline commit there as `source/`, write that commit ID to
`commit.txt`, then run:

```text
python3 consolidated/bench.py freeze ARTIFACT_ROOT --plan consolidated/plan.json
python3 consolidated/validate.py ARTIFACT_ROOT
python3 consolidated/bench.py build ARTIFACT_ROOT
python3 consolidated/bench.py references ARTIFACT_ROOT
python3 consolidated/bench.py preflight ARTIFACT_ROOT
python3 consolidated/memory.py ARTIFACT_ROOT
python3 consolidated/bench.py timings-focus ARTIFACT_ROOT
python3 consolidated/analyze.py ARTIFACT_ROOT --cohort focus
# Review and record the focused decision without changing the candidate.
python3 consolidated/bench.py timings-holdout ARTIFACT_ROOT
python3 consolidated/analyze.py ARTIFACT_ROOT --cohort both
```

Run from `.github/development/p2-hybrid/`. Keep builds, reference work and
diagnostic workloads separate from timings. The installed primecount executable
is required for the independent complete-count references.

## Recommendation for the first production version

Adopt this single conservative selector together with the shared initial
scheduling change after reviewing the documented first-call cost. Retain the
frozen versions and fixed multipliers as development controls. They do not need
to become production policy branches or public options.

The larger complete-call gains reproduce on new inputs, and the fallbacks avoid
the previously costly small parallel admissions. The tradeoffs include measured
first-call costs and missed opportunities. Expanding admission later
would be a separate, evidence-driven change; this candidate does not need a
larger fitted decision tree to serve as the first version.

In all 18 focused complete-call configurations, v2, v3 and the candidate choose
identical partitions and worker budgets. Their small timing differences cannot
establish a dispatch advantage. The meaningful choice among those policies
concerns the more general isolated-P2 contexts and how much first-call speed
we exchange for conservative warmed behavior. That favors one documented
policy with clear limits over permanently maintaining all three versions.
