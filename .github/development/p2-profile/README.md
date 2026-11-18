# Phi cache expansion: retained findings

The current Phi cache accepts arguments through 131,071, twice its original
65,535 endpoint. The maintainer accepted the measured memory and first-call
costs in exchange for faster warmed calls. This is a fixed bounded cache, not
an input-dependent policy or runtime calibration. See
[number-theory.md](../number-theory.md) for the current mathematical structure
and [P2 findings](../p2-hybrid/CONSOLIDATED.md) for the later streaming work.

## What the profiling established

Before hybrid P2, repeated Legendre-style Phi evaluations for outer quotients
were a major cost: at `1e13`, P2 accounted for about 60–67% of complete Meissel
time and 79–88% of Lehmer time in the measured serial/eight-thread cases.
These percentages describe that historical implementation, not today's P2.

Warmed eligible cache accesses already hit. The useful next question was reuse
among ineligible `(x,a)` states, rather than another synchronization change.
At `1e13`, 19.60%/28.70% of Meissel/Lehmer ineligible visits repeated a key;
at `1e12`, only 2.11%/3.70% repeated. A temporary replay predicted that admitting
`x <= 131071, a < 100` could avoid roughly 9.66%/14.34% of evaluator calls in
the sampled `1e13` computations. A full benchmark was still needed because
cross-call reuse and allocation costs were absent from that prediction.

The expansion preserves 100 rows, lazy 4,096-cell pages, atomic answers and
page publication. Tiny Phi handles `a <= 6` before cache access; monotonicity
and independent inclusion–exclusion give
`1 <= Phi(x,a) <= Phi(131071,6) = 25142 < 65536` for cacheable states.
Thus a 16-bit answer cell remains sufficient; zero continues to mean a miss.

## Measured benefit and costs

At `1e13`, complete native warmed calls improved by about 24–28% across all
three methods and one/eight requested threads. Every process pair improved.
Representative medians, in milliseconds:

| Method | Threads | First: original / doubled | Warmed: original / doubled |
| --- | ---: | ---: | ---: |
| Legendre | 1 | 447.02 / 419.14 | 386.07 / 278.26 |
| Meissel | 1 | 445.82 / 421.64 | 384.99 / 278.31 |
| Lehmer | 1 | 456.48 / 430.58 | 394.69 / 286.63 |
| Legendre | 8 | 122.85 / 116.77 | 108.00 / 80.70 |
| Meissel | 8 | 110.11 / 104.01 | 95.02 / 71.14 |
| Lehmer | 8 | 109.36 / 104.04 | 94.30 / 71.44 |

At `1e10`, first calls were consistently 14–20% slower, while warmed calls
were about 10–20% faster. At `1e12`, serial first calls were 2–5% slower;
eight-thread first-call changes were mixed, while warmed calls improved
16–23%. Small-input results did not establish a universal gain.
The repeated `1e14` confirmation improved both first and warmed calls in every
pair; a `1e15` screen had only one pair per method at eight threads and is
weaker evidence. No maximum-range run was made.

At tested inputs `1e12`–`1e14`, observed allocated cache payload grew from
533 pages / 4.1641 MiB to 1,454 pages / 11.3594 MiB. Including the cache object,
the increase was 7.2075 MiB. Maximum configured page payload grew from 12.5 to
25 MiB. These are cache components, **not process peak RSS**: tiny tables,
allocator bookkeeping, thread stacks, and per-call prime tables are excluded.

## Measurement and retained evidence

The main study covered `1e6`, `1e8`, `1e10`, `1e12`, and `1e13`, all three
methods, and thread budgets one/eight. Each configuration had three paired
baseline/candidate fresh processes, one first call and three warmed calls per
process, with randomized paired order. Complete native timings include table
construction/destruction and exclude startup and the R adapter. Timed builds
had no diagnostic instrumentation; builds and sanitizer checks did not overlap
timing. Warmed calls rebuild tables but retain the shared cache.

The larger confirmation used three pairs and one first/two warmed calls at
`1e14`; the `1e15` screen used one pair and one first/one warmed call.
Summaries use per-process medians and paired process ratios. Results are from
macOS arm64 and Apple Clang 17; they do not promise the same crossover or
percentage gain on other systems.

Retained files are the [main raw samples](results/cache-double.csv),
[larger raw samples](results/cache-double-large.csv), corresponding
[main](results/cache-double-summary.csv) and
[larger](results/cache-double-large-summary.csv) summaries,
[cache allocations](results/cache-double-memory.json), and compact
[methodology, source baseline and digests](results/provenance.json).
The current cache/exception/boundary checks are documented in
[testing.md](../testing.md); the old instrumentation is not a supported runner
for current production.

Full profiling traces, replay implementations and original reports remain
recoverable from commit `36303eb268288a2449e348ba6d5b911ffb2fac51`, under
`.github/development/p2-profile/`. The computational baseline for reproducing
that historical experiment is recorded separately in the provenance file.
