# P2 cleanup validation — November 11, 2026

This follow-up preserves the maintainer's formatting and signed-integer changes
and applies the agreed small cleanup. The selector thresholds, guard order,
quotient boundaries, phase barrier and R interfaces remain unchanged.

## Source changes

- `countPhiPrefix()` constructs its queue before choosing workers and caps
  launches by available jobs. A 64-query prefix has 33 jobs: 32 individual
  queries and one remaining 32-query job. An inherited budget of 36 now
  launches 33 workers. Budgets through 32 are unaffected by this cap.
- `streamWorkerBudget()` supplies the same formula to the selector and the
  streaming launcher. It retains the existing numeric-span grain and clamps.
- The stream class comment collects its input, seed, table-coverage and target
  ordering preconditions. The selector comment points to the retained
  consolidated benchmark rationale.

A comparison against the maintainer's staged source confirms that the entire
selector is identical after substituting the shared helper call for the
duplicated formula and removing comments/whitespace. The inherited Phi budget
is unchanged; only the actual number of prefix workers is additionally capped.

The signed queue and stream indices are covered by fresh validation. The tested
P2 SHA256 is
`b6553ed680d2445c1230d0b1f5e9522646a1f40643706c33a55b0bade50e2a8b`.
Both R and native snapshots match that source.

## Harness changes

`helpers.cpp` now uses signed tickets, including negative values and both signed
extremes. Comparisons with the historical unsigned queue explicitly convert
nonnegative sizes. The tests retain exhaustive queue lengths 0 through 130,
larger boundaries, concurrent coverage and intervals ending at `INT64_MAX`.

A test-only observer wraps the real `ThreadBatch`; its counters are accessed
only by the calling thread. It verifies actual launches and counts without
adding production hooks. Besides empty/singleton/chunk boundaries, it checks
the real selector geometry `x = 1e12`, `L = 1e6`, `a = 3304`, requested budget
64: the inherited Phi budget is 36 and its shortened 64-query prefix launches
33 workers. Shared stream-budget tests cover grain transitions, nonpositive
thread requests and large arithmetic inputs. Frozen-selector comparisons now
include budgets above 32.

The validator also retargets the retained independent private-stream oracle to
current production, checking repeated and odd/even endpoints, segment crossings,
large finishes and counts crossing `INT32_MAX`. This oracle deliberately uses
a synthetic seed to check increments; other tests verify absolute counts.

The source audit still checks three unchanged counting helpers against the
approved candidate. It explicitly records three intentionally revised bodies
and their behavioral coverage instead of claiming their tokens are identical.
Historical sources and results under `results/` were not overwritten.

## Validation results

A clean package copy was installed into a new temporary library. Its loaded
DLL was checked before running tests. All **864 expectations** passed across
eleven number-theory/error test files, with zero failures, errors, warnings or
skips. `semiPrimeCount.Rd` parsing/checking and installed examples also passed.
The package build had no compiler warnings. This used R 4.6.1 on arm64 macOS
15.7.4, with the existing R compiler configuration unchanged.

Both **ASan/UBSan and ThreadSanitizer** passed the following checks independently:

| Coverage | Checks per sanitizer |
| --- | ---: |
| Common independent/differential counting suite | 60,322 |
| Signed queue, concurrency and extreme intervals | 71,120 |
| Selector comparison with frozen candidate | 15,013 |
| Stream budgets and actual Phi launches | 739 |
| Independent private stream windows | 88 |

The 59,753 small exhaustive counts are included in the common suite's total.
Native builds explicitly used C++11 with `-Wall -Wextra -Wpedantic`; there were
no compiler warnings or sanitizer findings. `git diff --check` also passed.

## Focused worker-cap timing

After all builds and sanitizer tests finished, a small comparison used the
maintainer's staged source as the before version and this cleanup as the after
version. It timed only the 64-query Phi prefix above, excluding table setup.
Five randomized matched blocks each used a fresh process with one first call
and three warmed calls. All 80 counts from 20 processes matched.

| Budget | State | Before median (ms) | After median (ms) | Median paired change (ms) |
| ---: | --- | ---: | ---: | ---: |
| 8 | First | 0.533125 | 0.569625 | +0.046208 |
| 8 | Warmed | 0.109250 | 0.110791 | +0.002833 |
| 36 | First | 0.702583 | 0.623709 | -0.166958 |
| 36 | Warmed | 0.228166 | 0.216083 | -0.027958 |

Budget 8 is the unchanged-launch control; both of its comparisons varied in
direction across blocks. Budget 36 changes 36 launches to 33: its first-call
direction also varied, while warmed calls improved in every block (median
paired change -12.58%). These short measurements support removing unnecessary
launches in this case, not a general scaling claim or a new selector fit. The
normal R thread ceiling here is 10; this native private-helper diagnostic
deliberately tests a larger budget. Paired changes need not equal differences
between the independently reported medians.

## Reproduction and limits

[cleanup-results/](cleanup-results/) retains installation/R logs, source hashes,
the selector comparison, native commands/build/run logs, and raw focused timing
data with its reproduction script. The [R runner](r-tests.R) and updated
[native validator](validate.py) use the same isolated-directory preparation
described in [README.md](README.md). Use a new artifact directory: validation
refuses to overwrite an existing `native-validation/` directory. The focused
timing script also expects the maintainer's before source at
`before/staged-P2.cpp`, captured before applying the cleanup.

The temporary directory is recorded in `/tmp/rcppalgos-p2-cleanup-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-p2-cleanup-ntlp0qu5`.
The installed package in the normal library, working-tree build products and
Git index were not replaced. No commits were made.

Full `R CMD check`, other-platform testing, process peak-RSS measurement and a
new selector benchmark sweep were not run. This validation covers the cleanup
and reviewed integer changes; the previously accepted selector tradeoffs remain.
