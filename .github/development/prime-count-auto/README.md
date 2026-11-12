# Public R method benchmark for automatic prime counting

This study compares the complete public `primeCount()` calls using the cleaned
hybrid P2 implementation. It uses the isolated installation validated in
[the P2 cleanup](../p2-hybrid/production-integration/CLEANUP.md), before changing
automatic dispatch. All explicit counting methods, native sources, headers and
R wrappers matched that installation's source manifest at preparation time.

The final dispatcher is implemented and its fresh package, routing and help
checks are recorded in [VALIDATION.md](VALIDATION.md).

## Initial candidate and measurement

The initial rule was declared before timing:

- Keep the existing lookup/sieve below `1e5` and Legendre below `1e10`.
- Use Meissel from `1e10`, except choose Lehmer from `1e12` when more than one
  thread is usable.

The matrix has 43 input/thread cases: `1e8`, `1e9`, `3e9`, `1e10`, `3e10`,
`1e11`, `3e11`, `1e12`, `3e12`, and `1e13`, each at 1, 2, 4, and 8 threads;
plus `1e14` at 1, 4, and 8 threads. Three matched blocks randomize case order
and the three method orders within each case. Each method gets a fresh R
process with one first and three warmed calls. All **1,548 counts from 387
processes** matched independent local primecount 9.0 Gourdon and
Deleglise–Rivat results.

Unlike the earlier native phase harness, this measures the real installed R
wrapper, its argument conversion and dispatch, and the complete native call,
including table construction/destruction. Process/package loading is outside
the measured interval. Before counting, it warms only the small `auto` sieve
path (`n = 9973`), which does not evaluate Phi; the first measured call therefore
starts with a fresh Phi cache. The installed R code is byte-compiled, JIT is
disabled for the measurement sequence, and one garbage collection runs before
the sequence. Warm calls rebuild tables but retain Phi cache contents.

The timer is a test-only `std::chrono::steady_clock` DLL with a per-process
origin, called from R before and after the count. No timing code enters the
package, and no timer overhead is subtracted. No builds, sanitizers or other
CPU-heavy diagnostics ran alongside the timed matrix. The host's normal R
thread ceiling was 10, so all requested budgets were available.

An initial attempt with `microbenchmark::get_nanotime()` stopped on its
uptime/overflow precision warning after the first process; its four correct
counts and warning are retained separately and excluded from the declared
completed matrix. The full matrix restarted with the monotonic timer. A timer
build also initially stopped because `R CMD SHLIB` writes compiler/SDK
identification to stderr despite a clean exit; the runner now distinguishes
that informational output from compiler diagnostics. Both initial attempts
are preserved under `results/`; neither changed package code or the case order.

## Initial findings

The initial candidate switches methods in 31 cases. Relative to Legendre,
30/31 switched cases improve in every block for first calls, and 30/31 improve
in every block for warmed calls. The exceptions are different cases:

- At `3e10` with four threads, first-call Meissel is slower in all three
  blocks: median paired **+0.396 ms / +11.19%**. Warmed Meissel improves in
  every block: **-0.130 ms / -9.32%**.
- At `1e11` with eight threads, warmed Meissel has a small mixed change:
  **+0.0245 ms / +0.96%**. This is below the declared regression screen.

The screen requires an increase exceeding both 5% and 0.05 ms, with every
block slower. It identifies the one first-call case above. All raw observations
remain recorded; the screen is not a significance test or permission to ignore
smaller differences.

Representative complete R times are shown as **first / warmed milliseconds**:

| Input | Threads | Legendre | Meissel | Lehmer |
| ---: | ---: | ---: | ---: | ---: |
| `1e10` | 1 | 4.229 / 1.093 | 3.387 / 0.844 | 3.498 / 0.911 |
| `1e10` | 8 | 3.290 / 0.650 | 1.831 / 0.607 | 1.885 / 0.653 |
| `3e10` | 4 | 3.487 / 1.396 | 3.811 / 1.263 | 3.378 / 1.321 |
| `1e12` | 1 | 87.863 / 39.096 | 57.508 / 29.772 | 59.066 / 31.453 |
| `1e12` | 2 | 49.756 / 23.351 | 32.855 / 18.042 | 33.992 / 18.734 |
| `1e12` | 4 | 28.249 / 13.738 | 19.062 / 11.598 | 19.677 / 11.863 |
| `1e12` | 8 | 24.902 / 12.607 | 17.186 / 10.330 | 16.458 / 9.995 |
| `1e14` | 1 | 2532.681 / 2289.653 | 1395.914 / 1205.408 | 1484.597 / 1258.935 |
| `1e14` | 4 | 735.381 / 693.665 | 415.240 / 365.240 | 420.461 / 372.239 |
| `1e14` | 8 | 675.185 / 627.669 | 347.807 / 301.516 | 347.169 / 301.187 |

The method distinction needs a narrower interpretation than “parallel means
Lehmer.” Meissel consistently wins the serial larger cases. At two threads,
Lehmer is generally 2–5% slower warmed; at four threads, the methods are close,
usually favoring Meissel. Eight threads generally favor Lehmer by a few percent.
For example, at `1e12`, the paired warmed Lehmer changes versus Meissel are
+4.94%, +1.83%, and -3.24% at 2, 4, and 8 threads respectively. This initial
matrix does not measure the crossover at 5–7 threads.

## Focused follow-up

After reviewing the complete initial matrix, a separate confirmation matrix
was declared: `2e12`, `8e12`, and `8e13` at 4, 6, and 8 threads, plus `2e10`,
`3e10`, and `4e10` at four threads. It uses the same first/warmed methodology,
three matched blocks and a new seed. The new large inputs bracket the usable
thread decision; the smaller inputs independently repeat and contextualize the
first-call regression. Initial observations and analyses are not overwritten.

All **432 counts from 108 fresh processes** matched both independent reference
algorithms. The follow-up supports a more conservative Lehmer admission at six
usable threads, rather than the initial candidate's “more than one.” Paired
Lehmer changes versus Meissel are:

| Input | Threads | First change | Warmed change |
| ---: | ---: | ---: | ---: |
| `2e12` | 4 | +3.40% | +2.49% |
| `2e12` | 6 | -7.55% | -4.25% |
| `2e12` | 8 | -5.12% | -1.88% |
| `8e12` | 4 | +2.21% | +2.30% |
| `8e12` | 6 | -3.71% | -4.45% |
| `8e12` | 8 | -1.83% | -2.51% |
| `8e13` | 4 | -0.08% | +0.86% |
| `8e13` | 6 | -7.11% | -6.65% |
| `8e13` | 8 | -2.16% | -3.17% |

Every first-call block at six and eight threads favors Lehmer. Every warmed
block also favors it except `8e12`/six threads, whose directions are mixed.
At four threads Meissel wins consistently for `2e12` and `8e12`; `8e13` is
mixed. Five and seven threads were not measured, and these observations on one
host are not proof of an exact architectural crossover.

The cold-call cost is independently confirmed: Meissel at `3e10`/four threads
is **+0.321 ms / +9.32%** first and **-14.56%** warmed versus Legendre. At
`4e10`/four threads it is **+0.418 ms / +10.33%** first and **-7.72%** warmed.
Each of those directions holds in every block. At `2e10`/four threads Meissel
improves both first and warmed calls. This is a genuine first/warmed tradeoff
in the simple policy, not merely a single anomalous sample. The final dispatch
decision retains that limitation explicitly without introducing a special
input-specific exception around these observations.

[results-confirmation/](results-confirmation/) preserves this matrix separately.
Together the two completed matrices contain **1,980 matching counts from 495
fresh R processes**; `3e10`/four threads is intentionally measured in both.

## Final policy

The final rule keeps the proposed size boundaries and raises Lehmer's usable
thread requirement from more than one to at least six:

| Input | 1–5 usable threads | At least 6 usable threads |
| --- | --- | --- |
| `n < 1e5` | Existing lookup/sieve | Existing lookup/sieve |
| `1e5 <= n < 1e10` | Legendre | Legendre |
| `1e10 <= n < 1e12` | Meissel | Meissel |
| `n >= 1e12` | Meissel | Lehmer |

Usable threads means the validated request capped by the package maximum.
`nThreads = NULL` continues to mean one thread. The rule does not inspect cache
warmth, perform timing during a user's call, or change explicit method choices.
It adopts broad improvements over the historical policy while accepting the
documented 0.3–0.4 ms first-call costs at `3e10`–`4e10`/four threads. It does not
claim to select the fastest method for every input, platform or cache state.

[results/refined-analysis.json](results/refined-analysis.json) separately
evaluates this final rule from both raw matrices. Its independent reconstruction
verifies every expected count, input/thread combination, repeated-call index,
process identity, complete case matrix and timing-data hash. The initial
candidate analyses remain unchanged. The same first-call regression records
remain: one in the initial matrix and two in the follow-up (with `3e10`
intentionally repeated). This is a documented tradeoff, not an unreported
failure of the screen.

## Reading the results and reproducing

[results/](results/) contains the frozen suites, independent references, DLL
and source hashes, raw compressed timings, command metadata and summaries.
Warmed summaries first take the median of three calls within each process;
method comparisons then use paired block differences or ratios. These paired
changes need not equal differences or ratios of separately reported medians.
The comparison file includes both candidate-specific and fixed method pairs;
those roles can produce duplicate identical comparison rows.

Prepare a new artifact directory and point to an isolated validated library,
whose parent contains `package-source-manifest.json`. Its native/header/R
sources must match the current checkout. The original run used the cleanup
installation; after changing the adapter, use the new validated installation
for a new study. The old results and hashes still describe the original run:

```sh
python3 .github/development/prime-count-auto/benchmark.py prepare "$auto_root" \
    --library "$validated_root/library"
# Run timings only after all other CPU-heavy work has finished.
python3 .github/development/prime-count-auto/benchmark.py timings "$auto_root"
python3 .github/development/prime-count-auto/benchmark.py analyze "$auto_root"

python3 .github/development/prime-count-auto/confirmation.py prepare "$auto_root"
python3 .github/development/prime-count-auto/confirmation.py timings "$auto_root"
python3 .github/development/prime-count-auto/confirmation.py analyze "$auto_root"
python3 .github/development/prime-count-auto/refined-analysis.py "$auto_root"
```

The scripts refuse to overwrite existing measurement directories or timing
files and verify frozen runner/DLL hashes. References use the local
`/opt/homebrew/bin/primecount`; no package or tool was installed. The temporary
directory is recorded in `/tmp/rcppalgos-prime-count-auto-path`.

These are finite samples on one host with at most three matched blocks per
case. They compare first and repeated same-input calls separately; they do not
establish universal crossover points, behavior for every cache history, or
performance through the full supported numeric ceiling. Peak RSS was not
measured, and all methods still prepare the dense square-root pi table.
