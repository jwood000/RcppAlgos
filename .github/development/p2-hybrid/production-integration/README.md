# Production integration of the consolidated P2 policy

This page records the initial integration and its tested source snapshot.
The subsequent signed-index and worker-budget cleanup, with fresh R and
sanitizer validation, is recorded in [CLEANUP.md](CLEANUP.md). Historical
results under `results/` remain unchanged.

The maintainer accepted the tradeoffs in the
[consolidated comparison](../CONSOLIDATED.md). `src/P2.cpp` now implements that
policy. The historical experimental sources and results remain unchanged.

The implementation keeps the direct-table suffix, evaluates the expensive
prefix with shared Phi, and optionally counts a bounded lower quotient range
with a segmented sieve. Both Phi-only calls and hybrid prefixes use individual
jobs for the first 32 queries, then 32-query jobs from the same queue. The
reviewed cap schedule, guard order, worker budgets, and barrier between Phi
and streaming are preserved.

All helpers are private to `P2.cpp`. Production has one selector and a fixed
queue, with no experimental flags, diagnostic counters, alternative versions,
or new public headers. The chooser returns the first streamed index instead
of the experimental diagnostic structure. One constant supplies the stream
worker grain to both the chooser and launcher. No R signatures, registration,
default method, Phi cache, Legendre implementation, or language requirements
changed. The partial-sieve header and `semiPrimeCount` help now describe the
additional streaming buffers.

## R integration

A clean copy of the package was installed into a temporary library with
`R CMD INSTALL --preclean --no-multiarch`. The normal installed package and
working-tree object files were not replaced. The runner checks that its loaded
DLL is under the temporary library before testing.

All **864 expectations** passed across eleven files, with zero failures,
errors, warnings or skips:

| Coverage | Expectations |
| --- | ---: |
| `testPrimeCount.R` | 119 |
| `testPartialSieve.R` | 320 |
| `testErrors.R` | 207 |
| Eight other number-theory files | 218 |

The other files cover prime sieving, factorization and its sieve, primality,
Euler phi, divisors and their sieves. Prime-count tests run first in a fresh R
process to retain their existing cold-cache checks. Added tests cover:

- Short tails around a prime square, independently counted with a base-R pi
  table, plus singleton and empty cutoffs.
- Independent retained `8e10` semiprime counts at zero and Meissel cutoffs,
  serial/parallel evaluation, and integer/double result types.
- Known `1e12` Meissel and Lehmer counts with serial and parallel requests.

The changed Rd file passed parsing and checking, and its installed examples
ran successfully. The package build produced no warnings. This used R 4.6.1
on arm64 macOS 15.7.4; R's existing configuration selected GNU C++20. Separate
native checks explicitly used C++11. No build configuration was changed.

## Native preservation and concurrency

Both ASan/UBSan and ThreadSanitizer passed, each with:

- 60,322 unchanged common checks, including independent factor/prime-count
  oracles, larger differential cases, concurrent P2 calls, validation, and
  complete Meissel/Lehmer counts.
- 70,058 queue checks covering equivalence to the approved K=1 mapping,
  disjoint coverage, concurrency, overflow boundaries, and absence of an
  intermediate barrier.
- 9,313 selector checks comparing the production chooser with the frozen
  candidate at representative cap/admission boundaries and table/cutoff contexts.

Six preprocessed counting/stream helper bodies match the approved source after
removing whitespace, renaming the queue, and substituting the identical shared
worker-grain constant. The source audit separately verifies that constant's
value. The simplified queue and chooser are checked at runtime rather than
asserted to be textually identical.

One initial helper-test build warned about test variables that the compiler
could not prove initialized after a throwing check. Those test variables were
initialized explicitly; the original warning is retained. Production had no
warning, and the final helper builds are clean. The already-passing common
ASan suite was retained because neither its inputs nor production changed.
A later diagnostic-only audit assertion checked the shared constant before
normalization; no compiled input changed.

## Performance preservation

A fixed twelve-case replay compared the integrated production kernel directly
with the approved candidate, after all builds and other tests finished. It
includes early-cutoff Phi work, an enlarged-table fallback, the known small
first-call case, a small parallel admission, a short full stream, the 32L
promotion, and complete serial/parallel counters through `1e14`.

Three randomized matched blocks used one first and three warmed calls per
fresh process: **288 results from 72 processes**, all matching the retained
expectations. There were **zero primary regression flags** under the existing
screen (>5% and >0.05 ms paired median increase, every block slower).
Complete-call first changes ranged from -2.15% to +0.32%, and warmed changes
from -1.46% to -0.04%, relative to the approved candidate. At `1e14` with eight
threads, warmed Meissel/Lehmer totals differed by -0.15%/-0.57%.

Smaller systematic costs remain recorded: doubled-table Meissel-cutoff P2 at
`1e10` was +4.66% / +0.051251 ms first; Meissel-cutoff P2 at `1e12` was +1.03%
/ +0.053624 ms first and +1.78% / +0.062584 ms warmed; complete serial Lehmer
at `1e12` was +0.32% / +0.189835 ms first. These observations do not establish
their cause or statistical equivalence. No policy was retuned and the earlier
accepted first-call/fallback tradeoffs still apply.

Both timing builds describe their work after timing with the frozen selector.
The separate production-helper comparisons verify the production chooser;
identical diagnostic rows alone are not evidence of dispatch correctness.
Percentages and absolute changes are medians of paired block comparisons;
warmed values use the median of three calls within each process. A separate
audit recomputed all 144 paired distributions from the raw data.

## Artifacts and reproduction

[results/](results/) retains installation/R logs, test summaries, source hashes,
native validation logs, the initial harness warning, replay builds/commands,
raw compressed timings, and analysis. Binaries, object files and full temporary
checkouts remain outside the repository. The native comparison baseline is
commit `8b43ed790435b22eb9ec02627379082fe7deb103`.

Prepare a fresh artifact directory with that commit extracted as `baseline/`.
Copy current package sources into `package-source/`, excluding old object/shared
library files, and create `library/`. From the repository root:

```sh
MAKEFLAGS=-j4 R CMD INSTALL --preclean --no-multiarch \
    --library="$integration_root/library" "$integration_root/package-source"
Rscript --vanilla .github/development/p2-hybrid/production-integration/r-tests.R \
    "$integration_root" "$PWD"
python3 .github/development/p2-hybrid/production-integration/validate.py \
    "$integration_root"
python3 .github/development/p2-hybrid/production-integration/replay.py \
    prepare "$integration_root"
# Finish all other builds/tests before timing.
python3 .github/development/p2-hybrid/production-integration/replay.py \
    timings "$integration_root"
python3 .github/development/p2-hybrid/production-integration/replay.py \
    analyze "$integration_root"
```

The native validator and replay refuse to overwrite existing artifact
directories. Set `integration_root` to the prepared directory before running
these commands. The initial integration P2 hash is recorded in the audit and
build metadata. Its R and native snapshots matched the working source when
that validation ran.

All 780 inventoried tracked files outside the authorized change set and all
495 earlier experiment files are unchanged. No staging or commits were done.
Full `R CMD check`, other-platform runs, peak-RSS measurement and a new broad
policy benchmark were not run; this task integrates and checks the already
reviewed policy.
