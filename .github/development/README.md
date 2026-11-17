# Development notes and checks

This directory keeps the tools needed to check the current implementation and
the evidence behind decisions that are still in use. It is excluded from R
package distributions by `.Rbuildignore`.

| Start here | Purpose |
| --- | --- |
| [testing.md](testing.md) | Installation, focused R checks, native/sanitizer checks, failure injection, and registration |
| [number-theory.md](number-theory.md) | Source map, mathematical and concurrency invariants, memory costs, and decisions worth remembering |
| [DEFERRED_ISSUES.md](DEFERRED_ISSUES.md) | Outstanding investigations and their recorded disposition |

## Retained studies

| Study | Why it remains useful |
| --- | --- |
| [Phi cache and P2 profile](p2-profile/README.md) | Why repeated Phi work was expensive and why the cache argument limit was doubled |
| [P2 scheduling and streaming](p2-hybrid/CONSOLIDATED.md) | Current selector, first-query imbalance, final comparison, and known regressions |
| [Automatic method selection](prime-count-auto/README.md) | Evidence for the current size/thread dispatch rule |
| [Prime-sieve benchmarks](prime-sieve-benchmarks/README.md) | Sieve optimization measurements and their reproduction workflow |
| [nthPrime and count-only sieve](number-theory.md#sieving-and-nthprime-measurements) | Final native-sieve comparison, first/warmed calls, and memory tradeoffs |

The [sieve inspector](prime-sieve-inspector/README.md) is optional instrumentation
for investigating marking and extraction. [cache-inspector.md](cache-inspector.md)
describes its isolated-package workflow. Neither belongs in production builds.

The top-level `*-check.cpp`, `*-check.R`, and `prepare-*-check.py` files are
current correctness tools. Their build dependencies and uses are listed in
[testing.md](testing.md); keep them distinct from historical timing artifacts.

## Historical experiments

The November 2026 consolidation removed superseded selector variants, generated
source copies, repeated summaries, per-process command dumps, empty build logs,
and unsuccessful run directories. The findings, meaningful regressions, current
checks, and selected final raw measurements remain here. Archived observations
are measurements of particular source versions, not promises about every host
or the current checkout.

The complete pre-consolidation directory is in commit
`36303eb268288a2449e348ba6d5b911ffb2fac51`. To inspect an old report without
changing the working tree:

```sh
git show 36303eb:.github/development/p2-hybrid/SELECTOR_V3.md
```

To recover all original studies into a separate directory:

```sh
history_root=$(mktemp -d /tmp/rcppalgos-study-history.XXXXXX)
git archive 36303eb .github/development | tar -x -C "$history_root"
```

Earlier nthPrime calibration experiments were removed before that commit. Their
separate archive and its recorded checksum are described in
[number-theory.md](number-theory.md#decisions-not-to-repeat-without-new-evidence).

For future investigations, retain the question, source version, method, final
raw measurements, conclusion, and material limitations. Prefer updating the
relevant study to adding another parallel account. Disposable build products,
intermediate derived tables, and full process logs can stay outside the repo
unless they are needed to reproduce a failure.
