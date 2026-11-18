# Automatic prime-counting dispatch

The complete public-R benchmarks support this deliberately simple policy:

| Input | 1–5 usable threads | At least 6 usable threads |
| --- | --- | --- |
| `n < 1e5` | Existing lookup/sieve | Existing lookup/sieve |
| `1e5 <= n < 1e10` | Legendre | Legendre |
| `1e10 <= n < 1e12` | Meissel | Meissel |
| `n >= 1e12` | Meissel | Lehmer |

Usable threads means the validated request capped by the package maximum;
`nThreads = NULL` supplies one. Explicit methods remain selectable. Dispatch
neither times user calls nor inspects cache warmth. These are empirical choices,
not universal crossover points.

## Findings and accepted tradeoffs

Meissel substantially improves the larger serial cases over Legendre. At two
threads, Lehmer was generally 2–5% slower warmed than Meissel; four threads were
close and usually favored Meissel. The follow-up at 4/6/8 threads supported
admitting Lehmer at six. Five and seven threads were not measured.

Representative complete R times, **first / warmed milliseconds**:

| Input | Threads | Legendre | Meissel | Lehmer |
| ---: | ---: | ---: | ---: | ---: |
| `1e10` | 1 | 4.229 / 1.093 | 3.387 / 0.844 | 3.498 / 0.911 |
| `1e12` | 1 | 87.863 / 39.096 | 57.508 / 29.772 | 59.066 / 31.453 |
| `1e12` | 8 | 24.902 / 12.607 | 17.186 / 10.330 | 16.458 / 9.995 |
| `1e14` | 1 | 2532.681 / 2289.653 | 1395.914 / 1205.408 | 1484.597 / 1258.935 |
| `1e14` | 8 | 675.185 / 627.669 | 347.807 / 301.516 | 347.169 / 301.187 |

The final rule accepts real first/warmed tradeoffs. In the follow-up, Meissel
versus Legendre at four threads cost **0.321 ms / 9.32%** on first calls at
`3e10`, while warmed calls improved **14.56%**. At `4e10` the first-call cost
was **0.418 ms / 10.33%**, with a **7.72%** warmed improvement. Every matched
block agreed on these directions. The initial matrix independently found the
`3e10` first-call regression. At `1e11`/eight threads, the initial warmed
change was a small, mixed **+0.0245 ms / +0.96%**. No narrow input exception
was added for these cases.

## Measurement and provenance

The November 11, 2026 study measured **1,980 matching counts from 495 fresh R
processes** on one arm64 macOS 15.7.4 / R 4.6.1 host with ten available threads.
Each method/case/block used one first and three warmed calls. Calls include
the real R wrapper, conversions, tables, native computation and teardown;
process/package loading is excluded. The small sieve at `9973` warmed the
wrapper without filling the Phi cache. Warmed calls rebuilt tables but reused
Phi cache entries. Timing used a test-only monotonic-clock DLL, with JIT
disabled, one preceding garbage collection and no simultaneous heavy work.

The initial 43 cases span `1e8`–`1e14` at 1/2/4/8 threads, with three matched,
randomized blocks. The independent follow-up adds `2e12`, `8e12`, `8e13` at
4/6/8 threads and `2e10`, `3e10`, `4e10` at four threads. Every result matches
both local primecount 9.0 Gourdon and Deleglise–Rivat references. Warmed
comparisons first take each process's median, then paired block changes; those
changes need not equal ratios of the separately reported method medians.
The regression screen required more than 5% **and** 0.05 ms with every block
slower. It was a descriptive screen, not a significance test.

These are historical measurements of the cleaned hybrid P2 before automatic
dispatch changed, not timings of current HEAD. All 310 recorded native/header/R
source hashes match commit `e588cbd3453f1b92e9a29fe2db727d3b686dcf5b`, except
`src/P3.cpp`, which matches `0bd47cd07ae0c747b51f4696b3949994553f1e65`.
Each frozen `suite.json` retains the exact source and DLL hashes. Thus a
historical source reconstruction needs that P3 version as well as the base
commit. Commit `825bf1561b6c814c74e804c13cf88f970dde9ffd` implemented the final
six-thread policy. Later refactors have not rerun this matrix.

The integration's recorded checks passed 867 R expectations in eleven files,
253 sentinel routing/boundary checks, help/examples and explicit C++11 syntax
checking. The nine C++11 warnings also occurred in the previous adapter and
came from installed cpp11 headers. Compact validation summaries and source
hashes remain in `validation/`; the current routing fixture is
[`../prime-count-dispatch-check.R`](../prime-count-dispatch-check.R).
These historical checks do not replace current package validation.

## Retained evidence and replay

Each results directory keeps compressed raw timings, frozen suites, independent
reference values and compact evidence metadata. `results/r-environment.txt`
records the R environment. `SHA256.json` inventories this retained evidence.
Per-process command dumps, duplicate candidate summaries and unsuccessful
harness attempts remain in Git history. The original candidate in `suite.json`
used Lehmer with more than one thread; the final-policy analysis applies six.

Audit the retained matrices without installing or timing anything:

```sh
python3 .github/development/prime-count-auto/refined-analysis.py
# Optionally save the reconstructed detailed comparison to a new file:
python3 .github/development/prime-count-auto/refined-analysis.py --output /tmp/prime-count-auto-analysis.json
```

For a fresh study, use a new output directory and an isolated validated library
whose parent has `package-source-manifest.json` mapping package source paths to
SHA256 hashes. The installed sources must match the checkout. The runner uses
local `/opt/homebrew/bin/primecount`; adjust that path for another host. No tools
are downloaded. Run timings sequentially, after other heavy work finishes:

```sh
python3 .github/development/prime-count-auto/benchmark.py prepare "$auto_root" --library "$validated_root/library"
python3 .github/development/prime-count-auto/benchmark.py timings "$auto_root"
python3 .github/development/prime-count-auto/confirmation.py prepare "$auto_root"
python3 .github/development/prime-count-auto/confirmation.py timings "$auto_root"
python3 .github/development/prime-count-auto/refined-analysis.py "$auto_root" --output "$auto_root/final-analysis.json"
```

Preparation and timings refuse existing outputs and check runner/DLL hashes.
The historical `benchmark.py analyze` command evaluates the initial candidate;
use `refined-analysis.py` for the production rule. Original source and binary
hashes describe provenance; old temporary-library paths need not still exist
for the recorded-data audit.

Peak RSS, full-range counts, cross-platform performance and all cache histories
were not measured. All three methods still prepare a dense square-root pi
table. The evidence supports the stated tradeoffs on this host, without
promising the fastest method for every call.
