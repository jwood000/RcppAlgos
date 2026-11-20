# Optional native regression checks

Run these checks from the repository root when changing GMP conversion,
multiset counting/ranking/unranking, or native primality:

```sh
Rscript --vanilla .github/development/native-regressions/run.R
```

An optional argument selects a **new** output directory. The runner prints its
location and retains source/header snapshots, hashes, the temporary build,
session information, and successful check timings. An existing output path is
rejected so old build objects cannot be reused. Compilation and checks must both
succeed for the command to exit successfully.

The runner uses the existing cpp11 development toolchain and GMP dependency.
It compiles the actual production sources into an isolated R bridge; it does
not load an installed RcppAlgos package or modify package registration. The
temporary build requests C++14 because the existing counting factories use
`std::make_unique`; R's active compiler configuration may select a newer
standard. This does not change the package's compiler settings.
Existing user Makevars settings and library/include environment paths remain
available to the build. No test spawns worker threads.

These checks were moved from `src/test-*.cpp` so normal package installation
does not compile the Catch test framework. Assertions now throw ordinary
exceptions and remain active with `NDEBUG`. The original checks and exhaustive
case counts are retained:

| Fixture | Retained coverage |
| --- | --- |
| `gmp-convert.cpp` | Scalar/vector raw `bigz` conversion of zero, negative values, and signed values beyond `2^53`, with `negPoss = true` |
| `primality.cpp` | Native `IsPrime` at `INT_MAX`, a large prime square, and exceptions beyond its supported range |
| `multiset-gmp-buffers.cpp` | 5,000 seeded cases through four counters, independent multiplicity enumeration, poisoned/reused buffer contents and storage, and `ResetP2D` shape preservation |
| `multiset-gmp-internals.cpp` | 58,032 changing-multiplicity counting checks plus stored-multiplicity checks; 45,659 forced-GMP ranking/unranking cases compared with native-double implementations |

Public R tests continue to check user-facing behavior, independent enumeration,
and large GMP ranks. These fixtures additionally reach internal state and
dispatch paths that public calls do not expose. They are optional local checks:
`.Rbuildignore` excludes `.github`, so CRAN does not build or run them. Run them
explicitly during relevant development and before release; package test success
alone does not establish that these checks passed.
