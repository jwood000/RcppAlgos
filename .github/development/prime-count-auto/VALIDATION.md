# Automatic dispatch integration — November 11, 2026

The [matched benchmark](README.md) supports keeping the proposed input
boundaries while raising Lehmer's usable-thread requirement to six. The
implementation in `src/PrimeCount.cpp` selects:

- Existing lookup/sieve for `n < 1e5`.
- Legendre for `1e5 <= n < 1e10`.
- Meissel for larger inputs, except Lehmer when `n >= 1e12` and the usable
  thread budget is at least six.

The usable budget is `min(nThreads, maxThreads)` after the existing conversion
and validation; `NULL` still supplies one thread. The original request and
maximum are forwarded unchanged to the selected native driver. Explicit
methods, small-input thread handling, validation order, result conversion,
public formals and registration are unchanged. No numerical kernel, cache,
worker policy, compiler requirement or dependency changed in this integration.

Selection only changes which existing counting identity evaluates the same
validated input. Actual routing is tested independently of counting correctness:
equal answers alone would not reveal a dispatcher that always chose Legendre.

## Checks

A clean copy was installed in a new temporary library, and the test runner
verified the loaded DLL's location. All **867 expectations** passed across
eleven number-theory/error test files, with no failures, errors, warnings or
skips. The two original cache/batch tests now explicitly select Legendre so
they retain their intended coverage. Added automatic-boundary count/type
checks exercise both sides of `1e10` and `1e12` with 1, 2 and 6 threads on
this host.

The sentinel fixture passed **253 checks** through the real R wrapper,
converters, native adapter and cpp11 boundary. It verifies exact size
boundaries, five versus six usable threads, requests capped to 1/2/5/6,
unchanged raw argument forwarding, explicit overrides, aliases, maximum
exact input, validation, integer/double conversion and exception recovery.
Sentinel maximum-input checks do not allocate or execute a full-size count.

`primeCount.Rd` passed parsing and checking, and its installed examples ran.
The public formals remain `primeCount(n, nThreads = NULL, method = "auto")`;
the native registration still has four arguments. The existing partial-sieve
runner also checked `semiPrimeCount.Rd` and its examples. `git diff --check`
passed.

The clean package and sentinel builds had no compiler warnings. R 4.6.1 on
arm64 macOS 15.7.4 selected its existing `gnu++20` configuration, including
for the cpp11 fixture that requested CXX11; no build setting changed. A
separate explicit `clang++ -std=c++11 -fsyntax-only -Wall -Wextra -Wpedantic`
check succeeded. It reported nine warnings from installed cpp11 headers;
compiling the previous adapter with identical flags produced the same nine
warning locations/messages. Both outputs are retained.

An auxiliary formals check initially compared R's pairlist of formals against
an ordinary list and failed. The expected value was corrected to a pairlist;
the public function was unchanged, and the help/formals/registration check
then passed. The initial script and log are retained with the final check.

## Artifacts and limits

[validation/](validation/) retains installation, R-test, help, sentinel and
C++11 logs, commands, source hashes and the preservation summary. The tested
adapter SHA256 is
`b43945058d577eca9a4ba3b097f6ae90465b09a5829174a119e27a64d5080656`.
Its package snapshot matches the final working source. All counting kernels
retain their pre-dispatch hashes, including P2's
`b6553ed680d2445c1230d0b1f5e9522646a1f40643706c33a55b0bade50e2a8b`.

The temporary root is recorded by `/tmp/rcppalgos-prime-count-auto-path`:
`/var/folders/fd/v2_r1ssx14jd6yww77jfvnmm0000gn/T/rcppalgos-prime-count-auto-emmr8i33`.
For the clean installation use the established `R CMD INSTALL --preclean
--no-multiarch --library=...` workflow. The focused runner is
`../p2-hybrid/production-integration/r-tests.R`; the routing fixture is
`../prime-count-dispatch-check.R`. Run the latter from the repository root
with a fresh temporary output directory. Its counting kernels are sentinels,
so it verifies dispatch and the R/C++ boundary rather than worker execution.

The Git index and previously modified files outside this task were preserved;
no commits were made. The normal installed package and working-tree build
products were not replaced.

The benchmark contains 1,980 verified counts from 495 fresh R processes and
documents the first-call costs around `3e10`–`4e10` with four threads. It does
not establish a universally fastest policy. Full `R CMD check`, other-platform
testing, peak-RSS measurement and counts at the full numeric ceiling were not
run. No new sanitizer builds were needed for this adapter-only change: the
unchanged kernels' recent ASan/UBSan and TSan evidence remains recorded in the
[P2 cleanup validation](../p2-hybrid/production-integration/CLEANUP.md).
