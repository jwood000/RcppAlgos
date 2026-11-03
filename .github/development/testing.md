# Testing workflows

Run these commands from the repository root. Use a temporary library so an
investigation does not replace the normally installed package.

## Package installation and focused tests

The prime-count cache fix was built from a temporary copy of the package source,
without pre-existing object files. `package_source` below is that copy; using
the working tree instead will update its ignored build artifacts.

```sh
package_library=$(mktemp -d)
MAKEFLAGS=-j4 R CMD INSTALL --preclean --no-multiarch \
    --library="$package_library" "$package_source"
Rscript --vanilla - "$package_library" <<'RS'
library(RcppAlgos, lib.loc = commandArgs(TRUE)[1])
testthat::test_file("tests/testthat/testPrimeCount.R",
                    reporter = "summary", stop_on_failure = TRUE)
testthat::test_file("tests/testthat/testErrors.R",
                    reporter = "summary", stop_on_failure = TRUE)
RS
```

The parallel cases precede the serial cases in `testPrimeCount.R`. A standalone
run in a new R process therefore tests an initially empty phi cache. Running the
whole suite may populate it earlier. Count comparisons do not detect data races;
use ThreadSanitizer separately.

## Native prime-count sanitizer checks

This checks the production counting core and `Eratosthenes.cpp`, without the R
adapter. The extraction replaces only the R conversion include with the numeric
constants header and omits the registered R entry point. The cache, recursion,
and scheduling code are unchanged. The standalone checks call the internal
entry point directly; use supported inputs and realistic thread counts.

```sh
prime_count_check=$(mktemp -d)
python3 - "$prime_count_check" <<'PY'
from pathlib import Path
import sys
source = Path("src/PrimeCount.cpp").read_text()
assert source.count("[[cpp11::register]]") == 1
assert source.count('#include "CppConvert.h"') == 1
core = source.split("[[cpp11::register]]")[0].replace(
    '#include "CppConvert.h"', '#include "CppConvert/Constants.h"')
(Path(sys.argv[1]) / "prime-count-core.cpp").write_text(core)
PY
clang++ -std=c++11 -O1 -g -pthread -fsanitize=thread \
    -fno-omit-frame-pointer -Iinst/include -I"$prime_count_check" \
    .github/development/prime-count-cache-check.cpp src/Eratosthenes.cpp \
    -o "$prime_count_check/prime-count-tsan"
TSAN_OPTIONS=halt_on_error=1 \
    "$prime_count_check/prime-count-tsan" 4
TSAN_OPTIONS=halt_on_error=1 \
    "$prime_count_check/prime-count-tsan" 4 --allocation-failure
```

The checker includes the extracted core in its translation unit, so it can
inspect every populated cache entry without exposing production test hooks.
Its optional second argument replaces page allocations with failures in the
standalone executable. The package itself does not replace allocation functions.

The default mode runs known counts in increasing/decreasing order, compares
small counts and prime-square neighbors with an independent Eratosthenes sieve,
and compares every populated phi cell with independent partial-sieve counts.
The failure mode checks the result with no page allocations, verifies that no
page was published, then checks recovery, signs, zero entries, page boundaries,
and the final valid cache key.

AddressSanitizer and UndefinedBehaviorSanitizer use a separate executable:

```sh
clang++ -std=c++11 -O1 -g -pthread -fsanitize=address,undefined \
    -fno-omit-frame-pointer -Iinst/include -I"$prime_count_check" \
    .github/development/prime-count-cache-check.cpp src/Eratosthenes.cpp \
    -o "$prime_count_check/prime-count-asan-ubsan"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$prime_count_check/prime-count-asan-ubsan" 4
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    "$prime_count_check/prime-count-asan-ubsan" 4 --allocation-failure
```

The original implementation produced a ThreadSanitizer data-race report during
vector growth on a cold `10^10`, four-thread call. On the macOS toolchain that
report included
external-symbolizer warnings; `TSAN_OPTIONS=halt_on_error=1:symbolize=0` disables
symbolization without disabling race detection.

The integrated atomic cache passed both modes with ThreadSanitizer at 1, 2, 3,
4, and 10 threads, and with AddressSanitizer/UndefinedBehaviorSanitizer at four
threads. Each default-mode run checked 662 independent prime counts and
1,394,803 populated phi cells. Each failure-mode run rejected 646 page
allocations while retaining the correct result. There were no sanitizer reports.

## Atomic-cache contract

The shared cache preserves the original bounds (`a < 100`, `x <= 65535`) and
reuse across worker batches and calls. Each row owns up to sixteen fixed blocks
of 4,096 atomic cells, allocated on demand. Cells and pointers are explicitly
initialized for C++11. Published blocks remain alive and stationary until cache
destruction, which requires all workers to have joined.

The row mutex protects allocation only. Release/acquire publication makes block
initialization visible to readers. Cell loads/stores are relaxed: zero is a
miss, and each positive value is the complete phi answer. Duplicate computations
store the same value; a missed cache value only causes recomputation. Failed
page allocation skips the cache write and leaves the calculated answer intact.

This does not make the whole native entry point reentrant: `phiPrimes` and
`phiPi` remain shared mutable setup state. The recurrence, scheduling, and R
interface are unchanged. Separate scheduling findings from the initial audit
remain outside this cache fix.

On the tested arm64 macOS environment, both pointer and 16-bit atomics were
lock-free and a cell occupied two bytes. Correctness does not require lock-free
atomics; their size and performance may differ elsewhere.

## Native timing checks

Use the extracted core from above and the separate timing driver:

```sh
clang++ -std=c++11 -O3 -pthread -Iinst/include \
    .github/development/prime-count-check.cpp \
    "$prime_count_check/prime-count-core.cpp" src/Eratosthenes.cpp \
    -o "$prime_count_check/prime-count-timing"
"$prime_count_check/prime-count-timing" 1000000000000 10 4 37607912018
```

Arguments are `n threads repetitions expected_count`. Each invocation starts a
fresh process; repetition zero is cold, and later repetitions reuse the cache.
Output columns are `n,threads,repetition,count,milliseconds`. Timings include
native setup and cache allocation, but exclude R and process startup. Each
result is checked against the expected count. Sanitizer timings must not be
used as performance results. Keep before/after binaries, alternate their
execution order, and compare cold and warm calls separately.

The November 2, 2026 integration check used seven fresh processes per version
and thread count, each with one cold and three warm calls at `10^12`. Version
and thread-count order were shuffled with seed 20261102. All 168 results
matched. Median native times in milliseconds were:

| Workers | State | Original | Atomic cache |
| ---: | --- | ---: | ---: |
| 1 | cold | 85.135 | 85.724 |
| 1 | warm | 51.137 | 52.967 |
| 4 | cold | 102.453 | 28.635 |
| 4 | warm | 18.037 | 17.961 |
| 10 | cold | 161.335 | 24.433 |
| 10 | warm | 15.525 | 16.296 |

The original implementation is an unsafe timing baseline, not a correctness
reference. Cold parallel calls improve substantially; warmed single/ten-worker
calls retain a small overhead. Earlier prototype measurements also found a
small first-call allocation cost at `10^5` (about 19 to 41 microseconds).

The package rebuilt successfully under R 4.6.1, and both focused R test files
above passed. A full package check was not run for this localized internal
change. Validation used Apple clang 17 on arm64 macOS; other platforms have not
been checked in this investigation.

## Worker-reuse experiment (not retained)

On November 3, 2026, a per-call `PhiBatchPool` prototype reused workers for
parallel inputs above `1e10`, preserving the existing intervals and batch
barriers. Each worker received one interval per batch and slept between
batches using condition variables. A mutex protected dispatch and completion;
the counting work ran outside that mutex. Workers were joined at call exit.
Serial counting and the single-batch path at or below `1e10` were unchanged.

Decision: retain the original create/join-per-batch implementation. The small,
inconsistent timing differences below did not justify the extra synchronization
and lifecycle code. The prototype was removed; the atomic cache, 32-bit pi
table, setup improvements, and batch-boundary regression test were retained.
A fresh package build and both focused R test files passed again after
restoring the pre-pool source.

The prototype passed the existing native checker with
ThreadSanitizer at 2, 3, 4, 8, and 10 workers and AddressSanitizer/
UndefinedBehaviorSanitizer at four workers. Each run checked 662 independent
prime counts and 1,394,803 populated cache entries. A clean package build had
no compiler warnings; `testPrimeCount.R` (15 expectations) and `testErrors.R`
passed. The retained R regression crosses `1e10 - 1`, `1e10`, and `1e10 + 1` and
repeats the multiple-batch call with available worker counts from 2 through 4.

Two temporary diagnostics provided additional evidence:

- A compiled comparison replaced `phiForeman` with interval logging and
  compared the old and new schedules, sorting within each batch and excluding
  empty slots. All 3,398 batches and 13,027 nonempty intervals matched in 40
  cases: threshold neighbors and powers of ten through `1e15`, with 2, 3, 4,
  8, and 10 workers. These cases also had contiguous interval coverage.
- A harness extracting the pool class tested 804 batches, empty slots,
  once-only slot execution, worker exceptions followed by reuse, and idle
  shutdown. A modified diagnostic injected 100 constructor failures after
  zero through three successful launches. ThreadSanitizer reported no races.
  The fault injection was confined to the diagnostic, not production code.

These checks used Apple clang on arm64 macOS and C++11 for the native builds.
The experiment did not include a full package check or testing on other
platforms. It did not address the previously identified bounds problems at
excessive effective worker counts.

### Worker-reuse timings

Before/after native binaries used the same compiler and `-std=c++11 -O3
-pthread`. The baseline (`6cff583`) already contained the atomic cache, 32-bit pi table,
prime-vector move, and reuse of the setup count. Execution order was shuffled
with seed 20261103, with no other investigation builds or diagnostics running.
Calls included setup and pool destruction, but excluded process/R startup.

At `1e10`, `1e10 + 1`, `1e12`, and `1e13`, each version/worker combination
used seven fresh processes with one cold and three warm calls. The `1e14`
cases used three processes with one cold and two warm calls. The serial
`1e12` control used five processes with one cold and three warm calls.
All 524 counts matched. Median times in milliseconds were:

| Input | Workers | Before cold | Pool cold | Before warm | Pool warm |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `1e10` (control) | 4 | 2.524 | 2.349 | 0.944 | 0.889 |
| `1e10` (control) | 8 | 2.333 | 2.394 | 0.784 | 0.770 |
| `1e10 + 1` | 4 | 1.922 | 1.907 | 0.997 | 0.993 |
| `1e10 + 1` | 8 | 1.557 | 1.570 | 0.843 | 0.868 |
| `1e12` (serial control) | 1 | 90.300 | 88.509 | 55.134 | 54.359 |
| `1e12` | 4 | 28.685 | 27.915 | 18.804 | 18.385 |
| `1e12` | 8 | 25.131 | 25.194 | 16.900 | 16.552 |
| `1e13` | 4 | 169.418 | 165.382 | 145.663 | 144.065 |
| `1e13` | 8 | 141.419 | 138.838 | 126.845 | 127.900 |
| `1e14` | 4 | 1142.360 | 1103.988 | 1130.194 | 1095.499 |
| `1e14` | 8 | 935.681 | 941.031 | 934.904 | 944.292 |

The launch reduction was substantial: for example, `1e14` with eight workers
used eight threads instead of 875 launches across 110 batches. Elapsed-time
benefits were small and inconsistent, however. Four-worker `1e14` calls
improved about 3%, while eight-worker calls were about 1% slower. Unchanged
control paths also varied, so small differences should not be attributed
confidently to the pool. The prototype passed correctness checks, but did not
demonstrate a compelling performance advantage on this machine.

## Native registration

No registration changes are needed for internal prime-counting fixes. When a
registered interface changes, the verified workflow is
`cpp11::cpp_register(".")`; see [cache-inspector.md](cache-inspector.md)
for this repository's handling of the generated R wrappers.
