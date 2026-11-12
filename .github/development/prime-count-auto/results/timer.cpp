// Test-only monotonic timer. A per-process origin keeps elapsed doubles small
// regardless of host uptime. Both calls happen on R's calling thread.
#include <chrono>
#include <Rinternals.h>
#include <R_ext/Visibility.h>

extern "C" attribute_visible SEXP AutoDispatchMilliseconds() {
    using Clock = std::chrono::steady_clock;
    static const Clock::time_point origin = Clock::now();
    return Rf_ScalarReal(std::chrono::duration<double, std::milli>(
        Clock::now() - origin
    ).count());
}
