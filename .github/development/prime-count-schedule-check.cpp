// Scheduler bounds and coverage checks, without creating OS threads. Prepare
// prime-count-schedule-core.cpp with prepare-prime-count-schedule-check.py;
// compile with its directory and inst/include on the include path, linking
// src/Phi.cpp and src/Eratosthenes.cpp. See testing.md for commands.
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ScheduleCheck {
    struct Interval {
        int batch, slot;
        std::string region;
        std::int64_t lower, upper;
    };

    struct State {
        bool scheduleOnly = false;
        std::int64_t next = 6, a = -1, requestedA = -1;
        int effectiveThreads = 0, batch = -1, started = 0, joined = 0;
        const char* region = "serial";
        std::int64_t* sums = nullptr;
        std::size_t size = 0;
        std::vector<Interval> intervals;
    };
    State state;

    // Execute the production worker lambda immediately. The production
    // batch still owns joinable handles, catches worker exceptions, joins,
    // and propagates errors through its ordinary finish() path.
    class Thread {
        bool pending = true;
    public:
        template <typename Function>
        explicit Thread(Function&& function) {
            ++state.started;
            function();
        }
        Thread(const Thread&) = delete;
        Thread& operator=(const Thread&) = delete;
        Thread(Thread&& other) noexcept : pending(other.pending) {
            other.pending = false;
        }
        ~Thread() {
            // Mirror std::thread's lifetime rule, including early returns.
            if (pending) std::terminate();
        }
        bool joinable() const { return pending; }
        void join() {
            if (!pending) throw std::logic_error("joining a non-joinable worker");
            pending = false;
            ++state.joined;
        }
    };

    void BeforeSchedule(int threads) {
        // This ceiling is a fixture safeguard, not a production thread rule.
        // All retained worker requests fit below it after the real cap.
        if (threads < 1 || threads > 100000) {
            throw std::runtime_error("unsafe test allocation: unbounded worker request " +
                                     std::to_string(threads));
        }
        state.effectiveThreads = threads;
    }

    void BeginBatch(const char* region, std::vector<std::int64_t>& sums) {
        ++state.batch;
        state.region = region;
        state.sums = sums.data();
        state.size = sums.size();
    }

    bool BeforeForeman(std::int64_t& sum, std::int64_t lower,
                       std::int64_t upper, std::int64_t tableCount) {
        const std::int64_t a = state.requestedA < 0 ? tableCount : state.requestedA;
        state.a = a;
        if (a > tableCount || lower < 6 || lower > upper || upper > a ||
            lower != state.next) {
            std::ostringstream error;
            error << "invalid interval [" << lower << ',' << upper
                  << "), expected next=" << state.next << " within [6," << a
                  << "), batch=" << state.batch << " region=" << state.region;
            throw std::runtime_error(error.str());
        }
        int slot = -1;
        if (state.sums != nullptr) {
            // The instrumented foreman receives the production result slot.
            // Find it without subtracting potentially unrelated pointers.
            for (std::size_t i = 0; i < state.size; ++i) {
                if (&sum == &state.sums[i]) {
                    slot = static_cast<int>(i);
                    break;
                }
            }
            if (slot < 0) throw std::runtime_error("worker result is outside its batch");
        }
        state.intervals.push_back({state.batch, slot, state.region, lower, upper});
        state.next = upper;
        return state.scheduleOnly;
    }
}

#include "prime-count-schedule-core.cpp"

namespace {
    int cases = 0, realCounts = 0, phiCounts = 0;
    std::size_t intervals = 0;

    void Check(std::int64_t n, int threads, int maximum, bool scheduleOnly,
               std::int64_t expected = -1, bool dump = false,
               bool requireUsefulCap = true) {
        ScheduleCheck::state = ScheduleCheck::State();
        ScheduleCheck::state.scheduleOnly = scheduleOnly;
        std::int64_t result;
        try {
            result = PrimeCounting::PiLegendre(n, threads, maximum);
        } catch (const std::exception& error) {
            throw std::runtime_error("n=" + std::to_string(n) + " threads=" +
                                     std::to_string(threads) + " maxThreads=" +
                                     std::to_string(maximum) + ": " + error.what());
        }
        const auto& state = ScheduleCheck::state;
        if (state.a < 6 || state.next != state.a) {
            throw std::runtime_error("scheduler did not cover [6,a) completely");
        }
        const bool parallel = n >= 10000000 && threads > 1 && maximum > 1;
        if (parallel != (state.effectiveThreads != 0)) {
            throw std::runtime_error("parallel threshold or maxThreads handling changed");
        }
        if (state.started != state.joined ||
            (parallel && static_cast<std::size_t>(state.started) != state.intervals.size())) {
            throw std::runtime_error("workers were not all executed and joined");
        }
        if (requireUsefulCap && state.effectiveThreads > state.a - 6) {
            throw std::runtime_error("effective workers exceed available terms");
        }
        if (!scheduleOnly && expected >= 0 && result != expected) {
            throw std::runtime_error("incorrect prime count: " + std::to_string(result) +
                                     " instead of " + std::to_string(expected));
        }
        if (dump) {
            for (const auto& interval : state.intervals) {
                std::cout << n << ',' << threads << ',' << maximum << ','
                          << interval.batch << ',' << interval.region << ','
                          << interval.slot << ',' << interval.lower << ','
                          << interval.upper << '\n';
            }
        }
        ++cases;
        realCounts += !scheduleOnly;
        intervals += state.intervals.size();
    }

    void NormalSchedules(bool dump) {
        const std::int64_t inputs[] = {100000, 9999999, 10000000, 10000001,
                                      9999999999LL, 10000000000LL, 10000000001LL,
                                      1000000000000LL, 10000000000000LL};
        for (std::int64_t n : inputs) {
            for (int threads : {1, 2, 4, 8, 10}) {
                Check(n, threads, threads, true, -1, dump);
            }
        }
        for (int maximum : {1, 2, 4, 10}) {
            Check(10000000, 1000, maximum, true, -1, dump);
            Check(10000000001LL, 1000, maximum, true, -1, dump);
        }
    }

#if SCHEDULE_HAS_PHI_PARALLEL
    // Independent inclusion-exclusion for the short prefixes used below.
    // Each recursive branch adds one selected prime to the denominator.
    std::int64_t InclusionExclusion(std::int64_t x, std::int64_t first,
                                    std::int64_t a,
                                    const std::vector<std::int64_t>& primes) {
        std::int64_t result = x;
        for (std::int64_t i = first; i <= a && primes[i] <= x; ++i) {
            result -= InclusionExclusion(x / primes[i], i + 1, a, primes);
        }
        return result;
    }

    void CheckPhi(std::int64_t x, std::int64_t a, int threads,
                  const PrimeCounting::PhiTables& tables) {
        ScheduleCheck::state = ScheduleCheck::State();
        ScheduleCheck::state.requestedA = a;
        const std::int64_t expected = a <= 12 ?
            InclusionExclusion(x, 1, a, tables.primes()) :
            PrimeCounting::Phi(x, a, tables);
        const std::int64_t result = PrimeCounting::PhiParallel(x, a, threads, tables);
        const auto& state = ScheduleCheck::state;
        const bool delegated = a <= 6 || x <= tables.limit();
        const bool parallel = !delegated && x >= 10000000 && threads > 1;
        const bool coverage = delegated ? state.intervals.empty() :
            state.a == a && state.next == a;
        if (result != expected || !coverage ||
            parallel != (state.effectiveThreads != 0) ||
            state.started != state.joined ||
            (parallel && (state.effectiveThreads > a - 6 ||
                static_cast<std::size_t>(state.started) != state.intervals.size()))) {
            throw std::runtime_error("general Phi coverage/count failure: x=" +
                std::to_string(x) + " a=" + std::to_string(a) +
                " threads=" + std::to_string(threads));
        }
        ++cases;
        ++phiCounts;
        intervals += state.intervals.size();
    }

    void GeneralPhiSchedules() {
        const PrimeCounting::PhiTables tables(100000);
        const int huge = std::numeric_limits<int>::max();
        // The requested a is deliberately smaller than tables.primeCount().
        // a=7 leaves one term, and a=8 leaves two; excessive requests must
        // cap against this requested prefix rather than the prepared extent.
        for (std::int64_t x : {9999999LL, 10000000LL, 10000000001LL}) {
            for (std::int64_t a : {7, 8, 12, 64}) {
                for (int threads : {1, 4, huge}) CheckPhi(x, a, threads, tables);
            }
        }
        for (std::int64_t x : {0, 17, 100000}) {
            for (std::int64_t a : {0, 6, 7, 64}) CheckPhi(x, a, huge, tables);
        }
        for (int threads : {-1, 0, 1, 4, huge}) {
            CheckPhi(10000000001LL, tables.primeCount() - 1, threads, tables);
        }
        CheckPhi(10000000001LL, 6, huge, tables);
    }
#endif
}

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--dump-normal") {
            NormalSchedules(true);
            return 0;
        }
        if (argc == 5 && std::string(argv[1]) == "--case") {
            // Useful for checking a before-change snapshot. Bounds checks
            // precede real evaluation and the post-fix cap is not assumed.
            const std::int64_t n = std::stoll(argv[2]);
            const int threads = std::stoi(argv[3]);
            const bool scheduleOnly = std::string(argv[4]) == "schedule";
            Check(n, threads, threads, scheduleOnly, -1, false, false);
        } else if (argc == 1) {
            NormalSchedules(false);
            const int huge = std::numeric_limits<int>::max();
            for (int threads : {440, 441, 442, 1000, huge}) {
                Check(10000000, threads, threads, false, 664579);
            }
            for (std::int64_t n : {9999999LL, 10000000LL, 10000001LL}) {
                Check(n, huge, huge, false, 664579);
            }
            for (std::int64_t n : {9999999999LL, 10000000000LL, 10000000001LL}) {
                for (int threads : {1369, 1370, 2000, 10000, huge}) {
                    Check(n, threads, threads, true);
                }
            }
            for (int threads : {1370, 2000, 10000, huge}) {
                Check(10000000001LL, threads, threads, false, 455052511);
            }
            Check(10000000, huge, 1, false, 664579);
            Check(10000000, huge, 4, false, 664579);
#if SCHEDULE_HAS_PHI_PARALLEL
            GeneralPhiSchedules();
#endif
        } else {
            throw std::runtime_error("usage: checker [--dump-normal | --case n threads real|schedule]");
        }
        std::cout << "PASS: " << cases << " cases, " << realCounts
                  << " real prime counts, " << phiCounts << " general Phi counts, " << intervals
                  << " contiguous bounded intervals; no OS workers created\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
