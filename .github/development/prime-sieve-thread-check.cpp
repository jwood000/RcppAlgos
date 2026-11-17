// Prepare prime-sieve-thread-core.cpp with prepare-prime-sieve-thread-check.py.
// Compile this file alone with inst/include and the prepared directory on -I.
// The source copy substitutes only the thread wrapper and worker-entry hook;
// production code has no fault-injection hooks.
#include <array>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <new>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

namespace SieveThreadFaults {
    enum class Fault { None, Standard, NonStandard, BadAlloc, Ordered };
    struct WorkerError : std::runtime_error {
        int slot;
        explicit WorkerError(int value) : std::runtime_error(
            "injected sieve worker failure " + std::to_string(value)), slot(value) {}
    };
    struct LaunchError : std::system_error {
        LaunchError() : std::system_error(
            std::make_error_code(std::errc::resource_unavailable_try_again),
            "injected sieve thread launch failure") {}
    };
    constexpr int NONSTANDARD_ERROR = 8410;
    const std::chrono::seconds TIMEOUT(20);

    struct Scenario {
        Fault fault;
        int launchFailure, faultSlot;
        bool gated, combined;
        Scenario(Fault value = Fault::None, int launch = -1, int slot = 0,
                 bool gate = false, bool together = false)
            : fault(value), launchFailure(launch), faultSlot(slot),
              gated(gate), combined(together) {}
    };
    struct Counters {
        int attempts = 0, started = 0, bodies = 0, exited = 0;
        int joinAttempts = 0, joined = 0, gateArrivals = 0;
        bool launchFailed = false, returned = false;
        std::array<bool, 8> bodyDone = {{false}};
        std::vector<int> throwOrder;
        std::string protocolError;
    };
    std::mutex mutex;
    std::condition_variable changed;
    Scenario scenario;
    Counters state;
    bool release = false;
    thread_local int workerSlot = -1;

    void Configure(const Scenario& next = Scenario()) {
        std::lock_guard<std::mutex> lock(mutex);
        scenario = next;
        state = Counters();
        release = false;
    }
    Counters Inspect() {
        std::lock_guard<std::mutex> lock(mutex);
        return state;
    }
    int BeforeLaunch() {
        std::unique_lock<std::mutex> lock(mutex);
        const int slot = state.attempts++;
        if (slot >= 8) throw std::logic_error("fixture expects at most eight workers");
        if (slot == scenario.launchFailure) {
            if (scenario.combined && !changed.wait_for(lock, TIMEOUT, [] {
                    return state.bodyDone[0];
                })) state.protocolError = "worker exception did not finish before launch failure";
            state.launchFailed = true;
            changed.notify_all();
            throw LaunchError();
        }
        return slot;
    }
    void BeforeWorker() {
        std::unique_lock<std::mutex> lock(mutex);
        const int slot = workerSlot < 0 ? 0 : workerSlot;
        Fault selected = scenario.fault;
        if (selected == Fault::Ordered && slot < 2) {
            if (slot == 0) {
                if (!changed.wait_for(lock, TIMEOUT, [] { return state.bodyDone[1]; })) {
                    state.protocolError = "higher-slot exception did not finish first";
                }
                selected = Fault::Standard;
            } else selected = Fault::NonStandard;
        } else if (slot != scenario.faultSlot) return;
        if (selected == Fault::None || selected == Fault::Ordered) return;
        state.throwOrder.push_back(slot);
        lock.unlock();
        if (selected == Fault::Standard) throw WorkerError(slot);
        if (selected == Fault::BadAlloc) throw std::bad_alloc();
        throw NONSTANDARD_ERROR + slot;
    }
    void Started(int slot) {
        std::lock_guard<std::mutex> lock(mutex);
        workerSlot = slot;
        ++state.started;
        changed.notify_all();
    }
    void AfterBody(int slot) {
        std::unique_lock<std::mutex> lock(mutex);
        ++state.bodies;
        state.bodyDone[slot] = true;
        if (scenario.gated) {
            ++state.gateArrivals;
            changed.notify_all();
            if (!changed.wait_for(lock, TIMEOUT, [] { return release; })) {
                state.protocolError = "worker-exit gate timed out";
                release = true;
                changed.notify_all();
            }
        }
        ++state.exited;
        changed.notify_all();
    }
    // Actual std::threads retain their normal destruction and join semantics.
    class Thread {
        int slot;
        std::thread actual;
    public:
        template <typename Function>
        explicit Thread(Function work) : slot(BeforeLaunch()) {
            const int captured = slot;
            actual = std::thread([work, captured]() mutable {
                Started(captured);
                work();
                AfterBody(captured);
            });
        }
        Thread(Thread&&) = default;
        Thread& operator=(Thread&&) = default;
        Thread(const Thread&) = delete;
        Thread& operator=(const Thread&) = delete;
        bool joinable() const { return actual.joinable(); }
        void join() {
            {
                std::lock_guard<std::mutex> lock(mutex);
                ++state.joinAttempts;
                changed.notify_all();
            }
            actual.join();
            {
                std::lock_guard<std::mutex> lock(mutex);
                ++state.joined;
                changed.notify_all();
            }
        }
    };
    void Returned() {
        std::lock_guard<std::mutex> lock(mutex);
        state.returned = true;
        changed.notify_all();
    }
    bool ObserveBlockedJoin(int workers) {
        std::unique_lock<std::mutex> lock(mutex);
        const bool ready = changed.wait_for(lock, TIMEOUT, [&] {
            return state.returned || (state.started == workers &&
                state.bodies == workers && state.gateArrivals == workers &&
                state.joinAttempts > 0);
        });
        return ready && !state.returned && state.exited == 0 && state.joined == 0;
    }
    bool WaitReturned() {
        std::unique_lock<std::mutex> lock(mutex);
        return changed.wait_for(lock, TIMEOUT, [] { return state.returned; });
    }
    void Release() {
        std::lock_guard<std::mutex> lock(mutex);
        release = true;
        changed.notify_all();
    }
}

#include "prime-sieve-thread-core.cpp"

namespace SieveThreadChecks {
    void Check(bool value, const std::string& message) {
        if (!value) throw std::runtime_error(message);
    }
    void CheckJoined(const SieveThreadFaults::Counters& state) {
        Check(state.protocolError.empty(), state.protocolError);
        Check(state.started == state.bodies && state.started == state.exited &&
              state.started == state.joined, "caller observed unjoined workers");
    }
    template <typename T = std::int64_t>
    std::vector<T> Invoke(int maximum, int threads) {
        std::vector<std::vector<T>> chunks(threads);
        std::vector<T> serial;
        bool parallel = false;
        PrimeSieve::PrimeSieveMain(chunks, serial, 1, maximum, parallel,
                                   threads, threads, 2);
        if (!parallel) return serial;
        // Reduced budgets can leave unused middle vectors; preserve range order.
        for (const auto& chunk : chunks) serial.insert(serial.end(), chunk.begin(), chunk.end());
        return serial;
    }
    std::vector<std::int64_t> Reference(int maximum) {
        std::vector<bool> prime(maximum + 1, true);
        prime[0] = prime[1] = false;
        for (int p = 2; p <= maximum / p; ++p) {
            if (prime[p]) for (int multiple = p * p; multiple <= maximum; multiple += p) {
                prime[multiple] = false;
            }
        }
        std::vector<std::int64_t> result;
        for (int p = 2; p <= maximum; ++p) if (prime[p]) result.push_back(p);
        return result;
    }
    void Recover(int maximum, int threads, bool countOnly = false) {
        SieveThreadFaults::Configure();
        if (countOnly) {
            Check(PrimeSieve::CountPrimes(1, maximum, threads, threads, 2) ==
                      static_cast<std::int64_t>(Reference(maximum).size()),
                  "sieve count disagrees with independent reference");
        } else {
            Check(Invoke(maximum, threads) == Reference(maximum),
                  "sieve output disagrees with independent reference");
        }
        CheckJoined(SieveThreadFaults::Inspect());
    }
    void RunScenario(const SieveThreadFaults::Scenario& scenario,
                     int maximum, int threads, int workers, bool countOnly) {
        using namespace SieveThreadFaults;
        Configure(scenario);
        enum class Result { Success, Worker, NonStandard, BadAlloc, Launch, Unexpected };
        Result result = Result::Unexpected;
        int errorSlot = -1;
        std::string unexpected;
        std::thread coordinator([&] {
            try {
                if (countOnly) {
                    PrimeSieve::CountPrimes(1, maximum, threads, threads, 2);
                } else Invoke(maximum, threads);
                result = Result::Success;
            } catch (const LaunchError&) { result = Result::Launch;
            } catch (const WorkerError& error) {
                result = Result::Worker; errorSlot = error.slot;
            } catch (const std::bad_alloc&) { result = Result::BadAlloc;
            } catch (int error) {
                result = Result::NonStandard; errorSlot = error - NONSTANDARD_ERROR;
            } catch (const std::exception& error) { unexpected = error.what();
            } catch (...) { unexpected = "unknown exception"; }
            Returned();
        });
        const int expected = scenario.launchFailure >= 0 ? scenario.launchFailure : workers;
        const bool blocked = expected ? ObserveBlockedJoin(expected) : WaitReturned();
        Release();
        coordinator.join();
        const auto snapshot = Inspect();
        Check(blocked, "caller propagated before workers finished, or gate timed out");
        CheckJoined(snapshot);
        Check(snapshot.started == expected, "wrong number of launched workers");
        Check(snapshot.attempts == expected + (scenario.launchFailure >= 0),
              "unexpected launch after failure");
        Check(unexpected.empty(), unexpected);
        if (scenario.launchFailure >= 0) {
            Check(result == Result::Launch && snapshot.launchFailed,
                  "launch failure did not retain precedence");
            if (scenario.combined) Check(snapshot.throwOrder == std::vector<int>{0},
                                          "combined worker failure was not exercised");
        } else if (scenario.fault == Fault::BadAlloc) {
            Check(result == Result::BadAlloc, "allocation failure was not preserved");
        } else if (scenario.fault == Fault::NonStandard) {
            Check(result == Result::NonStandard && errorSlot == scenario.faultSlot,
                  "nonstandard exception was not preserved");
        } else {
            Check(result == Result::Worker && errorSlot == 0,
                  "standard worker exception was not preserved");
            if (scenario.fault == Fault::Ordered) {
                Check(snapshot.throwOrder == std::vector<int>({1, 0}),
                      "out-of-order worker failures were not exercised");
            }
        }
        Recover(maximum, threads, countOnly);
    }
}

int main() {
    try {
        using namespace SieveThreadFaults;
        using namespace SieveThreadChecks;
        int scenarios = 0, blocked = 0;
        for (int fixture = 0; fixture < 4; ++fixture) {
            const bool countOnly = fixture >= 2;
            const bool reduced = fixture % 2;
            const int maximum = reduced ? 450000 : 1000000;
            const int threads = reduced ? 8 : 4;
            const int workers = reduced ? 3 : 4;
            Recover(maximum, threads, countOnly);
            Check(Inspect().started == workers, "fixture no longer exercises intended budget");
            const std::vector<Scenario> cases = {
                Scenario(Fault::Standard, -1, 0, true),
                Scenario(Fault::NonStandard, -1, workers - 1, true),
                Scenario(Fault::BadAlloc, -1, 0, true),
                Scenario(Fault::Ordered, -1, 0, true),
                Scenario(Fault::None, 0, 0, true),
                Scenario(Fault::None, 1, 0, true),
                Scenario(Fault::None, workers - 1, 0, true),
                Scenario(Fault::Standard, 1, 0, true, true)
            };
            for (const auto& scenario : cases) {
                RunScenario(scenario, maximum, threads, workers, countOnly);
                ++scenarios;
                blocked += scenario.launchFailure != 0;
            }
        }
        Recover(1000000, 1);
        Check(Inspect().started == 0, "serial control unexpectedly launched workers");
        Recover(1000000, 1, true);
        Check(Inspect().started == 0, "serial count unexpectedly launched workers");
        const auto reference = Reference(450000);
        Configure();
        const auto integers = Invoke<int>(450000, 8);
        Check(integers.size() == reference.size() &&
              std::equal(integers.begin(), integers.end(), reference.begin()),
              "integer parallel output disagrees with reference");
        CheckJoined(Inspect());
        Configure();
        const auto doubles = Invoke<double>(450000, 8);
        Check(doubles.size() == reference.size() &&
              std::equal(doubles.begin(), doubles.end(), reference.begin()),
              "double parallel output disagrees with reference");
        CheckJoined(Inspect());
        // Exercise Big's hardware-derived sizing before allocating its segment.
        // Ratio 1024 used to overflow the int multiplication; INT_MAX also
        // exercises the segment-offset cap. The worker-entry hook throws after
        // real setup, avoiding hundreds of MiB for synthetic hardware budgets.
        for (int maximum : {1024, std::numeric_limits<int>::max()}) {
            Configure(Scenario(Fault::Standard));
            bool caught = false;
            try {
                PrimeSieve::CountPrimes(1000000000000000LL,
                    1000000000000031LL, 1, maximum, 1);
            } catch (const WorkerError&) { caught = true; }
            Check(caught, "synthetic Big sizing did not reach the worker hook");
            CheckJoined(Inspect());
            Check(Inspect().started == 0, "synthetic Big sizing launched a worker");
        }
        std::cout << "PASS scenarios=" << scenarios << " blocked_joins=" << blocked
                  << " recovery_checks=" << scenarios + 4
                  << " output_types=3 count_only=1 reduced_budget=3/8 sizing_checks=2\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
