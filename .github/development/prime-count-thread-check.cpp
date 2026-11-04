// Native exception/lifetime checks. Prepare the included copy with
// prepare-prime-count-thread-check.py, then compile with its output directory
// on the include path and link src/Phi.cpp and src/Eratosthenes.cpp. See testing.md.
#include <array>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <iostream>
#include <mutex>
#include <new>
#include <stdexcept>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

namespace ThreadFaults {
    enum class Region { First, Middle, Final, Small };
    enum class WorkerFault { None, Standard, NonStandard, BadAlloc, Ordered };

    const char* Name(Region region) {
        switch (region) {
            case Region::First: return "first";
            case Region::Middle: return "middle";
            case Region::Final: return "final";
            case Region::Small: return "small";
        }
        return "unknown";
    }

    struct WorkerError : std::runtime_error {
        int slot;
        explicit WorkerError(int value)
            : std::runtime_error("injected worker failure " + std::to_string(value)), slot(value) {}
    };
    struct LaunchError : std::system_error {
        LaunchError() : std::system_error(
            std::make_error_code(std::errc::resource_unavailable_try_again),
            "injected thread launch failure") {}
    };
    constexpr int NONSTANDARD_ERROR = 7410;

    struct Scenario {
        Region region;
        WorkerFault worker;
        int launchFailure;
        bool gated, simultaneous;
        int workerSlot;
        bool enabled;

        Scenario(Region r = Region::Small, WorkerFault w = WorkerFault::None,
                 int launch = -1, bool gate = true, bool together = false,
                 int slot = 0, bool active = true)
            : region(r), worker(w), launchFailure(launch), gated(gate),
              simultaneous(together), workerSlot(slot), enabled(active) {}
    };

    struct BatchRecord {
        Region region;
        int attempts = 0, started = 0, bodies = 0, exited = 0;
        int joinAttempts = 0, joined = 0;
        std::array<bool, 4> bodyDone = {{false, false, false, false}};
        explicit BatchRecord(Region value) : region(value) {}
    };

    struct State {
        std::mutex mutex;
        std::condition_variable changed;
        Scenario scenario;
        std::vector<BatchRecord> batches;
        std::vector<int> throwOrder;
        int target = -1, current = -1, gateArrivals = 0;
        bool release = false, returned = false, launchFailed = false;
        std::string protocolError;
    };
    State state;
    thread_local int workerBatch = -1;
    thread_local int workerSlot = -1;
    const std::chrono::seconds TIMEOUT(20);

    void Configure(const Scenario& scenario) {
        std::lock_guard<std::mutex> lock(state.mutex);
        state.scenario = scenario;
        state.batches.clear();
        state.throwOrder.clear();
        state.target = state.current = -1;
        state.gateArrivals = 0;
        state.release = state.returned = state.launchFailed = false;
        state.protocolError.clear();
    }

    void Observe() {
        Configure(Scenario(Region::Small, WorkerFault::None, -1, false, false, 0, false));
    }

    void BeginBatch(Region region) {
        std::lock_guard<std::mutex> lock(state.mutex);
        state.batches.emplace_back(region);
        state.current = static_cast<int>(state.batches.size()) - 1;
        if (state.scenario.enabled && state.target < 0 && region == state.scenario.region) {
            state.target = state.current;
        }
        state.changed.notify_all();
    }

    struct Ticket { int batch, slot; };

    Ticket BeforeLaunch() {
        std::unique_lock<std::mutex> lock(state.mutex);
        const int batch = state.current;
        if (batch < 0) throw std::logic_error("thread launch has no batch tag");
        const int slot = state.batches[batch].attempts++;
        if (slot >= 4) throw std::logic_error("test expects at most four workers per batch");
        if (batch == state.target && slot == state.scenario.launchFailure) {
            if (state.scenario.simultaneous) {
                // Let worker zero throw and return through the production catch
                // before simulating failure to construct the next thread.
                if (!state.changed.wait_for(lock, TIMEOUT, [&] {
                        return state.batches[batch].bodyDone[0];
                    })) {
                    state.protocolError = "worker did not finish before launch failure";
                }
            }
            state.launchFailed = true;
            state.changed.notify_all();
            throw LaunchError();
        }
        return {batch, slot};
    }

    void BeforeForeman() {
        std::unique_lock<std::mutex> lock(state.mutex);
        if (workerBatch != state.target || workerBatch < 0) return;
        const WorkerFault fault = state.scenario.worker;
        bool standard = false, nonstandard = false, badAllocation = false;
        if (fault == WorkerFault::Ordered && workerSlot < 2) {
            if (workerSlot == 0) {
                // Force the higher slot to be captured first. Production must
                // nevertheless propagate the lower slot's exception.
                if (!state.changed.wait_for(lock, TIMEOUT, [&] {
                        return state.batches[workerBatch].bodyDone[1];
                    })) {
                    state.protocolError = "higher-slot exception did not finish first";
                }
                standard = true;
            } else {
                nonstandard = true;
            }
        } else if (workerSlot == state.scenario.workerSlot) {
            standard = fault == WorkerFault::Standard;
            nonstandard = fault == WorkerFault::NonStandard;
            badAllocation = fault == WorkerFault::BadAlloc;
        }
        if (!standard && !nonstandard && !badAllocation) return;
        state.throwOrder.push_back(workerSlot);
        lock.unlock();
        if (badAllocation) throw std::bad_alloc();
        if (standard) throw WorkerError(workerSlot);
        throw NONSTANDARD_ERROR + workerSlot;
    }

    void WorkerStarted(const Ticket& ticket) {
        std::lock_guard<std::mutex> lock(state.mutex);
        ++state.batches[ticket.batch].started;
        workerBatch = ticket.batch;
        workerSlot = ticket.slot;
        state.changed.notify_all();
    }

    void AfterBody(const Ticket& ticket) {
        std::unique_lock<std::mutex> lock(state.mutex);
        auto& batch = state.batches[ticket.batch];
        batch.bodyDone[ticket.slot] = true;
        ++batch.bodies;
        if (ticket.batch == state.target && state.scenario.gated) {
            ++state.gateArrivals;
            state.changed.notify_all();
            if (!state.changed.wait_for(lock, TIMEOUT, [] { return state.release; })) {
                state.protocolError = "worker-exit gate timed out";
                state.release = true;
                state.changed.notify_all();
            }
        }
        ++state.batches[ticket.batch].exited;
        state.changed.notify_all();
    }

    // Only the generated test copy uses this wrapper; the actual worker is
    // still an ordinary std::thread, with its normal join/destruction rules.
    class Thread {
        Ticket ticket;
        std::thread actual;

    public:
        template <typename Callable>
        explicit Thread(Callable call) : ticket(BeforeLaunch()) {
            const Ticket captured = ticket;
            actual = std::thread([call, captured]() mutable {
                WorkerStarted(captured);
                call();
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
                std::lock_guard<std::mutex> lock(state.mutex);
                ++state.batches[ticket.batch].joinAttempts;
                state.changed.notify_all();
            }
            actual.join();
            {
                std::lock_guard<std::mutex> lock(state.mutex);
                ++state.batches[ticket.batch].joined;
                state.changed.notify_all();
            }
        }
    };

    void MarkReturned() {
        std::lock_guard<std::mutex> lock(state.mutex);
        state.returned = true;
        state.changed.notify_all();
    }

    struct Snapshot {
        std::vector<BatchRecord> batches;
        std::vector<int> throwOrder;
        int target;
        bool launchFailed;
        std::string protocolError;
    };
    Snapshot Inspect() {
        std::lock_guard<std::mutex> lock(state.mutex);
        return {state.batches, state.throwOrder, state.target,
                state.launchFailed, state.protocolError};
    }

    bool ObserveBlockedJoin(int expectedWorkers) {
        std::unique_lock<std::mutex> lock(state.mutex);
        const bool ready = state.changed.wait_for(lock, TIMEOUT, [&] {
            if (state.returned) return true;
            if (state.target < 0) return false;
            const auto& batch = state.batches[state.target];
            return batch.started == expectedWorkers && batch.bodies == expectedWorkers &&
                   state.gateArrivals == expectedWorkers && batch.joinAttempts > 0;
        });
        if (!ready || state.returned || state.target < 0) return false;
        const auto& batch = state.batches[state.target];
        return batch.exited == 0 && batch.joined == 0;
    }

    bool WaitReturned() {
        std::unique_lock<std::mutex> lock(state.mutex);
        return state.changed.wait_for(lock, TIMEOUT, [] { return state.returned; });
    }

    void ReleaseGate() {
        std::lock_guard<std::mutex> lock(state.mutex);
        state.release = true;
        state.changed.notify_all();
    }
}

#include "prime-count-thread-core.cpp"

namespace {
    void Check(bool condition, const std::string& message) {
        if (!condition) throw std::runtime_error(message);
    }

    void CheckAllJoined(const ThreadFaults::Snapshot& snapshot) {
        Check(snapshot.protocolError.empty(), snapshot.protocolError);
        for (const auto& batch : snapshot.batches) {
            Check(batch.started == batch.bodies && batch.started == batch.exited &&
                  batch.started == batch.joined, "started worker was not completed and joined");
        }
    }

    void Recover(std::int64_t n, std::int64_t expected) {
        ThreadFaults::Observe();
        Check(PrimeCounting::PiLegendre(n, 4, 4) == expected,
              "correct prime count did not recover after an injected failure");
        CheckAllJoined(ThreadFaults::Inspect());
    }

    void RunScenario(const ThreadFaults::Scenario& scenario, int batchSize,
                     std::int64_t n, std::int64_t expected) {
        using namespace ThreadFaults;
        Configure(scenario);
        enum class Result { Success, Worker, NonStandard, BadAlloc, Launch, Unexpected };
        Result result = Result::Unexpected;
        int errorSlot = -1;
        std::string unexpected;
        std::thread coordinator([&] {
            try {
                PrimeCounting::PiLegendre(n, 4, 4);
                result = Result::Success;
            } catch (const LaunchError&) {
                result = Result::Launch;
            } catch (const WorkerError& error) {
                result = Result::Worker;
                errorSlot = error.slot;
            } catch (const std::bad_alloc&) {
                result = Result::BadAlloc;
            } catch (int error) {
                result = Result::NonStandard;
                errorSlot = error - NONSTANDARD_ERROR;
            } catch (const std::exception& error) {
                unexpected = error.what();
            } catch (...) {
                unexpected = "unknown exception";
            }
            MarkReturned();
        });

        const int expectedWorkers = scenario.launchFailure >= 0 ? scenario.launchFailure : batchSize;
        const bool observed = expectedWorkers ? ObserveBlockedJoin(expectedWorkers) : WaitReturned();
        // Always unblock/join before making assertions, including on a failed test.
        ReleaseGate();
        coordinator.join();
        const auto snapshot = Inspect();
        const std::string label = std::string(Name(scenario.region)) + " batch: ";
        Check(observed, label + "propagation preceded worker completion, or gate timed out");
        CheckAllJoined(snapshot);
        Check(snapshot.target >= 0, label + "target region was not exercised");
        Check(snapshot.batches.size() == static_cast<std::size_t>(snapshot.target + 1),
              label + "another batch launched after the failed batch");
        const auto& batch = snapshot.batches[snapshot.target];
        Check(batch.started == expectedWorkers, label + "wrong number of started workers");
        Check(batch.attempts == expectedWorkers + (scenario.launchFailure >= 0 ? 1 : 0),
              label + "unexpected thread launch after failure");
        Check(unexpected.empty(), label + unexpected);
        if (scenario.launchFailure >= 0) {
            Check(result == Result::Launch && snapshot.launchFailed,
                  label + "launch failure did not retain precedence");
            if (scenario.simultaneous) {
                Check(snapshot.throwOrder == std::vector<int>{0},
                      label + "combined worker failure was not exercised");
            }
        } else if (scenario.worker == WorkerFault::NonStandard) {
            Check(result == Result::NonStandard && errorSlot == scenario.workerSlot,
                  label + "non-standard worker exception was not preserved");
        } else if (scenario.worker == WorkerFault::BadAlloc) {
            Check(result == Result::BadAlloc,
                  label + "allocation failure was not preserved");
        } else {
            Check(result == Result::Worker && errorSlot == 0,
                  label + "standard worker exception was not preserved");
            if (scenario.worker == WorkerFault::Ordered) {
                Check(snapshot.throwOrder == std::vector<int>({1, 0}),
                      label + "out-of-order failures were not exercised");
            }
        }
        Recover(n, expected);
    }
}

int main() {
    try {
        using namespace ThreadFaults;
        const std::int64_t largeN = 10000000001LL, smallN = 10000000000LL;
        const std::int64_t finalN = 999999999999LL, finalExpected = 37607912018LL;
        constexpr std::int64_t expected = 455052511;
        Recover(largeN, expected);
        const auto large = Inspect();
        int sizes[4] = {0, 0, 0, 0};
        for (const auto& batch : large.batches) {
            auto& size = sizes[static_cast<int>(batch.region)];
            if (!size) size = batch.started;
        }
        Recover(finalN, finalExpected);
        for (const auto& batch : Inspect().batches) {
            if (batch.region == Region::Final) sizes[static_cast<int>(Region::Final)] = batch.started;
        }
        Recover(smallN, expected);
        const auto small = Inspect();
        sizes[static_cast<int>(Region::Small)] = small.batches.front().started;
        Check(sizes[0] == 4 && sizes[1] == 4 && sizes[2] == 2 && sizes[3] == 4,
              "fixtures no longer cover all regions and a partial final batch: " +
              std::to_string(sizes[0]) + "," + std::to_string(sizes[1]) + "," +
              std::to_string(sizes[2]) + "," + std::to_string(sizes[3]));

        int scenarios = 0, blockedJoins = 0;
        for (Region region : {Region::First, Region::Middle, Region::Final, Region::Small}) {
            const int size = sizes[static_cast<int>(region)];
            const auto n = region == Region::Final ? finalN :
                           region == Region::Small ? smallN : largeN;
            const auto count = region == Region::Final ? finalExpected : expected;
            std::vector<Scenario> cases = {
                Scenario(region, WorkerFault::Standard),
                Scenario(region, WorkerFault::NonStandard, -1, true, false, size - 1),
                Scenario(region, WorkerFault::Ordered),
                Scenario(region, WorkerFault::None, 0),
                Scenario(region, WorkerFault::None, 1),
                Scenario(region, WorkerFault::Standard, 1, true, true)
            };
            if (size > 2) cases.emplace_back(region, WorkerFault::None, 2);
            if (region == Region::First) cases.emplace_back(region, WorkerFault::BadAlloc);
            for (const auto& scenario : cases) {
                RunScenario(scenario, size, n, count);
                ++scenarios;
                if (scenario.launchFailure != 0) ++blockedJoins;
            }
        }
        Observe();
        Check(PrimeCounting::PiLegendre(100000000, 1, 1) == 5761455,
              "serial count did not recover after fault checks");
        std::cout << "PASS regions=4 scenarios=" << scenarios
                  << " blocked_joins=" << blockedJoins
                  << " partial_final_workers=" << sizes[2]
                  << " recovery_counts=" << scenarios + 4 << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
