#pragma once

#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace CacheInspection {

    enum class Phase {
        Setup, Offsets, Buckets, DenseMark, BucketMark, Extract, Reset, Count
    };

    constexpr std::size_t PhaseCount = static_cast<std::size_t>(Phase::Count);
    extern const std::array<const char*, PhaseCount> PhaseNames;

    typedef std::chrono::steady_clock Clock;

    inline double SecondsSince(Clock::time_point start) {
        return std::chrono::duration<double>(Clock::now() - start).count();
    }

    struct KernelRecord {
        std::string kernel;
        int worker = 0;
        std::int_fast64_t lower = 0;
        std::int_fast64_t upper = 0;
        std::size_t segmentSpan = 0;
        std::size_t sieveBytes = 0;
        std::size_t offsetBytes = 0;
        std::size_t bucketEntries = 0;
        std::size_t bucketCapacityBytes = 0;
        std::size_t bufferBytes = 0;
        std::size_t primeCount = 0;
        std::size_t outputCapacityBytes = 0;
        std::size_t segments = 0;
        std::array<double, PhaseCount> seconds{};
        std::array<std::size_t, PhaseCount> calls{};
    };

    struct Inspector {
        // Zero retains the existing segment-size policy.
        std::size_t bigSegmentSpan = 0;
        int requestedThreads = 1;
        int effectiveThreads = 1;
        std::size_t sharedPrimeBytes = 0;
        std::size_t densePrimeBytes = 0;
        std::size_t bucketPrimeBytes = 0;
        std::size_t wheelBytes = 0;
        double baseSeconds = 0;
        double policySeconds = 0;
        double sieveSeconds = 0;
        double resultSeconds = 0;
        // Each worker writes only to its own collection.
        std::vector<std::vector<KernelRecord>> workers;
    };

    class PhaseTimer {
        KernelRecord* record;
        Phase phase;
        Clock::time_point start;

        void Record(Clock::time_point end);

    public:
        explicit PhaseTimer(KernelRecord* record_, Phase phase_ = Phase::Setup) :
            record(record_), phase(phase_) {
            if (record) start = Clock::now();
        }

        PhaseTimer(const PhaseTimer&) = delete;
        PhaseTimer& operator=(const PhaseTimer&) = delete;

        void next(Phase nextPhase) {
            if (record) {
                const Clock::time_point now = Clock::now();
                Record(now);
                phase = nextPhase;
                start = now;
            }
        }

        void stop() {
            if (record) {
                Record(Clock::now());
                record = nullptr;
            }
        }

        ~PhaseTimer() {
            stop();
        }
    };

    KernelRecord* BeginKernel(std::vector<KernelRecord>* records,
                              const char* kernel, std::int_fast64_t lower,
                              std::int_fast64_t upper);
}
