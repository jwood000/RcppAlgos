#pragma once

// Optional diagnostics for the experimental P2 implementation. Timing builds
// do not include this header or perform any of these writes. These capacities
// describe payloads, not allocator overhead, thread stacks, or process RSS.
#include <cstdint>
#include <vector>

namespace P2HybridMetrics {
    struct StreamStorage {
        bool active = false;
        bool phiSeed = false;
        std::int64_t start = 0, finish = 0, queries = 0;
        std::uint64_t sieveBytes = 0, nextBytes = 0;
    };

    struct Snapshot {
        std::int64_t x = 0, a = 0, tableLimit = 0, streamLimit = 0;
        std::int64_t phiQueries = 0, streamQueries = 0, directQueries = 0;
        std::uint64_t phiWorkers = 0, streamWorkers = 0;
        std::uint64_t phiSumBytes = 0, streamSumBytes = 0;
        std::vector<StreamStorage> intervals;
    };

    // Call P2 and read its snapshot on the same calling thread. Worker writes
    // target distinct, preallocated interval slots and finish before P2 returns.
    inline Snapshot& current() {
        static thread_local Snapshot value;
        return value;
    }
}
