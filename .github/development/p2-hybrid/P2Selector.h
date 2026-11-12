#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

// Frozen experimental v1 policy, chosen from the earlier fixed-multiplier
// anchors before its benchmark matrix is run. Small parallel calls retain
// Phi; larger calls gradually expand the streamed tail. Query-count guards
// avoid paying stream setup for a small tail after arbitrary cutoff/table
// choices. The middle-range parallel promotion reflects the earlier worker
// results. These are empirical policy choices, not mathematical thresholds.
// Do not retune this version after inspecting its results.
namespace P2Selector {
    enum class Reason {
        TooFewQueries, SmallParallel, Range8, Range16, Range32, Range64,
        ParallelPromotion, TooFewStreamQueries
    };

    struct Decision {
        std::int64_t limit;
        std::int64_t streamFirst;
        int multiplier; // Zero means Phi-only, with limit equal to L.
        Reason reason;
        std::int64_t nonTableQueries;
        std::int64_t phiQueries;
        std::int64_t streamQueries;
    };

    inline const char* reasonName(Reason reason) {
        switch (reason) {
            case Reason::TooFewQueries: return "too-few-queries";
            case Reason::SmallParallel: return "small-parallel";
            case Reason::Range8: return "range-8";
            case Reason::Range16: return "range-16";
            case Reason::Range32: return "range-32";
            case Reason::Range64: return "range-64";
            case Reason::ParallelPromotion: return "parallel-promotion";
            case Reason::TooFewStreamQueries: return "too-few-stream-queries";
        }
        return "unknown";
    }

    // Preconditions match P2 after input/coverage validation: x is positive,
    // L is the prepared table limit, primes is one-based, and [first,last)
    // is the non-table prime-index range. phiBudget >= 1 is the unchanged
    // production worker budget computed from that entire range. Empty ranges
    // are accepted without accessing primes[first]. No cache/timing state is
    // inspected, no storage is allocated, and the answer is independent of a
    // method label: a and table coverage affect it through the actual range.
    inline Decision choose(std::int64_t x, std::int64_t L,
                           std::int64_t first, std::int64_t last,
                           std::size_t phiBudget,
                           const std::vector<std::int64_t>& primes) {
        Decision result = {L, last, 0, Reason::TooFewQueries,
                           last - first, last - first, 0};
        if (result.nonTableQueries < 4096) return result;
        if (x < 100000000000LL && phiBudget > 1) {
            result.reason = Reason::SmallParallel;
            return result;
        }

        if (x < 100000000000LL) {
            result.multiplier = 8;
            result.reason = Reason::Range8;
        } else if (x < 3000000000000LL) {
            result.multiplier = 16;
            result.reason = Reason::Range16;
        } else if (x < 30000000000000LL) {
            result.multiplier = 32;
            result.reason = Reason::Range32;
        } else {
            result.multiplier = 64;
            result.reason = Reason::Range64;
        }

        if (result.multiplier == 16 && x >= 1000000000000LL && phiBudget >= 4) {
            result.multiplier = 32;
            result.reason = Reason::ParallelPromotion;
        }

        const std::int64_t largestQuotient = x / primes[first];
        // Preserve the fixed variants' capped, overflow-safe endpoint and
        // exact floor(x/p) <= Q iff p > floor(x/(Q+1)) partition.
        result.limit = L > largestQuotient / result.multiplier ?
            largestQuotient : L * result.multiplier;
        result.streamFirst = std::upper_bound(
            primes.begin() + first, primes.begin() + last,
            x / (result.limit + 1)) - primes.begin();
        result.streamQueries = last - result.streamFirst;
        if (result.streamQueries < 4096) {
            result.limit = L;
            result.streamFirst = last;
            result.multiplier = 0;
            result.reason = Reason::TooFewStreamQueries;
            result.streamQueries = 0;
        }
        result.phiQueries = result.streamFirst - first;
        return result;
    }
}
