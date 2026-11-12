#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

// Frozen experimental v3 policy. Retain v2 except that small parallel
// candidates passing its worker-budget test must scan at most 256 integers
// per streamed query. Short complete ranges still bypass these tests. The
// density bound is an empirical trial choice, not a mathematical crossover.
// Do not retune this version after inspecting its benchmark results.
namespace P2SelectorV3 {
    enum class Reason {
        AtMostOneQuery, ShortFullRange, TooFewQueries, TooFewStreamQueries,
        StreamWorkerDeficit, ParallelCoverage, Range8, Range16, Range32,
        Range64, ParallelPromotion, StreamTooSparse
    };

    struct Decision {
        std::int64_t limit;
        std::int64_t streamFirst;
        int multiplier; // Zero means Phi-only, with limit equal to L.
        Reason reason;
        std::int64_t nonTableQueries;
        std::int64_t phiQueries;
        std::int64_t streamQueries;

        // Preserve the proposed work even when a later guard rejects it.
        // The <=1-query guard does not evaluate a candidate: its limit is L,
        // first is last, and queries/span/budget are zero. candidateSpan ends
        // at the last actual target, which can be below candidateLimit.
        std::int64_t candidateLimit;
        std::int64_t candidateStreamFirst;
        std::int64_t candidateStreamQueries;
        std::int64_t candidateSpan;
        std::size_t candidateStreamBudget;
    };

    inline const char* reasonName(Reason reason) {
        switch (reason) {
            case Reason::AtMostOneQuery: return "at-most-one-query";
            case Reason::ShortFullRange: return "short-full-range";
            case Reason::TooFewQueries: return "too-few-queries";
            case Reason::TooFewStreamQueries: return "too-few-stream-queries";
            case Reason::StreamWorkerDeficit: return "stream-worker-deficit";
            case Reason::ParallelCoverage: return "parallel-coverage";
            case Reason::Range8: return "range-8";
            case Reason::Range16: return "range-16";
            case Reason::Range32: return "range-32";
            case Reason::Range64: return "range-64";
            case Reason::ParallelPromotion: return "parallel-promotion";
            case Reason::StreamTooSparse: return "stream-too-sparse";
        }
        return "unknown";
    }

    // Preconditions match P2 after input/coverage validation: positive x,
    // prepared limit L, one-based primes, and [first,last) containing exactly
    // the non-table queries. phiBudget >= 1 is computed from the full range
    // using production's 2048-term grain. nThreads is normalized here exactly
    // as in the stream launcher. Empty ranges never access primes[first].
    // Selection inspects no cache or timing state and allocates no storage.
    inline Decision choose(std::int64_t x, std::int64_t L,
                           std::int64_t first, std::int64_t last,
                           std::size_t phiBudget, int nThreads,
                           const std::vector<std::int64_t>& primes) {
        Decision result = {L, last, 0, Reason::AtMostOneQuery,
                           last - first, last - first, 0,
                           L, last, 0, 0, 0};
        if (result.nonTableQueries <= 1) return result;

        int multiplier;
        Reason reason;
        if (x < 100000000000LL) {
            multiplier = 8;
            reason = Reason::Range8;
        } else if (x < 3000000000000LL) {
            multiplier = 16;
            reason = Reason::Range16;
        } else if (x < 30000000000000LL) {
            multiplier = 32;
            reason = Reason::Range32;
        } else {
            multiplier = 64;
            reason = Reason::Range64;
        }

        if (multiplier == 16 && x >= 1000000000000LL && phiBudget >= 4) {
            multiplier = 32;
            reason = Reason::ParallelPromotion;
        }

        const bool smallParallel = x < 100000000000LL && phiBudget > 1;
        if (smallParallel) {
            multiplier = 16;
            reason = Reason::ParallelCoverage;
        }

        const std::int64_t largestQuotient = x / primes[first];
        // Use the existing overflow-safe cap and exact partition identity:
        // floor(x/p) <= Q iff p > floor(x/(Q+1)).
        result.candidateLimit = L > largestQuotient / multiplier ?
            largestQuotient : L * multiplier;
        result.candidateStreamFirst = std::upper_bound(
            primes.begin() + first, primes.begin() + last,
            x / (result.candidateLimit + 1)) - primes.begin();
        result.candidateStreamQueries = last - result.candidateStreamFirst;
        if (result.candidateStreamQueries > 0) {
            result.candidateSpan = x / primes[result.candidateStreamFirst] - L;
            // Match the current launcher's numeric-span grain. This is the
            // worker budget; empty quotient bands can launch fewer workers.
            constexpr std::int64_t MIN_WORK_PER_THREAD = 1048576;
            result.candidateStreamBudget = static_cast<std::size_t>(
                std::min<std::int64_t>(std::max(1, nThreads),
                    std::min(result.candidateStreamQueries,
                        std::max<std::int64_t>(
                            1, result.candidateSpan / MIN_WORK_PER_THREAD))));
        }

        // A short complete range avoids a separate Phi prefix and scans at
        // most two standard 32768-odd-number segments. It takes precedence
        // over the query-count guards and small-parallel budget comparison.
        if (result.candidateStreamFirst == first && result.candidateSpan <= 131072) {
            reason = Reason::ShortFullRange;
        } else if (result.nonTableQueries < 4096) {
            result.reason = Reason::TooFewQueries;
            return result;
        } else if (result.candidateStreamQueries < 4096) {
            result.reason = Reason::TooFewStreamQueries;
            return result;
        } else if (smallParallel && result.candidateStreamBudget < phiBudget) {
            result.reason = Reason::StreamWorkerDeficit;
            return result;
        } else if (smallParallel &&
                   result.candidateSpan / 256 + (result.candidateSpan % 256 != 0) >
                       result.candidateStreamQueries) {
            // Exactly span <= 256 * queries, without multiplying a count or
            // adding 255 to the span. Use the actual streamed endpoint, not
            // the cap. This guard changes only the small-parallel admission;
            // the short complete-range exception above still takes priority.
            result.reason = Reason::StreamTooSparse;
            return result;
        }

        result.limit = result.candidateLimit;
        result.streamFirst = result.candidateStreamFirst;
        result.multiplier = multiplier;
        result.reason = reason;
        result.phiQueries = result.streamFirst - first;
        result.streamQueries = result.candidateStreamQueries;
        return result;
    }
}
