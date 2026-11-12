#pragma once

// Counts-only Phi instrumentation for the matched-cutoff investigation. No
// clock reads, containers, shared diagnostic writes, or diagnostic locks occur
// in Phi. The caller owns one Stats object per worker and range, publishes it
// through Scope, and reads/merges it only after that worker has joined.
// Counts change execution cost and scheduling: never use this build for times.
#include <array>
#include <cstddef>
#include <cstdint>
#include <ostream>

namespace MatchedPhi {
    enum Counter {
        PublicPhi, PublicZero, PublicAZero, PublicOne, PublicTiny, PublicPi,
        KernelEntries, RecursiveEntries, KernelTiny, Expansions, TinyCalls,
        PiSuffixTerms, CacheGets, CacheEligible, CacheIneligible, CacheHits,
        CacheMisses, CacheNoPage, CacheEmptyCell, TopCacheGets, TopCacheEligible,
        TopCacheHits, TopCacheMisses, PutCalls, PutEligible, LockAttempts,
        PageAllocAttempts, PageAllocSuccess, PageAllocFailed, CacheStores,
        COUNTERS
    };

    inline const char* name(std::size_t i) {
        static const char* const names[COUNTERS] = {
            "public_phi", "public_zero", "public_a_zero", "public_one",
            "public_tiny", "public_pi", "kernel_entries", "recursive_entries",
            "kernel_tiny", "expansions", "tiny_calls", "pi_suffix_terms",
            "cache_gets", "cache_eligible", "cache_ineligible", "cache_hits",
            "cache_misses", "cache_no_page", "cache_empty_cell", "top_cache_gets",
            "top_cache_eligible", "top_cache_hits", "top_cache_misses", "put_calls",
            "put_eligible", "lock_attempts", "page_alloc_attempts",
            "page_alloc_success", "page_alloc_failed", "cache_stores"
        };
        return names[i];
    }

    struct Stats {
        std::array<std::uint64_t, COUNTERS> counts = {{0}};
        void add(const Stats& other) {
            for (std::size_t i = 0; i < COUNTERS; ++i) counts[i] += other.counts[i];
        }
        void writeJson(std::ostream& out) const {
            out << '{';
            for (std::size_t i = 0; i < COUNTERS; ++i) {
                if (i) out << ',';
                out << '"' << name(i) << "\":" << counts[i];
            }
            out << '}';
        }
        void writeCsv(std::ostream& out) const {
            for (std::size_t i = 0; i < COUNTERS; ++i) {
                if (i) out << ',';
                out << counts[i];
            }
        }
        // For successful complete Phi calls only; exceptions may interrupt
        // expansion/public-call identities. These are accounting invariants,
        // not independent checks of the mathematical Phi result.
        bool consistent() const {
            const auto& c = counts;
            return c[CacheGets] == c[CacheEligible] + c[CacheIneligible] &&
                c[CacheEligible] == c[CacheHits] + c[CacheMisses] &&
                c[CacheMisses] == c[CacheNoPage] + c[CacheEmptyCell] &&
                c[KernelEntries] == c[KernelTiny] + c[CacheGets] &&
                c[Expansions] == c[CacheMisses] + c[CacheIneligible] &&
                c[PutCalls] == c[Expansions] &&
                c[TinyCalls] == c[Expansions] + c[KernelTiny] + c[PublicTiny] &&
                c[KernelEntries] == c[RecursiveEntries] + c[PublicPhi] -
                    c[PublicZero] - c[PublicAZero] - c[PublicOne] -
                    c[PublicTiny] - c[PublicPi] &&
                c[TopCacheEligible] == c[TopCacheHits] + c[TopCacheMisses] &&
                c[PageAllocAttempts] == c[PageAllocSuccess] + c[PageAllocFailed] &&
                c[CacheStores] <= c[PutEligible] &&
                c[PutEligible] <= c[PutCalls];
        }
    };

    inline void writeCsvHeader(std::ostream& out) {
        for (std::size_t i = 0; i < COUNTERS; ++i) {
            if (i) out << ',';
            out << name(i);
        }
    }

    // Exactly one definition lives in the generated Phi.cpp. Scope can wrap
    // one public call or several calls in the same phase. Do not redirect it
    // during a recursive evaluation or share a Stats between running workers.
    extern thread_local Stats* active;
    extern thread_local unsigned depth;

    struct Scope {
        Stats* previous;
        unsigned previousDepth;
        explicit Scope(Stats& destination)
            : previous(active), previousDepth(depth) {
            active = &destination;
            depth = 0;
        }
        ~Scope() { active = previous; depth = previousDepth; }
        Scope(const Scope&) = delete;
        Scope& operator=(const Scope&) = delete;
    };

    inline void add(Counter counter, std::uint64_t count = 1) {
        if (active) active->counts[counter] += count;
    }
    struct Evaluation {
        bool enabled;
        Evaluation() : enabled(active != nullptr) {
            if (!enabled) return;
            add(KernelEntries);
            if (depth++) add(RecursiveEntries);
        }
        ~Evaluation() { if (enabled) --depth; }
    };
    inline void cacheGet(bool eligible, bool hit, bool page) {
        if (!active) return;
        add(CacheGets);
        const bool top = depth == 1;
        if (top) add(TopCacheGets);
        if (!eligible) { add(CacheIneligible); return; }
        add(CacheEligible);
        if (top) add(TopCacheEligible);
        if (hit) {
            add(CacheHits);
            if (top) add(TopCacheHits);
        } else {
            add(CacheMisses);
            add(page ? CacheEmptyCell : CacheNoPage);
            if (top) add(TopCacheMisses);
        }
    }
}
