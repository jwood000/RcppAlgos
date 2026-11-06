#pragma once

// Diagnostic-only collection. Production translation units do not include this
// header. Timing mode deliberately removes recursive counters; compare its
// total time with an uninstrumented build to assess clock/probe overhead.
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <vector>

namespace P2ProfileMetrics {
    constexpr std::size_t BINS = 8;
    enum Counter {
        Queries, CompletedQueries, PublicPhi, KernelEntries, RecursiveEntries,
        Expansions, TinyCalls, ExitZero, ExitAZero, ExitOne, ExitTiny, ExitPi,
        PiSuffixTerms, CacheGets, CacheEligible, CacheIneligible, CacheHits,
        CacheMisses, CacheNoPage, CacheEmptyCell, TopCacheGets, TopCacheEligible,
        TopCacheHits, TopCacheMisses, PutCalls, PutEligible, LockAttempts,
        PageAllocAttempts, PageAllocSuccess, PageAllocFailed, CacheStores, COUNTERS
    };
    using Clock = std::chrono::steady_clock;
    using Counts = std::array<std::uint64_t, COUNTERS>;

    struct Bin {
        Counts counts = {{0}};
        std::uint64_t rootLookupNs = 0, phiNs = 0;
        void add(const Bin& other) {
            for (std::size_t i = 0; i < COUNTERS; ++i) counts[i] += other.counts[i];
            rootLookupNs += other.rootLookupNs;
            phiNs += other.phiNs;
        }
    };
    struct WorkerData {
        std::array<Bin, BINS> bins;
        std::uint64_t calls = 0, chunks = 0, wallNs = 0, maxChunkNs = 0;
        std::int64_t maxChunkFirst = 0, maxChunkLast = 0;
        void add(const WorkerData& other) {
            for (std::size_t i = 0; i < BINS; ++i) bins[i].add(other.bins[i]);
            calls += other.calls; chunks += other.chunks; wallNs += other.wallNs;
            if (other.maxChunkNs > maxChunkNs) {
                maxChunkNs = other.maxChunkNs;
                maxChunkFirst = other.maxChunkFirst;
                maxChunkLast = other.maxChunkLast;
            }
        }
        Bin total() const {
            Bin result;
            for (const auto& bin : bins) result.add(bin);
            return result;
        }
    };
    struct ThreadState {
        bool active = false;
        std::size_t bin = 0;
        unsigned depth = 0;
        std::int64_t first = 0, last = 0;
        WorkerData data;
    };
    struct Session {
        bool enabled = false, mixedCases = false;
        std::int64_t x = 0, a = 0, first = 0, last = 0;
        int requested = 0;
        std::uint64_t calls = 0, directTerms = 0, directNs = 0;
        std::vector<WorkerData> workers;
    };

    // Function-local statics in these external inline functions are shared
    // across translation units, but each thread owns its own separate state.
    inline ThreadState& local() { static thread_local ThreadState value; return value; }
    inline Session& session() { static thread_local Session value; return value; }
    inline std::uint64_t elapsed(Clock::time_point start, Clock::time_point end) {
        return static_cast<std::uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());
    }
    inline void begin() { session() = Session(); session().enabled = true; }
    inline void end() { session().enabled = false; }

    struct Call {
        Session* output = nullptr;
        std::int64_t first = 0, last = 0;
        Clock::time_point directStart;
        Call(std::int64_t x, std::int64_t a, int threads) {
            if (!session().enabled) return;
            output = &session();
            if (output->calls && (output->x != x || output->a != a ||
                                  output->requested != threads)) output->mixedCases = true;
            output->x = x; output->a = a; output->requested = threads;
            ++output->calls;
        }
        void range(std::int64_t begin, std::int64_t end, std::int64_t b) {
            first = begin; last = end;
            if (!output) return;
            output->first = first; output->last = last;
            output->directTerms += static_cast<std::uint64_t>(b - last + 1);
        }
        void prepare(std::size_t workers) {
            if (output && output->workers.size() < workers) output->workers.resize(workers);
        }
        void start_direct() {
#ifdef P2_PROFILE_TIMING
            if (output) directStart = Clock::now();
#endif
        }
        void end_direct() {
#ifdef P2_PROFILE_TIMING
            if (output) output->directNs += elapsed(directStart, Clock::now());
#endif
        }
    };

    struct Worker {
        WorkerData* destination = nullptr;
        Clock::time_point start;
        Worker(Call& call, std::size_t slot) {
            if (!call.output) return;
            destination = &call.output->workers[slot];
            local() = ThreadState();
            local().active = true;
            local().first = call.first; local().last = call.last;
            local().data.calls = 1;
#ifdef P2_PROFILE_TIMING
            start = Clock::now();
#endif
        }
        ~Worker() {
            if (!destination) return;
#ifdef P2_PROFILE_TIMING
            local().data.wallNs = elapsed(start, Clock::now());
#endif
            // Each completed worker publishes once into its own fixed slot.
            // The vector is prepared before launch; ThreadBatch joins before
            // the caller reads it or starts another P2 invocation.
            destination->add(local().data);
            local().active = false;
        }
    };
    struct Chunk {
        bool active;
        std::int64_t first, last;
        Clock::time_point start;
        Chunk(std::int64_t begin, std::int64_t end)
            : active(local().active), first(begin), last(end) {
            if (!active) return;
            ++local().data.chunks;
#ifdef P2_PROFILE_TIMING
            start = Clock::now();
#endif
        }
        ~Chunk() {
#ifdef P2_PROFILE_TIMING
            if (!active) return;
            const auto ns = elapsed(start, Clock::now());
            if (ns > local().data.maxChunkNs) {
                local().data.maxChunkNs = ns;
                local().data.maxChunkFirst = first;
                local().data.maxChunkLast = last;
            }
#endif
        }
    };
    struct Query {
        bool active;
        Clock::time_point start, afterRoot;
        explicit Query(std::int64_t index) : active(local().active) {
            if (!active) return;
            local().bin = static_cast<std::size_t>(
                (index - local().first) * BINS / (local().last - local().first));
            ++local().data.bins[local().bin].counts[Queries];
#ifdef P2_PROFILE_TIMING
            start = Clock::now();
#endif
        }
        void root_done() {
#ifdef P2_PROFILE_TIMING
            if (active) afterRoot = Clock::now();
#endif
        }
        void phi_done() {
            if (!active) return;
            auto& bin = local().data.bins[local().bin];
            ++bin.counts[CompletedQueries];
#ifdef P2_PROFILE_TIMING
            const auto finish = Clock::now();
            bin.rootLookupNs += elapsed(start, afterRoot);
            bin.phiNs += elapsed(afterRoot, finish);
#endif
        }
    };

    inline void add(Counter counter, std::uint64_t value = 1) {
#ifndef P2_PROFILE_TIMING
        if (local().active) local().data.bins[local().bin].counts[counter] += value;
#else
        (void)counter; (void)value;
#endif
    }
    struct Evaluation {
#ifndef P2_PROFILE_TIMING
        bool active;
        Evaluation() : active(local().active) {
            if (!active) return;
            add(KernelEntries);
            if (local().depth++) add(RecursiveEntries);
        }
        ~Evaluation() { if (active) --local().depth; }
#else
        Evaluation() {}
        ~Evaluation() {}
#endif
    };
    inline void cache_get(bool eligible, bool hit, bool page) {
#ifndef P2_PROFILE_TIMING
        if (!local().active) return;
        add(CacheGets);
        const bool top = local().depth == 1;
        if (top) add(TopCacheGets);
        if (!eligible) { add(CacheIneligible); return; }
        add(CacheEligible);
        if (top) add(TopCacheEligible);
        if (hit) { add(CacheHits); if (top) add(TopCacheHits); }
        else {
            add(CacheMisses);
            add(page ? CacheEmptyCell : CacheNoPage);
            if (top) add(TopCacheMisses);
        }
#else
        (void)eligible; (void)hit; (void)page;
#endif
    }

    inline void write_bin(std::ostream& out, const Bin& bin) {
        static const char* names[COUNTERS] = {
            "queries", "completed_queries", "public_phi", "kernel_entries", "recursive_entries",
            "expansions", "tiny_calls", "exit_zero", "exit_a_zero", "exit_one", "exit_tiny",
            "public_pi_exit", "pi_suffix_terms", "cache_gets", "cache_eligible", "cache_ineligible",
            "cache_hits", "cache_misses", "cache_no_page", "cache_empty_cell", "top_cache_gets",
            "top_cache_eligible", "top_cache_hits", "top_cache_misses", "put_calls", "put_eligible",
            "lock_attempts", "page_alloc_attempts", "page_alloc_success", "page_alloc_failed", "cache_stores"
        };
        out << ",\"root_lookup_ns\":" << bin.rootLookupNs << ",\"phi_ns\":" << bin.phiNs
            << ",\"counters\":{";
        for (std::size_t i = 0; i < COUNTERS; ++i) {
            if (i) out << ',';
            out << '\"' << names[i] << "\":" << bin.counts[i];
        }
        out << '}';
    }
    inline void report(std::ostream& out) {
        const auto& result = session();
        Bin total;
        for (const auto& worker : result.workers) total.add(worker.total());
        out << "{\"type\":\"p2_profile\",\"kind\":\"summary\",\"mode\":\""
#ifdef P2_PROFILE_TIMING
            << "timing"
#else
            << "counts"
#endif
            << "\",\"x\":" << result.x << ",\"a\":" << result.a
            << ",\"requested\":" << result.requested << ",\"calls\":" << result.calls
            << ",\"first\":" << result.first << ",\"last\":" << result.last
            << ",\"direct_terms\":" << result.directTerms << ",\"direct_ns\":" << result.directNs
            << ",\"worker_slots\":" << result.workers.size()
            << ",\"mixed_cases\":" << (result.mixedCases ? "true" : "false");
        write_bin(out, total); out << "}\n";
        for (std::size_t slot = 0; slot < result.workers.size(); ++slot) {
            const auto& worker = result.workers[slot];
            out << "{\"type\":\"p2_profile\",\"kind\":\"worker\",\"slot\":" << slot
                << ",\"calls\":" << worker.calls << ",\"chunks\":" << worker.chunks
                << ",\"wall_ns\":" << worker.wallNs << ",\"max_chunk_ns\":" << worker.maxChunkNs
                << ",\"max_chunk_first\":" << worker.maxChunkFirst
                << ",\"max_chunk_last\":" << worker.maxChunkLast;
            write_bin(out, worker.total()); out << "}\n";
            for (std::size_t bin = 0; bin < BINS; ++bin) {
                if (!worker.bins[bin].counts[Queries]) continue;
                out << "{\"type\":\"p2_profile\",\"kind\":\"quantile\",\"slot\":" << slot
                    << ",\"bin\":" << bin;
                write_bin(out, worker.bins[bin]); out << "}\n";
            }
        }
    }
}
