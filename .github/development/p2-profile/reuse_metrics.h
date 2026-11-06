#pragma once

// Serial diagnostic tracing only. No memoization, timers, result substitution,
// or production hooks. Inclusive costs cover all evaluator calls/direct suffix
// terms, including descendants that are not themselves eligible for recording.
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace P2ReuseMetrics {
    constexpr std::uint32_t UNSET = std::numeric_limits<std::uint32_t>::max();

    struct Record {
        std::uint64_t x;
        std::uint32_t a, query, end, kernel, suffix, keyId;
    };
    static_assert(sizeof(Record) == 32, "reuse trace records must occupy 32 bytes");

    struct State {
        std::vector<Record> records;
        std::uint64_t totalKernels = 0, totalSuffix = 0, queries = 0;
        std::int64_t first = 0, last = 0;
        std::size_t maxEvents = 0;
        std::uint32_t currentQuery = UNSET;
        bool enabled = false, hasRange = false;
        const char* error = nullptr;
    };

    // External inline function-local static: one shared object across the
    // diagnostic translation units. The instrumented P2 rejects worker use
    // while enabled; the driver must keep the entire trace invocation serial.
    inline State& state() { static State value; return value; }
    inline void fail(const char* message) noexcept {
        if (!state().error) state().error = message;
    }
    inline void begin(std::size_t maxEvents) {
        if (maxEvents == 0 || maxEvents > UNSET) {
            throw std::out_of_range("reuse maxEvents must fit a positive uint32_t");
        }
        state() = State();
        state().maxEvents = maxEvents;
        state().records.reserve(maxEvents);
        state().enabled = true;
    }
    inline void end() {
        auto& value = state();
        value.enabled = false;
        if (value.hasRange && value.queries !=
            static_cast<std::uint64_t>(value.last - value.first)) {
            fail("reuse trace did not visit the complete P2 query range");
        }
        if (value.error) throw std::runtime_error(value.error);
    }
    inline void require_serial(int nThreads) {
        if (state().enabled && nThreads > 1) {
            fail("reuse tracing requires serial P2");
            throw std::invalid_argument(state().error);
        }
    }
    inline void range(std::int64_t first, std::int64_t last) noexcept {
        auto& value = state();
        if (!value.enabled) return;
        if (value.hasRange || first < 0 || last < first ||
            static_cast<std::uint64_t>(last) > UNSET) {
            fail("reuse trace requires one bounded serial P2 range");
            return;
        }
        value.hasRange = true;
        value.first = first; value.last = last;
    }
    inline void query(std::int64_t index) noexcept {
        auto& value = state();
        if (!value.enabled) return;
        if (!value.hasRange || index < value.first || index >= value.last ||
            index != value.first + static_cast<std::int64_t>(value.queries)) {
            fail("reuse P2 queries are not a contiguous increasing range");
            return;
        }
        value.currentQuery = static_cast<std::uint32_t>(index);
        ++value.queries;
    }
    inline void suffix_terms(std::int64_t count) noexcept {
        auto& value = state();
        if (!value.enabled) return;
        if (count < 0 || static_cast<std::uint64_t>(count) >
            std::numeric_limits<std::uint64_t>::max() - value.totalSuffix) {
            fail("reuse suffix counter overflow");
            return;
        }
        value.totalSuffix += static_cast<std::uint64_t>(count);
    }

    class Evaluation {
        std::uint32_t recordIndex = UNSET;
        std::uint64_t kernelBefore = 0, suffixBefore = 0;

    public:
        Evaluation(std::int64_t x, std::int64_t a) {
            auto& value = state();
            if (!value.enabled) return;
            kernelBefore = value.totalKernels;
            suffixBefore = value.totalSuffix;
            if (value.totalKernels == std::numeric_limits<std::uint64_t>::max()) {
                fail("reuse evaluator counter overflow");
                return;
            }
            ++value.totalKernels;

            // a == 6 returns the tiny-table result. Every other ineligible
            // cache request necessarily expands; record before its children.
            // Cache eligibility constants are checked by instrument_reuse.py.
            if (a <= 6 || (a < 100 && x <= 65535) || value.error) return;
            if (x < 0 || a < 0 || static_cast<std::uint64_t>(a) > UNSET ||
                value.currentQuery == UNSET) {
                fail("reuse node fields exceed the trace representation");
                return;
            }
            if (value.records.size() >= value.maxEvents) {
                // Keep computing/counting, but end() rejects this incomplete
                // trace. No new allocation or exception occurs in the kernel.
                fail("reuse trace exceeded maxEvents");
                return;
            }
            recordIndex = static_cast<std::uint32_t>(value.records.size());
            value.records.push_back({static_cast<std::uint64_t>(x),
                static_cast<std::uint32_t>(a), value.currentQuery, 0, 0, 0, UNSET});
        }
        Evaluation(const Evaluation&) = delete;
        Evaluation& operator=(const Evaluation&) = delete;
        ~Evaluation() noexcept {
            if (recordIndex == UNSET) return;
            auto& value = state();
            const auto kernel = value.totalKernels - kernelBefore;
            const auto suffix = value.totalSuffix - suffixBefore;
            if (value.records.size() > UNSET || kernel > UNSET || suffix > UNSET ||
                recordIndex >= value.records.size()) {
                fail("reuse subtree fields exceed uint32_t");
                return;
            }
            // Save an index, never a pointer/reference that a vector growth
            // could invalidate. end is the exclusive preorder subtree end.
            auto& record = value.records[recordIndex];
            record.end = static_cast<std::uint32_t>(value.records.size());
            record.kernel = static_cast<std::uint32_t>(kernel);
            record.suffix = static_cast<std::uint32_t>(suffix);
        }
    };
}
