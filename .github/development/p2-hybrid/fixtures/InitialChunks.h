#pragma once

// Temporary P2 experiment. One ticket queue covers the entire non-table range;
// completing the initial small chunks is NOT a prerequisite for later claims.
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <stdexcept>

namespace InitialChunks {
    template <std::int64_t INITIAL>
    class Queue {
        static_assert(INITIAL == 1 || INITIAL == 2 || INITIAL == 4 ||
                      INITIAL == 8 || INITIAL == 32, "unsupported initial chunk");
        const std::int64_t first_, last_, headTerms_;
        const std::uint64_t headTasks_, tasks_;
        std::atomic<std::uint64_t> next_;

        static std::int64_t length(std::int64_t first, std::int64_t last) {
            if (first < 0 || last < first)
                throw std::out_of_range("invalid initial-chunk interval");
            return last - first;
        }
        static std::uint64_t tasks(std::int64_t terms, std::int64_t grain) {
            // Division followed by a remainder test avoids terms+grain-1.
            return static_cast<std::uint64_t>(terms / grain + (terms % grain != 0));
        }

    public:
        Queue(std::int64_t first, std::int64_t last)
            : first_(first), last_(last),
              headTerms_(std::min<std::int64_t>(32, length(first, last))),
              headTasks_(tasks(headTerms_, INITIAL)),
              tasks_(headTasks_ + tasks(last - first - headTerms_, 32)), next_(0) {}

        std::uint64_t size() const { return tasks_; }

        // Pure ticket-to-range mapping, also exposed for allocation-free edge
        // checks. A valid offset is strictly below last-first, so adding first
        // cannot overflow even for a range ending at INT64_MAX.
        bool range(std::uint64_t ticket, std::int64_t& begin,
                   std::int64_t& end) const {
            if (ticket >= tasks_) return false;
            if (ticket < headTasks_) {
                const auto offset = static_cast<std::int64_t>(ticket) * INITIAL;
                begin = first_ + offset;
                end = begin + std::min(INITIAL, headTerms_ - offset);
            } else {
                const auto offset = headTerms_ +
                    static_cast<std::int64_t>(ticket - headTasks_) * 32;
                begin = first_ + offset;
                end = begin + std::min<std::int64_t>(32, last_ - begin);
            }
            return true;
        }

        bool claim(std::int64_t& begin, std::int64_t& end) {
            // This counter assigns unique work only; publication and result
            // visibility still come from ThreadBatch launch/join. The kernel
            // makes one failed terminal claim per worker at most.
            return range(next_.fetch_add(1, std::memory_order_relaxed), begin, end);
        }
    };
}
