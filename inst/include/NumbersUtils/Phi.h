#pragma once

#include <cstdint>
#include <vector>

namespace PrimeCounting {

    // Canonical primes and pi values through an inclusive limit. Preparation
    // is separate from evaluation so callers can reuse these tables.
    class PhiTables {
        std::vector<std::int64_t> primes_;
        std::vector<std::int32_t> pi_;

        friend std::int64_t Phi(std::int64_t x, std::int64_t a,
                                const PhiTables& tables);

    public:
        // Initial preparation range: 0 <= limit <= 94906265, matching the
        // largest root table needed for the package's 53-bit input range.
        explicit PhiTables(std::int64_t limit);

        std::int64_t limit() const;
        std::int64_t primeCount() const;

        // Read-only, one-based prime sequence; element zero is a sentinel.
        const std::vector<std::int64_t>& primes() const { return primes_; }
    };

    // Count integers in [1, x] not divisible by the first a primes.
    // Requires 0 <= x <= 2^53 - 1 and 0 <= a <= tables.primeCount().
    // x may exceed tables.limit(); smaller tables only limit the available a
    // and lookup shortcuts. Invalid inputs throw std::out_of_range.
    // Evaluation is serial. Concurrent calls may share immutable tables;
    // the bounded cache is shared safely across calls and table instances.
    std::int64_t Phi(std::int64_t x, std::int64_t a, const PhiTables& tables);
}
