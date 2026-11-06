#pragma once

#include "NumbersUtils/Phi.h"
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace PrimeCounting {
namespace PartialSieveDetail {

    constexpr std::int64_t MAX_X = (std::int64_t(1) << 53) - 1;

    inline void Validate(std::int64_t x, std::int64_t a,
                         const PhiTables& tables) {
        if (x < 0 || x > MAX_X) {
            throw std::out_of_range("Partial sieve x must be between 0 and 2^53 - 1");
        }
        if (a < 0 || a > tables.primeCount()) {
            throw std::out_of_range("Partial sieve a must be between 0 and the prepared prime count");
        }
    }

    // Floating-point estimates are corrected by exact comparisons. Successive
    // divisions avoid overflowing a square or cube at a representation boundary.
    // Call only with nonnegative inputs in the validated 53-bit domain.
    inline std::int64_t IntegerSqrt(std::int64_t x) {
        if (x == 0) return 0;
        std::int64_t root = static_cast<std::int64_t>(std::sqrt(static_cast<double>(x)));
        while (root > x / root) --root;
        while (root + 1 <= x / (root + 1)) ++root;
        return root;
    }

    inline std::int64_t IntegerCbrt(std::int64_t x) {
        if (x == 0) return 0;
        std::int64_t root = static_cast<std::int64_t>(std::cbrt(static_cast<double>(x)));
        if (root < 1) root = 1;
        while (root > x / root / root) --root;
        while (root + 1 <= x / (root + 1) / (root + 1)) ++root;
        return root;
    }
}
}
