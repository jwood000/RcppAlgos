#include "NumbersUtils/NthPrime.h"
#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "CppConvert/Constants.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace PrimeCounting {
    namespace {

        using PartialSieveDetail::IntegerSqrt;

        constexpr std::int64_t SIEVE_BLOCK_SIZE = 65536;

        void ensureTables(
            std::int64_t maximum, std::unique_ptr<PhiTables>& tables
        ) {

            const std::int64_t root = IntegerSqrt(maximum);

            if (!tables || tables->limit() < root) {
                // No evaluator or sieve holds references across this call. Release
                // the old table first to avoid holding two dense pi arrays at once.
                tables.reset();
                tables.reset(new PhiTables(root));
            }
        }

        std::int64_t clampEstimate(double estimate) {
            // Clamp before conversion. Floating-point estimates only choose search
            // locations; exact counts and sieving determine the result.
            return static_cast<std::int64_t>(
                std::max(
                    2.0,
                    std::min(
                        Significand53,
                        std::round(estimate)
                    )
                )
            );
        }

        // anchor is an integer in [2, Significand53Int], count == pi(anchor), and n is a
        // supported rank. Forward selection excludes anchor; backward selection
        // includes it. In particular, count == n still requires finding the last
        // prime <= anchor: the anchor itself may be composite.
        std::int64_t sieveNthPrime(std::int64_t n, std::int64_t anchor,
                                   std::int64_t count,
                                   std::unique_ptr<PhiTables>& tables) {

            const bool forward = count < n;
            std::int64_t remaining = forward ? n - count : count - n + 1;
            std::int64_t edge = forward ? anchor + 1 : anchor;
            std::vector<unsigned char> sieve(SIEVE_BLOCK_SIZE);

            while (edge >= 2 && edge <= Significand53Int) {
                const std::int64_t low = forward ? edge :
                    std::max<std::int64_t>(2, edge - SIEVE_BLOCK_SIZE + 1);

                const std::int64_t high = forward ?
                    std::min(Significand53Int, edge + SIEVE_BLOCK_SIZE - 1) : edge;

                const std::int64_t length = high - low + 1;
                ensureTables(high, tables);
                std::fill(sieve.begin(), sieve.begin() + length, 1);

                // Simple bounded segmented sieve. Starting at p^2 preserves p
                // itself; primes through floor(sqrt(high)) mark every composite.
                const auto& primes = tables->primes();

                for (std::size_t i = 1; i < primes.size(); ++i) {
                    const std::int64_t p = primes[i];

                    if (p > high / p) {
                        break;
                    }

                    const std::int64_t first = std::max(
                        p * p,
                        low + (p - low % p) % p
                    );

                    for (std::int64_t value = first; value <= high; value += p) {
                        sieve[value - low] = 0;
                    }
                }

                if (forward) {
                    for (std::int64_t i = 0; i < length; ++i) {
                        if (sieve[i] && --remaining == 0) {
                            return low + i;
                        }
                    }

                    edge = high + 1;
                } else {
                    for (std::int64_t i = length - 1; i >= 0; --i) {
                        if (sieve[i] && --remaining == 0) {
                            return low + i;
                        }
                    }

                    edge = low - 1;
                }
            }

            // Valid ranks and an exact anchor count must reach their target.
            throw std::logic_error(
                "nthPrime exhausted the supported search range"
            );
        }

    } // namespace

    std::int64_t NthPrime(std::int64_t n, int nThreads) {
        if (n < 1 || n > MAX_NTH_PRIME_INDEX) {
            throw std::out_of_range(
                "nthPrime n must be between 1 and 252252704148404"
            );
        }

        const std::int64_t smallPrimes[] = {
            0, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
            41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97
        };

        if (n < static_cast<std::int64_t>(sizeof(smallPrimes) / sizeof(*smallPrimes))) {
            return smallPrimes[n];
        }

        const double logN = std::log(static_cast<double>(n));
        const double logLogN = std::log(logN);

        std::int64_t anchor = clampEstimate(
            n * (logN + logLogN - 1 + (logLogN - 2) / logN)
        );

        // n(log n + log log n) plans table capacity. Although it is a known
        // upper bound for n >= 6 in real arithmetic, correctness does not rely
        // on its floating-point evaluation: ensureTables checks every endpoint.
        std::unique_ptr<PhiTables> tables;
        ensureTables(
            std::max(anchor, clampEstimate(n * (logN + logLogN))), tables
        );

        std::int64_t lower = 2;
        std::int64_t upper = Significand53Int;
        std::int64_t count;
        nThreads = std::max(1, nThreads);

        for (int correction = 0; ; ++correction) {
            ensureTables(anchor, tables);
            count = PiAuto(anchor, nThreads, *tables);

            // Density near x is approximately 1/log(x). This estimates the
            // remaining sieve distance, not a guaranteed error bound. Permit
            // at most two corrections; the exact sieve handles any remainder.
            const double distance = static_cast<double>(n - count) *
                std::log(static_cast<double>(anchor)
            );

            if (std::abs(distance) <= SIEVE_BLOCK_SIZE || correction == 2) {
                break;
            }

            // Exact counts maintain lower <= p_n <= upper. These clamps keep
            // a heuristic correction within the remaining search interval.
            if (count < n) {
                lower = anchor + 1;
            } else {
                upper = anchor;
            }

            const std::int64_t next = std::max(
                lower,
                std::min(
                    upper,
                    clampEstimate(anchor + distance)
                )
            );

            if (next == anchor) {
                break;
            }

            anchor = next;
        }

        return sieveNthPrime(n, anchor, count, tables);
    }
}
