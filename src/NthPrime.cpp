#include "NumbersUtils/NthPrime.h"
#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/Eratosthenes.h"
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

        // Bound retained prime values in the final selection, independently of
        // the larger count-only intervals used for corrections.
        constexpr std::int64_t FINAL_SIEVE_LIMIT = 1000000;

        std::int64_t sieveLimit(std::int64_t x) {
            // The search location, not the requested rank, sets the allowance.
            // These broad bands are a work policy, not an error bound.
            if (x < 10000000000LL) return 1000000;
            if (x < 1000000000000LL) return 10000000;
            if (x < 100000000000000LL) return 100000000;
            return 1000000000;
        }

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

        double estimatedDistance(std::int64_t n, std::int64_t anchor,
                                 std::int64_t count) {
            return static_cast<double>(n - count) *
                std::log(static_cast<double>(anchor));
        }

        std::int64_t nextAnchor(std::int64_t n, std::int64_t anchor,
                                std::int64_t count, std::int64_t& lower,
                                std::int64_t& upper) {
            // Exact counts maintain lower <= p_n <= upper. The density
            // estimate only chooses a new location inside that bracket.
            if (count < n) {
                lower = anchor + 1;
            } else {
                upper = anchor;
            }

            return std::max(
                lower,
                std::min(
                    upper,
                    clampEstimate(anchor + estimatedDistance(n, anchor, count))
                )
            );
        }

        std::int64_t countAt(std::int64_t next, std::int64_t anchor,
                             std::int64_t count, int nThreads,
                             int maxThreads, int maxCores) {

            // count == pi(anchor). The disjoint interval supplies pi(next)
            // exactly, including when either endpoint is itself prime.
            if (next > anchor) {
                return count + PiPrimeRange(anchor + 1, next,
                                            nThreads, maxThreads, maxCores);
            }

            if (next < anchor) {
                return count - PiPrimeRange(next + 1, anchor,
                                            nThreads, maxThreads, maxCores);
            }

            return count;
        }

        // Forward selection excludes anchor; backward selection includes it.
        // Equal counts still require finding the last prime <= anchor, since
        // the anchor may lie on a composite plateau.
        std::int64_t sieveNthPrime(std::int64_t n, std::int64_t anchor,
                                   std::int64_t count, int nThreads,
                                   int maxThreads, int maxCores) {

            const bool forward = count < n;
            std::int64_t remaining = forward ? n - count : count - n + 1;
            std::int64_t edge = forward ? anchor + 1 : anchor;

            while (edge >= 2 && edge <= Significand53Int) {
                // Leave room around the density estimate, but retain at most
                // one bounded window. A short estimate only causes another
                // adjacent window; it cannot change the selected rank.
                const double estimate = 2.0 * static_cast<double>(remaining) *
                    std::log(static_cast<double>(edge));

                const std::int64_t width = static_cast<std::int64_t>(
                    std::max(
                        1.0,
                        std::min(
                            static_cast<double>(FINAL_SIEVE_LIMIT),
                            std::ceil(estimate)
                        )
                    )
                );

                const std::int64_t low = forward ? edge :
                    std::max<std::int64_t>(2, edge - width + 1);

                const std::int64_t high = forward ?
                    std::min(Significand53Int, edge + width - 1) : edge;

                std::vector<std::vector<std::int64_t>> primeList;
                std::vector<std::int64_t> primes;
                bool parallel = false;

                PrimeSieve::PrimeSieveMain(
                    primeList, primes, low, high,
                    parallel, nThreads, maxThreads, maxCores
                );

                // Rows are ordered intervals. Walk them directly instead of
                // flattening or storing primes from the correction interval.
                if (!parallel) {
                    const std::int64_t size = primes.size();

                    if (remaining <= size) {
                        return primes[forward ? remaining - 1 : size - remaining];
                    }

                    remaining -= size;
                } else if (forward) {
                    for (const auto& row : primeList) {
                        const std::int64_t size = row.size();

                        if (remaining <= size) {
                            return row[remaining - 1];
                        }

                        remaining -= size;
                    }
                } else {
                    for (auto row = primeList.crbegin(); row != primeList.crend(); ++row) {
                        const std::int64_t size = row->size();

                        if (remaining <= size) {
                            return (*row)[size - remaining];
                        }

                        remaining -= size;
                    }
                }

                edge = forward ? high + 1 : low - 1;
            }

            // Valid ranks and an exact anchor count must reach their target.
            throw std::logic_error("nthPrime exhausted the supported search range");
        }

    } // namespace

    std::int64_t NthPrime(std::int64_t n, int nThreads,
                           int maxThreads, int maxCores) {
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

        maxThreads = std::max(1, maxThreads);
        maxCores = std::max(1, maxCores);
        nThreads = std::max(1, std::min(nThreads, maxThreads));

        // At most two whole-prefix counts. Stop after the first when the
        // estimated correction fits the allowance for this search location.
        for (int fullCounts = 0; ; ++fullCounts) {
            ensureTables(anchor, tables);
            count = PiAuto(anchor, nThreads, *tables);

            if (fullCounts == 1 ||
                std::abs(estimatedDistance(n, anchor, count)) <= sieveLimit(anchor)) {
                break;
            }

            const std::int64_t next = nextAnchor(n, anchor, count, lower, upper);

            if (next == anchor) {
                break;
            }

            anchor = next;
        }

        // Native sieving prepares its own base primes. No worker or evaluator
        // still refers to these dense tables after PiAuto returns.
        tables.reset();

        // Refine large residuals with interval counts, not new prefix counts.
        // Limit repeated sieve preparation; exact selection handles any rest.
        for (int correction = 0; correction < 2; ++correction) {
            if (std::abs(estimatedDistance(n, anchor, count)) <= FINAL_SIEVE_LIMIT) {
                break;
            }

            const std::int64_t next = nextAnchor(n, anchor, count, lower, upper);

            if (next == anchor) {
                break;
            }

            count = countAt(next, anchor, count, nThreads, maxThreads, maxCores);
            anchor = next;
        }

        return sieveNthPrime(n, anchor, count, nThreads, maxThreads, maxCores);
    }

    std::int64_t NthPrime(std::int64_t n, int nThreads) {
        // Existing native callers may supply only a usable worker budget.
        // R supplies both hardware limits through the full overload.
        const int budget = std::max(1, nThreads);
        return NthPrime(n, budget, budget, budget);
    }
}
