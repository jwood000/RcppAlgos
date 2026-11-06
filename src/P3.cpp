#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/ThreadBatch.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace PrimeCounting {
namespace {

    std::int64_t countP3Row(std::int64_t quotient, std::int64_t first,
                            std::int64_t last,
                            const std::vector<std::int64_t>& primes,
                            const std::vector<std::int32_t>& pi) {
        std::int64_t sum = 0;

        for (std::int64_t j = first; j <= last; ++j) {
            // k >= j and p[k] <= floor(x / p[i] / p[j]). The subtraction
            // includes equal prime factors and counts each nondecreasing triple
            // i <= j <= k exactly once.
            sum += static_cast<std::int64_t>(pi[quotient / primes[j]]) - j + 1;
        }

        return sum;
    }

    std::int64_t countP3Serial(std::int64_t x, std::int64_t a,
                               std::int64_t b, const PhiTables& tables) {
        const auto& primes = tables.primes();
        const auto& pi = tables.pi();
        std::int64_t sum = 0;

        for (std::int64_t i = a + 1; i <= b; ++i) {
            const std::int64_t quotient = x / primes[i];
            const std::int64_t last = pi[PartialSieveDetail::IntegerSqrt(quotient)];
            sum += countP3Row(quotient, i, last, primes, pi);
        }

        return sum;
    }

    struct P3Row {
        std::int64_t quotient;
        std::int64_t last;
        std::int64_t cumulativeWork;
    };

    std::int64_t countP3Rows(const std::vector<P3Row>& rows,
                             std::size_t begin, std::size_t end,
                             std::int64_t a, const PhiTables& tables) {
        const auto& primes = tables.primes();
        const auto& pi = tables.pi();
        std::int64_t sum = 0;

        for (std::size_t row = begin; row < end; ++row) {
            const std::int64_t first = a + 1 + static_cast<std::int64_t>(row);
            sum += countP3Row(rows[row].quotient, first, rows[row].last, primes, pi);
        }

        return sum;
    }

} // namespace

    std::int64_t P3(std::int64_t x, std::int64_t a, int nThreads,
                    const PhiTables& tables) {
        PartialSieveDetail::Validate(x, a, tables);
        if (x < 8) return 0;

        // For i <= j <= k, p[i]^3 <= x. The caller must prepare enough
        // primes to determine this outer bound, even when the sum is empty.
        const std::int64_t cubeRoot = PartialSieveDetail::IntegerCbrt(x);
        if (tables.limit() < cubeRoot) {
            throw std::out_of_range("P3 requires tables through floor(cuberoot(x))");
        }

        const std::int64_t b = tables.pi()[cubeRoot];
        if (a >= b) return 0;

        // The largest pi lookup occurs at i = j = a + 1. Divide in stages
        // to avoid forming a product. This coverage also supplies every
        // upper j bound: floor(sqrt(x / p[i])) <= floor(x / p[a+1]^2).
        // P3 borrows the caller's dense table; it never builds a larger one.
        const std::int64_t firstPrime = tables.primes()[a + 1];
        const std::int64_t maxLookup = x / firstPrime / firstPrime;
        if (tables.limit() < maxLookup) {
            throw std::out_of_range("P3 requires tables through floor(x / p[a+1]^2)");
        }

        if (nThreads <= 1 || b - a == 1) return countP3Serial(x, a, b, tables);

        // Small sums do not amortize worker launch costs. Row lengths decrease
        // with i, so the first row times the number of rows bounds all work.
        // Use division for this conservative check, avoiding both a product
        // and row-metadata allocation when the sum is certainly small.
        constexpr std::int64_t minParallelWork = 131072;
        const std::int64_t maxRowWork = tables.pi()[
            PartialSieveDetail::IntegerSqrt(x / firstPrime)
        ] - a;
        if (maxRowWork <= (minParallelWork - 1) / (b - a)) {
            return countP3Serial(x, a, b, tables);
        }

        // Rows have different lengths. Count their actual inner iterations
        // once, then partition contiguous rows by this cost rather than by
        // the number of outer primes. Every iteration contributes at least
        // one triple, so cumulativeWork <= P3(x, a) <= x fits in int64_t.
        std::vector<P3Row> rows;
        rows.reserve(static_cast<std::size_t>(b - a));
        std::int64_t totalWork = 0;

        for (std::int64_t i = a + 1; i <= b; ++i) {
            const std::int64_t quotient = x / tables.primes()[i];
            const std::int64_t last = tables.pi()[
                PartialSieveDetail::IntegerSqrt(quotient)
            ];
            totalWork += last - i + 1;
            rows.push_back({quotient, last, totalWork});
        }

        // Conservative initial grain: tiny sums stay serial. This controls
        // launch overhead; it is a tuning choice, not a mathematical limit.
        constexpr std::int64_t minWorkPerThread = 32768;
        const std::size_t workers = static_cast<std::size_t>(
            std::min<std::int64_t>(nThreads, std::min<std::int64_t>(
                b - a, totalWork / minWorkPerThread
            ))
        );

        if (totalWork < minParallelWork || workers < 2) {
            return countP3Rows(rows, 0, rows.size(), a, tables);
        }

        // Declare all referenced storage before the batch. Its destructor
        // joins existing workers if a subsequent launch or worker fails.
        std::vector<std::int64_t> sums(workers, 0);
        ThreadBatch batch(workers);
        std::size_t begin = 0;

        for (std::size_t worker = 0; worker < workers; ++worker) {
            std::size_t end = rows.size();

            if (worker + 1 < workers) {
                const std::int64_t done = begin ? rows[begin - 1].cumulativeWork : 0;
                const std::int64_t target = done + (totalWork - done) /
                    static_cast<std::int64_t>(workers - worker);
                const auto cut = std::lower_bound(
                    rows.begin() + begin, rows.end(), target,
                    [](const P3Row& row, std::int64_t work) {
                        return row.cumulativeWork < work;
                    }
                );
                end = static_cast<std::size_t>(cut - rows.begin()) + 1;
                // Retain at least one complete row for each later worker.
                end = std::min(end, rows.size() - (workers - worker - 1));
            }

            batch.launch([&, begin, end, worker] {
                sums[worker] = countP3Rows(rows, begin, end, a, tables);
            });
            begin = end;
        }

        batch.finish();
        return std::accumulate(sums.begin(), sums.end(), std::int64_t(0));
    }
}
