#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/Phi.h"
#include "CppConvert/Constants.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <numeric>
#include <thread>
#include <vector>

// "PiLegendre" based off of the highly optimized
// "pi_legendre.cpp" algorithm by Kim Walisch, which calculates
// the numbers of primes less than n using Legendre's formula.
// Kim Walisch's official github repo for "pi_legendre" is:
//                      https://github.com/kimwalisch/primecount

namespace PrimeCounting {
namespace {

    void phiForeman(std::int64_t &mySum, std::int64_t lowerBound,
                    std::int64_t upperBound, std::int64_t x,
                    const PhiTables& tables) {

        const auto& primes = tables.primes();
        for (std::int64_t i = lowerBound; i < upperBound; ++i) {
            mySum -= Phi(x / primes[i + 1], i, tables);
        }
    }

    // One batch owns its workers and their exception slots. Result storage
    // must outlive the batch, including when a later thread fails to start.
    class PhiThreadBatch {
        const PhiTables& tables;
        std::vector<std::exception_ptr> errors;
        std::vector<std::thread> threads;

        void join() {
            for (auto& thread : threads) {
                if (thread.joinable()) thread.join();
            }
        }

    public:
        PhiThreadBatch(int capacity, const PhiTables& preparedTables)
            : tables(preparedTables), errors(capacity) {
            threads.reserve(capacity);
        }

        PhiThreadBatch(const PhiThreadBatch&) = delete;
        PhiThreadBatch& operator=(const PhiThreadBatch&) = delete;

        ~PhiThreadBatch() {
            // Also runs if a launch fails after earlier workers have started.
            join();
        }

        void launch(std::int64_t& sum, std::int64_t lower,
                    std::int64_t upper, std::int64_t x) {
            const std::size_t slot = threads.size();
            threads.emplace_back([this, &sum, lower, upper, x, slot] {
                try {
                    phiForeman(sum, lower, upper, x, tables);
                } catch (...) {
                    // No R API calls here. Each worker writes only its slot.
                    errors[slot] = std::current_exception();
                }
            });
        }

        void finish() {
            join();
            // Joining makes every worker's writes visible. Propagate on the
            // calling thread before any partial sums are used.
            for (const auto& error : errors) {
                if (error) std::rethrow_exception(error);
            }
        }
    };

    const double getChunkFactor(std::int64_t x) {
        static constexpr std::array<double, 9> nums = {{1e10, 1e12, 2e13, 5e13, 8e13, 1e14, 5e14, 1e15, 1e16}};
        static constexpr std::array<double, 9> factor = {{1.3, 1.2, 1.1, 1.07, 1.05, 1.01, 1.007, 1.006, 1.005}};

        const auto it = std::upper_bound(nums.cbegin(), nums.cend(), static_cast<double>(x));
        return std::log(factor[it - nums.cbegin()]);
    }

    std::int64_t phiMain(std::int64_t x, std::int64_t a, int nThreads,
                         bool Parallel, const PhiTables& tables) {

        // Legendre supplies a as its prepared root table's prime count.
        // For x >= 100000, all six tiny-table primes are below sqrt(x),
        // so strt is always 6.
        const std::int64_t piSqrtx = a;
        const std::int64_t strt = 6;
        std::int64_t mySum = Phi(x, strt, tables) + piSqrtx - a;

        if (Parallel) {
            constexpr double divLim = 10000000000.0;    //  1e10

            // The summation covers [strt, piSqrtx). More workers than terms
            // can make the small-input chunk zero and move lower backwards.
            nThreads = static_cast<int>(std::min<std::int64_t>(
                nThreads, piSqrtx - strt
            ));

            // We know x, t, and we need to solve for n.
            // x / (1.5^t)^n < 1e10 -->>> 1e10 * (1.5^t)^n = x --->>> x / 1e10 = (1.5^t)^n
            // log(x / 1e10) = n * log(1.5^t) --->> log(x / 1e10) / log(1.5^t) :
            //                 log(x / 1e10) / (t * log(1.5))
            std::int64_t myRange = (piSqrtx - strt) + 1;
            std::int64_t lower = strt;
            std::int64_t upper;

            if (x > divLim) {
                const double chunkFactor = getChunkFactor(x);
                const int nLoops = 1 + std::ceil(std::log(x / divLim) / (nThreads * chunkFactor));
                const double multOne = std::exp(std::log(myRange) / (nThreads * nLoops));
                const int firstStep = std::ceil(std::log10(Significand53 / x)) + 1;
                std::int64_t base = lower;

                //*********************** Begin of firstThreads *****************************
                // Typically, the size of multOne will be very small resulting in not very
                // meaningful chunk sizes in the first few loops. Because of this, we
                // gurantee that each of the first threads will have at least firstStep
                // iterations and that the threads following will have at least firstStep + 1

                // Here we have m = multOne, p = power, and s = firstStep:
                // m^(p + 1) - m^p = m^p * (m - 1) > (s + 1)  -->>  m^p > (s + 1) / (m - 1)  -->>
                // p * log(m) > log((s + 1) / (m - 1))  -->>  p > log((s + 1) / (m - 1)) / log(m)
                double power = std::ceil(std::log((firstStep + 1) / (multOne - 1)) / std::log(multOne)) + 1;

                // Same as above: lower + (upper - lower) / s   and given:   m^p + lower = upper
                //   -->>  (lower + m^p / s)  -->>  firstLoops = ceil(((m^p + lower) / s) / nThreads)
                const int firstLoops = std::ceil(((std::pow(multOne, power) +
                                                 base) / firstStep) / nThreads);
                upper = lower + firstStep;
                std::vector<std::int64_t> mySums;

                for (int i = 0; i < firstLoops && lower < piSqrtx; ++i) {
                    std::vector<std::int64_t> intermediate(nThreads, 0);
                    PhiThreadBatch batch(nThreads, tables);

                    for (int j = 0; j < nThreads && lower < piSqrtx;
                         lower = upper, upper += firstStep, ++j) {
                        upper = std::min(upper, piSqrtx);
                        batch.launch(intermediate[j], lower, upper, x);
                    }

                    batch.finish();

                    mySums.insert(mySums.end(), intermediate.cbegin(),
                                  intermediate.cend());
                }

                // A large worker request can finish all terms in these first
                // batches. Every launched worker has already been joined.
                if (lower == piSqrtx) {
                    return mySum + std::accumulate(mySums.cbegin(), mySums.cend(),
                                                  static_cast<std::int64_t>(0));
                }

                //*********************** End of firstThreads *****************************
                //*************************************************************************
                //*********************** Begin of midThreads *****************************

                while (std::pow(multOne, power) < (upper - base))
                    ++power;

                double dblChunk = std::pow(multOne, power);
                upper = static_cast<std::int64_t>(dblChunk) + base;

                while (lower < piSqrtx &&
                       (dblChunk * std::pow(multOne, nThreads - 1) + base) < piSqrtx) {
                    std::vector<std::int64_t> intermediate(nThreads, 0);
                    PhiThreadBatch batch(nThreads, tables);

                    for (int j = 0; j < nThreads && lower < piSqrtx; lower = upper, ++j,
                            dblChunk *= multOne, upper = static_cast<std::int64_t>(dblChunk) + base) {
                        // Repeated multiplication can round differently from
                        // the pow estimate used to admit the whole batch.
                        upper = std::min(upper, piSqrtx);
                        batch.launch(intermediate[j], lower, upper, x);
                    }

                    batch.finish();

                    mySums.insert(mySums.end(), intermediate.cbegin(),
                                  intermediate.cend());
                }

                std::vector<std::int64_t> intermediate(nThreads, 0);
                PhiThreadBatch batch(nThreads, tables);

                for (int j = 0; j < (nThreads - 1) && upper < piSqrtx; lower = upper, ++j,
                        dblChunk *= multOne, upper = static_cast<std::int64_t>(dblChunk) + base) {
                    batch.launch(intermediate[j], lower, upper, x);
                }

                batch.launch(intermediate.back(), lower, piSqrtx, x);

                batch.finish();

                mySums.insert(mySums.end(), intermediate.cbegin(),
                              intermediate.cend());

                mySum += std::accumulate(mySums.cbegin(), mySums.cend(),
                                         static_cast<std::int64_t>(0));
            } else {
                std::vector<std::int64_t> mySums(nThreads, 0);
                PhiThreadBatch batch(nThreads, tables);
                std::int64_t chunk = myRange / nThreads;
                upper = lower + chunk - 1;

                for (int j = 0; j < (nThreads - 1); lower = upper, upper += chunk, ++j) {
                    batch.launch(mySums[j], lower, upper, x);
                }

                batch.launch(mySums.back(), lower, piSqrtx, x);

                batch.finish();

                mySum += std::accumulate(mySums.cbegin(), mySums.cend(),
                                         static_cast<std::int64_t>(0));
            }
        } else {
            phiForeman(mySum, strt, piSqrtx, x, tables);
        }

        return mySum;
    }

} // namespace

    // Historical timings before shared-Phi integration (nThreads = 8).
    // Counts verified by Kim Walisch's primecount library.
    //  10^9 -->>          50,847,534   -->>   167.5 microseconds
    // 10^10 -->>         455,052,511   -->>   677.6 microseconds
    // 10^11 -->>       4,118,054,813   -->>   2.623 milliseconds
    // 10^12 -->>      37,607,912,018   -->>   14.40 milliseconds
    // 10^13 -->>     346,065,536,839   -->>   103.3 milliseconds
    // 10^14 -->>   3,204,941,750,802   -->>   0.689 seconds
    // 10^15 -->>  29,844,570,422,669   -->>   5.096 seconds
    // MAX VALUE (2^53 - 1) -->>
    //            252,252,704,148,404   -->>  37.234 seconds

    std::int64_t PiLegendre(std::int64_t n, int nThreads, int maxThreads) {

        const std::int64_t sqrtBound = std::sqrt(static_cast<double>(n));
        const PhiTables tables(sqrtBound);

        bool Parallel = false;

        if (nThreads > 1 && maxThreads > 1) {
            Parallel = true;
            if (nThreads > maxThreads) nThreads = maxThreads;
            if (n < 1e7) Parallel = false;
        }

        const std::int64_t piSqrt = tables.primeCount();
        const std::int64_t phiSqrt = phiMain(n, piSqrt, nThreads, Parallel, tables);
        const std::int64_t int64result = piSqrt + phiSqrt - 1;

        return int64result;
    }
}
