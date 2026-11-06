#include "NumbersUtils/Phi.h"
#include "NumbersUtils/ThreadBatch.h"
#include "CppConvert/Constants.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

// Parallel outer Phi summation extracted from PiLegendre.cpp. Its interval
// formulas and batch barriers retain the existing Legendre scheduling policy.
// This computes general Phi(x, a), but the chunk sizes and thresholds have
// only been tuned for Legendre's cutoff. Reassess them for other workloads.

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

    void launchPhi(ThreadBatch& batch, std::int64_t& sum,
                   std::int64_t lower, std::int64_t upper, std::int64_t x,
                   const PhiTables& tables) {
        // Bounds are captured by value because scheduling advances them while
        // workers run. Tables and each distinct result slot outlive the batch.
        batch.launch([&sum, lower, upper, x, &tables] {
            phiForeman(sum, lower, upper, x, tables);
        });
    }

    double getChunkFactor(std::int64_t x) {
        static constexpr std::array<double, 9> nums = {
            {1e10, 1e12, 2e13, 5e13, 8e13, 1e14, 5e14, 1e15, 1e16}
        };

        static constexpr std::array<double, 9> factor = {
            {1.3, 1.2, 1.1, 1.07, 1.05, 1.01, 1.007, 1.006, 1.005}
        };

        const auto it = std::upper_bound(
            nums.cbegin(), nums.cend(), static_cast<double>(x)
        );

        return std::log(factor[it - nums.cbegin()]);
    }

    std::int64_t phiMain(std::int64_t x, std::int64_t a, int nThreads,
                         bool Parallel, const PhiTables& tables) {

        // Repeatedly apply Phi(x, i + 1) = Phi(x, i) - Phi(x / p[i + 1], i).
        // The caller has handled a <= 6; the remaining sum covers [6, a).
        // a is independent of the prepared table's extent.
        const std::int64_t strt = 6;
        std::int64_t mySum = Phi(x, strt, tables);

        if (Parallel) {
            constexpr double divLim = 10000000000.0;    //  1e10

            // The summation covers [strt, a). More workers than terms
            // can make the small-input chunk zero and move lower backwards.
            nThreads = static_cast<int>(std::min<std::int64_t>(
                nThreads, a - strt
            ));

            // We know x, t, and we need to solve for n.
            // x / (1.5^t)^n < 1e10 -->>> 1e10 * (1.5^t)^n = x --->>> x / 1e10 = (1.5^t)^n
            // log(x / 1e10) = n * log(1.5^t) --->> log(x / 1e10) / log(1.5^t) :
            //                 log(x / 1e10) / (t * log(1.5))
            std::int64_t myRange = (a - strt) + 1;
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
                double power = std::ceil(
                    std::log((firstStep + 1) / (multOne - 1)) / std::log(multOne)
                ) + 1;

                // Same as above: lower + (upper - lower) / s   and given:   m^p + lower = upper
                //   -->>  (lower + m^p / s)  -->>  firstLoops = ceil(((m^p + lower) / s) / nThreads)
                const int firstLoops = std::ceil(
                    ((std::pow(multOne, power) + base) / firstStep) / nThreads
                );

                upper = lower + firstStep;
                std::vector<std::int64_t> mySums;

                for (int i = 0; i < firstLoops && lower < a; ++i) {
                    std::vector<std::int64_t> intermediate(nThreads, 0);
                    ThreadBatch batch(nThreads);

                    for (int j = 0; j < nThreads && lower < a;
                         lower = upper, upper += firstStep, ++j) {
                        upper = std::min(upper, a);
                        launchPhi(batch, intermediate[j], lower, upper, x, tables);
                    }

                    batch.finish();
                    mySums.insert(mySums.end(), intermediate.cbegin(),
                                  intermediate.cend());
                }

                // A large worker request can finish all terms in these first
                // batches. Every launched worker has already been joined.
                if (lower == a) {
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

                while (lower < a &&
                       (dblChunk * std::pow(multOne, nThreads - 1) + base) < a) {
                    std::vector<std::int64_t> intermediate(nThreads, 0);
                    ThreadBatch batch(nThreads);

                    for (int j = 0; j < nThreads && lower < a; lower = upper, ++j,
                            dblChunk *= multOne, upper = static_cast<std::int64_t>(dblChunk) + base) {
                        // Repeated multiplication can round differently from
                        // the pow estimate used to admit the whole batch.
                        upper = std::min(upper, a);
                        launchPhi(batch, intermediate[j], lower, upper, x, tables);
                    }

                    batch.finish();

                    mySums.insert(mySums.end(), intermediate.cbegin(),
                                  intermediate.cend());
                }

                std::vector<std::int64_t> intermediate(nThreads, 0);
                ThreadBatch batch(nThreads);

                for (int j = 0; j < (nThreads - 1) && upper < a; lower = upper, ++j,
                        dblChunk *= multOne, upper = static_cast<std::int64_t>(dblChunk) + base) {
                    launchPhi(batch, intermediate[j], lower, upper, x, tables);
                }

                launchPhi(batch, intermediate.back(), lower, a, x, tables);

                batch.finish();

                mySums.insert(mySums.end(), intermediate.cbegin(),
                              intermediate.cend());

                mySum += std::accumulate(mySums.cbegin(), mySums.cend(),
                                         static_cast<std::int64_t>(0));
            } else {
                std::vector<std::int64_t> mySums(nThreads, 0);
                ThreadBatch batch(nThreads);
                std::int64_t chunk = myRange / nThreads;
                upper = lower + chunk - 1;

                for (int j = 0; j < (nThreads - 1); lower = upper, upper += chunk, ++j) {
                    launchPhi(batch, mySums[j], lower, upper, x, tables);
                }

                launchPhi(batch, mySums.back(), lower, a, x, tables);

                batch.finish();

                mySum += std::accumulate(mySums.cbegin(), mySums.cend(),
                                         static_cast<std::int64_t>(0));
            }
        } else {
            phiForeman(mySum, strt, a, x, tables);
        }

        return mySum;
    }

} // namespace

    std::int64_t PhiParallel(std::int64_t x, std::int64_t a, int nThreads,
                            const PhiTables& tables) {
        // Delegate invalid inputs to Phi so validation and error messages stay
        // identical, including invalid a when x == 0. Do this before indexing
        // primes or entering the scheduler's logarithms and worker arithmetic.
        if (x < 0 || x > static_cast<std::int64_t>(Significand53) ||
            a < 0 || a > tables.primeCount()) {
            return Phi(x, a, tables);
        }

        // No outer sum is needed for tiny a. Within the prepared table, serial
        // Phi can use its boundary/direct-pi shortcuts; in particular this
        // covers x == 0 and x <= p[a], without scheduling an already-known sum.
        if (a <= 6 || x <= tables.limit()) return Phi(x, a, tables);

        const bool Parallel = nThreads > 1 && x >= 10000000;
        return phiMain(x, a, nThreads, Parallel, tables);
    }
}
