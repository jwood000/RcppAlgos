// Independent native checks for P2/P3. Compile with -std=c++11 -pthread
// -Iinst/include and link src/P2.cpp, src/P3.cpp, src/Phi.cpp, and
// src/Eratosthenes.cpp. Use separate ASan/UBSan and TSan executables.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/Phi.h"
#include "NumbersUtils/ThreadBatch.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
    struct Reference {
        std::vector<int> primes, pi, smallest, omega;

        Reference() : pi(9000001, 0), smallest(100001, 0), omega(100001, 0) {
            const int maximum = static_cast<int>(pi.size()) - 1;
            std::vector<char> prime(maximum + 1, 1);
            prime[0] = prime[1] = 0;
            for (int p = 2; p * p <= maximum; ++p) if (prime[p]) {
                for (int n = p * p; n <= maximum; n += p) prime[n] = 0;
            }
            for (int n = 1; n <= maximum; ++n) {
                pi[n] = pi[n - 1] + prime[n];
                if (prime[n]) primes.push_back(n);
            }

            // A separate smallest-prime-factor sieve supplies the small oracle.
            for (int p = 2; p < static_cast<int>(smallest.size()); ++p) {
                if (smallest[p]) continue;
                for (int n = p; n < static_cast<int>(smallest.size()); n += p) {
                    if (!smallest[n]) smallest[n] = p;
                }
            }
            for (int n = 2; n < static_cast<int>(omega.size()); ++n) {
                omega[n] = omega[n / smallest[n]] + 1;
            }
        }

        std::int64_t factors(int x, int a, int multiplicity) const {
            const int cutoff = a ? primes[a - 1] : 1;
            std::int64_t count = 0;
            for (int n = 2; n <= x; ++n) {
                if (omega[n] == multiplicity && smallest[n] > cutoff) ++count;
            }
            return count;
        }

        std::int64_t pairs(std::int64_t x, int a) const {
            std::int64_t count = 0;
            for (std::size_t i = a; primes[i] <= x / primes[i]; ++i) {
                const auto quotient = x / primes[i];
                if (quotient >= static_cast<std::int64_t>(pi.size())) {
                    throw std::runtime_error("P2 reference table too short");
                }
                count += pi[quotient] - static_cast<std::int64_t>(i);
            }
            return count;
        }

        std::int64_t triples(std::int64_t x, int a, bool enumerate) const {
            std::int64_t count = 0;
            for (std::size_t i = a; primes[i] <= x / primes[i] / primes[i]; ++i) {
                const auto quotient = x / primes[i];
                for (std::size_t j = i; primes[j] <= quotient / primes[j]; ++j) {
                    const auto upper = quotient / primes[j];
                    if (upper >= static_cast<std::int64_t>(pi.size())) {
                        throw std::runtime_error("P3 reference table too short");
                    }
                    if (enumerate) {
                        // Explicitly visit every prime-factor multiset i<=j<=k.
                        for (std::size_t k = j; primes[k] <= upper; ++k) ++count;
                    } else {
                        count += pi[upper] - static_cast<std::int64_t>(j);
                    }
                }
            }
            return count;
        }
    };

    int IntegerRoot(std::int64_t x, int power) {
        int lower = 0, upper = static_cast<int>(std::min<std::int64_t>(x, 100000000));
        while (lower < upper) {
            const int middle = lower + (upper - lower + 1) / 2;
            std::int64_t product = 1;
            bool fits = true;
            for (int i = 0; i < power; ++i) {
                if (product > x / middle) { fits = false; break; }
                product *= middle;
            }
            if (fits) lower = middle;
            else upper = middle - 1;
        }
        return lower;
    }

    void Check(std::int64_t actual, std::int64_t expected, const char* label,
               std::int64_t x, std::int64_t a, int workers) {
        if (actual != expected) {
            throw std::runtime_error(std::string(label) + " mismatch at x=" +
                std::to_string(x) + " a=" + std::to_string(a) + " workers=" +
                std::to_string(workers) + ": actual=" + std::to_string(actual) +
                " expected=" + std::to_string(expected));
        }
    }

    template <typename Callable>
    void ExpectRange(Callable call, const char* label) {
        try { call(); }
        catch (const std::out_of_range&) { return; }
        throw std::runtime_error(std::string("missing out_of_range: ") + label);
    }

    struct ReleaseOnExit {
        std::atomic<bool>& flag;
        explicit ReleaseOnExit(std::atomic<bool>& value) : flag(value) {}
        ~ReleaseOnExit() { flag.store(true, std::memory_order_release); }
    };
}

int main() {
    try {
        using PrimeCounting::P2;
        using PrimeCounting::P3;
        using PrimeCounting::Phi;
        using PrimeCounting::PhiTables;
        const std::int64_t maxX = 9007199254740991LL;
        std::size_t rootChecks = 0;
        const auto checkRoots = [&](std::int64_t x) {
            Check(PrimeCounting::PartialSieveDetail::IntegerSqrt(x), IntegerRoot(x, 2),
                  "integer square root", x, 0, 0);
            Check(PrimeCounting::PartialSieveDetail::IntegerCbrt(x), IntegerRoot(x, 3),
                  "integer cube root", x, 0, 0);
            rootChecks += 2;
        };
        for (std::int64_t x = 0; x <= 64; ++x) checkRoots(x);
        // The binary-search oracle uses division only, independently of the
        // floating-point estimates and correction loops under test.
        const std::int64_t squareRoots[] = {3, 17, 1000000, 67108863, 94906265};
        const std::int64_t cubeRoots[] = {3, 17, 1000, 65535, 208063};
        for (std::int64_t root : squareRoots) {
            const std::int64_t square = root * root;
            for (int delta = -1; delta <= 1; ++delta) checkRoots(square + delta);
        }
        for (std::int64_t root : cubeRoots) {
            const std::int64_t cube = root * root * root;
            for (int delta = -1; delta <= 1; ++delta) checkRoots(cube + delta);
        }
        for (int delta = -2; delta <= 0; ++delta) checkRoots(maxX + delta);

        const Reference reference;
        std::map<int, std::unique_ptr<PhiTables>> prepared;
        const auto table = [&](int limit) -> const PhiTables& {
            auto found = prepared.find(limit);
            if (found == prepared.end()) {
                found = prepared.emplace(limit, std::unique_ptr<PhiTables>(new PhiTables(limit))).first;
            }
            return *found->second;
        };

        const int budgets[] = {-2, 0, 1, 2, 3, 4, 8};
        std::size_t small = 0;
        const auto& smallTable = table(512);
        for (int a = 0; a <= reference.pi[512]; ++a) {
            const int cutoff = a ? reference.primes[a - 1] : 1;
            std::int64_t count2 = 0, count3 = 0;
            for (int x = 0; x <= 512; ++x) {
                if (x >= 2 && reference.smallest[x] > cutoff) {
                    count2 += reference.omega[x] == 2;
                    count3 += reference.omega[x] == 3;
                }
                const int workers = budgets[(x + a) % 7];
                Check(P2(x, a, workers, smallTable), count2, "small P2", x, a, workers);
                Check(P3(x, a, workers, smallTable), count3, "small P3", x, a, workers);
                small += 2;
            }
        }

        std::size_t factorBoundaries = 0;
        const int products[] = {37 * 37, 37 * 37 * 37, 37 * 37 * 41, 37 * 41 * 41};
        for (int product : products) for (int delta = -1; delta <= 1; ++delta) {
            const int x = product + delta;
            for (int a = reference.pi[37] - 1; a <= reference.pi[37] + 1; ++a) {
                Check(P2(x, a, 4, table(1000)), reference.factors(x, a, 2),
                      "factor boundary P2", x, a, 4);
                Check(P3(x, a, 4, table(1000)), reference.factors(x, a, 3),
                      "factor boundary P3", x, a, 4);
                factorBoundaries += 2;
            }
        }

        // Retain endpoints around the former stream's block boundary as pi
        // query regressions, including the prime square 257^2=66049.
        std::size_t queryChecks = 0;
        for (int end : {66046, 66047, 66048, 66049, 66050}) {
            const std::int64_t x = 2 * end;
            Check(P2(x, 0, 4, table(512)), reference.pairs(x, 0),
                  "quotient boundary P2", x, 0, 4);
            ++queryChecks;
        }
        // Former worker-grain boundaries remain numerical regression inputs;
        // the current scheduler instead partitions prime-index chunks.
        for (int delta = -1; delta <= 1; ++delta) {
            const std::int64_t x = 2LL * (4096 + 2 * 1048576 + delta);
            Check(P2(x, 0, 4, table(4096)), reference.pairs(x, 0),
                  "quotient threshold P2", x, 0, 4);
            ++queryChecks;
        }
        // Keep this sparse quotient sequence near x/5, x/3, and x/2. The
        // independent oracle remains bounded at nine million.
        Check(P2(18000000, 0, 8, table(IntegerRoot(18000000, 2))),
              reference.pairs(18000000, 0), "sparse quotients P2", 18000000, 0, 8);
        ++queryChecks;

        // At p_b^2, only the last quotient fits this root-sized table.
        // Leave 4097/4096/4095 computed queries to straddle the two-worker
        // threshold and exact/partial chunks of 32 prime indices.
        const int boundaryB = 10000;
        const int boundaryPrime = reference.primes[boundaryB - 1];
        const std::int64_t boundaryX = std::int64_t(boundaryPrime) * boundaryPrime;
        for (int delta = -1; delta <= 1; ++delta) {
            const int a = boundaryB - 4097 + delta;
            Check(P2(boundaryX, a, 4, table(boundaryPrime)),
                  reference.pairs(boundaryX, a), "query chunk boundary P2",
                  boundaryX, a, 4);
            ++queryChecks;
        }

        const std::int64_t pairX = 10000000000LL;
        const int pairA = reference.pi[IntegerRoot(pairX, 3)];
        const auto pairExpected = reference.pairs(pairX, pairA);
        for (int workers : {1, 2, 3, 4, 8}) {
            Check(P2(pairX, pairA, workers, table(100000)), pairExpected,
                  "large independent P2", pairX, pairA, workers);
            ++queryChecks;
        }
        Check(P2(pairX, pairA, 4, table(1000000)), pairExpected,
              "oversized table P2", pairX, pairA, 4);
        ++queryChecks;

        const std::int64_t enumerationX = 100000000;
        const int enumerationA = reference.pi[IntegerRoot(enumerationX, 4)];
        const auto enumerated = reference.triples(enumerationX, enumerationA, true);
        Check(P3(enumerationX, enumerationA, 4, table(10000)), enumerated,
              "explicit triples P3", enumerationX, enumerationA, 4);

        const std::int64_t tripleX = 100000000000LL;
        const int tripleA = reference.pi[IntegerRoot(tripleX, 4)];
        const int firstPrime = reference.primes[tripleA];
        const int tripleLimit = static_cast<int>(tripleX / firstPrime / firstPrime);
        const auto tripleExpected = reference.triples(tripleX, tripleA, false);
        std::size_t tripleChecks = 1;
        for (int workers : {1, 2, 3, 4, 8}) {
            Check(P3(tripleX, tripleA, workers, table(tripleLimit)), tripleExpected,
                  "large independent P3", tripleX, tripleA, workers);
            ++tripleChecks;
        }
        Check(P3(tripleX, tripleA, 3, table(IntegerRoot(tripleX, 2))), tripleExpected,
              "oversized table P3", tripleX, tripleA, 3);
        ++tripleChecks;

        std::size_t coverage = 0;
        for (int x = 0; x < 4; ++x) {
            Check(P2(x, 0, 8, table(0)), 0, "tiny P2", x, 0, 8);
            ++coverage;
        }
        for (int x = 0; x < 8; ++x) {
            Check(P3(x, 0, 8, table(0)), 0, "tiny P3", x, 0, 8);
            ++coverage;
        }
        Check(P2(120, 4, 8, table(10)), 0, "empty square cutoff", 120, 4, 8);
        Check(P2(121, 4, 8, table(11)), 1, "exact square table", 121, 4, 8);
        Check(P2(121, 5, 8, table(11)), 0, "all prepared primes excluded", 121, 5, 8);
        Check(P3(342, 3, 8, table(6)), 0, "empty cube cutoff", 342, 3, 8);
        Check(P3(343, 3, 8, table(7)), 1, "exact cube table", 343, 3, 8);
        Check(P3(344, 3, 8, table(7)), 1, "cube neighbor", 344, 3, 8);
        Check(P3(343, 4, 8, table(7)), 0, "empty final prime", 343, 4, 8);
        Check(P2(1000, 0, 3, table(31)), reference.factors(1000, 0, 2),
              "minimal root table P2", 1000, 0, 3);
        Check(P3(1000, 0, 3, table(250)), reference.factors(1000, 0, 3),
              "minimal quotient table P3", 1000, 0, 3);
        coverage += 9;
        // This is a valid maximal input with an empty triple sum. Preparing
        // only its cube-root table avoids an enormous square-root allocation.
        const int maximumCubeRoot = IntegerRoot(maxX, 3);
        const int maximumCubeA = reference.pi[maximumCubeRoot];
        Check(P3(maxX, maximumCubeA, 8, table(maximumCubeRoot)), 0,
              "maximum x empty P3", maxX, maximumCubeA, 8);
        ++coverage;

        std::size_t invalid = 0;
        ExpectRange([&] { P2(121, 4, 8, table(10)); }, "P2 square coverage");
        ExpectRange([&] { P2(1000, 0, 8, table(30)); }, "P2 root coverage");
        ExpectRange([&] { P3(343, 3, 8, table(6)); }, "P3 cube coverage");
        ExpectRange([&] { P3(1000, 0, 8, table(9)); }, "P3 root coverage");
        ExpectRange([&] { P3(1000, 0, 8, table(249)); }, "P3 quotient coverage");
        ExpectRange([&] { P3(tripleX, tripleA, 8, table(tripleLimit - 1)); },
                    "P3 large quotient coverage");
        invalid += 6;
        const std::int64_t invalidX[] = {-1, maxX + 1, std::numeric_limits<std::int64_t>::max()};
        for (int limit : {0, 17, 512}) {
            const auto& context = table(limit);
            for (std::int64_t x : invalidX) {
                ExpectRange([&] { P2(x, 0, 8, context); }, "invalid P2 x");
                ExpectRange([&] { P3(x, 0, 8, context); }, "invalid P3 x");
                invalid += 2;
            }
            const std::int64_t invalidA[] = {-1, context.primeCount() + 1,
                                            std::numeric_limits<std::int64_t>::max()};
            for (std::int64_t a : invalidA) {
                ExpectRange([&] { P2(0, a, 8, context); }, "invalid P2 a at zero");
                ExpectRange([&] { P3(0, a, 8, context); }, "invalid P3 a at zero");
                invalid += 2;
            }
        }
        ExpectRange([&] { P2(maxX, 0, 8, table(512)); }, "maximum x P2 coverage");
        ExpectRange([&] { P3(maxX, 0, 8, table(512)); }, "maximum x P3 coverage");
        invalid += 2;

        {
            const auto& pairTables = table(100000);
            const auto& tripleTables = table(tripleLimit);
            std::int64_t results[2] = {0, 0};
            std::atomic<bool> start(false);
            PrimeCounting::ThreadBatch callers(2);
            ReleaseOnExit releaseOnExit(start);
            callers.launch([&] {
                while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                results[0] = P2(pairX, pairA, 2, pairTables);
            });
            callers.launch([&] {
                while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                results[1] = P3(tripleX, tripleA, 3, tripleTables);
            });
            start.store(true, std::memory_order_release);
            callers.finish();
            Check(results[0], pairExpected, "concurrent P2", pairX, pairA, 2);
            Check(results[1], tripleExpected, "concurrent P3", tripleX, tripleA, 3);
        }

        // Test-only identities; these are not new prime-counting drivers.
        struct KnownCount { std::int64_t x, expected; };
        const KnownCount counts[] = {
            {1000000, 78498}, {10000000, 664579}, {100000000, 5761455},
            {1000000000, 50847534}, {10000000000LL, 455052511},
            {100000000000LL, 4118054813LL}, {1000000000000LL, 37607912018LL}
        };
        std::size_t identities = 0;
        for (const auto& count : counts) {
            const auto& context = table(IntegerRoot(count.x, 2));
            const int cubeA = reference.pi[IntegerRoot(count.x, 3)];
            const auto meissel = Phi(count.x, cubeA, context) + cubeA - 1 -
                                 P2(count.x, cubeA, 4, context);
            Check(meissel, count.expected, "Meissel identity", count.x, cubeA, 4);
            ++identities;
            if (count.x <= 10000000000LL) {
                const int fourthA = reference.pi[IntegerRoot(count.x, 4)];
                const auto lehmer = Phi(count.x, fourthA, context) + fourthA - 1 -
                    P2(count.x, fourthA, 4, context) - P3(count.x, fourthA, 4, context);
                Check(lehmer, count.expected, "Lehmer identity", count.x, fourthA, 4);
                ++identities;
            }
        }
        std::cout << "PASS roots=" << rootChecks << " small=" << small
                  << " factor_boundaries=" << factorBoundaries
                  << " query_p2=" << queryChecks << " independent_p3=" << tripleChecks
                  << " coverage=" << coverage << " invalid=" << invalid
                  << " concurrent=2 identities=" << identities << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
