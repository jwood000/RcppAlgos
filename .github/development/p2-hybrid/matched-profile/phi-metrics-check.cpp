// Native checker for the temporary Phi counters. The expected partial-sieve
// counts come from a direct survivor sieve, independent of the Phi recurrence.
#include "NumbersUtils/Phi.h"
#include "phi_metrics.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {
    void require(bool condition) {
        if (!condition) throw std::runtime_error("Phi metrics check failed");
    }
}

int main() {
    using PrimeCounting::Phi;
    using PrimeCounting::PhiTables;
    constexpr int LIMIT = 131071;
    const PhiTables tables(1000);
    std::vector<int> primes;
    for (int n = 2; n <= 1000; ++n) {
        bool prime = true;
        for (int p : primes) {
            if (p > n / p) break;
            if (n % p == 0) { prime = false; break; }
        }
        if (prime) primes.push_back(n);
    }
    require(static_cast<std::size_t>(tables.primeCount()) == primes.size());
    std::vector<unsigned char> survives(LIMIT + 1, 1);
    survives[0] = 0;
    std::vector<std::vector<int>> expected(primes.size() + 1);
    for (std::size_t a = 0; a <= primes.size(); ++a) {
        if (a) for (int n = primes[a - 1]; n <= LIMIT; n += primes[a - 1]) survives[n] = 0;
        expected[a].resize(LIMIT + 1);
        for (int n = 1; n <= LIMIT; ++n) expected[a][n] = expected[a][n - 1] + survives[n];
    }
    std::uint64_t checks = 0;
    MatchedPhi::Stats stats;
    {
        MatchedPhi::Scope scope(stats);
        for (int a = 0; a <= 30; ++a) {
            for (int n = 0; n <= 3000; ++n) {
                require(Phi(n, a, tables) == expected[a][n]);
                require(stats.consistent());
                ++checks;
            }
        }
        for (int a : {6, 7, 25, 80, 99, 100, 168}) {
            for (int n : {4095, 4096, 4097, 65535, 65536, 65537, 131070, 131071}) {
                require(Phi(n, a, tables) == expected[a][n]);
                require(stats.consistent());
                ++checks;
            }
        }
        MatchedPhi::Stats inner;
        {
            MatchedPhi::Scope nested(inner);
            require(Phi(LIMIT, 7, tables) == expected[7][LIMIT]);
            require(inner.consistent());
            ++checks;
        }
        require(MatchedPhi::active == &stats && MatchedPhi::depth == 0);
    }
    require(MatchedPhi::active == nullptr && MatchedPhi::depth == 0);
    const auto before = stats.counts;
    require(Phi(LIMIT, 7, tables) == expected[7][LIMIT]);
    require(stats.counts == before);
    ++checks;
    // Each worker owns its Stats. The parent only reads the slots after join.
    MatchedPhi::Stats workers[2];
    bool okay[2] = {true, true};
    std::thread threads[2];
    for (int slot = 0; slot < 2; ++slot) {
        threads[slot] = std::thread([&, slot] {
            MatchedPhi::Scope scope(workers[slot]);
            for (int repeat = 0; repeat < 3; ++repeat) {
                for (int a = 7; a <= 100; ++a) {
                    const int n = LIMIT - (a * 13 + slot * 17 + repeat);
                    if (Phi(n, a, tables) != expected[a][n]) okay[slot] = false;
                }
            }
        });
    }
    for (auto& thread : threads) thread.join();
    for (int slot = 0; slot < 2; ++slot) {
        require(okay[slot] && workers[slot].consistent());
        stats.add(workers[slot]);
        checks += 3 * 94;
    }
    require(stats.consistent());
    require(stats.counts[MatchedPhi::CacheHits] > 0);
    require(stats.counts[MatchedPhi::CacheMisses] > 0);
    require(stats.counts[MatchedPhi::PageAllocSuccess] > 0);
    require(stats.counts[MatchedPhi::CacheIneligible] > 0);
    std::cout << "{\"checks\":" << checks << ",\"counts\":";
    stats.writeJson(std::cout);
    std::cout << "}\n";
}
