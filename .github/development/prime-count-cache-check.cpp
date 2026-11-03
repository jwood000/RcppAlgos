// Standalone cache checks, including fault injection; see testing.md.
#include <atomic>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>
#include <string>
#include "prime-count-core.cpp"

namespace {
    std::atomic<bool> failPageAllocation(false);
    std::atomic<std::size_t> allocationFailures(0);

    void Check(bool condition, const char* message) {
        if (!condition) throw std::runtime_error(message);
    }

    std::vector<int> PhiReference(const std::vector<int>& primes, int a,
                                  int maximum = 65535) {
        std::vector<char> survives(maximum + 1, 1);
        survives[0] = 0;
        for (int i = 0; i < a; ++i) {
            for (int x = primes[i]; x <= maximum; x += primes[i]) survives[x] = 0;
        }
        std::vector<int> counts(maximum + 1, 0);
        for (int x = 1; x <= maximum; ++x) counts[x] = counts[x-1] + survives[x];
        return counts;
    }

    int CheckPhiWorker(const std::vector<int>& primes, int a,
                       const std::vector<int>& inputs) {
        const auto phi = PhiReference(primes, a,
                                      *std::max_element(inputs.begin(), inputs.end()));
        for (int x : inputs) {
            const auto positive = PrimeCounting::phiWorker<1>(x, a);
            const auto negative = PrimeCounting::phiWorker<-1>(x, a);
            if (positive != phi[x] || negative != -phi[x]) {
                throw std::runtime_error("independent phi mismatch at x=" +
                    std::to_string(x) + ", a=" + std::to_string(a));
            }
        }
        return 2 * inputs.size();
    }
}

// Only the standalone executable replaces this allocation overload. Other
// allocations proceed normally, and the production code has no test hook.
void* operator new(std::size_t bytes, const std::nothrow_t&) noexcept {
    if (bytes == sizeof(PrimeCounting::PhiPage) &&
        failPageAllocation.load(std::memory_order_relaxed)) {
        allocationFailures.fetch_add(1, std::memory_order_relaxed);
        return nullptr;
    }
    try {
        return ::operator new(bytes);
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}

int main(int argc, char** argv) {
    try {
        const int threads = argc > 1 ? std::atoi(argv[1]) : 4;
        Check(threads > 0, "positive thread count required");

        const int limit = 2000000;
        std::vector<char> prime(limit + 1, 1);
        prime[0] = prime[1] = 0;
        for (int p = 2; p * p <= limit; ++p) {
            if (prime[p]) {
                for (int x = p * p; x <= limit; x += p) prime[x] = 0;
            }
        }
        std::vector<int> pi(limit + 1, 0), primes;
        for (int x = 1; x <= limit; ++x) {
            pi[x] = pi[x-1] + prime[x];
            if (prime[x]) primes.push_back(x);
        }

        if (argc > 2 && std::string(argv[2]) == "--allocation-failure") {
            failPageAllocation.store(true, std::memory_order_relaxed);
            Check(PrimeCounting::MainPrimeCount(10000000, threads, threads) == 664579,
                  "prime count changed when page allocations failed");
            Check(allocationFailures.load(std::memory_order_relaxed) > 0,
                  "allocation failure was not exercised");
            for (const auto& row : PrimeCounting::phiCache) {
                for (const auto& page : row.pages) {
                    Check(page.value.load(std::memory_order_relaxed) == nullptr,
                          "failed allocation published a page");
                }
            }
            failPageAllocation.store(false, std::memory_order_relaxed);

            // Recovery, unallocated blocks, zero entries, signs, and both ends
            // of the cache's range, using independently calculated phi values.
            for (int a : {7, 99}) {
                const auto phi = PhiReference(primes, a);
                PrimeCounting::updateCache(4096, a, -phi[4096]);
                Check(PrimeCounting::getCachedPhi(4096, a) == phi[4096],
                      "cache did not recover after allocation failure");
                Check(PrimeCounting::getCachedPhi(4095, a) == 0,
                      "unallocated block was not a miss");
                Check(PrimeCounting::getCachedPhi(4097, a) == 0,
                      "new cells were not initialized to zero");
                for (int x : {1, 4095, 4097, 65535}) {
                    PrimeCounting::updateCache(x, a, phi[x]);
                    Check(PrimeCounting::getCachedPhi(x, a) == phi[x],
                          "cache boundary mismatch");
                }
            }
            PrimeCounting::updateCache(65536, 7, 1);
            PrimeCounting::updateCache(1, 100, 1);
            Check(PrimeCounting::getCachedPhi(65536, 7) == 0 &&
                  PrimeCounting::getCachedPhi(1, 100) == 0,
                  "out-of-range cache key was accepted");
            Check(PrimeCounting::MainPrimeCount(10000000, threads, threads) == 664579,
                  "prime count changed after cache recovery");
            std::cout << "PASS allocation failure and recovery; rejected allocations="
                      << allocationFailures.load(std::memory_order_relaxed) << "\n";
            return 0;
        }

        // Match the tables for parent n=10^8, using only the independent sieve.
        // Disable cache publication so both signs exercise the worker itself.
        const int parentRoot = 10000;
        PrimeCounting::phiPi.assign(pi.begin(), pi.begin() + parentRoot + 1);
        PrimeCounting::phiPrimes.assign(1, 0);
        for (int p : primes) {
            if (p > parentRoot) break;
            PrimeCounting::phiPrimes.push_back(p);
        }
        failPageAllocation.store(true, std::memory_order_relaxed);
        int phiComparisons = 0;
        for (int p : {41, 43, 47}) {
            const int cube = p * p * p;
            for (int delta = -1; delta <= 1; ++delta) {
                // Stop before, at, and after the prime whose cube is crossed.
                phiComparisons += CheckPhiWorker(primes, pi[p] + delta,
                                                {cube - 1, cube, cube + 1});
            }
        }
        const int tableBoundary = 127 * static_cast<int>(PrimeCounting::phiPi.size());
        for (int delta = -1; delta <= 1; ++delta) {
            // Here the square condition already holds; table membership decides.
            phiComparisons += CheckPhiWorker(primes, pi[127] + delta,
                {tableBoundary - 1, tableBoundary, tableBoundary + 1});
        }
        phiComparisons += CheckPhiWorker(primes, 8, {4000}); // Entire loop is direct.
        // An empty recursive loop is preempted by the initial direct shortcut
        // for valid tables; exercise that degenerate input without artificial state.
        phiComparisons += CheckPhiWorker(primes, 7, {169});
        failPageAllocation.store(false, std::memory_order_relaxed);

        struct Case { std::int64_t n, expected; };
        const Case cases[] = {{10000000000LL, 455052511LL}, {100000000, 5761455},
            {1000000000000LL, 37607912018LL}, {10000000, 664579},
            {10000000000LL, 455052511LL}};
        for (const auto& c : cases) {
            Check(PrimeCounting::MainPrimeCount(c.n, threads, threads) == c.expected,
                  "known prime count mismatch");
        }

        int comparisons = 0;
        for (int n = 100000; n <= limit; n += 9973) {
            Check(PrimeCounting::MainPrimeCount(n, threads, threads) == pi[n],
                  "independent sieve mismatch");
            ++comparisons;
        }
        for (int p : primes) {
            if (p < 317) continue;
            if (p > 1400) break;
            for (int delta = -1; delta <= 1; ++delta) {
                const int n = p * p + delta;
                Check(PrimeCounting::MainPrimeCount(n, threads, threads) == pi[n],
                      "prime-square boundary mismatch");
                ++comparisons;
            }
        }

        std::size_t populated = 0;
        for (int a = 1; a < PrimeCounting::MAX_A; ++a) {
            const auto phi = PhiReference(primes, a);
            for (int x = 1; x < 65536; ++x) {
                const auto value = PrimeCounting::getCachedPhi(x, a);
                if (value) {
                    Check(value == phi[x], "cached phi mismatch");
                    ++populated;
                }
            }
        }
        PrimeCounting::PhiCell cell;
        PrimeCounting::PublishedPhiPage page;
        std::cout << "PASS threads=" << threads << " independent_pi=" << comparisons
                  << " independent_phi=" << phiComparisons
                  << " valid_cache_cells=" << populated
                  << " atomic16_bytes=" << sizeof(cell.value)
                  << " atomic16_lock_free=" << cell.value.is_lock_free()
                  << " pointer_lock_free=" << page.value.is_lock_free() << "\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
