// Candidate-only cache-range checks. Define PHI_SOURCE as the temporary
// doubled-cache Phi.cpp path; link Eratosthenes.cpp, not a second Phi.cpp.
// Existing general-Phi and prime-count checkers supply broader regressions.
#include "NumbersUtils/Phi.h"
#include "NumbersUtils/ThreadBatch.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <new>
#include <stdexcept>
#include <thread>
#include <vector>
#ifndef PHI_SOURCE
#error Define PHI_SOURCE to the temporary doubled-cache Phi.cpp
#endif
#include PHI_SOURCE

namespace {
    std::atomic<bool> rejectPages(false);
    std::atomic<std::size_t> rejectedPages(0), allocatedPages(0);

    void Check(bool condition, const char* message) {
        if (!condition) throw std::runtime_error(message);
    }

    std::vector<int> PartialCounts(const std::vector<int>& primes, int a, int maximum) {
        std::vector<char> survives(maximum + 1, 1);
        survives[0] = 0;
        for (int i = 0; i < a; ++i) {
            for (int n = primes[i]; n <= maximum; n += primes[i]) survives[n] = 0;
        }
        std::vector<int> counts(maximum + 1, 0);
        for (int n = 1; n <= maximum; ++n) counts[n] = counts[n - 1] + survives[n];
        return counts;
    }

    struct ReleaseOnExit {
        std::atomic<bool>& flag;
        explicit ReleaseOnExit(std::atomic<bool>& value) : flag(value) {}
        ~ReleaseOnExit() { flag.store(true, std::memory_order_release); }
    };
}

void* operator new(std::size_t bytes, const std::nothrow_t&) noexcept {
    if (bytes == sizeof(PrimeCounting::PhiPage) && rejectPages.load(std::memory_order_relaxed)) {
        rejectedPages.fetch_add(1, std::memory_order_relaxed);
        return nullptr;
    }
    try {
        void* pointer = ::operator new(bytes);
        if (bytes == sizeof(PrimeCounting::PhiPage)) allocatedPages.fetch_add(1, std::memory_order_relaxed);
        return pointer;
    } catch (const std::bad_alloc&) { return nullptr; }
}

void operator delete(void* pointer, const std::nothrow_t&) noexcept { ::operator delete(pointer); }

int main() {
    try {
        using namespace PrimeCounting;
        Check(CACHE_ROWS == 100 && PAGE_SIZE == 4096 && PAGE_COUNT == 32,
              "expected doubled cache geometry");
        Check(sizeof(PhiCell::value) == sizeof(std::atomic<std::uint16_t>),
              "cache cell is no longer atomic uint16");

        std::vector<char> prime(10001, 1);
        prime[0] = prime[1] = 0;
        for (int p = 2; p * p <= 10000; ++p) if (prime[p]) {
            for (int n = p * p; n <= 10000; n += p) prime[n] = 0;
        }
        std::vector<int> primes;
        std::vector<std::int32_t> pi(10001, 0);
        for (int n = 1; n <= 10000; ++n) {
            pi[n] = pi[n - 1] + prime[n];
            if (prime[n]) primes.push_back(n);
        }
        std::map<int, std::vector<int>> reference;
        for (int a : {0, 1, 6, 7, 99, 100}) reference.emplace(a, PartialCounts(primes, a, 131073));
        Check(reference.at(6)[131071] == 25142 &&
              reference.at(6)[131071] < std::numeric_limits<std::uint16_t>::max(),
              "independent six-prime upper bound failed");
        Check(reference.at(1)[131071] == 65536, "large public base-case fixture failed");
        for (int a : {7, 99}) {
            for (int x = 1; x <= 131071; ++x) {
                Check(reference.at(a)[x] > 0 && reference.at(a)[x] <= reference.at(6)[131071],
                      "positive cached values exceed uint16 bound");
            }
        }

        const std::vector<int> storedInputs = {
            65535, 65536, 65537, 69631, 69632, 126975, 126976, 131070, 131071
        };
        {
            PhiCache cache;
            rejectPages.store(true, std::memory_order_relaxed);
            cache.put(65536, 7, reference.at(7)[65536]);
            cache.put(131071, 99, reference.at(99)[131071]);
            Check(cache.get(65536, 7) == 0 && cache.get(131071, 99) == 0 &&
                  rejectedPages.load() == 2 && allocatedPages.load() == 0,
                  "failed new-half allocations published a value");
            rejectPages.store(false, std::memory_order_relaxed);
            cache.put(65536, 7, -reference.at(7)[65536]);
            Check(cache.get(65536, 7) == reference.at(7)[65536] &&
                  cache.get(65537, 7) == 0 && cache.get(65535, 7) == 0,
                  "new-page recovery, zero initialization or neighboring miss failed");

            std::atomic<bool> start(false);
            ThreadBatch workers(4);
            ReleaseOnExit release(start);
            for (int worker = 0; worker < 4; ++worker) workers.launch([&, worker] {
                while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                for (int repeat = 0; repeat < 4; ++repeat) for (int a : {7, 99}) {
                    for (int x : storedInputs) {
                        const auto expected = reference.at(a)[x];
                        cache.put(x, a, worker % 2 ? expected : -expected);
                        Check(cache.get(x, a) == expected, "concurrent new-page publication mismatch");
                    }
                }
            });
            start.store(true, std::memory_order_release);
            workers.finish();
            Check(allocatedPages.load() == 10, "duplicate writers allocated unexpected page count");
            cache.put(131072, 7, reference.at(7)[131072]);
            cache.put(131071, 100, reference.at(100)[131071]);
            Check(cache.get(131072, 7) == 0 && cache.get(131071, 100) == 0 &&
                  allocatedPages.load() == 10, "out-of-range key allocated or returned a cache cell");
        }

        const PhiTables small(541), large(10000);
        Check(large.pi() == pi, "prepared pi disagrees with independent sieve");
        std::vector<std::int64_t> expectedPrimes(1, 0);
        expectedPrimes.insert(expectedPrimes.end(), primes.begin(), primes.end());
        Check(large.primes() == expectedPrimes, "prepared primes disagree with independent sieve");
        std::vector<int> inputs = storedInputs;
        inputs.push_back(131072); inputs.push_back(131073);

        // Start public evaluation concurrently while the shared cache's newly
        // eligible pages are still cold. The earlier cache was a local object.
        std::atomic<bool> start(false);
        ThreadBatch workers(4);
        ReleaseOnExit release(start);
        for (int worker = 0; worker < 4; ++worker) workers.launch([&, worker] {
            while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
            const auto& tables = worker % 2 ? small : large;
            for (int repeat = 0; repeat < 2; ++repeat) for (int a : {0, 1, 6, 7, 99, 100}) {
                for (int x : inputs) {
                    Check(Phi(x, a, tables) == reference.at(a)[x], "concurrent public Phi mismatch");
                }
            }
        });
        start.store(true, std::memory_order_release);
        workers.finish();

        PhiEvaluator evaluator(large.primes(), pi);
        std::size_t signedChecks = 0;
        for (int a : {6, 7, 99, 100}) for (int x : inputs) {
            // x exceeds both the prepared table and p[a], so these are valid
            // private entries after the public direct shortcut has failed.
            Check(evaluator.evaluate<1>(x, a) == reference.at(a)[x] &&
                  evaluator.evaluate<-1>(x, a) == -reference.at(a)[x], "signed boundary Phi mismatch");
            signedChecks += 2;
            if (a == 6 || a == 100 || x > 131071) {
                Check(sharedCache().get(x, a) == 0, "bypassed/out-of-range shared key was populated");
            } else {
                Check(sharedCache().get(x, a) == reference.at(a)[x], "new-half shared cell mismatch");
            }
        }
        std::cout << "PASS boundary_inputs=" << inputs.size() << " positive_bounds=262142"
                  << " duplicate_publications=288 public_concurrent=528 signed=" << signedChecks
                  << " rejected_allocations=" << rejectedPages.load() << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
