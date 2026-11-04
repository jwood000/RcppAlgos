// Standalone checks for the general Phi prototype; no R runtime is required.
/* Run from the repository root:
phi_check=$(mktemp -d)
clang++ -std=c++11 -O1 -g -pthread -fsanitize=address,undefined \
    -fno-omit-frame-pointer -Iinst/include .github/development/phi-check.cpp \
    src/Phi.cpp src/Eratosthenes.cpp -o "$phi_check/phi-asan-ubsan"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$phi_check/phi-asan-ubsan"
ASAN_OPTIONS=halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 "$phi_check/phi-asan-ubsan" --no-cache
clang++ -std=c++11 -O1 -g -pthread -fsanitize=thread \
    -fno-omit-frame-pointer -Iinst/include .github/development/phi-check.cpp \
    src/Phi.cpp src/Eratosthenes.cpp -o "$phi_check/phi-tsan"
TSAN_OPTIONS=halt_on_error=1 "$phi_check/phi-tsan"
TSAN_OPTIONS=halt_on_error=1 "$phi_check/phi-tsan" --no-cache
*/
#include "NumbersUtils/Phi.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {
    std::atomic<bool> rejectCacheAllocations(false);
    std::atomic<std::size_t> rejectedAllocations(0);

    struct IndependentSieve {
        std::vector<int> primes;
        std::vector<int> pi;

        explicit IndependentSieve(int maximum) : pi(maximum + 1, 0) {
            std::vector<char> prime(maximum + 1, 1);
            prime[0] = prime[1] = 0;
            for (int p = 2; p * p <= maximum; ++p) {
                if (prime[p]) {
                    for (int x = p * p; x <= maximum; x += p) prime[x] = 0;
                }
            }
            for (int x = 1; x <= maximum; ++x) {
                pi[x] = pi[x - 1] + prime[x];
                if (prime[x]) primes.push_back(x);
            }
        }
    };

    std::vector<int> PartialCounts(const std::vector<int>& primes, int a,
                                   int maximum) {
        std::vector<char> survives(maximum + 1, 1);
        survives[0] = 0;
        for (int i = 0; i < a; ++i) {
            for (int x = primes[i]; x <= maximum; x += primes[i]) survives[x] = 0;
        }
        std::vector<int> counts(maximum + 1, 0);
        for (int x = 1; x <= maximum; ++x) counts[x] = counts[x - 1] + survives[x];
        return counts;
    }

    void CheckPhi(std::int64_t x, std::int64_t a,
                  const PrimeCounting::PhiTables& tables,
                  std::int64_t expected, const char* category) {
        const auto actual = PrimeCounting::Phi(x, a, tables);
        if (actual != expected) {
            throw std::runtime_error(std::string(category) + ": Phi(" +
                std::to_string(x) + ", " + std::to_string(a) + ") with limit=" +
                std::to_string(tables.limit()) + " returned " + std::to_string(actual) +
                ", expected " + std::to_string(expected));
        }
    }

    std::size_t CheckInputs(const IndependentSieve& sieve,
                           const PrimeCounting::PhiTables& tables, int a,
                           const std::vector<int>& inputs) {
        const auto reference = PartialCounts(sieve.primes, a,
            *std::max_element(inputs.begin(), inputs.end()));
        for (int x : inputs) CheckPhi(x, a, tables, reference[x], "targeted sieve");
        return inputs.size();
    }

    template <typename Callable>
    void ExpectOutOfRange(Callable call, const char* label) {
        try {
            call();
        } catch (const std::out_of_range&) {
            return;
        }
        throw std::runtime_error(std::string("missing out_of_range: ") + label);
    }

    std::int64_t InclusionExclusion(std::int64_t x, int a,
                                    const std::vector<int>& primes) {
        // Used only for a <= 8: all products and the signed running sum fit int64.
        std::int64_t result = x;
        for (unsigned mask = 1; mask < (1u << a); ++mask) {
            std::int64_t product = 1;
            int selected = 0;
            for (int i = 0; i < a; ++i) {
                if (mask & (1u << i)) {
                    product *= primes[i];
                    ++selected;
                }
            }
            result += (selected % 2 ? -1 : 1) * (x / product);
        }
        return result;
    }

    using Contexts = std::map<int, std::unique_ptr<PrimeCounting::PhiTables>>;

    std::size_t CheckConcurrent(const IndependentSieve& sieve, const Contexts& contexts) {
        struct Configuration { int limit, a; };
        const Configuration configurations[] = {
            {17, 7}, {10000, 7}, {19, 8}, {10000, 8},
            {97, 24}, {100, 24}, {10000, 24}, {97, 25},
            {100, 25}, {10000, 25}, {10000, 31}, {10000, 64}, {10000, 99}
        };
        struct Job {
            const PrimeCounting::PhiTables* tables;
            int x, a, expected;
        };
        std::vector<Job> jobs;
        for (const auto& configuration : configurations) {
            const auto reference = PartialCounts(sieve.primes, configuration.a, 65535);
            for (int i = 0; i < 32; ++i) {
                const int x = 30000 + (i * 977) % 35000;
                jobs.push_back({contexts.at(configuration.limit).get(), x,
                                configuration.a, reference[x]});
            }
        }

        // This runs before other evaluations, keeping these cache pages cold.
        // Workers share contexts and duplicate keys across different table sizes.
        std::atomic<bool> start(false);
        std::exception_ptr failures[4];
        std::vector<std::thread> workers;
        try {
            for (std::size_t worker = 0; worker < 4; ++worker) {
                workers.emplace_back([&, worker] {
                    while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                    try {
                        for (int repeat = 0; repeat < 2; ++repeat) {
                            for (std::size_t i = 0; i < jobs.size(); ++i) {
                                const auto& job = jobs[(i + worker * 17) % jobs.size()];
                                CheckPhi(job.x, job.a, *job.tables, job.expected, "concurrent");
                            }
                        }
                    } catch (...) {
                        failures[worker] = std::current_exception();
                    }
                });
            }
        } catch (...) {
            start.store(true, std::memory_order_release);
            for (auto& worker : workers) worker.join();
            throw;
        }
        start.store(true, std::memory_order_release);
        for (auto& worker : workers) worker.join();
        for (const auto& failure : failures) {
            if (failure) std::rethrow_exception(failure);
        }
        return jobs.size() * 8;
    }
}

// Only this standalone executable replaces nothrow allocation. Ordinary
// allocations for tables and independent references still use the normal heap.
void* operator new(std::size_t bytes, const std::nothrow_t&) noexcept {
    if (rejectCacheAllocations.load(std::memory_order_relaxed)) {
        rejectedAllocations.fetch_add(1, std::memory_order_relaxed);
        return nullptr;
    }
    try {
        return ::operator new(bytes);
    } catch (const std::bad_alloc&) {
        return nullptr;
    }
}

void operator delete(void* address, const std::nothrow_t&) noexcept {
    ::operator delete(address);
}

int main(int argc, char** argv) {
    try {
        if (argc > 2 || (argc == 2 && std::string(argv[1]) != "--no-cache")) return 2;
        const bool noCache = argc == 2;
        rejectCacheAllocations.store(noCache, std::memory_order_relaxed);
        const IndependentSieve sieve(1300000);
        Contexts contexts;
        for (int limit : {0, 1, 2, 13, 17, 19, 97, 100, 10000}) {
            auto tables = std::unique_ptr<PrimeCounting::PhiTables>(new PrimeCounting::PhiTables(limit));
            if (tables->limit() != limit || tables->primeCount() != sieve.pi[limit]) {
                throw std::runtime_error("incorrect PhiTables metadata at limit=" + std::to_string(limit));
            }
            contexts.emplace(limit, std::move(tables));
        }

        const auto concurrent = CheckConcurrent(sieve, contexts);
        std::size_t exhaustive = 0;
        for (const auto& context : contexts) {
            const auto& tables = *context.second;
            std::vector<char> survives(2001, 1);
            survives[0] = 0;
            for (int a = 0; a <= tables.primeCount(); ++a) {
                if (a) {
                    const int p = sieve.primes[a - 1];
                    for (int x = p; x <= 2000; x += p) survives[x] = 0;
                }
                int expected = 0;
                for (int x = 0; x <= 2000; ++x) {
                    expected += survives[x];
                    CheckPhi(x, a, tables, expected, "exhaustive sieve");
                    ++exhaustive;
                }
            }
        }

        // Use the Legendre identity without calling the existing implementation.
        // Keep smaller contexts alive while larger temporary tables are released.
        struct KnownCount { std::int64_t n; int limit; std::int64_t expected; };
        const KnownCount knownCounts[] = {
            {1000000, 1000, 78498}, {100000000, 10000, 5761455},
            {10000000000LL, 100000, 455052511},
            {1000000000000LL, 1000000, 37607912018LL}
        };
        std::size_t known = 0;
        for (const auto& count : knownCounts) {
            if (noCache && count.n > 100000000) continue;
            const PrimeCounting::PhiTables tables(count.limit);
            const auto a = tables.primeCount();
            if (a != sieve.pi[count.limit]) {
                throw std::runtime_error("table prime count mismatch at limit=" +
                                         std::to_string(count.limit));
            }
            const auto actual = PrimeCounting::Phi(count.n, a, tables) + a - 1;
            if (actual != count.expected) {
                throw std::runtime_error("Legendre identity mismatch at n=" +
                    std::to_string(count.n) + ": returned " + std::to_string(actual) +
                    ", expected " + std::to_string(count.expected));
            }
            ++known;
        }

        std::size_t targeted = 0;
        for (int limit : {2, 13, 17, 19, 97, 100}) {
            const int a = sieve.pi[limit];
            const int p = sieve.primes[a - 1];
            targeted += CheckInputs(sieve, *contexts.at(limit), a,
                {p - 1, p, p + 1, p * p - 1, p * p, p * p + 1});
        }
        for (int p : {41, 43, 47}) {
            const int cube = p * p * p;
            for (int a = sieve.pi[p] - 1; a <= sieve.pi[p] + 1; ++a) {
                // Interleave contexts with identical canonical primes but different
                // lookup extents, including values beyond each small table.
                for (int limit : {10000, 97, 100, 97}) {
                    targeted += CheckInputs(sieve, *contexts.at(limit), a,
                                           {cube - 1, cube, cube + 1});
                }
            }
        }
        for (int a : {30, 31, 32}) {
            targeted += CheckInputs(sieve, *contexts.at(10000), a,
                                   {1270126, 1270127, 1270128});
        }
        for (int a : {1228, 1229}) {
            targeted += CheckInputs(sieve, *contexts.at(10000), a,
                                   {9972, 9973, 9974, 9999, 10000, 10001});
        }

        std::size_t large = 0;
        const std::int64_t maxX = 9007199254740991LL;
        const std::int64_t largeInputs[] = {
            maxX, maxX - 1, 4503599627370495LL, 1000000000000000LL
        };
        for (std::int64_t x : largeInputs) {
            for (const auto& context : contexts) {
                const auto& tables = *context.second;
                for (int a = 0; a <= std::min<std::int64_t>(8, tables.primeCount()); ++a) {
                    CheckPhi(x, a, tables, InclusionExclusion(x, a, sieve.primes),
                             "large inclusion-exclusion");
                    ++large;
                }
            }
        }

        std::size_t invalid = 0;
        const std::int64_t invalidLimits[] = {
            -1, 94906266, std::numeric_limits<std::int64_t>::max()
        };
        const std::int64_t invalidInputs[] = {
            -1, maxX + 1, std::numeric_limits<std::int64_t>::max()
        };
        for (std::int64_t limit : invalidLimits) {
            ExpectOutOfRange([&] { PrimeCounting::PhiTables tables(limit); }, "table limit");
            ++invalid;
        }
        for (const auto& context : contexts) {
            const auto& tables = *context.second;
            for (std::int64_t x : invalidInputs) {
                ExpectOutOfRange([&] { PrimeCounting::Phi(x, 0, tables); }, "x");
                ++invalid;
            }
            const std::int64_t invalidA[] = {
                -1, tables.primeCount() + 1, std::numeric_limits<std::int64_t>::max()
            };
            for (std::int64_t a : invalidA) {
                ExpectOutOfRange([&] { PrimeCounting::Phi(0, a, tables); }, "a");
                ++invalid;
            }
        }
        if (noCache && !rejectedAllocations.load(std::memory_order_relaxed)) {
            throw std::runtime_error("--no-cache did not exercise allocation failure");
        }
        std::cout << "PASS mode=" << (noCache ? "no-cache" : "cache")
                  << " exhaustive=" << exhaustive << " targeted=" << targeted
                  << " large=" << large << " concurrent=" << concurrent
                  << " known_pi=" << known
                  << " invalid=" << invalid << " rejected_allocations="
                  << rejectedAllocations.load(std::memory_order_relaxed) << "\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
