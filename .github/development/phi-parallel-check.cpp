// Focused native checks for PhiParallel and its generic batch helper.
// Compile with -std=c++11 -pthread -Iinst/include and link src/Phi.cpp,
// src/PhiParallel.cpp, and src/Eratosthenes.cpp. Use separate ASan/UBSan and
// TSan executables; this checker neither needs R nor benchmarks performance.
#include "NumbersUtils/Phi.h"
#include "NumbersUtils/ThreadBatch.h"
#include <algorithm>
#include <array>
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
    void Check(bool condition, const char* message) {
        if (!condition) throw std::runtime_error(message);
    }

    struct IndependentSieve {
        std::vector<int> primes, pi;
        explicit IndependentSieve(int maximum) : pi(maximum + 1, 0) {
            std::vector<char> prime(maximum + 1, 1);
            prime[0] = prime[1] = 0;
            for (int p = 2; p * p <= maximum; ++p) if (prime[p]) {
                for (int x = p * p; x <= maximum; x += p) prime[x] = 0;
            }
            for (int x = 1; x <= maximum; ++x) {
                pi[x] = pi[x - 1] + prime[x];
                if (prime[x]) primes.push_back(x);
            }
        }
    };

    void CheckValue(std::int64_t x, std::int64_t a, int workers,
                    const PrimeCounting::PhiTables& tables, std::int64_t expected) {
        const auto actual = PrimeCounting::PhiParallel(x, a, workers, tables);
        if (actual != expected) {
            throw std::runtime_error("PhiParallel mismatch: x=" + std::to_string(x) +
                " a=" + std::to_string(a) + " workers=" + std::to_string(workers) +
                " limit=" + std::to_string(tables.limit()) + " actual=" +
                std::to_string(actual) + " expected=" + std::to_string(expected));
        }
    }

    std::array<std::int64_t, 3> PartialReference(const std::vector<int>& primes, int a) {
        const int maximum = 10000001;
        std::vector<char> survives(maximum + 1, 1);
        survives[0] = 0;
        for (int i = 0; i < a; ++i) {
            for (int x = primes[i]; x <= maximum; x += primes[i]) survives[x] = 0;
        }
        std::int64_t count = 0;
        std::array<std::int64_t, 3> result;
        for (int x = 1; x <= maximum; ++x) {
            count += survives[x];
            if (x >= maximum - 2) result[x - (maximum - 2)] = count;
        }
        return result;
    }

    std::int64_t InclusionExclusion(std::int64_t x, int a,
                                    const std::vector<int>& primes) {
        // At most eight primes: products and the signed running sum fit int64.
        std::int64_t result = x;
        for (unsigned mask = 1; mask < (1u << a); ++mask) {
            std::int64_t product = 1;
            int selected = 0;
            for (int i = 0; i < a; ++i) if (mask & (1u << i)) {
                product *= primes[i];
                ++selected;
            }
            result += (selected % 2 ? -1 : 1) * (x / product);
        }
        return result;
    }

    int IntegerRoot(std::int64_t n, int power) {
        int lower = 0, upper = 100000;
        while (lower < upper) {
            const int middle = lower + (upper - lower + 1) / 2;
            std::int64_t product = 1;
            bool fits = true;
            for (int i = 0; i < power; ++i) {
                if (product > n / middle) { fits = false; break; }
                product *= middle;
            }
            if (fits) lower = middle;
            else upper = middle - 1;
        }
        return lower;
    }

    template <typename Callable>
    void ExpectOutOfRange(Callable call) {
        try { call(); }
        catch (const std::out_of_range&) { return; }
        throw std::runtime_error("invalid PhiParallel input did not throw out_of_range");
    }

    struct ReleaseOnExit {
        std::atomic<bool>& flag;
        explicit ReleaseOnExit(std::atomic<bool>& value) : flag(value) {}
        ~ReleaseOnExit() { flag.store(true, std::memory_order_release); }
    };

    void CheckGenericBatch() {
        using PrimeCounting::ThreadBatch;
        {
            ThreadBatch empty(0);
            empty.finish();
        }
        {
            std::string text;
            std::vector<int> squares(4, 0);
            ThreadBatch batch(2);
            batch.launch([&] { text = "ordinary callable"; });
            batch.launch([&] {
                for (int i = 0; i < 4; ++i) squares[i] = (i + 1) * (i + 1);
            });
            batch.finish();
            Check(text == "ordinary callable" && squares == std::vector<int>({1, 4, 9, 16}),
                  "generic batch result mismatch");
        }
        {
            std::atomic<bool> release(false), completed(false);
            bool caught = false;
            try {
                ThreadBatch batch(1);
                // This guard releases the worker before batch destruction on
                // overflow, so the catch can verify its joined completion.
                ReleaseOnExit releaseOnExit(release);
                batch.launch([&] {
                    while (!release.load(std::memory_order_acquire)) std::this_thread::yield();
                    completed.store(true, std::memory_order_release);
                });
                batch.launch([] {});
            } catch (const std::length_error&) {
                caught = true;
                Check(completed.load(std::memory_order_acquire),
                      "overflow escaped before the started worker completed");
            }
            Check(caught, "generic batch accepted work beyond its capacity");
        }
        {
            int result = 0;
            bool caught = false;
            ThreadBatch batch(2);
            batch.launch([] { throw std::runtime_error("generic worker failure"); });
            batch.launch([&] { result = 42; });
            try { batch.finish(); }
            catch (const std::runtime_error& error) {
                caught = std::string(error.what()) == "generic worker failure";
            }
            Check(caught && result == 42, "generic worker exception/result propagation failed");
        }
    }
}

int main() {
    try {
        using PrimeCounting::Phi;
        using PrimeCounting::PhiParallel;
        using PrimeCounting::PhiTables;
        const IndependentSieve sieve(100000);
        std::map<int, std::unique_ptr<PhiTables>> tables;
        for (int limit : {0, 1, 2, 13, 17, 19, 100, 10000, 100000}) {
            tables.emplace(limit, std::unique_ptr<PhiTables>(new PhiTables(limit)));
        }

        std::size_t small = 0;
        const int budgets[] = {-3, 0, 1, 2, 3, 4, 8};
        for (int limit : {0, 1, 2, 13, 17, 19, 100}) {
            std::vector<char> survives(129, 1);
            survives[0] = 0;
            for (int a = 0; a <= sieve.pi[limit]; ++a) {
                if (a) for (int x = sieve.primes[a - 1]; x <= 128; x += sieve.primes[a - 1]) {
                    survives[x] = 0;
                }
                int expected = 0;
                for (int x = 0; x <= 128; ++x) {
                    expected += survives[x];
                    CheckValue(x, a, budgets[(x + a) % 7], *tables.at(limit), expected);
                    ++small;
                }
            }
        }

        const int fullA = sieve.pi[3162];
        std::map<int, std::array<std::int64_t, 3>> nearReference;
        for (int a : {7, 8, 25, 64, fullA}) {
            nearReference.emplace(a, PartialReference(sieve.primes, a));
        }
        struct NearCase { int x, a, workers, limit; };
        const NearCase nearCases[] = {
            {9999999, 25, 4, 10000}, {10000000, 7, 2, 17},
            {10000001, 8, 3, 19}, {10000000, 25, 4, 10000},
            {10000001, 64, 8, 10000}, {10000000, fullA, 2, 10000},
            {10000001, 25, 0, 10000}, {10000001, 25, -3, 10000},
            {10000001, 25, 1, 10000}, {10000001, 7, 8, 10000}
        };
        for (const auto& c : nearCases) {
            CheckValue(c.x, c.a, c.workers, *tables.at(c.limit),
                       nearReference.at(c.a)[c.x - 9999999]);
        }

        struct LargeCase { std::int64_t x; int a, workers; };
        const std::int64_t maxX = 9007199254740991LL;
        const LargeCase largeCases[] = {
            {maxX, 0, 8}, {maxX, 1, 8}, {maxX, 6, 8}, {maxX, 7, 8},
            {maxX, 8, 8}, {maxX - 1, 8, 3},
            {10000000001LL, 8, 4}, {10000000001LL, 7, 2}
        };
        for (const auto& c : largeCases) {
            CheckValue(c.x, c.a, c.workers, *tables.at(10000),
                       InclusionExclusion(c.x, c.a, sieve.primes));
        }

        // Requested cutoffs are much smaller than the prepared prime count.
        // A small fixed set limits actual worker creation, unlike a full grid.
        std::size_t differential = 0;
        const std::int64_t differentialInputs[] = {100000000, 10000000000LL};
        for (std::int64_t x : differentialInputs) for (int power : {3, 4}) {
            const int center = sieve.pi[IntegerRoot(x, power)];
            for (int delta = -1; delta <= 1; ++delta) {
                const int a = center + delta;
                CheckValue(x, a, delta + 3, *tables.at(100000),
                           Phi(x, a, *tables.at(100000)));
                ++differential;
            }
        }
        for (int power : {3, 4}) {
            const std::int64_t x = 100000000000LL;
            const int a = sieve.pi[IntegerRoot(x, power)];
            CheckValue(x, a, 3, *tables.at(100000), Phi(x, a, *tables.at(100000)));
            ++differential;
        }

        {
            std::int64_t results[2] = {0, 0};
            std::atomic<bool> start(false);
            PrimeCounting::ThreadBatch callers(2);
            ReleaseOnExit releaseOnExit(start);
            callers.launch([&] {
                while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                results[0] = PhiParallel(10000000, 8, 2, *tables.at(19));
            });
            callers.launch([&] {
                while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                results[1] = PhiParallel(10000001, 25, 3, *tables.at(10000));
            });
            start.store(true, std::memory_order_release);
            callers.finish();
            Check(results[0] == nearReference.at(8)[1] && results[1] == nearReference.at(25)[2],
                  "concurrent PhiParallel caller mismatch");
        }

        std::size_t invalid = 0;
        const std::int64_t invalidX[] = {-1, maxX + 1, std::numeric_limits<std::int64_t>::max()};
        for (int limit : {0, 17, 100000}) {
            const auto& prepared = *tables.at(limit);
            for (std::int64_t x : invalidX) {
                ExpectOutOfRange([&] { PhiParallel(x, 0, 8, prepared); });
                ++invalid;
            }
            const std::int64_t invalidA[] = {-1, prepared.primeCount() + 1,
                                            std::numeric_limits<std::int64_t>::max()};
            for (std::int64_t a : invalidA) {
                ExpectOutOfRange([&] { PhiParallel(0, a, 8, prepared); });
                ++invalid;
            }
        }
        CheckGenericBatch();
        std::cout << "PASS small=" << small << " near_threshold="
                  << sizeof(nearCases) / sizeof(nearCases[0]) << " large_exact="
                  << sizeof(largeCases) / sizeof(largeCases[0])
                  << " differential=" << differential
                  << " concurrent=2 invalid=" << invalid << " generic_batch=4\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
