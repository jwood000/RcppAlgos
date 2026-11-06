// Independent native checks for the Meissel and Lehmer drivers.
// Compile with -std=c++11 -pthread -Iinst/include and link src/PiMeissel.cpp,
// src/PiLehmer.cpp, src/PiLegendre.cpp, src/PhiParallel.cpp, src/Phi.cpp,
// src/P2.cpp, src/P3.cpp, and src/Eratosthenes.cpp. Build separate executables
// with -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined and
// -fsanitize=thread, then run each without arguments.
#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/ThreadBatch.h"
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
    using CountFunction = std::int64_t (*)(std::int64_t, int, int);
    struct Method {
        const char* name;
        CountFunction count;
    };

    // An ordinary segmented Eratosthenes sieve supplies every selected small
    // and cutoff-boundary answer. It uses no production tables or Phi/P2/P3
    // identities. Only the requested checkpoints and one byte block persist.
    class Reference {
        std::map<std::int64_t, std::int64_t> counts;

    public:
        explicit Reference(std::vector<std::int64_t> checkpoints) {
            std::sort(checkpoints.begin(), checkpoints.end());
            checkpoints.erase(std::unique(checkpoints.begin(), checkpoints.end()),
                              checkpoints.end());
            const std::int64_t maximum = checkpoints.back();
            int root = 0;
            while (root + 1 <= maximum / (root + 1)) ++root;
            std::vector<unsigned char> composite(root + 1, 0);
            std::vector<int> primes;
            for (int p = 2; p <= root; ++p) {
                if (composite[p]) continue;
                primes.push_back(p);
                if (p <= root / p) {
                    for (int n = p * p; n <= root; n += p) composite[n] = 1;
                }
            }

            constexpr std::int64_t blockSize = 65536;
            std::size_t checkpoint = 0;
            std::int64_t count = 0;
            for (std::int64_t low = 0; low <= maximum; low += blockSize) {
                const std::int64_t high = std::min(maximum, low + blockSize - 1);
                composite.assign(static_cast<std::size_t>(high - low + 1), 0);
                for (int p : primes) {
                    if (std::int64_t(p) * p > high) break;
                    const std::int64_t first = std::max(std::int64_t(p) * p,
                        ((low + p - 1) / p) * p);
                    for (std::int64_t n = first; n <= high; n += p) {
                        composite[static_cast<std::size_t>(n - low)] = 1;
                    }
                }
                for (std::int64_t n = low; n <= high; ++n) {
                    if (n >= 2 && !composite[static_cast<std::size_t>(n - low)]) ++count;
                    if (checkpoint < checkpoints.size() && n == checkpoints[checkpoint]) {
                        counts.emplace(n, count);
                        ++checkpoint;
                    }
                }
            }
        }

        std::int64_t operator()(std::int64_t n) const { return counts.at(n); }
    };

    void Check(const Method& method, std::int64_t n, int requested, int maximum,
               std::int64_t expected, const char* category) {
        const auto actual = method.count(n, requested, maximum);
        if (actual != expected) {
            throw std::runtime_error(std::string(category) + " " + method.name +
                " at n=" + std::to_string(n) + " requested=" +
                std::to_string(requested) + " maximum=" + std::to_string(maximum) +
                ": actual=" + std::to_string(actual) +
                " expected=" + std::to_string(expected));
        }
    }

    struct ReleaseOnExit {
        std::atomic<bool>& flag;
        explicit ReleaseOnExit(std::atomic<bool>& value) : flag(value) {}
        ~ReleaseOnExit() { flag.store(true, std::memory_order_release); }
    };
}

int main() {
    try {
        const Method methods[] = {
            {"Meissel", PrimeCounting::PiMeissel},
            {"Lehmer", PrimeCounting::PiLehmer}
        };
        std::vector<std::int64_t> checkpoints, samples, boundaries;
        for (int n = 0; n <= 3000; ++n) checkpoints.push_back(n);
        for (int n : {100000, 100001, 999999, 1000000, 1234567, 2000000}) {
            samples.push_back(n);
        }
        std::uint32_t state = 0x73516u;
        for (int i = 0; i < 64; ++i) {
            state = state * 1664525u + 1013904223u;
            samples.push_back(3001 + state % 1997000u);
        }
        for (int p : {2, 3, 17, 97}) {
            std::int64_t power = p;
            for (int exponent = 2; exponent <= 4; ++exponent) {
                power *= p;
                for (int delta = -1; delta <= 1; ++delta) boundaries.push_back(power + delta);
            }
        }
        checkpoints.insert(checkpoints.end(), samples.begin(), samples.end());
        checkpoints.insert(checkpoints.end(), boundaries.begin(), boundaries.end());
        const Reference reference(checkpoints);

        std::size_t exhaustive = 0, sampled = 0, cutoff = 0;
        for (const auto& method : methods) {
            for (int n = 0; n <= 3000; ++n) {
                Check(method, n, 1, 1, reference(n), "exhaustive");
                ++exhaustive;
            }
            for (auto n : samples) {
                Check(method, n, 4, 4, reference(n), "sampled");
                ++sampled;
            }
            const int budgets[] = {1, 2, 4, 8};
            for (std::size_t i = 0; i < boundaries.size(); ++i) {
                const int threads = budgets[i % 4];
                Check(method, boundaries[i], threads, threads,
                      reference(boundaries[i]), "root cutoff");
                ++cutoff;
            }
        }

        struct Known { std::int64_t n, count; int threads; };
        const Known known[] = {
            {1000000LL, 78498LL, 1},
            {10000000LL, 664579LL, 2},
            {100000000LL, 5761455LL, 4},
            {1000000000LL, 50847534LL, 8},
            {10000000000LL, 455052511LL, 2},
            {100000000000LL, 4118054813LL, 4},
            {1000000000000LL, 37607912018LL, 8}
        };
        std::size_t knownChecks = 0;
        for (const auto& item : known) {
            for (const auto& method : methods) {
                Check(method, item.n, item.threads, item.threads, item.count, "known pi");
                ++knownChecks;
            }
        }

        // Exercise normalization of both controls without permitting an
        // unbounded OS-thread request. INT_MAX is always capped by four.
        const int threadBudgets[][2] = {
            {-3, 4}, {0, 4}, {1, 4}, {2, 4}, {4, 2}, {8, 4},
            {std::numeric_limits<int>::max(), 4}, {4, 0}, {4, -2}, {4, 1}
        };
        std::size_t budgetChecks = 0;
        for (const auto& method : methods) for (const auto& budget : threadBudgets) {
            Check(method, 10000000000LL, budget[0], budget[1], 455052511LL,
                  "thread clamp");
            ++budgetChecks;
        }

        // Legendre is an additional comparison inside its native domain;
        // each expected result was independently checked above.
        const Method legendre = {"Legendre", PrimeCounting::PiLegendre};
        std::size_t legendreChecks = 0;
        for (std::int64_t n : {100000LL, 100001LL, 1234567LL}) {
            Check(legendre, n, 3, 3, reference(n), "supplemental comparison");
            ++legendreChecks;
        }
        for (int index : {0, 2, 4}) {
            const auto& item = known[index];
            Check(legendre, item.n, 4, 4, item.count, "supplemental comparison");
            ++legendreChecks;
        }

        const std::int64_t maxX = 9007199254740991LL;
        const std::int64_t invalidInputs[] = {
            std::numeric_limits<std::int64_t>::min(), -1, maxX + 1,
            std::numeric_limits<std::int64_t>::max()
        };
        std::size_t invalid = 0, recovery = 0;
        for (const auto& method : methods) for (auto n : invalidInputs) {
            bool caught = false;
            try { method.count(n, 4, 4); }
            catch (const std::out_of_range&) { caught = true; }
            if (!caught) {
                throw std::runtime_error(std::string("missing out_of_range for ") +
                    method.name + " at n=" + std::to_string(n));
            }
            ++invalid;
            Check(method, 1000000, 2, 2, 78498, "invalid-input recovery");
            ++recovery;
        }

        // Each call prepares a different context while sharing only the Phi
        // cache. The outer gate encourages overlap; inner worker counts stay
        // small. ReleaseOnExit prevents a partial launch failure from leaving
        // an earlier worker waiting while the batch destructor joins it.
        std::atomic<bool> start(false);
        std::int64_t concurrentResults[2] = {0, 0};
        PrimeCounting::ThreadBatch batch(2);
        ReleaseOnExit release(start);
        batch.launch([&] {
            while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
            concurrentResults[0] = PrimeCounting::PiMeissel(1000000000LL, 2, 2);
        });
        batch.launch([&] {
            while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
            concurrentResults[1] = PrimeCounting::PiLehmer(10000000000LL, 3, 3);
        });
        start.store(true, std::memory_order_release);
        batch.finish();
        if (concurrentResults[0] != 50847534 || concurrentResults[1] != 455052511) {
            throw std::runtime_error("concurrent methods mismatch");
        }

        std::cout << "PASS exhaustive=" << exhaustive << " samples=" << sampled
                  << " cutoff_boundaries=" << cutoff << " known=" << knownChecks
                  << " thread_budgets=" << budgetChecks << " legendre=" << legendreChecks
                  << " invalid=" << invalid << " recovery=" << recovery
                  << " concurrent=2\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
