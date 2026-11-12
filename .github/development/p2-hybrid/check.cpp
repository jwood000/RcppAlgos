// Independent checks for the temporary bounded-stream P2 candidates.
// Link production P2.cpp once with -DP2=P2Reference, then link either the
// unchanged P2.cpp or P2Hybrid.cpp normally. No production hooks are needed.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/Phi.h"
#include "NumbersUtils/PrimeCount.h"
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

#ifdef P2_HYBRID_STREAM_CHECK
// Test-only inclusion exposes the private stream within this translation unit.
// No production hook or full P2 sweep through these large endpoints is needed.
#include "P2Hybrid.cpp"
#endif

namespace PrimeCounting {
    std::int64_t P2Reference(std::int64_t, std::int64_t, int, const PhiTables&);
}

namespace {
    std::size_t checks = 0;

    void Check(std::int64_t actual, std::int64_t expected, const char* label,
               std::int64_t x, std::int64_t a, int workers) {
        ++checks;
        if (actual != expected) {
            throw std::runtime_error(std::string(label) + " mismatch at x=" +
                std::to_string(x) + " a=" + std::to_string(a) + " workers=" +
                std::to_string(workers) + ": actual=" + std::to_string(actual) +
                " expected=" + std::to_string(expected));
        }
    }

    // The oracle root uses binary search and exact division, independently of
    // the production floating-point estimate and correction loop.
    int Root(std::int64_t x) {
        int low = 0, high = static_cast<int>(std::min<std::int64_t>(x, 94906265));
        while (low < high) {
            const int middle = low + (high - low + 1) / 2;
            if (middle <= x / middle) low = middle;
            else high = middle - 1;
        }
        return low;
    }

    struct Reference {
        std::vector<int> primes, pi;
        std::vector<std::vector<int>> smallCounts;

        Reference() : pi(9000001, 0) {
            const int maximum = static_cast<int>(pi.size()) - 1;
            std::vector<unsigned char> sieve(maximum + 1, 1);
            sieve[0] = sieve[1] = 0;
            for (int p = 2; p <= maximum / p; ++p) if (sieve[p]) {
                for (int n = p * p; n <= maximum; n += p) sieve[n] = 0;
            }
            for (int n = 1; n <= maximum; ++n) {
                pi[n] = pi[n - 1] + sieve[n];
                if (sieve[n]) primes.push_back(n);
            }

            // Trial division supplies multiplicities and least factors for
            // the exhaustive oracle. It does not use the native Phi/P2 code
            // or the prime-count summation used by the larger reference.
            std::vector<int> omega(513, 0), least(513, 0);
            for (int n = 2; n <= 512; ++n) {
                int remaining = n;
                for (int divisor = 2; divisor <= remaining / divisor; ++divisor) {
                    while (remaining % divisor == 0) {
                        if (!least[n]) least[n] = divisor;
                        ++omega[n];
                        remaining /= divisor;
                    }
                }
                if (remaining > 1) {
                    if (!least[n]) least[n] = remaining;
                    ++omega[n];
                }
            }
            smallCounts.resize(pi[512] + 1, std::vector<int>(513, 0));
            for (int a = 0; a <= pi[512]; ++a) {
                const int cutoff = a ? primes[a - 1] : 1;
                for (int n = 1; n <= 512; ++n) {
                    smallCounts[a][n] = smallCounts[a][n - 1] +
                        (omega[n] == 2 && least[n] > cutoff);
                }
            }
        }

        std::int64_t pairs(std::int64_t x, int a) const {
            std::int64_t sum = 0;
            for (std::size_t i = static_cast<std::size_t>(a);
                 i < primes.size() && primes[i] <= x / primes[i]; ++i) {
                const auto q = x / primes[i];
                if (q >= static_cast<std::int64_t>(pi.size())) {
                    throw std::runtime_error("independent pi oracle too short");
                }
                sum += pi[q] - static_cast<std::int64_t>(i);
            }
            return sum;
        }
    };

    template <typename Callable>
    void ExpectRange(Callable call, const char* label) {
        ++checks;
        try { call(); }
        catch (const std::out_of_range&) { return; }
        throw std::runtime_error(std::string("missing out_of_range: ") + label);
    }

    struct ReleaseOnExit {
        std::atomic<bool>& flag;
        explicit ReleaseOnExit(std::atomic<bool>& value) : flag(value) {}
        ~ReleaseOnExit() { flag.store(true, std::memory_order_release); }
    };

#ifdef P2_HYBRID_STREAM_CHECK
    void CheckStreamWindows(const Reference& reference) {
        const std::int64_t finishes[] = {
            32769LL * 32769 - 1, 32769LL * 32769,
            2147483648LL + 4096, 128LL * 94906265
        };
        for (const auto finish : finishes) {
            const int root = Root(finish);
            const std::int64_t buffer = std::max(32768, root);
            const std::int64_t start = finish - 2 * buffer - 7;
            const PrimeCounting::PhiTables tables(root);
            // Only the increment is relevant to this private stream test.
            // A synthetic base also makes the count cross INT32_MAX.
            const std::int64_t seed = 2147483647;
            PrimeCounting::PrimeCountStream stream(start, finish, seed, tables);
            const std::int64_t firstOdd = start + (start % 2 == 0 ? 1 : 2);
            const std::int64_t nextBlock = firstOdd + 2 * buffer;
            std::vector<std::int64_t> targets = {
                start, start + 1, start + 2, nextBlock - 2, nextBlock - 1,
                nextBlock, nextBlock + 1, nextBlock + 2, finish - 1, finish
            };
            for (std::int64_t value : {2147483646LL, 2147483647LL,
                                       2147483648LL, 2147483649LL}) {
                if (value >= start && value <= finish) targets.push_back(value);
            }
            std::sort(targets.begin(), targets.end());
            targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
            std::int64_t next = start + 1, expected = seed;
            for (const auto target : targets) {
                for (; next <= target; ++next) {
                    bool prime = next >= 2;
                    for (const auto divisor : reference.primes) {
                        if (divisor > next / divisor) break;
                        if (next % divisor == 0) { prime = false; break; }
                    }
                    expected += prime;
                }
                Check(stream.countTo(target), expected, "large stream window",
                      target, 0, 1);
                Check(stream.countTo(target), expected, "repeated stream endpoint",
                      target, 0, 1);
            }
        }
    }
#endif
}

int main() {
    try {
        using PrimeCounting::P2;
        using PrimeCounting::P2Reference;
        using PrimeCounting::PhiTables;
        const Reference reference;
#ifdef P2_HYBRID_STREAM_CHECK
        CheckStreamWindows(reference);
        std::cout << "PASS stream_checks=" << checks << '\n';
#else
        const int budgets[] = {-2, 0, 1, 2, 8};
        std::map<int, std::unique_ptr<PhiTables>> smallTables;
        const auto smallTable = [&](int limit) -> const PhiTables& {
            auto found = smallTables.find(limit);
            if (found == smallTables.end()) {
                found = smallTables.emplace(limit,
                    std::unique_ptr<PhiTables>(new PhiTables(limit))).first;
            }
            return *found->second;
        };

        // Minimal contexts exercise the streaming path; a fixed table through
        // 512 alone would answer every small quotient directly and miss it.
        for (int x = 0; x <= 512; ++x) {
            const int limits[] = {Root(x), Root(x) + 17, 512};
            for (int limit : limits) {
                const auto& tables = smallTable(limit);
                for (int a = 0; a <= tables.primeCount(); ++a) {
                    const int workers = budgets[(x + a + limit) % 5];
                    Check(P2(x, a, workers, tables), reference.smallCounts[a][x],
                          "small factors", x, a, workers);
                }
            }
        }
        std::cout << "small_checks=" << checks << '\n';

        // Exercise exact q=L and q=mL split decisions, including both ends
        // of each floor-division plateau. Each n fits the stated table limit.
        struct Endpoint { int limit, prime; };
        const Endpoint endpoints[] = {{128, 11}, {256, 17}, {1000, 113}};
        for (const auto& item : endpoints) {
            const PhiTables tables(item.limit);
            const int index = reference.pi[item.prime];
            for (int multiplier : {1, 2, 4, 8}) {
                for (int delta = -1; delta <= 1; ++delta) {
                    for (int remainder : {0, item.prime - 1}) {
                        const auto x = std::int64_t(item.prime) *
                            (multiplier * item.limit + delta) + remainder;
                        for (int a : {0, index - 1, index}) {
                            Check(P2(x, a, 8, tables), reference.pairs(x, a),
                                  "quotient split", x, a, 8);
                        }
                    }
                }
            }
        }

        // Wider experimental cutoffs need larger prepared limits: requiring
        // L > 32*p keeps these n=p*(mL+delta)+remainder cases below L^2.
        // Even the a=0 quotient stays within the independent 9e6 pi oracle.
        const Endpoint wideEndpoints[] = {{512, 11}, {1024, 17}, {4096, 113}};
        for (const auto& item : wideEndpoints) {
            const PhiTables tables(item.limit);
            const int index = reference.pi[item.prime];
            for (int multiplier : {16, 32}) {
                for (int delta = -1; delta <= 1; ++delta) {
                    for (int remainder : {0, item.prime - 1}) {
                        const auto x = std::int64_t(item.prime) *
                            (multiplier * item.limit + delta) + remainder;
                        for (int a : {0, index - 1, index}) {
                            Check(P2(x, a, 8, tables), reference.pairs(x, a),
                                  "wide quotient split", x, a, 8);
                        }
                    }
                }
            }
        }

        // Preserve the earlier boundaries and add valid 64L/128L cases.
        // These choices have L > 128*p, while x/2 remains below 9e6.
        const Endpoint extendedEndpoints[] = {{2048, 11}, {4096, 17}, {4096, 31}};
        for (const auto& item : extendedEndpoints) {
            const PhiTables tables(item.limit);
            const int index = reference.pi[item.prime];
            for (int multiplier : {64, 128}) {
                for (int delta = -1; delta <= 1; ++delta) {
                    for (int remainder : {0, item.prime - 1}) {
                        const auto x = std::int64_t(item.prime) *
                            (multiplier * item.limit + delta) + remainder;
                        for (int a : {0, index - 1, index}) {
                            Check(P2(x, a, 8, tables), reference.pairs(x, a),
                                  "extended quotient split", x, a, 8);
                        }
                    }
                }
            }
        }

        // The initial 32768-odd block ends near L+65536. Targets include an
        // even endpoint and the odd endpoint that requires a second block.
        {
            const int limit = 100000, prime = 10007;
            const PhiTables tables(limit);
            const int a = reference.pi[prime] - 1;
            for (int delta = -2; delta <= 2; ++delta) {
                const auto x = std::int64_t(prime) * (limit + 65536 + delta);
                Check(P2(x, a, 8, tables), reference.pairs(x, a),
                      "sieve block endpoint", x, a, 8);
            }
        }

        for (int prime : {17, 37, 127, 257}) {
            for (int power = 2; power <= 3; ++power) {
                std::int64_t value = std::int64_t(prime) * prime;
                if (power == 3) value *= prime;
                for (int delta = -1; delta <= 1; ++delta) {
                    const auto x = value + delta;
                    const PhiTables tables(Root(x));
                    const int a = reference.pi[prime] - 1;
                    Check(P2(x, a, 8, tables), reference.pairs(x, a),
                          "prime power", x, a, 8);
                }
            }
        }

        // Larger independent pi counts exercise real stream workers for all
        // multipliers. The cutoff keeps every quotient within the oracle.
        {
            const std::int64_t x = 10000000000000LL;
            const int a = reference.pi[1200000];
            const PhiTables tables(Root(x));
            const auto expected = reference.pairs(x, a);
            for (int workers : budgets) {
                Check(P2(x, a, workers, tables), expected,
                      "large independent stream", x, a, workers);
            }
            Check(P2Reference(x, a, 1, tables), expected,
                  "baseline independent stream", x, a, 1);
        }

        // Both a larger and a smaller cutoff, and oversized prepared tables,
        // test the combination of unchanged Phi-prefix and stream-tail work.
        for (const std::int64_t x : {100000000LL, 10000000000LL, 1000000000000LL}) {
            const int root = Root(x);
            for (int limit : {root, root * 2}) {
                const PhiTables tables(limit);
                const int cutoffs[] = {
                    tables.pi()[PrimeCounting::PartialSieveDetail::IntegerCbrt(x)],
                    tables.pi()[Root(root)], tables.pi()[root] - 1,
                    tables.pi()[root]
                };
                for (int a : cutoffs) {
                    const auto expected = P2Reference(x, a, 1, tables);
                    for (int workers : {1, 2, 8}) {
                        Check(P2(x, a, workers, tables), expected,
                              "baseline differential", x, a, workers);
                    }
                }
            }
        }

        // Two independent calls share immutable tables and the process-wide
        // Phi cache, while their stream buffers and subtotal slots stay local.
        {
            const std::int64_t x = 10000000000LL;
            const PhiTables tables(Root(x));
            const int a = reference.pi[2154];
            const auto expected = P2Reference(x, a, 1, tables);
            std::int64_t sums[2] = {0, 0};
            std::atomic<bool> start(false);
            PrimeCounting::ThreadBatch callers(2);
            ReleaseOnExit release(start);
            for (int slot = 0; slot < 2; ++slot) {
                callers.launch([&, slot] {
                    while (!start.load(std::memory_order_acquire)) std::this_thread::yield();
                    sums[slot] = P2(x, a, 2, tables);
                });
            }
            start.store(true, std::memory_order_release);
            callers.finish();
            Check(sums[0], expected, "concurrent call 0", x, a, 2);
            Check(sums[1], expected, "concurrent call 1", x, a, 2);
        }

        for (const int limit : {0, 17, 512}) {
            const PhiTables tables(limit);
            for (std::int64_t x : {-1LL, 9007199254740992LL,
                                   std::numeric_limits<std::int64_t>::max()}) {
                ExpectRange([&] { P2(x, 0, 8, tables); }, "invalid x");
            }
            for (std::int64_t a : {-1LL, tables.primeCount() + 1}) {
                ExpectRange([&] { P2(0, a, 8, tables); }, "invalid a at zero");
            }
        }
        const PhiTables shortTable(10);
        ExpectRange([&] { P2(121, 4, 8, shortTable); }, "insufficient root table");

        struct Known { std::int64_t x, pi; };
        const Known known[] = {{1000000, 78498}, {10000000000LL, 455052511},
                               {1000000000000LL, 37607912018LL}};
        for (const auto& value : known) {
            for (int workers : {1, 8}) {
                Check(PrimeCounting::PiMeissel(value.x, workers, workers), value.pi,
                      "complete Meissel", value.x, 0, workers);
                Check(PrimeCounting::PiLehmer(value.x, workers, workers), value.pi,
                      "complete Lehmer", value.x, 0, workers);
            }
        }
        std::cout << "PASS checks=" << checks << '\n';
#endif
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
