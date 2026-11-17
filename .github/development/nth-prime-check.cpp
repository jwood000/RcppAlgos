// Native checks for nth-prime selection, exact range updates and reused tables.
// Include the implementation to exercise private sieve anchors without adding
// a testing API to the package. Do not also link src/NthPrime.cpp.
// Link PiAuto.cpp, PiLegendre.cpp, PiMeissel.cpp, PiLehmer.cpp, PiPrime.cpp,
// Phi.cpp, PhiParallel.cpp, P2.cpp, P3.cpp and Eratosthenes.cpp with
// -std=c++11 -pthread.
#include "CppConvert/Constants.h"
#include "../../src/NthPrime.cpp"
#include <iostream>
#include <limits>
#include <set>
#include <string>

namespace {
    struct Reference {
        std::vector<std::int64_t> primes;
        std::vector<int> pi;

        explicit Reference(int maximum) : primes(1, 0), pi(maximum + 1, 0) {
            std::vector<unsigned char> prime(maximum + 1, 1);
            prime[0] = prime[1] = 0;
            for (int p = 2; p <= maximum / p; ++p) {
                if (!prime[p]) continue;
                for (int multiple = p * p; multiple <= maximum; multiple += p) {
                    prime[multiple] = 0;
                }
            }
            for (int n = 1; n <= maximum; ++n) {
                pi[n] = pi[n - 1] + prime[n];
                if (prime[n]) primes.push_back(n);
            }
        }
    };

    std::size_t checks = 0;

    void Check(std::int64_t actual, std::int64_t expected,
               const std::string& label) {
        ++checks;
        if (actual != expected) {
            throw std::runtime_error(label + ": got " + std::to_string(actual) +
                                     ", expected " + std::to_string(expected));
        }
    }

    template <typename Callable>
    void ExpectRange(Callable call, const char* label) {
        ++checks;
        try {
            call();
        } catch (const std::out_of_range&) {
            return;
        }
        throw std::runtime_error(std::string("Missing out_of_range: ") + label);
    }

    // Deliberately independent of the production floating-estimate root.
    int Root(int x) {
        int low = 0, high = x;
        while (low < high) {
            const int middle = low + (high - low + 1) / 2;
            if (middle <= x / middle) low = middle;
            else high = middle - 1;
        }
        return low;
    }
}

int main() {
    try {
        using namespace PrimeCounting;
        const int referenceMaximum = 3000000;
        const Reference reference(referenceMaximum);

        // Complete small rank range exercises lookup and estimate/count paths.
        for (std::int64_t rank = 1; rank <= 1000; ++rank) {
            Check(NthPrime(rank, 1), reference.primes[rank], "small nth prime");
        }
        for (std::int64_t rank : {1001, 1234, 10000, 50000, 78000}) {
            for (int workers : {1, 2, 4}) {
                Check(NthPrime(rank, workers), reference.primes[rank],
                      "larger nth prime");
            }
        }

        // Hardware inputs are fixed rather than inferred from the test host.
        // An oversized request must be capped before native sieve launch.
        struct Budget { int requested, maximum, cores; };
        const Budget budgets[] = {
            {1, 1, 1}, {4, 4, 2}, {8, 2, 1},
            {std::numeric_limits<int>::max(), 4, 1}
        };
        for (const auto& budget : budgets) {
            Check(NthPrime(78000, budget.requested, budget.maximum, budget.cores),
                  reference.primes[78000], "fixed native budgets");
        }

        // The magnitude policy is separate from numerical correctness. Check
        // every transition directly without sieving enormous intervals.
        struct LimitCase { std::int64_t x, expected; };
        const LimitCase limits[] = {
            {2, 1000000},
            {9999999999LL, 1000000},
            {10000000000LL, 10000000},
            {10000000001LL, 10000000},
            {999999999999LL, 10000000},
            {1000000000000LL, 100000000},
            {1000000000001LL, 100000000},
            {99999999999999LL, 100000000},
            {100000000000000LL, 1000000000},
            {100000000000001LL, 1000000000},
            {Significand53Int, 1000000000}
        };
        for (const auto& item : limits) {
            Check(sieveLimit(item.x), item.expected, "magnitude sieve limit");
        }

        // Every update starts from a known exact prefix count. Prime and
        // composite endpoints, a repeated anchor, and direction changes catch
        // inclusion mistakes in (anchor, next] and (next, anchor].
        const std::int64_t nextAnchors[] = {
            2, 3, 4, 2, 5, 7, 8, 97, 96, 100, 31, 30, 29, 31,
            65536, 65535, 65537, 999983, 999982, 1000000, referenceMaximum, 2
        };
        for (int workers : {1, 4}) {
            std::int64_t anchor = 2;
            std::int64_t count = reference.pi[anchor];
            for (std::int64_t next : nextAnchors) {
                count = countAt(next, anchor, count, workers, 4, 2);
                Check(count, reference.pi[next], "exact interval update");
                anchor = next;
            }
        }
        Check(countAt(Significand53Int, Significand53Int, MAX_NTH_PRIME_INDEX,
                      1, 1, 1), MAX_NTH_PRIME_INDEX,
              "unchanged maximum anchor needs no interval");

        // Force both directions, equality plateaus, and a final answer of 2.
        std::set<int> anchors;
        for (int x = 2; x <= 200; ++x) anchors.insert(x);
        for (int center : {31 * 31, 97 * 97, 257 * 257,
                           65536, 131072, 196608, 999983}) {
            for (int offset = -2; offset <= 2; ++offset) {
                anchors.insert(center + offset);
            }
        }
        for (int anchor : anchors) {
            const int count = reference.pi[anchor];
            std::set<int> ranks = {1, count - 1, count, count + 1, count + 2};
            ranks.insert(reference.pi[std::min(anchor + 131075, 1000000)]);
            for (int rank : ranks) {
                if (rank < 1 || rank >= static_cast<int>(reference.primes.size())) {
                    continue;
                }
                Check(sieveNthPrime(rank, anchor, count, 1, 1, 1),
                      reference.primes[rank],
                      "sieve rank=" + std::to_string(rank) +
                      " anchor=" + std::to_string(anchor));
            }
        }
        // These displacements exceed twice the final output-window cap.
        // Independent endpoint values verify several contiguous native-sieve
        // windows, in both directions and with actual parallel-size work.
        const std::int64_t finalRank = reference.primes.size() - 1;
        Check(reference.primes[finalRank] > 2 * FINAL_SIEVE_LIMIT, true,
              "multi-window fixture exceeds twice the cap");
        for (int workers : {1, 4}) {
            Check(sieveNthPrime(finalRank, 2, 1, workers, 4, 2),
                  reference.primes[finalRank], "multi-window forward");
            Check(sieveNthPrime(1, referenceMaximum, reference.pi[referenceMaximum],
                               workers, 4, 2), 2, "multi-window backward");
        }

        // Prefix counting can still retain a sufficiently large table.
        std::unique_ptr<PhiTables> tables;
        ensureTables(100, tables);
        const PhiTables* retained = tables.get();
        ensureTables(100, tables);
        Check(tables.get() == retained, true, "sufficient tables retained");
        std::unique_ptr<PhiTables> growing;
        ensureTables(100, growing);
        Check(growing->limit() >= 10, true, "initial table coverage");
        ensureTables(1000000, growing);
        Check(growing->limit() >= 1000, true, "grown table coverage");

        // A larger context must not change the mathematical cutoff, especially
        // in Legendre where table.primeCount() would be wrong for smaller x.
        const PhiTables oversized(2000);
        std::vector<int> probes;
        for (int x = 0; x <= 200; ++x) probes.push_back(x);
        for (int center : {31 * 31, 97 * 97, 257 * 257, 100000, 999983}) {
            for (int offset = -1; offset <= 1; ++offset) {
                probes.push_back(center + offset);
            }
        }
        for (int x : probes) {
            const PhiTables minimum(Root(x));
            for (int workers : {-2, 0, 1, 4}) {
                Check(PiLegendre(x, workers, minimum), reference.pi[x],
                      "minimal Legendre");
                Check(PiLegendre(x, workers, oversized), reference.pi[x],
                      "oversized Legendre");
                Check(PiMeissel(x, workers, minimum), reference.pi[x],
                      "minimal Meissel");
                Check(PiMeissel(x, workers, oversized), reference.pi[x],
                      "oversized Meissel");
                Check(PiLehmer(x, workers, minimum), reference.pi[x],
                      "minimal Lehmer");
                Check(PiLehmer(x, workers, oversized), reference.pi[x],
                      "oversized Lehmer");
                Check(PiAuto(x, workers, oversized), reference.pi[x],
                      "oversized auto");
            }
        }

        const PhiTables shortTables(30);
        ExpectRange([&] { PiLegendre(31 * 31, 1, shortTables); }, "short Legendre");
        ExpectRange([&] { PiMeissel(31 * 31, 1, shortTables); }, "short Meissel");
        ExpectRange([&] { PiLehmer(31 * 31, 1, shortTables); }, "short Lehmer");
        ExpectRange([&] { PiAuto(31 * 31, 1, shortTables); }, "short auto");
        for (std::int64_t bad : {std::int64_t(-1), Significand53Int + 1}) {
            ExpectRange([&] { PiLegendre(bad, 1, oversized); }, "invalid Legendre x");
            ExpectRange([&] { PiMeissel(bad, 1, oversized); }, "invalid Meissel x");
            ExpectRange([&] { PiLehmer(bad, 1, oversized); }, "invalid Lehmer x");
            ExpectRange([&] { PiAuto(bad, 1, oversized); }, "invalid auto x");
        }
        const PhiTables largeTables(1010000);
        for (int workers : {1, 4}) {
            Check(PiLegendre(100000000, workers, largeTables), 5761455,
                  "known large Legendre");
            Check(PiMeissel(100000000, workers, largeTables), 5761455,
                  "known large Meissel");
            Check(PiLehmer(100000000, workers, largeTables), 5761455,
                  "known large Lehmer");
        }
        Check(PiAuto(10000000000LL, 2, largeTables), 455052511,
              "known auto Meissel");
        Check(PiAuto(1000000000000LL, 6, largeTables), 37607912018LL,
              "known auto Lehmer");
        ExpectRange([&] { NthPrime(0, 1); }, "zero rank");
        ExpectRange([&] { NthPrime(-1, 1); }, "negative rank");
        ExpectRange([&] { NthPrime(MAX_NTH_PRIME_INDEX + 1, 1); }, "excess rank");

        std::cout << checks << " native nth-prime checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
