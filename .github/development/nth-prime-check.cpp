// Native checks for nth-prime selection and reusable prime-count tables.
// Include the implementation to exercise private sieve anchors without adding
// a testing API to the package. Do not also link src/NthPrime.cpp.
// Link PiAuto.cpp, PiLegendre.cpp, PiMeissel.cpp, PiLehmer.cpp, PiPrime.cpp,
// Phi.cpp, PhiParallel.cpp, P2.cpp, P3.cpp and Eratosthenes.cpp with
// -std=c++11 -pthread.
#include "CppConvert/Constants.h"
#include "../../src/NthPrime.cpp"
#include <iostream>
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
        const Reference reference(1000000);

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

        // Force both directions, equality plateaus, and a final answer of 2.
        // Reuse one table owner between calls to also exercise retained setup.
        std::unique_ptr<PhiTables> tables;
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
                Check(sieveNthPrime(rank, anchor, count, tables),
                      reference.primes[rank],
                      "sieve rank=" + std::to_string(rank) +
                      " anchor=" + std::to_string(anchor));
            }
        }
        // Large displacements force several whole sieve blocks in each direction.
        Check(sieveNthPrime(78000, 2, 1, tables), reference.primes[78000],
              "multi-block forward");
        Check(sieveNthPrime(1, 1000000, reference.pi[1000000], tables), 2,
              "multi-block backward");

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
