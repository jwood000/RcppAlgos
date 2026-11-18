// Standalone C++11 regression checks for the shared counting/generation sieve.
// Compile with -Iinst/include and src/PiPrime.cpp; Eratosthenes.cpp is included
// here so kernel transitions can be exercised without production test hooks.
// Add --wide for the >1e9 interval-window check (omit in sanitizer runs).
#include "NumbersUtils/PrimeCount.h"
#include "NumbersUtils/Eratosthenes.h"
#include "../../src/Eratosthenes.cpp"
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace {
    using Integer = std::int64_t;
    constexpr Integer MAX_VALUE = (Integer(1) << 53) - 1;
    // Keep an independent expression to catch changes to the shared bound.
    static_assert(Significand53Int == MAX_VALUE, "53-bit integer bound changed");
    std::size_t assertions = 0;

    void Check(bool condition, const std::string& message) {
        ++assertions;
        if (!condition) throw std::runtime_error(message);
    }

    std::string Interval(Integer low, Integer high) {
        return "[" + std::to_string(low) + ", " + std::to_string(high) + "]";
    }

    // Independent ordinary Eratosthenes sieve, followed by a byte-per-integer
    // segmented sieve. No wheels, package prime tables, or shared kernels.
    class Reference {
        std::vector<int> primes;
        int limit;
    public:
        explicit Reference(int maximum) : limit(maximum) {
            std::vector<unsigned char> sieve(maximum + 1, 1);
            sieve[0] = sieve[1] = 0;
            for (int p = 2; p <= maximum / p; ++p) {
                if (sieve[p]) for (int j = p * p; j <= maximum; j += p) sieve[j] = 0;
            }
            for (int p = 2; p <= maximum; ++p) if (sieve[p]) primes.push_back(p);
        }

        std::vector<Integer> range(Integer low, Integer high) const {
            if (high < low || high < 2) return {};
            low = std::max<Integer>(2, low);
            Check(high / limit < limit, "reference sieve lacks base-prime coverage");
            std::vector<unsigned char> sieve(high - low + 1, 1);
            for (Integer p : primes) {
                if (p > high / p) break;
                const Integer start = std::max(p * p, ((low + p - 1) / p) * p);
                for (Integer j = start; j <= high; j += p) sieve[j - low] = 0;
            }
            std::vector<Integer> result;
            for (Integer j = low; j <= high; ++j) if (sieve[j - low]) result.push_back(j);
            return result;
        }

        std::vector<int> base(Integer maximum) const {
            std::vector<int> result;
            for (int p : primes) {
                if (p == 2) continue;
                result.push_back(p);
                if (p > maximum / p) return result;
            }
            throw std::runtime_error("reference lacks extra base prime");
        }
    };

    template <typename T>
    std::vector<T> Generate(Integer low, Integer high, int threads = 1,
                            int maxThreads = 1, int maxCores = 1,
                            bool* usedParallel = nullptr) {
        std::vector<std::vector<T>> chunks(threads);
        std::vector<T> result;
        bool parallel = false;
        PrimeSieve::PrimeSieveMain(chunks, result, low, high, parallel,
                                   threads, maxThreads, maxCores);
        if (parallel) for (const auto& chunk : chunks) {
            result.insert(result.end(), chunk.begin(), chunk.end());
        }
        if (usedParallel) *usedParallel = parallel;
        Check(std::is_sorted(result.begin(), result.end()), "generation is not sorted");
        Check(std::adjacent_find(result.begin(), result.end()) == result.end(),
              "generation contains duplicate primes");
        return result;
    }

    template <typename T>
    void CheckGenerated(const std::vector<T>& actual,
                        const std::vector<Integer>& expected) {
        Check(actual.size() == expected.size() &&
              std::equal(actual.begin(), actual.end(), expected.begin()),
              "generated values differ from independent sieve");
    }

    void CheckRange(const Reference& reference, Integer low, Integer high,
                    bool generate = false) {
        const auto expected = reference.range(low, high);
        Check(PrimeCounting::PiPrimeRange(low, high) == Integer(expected.size()),
              "range count differs for " + Interval(low, high));
        // This bypasses the precomputed-table shortcut in PiPrimeRange.
        if (low >= 1 && low <= high) {
            Check(PrimeSieve::CountPrimes(low, high) == Integer(expected.size()),
                  "direct count differs for " + Interval(low, high));
        }
        if (generate && low >= 1 && low <= high) {
            CheckGenerated(Generate<Integer>(low, high), expected);
        }
    }

    template <typename Function>
    void CheckOutOfRange(Function work) {
        bool caught = false;
        try { work(); } catch (const std::out_of_range&) { caught = true; }
        Check(caught, "invalid bounds did not throw std::out_of_range");
    }

    void SmallAndBoundaryChecks(const Reference& reference) {
        const auto primes = reference.range(0, 4096);
        for (Integer n = 0; n <= 4096; ++n) {
            const Integer expected = std::upper_bound(primes.begin(), primes.end(), n) - primes.begin();
            Check(PrimeCounting::PiPrime(n) == expected, "small prefix mismatch");
            CheckRange(reference, n, n);
        }
        for (Integer center : {2310, 30030, 31627, 32340, 64680, 99999, 100000}) {
            for (int offset = -2; offset <= 2; ++offset) {
                CheckRange(reference, center + offset, center + offset, true);
                CheckRange(reference, center - 37, center + offset, true);
                CheckRange(reference, center + offset, center + 73, true);
            }
            Check(PrimeCounting::PiPrime(center) == Integer(reference.range(0, center).size()),
                  "prefix differs across table/segment boundary");
        }
        Check(PrimeCounting::PiPrime(1000000) == 78498, "known million prefix mismatch");
        Check(PrimeCounting::PiPrime(10000000) == 664579, "known ten-million prefix mismatch");
        for (Integer prime : {2, 3, 13, 179, 997, 31627, 1000003}) {
            for (int offset = -1; offset <= 1; ++offset) {
                const Integer square = prime * prime;
                CheckRange(reference, std::max<Integer>(0, square - 23), square + offset, true);
            }
        }
        const Integer smallCut = Integer(Almost2310L1Cache - 100) * (Almost2310L1Cache - 100);
        CheckRange(reference, smallCut - 1001, smallCut + 1001, true);
        std::mt19937 random(9071);
        for (int i = 0; i < 120; ++i) {
            const Integer low = random() % 100000000;
            CheckRange(reference, low, low + random() % 6000, true);
        }
        Check(PrimeCounting::PiPrimeRange(100, 0) == 0, "reversed range is not empty");
        Check(PrimeCounting::PiPrimeRange(MAX_VALUE, MAX_VALUE - 1) == 0,
              "reversed upper-bound range is not empty");
        CheckOutOfRange([] { PrimeCounting::PiPrime(-1); });
        CheckOutOfRange([] { PrimeCounting::PiPrime(MAX_VALUE + 1); });
        CheckOutOfRange([] { PrimeCounting::PiPrimeRange(-1, 0); });
        CheckOutOfRange([] { PrimeCounting::PiPrimeRange(0, -1); });
        CheckOutOfRange([] { PrimeCounting::PiPrimeRange(MAX_VALUE + 1, 0); });
        CheckOutOfRange([] { PrimeCounting::PiPrimeRange(0, MAX_VALUE + 1); });
    }

    struct Prepared {
        std::vector<int> main, first, second;
        std::vector<char> wheel;
        explicit Prepared(const Reference& reference, Integer high, std::size_t segments)
            : main(reference.base(high)), wheel(NUM30030, 1) {
            const auto split = std::lower_bound(main.begin(), main.end(), segments * NUM30030 / 2);
            first.assign(main.begin(), split);
            second.assign(split, main.end());
            std::size_t index = 0;
            for (auto gap : ARR_WHEEL30030) { wheel[index] = 0; index += gap; }
        }
    };

    void KernelChecks(const Reference& reference) {
        // Force every kernel and both mixed dispatch boundaries. Big's final
        // fixture has both densely marked primes and scheduled bucket primes.
        for (Integer low : {100000000LL, 1000000000000LL}) {
            const Integer high = low + 120137;
            const auto expected = reference.range(low, high);
            const Prepared inputs(reference, high, 1);
            Integer count = 0;
            std::vector<Integer> generated;
            if (low < 1000000000) {
                PrimeSieve::PrimeSieveSmall(inputs.main, count, low, high);
                PrimeSieve::PrimeSieveSmall(inputs.main, generated, low, high);
                Check(count == Integer(expected.size()), "forced Small count mismatch");
                CheckGenerated(generated, expected);
            }
            count = 0;
            generated.clear();
            PrimeSieve::PrimeSieveMedium(inputs.main, count, low, high);
            PrimeSieve::PrimeSieveMedium(inputs.main, generated, low, high);
            Check(count == Integer(expected.size()), "forced Medium count mismatch");
            CheckGenerated(generated, expected);
            count = 0;
            generated.clear();
            PrimeSieve::PrimeSieveBig(inputs.first, inputs.second, inputs.wheel,
                                      count, 1, low, high);
            PrimeSieve::PrimeSieveBig(inputs.first, inputs.second, inputs.wheel,
                                      generated, 1, low, high);
            Check(count == Integer(expected.size()), "forced Big count mismatch");
            CheckGenerated(generated, expected);
        }
        const Integer smallCut = 100000000, mediumCut = 100002000;
        for (Integer low : {smallCut - 101, smallCut, mediumCut - 101, mediumCut}) {
            const Integer high = mediumCut + 107;
            const Prepared inputs(reference, high, 1);
            const auto expected = reference.range(low, high);
            Integer count = 0;
            std::vector<Integer> generated;
            PrimeSieve::PrimeWorker(inputs.main, inputs.first, inputs.second,
                inputs.wheel, count, low, high, smallCut, mediumCut, 1);
            PrimeSieve::PrimeWorker(inputs.main, inputs.first, inputs.second,
                inputs.wheel, generated, low, high, smallCut, mediumCut, 1);
            Check(count == Integer(expected.size()), "mixed kernel count mismatch");
            CheckGenerated(generated, expected);
        }
    }

    void ParallelAndHardwareChecks(const Reference& reference) {
        const Integer low = 100003, high = 2000033;
        const auto expected = reference.range(low, high);
        for (int threads : {1, 2, 3, 8}) {
            Check(PrimeCounting::PiPrimeRange(low, high, threads, 10, 10) == Integer(expected.size()),
                  "parallel Small count mismatch");
            CheckGenerated(Generate<int>(low, high, threads, 10, 10), expected);
            CheckGenerated(Generate<double>(low, high, threads, 10, 10), expected);
            CheckGenerated(Generate<Integer>(low, high, threads, 10, 10), expected);
        }
        for (int threads : {-3, 0, 1, 2, 19}) {
            for (int maximum : {-2, 0, 1, 3}) {
                Check(PrimeCounting::PiPrimeRange(low, high, threads, maximum, 0) == Integer(expected.size()),
                      "count budget normalization mismatch");
            }
        }
        // These values overflowed the original int arithmetic before the
        // result was assigned to a wider type. Only a small budget is used.
        for (int cores : {342, std::numeric_limits<int>::max()}) {
            Check(PrimeCounting::PiPrimeRange(low, high, 2, 1024, cores) == Integer(expected.size()),
                  "large hardware limit count mismatch");
            CheckGenerated(Generate<Integer>(low, high, 2, 1024, cores), expected);
        }
        // Independent native budgets can request far more workers than cores.
        // That formerly inverted the Small/Medium cutoffs, sending Small past
        // the prepared prime-table coverage even for this short interval.
        const int unlimited = std::numeric_limits<int>::max();
        Check(PrimeCounting::PiPrimeRange(2, 50000, unlimited, unlimited, 1) == 5133,
              "independent native budgets inverted the kernel cutoffs");
        const Integer usefulBoundary = 2 + 8 * Almost2310L1Cache;
        for (int offset = -1; offset <= 1; ++offset) {
            const Integer finish = usefulBoundary + offset;
            const auto expectedRange = reference.range(2, finish);
            Check(PrimeCounting::PiPrimeRange(2, finish, unlimited, unlimited, 1) ==
                      Integer(expectedRange.size()), "large request reduction count mismatch");
            std::vector<Integer> slots;
            Integer serialCount = 0;
            bool parallel = false;
            PrimeSieve::SieveMain(slots, serialCount, 2, finish, parallel,
                                   unlimited, unlimited, 1);
            Check(parallel == (offset >= 0), "useful-worker boundary changed");
            Check(slots.size() == (parallel ? 2u : 0u),
                  "count slots allocated before useful-worker reduction");
            Integer actual = serialCount;
            for (Integer count : slots) actual += count;
            Check(actual == Integer(expectedRange.size()), "reduced slot total mismatch");
        }
        const Integer smallCut = Integer(Almost2310L1Cache - 100) * (Almost2310L1Cache - 100);
        Check(PrimeCounting::PiPrimeRange(smallCut - 101, smallCut + 103,
                  unlimited, unlimited, 1) ==
                  Integer(reference.range(smallCut - 101, smallCut + 103).size()),
              "equal Small/Medium cutoffs lost or duplicated primes");
    }

    void LargeChecks() {
        // Independent inclusive interval references from primesieve 12.17:
        // primesieve LOW HIGH --count --quiet --threads=1
        struct Case { Integer low, high, expected; bool parallel; };
        const Case cases[] = {
            {1000000000000LL, 1000010000031LL, 361727, true},
            {1000000000000000LL, 1000000010000031LL, 289394, false},
            {1000000000000000LL, 1000000032000031LL, 926210, true},
            {MAX_VALUE - 1000, MAX_VALUE, 25, false}
        };
        for (const auto& item : cases) {
            Check(PrimeCounting::PiPrimeRange(item.low, item.high) == item.expected,
                  "large serial count mismatch for " + Interval(item.low, item.high));
            Check(PrimeCounting::PiPrimeRange(item.low, item.high, 8, 10, 10) == item.expected,
                  "large parallel count mismatch");
            bool parallel = false;
            const auto generated = Generate<Integer>(item.low, item.high, 8, 10, 10, &parallel);
            Check(Integer(generated.size()) == item.expected, "large generation mismatch");
            Check(parallel == item.parallel, "fixture did not exercise intended parallel path");
            if (!generated.empty()) {
                Check(generated.front() >= item.low && generated.back() <= item.high,
                      "generated prime escaped interval");
            }
        }
        Check(PrimeCounting::PiPrimeRange(9007199254740881LL, 9007199254740881LL) == 1,
              "last supported prime singleton mismatch");
        Check(PrimeCounting::PiPrimeRange(MAX_VALUE, MAX_VALUE) == 0,
              "maximum endpoint singleton mismatch");
        const auto largeDoubles = Generate<double>(MAX_VALUE - 1000, MAX_VALUE);
        Check(largeDoubles.size() == 25 && largeDoubles.back() == 9007199254740881.0,
              "double generation at the exact-integer limit mismatch");
    }
}

int main(int argc, char** argv) {
    try {
        const bool wide = argc == 2 && std::string(argv[1]) == "--wide";
        Check(argc == 1 || wide, "usage: pi-prime-check [--wide]");
        const Reference reference(1100000);
        SmallAndBoundaryChecks(reference);
        std::cout << "PASS small, boundary, domain checks\n" << std::flush;
        KernelChecks(reference);
        ParallelAndHardwareChecks(reference);
        std::cout << "PASS kernels, transitions, budgets, output types\n" << std::flush;
        LargeChecks();
        if (wide) {
            Check(PrimeCounting::PiPrimeRange(1000000000000LL, 1001000000123LL,
                                               8, 10, 10) == 36190994,
                  "multi-window count mismatch");
        }
        std::cout << "PASS assertions=" << assertions << " wide=" << wide << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
