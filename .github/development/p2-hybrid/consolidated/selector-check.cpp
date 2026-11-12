// Focused checks of the consolidated policy's new admission boundaries.
// Independent prime sieve and integer roots; no Phi calls or benchmark timers.
#include "P2SelectorCandidate.h"
#include "P2SelectorV2.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
std::uint64_t checks = 0;
void check(bool condition) {
    ++checks;
    if (!condition) throw std::runtime_error("consolidated selector check failed");
}
std::int64_t root(std::int64_t x) {
    std::int64_t low = 0, high = x;
    while (low < high) {
        const auto middle = low + (high - low + 1) / 2;
        if (middle <= x / middle) low = middle;
        else high = middle - 1;
    }
    return low;
}
struct Oracle {
    std::vector<std::int64_t> primes{0};
    Oracle() {
        const int maximum = 2000001;
        std::vector<bool> sieve(maximum + 1, true);
        sieve[0] = sieve[1] = false;
        for (int p = 2; p <= maximum / p; ++p) if (sieve[p])
            for (int n = p * p; n <= maximum; n += p) sieve[n] = false;
        for (int n = 2; n <= maximum; ++n) if (sieve[n]) primes.push_back(n);
    }
    std::int64_t pi(std::int64_t x) const {
        return std::upper_bound(primes.begin(), primes.end(), x) - primes.begin() - 1;
    }
};

P2SelectorCandidate::Decision choose(const Oracle& oracle, std::int64_t x,
        std::int64_t L, std::int64_t a, int requested) {
    const auto& primes = oracle.primes;
    const auto first = a + 1, b = oracle.pi(root(x));
    const std::int64_t last = std::upper_bound(primes.begin() + first,
        primes.begin() + b + 1, x / (L + 1)) - primes.begin();
    const auto budget = static_cast<std::size_t>(std::min<std::int64_t>(
        std::max(1, requested), std::max<std::int64_t>(1, (last - first) / 2048)));
    const auto result = P2SelectorCandidate::choose(x, L, first, last,
                                                  budget, requested, primes);
    const auto v2 = P2SelectorV2::choose(x, L, first, last, budget, requested, primes);
    check(result.nonTableQueries == last - first);
    check(result.phiQueries + result.streamQueries == result.nonTableQueries);
    check(result.streamFirst >= first && result.streamFirst <= last);
    check(result.candidateLimit == v2.candidateLimit);
    check(result.candidateStreamFirst == v2.candidateStreamFirst);
    check(result.candidateSpan == v2.candidateSpan);
    check(result.candidateStreamBudget == v2.candidateStreamBudget);
    if (result.multiplier == 0) {
        check(result.limit == L && result.streamFirst == last && result.streamQueries == 0);
    } else {
        check(result.limit == v2.limit && result.streamFirst == v2.streamFirst);
        check(result.multiplier == v2.multiplier);
        for (auto i = first; i < last; ++i)
            check((x / primes[i] <= result.limit) == (i >= result.streamFirst));
    }
    const bool addedGate = result.reason == P2SelectorCandidate::Reason::EnlargedTable ||
                          result.reason == P2SelectorCandidate::Reason::EarlyPrefix;
    if (addedGate) {
        check(x < 100000000000LL && budget > 1 && v2.multiplier > 0);
        check(result.reason == P2SelectorCandidate::Reason::EnlargedTable ?
              L != root(x) : (L == root(x) && primes[first] <= root(root(x))));
    } else {
        check(result.limit == v2.limit && result.streamFirst == v2.streamFirst);
        check(result.multiplier == v2.multiplier);
    }
    return result;
}
}

int main() {
    try {
        const Oracle oracle;
        using Reason = P2SelectorCandidate::Reason;
        // Exact fourth-power boundary: 503 itself is a prime. At 503^4-1,
        // first=503 is above the fourth root; at/above 503^4 it is not.
        const std::int64_t p = 503, power = p * p * p * p;
        const auto at = oracle.pi(p) - 1, after = at + 1;
        for (int delta : {-1, 0, 1}) {
            const auto x = power + delta, L = root(x);
            check(choose(oracle, x, L, at, 2).reason ==
                  (delta < 0 ? Reason::ParallelCoverage : Reason::EarlyPrefix));
            check(choose(oracle, x, L, after, 2).reason == Reason::ParallelCoverage);
            check(choose(oracle, x, L + 1, after, 2).reason == Reason::EnlargedTable);
        }
        // The short-whole-range exception precedes the enlarged-table guard,
        // even with enough terms for a genuine two-worker Phi budget.
        const auto shortFull = choose(oracle, 10000000000LL, 100001, oracle.pi(50000), 2);
        check(shortFull.reason == Reason::ShortFullRange);
        check(shortFull.nonTableQueries >= 4096 && shortFull.candidateSpan <= 131072);
        // Empty and singleton ranges must never inspect a nonexistent prime.
        const std::vector<std::int64_t> empty;
        check(P2SelectorCandidate::choose(1, 1, 0, 0, 1, 1, empty).reason == Reason::AtMostOneQuery);
        check(P2SelectorCandidate::choose(1, 1, 0, 1, 1, 1, empty).reason == Reason::AtMostOneQuery);
        // Existing caps and guards are preserved outside the new two gates.
        for (const auto x : {100000000LL, 10000000000LL, 80000000000LL,
                            100000000000LL, 1000000000000LL}) {
            const auto r = root(x), b = oracle.pi(r);
            for (const auto L : {r, r + 1, 2 * r})
                for (const auto a : {std::int64_t(0), oracle.pi(root(r)), b / 2, b - 2, b - 1, b})
                    for (int threads : {-2, 0, 1, 2, 4, 8}) choose(oracle, x, L, a, threads);
        }
        std::cout << "PASS selector_checks=" << checks << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
