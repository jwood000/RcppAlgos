// Replay prior selector contexts using only independent primes and the actual
// frozen policy headers. No production helper, Phi table/cache, or P2 call.
#include "P2SelectorV2.h"
#include "P2SelectorV3.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
    using I = std::int64_t;

    I integerSqrt(I x) {
        I low = 0, high = x;
        while (low < high) {
            const I middle = low + (high - low + 1) / 2;
            if (middle <= x / middle) low = middle;
            else high = middle - 1;
        }
        return low;
    }

    std::vector<I> independentPrimes(I limit) {
        std::vector<unsigned char> prime(static_cast<std::size_t>(limit) + 1, 1);
        prime[0] = prime[1] = 0;
        for (I p = 2; p <= limit / p; ++p) if (prime[p]) {
            for (I n = p * p; n <= limit; n += p) prime[n] = 0;
        }
        std::vector<I> result(1, 0); // Match the policies' one-based indexing.
        for (I n = 2; n <= limit; ++n) if (prime[n]) result.push_back(n);
        return result;
    }

    struct Context {
        std::string id, reason;
        I x, a, limit, b, nonTable, direct, phiBudget, selectedQ;
        I phiQueries, streamQueries, candidateQ, candidateFirst;
        I candidateQueries, candidateSpan, candidateBudget;
        int threads, multiplier;

        bool read(std::istream& input) {
            return static_cast<bool>(input >> id >> x >> a >> limit >> threads
                >> b >> nonTable >> direct >> phiBudget >> selectedQ >> multiplier
                >> reason >> phiQueries >> streamQueries >> candidateQ
                >> candidateFirst >> candidateQueries >> candidateSpan >> candidateBudget);
        }
    };

    void require(bool condition, const Context& c, const char* message) {
        if (!condition) throw std::runtime_error(c.id + ": " + message);
    }

    template <typename A, typename B>
    bool sameCandidate(const A& a, const B& b) {
        return a.candidateLimit == b.candidateLimit &&
            a.candidateStreamFirst == b.candidateStreamFirst &&
            a.candidateStreamQueries == b.candidateStreamQueries &&
            a.candidateSpan == b.candidateSpan &&
            a.candidateStreamBudget == b.candidateStreamBudget;
    }

    template <typename A, typename B>
    bool sameSelection(const A& a, const B& b) {
        return a.limit == b.limit && a.streamFirst == b.streamFirst &&
            a.multiplier == b.multiplier && a.nonTableQueries == b.nonTableQueries &&
            a.phiQueries == b.phiQueries && a.streamQueries == b.streamQueries;
    }
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("usage: replay contexts.tsv");
        std::ifstream input(argv[1]);
        std::size_t expected;
        I maximum;
        if (!(input >> expected >> maximum) || maximum < 2) {
            throw std::runtime_error("invalid context header");
        }
        const auto primes = independentPrimes(maximum);
        std::size_t contexts = 0, bypasses = 0, evaluated = 0;
        std::vector<std::string> changed;
        Context c;
        while (c.read(input)) {
            ++contexts;
            const I root = integerSqrt(c.x);
            require(c.limit >= root && c.limit <= maximum, c, "table coverage");
            const I b = std::upper_bound(primes.begin(), primes.end(), root) -
                primes.begin() - 1;
            require(b == c.b, c, "independent square-root prime count");
            if (c.x < 4 || c.a >= b) {
                require(c.reason == "empty" && c.nonTable == 0 && c.direct == 0 &&
                        c.phiQueries == 0 && c.streamQueries == 0, c, "empty bypass");
                ++bypasses;
                continue;
            }
            const I first = c.a + 1;
            const I last = std::upper_bound(primes.begin() + first,
                primes.begin() + b + 1, c.x / (c.limit + 1)) - primes.begin();
            require(last - first == c.nonTable && b + 1 - last == c.direct,
                    c, "independent table partition");
            if (last == first) {
                require(c.reason == "direct-only" && c.phiQueries == 0 &&
                        c.streamQueries == 0, c, "direct-table bypass");
                ++bypasses;
                continue;
            }
            const I budget = std::min<I>(std::max(1, c.threads),
                std::max<I>(1, (last - first) / 2048));
            require(budget == c.phiBudget, c, "Phi budget");
            const auto v2 = P2SelectorV2::choose(c.x, c.limit, first, last,
                static_cast<std::size_t>(budget), c.threads, primes);
            const auto v3 = P2SelectorV3::choose(c.x, c.limit, first, last,
                static_cast<std::size_t>(budget), c.threads, primes);
            ++evaluated;
            require(v2.limit == c.selectedQ && v2.multiplier == c.multiplier &&
                    P2SelectorV2::reasonName(v2.reason) == c.reason &&
                    v2.nonTableQueries == c.nonTable && v2.phiQueries == c.phiQueries &&
                    v2.streamQueries == c.streamQueries &&
                    v2.streamFirst == first + c.phiQueries, c, "recorded v2 selection");
            require(v2.candidateLimit == c.candidateQ &&
                    v2.candidateStreamFirst == c.candidateFirst &&
                    v2.candidateStreamQueries == c.candidateQueries &&
                    v2.candidateSpan == c.candidateSpan &&
                    v2.candidateStreamBudget == static_cast<std::size_t>(c.candidateBudget),
                    c, "recorded v2 candidate");
            require(sameCandidate(v2, v3), c, "v2/v3 candidate equality");

            // Test the intended scope with a different exact formulation:
            // the product is bounded by 256*pi(10^7) in this replay.
            const bool shouldChange = v2.reason == P2SelectorV2::Reason::ParallelCoverage &&
                v2.candidateSpan > 256 * v2.candidateStreamQueries;
            if (shouldChange) {
                require(v3.limit == c.limit && v3.streamFirst == last &&
                        v3.multiplier == 0 && v3.reason == P2SelectorV3::Reason::StreamTooSparse &&
                        v3.nonTableQueries == c.nonTable && v3.phiQueries == c.nonTable &&
                        v3.streamQueries == 0, c, "expected sparse-stream fallback");
                changed.push_back(c.id);
            } else {
                require(sameSelection(v2, v3) &&
                        std::string(P2SelectorV2::reasonName(v2.reason)) ==
                            P2SelectorV3::reasonName(v3.reason),
                        c, "unexpected selection change");
            }
        }
        if (!input.eof() || contexts != expected) throw std::runtime_error("incomplete input");
        const std::vector<std::string> expectedChanges = {
            "p2-meissel-10000000000-t4-L2", "p2-meissel-70000000000-t4-L2",
            "p2-zero-10000000000-t4-L2", "p2-zero-70000000000-t4-L2"};
        std::sort(changed.begin(), changed.end());
        if (changed != expectedChanges) throw std::runtime_error("unexpected changed context IDs");
        std::cout << "{\"contexts\":" << contexts << ",\"maximum_prime_limit\":" << maximum
            << ",\"independent_primes\":" << primes.size() - 1
            << ",\"kernel_bypasses\":" << bypasses << ",\"evaluated_policy_pairs\":" << evaluated
            << ",\"candidate_equality_pairs\":" << evaluated
            << ",\"unchanged_contexts\":" << contexts - changed.size()
            << ",\"changed_contexts\":" << changed.size() << ",\"changed_ids\":[";
        for (std::size_t i = 0; i < changed.size(); ++i) {
            if (i) std::cout << ',';
            std::cout << '"' << changed[i] << '"';
        }
        std::cout << "],\"phi_or_p2_calls\":0,\"passed\":true}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
