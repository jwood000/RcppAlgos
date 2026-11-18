// Test-only inclusion exposes private production helpers without adding hooks
// or exported interfaces to the package. The runner supplies a frozen copy.
#include "NumbersUtils/ThreadBatch.h"
#include <iostream>
#include <limits>
#include <thread>
#include <utility>

namespace LaunchObservation {
std::vector<std::size_t> capacities;
std::size_t launches = 0;
void reset() { capacities.clear(); launches = 0; }
}

namespace PrimeCounting {
// Keep the real worker launcher and its exception/join behavior. Only record
// caller-thread construction/launches; worker bodies do not touch these logs.
class ObservedThreadBatch {
    ThreadBatch batch;
public:
    explicit ObservedThreadBatch(std::size_t capacity) : batch(capacity) {
        LaunchObservation::capacities.push_back(capacity);
    }
    template <typename Function>
    void launch(Function work) {
        batch.launch(work);
        ++LaunchObservation::launches;
    }
    void finish() { batch.finish(); }
};
}

#define ThreadBatch ObservedThreadBatch
#include "ProductionP2.cpp"
#undef ThreadBatch
#include "P2SelectorCandidate.h"
#include "InitialChunks.h"

namespace {
std::uint64_t queueChecks = 0, selectorChecks = 0, workerChecks = 0;
void check(bool value, std::uint64_t& count, const char* message) {
    ++count;
    if (!value) throw std::runtime_error(message);
}
using Range = std::pair<std::int64_t, std::int64_t>;

std::vector<Range> expected(std::int64_t first, std::int64_t last) {
    std::vector<Range> result;
    const auto head = first + std::min<std::int64_t>(32, last - first);
    while (first < head) {
        result.emplace_back(first, first + 1);
        ++first;
    }
    while (first < last) {
        const auto end = first + std::min<std::int64_t>(32, last - first);
        result.emplace_back(first, end);
        first = end;
    }
    return result;
}

void queueRange(std::int64_t first, std::int64_t last) {
    PrimeCounting::PhiWorkQueue actual(first, last);
    InitialChunks::Queue<1> frozen(first, last);
    const auto reference = expected(first, last);
    check(static_cast<std::size_t>(actual.size()) == reference.size() &&
          static_cast<std::uint64_t>(actual.size()) == frozen.size(),
          queueChecks, "queue size");
    for (std::int64_t ticket = 0; ticket < actual.size(); ++ticket) {
        std::int64_t begin = -1, end = -1, frozenBegin = -1, frozenEnd = -1;
        check(actual.range(ticket, begin, end) && frozen.range(ticket, frozenBegin, frozenEnd),
              queueChecks, "queue range");
        check(Range(begin, end) == reference[ticket] && begin == frozenBegin && end == frozenEnd,
              queueChecks, "queue coverage");
        begin = end = frozenBegin = frozenEnd = -1;
        check(actual.claim(begin, end) && frozen.claim(frozenBegin, frozenEnd),
              queueChecks, "queue claim");
        check(Range(begin, end) == reference[ticket] && begin == frozenBegin && end == frozenEnd,
              queueChecks, "queue claim mapping");
    }
    for (auto ticket : {std::numeric_limits<std::int64_t>::min(), std::int64_t(-1),
                        actual.size(), actual.size() + 1,
                        std::numeric_limits<std::int64_t>::max()}) {
        std::int64_t begin = -1, end = -1;
        check(!actual.range(ticket, begin, end) && begin == -1 && end == -1,
              queueChecks, "queue terminal range");
    }
    std::int64_t begin = -1, end = -1;
    check(!actual.claim(begin, end), queueChecks, "queue terminal claim");
}

void concurrent(std::int64_t first, std::int64_t last, int workers) {
    PrimeCounting::PhiWorkQueue queue(first, last);
    std::vector<std::vector<Range>> results(workers);
    std::vector<std::thread> threads;
    for (int slot = 0; slot < workers; ++slot) threads.emplace_back([&, slot] {
        std::int64_t begin, end;
        while (queue.claim(begin, end)) results[slot].emplace_back(begin, end);
    });
    for (auto& thread : threads) thread.join();
    std::vector<Range> combined;
    for (const auto& result : results) combined.insert(combined.end(), result.begin(), result.end());
    std::sort(combined.begin(), combined.end());
    check(combined == expected(first, last), queueChecks, "concurrent queue coverage");
}

void checkQueue() {
    const auto maximum = std::numeric_limits<std::int64_t>::max();
    for (std::int64_t length = 0; length <= 130; ++length)
        for (auto first : {std::int64_t(0), std::int64_t(1), std::int64_t(590), maximum - 130})
            queueRange(first, first + length);
    for (auto length : {2047, 2048, 4095, 4096, 8191, 8192, 8193}) queueRange(591, 591 + length);
    for (auto length : {0, 1, 2, 31, 32, 33, 63, 64, 65, 2047, 4095, 4096, 8193})
        for (auto workers : {1, 2, 3, 4, 8}) concurrent(maximum - 8193, maximum - 8193 + length, workers);
    for (auto first : {std::int64_t(0), std::int64_t(1), std::int64_t(31)}) {
        PrimeCounting::PhiWorkQueue actual(first, maximum);
        InitialChunks::Queue<1> frozen(first, maximum);
        check(static_cast<std::uint64_t>(actual.size()) == frozen.size(), queueChecks, "extreme queue size");
        for (auto ticket : {std::int64_t(0), std::int64_t(31), std::int64_t(32), actual.size() - 1}) {
            std::int64_t begin = -1, end = -1, frozenBegin = -1, frozenEnd = -1;
            check(actual.range(ticket, begin, end) && frozen.range(ticket, frozenBegin, frozenEnd),
                  queueChecks, "extreme queue range");
            check(begin == frozenBegin && end == frozenEnd && first <= begin && begin < end && end <= maximum,
                  queueChecks, "extreme queue mapping");
            if (ticket == actual.size() - 1) check(end == maximum, queueChecks, "extreme queue last");
        }
    }
    for (auto range : {Range(-1, 0), Range(2, 1), Range(-2, -1)}) {
        bool caught = false;
        try { PrimeCounting::PhiWorkQueue queue(range.first, range.second); }
        catch (const std::out_of_range&) { caught = true; }
        check(caught, queueChecks, "invalid queue interval");
    }
    // Hold the first task unfinished while another thread reaches later jobs.
    PrimeCounting::PhiWorkQueue queue(1, 98);
    std::atomic<bool> held(false), reached(false);
    bool firstOK = false, tailOK = false;
    std::thread first([&] {
        std::int64_t begin, end;
        firstOK = queue.claim(begin, end) && begin == 1;
        held.store(true, std::memory_order_release);
        while (!reached.load(std::memory_order_acquire)) std::this_thread::yield();
    });
    std::thread second([&] {
        while (!held.load(std::memory_order_acquire)) std::this_thread::yield();
        std::int64_t begin, end;
        while (queue.claim(begin, end)) if (begin >= 33) { tailOK = true; break; }
        reached.store(true, std::memory_order_release);
    });
    first.join();
    second.join();
    check(firstOK && tailOK, queueChecks, "queue added a barrier");
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
        constexpr int maximum = 10000001;
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

void checkWorkers(const Oracle& oracle) {
    // Grain boundaries and large arithmetic inputs, without launching large
    // thread counts or allocating tables for those inputs.
    constexpr std::int64_t grain = 1048576;
    for (auto queries : std::vector<std::int64_t>{1, 2, 31, 32, 33, 64, 4096,
                                                std::numeric_limits<std::int64_t>::max()}) {
        for (auto span : std::vector<std::int64_t>{1, grain - 1, grain, 2 * grain - 1, 2 * grain,
                                                 33 * grain, 64 * grain, std::numeric_limits<std::int64_t>::max()}) {
            for (int threads : {-2, 0, 1, 2, 8, 32, 33, 36, 64,
                                std::numeric_limits<int>::max()}) {
                std::int64_t expected = threads > 0 ? threads : 1;
                if (expected > queries) expected = queries;
                const auto bySpan = span < grain ? 1 : span / grain;
                if (expected > bySpan) expected = bySpan;
                check(PrimeCounting::streamWorkerBudget(queries, span, threads) == expected,
                      workerChecks, "stream worker budget");
            }
        }
    }

    const std::int64_t x = 100000000, first = 10;
    const PrimeCounting::PhiTables tables(root(x));
    for (std::int64_t terms : {0, 1, 31, 32, 33, 63, 64, 65}) {
        const auto ranges = expected(first, first + terms);
        std::int64_t reference = 0;
        for (auto i = first; i < first + terms; ++i)
            reference += oracle.pi(x / tables.primes()[i]) - (i - 1);
        for (std::size_t budget : {1, 8, 32, 33, 36, 64}) {
            LaunchObservation::reset();
            const auto actual = PrimeCounting::countPhiPrefix(x, first, first + terms, budget, tables);
            check(actual == reference, workerChecks, "prefix independent count");
            const auto workers = std::min(budget, ranges.size());
            if (workers <= 1) {
                check(LaunchObservation::capacities.empty() && LaunchObservation::launches == 0,
                      workerChecks, "empty/serial prefix launched workers");
            } else {
                check(LaunchObservation::capacities == std::vector<std::size_t>{workers} &&
                      LaunchObservation::launches == workers,
                      workerChecks, "prefix launched more workers than jobs");
            }
        }
    }

    // This public-call geometry really inherits more workers than the
    // shortened prefix has jobs; it is not only an artificial helper budget.
    const std::int64_t largeX = 1000000000000LL, a = 3304, L = 1000000;
    const auto b = oracle.pi(L), begin = a + 1;
    const auto last = std::upper_bound(oracle.primes.begin() + begin,
        oracle.primes.begin() + b + 1, largeX / (L + 1)) - oracle.primes.begin();
    const auto budget = static_cast<std::size_t>(std::min<std::int64_t>(64,
        std::max<std::int64_t>(1, (last - begin) / 2048)));
    const auto end = PrimeCounting::chooseStreamFirst(largeX, L, begin, last, budget, 64, oracle.primes);
    check(end - begin == 64 && budget == 36, workerChecks, "shortened prefix regression geometry");
    const PrimeCounting::PhiTables largeTables(L);
    const auto reference = PrimeCounting::countP2Range(largeX, begin, end, largeTables);
    LaunchObservation::reset();
    check(PrimeCounting::countPhiPrefix(largeX, begin, end, budget, largeTables) == reference,
          workerChecks, "shortened prefix serial agreement");
    check(LaunchObservation::capacities == std::vector<std::size_t>{33} && LaunchObservation::launches == 33,
          workerChecks, "shortened prefix actual launches");
}

void selection(const Oracle& oracle, std::int64_t x, std::int64_t L,
               std::int64_t a, int threads) {
    const auto& primes = oracle.primes;
    const auto first = a + 1, b = oracle.pi(root(x));
    const std::int64_t last = std::upper_bound(primes.begin() + first,
        primes.begin() + b + 1, x / (L + 1)) - primes.begin();
    const auto budget = static_cast<std::size_t>(std::min<std::int64_t>(std::max(1, threads),
        std::max<std::int64_t>(1, (last - first) / 2048)));
    const auto frozen = P2SelectorCandidate::choose(x, L, first, last, budget, threads, primes);
    const auto actual = PrimeCounting::chooseStreamFirst(x, L, first, last, budget, threads, primes);
    check(actual == frozen.streamFirst, selectorChecks, "production selector differs from frozen candidate");
    check(first <= actual && actual <= last, selectorChecks, "selector range");
    if (actual < last) {
        check(x / primes[actual] <= frozen.limit, selectorChecks, "stream first quotient");
        if (actual > first) check(x / primes[actual - 1] > frozen.limit, selectorChecks, "Phi last quotient");
    }
}

void checkSelector() {
    const Oracle oracle;
    checkWorkers(oracle);
    const std::int64_t prime = 503, power = prime * prime * prime * prime;
    const auto before = oracle.pi(prime) - 1, after = before + 1;
    for (int delta : {-1, 0, 1}) {
        const auto x = power + delta, L = root(x);
        for (auto a : {before, after})
            for (auto limit : {L, L + 1}) selection(oracle, x, limit, a, 2);
    }
    selection(oracle, 10000000000LL, 100001, oracle.pi(50000), 2);
    const std::vector<std::int64_t> empty;
    for (std::int64_t last : {0, 1})
        check(PrimeCounting::chooseStreamFirst(1, 1, 0, last, 1, 1, empty) == last,
              selectorChecks, "empty/singleton selector access");
    // All cap transitions, parallel promotion and useful-worker changes.
    for (const auto x : {100000000LL, 10000000000LL, 80000000000LL,
                        99999999999LL, 100000000000LL, 999999999999LL, 1000000000000LL,
                        2999999999999LL, 3000000000000LL, 29999999999999LL,
                        30000000000000LL, 100000000000000LL}) {
        const auto r = root(x), b = oracle.pi(r);
        for (auto L : {r, r + 1, 2 * r}) {
            std::vector<std::int64_t> cutoffs = {0, oracle.pi(root(r)), b / 2, b - 2, b - 1, b};
            for (auto tail : {2047, 2048, 4095, 4096, 4097, 8192})
                cutoffs.push_back(std::max<std::int64_t>(0, b - tail));
            for (auto a : cutoffs)
                for (int threads : {-2, 0, 1, 2, 3, 4, 6, 8, 32, 33, 36, 64, 128})
                    selection(oracle, x, L, a, threads);
        }
    }
}
}

int main() {
    try {
        checkQueue();
        checkSelector();
        std::cout << "PASS queue_checks=" << queueChecks << " selector_checks=" << selectorChecks
                  << " worker_checks=" << workerChecks << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
