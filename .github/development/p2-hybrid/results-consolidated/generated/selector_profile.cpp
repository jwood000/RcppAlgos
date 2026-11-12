// Generated consolidated comparison; diagnostic decisions remain outside timers.
// Standalone selector v3 comparison experiment. All diagnostic decisions and capacity reads
// occur AFTER the timed kernels. No production file includes this harness.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/Phi.h"
#include "selector_v3_cases.h"
#include "P2Selector.h"
#ifdef P2_HYBRID_SELECTOR_V2
#include "P2SelectorV2.h"
#endif
#ifdef P2_HYBRID_SELECTOR_V3
#include "P2SelectorV3.h"
#endif
#ifdef P2_HYBRID_SELECTOR_CANDIDATE
#include "P2SelectorCandidate.h"
#endif
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef P2_HYBRID_MEMORY
#include "hybrid_metrics.h"
#endif

#if (defined(P2_SELECTOR_BASELINE) + defined(P2_HYBRID_SELECTOR) + defined(P2_HYBRID_SELECTOR_V2) + defined(P2_HYBRID_SELECTOR_V3) + defined(P2_HYBRID_SELECTOR_CANDIDATE) + defined(P2_HYBRID_MULTIPLIER)) != 1
#error "Select exactly one baseline or hybrid policy build"
#endif

namespace {
    using Clock = std::chrono::steady_clock;
    using PrimeCounting::PhiTables;

    double Milliseconds(Clock::time_point first, Clock::time_point last) {
        return std::chrono::duration<double, std::milli>(last - first).count();
    }

    std::int64_t Integer(const char* text) {
        std::size_t used = 0;
        const std::string value(text);
        const auto result = std::stoll(value, &used);
        if (used != value.size()) throw std::runtime_error("invalid integer: " + value);
        return result;
    }

    const char* BuildKind() {
#if defined(P2_SELECTOR_BASELINE) && defined(P2_INITIAL_CHUNK)
        return "phi-initial-1";
#elif defined(P2_SELECTOR_BASELINE)
        return "baseline";
#elif defined(P2_HYBRID_SELECTOR_CANDIDATE)
        return "candidate";
#elif defined(P2_HYBRID_SELECTOR_V3)
        return "selector-v3";
#elif defined(P2_HYBRID_SELECTOR_V2)
        return "selector-v2";
#elif defined(P2_HYBRID_SELECTOR)
        return "selector-v1";
#else
        return "fixed";
#endif
    }

    struct Options {
        std::string mode, method;
        std::int64_t x, repeats, expected, limit;
        int threads, tableNumerator, tableDenominator;

        Options(int argc, char** argv) {
            if (argc < 7 || argc > 9) throw std::runtime_error(
                "usage: selector-profile phases|p2 cutoff x threads repeats expected [table_numerator [table_denominator]]\n"
                "       selector-profile reference cutoff x threads 1 -1 [table_numerator [table_denominator]]\n"
                "       selector-profile oracle cutoff x 1 1 -1 [table_numerator [table_denominator]]");
            mode = argv[1];
            method = argv[2];
            x = Integer(argv[3]);
            const auto requested = Integer(argv[4]);
            repeats = Integer(argv[5]);
            expected = Integer(argv[6]);
            const auto numerator = argc >= 8 ? Integer(argv[7]) : 1;
            const auto denominator = argc == 9 ? Integer(argv[8]) : 1;
            const bool timed = mode == "phases" || mode == "p2";
            if ((!timed && mode != "reference" && mode != "oracle") ||
                !P2SelectorV3Cases::Valid(method) ||
                x < 2 || x > PrimeCounting::PartialSieveDetail::MAX_X ||
                requested < std::numeric_limits<int>::min() ||
                requested > std::numeric_limits<int>::max() ||
                repeats < 1 || repeats > std::numeric_limits<int>::max() ||
                denominator < 1 || numerator < denominator ||
                numerator > std::numeric_limits<int>::max() ||
                (timed && (expected < 0 || expected > x)) ||
                (!timed && (expected != -1 || repeats != 1)) ||
                (mode == "phases" && (numerator != denominator ||
                    (method != "meissel" && method != "lehmer"))) ||
                (mode == "oracle" && (x > 100000000 || requested != 1))) {
                throw std::runtime_error("invalid selector profile arguments");
            }
#ifndef P2_SELECTOR_BASELINE
            if (mode == "reference") throw std::runtime_error(
                "reference mode requires a production baseline build");
#endif
            threads = static_cast<int>(requested);
            tableNumerator = static_cast<int>(numerator);
            tableDenominator = static_cast<int>(denominator);
            const auto root = PrimeCounting::PartialSieveDetail::IntegerSqrt(x);
            if (root > std::numeric_limits<std::int64_t>::max() / numerator) {
                throw std::runtime_error("table ratio product exceeds int64_t");
            }
            // Divide after multiplying: truncating the ratio or root first
            // would change floor(root * numerator / denominator).
            limit = root * numerator / denominator;
            if (limit > PrimeCounting::PartialSieveDetail::IntegerSqrt(
                    PrimeCounting::PartialSieveDetail::MAX_X)) {
                throw std::runtime_error("table ratio exceeds the prepared-table limit");
            }
        }
    };

    struct TableInfo {
        std::int64_t a = 0, limit = 0, primeCount = 0;
        P2SelectorV3Cases::Values cutoffs{0, 0, 0, 0};
        std::uint64_t piSize = 0, piCapacity = 0, primesSize = 0, primesCapacity = 0;
        bool hasPayload = false;

        TableInfo() = default;
        TableInfo(const Options& options, std::int64_t cutoff, const PhiTables& tables)
            : a(cutoff), limit(tables.limit()), primeCount(tables.primeCount()),
              cutoffs(P2SelectorV3Cases::Describe(options.x, tables)),
              piSize(tables.pi().size()), piCapacity(tables.pi().capacity()),
              primesSize(tables.primes().size()), primesCapacity(tables.primes().capacity()),
              hasPayload(true) {}
    };

    struct Trace {
        std::int64_t limit = 0, largest = 0, nonTable = 0;
        std::int64_t phiQueries = 0, streamQueries = 0, directQueries = 0;
        std::int64_t streamSpan = 0;
        std::uint64_t phiBudget = 0, phiWorkers = 0;
        std::uint64_t streamBudget = 0, streamWorkers = 0, phiSeeds = 0;
        bool hasCandidate = false;
        std::int64_t candidateLimit = 0, candidateStreamFirst = 0;
        std::int64_t candidateStreamQueries = 0, candidateSpan = 0;
        std::size_t candidateStreamBudget = 0;
        int multiplier = 0;
        const char* reason = "empty";
    };

    // Mirrors the existing scheduling bounds without evaluating Phi, launching
    // workers, or modifying cache state. It is never called inside a timer.
    Trace DescribeChoice(const Options& options, std::int64_t a,
                         const PhiTables& tables) {
        Trace trace;
        trace.limit = tables.limit();
        const auto root = PrimeCounting::PartialSieveDetail::IntegerSqrt(options.x);
        const std::int64_t b = tables.pi()[root];
        if (options.x < 4 || a >= b) return trace;
        const auto& primes = tables.primes();
        const std::int64_t first = a + 1;
        const std::int64_t last = std::upper_bound(
            primes.begin() + first, primes.begin() + b + 1,
            options.x / (tables.limit() + 1)) - primes.begin();
        trace.largest = options.x / primes[first];
        trace.directQueries = b + 1 - last;
        trace.nonTable = last - first;
        if (first == last) {
            trace.reason = "direct-only";
            return trace;
        }
        trace.phiBudget = static_cast<std::uint64_t>(std::min<std::int64_t>(
            std::max(1, options.threads), std::max<std::int64_t>(1, (last - first) / 2048)));
        std::int64_t streamFirst = last;
#ifdef P2_HYBRID_SELECTOR_CANDIDATE
        const auto choice = P2SelectorCandidate::choose(options.x, tables.limit(), first, last,
            trace.phiBudget, options.threads, primes);
        trace.limit = choice.limit;
        trace.multiplier = choice.multiplier;
        trace.reason = P2SelectorCandidate::reasonName(choice.reason);
        streamFirst = choice.streamFirst;
        trace.hasCandidate = true;
        trace.candidateLimit = choice.candidateLimit;
        trace.candidateStreamFirst = choice.candidateStreamFirst;
        trace.candidateStreamQueries = choice.candidateStreamQueries;
        trace.candidateSpan = choice.candidateSpan;
        trace.candidateStreamBudget = choice.candidateStreamBudget;
#elif defined(P2_HYBRID_SELECTOR_V3)
        const auto choice = P2SelectorV3::choose(options.x, tables.limit(), first, last,
            trace.phiBudget, options.threads, primes);
        trace.limit = choice.limit;
        trace.multiplier = choice.multiplier;
        trace.reason = P2SelectorV3::reasonName(choice.reason);
        streamFirst = choice.streamFirst;
        trace.hasCandidate = true;
        trace.candidateLimit = choice.candidateLimit;
        trace.candidateStreamFirst = choice.candidateStreamFirst;
        trace.candidateStreamQueries = choice.candidateStreamQueries;
        trace.candidateSpan = choice.candidateSpan;
        trace.candidateStreamBudget = choice.candidateStreamBudget;
#elif defined(P2_HYBRID_SELECTOR_V2)
        const auto choice = P2SelectorV2::choose(options.x, tables.limit(), first, last,
            trace.phiBudget, options.threads, primes);
        trace.limit = choice.limit;
        trace.multiplier = choice.multiplier;
        trace.reason = P2SelectorV2::reasonName(choice.reason);
        streamFirst = choice.streamFirst;
        trace.hasCandidate = true;
        trace.candidateLimit = choice.candidateLimit;
        trace.candidateStreamFirst = choice.candidateStreamFirst;
        trace.candidateStreamQueries = choice.candidateStreamQueries;
        trace.candidateSpan = choice.candidateSpan;
        trace.candidateStreamBudget = choice.candidateStreamBudget;
#elif defined(P2_HYBRID_SELECTOR)
        const auto choice = P2Selector::choose(options.x, tables.limit(), first, last,
                                               trace.phiBudget, primes);
        trace.limit = choice.limit;
        trace.multiplier = choice.multiplier;
        trace.reason = P2Selector::reasonName(choice.reason);
        streamFirst = choice.streamFirst;
#elif defined(P2_HYBRID_MULTIPLIER)
        trace.multiplier = P2_HYBRID_MULTIPLIER;
        trace.limit = tables.limit() > trace.largest / P2_HYBRID_MULTIPLIER ?
            trace.largest : tables.limit() * P2_HYBRID_MULTIPLIER;
        streamFirst = std::upper_bound(primes.begin() + first, primes.begin() + last,
            options.x / (trace.limit + 1)) - primes.begin();
        trace.reason = "fixed-multiplier";
#else
        trace.reason = "production-phi";
#endif
        trace.phiQueries = streamFirst - first;
        trace.streamQueries = last - streamFirst;
        trace.phiWorkers = std::min<std::uint64_t>(trace.phiBudget, trace.phiQueries);
        if (streamFirst == last) return trace;
        const std::int64_t span = options.x / primes[streamFirst] - tables.limit();
        trace.streamSpan = span;
        const auto workers = std::min<std::int64_t>(std::max(1, options.threads),
            std::min(trace.streamQueries, std::max<std::int64_t>(1, span / 1048576)));
        trace.streamBudget = static_cast<std::uint64_t>(workers);
        if (workers == 1) {
            trace.streamWorkers = 1;
            return trace;
        }
        const auto step = span / workers, extra = span % workers;
        std::int64_t lower = tables.limit(), remaining = last;
        for (std::int64_t slot = 0; slot < workers && remaining > streamFirst; ++slot) {
            const auto upper = lower + step + (slot < extra);
            const std::int64_t begin = std::upper_bound(
                primes.begin() + streamFirst, primes.begin() + remaining,
                options.x / (upper + 1)) - primes.begin();
            if (begin < remaining) {
                ++trace.streamWorkers;
                trace.phiSeeds += lower > tables.limit();
            }
            remaining = begin;
            lower = upper;
        }
        return trace;
    }

    struct Storage {
        bool present = false;
        std::uint64_t sieve = 0, next = 0, phiSums = 0, streamSums = 0, scratch = 0;
    };

    Storage DescribeStorage(const Trace& trace) {
        Storage result;
#if defined(P2_HYBRID_MEMORY) && !defined(P2_SELECTOR_BASELINE)
        const auto& snapshot = P2HybridMetrics::current();
        std::uint64_t workers = 0, phiSeeds = 0;
        std::int64_t queries = 0;
        for (const auto& interval : snapshot.intervals) if (interval.active) {
            ++workers;
            phiSeeds += interval.phiSeed;
            queries += interval.queries;
            result.sieve += interval.sieveBytes;
            result.next += interval.nextBytes;
        }
        if (snapshot.phiQueries != trace.phiQueries ||
            snapshot.streamQueries != trace.streamQueries ||
            snapshot.directQueries != trace.directQueries ||
            snapshot.phiWorkers != trace.phiWorkers ||
            snapshot.streamWorkers != trace.streamWorkers ||
            workers != trace.streamWorkers || phiSeeds != trace.phiSeeds ||
            queries != trace.streamQueries ||
            (trace.nonTable && snapshot.streamLimit != trace.limit)) {
            throw std::runtime_error("selection trace disagrees with instrumented kernel");
        }
        result.phiSums = snapshot.phiSumBytes;
        result.streamSums = snapshot.streamSumBytes;
        result.scratch = std::max(result.phiSums, result.streamSums + result.sieve + result.next);
        result.present = true;
#else
        (void) trace;
#endif
        return result;
    }

    struct Row {
        std::int64_t repeat = 0, phi = 0, p2 = 0, p3 = 0, primeCount = 0;
        std::int64_t oracleLimit = 0;
        bool timed = false, hasPrimeCount = false, hasTrace = false;
        double setup = 0, phiMs = 0, p2Ms = 0, p3Ms = 0, teardown = 0;
        TableInfo table;
        Trace trace;
        Storage storage;
    };

    void Header() {
        std::cout << "mode,method,x,threads,repeat,iterations,cache_state,priming,"
            "setup_ms,phi_ms,p2_ms,p3_ms,warmup_ms,teardown_ms,total_ms,"
            "phi,p2,p3,prime_count,a,table_limit,table_prime_count,"
            "pi_size,pi_capacity,pi_size_bytes,pi_capacity_bytes,"
            "primes_size,primes_capacity,primes_size_bytes,primes_capacity_bytes,"
            "table_size_bytes,table_capacity_bytes,"
            "table_multiplier,root,b,meissel_a,lehmer_a,zero_a,one_a,half_meissel_a,"
            "double_meissel_a,half_root_a,near_root_a,last_a,empty_a,tail4095_a,tail4096_a,cutoff_definition,"
            "build_kind,build_multiplier,measurement,verification,expected,result,timing_definition,"
            "selected_q,selected_k,selection_reason,largest_quotient,non_table_queries,"
            "phi_queries,stream_queries,direct_queries,phi_worker_budget,phi_workers,"
            "stream_worker_budget,stream_workers,seed_queries,phi_seed_queries,"
            "sieve_capacity_bytes,next_capacity_bytes,phi_sum_bytes,stream_sum_bytes,"
            "conservative_scratch_bytes,oracle_limit,stream_span,integers_per_stream_query,"
            "candidate_q,candidate_stream_first,candidate_stream_queries,candidate_stream_span,"
            "candidate_stream_worker_budget,quarter_root_a,three_quarter_root_a,"
            "tail2_a,tail4_a,tail8_a,tail16_a,tail31_a,tail64_a,tail256_a,tail1024_a,"
            "tail2048_a,tail4094_a,tail4097_a,tail8192_a,table_numerator,table_denominator\n";
    }

    void Print(const Options& options, const Row& row) {
        const auto& table = row.table;
        const bool reference = options.mode == "reference";
        const bool oracle = options.mode == "oracle";
        const char* measurement = reference ? "reference-only" : oracle ? "independent-oracle" :
#ifdef P2_HYBRID_MEMORY
            "memory-probe";
#else
            "timing";
#endif
        std::cout << std::setprecision(12) << options.mode << ',' << options.method << ','
            << options.x << ',' << options.threads << ',' << row.repeat << ",1,"
            << (row.repeat ? "warm" : "first") << ",plain,";
        if (row.timed) {
            std::cout << row.setup << ',' << row.phiMs << ',' << row.p2Ms << ',' << row.p3Ms
                << ",0," << row.teardown << ','
                << row.setup + row.phiMs + row.p2Ms + row.p3Ms + row.teardown;
        } else std::cout << ",,,,,,";
        std::cout << ',';
        if (row.hasPrimeCount) std::cout << row.phi;
        std::cout << ',' << row.p2 << ',';
        if (row.hasPrimeCount) std::cout << row.p3;
        std::cout << ',';
        if (row.hasPrimeCount) std::cout << row.primeCount;
        std::cout << ',' << table.a << ',' << table.limit << ',' << table.primeCount << ',';
        if (table.hasPayload) {
            const auto pb = table.piSize * sizeof(std::int32_t);
            const auto pcb = table.piCapacity * sizeof(std::int32_t);
            const auto qb = table.primesSize * sizeof(std::int64_t);
            const auto qcb = table.primesCapacity * sizeof(std::int64_t);
            std::cout << table.piSize << ',' << table.piCapacity << ',' << pb << ',' << pcb
                << ',' << table.primesSize << ',' << table.primesCapacity << ',' << qb << ','
                << qcb << ',' << pb + qb << ',' << pcb + qcb;
        } else std::cout << ",,,,,,,,,";
        const auto& c = table.cutoffs;
        // Keep the old table_multiplier field for denominator-one callers.
        // New consumers use the explicit numerator/denominator pair below.
        std::cout << ',' << options.tableNumerator << ',' << c.root << ',' << c.b << ','
            << c.meissel << ',' << c.lehmer << ",0,1," << c.get("half-meissel") << ','
            << c.get("double-meissel") << ',' << c.get("half-root") << ','
            << c.get("near-root") << ',' << c.get("last") << ',' << c.b << ','
            << c.get("tail4095") << ',' << c.get("tail4096") << ','
            << P2SelectorV3Cases::Definition(options.method) << ',' << BuildKind() << ',';
#ifdef P2_HYBRID_MULTIPLIER
        std::cout << P2_HYBRID_MULTIPLIER;
#else
        std::cout << 0;
#endif
        std::cout << ',' << measurement << ','
            << (reference ? "unchecked-production-reference" : oracle ? "independent-dense-pi" : "expected-match")
            << ',' << options.expected << ',' << (row.hasPrimeCount ? row.primeCount : row.p2)
            << ",phase-sum-excludes-diagnostics,";
        if (row.hasTrace) {
            const auto& t = row.trace;
            std::cout << t.limit << ',' << t.multiplier << ',' << t.reason << ',' << t.largest
                << ',' << t.nonTable << ',' << t.phiQueries << ',' << t.streamQueries << ','
                << t.directQueries << ',' << t.phiBudget << ',' << t.phiWorkers << ','
                << t.streamBudget << ',' << t.streamWorkers << ',' << t.streamWorkers << ',' << t.phiSeeds;
        } else std::cout << ",,,,,,,,,,,,,";
        std::cout << ',';
        if (row.storage.present) {
            const auto& s = row.storage;
            std::cout << s.sieve << ',' << s.next << ',' << s.phiSums << ',' << s.streamSums << ',' << s.scratch;
        } else std::cout << ",,,,";
        std::cout << ',';
        if (oracle) std::cout << row.oracleLimit;
        std::cout << ',';
        if (row.hasTrace) std::cout << row.trace.streamSpan;
        std::cout << ',';
        if (row.hasTrace && row.trace.streamQueries) {
            std::cout << static_cast<double>(row.trace.streamSpan) / row.trace.streamQueries;
        }
        std::cout << ',';
        if (row.trace.hasCandidate) {
            const auto& t = row.trace;
            std::cout << t.candidateLimit << ',' << t.candidateStreamFirst << ','
                << t.candidateStreamQueries << ',' << t.candidateSpan << ',' << t.candidateStreamBudget;
        } else std::cout << ",,,,";
        std::cout << ',' << c.get("quarter-root") << ',' << c.get("three-quarter-root");
        for (const char* label : {"tail2", "tail4", "tail8", "tail16", "tail31", "tail64",
                                  "tail256", "tail1024", "tail2048", "tail4094", "tail4097", "tail8192"}) {
            std::cout << ',' << c.get(label);
        }
        std::cout << ',' << options.tableNumerator << ',' << options.tableDenominator << '\n';
    }

    void Check(std::int64_t actual, std::int64_t expected) {
        if (actual != expected) throw std::runtime_error("count mismatch: actual=" +
            std::to_string(actual) + " expected=" + std::to_string(expected));
    }

    void Describe(Row& row, const Options& options, std::int64_t a, const PhiTables& tables) {
        row.table = TableInfo(options, a, tables);
        row.trace = DescribeChoice(options, a, tables);
        row.hasTrace = true;
        row.storage = DescribeStorage(row.trace);
    }

    void Phases(const Options& options) {
        for (std::int64_t repeat = 0; repeat < options.repeats; ++repeat) {
            Row row;
            row.repeat = repeat;
#ifndef P2_HYBRID_MEMORY
            row.timed = true;
#endif
            const auto setupStart = Clock::now();
            auto tables = std::unique_ptr<PhiTables>(new PhiTables(options.limit));
            const auto a = P2SelectorV3Cases::Cutoff(options.method, options.x, *tables);
            row.setup = Milliseconds(setupStart, Clock::now());
            const auto phiStart = Clock::now();
            row.phi = PrimeCounting::PhiParallel(options.x, a, options.threads, *tables);
            row.phiMs = Milliseconds(phiStart, Clock::now());
            const auto p2Start = Clock::now();
            row.p2 = PrimeCounting::P2(options.x, a, options.threads, *tables);
            row.p2Ms = Milliseconds(p2Start, Clock::now());
            if (options.method == "lehmer") {
                const auto p3Start = Clock::now();
                row.p3 = PrimeCounting::P3(options.x, a, options.threads, *tables);
                row.p3Ms = Milliseconds(p3Start, Clock::now());
            }
            row.primeCount = row.phi + a - 1 - row.p2 - row.p3;
            row.hasPrimeCount = true;
            Describe(row, options, a, *tables);
            const auto teardownStart = Clock::now();
            tables.reset();
            row.teardown = Milliseconds(teardownStart, Clock::now());
            Check(row.primeCount, options.expected);
            Print(options, row);
        }
    }

    void Isolated(const Options& options) {
        const auto setupStart = Clock::now();
        auto tables = std::unique_ptr<PhiTables>(new PhiTables(options.limit));
        const auto a = P2SelectorV3Cases::Cutoff(options.method, options.x, *tables);
        const auto setupMs = Milliseconds(setupStart, Clock::now());
        for (std::int64_t repeat = 0; repeat < options.repeats; ++repeat) {
            Row row;
            row.repeat = repeat;
#ifndef P2_HYBRID_MEMORY
            row.timed = options.mode == "p2";
#endif
            row.setup = repeat ? 0 : setupMs;
            const auto p2Start = Clock::now();
            row.p2 = PrimeCounting::P2(options.x, a, options.threads, *tables);
            row.p2Ms = Milliseconds(p2Start, Clock::now());
            Describe(row, options, a, *tables);
            if (repeat + 1 == options.repeats) {
                const auto teardownStart = Clock::now();
                tables.reset();
                row.teardown = Milliseconds(teardownStart, Clock::now());
            }
            if (options.mode == "p2") Check(row.p2, options.expected);
            Print(options, row);
        }
    }

    // Independent exact roots for the bounded oracle; no native root helper.
    std::int64_t OracleRoot(std::int64_t x, int power) {
        std::int64_t low = 0, high = x;
        while (low < high) {
            const auto middle = low + (high - low + 1) / 2;
            auto quotient = x;
            for (int i = 1; i < power; ++i) quotient /= middle;
            if (middle <= quotient) low = middle;
            else high = middle - 1;
        }
        return low;
    }

    void Oracle(const Options& options) {
        const auto root = OracleRoot(options.x, 2);
        const auto limit = root * options.tableNumerator / options.tableDenominator;
        const auto maximum = std::max(options.x / 2, limit);
        std::vector<unsigned char> prime(static_cast<std::size_t>(maximum) + 1, 1);
        prime[0] = prime[1] = 0;
        for (std::int64_t p = 2; p <= maximum / p; ++p) if (prime[p]) {
            for (auto multiple = p * p; multiple <= maximum; multiple += p) prime[multiple] = 0;
        }
        std::vector<std::int32_t> pi(static_cast<std::size_t>(maximum) + 1, 0);
        std::vector<std::int64_t> primes;
        for (std::int64_t n = 1; n <= maximum; ++n) {
            pi[n] = pi[n - 1] + prime[n];
            if (n <= root && prime[n]) primes.push_back(n);
        }
        const P2SelectorV3Cases::Values values(root, pi[root], pi[OracleRoot(options.x, 3)],
                                              pi[OracleRoot(options.x, 4)]);
        const auto a = values.get(options.method);
        if (a > pi[limit]) throw std::runtime_error("cutoff exceeds the prepared prime count");
        std::int64_t count = 0;
        for (std::int64_t index = a; index < static_cast<std::int64_t>(primes.size()); ++index) {
            count += pi[options.x / primes[index]] - index;
        }
        Row row;
        row.p2 = count;
        row.table.a = a;
        row.table.limit = limit;
        row.table.primeCount = pi[limit];
        row.table.cutoffs = values;
        row.oracleLimit = maximum;
        Print(options, row);
    }
}

int main(int argc, char** argv) {
    try {
        const Options options(argc, argv);
        Header();
        if (options.mode == "phases") Phases(options);
        else if (options.mode == "oracle") Oracle(options);
        else Isolated(options);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
