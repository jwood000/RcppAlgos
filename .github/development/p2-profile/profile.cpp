// Standalone phase/P2 profiler. Link production P2.cpp, P3.cpp, PhiParallel.cpp,
// Phi.cpp and Eratosthenes.cpp; use -std=c++11 -pthread -Iinst/include.
// See README.md for coordinated runs and the optional diagnostic build.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/Phi.h"
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#ifdef P2_PROFILE_METRICS
#include "profile_metrics.h"
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

    struct Options {
        std::string mode, method, priming;
        std::int64_t x, amount, expected;
        int threads;

        Options(int argc, char** argv) {
            if (argc != 7 && argc != 8) throw std::runtime_error(
                "usage: profile phases meissel|lehmer x threads repeats expected_pi\n"
                "       profile p2 meissel|lehmer x threads repeats expected_p2 [plain|prime-phi]\n"
                "       profile sample meissel|lehmer x threads duration_ms expected_p2 [plain|prime-phi]"
            );
            mode = argv[1];
            method = argv[2];
            x = Integer(argv[3]);
            const auto requested = Integer(argv[4]);
            amount = Integer(argv[5]);
            expected = Integer(argv[6]);
            priming = argc == 8 ? argv[7] : "plain";
            if ((mode != "phases" && mode != "p2" && mode != "sample") ||
                (method != "meissel" && method != "lehmer") ||
                (priming != "plain" && priming != "prime-phi") ||
                (mode == "phases" && argc != 7) ||
                x < 2 || x > PrimeCounting::PartialSieveDetail::MAX_X ||
                requested < 1 || requested > std::numeric_limits<int>::max() ||
                amount < 1 || amount > std::numeric_limits<int>::max() ||
                expected < 0 || expected > x) {
                throw std::runtime_error("invalid profile arguments");
            }
            threads = static_cast<int>(requested);
        }
    };

    struct TableInfo {
        std::int64_t a = 0, limit = 0, primeCount = 0;
        std::uint64_t piSize = 0, piCapacity = 0, primesSize = 0, primesCapacity = 0;

        TableInfo() = default;
        TableInfo(const Options& options, const PhiTables& tables) {
            limit = tables.limit();
            const auto cutoff = options.method == "meissel" ?
                PrimeCounting::PartialSieveDetail::IntegerCbrt(options.x) :
                PrimeCounting::PartialSieveDetail::IntegerSqrt(limit);
            a = tables.pi()[cutoff];
            primeCount = tables.primeCount();
            piSize = tables.pi().size();
            piCapacity = tables.pi().capacity();
            primesSize = tables.primes().size();
            primesCapacity = tables.primes().capacity();
        }
    };

    struct Row {
        std::int64_t repeat = 0;
        std::uint64_t iterations = 1;
        std::string cacheState;
        double setup = 0, phiMs = 0, p2Ms = 0, p3Ms = 0;
        double warmup = 0, teardown = 0, total = 0;
        std::int64_t phi = 0, p2 = 0, p3 = 0, primeCount = 0;
        bool hasPhi = false, hasP3 = false, hasPrimeCount = false;
        TableInfo tables;
    };

    void Header() {
        std::cout << "mode,method,x,threads,repeat,iterations,cache_state,priming,"
            "setup_ms,phi_ms,p2_ms,p3_ms,warmup_ms,teardown_ms,total_ms,"
            "phi,p2,p3,prime_count,a,table_limit,table_prime_count,"
            "pi_size,pi_capacity,pi_size_bytes,pi_capacity_bytes,"
            "primes_size,primes_capacity,primes_size_bytes,primes_capacity_bytes,"
            "table_size_bytes,table_capacity_bytes\n";
    }

    void Print(const Options& options, const Row& row) {
        const auto& table = row.tables;
        const auto piBytes = table.piSize * sizeof(std::int32_t);
        const auto piCapacityBytes = table.piCapacity * sizeof(std::int32_t);
        const auto primeBytes = table.primesSize * sizeof(std::int64_t);
        const auto primeCapacityBytes = table.primesCapacity * sizeof(std::int64_t);
        std::cout << options.mode << ',' << options.method << ',' << options.x << ','
                  << options.threads << ',' << row.repeat << ',' << row.iterations << ','
                  << row.cacheState << ',' << options.priming << ',' << std::setprecision(12)
                  << row.setup << ',' << row.phiMs << ',' << row.p2Ms << ',' << row.p3Ms << ','
                  << row.warmup << ',' << row.teardown << ',' << row.total << ',';
        if (row.hasPhi) std::cout << row.phi;
        std::cout << ',' << row.p2 << ',';
        if (row.hasP3) std::cout << row.p3;
        std::cout << ',';
        if (row.hasPrimeCount) std::cout << row.primeCount;
        std::cout << ',' << table.a << ',' << table.limit << ',' << table.primeCount << ','
                  << table.piSize << ',' << table.piCapacity << ','
                  << piBytes << ',' << piCapacityBytes << ','
                  << table.primesSize << ',' << table.primesCapacity << ','
                  << primeBytes << ',' << primeCapacityBytes << ','
                  << piBytes + primeBytes << ',' << piCapacityBytes + primeCapacityBytes << '\n';
    }

    void Check(std::int64_t actual, std::int64_t expected, const char* label) {
        if (actual != expected) throw std::runtime_error(std::string(label) +
            " mismatch: actual=" + std::to_string(actual) +
            " expected=" + std::to_string(expected));
    }

    void BeginP2() {
#ifdef P2_PROFILE_METRICS
        P2ProfileMetrics::begin();
#endif
    }

    void EndP2() {
#ifdef P2_PROFILE_METRICS
        P2ProfileMetrics::end();
#endif
    }

    void ReportP2() {
#ifdef P2_PROFILE_METRICS
        P2ProfileMetrics::report(std::cout);
#endif
    }

    void Phases(const Options& options) {
        for (std::int64_t repeat = 0; repeat < options.amount; ++repeat) {
            Row row;
            row.repeat = repeat;
            row.cacheState = repeat ? "process-warm" : "process-first";
            const auto totalStart = Clock::now();
            auto tables = std::unique_ptr<PhiTables>(new PhiTables(
                PrimeCounting::PartialSieveDetail::IntegerSqrt(options.x)));
            row.tables = TableInfo(options, *tables);
            const auto setupEnd = Clock::now();
            row.setup = Milliseconds(totalStart, setupEnd);
            row.phi = PrimeCounting::PhiParallel(options.x, row.tables.a, options.threads, *tables);
            const auto phiEnd = Clock::now();
            row.phiMs = Milliseconds(setupEnd, phiEnd);
            row.hasPhi = true;
            BeginP2();
            const auto p2Start = Clock::now();
            row.p2 = PrimeCounting::P2(options.x, row.tables.a, options.threads, *tables);
            row.p2Ms = Milliseconds(p2Start, Clock::now());
            EndP2();
            if (options.method == "lehmer") {
                const auto p3Start = Clock::now();
                row.p3 = PrimeCounting::P3(options.x, row.tables.a, options.threads, *tables);
                row.p3Ms = Milliseconds(p3Start, Clock::now());
            }
            row.hasP3 = true; // At the Meissel cutoff, P3 is identically zero.
            row.primeCount = row.phi + row.tables.a - 1 - row.p2 - row.p3;
            row.hasPrimeCount = true;
            const auto teardownStart = Clock::now();
            tables.reset();
            const auto totalEnd = Clock::now();
            row.teardown = Milliseconds(teardownStart, totalEnd);
            row.total = Milliseconds(totalStart, totalEnd);
            Check(row.primeCount, options.expected, "prime count");
            Print(options, row);
            ReportP2();
        }
    }

    void Isolated(const Options& options) {
        const auto totalStart = Clock::now();
        auto tables = std::unique_ptr<PhiTables>(new PhiTables(
            PrimeCounting::PartialSieveDetail::IntegerSqrt(options.x)));
        const TableInfo info(options, *tables);
        const auto setupEnd = Clock::now();
        const auto setupMs = Milliseconds(totalStart, setupEnd);
        std::int64_t phi = 0;
        double phiMs = 0;
        if (options.priming == "prime-phi") {
            const auto phiStart = Clock::now();
            phi = PrimeCounting::PhiParallel(options.x, info.a, options.threads, *tables);
            phiMs = Milliseconds(phiStart, Clock::now());
        }

        if (options.mode == "sample") {
            Row row;
            row.tables = info;
            row.setup = setupMs;
            row.phi = phi;
            row.phiMs = phiMs;
            row.hasPhi = options.priming == "prime-phi";
            row.cacheState = "warm-sample";
            const auto warmupStart = Clock::now();
            row.p2 = PrimeCounting::P2(options.x, info.a, options.threads, *tables);
            row.warmup = Milliseconds(warmupStart, Clock::now());
            Check(row.p2, options.expected, "P2 warm-up");
            row.iterations = 0;
            BeginP2();
            const auto sampleStart = Clock::now();
            do {
                const auto p2Start = Clock::now();
                row.p2 = PrimeCounting::P2(options.x, info.a, options.threads, *tables);
                row.p2Ms += Milliseconds(p2Start, Clock::now());
                Check(row.p2, options.expected, "P2 sample");
                ++row.iterations;
            } while (Milliseconds(sampleStart, Clock::now()) < options.amount);
            EndP2();
            const auto teardownStart = Clock::now();
            tables.reset();
            const auto totalEnd = Clock::now();
            row.teardown = Milliseconds(teardownStart, totalEnd);
            row.total = Milliseconds(totalStart, totalEnd);
            Print(options, row);
            ReportP2();
            return;
        }

        for (std::int64_t repeat = 0; repeat < options.amount; ++repeat) {
            Row row;
            row.repeat = repeat;
            row.tables = info;
            row.phi = phi;
            row.hasPhi = options.priming == "prime-phi";
            row.cacheState = repeat ? "p2-warm" :
                (row.hasPhi ? "phi-primed-p2-first" : "p2-first");
            row.setup = repeat ? 0 : setupMs;
            row.phiMs = repeat ? 0 : phiMs;
            const auto iterationStart = repeat ? Clock::now() : totalStart;
            BeginP2();
            const auto p2Start = Clock::now();
            row.p2 = PrimeCounting::P2(options.x, info.a, options.threads, *tables);
            row.p2Ms = Milliseconds(p2Start, Clock::now());
            EndP2();
            if (repeat + 1 == options.amount) {
                const auto teardownStart = Clock::now();
                tables.reset();
                row.teardown = Milliseconds(teardownStart, Clock::now());
            }
            row.total = Milliseconds(iterationStart, Clock::now());
            Check(row.p2, options.expected, "P2");
            Print(options, row);
            ReportP2();
        }
    }
}

int main(int argc, char** argv) {
    try {
        const Options options(argc, argv);
        Header();
        if (options.mode == "phases") Phases(options);
        else Isolated(options);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
