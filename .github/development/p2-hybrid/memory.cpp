// Memory/work probe for P2Hybrid.cpp. Compile separately with
// -DP2_HYBRID_MEMORY -DP2_HYBRID_MULTIPLIER=2,4,8. This is not a timing binary.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include "NumbersUtils/Phi.h"
#include "hybrid_metrics.h"
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
    std::int64_t integer(const char* text) {
        std::size_t used = 0;
        const std::string value(text);
        const std::int64_t result = std::stoll(value, &used);
        if (used != value.size()) throw std::runtime_error("invalid integer argument");
        return result;
    }
}

int main(int argc, char** argv) {
    try {
        if (argc != 6) throw std::runtime_error(
            "usage: memory meissel|lehmer x threads repeats expected_pi");
        const std::string method(argv[1]);
        const std::int64_t x = integer(argv[2]);
        const std::int64_t requested = integer(argv[3]);
        const std::int64_t repeats = integer(argv[4]);
        const std::int64_t expected = integer(argv[5]);
        using namespace PrimeCounting;
        if ((method != "meissel" && method != "lehmer") ||
            x < 2 || x > PartialSieveDetail::MAX_X ||
            requested < 1 || requested > std::numeric_limits<int>::max() ||
            repeats < 1 || repeats > std::numeric_limits<int>::max() ||
            expected < 0 || expected > x) {
            throw std::runtime_error("invalid probe arguments");
        }
        const int threads = static_cast<int>(requested);
        std::cout << "multiplier,method,x,threads,repeat,count,table_limit,stream_limit,"
                     "phi_queries,stream_queries,direct_queries,phi_workers,stream_workers,"
                     "seed_queries,phi_seed_queries,sieve_capacity_bytes,next_capacity_bytes,"
                     "phi_sum_bytes,stream_sum_bytes,conservative_scratch_bytes,"
                     "table_capacity_bytes\n";
        for (std::int64_t repeat = 0; repeat < repeats; ++repeat) {
            const std::int64_t root = PartialSieveDetail::IntegerSqrt(x);
            const PhiTables tables(root);
            const std::int64_t cutoff = method == "meissel" ?
                PartialSieveDetail::IntegerCbrt(x) : PartialSieveDetail::IntegerSqrt(root);
            const std::int64_t a = tables.pi()[cutoff];
            const std::int64_t phi = PhiParallel(x, a, threads, tables);
            const std::int64_t p2 = P2(x, a, threads, tables);
            const auto& metrics = P2HybridMetrics::current();
            const std::int64_t p3 = method == "lehmer" ? P3(x, a, threads, tables) : 0;
            const std::int64_t count = phi + a - 1 - p2 - p3;
            if (count != expected) throw std::runtime_error("incorrect complete prime count");

            std::uint64_t seeds = 0, phiSeeds = 0, sieve = 0, next = 0;
            std::int64_t streamedQueries = 0;
            for (const auto& interval : metrics.intervals) if (interval.active) {
                ++seeds;
                phiSeeds += interval.phiSeed;
                sieve += interval.sieveBytes;
                next += interval.nextBytes;
                streamedQueries += interval.queries;
            }
            if (seeds != metrics.streamWorkers || streamedQueries != metrics.streamQueries ||
                metrics.phiQueries + metrics.streamQueries + metrics.directQueries !=
                    std::max<std::int64_t>(0, tables.pi()[root] - a)) {
                throw std::runtime_error("inconsistent hybrid work counters");
            }
            // Prefix and streaming phases are sequential. This is a capacity
            // upper bound for tracked scratch payload, since worker lifetimes
            // need not overlap completely. It excludes the diagnostic snapshot
            // itself, ThreadBatch storage, thread stacks, allocator bookkeeping,
            // shared Phi cache, table construction scratch, and process RSS.
            const std::uint64_t peak = std::max(metrics.phiSumBytes,
                metrics.streamSumBytes + sieve + next);
            const std::uint64_t tableBytes = tables.pi().capacity() * sizeof(std::int32_t) +
                tables.primes().capacity() * sizeof(std::int64_t);
            std::cout << P2_HYBRID_MULTIPLIER << ',' << method << ',' << x << ','
                      << threads << ',' << repeat << ',' << count << ','
                      << metrics.tableLimit << ',' << metrics.streamLimit << ','
                      << metrics.phiQueries << ',' << metrics.streamQueries << ','
                      << metrics.directQueries << ',' << metrics.phiWorkers << ','
                      << metrics.streamWorkers << ',' << seeds << ',' << phiSeeds << ','
                      << sieve << ',' << next << ',' << metrics.phiSumBytes << ','
                      << metrics.streamSumBytes << ',' << peak << ',' << tableBytes << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
