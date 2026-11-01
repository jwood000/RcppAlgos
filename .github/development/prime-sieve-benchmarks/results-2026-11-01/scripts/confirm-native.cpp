// Compiled by confirm-native.py against its retained, minimally instrumented clone.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "sieve-clone.cpp"

struct Configuration {
    std::string name;
    std::size_t bytes;
};

struct Run {
    std::vector<double> primes;
    double milliseconds;
    int actualThreads;
    std::size_t segmentSpan;
};

Run Execute(std::int_fast64_t lower, std::int_fast64_t upper, int requested,
            int maxThreads, int maxCores, std::size_t bytes, bool inspect) {
    std::vector<std::vector<double>> chunks(requested);
    Run result;
    bool parallel = false;
    CacheInspection::Inspector inspector;
    inspector.bigSegmentSpan = (bytes * 8 / 30030) * 30030;

    const auto start = std::chrono::steady_clock::now();
    PrimeSieve::PrimeSieveMain(chunks, result.primes, lower, upper, parallel,
        requested, maxThreads, maxCores, inspect ? &inspector : nullptr);
    result.milliseconds = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();

    // Allocation, concatenation, and exact-output verification are outside timing.
    if (parallel) {
        std::size_t count = 0;
        for (const auto& chunk : chunks) count += chunk.size();
        result.primes.reserve(count);
        for (const auto& chunk : chunks) {
            result.primes.insert(result.primes.end(), chunk.begin(), chunk.end());
        }
    }
    result.actualThreads = inspect ? inspector.effectiveThreads : (parallel ? -1 : 1);
    // Default sizing uses requested threads before any workload-driven fallback.
    const int policyThreads = (requested > 1 && maxThreads > 1) ?
        std::min(requested, maxThreads) : requested;
    result.segmentSpan = bytes ? inspector.bigSegmentSpan :
        std::max<std::size_t>(1, (maxThreads / policyThreads) * 64) * 30030;
    return result;
}

int main(int argc, char** argv) {
    try {
        if (argc != 12) throw std::runtime_error("expected output, repeats, seed, policy, and three configurations");
        const int repeats = std::stoi(argv[2]);
        const unsigned seed = static_cast<unsigned>(std::stoul(argv[3]));
        const int maxThreads = std::stoi(argv[4]);
        const int maxCores = std::stoi(argv[5]);
        std::vector<Configuration> configurations;
        for (int i = 6; i < argc; i += 2) {
            configurations.push_back({argv[i], static_cast<std::size_t>(std::stoull(argv[i + 1]))});
        }
        std::ofstream output(argv[1]);
        if (!output) throw std::runtime_error("cannot open output CSV");
        output << "base,lower,upper,width,requested,actual,config,budget_bytes,segment_span,"
                  "round,order,native_ms,prime_count,matches_reference\n";
        output << std::setprecision(17);
        std::mt19937 random(seed);
        std::vector<std::size_t> order(configurations.size());
        std::iota(order.begin(), order.end(), 0);

        for (std::int_fast64_t base : {4000000000000000LL, 8000000000000000LL}) {
            const std::int_fast64_t lower = base + 1234567;
            for (std::int_fast64_t width : {1000LL, 1000000LL, 100000000LL}) {
                const std::int_fast64_t upper = lower + width;
                for (int requested : {1, 4}) {
                    std::cerr << "Native confirmation: lower=" << lower << " width=" << width
                              << " requested=" << requested << '\n';
                    const Run reference = Execute(lower, upper, requested,
                        maxThreads, maxCores, 0, false);
                    std::shuffle(order.begin(), order.end(), random);
                    for (std::size_t i : order) {
                        const Run warmup = Execute(lower, upper, requested,
                            maxThreads, maxCores, configurations[i].bytes, true);
                        if (warmup.primes != reference.primes) {
                            throw std::runtime_error("warmup output mismatch: " + configurations[i].name);
                        }
                    }
                    for (int round = 1; round <= repeats; ++round) {
                        std::shuffle(order.begin(), order.end(), random);
                        for (std::size_t position = 0; position < order.size(); ++position) {
                            const Configuration& config = configurations[order[position]];
                            const Run result = Execute(lower, upper, requested,
                                maxThreads, maxCores, config.bytes, true);
                            const bool matches = result.primes == reference.primes;
                            output << base << ',' << lower << ',' << upper << ',' << width << ','
                                   << requested << ',' << result.actualThreads << ',' << config.name << ',';
                            if (config.bytes) output << config.bytes;
                            else output << "NA";
                            output << ',' << result.segmentSpan << ',' << round << ',' << position + 1
                                   << ',' << result.milliseconds << ',' << result.primes.size() << ','
                                   << (matches ? "TRUE" : "FALSE") << '\n';
                            if (!matches) {
                                output.flush();
                                throw std::runtime_error("timed output mismatch: " + config.name);
                            }
                        }
                        output.flush();
                        if (!output) throw std::runtime_error("failed writing native confirmation CSV");
                    }
                }
            }
        }
        std::cerr << "All native confirmation outputs matched.\n";
    } catch (const std::exception& error) {
        std::cerr << "Native confirmation failed: " << error.what() << '\n';
        return 1;
    }
}
