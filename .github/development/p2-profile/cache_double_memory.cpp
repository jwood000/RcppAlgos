// Separate memory diagnostic; never use this executable for timing claims.
// Link one memory-probe Phi.cpp and the usual complete native counting sources.
#include "NumbersUtils/PrimeCount.h"
#include "cache_memory_probe.h"
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        if (argc != 6) throw std::runtime_error(
            "usage: cache-memory legendre|meissel|lehmer n threads repeats expected"
        );
        const std::string method(argv[1]);
        using Counter = std::int64_t (*)(std::int64_t, int, int);
        Counter count;
        if (method == "legendre") count = PrimeCounting::PiLegendre;
        else if (method == "meissel") count = PrimeCounting::PiMeissel;
        else if (method == "lehmer") count = PrimeCounting::PiLehmer;
        else throw std::runtime_error("unknown prime-count method");

        const auto n = std::stoll(argv[2]);
        const int threads = std::stoi(argv[3]);
        const int repeats = std::stoi(argv[4]);
        const auto expected = std::stoll(argv[5]);
        constexpr std::int64_t maxInput = (std::int64_t(1) << 53) - 1;
        if (n < (method == "legendre" ? 100000 : 0) || n > maxInput ||
                threads < 1 || repeats < 1 || expected < 0) {
            throw std::runtime_error("invalid memory-probe arguments");
        }

        for (int repeat = 0; repeat < repeats; ++repeat) {
            const auto result = count(n, threads, threads);
            if (result != expected) throw std::runtime_error("prime count mismatch");
            // Native drivers return only after joining workers and destroying
            // their per-call tables. The shared cache remains alive to inspect.
            const auto stats = PrimeCounting::PhiCacheStats();
            if (stats.allocatedPages > stats.possiblePages ||
                    stats.allocatedPageBytes != stats.allocatedPages * stats.pageBytes ||
                    stats.possiblePageBytes != stats.possiblePages * stats.pageBytes) {
                throw std::runtime_error("inconsistent cache memory statistics");
            }
            std::cout << "{\"method\":\"" << method << "\",\"n\":" << n
                      << ",\"threads\":" << threads << ",\"repeat\":" << repeat
                      << ",\"count\":" << result
                      << ",\"allocated_pages\":" << stats.allocatedPages
                      << ",\"page_bytes\":" << stats.pageBytes
                      << ",\"allocated_page_bytes\":" << stats.allocatedPageBytes
                      << ",\"cache_object_bytes\":" << stats.cacheObjectBytes
                      << ",\"possible_pages\":" << stats.possiblePages
                      << ",\"possible_page_bytes\":" << stats.possiblePageBytes << "}\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
