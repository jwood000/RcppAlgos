// Complete native driver timings, including table construction/destruction.
// CLI: benchmark legendre|meissel|lehmer n threads repeats expected
// CSV: method,n,threads,repeat,count,total_ms. Repeat zero is the first call;
// later calls reuse the process-wide Phi cache but rebuild their own tables.
// Use a fresh process for each method/input/thread configuration.
#include "NumbersUtils/PrimeCount.h"
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    try {
        if (argc != 6) throw std::runtime_error(
            "usage: benchmark legendre|meissel|lehmer n threads repeats expected"
        );

        const std::string method(argv[1]);
        using Counter = std::int64_t (*)(std::int64_t, int, int);
        Counter count;
        if (method == "legendre") count = PrimeCounting::PiLegendre;
        else if (method == "meissel") count = PrimeCounting::PiMeissel;
        else if (method == "lehmer") count = PrimeCounting::PiLehmer;
        else throw std::runtime_error("unknown prime-count method");

        const std::int64_t n = std::stoll(argv[2]);
        const int threads = std::stoi(argv[3]);
        const int repeats = std::stoi(argv[4]);
        const std::int64_t expected = std::stoll(argv[5]);
        constexpr std::int64_t maxInput = (std::int64_t(1) << 53) - 1;
        if (n < (method == "legendre" ? 100000 : 0) || n > maxInput ||
            threads < 1 || repeats < 1 || expected < 0) {
            throw std::runtime_error("invalid benchmark arguments");
        }

        using Clock = std::chrono::steady_clock;
        for (int repeat = 0; repeat < repeats; ++repeat) {
            const auto start = Clock::now();
            const std::int64_t result = count(n, threads, threads);
            const double milliseconds = std::chrono::duration<double, std::milli>(
                Clock::now() - start
            ).count();
            if (result != expected) {
                throw std::runtime_error("prime count mismatch: expected " +
                    std::to_string(expected) + ", got " + std::to_string(result));
            }
            std::cout << method << ',' << n << ',' << threads << ',' << repeat
                      << ',' << result << ',' << std::setprecision(9)
                      << milliseconds << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
