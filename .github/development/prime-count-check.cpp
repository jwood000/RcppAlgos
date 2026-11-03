// Native driver for the extracted prime-counting core; see testing.md.
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdint>

namespace PrimeCounting {
    std::int64_t MainPrimeCount(std::int64_t, int, int);
}

int main(int argc, char** argv) {
    if (argc != 5) return 2;
    const std::int64_t n = std::strtoll(argv[1], nullptr, 10);
    const int threads = std::atoi(argv[2]);
    const int repeats = std::atoi(argv[3]);
    const std::int64_t expected = std::strtoll(argv[4], nullptr, 10);
    if (n < 100000 || threads < 1 || repeats < 1) return 2;

    for (int i = 0; i < repeats; ++i) {
        const auto start = std::chrono::steady_clock::now();
        const auto result = PrimeCounting::MainPrimeCount(n, threads, threads);
        const double ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
        std::cout << n << "," << threads << "," << i << ","
                  << result << "," << std::setprecision(9) << ms << "\n";
        if (result != expected) return 1;
    }
}
