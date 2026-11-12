#include "VariantP2.cpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    const std::size_t budget = std::stoul(argv[1]);
    const PrimeCounting::PhiTables tables(1000000);
    for (int repeat = 0; repeat < 4; ++repeat) {
        const auto start = std::chrono::steady_clock::now();
        const auto result = PrimeCounting::countPhiPrefix(1000000000000LL, 3305, 3369, budget, tables);
        const auto end = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(end-start).count();
        std::cout << repeat << ',' << result << ',' << std::setprecision(12) << ms << '\n';
    }
}
