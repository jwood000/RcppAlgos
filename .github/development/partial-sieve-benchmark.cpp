// Native timing driver; prepared-table setup is reported separately.
// Link P2.cpp, P3.cpp, Phi.cpp and Eratosthenes.cpp. See testing.md.
#include "NumbersUtils/PartialSieve.h"
#include "NumbersUtils/PartialSieveUtils.h"
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
    using Clock = std::chrono::steady_clock;
    double Milliseconds(Clock::time_point start) {
        return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    }

    // Straight implementation of the agreed P2 formula. Every out-of-table
    // quotient gets a separate serial Legendre/Phi evaluation, reusing tables.
    // This is a diagnostic baseline, not another production counting method.
    std::int64_t P2Reference(std::int64_t x, std::int64_t a,
                             const PrimeCounting::PhiTables& tables) {
        using namespace PrimeCounting;
        std::int64_t sum = 0;
        const auto& primes = tables.primes();
        const auto& pi = tables.pi();
        const std::int64_t b = pi[PartialSieveDetail::IntegerSqrt(x)];
        for (std::int64_t i = a + 1; i <= b; ++i) {
            const auto quotient = x / primes[i];
            const auto c = pi[PartialSieveDetail::IntegerSqrt(quotient)];
            const auto count = quotient <= tables.limit() ? pi[quotient] :
                Phi(quotient, c, tables) + c - 1;
            sum += count - (i - 1);
        }
        return sum;
    }
}

int main(int argc, char** argv) {
    try {
        if (argc != 7) throw std::runtime_error(
            "usage: benchmark p2|p3|reference|verify x root_power threads repeats expected"
        );
        const std::string method(argv[1]);
        const std::int64_t x = std::stoll(argv[2]);
        const int power = std::stoi(argv[3]);
        const int threads = std::stoi(argv[4]);
        const int repeats = std::stoi(argv[5]);
        const std::int64_t expected = std::stoll(argv[6]);
        if (x < 8 || x > PrimeCounting::PartialSieveDetail::MAX_X ||
            (power != 3 && power != 4) || threads < 1 || repeats < 1) {
            throw std::runtime_error("invalid benchmark arguments");
        }
        using namespace PrimeCounting;
        const auto setupStart = Clock::now();
        const PhiTables tables(PartialSieveDetail::IntegerSqrt(x));
        const auto cutoff = power == 3 ? PartialSieveDetail::IntegerCbrt(x) :
            PartialSieveDetail::IntegerSqrt(PartialSieveDetail::IntegerSqrt(x));
        const auto a = tables.pi()[cutoff];
        const auto setup = Milliseconds(setupStart);
        if (method == "verify") {
            const auto p2 = P2(x, a, threads, tables);
            const auto p3 = P3(x, a, threads, tables);
            const auto phi = Phi(x, a, tables);
            if (phi + a - 1 - p2 - p3 != expected) {
                throw std::runtime_error("known prime count identity failed");
            }
            std::cout << p2 << ',' << p3 << ',' << phi << ',' << a << '\n';
            return 0;
        }
        for (int repeat = 0; repeat < repeats; ++repeat) {
            const auto start = Clock::now();
            std::int64_t count;
            if (method == "p2") count = P2(x, a, threads, tables);
            else if (method == "p3") count = P3(x, a, threads, tables);
            else if (method == "reference") count = P2Reference(x, a, tables);
            else throw std::runtime_error("unknown benchmark method");
            const auto ms = Milliseconds(start);
            if (count != expected) throw std::runtime_error("partial count mismatch");
            std::cout << method << ',' << x << ',' << power << ',' << threads << ','
                      << repeat << ',' << count << ',' << std::setprecision(9)
                      << setup << ',' << ms << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
