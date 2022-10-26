#pragma once

#include <cstdint>
#include <cstdlib>
#include <bitset>
#include <limits>

constexpr int L1_CACHE_SIZE = 32768;
constexpr double dblIntMax = std::numeric_limits<int>::max();

// std::pow(2, 53) - 1 is not constexpr as of C++11
constexpr double Significand53 = 9007199254740991.0;

// Based off the internal limitations of sample, we
// cannot utilize the full range of 53-bit significand
// precision. Here is the condition from `do_dample2`:
//     if (!R_FINITE(dn) || dn < 0 || dn > 4.5e15 || (k > 0 && dn == 0))
// Here is the source (line ~2100):
//     https://github.com/wch/r-source/blob/trunk/src/main/unique.c
constexpr double SampleLimit = 4500000000000000.0;

// in R console: print(sqrt(.Machine$double.eps), digits = 16)
// [1] 0.00000001490116119384766
// Which is also 2^(-26)
constexpr double DefaultTolerance = std::numeric_limits<float>::epsilon() / 8.0;

// 2^63 -->> 922337203685477580
constexpr double my63Max = 922337203685477580.0;

// std::sqrt(std::pow(2, 63)) -->> 3037000499.97605
constexpr std::int64_t Sqrt63Max = static_cast<std::int64_t>(3037000499.97);

constexpr int maxVecSize = std::numeric_limits<int>::max() / 2;
constexpr std::size_t intSize = sizeof(int);
constexpr std::size_t numb = 8 * intSize;
constexpr unsigned long int oneThousand = 1000u;
constexpr size_t wordSize = sizeof(std::bitset<1>) * 8;

enum class VecType {
    Integer   = 1,
    Numeric   = 2,
    Logical   = 3,
    Character = 4,
    Complex   = 5,
    Raw       = 6
};
