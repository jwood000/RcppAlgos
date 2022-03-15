#include "Constraints/ConstraintsUtils.h"

// Consider finding partitions of 50 of length 5 with no zeros.
// The first few results are:
//                  1 1 1 1 46
//                  1 1 1 2 45
//                  1 1 1 3 44
//
// If we cont. until the num in the 3rd col changes, we have:
//              [22,] 1 1 1 22 25
//              [23,] 1 1 1 23 24
//              [24,] 1 1 2  2 44
//
// That's 23 results, which is expected as we are decrementing
// the right most number and incrementing the number in the
// 4th column. Thus we have:
//                   (46 - n) > (1 + n)
//              ==>        45 > 2n
//              ==>      22.5 > n
//              ==>         n = 22
//
// Including the first result (i.e. 1 46), we have n + 1 = 23.
// In general, it can be show that when we are presented with
// a situaiton of the form: 1 1 ... 1 m, then the number of
// results for the "1st" group (i.e. while the 2nd to the last
// column stays the same), is given by floor((1 + m) / 2).
//
// Continuing on with the next group, we have:
//                   1 1 2 2 44
//
// We need to find the number of result of (2 44). This is
// equivalent to (1 43). The next few groups are (3 42), (4, 40),
// (5, 38), and (6, 36). These all get converted to (1, 40),
// (1, 37), (1, 34), and (1, 31). The number of results for
// each are given by (assume integer division and don't
// forget to add 1):
//          44 / 2, 41 / 2, 38 / 2, 35 / 2, 32 / 2
//
// The total number of result for all groups is given by:
//          47 / 2 + 44 / 2 + ... + 5 / 2 + 2 / 2
//
// This can be converted to:
// {(47 - 0) + (47 - 3) + ... + (47 - 42) + (47 - 45)} / 2
//
// The 8 below comes from adding 1 to every other term
// that needs to be removed (leave as a proof to the reader)
//   ==>>  (47 * 16 - 3 * (0 + 1 + 2 + .... + 14 + 15) - 8) / 2
//   ==>>    (47 * 16 - 3 * (16 * 15 / 2) - 8) / 2

std::uint64_t SumSection(std::uint64_t n) {
    const std::uint64_t nIter = n / 3;
    const std::uint64_t sumOne = nIter * (n - 1);
    const std::uint64_t sumTwo = 3 * (nIter * (nIter - 1)) / 2;
    return (sumOne - sumTwo - nIter / 2) / 2;
}

int GetMaxWidth(double target) {
    const double discriminant = 1.0 + 8.0 * target;
    int max_width = (-1 + std::sqrt(discriminant)) / 2;
    return max_width;
}

void CheckMultIsInt(double x, double y) {
    if ((x * y) > dblIntMax) {
        Rf_error("Sorry, this case is too large!");
    }
}
