#include "Combinations/NthCombination.h"
#include "Combinations/BigComboCount.h"
#include "Combinations/ComboCount.h"
#include <algorithm> // std::sort
#include <numeric>   // std::iota
#include "CppConvert/GmpxxCopy.h"
#include <cstdint>
#include <limits>
#include <vector>
#include <cmath>

// ******* Overview of the Crucial Part of the Algorithm *******
// -------------------------------------------------------------
//
// Initial Setup:
//
// idx1 = (r - 1) * grpSize - 1      Penultimate upper bound (inclusive)
//
// idx2 = z.size() - 1               Last upper bound (inclusive)
//
// low_one = (r - 2) * grpSize + 1   One plus the lower bound in the section
//                                   where idx1 is present
//
// low_one is one plus the lower bound in the idx1 section, so to obtain
// the current upper bound, we must first add the size of a section
// (i.e. grpSize) and subtract one. We can now compute the length we need
// to reset v by subtracting idx1. E.g.
//
// Given a portion of v w/ low_one = 9, grpSize = 4, idx1 = 9, 6 groups
// (24 subjects) and base 0:
//
//              prev sections   lower bound (index = 8)
//                  /  \        |
//            ............. 8 | 9 12 23 24 | 10 20 21 22 | 11 ...
//                                 |
//                               idx1 (equal to low_one, in this case)
//
// Sort v past idx1:
//                              size of left range (len_rng - 1)
//                                    /  \
//                                   |    |
//                      ... 8 | 9 12 10 11 | 13 14 15 16 | 17...
//
// Initially idx3 is set to (idx1 + 1), so for this example idx3 = 10. This
// together with low_one = 9 and grpSize = 4, we have all of the ingredients to
// obtain len_rng. We can calculate len_rng = low_one + grpSize - idx3 = 3.
// Once we set idx3 below, we need to rotate at the index just after idx3, thus
// the distance (zbeg + idx3 + len_rng) - (zbeg + idx3 + 1) = 2.
//
// Determine the index, idx3, such that v[idx3] > v[idx1]
//
//                      ... 8 | 9 12 10 11 | 13 14 15 16 | 17 ...
//                                 |          |
//                               idx1       idx3
//
// Swap idx1 and idx3:
//                      ... 8 | 9 13 10 11 | 12 14 15 16 | 17...
//
// Move enough indices after idx1 to fill that specific group:
//
//                      ... 8 | 9 13 __ __ | 10 11 12 14 | 15 16 ...
//
// Identify and move indices that are successively incrementing values of
// v past idx1:
//                      ... 8 | 9 13 14 15 | 10 11 12 16 | 17 ...
//
// The last two steps are accomplished with std::rotate.
//
// From https://en.cppreference.com/w/cpp/algorithm/rotate
//
// std::rotate (defined in header <algorithm>)
// template< class ForwardIt >
// ForwardIt rotate( ForwardIt first, ForwardIt n_first, ForwardIt last );
//
// 1) Performs a left rotation on a range of elements.
//    Specifically, std::rotate swaps the elements in the range [first, last)
//    in such a way that the element n_first becomes the first element of the
//    new range and n_first - 1 becomes the last element.
//
// This completes the algorithm.

bool nextComboGroup(std::vector<int> &z, int grpSize,
                    int idx1, int idx2, int low_one) {

    while (idx2 > idx1 && z[idx2] > z[idx1]) {
        --idx2;
    }

    if ((idx2 + 1) < static_cast<int>(z.size())) {
        // Formally, we had the following conditional:
        //
        // if (z[idx2 + 1] > z[idx1]) {std::swap(z[idx1], z[idx2 + 1]);}
        //
        // This is unnecessary because we know that idx2 always starts out as
        // z.size() - 1 and if the while loop is never entered, idx2 + 1 will
        // be equal to z.size(), meaning:
        //
        //       (idx2 + 1) < static_cast<int>(z.size()) -->> false
        //
        // and this section will not be entered.
        //
        // Otherwise, both of the checks are true at least once. That is,
        // idx2 > idx1 and z[idx2] > z[idx1] for at least one iteration. When
        // the check becomes false, we know that the previous iteration was
        // true, thus z[idx2 + 1] > z[idx1] is true. Thus we can safely
        // proceed with our swap.

        std::swap(z[idx1], z[idx2 + 1]);
        return true;
    }

    const auto zbeg = z.begin();

    while (idx1 > 0) {
        const int tipPnt = z[idx2];

        while (idx1 > low_one && tipPnt < z[idx1]) {
            --idx1;
        }

        if (tipPnt > z[idx1]) { // **Crucial Part**
            int idx3 = idx1 + 1;
            std::sort(zbeg + idx3, z.end());

            // length of left range plus one. The plus one is needed as we are
            // rotating at a pivot just to the right (i.e. plus one).
            const int len_rng = low_one + grpSize - idx3;

            while (z[idx3] < z[idx1]) {
                ++idx3;
            }

            std::swap(z[idx3], z[idx1]);
            std::rotate(zbeg + idx1 + 1,
                        zbeg + idx3 + 1, zbeg + idx3 + len_rng);
            return true;
        } else {
            idx1    -= 2;
            idx2    -= grpSize;
            low_one -= grpSize;
        }
    }

    return false;
}

// The general algorithm where the length of each group can vary. The skeleton
// of the algorithm is exactly the same as the special case. The only
// difference is that we are aware of what group we are in as well as the size
// of each group, hence the vector grpSize (above, grpSize is a constant).
bool nextComboGroup(std::vector<int> &z,
                    const std::vector<int> &grpSize,
                    int idx1, int idx2, int lbound) {

    while (idx2 > idx1 && z[idx2] > z[idx1]) {
        --idx2;
    }

    if ((idx2 + 1) < static_cast<int>(z.size())) {
        std::swap(z[idx1], z[idx2 + 1]);
        return true;
    }

    const auto zbeg = z.begin();

    // Start at penultimate group
    for (int g_idx = grpSize.size() - 2; g_idx >= 0; --g_idx) {

        const int tipPnt = z[idx2];

        while (idx1 > lbound && tipPnt < z[idx1]) {
            --idx1;
        }

        if (z[idx1] < tipPnt) {
            int idx3 = idx1 + 1;
            std::sort(zbeg + idx3, z.end());

            const int len_rng = lbound + grpSize[g_idx] - idx1;

            while (z[idx3] < z[idx1]) {
                ++idx3;
            }

            std::swap(z[idx3], z[idx1]);
            std::rotate(zbeg + idx1 + 1,
                        zbeg + idx3 + 1, zbeg + idx3 + len_rng);
            return true;
        } else if (g_idx > 0) {
            idx2   -= grpSize[g_idx + 1];
            lbound -= grpSize[g_idx - 1];
            --idx1;
        }
    }

    return false;
}

double numGroupCombs(int n, int numGroups, int grpSize) {

    double result = 1;

    for (double i = n; i > numGroups; --i) {
        result *= i;
    }

    if (result < std::numeric_limits<double>::max()) {
        double myDiv = 1;

        for (double i = 2; i <= grpSize; ++i) {
            myDiv *= i;
        }

        result /= std::pow(myDiv, numGroups);
        return std::round(result);
    } else {
        return std::numeric_limits<double>::infinity();
    }
}

void numGroupCombsGmp(mpz_class &result, int n,
                      int numGroups, int grpSize) {

    for (int i = n; i > numGroups; --i) {
        result *= i;
    }

    mpz_class myDiv(1);

    for (int i = 2; i <= grpSize; ++i) {
        myDiv *= i;
    }

    mpz_pow_ui(myDiv.get_mpz_t(), myDiv.get_mpz_t(), numGroups);
    mpz_divexact(result.get_mpz_t(), result.get_mpz_t(), myDiv.get_mpz_t());
}

std::vector<int> nthComboGroup(int n, int grpSize, int r,
                               double myIndex, double total) {

    int s = n - 1;
    const int g = grpSize - 1;
    std::int64_t temp   = static_cast<std::int64_t>(nChooseK(s, g));
    std::int64_t secLen = static_cast<std::int64_t>(total) / temp;

    std::vector<int>  res(n, 0);
    std::vector<char> idx_used(n, 0);
    std::vector<int>  v(s);
    std::iota(v.begin(), v.end(), 1);

    int myMin = 0;
    std::int64_t ind1 = myIndex;
    std::int64_t ind2 = myIndex;

    mpz_class mpzDefault;

    for (int j = 0; j < (r - 1); ++j) {
        ind2 = ind2 / secLen;
        res[j * grpSize]  = myMin;
        idx_used[myMin] = 1;
        const std::vector<int> comb = (g == 1) ? std::vector<int>(1, ind2) :
            nthComb(s, g, ind2, mpzDefault, v);

        for (int k = j * grpSize + 1, i = 0; i < g; ++k, ++i) {
            res[k] = v[comb[i]];
            idx_used[res[k]] = 1;
        }

        v.clear();

        for (int i = 1; i < n; ++i) {
            if (!idx_used[i]) {
                v.push_back(i);
            }
        }

        myMin = v.front();
        v.erase(v.begin());
        ind1 -= ind2 * secLen;
        ind2 = ind1;

        s -= grpSize;
        temp = nChooseK(s, g);
        secLen /= temp;
    }

    res[(r - 1) * grpSize] = myMin;

    for (int k = (r - 1) * grpSize + 1, i = 0; i < g; ++k, ++i) {
        res[k] = v[i];
    }

    return res;
}

std::vector<int> nthComboGroupGmp(int n, int grpSize, int r,
                                  const mpz_class &lowerMpz,
                                  const mpz_class &computedRowMpz) {
    mpz_class ind1(lowerMpz);
    mpz_class ind2(lowerMpz);

    int s = n - 1;
    const int g = grpSize - 1;

    mpz_class temp(1);
    mpz_class secLen(1);

    nChooseKGmp(temp, s, g);
    mpz_divexact(secLen.get_mpz_t(), computedRowMpz.get_mpz_t(),
                 temp.get_mpz_t());

    std::vector<int>  res(n, 0);
    std::vector<char> idx_used(n, 0);
    std::vector<int>  v(s);
    std::iota(v.begin(), v.end(), 1);

    int myMin = 0;
    constexpr double dblDefault = 0;

    for (int j = 0; j < (r - 1); ++j) {
        ind2 /= secLen;
        res[j * grpSize] = myMin;
        idx_used[myMin] = 1;

        const std::vector<int> comb = (g == 1) ?
            std::vector<int>(1, ind2.get_si()) :
            nthCombGmp(s, g, dblDefault, ind2, v);

        for (int k = j * grpSize + 1, i = 0; i < g; ++k, ++i) {
            res[k] = v[comb[i]];
            idx_used[res[k]] = 1;
        }

        v.clear();

        for (int i = 1; i < n; ++i) {
            if (!idx_used[i]) {
                v.push_back(i);
            }
        }

        myMin = v.front();
        v.erase(v.begin());
        temp = ind2 * secLen;
        ind1 -= temp;
        ind2 = ind1;

        s -= grpSize;
        nChooseKGmp(temp, s, g);
        mpz_divexact(secLen.get_mpz_t(), secLen.get_mpz_t(), temp.get_mpz_t());
    }

    res[(r - 1) * grpSize] = myMin;

    for (int k = (r - 1) * grpSize + 1, i = 0; i < g; ++k, ++i) {
        res[k] = v[i];
    }

    return res;
}
