#include <algorithm>
#include <numeric>
#include <vector>

template <int one_or_zero>
void NextCompositionRep(std::vector<int> &z, int lastCol) {

    if (z[lastCol] != one_or_zero) {
        --z[lastCol];
        ++z[lastCol - 1];
    } else {
        int j = lastCol - 1;

        while (j > 0 && z[j] == one_or_zero) {
            --j;
        }

        if (j > 0) {
            ++z[j - 1];
            std::reverse(z.begin() + j, z.end());
            --z[lastCol];
        }
    }
}

int NextDistinctBlock(
    const std::vector<int> &v, std::vector<int> &idx, int target, int m
) {

    int lenV = v.size();
    const int testMax = std::accumulate(v.cend() - m, v.cend(), 0);

    // The length is too small
    if (testMax < target) {
        return -2;
    }

    const int testMin = std::accumulate(v.cbegin(), v.cbegin() + m, 0);

    // The length is too long
    if (testMin > target) {
        return -1;
    }

    std::iota(idx.begin(), idx.end(), 0);

    auto GetSum = [](const std::vector<int> &_v,
                     const std::vector<int> &_idx) {
        return std::accumulate(
            _idx.begin(), _idx.end(), 0, [&](int sum, int i) {
                return sum + _v[i];
            }
        );
    };

    int tempSum = GetSum(v, idx);
    bool keepGoing = tempSum != target;
    int partial = target - (tempSum - v[idx.back()]);
    int j = m - 1;

    while (keepGoing) {
        auto lower = std::lower_bound(v.begin() + j, v.end(), partial);

        if (lower != v.end() && *lower == partial) {
            idx[m - 1] = std::distance(v.begin(), lower);
            return 1;
        }

        int k = m - 2;
        int g = lenV - 2;

        while (k > 0 && idx[k] == g) {
            --k;
            --g;
        }

        if (k == 0 && idx.front() == (lenV - m)) {
            return 0;
        }

        bool impossible = true;

        while (impossible && k >= 0) {
            ++idx[k];

            for (int i = k + 1; i < m; ++i) {
                idx[i] = idx[i - 1] + 1;
            }

            tempSum = GetSum(v, idx);
            int maxSum = tempSum;

            for (int i = k; i < m; ++i) {
                maxSum -= v[idx[i]];
                maxSum += v[lenV - (m - i)];
            }

            impossible = (tempSum > target) ||
                (maxSum < target);
            --k;
        }

        if (impossible) {
            return 0;
        }

        keepGoing = tempSum != target;
        partial = target - (tempSum - v[idx.back()]);
        j = idx.back();
    }

    return 1;
}

bool NextRoutine(
    std::vector<int> &z, std::vector<int> &complement,
    int &i1, int &i2, int lastCol, int lastIdx
) {

    auto&& ref_one = z[lastCol - 1];
    auto&& ref_two = z[lastCol];

    // See commentary in else block of res2 == 1 conditional.
    if (i1 < 0) {
        std::swap(ref_one, ref_two);
        return true;
    }

    // This will be used to see if we need to swap the last two elements of z
    bool less_flag = ref_one < ref_two;

    const int target = ref_one + ref_two;
    int check = complement[i1] + complement[i2];

    while (check != target) {
        if (check < target) {
            ++i1;
        } else {
            --i2;
        }

        check = complement[i1] + complement[i2];
    }

    if (less_flag) {
        if (i1 == i2) {
            std::swap(ref_one, ref_two);
            ++i1;
            --i2;
            return true;
        } else if (ref_two - ref_one < complement[i1] - complement[i2]) {
            // In this case we already have figured out what i1 and i2 will
            // be for the iteration, so no need to do anything to them.
            std::swap(ref_one, ref_two);
            return true;
        }
    }

    std::swap(ref_one, complement[i1]);
    std::swap(ref_two, complement[i2]);
    std::sort(complement.begin(), complement.end());

    if (i1 < lastIdx) ++i1; else return false;
    if (i2 > 0) --i2; else return false;
    return true;
}

// The algorithm for NextCompositionDistinct is similar in principle to the
// traditional NextCompositionRep algorithm, however the details are far more
// complicated. For example, in the traditional algo, the idea for each
// iteration is to first check if the last two elements can be altered to
// produce another composition. This is fairly straightforward when repetition
// is allowed. We simply check to see if the last element is minimal. If it
// isn't, we increment the penultimate index and decrement the final index to
// obtain the next composition. If it is minimal, we back track from the
// penultimate index until we reach an index that isn't maximized. Once we do,
// we increment that index and adjust the remaining indices. The algorithm
// continues on until all indices are maximized starting with the first index.
//
// We can call the first part of the traditional algorithm the routine method
// and the second part is setting up the next block so that we can continue
// calling the routine method. We have broken up the algorithm below to reflect
// these ideas.
//
// You will note that for the 1st part, we can't check the last element as we
// were with the traditional algorithm as the final element for each block
// is not constant in the distinct case. Instead, we determine the maximum
// that the penultimate index can be for each block and use that as our check.
// When that index is maximized, we know that we have exhausted the
// compositions in this block and thus need to find the next block.
//
// For the 2nd part, we are following the same pattern as with the traditional
// algorithm, there are just several pieces that are quite evasive. For
// example, we can't just check to see if an index has reached some maximal
// constant value. Each index in the distinct case has its own particular
// maximal value at any given state. A few cases below:
//
// CASE 1
// target: 20; m = 4; global maximum: 14; current iteration: 8 9 2 1
//
// We can't simply increment the 9 above to 10 as there would be no solutions.
// The maximal value for the 2nd index is in fact 9 and we must move to the
// 1st index and increment from 8 to 9.
//
// CASE 2
// target: 20; m = 4; global maximum: 14; current iteration: 3 12 4 1
//
// Again, we can't increment the 12 to 13 as this would lead to no solutions as
// (3 + 13) = 16 ==>> we need to find two elements from:
//
//                 1 2 4 5 6 7 8 9 10 11 13 14
//
// that sum to 4. There is no solution. The correct course is to increment 12
// to 14. These are just a few cases, but hopefully the point is made.

void NextCompositionDistinct(
    std::vector<int> &z, std::vector<int> &complement, int &i1,
    int &i2, int &myMax, int lastCol, int lastIdx, int target
) {

    if (z[lastCol - 1] < myMax) {
        NextRoutine(z, complement, i1, i2, lastCol, lastIdx);
    } else if (lastCol > 1) {
        int m = 2;
        int j = lastCol - m;

        complement.insert(complement.end(), z.end() - m, z.end());
        std::sort(complement.begin(), complement.end());
        bool keepGoing = true;

        while (j >= 0 && keepGoing) {
            int res = 0;
            std::vector<int> idx(m);

            while (res == 0) {
                // We first check to see if there is an element in complement
                // that is larger than z[j]
                auto upper = std::upper_bound(
                    complement.begin(), complement.end(), z[j]
                );

                if (upper != complement.end()) {
                    // We found one! Now, we swap and calculate the partial sum
                    // that we need to find in the current complement.
                    std::swap(*upper, z[j]);
                    int partial = target -
                        std::accumulate(z.cbegin(), z.cend() - m, 0);
                    res = NextDistinctBlock(complement, idx, partial, m);
                } else {
                    res = 100;
                }
            }

            if (res == 1) {
                // We found a viable block.
                keepGoing = false;

                // Set the remaining indices of z and remove these elements
                // from complement. Note, we must do this from the end as doing
                // it from the beginning would cause issues with indexing. E.g.
                //
                //                             0, 1, 2, 3, 4,  5,  6
                // idx = {2, 6}; complement = {3, 5, 6, 7, 8, 10, 11}
                //                   remove these    ^            ^
                //
                // remove idx = 2:
                // complement = {3, 5, 7, 8, 10, 11}
                //
                // remove idx = 6:
                // complement = {3, 5, 7, 8, 10, 11} <<-- error: out of bounds
                //
                // ********************** CORRECT WAY *************************
                //
                // remove idx = 6:
                // complement = {3, 5, 6, 7, 8, 10}
                //
                // remove idx = 2:
                // complement = {3, 5, 7, 8, 10}

                for (int i = m - 1, k = lastCol; i >= 0; --i, --k) {
                    z[k] = complement[idx[i]];
                    complement.erase(complement.begin() + idx[i]);
                }

                // Reset myMax to the current last element of z
                myMax = z[lastCol];

                // Ensure that the current complement affords new solutions
                idx.resize(2);
                int lastTwo = z[lastCol - 1] + z[lastCol];
                int res2 = NextDistinctBlock(complement, idx, lastTwo, 2);

                if (res2 == 1) {
                    i1 = idx.front();
                    i2 = idx.back();
                } else {
                    // This means that there are no further solutions in
                    // complement. The next iteration we simply will swap the
                    // last two elements. By doing this, we will also ensure
                    // that z[lastCol - 1] will be equal to myMax which will
                    // force a new block calculation.
                    i1 = -1;
                }
            } else if (res == 100) {
                // The above did not produce any new fruitful solutions so we
                // must back out further and start finding an even larger
                // block to start testing.
                complement.push_back(z[j]);
                --j;
                ++m;
            }
        }
    }
}

template void NextCompositionRep<0>(std::vector<int>&, int);
template void NextCompositionRep<1>(std::vector<int>&, int);
