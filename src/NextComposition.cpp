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

int GetSum(const std::vector<int>& v, const std::vector<int>& idx, int m) {

    int sum = 0;

    for (int i = 0; i < m; ++i) {
        sum += v[idx[i]];
    }

    return sum;
}

bool LowerBound(const std::vector<int> &v, int target,
                int partVal, int &idx, int low) {

    const int bound = target - partVal;

    if (v[idx] <= bound) {
        return false;
    } else if (v[low] < bound) {
        auto lower = std::find_if(
            v.cbegin() + low, v.cbegin() + idx, [=](int v_i) {
                return v_i >= bound;
            }
        );

        idx = std::distance(v.cbegin(), lower);
        return v[idx] > bound;
    } else {
        idx = low;
        return false;
    }
}

void LowerBoundLast(const std::vector<int> &v, int target,
                    int partVal, int &idx, int low) {

    const int bound = target - partVal;

    if (v[idx] > bound && v[low] < bound) {
        while (idx > low && v[idx] > bound) {
            --idx;
        }
    } else {
        idx = low;
    }
}

int GetLowerBound(const std::vector<int> &v, std::vector<int> &z,
                  int n, int m, int strt, int target) {

    int currPartial = 0;
    const int lastCol = m - 1;
    std::vector<int> vPass(m);
    vPass.assign(v.crbegin(), v.crbegin() + m);
    int partVal = std::accumulate(vPass.cbegin(), vPass.cbegin() + m - 1, 0);

    if (strt == 0) {
        const int testMax = partVal + vPass.back();

        if (testMax < target) {
            return 0;
        }
    }

    int currPos = n - m;

    if (strt) {
        for (int i = 0; i < strt; ++i) {
            vPass[i] = v[z[i]];
            partVal = partVal + vPass[i];
            ++currPos;
            partVal -= v[currPos];
        }

        currPartial = std::accumulate(
            vPass.cbegin(), vPass.cbegin() + strt, 0
        );

        for (int i = strt, j = 1; i < m; ++i, ++j) {
            vPass[i] = v[z[strt - 1] + j];
        }
    } else {
        vPass.assign(v.cbegin(), v.cbegin() + m);
    }

    const int testMin = std::accumulate(vPass.cbegin(), vPass.cbegin() + m, 0);

    if (testMin > target) {
        return 0;
    }

    int idx = n - m + strt;
    int lowBnd = (strt) ? z[strt - 1] + 1 : 0;

    for (int i = strt; i < lastCol; ++i) {
        if (LowerBound(v, target, partVal, idx, lowBnd)) {
            if (idx > lowBnd) {
                const int numIterLeft = m - i;

                for (int j = 0, k = idx; j < numIterLeft; ++j, ++k) {
                    vPass[j] = v[k];
                }

                const int minRemaining = std::accumulate(
                    vPass.cbegin(), vPass.cbegin() + numIterLeft, 0
                );
                const int currMin = minRemaining + currPartial;

                if (currMin > target) {
                    --idx;
                }
            }
        }

        z[i] = idx;
        partVal += v[idx];
        currPartial += v[idx];

        ++idx;
        ++currPos;

        lowBnd = idx;
        idx = currPos;
        partVal -= v[currPos];
    }

    LowerBoundLast(v, target, partVal, idx, lowBnd);
    z[lastCol] = idx;
    return 1;
}

void NextSection(const std::vector<int> &v, std::vector<int> &z,
                 bool &check, int &test, int m, int n, int target) {

    for (int i = m - 2, nMinusM = n - m; i >= 0 && !check; --i) {
        if (z[i] != (nMinusM + i)) {
            ++z[i];
            GetLowerBound(v, z, n, m, i + 1, target);
            test  = GetSum(v, z, m);
            check = test <= target;
        }
    }
}

int FilterProspects(const std::vector<int> &v, std::vector<int> &z,
                    bool &check, int &test, int m, int m1, int n, int target) {

    while (check) {
        if (test == target) {
            return 1;
        }

        check = z[m1] != (n - 1);

        if (check) {
            ++z[m1];
            test  = GetSum(v, z, m);
            check = test <= target;
        }
    }

    return 0;
}

int NextDistinctBlock(
    const std::vector<int> &v, std::vector<int> &idx, int target, int m
) {

    const int n = v.size();

    std::iota(idx.begin(), idx.end(), 0);
    GetLowerBound(v, idx, n, m, 0, target);

    int test   = GetSum(v, idx, m);
    int res    = test == target;
    bool check = test <= target;

    while (res == 0 && check) {
        res = FilterProspects(v, idx, check, test, m, m - 1, n, target);
        NextSection(v, idx, check, test, m, n, target);
    }

    return res;
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

    // We need a proof that complement is guaranteed to be sorted after
    // the swaps. Intuitively, it makes sense, however rigor is required.
    // Sorting is expensive based off of empirical tests. If we end up
    // proving that we need to sort, we can do something like the below:
    //
    // sort the small affected region [i1, i2 + 1]
    // if (i1 > i2) {
    //     std::sort(complement.begin() + i2, complement.begin() + i1 + 1);
    // } else {
    //     std::sort(complement.begin() + i1, complement.begin() + i2 + 1);
    // }

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
    std::vector<int> &z, std::vector<int> &complement, std::vector<int> &idx,
    int &i1, int &i2, int &myMax, int lastCol, int lastIdx, int target
) {

    if (z[lastCol - 1] < myMax) {
        NextRoutine(z, complement, i1, i2, lastCol, lastIdx);
    } else if (lastCol > 1) {
        int m = 2;
        int j = lastCol - m;

        complement.insert(
            std::lower_bound(
                complement.begin(), complement.end(), z[lastCol]
            ),
            z[lastCol]
        );
        complement.insert(
            std::lower_bound(
                complement.begin(), complement.end(), z[lastCol - 1]
            ),
            z[lastCol - 1]
        );

        bool keepGoing = true;

        while (j >= 0 && keepGoing) {
            int res = 0;

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
                    idx.resize(m);
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
