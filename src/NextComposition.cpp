#include <algorithm>
#include <numeric>
#include <vector>

std::vector<int> PrepareComplement(std::vector<int> z, int target) {

    const int z_size = z.size() - std::count(z.cbegin(), z.cend(), 0);

    // Here we are trying to find the maximum possible value of z. We do this
    // by assuming that the first n - 1 elements are minimized (n is the size
    // of z). That is, the first n - 1 elements are 1, 2, ..., n - 1. This
    // will give us a sum of (n * (n - 1) / 2). Finally, we take the target
    // and subtract the triangle number from it. E.g. if the size of z is 5
    // and the target is 30, the absolute maximal element will be:
    //
    //   30 - (1 + 2 + 3 + 4) = 30 - (5 * 4 / 2) = 30 - 10 = 20.
    //
    // This will be the starter for determining the elements that could make
    // up the complement.

    int myMax = target - static_cast<int>((z_size * (z_size - 1)) / 2);
    std::vector<int> complement;

    if (myMax > z_size) {
        std::sort(z.begin(), z.end());
        std::vector<int> myRange(myMax);
        std::iota(myRange.begin(), myRange.end(), 1);
        std::set_difference(myRange.begin(), myRange.end(), z.begin(), z.end(),
                            std::inserter(complement, complement.begin()));
    }

    return complement;
}

int GetMax(const std::vector<int> &z, const std::vector<int> &complement) {

    int res = z.back() > z[z.size() - 2] ? z.back() : z[z.size() - 2];
    const int last_two = std::accumulate(z.end() - 2, z.end(), 0);

    int target = last_two - complement.front();

    auto lower = std::lower_bound(
        complement.cbegin(), complement.cend(), target
    );

    bool foundTar = lower != complement.end() && *lower == target;
    const int comp_size = complement.size();

    // We've already checked the first element above with front()
    int j = 1;

    while (!foundTar && j < comp_size && target > res) {
        target = last_two - complement[j];

        lower = std::lower_bound(
            complement.cbegin(), complement.cend(), target
        );

        foundTar = lower != complement.end() && *lower == target;
        ++j;
    }

    if (foundTar) {
        --j;

        if (complement[j] > target && complement[j] > res) {
            res = complement[j];
        } else if (target > res) {
            res = target;
        }
    }

    return res;
}

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

int RangeSum(const std::vector<int>& v, int low, int high) {
    return std::accumulate(v.cbegin() + low, v.cbegin() + high, 0);
}

int FindBacktrackIndex(const std::vector<int>& idx, int k, int g) {

    // We are finding the first index that isn't maximized
    while (k > 0 && idx[k] == g) {
        --k;
        --g;
    }

    return k;
}

/* ---------------------------------------------------------------------------
 Overview of How These Routines Work Together
 ---------------------------------------------------------------------------

 These routines cooperate to find pairs of indices in the sorted vector "v"
 (the complement) whose values sum to a required target. For m = 2 this
 means finding distinct indices i < j with v[i] + v[j] == target.

 ---------------------------------------------------------------------------
 1. NextDistinctBlock(v, idx, tailSum, target, m)
 ---------------------------------------------------------------------------

 Purpose:
 Find the FIRST valid m-tuple of indices whose values sum to "target".
 The search is lexicographic, so for m = 2 this yields the pair with the
 smallest possible idx[0], and among those the smallest idx[1].

 Behavior:
 * Initializes idx to {0,1,...}, the smallest lexicographic combination.
 * Iteratively advances indices without ever decreasing them.
 * Uses pruning (via tailSum and partial sums) to eliminate branches that
 cannot contain ANY valid solution.
 * Because pruning only removes globally impossible branches, earlier
 possible pairs are never skipped.  Thus the first discovered solution
 is the lexicographically smallest valid one.
 * Returns 1 with idx filled on success, or 0/-1/-2 when impossible.

 Role:
 "Find the lexicographically earliest valid pair, or report none."

 ---------------------------------------------------------------------------
 2. NextDistinctBlock2(v, idx, target, maxLast)
 ---------------------------------------------------------------------------

 Purpose:
 Search for the FIRST valid pair satisfying BOTH:
 - v[idx[0]] + v[idx[1]] == target
 - v[idx[k]] < maxLast for both k = 0,1.

 Behavior:
 * Scans idx[0] upward and performs lower_bound on the remaining tail.
 * Returns the earliest lexicographic pair satisfying the maxLast filter.

 Role:
 "Find the earliest valid pair under the added maxLast constraint."

 ---------------------------------------------------------------------------
 3. CompsDistinctSetup(...)
 ---------------------------------------------------------------------------

 Purpose:
 High-level driver that coordinates all steps needed to choose the pair
 that interacts correctly with the current composition z.

 Workflow:
 1. Build complement of z and compute lastTwo = z[last] + z[last-1].
 2. Call NextDistinctBlock to obtain the FIRST possible pair.  If none
 exists, report idx_1 = -1 so the caller can restructure z.
 3. If a solution exists, call NextDistinctBlock2 to find the FIRST pair
 that also satisfies v[i] < z.back().  If found, use it.
 4. If the restricted search fails, fall back to a manual sweep:
 - Begin at the first v[idx_1] > z.back().
 - For each idx_1, try idx_2 from 0 upward.
 - Because a solution is guaranteed, this finds the FIRST pair with
 idx_1 > z.back().
 5. Compute myMax for later steps.

 Role:
 "Select the correct next complementary pair, preferring restricted
 lexicographically minimal pairs when possible."
 --------------------------------------------------------------------------- */

int NextDistinctBlock2(
    const std::vector<int> &v, std::vector<int> &idx, int target, int maxLast
) {

    const int n = v.size();
    if (n < 2) return 0;

    // Quick bounding check: smallest possible 2-sum is too large
    const int minSum = v[0] + v[1];
    if (minSum > target) return -1;

    // Largest possible 2-sum is too small → no solution
    const int maxSum = v[n - 2] + v[n - 1];
    if (maxSum < target) return -2;

    int i = 0;
    // For each v[i], we search for v[j] = target - v[i]
    int partial = target - v.front();

    // Scan i until v[i] hits the cutoff maxLast
    while (i < (n - 1) && v[i] < maxLast) {
        // Search in the strictly larger index range [i+1, end)
        auto lower = std::lower_bound(
            v.cbegin() + (i + 1), v.cend(), partial
        );

        // Check if we found exact match and the second term respects < maxLast
        if (lower != v.cend() && *lower == partial && partial < maxLast) {
            idx[0] = i;
            idx[1] = std::distance(v.cbegin(), lower);
            return 1; // Found valid pair
        }

        // Otherwise advance i and recompute the needed complement
        ++i;
        partial = target - v[i];
    }

    return 0; // No valid pair satisfying the maxLast constraint
}

int NextDistinctBlock(const std::vector<int> &v, std::vector<int> &idx,
                      std::vector<int> &tailSum, int target, int m) {

    int n = v.size();
    // Largest m-sum available from v (rightmost m values)
    const int testMax = std::accumulate(v.cend() - m, v.cend(), 0);

    // If even the largest m elements can't reach target → impossible
    if (testMax < target) {
        return -2;
    }

    // Smallest possible m-sum from v (first m values)
    const int testMin = std::accumulate(v.cbegin(), v.cbegin() + m, 0);

    // If minimum m-sum exceeds target → impossible
    if (testMin > target) {
        return -1;
    }

    // Initialize idx = {0,1,...,m-1}
    std::iota(idx.begin(), idx.end(), 0);

    // Sum of v[idx[0]] + ... + v[idx[m-1]]
    int tempSum = GetSum(v, idx, m);

    // partial = required value for the last index position
    int partial = target - (tempSum - v[idx.back()]);

    // tailSum[k] = sum of the last (k+1) values of v
    // Used to prune: if even the largest possible tail can't reach target
    tailSum.resize(m);
    std::partial_sum(v.crbegin(), v.crbegin() + m, tailSum.begin());

    // Main search loop: increment idx lexicographically until sum == target
    while (tempSum != target) {

        // Try to place the last index via binary search
        auto lower = std::lower_bound(
            v.cbegin() + idx.back(), v.cend(), partial
        );

        if (lower != v.cend() && *lower == partial) {
            // Found correct last component → full solution
            idx.back() = std::distance(v.cbegin(), lower);
            return 1;
        }

        // Need to backtrack to an earlier index position k
        // Find the earliest index we can increment
        int k = FindBacktrackIndex(idx, m - 2, n - 2);

        bool impossible = true;

        // front_partial = sum of v[idx[0..k]]
        int front_partial = GetSum(v, idx, k);

        // Attempt to find a valid backtracking point
        for (; impossible && k >= 0; --k) {

            // tempSum is the total constructed using:
            // - fixed prefix idx[0..k]
            // - the next (m - k) values chosen as the smallest possible after idx[k]
            tempSum = front_partial +
                RangeSum(v, idx[k] + 1, idx[k] + m - k + 1);

            // Pruning conditions:
            // 1) tempSum > target → overshoot
            // 2) even using the largest tail (tailSum) we cannot reach target
            impossible = (tempSum > target) ||
                ((front_partial + tailSum[m - k - 1]) < target);

            // Prepare next prefix sum for k-1
            if (k > 0) {
                front_partial -= v[idx[k - 1]];
            }
        }

        if (impossible) {
            // No valid backtracking point → no solution in this block
            return 0;
        }

        // Rebuild the tail of the index vector lexicographically
        std::iota(idx.begin() + k + 1, idx.end(), idx[k + 1] + 1);

        // Recompute new required last value
        partial = target - (tempSum - v[idx.back()]);
    }

    return 1; // Found exact m-sum
}

int CompsDistinctSetup(
    const std::vector<int> &z, std::vector<int> &complement,
    int &tar, int &idx_1, int &idx_2, int &myMax
) {

    // tar = total sum of z
    tar = std::accumulate(z.cbegin(), z.cend(), 0);

    // Build complement = {0..tar} \ z, sorted
    complement = PrepareComplement(z, tar);

    if (complement.empty()) {
        return 0;
    }

    // Sum of last two elements of z; target for 2-sum search
    int lastTwo = std::accumulate(z.rbegin(), z.rbegin() + 2, 0);

    std::vector<int> idx(2);
    std::vector<int> tailSum;

    // First: check if ANY 2-sum exists in complement equal to lastTwo
    int sol_exist = NextDistinctBlock(
        complement, idx, tailSum, lastTwo, 2
    );

    if (sol_exist == 1) {

        // Try finding the *best* solution — one respecting maxLast = z.back()
        int best_sol = NextDistinctBlock2(
            complement, idx, lastTwo, z.back()
        );

        if (best_sol == 1) {
            // Found an optimal pair directly
            idx_1 = idx.front();
            idx_2 = idx.back();
        } else {
            // Need fallback: manually find first idx1 > z.back()

            auto it = std::upper_bound(
                complement.cbegin(), complement.cend(), z.back()
            );

            idx_1 = std::distance(complement.cbegin(), it);
            idx_2 = 0;

            int testSum = complement[idx_1] + complement[idx_2];

            // Monotone two-pointer search with guaranteed existence of solution
            while (testSum != lastTwo) {

                // Increase idx_2 while sum is too small
                while (testSum <= lastTwo && idx_2 < idx_1) {
                    ++idx_2;
                    testSum = complement[idx_1] + complement[idx_2];
                }

                // If not solved, advance idx_1 and reset idx_2
                if (testSum != lastTwo) {
                    ++idx_1;
                    idx_2 = 0;
                    testSum = complement[idx_1] + complement[idx_2];
                }
            }
        }

    } else {
        // No further 2-sum solutions exist in complement.
        // Next iteration will swap last two entries of z,
        // causing a new block to be computed on next call.

        idx_1 = -1;

        // For idx_2 choose the largest complement value < z.back()
        auto it_last = std::lower_bound(
            complement.cbegin(), complement.cend(), z.back()
        );

        idx_2 = (it_last != complement.cend()) ?
        std::distance(complement.cbegin(), --it_last) :
            complement.size() - 1;
    }

    // myMax = maximum allowed value for next step (z-dependent)
    myMax = GetMax(z, complement);
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

    // We need a proof that complement is guaranteed to be sorted after
    // the swaps. Intuitively, it makes sense, however rigor is required.
    // Sorting is expensive based off of empirical tests. Until we have
    // a proof, we will sort the affected area.
    //
    // Sort the small affected region [i1, i2 + 1]
    if (i1 > i2) {
        std::sort(complement.begin() + i2, complement.begin() + i1 + 1);
    } else {
        std::sort(complement.begin() + i1, complement.begin() + i2 + 1);
    }

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
    std::vector<int> &tailSum, int &i1, int &i2, int &myMax, int lastCol,
    int lastIdx, int target
) {

    if (z[lastCol - 1] < myMax) {
        NextRoutine(z, complement, i1, i2, lastCol, lastIdx);
    } else if (lastCol > 1) {
        int m = 2;
        int j = lastCol - m;

        for (int i = lastCol - 1; i <= lastCol; ++i) {
            complement.insert(
                std::lower_bound(complement.begin(), complement.end(), z[i]),
                z[i]
            );
        }

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
                    res = NextDistinctBlock(
                        complement, idx, tailSum, partial, m
                    );
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

                int res2 = NextDistinctBlock2(
                    complement, idx, lastTwo, z.back()
                );

                if (res2 == 1) {
                    i1 = idx.front();
                    i2 = idx.back();
                } else {
                    // This means that there are no further solutions in
                    // complement. The next iteration we simply will swap the
                    // last two elements of z. By doing this, we will also
                    // ensure that z[lastCol - 1] will be equal to myMax which
                    // will force a new block calculation.
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
