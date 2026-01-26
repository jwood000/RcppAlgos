#include <algorithm> // std::set_difference; std::count; std::min; std::sort;
                     // std::accumulate; std::lower_bound; std::partial_sum
#include <iterator>  // std::inserter; std::distance
#include <numeric>   // std::iota
#include <cstdint>
#include <vector>

static inline int minTailSum(int t) {
    // 1 + 2 + ... + t
    return (t <= 0) ? 0 : (t * (t + 1)) / 2;
}

// Returns true iff the vector is "maximized at each non-zero position"
// under: fixed total sum, positive parts, strictly decreasing after first
// non-zero, and allowing leading zeros only.
bool IsMaximizedGreedySuffix(const std::vector<int>& v, int target, int nz) {

    const int n = static_cast<int>(v.size());

    // Now check greedy-max condition at each position. We guarantee nz is
    // within bounds in the calling code.
    int prefixSum = v[nz];

    for (int i = nz + 1; i < n; ++i) {
        const int t = n - i - 1;             // slots after i
        const int R = target - prefixSum;    // sum remaining INCLUDING a[i..end]
        const int reserve = minTailSum(t);   // must leave at least this for tail
        const int upper = v[i - 1] - 1;      // strict decrease constraint

        // Candidate from sum constraint:
        // v[i] <= R - reserve (so that tail can be at least 1..t)
        const int cand = R - reserve;
        const int expected = std::min(upper, cand);

        if (v[i] != expected) return false;
        prefixSum += v[i];
    }

    return true;
}

void BinaryNextElem(
    int &uppBnd, int &lowBnd, int &ind, int lastElem,
    int target, int partial, const std::vector<int> &v
) {

    int dist = target - (partial + v[ind]);

    while ((uppBnd - lowBnd) > 1 && dist != 0) {
        const int mid = (uppBnd - lowBnd) / 2;
        ind = lowBnd + mid;
        dist = target - (partial + v[ind]);

        if (dist > 0) {
            lowBnd = ind;
        } else {
            uppBnd = ind;
        }
    }

    if (dist < 0) {
        ind = lowBnd;
        dist = target - (partial + v[ind]);
    }

    if (dist > 0 && ind < lastElem) {
        ++ind;
    }
}

int GetFirstPartitionDistinct(const std::vector<int> &v, std::vector<int> &z,
                              int target, int m, int lenV) {

    // Compute maximum possible sum: largest m elements
    int testMax = 0;

    for (int i = lenV - m; i < lenV; ++i) {
        testMax += v[i];
    }

    if (testMax < target) {
        return -2; // length too small
    }

    // Compute minimum possible sum: smallest m elements
    int testMin = 0;

    for (int i = 0; i < m; ++i) {
        testMin += v[i];
    }

    if (testMin > target) {
        return -1; // length too large
    }

    const int lastElem = lenV - 1;

    // Initial positions for the first element
    int currPos = lenV - m;
    int mid = currPos / 2;

    int partial = testMax - v[currPos];
    int dist = target - (partial + v[mid]);

    int lowBnd = (dist > 0) ? mid : 0;
    int uppBnd = (dist > 0) ? currPos : mid;
    int ind = mid;

    // Fill first m-1 columns
    for (int i = 0; i < m - 1; ++i) {

        BinaryNextElem(uppBnd, lowBnd, ind,
                       lastElem, target, partial, v);

        z[i] = ind;
        partial += v[ind];

        // Distinct-case advancement
        ++ind;
        ++currPos;

        lowBnd = ind;
        uppBnd = currPos;

        int span = uppBnd - lowBnd;
        mid = span / 2;

        ind = lowBnd + mid;
        partial -= v[currPos];
    }

    // Last column
    BinaryNextElem(uppBnd, lowBnd, ind,
                   lastElem, target, partial, v);

    z[m - 1] = ind;

    // Final verification
    int finalCheck = 0;

    for (int i = 0; i < m; ++i) {
        finalCheck += v[z[i]];
    }

    if (finalCheck != target) {
        return 0;
    }

    return 1;
}

bool IsComplementZeroBased(bool firstZero, bool isWeak, bool IsGen) {

    // N.B. Only mapped cases with `v[0] == 0` are considered.
    // If zero appears elsewhere in `v` (e.g., due to a shift placing zero in
    // the interior), the intended semantics are unclear and such cases are
    // deliberately excluded.

    if (isWeak) {
        return true;
    }

    if (firstZero && !isWeak) {
        // Non-weak case with v containing a leading zero.
        //
        // When execution reaches this point, z already contains exactly the
        // elements needed for the remainder of the algorithm. There are two
        // possibilities:
        //
        // 1) z contains zeros:
        //    Those zeros map to leading zeros in every result. In the non-weak
        //    setting, we do not allow zero to appear in the composition, so
        //    zero must not appear in the complement.
        //
        // 2) z contains no zeros:
        //    We have reached the region where no mapped zeros appear in any
        //    composition. To preserve this in the non-weak case, zero must
        //    again be excluded from the complement.
        //
        // In both cases, zeros in the complement are disallowed.
        return false;
    }

    if (!isWeak && IsGen) {
        // In the non-weak setting, once zeros are excluded from v, we must
        // allow zero in the complement even if the current indexing vector z
        // does not contain it. Although some intermediate results may be
        // formed without zero (i.e. z may omit zero at this stage), later
        // lexicographic results may require it. For example, the final
        // lexicographic result corresponds to a reversal of the initial
        // indexing, implying z.back() -> 0.
        return true;
    }

    return false;
}

std::vector<int> PrepareComplement(std::vector<int> z, int target, int idx_max,
                                   bool startAtZero, int zeroBudget) {

    int z_size = static_cast<int>(z.size()) - zeroBudget;
    if (startAtZero) z_size = std::max(0, z_size - 1);

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

    const int tri = (z_size <= 1) ? 0 : (z_size * (z_size - 1)) / 2;
    const int myMax = std::min(target - tri, idx_max);
    std::vector<int> complement;

    if (myMax > z_size) {
        std::sort(z.begin(), z.end());
        std::vector<int> myRange(myMax + static_cast<int>(startAtZero));
        std::iota(myRange.begin(), myRange.end(),
                  startAtZero ? 0 : 1);

        std::set_difference(
            myRange.begin(), myRange.end(), z.begin(), z.end(),
            std::inserter(complement, complement.begin())
        );

        const int z_zeros = std::count(z.cbegin(), z.cend(), 0);
        const int cmp_zeros = std::count(
            complement.cbegin(), complement.cend(), 0
        );
        const int num_zeros_needed = zeroBudget - cmp_zeros - z_zeros;

        if (num_zeros_needed > 0) {
            complement.insert(complement.begin(), num_zeros_needed, 0);
        }
    }

    return complement;
}

// Contract:
//   - z.size() >= 2
//   - complement is non-empty
// These invariants are guaranteed by the calling code.
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

//  ---------------------------------------------------------------------------
//  Overview of How These Routines Work Together
//  ---------------------------------------------------------------------------
//
//  These routines select a pair of indices in the sorted vector "v"
//  (the complement) whose values sum to a required target. For m = 2 this
//  means finding distinct indices i != j with v[i] + v[j] == target.
//
//  ---------------------------------------------------------------------------
//  1. NextDistinctBlock2(v, idx, target, maxLast)
//  ---------------------------------------------------------------------------
//
//  Purpose:
//  Search for the FIRST valid pair satisfying BOTH:
//  - v[idx[0]] + v[idx[1]] == target
//  - v[idx[0]] < maxLast and v[idx[1]] < maxLast
//
//  Behavior:
//  * Considers candidates strictly less than maxLast.
//  * For each candidate, uses lower_bound to search for the complementary
//    value needed to reach target.
//  * Returns 1 with idx filled on success, or 0 when no such constrained
//    pair exists (and -1/-2 for quick global impossibility checks).
//
//  Role:
//  "Find the earliest 2-sum pair under the maxLast constraint."
//
//  ---------------------------------------------------------------------------
//  2. CompsDistinctSetup(...)
//  ---------------------------------------------------------------------------
//
//  Purpose:
//  Driver that builds the complement of z and attempts to select the next
//  complementary pair consistent with the current composition state. This
//  routine only selects pairs with both values < z.back(). If none exist, it
//  signals idx_1 = -1 and the caller must restructure z. It does not attempt
//  an alternate pair selection within the same block.
//
//  Workflow:
//  1. Build complement of z and compute lastTwo = z[last] + z[last-1].
//
//  2. Call NextDistinctBlock2 to find the earliest 2-sum pair with both
//     values < z.back().
//
//  3. If found, return its indices in idx_1/idx_2.
//     Otherwise set idx_1 = -1 so the caller can restructure z.
//
//  4. Compute myMax for later steps.
//
//  Role:
//  "Select the next complementary pair when available under the z.back()
//   constraint; otherwise signal the caller to restructure."
//  ---------------------------------------------------------------------------

int NextDistinctBlock2(const std::vector<int> &v, std::vector<int> &idx,
                       int target, int maxLast) {

    const int n = v.size();
    if (n < 2) return 0;

    // Quick bounding check: smallest possible 2-sum is too large
    const int minSum = v[0] + v[1];
    if (minSum > target) return -1;

    // Largest possible 2-sum is too small → no solution
    const int maxSum = v[n - 2] + v[n - 1];
    if (maxSum < target) return -2;

    // Get iterator to the first element >= maxLast
    auto upper = std::lower_bound(v.cbegin(), v.cend(), maxLast);

    // Move back to get largest element < maxLast (if any)
    while (upper != v.cbegin()) {
        --upper;

        int partial = target - *upper;

        // ---- Early exit: partial too large to find ----
        if (partial > v[n - 1]) break;

        auto lower = std::lower_bound(v.cbegin(), v.cend(), partial);

        if (lower != v.cend() && *lower == partial && lower != upper) {
            idx[0] = std::distance(v.cbegin(), lower);
            idx[1] = std::distance(v.cbegin(), upper);
            return 1; // Found valid pair
        }
    }

    return 0; // No valid pair satisfying the maxLast constraint
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

    // Initialize idx = {0, 1, ..., m - 1}
    std::iota(idx.begin(), idx.end(), 0);

    // Sum of v[idx[0]] + ... + v[idx[m - 1]]
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

        // Start backtracking from the rightmost pivot position (m - 2).
        // The last index (m - 1) is handled separately via binary search,
        // so the earliest position that can be incremented is m - 2.
        // Earlier versions attempted to skip already-maximized suffixes,
        // but in practice the pruning logic below handles those cases,
        // so starting from m - 2 is always correct and simpler.
        int k = m - 2;

        bool impossible = true;

        // front_partial = sum of v[idx[0..(k - 1)]]
        int front_partial = GetSum(v, idx, k);

        // Attempt to find a valid backtracking point
        for (; impossible && k >= 0; --k) {

            // tempSum is the minimal total achievable after deciding to
            // backtrack at k:
            //
            // - The prefix idx[0..k-1] is fixed.
            // - idx[k] itself has already been tested and cannot lead to a
            //   solution, so the tail must start strictly after it.
            // - Therefore we fill the remaining (m - k) positions with the
            //   smallest possible values starting at idx[k] + 1.
            //
            // The assignment to tempSum below computes the smallest sum
            // attainable for this backtrack which is used to decide whether
            // further search at this k is viable.

            const int low  = idx[k] + 1;
            const int high = idx[k] + m - k + 1;

            // if we can't fit the minimal tail, this k is invalid => force
            // impossible.
            if (high > n) {
                impossible = true;
            } else {
                tempSum = front_partial + RangeSum(v, low, high);

                // Pruning conditions:
                // 1) tempSum > target → overshoot
                // 2) Even using the largest tail (from tailSum) we cannot
                //    reach target.
                impossible = (tempSum > target) ||
                    ((front_partial + tailSum[m - k - 1]) < target);
            }

            // Prepare next prefix sum for k - 1
            if (k > 0) {
                front_partial -= v[idx[k - 1]];
            }
        }

        if (impossible) {
            // No valid backtracking point → no solution in this block
            return 0;
        }

        // If we find a viable pivot at K, the for-loop will exit after
        // decrementing k, so on exit k == K - 1 and pivot == k + 1.
        const int pivot = k + 1;

        // Rebuild the tail lexicographically by advancing idx[pivot] and
        // resetting idx[pivot+1..] to the minimal increasing continuation.
        std::iota(idx.begin() + pivot, idx.end(), idx[pivot] + 1);


        // Recompute new required last value
        partial = target - (tempSum - v[idx.back()]);
    }

    return 1; // Found exact m-sum
}

int CompsDistinctSetup(
    const std::vector<int> &z, std::vector<int> &complement,
    int &tar, int &idx_1, int &idx_2, int &myMax, int idx_max,
    bool startAtZero, int zeroBudget
) {

    // tar = total sum of z
    tar = std::accumulate(z.cbegin(), z.cend(), 0);

    // Build complement = {0..tar} \ z, sorted
    complement = PrepareComplement(z, tar, idx_max, startAtZero, zeroBudget);

    if (complement.empty()) {
        return 0;
    }

    // Sum of last two elements of z; target for 2-sum search
    int lastTwo = std::accumulate(z.rbegin(), z.rbegin() + 2, 0);

    std::vector<int> idx(2);
    std::vector<int> tailSum;

    int best_sol = NextDistinctBlock2(
        complement, idx, lastTwo, z.back()
    );

    if (best_sol == 1) {
        // Found an optimal pair directly
        idx_1 = idx.front();
        idx_2 = idx.back();
    } else {
        // No 2-sum exists satisfying the < maxLast constraint. Next iteration
        // will swap last two entries of z, causing a new block to be computed
        // on the next call.

        idx_1 = -1;

        // For idx_2 choose the largest complement value < z.back()
        auto it_last = std::lower_bound(
            complement.cbegin(), complement.cend(), z.back()
        );

        idx_2 = (it_last != complement.cbegin()) ?
            std::distance(complement.cbegin(), --it_last) :
                complement.size() - 1;
    }

    // myMax = maximum allowed value for next step (z-dependent)
    myMax = GetMax(z, complement);
    return 1;
}
