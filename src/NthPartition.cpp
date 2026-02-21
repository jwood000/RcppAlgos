#include "cpp11/R.hpp"
#include "cpp11/protect.hpp"

#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsUtils.h"
#include "Partitions/NthPartition.h"
#include "CppConvert/Constants.h"
#include <numeric>  // std::accumulate

std::vector<int> nthLengthOne(int n, int m, int cap, int k,
                              double dblIdx, const mpz_class &mpzIdx) {

    std::vector<int> res(m, n - 1);
    return res;
}

//*********************** Compositions Functions **************************//

std::vector<int> nthCompsRep(int n, int m, int cap, int k,
                             double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0) {
        for (double temp = CountCompsRepLen(n, m); temp <= dblIdx; ++j) {
            --n;
            dblIdx -= temp;
            temp = CountCompsRepLen(n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthCompsRepZero(int n, int m, int cap, int k,
                                 double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, j = incr_j) {
        double temp = incr_j ? CountCompsRepLen(n, m) :
                    CountCompsRepZNotWk(n, m);

        for (; temp <= dblIdx; ++j) {
            incr_j = true;
            --n;
            dblIdx -= temp;
            temp = CountCompsRepLen(n, m);
        }

        if (incr_j) --n;
        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

// nthCompsDistinct
// ----------------
// Unranks the dblIdx-th distinct composition of n into width parts,
// using positive integers, with each part distinct, and bounded by cap.
//
// Output ordering:
//   - Lexicographic order on the resulting vector:
//
//                   (v[0], v[1], ...,v[width - 1])
//
//   - Because parts are distinct, each value can appear at most once.
//
// High-level strategy (classic unranking by "block subtraction"):
//
//   At each position i (left to right), we consider fixing the next part to
//   each feasible candidate value in increasing order. For each candidate
//   choice, we count how many full compositions exist with that fixed prefix.
//   Those counts form contiguous "blocks" in lexicographic order.
//
//   If dblIdx is larger than the size of the first block, we subtract the
//   block size and advance to the next candidate. Repeat until dblIdx falls
//   inside the current candidate's block. That candidate is the correct value
//   for position i.
//
// Key data structures:
//
//   - mask[v]:
//       1 if value v has been used in the prefix, else 0. Used to enforce
//       distinctness and to find the next unused candidate fast.
//
//   - allowed:
//       A compact list of remaining feasible values for the suffix.
//       It is continuously rebuilt (bounded) via UpdateAllowed so the DP
//       counting stays small and fast.
//
//   - partial_sum:
//       Sum of the fixed prefix so far (including the current chosen value).
//
//   - m (mutated inside the function):
//       Number of parts remaining to fill *after* the current position.
//       This is why the code does --m up front and then --m each iteration.
//
// Important feasibility bound:
//
//   The smallest possible sum of width distinct positive integers is
//     1 + 2 + ... + width = width * (width + 1) / 2.
//   If we fix a prefix, the minimal possible suffix is the sum of the smallest
//   unused values. UpdateAllowed uses this to restrict candidate values.
//
// Notes about k/mpzIdx:
//   In this double-based version, mpzIdx is unused; it exists to mirror
//   a big-integer path elsewhere.
//
std::vector<int> nthCompsDistinct(int n, int m, int cap, int k,
                                  double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_val = std::min(cap, n - (width * (width - 1)) / 2);

    std::vector<char> mask(max_val + 1, 0);
    const int mask_size = mask.size();
    std::vector<int> res(width, 0);
    --m;

    mask[1] = 1;
    std::vector<int> allowed(max_val - 1);
    std::iota(allowed.begin(), allowed.end(), 2);

    int partial_sum = 1;
    int cur_val = 1;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        double temp = CountCompDistLenRstrctd(n - partial_sum, m, allowed);

        for (; temp <= dblIdx; cur_val = j) {
            while ((j + 1) < mask_size && mask[j + 1]) {
                ++j;
            }

            partial_sum += (j + 1 - cur_val);
            UpdateAllowed(mask, allowed, i, j + 1, width,
                          n, cur_val, partial_sum);

            dblIdx -= temp;
            temp = CountCompDistLenRstrctd(n - partial_sum, m, allowed);

            if (temp <= dblIdx) {
                ++j;
            }
        }

        res[i] = j;
        j = 0;

        while ((j + 1) < mask_size && mask[j + 1]) {
            ++j;
        }

        cur_val = j + 1;
        partial_sum += (j + 1);
        UpdateAllowed(mask, allowed, i + 1, j + 1, width,
                      n, cur_val, partial_sum);
    }

    res[width - 1] = n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthCompsDistinctMZ(int n, int m, int cap, int k,
                                    double dblIdx, const mpz_class &mpzIdx) {

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    double temp = (cap == n) ?
        CountCompsDistinctLen(n, k) :
        CountCompDistLenRstrctd(n, k, allowed);

    while (dblIdx >= temp && k < m) {
        dblIdx -= temp;
        ++k;
        temp = (cap == n) ?
            CountCompsDistinctLen(n, k) :
            CountCompDistLenRstrctd(n, k, allowed);
    }

    std::vector<int> res = nthCompsDistinct(n, k, cap, k, dblIdx, mpzIdx);

    for (auto& z_i: res) {
        ++z_i;
    }

    if (m > k) res.insert(res.begin(), m - k, 0);
    return res;
}

std::vector<int> nthCompsDistinctMZWeak(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_val = std::min(cap, n - (k * (k - 1)) / 2);
    int zeros_remaining = width - k;

    std::vector<char> mask(max_val + 1, 0);
    std::vector<int> res(width, 0);

    --zeros_remaining;
    --m;

    mask[0] = 1;
    std::vector<int> allowed(max_val);
    std::iota(allowed.begin(), allowed.end(), 1);

    int partial_sum = 0;
    int cur_val = 0;

    bool anyZeros = true;
    const int mask_size = mask.size();

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        int strtLen = std::max(1, m - zeros_remaining);
        if (anyZeros) ++zeros_remaining;

        double temp = CountCompsDistinctRstrctdMZWeak(
            n - partial_sum, m, allowed, strtLen
        );

        for (; temp <= dblIdx && n > partial_sum; cur_val = j) {
            ++j;

            while (j < mask_size && mask[j]) {
                ++j;
            }

            partial_sum += (j - cur_val);
            UpdateAllowed(mask, allowed, i, j, width - zeros_remaining,
                          n, cur_val, partial_sum);

            dblIdx -= temp;
            strtLen = std::max(1, m - zeros_remaining);
            temp = CountCompsDistinctRstrctdMZWeak(
                n - partial_sum, m, allowed, strtLen
            );
        }

        if (j == 0) --zeros_remaining;
        res[i] = j;

        anyZeros = zeros_remaining > 0;
        j = anyZeros ? 0 : 1;

        while (!anyZeros && j < mask_size && mask[j]) {
            ++j;
        }

        if (anyZeros) --zeros_remaining;
        cur_val = j;
        partial_sum += j;

        UpdateAllowed(mask, allowed, i + 1, j, width - zeros_remaining,
                      n, cur_val, partial_sum);
    }

    res[width - 1] = n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

//************************* Partition Functions ***************************//

std::vector<int> nthPartsRepLen(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        for (double temp = CountPartsRepLen(n, m); temp <= dblIdx; ++j) {
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsRepLen(n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsRepShort(int n, int m, int cap, int k,
                                  double dblIdx, const mpz_class &mpzIdx) {

    return nthPartsRepLen(n, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRep(int n, int m, int cap, int k,
                             double dblIdx, const mpz_class &mpzIdx) {

    return nthPartsRepLen(n * 2, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRepCap(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        for (double temp = CountPartsRepLenRstrctd(n, m, allowed);
             temp <= dblIdx; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            dblIdx -= temp;
            temp = CountPartsRepLenRstrctd(n, m, allowed);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctLen(int n, int m, int cap, int k,
                                     double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        for (double temp = CountPartsDistinctLen(n, m); temp <= dblIdx; ++j) {
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsDistinctLen(n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctMultiZero(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    std::vector<int> empty_allowed;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        double temp = (incr_j || i >= (width - k)) ?
                      CountPartsDistinctLen(n, m) :
                      CountPartsDistinctMultiZero(n, m, empty_allowed, k);

        for (; temp <= dblIdx; ++j) {
            incr_j = true;
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsDistinctLen(n, m);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctCap(int n, int m, int cap, int k,
                                     double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --cap;
    --m;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        for (double temp = CountPartsDistLenRstrctd(n, m, allowed);
             temp <= dblIdx; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            dblIdx -= temp;
            temp = CountPartsDistLenRstrctd(n, m, allowed);
        }

        res[i] = j;
        allowed.pop_back();
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctCapMZ(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        double temp = (incr_j || i >= (width - k)) ?
                CountPartsDistLenRstrctd(n, m, allowed) :
                CountPartsDistinctRstrctdMZ(n, m, allowed, k);

        for (; temp <= dblIdx; ++j) {
            incr_j = true;
            n -= (m + 1);
            allowed.pop_back();
            dblIdx -= temp;
            temp = CountPartsDistLenRstrctd(n, m, allowed);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            allowed.pop_back();
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

//*********************** Starting Gmp Functions **************************//

std::vector<int> nthCompsRepGmp(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    std::unique_ptr<CountClass> Counter = MakeCount(
        PartitionType::CompRepNoZero
    );

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0) {
        Counter->GetCount(temp, n, m);

        for (; cmp(temp, index) <= 0; ++j) {
            --n;
            index -= temp;
            Counter->GetCount(temp, n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthCompsRepZeroGmp(int n, int m, int cap, int k,
                                    double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);
    std::unique_ptr<CountClass> Counter = MakeCount(
        PartitionType::CmpRpZroNotWk
    );

    std::vector<int> empty_allowed;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, j = incr_j) {
        Counter->GetCount(temp, n, m, empty_allowed, k, !incr_j);

        for (; cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            --n;
            index -= temp;
            Counter->GetCount(temp, n, m, empty_allowed, k, false);
        }

        temp = 0;
        if (incr_j) --n;
        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthCompsDistinctGmp(int n, int m, int cap, int k,
                                     double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_val = std::min(cap, n - (width * (width - 1)) / 2);

    std::vector<char> mask(max_val + 1, 0);
    const int mask_size = mask.size();
    std::vector<int> res(width, 0);
    --m;

    mask[1] = 1;
    std::vector<int> allowed(max_val - 1);
    std::iota(allowed.begin(), allowed.end(), 2);

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::PrmDstPrtCap;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    int partial_sum = 1;
    int cur_val = 1;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        Counter->GetCount(temp, n - partial_sum, m, allowed);

        for (; cmp(temp, index) <= 0; cur_val = j) {
            while ((j + 1) < mask_size && mask[j + 1]) {
                ++j;
            }

            partial_sum += (j + 1 - cur_val);
            UpdateAllowed(mask, allowed, i, j + 1, width,
                          n, cur_val, partial_sum);

            index -= temp;
            Counter->GetCount(temp, n - partial_sum, m, allowed);

            if (cmp(temp, index) <= 0) {
                ++j;
            }
        }

        res[i] = j;
        j = 0;

        while ((j + 1) < mask_size && mask[j + 1]) {
            ++j;
        }

        cur_val = j + 1;
        partial_sum += (j + 1);
        UpdateAllowed(mask, allowed, i + 1, j + 1, width,
                      n, cur_val, partial_sum);
    }

    res[width - 1] = n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthCompsDistinctMZGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    mpz_class temp;
    mpz_class index(mpzIdx);

    if (cap == n) {
        const PartitionType ptype = PartitionType::CmpDstctNoZero;
        std::unique_ptr<CountClass> Counter = MakeCount(ptype);

        Counter->SetArrSize(ptype, n, m);
        Counter->InitializeMpz();
        Counter->GetCount(temp, n, k);

        while (cmp(index, temp) >= 0 && k < m) {
            index -= temp;
            ++k;
            Counter->GetCount(temp, n, k);
        }
    } else {
        std::vector<int> allowed(cap);
        std::iota(allowed.begin(), allowed.end(), 1);

        const PartitionType ptype = PartitionType::CmpDstctCapped;
        std::unique_ptr<CountClass> Counter = MakeCount(ptype);

        Counter->SetArrSize(ptype, n, m);
        Counter->InitializeMpz();
        Counter->GetCount(temp, n, k, allowed);

        while (cmp(index, temp) >= 0 && k < m) {
            index -= temp;
            ++k;
            Counter->GetCount(temp, n, k, allowed);
        }
    }

    std::vector<int> res;

    if (cmp(index, Significand53) > 0) {
        res = nthCompsDistinctGmp(n, k, cap, k, dblIdx, index);
    } else {
        dblIdx = index.get_d();
        res = nthCompsDistinct(n, k, cap, k, dblIdx, index);
    }

    for (auto& z_i: res) {
        ++z_i;
    }

    if (m > k) res.insert(res.begin(), m - k, 0);
    return res;
}

std::vector<int> nthCompsDistinctMZWeakGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_val = std::min(cap, n - (k * (k - 1)) / 2);
    int zeros_remaining = width - k;

    std::vector<char> mask(max_val + 1, 0);
    std::vector<int> res(width, 0);

    --zeros_remaining;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    mask[0] = 1;
    std::vector<int> allowed(max_val);
    std::iota(allowed.begin(), allowed.end(), 1);

    int partial_sum = 0;
    int cur_val = 0;

    bool anyZeros = true;
    const int mask_size = mask.size();

    const PartitionType ptype = PartitionType::CmpDstCapMZWeak;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        int strtLen = std::max(1, m - zeros_remaining);
        if (anyZeros) ++zeros_remaining;

        Counter->GetCount(
            temp, n - partial_sum, m, allowed, strtLen, true
        );

        for (; cmp(temp, index) <= 0 && n > partial_sum; cur_val = j) {
            ++j;

            while (j < mask_size && mask[j]) {
                ++j;
            }

            partial_sum += (j - cur_val);
            UpdateAllowed(mask, allowed, i, j, width - zeros_remaining,
                          n, cur_val, partial_sum);

            index -= temp;
            strtLen = std::max(1, m - zeros_remaining);
            Counter->GetCount(
                temp, n - partial_sum, m, allowed, strtLen, true
            );
        }

        if (j == 0) --zeros_remaining;
        res[i] = j;

        anyZeros = zeros_remaining > 0;
        j = anyZeros ? 0 : 1;

        while (!anyZeros && j < mask_size && mask[j]) {
            ++j;
        }

        if (anyZeros) --zeros_remaining;
        cur_val = j;
        partial_sum += j;

        UpdateAllowed(mask, allowed, i + 1, j, width - zeros_remaining,
                      n, cur_val, partial_sum);
    }

    res[width - 1] = n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsRepLenGmp(int n, int m, int cap, int k,
                                   double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::RepShort;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        Counter->GetCount(temp, n, m);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            index -= temp;
            Counter->GetCount(temp, n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsRepShortGmp(int n, int m, int cap, int k,
                                     double dblIdx, const mpz_class &mpzIdx) {

    return nthPartsRepLenGmp(n, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRepGmp(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    return nthPartsRepLenGmp(n * 2, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRepCapGmp(int n, int m, int cap, int k,
                                   double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::RepCapped;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        Counter->GetCount(temp, n, m, allowed);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            index -= temp;
            Counter->GetCount(temp, n, m, allowed);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctLenGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DstctNoZero;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        Counter->GetCount(temp, n, m);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            index -= temp;
            Counter->GetCount(temp, n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctMultiZeroGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DstctMultiZero;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> empty_allowed;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        Counter->GetCount(temp, n, m, empty_allowed, k, bLiteral);

        for (; cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            n -= (m + 1);
            index -= temp;
            Counter->GetCount(temp, n, m, empty_allowed, k, false);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctCapGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --cap;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DstctCapped;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        Counter->GetCount(temp, n, m, allowed, k);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            allowed.pop_back();
            index -= temp;
            Counter->GetCount(temp, n, m, allowed, k);
        }

        res[i] = j;
        allowed.pop_back();
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctCapMZGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DstctCappedMZ;
    std::unique_ptr<CountClass> Counter = MakeCount(ptype);

    Counter->SetArrSize(ptype, n, m);
    Counter->InitializeMpz();

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        Counter->GetCount(temp, n, m, allowed, k, bLiteral);

        for (; cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            n -= (m + 1);
            allowed.pop_back();
            index -= temp;
            Counter->GetCount(temp, n, m, allowed, k, false);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            allowed.pop_back();
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> EmptyReturn(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {
    std::vector<int> res(m);
    return(res);
}

nthPartsPtr GetNthPartsFunc(PartitionType ptype, bool IsGmp) {

    if (IsGmp) {
        switch (ptype) {
            case PartitionType::LengthOne: {
                return(nthPartsPtr(nthLengthOne));
            } case PartitionType::DstctCapped: {
                return(nthPartsPtr(nthPartsDistinctCapGmp));
            } case PartitionType::DstctCappedMZ: {
                return(nthPartsPtr(nthPartsDistinctCapMZGmp));
            } case PartitionType::DstctNoZero: {
                return(nthPartsPtr(nthPartsDistinctLenGmp));
            } case PartitionType::DstctOneZero: {
                return(nthPartsPtr(nthPartsDistinctLenGmp));
            } case PartitionType::DstctMultiZero: {
                return(nthPartsPtr(nthPartsDistinctMultiZeroGmp));
            } case PartitionType::DstctStdAll: {
                return(nthPartsPtr(nthPartsDistinctMultiZeroGmp));
            } case PartitionType::RepCapped: {
                return(nthPartsPtr(nthPartsRepCapGmp));
            } case PartitionType::RepNoZero: {
                return(nthPartsPtr(nthPartsRepLenGmp));
            } case PartitionType::RepShort: {
                return(nthPartsPtr(nthPartsRepShortGmp));
            } case PartitionType::RepStdAll: {
                return(nthPartsPtr(nthPartsRepGmp));
            } case PartitionType::CompRepNoZero: {
                return(nthPartsPtr(nthCompsRepGmp));
            } case PartitionType::CompRepWeak: {
                return(nthPartsPtr(nthCompsRepGmp));
            } case PartitionType::CmpRpZroNotWk: {
                return(nthPartsPtr(nthCompsRepZeroGmp));
            } case PartitionType::CmpDstctNoZero: {
                return(nthPartsPtr(nthCompsDistinctGmp));
            } case PartitionType::CmpDstctCapped: {
                return(nthPartsPtr(nthCompsDistinctGmp));
            } case PartitionType::CmpDstctWeak: {
                return(nthPartsPtr(nthCompsDistinctGmp));
            } case PartitionType::CmpDstCapWeak: {
                return(nthPartsPtr(nthCompsDistinctGmp));
            } case PartitionType::CmpDstctZNotWk: {
                return(nthPartsPtr(nthCompsDistinctMZGmp));
            } case PartitionType::CmpDstCapMZNotWk: {
                return(nthPartsPtr(nthCompsDistinctMZGmp));
            } case PartitionType::CmpDstctMZWeak: {
                return(nthPartsPtr(nthCompsDistinctMZWeakGmp));
            } case PartitionType::CmpDstCapMZWeak: {
                return(nthPartsPtr(nthCompsDistinctMZWeakGmp));
            } case PartitionType::NoSolution: {
                return(nthPartsPtr(EmptyReturn));
            } default : {
                return nullptr;
            }
        }
    } else {
        switch (ptype) {
            case PartitionType::LengthOne: {
                return(nthPartsPtr(nthLengthOne));
            } case PartitionType::DstctCapped: {
                return(nthPartsPtr(nthPartsDistinctCap));
            } case PartitionType::DstctCappedMZ: {
                return(nthPartsPtr(nthPartsDistinctCapMZ));
            } case PartitionType::DstctNoZero: {
                return(nthPartsPtr(nthPartsDistinctLen));
            } case PartitionType::DstctOneZero: {
                return(nthPartsPtr(nthPartsDistinctLen));
            } case PartitionType::DstctMultiZero: {
                return(nthPartsPtr(nthPartsDistinctMultiZero));
            } case PartitionType::DstctStdAll: {
                return(nthPartsPtr(nthPartsDistinctMultiZero));
            } case PartitionType::RepCapped: {
                return(nthPartsPtr(nthPartsRepCap));
            } case PartitionType::RepNoZero: {
                return(nthPartsPtr(nthPartsRepLen));
            } case PartitionType::RepShort: {
                return(nthPartsPtr(nthPartsRepShort));
            } case PartitionType::RepStdAll: {
                return(nthPartsPtr(nthPartsRep));
            } case PartitionType::CompRepNoZero: {
                return(nthPartsPtr(nthCompsRep));
            } case PartitionType::CompRepWeak: {
                return(nthPartsPtr(nthCompsRep));
            } case PartitionType::CmpRpZroNotWk: {
                return(nthPartsPtr(nthCompsRepZero));
            } case PartitionType::CmpDstctNoZero: {
                return(nthPartsPtr(nthCompsDistinct));
            } case PartitionType::CmpDstctCapped: {
                return(nthPartsPtr(nthCompsDistinct));
            } case PartitionType::CmpDstctWeak: {
                return(nthPartsPtr(nthCompsDistinct));
            } case PartitionType::CmpDstCapWeak: {
                return(nthPartsPtr(nthCompsDistinct));
            } case PartitionType::CmpDstctZNotWk: {
                return(nthPartsPtr(nthCompsDistinctMZ));
            } case PartitionType::CmpDstCapMZNotWk: {
                return(nthPartsPtr(nthCompsDistinctMZ));
            } case PartitionType::CmpDstctMZWeak: {
                return(nthPartsPtr(nthCompsDistinctMZWeak));
            } case PartitionType::CmpDstCapMZWeak: {
                return(nthPartsPtr(nthCompsDistinctMZWeak));
            } case PartitionType::NoSolution: {
                return(nthPartsPtr(EmptyReturn));
            } default : {
                return nullptr;
            }
        }
    }
}

nthPartsPtr GetNthPartsFuncOrStop(PartitionType ptype, bool IsGmp) {
    if (auto res = GetNthPartsFunc(ptype, IsGmp)) {
        return res;
    }

    cpp11::stop(
        "No algorithm available for PartitionType = " + GetPTypeName(ptype)
    );

    return nullptr;
}
