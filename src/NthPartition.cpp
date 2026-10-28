#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/MultisetCountEngine.h"
#include "Partitions/MultisetCountClass.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/PartitionsCount.h"
#include "Partitions/NthPartition.h"
#include <memory>
#include <numeric>
#include <vector>

// See NthPartition.h for block-scan and offset-representation conventions.

std::vector<int> nthPartsRepLen(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

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

    // See OFFSET (POSITIVE) REPRESENTATION NOTE
    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsRepLenGmp(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

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

    // See OFFSET (POSITIVE) REPRESENTATION NOTE
    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsRepShort(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> &Reps
) {

    return nthPartsRepLen(n, m, cap, k, dblIdx, mpzIdx, Reps);
}

std::vector<int> nthPartsRepShortGmp(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> &Reps
) {

    return nthPartsRepLenGmp(n, m, cap, k, dblIdx, mpzIdx, Reps);
}

std::vector<int> nthPartsRep(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> &Reps
) {

    return nthPartsRepLen(n * 2, m, cap, k, dblIdx, mpzIdx, Reps);
}

std::vector<int> nthPartsRepGmp(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> &Reps
) {

    return nthPartsRepLenGmp(n * 2, m, cap, k, dblIdx, mpzIdx, Reps);
}

std::vector<int> nthPartsRepCap(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

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

    // See OFFSET (POSITIVE) REPRESENTATION NOTE
    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsRepCapGmp(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

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

    // See OFFSET (POSITIVE) REPRESENTATION NOTE
    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctLen(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

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

    // See OFFSET (POSITIVE) REPRESENTATION NOTE
    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctLenGmp(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DistinctNoZero;
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

    // See OFFSET (POSITIVE) REPRESENTATION NOTE
    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctMZ(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
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
                      CountPartsDistinctMZ(n, m, empty_allowed, k);

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

    res.back() = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctMZGmp(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DistinctMZ;
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

    res.back() = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctCap(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

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

    // See OFFSET (POSITIVE) REPRESENTATION NOTE
    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctCapGmp(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --cap;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DistinctCapped;
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

    // See OFFSET (POSITIVE) REPRESENTATION NOTE
    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctCapMZ(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
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

    res.back() = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsDistinctCapMZGmp(
    int n, int m, int cap, int k, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> & /* Reps */
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::DistinctCappedMZ;
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

    res.back() = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

// Both multiset unrankers require positive input multiplicities. Only the
// front bucket is consumed; later buckets stay positive until reached.
// j is the original bucket offset, while allowed stays normalized to 1..k.
// Dropping the front multiplicity therefore removes allowed.back(). Advancing
// j shifts the candidate and its m remaining parts down by one, so n -= m + 1.
std::vector<int> nthPartsMultiset(
    int n, int m, int cap, int /* k */, double dblIdx,
    const mpz_class &mpzIdx, const std::vector<int> &Reps
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    std::vector<int> lclReps(Reps.begin(), Reps.end());
    std::vector<int> res(width, 0);
    --n;
    --m;

    MultisetCountEngine<false> counter;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        --lclReps.front();

        if (lclReps.front() < 0) {
            lclReps.erase(lclReps.begin());
            allowed.pop_back();
            n -= (m + 1);
            ++j;

            // Reserve the candidate from the new bucket before counting.
            --lclReps.front();
        }

        double temp = counter.Count(n, m, allowed, lclReps);

        for (; temp <= dblIdx; ++j) {
            lclReps.erase(lclReps.begin());
            allowed.pop_back();

            // Every new candidate consumes one occurrence of its own value.
            --lclReps.front();

            n -= (m + 1);
            dblIdx -= temp;
            temp = counter.Count(n, m, allowed, lclReps);
        }

        res[i] = j;
    }

    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

std::vector<int> nthPartsMultisetGmp(
    int n, int m, int cap, int /* k */, double /* dblIdx */,
    const mpz_class &mpzIdx, const std::vector<int> &Reps
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> allowed(cap);
    std::iota(allowed.begin(), allowed.end(), 1);

    std::vector<int> lclReps(Reps.begin(), Reps.end());
    std::vector<int> res(width, 0);
    --n;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);
    PartsMultiset counter(Reps);
    counter.SetArrSize(PartitionType::Multiset, n, m);
    counter.InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        --lclReps.front();

        if (lclReps.front() < 0) {
            lclReps.erase(lclReps.begin());
            allowed.pop_back();
            n -= (m + 1);
            ++j;

            // Reserve the candidate from the new bucket before counting.
            --lclReps.front();
        }

        counter.GetCount(temp, n, m, allowed, lclReps);

        for (; cmp(temp, index) <= 0; ++j) {
            lclReps.erase(lclReps.begin());
            allowed.pop_back();

            // Every new candidate consumes one occurrence of its own value.
            --lclReps.front();

            n -= (m + 1);
            index -= temp;
            counter.GetCount(temp, n, m, allowed, lclReps);
        }

        res[i] = j;
    }

    res.back() = (max_n - width) - std::accumulate(res.begin(), res.end(), 0);
    return res;
}
