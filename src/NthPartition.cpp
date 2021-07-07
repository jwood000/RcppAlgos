#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsTypes.h"
#include <numeric>  // std::accumulate

#define R_NO_REMAP
#include <Rinternals.h>
#include <R.h>

// The variable k is either strtLen or cap
using nthPartsPtr = std::vector<int> (*const)(int n, int m, int k,
                                              double dblIdx, mpz_t mpzIdx);

std::vector<int> nthPartsRepLen(int n, int m, int k,
                                double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        double temp = CountPartsRepLen(n, m);

        for (; temp <= dblIdx; ++j) {
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsRepLen(n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsRepShort(int n, int m, int k,
                                  double dblIdx, mpz_t mpzIdx) {

    return nthPartsRepLen(n + m, m, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRep(int n, int m, int k,
                             double dblIdx, mpz_t mpzIdx) {

    return nthPartsRepLen(n * 2, n, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRepLenCap(int n, int m, int k,
                                   double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        double temp = CountPartsRepLenCap(n, m, k);

        for (; temp <= dblIdx; ++j) {
            n -= (m + 1);
            --k;
            dblIdx -= temp;
            temp = CountPartsRepLenCap(n, m, k);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctLen(int n, int m, int k,
                                     double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        double temp = CountPartsDistinctLen(n, m);

        for (; temp <= dblIdx; ++j) {
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsDistinctLen(n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctOneZero(int n, int m, int k,
                                         double dblIdx, mpz_t mpzIdx) {

    return nthPartsDistinctLen(n + m, m, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsDistinctMultiZero(int n, int m, int k,
                                           double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        double temp = (incr_j || i >= (width - k)) ?
                      CountPartsDistinctLen(n, m) :
                      CountPartsDistinctMultiZero(n, m, k);

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

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctLenCap(int n, int m, int k,
                                        double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --m;
    --k;

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, --k) {
        double temp = CountPartsDistinctLenCap(n, m, k);

        for (; temp <= dblIdx; ++j) {
            n -= (m + 1);
            --k;
            dblIdx -= temp;
            temp = CountPartsDistinctLenCap(n, m, k);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsRepLenGmp(int n, int m, int k,
                                   double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    mpz_t temp;
    mpz_t index;
    
    mpz_init(temp);
    mpz_init(index);
    
    mpz_set(index, mpzIdx);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        CountPartsRepLen(temp, n, m);

        for (; mpz_cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            mpz_sub(index, index, temp);
            CountPartsRepLen(temp, n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    
    mpz_clear(temp);
    mpz_clear(index);

    return res;
}

std::vector<int> nthPartsRepShortGmp(int n, int m, int k,
                                     double dblIdx, mpz_t mpzIdx) {

    return nthPartsRepLenGmp(n + m, m, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRepGmp(int n, int m, int k,
                                double dblIdx, mpz_t mpzIdx) {

    return nthPartsRepLenGmp(n * 2, n, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsRepLenCapGmp(int n, int m, int k,
                                      double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    mpz_t temp;
    mpz_t index;
    
    mpz_init(temp);
    mpz_init(index);

    mpz_set(index, mpzIdx);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        CountPartsRepLenCap(temp, n, m, k);

        for (; mpz_cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            --k;
            mpz_sub(index, index, temp);
            CountPartsRepLenCap(temp, n, m, k);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);

    mpz_clear(temp);
    mpz_clear(index);

    return res;
}

std::vector<int> nthPartsDistinctLenGmp(int n, int m, int k,
                                        double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --m;

    mpz_t temp;
    mpz_t index;
    
    mpz_init(temp);
    mpz_init(index);

    mpz_set(index, mpzIdx);

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        CountPartsDistinctLen(temp, n, m);

        for (; mpz_cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            mpz_sub(index, index, temp);
            CountPartsDistinctLen(temp, n, m);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    
    mpz_clear(temp);
    mpz_clear(index);

    return res;
}

std::vector<int> nthPartsDistinctOneZeroGmp(int n, int m, int k,
                                            double dblIdx, mpz_t mpzIdx) {

    return nthPartsDistinctLenGmp(n + m, m, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsDistinctMultiZeroGmp(int n, int m, int k,
                                              double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    mpz_t temp;
    mpz_t index;
    
    mpz_init(temp);
    mpz_init(index);

    mpz_set(index, mpzIdx);

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        if (incr_j || i >= (width - k)) {
            CountPartsDistinctLen(temp, n, m);
        } else {
            CountPartsDistinctMultiZero(temp, n, m, k);
        }

        for (; mpz_cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            n -= (m + 1);
            mpz_sub(index, index, temp);
            CountPartsDistinctLen(temp, n, m);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);

    mpz_clear(temp);
    mpz_clear(index);

    return res;
}

std::vector<int> nthPartsDistinctLenCapGmp(int n, int m, int k,
                                           double dblIdx, mpz_t mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    n -= m;
    --m;
    --k;

    mpz_t temp;
    mpz_t index;
    
    mpz_init(temp);
    mpz_init(index);

    mpz_set(index, mpzIdx);

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, --k) {
        CountPartsDistinctLenCap(temp, n, m, k);

        for (; mpz_cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            --k;
            mpz_sub(index, index, temp);
            CountPartsDistinctLenCap(temp, n, m, k);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    
    mpz_clear(temp);
    mpz_clear(index);

    return res;
}

nthPartsPtr GetNthPartsFunc(PartitionType ptype, bool IsGmp) {

    if (IsGmp) {
        switch (ptype) {
            case PartitionType::DstctCapped: {
                return(nthPartsPtr(nthPartsDistinctLenCapGmp));
            } case PartitionType::DstctNoZero : {
                return(nthPartsPtr(nthPartsDistinctLenGmp));
            } case PartitionType::DstctOneZero: {
                return(nthPartsPtr(nthPartsDistinctOneZeroGmp));
            } case PartitionType::DstctSpecial : {
                return(nthPartsPtr(nthPartsDistinctMultiZeroGmp));
            } case PartitionType::DstctStdAll: {
                return(nthPartsPtr(nthPartsDistinctMultiZeroGmp));
            } case PartitionType::RepCapped : {
                return(nthPartsPtr(nthPartsRepLenCapGmp));
            } case PartitionType::RepNoZero: {
                return(nthPartsPtr(nthPartsRepLenGmp));
            } case PartitionType::RepShort : {
                return(nthPartsPtr(nthPartsRepShortGmp));
            } case PartitionType::RepStdAll : {
                return(nthPartsPtr(nthPartsRepGmp));
            } case PartitionType::Multiset : {
                Rf_error("Investigate multiset algo later");
            } case PartitionType::CoarseGrained : {
                Rf_error("No algo available");
            } case PartitionType::NotPartition : {
                Rf_error("Error... Not partition! This should not happen");
            }
        }
    } else {
        switch (ptype) {
            case PartitionType::DstctCapped: {
                return(nthPartsPtr(nthPartsDistinctLenCap));
            } case PartitionType::DstctNoZero : {
                return(nthPartsPtr(nthPartsDistinctLen));
            } case PartitionType::DstctOneZero: {
                return(nthPartsPtr(nthPartsDistinctOneZero));
            } case PartitionType::DstctSpecial : {
                return(nthPartsPtr(nthPartsDistinctMultiZero));
            } case PartitionType::DstctStdAll: {
                return(nthPartsPtr(nthPartsDistinctMultiZero));
            } case PartitionType::RepCapped : {
                return(nthPartsPtr(nthPartsRepLenCap));
            } case PartitionType::RepNoZero: {
                return(nthPartsPtr(nthPartsRepLen));
            } case PartitionType::RepShort : {
                return(nthPartsPtr(nthPartsRepShort));
            } case PartitionType::RepStdAll : {
                return(nthPartsPtr(nthPartsRep));
            } case PartitionType::Multiset : {
                Rf_error("Investigate multiset algo later");
            } case PartitionType::CoarseGrained : {
                Rf_error("No algo available");
            } case PartitionType::NotPartition : {
                Rf_error("Error... Not partition! This should not happen");
            }
        }
    }
}
