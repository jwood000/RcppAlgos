#include "Partitions/PartitionsCountMultiset.h"
#include "Partitions/PartitionsCountDistinct.h"
#include "Partitions/BigPartsCountDistinct.h"
#include "Partitions/PartitionsCountRep.h"
#include "Partitions/BigPartsCountRep.h"
#include "Partitions/PartitionsTypes.h"
#include "Partitions/PartitionsCount.h"
#include <numeric>  // std::accumulate

#include "cpp11/R.hpp"
#include "cpp11/protect.hpp"

// The variable k is strtLen
using nthPartsPtr = std::vector<int> (*const)(int n, int m, int cap, int k,
                                      double dblIdx, const mpz_class &mpzIdx);

//*********************** Compositions Funcitons **************************//

std::vector<int> nthCompsRep(int n, int m, int cap, int k,
                             double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0) {
        for (double temp = CountCompsRepLen(n, m, cap, k);
             temp <= dblIdx; ++j) {
            --n;
            dblIdx -= temp;
            temp = CountCompsRepLen(n, m, cap, k);
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
        double temp = incr_j ? CountCompsRepLen(n, m, cap, k) :
                               CountCompsRepZero(n, m, cap, k);

        for (; temp <= dblIdx; ++j) {
            incr_j = true;
            --n;
            dblIdx -= temp;
            temp = CountCompsRepLen(n, m, cap, k);
        }

        if (incr_j) --n;
        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

//************************* Paritions Funcitons ***************************//

std::vector<int> nthPartsRepLen(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        for (double temp = CountPartsRepLen(n, m, cap, k);
             temp <= dblIdx; ++j) {
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsRepLen(n, m, cap, k);
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

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        for (double temp = CountPartsRepLenCap(n, m, cap, k);
             temp <= dblIdx; ++j) {
            n -= (m + 1);
            --cap;
            dblIdx -= temp;
            temp = CountPartsRepLenCap(n, m, cap, k);
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
        for (double temp = CountPartsDistinctLen(n, m, cap, k);
             temp <= dblIdx; ++j) {
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsDistinctLen(n, m, cap, k);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctOneZero(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    return nthPartsDistinctLen(n, m, cap, k, dblIdx, mpzIdx);
}

std::vector<int> nthPartsDistinctMultiZero(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    bool incr_j = false;
    --m;

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        double temp = (incr_j || i >= (width - k)) ?
                      CountPartsDistinctLen(n, m, cap, k) :
                      CountPartsDistinctMultiZero(n, m, cap, k);

        for (; temp <= dblIdx; ++j) {
            incr_j = true;
            n -= (m + 1);
            dblIdx -= temp;
            temp = CountPartsDistinctLen(n, m, cap, k);
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

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, --cap) {
        for (double temp = CountPartsDistinctLenCap(n, m, cap, k);
             temp <= dblIdx; ++j) {
            n -= (m + 1);
            --cap;
            dblIdx -= temp;
            temp = CountPartsDistinctLenCap(n, m, cap, k);
        }

        res[i] = j;
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

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        double temp = (incr_j || i >= (width - k)) ?
                      CountPartsDistinctLenCap(n, m, cap, k) :
                      CountPartsDistinctCapMZ(n, m, cap, k);

        for (; temp <= dblIdx; ++j) {
            incr_j = true;
            n -= (m + 1);
            --cap;
            dblIdx -= temp;
            temp = CountPartsDistinctLenCap(n, m, cap, k);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            --cap;
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

//*********************** Starting Gmp Funcitons **************************//

std::vector<int> nthCompsRepGmp(int n, int m, int cap, int k,
                                double dblIdx, const mpz_class &mpzIdx) {

    const int width = m;
    const int max_n = n;

    std::vector<int> res(width);
    --n;
    --m;

    mpz_class temp;
    mpz_class index(mpzIdx);

    const PartitionType ptype = PartitionType::RepNoZero;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype, true);

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m, j = 0) {
        myClass->GetCount(temp, n, m, cap, k);

        for (; cmp(temp, index) <= 0; ++j) {
            --n;
            index -= temp;
            myClass->GetCount(temp, n, m, cap, k);
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

    const PartitionType ptype = PartitionType::RepShort;
    std::unique_ptr<CountClass> myClass = MakeCount(ptype, true);

    for (int i = 0, j = 0; i < (width - 1); ++i, --m, j = incr_j) {
        myClass->GetCount(temp, n, m, cap, k, !incr_j);

        for (; cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            --n;
            index -= temp;
            myClass->GetCount(temp, n, m, cap, k, false);
        }

        temp = 0;
        if (incr_j) --n;
        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
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
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        myClass->GetCount(temp, n, m, cap, k);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            index -= temp;
            myClass->GetCount(temp, n, m, cap, k);
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
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --n, --m) {
        myClass->GetCount(temp, n, m, cap, k);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            --cap;
            index -= temp;
            myClass->GetCount(temp, n, m, cap, k);
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
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j) {
        myClass->GetCount(temp, n, m, cap, k);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            index -= temp;
            myClass->GetCount(temp, n, m, cap, k);
        }

        res[i] = j;
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), width);
    return res;
}

std::vector<int> nthPartsDistinctOneZeroGmp(
    int n, int m, int cap, int k, double dblIdx, const mpz_class &mpzIdx
) {

    return nthPartsDistinctLenGmp(n, m, cap, k, dblIdx, mpzIdx);
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
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        myClass->GetCount(temp, n, m, cap, k, bLiteral);

        for (; cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            n -= (m + 1);
            index -= temp;
            myClass->GetCount(temp, n, m, cap, k, false);
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
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, n -= m, --m, ++j, --cap) {
        myClass->GetCount(temp, n, m, cap, k);

        for (; cmp(temp, index) <= 0; ++j) {
            n -= (m + 1);
            --cap;
            index -= temp;
            myClass->GetCount(temp, n, m, cap, k);
        }

        res[i] = j;
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
    std::unique_ptr<CountClass> myClass = MakeCount(ptype);

    myClass->SetArrSize(ptype, n, m, cap);
    myClass->InitializeMpz();

    for (int i = 0, j = 0; i < (width - 1); ++i, --m) {
        const bool bLiteral = !(incr_j || i >= (width - k));
        myClass->GetCount(temp, n, m, cap, k, bLiteral);

        for (; cmp(temp, index) <= 0; ++j) {
            incr_j = true;
            n -= (m + 1);
            --cap;
            index -= temp;
            myClass->GetCount(temp, n, m, cap, k, false);
        }

        res[i] = j;

        if (incr_j || (i + 1) >= (width - k)) {
            ++j;
            n -= m;
            --cap;
        }
    }

    res[width - 1] = max_n - std::accumulate(res.begin(), res.end(), 0);
    return res;
}

nthPartsPtr GetNthPartsFunc(PartitionType ptype, bool IsGmp, bool IsComp) {

    if (IsComp && IsGmp) {
        switch (ptype) {
            case PartitionType::RepNoZero : {
                return(nthPartsPtr(nthCompsRepGmp));
            } case PartitionType::RepShort : {
                return(nthPartsPtr(nthCompsRepZeroGmp));
            } case PartitionType::RepStdAll : {
                return(nthPartsPtr(nthCompsRepZeroGmp));
            }default : {
                cpp11::stop("No algorithm available");
            }
        }
    } else if (IsComp) {
        switch (ptype) {
            case PartitionType::RepNoZero : {
                return(nthPartsPtr(nthCompsRep));
            } case PartitionType::RepShort : {
                return(nthPartsPtr(nthCompsRepZero));
            } case PartitionType::RepStdAll : {
                return(nthPartsPtr(nthCompsRepZero));
            }default : {
                cpp11::stop("No algorithm available");
            }
        }
    } else if (IsGmp) {
        switch (ptype) {
            case PartitionType::DstctCapped : {
                return(nthPartsPtr(nthPartsDistinctCapGmp));
            } case PartitionType::DstctCappedMZ : {
                return(nthPartsPtr(nthPartsDistinctCapMZGmp));
            } case PartitionType::DstctNoZero : {
                return(nthPartsPtr(nthPartsDistinctLenGmp));
            } case PartitionType::DstctOneZero : {
                return(nthPartsPtr(nthPartsDistinctOneZeroGmp));
            } case PartitionType::DstctMultiZero : {
                return(nthPartsPtr(nthPartsDistinctMultiZeroGmp));
            } case PartitionType::DstctStdAll : {
                return(nthPartsPtr(nthPartsDistinctMultiZeroGmp));
            } case PartitionType::RepCapped : {
                return(nthPartsPtr(nthPartsRepCapGmp));
            } case PartitionType::RepNoZero : {
                return(nthPartsPtr(nthPartsRepLenGmp));
            } case PartitionType::RepShort : {
                return(nthPartsPtr(nthPartsRepShortGmp));
            } case PartitionType::RepStdAll : {
                return(nthPartsPtr(nthPartsRepGmp));
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    } else {
        switch (ptype) {
            case PartitionType::DstctCapped : {
                return(nthPartsPtr(nthPartsDistinctCap));
            } case PartitionType::DstctCappedMZ : {
                return(nthPartsPtr(nthPartsDistinctCapMZ));
            } case PartitionType::DstctNoZero : {
                return(nthPartsPtr(nthPartsDistinctLen));
            } case PartitionType::DstctOneZero : {
                return(nthPartsPtr(nthPartsDistinctOneZero));
            } case PartitionType::DstctMultiZero : {
                return(nthPartsPtr(nthPartsDistinctMultiZero));
            } case PartitionType::DstctStdAll : {
                return(nthPartsPtr(nthPartsDistinctMultiZero));
            } case PartitionType::RepCapped : {
                return(nthPartsPtr(nthPartsRepCap));
            } case PartitionType::RepNoZero : {
                return(nthPartsPtr(nthPartsRepLen));
            } case PartitionType::RepShort : {
                return(nthPartsPtr(nthPartsRepShort));
            } case PartitionType::RepStdAll : {
                return(nthPartsPtr(nthPartsRep));
            } default : {
                cpp11::stop("No algorithm available");
            }
        }
    }
}
