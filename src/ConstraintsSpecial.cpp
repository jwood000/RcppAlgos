#include "Constraints/UserConstraintFuns.h"
#include "Constraints/ConstraintsUtils.h"
#include "ClassUtils/NextCombinatorics.h"
#include "NthResult.h"
#include <thread>

template <typename T>
void CnstrntLowerWorker(
    const std::vector<T> &v, const std::vector<T> &targetVals,
    const std::vector<int> &freqs, const std::vector<std::string> &compVec,
    std::vector<T> &cnstrntVec, std::vector<T> &resVec, std::vector<int> &z,
    nextIterPtr nextIter, funcPtr<T> fun, compPtr<T> compOne, int m, int n1,
    int m1, int maxRows, bool xtraCol
) {

    int count = 0;
    std::vector<T> testVec(m);

    if (compVec.size() == 1) {
        do {
            for (int j = 0; j < m; ++j) {
                testVec[j] = v[z[j]];
            }

            const T testVal = fun(testVec, m);

            if (compOne(testVal, targetVals)) {
                cnstrntVec.insert(cnstrntVec.end(),
                                  testVec.begin(), testVec.end());
                if (xtraCol) resVec.push_back(testVal);
            }

            ++count;
        } while (count < maxRows && nextIter(freqs, z, n1, m1));
    } else {
        compPtr<T> compTwo = GetCompPtr<T>(compVec.back());
        std::vector<T> targetVals2(1, targetVals.back());

        do {
            for (int j = 0; j < m; ++j) {
                testVec[j] = v[z[j]];
            }

            const T testVal = fun(testVec, m);

            if (compOne(testVal, targetVals) ||
                compTwo(testVal, targetVals2)) {
                cnstrntVec.insert(cnstrntVec.end(),
                                  testVec.begin(), testVec.end());
                if (xtraCol) resVec.push_back(testVal);
            }

            ++count;
        } while (count < maxRows && nextIter(freqs, z, n1, m1));
    }
}

template <typename T>
void CnstrntSpcWorker(
    const std::vector<T> &v, const std::vector<T> &targetVals,
    const std::vector<int> &freqs, const std::vector<std::string> &compVec,
    std::vector<T> &cnstrntVec, std::vector<T> &resVec, std::vector<int> &z,
    nextIterPtr nextIter, funcPtr<T> fun, compPtr<T> compOne, int m, int n1,
    int m1, int maxRows, bool xtraCol
) {

    int count = 0;
    std::vector<T> testVec(m);

    if (compVec.size() == 1) {
        do {
            for (int j = 0; j < m; ++j) {
                testVec[j] = v[z[j]];
            }

            const T testVal = fun(testVec, m);

            if (compOne(testVal, targetVals)) {
                cnstrntVec.insert(cnstrntVec.end(),
                                  testVec.begin(), testVec.end());
                if (xtraCol) resVec.push_back(testVal);
                ++count;
            }
        } while (count < maxRows && nextIter(freqs, z, n1, m1));
    } else {
        compPtr<T> compTwo = GetCompPtr<T>(compVec.back());
        std::vector<T> targetVals2(1, targetVals.back());

        do {
            for (int j = 0; j < m; ++j) {
                testVec[j] = v[z[j]];
            }

            const T testVal = fun(testVec, m);

            if (compOne(testVal, targetVals) ||
                compTwo(testVal, targetVals2)) {
                cnstrntVec.insert(cnstrntVec.end(),
                                  testVec.begin(), testVec.end());
                if (xtraCol) resVec.push_back(testVal);
                ++count;
            }
        } while (count < maxRows && nextIter(freqs, z, n1, m1));
    }
}

// This is called when we can't easily produce a (loose) monotonic sequence overall,
// and we must generate and test every possible combination/permutation. This occurs
// when we are using "prod" and we have negative numbers involved. We also call this
// when lower is invoked implying that we are testing a specific range.
template <typename T>
void ConstraintsSpecial(
    const std::vector<T> &v, const std::vector<T> &targetVals,
    const std::vector<std::string> &compVec, const std::vector<int> &myRep,
    std::vector<int> freqs, std::vector<T> &cnstrntVec,
    std::vector<T> &resVec, const std::string &mainFun, std::vector<int> &z,
    double lower, mpz_t lowerMpz, int n, int m, int maxRows, int nThreads,
    bool IsRep, bool xtraCol, bool IsComb, bool IsMult, bool IsGmp
) {

    // Needed to determine if nextFullPerm or nextPerm will be called
    const bool IsFullPerm = (IsComb || IsRep) ? false :
                            (m == n || m == static_cast<int>(freqs.size()));

    if (!IsComb && freqs.size() == 0) {
        freqs.resize(v.size());
        std::iota(freqs.begin(), freqs.end(), 0);
    }

    const funcPtr<T> fun = GetFuncPtr<T>(mainFun);
    const compPtr<T> compOne = GetCompPtr<T>(compVec.front());
    const nextIterPtr nextIter = GetNextIterPtr(IsComb, IsMult,
                                                IsRep, IsFullPerm);

    const int n1 = IsComb ? n - 1 : (IsMult ? freqs.size() - 1 : n - 1);
    const int m1 = m - 1;

    if (nThreads > 1) {
        std::vector<std::thread> threads;

        const int stepSize = maxRows / nThreads;
        int nextStep = stepSize;
        int step = 0;

        const nthResultPtr nthResFun = GetNthResultFunc(IsComb, IsMult,
                                                        IsRep, IsGmp);
        std::vector<std::vector<int>> zs(nThreads, z);
        std::vector<std::vector<T>> resThrd(nThreads);
        std::vector<std::vector<T>> cnstrThrd(nThreads);

        for (int j = 0; j < (nThreads - 1);
             ++j, step += stepSize, nextStep += stepSize) {

            threads.emplace_back(std::cref(CnstrntLowerWorker<T>),
                                 std::cref(v), std::cref(targetVals),
                                 std::cref(freqs), std::cref(compVec),
                                 std::ref(cnstrThrd[j]),
                                 std::ref(resThrd[j]), std::ref(zs[j]),
                                 nextIter, fun, compOne, m, n1, m1,
                                 stepSize, xtraCol);

            SetNextIter(myRep, zs[j + 1], nthResFun, lower, lowerMpz,
                        stepSize, n, m, IsGmp, IsComb, IsRep, IsMult);

        }

        const int leftOver = maxRows - ((nThreads - 1) * stepSize);

        threads.emplace_back(std::cref(CnstrntLowerWorker<T>),
                             std::cref(v), std::cref(targetVals),
                             std::cref(freqs), std::cref(compVec),
                             std::ref(cnstrThrd.back()),
                             std::ref(resThrd.back()), std::ref(zs.back()),
                             nextIter, fun, compOne, m, n1, m1,
                             leftOver, xtraCol);

        for (auto& thr: threads) {
            thr.join();
        }

        for (int i = 0; i < nThreads; ++i) {
            cnstrntVec.insert(cnstrntVec.end(), cnstrThrd[i].begin(),
                              cnstrThrd[i].end());
            resVec.insert(resVec.end(), resThrd[i].begin(),
                          resThrd[i].end());
        }
    } else {
        if (lower > 0) {
            CnstrntLowerWorker(v, targetVals, freqs, compVec,
                               cnstrntVec, resVec, z, nextIter, fun,
                               compOne, m, n1, m1, maxRows, xtraCol);
        } else {
            CnstrntSpcWorker(v, targetVals, freqs, compVec,
                             cnstrntVec, resVec, z, nextIter, fun,
                             compOne, m, n1, m1, maxRows, xtraCol);
        }
    }
}

template void ConstraintsSpecial(
    const std::vector<int>&, const std::vector<int>&,
    const std::vector<std::string>&, const std::vector<int>&,
    std::vector<int>, std::vector<int>&, std::vector<int>&,
    const std::string&, std::vector<int>&, double, mpz_t,
    int, int, int, int, bool, bool, bool, bool, bool
);

template void ConstraintsSpecial(
    const std::vector<double>&, const std::vector<double>&,
    const std::vector<std::string>&, const std::vector<int>&,
    std::vector<int>, std::vector<double>&, std::vector<double>&,
    const std::string&, std::vector<int>&, double, mpz_t,
    int, int, int, int, bool, bool, bool, bool, bool
);
