#include "Constraints/ConstraintsUtils.h"
#include "Constraints/NextGeneralRes.h"

// This function applys a constraint function to a vector v with respect
// to a constraint value "target". The main idea is that combinations are
// added successively, until a particular combination exceeds the given
// constraint value for a given constraint function. After this point, we
// can safely skip several combinations knowing that they will exceed the
// given constraint value.

template <typename T>
void ConstraintsGeneral(std::vector<T> &v, std::vector<int> &Reps,
                        const std::vector<std::string> &comparison,
                        std::vector<T> &cnstrntVec,
                        std::vector<T> &resultsVec,
                        std::vector<T> &targetVals,
                        const std::string &myFun, double numRows,
                        int n, int m, bool IsRep, bool IsComb,
                        bool IsMult, bool bUserRows, bool xtraCol) {

    // myFun is one of the following general functions: "prod", "sum",
    // "mean", "min", or "max"; The comparison vector contains up to 2 of the
    // following comparison operator:
    //           "<", "<=", ">", ">=", "==", ">,<", ">=,<", ">,<=", ">=,<="

    T testVal;
    int count = 0;
    const int maxRows = std::min(dblIntMax, numRows);

    if (bUserRows) {
        cnstrntVec.reserve(m * maxRows);
        resultsVec.reserve(maxRows);
    }

    const funcPtr<T> fun = GetFuncPtr<T>(myFun);
    const partialPtr<T> partial = GetPartialPtr<T>(myFun);
    const nextCnstrtPtr<T> nextCnstrnt = GetCnstrtPtr<T>(IsMult, IsRep);

    for (std::size_t nC = 0; nC < comparison.size(); ++nC) {

        const compPtr<T> compOne = GetCompPtr<T>(comparison[nC]);
        compPtr<T> compTwo = compOne;

        if (comparison[nC] == ">" || comparison[nC] == ">=") {
            if (IsMult) {
                for (int i = 0; i < (n - 1); ++i) {
                    for (int j = i + 1; j < n; ++j) {
                        if (v[i] < v[j]) {
                            std::swap(v[i], v[j]);
                            std::swap(Reps[i], Reps[j]);
                        }
                    }
                }
            } else {
                std::sort(v.begin(), v.end(), std::greater<double>());
            }
        } else {
            if (IsMult) {
                for (int i = 0; i < (n - 1); ++i) {
                    for (int j = i + 1; j < n; ++j) {
                        if (v[i] > v[j]) {
                            std::swap(v[i], v[j]);
                            std::swap(Reps[i], Reps[j]);
                        }
                    }
                }
            } else {
                std::sort(v.begin(), v.end());
            }

            const auto itComp = std::find(compSpecial.cbegin(),
                                          compSpecial.cend(), comparison[nC]);

            if (itComp != compSpecial.end()) {
                int myIndex = std::distance(compSpecial.cbegin(), itComp);
                compTwo = GetCompPtr<T>(compHelper[myIndex]);
            }
        }

        std::vector<int> z(m);
        std::vector<T> testVec(m);

        bool check_0 = true;
        bool check_1 = true;

        int maxZ = n - 1;
        const int m1 = m - 1;
        const int m2 = m - 2;
        const int nMinusM = (n - m);

        if (m == 1) {
            int ind = 0;
            testVal = v[ind];
            check_0 = compTwo(testVal, targetVals);

            while (check_0 && check_1) {
                if (compOne(testVal, targetVals)) {
                    for (int k = 0; k < m; ++k) {
                        cnstrntVec.push_back(v[ind]);
                    }

                    ++count;

                    if (xtraCol) {
                        resultsVec.push_back(testVal);
                    }

                    check_1 =  (count < maxRows);
                }

                check_0 = ind != maxZ;

                if (check_0) {
                    ++ind;
                    testVal = v[ind];
                    check_0 = compTwo(testVal, targetVals);
                }
            }
        } else if (IsMult) {
            int freqsSize = std::accumulate(Reps.cbegin(), Reps.cend(), 0);
            std::vector<int> freqs;
            std::vector<int> zIndex;
            const int pentExtreme = freqsSize - m;

            for (int i = 0, k = 0; i < n; ++i) {
                zIndex.push_back(k);

                for (int j = 0; j < Reps[i]; ++j, ++k) {
                    freqs.push_back(i);
                }
            }

            z.assign(freqs.cbegin(), freqs.cbegin() + m);

            while (check_1) {
                SectionOne(v, testVec, z, targetVals, cnstrntVec,
                           resultsVec, check_0, check_1, count, partial,
                           fun, compOne, compTwo, m, m1,
                           maxRows, maxZ, IsComb, xtraCol);

                nextCnstrnt(v, targetVals, freqs, zIndex, testVec, z,
                            fun, compTwo, m, m1, m2, 0, 0, pentExtreme,
                            check_0, check_1);
            }
        } else if (IsRep) {
            std::vector<int> emptyVec;
            v.erase(std::unique(v.begin(), v.end()), v.end());
            maxZ = static_cast<int>(v.size()) - 1;
            z.assign(m, 0);

            while (check_1) {
                SectionOne(v, testVec, z, targetVals, cnstrntVec,
                           resultsVec, check_0, check_1, count, partial,
                           fun, compOne, compTwo, m, m1,
                           maxRows, maxZ, IsComb, xtraCol);

                nextCnstrnt(v, targetVals, emptyVec, emptyVec, testVec,
                            z, fun, compTwo, m, m1, m2, 0, maxZ, 0,
                            check_0, check_1);
            }
        } else {
            std::vector<int> emptyVec;
            std::iota(z.begin(), z.end(), 0);

            while (check_1) {
                SectionOne(v, testVec, z, targetVals, cnstrntVec,
                           resultsVec, check_0, check_1, count, partial,
                           fun, compOne, compTwo, m, m1,
                           maxRows, maxZ, IsComb, xtraCol);

                nextCnstrnt(v, targetVals, emptyVec, emptyVec, testVec,
                            z, fun, compTwo, m, m1, m2, nMinusM, 0, 0,
                            check_0, check_1);
            }
        }

        targetVals.erase(targetVals.begin());
    }
}

template void ConstraintsGeneral(std::vector<int>&, std::vector<int>&,
                                 const std::vector<std::string>&,
                                 std::vector<int>&, std::vector<int>&,
                                 std::vector<int>&, const std::string&,
                                 double, int, int, bool, bool,
                                 bool, bool, bool);

template void ConstraintsGeneral(std::vector<double>&, std::vector<int>&,
                                 const std::vector<std::string>&,
                                 std::vector<double>&, std::vector<double>&,
                                 std::vector<double>&, const std::string&,
                                 double, int, int, bool, bool,
                                 bool, bool, bool);
