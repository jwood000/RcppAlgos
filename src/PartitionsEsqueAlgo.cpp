#include "Constraints/ConstraintsUtils.h"
#include "Constraints/GetLowerBound.h"

// This algorithm is similar to the ConstraintsGeneral algorithm. The main
// difference is in getting the next section of combinations. With the
// constraints algo, we rely on an iterative approach and with this algo we
// make heavy use of the respective "GetLowerBound" algos. Given a particular
// starting index vector as well as the specific index we are trying to set,
// these algorithms return the greatest lexicographical combination such that
// when the fun (e.g. "sum") is applied, the value doesn't exceed
// the minimum target value. This algo is particularly effective when we need
// to find combinations of a vector such that the value when the
// fun is applied is between a range (e.g. comparisonFun = "==" and
// tolerance = 100, or comparisonFun = c(">", "<") and
// limitConstraints = c(60, 62))

template <typename T>
void PartitionsEsqueAlgo(std::vector<T> &v,
                         const std::vector<T> &targetVals,
                         const std::vector<int> &Reps,
                         const std::string &myFun,
                         const std::string &comparison,
                         std::vector<T> &cnstrntVec,
                         std::vector<T> &resultsVec, int maxRows,
                         int n, int m, bool IsRep, bool IsComb,
                         bool xtraCol, bool IsMult, bool bUserRows) {

    // N.B. The range of variables below are different in ConstraintsGeneral
    // myFun is one of the following general functions: "prod", "sum", or "mean";
    // comparison is one of the comparison operator:
    //                             "==", ">,<", ">=,<", ">,<=", ">=,<=";

    T testVal;
    int count = 0;

    if (bUserRows) {
        cnstrntVec.reserve(m * maxRows);
        resultsVec.reserve(maxRows);
    }

    const funcPtr<T> fun = GetFuncPtr<T>(myFun);
    const partialPtr<T> partial = GetPartialPtr<T>(myFun);
    const partialReducePtr<T> reduce = GetPartialReducePtr<T>(myFun);
    const compPtr<T> compFunOne = GetCompPtr<T>(comparison);

    const auto itComp = std::find(compSpecial.cbegin(),
                                  compSpecial.cend(), comparison);

    if (itComp == compSpecial.cend()) {
        // This should not happen. As we have determined that the PartType is
        // PartitionEsque, the comparison must be "==" or one of the
        // "between" operators (i.e. those found in compSpecial)
        Rf_error("We have encountered an error!!");
    }

    const int myIndex = std::distance(compSpecial.cbegin(), itComp);
    const compPtr<T> compFunTwo = GetCompPtr<T>(compHelper[myIndex]);

    std::vector<int> z(m);
    std::vector<T> testVec(m);

    bool check_0 = true;
    bool check_1 = true;

    int maxZ = n - 1;
    const int m1 = m - 1;
    const int m2 = m - 2;
    const T currPartial = (myFun == "prod") ? 1 : 0;

    const T tarMin = *std::min_element(targetVals.cbegin(),
                                          targetVals.cend());
    const T tarMax = *std::max_element(targetVals.cbegin(),
                                          targetVals.cend());

    if (IsMult) {
        int freqsSize = std::accumulate(Reps.cbegin(), Reps.cend(), 0);
        std::vector<int> freqs, zIndex;
        const int pentExtreme = freqsSize - m;

        for (int i = 0, k = 0; i < n; ++i) {
            zIndex.push_back(k);

            for (int j = 0; j < Reps[i]; ++j, ++k) {
                freqs.push_back(i);
            }
        }

        check_1 = GetLowerBoundMulti(freqs, Reps, v, z, fun, reduce, partial,
                                     currPartial, tarMin, tarMax, n, m);

        while (check_1) {
            SectionOne(v, testVec, z, targetVals, cnstrntVec,
                       resultsVec, check_0, check_1, count, partial,
                       fun, compFunOne, compFunTwo, m, m1, maxRows,
                       maxZ, IsComb, xtraCol);

            if (check_1) {
                bool noChange = true;

                for (int i = m2; i >= 0; --i) {
                    if (z[i] != freqs[pentExtreme + i]) {
                        ++z[i];
                        testVec[i] = v[z[i]];

                        GetLowerBoundMulti(freqs, Reps, v, z, fun, reduce,
                                           partial, currPartial, tarMin,
                                           tarMax, n, m, i + 1);

                        for (int j = i + 1, k = zIndex[z[i]] + 1;
                             j <= m1; ++j, ++k) {
                            testVec[j] = v[freqs[k]];
                        }

                        testVal = fun(testVec, m);
                        check_0 = compFunTwo(testVal, targetVals);
                        noChange = false;

                        if (check_0) {
                            break;
                        }
                    }
                }

                check_1 = (!noChange && check_0);
            }
        }
    } else if (IsRep) {

        v.erase(std::unique(v.begin(), v.end()), v.end());
        maxZ = static_cast<int>(v.size()) - 1;

        check_1 = GetLowerBoundRep(v, z, fun, reduce, partial,
                                   currPartial, tarMin, tarMax, n, m);

        while (check_1) {
            SectionOne(v, testVec, z, targetVals, cnstrntVec,
                       resultsVec, check_0, check_1, count, partial,
                       fun, compFunOne, compFunTwo, m, m1, maxRows,
                       maxZ, IsComb, xtraCol);

            if (check_1) {
                bool noChange = true;

                for (int i = m2; i >= 0; --i) {
                    if (z[i] != maxZ) {
                        ++z[i];
                        testVec[i] = v[z[i]];

                        GetLowerBoundRep(v, z, fun, reduce, partial,
                                         currPartial, tarMin, tarMax,
                                         n, m, i + 1);

                        for (int k = i + 1; k < m; ++k) {
                            testVec[k] = v[z[k]];
                        }

                        testVal = fun(testVec, m);
                        check_0 = compFunTwo(testVal, targetVals);
                        noChange = false;

                        if (check_0) {
                            break;
                        }
                    }
                }

                check_1 = (!noChange && check_0);
            }
        }
    } else {

        const int nMinusM = (n - m);
        check_1 = GetLowerBoundNoRep(v, z, fun, reduce, partial,
                                     currPartial, tarMin, tarMax, n, m);

        while (check_1) {
            SectionOne(v, testVec, z, targetVals, cnstrntVec,
                       resultsVec, check_0, check_1, count, partial,
                       fun, compFunOne, compFunTwo, m, m1, maxRows,
                       maxZ, IsComb, xtraCol);

            if (check_1) {
                bool noChange = true;

                for (int i = m2; i >= 0; --i) {
                    if (z[i] != (nMinusM + i)) {
                        ++z[i];
                        testVec[i] = v[z[i]];

                        GetLowerBoundNoRep(v, z, fun, reduce, partial,
                                           currPartial, tarMin, tarMax,
                                           n, m, i + 1);

                        for (int k = (i + 1); k < m; ++k) {
                            testVec[k] = v[z[k]];
                        }

                        testVal = fun(testVec, m);
                        check_0 = compFunTwo(testVal, targetVals);
                        noChange = false;

                        if (check_0) {
                            break;
                        }
                    }
                }

                check_1 = (!noChange && check_0);
            }
        }
    }
}

template void PartitionsEsqueAlgo(std::vector<int>&,
                                  const std::vector<int>&,
                                  const std::vector<int>&,
                                  const std::string&, const std::string&,
                                  std::vector<int>&, std::vector<int>&,
                                  int, int, int, bool, bool, bool,
                                  bool, bool);

template void PartitionsEsqueAlgo(std::vector<double>&,
                                  const std::vector<double>&,
                                  const std::vector<int>&,
                                  const std::string&, const std::string&,
                                  std::vector<double>&, std::vector<double>&,
                                  int, int, int, bool, bool, bool,
                                  bool, bool);
