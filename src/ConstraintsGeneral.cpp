#include "Constraints/ConstraintsUtils.h"
#include "Constraints/NextGeneralRes.h"

// This function applys a constraint function to a vector v with respect
// to a constraint value "target". The main idea is that combinations are
// added successively, until a particular combination exceeds the given
// constraint value for a given constraint function. After this point, we
// can safely skip several combinations knowing that they will exceed the
// given constraint value.

template <typename T>
void ConstraintsGeneral(std::vector<T> v, std::vector<int> Reps,
                        const std::vector<std::string> comparison,
                        std::vector<T> combinatoricsVec,
                        std::vector<T> resultsVec,
                        std::vector<T> targetVals,
                        const std::string myFun, double numRows,
                        int n, int m, bool IsRep, bool IsComb,
                        bool IsMult, bool bUserRows, bool xtraCol) {
    
    // myFun is one of the following general functions: "prod", "sum",
    // "mean", "min", or "max"; The comparison vector contains up to 2 of the
    // following comparison operator: 
    //           "<", "<=", ">", ">=", "==", ">,<", ">=,<", ">,<=", ">=,<="
    
    T testVal;
    int count = 0;
    
    constexpr double dblIntMax = std::numeric_limits<int>::max();
    const int maxRows = std::min(dblIntMax, numRows);

    if (bUserRows) {
        combinatoricsVec.reserve(m * maxRows);
        resultsVec.reserve(maxRows);
    }
    
    const funcPtr<T> constraintFun = GetFuncPtr<T>(myFun);
    const partialPtr<T> partialFun = GetPartialPtr<T>(myFun);
    
    for (std::size_t nC = 0; nC < comparison.size(); ++nC) {
        
        const compPtr<T> compFunOne = GetCompPtr<T>(comparison[nC]);
        compPtr<T> compFunTwo = compFunOne;
        
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
                compFunTwo = GetCompPtr<T>(compHelper[myIndex]);
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
            check_0 = compFunTwo(testVal, targetVals);
            
            while (check_0 && check_1) {
                if (compFunOne(testVal, targetVals)) {
                    for (int k = 0; k < m; ++k) {
                        combinatoricsVec.push_back(v[ind]);
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
                    check_0 = compFunTwo(testVal, targetVals);
                }
            }
        } else if (IsMult) {
            int freqsSize = std::accumulate(Reps.cbegin(), Reps.cend(), 0);
            std::vector<int> freqs, zIndex;
            const int pentExtreme = freqsSize - m;
            
            for (int i = 0, k = 0; i < n; ++i) {
                zIndex.push_back(k);
                
                for (int j = 0; j < Reps[i]; ++j, ++k) {
                    freqs.push_back(i);
                }
            }
            
            z.assign(freqs.cbegin(), freqs.cbegin() + m);
            auto check_poincheck_1 = std::chrono::steady_clock::now();
            
            while (check_1) {
                SectionOne(v, testVec, z, targetVals, combinatoricsVec,
                           resultsVec, check_0, check_1, count, partialFun,
                           constraintFun, compFunOne, compFunTwo, m, m1,
                           maxRows, maxZ, IsComb, xtraCol);
                
                NextCnstrntMulti(v, targetVals, freqs, zIndex, testVec,
                                 z, constraintFun, compFunTwo, m, m1, m2,
                                 pentExtreme, check_0, check_1);
                
                const auto check_point_2 = std::chrono::steady_clock::now();
                
                if (check_point_2 - check_poincheck_1 > timeout) {
                    // RcppThread::checkUserInterrupt();
                    check_poincheck_1 = std::chrono::steady_clock::now();
                }
            }
        } else if (IsRep) {
            v.erase(std::unique(v.begin(), v.end()), v.end());
            maxZ = static_cast<int>(v.size()) - 1;
            z.assign(m, 0);
            
            auto check_poincheck_1 = std::chrono::steady_clock::now();
            
            while (check_1) {
                SectionOne(v, testVec, z, targetVals, combinatoricsVec,
                           resultsVec, check_0, check_1, count, partialFun,
                           constraintFun, compFunOne, compFunTwo, m, m1,
                           maxRows, maxZ, IsComb, xtraCol);

                NextCnstrntRep(v, targetVals, testVec, z, constraintFun,
                               compFunTwo, m, m2, maxZ, check_0, check_1);
                
                const auto check_point_2 = std::chrono::steady_clock::now();
                
                if (check_point_2 - check_poincheck_1 > timeout) {
                    // RcppThread::checkUserInterrupt();
                    check_poincheck_1 = std::chrono::steady_clock::now();
                }
            }
        } else {
            std::iota(z.begin(), z.end(), 0);
            auto check_poincheck_1 = std::chrono::steady_clock::now();
            
            while (check_1) {
                SectionOne(v, testVec, z, targetVals, combinatoricsVec,
                           resultsVec, check_0, check_1, count, partialFun,
                           constraintFun, compFunOne, compFunTwo, m, m1,
                           maxRows, maxZ, IsComb, xtraCol);

                NextCnstrntDistinct(v, targetVals, testVec, z,
                                    constraintFun, compFunTwo, m,
                                    m2, nMinusM, check_0, check_1);
                
                const auto check_point_2 = std::chrono::steady_clock::now();
                
                if (check_point_2 - check_poincheck_1 > timeout) {
                    // RcppThread::checkUserInterrupt();
                    check_poincheck_1 = std::chrono::steady_clock::now();
                }
            }
        }
        
        targetVals.erase(targetVals.begin());
    }
}

template void ConstraintsGeneral(std::vector<int>, std::vector<int>, 
                                 const std::vector<std::string>,
                                 std::vector<int>, std::vector<int>,
                                 std::vector<int>, const std::string,
                                 double, int, int, bool, bool,
                                 bool, bool, bool);

template void ConstraintsGeneral(std::vector<double>, std::vector<int>, 
                                 const std::vector<std::string>,
                                 std::vector<double>, std::vector<double>,
                                 std::vector<double>, const std::string,
                                 double, int, int, bool, bool,
                                 bool, bool, bool);
