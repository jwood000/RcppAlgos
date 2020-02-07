#include "ConstraintsUtils.h"

// This function applys a constraint function to a vector v with respect
// to a constraint value "target". The main idea is that combinations are added
// successively, until a particular combination exceeds the given constraint
// value for a given constraint function. After this point, we can safely skip
// several combinations knowing that they will exceed the given constraint value.
template <typename typeRcpp, typename typeVector>
typeRcpp CombinatoricsConstraints(int n, int m, std::vector<typeVector> &v, bool isRep, 
                                  std::string myFun, std::vector<std::string> comparison,
                                  std::vector<typeVector> targetVals, double numRows, 
                                  bool IsComb, bool xtraCol, std::vector<int> &Reps, 
                                  bool IsMult, bool bUserRows) {
    
    // myFun is one of the following general functions: "prod", "sum", "mean", "min", or "max";
    // The comparison vector contains up to 2 of the following comparison operator: 
    //                             "<", "<=", ">", ">=", "==", ">,<", ">=,<", ">,<=", ">=,<=";
    
    typeVector testVal;
    int count = 0;
    const int maxRows = std::min(static_cast<double>(std::numeric_limits<int>::max()), numRows);
    
    std::vector<typeVector> combinatoricsVec;
    std::vector<typeVector> resultsVec;
    
    if (bUserRows) {
        combinatoricsVec.reserve(m * maxRows);
        resultsVec.reserve(maxRows);
    }
    
    const Rcpp::XPtr<funcPtr<typeVector>> xpFun = putFunPtrInXPtr<typeVector>(myFun);
    const funcPtr<typeVector> constraintFun = *xpFun;
    
    const Rcpp::XPtr<partialPtr<typeVector>> xpPartial = putPartialPtrInXPtr<typeVector>(myFun);
    const partialPtr<typeVector> partialFun = *xpPartial;
    
    for (std::size_t nC = 0; nC < comparison.size(); ++nC) {
        
        Rcpp::XPtr<compPtr<typeVector>> xpCompOne = putCompPtrInXPtr<typeVector>(comparison[nC]);
        compPtr<typeVector> compFunOne = *xpCompOne;
        
        Rcpp::XPtr<compPtr<typeVector>> xpCompTwo = xpCompOne;
        compPtr<typeVector> compFunTwo;
        
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
            compFunTwo = *xpCompOne;
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
            
            const auto itComp = std::find(compSpecial.cbegin(), compSpecial.cend(), comparison[nC]);
            
            if (itComp != compSpecial.end()) {
                int myIndex = std::distance(compSpecial.cbegin(), itComp);
                Rcpp::XPtr<compPtr<typeVector>> xpCompThree = putCompPtrInXPtr<typeVector>(compHelper[myIndex]);
                compFunTwo = *xpCompThree;
            } else {
                compFunTwo = *xpCompOne;
            }
        }
        
        std::vector<int> z(m);
        std::vector<typeVector> testVec(m);
        
        bool t_0 = true;
        bool t_1 = true;
        
        int maxZ = n - 1;
        const int m1 = m - 1;
        const int m2 = m - 2;
        
        if (m == 1) {
            
            int ind = 0;
            testVal = v[ind];
            t_0 = compFunTwo(testVal, targetVals);
            
            while (t_0 && t_1) {
                if (compFunOne(testVal, targetVals)) {
                    for (int k = 0; k < m; ++k)
                        combinatoricsVec.push_back(v[ind]);
                    
                    ++count;
                    
                    if (xtraCol)
                        resultsVec.push_back(testVal);
                    
                    t_1 =  (count < maxRows);
                }
                
                t_0 = ind != maxZ;
                
                if (t_0) {
                    ++ind;
                    testVal = v[ind];
                    t_0 = compFunTwo(testVal, targetVals);
                }
            }
        } else {
            
            if (IsMult) {
                int freqsSize = std::accumulate(Reps.cbegin(), Reps.cend(), 0);
                std::vector<int> freqs, zIndex;
                const int pentExtreme = freqsSize - m;
                
                for (int i = 0, k = 0; i < n; ++i) {
                    zIndex.push_back(k);
                    
                    for (int j = 0; j < Reps[i]; ++j, ++k)
                        freqs.push_back(i);
                }
                
                z.assign(freqs.cbegin(), freqs.cbegin() + m);
                auto check_point_1 = std::chrono::steady_clock::now();
                
                while (t_1) {
                    SectionOne(v, testVec, z, targetVals, combinatoricsVec,
                               resultsVec, t_0, t_1, count, partialFun, constraintFun,
                               compFunOne, compFunTwo, m, m1, maxRows, maxZ, IsComb, xtraCol);
                    
                    if (t_1) {
                        bool noChange = true;
                        
                        for (int i = m2; i >= 0; --i) {
                            if (z[i] != freqs[pentExtreme + i]) {
                                ++z[i];
                                testVec[i] = v[z[i]];
                                
                                for (int j = i + 1, k = zIndex[z[i]] + 1; j <= m1; ++j, ++k) {
                                    z[j] = freqs[k];
                                    testVec[j] = v[z[j]];
                                }
                                
                                testVal = constraintFun(testVec, m);
                                t_0 = compFunTwo(testVal, targetVals);
                                noChange = false;
                                if (t_0) {break;}
                            }
                        }
                        
                        t_1 = (!noChange && t_0);
                    }
                    
                    const auto check_point_2 = std::chrono::steady_clock::now();
                    
                    if (check_point_2 - check_point_1 > timeout) {
                        Rcpp::checkUserInterrupt();
                        check_point_1 = std::chrono::steady_clock::now();
                    }
                }
                
            } else if (isRep) {
                
                v.erase(std::unique(v.begin(), v.end()), v.end());
                maxZ = static_cast<int>(v.size()) - 1;
                z.assign(m, 0);
                
                auto check_point_1 = std::chrono::steady_clock::now();
                
                while (t_1) {
                    SectionOne(v, testVec, z, targetVals, combinatoricsVec,
                               resultsVec, t_0, t_1, count, partialFun, constraintFun,
                               compFunOne, compFunTwo, m, m1, maxRows, maxZ, IsComb, xtraCol);
                    
                    if (t_1) {
                        bool noChange = true;
                        
                        for (int i = m2; i >= 0; --i) {
                            if (z[i] != maxZ) {
                                ++z[i];
                                testVec[i] = v[z[i]];
                                
                                for (int k = i + 1; k < m; ++k) {
                                    z[k] = z[k - 1];
                                    testVec[k] = v[z[k]];
                                }
                                
                                testVal = constraintFun(testVec, m);
                                t_0 = compFunTwo(testVal, targetVals);
                                noChange = false;
                                if (t_0) {break;}
                            }
                        }
                        
                        t_1 = (!noChange && t_0);
                    }
                    
                    const auto check_point_2 = std::chrono::steady_clock::now();
                    
                    if (check_point_2 - check_point_1 > timeout) {
                        Rcpp::checkUserInterrupt();
                        check_point_1 = std::chrono::steady_clock::now();
                    }
                }
                
            } else {
                
                const int nMinusM = (n - m);
                std::iota(z.begin(), z.end(), 0);
                auto check_point_1 = std::chrono::steady_clock::now();
                
                while (t_1) {
                    SectionOne(v, testVec, z, targetVals, combinatoricsVec,
                               resultsVec, t_0, t_1, count, partialFun, constraintFun,
                               compFunOne, compFunTwo, m, m1, maxRows, maxZ, IsComb, xtraCol);
                    
                    if (t_1) {
                        bool noChange = true;
                        
                        for (int i = m2; i >= 0; --i) {
                            if (z[i] != (nMinusM + i)) {
                                ++z[i];
                                testVec[i] = v[z[i]];
                                
                                for (int k = i + 1; k < m; ++k) {
                                    z[k] = z[k - 1] + 1;
                                    testVec[k] = v[z[k]];
                                }
                                
                                testVal = constraintFun(testVec, m);
                                t_0 = compFunTwo(testVal, targetVals);
                                noChange = false;
                                if (t_0) {break;}
                            }
                        }
                        
                        t_1 = (!noChange && t_0);
                    }
                    
                    const auto check_point_2 = std::chrono::steady_clock::now();
                    
                    if (check_point_2 - check_point_1 > timeout) {
                        Rcpp::checkUserInterrupt();
                        check_point_1 = std::chrono::steady_clock::now();
                    }
                }
            }
            
            targetVals.erase(targetVals.begin());
        }
    }
    const int numCols = xtraCol ? (m + 1) : m;
    typeRcpp combinatoricsMatrix = Rcpp::no_init_matrix(count, numCols);
    
    for (int i = 0, k = 0; i < count; ++i)
        for (int j = 0; j < m; ++j, ++k)
            combinatoricsMatrix(i, j) = combinatoricsVec[k];
    
    if (xtraCol)
        for (int i = 0; i < count; ++i)
            combinatoricsMatrix(i, m) = resultsVec[i];
    
    if (count > std::numeric_limits<int>::max()) {
        Rcpp::warning("The algorithm terminated early as the number of "
                          "results meeting the criteria exceeds 2^31 - 1.");
    }
    
    return combinatoricsMatrix;
}

template Rcpp::IntegerMatrix CombinatoricsConstraints(int, int, std::vector<int>&, bool, std::string myFun,
                                                      std::vector<std::string>, std::vector<int>, double,
                                                      bool, bool, std::vector<int>&, bool, bool);

template Rcpp::NumericMatrix CombinatoricsConstraints(int, int, std::vector<double>&, bool, std::string myFun,
                                                      std::vector<std::string>, std::vector<double>, double,
                                                      bool, bool, std::vector<int>&, bool, bool);


