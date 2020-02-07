#include "GetLowerBound.h"

// This algorithm is similar to the ConstraintsGeneral algorithm. The main difference
// is in getting the next section of combinations. With the constraints algo, we rely
// on an iterative approach and with this algo we make heavy use of the respective
// "GetLowerBound" algos. Given a particular starting index vector as well as the
// specific index we are trying to set, these algorithms return the greatest lexico-
// graphical combination such that when the constraintFun (e.g. "sum") is applied,
// the value doesn't exceed the minimum target value. This algo is particularly
// effective when we need to find combinations of a vector such that the value when
// the constrainfFun is applied is between a range (e.g. comparisonFun = "==" and 
// tolerance = 100, or comparisonFun = c(">", "<") and limitConstraints = c(60, 62))
template <typename typeRcpp, typename typeVector>
typeRcpp PartitionEsqueAlgo(int n, int m, std::vector<typeVector> &v, bool isRep, 
                            std::string myFun, const std::string &comparison,
                            std::vector<typeVector> targetVals, double numRows, 
                            bool IsComb, bool xtraCol, std::vector<int> &Reps, 
                            bool IsMult, bool bUserRows) {
    
    // N.B. The range of variables below are different in ConstraintsGeneral
    // myFun is one of the following general functions: "prod", "sum", or "mean";
    // comparison is one of the comparison operator: 
    //                             "==", ">,<", ">=,<", ">,<=", ">=,<=";
    
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
    
    const Rcpp::XPtr<partialReducePtr<typeVector>> xpPartialRed = putPartialReduceInXPtr<typeVector>(myFun);
    const partialReducePtr<typeVector> partialReduce = *xpPartialRed;
        
    const Rcpp::XPtr<compPtr<typeVector>> xpCompOne = putCompPtrInXPtr<typeVector>(comparison);
    const compPtr<typeVector> compFunOne = *xpCompOne;
    
    if (IsMult) {
        for (int i = 0; i < (n-1); ++i) {
            for (int j = (i+1); j < n; ++j) {
                if (v[i] > v[j]) {
                    std::swap(v[i], v[j]);
                    std::swap(Reps[i], Reps[j]);
                }
            }
        }
    } else {
        std::sort(v.begin(), v.end());
    }
    
    const auto itComp = std::find(compSpecial.cbegin(), compSpecial.cend(), comparison);
    
    if (itComp == compSpecial.cend()) {
        // This should not happen. As we have determined that the PartType is PartitionEsque, the
        // comparison must be "==" or one of the "between" operators (i.e. those found in compSpecial)
        Rcpp::stop("We have encountered an error!!");
    }
    
    int myIndex = std::distance(compSpecial.cbegin(), itComp);
    const Rcpp::XPtr<compPtr<typeVector>> xpCompTwo = putCompPtrInXPtr<typeVector>(compHelper[myIndex]);
    const compPtr<typeVector> compFunTwo = *xpCompTwo;
    
    std::vector<int> z(m);
    std::vector<typeVector> testVec(m);
    
    bool t_0 = true;
    bool t_1 = true;
    
    int maxZ = n - 1;
    const int m1 = m - 1;
    const int m2 = m - 2;
    const typeVector currPartial = (myFun == "prod") ? 1 : 0;
    
    const typeVector targetMin = *std::min_element(targetVals.cbegin(), targetVals.cend());
    const typeVector targetMax = *std::max_element(targetVals.cbegin(), targetVals.cend());
    
    if (IsMult) {
        int freqsSize = std::accumulate(Reps.cbegin(), Reps.cend(), 0);
        std::vector<int> freqs, zIndex;
        const int pentExtreme = freqsSize - m;
        
        for (int i = 0, k = 0; i < n; ++i) {
            zIndex.push_back(k);
            
            for (int j = 0; j < Reps[i]; ++j, ++k)
                freqs.push_back(i);
        }
        
        t_1 = GetLowerBoundMulti<typeVector>(n, m, v, z, freqs, targetMin,
                                             targetMax, Reps, constraintFun,
                                             partialReduce, currPartial, partialFun);
        
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
                        
                        GetLowerBoundMulti<typeVector>(n, m, v, z, freqs, 
                                                       targetMin, targetMax, Reps,
                                                       constraintFun, partialReduce,
                                                       currPartial, partialFun, i + 1);
                        
                        for (int j = i + 1, k = zIndex[z[i]] + 1; j <= m1; ++j, ++k)
                            testVec[j] = v[freqs[k]];
                        
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
        
        t_1 = GetLowerBoundRep<typeVector>(n, m, v, z, targetMin, targetMax,
                                           constraintFun, partialReduce,
                                           currPartial, partialFun);
        
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
                        
                        GetLowerBoundRep<typeVector>(n, m, v, z, targetMin, targetMax,
                                                     constraintFun, partialReduce,
                                                     currPartial, partialFun, i + 1);
                        
                        for (int k = i + 1; k < m; ++k)
                            testVec[k] = v[z[k]];
                        
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
        t_1 = GetLowerBoundNoRep<typeVector>(n, m, v, z, targetMin, targetMax,
                                             constraintFun, partialReduce,
                                             currPartial, partialFun);
        
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
                        
                        GetLowerBoundNoRep<typeVector>(n, m, v, z, targetMin, targetMax,
                                                       constraintFun, partialReduce,
                                                       currPartial, partialFun, i + 1);

                        for (int k = (i + 1); k < m; ++k)
                            testVec[k] = v[z[k]];

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

template Rcpp::IntegerMatrix PartitionEsqueAlgo(int, int, std::vector<int>&, bool, std::string myFun,
                                                const std::string&, std::vector<int>, double,
                                                bool, bool, std::vector<int>&, bool, bool);

template Rcpp::NumericMatrix PartitionEsqueAlgo(int, int, std::vector<double>&, bool, std::string myFun,
                                                const std::string&, std::vector<double>, double,
                                                bool, bool, std::vector<int>&, bool, bool);


