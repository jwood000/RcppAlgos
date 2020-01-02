#ifndef CONSTRAINTS_UTILS_H
#define CONSTRAINTS_UTILS_H

#include <Rcpp.h>
#include "UserConstraintFuns.h"
#include "GmpCombPermUtils.h"

const std::vector<std::string> compForms = {"<", ">", "<=", ">=", "==", "=<", "=>"};
const std::vector<std::string> compSpecial = {"==", ">,<", ">=,<", ">,<=", ">=,<="};
const std::vector<std::string> compHelper = {"<=", "<", "<", "<=", "<="};

// in R console: print(sqrt(.Machine$double.eps), digits = 16)
// [1] 0.00000001490116119384766
// Which is also 2^(-26)
constexpr double defaultTolerance = 0.00000001490116119384766;

struct distinctType {
    int limit = 0;
    bool getAll = false;
};

distinctType DistinctAttr(int lenV, int m, bool IsRep, bool IsMult, int64_t target,
                          const std::vector<int> &Reps, bool IncludeZero);

bool CheckIsInteger(const std::string &funPass, int n, int m,
                    const std::vector<double> &vNum, const std::vector<double> &targetVals,
                    funcPtr<double> myFunDbl, bool checkLim = false);

void SetStartPartitionZ(PartitionType PartType, distinctType distinctTest,
                        std::vector<int> &z, const std::vector<int> &Reps,
                        int target, int lenV, int m, bool IncludeZero);

void ConstraintSetup(std::vector<std::string> &compFunVec,
                     std::vector<double> &targetVals, bool &IsBetweenComp);

void AdjustTargetVals(int n, VecType myType, std::vector<double> &targetVals,
                      std::vector<int> &targetIntVals, const SEXP &Rtolerance,
                      std::vector<std::string> &compFunVec, double &tolerance,
                      const std::string &mainFun, const std::vector<double> &vNum);

template <typename typeVector>
inline void BruteNextElem(int &ind, int lowBnd, typeVector targetMin,
                          typeVector partial, int m, const std::vector<typeVector> &v,
                          partialPtr<typeVector> partialFun) {
    
    typeVector dist = targetMin - partialFun(partial, v[ind], m);
    const int origInd = ind;
    
    while (ind > lowBnd && dist < 0) {
        --ind;
        dist = targetMin - partialFun(partial, v[ind], m);
    }
    
    // The main idea of this algorithm is to ensure we are finding the largest 
    // lexicographical combination that gets us as close as possible to the target
    // value without exceeding that value. As we are only implementing this in a
    // in a semi-monotic manner, we are relying on the iteratively increasing
    // algos of the main combinatoric constraint subroutine to find our actual
    // soln's. This means that if we ever have a situation where dist > 0, we must
    // increase ind in order to make dist < 0.
    //
    // The only exceptions are if ind == (lowBnd + 1) or if ind is unchanged 
    // (i.e. origInd == ind).
    //
    // The latter case is easier to explain. Basically, origInd can be viewed as
    // an upperbound and thus should the while loop never execute, we don't want
    // to increment ind to an impossible value.
    //
    // For the first case, if the loop above executed until ind == (lowBnd + 1),
    // this would mean that dist < 0 for lowBnd + 2 and dist >= 0 for lowBnd + 1.
    // If we were to increment ind to lowBnd + 2, the upperbound for the next
    // iteration will equal the lowerbound and the dist will forever be less than
    // zero moving forward. This will make finding a solution impossible. We have
    // effectively generated a combination that is greater than the combination
    // that would result in the minimum (which is not good).
    if (dist > 0 && (lowBnd + 1) != ind && origInd != ind)
        ++ind;
}

template <typename typeVector>
typeVector PartialReduce(int m, typeVector partial, 
                         typeVector w, std::string myFun) {
    
    if (myFun == "prod") {
        partial /= w;
    } else if (myFun == "sum") {
        partial -= w;
    } else if (myFun == "mean") {
        partial = (partial * static_cast<double>(m) - w) / static_cast<double>(m - 1);
    }
    
    return partial;
}

template <typename typeVector>
int GetLowerBound(int n, int m, const std::vector<typeVector> &v, bool IsRep,
                  bool IsMult, std::vector<int> &z, const std::vector<int> &freqs, 
                  const std::vector<typeVector> &targetVals, const std::vector<int> &Reps,
                  funcPtr<typeVector> constraintFun, partialPtr<typeVector> partialFun,
                  const std::string &myFun) {
    
    const int lastElem = n - 1;
    const int lastCol = m - 1;
    
    std::vector<typeVector> vPass(m);
    const typeVector targetMin = *std::min_element(targetVals.cbegin(), targetVals.cend());
    const typeVector targetMax = *std::max_element(targetVals.cbegin(), targetVals.cend());
    
    if (IsRep) {
        std::fill(vPass.begin(), vPass.end(), v.back());
    } else if (IsMult) {
        const int lenMinusM = freqs.size() - m;
        
        for (int i = freqs.size() - 1, j = 0; i >= lenMinusM; --i, ++j)
            vPass[j] = v[freqs[i]];
    } else {
        vPass.assign(v.crbegin(), v.crbegin() + m);
    }
    
    typeVector partial = constraintFun(vPass, m - 1);
    const typeVector testMax = partialFun(partial, vPass.back(), m);
    if (testMax < targetMin)  {return 0;}
    
    if (IsRep) {
        std::fill(vPass.begin(), vPass.end(), v[0]);
    } else if (IsMult) {
        for (int i = 0; i < m; ++i)
            vPass[i] = v[freqs[i]];
    } else {
        vPass.assign(v.cbegin(), v.cbegin() + m);
    }
    
    const typeVector testMin = constraintFun(vPass, m);
    if (testMin > targetMax)  {return 0;}
    
    int zExpCurrPos = IsMult ? freqs.size() - m : 0;
    int currPos = IsMult ? freqs[zExpCurrPos] : (IsRep ? lastElem: (n - m));
    
    int ind = currPos;
    int lowBnd = 0;
    std::vector<int> repsCounter;
    
    if (IsMult)
        repsCounter.assign(Reps.cbegin(), Reps.cend());
    
    for (int i = 0; i < lastCol; ++i) {
        BruteNextElem(ind, lowBnd, targetMin, partial, m, v, partialFun);
        z[i] = ind;
        partial = partialFun(partial, v[ind], m);
        
        if (IsMult) {
            --repsCounter[ind];
            
            if (repsCounter[ind] == 0)
                ++ind;
            
            ++zExpCurrPos;
            currPos = freqs[zExpCurrPos];
        } else if (!IsRep) {
            ++ind;
            ++currPos;
        }
        
        lowBnd = ind;
        ind = currPos;
        partial = PartialReduce(m, partial, v[currPos], myFun);
    }
    
    BruteNextElem(ind, lowBnd, targetMin, partial, m, v, partialFun);
    z[lastCol] = ind;
    return 1;
}

template <typename typeVector>
void GetPartitionCase(const std::vector<std::string> &compFunVec, std::vector<typeVector> &v,
                      const std::string &mainFun, typeVector target, PartitionType &PartType,
                      distinctType &distinctTest, const SEXP &Rlow, std::vector<int> &Reps,
                      int lenV, int &m, double tolerance, bool IsMult, bool IsRep) {
    
    bool PartitionCase = false;
    PartType = PartitionType::NotPartition;
    bool bLower = false;
    
    if (!Rf_isNull(Rlow)) {
        auto tempLower = FromCpp14::make_unique<mpz_t[]>(1);
        mpz_init(tempLower[0]);
        
        createMPZArray(Rlow, tempLower.get(), 1, "lower");
        bLower = mpz_cmp_si(tempLower[0], 1) > 0;
    }
    
    if (!compFunVec.empty() && !bLower) {
        if (IsMult) {
            for (int i = 0; i < (lenV - 1); ++i) {
                for (int j = (i + 1); j < lenV; ++j) {
                    if (v[i] > v[j]) {
                        std::swap(v[i], v[j]);
                        std::swap(Reps[i], Reps[j]);
                    }
                }
            }
        } else {
            std::sort(v.begin(), v.end());
        }

        int64_t tarTest = 0;

        // We don't need to check the edge case when lenV == target and there is a
        // zero. Remember, lenV is the length of the vector v, so we could have a
        // situation where v = c(0, 2, 3, 4, 5) -->> length(v) = 5. This would
        // cause a problem if we were to allow this through, however, in the
        // calling code (i.e. Combinatorics.cpp), we ensure that the distance
        // between each element is the same. This means for the example we gave,
        // we would have length(unique(diff(v))) > 1, which means PartitionCase
        // (see Combinatorics.cpp) would be false, and thus the general
        // algorithm would be executed.
        //
        // We do have to ensure that the smallest element is non-negative, othe-
        // rwise, cases like v = seq(-8, 10, 2), m = 7, rep = TRUE, & limit = 10
        // would pass as v = 0:9, m = 7, rep = TRUE, & limit = 9, --or--
        // v = 1:10, m = 7, rep = TRUE, & limit = 10

        if (compFunVec[0] == "==" && mainFun == "sum" && lenV > 1 && m > 1) {
            std::vector<typeVector> pTest(v.cbegin(), v.cend());
            std::sort(pTest.begin(), pTest.end());
            const typeVector tarDiff = pTest[1] - pTest[0];

            if (static_cast<int64_t>(pTest[0]) == pTest[0]) {
                PartitionCase = true;

                for (int i = 1; i < lenV; ++i) {
                    const typeVector testDiff = pTest[i] - pTest[i - 1];

                    if (std::abs(testDiff - tarDiff) > tolerance
                            || static_cast<int64_t>(pTest[i]) != pTest[i]) {
                        PartitionCase = false;
                        break;
                    }
                }
            }

            tarTest = static_cast<int64_t>(target);

            if (PartitionCase)
                PartitionCase = (tarTest == target);
        }
        
        if (PartitionCase) {
            const typeVector myMax = v.back();
            const bool IncludeZero = (v.front() == 0);
            PartType = PartitionType::PartGeneral;
            
            if (myMax == tarTest && (lenV == tarTest || lenV == (tarTest + 1)) && v.front() >= 0) {
                distinctTest = DistinctAttr(lenV, m, IsRep, IsMult, tarTest, Reps, IncludeZero);
                
                if (distinctTest.limit > 0) {
                    
                    m = distinctTest.limit;
                    
                    if (IncludeZero) {
                        if (IsMult) {
                            if (distinctTest.getAll) {
                                PartType = PartitionType::PartDstctStdAll;
                            } else if (Reps[0] >= (m - 1)) {
                                PartType = PartitionType::PartDstctShort;
                            } else {
                                PartType = PartitionType::PartDstctSpecial;
                            }
                        } else {
                            PartType = PartitionType::PartDstctOneZero;
                        }
                    } else {
                        PartType = PartitionType::PartDstctNoZero;
                    }
                } else if (IsRep) {
                    if (m >= lenV) { 
                        m = (IncludeZero) ? lenV - 1 : lenV;
                    }
                    
                    if (IncludeZero) {
                        PartType = PartitionType::PartTraditional;
                    } else {
                        PartType = PartitionType::PartTradNoZero;
                    }
                }
            }
        }
    }
}

template <typename typeVector>
bool CheckSpecialCase(int n, bool bLower, const std::string &mainFun,
                      const std::vector<typeVector> &vNum) {
    
    bool result = false;
    
    // If bLower, the user is looking to test a particular range. Otherwise, the constraint algo
    // will simply return (upper - lower) # of combinations/permutations that meet the criteria
    if (bLower) {
        result = true;
    } else if (mainFun == "prod") {
        for (int i = 0; i < n; ++i) {
            if (vNum[i] < 0) {
                result = true;
                break;
            }
        }
    }
    
    return result;
}

#endif
