#include "ConstraintsUtils.h"

// The main idea of this algorithm is to ensure we are finding the largest 
// lexicographical combination that gets us as close as possible to the target
// value without exceeding that value. We are only implementing this in a
// semi-monotic manner as we are relying on the iteratively increasing
// algos of the main combinatoric constraint subroutine to find our actual
// soln's. This means that if dist > 0, we increase ind in order to force
// dist < 0. If we did not increment ind, we would possibly miss out on the
// maximal lexicographic combination. The only time we don't want to increment
// is when the loop never gets entered (i.e. ind doesn't change) as if we were
// to increase, we would be incrementing pass the upper limit -or- when we are
// on the last index as we have no further opportunities to decrease the value.
//
// N.B. We are returning a boolean. When true is returned, this tells the
// calling function that we currently have dist > 0. In this situation we need
// to ensure that the remaining indices can be reduced enough to obtain a value
// less than the targetMin
template <typename typeVector>
inline bool BruteNextElem(int &ind, int lowBnd, typeVector targetMin,
                          typeVector partial, int m, const std::vector<typeVector> &v,
                          partialPtr<typeVector> partialFun, bool isLast = false) {
    
    typeVector dist = targetMin - partialFun(partial, v[ind], m);
    const int origInd = ind;
    
    while (ind > lowBnd && dist < 0) {
        --ind;
        dist = targetMin - partialFun(partial, v[ind], m);
    }
    
    if (dist > 0 && ind != origInd && !isLast) {
        ++ind;
        return true;
    } else {
        return false;
    }
}

template <typename typeVector>
int GetLowerBoundNoRep(int n, int m, const std::vector<typeVector> &v, std::vector<int> &z,
                       typeVector targetMin, typeVector targetMax, funcPtr<typeVector> constraintFun,
                       partialReducePtr<typeVector> partialReduce, typeVector currPartial,
                       partialPtr<typeVector> partialFun, int strt) {
    
    const int lastCol = m - 1;
    std::vector<typeVector> vPass(m);
    vPass.assign(v.crbegin(), v.crbegin() + m);
    typeVector partial = constraintFun(vPass, m - 1);
    
    if (strt == 0) {
        const typeVector testMax = partialFun(partial, vPass.back(), m);
        if (testMax < targetMin)  {return 0;}
    }
    
    int currPos = n - m;
    
    if (strt) {
        for (int i = 0; i < strt; ++i) {
            vPass[i] = v[z[i]];
            partial = partialFun(partial, vPass[i], m);
            ++currPos;
            partialReduce(m, partial, v[currPos]);
        }
        
        currPartial = constraintFun(vPass, strt);
        
        for (int i = strt, j = 1; i < m; ++i, ++j)
            vPass[i] = v[z[strt - 1] + j];
    } else {
        vPass.assign(v.cbegin(), v.cbegin() + m);
    }
    
    const typeVector testMin = constraintFun(vPass, m);
    if (testMin > targetMax)  {return 0;}
    
    int ind = n - m + strt;
    int lowBnd = (strt) ? z[strt - 1] + 1 : 0;
    
    for (int i = strt; i < lastCol; ++i) {
        if (BruteNextElem(ind, lowBnd, targetMin, partial, m, v, partialFun)) {
            if (ind > lowBnd) {
                const int numIterLeft = m - i;
                
                for (int j = 0, k = ind; j < numIterLeft; ++j, ++k)
                    vPass[j] = v[k];
                
                const typeVector minRemaining = constraintFun(vPass, numIterLeft);
                const typeVector currMin = partialFun(minRemaining, currPartial, m);
                
                if (currMin > targetMin) {
                    --ind;
                }
            }
        }
        
        z[i] = ind;
        partial = partialFun(partial, v[ind], m);
        currPartial = partialFun(currPartial, v[ind], m);
        
        ++ind;
        ++currPos;
        
        lowBnd = ind;
        ind = currPos;
        partialReduce(m, partial, v[currPos]);
    }
    
    BruteNextElem(ind, lowBnd, targetMin, partial, m, v, partialFun, true);
    z[lastCol] = ind;
    return 1;
}

template <typename typeVector>
int GetLowerBoundRep(int n, int m, const std::vector<typeVector> &v, std::vector<int> &z,
                     typeVector targetMin, typeVector targetMax, funcPtr<typeVector> constraintFun,
                     partialReducePtr<typeVector> partialReduce, typeVector currPartial,
                     partialPtr<typeVector> partialFun, int strt) {
    
    const int lastElem = n - 1;
    const int lastCol = m - 1;
    
    std::vector<typeVector> vPass(m);
    std::fill(vPass.begin(), vPass.end(), v.back());
    typeVector partial = constraintFun(vPass, m - 1);
    
    if (strt == 0) {
        const typeVector testMax = partialFun(partial, vPass.back(), m);
        if (testMax < targetMin)  {return 0;}
    }
    
    if (strt) {
        for (int i = 0; i < strt; ++i) {
            vPass[i] = v[z[i]];
            partial = partialFun(partial, vPass[i], m);
            partialReduce(m, partial, v[lastElem]);
        }
        
        currPartial = constraintFun(vPass, strt);
        
        for (int i = strt; i < m; ++i)
            vPass[i] = v[z[strt - 1]];
    } else {
        std::fill(vPass.begin(), vPass.end(), v[0]);
    }
    
    const typeVector testMin = constraintFun(vPass, m);
    if (testMin > targetMax)  {return 0;}
    
    int ind = lastElem;
    int lowBnd = (strt) ? z[strt - 1] : 0;
    
    for (int i = strt; i < lastCol; ++i) {
        if (BruteNextElem(ind, lowBnd, targetMin, partial, m, v, partialFun)) {
            if (ind > lowBnd) {
                const int numIterLeft = m - i;
                
                for (int j = 0; j < numIterLeft; ++j)
                    vPass[j] = v[ind];
                
                const typeVector minRemaining = constraintFun(vPass, numIterLeft);
                const typeVector currMin = partialFun(minRemaining, currPartial, m);
                
                if (currMin > targetMin) {
                    --ind;
                }
            }
        }
        
        z[i] = ind;
        partial = partialFun(partial, v[ind], m);
        currPartial = partialFun(currPartial, v[ind], m);
        
        lowBnd = ind;
        ind = lastElem;
        partialReduce(m, partial, v[lastElem]);
    }
    
    BruteNextElem(ind, lowBnd, targetMin, partial, m, v, partialFun, true);
    z[lastCol] = ind;
    return 1;
}

template <typename typeVector>
int GetLowerBoundMulti(int n, int m, const std::vector<typeVector> &v, std::vector<int> &z,
                       const std::vector<int> &freqs, typeVector targetMin, typeVector targetMax,
                       const std::vector<int> &Reps, funcPtr<typeVector> constraintFun,
                       partialReducePtr<typeVector> partialReduce, typeVector currPartial,
                       partialPtr<typeVector> partialFun, int strt) {
    
    const int lastCol = m - 1;
    const int lenMinusM = freqs.size() - m;
    
    std::vector<typeVector> vPass(m);
    
    for (int i = freqs.size() - 1, j = 0; i >= lenMinusM; --i, ++j)
        vPass[j] = v[freqs[i]];
    
    typeVector partial = constraintFun(vPass, m - 1);
    
    if (strt == 0) {
        const typeVector testMax = partialFun(partial, vPass.back(), m);
        if (testMax < targetMin)  {return 0;}
    }
    
    int zExpCurrPos = freqs.size() - m;
    std::vector<int> repsCounter(Reps.cbegin(), Reps.cend());
    
   if (strt) {
       for (int i = 0; i < strt; ++i) {
           vPass[i] = v[z[i]];
           partial = partialFun(partial, vPass[i], m);
           --repsCounter[z[i]];
           ++zExpCurrPos;
           partialReduce(m, partial, v[freqs[zExpCurrPos]]);
       }
       
       currPartial = constraintFun(vPass, strt);
    
       if (z[strt - 1] != freqs.back()) {
           const auto it = std::find(freqs.begin(), freqs.end(), z[strt - 1] + 1);
    
           // Find the first index in freqs that equals z[strt - 1] + 1
           // We want to get the next index after z[strt - 1], so we must
           // take into account repsCounter, which keeps track of how many
           // of each index is left.
           const int myInd = std::distance(freqs.begin(), it);
           const int freqsStrt = myInd - repsCounter[z[strt - 1]];
    
           for (int i = strt, j = freqsStrt; i < m; ++i, ++j)
               vPass[i] = v[freqs[j]];
       } else {
          for (int i = strt; i < m; ++i)
              vPass[i] = v[freqs.back()];
       }
   } else {
        for (int i = 0; i < m; ++i)
            vPass[i] = v[freqs[i]];
   }
    
    const typeVector testMin = constraintFun(vPass, m);
    if (testMin > targetMax) {return 0;}
    
    int ind = freqs[freqs.size() - m + strt];
    int lowBnd = 0;
    
    if (strt) {
        lowBnd = repsCounter[z[strt - 1]] ? z[strt - 1] : z[strt - 1] + 1;
    }
    
    for (int i = strt; i < lastCol; ++i) {
        if (BruteNextElem(ind, lowBnd, targetMin, partial, m, v, partialFun)) {
            if (ind > lowBnd && repsCounter[ind - 1]) {
                const int numIterLeft = m - i;
                const auto it = std::find(freqs.begin(), freqs.end(), ind + 1);
                const int myInd = std::distance(freqs.begin(), it);
                const int freqsStrt = myInd - repsCounter[ind];
    
                for (int j = 0, k = freqsStrt; j < numIterLeft; ++j, ++k)
                    vPass[j] = v[freqs[k]];
    
                const typeVector minRemaining = constraintFun(vPass, numIterLeft);
                const typeVector currMin = partialFun(minRemaining, currPartial, m);
    
                if (currMin > targetMin) {
                    --ind;
                }
            }
        }
        
        z[i] = ind;
        partial = partialFun(partial, v[ind], m);
        currPartial = partialFun(currPartial, v[ind], m);
        
        --repsCounter[ind];
        
        if (repsCounter[ind] == 0)
            ++ind;
        
        ++zExpCurrPos;
        lowBnd = ind;
        ind = freqs[zExpCurrPos];
        partialReduce(m, partial, v[ind]);
    }
    
    BruteNextElem(ind, lowBnd, targetMin, partial, m, v, partialFun, true);
    z[lastCol] = ind;
    return 1;
}


template int GetLowerBoundNoRep(int, int, const std::vector<int>&, std::vector<int>&,
                                int, int, funcPtr<int>, partialReducePtr<int>, int,
                                partialPtr<int>, int);
template int GetLowerBoundNoRep(int, int, const std::vector<double>&, std::vector<int>&,
                                double, double, funcPtr<double>, partialReducePtr<double>, double,
                                partialPtr<double>, int);

template int GetLowerBoundRep(int, int, const std::vector<int>&, std::vector<int>&,
                              int, int, funcPtr<int>, partialReducePtr<int>, int,
                              partialPtr<int>, int);
template int GetLowerBoundRep(int, int, const std::vector<double>&, std::vector<int>&,
                              double, double, funcPtr<double>, partialReducePtr<double>, double,
                              partialPtr<double>, int);

template int GetLowerBoundMulti(int, int, const std::vector<int>&, std::vector<int>&,
                                const std::vector<int>&, int, int,
                                const std::vector<int>&, funcPtr<int>,
                                partialReducePtr<int>, int,
                                partialPtr<int>, int);
template int GetLowerBoundMulti(int, int, const std::vector<double>&, std::vector<int>&,
                                const std::vector<int>&, double, double,
                                const std::vector<int>&, funcPtr<double>,
                                partialReducePtr<double>, double,
                                partialPtr<double>, int);

