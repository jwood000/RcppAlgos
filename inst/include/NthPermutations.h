#ifndef RcppAlgos_NthPermutation_h
#define RcppAlgos_NthPermutation_h

#include <CombPermUtility.h>

template <typename typeStd>
std::vector<typeStd> nthPermutation(int n, int r, double myIndex, 
                                    std::vector<typeStd> v, bool isRep);

template <typename typeStd>
std::vector<typeStd> nthPermutation(int n, int r, double myIndex, 
                                    std::vector<typeStd> v, bool isRep) {
    
    int j = 0, n1 = n;
    double temp, index1 = myIndex;
    std::vector<typeStd> res(r);
    
    if (isRep) {
        
        temp = std::pow((double) n, (double) r);
        
        for (int k = 0; k < r; k++) {
            temp /= n;
            j = (int) std::trunc(index1 / temp);
            res[k] = v[j];
            index1 -= (temp * (double) j);
        }
        
    } else {
        
        temp = NumPermsNoRep(n, r);
        std::vector<int> indexVec(n);
        std::iota(indexVec.begin(), indexVec.end(), 0);
        
        for (int k = 0; k < r; k++, n1--) {
            temp /= n1;
            j = (int) std::trunc(index1 / temp);
            res[k] = v[indexVec[j]];
            index1 -= (temp * (double) j);
            indexVec.erase(indexVec.begin() + j);
        }
    }
    
    return res;
}

#endif
