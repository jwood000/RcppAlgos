#ifndef RcppAlgos_NthCombination_h
#define RcppAlgos_NthCombination_h

#include <CombPermUtility.h>

template <typename typeStd>
std::vector<typeStd> nthCombination(int n, int r, double myIndex, 
                                    std::vector<typeStd> v, bool isRep);

template <typename typeStd>
std::vector<typeStd> nthCombination(int n, int r, double myIndex, 
                                    std::vector<typeStd> v, bool isRep) {
    
    int j = 0, n1 = n, r1 = r - 1;
    double test, index1 = myIndex, index2 = myIndex;
    std::vector<typeStd> res(r);
    
    if (isRep) {
        for (int k = 0; k < r; k++, r1--) {
            test = NumCombsWithRep(n1, r1);
            while (test <= index1) {
                index2 -= NumCombsWithRep(n1, r1);
                n1--;
                j++;
                test += NumCombsWithRep(n1, r1);
            }
            res[k] = v[j];
            index1 = index2;
        }
    } else {
        n1--;
        
        for (int k = 0; k < r; k++, n1--, r1--, j++) {
            test = nChooseK(n1, r1);
            while (test <= index1) {
                index2 -= nChooseK(n1, r1);
                n1--;
                j++;
                test += nChooseK(n1, r1);
            }
            res[k] = v[j];
            index1 = index2;
        }
    }
    
    return res;
}

#endif
