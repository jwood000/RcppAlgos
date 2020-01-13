#ifndef COMBINATION_RESULTS_H
#define COMBINATION_RESULTS_H

#include "UserConstraintFuns.h"

template <typename typeMatrix, typename typeVector>
void ComboGenResNoRep(typeMatrix &matRcpp, const std::vector<typeVector> &v,
                      std::vector<int> z, int n, int m, int strt, int nRows,
                      const std::vector<int> &freqs, funcPtr<typeVector> myFun);
    
template <typename typeMatrix, typename typeVector>
void ComboGenResRep(typeMatrix &matRcpp, const std::vector<typeVector> &v, 
                    std::vector<int> z, int n, int m, int strt, int nRows,
                    const std::vector<int> &freqs, funcPtr<typeVector> myFun);

template <typename typeMatrix, typename typeVector>
void MultisetComboResult(typeMatrix &matRcpp, const std::vector<typeVector> &v,
                         std::vector<int> z, int n, int m, int strt, int nRows,
                         const std::vector<int> &freqs, funcPtr<typeVector> myFun);

#endif
