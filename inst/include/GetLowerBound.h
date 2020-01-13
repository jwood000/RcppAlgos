#ifndef GET_LOWER_BOUND_H
#define GET_LOWER_BOUND_H

#include "ConstraintsUtils.h"

template <typename typeVector>
int GetLowerBoundNoRep(int n, int m, const std::vector<typeVector> &v, std::vector<int> &z,
                       typeVector targetMin, typeVector targetMax, funcPtr<typeVector> constraintFun,
                       partialReducePtr<typeVector> partialReduce, typeVector,
                       partialPtr<typeVector> partialFun, int strt = 0);

template <typename typeVector>
int GetLowerBoundRep(int n, int m, const std::vector<typeVector> &v, std::vector<int> &z,
                     typeVector targetMin, typeVector targetMax, funcPtr<typeVector> constraintFun,
                     partialReducePtr<typeVector> partialReduce, typeVector,
                     partialPtr<typeVector> partialFun, int strt = 0);

template <typename typeVector>
int GetLowerBoundMulti(int n, int m, const std::vector<typeVector> &v, std::vector<int> &z,
                       const std::vector<int> &freqs, typeVector targetMin, typeVector targetMax,
                       const std::vector<int> &Reps, funcPtr<typeVector> constraintFun,
                       partialReducePtr<typeVector> partialReduce, typeVector,
                       partialPtr<typeVector> partialFun, int strt = 0);

#endif
