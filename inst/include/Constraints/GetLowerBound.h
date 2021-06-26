#ifndef GET_LOWER_BOUND_H
#define GET_LOWER_BOUND_H

#include "ConstraintsUtils.h"

template <typename T>
int GetLowerBoundNoRep(int n, int m, const std::vector<T> &v,
                       std::vector<int> &z, T targetMin, T targetMax,
                       funcPtr<T> constraintFun,
                       partialReducePtr<T> partialReduce, T,
                       partialPtr<T> partialFun, int strt = 0);

template <typename T>
int GetLowerBoundRep(int n, int m, const std::vector<T> &v,
                     std::vector<int> &z, T targetMin, T targetMax,
                     funcPtr<T> constraintFun,
                     partialReducePtr<T> partialReduce, T,
                     partialPtr<T> partialFun, int strt = 0);

template <typename T>
int GetLowerBoundMulti(int n, int m, const std::vector<T> &v,
                       std::vector<int> &z, const std::vector<int> &freqs,
                       T targetMin, T targetMax,
                       const std::vector<int> &Reps,
                       funcPtr<T> constraintFun,
                       partialReducePtr<T> partialReduce, T,
                       partialPtr<T> partialFun, int strt = 0);

#endif