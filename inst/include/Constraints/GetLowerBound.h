#ifndef GET_LOWER_BOUND_H
#define GET_LOWER_BOUND_H

#include "ConstraintsUtils.h"

template <typename T>
int GetLowerBoundNoRep(const std::vector<T> &v, std::vector<int> &z,
                       funcPtr<T> fun, partialReducePtr<T> reduce,
                       partialPtr<T> partial, T currPartial, T tarMin,
                       T tarMax, int n, int m, int strt = 0);

template <typename T>
int GetLowerBoundRep(const std::vector<T> &v, std::vector<int> &z,
                     funcPtr<T> fun, partialReducePtr<T> reduce,
                     partialPtr<T> partial, T currPartial, T tarMin,
                     T tarMax, int n, int m, int strt = 0);

template <typename T>
int GetLowerBoundMulti(const std::vector<int> &freqs,
                       const std::vector<int> &Reps,
                       const std::vector<T> &v, std::vector<int> &z,
                       funcPtr<T> fun, partialReducePtr<T> reduce,
                       partialPtr<T> partial, T currPartial, T tarMin,
                       T tarMax, int n, int m, int strt = 0);

#endif