#ifndef PARTIAL_REDUCE_H
#define PARTIAL_REDUCE_H

#include <string>

template <typename T>
using partialReducePtr = void (*const)(int m, T &partial, T w);

template <typename T>
partialReducePtr<T> GetPartialReducePtr(const std::string &myFun);

#endif