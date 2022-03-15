#ifndef USER_CONSTRAINT_FUNS_H
#define USER_CONSTRAINT_FUNS_H

#include <vector>
#include <numeric>     // std::accumulate
#include <algorithm>   // std::min, std::min/max_element, std::find
#include <string>
#include <type_traits> // std::is_integral
#include <iterator>    // std::distance
#include <array>

enum myComps {
    LT   = 0,
    GT   = 1,
    LE   = 2,
    GE   = 3,
    EQ   = 4,
    GTLT = 5,
    GELT = 6,
    GTLE = 7,
    GELE = 8
};

static const std::array<std::string, 5> mainFunSet = {{
    "prod", "sum", "mean", "min", "max"
}};

static const std::array<std::string, 9> compVec = {{
       "<", ">",
      "<=", ">=",
      "==",
     ">,<", ">=,<",
    ">,<=", ">=,<="
}};

template <typename T>
using funcPtr = T (*)(const std::vector<T> &v, int mySize);

template <typename T>
using partialPtr = T (*)(T partial, T w, int mySize);

template <typename T>
using compPtr = bool (*)(T x, const std::vector<T> &y);

template <typename T>
using reducePtr = void (*const)(int m, T &partial, T w);

template <typename T>
funcPtr<T> GetFuncPtr(const std::string &fstr);

template <typename T>
partialPtr<T> GetPartialPtr(const std::string &fstr);

template <typename T>
compPtr<T> GetCompPtr(const std::string &fstr);

template <typename T>
reducePtr<T> GetReducePtr(const std::string &fstr);

#endif
