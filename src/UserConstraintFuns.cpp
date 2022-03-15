#include "Constraints/UserConstraintFuns.h"

// Below, we define five main functions that will be utilized
// as constraint functions. We also define five comparison
// operations (<, <=, >, >=, ==).

template <typename T>
T prod(const std::vector<T> &v, int mySize) {
    T myProduct = 1;
    for (int i = 0; i < mySize; ++i) {myProduct *= v[i];}
    return myProduct;
}

template <typename T>
T sum(const std::vector<T> &v, int mySize) {
    return (std::accumulate(v.cbegin(), v.cbegin() + mySize, static_cast<T>(0)));
}

template <typename T>
T mean(const std::vector<T> &v, int mySize) {
    const T mySum = sum(v, mySize);
    return (mySum / static_cast<double>(mySize));
}

template <typename T>
T max(const std::vector<T> &v, int mySize) {
    return (*std::max_element(v.cbegin(), v.cbegin() + mySize));
}

template <typename T>
T min(const std::vector<T> &v, int mySize) {
    return (*std::min_element(v.cbegin(), v.cbegin() + mySize));
}

// Helper functions to the above. They achieve the same result when
// itermediates are known. They are simpler and more efficienct as
// they are only performing a couple of operations as opposed to
// traversing the entire vector.
template <typename T>
T prodPartial(T partial, T w, int mySize) {
    return (partial * w);
}

template <typename T>
T sumPartial(T partial, T w, int mySize) {
    return (partial + w);
}

template <typename T>
T meanPartial(T partial, T w, int mySize) {
    return (partial + (w - partial) / static_cast<double>(mySize));
}

template <typename T>
T maxPartial(T partial, T w, int mySize) {
    return std::max(partial, w);
}

template <typename T>
T minPartial(T partial, T w, int mySize) {
    return std::min(partial, w);
}

// Standard comparison functions
template <typename T>
bool less(T x, const std::vector<T> &y) {return x < y[0];}

template <typename T>
bool greater(T x, const std::vector<T> &y) {return x > y[0];}

template <typename T>
bool lessEqual(T x, const std::vector<T> &y) {return x <= y[0];}

template <typename T>
bool greaterEqual(T x, const std::vector<T> &y) {return x >= y[0];}

template <typename T>
bool equalInt(T x, const std::vector<T> &y) {return x == y[0];}

// Compound comparison functions for finding values between a given range
template <typename T>
bool greaterLess(T x, const std::vector<T> &y) {return (x < y[0]) && (x > y[1]);}

template <typename T>
bool greaterEqlLess(T x, const std::vector<T> &y) {return (x < y[0]) && (x >= y[1]);}

template <typename T>
bool greaterLessEql(T x, const std::vector<T> &y) {return x <= y[0] && x > y[1];}

template <typename T>
bool greaterEqlLessEql(T x, const std::vector<T> &y) {return x <= y[0] && x >= y[1];}

template <typename T>
funcPtr<T> GetFuncPtr(const std::string &fstr) {
    if (fstr == "prod") {
        return(funcPtr<T>(prod));
    } else if (fstr == "sum") {
        return(funcPtr<T>(sum));
    } else if (fstr == "mean") {
        return(funcPtr<T>(mean));
    } else if (fstr == "max") {
        return(funcPtr<T>(max));
    } else {
        return(funcPtr<T>(min));
    }
}

template <typename T>
partialPtr<T> GetPartialPtr(const std::string &fstr) {
    if (fstr == "prod") {
        return(partialPtr<T>(prodPartial));
    } else if (fstr == "sum") {
        return(partialPtr<T>(sumPartial));
    } else if (fstr == "mean") {
        return(partialPtr<T>(meanPartial));
    } else if (fstr == "max") {
        return(partialPtr<T>(maxPartial));
    } else {
        return(partialPtr<T>(minPartial));
    }
}

// N.B. With equality check for double data type we must call greaterEqlLessEql
// function with y being altered in the calling function to give a range (y - e, y + e)
template <typename T>
compPtr<T> GetCompPtr(const std::string &fstr) {

    const auto it = std::find(compVec.cbegin(), compVec.cend(), fstr);
    const int myIndex = std::distance(compVec.cbegin(), it);

    switch(myIndex) {
        case LT: {
            return(compPtr<T>(less));
        } case GT: {
            return(compPtr<T>(greater));
        } case LE: {
            return(compPtr<T>(lessEqual));
        } case GE: {
            return(compPtr<T>(greaterEqual));
        } case EQ: {
            if (std::is_integral<T>::value) {
                return(compPtr<T>(equalInt));
            } else {
                return(compPtr<T>(greaterEqlLessEql));
            }
        } case GTLT: {
            return(compPtr<T>(greaterLess));
        } case GELT: {
            return(compPtr<T>(greaterEqlLess));
        } case GTLE: {
            return(compPtr<T>(greaterLessEql));
        } default: {
            return(compPtr<T>(greaterEqlLessEql));
        }
    }
}

template <typename T>
void ReduceProd(int m, T &partial, T w) {
    partial /= w;
}

template <typename T>
void ReduceSum(int m, T &partial, T w) {
    partial -= w;
}

template <typename T>
void ReduceMean(int m, T& partial, T w) {
    partial = (partial * static_cast<double>(m) - w) / static_cast<double>(m - 1);
}

template <typename T>
reducePtr<T> GetReducePtr(const std::string &fstr) {

    if (fstr == "prod") {
        return(reducePtr<T>(ReduceProd));
    } else if (fstr == "sum") {
        return(reducePtr<T>(ReduceSum));
    } else {
        return(reducePtr<T>(ReduceMean));
    }
}

template compPtr<int> GetCompPtr(const std::string&);
template compPtr<double> GetCompPtr(const std::string&);

template partialPtr<int> GetPartialPtr(const std::string&);
template partialPtr<double> GetPartialPtr(const std::string&);

template funcPtr<int> GetFuncPtr(const std::string&);
template funcPtr<double> GetFuncPtr(const std::string&);

template reducePtr<int> GetReducePtr(const std::string&);
template reducePtr<double> GetReducePtr(const std::string&);
