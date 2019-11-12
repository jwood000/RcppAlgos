#ifndef USER_CONSTRAINT_FUNS_H
#define USER_CONSTRAINT_FUNS_H

#include <Rcpp.h>

template <typename stdType>
using funcPtr = stdType (*)(const std::vector<stdType> &v, int mySize);

template <typename stdType>
using partialPtr = stdType (*)(stdType partial, stdType w, int mySize);

template <typename stdType>
using compPtr = bool (*)(stdType x, const std::vector<stdType> &y);

// Below, we define five main functions that will be utilized
// as constraint functions. We also define five comparison
// operations (<, <=, >, >=, ==). The framework is based on
// the information posted by Dirk Eddelbuettel from
// this Rcpp Gallery (http://gallery.rcpp.org/articles/passing-cpp-function-pointers/)

template <typename stdType>
stdType prod(const std::vector<stdType> &v, int mySize) {
    stdType myProduct = 1;
    for (int i = 0; i < mySize; ++i) {myProduct *= v[i];}
    return (myProduct);
}

template <typename stdType>
stdType sum(const std::vector<stdType> &v, int mySize) {
    return (std::accumulate(v.cbegin(), v.cbegin() + mySize, static_cast<stdType>(0)));
}

template <typename stdType>
stdType mean(const std::vector<stdType> &v, int mySize) {
    const stdType mySum = sum(v, mySize);
    return (mySum / static_cast<double>(mySize));
}

template <typename stdType>
stdType max(const std::vector<stdType> &v, int mySize) {
    return (*std::max_element(v.cbegin(), v.cbegin() + mySize));
}

template <typename stdType>
stdType min(const std::vector<stdType> &v, int mySize) {
    return (*std::min_element(v.cbegin(), v.cbegin() + mySize));
}

// Helper functions to the above. They achieve the same result when
// itermediates are known. They are simpler and more efficienct as
// they are only performing a couple of operations as opposed to
// traversing the entire vector.
template <typename stdType>
stdType prodPartial(stdType partial, stdType w, int mySize) {
    return (partial * w);
}

template <typename stdType>
stdType sumPartial(stdType partial, stdType w, int mySize) {
    return (partial + w);
}

template <typename stdType>
stdType meanPartial(stdType partial, stdType w, int mySize) {
    return (partial + (w - partial) / static_cast<double>(mySize));
}

template <typename stdType>
stdType maxPartial(stdType partial, stdType w, int mySize) {
    return std::max(partial, w);
}

template <typename stdType>
stdType minPartial(stdType partial, stdType w, int mySize) {
    return std::min(partial, w);
}

// Standard comparison functions
template <typename stdType>
bool less(stdType x, const std::vector<stdType> &y) {return x < y[0];}

template <typename stdType>
bool greater(stdType x, const std::vector<stdType> &y) {return x > y[0];}

template <typename stdType>
bool lessEqual(stdType x, const std::vector<stdType> &y) {return x <= y[0];}

template <typename stdType>
bool greaterEqual(stdType x, const std::vector<stdType> &y) {return x >= y[0];}

template <typename stdType>
bool equalInt(stdType x, const std::vector<stdType> &y) {return x == y[0];}


// Compound comparison functions for finding values between a given range
template <typename stdType>
bool greaterLess(stdType x, const std::vector<stdType> &y) {return (x < y[0]) && (x > y[1]);}

template <typename stdType>
bool greaterEqlLess(stdType x, const std::vector<stdType> &y) {return (x < y[0]) && (x >= y[1]);}

template <typename stdType>
bool greaterLessEql(stdType x, const std::vector<stdType> &y) {return x <= y[0] && x > y[1];}

template <typename stdType>
bool greaterEqlLessEql(stdType x, const std::vector<stdType> &y) {return x <= y[0] && x >= y[1];}


template <typename stdType>
Rcpp::XPtr<funcPtr<stdType>> putFunPtrInXPtr(const std::string &fstr) {
    if (fstr == "prod")
        return(Rcpp::XPtr<funcPtr<stdType>>(new funcPtr<stdType>(&prod)));
    else if (fstr == "sum")
        return(Rcpp::XPtr<funcPtr<stdType>>(new funcPtr<stdType>(&sum)));
    else if (fstr == "mean")
        return(Rcpp::XPtr<funcPtr<stdType>>(new funcPtr<stdType>(&mean)));
    else if (fstr == "max")
        return(Rcpp::XPtr<funcPtr<stdType>>(new funcPtr<stdType>(&max)));
    else if (fstr == "min")
        return(Rcpp::XPtr<funcPtr<stdType>>(new funcPtr<stdType>(&min)));
    else
        return Rcpp::XPtr<funcPtr<stdType>>(R_NilValue); // runtime error as NULL no XPtr
}

template <typename stdType>
Rcpp::XPtr<partialPtr<stdType>> putPartialPtrInXPtr(const std::string &fstr) {
    if (fstr == "prod")
        return(Rcpp::XPtr<partialPtr<stdType>>(new partialPtr<stdType>(&prodPartial)));
    else if (fstr == "sum")
        return(Rcpp::XPtr<partialPtr<stdType>>(new partialPtr<stdType>(&sumPartial)));
    else if (fstr == "mean")
        return(Rcpp::XPtr<partialPtr<stdType>>(new partialPtr<stdType>(&meanPartial)));
    else if (fstr == "max")
        return(Rcpp::XPtr<partialPtr<stdType>>(new partialPtr<stdType>(&maxPartial)));
    else if (fstr == "min")
        return(Rcpp::XPtr<partialPtr<stdType>>(new partialPtr<stdType>(&minPartial)));
    else
        return Rcpp::XPtr<partialPtr<stdType>>(R_NilValue); // runtime error as NULL no XPtr
}

const std::vector<std::string> compVec = {"<", ">", "<=", ">=", "==",
                                          ">,<", ">=,<", ">,<=", ">=,<="};
enum myComps {
    LT = 0,
    GT = 1,
    LE = 2,
    GE = 3,
    EQ = 4,
    GTLT = 5,
    GELT = 6,
    GTLE = 7,
    GELE = 8
};

// N.B. With equality check for double data type we must call greaterEqlLessEql
// function with y being altered in the calling function to give a range (y - e, y + e)
template <typename stdType>
Rcpp::XPtr<compPtr<stdType>> putCompPtrInXPtr(std::string fstr) {
    
    auto it = std::find(compVec.cbegin(), compVec.cend(), fstr);
    const int myIndex = std::distance(compVec.cbegin(), it);
    
    switch(myIndex) {
        case LT:
            return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&less)));
        case GT:
            return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&greater)));
        case LE:
            return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&lessEqual)));
        case GE:
            return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&greaterEqual)));
        case EQ:
            if (std::is_integral<stdType>::value)
                return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&equalInt)));
            else
                return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&greaterEqlLessEql)));
        case GTLT:
            return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&greaterLess)));
        case GELT:
            return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&greaterEqlLess)));    
        case GTLE:
            return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&greaterLessEql)));
        case GELE:
            return(Rcpp::XPtr<compPtr<stdType>>(new compPtr<stdType>(&greaterEqlLessEql)));
        default:
            return Rcpp::XPtr<compPtr<stdType>>(R_NilValue); // runtime error as NULL no XPtr
    }
}

#endif
