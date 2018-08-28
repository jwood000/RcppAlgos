#ifndef CONSTRAINTS_UTILS_H
#define CONSTRAINTS_UTILS_H

template <typename stdType>
using funcPtr = stdType (*)(std::vector<stdType> &v, unsigned long int &mySize);

template <typename stdType>
using compPtr = bool (*)(stdType &x, std::vector<stdType> &y);

// Below, we define five functions that will be utilized
// as constraint functions. We also define five comparison
// operations (<, <=, >, >=, ==). The framework is based on
// the information posted by Dirk Eddelbuettel from
// this Rcpp Gallery (http://gallery.rcpp.org/articles/passing-cpp-function-pointers/)

template <typename stdType>
stdType prod(std::vector<stdType> &v, unsigned long int &mySize) {
    stdType myProduct = 1;
    for (std::size_t i = 0; i < mySize; ++i) {myProduct *= v[i];}
    return (myProduct);
}

template <typename stdType>
stdType sum(std::vector<stdType> &v, unsigned long int &mySize) {
    stdType mySum = 0;
    for (std::size_t i = 0; i < mySize; ++i) {mySum += v[i];}
    return (mySum);
}

template <typename stdType>
stdType mean(std::vector<stdType> &v, unsigned long int &mySize) {
    double mySum = sum(v, mySize);
    return (mySum / mySize);
}

template <typename stdType>
stdType max(std::vector<stdType> &v, unsigned long int &mySize) {
    return (*std::max_element(v.begin(), v.end()));
}

template <typename stdType>
stdType min(std::vector<stdType> &v, unsigned long int &mySize) {
    return (*std::min_element(v.begin(), v.end()));
}

// Standard comparison functions
template <typename stdType>
bool less(stdType &x, std::vector<stdType> &y) {return x < y[0];}

template <typename stdType>
bool greater(stdType &x, std::vector<stdType> &y) {return x > y[0];}

template <typename stdType>
bool lessEqual(stdType &x, std::vector<stdType> &y) {return x <= y[0];}

template <typename stdType>
bool greaterEqual(stdType &x, std::vector<stdType> &y) {return x >= y[0];}

template <typename stdType>
bool equalDbl(stdType &x, std::vector<stdType> &y) {
    return std::abs(x - y[0]) <= std::numeric_limits<double>::epsilon();
}

template <typename stdType>
bool equalInt(stdType &x, std::vector<stdType> &y) {return x == y[0];}


// Compound comparison functions for finding values between a given range
template <typename stdType>
bool greaterLess(stdType &x, std::vector<stdType> &y) {return (x < y[0]) && (x > y[1]);}

template <typename stdType>
bool greaterEqlLess(stdType &x, std::vector<stdType> &y) {return (x < y[0]) && (x >= y[1]);}

template <typename stdType>
bool greaterLessEql(stdType &x, std::vector<stdType> &y) {return x <= y[0] && x > y[1];}

template <typename stdType>
bool greaterEqlLessEql(stdType &x, std::vector<stdType> &y) {return x <= y[0] && x >= y[1];}


template <typename stdType>
Rcpp::XPtr<funcPtr<stdType> > putFunPtrInXPtr(std::string fstr) {
    if (fstr == "prod")
        return(Rcpp::XPtr<funcPtr<stdType> >(new funcPtr<stdType>(&prod)));
    else if (fstr == "sum")
        return(Rcpp::XPtr<funcPtr<stdType> >(new funcPtr<stdType>(&sum)));
    else if (fstr == "mean")
        return(Rcpp::XPtr<funcPtr<stdType> >(new funcPtr<stdType>(&mean)));
    else if (fstr == "max")
        return(Rcpp::XPtr<funcPtr<stdType> >(new funcPtr<stdType>(&max)));
    else if (fstr == "min")
        return(Rcpp::XPtr<funcPtr<stdType> >(new funcPtr<stdType>(&min)));
    else
        return Rcpp::XPtr<funcPtr<stdType> >(R_NilValue); // runtime error as NULL no XPtr
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

template <typename stdType>
Rcpp::XPtr<compPtr<stdType> > putCompPtrInXPtr(std::string fstr) {
    
    std::vector<std::string>::const_iterator it = std::find(compVec.begin(), 
                                                            compVec.end(), fstr);
    int myIndex = std::distance(compVec.begin(), it);
    
    switch(myIndex) {
        case LT:
            return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&less)));
        case GT:
            return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&greater)));
        case LE:
            return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&lessEqual)));
        case GE:
            return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&greaterEqual)));
        case EQ:
            if (std::is_integral<stdType>::value)
                return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&equalInt)));
            else
                return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&equalDbl)));
        case GTLT:
            return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&greaterLess)));
        case GELT:
            return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&greaterEqlLess)));    
        case GTLE:
            return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&greaterLessEql)));
        case GELE:
            return(Rcpp::XPtr<compPtr<stdType> >(new compPtr<stdType>(&greaterEqlLessEql)));
        default:
            return Rcpp::XPtr<compPtr<stdType> >(R_NilValue); // runtime error as NULL no XPtr
    }
}

#endif
