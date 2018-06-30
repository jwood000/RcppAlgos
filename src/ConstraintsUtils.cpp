#include <Rcpp.h>

typedef double (*funcPtr)(std::vector<double> &x);
typedef bool (*compPtr)(double &x, std::vector<double> &y);

// Below, we define five functions that will be utilized
// as constraint functions. We also define five comparison
// operations (<, <=, >, >=, ==). The framework is based on
// the information posted by Dirk Eddelbuettel from
// this Rcpp Gallery (http://gallery.rcpp.org/articles/passing-cpp-function-pointers/)

double prodCpp(std::vector<double>& v) {
    std::vector<double>::iterator it, vEnd = v.end();
    double myProduct = 1.0;
    for (it = v.begin(); it < vEnd; ++it) {myProduct *= *it;}
    return (myProduct);
}

double sumCpp(std::vector<double>& v) {
    std::vector<double>::iterator it, vEnd = v.end();
    double mySum = 0.0;
    for (it = v.begin(); it < vEnd; ++it) {mySum += *it;}
    return (mySum);
}

double meanCpp(std::vector<double>& v) {
    double s = v.size();
    double mySum = sumCpp(v);
    return (mySum / s);
}

double maxCpp(std::vector<double>& v) {
    return (*std::max_element(v.begin(), v.end()));
}

double minCpp(std::vector<double>& v) {
    return (*std::min_element(v.begin(), v.end()));
}

// Standard comparison functions
bool lessCpp(double &x, std::vector<double> &y) {return x < y[0];}
bool greaterCpp(double &x, std::vector<double> &y) {return x > y[0];}
bool lessEqualCpp(double &x, std::vector<double> &y) {return x <= y[0];}
bool greaterEqualCpp(double &x, std::vector<double> &y) {return x >= y[0];}
bool equalCpp(double &x, std::vector<double> &y) {return x == y[0];}

// Compound comparison functions for finding values between a given range
bool greaterLessCpp(double &x, std::vector<double> &y) {return x < y[0] && x > y[1];}
bool greaterEqlLessCpp(double &x, std::vector<double> &y) {return x < y[0] && x >= y[1];}
bool greaterLessEqlCpp(double &x, std::vector<double> &y) {return x <= y[0] && x > y[1];}
bool greaterEqlLessEqlCpp(double &x, std::vector<double> &y) {return x <= y[0] && x >= y[1];}

Rcpp::XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
    if (fstr == "prod")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&prodCpp)));
    else if (fstr == "sum")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&sumCpp)));
    else if (fstr == "mean")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&meanCpp)));
    else if (fstr == "max")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&maxCpp)));
    else if (fstr == "min")
        return(Rcpp::XPtr<funcPtr>(new funcPtr(&minCpp)));
    else
        return Rcpp::XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}

const std::vector<std::string> compVec = {"<", ">", "<=", ">=", "==",
                                          ">,<", ">=,<", ">,<=", ">=,<="};
enum myComps {
    less = 0,
    greater = 1,
    lessEql = 2,
    greaterEql = 3,
    eql = 4,
    greaterLess = 5,
    greaterEqlLess = 6,
    greaterLessEql = 7,
    greaterEqlLessEql = 8
};

Rcpp::XPtr<compPtr> putCompPtrInXPtr(std::string fstr) {
    
    std::vector<std::string>::const_iterator it = std::find(compVec.begin(), 
                                                            compVec.end(), fstr);
    int myIndex = std::distance(compVec.begin(), it);
    
    switch(myIndex) {
        case less:
            return(Rcpp::XPtr<compPtr>(new compPtr(&lessCpp)));
        case greater:
            return(Rcpp::XPtr<compPtr>(new compPtr(&greaterCpp)));
        case lessEql:
            return(Rcpp::XPtr<compPtr>(new compPtr(&lessEqualCpp)));
        case greaterEql:
            return(Rcpp::XPtr<compPtr>(new compPtr(&greaterEqualCpp)));
        case eql:
            return(Rcpp::XPtr<compPtr>(new compPtr(&equalCpp)));
        case greaterLess:
            return(Rcpp::XPtr<compPtr>(new compPtr(&greaterLessCpp)));
        case greaterEqlLess:
            return(Rcpp::XPtr<compPtr>(new compPtr(&greaterEqlLessCpp)));
        case greaterLessEql:
            return(Rcpp::XPtr<compPtr>(new compPtr(&greaterLessEqlCpp)));
        case greaterEqlLessEql:
            return(Rcpp::XPtr<compPtr>(new compPtr(&greaterEqlLessEqlCpp)));
        default:
            return Rcpp::XPtr<compPtr>(R_NilValue); // runtime error as NULL no XPtr
    }
}
