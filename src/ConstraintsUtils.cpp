#include <Rcpp.h>

typedef double (*funcPtr)(std::vector<double>& x);
typedef bool (*compPtr)(double& x, double& y);

// Below, we define five functions that will be utilized
// as constraint functions. We also define five comparison
// operations (<, <=, >, >=, ==). The framework is based on
// the information posted by Dirk Eddelbuettel from
// this Rcpp Gallery (http://gallery.rcpp.org/articles/passing-cpp-function-pointers/)

double prodCpp(std::vector<double>& v) {
    std::vector<double>::iterator it, vEnd = v.end();
    double myProduct = 1.0;
    for (it = v.begin(); it < vEnd; it++) {myProduct *= *it;}
    return(myProduct);
}

double sumCpp(std::vector<double>& v) {
    std::vector<double>::iterator it, vEnd = v.end();
    double mySum = 0.0;
    for (it = v.begin(); it < vEnd; it++) {mySum += *it;}
    return(mySum);
}

double meanCpp(std::vector<double>& v){
    double s = v.size();
    double mySum = sumCpp(v);
    return (mySum/s);
}

double maxCpp(std::vector<double>& v) {
    std::vector<double>::iterator y = std::max_element(v.begin(), v.end());
    return v[std::distance(v.begin(), y)];
}

double minCpp(std::vector<double>& v) {
    std::vector<double>::iterator y = std::min_element(v.begin(), v.end());
    return v[std::distance(v.begin(), y)];
}

// Comparison functions
bool lessCpp(double& x, double& y) {return x < y;}
bool greaterCpp(double& x, double& y) {return x > y;}
bool lessEqualCpp(double& x, double& y) {return x <= y;}
bool greaterEqualCpp(double& x, double& y) {return x >= y;}
bool equalCpp(double& x, double& y) {return x == y;}

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

Rcpp::XPtr<compPtr> putCompPtrInXPtr(std::string fstr) {
    if (fstr == "<")
        return(Rcpp::XPtr<compPtr>(new compPtr(&lessCpp)));
    else if (fstr == ">")
        return(Rcpp::XPtr<compPtr>(new compPtr(&greaterCpp)));
    else if (fstr == "<=")
        return(Rcpp::XPtr<compPtr>(new compPtr(&lessEqualCpp)));
    else if (fstr == ">=")
        return(Rcpp::XPtr<compPtr>(new compPtr(&greaterEqualCpp)));
    else if (fstr == "==")
        return(Rcpp::XPtr<compPtr>(new compPtr(&equalCpp)));
    else
        return Rcpp::XPtr<compPtr>(R_NilValue); // runtime error as NULL no XPtr
}
